from __future__ import print_function

import os
import functools
import inspect
import types
import subprocess
from glob import glob
from collections import OrderedDict
import textwrap

from .utils import download_chandra, ChangeDir
from . import database


class ExternalBaseWrapper(object):
    '''Wrap a function that is run as part of a test.

    The wrapped function will have two important methods: ``__call__``
    and ``source``, which returns the source as run in a test.
    For example, the function might dynamically insert a filename
    and path into a string template to generate a set of `CIAO`_ commands.
    For display purposes, we do not want to see that template,
    but the complete commands run as part of a test,
    e.g. the CIAO commands to extract a spectrum.

    The `ExternalBaseWrapper` is the base class to define such wrappers.
    '''

    interpreter = 'none'
    '''Used to display the source of the wrapped function with a specific
    syntax highlighting.
    '''

    program = 'unknown'
    '''Used to generate header in the source listing page.'''

    def __init__(self, f):
        self.f = f
        functools.update_wrapper(self, f)

    def source(self):
        raise NotImplementedError

    @property
    def title(self):
        '''First line of function doctring

        If the function has no docstring, return the function name instead.
        '''
        if self.f.__doc__ is None:
            return self.f.__name__
        else:
            return self.f.__doc__.split('\n')[0]

    @property
    def description(self):
        '''Function docstring except for the first line

        The first line of the docstring is considered the title and is removed,
        the remaining lines are dedented.
        If the docstring does not have any lines beyond the first, return
        an empty string.
        '''
        if self.f.__doc__ is None:
            return ''
        else:
            lines = self.f.__doc__.split('\n')
            if len(lines) > 1:
                return textwrap.dedent('\n'.join(lines[1:]))
            else:
                return ''

    def __call__(self, obj, conf):
        '''
        Parameters
        ----------
        conf : `~ConfigParser.ConfigParser` instance
            The configures the path to external programs like marx.
        '''
        self.f(obj)

    def __get__(self, instance, cls=None):
        '''This implements the descriptor protocol and allows this
        callable class to be used as a bound method.
        See https://docs.python.org/2/howto/descriptor.html#functions-and-methods
        '''
        self.instance = instance
        return types.MethodType(self, instance, cls)


class Python(ExternalBaseWrapper):
    '''Wrap a python test function.

    Wrap a function written in python. ``source`` attempts to parse the python
    code and identify other python functions that are called by the wrapped
    function. These are included in the output so that the reader can see all
    important steps of the test.

    The wrapped function is just called normally. It should be defined like any
    other python function.
    '''
    interpreter = "python"
    program = 'Python'

    def source(self):
        fsource = inspect.getsourcelines(self.f)[0][2:]
        fsource = textwrap.dedent(''.join(fsource))
        # find other functions defined in same module (e.g. plotting helpers)
        for k, v in self.f.func_globals.iteritems():
            if hasattr(v, '__module__') and v.__module__ == self.f.__module__:
                # and if they are mentioned in the source code, list them, too.
                if k in fsource:
                    fsource = textwrap.dedent(''.join(inspect.getsourcelines(v)[0])) + '\n' + fsource
        preamble = '# Note that this code might not run if you directly copy and paste it:\n# - Not all import statements are shown here\n# - `self` is a reference to a test instance, which allows access to\n#   parameters such as the directory where the test is run etc.\n\n'

        return preamble + fsource

    def __call__(self, obj, conf):
        return self.f(obj)


class Sherpa(ExternalBaseWrapper):
    '''Wrap a Sherpa script.

    The wrapped function should return a single string
    (possibly with ``\\n`` in there). In order to execute that, it is written
    to a temporary file and executed in a subshell after proper CIAO
    initialization.
    '''
    interpreter = "python"
    program = 'Sherpa'

    def source(self):
        return self.f(self.instance).replace(self.instance.basepath + '/', '')

    def __call__(self, obj, conf):
        with open('sherpa_script.py', 'w') as f:
            f.write(self.f(obj))
        subprocess.call([conf.get('CIAO', 'setup') + '\nsherpa -b sherpa_script.py'],
                         cwd=obj.basepath, shell=True)


class Ciao(ExternalBaseWrapper):
    '''Wrap functions that generate `CIAO`_ shell commands.

    The output format of the wrapped function is a list of strings.

    On display with ``source()`` absolute path in the CIAO commands will be
    replaced with relative path, which are typically shorter in print.
    '''
    interpreter = "shell"
    program = 'CIAO'

    def source(self):
        return '\n'.join(self.f(self.instance)).replace(self.instance.basepath + '/', '')

    def __call__(self, obj, conf):
        commands = self.f(obj)
        commands.insert(0, conf.get(self.program, 'setup'))
        print('\n'.join(commands))
        subprocess.call(['\n'.join(commands)], shell=True, cwd=obj.basepath)


class Marx(ExternalBaseWrapper):
    '''Wrap a function that generates MARX input parameters.

    The output format of the wrapped function is either

    - a dictionary of MARX parameters.
    - or a list of such dictionaries to run a batch of marx simulations.

    All parameters that are not set in this dictionary stay
    at the default set in ``marx.par``.

    Unless there is a ``marx.par`` file in the current directory, marx is
    called with the file defined in the configuration file.
    '''
    interpreter = "shell"
    program = 'marx'

    def source(self):
        '''Assemble string for marx assuming everything is in the path.'''
        par = self.f(self.instance)
        if isinstance(par, dict):
            par = [par]

        calls = []
        for p in par:
            marxcall = ['{0}={1}'.format(k, v) for k, v in p.iteritems()]
            marxcall.insert(0, self.program)
            calls.append(' '.join(marxcall).replace(self.instance.basepath + '/', ''))
        return '\n'.join(calls)

    def __call__(self, obj, conf):
        '''Assemble string for marxcall using path information from env.'''
        par = self.f(self.instance)
        if isinstance(par, dict):
            par = [par]

        for p in par:
            marxcall = ['{0}={1}'.format(k, v) for k, v in p.iteritems()]
            marxcall.insert(0, os.path.join(conf.get('marx', 'binpath'),
                                            self.program))
            if not os.path.exists('{0}.par'.format(self.program)):
                marxcall.insert(1, '@@{0}'.format(conf.get('marx', self.program + 'parfile')))
            subprocess.call([' '.join(marxcall)], shell=True, cwd=obj.basepath)


class Marxasp(Marx):
    '''Wrap a function that generates :marxtool:`marxasp` input parameters.

    The output format of the wrapped function is either

    - a dictionary of :marxtool:`marxasp` parameters.
    - or a list of such dictionaries to run a batch of :marxtool:`marxasp` commands.

    All parameters that are not set in this dictionary stay
    at the default set in ``marxasp.par``.

    Unless there is a ``marxasp.par`` file in the current directory, marx is
    called with the file defined in the configuration file.
    '''
    program = 'marxasp'


class Marx2fits(ExternalBaseWrapper):
    '''Wrap a function that generates :marxtool:`marx2fits` input parameters.

    The output format of the wrapped function is either

    - a tuple of 3 strings, which set (in this order) the parameter (e.g.
      ``--pixadj=EDSER``), the marx directory and the output file name.
      See :marxtool:`marx2fits` for details.
    - a tuple of three arrays, where each array is a list of strings. In this case
      a batch of :marxtool:`marx2fits` commands is generated.
    '''
    interpreter = "shell"
    program = 'marx2fits'

    def args2list(self, args):
        out = []
        for a in args:
            if isinstance(a, list):
                out.append(a)
            else:
                out.append([a])
        return out

    def source(self):
        '''Assemble string for marx2fits call assuming PATH is set.'''
        options, marxdir, outfile = self.args2list(self.f(self.instance))
        source = []
        for op, md, out in zip(options, marxdir, outfile):
            source.append(' '.join([self.program, op, md, out]).replace(self.instance.basepath + '/', ''))
        return '\n'.join(source)

    def __call__(self, obj, conf):
        options, marxdir, outfile = self.args2list(self.f(self.instance))
        exefile = os.path.join(conf.get('marx', 'binpath'), self.program)
        for op, md, out in zip(options, marxdir, outfile):
            marx2fitscall = ' '.join([exefile, op, md, out])
            subprocess.call([marx2fitscall], shell=True, cwd=obj.basepath)


class SAOTraceLua(ExternalBaseWrapper):
    '''Wrap a function that writes a lua input file for `SAOTrace`_.

    The output of that function is just a string. This string is written to a
    file called ``saotrace_source.lua`` in the appropriate directory.
    '''
    interpreter = "lua"
    program = "Lua input for SAOTrace"

    def source(self):
        return self.f(self.instance)

    def __call__(self, obj, conf):
        with open('saotrace_source.lua', 'w') as f:
            for line in self.f(obj):
                f.write(line)


class SAOTrace(Ciao):
    '''Wrap a function that return shall commands to run `SAOTrace`_.

    The output format of the wrapped function is a list of strings. This
    wrapper will set up the appropriate shell envisonment to run the
    commands.
    '''
    interpreter = "shell"
    program = "SAOTrace"


class MarxTest(object):
    '''Base class to specify a |marx| test.

    Each test for |marx| is written as class that is derived from this
    class. In this way, the core functionality to execute a test can be
    provided by this base class, but it can be extended or overwritten
    with no limits by writing arbitrary python code.

    The basic requirements for a |marx| test are:

    - Set up a test by making a directory structure to save temporary and
      final products, download Chandra data from the archive that can
      be compared to simulations etc.

    - Run a one or more pieces of code to execute a full test. This
      can include, but is not limited to:

      - Python code
      - Shell code (bash by default)
      - `CIAO`_ commands (executed in a bash shell)
      - `Sherpa`_ code (executed in a subshell, because it
        requires its own setup independent from the system python)
      - `SAOTrace`_ commands (executed in a bash shell)
      - |marx| or any of its tools like :marxtool:`marx2fits`.

      Individual shell sessions need to be independent, because they require
      different and sometimes mutally exclusive setup commands.

    - While running, produce output that indicates the result of a test, e.g.
      figures or numbers.

    - Document how the test is run for accountability and to serve as
      example for users.

    This base class provides hooks to fill all these requirements, some are
    implemented here, some need to the specified in the subclasses in one of
    two ways:

    - Class attributes. Derived classes can overwrite these attributes, see
      the docstring of eahc attribute for details. For example,

    Parameters
    ----------
    conf : `~ConfigParser.ConfigParser`
        Configuration Parser with appropriate configuration file loaded
    '''

    obsid = None
    '''Set to a Chandra ObsID in a derived class'''

    download_all = False
    '''Set to ``True`` to download all files in the Chandra archive.
    If ``False`` files will be downloaded one by one when they are needed.
    '''

    title = ''
    '''String for display of these test in rendered docs.

    Should not be longer than a few words.'''

    figures = OrderedDict()
    '''Speficfy any figures that this test generates.

    The format is an `~collection.OrderedDict` where the keys are the
    file names of the figures (no file extension - the code performs some
    name mangling adding the name of the class to avoid name conflicts on
    top of the given name) and the values are again dictionaries with entries
    that set the rst figure properties (e.g. caption). An simple example is:
    ``figures = OrderedDict([('fig1', {'caption': 'Figure caption 1})])``.
    '''

    summary = ''
    '''String to display in the rendered docs to summarize the test result.'''

    version = 1
    '''If the definition of a test changes so that results cannot be
    compared to previous runs any longer, the version needs to be incremented.
    However, usually it is better to rename the test.
    '''

    expresults = []
    '''Expected test results that will be stored in the test database.

    The format is a list, where each element of the list is a dictionary
    with keys for the results database, e.g.
    ``expresults = [{'name': 't1', 'title': 'max 20 characters', 'description': 'text', 'value': 1}]``.
    See the description of the database for a complete list of fields and
    their meaning.
    '''

    def __init__(self, conf):
        self.conf = conf
        self.outpath = os.path.abspath(self.conf.get('Output', 'outpath'))
        self.basepath = os.path.abspath(os.path.join(self.conf.get('Output','temppath'),
                                                     self.name))
        if hasattr(self, 'obsid') and self.download_all:
            download_chandra(self.obsid, self.datapath)

        for t in self.expresults:
            database.insert_expected_result(self.conf, self.name, t['name'],
                                            self.version, t['value'],
                                            t['title'],
                                            t.get('description', ''),
                                            t.get('unit', ''),
                                            t.get('acceptable', None))

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def doc(self):
        '''First line of class doc has no indentation. Bring all lines
        to the same level.
        '''
        t = self.__doc__.split('\n')
        return t[0] + '\n' + textwrap.dedent('\n'.join(t[1:]))

    @property
    def datapath(self):
        '''Path to downloaded Chandra data.'''
        return os.path.join(self.basepath, 'download')

    def figpath(self, name):
        '''Return path and filename for a figure with name.

        If the target directory does not exist, it will be created.

        Parameters
        ----------
        name : string
            name of figure (without extension)

        Returns
        -------
        pathname : string
            File path and complete filename. The filename is composed of the
            name of the test, the input ``name`` and the file format defined in
            the configuration file.
        '''
        fileformat = self.conf.get('Output', 'plotformat')
        if not os.path.exists(os.path.join(self.outpath, 'figures')):
            os.makedirs(os.path.join(self.outpath, 'figures'))

        return os.path.join(self.outpath, 'figures',
                            '{0}_{1}.{2}'.format(self.name, name,
                                                 fileformat))

    @property
    def steplist(self):
        '''Get a sorted list of methods called ``step_*``.
        '''
        steplist = [method for method in dir(self) if (method[:5] == 'step_')]
        steplist.sort(key=lambda s: int(s.split('_')[1]))
        return [getattr(self, s) for s in steplist]

    def run(self):
        '''Setup test, make output directories and run all steps.'''
        self.pkg_data = os.path.join(self.conf.get('tests', 'path'),
                                     self.conf.get('tests', 'modulename'),
                                     'data')
        if not os.path.exists(self.outpath):
                os.makedirs(self.outpath)

        with ChangeDir(self.basepath):
            for step in self.steplist:
                step(self.conf)

    def get_data_file(self, filetype, download=True):
        '''Retrieve filename for a Chandra data file.

        Parameters
        ----------
        filetype : string
            A string that includes part of the filename, e.g. ``evt2``.
        download : bool
            If ``True`` a file that is not found on disk, will be downloaded
            from the Chandra ftp archive, is possible.

        Returns
        -------
        filename : string or None
            Return filename and full path in the local directory structure.
            If a file is not found locally and download if ``False`` (or the
            file does not exist in the ftp archive either), the filename is
            set to ``None``. Check the spelling of ``filetype``!
        '''
        filename = glob(os.path.join(self.datapath,
                                     '*',  # primary or secondary
                                     '*{0}*'.format(filetype)))
        if len(filename) == 0:
            if download:
                download_chandra(self.obsid, self.datapath, [filetype])
                return self.get_data_file(filetype, download=False)
            else:
                return None

        elif len(filename) == 1:
            return filename[0]

        else:
            raise ValueError('Specification not unique. Found {0}'.format(filename))

    def save_test_result(self, testname, value, sigma=None, sigma2=None):
        '''Save a test result to the test database.

        Parameters
        ----------
        testname : string
            Name of this test. Needs to match on of the names in
            expresults.
        value : typically a float
            results of the test
        sigma : typically a float
            Uncertainty on ``value`` e.g. for fit results. Can be ``NONE``
            if there is no uncertainty associated with the test result.
        sigma2 : typically a float
            In some cases a single value for ``sigma`` might not be sufficient,
            e.g. if error up and down a different. In this case, ``sigma2`` hold the
            second number.
        '''
        database.insert_test_run(self.conf, self.name, testname, self.version,
                                 value, sigma, sigma2)
