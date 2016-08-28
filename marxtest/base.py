import os
import functools
import inspect
import types
import subprocess
from glob import glob
from collections import OrderedDict
import textwrap

from .utils import download_chandra, ChangeDir
from .run_external import external_settings


class OutputNumber(object):
    '''
    '''
    max_title_length = 20

    def __init__(self, testname, value, title, unit='',
                 error=None, expected=None, description=None):
        self.testname = testname
        self.value = value
        self._title = title
        self.error = error
        self.expected = expected
        self.unit = unit
        self._description = description

    @property
    def title(self):
        return self._title

    @title.setter
    def title(self, title):
        if len(title) > self.max_title_length:
            raise ValueError('Title can only be {0} characters long'.format(self.max_title_length))
        else:
            self._title = title

    @property
    def description(self):
        if self._description is not None:
            return self._description
        else:
            return self.title

    @description.setter
    def description(self, description):
        self._description = description

    def to_dict(self):
        return {'testname': self.testname,
                'value': self.value,
                'title': self.title,
                'error': self.error,
                'expected': self.expected,
                'description': self.description,
                'unit': self.unit
                }

def run_external(cmdlist, setup=None, **kwargs):
    '''
    Parameters
    ----------
    cmdlist : list
        List of shell commands to be executed
    setup : string or None
        If this is a string, then ``cmdlist`` is prefixed with the
        corresponding entry in ``external_setup``.
    '''
    subprocess.call([external_settings[setup] + '\n' + '\n'.join(cmdlist)],
                    shell=True, **kwargs)


class ExternalBaseWrapper(object):
    '''Wrap a function that is run as part of a test.

    The wrapped function will have two important methods: ``__call__``
    and ``source``, which returns the source as run in a test.
    In the simplest case, the function might dynamically insert a filename
    and path into a string template. For display purposes, we do not want
    to see that template, but the complete commands run as part of a test,
    e.g. the CIAO commands to extract a spectrum.
    '''
    interpreter = 'undefined'
    program = 'unknown'

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

    def __call__(self, obj):
        self.f(obj)

    def __get__(self, instance, cls=None):
        ''' This implements the descriptor protocol and allows this
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
        return fsource

    def __call__(self, obj):
        return self.f(obj)


class Sherpa(ExternalBaseWrapper):
    '''Wrap a Sherpa script.

    The wrapped function should return a single string
    (possibly with \n in there). In order to execute that, it is written
    to a temporary file and executed in a subshell after proper CIAO
    initialization.
    '''
    interpreter = "python"
    program = 'Sherpa'

    def source(self):
        return self.f(self.instance).replace(self.instance.basepath + '/', '')

    def __call__(self, obj):
        with open('sherpa_script.py', 'w') as f:
            f.write(self.f(obj))
        run_external(['sherpa -b sherpa_script.py'], setup='CIAO', cwd=obj.basepath)


class Ciao(ExternalBaseWrapper):
    '''Wrap functions that generate CIAO shell commands.

    The output format of the wrapped function is a list of strings.

    On display with ``source`` absolute path in the CIAO commands will be
    replaced with relative path, which are typically shorter in print.
    '''
    interpreter = "shell"
    program = 'CIAO'

    def source(self):
        return '\n'.join(self.f(self.instance)).replace(self.instance.basepath + '/', '')

    def __call__(self, obj):
        run_external(self.f(obj), setup=self.program, cwd=obj.basepath)


class Marx(ExternalBaseWrapper):
    '''Wrap a function that generates MARX input parameters.

    The output format of the wrapped function is a dictionary of MARX
    parameters. All parameters that are not set in this dictionary stay
    at the default set in marx.par in the marx installation.
    '''
    interpreter = "shell"
    program = 'marx'

    def source(self):
        '''Assemble string for marx assuming everything is in the path.'''
        par = self.f(self.instance)
        marxcall = ['{0}={1}'.format(k, v) for k, v in par.iteritems()]
        marxcall.insert(0, self.program)

        return ' '.join(marxcall).replace(self.instance.basepath + '/', '')

    def __call__(self, obj):
        '''Assemble string for marxcall using path information from env.'''
        par = self.f(self.instance)
        marxcall = ['{0}={1}'.format(k, v) for k, v in par.iteritems()]
        marxcall.insert(0, os.path.join(external_settings['marxpath'],
                                        self.program))
        marxcall.insert(1, '@@{0}'.format(external_settings[self.program + 'parfile']))

        run_external([' '.join(marxcall)], cwd=obj.basepath)


class Marxasp(Marx):
    program = 'marxasp'


class Marx2fits(ExternalBaseWrapper):
    interpreter = "shell"
    program = 'marx2fits'

    def source(self):
        '''Assemble string for marx2fits call assuming PATH is set.'''
        options, marxdir, outfile = self.f(self.instance)
        return ' '.join([self.program, options, marxdir, outfile]).replace(self.instance.basepath + '/', '')

    def __call__(self, obj):
        options, marxdir, outfile = self.f(self.instance)
        marx2fitscall = ' '.join([os.path.join(external_settings['marxpath'], self.program),
                                  options, marxdir, outfile])
        run_external([marx2fitscall], cwd=obj.basepath)


class SAOTraceLua(ExternalBaseWrapper):
    interpreter = "lua"
    program = "Lua input for SAOTrace"

    def source(self):
        return self.f(self.instance)

    def __call__(self, obj):
        with open('saotrace_source.lua', 'w') as f:
            for line in self.f(obj):
                f.write(line)


class SAOTrace(Ciao):
    interpreter = "shell"
    program = "SAOTrace"


class MarxTest(object):
    '''

    Parameters
    ----------
    env : dict
        This dictionary holds the configuration of the test environment.
        The following keys are required:

        - temppath: path to a directory used for temporary files (e.g. data download,
          CIAO output etc.). The tests will create appropriate sub-directories.
          To aid debugging, this directory is **not** automatically deleted
          when the tests are done. However, removing it will not impact the
          webpages that summarize the test results.
        - outpath: path to a directory used for permanent output (final diagnostic
          figures, goodness-of-fit tables, html files etc.)
        - plotformat: defaults to 'png' is not given.
    '''

    obsid = None

    download_all = False
    '''Set to ``True`` to download all files in the Chandra archive.
    If ``False`` files will be downloaded one by one when they are needed.
    '''

    figures = OrderedDict()

    def __init__(self, env):
        self.env = env
        self.env['outpath'] = os.path.abspath(env['outpath'])
        self.basepath = os.path.abspath(os.path.join(self.env['temppath'],
                                                     self.name))
        if hasattr(self, 'obsid') and self.download_all:
            download_chandra(self.obsid, self.datapath)

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def datapath(self):
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
            File and complete filename. The filename is composed of the name
            of the test, the input ``name`` and the file format defined in
            ``env['plotformat']``.
        '''
        if not os.path.exists(os.path.join(self.env['outpath'], 'figures')):
            os.makedirs(os.path.join(self.env['outpath'], 'figures'))

        return os.path.join(self.env['outpath'], 'figures',
                            '{0}_{1}.{2}'.format(self.name, name,
                                                  self.env.get('plotformat', 'png')))

    @property
    def steplist(self):
        steplist = [method for method in dir(self) if (method[:5] == 'step_')]
        steplist.sort(key=lambda s: int(s.split('_')[1]))
        return [getattr(self, s) for s in steplist]

    def run(self):
        self.pkg_data = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                     'tests', 'data'))
        if not os.path.exists(self.env['outpath']):
                os.makedirs(self.env['outpath'])

        with ChangeDir(self.basepath):
            for step in self.steplist:
                step()

    def get_data_file(self, filetype):
        filename = glob(os.path.join(self.datapath,
                                     '*',  # primary or secondary
                                     '*{0}*'.format(filetype)))
        if len(filename) == 0:
            download_chandra(self.obsid, self.datapath, [filetype])
            return self.get_data_file(filetype)

        elif len(filename) == 1:
            return filename[0]

        else:
            raise ValueError('Specification not unique. Found {0}'.format(filename))
