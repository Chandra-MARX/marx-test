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

    def __call__(self, obj, conf):
        '''
        Parameters
        ----------
        conf : `~ConfigParser.ConfigParser` instance
            The configures the path to external programs like marx.
        '''
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
    '''Wrap functions that generate CIAO shell commands.

    The output format of the wrapped function is a list of strings.

    On display with ``source`` absolute path in the CIAO commands will be
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
    at the default set in marx.par in the marx installation.

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
    program = 'marxasp'


class Marx2fits(ExternalBaseWrapper):
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
    interpreter = "lua"
    program = "Lua input for SAOTrace"

    def source(self):
        return self.f(self.instance)

    def __call__(self, obj, conf):
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
    conf : `~ConfigParser.ConfigParser`
        Configuration Parser with appropriate configuration file loaded
    '''

    obsid = None

    download_all = False
    '''Set to ``True`` to download all files in the Chandra archive.
    If ``False`` files will be downloaded one by one when they are needed.
    '''

    figures = OrderedDict()

    version = 1
    '''If the definition of a test changes so that results cannot be
    compared to previous runs any longer, the version needs to be incremented.
    However, usually it is better to rename the test.
    '''

    expresults = []

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
        steplist = [method for method in dir(self) if (method[:5] == 'step_')]
        steplist.sort(key=lambda s: int(s.split('_')[1]))
        return [getattr(self, s) for s in steplist]

    def run(self):
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
        database.insert_test_run(self.conf, self.name, testname, self.version,
                                 value, sigma, sigma2)
