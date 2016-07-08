import os
import functools
import inspect
import types
import subprocess
from glob import glob

from .utils import download_chandra, ChangeDir
from .run_external import external_settings

from astropy.io import fits
from astropy.coordinates import SkyCoord


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

    interpreter = 'undefined'
    program = 'unknown'

    def __init__(self, f):
        self.f = f
        functools.update_wrapper(self, f)

    def source(self):
        raise NotImplementedError

    def html(self):
        print 'test text' + self.source() + 'test text 2'

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

    interpreter = "python"
    program = 'Python'

    def source(self):
        return inspect.getsourcelines(self.f)[2:]

    def __call__(self, obj):
        self.f(obj)


class Ciao(ExternalBaseWrapper):

    interpreter = "shell"
    program = 'CIAO'

    def source(self):
        return self.f(self.instance)

    def __call__(self, obj):
        run_external(self.f(obj), setup=self.program, cwd=obj.basepath)


class Marx(ExternalBaseWrapper):

    interpreter = "shell"
    program = 'marx'

    def source(self):
        '''Assemble string for marx call.'''
        par = self.f(self.instance)
        marxcall = ['{0}={1}'.format(k, v) for k, v in par.iteritems()]
        marxcall.insert(0, os.path.join(external_settings['marxpath'],
                                        self.program))
        marxcall.insert(1, '@@{0}'.format(external_settings[self.program + 'parfile']))

        return marxcall

    def __call__(self, obj):
        run_external(self.source(), cwd=obj.basepath)


class Marx2asp(Marx):
    program = 'marx2asp'


class Marx2fits(ExternalBaseWrapper):
    interpreter = "shell"
    program = 'marx2fits'

    def source(self):
        '''Assemble string for marx2fits call'''
        options, marxdir, outfile = self.f(self.instance)
        return ' '.join(os.path.join(external_settings['marxpath'], self.program),
                        options, marxdir, outfile)

    def __call__(self, obj):
        run_external(self.source(), cwd=obj.basepath)


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
    '''

    obsid = None

    download_all = False
    '''Set to ``True`` to download all files in the Chandra archive.
    If ``False`` files will be downloaded one by one when they are needed.
    '''

    def __init__(self, env):
        self.env = env
        if hasattr(self, 'obsid') and self.download_all:
            download_chandra(self.obsid, self.datapath)

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def basepath(self):
        return os.path.join(self.env['temppath'], self.name)

    @property
    def datapath(self):
        return os.path.join(self.basepath, 'download')

    def run(self, label, marxversion, ciaoversion):
        with ChangeDir(self.basepath):
            steplist = [method for method in dir(self) if (method[:5] == 'step_')]
            steplist.sort(key=lambda s: int(s.split('_')[1]))
            for step in steplist:
                getattr(self, step)()

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

    def marxpars_from_asol(self, asolfile):

        asol = fits.getheader(asolfile)

        skyco = SkyCoord.from_name(asol.meta['OBJECT'].split('/')[0])
        # Set sensible default parameters
        marx_pars = {
                     # keep things simple so that sky and det coos are aligned
                     'SourceRA': skyco.ra.value,
                     'SourceDEC': skyco.dec.value,
                     'RA_Nom': asol['RA_NOM'],
                     'Dec_Nom': asol['DEC_NOM'],
                     'Roll_Nom': asol['ROLL_NOM'],
                     'GratingType': asol['GRATING'],
                     'ExposureTime': asol['TSTOP'] - asol['TSTART'],
                     'DitherModel': 'FILE',
                     'DitherFile': asolfile,
                     'TStart': asol['TSTART'],
                    }
        return marx_pars
