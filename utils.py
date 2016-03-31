from __future__ import print_function

import sys
import os
import subprocess
import copy


class ChangeDir(object):
    """ Context manager to step into a directory temporarily.

    See: https://pythonadventures.wordpress.com/2013/12/15/chdir-a-context-manager-for-switching-working-directories/
    """
    def __init__(self, directory):
        self.old_dir = os.getcwd()
        self.new_dir = directory

    def __enter__(self):
        os.chdir(self.new_dir)

    def __exit__(self, *args):
        os.chdir(self.old_dir)


def get_immidiate_subdirs():
    return next(os.walk('.'))[1]


def checkout_hash(directory, githash, j=1):
    with open(os.devnull, 'wb', 0) as DEVNULL:
        with ChangeDir(os.path.join(directory, 'dist')):
            print('{0}: git fetch'.format(directory), end=', ')
            sys.stdout.flush()
            subprocess.call(['git', 'fetch'])
                            #stdout=DEVNULL, stderr=subprocess.STDOUT)
            print('git checkout', end=', ')
            sys.stdout.flush()
            subprocess.call(['git', 'checkout', githash],
                            stdout=DEVNULL, stderr=subprocess.STDOUT)
            print('make', end=', ')
            sys.stdout.flush()
            subprocess.call(['make', '-j{0}'.format(j)],
                            stdout=DEVNULL, stderr=subprocess.STDOUT)
            print('make install')
            sys.stdout.flush()
            subprocess.call(['make', 'install'],
                            stdout=DEVNULL, stderr=subprocess.STDOUT)


def call_marx(tool='marx', marxdir='', parfile=None, **kwargs):
    '''Call marx from python

    This will raise an exception if marx does not terminate successfully.

    Parameters
    ----------
    tool : string
        'marx' or 'marx2fits' or any other of the compiled marx tools
    marxdir : string
        path to marx executable. Can be empty if marx is on the PATH.
    parfile : string or None
        marx parameter file. Set to ``None`` to use marx.par in the
        current directory.
    kwargs : various
        Any marx argument, e.g. ``SourceFlux=0.001``.
    '''
    marxcall = ['{0}={1}'.format(k, v) for k, v in kwargs.iteritems()]
    marxcall.insert(0, os.path.join(marxdir, tool))
    if parfile is not None:
        marxcall.insert(1, '@@{0}'.format(parfile))

    retcode = subprocess.check_call(marxcall)
