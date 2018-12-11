from __future__ import print_function

import sys
import os
import subprocess
import re
import shutil
import gzip
import numpy as np
import ftputil


class ChangeDir(object):
    """ Context manager to step into a directory temporarily.

    Parameters
    ----------
    directory : string
        Directory to step into.

    make_dir : boolean
        If True, the directory is created should it not exist yet.

    See: https://pythonadventures.wordpress.com/2013/12/15/chdir-a-context-manager-for-switching-working-directories/
    """
    def __init__(self, directory, make_dir=True):
        self.old_dir = os.getcwd()
        self.new_dir = directory
        self.make_dir = make_dir

    def __enter__(self):
        if self.make_dir:
            if not os.path.exists(self.new_dir):
                os.makedirs(self.new_dir)

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


def _name_matches(name, products):
    '''Check is name is match by any of the patters in the list products'''
    if products is None:
        return True
    for p in products:
        if re.search(p, name):
            return True
    return False


def _copy_into_dir(host, t_dir, products):
    '''copy products from current host dir into t_dir

    Unzip all files ending on .gz. This routine skips files that already exist
    in the target directory with the same filename.
    '''
    names = host.listdir(host.curdir)
    for name in names:
        if host.path.isfile(name) and _name_matches(name, products):
            if not os.path.exists(t_dir):
                os.makedirs(t_dir)
            # Remote name, local name, binary mode
            outname = os.path.join(t_dir, name)
            if not os.path.exists(outname[:-3]):
                host.download(name, outname)
                if name[-3:] == '.gz':
                    with gzip.open(outname, 'rb') as f_in, open(outname[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    os.remove(outname)


def download_chandra(obsid, targetdir, products=None):
    '''Download Chandra data.

    Files are retrieved from the Chandra ftp archive and unzipped.
    This routine skips files that already exist
    in the target directory with the same filename.

    Parameters
    ----------
    obsid : int or str
        ObsID to download
    targetdir : str
        Directory to save products in
    products : None or list of strings
        If not ``None`` then this should be a list of strings that will be
        searched for in the filenames, e.g. ``products=["fits"]`` will only
        download fits files and ``products=["evt2", "asol"]`` will only
        download event files and the asol file.
    '''
    obsid = str(obsid)
    dataroot = 'pub/byobsid/{0}/{1}/'.format(obsid[-1], obsid)
    # Download some files from the login directory.
    with ftputil.FTPHost('cda.cfa.harvard.edu', 'anonymous', 'marx-help@space.mit.edu') as host:
        host.chdir(dataroot)
        _copy_into_dir(host, targetdir, products)
        host.chdir('primary')
        _copy_into_dir(host, os.path.join(targetdir, 'primary'), products)
        host.chdir('../secondary')
        _copy_into_dir(host, os.path.join(targetdir, 'secondary'), products)


def find_centroid(x, x0, width, func=np.mean):
    '''Find centroid of a distribution

    Parameters
    ----------
    x : np.array
    x0 : float
        starting guess
    width : float
        use only data values in ``x0 - width, x0 + width``
    func : callable
    '''
    ind = np.abs(x - x0) < width
    return func(x[ind])


def colname_case(table, name):
    '''Extract a column from a Table case insensitive to the column name.

    If several columns in the table match ``name`` in a case-insensitive way,
    this just returns the first match.

    Parameters
    ----------
    table : `astropy.table.Table`
        Table with columns
    name : string
        column name, case insensitive.

    Returns
    -------
    col : `astropy.table.Column`
        Table column
    '''
    for n in table.colnames:
        if name.lower() == n.lower():
            return table[n]
    else:
        raise KeyError('{0} not found in {1}'.format(name, table.colnames))
