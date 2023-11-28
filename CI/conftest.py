# Licensed under GPL 2 - see LICENSE file
import os
import subprocess
from contextlib import chdir
import shutil
import glob
import gzip

import pytest
from ciao_contrib.cda.data import download_chandra_obsids

from test.utils import check_no_warnings


@pytest.fixture(scope="session", autouse=True)
def set_env():
    out = subprocess.run(['marx'], capture_output=True)
    if 'Unable to locate a marx.par parameter file' in out.stderr.decode('utf-8'):
        os.environ['PFILES'] = out.stderr.decode('utf-8').split('\n')[4].strip()


@pytest.fixture(scope="session")
def obsid11005(tmp_path_factory):
    """Provide one OBSID downloaded from the archive"""
    fn = tmp_path_factory.mktemp("data")
    with chdir(fn):
        assert download_chandra_obsids([11005])
        # Marx needs the asol file to be unzipped
        # So we make a copy
        asol = glob.glob('*/primary/*asol1.fits.gz')[0]
        with gzip.open(asol, 'rb') as f_in:
            with open(asol.replace('.gz', ''), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    return fn / "11005"


@pytest.fixture(scope="session")
def default_marx_sim(tmp_path_factory, marxpar):
    """Provide output for a single marx simulation"""
    fn = tmp_path_factory.mktemp("sim")
    for m in ['marx', 'marxasp', 'marxpileup']:
        shutil.copyfile(marxpar / f'{m}.par', fn / f'{m}.par')
    with chdir(fn):
        out = subprocess.run(['marx'],
                    check=True, capture_output=True)
        check_no_warnings(out)
        out = subprocess.run(['marx2fits', '--pixadj=EDSER',
                    'point',  'point.fits'],
                    check=True, capture_output=True)
        check_no_warnings(out)
        out = subprocess.run(['marxasp', "MarxDir=point",
                     'OutputFile=asol1.fits'],
                      check=True, capture_output=True)
        check_no_warnings(out)
    return fn