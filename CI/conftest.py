import os
import subprocess
from contextlib import chdir
import shutil

import pytest
from ciao_contrib.cda.data import download_chandra_obsids

@pytest.fixture(scope="session")
def obsid11005(tmp_path_factory):
    """Provide one OBSID downloaded from the archive"""
    fn = tmp_path_factory.mktemp("data")
    with chdir(fn):
        assert download_chandra_obsids([11005])
    return fn / "11005"


@pytest.fixture(scope="session")
def marxpar(tmp_path_factory):
    """Directory with MARX parameter files"""
    env = os.environ.copy()
    if 'PFILES' in env:
        del env['PFILES']
    if 'UPARM' in env:
        del env['UPARM']
    outprocess = subprocess.run(['marx'], capture_output=True, env=env)
    parpath = outprocess.stderr.decode('utf-8').split('\n')[4].strip()
    outpath = tmp_path_factory.mktemp("pars")
    for m in ['marx', 'marxasp', 'marxpileup']:
        shutil.copyfile(parpath + f'{m}.par', outpath / f'{m}.par')
    return outpath


@pytest.fixture(scope="session")
def default_marx_sim(tmp_path_factory, marxpar):
    """Provide output for a single marx simulation"""
    fn = tmp_path_factory.mktemp("sim")
    for m in ['marx', 'marxasp', 'marxpileup']:
        shutil.copyfile(marxpar / f'{m}.par', fn / f'{m}.par')
    with chdir(fn):
        subprocess.run(['marx'],
                    check=True)
        subprocess.run(['marx2fits', '--pixadj=EDSER',
                    'point',  'point.fits'],
                   check=True)
        subprocess.run(['marxasp', "MarxDir=point",
                     'OutputFile=asol1.fits'],
                    check=True)
    return fn