import os
import subprocess
from contextlib import chdir

def test_version_message():
    outprocess = subprocess.run(['marx', '--version'], capture_output=True, check=True)
    out = outprocess.stderr.decode('utf-8')
    assert 'MARX version' in out
    assert 'Supported Sources: POINT, GAUSS, LINE, BETA, RAYFILE, DISK, USER, SAOSAC, IMAGE, SIMPUT' in out
    assert 'Supported Detectors: ACIS-S, ACIS-I, HRC-S, HRC-I' in out
    for n in ['HRMA Pitch/Yaw', 'Wfold Scattering', 'Drake Flat', 'Dynamic Linking', 'ACIS Streak', 'ACIS FEF Support', 'Dither Support']:
        assert f'{n} : yes' in out

def test_no_par_file_noPFILE(tmp_path):
    """Run marx in empty directory and check that error is raised for no par file
    """
    env = os.environ.copy()
    if 'PFILES' in env:
        del env['PFILES']
    if 'UPARM' in env:
        del env['UPARM']
    with chdir(tmp_path):
        outprocess = subprocess.run(['marx'], capture_output=True, env=env)
    assert outprocess.returncode == 1
    out = outprocess.stderr.decode('utf-8')
    assert "Unable to open parameter file marx.par." in out
