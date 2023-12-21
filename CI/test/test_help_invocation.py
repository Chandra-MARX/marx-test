# Licensed under GPL 2 - see LICENSE file
import os
import subprocess
from contextlib import chdir

from .utils import check_no_warnings

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

def test_warning_asol_end(tmp_path, obsid11005):
    """Test output related to asol files

    - warning is raised if the asol file does not cover the entire exposure
    - Message when trying to general an asol file for an observations that used FILE dither
    """
    with chdir(tmp_path):
        out = subprocess.run(['marx', 'DitherModel=FILE',
                              f'DitherFile={obsid11005}/primary/pcadf11005_000N001_asol1.fits',
                              # Just so we use different detectors for different tests
                              'DetectorType=HRC-I',
                              'GratingType=LETG',
                              # Make the simulation run fast, don't need many photons
                              'SourceFlux=0.001',
                              ],
                             check=True, capture_output=True)
        outasp = subprocess.run(['marxasp'],
                                capture_output=True)
    check_no_warnings(out)
    assert 'Simulation stopped early because end of ASPSOL file was reached' in out.stdout.decode('utf-8')
    assert outasp.returncode == 255
    assert 'There is no need to run marxasp to generate a new ASPSOL file.' in outasp.stderr.decode('utf-8')


def test_verbosity(tmp_path):
    """Check verbosity setting does something
    """
    out = []
    with chdir(tmp_path):
        for i in range(3):
            out.append(subprocess.run(['marx',
                                       f'Verbose={i}',
                              # Just so we use different detectors for different tests
                              'DetectorType=HRC-S',
                              'GratingType=HETG',
                              # Make the simulation run fast, don't need many photons
                              'SourceFlux=0.001',
                              ],
                             check=True, capture_output=True))
    for o in out:
        check_no_warnings(o)
    # Check that there is no errand output
    # (e.g. print statements left behind from debugging)
    assert len(out[0].stdout) == 0
    # check that higher verbosity setting do something
    assert len(out[1].stdout) > len(out[0].stdout)
    assert len(out[2].stdout) > len(out[1].stdout)
