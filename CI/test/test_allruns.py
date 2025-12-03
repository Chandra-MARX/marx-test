# Licensed under GPL 2 - see LICENSE file
"""Test that all detectors and gratings run without errors.

This does not test that the output is correct, only that no errors
are raised during the simulation.
"""
import subprocess
from contextlib import chdir

from astropy.table import Table
import numpy as np
import pytest

from .utils import check_no_warnings

@pytest.mark.parametrize("detector", ["ACIS-I", "ACIS-S", "HRC-I", "HRC-S"])
@pytest.mark.parametrize("grating", ["NONE", "LETG", "HETG"])
def test_all_run(tmp_path, detector, grating):
    """
    """
    with chdir(tmp_path):
        out = subprocess.run(['marx', f'DetectorType={detector}',
                               f'GratingType={grating}',
                              # Make the simulation run fast, don't need many photons
                              'SourceFlux=0.001',
                              ],
                             check=True, capture_output=True)
        assert out.returncode == 0
        check_no_warnings(out)
        outasp = subprocess.run(['marxasp'],
                                capture_output=True)
        assert outasp.returncode == 0
        check_no_warnings(outasp)
        out2fits = subprocess.run(['marx2fits', 'point', 'point.fits'],
                                     check=True, capture_output=True)
        assert out2fits.returncode == 0
        # Some generic checks that the output file looks reasonable
        evt = Table.read('point.fits', hdu=1)
        assert len(evt) > 100  # at least 100 events

        if grating == 'NONE':
            assert np.std(evt['X']) > 10  # not all events at same location
            assert np.std(evt['Y']) > 10  # not all events at same location
            assert np.std(evt['X']) < 500 # but not as spread out as a grating would
        else:
            assert 'ORDER' in evt.colnames  # grating orders present
            np.std(evt['X']) > 5000
