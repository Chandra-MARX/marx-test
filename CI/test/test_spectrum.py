# Licensed under GPL 2 - see LICENSE file
"""Test the energy of photons that are generated for different input spectra"""
import pycrates
import pytest

def test_default_energy(default_marx_sim):
    """Check that the default energy is 1.775 keV"""
    r = pycrates.read_file(str(default_marx_sim / 'point.fits'))
    # MARX_ENERGY is in keV
    assert r.get_column('MARX_ENERGY').values.min() == pytest.approx(1.775)
    assert r.get_column('MARX_ENERGY').values.max() == pytest.approx(1.775)
    # ENERGY is in eV
    assert r.get_column('ENERGY').values.mean() == pytest.approx(1775, abs=10, rel=0.01)
    assert r.get_column('ENERGY').values.std() == pytest.approx(100, rel=0.1)
    assert r.get_column('PHA').values.mean() == pytest.approx(369, abs=10, rel=0.01)
    assert r.get_column('PHA').values.std() == pytest.approx(26, rel=0.2)
    assert r.get_column('PI').values.mean() == pytest.approx(121.5, abs=10, rel=0.01)
    assert r.get_column('PI').values.std() == pytest.approx(7.1, rel=0.1)
    assert r.get_column('GRADE').values.min() == 0
    assert r.get_column('GRADE').values.max() == 6
    assert r.get_column('FLTGRADE').values.min() == 0
    assert r.get_column('FLTGRADE').values.max() == 209
    # And some other checks that are not energy related
    # but that we might do as well while the file is read in already
    assert r.get_column('CCD_ID').values.min() == 4
    assert r.get_column('CCD_ID').values.max() == 9
    assert r.get_column('NODE_ID').values.min() == 0
    assert r.get_column('NODE_ID').values.max() == 3
    assert r.get_column('SHELL').values.min() == 0
    assert r.get_column('SHELL').values.max() == 3
    # Exact numbers depends on Poisson draw
    assert r.get_column('EXPNO').values.min() < 3
    assert r.get_column('EXPNO').values.max() > 3000
    
    