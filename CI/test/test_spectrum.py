# Licensed under GPL 2 - see LICENSE file
"""Test the energy of photons that are generated for different input spectra"""
from contextlib import chdir
import os
import subprocess

import pycrates
import pytest
from ciao_contrib import runtool as rt
from sherpa.astro import ui

from .utils import check_no_warnings

def test_default_energy(default_marx_sim):
    """Check that the default energy is 1.775 keV

    This implicitly tests that the ACIS module (at the time of writing using the FEF formalism)
    converts marx energy to PI, PHA, and energy correctly.
    """
    r = pycrates.read_file(str(default_marx_sim / 'point.fits'))
    # MARX_ENERGY is in keV
    assert r.get_column('MARX_ENERGY').values.min() == pytest.approx(1.775)
    assert r.get_column('MARX_ENERGY').values.max() == pytest.approx(1.775)
    # ENERGY is in eV
    assert r.get_column('ENERGY').values.mean() == pytest.approx(1775, abs=10, rel=0.01)
    assert r.get_column('ENERGY').values.std() == pytest.approx(100, rel=0.3)
    assert r.get_column('PHA').values.mean() == pytest.approx(369, abs=10, rel=0.01)
    assert r.get_column('PHA').values.std() == pytest.approx(26, rel=0.3)
    assert r.get_column('PI').values.mean() == pytest.approx(121.5, abs=10, rel=0.01)
    assert r.get_column('PI').values.std() == pytest.approx(7.1, rel=0.3)
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
    

def test_absorbed_powerlaw_spectrum(tmp_path):
    """Simulate an absorbed powerlaw spectrum for a point source

    Then use CIAO to extract and fit back that parameters.
    This tests energy conversion, spectral shape, contamination,
    and output formats so it's pretty comprehensive.
    On the other hand, small changes (say, a slip by one bin)
    would be missed because a powerlaw is smooth and the
    assert statements are not too stringent to prevent test
    failures for numerical reasons or for minor changes in the CALDB.
    """
    data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

    with chdir(tmp_path):
        out = subprocess.run(['marx ExposureTime=100000 TStart=2015.5 GratingType=NONE ' +
                        'SourceFlux=-1 SpectrumType="FILE" ' +
                        f'SpectrumFile={data_path}/abspowerlaw.spec ' +
                        'RA_Nom=30. Dec_Nom=40. SourceRA=30 SourceDEC=40',
        ],
                             shell=True,
                              check=True, capture_output=True)
        check_no_warnings(out)
        out = subprocess.run(['marx2fits', '--pixadj=EDSER',
                              'point',  'point.fits'],
                             check=True, capture_output=True)
        check_no_warnings(out)
        out = subprocess.run(['marxasp', 'MarxDir=point',
                     'OutputFile=asol1.fits'],
                      check=True, capture_output=True)
        check_no_warnings(out)

        phagrid = "pi=1:1024:1"
        rt.asphist(infile="asol1.fits", outfile="asp.fits", evtfile="point.fits")
        rt.dmextract(infile=f"point.fits[sky=circle(4096.5,4096.5,20)][bin {phagrid}]",
                     outfile="pha.fits")

        # For ACIS-I, use engrid="0.3:11.0:0.003". This reflects a limitation of mkrmf.
        engrid="0.3:12.0:0.003"

        rt.mkarf(mirror="hrma", detsubsys="ACIS-7;UNIFORM;bpmask=0", grating="NONE",
                 outfile="arf.fits", obsfile="point.fits",
                 engrid=engrid, asphistfile="asp.fits",
                 sourcepixelx=4096.5, sourcepixely=4096.5,
                 maskfile=None, pbkfile=None, dafile=None, verbose=0)


        # can get those numbers with dmcoords or ciao_contrib.cords.chandra.sky_to_chandra
        # but for this test, we hardcode all locations, so can hardcode this, too.
        chipx=220.7
        chipy=531.8
        fef=rt.acis_fef_lookup('point.fits', 7, chipx, chipy)
        rt.mkrmf(infile=fef,
                 outfile="rmf.fits", axis1=f"energy={engrid}", axis2=phagrid,
                 verbose=0)

        ui.load_pha('pha.fits')
        ui.load_arf('arf.fits')
        ui.load_rmf('rmf.fits')
        abs1 = ui.xsphabs(name='abs1')
        pl = ui.xspowerlaw(name='pl')
        ui.set_source(abs1 * pl)
        ui.group_counts(1, 10)
        ui.fit()

        assert abs1.nH.val == pytest.approx(1.0, rel=0.2)
        assert pl.PhoIndex.val == pytest.approx(1.8, rel=0.2)
        assert pl.norm.val == pytest.approx(0.001, rel=0.2)

        # There is some ambiguity between nH and norm. To get a more stringent test,
        # we can fix nH and check the norm.
        abs1.nH.frozen = True
        abs1.nH.val = 1.0
        ui.fit()
        assert pl.norm.val == pytest.approx(0.001, rel=0.05)
        assert pl.PhoIndex.val == pytest.approx(1.8, rel=0.05)


