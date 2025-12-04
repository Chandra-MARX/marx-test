# Licensed under GPL 2 - see LICENSE file
"""Test for different build-in marx sources.

As with any Monte-Carlo code it is difficult to verify correctness
with limited run time. Also, output changes as the CALDB information
that is used by marx is updated. To avoid having to update the tests
every time the CALDB is updated, many of the numerical comparisons
are relatively loose.
"""
from contextlib import chdir
import gzip
import shutil
import glob
import os
import re
import subprocess

from astropy.io import fits
import numpy as np
import pytest
import pycrates
from coords.chandra import sky_to_chandra
from ciao_contrib.cda.data import download_chandra_obsids
from ciao_contrib import runtool as rt

from .utils import check_no_warnings


nphotons = re.compile(r"Total photons: (?P<nphotons>[0-9]+), " +
                       "Total Photons detected: (?P<ndetected>[0-9]+)")


def test_point_source_no_grating(tmp_path):
    """Some regression tests of the default setup."""
    with chdir(tmp_path):
        out = subprocess.run(['marx', 'GratingType=NONE', "Verbose=2"],
                              capture_output=True, check=True)
        assert 'hrma/corr_3.dat' in out.stdout.decode('utf-8')
        check_no_warnings(out)

        out = subprocess.run(['marx2fits', '--pixadj=EDSER', 'point',  'point.fits'],
                              check=True, capture_output=True)
        check_no_warnings(out)

        r = pycrates.read_file('point.fits')
        x = r.get_column('x').values
        y = r.get_column('y').values
        # Very rough number so that it still passes when ACIS contamination
        # increases.
        assert 5000 < len(x) < 20000
        assert 4015 < np.median(x) < 4020
        # Loose number in case scattering or dither uncertainties change.
        assert 20 < x.std() < 50
        assert 4015 < np.median(x) < 4020
        assert 20 < y.std() < 50


def test_S_IMAGE_position(tmp_path):
    """Check the position of an image source.
    
    This test uses an input image where only a single pixel is non-zero.

    See https://github.com/Chandra-MARX/marx/issues/50
    """
    data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

    with chdir(tmp_path):
        assert download_chandra_obsids([7901], filetypes='asol')
        
        # marx does not read gz files, so unzip
        asol = glob.glob('*/primary/*asol1.fits.gz')[0]
        with gzip.open(asol, 'rb') as f_in:
            with open(asol.replace('.gz', ''), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        out = subprocess.run(['marx',
                        'ExposureTime=37144.58925002813',
                        'OutputDir=delta_out',
                        'GratingType=NONE',
                        'SourceFlux=0.01',
                        'MinEnergy=1',
                        'MaxEnergy=1',
                        'SourceRA=195.9341611742938',
                        'SourceDEC=-24.22948889322631',
                        'SourceType=IMAGE',
                        f'S-ImageFile={data_path}/one_pixel.img',
                        'DitherModel=FILE',
                        f'DitherFile={asol.removesuffix(".gz")}',
                        'DetOffsetX=0.001444942264670734',
                        'DetOffsetZ=-0.01005726120280315',
                        'RA_Nom=195.93423610867',
                        'Dec_Nom=-24.229557226578',
                        'Roll_Nom=80.721746267494',
                        'ACIS_Frame_Transfer_Time=0',  # No read-out streak that would confuse the mean position],
        ],
                              check=True, capture_output=True)
        check_no_warnings(out)
        out = subprocess.run(['marx2fits', '--pixadj=EDSER',
                              'delta_out',  'delta_out.fits'],
                             check=True, capture_output=True)
        check_no_warnings(out)

        r = pycrates.read_file('delta_out.fits[cols x,y]')

    rt.dmstat(f'{data_path}/one_pixel.img')
    x_pix, y_pix = rt.dmstat.out_cntrd_phys.split(',')
    x_pix = float(x_pix)
    y_pix = float(y_pix)
    assert pytest.approx(x_pix, abs=0.1) == r.get_column('x').values.mean()
    assert pytest.approx(y_pix, abs=0.1) == r.get_column('y').values.mean()


def test_image_rotation(tmp_path):
    data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

    with chdir(tmp_path):
        assert download_chandra_obsids([7901], filetypes='asol')

        # marx does not read gz files, so unzip
        asol = glob.glob('*/primary/*asol1.fits.gz')[0]
        with gzip.open(asol, 'rb') as f_in:
            with open(asol.replace('.gz', ''), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        out = subprocess.run(['marx',
                        'ExposureTime=37144.58925002813',
                        'OutputDir=asym_shape',
                        'GratingType=NONE',
                        'SourceFlux=0.01',
                        'MinEnergy=1',
                        'MaxEnergy=1',
                        'SourceRA=195.9341611742938',
                        'SourceDEC=-24.22948889322631',
                        'SourceType=IMAGE',
                        f'S-ImageFile={data_path}/asym_shape.img',
                        'DitherModel=FILE',
                        f'DitherFile={asol.removesuffix(".gz")}',
                        # Note:
                        # Compared to test_S_IMAGE_position we do no give
                        # DetOffsetX/Z and RA/Dec/roll nom here.
                        # These are taken from the asol file directly.
                        'ACIS_Frame_Transfer_Time=0',  # No read-out streak that would confuse the mean position],
        ],
                              check=True, capture_output=True)
        check_no_warnings(out)
        out = subprocess.run(['marx2fits', '--pixadj=EDSER',
                              'asym_shape',  'asym_shape.fits'],
                             check=True, capture_output=True)
        check_no_warnings(out)
        # Bin with the same specification that were used to generate the input image
        rt.dmcopy('asym_shape.fits[bin sky=4046.5:4146.5:1]', 'asym_shape_cut.fits', clobber=True, opt='image')
        # Threshold such that the image only has 1 and 0, like the input image
        rt.dmimgthresh('asym_shape_cut.fits', 'new.img', cut="50%")
        rt.dmimgthresh('new.img', 'new2.img', cut=":50%", value=1)

        # compare pixel by pixel
        orig = pycrates.read_file(f'{data_path}/asym_shape.img')
        new = pycrates.read_file('new2.img')
        assert orig.get_image().values == pytest.approx(new.get_image().values)


@pytest.mark.parametrize("SourceType", ["POINT", "IMAGE"])
def test_error_outside_of_FOV(tmp_path, SourceType):
    """Check the position of an image source.

    This test uses an input image where only a single pixel is non-zero.

    See https://github.com/Chandra-MARX/marx/issues/50
    """
    with chdir(tmp_path):
        out = subprocess.run(['marx',
                              'ExposureTime=1000',
                              'SourceRA=195.9',
                              'SourceDEC=-24.2',
                              f'SourceType={SourceType}',
                              'RA_Nom=25.',
                              'Dec_Nom=71.',

        ],
                              capture_output=True)
        assert out.returncode == 1
        assert "Source is located at RA=195.900000, Dec=-24.200000 but nominal pointing is RA=25.000000 Dec=71.000000." in out.stderr.decode('utf-8')
        assert "Rays will not hit the telescope." in out.stderr.decode('utf-8')


def test_rayfiles(tmp_path):
    """Generate a rayfile and then use it with different detectors.

    """
    with chdir(tmp_path):

        outdump = subprocess.run(['marx',
                                  'DumpToRayFile=yes',
                                  'SourceType=POINT',
                                  'MinEnergy=0.5',
                                  'MaxEnergy=2.5',
                                  # Low flux so it's unlikely to have more than one
                                  # photon per read-out
                                  'SourceFlux=0.00001',],
                                  capture_output=True, check=True)
        # https://github.com/Chandra-MARX/marx/issues/30
        assert 'Diffracting from HETG' not in outdump.stdout.decode('utf-8')
        assert 'Detecting with ACIS-S' not in outdump.stdout.decode('utf-8')
        assert 'Diffracting from HETG' not in outdump.stdout.decode('utf-8')
        assert '(efficiency: 1.000000)' in outdump.stdout.decode('utf-8')
        n_dump = nphotons.search(outdump.stdout.decode('utf-8'))
        n_generated = int(n_dump.groupdict()['nphotons'])
        assert n_generated == int(n_dump.groupdict()['ndetected'])
        check_no_warnings(outdump)

        # Switch off most sources of randomness
        # However, there are still some effects that make that
        # the simulations don't match exactly.
        # For example, there is the scatter in the HRMA, etc.
        # I did not actually trace down all effects because that's not
        # needed for this test.
        # A significant fraction (~70%) of the photons in the rayfile never
        # hits the detector, because they never hit the mirror shells in the HRMA.
        # Thus, this is a small subset of photons that are available for detection and
        # many of those overlap between the two runs.
        marx_args = ['marx',
                     'GratingType=NONE',
                     'SourceType=RAYFILE',
                     'DetectorType=ACIS-I',
                     'DetIdeal=yes',
                     'HRMA_Ideal=yes',
                     'HRMA_Use_Blur=no']
        outi = subprocess.run(marx_args + ['OutputDir=acisi'],
                               capture_output=True, check=True)
        assert 'Detecting with ACIS-I' in outi.stdout.decode('utf-8')
        check_no_warnings(outi)
        n_outi = nphotons.search(outi.stdout.decode('utf-8'))
        n_outi = int(n_outi.groupdict()['ndetected'])
        assert n_generated > n_outi

        # Run again with a different random seed
        outi2 = subprocess.run(marx_args + [
                               'OutputDir=acisi2',
                               'RandomSeed=12345'],
                               capture_output=True, check=True)
        assert 'Detecting with ACIS-I' in outi2.stdout.decode('utf-8')
        check_no_warnings(outi2)
        n_outi2 = nphotons.search(outi2.stdout.decode('utf-8'))
        n_outi2 = int(n_outi2.groupdict()['ndetected'])
        assert  n_generated > n_outi2

        outg = subprocess.run(['marx',
                               'GratingType=LETG',
                               'SourceType=RAYFILE',
                               'OutputDir=letg',
                               'DetectorType=ACIS-S'],
                               capture_output=True, check=True)
        assert 'Detecting with ACIS-S' in outg.stdout.decode('utf-8')
        check_no_warnings(outg)
        n_outg = nphotons.search(outg.stdout.decode('utf-8'))
        n_outg = int(n_outg.groupdict()['ndetected'])
        assert  n_outg < n_outi
        for dir in ['acisi', 'acisi2', 'letg']:
            out = subprocess.run(['marx2fits', '--pixadj=EXACT', dir,
                            f'{dir}.fits'], check=True, capture_output=True)
            check_no_warnings(out)
            assert n_generated > int(nphotons.search(outi.stdout.decode('utf-8')).group('ndetected'))

        acisi = pycrates.read_file('acisi.fits')
        acisi2 = pycrates.read_file('acisi2.fits')
        letg = pycrates.read_file('letg.fits')

        # Photons are all in the same location, but not all photons might be detected in both cases
        # So, we need to select carefully which ones to compare
        t = acisi.get_column('time').values
        t2 = acisi2.get_column('time').values
        unique, ind1, ind2 = np.intersect1d(t, t2, return_indices=True)

        # Sometimes, there are two photons with the same time stamp
        # and they might have different energies.
        # So, to make the test robust, allow a few mismatches.
        delta_E = acisi.get_column('MARX_ENERGY').values[ind1] - acisi2.get_column('MARX_ENERGY').values[ind2]
        assert (delta_E == 0).sum() > 0.9 * len(delta_E)

        # Due to HRMA scatter etc., the positions might not match exactly
        delta_x = acisi.get_column('chipx').values[ind1] - acisi2.get_column('chipx').values[ind2]
        assert (delta_x == 0).sum() > 0.3 * len(delta_x)
        assert (np.abs(delta_x <2)).sum() > 0.9 * len(delta_x)

        # some PHA might match by accident, but usually they should be different with a
        # different random seed
        assert (acisi.get_column('PHA').values[ind1] == acisi2.get_column('PHA').values[ind2]).sum() < 0.9 * n_outi

        # grating spreads out photons, this is just a simple test to ensure the grating is in fact simulated
        assert letg.get_column('x').values.std() > acisi.get_column('x').values.std()
        assert letg.get_column('y').values.std() > acisi.get_column('y').values.std()

        # Should not issue a warning message, even if SourceRA/DEC differs from nominal pointing.
        ra_nom = 25.
        dec_nom = 77.
        out = subprocess.run(['marx',
                               'GratingType=LETG',
                               'SourceType=RAYFILE',
                               'OutputDir=pointing',
                               'DetectorType=ACIS-S',
                               f'SourceRA={ra_nom+180}',
                               f'SourceDEC={-dec_nom}',
                               f'RA_Nom={ra_nom}',
                               f'Dec_Nom={dec_nom}',],
                               capture_output=True, check=True)
        check_no_warnings(out)
        subprocess.run(['marx2fits', '--pixadj=EXACT', 'pointing',
                        'pointing.fits'], check=True, capture_output=True)

        # Positions should be very close to the nominal pointing case
        pointing = pycrates.read_file('pointing.fits')
        coos = sky_to_chandra(fits.getheader('pointing.fits', 1),
                                pointing.get_column('x').values,
                                pointing.get_column('y').values)
        assert [ra_nom, dec_nom] == pytest.approx([np.mean(coos['ra']),
                                                   np.mean(coos['dec'])], abs=0.2)