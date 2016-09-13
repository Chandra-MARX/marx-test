'''
|marx| offers several different source shapes. Tests in this module exercise
those sources (except ``SAOSAC``, which is heavily used in
:ref:`sect-PSF` already).


'''
import shutil
import subprocess
import os
from marxtest import base

title = 'Sources in |marx|'

tests = ['GeometricSources', 'ImageSource', 'UserSource']


class GeometricSources(base.MarxTest):
    '''This test exercises build-in |marx| sources with different geometric
    shapes.
    Most source types have parameters, and not all parameters are tested here.
    See :ref:`sect-sourcemodels` for a detailed description of source
    parameters.
    '''

    title = 'Build-in geometric sources'

    figures = OrderedDict([('ds9', {'alternative': 'Six PSFs.',
                                    'caption': '`ds9`_ image of the simulated PSFs in alphabeticcal order (beta distribution, disk, disk with hole, Gauss, line, and point.'})
                       ])

    @base.Marx
    def step_10(self):
        return [{'OutputDir': 'point'},
                {'SourceType': 'GAUSS', 'S-GaussSigma': 20,
                 'OutputDir': 'gauss'},
                {'SourceType': 'BETA', 'S-BetaCoreRadius': 10,
                 'S-BetaBeta': 0.6, 'OutputDir': 'beta'},
                {'SourceType': 'DISK',
                 'S-DiskTheta0': 0, 'S-DiskTheta1': 20,
                 'OutputDir': 'disk'},
                {'SourceType': 'DISK',
                 'S-DiskTheta0': 10, 'S-DiskTheta1': 20,
                 'OutputDir': 'diskhole'},
                {'SourceType': 'LINE', 'S-LinePhi': 45, 'S-LineTheta': 30,
                 'OutputDir': 'line'},
                ]

    # more to come for SAOSAC, RAYFILE, SIMPUT, USER
    # but first make something work here
    @base.Marx2fits
    def step_20(self):
        dirs = ['point', 'gauss', 'beta', 'disk', 'diskhole',
                'line']
        return ['--pixadj=EDSER'] * len(dirs), dirs, [d + '.fits' for d in dirs]

    @base.Ciao
    def step_30(self):
        '''ds9 images of the PSF'''
        return ['''ds9 -width 800 -height 500 -log -cmap heat *.fits -pan to 4018 4141 physical -match frame wcs -saveimage {0} -exit'''.format(self.figpath(self.figures.keys()[0]))]


class ImageSource(base.MarxTest):
    '''An image can be used as |marx| input. In this case, the intensity of the
    X-ray radiation on that sky is taken to be proportional to the value of the
    image at that point.
    '''

    title = 'Image as source'

    figures = OrderedDict([('ds9', {'alternative': 'The simulated events generally follow the input image, but with significant noise because of the short observation time.',
                                    'caption': '`ds9`_ shows the input image (left) and the simulated event list (right).'})
                       ])

    @base.Python
    def step_0(self):
        '''Make input image

        In this example we use python to make a simple image as input.
        We setup a 3-d box and fill it with an emitting shell. We then
        integrate along one dimension to obtain a collapsed image.
        Physically, this represents the thin shell of a supernova
        explosion.
        '''
        import numpy as np
        from astropy.wcs import WCS
        from astropy.io import fits

        # Actually to make this run faster, we'll do only one quadrant here
        cube = np.zeros((201, 201, 201))
        mg = np.mgrid[0: 201., 0:201, 0:201 ]
        d = np.sqrt(mg[0, :, :, :]**2 + mg[1, :, :, :]**2 + mg[2, :, :, :]**2)
        cube[(d > 160.) & (d < 170)] = 1
        im = cube.sum(axis=0)
        # Now rotate and put the four quarters together
        image = np.zeros((401, 401))
        image[:201, :201] = np.fliplr(np.flipud(im))
        image[:201, 200:] = np.flipud(im)
        image[200:, :201] = np.fliplr(im)
        image[200:, 200:] = im

        # Create a new WCS object.
        w = WCS(naxis=2)
        w.wcs.crpix = [-234.75, 8.3393]
        # Pixel size of our image shall be 1 arcsec
        w.wcs.cdelt = [1. / 3600., 1. / 3600.]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        # Now, write out the WCS object as a FITS header
        header = w.to_header()

        # header is an astropy.io.fits.Header object.  We can use it to create a new
        # PrimaryHDU and write it to a file.
        hdu = fits.PrimaryHDU(header=header, data=image)

        # Save to FITS file
        hdu.writeto(os.path.join(self.basepath, 'input_image.fits'), clobber=True)

    @base.Marx
    def step_1(self):
        '''Run |marx|.

        We run a monoenergetic simulation here for the Si XIII line at 6.65 Ang.
        '''
        return {'SourceType': "IMAGE", 'S-ImageFile': 'input_image.fits',
                'MinEnergy': 1.9, 'MaxEnergy': 1.9,
                'OutputDir': 'image'}

    @base.Marx2fits
    def step_2(self):
        return '--pixadj=EDSER', 'image', 'image.fits'

    @base.Ciao
    def step_30(self):
        '''ds9 images of the PSF'''
        return ['''ds9 -width 800 -height 500 -log -cmap heat input_image.fits image.fits -pan to 4018 4141 physical -saveimage {0} -exit'''.format(self.figpath(self.figures.keys()[0]))]


class RayfileSource(base.MarxTest):
    '''|marx| is a Monte-Carlo code, thus the exact distribution of photons
    on the sky will be different every time the code is run. Sometimes it
    can be useful to generate a list of photons with position, time and
    energy from the source on the sky and then "observe" the exact same list
    with different instrument configurations so that any differences in the
    result are only due to the different configuration and not to random
    fluctuations in the source.

    In this example, we look at a relatively large, diffuse emission region
    with a very soft spectrum (for simplicity we are using a flat spectrum).
    We compare a simulation using ACIS-S and ACIS-I. ACIS-S has a better
    response to soft photons, but some parts of the source may not be in the
    field-of-view; ACIS-I is less efficient for soft photons, but has a
    larger field-of-view.
    '''

    title = 'Using a RAYFILE source'

    figures = OrderedDict([('ds9', {'alternative': 'As described above, ACIS-S shows more photons, but ACIS-I does include more the wings of the GAussian source distribution',
                                    'caption': '`ds9`_ shows the ACIS-I (left) and ACIS-S image (right). Both sources are generated from the same photon list. Sometimes the same pattern of photons can be seen in both images, but with a few events missing on ACIS-I due to the lower soft response.'})
                       ])

    @base.Marx
    def step_1(self):
        '''Write ray file

        '''
        return {'SourceType': 'GAUSS', 'S-GaussSigma': 300,
                'DumpToRayFile': 'yes', 'MinEnergy': 0.3, 'MaxEnergy': 0.5}

    @base.Marx
    def step_2(self):
        '''ACIS-S'''
        return {'SourceType': 'RAYFILE', 'RayFile': 'marx.output',
                'OutputDir': 'aciss', 'DetectorType': 'ACIS-S'}

    @base.Marx
    def step_3(self):
        '''ACIS-I'''
        return {'SourceType': 'RAYFILE', 'RayFile': 'marx.output',
                'OutputDir': 'acisi', 'DetectorType': 'ACIS-I'}

    @base.Marx2fits
    def step_4(self):
        '''Turn into fits files

        We use the ``EXACT`` setting here to make the comparison simpler.
        The default EDSER (energy-dependent sub-pixel event repositioning)
        shifts photons of the same energy by a different amount for ACIS-S and
        ACIS-I, which would make it harder to compare the resulting images.
        '''
        return [ '--pixadj=EXACT',  '--pixadj=EXACT'], ['acisi', 'aciss'], ['i.fits', 's.fits']

    @base.Ciao
    def step_30(self):
        '''ds9 images of the PSF'''
        return ['''ds9 -width 800 -height 500 -log -cmap heat i.fits s.fits -pan to 4018 4141 physical -match frame wcs -saveimage {0} -exit'''.format(self.figpath(self.figures.keys()[0]))]


class SimputSource(base.MarxTest):
    pass


class UserSource(base.MarxTest):
    '''Run an example for a USER source.

    |marx| comes with several examples for user written source in C.
    These can be compiled as shared objects and dynamically linked into |marx|
    at run time.
    To test this, we copy one of the source files from the installed |marx|
    version and compile it with gcc. This particular case is not very useful,
    because |marx| already has a point source with the same properties
    build-in. The purpose of this test is only to have an automatic check that
    the dynamic linking works.
    '''

    title = 'Compiling a USER source'

    figures = OrderedDict([('ds9', {'alternative': 'A point source',
                                    'caption': '`ds9`_ shows that the distribution of source is indeed a point source.'})
                       ])

    @base.Python
    def step_1(self):
        '''compile USER code

        |marx| ships with a few examples of user sources. We pick one
        of them, copy them to the right directory and compile it with gcc.
        '''
        marxpath = self.conf.get('marx', 'path')
        src = os.path.join(marxpath,
                           'share', 'doc', 'marx', 'examples', 'user-source')
        for f in ['point.c', 'user.h']:
            shutil.copy(os.path.join(src, f),
                        os.path.join(self.basepath, f))

        jdmath_h = os.path.join(marxpath, 'include')
        jdmath_a = os.path.join(marxpath, 'lib', 'libjdmath.a')

        subprocess.call(['gcc', '-I' + jdmath_h, jdmath_a,
                         '-shared', 'point.c', '-o', 'point.so'])

    @base.Marx
    def step_2(self):
        '''run USER source'''
        return {'SourceType': 'USER',
                'UserSourceFile': os.path.join(self.basepath, 'point.so')}

    @base.Marx2fits
    def step_3(self):
        'turn into fits file'
        return '--pixadj=EDSER', 'point', 'point.fits'

    @base.Ciao
    def step_30(self):
        '''ds9 images of the PSF'''
        return ['''ds9 -width 800 -height 500 -log -cmap heat point.fits -pan to 4018 4141 physical -zoom 8 -saveimage {0} -exit'''.format(self.figpath(self.figures.keys()[0]))]
