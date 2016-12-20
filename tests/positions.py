'''
In a |marx| simulation, a source is placed in some sky position. |marx| simulates
photons coming from that position and places them on the chip. With a known
aspect solution, chip coordinates can then be transformed back to sky
coordinates. In general, this will not recover the exact sky position where a
photons started out. A big part of that is scatter in the mirrors, which blurs
the image (see :ref:`sect-tests.PSF` for tests of the PSF that |marx|
simulates).
However, with a large number of photons, we can fit the average position which
should be close to the real sky position.

In real observations, other factors contribute, such as the finite
resolution of the detectors (|marx| usually takes that into account, but it can
be switched of through the ``--pixadj="EXACT"`` switch in :marxtool:`marx2fits`)
and the uncertainty of the aspect solution.

Within a single observation, positions will be less certain for fainter sources
due to Poisson statistics and for source at a larger off-axis angles due to the
larger PSF.
'''
import shutil
import subprocess
import os
from collections import OrderedDict
from marxtest import base

title = 'Coordinates on the sky and the chip'

tests = ['RegularGrid', 'RegularGridHRCI']


class RegularGrid(base.MarxTest):
    '''In this example we place a radial grid of sources on the sky. Each source
    emits an equal number of photons (exactly, no Poisson statistics) so that
    we can compare the accuracy of the position we recover. Note that the
    *detected* number of photons will be smaller for off-axis photons because
    of vignetting!

    We write a short C code that generates the photons in this manner, compile
    it, and call is as a ``USER`` source (see :ref:`sect-usersource`).
    '''

    DetectorType = 'ACIS-I'

    title = 'Regular Grid (ACIS)'

    figures = OrderedDict([('ds9', {'alternative': 'Sources positioned like knots in a spider web.',
                                    'caption': '`ds9`_ image of the simulation. The size of the PSF increases further away from the aimpoint.'}),
                           ('hist', {'alternative': 'Plot is described in the caption.',
                                     'caption': '*left*: The error in the position (measured radially to the optical axis) increases with the distance to the optical axis. One part of this is just that the effective area and thus the number of counts decreases. There is also a systematic trend where sources at larger off-acis angle are systematically fitted too close to the center. Further investigation is necessary to check if this is |marx| related to just due to the use of :ciao:`celldetect`. In any case, the typical offset is below 0.2 arcsec, which is less then half a pixel in ACIS. *right*: Difference in position angle between input and fit. (Outliers beyond the plot range are not shown.)'})
    ])

    summary = 'The input position is typically recovered to much better than one pixel for sources with a few hundred counts. There is a small systematic trend that needs to be studied further.'

    @base.CCode
    def step_5(self):
        '''C code for a grid of sources.

        (``user.h`` is shipped with |marx|.)'''
        ccode='''
#include <stdio.h>
#include <math.h>
#include "user.h"

static double Source_CosX;
static double Source_CosY;
static double Source_CosZ;

int user_open_source (char **argv, int argc, double area,
		      double cosx, double cosy, double cosz)
{
   return 0;
}

void user_close_source (void)
{
}

static double To_Radians = (M_PI / 180.0 / 3600.0);
#define ARC_SECONDS_PER_CELL 50
#define ANGULAR_STEPS 16

int user_create_ray (double *delta_t, double *energy,
		     double *cosx, double *cosy, double *cosz)
{
   static int last_i = 0;
   static int last_j = 0;
   double theta, phi;
   double cos_theta, sin_theta;

   if (last_j == ANGULAR_STEPS){
        last_j = 0;
        last_i++;
   }
   if (last_i == 20) last_i = 0;

   theta = To_Radians * last_i * ARC_SECONDS_PER_CELL;
   phi = (10. /180 * M_PI) + last_j * 2 * M_PI / ANGULAR_STEPS;

   sin_theta = sin(theta);

   *cosx = -cos (theta);
   *cosy = sin_theta * cos (phi);
   *cosz = sin_theta * sin (phi);

   *delta_t = -1.0;
   *energy = -1.0;

   if (last_i ==0){
     last_i++;
        }
   else {
     last_j++;
        }

   return 0;
}

int main (int a, char **b)
{
   (void) a;
   (void) b;
   return 1;
}'''
        return 'radialgrid.c', ccode

    @base.Python
    def step_6(self):
        '''compile USER code'''
        marxpath = self.conf.get('marx', 'path')
        src = os.path.join(marxpath, 'share', 'doc', 'marx', 'examples',
                           'user-source', 'user.h')
        shutil.copy(os.path.join(src),
                    os.path.join(self.basepath, 'user.h'))

        subprocess.call(['gcc', '-lm', '-fPIC',
                         '-shared', 'radialgrid.c', '-o', 'radialgrid.so'])

    @base.Marx
    def step_7(self):
        '''run USER source'''
        return {'SourceType': 'USER', 'OutputDir': 'points',
                'GratingType': 'NONE',
                'SourceRA': 90., 'SourceDEC': 0.,
                'RA_Nom': 90., 'Dec_Nom': 0, 'Roll_Nom': 0,
                'DetectorType': self.DetectorType,
                'UserSourceFile': os.path.join(self.basepath, 'radialgrid.so'),
                'NumRays': -100000, 'ExposureTime': 0}

    @base.Marx2fits
    def step_8(self):
        'turn into fits file'
        return '--pixadj=EDSER', 'points', 'points.fits'

    @base.Ciao
    def step_10(self):
        '''ds9 image of the PSF'''
        return ['''ds9 -width 500 -height 500 -log -cmap heat points.fits -pan to 4097 4097 physical -zoom 0.5 -bin factor 2 -grid -saveimage {0} -exit'''.format(self.figpath(self.figures.keys()[0]))]

    @base.Ciao
    def step_11(self):
        '''Source detection'''
        out = ['dmcopy "points.fits[EVENTS][bin x=3000:5100:2,y=3000:5100:2]" im.fits  option=image clobber=yes',
          'celldetect im.fits src.fits clobber=yes'
               ]
        return out

    @base.Python
    def step_15(self):
        '''Check position of detected sources'''
        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.table import Table
        from astropy.coordinates import SkyCoord
        src = Table.read('src.fits')
        # Find distance from input position.
        src['RA_INPUT'] = src['RA'] - (src['RA'] // (360./16.)) * (360./16.) - 10.
        # Problem: Might expect source at 1.0,
        # but measure at 0.99. In this case distance to next lower source
        # is 0.99. Thus shift input by 0.005 (about 50 arcsec / 2)
        # before integer devision
        src['DEC_INPUT'] = src['DEC'] - ((0.005 + src['DEC']) // (50./3600.)) * (50./3600.)

        cen = SkyCoord(90., 0, unit='deg')
        det = SkyCoord(src['RA'], src['DEC'], unit='deg')
        d = cen.separation(det).arcsec
        d_err = np.mod(d + 10, 50.) - 10

        ang = cen.position_angle(det).degree
        # Subtract offset that we placed in the C code to avoid 0./360. ambiguity
        # Step width is 360./16 = 22.5 deg
        # Offset is 10 deg. Complement we find here is 12.5 deg.
        ang = ang - 12.5

        ang_err = np.mod(ang + 2, 360. / 16.) - 2
        ind = d > 10

        fig = plt.figure(figsize=(8, 4))
        ax1 = plt.subplot(121)
        scat1 = ax1.scatter(d, d_err, c=src['NET_COUNTS'], lw=1)
        ax1.set_xlabel('distance [arcsec]')
        ax1.set_ylabel('distance error [arcsec]')
        ax1.set_xlim([-10, 620])
        ax1.set_ylim([-1, 0.5])
        ax2 = plt.subplot(122)
        scat2 = ax2.scatter(ang, ang_err, c=src['NET_COUNTS'], lw=1)
        ax2.set_xlabel('pos ang [deg]')
        ax2.set_ylabel('pos ang error [deg]')
        ax2.set_xlim([-5, 350])
        ax2.set_ylim([-0.3, 0.3])
        cbar2 = fig.colorbar(scat2, ax=ax2)
        cbar2.set_label('net counts per source')

        fig.savefig(self.figpath(self.figures.keys()[1]))


class RegularGridHRCI(RegularGrid):
    '''Same as above, but with HRC-I as a detector.

    The field-of-view for the HRC-I is larger for than for ACIS-I, but the PSF becomes
    very large at large off-axis angles and thus the positional uncertainty
    will be so large that a comparison to |marx| is no longer helpful to test
    the accuracy of the |marx| simulations.
    '''

    figures = OrderedDict([('ds9', {'alternative': 'Sources positioned like knots in a spider web.',
                                    'caption': 'See previous example'}),
                           ('hist', {'alternative': 'Plot is described in the caption.',
                                     'caption': 'See previous example. The same trends are visible with a slightly larger scatter.'})
    ])

    summary = 'The input position is typically recovered to better than 0.2 pixels for sources with a few hundred counts.'

    DetectorType = 'HRC-I'
    title = 'Regular grod (HRC)'

    @base.Ciao
    def step_10(self):
        '''ds9 image of the PSF'''
        return ['''ds9 -width 500 -height 500 -log -cmap heat points.fits -pan to 16392 16392 physical -bin factor 16 -grid -saveimage {0} -exit'''.format(self.figpath(self.figures.keys()[0]))]

    @base.Ciao
    def step_11(self):
        '''Source detection'''
        out = ['dmcopy "points.fits[EVENTS][bin x=8500:24500:8,y=8500:24500:8]" im.fits  option=image clobber=yes',
          'celldetect im.fits src.fits clobber=yes'
               ]
        return out
