'''
In every |marx| simulation, one or more sources are placed at some sky position.
|marx| simulates photons coming from that position, traces them through the
mirror and gratings and finally places them on the chip. With a known
aspect solution, chip coordinates can then be transformed back to sky
coordinates. In general, this will not recover the exact sky position where a
photon started out. A big part of that is scatter in the mirrors, which blurs
the image (see :ref:`sect-tests.PSF` for tests of the PSF).
However, with a large number of photons, we can fit the average position which
should be close to the real sky position.

In real observations, other factors contribute, such as the finite
resolution of the detectors (|marx| usually takes that into account, but it can
be switched of through the ``--pixadj="EXACT"`` switch in :marxtool:`marx2fits`)
and the uncertainty of the aspect solution.

Within a single observation, positions will be less certain for fainter sources
(due to Poisson statistics) and for sources at a larger off-axis angles (due to the
larger PSF).
'''
import shutil
import subprocess
import os
from collections import OrderedDict
from marxtest import base
from marxtest.process_utils import marxpars_from_asol

title = 'Coordinates on the sky and the chip'

tests = ['ONC', 'RegularGrid', 'RegularGridHRCI']


class ONC(base.MarxTest):
    '''The `Orion Nebula Cluster (ONC) <http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=onc>`_
    is a dense star forming region with about 1600 X-ray sources observed
    in the COUP survey by
    `Getman et al (2005) <http://adsabs.harvard.edu/abs/2005ApJS..160..319G>`_ .
    We simulate this field with |marx| and then run a source detection to check
    how well we recover the input coordinates. This will depend on the number
    of counts detected and the position in the field.
    To simplify the simulation input, we assume that all sources have flat
    lightcurves and are
    monoenergetic at the observed mean energy (the energy matters because
    the effective area is energy dependent and so is the PSF).
    We write a short C code that reads an input coordiante list and generates
    the photons in this manner. We compile
    the code, and call it as a :ref:`sect-usersource`.
    '''

    title = 'Chandra Orion Ultradeep project'
    obsid = 3744

    figures = OrderedDict([('ds9', {'alternative': '',
                                    'caption': '`ds9`_ image of the observed data (left) and simulation (right). The sources detected in the simulation are overlayed. There are few cases where the read-out streak is identified as source or where two close sources are detected as one larger resolved source. The COUP catalog used as input is based on much longer merged observations and has been checked against optical and IR observations to remove such spurious detections.'}),
                           ('dist', {'alternative': 'Scatter plot with distance from aimpoint vs coordinate error in the fit.',
                                     'caption': 'Apart from a few outliers close to the aimpoint  (mostly confused sources, see above), the distribution of coordinate errors follows spreads out with increasing distance, i.e. size of the PSF.'})
    ])

    summary='For this field, we know the true input coordinates so we can check how well |marx| reproduces those. In the center of the field (about one armin) the coordiante error is less than the size of an ACIS pixel for all sources and the average error never grows much beyond 1 ACIS pixel even for far off-axis source. The upper envelope of the distribution of errors is approximate linear and reaches 1 arcsec at a distance of 200 arcsec. No strong correlation of coordiante error and count rate of the source is apparent, indicating that the dominant error is not just due to Poisson counting statistics.'

    @base.Python
    def step_2(self):
        '''Make input coordinate table

        Coordinates are relative to pointing direction in arcmin'''
        import os
        from astropy.table import Table
        from astropy.io import fits

        asolfile = self.get_data_file('asol')
        asol = fits.getheader(asolfile, 1)
        coup = Table.read(os.path.join(self.pkg_data, 'COUP.tsv'),
                          format='ascii.fast_tab')
        tab = Table()
        tab['RA'] = (coup['RAJ2000'] - asol['RA_NOM']) * 60
        tab['DEC'] = (coup['DEJ2000'] - asol['DEC_NOM']) * 60
        tab['weight'] = 10**(coup['Lt'] - 27)
        tab['energy'] = coup['<E>']
        tab.write('coup.marxin', format='ascii.no_header', overwrite=True)

    @base.CCode
    def step_5(self):
        '''C code for a grid of sources.

        (``user.h`` and ``jdmath.h`` are shipped with |marx|.)'''
        ccode=r'''
#include <stdio.h>
#include <stdlib.h>
#include <jdmath.h>
#include "user.h"

/* This user source implements many point sources via a file that
 * specifies the source positions and energies.  The current implementation
 * assumes the format:
 *  RA  Dec weight energy
 * Here RA, Dec specifiy the source position, weight specifies the strength
 * of the source in relation to the others.
 */
typedef struct
{
   double cosx, cosy, cosz;
   double weight;
   double energy;
}
Point_Source_Type;

static unsigned int Num_Points;
static Point_Source_Type *Point_Sources;
static unsigned int Max_Num_Points;

static char *do_realloc (char *p, unsigned int len)
{
   if (p == NULL)
     p = malloc (len);
   else
     p = realloc (p, len);

   if (p == NULL)
     fprintf (stderr, "Not enough memory\n");

   return p;
}

static void free_sources (void)
{
   if (Point_Sources == NULL)
     return;

   free ((char *) Point_Sources);
   Point_Sources = NULL;
}

static int add_source (double ra, double dec, double weight, double energy)
{
   Point_Source_Type *p;
   double cosx, cosy, cosz;

   /* Convert to God's units from arc-min */
   ra = ra * (PI/(180.0 * 60.0));
   dec = dec * (PI/(180.0 * 60.0));

   if (Max_Num_Points == Num_Points)
     {
        Max_Num_Points += 32;
        p = (Point_Source_Type *)do_realloc ((char *)Point_Sources, Max_Num_Points * sizeof (Point_Source_Type));
        if (p == NULL)
          {
             free_sources ();
             return -1;
          }
        Point_Sources = p;
     }

   p = Point_Sources + Num_Points;
   /* Note the the minus sign is to generate a vector pointing from the
    * source to the origin
    */
   p->cosx = -cos (dec) * cos (ra);
   p->cosy = -cos (dec) * sin(ra);
   p->cosz = -sin (dec);

   p->weight = weight;
   p->energy = energy;
   Num_Points += 1;

   return 0;
}

static void normalize_sources (void)
{
   double total;
   unsigned int i;

   total = 0;
   for (i = 0; i < Num_Points; i++)
     {
	Point_Sources[i].weight += total;
	total = Point_Sources[i].weight;
     }

   for (i = 0; i < Num_Points; i++)
     Point_Sources[i].weight /= total;

   /* Make sure no round-off error affects the weight of the last point */
   Point_Sources[Num_Points - 1].weight = 1.0;
}

int user_open_source (char **argv, int argc, double area,
		      double cosx, double cosy, double cosz)
{
   FILE *fp;
   char line[1024];
   char *file;
   unsigned int linenum;

   file = argv[0];
   if (file == NULL)
     {
	fprintf (stderr, "UserSource Model requires FILE as argument\n");
	return -1;
     }

   fp = fopen (file, "r");
   if (fp == NULL)
     {
	fprintf (stderr, "Unable to open %s\n", file);
	return -1;
     }

   linenum = 0;
   while (NULL != fgets (line, sizeof (line), fp))
     {
	double ra, dec, weight, energy;

	linenum++;
	if (4 != sscanf (line, "%lf %lf %lf %lf", &ra, &dec, &weight, &energy))
	  continue;

	if (weight <= 0.0)
	  {
	     fprintf (stderr, "weight on line %d of %s must be positive\n",
		      linenum, file);
	     free_sources ();
	     return -1;
	  }

	if (-1 == add_source (ra, dec, weight, energy))
	  {
	     fclose (fp);
	     return -1;
	  }
     }

   fclose (fp);
   if (Num_Points == 0)
     {
	fprintf (stderr, "%s contains no sources\n", file);
	return -1;
     }

   normalize_sources ();
   return 0;
}

void user_close_source (void)
{
   free_sources ();
}


int user_create_ray (double *delta_t, double *energy,
		     double *cosx, double *cosy, double *cosz)
{
   double r;
   Point_Source_Type *p;

   p = Point_Sources;

   r = JDMrandom ();
   while (r > p->weight)
     p++;

   *delta_t = -1.0;
   *energy = p->energy;
   *cosx = p->cosx;
   *cosy = p->cosy;
   *cosz = p->cosz;

   return 0;
}

int main (int a, char **b)
{
   (void) a;
   (void) b;
   return 1;
}
'''
        return 'pnts.c', ccode

    @base.Python
    def step_6(self):
        '''compile USER code

        |marx| ships with a few examples of user sources. We pick one
        of them, copy them to the right directory and compile it with gcc.
        '''
        marxpath = self.conf.get('marx', 'path')
        src = os.path.join(marxpath,
                           'share', 'doc', 'marx', 'examples', 'user-source')
        shutil.copy(os.path.join(src, 'user.h'),
                    os.path.join(self.basepath, 'user.h'))

        jdmath_h = os.path.join(marxpath, 'include')
        jdmath_a = os.path.join(marxpath, 'lib', 'libjdmath.a')

        subprocess.call(['gcc',
                         '-shared', 'pnts.c', '-o', 'pnts.so', '-fPIC',
                         '-I' + jdmath_h, jdmath_a])

    @base.Shell
    def step_7(self):
        '''Unzip fits file.

        MARX cannot read zipped fits files, so we need to unzip the .fits.gz asol
        files that we downloaded from the archive. On the other hand, `CIAO`_
        tools work on both zipped or unzipped files, so there is no need to
        unzip all of them, just the files that MARX reads as input.
        '''
        asol = self.get_data_file('asol')
        return [f'gunzip -f {asol}']

    @base.Marx
    def step_8(self):
        '''run marx USER source matching observation'''
        asol = self.get_data_file('asol')
        evt = self.get_data_file('evt2')
        pars = marxpars_from_asol(self.conf, asol, evt)
        pars['OutputDir'] = 'COUP'
        pars['SourceType'] = 'USER'
        pars['UserSourceFile'] = os.path.join(self.basepath, 'pnts.so')
        pars['UserSourceArgs'] = os.path.join(self.basepath, 'coup.marxin')

        return pars

    @base.Marx2fits
    def step_9(self):
        'turn into fits file'
        return '--pixadj=EDSER', 'COUP', 'COUP.fits'

    @base.Ciao
    def step_10(self):
        '''ds9 image of the PSF

        In the observation, the brightest sources are piled-up. We don't bother
        simulating this here, so we just set the scaling limits to bring out
        the fainter details and ignore the bright peaks.
        '''
        return ['''ds9 -log -cmap heat {0} COUP.fits -scale limits 0 2000 -frame 1 -regions command 'text 5:35:15 -5:22:09 # text=Observation font="helvetica 24"' -frame 2 -regions command 'text 5:35:15 -5:22:09 # text=MARX font="helvetica 24"' -region load src.fits -saveimage {1} -exit'''.format(self.get_data_file("evt2"), self.figpath(list(self.figures.keys())[0]))]

    @base.Ciao
    def step_11(self):
        '''Source detection'''
        out = ['dmcopy "COUP.fits[EVENTS][bin x=2500:5500:2,y=2500:5500:2]" im.fits  option=image clobber=yes',
             'mkpsfmap im.fits psf.map 1.4 ecf=0.5',
          'celldetect im.fits src.fits psffile=psf.map clobber=yes'
               ]
        return out

    @base.Python
    def step_15(self):
        '''Check position of detected sources'''
        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.table import Table
        from astropy.coordinates import SkyCoord
        from astropy.io import fits
        src = Table.read('src.fits')
        srcin = Table.read(os.path.join(self.pkg_data, 'COUP.tsv'),
                           format='ascii.fast_tab')

        src_co = SkyCoord(src['RA'], src['DEC'], unit='deg')
        srcin_co = SkyCoord(srcin['RAJ2000'], srcin['DEJ2000'], unit='deg')
        idx, d2d, d3d = src_co.match_to_catalog_sky(srcin_co)

        asolfile = self.get_data_file('asol')
        asol = fits.getheader(asolfile, 1)
        cen = SkyCoord(asol['RA_NOM'], asol['DEC_NOM'], unit='deg')
        d = cen.separation(src_co).arcsec

        fig = plt.figure()
        ax1 = plt.subplot(111)
        scat1 = ax1.scatter(d, d2d.arcsec, c=np.log10(src['NET_COUNTS']), lw=1)
        ax1.set_xlabel('distance from aimpoint [arcsec]')
        ax1.set_ylabel('coordinate error [arcsec]')
        ax1.set_xlim([0, 350])
        ax1.set_ylim([0, 2])
        cbar1 = fig.colorbar(scat1, ax=ax1)
        cbar1.set_label('log(net counts per source)')

        fig.savefig(self.figpath(list(self.figures.keys())[1]))


class RegularGrid(base.MarxTest):
    '''In this example we place a radial grid of sources on the sky. Each source
    emits an equal number of photons (exactly, no Poisson statistics) so that
    we can compare the accuracy of the position we recover. Note that the
    *detected* number of photons will be smaller for off-axis photons because
    of vignetting!

    We write a short C code that generates the photons in this manner, compile
    it, and call is as a :ref:`sect-usersource`.
    '''

    DetectorType = 'ACIS-I'

    title = 'Regular Grid (ACIS)'

    figures = OrderedDict([('ds9', {'alternative': 'Sources positioned like knots in a spider web.',
                                    'caption': '`ds9`_ image of the simulation. The size of the PSF increases further away from the aimpoint.'}),
                           ('hist', {'alternative': 'Plot is described in the caption.',
                                     'caption': '*left*: The error in the position (measured radially to the optical axis) increases with the distance to the optical axis. One part of this is just that the effective area and thus the number of counts decreases. There is also a systematic trend where sources at larger off-acis angle are systematically fitted too close to the center. Further investigation is necessary to check if this is a problem of |marx| related or :ciao:`celldetect`. In any case, the typical offset is below 0.2 arcsec, which is less then half a pixel in ACIS. *right*: Difference in position angle between input and fit. (Outliers beyond the plot range are not shown.)'})
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
        return ['''ds9 -width 500 -height 500 -log -cmap heat points.fits -pan to 4097 4097 physical -zoom 0.5 -bin factor 2 -grid -saveimage {0} -exit'''.format(self.figpath(list(self.figures.keys())[0]))]

    @base.Ciao
    def step_11(self):
        '''Source detection'''
        out = ['dmcopy "points.fits[EVENTS][bin x=3000:5100:2,y=3000:5100:2]" im.fits  option=image clobber=yes',
                 'mkpsfmap im.fits psf.map 1.4 ecf=0.5 clobber=yes',
          'celldetect im.fits src.fits psffile=psf.map clobber=yes'
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

        fig.savefig(self.figpath(list(self.figures.keys())[1]))


class RegularGridHRCI(RegularGrid):
    '''Same as above, but with HRC-I as a detector.

    The field-of-view for the HRC-I is larger for than for ACIS-I, but the PSF becomes
    very large at large off-axis angles and thus the positional uncertainty
    will be so large that a comparison to |marx| is no longer helpful to test
    the accuracy of the |marx| simulations.
    '''

    figures = OrderedDict([('ds9', {'alternative': 'Sources positioned like knots in a spider web. The image is very similar to the previous ACIS example.',
                                    'caption': '`ds9`_ image of the simulation. The size of the PSF increases further away from the aimpoint.'}),
                           ('hist', {'alternative': 'Plot is described in the caption.',
                                     'caption': 'See previous example. The same trends are visible with a slightly larger scatter.'})
    ])

    summary = 'In the central few arcmin the input position is typically recovered to better than 0.2 pixels for sources with a few hundred counts.'

    DetectorType = 'HRC-I'
    title = 'Regular grid (HRC)'

    @base.Ciao
    def step_10(self):
        '''ds9 image of the PSF'''
        return ['''ds9 -width 500 -height 500 -log -cmap heat points.fits -pan to 16392 16392 physical -bin factor 16 -grid -saveimage {0} -exit'''.format(self.figpath(list(self.figures.keys())[0]))]

    @base.Ciao
    def step_11(self):
        '''Source detection'''
        out = ['dmcopy "points.fits[EVENTS][bin x=8500:24500:8,y=8500:24500:8]" im.fits  option=image clobber=yes',
                 'mkpsfmap im.fits psf.map 1.4 ecf=0.5 clobber=yes',
          'celldetect im.fits src.fits psffile=psf.map clobber=yes'
               ]
        return out
