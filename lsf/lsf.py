import os
import subprocess

from .. import utils

'''
These are tests that are run for one architecture only.
'''
marxdir = '/nfs/melkor/d1/guenther/marx/installed/dev/bin/'
outputdir = '/nfs/melkor/d1/guenther/marx/test/lsffix/'

parfile = os.path.join(marxdir, '..', 'share', 'marx', 'pfiles', 'marx.par')
energylist = [1.2]
energylist = np.logspace(np.log10(.2), 1., 50)

with utils.ChangeDir(outputdir):
    for e in energylist:
        name = 'energy{0:.2f}'.format(e)
        utils.call_marx(marxdir=marxdir, parfile=parfile, OutputDir=name,
                        MaxEnergy=e, MinEnergy=e, SpectrumType="FLAT",
                        # keep things simple so that sky and det coos are aligned
                        SourceRA=0., SourceDEC=0.,
                        RA_Nom=0., Dec_Nom=0., Roll_Nom=0.,
                        # Set number of photons
                        numRays=-10000000, ExposureTime=0,
                        # Now some parameters to speed this up
                        SourceFlux=0.1, dNumRays=100000)
        subprocess.check_call([os.path.join(marxdir, 'marx2fits'),
                               '--pixadj=RANDOMIZE',
                               name, name + '.fits'])



import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy import constants as c
dat = Table.read("/melkor/d1/guenther/marx/test/lsf/energy1.20RANDOMIZE.fits", hdu=1)

'''Yes, CIAO tools can do most of what I do below, but I want an independent
implementation if practical.
So, I will try to do simple things in python, but use CIAO tools for
complicated stuff.
Also: Cannot use CIAO and anaconda in the same shell.'''




# http://cxc.harvard.edu/proposer/POG/html/chap8.html
ang_meg = np.deg2rad(4.725)
ang_heg = np.deg2rad(-5.235)

# The number 400 depends on energy...

x0 = utils.find_centroid(dat['X'].data, 4000, 100)
y0 = utils.find_centroid(dat['Y'].data, 4100, 100)


disp3 = (dat['X'] - x0) * np.cos(ang_meg) - (dat['Y'] - y0) * np.sin(ang_meg)
cdisp3 = (dat['Y'] - y0) * np.cos(ang_meg) + (dat['X'] - x0) * np.sin(ang_meg)

p_heg = 2000.81 * u.Angstrom
p_meg = 4001.95 * u.Angstrom

pix_size_acis = 23.985 * u.micron

f = 8632.65 * u.mm


def estimate_order_pos(order, energy, p, pix_size=pix_size_acis):
    '''estimate ignoring that ACIS-S follows a curve etc.

    Return
    ------
    d : float
        distance to the 0th order in pixels
    '''
    wave = (energy * u.keV).to(u.m, equivalencies=u.spectral())
    return (order * f * wave / p  / pix_size).decompose().value


xpos = find_centroid(disp3, estimate_order_pos(0, 1.2, p_heg), 50)

ind3 = (np.abs(disp3 - xpos) < 20) & (np.abs(cdisp3) < 20)
hist3, edges3 = np.histogram(disp3[ind3], bins=1000)
h = hist3 / hist3.sum(dtype=np.float)
x = 0.5* (edges3[:-1] + edges3[1:])

from astropy.modeling import models, fitting
g1 = models.Gaussian1D(amplitude=.01, mean=xpos, stddev=0.5)
g2 = models.Gaussian1D(amplitude=.01, mean=xpos, stddev=1.5)
l1 = models.Lorentz1D(amplitude = .01, x_0=xpos, fwhm=0.5)
l2 = models.Lorentz1D(amplitude = .01, x_0=xpos, fwhm=1.5)
m = g1+g2+l1+l2

fit_t = fitting.LevMarLSQFitter()
t = fit_t(m, x, h)

plt.plot(x, h)
plt.plot(x, t(x))
for submodel in t:
    plt.plot(x, submodel(x))
