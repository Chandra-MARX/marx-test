from copy import deepcopy
import os
import subprocess
from glob import glob

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from utils import ChangeDir, download_chandra, call_marx, find_centroid

marxdir = '/nfs/melkor/d1/guenther/marx/installed/dev/bin/'
outputdir = '/nfs/melkor/d1/guenther/marx/test/chipxy/'

parfile = os.path.join(marxdir, '..', 'share', 'marx', 'pfiles', 'marx.par')

marx_pars = {'marxdir': marxdir, 'parfile': parfile,
             'SpectrumType': "FLAT",
             # keep things simple so that sky and det coos are aligned
             'SourceRA': 0., 'SourceDEC': 0.,
             'RA_Nom': 0., 'Dec_Nom': 0., 'Roll_Nom': 0.,
             # Set number of photons
             'numRays': 0, 'ExposureTime': 0
             }
obsid = '13683'
obsid = '2074'

with ChangeDir(outputdir):

    download_chandra(obsid, obsid, ['evt2', 'asol'])

    asolfile = glob(os.path.join(obsid, 'primary', '*asol*'))[0]
    asol = Table.read(asolfile)
    # object in 2074 is 'NGC7619/NGC7626'
    skyco = SkyCoord.from_name(asol.meta['OBJECT'].split('/')[0])
    # Set sensible default parameters
    marx_pars = {'marxdir': marxdir, 'parfile': parfile,
                 'SpectrumType': "FLAT",
                 # keep things simple so that sky and det coos are aligned
                 'SourceRA': skyco.ra.value,
                 'SourceDEC': skyco.dec.value,
                 'RA_Nom': asol.meta['RA_NOM'],
                 'Dec_Nom': asol.meta['DEC_NOM'],
                 'Roll_Nom': asol.meta['ROLL_NOM'],
                 'GratingType': 'NONE',
                 'ExposureTime': asol.meta['TSTOP'] - asol.meta['TSTART']
                 }

    name = 'asol'
    mpars = deepcopy(marx_pars)
    mpars['OutputDir'] = name
    mpars['DitherModel'] = 'FILE'
    mpars['DitherFile'] = asolfile
    mpars['TStart'] = asol.meta['TSTART']
    call_marx(**mpars)
    subprocess.check_call([os.path.join(marxdir, 'marx2fits'),
                           '--pixadj=RANDOMIZE',
                           name, name + '.fits'])

    asol = Table.read('asol.fits')
    evt2file = glob(os.path.join(obsid, 'primary', '*evt2*'))[0]
    obs = Table.read(evt2file)
    # Simulation has a lot more counts than observation and no background
    # So comparison metric is not very easy.
    ind = (obs['energy'] > 300) & (obs['energy'] < 3000)
    # reduce background by cutting in the likely region
    indposx = (obs['chipx'] > 100) & (obs['chipx'] < 300)
    indposy = (obs['chipy'] > 400) & (obs['chipy'] < 600)
    obsfiltered = obs[ind & indposx & indposy]

    '''

    The numerical test bins the image up in CHIPX and CHIPY and tests if
    a cross correlation shows no shift greater than +- 2 pixel.
    '''
    status = 'pass'

    for c in ['chipx', 'chipy']:
        cent = find_centroid(asol[c.upper()], 200, 1e4)
        r = [cent - 50, cent + 50]
        hs, x = np.histogram(asol[c.upper()], range=r, bins=50)
        ho, x = np.histogram(obsfiltered[c], range=r, bins=50)
        hs = 1.0 * hs / hs.sum()
        ho = 1.0 * ho / ho.sum()

        corr = np.correlate(hs, ho, mode="same")
        if np.max(corr) < 0.2 or np.abs(np.argmax(corr) - 25) > 2:
            status = 'fail'
        elif np.abs(np.argmax(corr) - 25) > 1:
            status = 'warning'

    if (plotmode == 'all') or status != 'pass':
        import matplotlib.pyplot as plt
        plt.clf()
        n = len(asol)
        plt.plot(asol['CHIPX'] + np.random.rand(n) - 0.5,
                 asol['CHIPY'] + np.random.rand(n) - 0.5, 'b.', label='sim')
        n = len(obsfiltered)
        plt.plot(obsfiltered['chipx'] + np.random.rand(n) - 0.5,
                 obsfiltered['chipy'] + np.random.rand(n) - 0.5, 'm.',
                 label='obs')
        plt.legend()
        plt.xlabel('CHIPX')
        plt.ylabel('CHIPY')
        plt.title('Do CHIP coordinates agree in obs and asol based sim?')
        plt.xlim([150, 260])
        plt.ylim([450, 560])
        plt.plot(os.path.join(plotpath, 'chipxy.png'))
