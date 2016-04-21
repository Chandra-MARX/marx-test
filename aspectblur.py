from copy import deepcopy
import os
import subprocess
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.table import Table

from utils import ChangeDir, download_chandra, call_marx, find_centroid

marxdir = '/nfs/melkor/d1/guenther/marx/installed/dev/bin/'
outputdir = '/nfs/melkor/d1/guenther/marx/test/aspectblur/'

parfile = os.path.join(marxdir, '..', 'share', 'marx', 'pfiles', 'marx.par')

marx_pars = {'marxdir': marxdir, 'parfile': parfile,
             'SpectrumType': "FLAT",
             # keep things simple so that sky and det coos are aligned
             'SourceRA': 0., 'SourceDEC': 0.,
             'RA_Nom': 0., 'Dec_Nom': 0., 'Roll_Nom': 0.,
             # Set number of photons
             'numRays': 0, 'ExposureTime': 0
             }

obsid = '13182'  # AR Lac  HRC-I

# this is an open cluster - not a point source!
obsid = '4469'   # tau Canis Major - ACIS-I

''' TYC
there is a little pile-up in the central pixel!

sherpa> fluxdensity = get_source_plot(lo=0.3,hi=3.0).y
sherpa> energy = get_source_plot(lo=0.3,hi=3.0).xhi
sherpa> save_arrays("Obsid15713_source_flux_marx.dat",[energy,fluxdensity],["keV","photons/s/cm**2/keV"],ascii=True)
'''
obsid = '15713'

with ChangeDir(outputdir):

    ### Uncomment the following line###
    download_chandra(obsid, obsid)

    asolfile = glob(os.path.join(obsid, 'primary', '*asol*'))[0]
    asol = Table.read(asolfile)

    ## Uncomment the following line ###
    # skyco = SkyCoord.from_name(asol.meta['OBJECT'].split('/')[0])
    # Set sensible default parameters
    marx_pars = {'marxdir': marxdir, 'parfile': parfile,
                 'SpectrumType': "FILE",
                 'SpectrumFile': "/melkor/d1/guenther/marx/test/marx-test/Obsid15713_source_flux_marx.dat",
                 'SourceFlux': -1,
                 # keep things simple so that sky and det coos are aligned
                 'SourceRA': skyco.ra.value,
                 'SourceDEC': skyco.dec.value,
                 'RA_Nom': asol.meta['RA_NOM'],
                 'Dec_Nom': asol.meta['DEC_NOM'],
                 'Roll_Nom': asol.meta['ROLL_NOM'],
                 'GratingType': 'NONE',
                 'ExposureTime': asol.meta['TSTOP'] - asol.meta['TSTART'],
                 'DitherModel': 'FILE',
                 'DitherFile': asolfile,
                 'TStart': asol.meta['TSTART'],
                 }

    for i in range(10):
        name = 'marxrun{0}'.format(i)
        mpars = deepcopy(marx_pars)
        mpars['OutputDir'] = name
        call_marx(**mpars)
        # subprocess.check_call([os.path.join(marxdir, 'marxpileup'),
        #                        'MarxOutputDir={0}'.format(name),
        #                        'FrameTime=0.4'])
        subprocess.check_call([os.path.join(marxdir, 'marx2fits'),
                               '--pixadj=EDSER',
                               os.path.join(name),
                               name + '.fits'])
        subprocess.check_call([os.path.join(marxdir, 'marx2fits'),
                               '--pixadj=RANDOMIZE',
                               name,
                               name + 'rand.fits'])

with ChangeDir(outputdir):

    evt2file = glob(os.path.join(obsid, 'primary', '*evt2*'))[0]
    obs = Table.read(evt2file)
    # In this energy range we have the most counts and we limited
    # the simulations to the same range.
    ind = (obs['energy'] > 300) & (obs['energy'] < 3000)
    # reduce background by cutting in the likely region
    indposx = (obs['x'] > 4000) & (obs['x'] < 4150)
    indposy = (obs['y'] > 4000) & (obs['y'] < 4150)
    obsfiltered = obs[ind & indposx & indposy]
    data = obsfiltered
    centx = find_centroid(data['x'], 4096, 10)
    centy = find_centroid(data['y'], 4074, 10)

    xy = np.empty((2, len(data)))
    xy[0, :] = data['x'] - centx
    xy[1, :] = data['y'] - centy
    r = np.linalg.norm(xy, axis=0)

    val, edges = np.histogram(r, range=[0, 5], bins=25)

    plt.plot(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'k', lw=3, label='Obs (EDSER)')
    #plt.plot(0.1 + edges[:-1], val, 'k', lw=3, label='Obs (EDSER)')


    evt2file = glob(os.path.join(obsid, 'randomize', '*evt2*'))[0]
    obs = Table.read(evt2file)
    # In this energy range we have the most counts and we limited
    # the simulations to the same range.
    ind = (obs['energy'] > 300) & (obs['energy'] < 3000)
    # reduce background by cutting in the likely region
    indposx = (obs['x'] > 4000) & (obs['x'] < 4150)
    indposy = (obs['y'] > 4000) & (obs['y'] < 4150)
    obsfiltered = obs[ind & indposx & indposy]
    data = obsfiltered
    centx = find_centroid(data['x'], 4096, 10)
    centy = find_centroid(data['y'], 4074, 10)

    xy = np.empty((2, len(data)))
    xy[0, :] = data['x'] - centx
    xy[1, :] = data['y'] - centy
    r = np.linalg.norm(xy, axis=0)

    val, edges = np.histogram(r, range=[0,5], bins=25)

    plt.plot(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'g', lw=3, label='Obs (RANDOMIZE)')

    for i in range(10):

        for opt, c in zip(['', 'rand'], 'rb'):
            name = 'marxrun{0}{1}'.format(i, opt)
            #  data = asol[randinds[i, :]]
            data = Table.read('{0}.fits'.format(name))
            # pile-up spectrum could have other energies
            ind = (data['ENERGY'] > 300) & (data['ENERGY'] < 3000)
            data = data[ind]
            centx = find_centroid(data['X'], 4096, 10)
            centy = find_centroid(data['Y'], 4074, 10)
            xy = np.empty((2, len(data)))
            xy[0, :] = data['X'] - centx
            xy[1, :] = data['Y'] - centy
            r = np.linalg.norm(xy, axis=0)
            val, edges = np.histogram(r, range=[0, 5], bins=25)
            if i == 0:
                plotargs = {'label': 'sim {0}'.format(opt), 'color': c}
            else:
                plotargs = {'color': c}
            plt.plot(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), **plotargs)
            # plt.plot(0.1 + edges[:-1], val, **plotargs)
