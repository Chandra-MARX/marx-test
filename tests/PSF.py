'''
The `point-spread function (PSF) for Chandra <http://cxc.harvard.edu/ciao/PSFs/psf_central.html>`_ describes how the light from a point source is spread over a larger area on the detector. Several effects contribute to this, e.g. the uncertainty in the pointing, the fact that the detectors are flat, while the focal plane of the mirror is curved (specifically for large off-axis angles) and the pixalization of data on detector read-out.

The following tests compare |marx| simulations, `SAOTrace`_ simulations, and data to look at different aspects of the Chandra PSF.
'''
from __future__ import division

import os
import numpy as np
from collections import OrderedDict
from astropy.table import Table
from astropy.io import fits
import astropy
import astropy.stats

from marxtest import base
from marxtest.utils import colname_case as cc

from marxtest.process_utils import (marxpars_from_asol,
                                    spectrum_from_fluxcorrection)
# on import, this registers the plotting scale "power"
from marxtest import plot_utils

tests = ['ACISSPSF', 'ACISIPSF', 'HRCIPSF', 'OffAxisPSF',
         'Wings',
         'CompMARXSAOTraceenergies', 'CompMARXSAOTraceoffaxis']

title = 'Point Spread Function (PSF)'


def plot_ecf(ax, files, filterargs, bgfilterargs):
    '''
    Parameters
    ----------
    ax : matplotlib axes
        Axes where the plot will be placed
    file : list
        List of three strings with filename for observations,
        Marx simulation, and SAOTRace + marx simulation
    filterargs : dict
        See ``filter_events`` for details
    bgfilterargs : dict
        Select a background region (same size, no automatic scaling)

    Returns
    -------
    err_marx : float
        Flux error if the MARX simulated PSF is assumed to be right.
        This is calculated as follows: We find the radius that encircles 90% of
        all counts in the MARX simulation. Then, we extract the observed counts
        using that radius. If the simulation is correct, that radius should
        contain 90% of the observed counts, too. If, e.g. the simulated radius
        is too small, we may extract only 80 % of the counts. The ratio between
        the two would be 8/9=0.88, meaning that all fluxes extracted using
        this radius are 12 % too small.
    err_sao: float
        Same for the simulation that used SAOTrace + MARX.
    '''
    # In this energy range we have the most counts and we limited
    # the simulations to the same range.
    obs = filter_events(Table.read(files[0]), **filterargs)
    bkg = filter_events(Table.read(files[0]), **bgfilterargs)
    r_obs = radial_distribution(obs['x'], obs['y'])
    r_obs.sort()
    ax.plot(r_obs, np.linspace(0, 1 + len(bkg) / len(r_obs), len(r_obs)),
            label='Observation')

    simmarx = filter_events(Table.read(files[1]), **filterargs)
    r_marx = radial_distribution(simmarx['X'], simmarx['Y'])
    r_marx.sort()
    ax.plot(r_marx, np.linspace(0, 1, len(r_marx)), label='MARX')

    simsao = filter_events(Table.read(files[2]), **filterargs)
    r_sao = radial_distribution(simsao['X'], simsao['Y'])
    r_sao.sort()
    ax.plot(r_sao, np.linspace(0, 1, len(r_sao)), label='SAOTrace + MARX')

    ax.set_xscale('power', power=0.5)
    ax.set_xlabel('radius [pixel]')
    ax.set_ylabel('enclosed count fraction')
    ax.legend(loc='lower right')
    ax.set_xticks([.1, .4, .7, 1, 2, 3, 4, 5])
    ax.set_xlim([0.1, 5])
    ax.set_ylim([0, 1.])

    psf90_marx = np.percentile(r_marx, 90)
    ecf_at_that_rad = np.argmin(np.abs(r_obs - psf90_marx)) / len(r_obs)
    err_marx = ecf_at_that_rad / 0.9 - 1

    psf90_sao = np.percentile(r_sao, 90)
    ecf_at_that_rad = np.argmin(np.abs(r_obs - psf90_sao)) / len(r_obs)
    err_sao = ecf_at_that_rad / 0.9 - 1

    return err_marx, err_sao


def filter_events(evt, circle=None, energy=None, time=None):
    '''A simple event filter function.

    The parameters offer several filters, it set to ``None`` a filter is not
    applied.
    This simple function does not duplicate the entire :ciao:`dm` syntax, only
    a few simple filters useful for plotting in pure-python functions.

    Parameters
    ----------
    evt : `astropy.table.Table`
        Event table
    circle : tuple or ``None``
        Tuple of the form (x, y, r) where (x, y) is the center of
        a circle and r the radius. All number are in detector coordinates.
    energy : tuple or ``None``
        Tuple of the form (lower, upper) with all energies in eV.
    time : tuple or ``None``
        Tuple of the form (star, end) with all times in seconds.

    Returns
    -------
    evt_filt : ~astropy.table.Table`
        Filtered event list
    '''
    if circle is not None:
        d = np.sqrt((cc(evt, 'x') - circle[0])**2 + (cc(evt, 'y') - circle[1])**2)
        indcirc = d < circle[2]
    else:
        indcirc = np.ones(len(evt), dtype=bool)

    if energy is not None:
        en = cc(evt, 'energy')
        inden = (en > energy[0]) & (en < energy[1])
    else:
        inden = np.ones_like(indcirc)

    if time is not None:
        t = cc(evt, 'time')
        indt = (t > time[0]) & (t < time[1])
    else:
        indt = np.ones_like(indcirc)

    return evt[indcirc & inden & indt]


def radial_distribution(x, y):
    '''Calculate the radial distance from the mean position for events.

    Parameters
    ----------
    x, y : np.array
        x and y coordinates of events

    Returns
    -------
    r : np.array
        radial distance from mean position of the events.
    '''
    centx = x - astropy.stats.sigma_clipped_stats(x, sigma=1.5, iters=20)[0]
    centy = y - astropy.stats.sigma_clipped_stats(y, sigma=1.5, iters=20)[0]

    xy = np.vstack([centx, centy])
    return np.linalg.norm(xy, axis=0)


class HRCIPSF(base.MarxTest):
    '''Same as above for an HRC-I observations of AR Lac.
    '''

    title = 'On-axis PSF for an HRC-I observation'
    obsid = 13182
    download_all = True

    figures = OrderedDict([('ECF', {'alternative': 'TBD',
                                    'caption': 'Enclosed count fraction for observation and simulations.'})
                       ])

    summary = '''In contrast to the ACIS simulations, in the HRC-I the observed PSF is in fact narrower than all simulations in the range shown here. The difference is most prominent between about one and two HRC resolution elements where the |marx| simulation is about midway between the true observed PSF shape and the `SAOTrace`_  + |marx| predicted shape.'''

    source = {'x': 16405,
              'y': 16500,
              'r': 15,
              'bg_x': 16043,
              'bg_y': 16700
    }
    '''Source position and extraction radius.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''

    expresults = [{'name': 'marx90', 'title': 'Flux err (MARX)',
                   'description': 'Flux error if the MARX simulated PSF is assumed to be right. This is calculated as follows: We find the radius that encircles 90% of all counts in the MARX simulation. Then, we extract the observed counts using that radius. If the simulation is correct, that radius should contain 90% of the observed counts, too. If, e.g. the simulated radius is too small, we may extract only 80 % of the counts. The ratio between the two would be 8/9=0.88, meaning that all fluxes extracted using this radius are 12 % too small.',
                   'value': 0},
                  {'name': 'saotrace90', 'title': 'Fluxerr(SAOTr+MARX)',
                   'description': 'Same, but using a model of SAOTrace + MARX',
                   'value': 1}]

    @property
    def source_reg(self):
        return "circle({0}, {1}, {2})".format(self.source['x'],
                                              self.source['y'],
                                              self.source['r'])

    @base.Python
    def step_0(self):
        '''Copy input spectra

        Since the HRC has no intrisinc energy resolution, the source spectrum
        can not be taken from the observed data. Instead, we fitted a model
        to a grating spectrum of the same source and saved that model
        spectrum to a file. AR Lac is known to be time variable and to show
        stellar flares, so this spectrum is only an approximation, but it
        is the best we can do here.

        Spectrum files are part of the marx-test distribution.
        Copy them from the source code into the right directory.
        '''
        import shutil
        shutil.copy(os.path.join(self.pkg_data, 'ARLac_input_spec_marx.tbl'),
                    os.path.join(self.basepath, 'input_spec_marx.tbl'))
        shutil.copy(os.path.join(self.pkg_data, 'ARLac_input_spec_saotrace.rdb'),
                    os.path.join(self.basepath, 'input_spec_saotrace.rdb'))

    @base.Marx
    def step_1(self):
        '''Set marx parameters appropriate for observation'''
        asol = self.get_data_file('asol')
        evt = self.get_data_file('evt2')
        pars = marxpars_from_asol(self.conf, asol, evt)
        pars['SpectrumType'] = 'FILE'
        pars['SpectrumFile'] = 'input_spec_marx.tbl'
        pars['OutputDir'] = 'marx_only'
        pars['SourceFlux'] = -1
        return pars

    @base.Marx2fits
    def step_2(self):
        '''No EDSER is available for HRC data'''
        return ('--pixadj=NONE', 'marx_only', 'marx_only.fits')

    @base.SAOTraceLua
    def step_3(self):
        '''`SAOTrace`_ input matching observation'''
        asol = self.get_data_file('asol')
        asolh = fits.getheader(asol, 1)
        evt = fits.getheader(self.get_data_file('evt2'), 1)
        return '''
ra_pnt = {rapnt}
dec_pnt = {decpnt}
roll_pnt = {rollpnt}

dither_asol_chandra{{ file = "{asol}",
                     ra = ra_pnt, dec = dec_pnt, roll = roll_pnt }}

point{{ position = {{ ra = {ra},
          dec = {dec},
          ra_aimpt = ra_pnt,
          dec_aimpt = dec_pnt,
       }},
       spectrum = {{ {{ file = "input_spec_saotrace.rdb",
                      units = "photons/s/cm2",
                      scale = 1,
                      format = "rdb",
                      emin = "ENERG_LO",
                      emax = "ENERG_HI",
                      flux = "FLUX"}} }}
    }}
'''.format(asol=asol, rapnt=asolh['RA_NOM'], decpnt=asolh['DEC_NOM'],
           rollpnt=asolh['ROLL_NOM'],
           ra=evt['RA_TARG'], dec=evt['DEC_TARG'])

    @base.SAOTrace
    def step_4(self):
        '''
        CXO time numbers are large and round-off error can appear
        which make `SAOTrace`_ fail.
        Therefore, shorten all times by about 0.01 sec to make sure.
        '''
        asol = self.get_data_file('asol')
        asolh = fits.getheader(asol, 1)

        limit = asolh['TSTOP'] - asolh['TSTART']
        limit = limit - 0.02
        return ['trace-nest tag=saotrace srcpars=saotrace_source.lua tstart={tstart}  limit={limit} limit_type=sec'.format(tstart=asolh['TSTART'] + 0.01, limit=limit)]

    @base.Marx
    def step_5(self):
        '''Run |marx| with `SAOTrace`_ ray file as input'''
        asol = self.get_data_file('asol')
        evt = self.get_data_file('evt2')
        pars = marxpars_from_asol(self.conf, asol, evt)
        pars['SpectrumType'] = 'FILE'
        pars['SpectrumFile'] = 'input_spec_marx.tbl'
        pars['OutputDir'] = 'marx_saotrace'
        pars['SourceType'] = 'SAOSAC'
        pars['SAOSACFile'] = 'saotrace.fits'
        pars['SourceFlux'] = -1
        return pars

    @base.Marx2fits
    def step_6(self):
        '''Same settings as the marx2fits run above'''
        return ('--pixadj=NONE', 'marx_saotrace', 'marx_saotrace.fits')

    @base.Python
    def step_20(self):
        '''Extract radial count distribution'''
        import matplotlib.pyplot as plt
        filterargs = {'circle': (self.source['x'], self.source['y'],
                                 5 * self.source['r'])
                      }

        bgfilterargs =  {'circle': (self.source['bg_x'], self.source['bg_y'],
                                 5 * self.source['r'])
                        }

        evt2file = self.get_data_file('evt2')
        files = [evt2file, 'marx_only.fits', 'marx_saotrace.fits']

        fig = plt.figure()
        ax = fig.add_subplot(111)

        err_marx, err_sao = plot_ecf(ax, files, filterargs, bgfilterargs)
        ax.set_xlim([0, 20.])

        fig.savefig(self.figpath(self.figures.keys()[0]))

        self.save_test_result('marx90', err_marx)
        self.save_test_result('saotrace90', err_sao)


class ACISSPSF(HRCIPSF):
    '''
    The PSF depends on many things, some of which are common to all observations
    like the shape of the mirror, and some are due to detector effects.
    For ACIS detectors, the
    `sub-pixel event repositioning (EDSER) <http://cxc.harvard.edu/ciao/why/acissubpix.html>`_
    can improve
    the quality of an image by repositioning events based on the event grade.
    This correction depends on the type pf chip (FI or BI). This test compares
    the simulation of a point source on a BI ACIS-S chip to an observation.
    The observed object is TYC 8241 2652 1, a young star, and was observed in
    1/8 sub-array mode to reduce pile-up. The pile-up fraction in the data is
    about 5% in the brightest pixel.
    '''

    title = 'On-axis PSF on an ACIS-BI chip'
    obsid = 15713
    download_all = True

    figures = OrderedDict([('ECF', {'alternative': 'TBD',
                                    'caption': 'Enclosed count fraction for observation and simulations. The simulations are run with a count number similar to the (fairly short) observation and thus there is some wiggling due to Poisson noise.'})
                       ])

    summary = '''For BI chips (here ACIS-S3) the |marx| simulated PSF is too narrow. This effect is most pronounced at small radii around 1 pixel. At larger radii the |marx| simulation comes closer to the observed distribution. Running the simulation with a combination of `SAOTrace`_ and |marx| gives PSF distributions that are closer to the observed numbers. However, in the range around 1 pix, the simulations are still too wide. Depending on the way sub-pixel information is handeled, there is a notable difference in the size of the effect. Using energy-dependent sub-pixel event repositioning (EDSER) requires not only a good mirror model, but also a realistic treatment of the flight grades assigned on the detector. This is where the difference between simulation and observations is largest.

The figures show that there are at least two factors contributing to the difference: The mirror model and problems in the simulation of the EDSER algorithm.

The put those numbers into perspective: If I used the |marx| simulation to determine the size of the extraction region that encloses 80% of all source counts, I would find a radius close to 0.9 pixel. However, in reality, such a region will only contain about 65% of all flux, causing me to underestimate the total X-ray flux in this source by 15%. (The exact number depends on the spectrum of the source in question, but for most on-axis sources they will be similar.)
'''

    source = {'x': 4096,
              'y': 4074,
              'r': 4,
              'bg_x': 4044,
              'bg_y': 4123}
    '''Source position and extraction radius.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''

    @base.Ciao
    def step_0(self):
        '''Obtain spectrum from observed data'''
        asolfile = self.get_data_file('asol')
        evtfile = self.get_data_file('evt2')
        commands = spectrum_from_fluxcorrection(self.conf,
                                                asolfile, evtfile,
                                                self.source['x'],
                                                self.source['y'],
                                                self.source_reg)
        return commands

    @base.Marx2fits
    def step_2(self):
        '''Use the EDSER subpixel algorithm'''
        return ('--pixadj=EDSER', 'marx_only', 'marx_only.fits')

    @base.Marx2fits
    def step_6(self):
        '''Same setting as above for comparison'''
        return ('--pixadj=EDSER', 'marx_saotrace', 'marx_saotrace.fits')

    @base.Marx2fits
    def step_10(self):
        '''USE RANDOMIZE as an alternative to EDSER'''
        return ('--pixadj=RANDOMIZE', 'marx_only', 'marx_only_rand.fits')

    @base.Marx2fits
    def step_11(self):
        '''Again, same setting as above'''
        return ('--pixadj=RANDOMIZE', 'marx_saotrace', 'marx_saotrace_rand.fits')

    @base.Ciao
    def step_12(self):
        '''Reprocess the observation with RANDOMIZE'''
        evt2 = self.get_data_file('evt2')
        asol = self.get_data_file('asol')
        return ['acis_process_events infile={evt} outfile=obs_rand.fits acaofffile={asol} doevtgrade=no calculate_pi=no pix_adj=RANDOMIZE clobber=yes'.format(evt=evt2, asol=asol)]

    @base.Python
    def step_20(self):
        '''Compare radial event distributions'''
        import matplotlib.pyplot as plt

        filterargs = {'energy': [300, 3000],
                      'circle': (self.source['x'], self.source['y'],
                                 5 * self.source['r'])
                     }
        bgfilterargs = {'energy': [300, 3000],
                        'circle': (self.source['bg_x'], self.source['bg_y'],
                                 5 * self.source['r'])
        }

        evt2file = self.get_data_file('evt2')

        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        for ax, files, title in zip([ax1, ax2],
                                    [[evt2file, 'marx_only.fits', 'marx_saotrace.fits'],
                                     ['obs_rand.fits', 'marx_only_rand.fits', 'marx_saotrace_rand.fits']],
                                    ['EDSER', 'RANDOMIZE']):
            err_marx, err_sao = plot_ecf(ax, files, filterargs, bgfilterargs)

            ax.set_title(title)

            if title == 'EDSER':
                self.save_test_result('marx90', err_marx)
                self.save_test_result('saotrace90', err_sao)

        fig.savefig(self.figpath(self.figures.keys()[0]))


class ACISIPSF(ACISSPSF):
    '''Same as above for a FI chip. The target of this observations is
    tau Canis Majores.
    '''

    title = 'On-axis PSF on an ACIS-FI chip'
    obsid = 4469

    figures = OrderedDict([('ECF', {'alternative': 'TBD',
                                    'caption': 'Enclosed count fraction for observation and simulations.'})
                       ])

    summary = '''As for the BI chip, |marx| simulations indicate a PSF that is narrower than the observed distribution. Using `SAOTrace`_ as a mirror model, gives better results but they still differ significantly from the observed distribution.'''

    source = {'x': 4140,
              'y': 4102,
              'r': 4,
              'bg_x': 4182,
              'bg_y': 4127}


class OffAxisPSF(HRCIPSF):
    '''
    For an off-axis source the differences between detector pixel size and different
    event repositioning algorithms is not important any longer. Instead, the size of
    the PSF is dominated by optics errors because the detector plane deviates from
    the curved focal plane and because a Wolter type I optic has fundamental limitations
    for off-axis sources.

    The PSF is much larger and it is no longer round. As this test shows, the shadows of
    the support struts become visible. The code to run this test is very similar to
    :ref:`sect-ex-simobs`, where the individual steps are explained in greater detail.
    '''

    title = 'Off-axis PSF'
    obsid = 1068
    download_all = True

    figures = OrderedDict([('ds9', {'alternative': 'Three very similar PSFs.',
                                    'caption': '`ds9`_ image of the PSF in the observation (top left), the simulation using only |marx| (top right), and the simulation using `SAOTrace`_ to trace the mirror and |marx| as the instrument model (bottom left)'})
                            ])

    summary = '''The structure of the PSFs is very similar, emphasizing how good both mirror
models are. On closer inspection, there is a small shadow just above and to the
right of the point where the support strut shadows meet. This feature is a
little smaller in |marx| than in `SAOTrace`_ or the real data due to the
simplification that the |marx| mirror model makes.'''

    source = {'x': 5246,
              'y': 6890,
              'r': 200}

    @base.Ciao
    def step_0(self):
        '''Get spectrum from flux correction

        Use one of the asol files in only.
        Good enough for this example and simpler than stack syntax
        '''
        asolfile = self.get_data_file('asol')
        evtfile = self.get_data_file('evt2')
        commands = spectrum_from_fluxcorrection(self.conf,
                                                asolfile, evtfile,
                                                self.source['x'],
                                                self.source['y'],
                                                self.source_reg)
        return commands

    @base.Ciao
    def step_20(self):
        '''ds9 images of the PSF'''
        return ['''ds9 -width 800 -height 500 -log -cmap heat {0} marx_only.fits marx_saotrace.fits -pan to 5256 6890 physical -bin about 5256 6890 -match frame wcs -match bin -frame 1 -regions command 'text 5:39:27.987 -69:43:52.31 # text=Observation font="helvetica 24"' -frame 2 -regions command 'text 5:39:27.987 -69:43:52.31 # text="only MARX" font="helvetica 24"' -frame 3 -regions command 'text 5:39:27.987 -69:43:52.31 # text=SAOTrace font="helvetica 24"' -regions command 'text 5:39:26.691 -69:44:09.93 # text="+ MARX" font="helvetica 24"' -saveimage {1} -exit'''.format(self.get_data_file('evt2'), self.figpath(self.figures.keys()[0]))]


class Wings(HRCIPSF):
    '''The Chandra PSF falls in two regimes: The shape of the core is set by quasi-specular
    reflection from low-frequency surface errors. This component is strongly peaked and dominates
    in the innermost few arcseconds (for and on-axis source). The out part of the PSF has the shape of a powerlaw with an spectral index close to two and is caused by micro-roughness on the mirror surfaces. The `Proposer's Observatory Guide <http://cxc.harvard.edu/proposer/POG/html/chap4.html#tth_sEc4.2.3>`_ and `a memo by T.J. Gaetz <http://cxc.harvard.edu/cal/Hrma/rsrc/Publish/Optics/PSFWings/wing_analysis_rev1b.pdf>`_ discuss this in a lot more detail. The memo is based on a very detailed analysis of a deep Her-X1 observation. In this test we replicate that observation in |marx| and `SAOTrace`_ simulations.

    In the memo, T.J. Gaetz goes to great lengths to separate the effect of the intrinsic Chandra PSF on the data from other influences such as the ACIS contamination. When we compare data and |marx| simulation below, we chose a somewhat simpler approach. We perform steps to mitigate the impact of other astrophysical sources (which are not part of the |marx| simulation) and background (which is not part of the |marx| simulation). On the other hand, we just compare the observed count rates with no correction for the ACIS contamination, dithering near chip edges etc. because those effects are included in both the observed and the simulated data.

    Read the code for this test to see all details of the data extraction, but here are the most
    important points:

    - The source is heavily piled-up, thus the source spectrum is taken from the read-out streak.
    - Because of pile-up in the data, all fitting is restricted to radii above 15 arcsec.
    '''

    title = 'Wings of the Chandra PSF'
    obsid = 3662
    download_all = True

    figures = OrderedDict([('PowerLaws', {'alternative': 'Four plots show that the flux in the scattering halo drops as a powerlaw with increasing distance from the source. See text for a comparison of data and observations.',
                                          'caption': "The flux in the scattering halo drops as a powerlaw with increasing distance from the source. The plots show the data (with purely statistical uncertainties) and powerlaw fits. Due to the background that is part of the observations but not the simulated data, the powerlaw flattens out for the observed data in some cases. We add a constant to the model to account for this. `SAOTrace`_ simulations produce a fairly smooth powerlaw similar to the observations. The pure |marx| simulations have much larger deviations due to |marx|'s simplified mirror model."}),
                           ('PowerSlope', {'alternative': 'See caption and summary text for a description.',
                                    'caption': 'The surface brightness falls off with increasing distance from the source position as a powerlaw. The exponent of this powerlaw depends on the photon energy as shown in the figure. If spectra are extracted in the scattering halo the observed spectrum will thus change with distance from the source. For the observed data, the low energy end (below about 1 keV) has to be taken with a grain of salt. We applied a correction for the ACIS contamination in the data reduction, but the distribution of the contaminant over the detector area is not well known, leading to systematic uncertainties. For energies above 5 keV, the outer mirror shell does not contribute much effective area any longer. This shell has the roughest surface, which explains why we observe steeper powerlaw slopes for higher eneergies. Both the |marx| and the `SAOTrace` mirror model also have energy dependent exponents.'})
                       ])

    summary = '''Both |marx| and `SAOTrace`_ implement a mirror model that has wide scattering wings which are also seen in actual Chandra data. These wings can be traced out to at least 8 arcmin. On the other hand, these scattering wings are very weak (the flux per surface area drops by about four orders of magnitute between 10 arcsec and 500 arcsec for an on-axis source) and thus they are not not relevant in any but the very brightest objects. The `SAOTrace` mirror model produces a smoother distribution that is more similar to the observed data tan |marx|. In both mirror models the slope depends on the energy. `SAOTrace`_ reproduces the shape of the observed data very well, but the exponent is significantly too steep everywhere. Thus, spectra extracted in the halo will be compatible with observed data, but the flux in the simulation will be too low. On the other hand, |marx| produces more realistic exponents above 1 keV, so it will match the observed surface brightness better, but not the spectral properties. '''

    source = {'x': 4174,
              'y': 4043}
    '''Source position.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''

    source_reg = 'rotbox(3922.6168,3556.9146,12,512,332.71498)'
    '''Source is heavily piled-up so we extract the source spectrum from the read-out streak.'''

    expresults = [{'name': 'marx90', 'title': 'Flux err (MARX)',
                   'description': 'Flux error if the MARX simulated PSF is assumed to be right. This is calculated as follows: We find the radius that encircles 90% of all counts in the MARX simulation. Then, we extract the observed counts using that radius. If the simulation is correct, that radius should contain 90% of the observed counts, too. If, e.g. the simulated radius is too small, we may extract only 80 % of the counts. The ratio between the two would be 8/9=0.88, meaning that all fluxes extracted using this radius are 12 % too small.',
                   'value': 0},
                  {'name': 'saotrace90', 'title': 'Fluxerr(SAOTr+MARX)',
                   'description': 'Same, but using a model of SAOTrace + MARX',
                   'value': 1}]

    energy_bins = np.array([.4, .5, .6, .7, .9, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.,
                    3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.])
    '''The analysis makes narrow band count rate images.
    Bands are defined here. If they are too wide, then the effective area and spectrum
    within one band cannot be taken as constant any longer, but if they are too narrow
    the test will need a very long time to run.
    '''

    @base.Ciao
    def step_0(self):
        '''Obtain spectrum from observed data'''
        asolfile = self.get_data_file('asol')
        evtfile = self.get_data_file('evt2')
        commands = spectrum_from_fluxcorrection(self.conf,
                                                asolfile, evtfile,
                                                self.source['x'],
                                                self.source['y'],
                                                self.source_reg,
                                                # Rough estimate of normalization
                                                # Read-out streak region is 512 pix long
                                                0.00004*512/3.2)
        return commands

    @base.Marx2fits
    def step_2(self):
        '''Use the EDSER subpixel algorithm'''
        return ('--pixadj=EDSER', 'marx_only', 'marx_only.fits')

    @base.Marx2fits
    def step_6(self):
        '''Same setting as above for comparison'''
        return ('--pixadj=EDSER', 'marx_saotrace', 'marx_saotrace.fits')

    @base.Ciao
    def step_14(self):
        '''Make image and run cell detect to find exclusion regions for extraction.'''
        evtfile = self.get_data_file('evt2')
        asolfile = self.get_data_file('asol')
        return ['dmcopy "{0}[bin x=3700:4300:1,y=3700:4300:1][opt type=i4]" im.fits option=image clobber=yes'.format(evtfile),
               'celldetect im.fits im_src.fits clobber=yes',
               'asphist {0} asphist_7.fits evtfile="{1}[ccd_id=7]" clobber=yes'.format(asolfile, evtfile)
               ]

    @base.Ciao
    def step_15(self):
        '''Create exposure maps

        Also, make exposure map for S-3 chip to correct for dither etc. Mirror effective
        area is not included in that calculation, because all the photons we care for come
        from the same source, not from different positions on the sky, thus the mirror
        effective area is the constant for the source we look at.

        The instrument map depends strongly on the energy. We will thus pass in the
        source spectrum. This gives us one map appropriately scaled for the whole spectrum.
        However, in later analysis, we want to split the data in narrower bands.
        When looking at that accurately, we need to have one map for each band.

        Also, the observed data contains bad pixels, while the marx simulated data does not.
        This means that we need to make a second set of maps with no bad pixels in them.
        '''
        evtfile = self.get_data_file('evt2')
        out = []
        for i in range(len(self.energy_bins) - 1):
            e = self.energy_bins[i]
            e_mid = 0.5 * (self.energy_bins[i] + self.energy_bins[i + 1])
            for name, bpix in zip(['', '_marx'], ['', ';BPMASK=0x00000']):
                out.append('mkinstmap instmap{4}{1}.fits monoenergy={3} pixelgrid="1:1024:#1024,1:1024:#1024" obsfile="{0}[EVENTS]" dafile=NONE detsubsys="ACIS-S3{2}" mirror="HRMA;AREA=1" grating=NONE maskfile=NONE spectrumfil=NONE clobber=yes'.format(evtfile, name, bpix, e_mid, '_{0:3.1f}'.format(e)))
                out.append('mkexpmap asphist_7.fits expmap{1}{0}.fits instmap{1}{0}.fits xygrid="2975:4400:#512,3250:4650:#512" clobber=yes useavgaspect=no'.format(name, '_{0:3.1f}'.format(e)))

        return out

    @base.Python
    def step_16(self):
        '''Clean bogus sources from src file and write regions

        celldetect always detects plenty of sources on the read-out streak and
        around the piled-up regions.
        For processing speed and clarity we want to remove those bogus detections
        from the source list.
        '''
        import numpy as np
        from astropy.table import Table
        src = Table.read(os.path.join(self.basepath, 'im_src.fits'), hdu=1)
        # Get sources within 50 pix of Her X-1. We don't look there anyway because of pile-up
        src['d_herx1'] = np.sqrt((src['X'] - self.source['x'])**2 +
                                 (src['Y'] - self.source['y'])**2)
        d_src = src['d_herx1'] < 50.
        # Get sources on the read-out streak
        # 2 points on streak
        p1 = np.array([3772., 3264.])
        p2 = np.array([4251., 4193.])
        # n is normal to read-out streak
        n = np.array([1., - (p2[0] - p1[0]) / (p2[1] - p1[1])])
        n = n / np.linalg.norm(n)
        x = np.vstack([src['X'], src['Y']])
        d_line = np.abs(np.dot(n, (p1 - x.T).T)) < 2.
        src = src[~d_src & ~d_line]

        # box for readout streak
        streak = 'rotbox(4010.0285,3725.3416,18.240189,1085.0354,332.74569)'
        # Assemble regions for extraction:
        n_rings = 50
        r_rings = np.logspace(np.log10(7.), np.log10(500.), n_rings + 1) / 0.492
        with open(os.path.join(self.basepath, 'annuli.stk'), 'w') as f:
            for i in range(n_rings):
                reg = 'annulus({0}, {1}, {2}, {3})'.format(self.source['x'], self.source['y'],
                                                           r_rings[i], r_rings[i + 1])
                reg += '-' + streak
                # find sources that overlap this ring
                ind1 = src['d_herx1'] + src['R'][:, 0] > r_rings[i]
                ind2 = src['d_herx1'] - src['R'][:, 0] < r_rings[i + 1]
                src_intersect = src[ind1 & ind2]
                for s in src_intersect:
                    reg += '-ellipse({0}, {1}, {2}, {3}, {4})'.format(s['X'], s['Y'],
                                                                        s['R'][0], s['R'][1],
                                                                        s['ROTANG'])
                f.write(reg + '\n')

    @base.Ciao
    def step_17(self):
        '''Extract data in energy bins'''
        evtfile = self.get_data_file('evt2')

        out = []
        for efile, name, outname in zip([evtfile, 'marx_only.fits', 'marx_saotrace.fits'],
                                        ['', '_marx', '_marx'],
                                        ['_obs', '_marx', '_saotrace']):
            for i in range(len(self.energy_bins) - 1):
                e = self.energy_bins[i]
                e_upper = self.energy_bins[i + 1]
                # Energy in fits file in in eV not, keV, thus a factor of 1000 here
                out.append('dmextract "{efile}[energy={elower}:{eupper}][bin sky=@annuli.stk]" profile{estr}{outname}.fits exp=expmap{estr}{name}.fits opt=generic2 clobber=yes'.format(efile=efile, elower=e*1000, eupper=e_upper*1000, estr='_{0:3.1f}'.format(e), name=name, outname=outname))
        return out

    @base.Python
    def step_20(self):
        '''Analyze and plot results - Plot 1'''
        import numpy as np
        from matplotlib import pyplot as plt
        from astropy.table import Table
        from astropy.modeling import models
        import saba
        powconst = models.PowerLaw1D + models.Const1D
        fit_g = saba.SherpaFitter(statistic='chi2', optimizer='levmar')

        fig = plt.figure(figsize=(10, 8))
        energy_index = [1, 8, 18, 24]
        fullname = {'obs': 'Observation', 'marx': 'MARX',
                    'saotrace': 'SAOTrace + MARX'}
        for i in range(len(energy_index)):
            e = self.energy_bins[energy_index[i]]
            e_upper = self.energy_bins[energy_index[i] + 1]
            ax = fig.add_subplot(2, 2, i + 1)
            for j, name in enumerate(['obs', 'marx', 'saotrace']):
                tab = Table.read('profile_{0:3.1f}_{1}.fits'.format(e, name))
                r = np.mean(tab['R'], axis=1)
                r_arcsec = r * 0.492
                indfit = r_arcsec > 15.
                plots = ax.errorbar(r_arcsec, tab['SUR_FLUX'],
                                    yerr=tab['SUR_FLUX_ERR'],
                                    label='__no_legend__')

                if name == 'obs':
                    mymod = powconst(amplitude_1=tab['SUR_FLUX'][-1],
                                     bounds={'amplitude_1': (0, None)})
                else:
                    mymod = models.PowerLaw1D()
                # Multiply data values with 1e6 for numerical stability
                g = fit_g(mymod, r[indfit], tab['SUR_FLUX'][indfit] * 1e6,
                          err=tab['SUR_FLUX_ERR'][indfit] * 1e6)
                ax.plot(r_arcsec, g(r) / 1e6, label=fullname[name],
                        color=plots[0].get_color())
                ax.set_title('{0:3.1f} - {1:3.1f} keV'.format(e, e_upper))
                ax.loglog()
                ax.set_xlim(6, 550)
                if i in [0, 2]:
                    ax.set_ylabel('flux per area')
                if i in [2, 3]:
                    ax.set_xlabel('radius [arcsec]')

        ax.legend(fontsize='small', loc='lower left')
        #fig.subplots_adjust(top=.95, right=.99, wspace=.2, hspace=.25)
        fig.savefig(self.figpath(self.figures.keys()[0]))

    @base.Python
    def step_21(self):
        '''Analyze and plot results - Plot 2'''
        import numpy as np
        from matplotlib import pyplot as plt
        from astropy.table import Table
        from astropy.modeling import models, fitting
        import saba

        powconst = models.PowerLaw1D + models.Const1D
        fit_g = saba.SherpaFitter(statistic='chi2', optimizer='levmar')
        slopes = np.zeros((len(self.energy_bins) - 1, 3))
        for i in range(len(self.energy_bins) - 1):
            e = self.energy_bins[i]
            for j, name in enumerate(['obs', 'marx', 'saotrace']):
                tab = Table.read('profile_{0:3.1f}_{1}.fits'.format(e, name))
                r = np.mean(tab['R'], axis=1)
                r_arcsec = r * 0.492
                indfit = r_arcsec > 15.
                if name == 'obs':
                    mymod = powconst(amplitude_1=tab['SUR_FLUX'][-1],
                                     bounds={'amplitude_1': (0, None)})
                else:
                    mymod = models.PowerLaw1D()
                # Multiply data values with 1e6 for numerical stability
                g = fit_g(mymod, r[indfit], tab['SUR_FLUX'][indfit] * 1e6,
                          err=tab['SUR_FLUX_ERR'][indfit] * 1e6)
                if name == 'obs':
                    slopes[i, j] = g.alpha_0.value
                else:
                    slopes[i, j] = g.alpha.value
                print g

        fig = plt.figure()
        ax = fig.add_subplot(111)
        e_mid = 0.5 * (self.energy_bins[:-1] + self.energy_bins[1:])
        for i, name in enumerate(['Observation', 'MARX', 'SAOTrace + MARX']):
            ax.plot(e_mid, slopes[:, i], label=name)
            ax.legend()
            ax.set_xlabel('energy [keV]')
            ax.set_ylabel('slope of powerlaw')
        fig.savefig(self.figpath(self.figures.keys()[1]))


class CompMARXSAOTraceenergies(base.MarxTest):
    '''The gold standard to test the fidelity of |marx| and `SAOTrace`_
    obviously is to compare simulations to observations.
    However, it is also instructive to look at a few idealized cases with no
    observational counterpart so we can simulate high fidelity PSFs with a
    large number of counts without worrying about background or pile-up.
    The `Chandra Proposers Observatory Guide <http://cxc.harvard.edu/proposer/POG/html/chap4.html#tth_sEc4.2.3>`_
    contains a long section on the PSF and encircled energy based on `SAOTrace`_
    simulations. Some of those simulations are repeated here to compare them
    to pure |marx| simulations.
    '''

    title = 'On-axis PSF at different energies'

    summary = "|marx| and `SAOTrace`_ simulations predict a very similar PSF shape, but for most energies the `SAOTrace`_ model predicts a slightly broader PSF."

    figures = OrderedDict([('PSFimages', {'alternative': 'Gallery of PSf images. See caption for a description.',
                                          'caption': 'PSF images for different discrete energies for pure |marx| and `SAOTrace`_ + |marx| simulations. The color scale is linear, but the absolute scaling is different for different images, because the effective area and thus the number of detected photons is lower at higher energies. Contour lines mark flux levels at 30%, 60%, and 90% of the peak flux level. At higher energies, the PSF becomes asymmetric, because mirror pair 6, which is most important at high energies, is slightly tilted with respect to the nominal aimpoint. At the same time, the scatter, and thus the size of the PSF, increase at higher energies.'}),
                           ('ECF', {'alternative': 'The |marx| PSF is generally narrower than the  `SAOTrace`_ + |marx| PSF. The agreement is best for energies around 2-4 keV.',
                                    'caption': 'A different way to present the width of the PSF is the encircled count fraction. **solid line**: |marx| only, **dotted lines**: `SAOTrace`_ + |marx|. For each photon, the radial distance is calculated from the nominal source position. Because the center is offset for hard photons, the PSF appears wider at those energies.'})])

    parameter = [0.25, 0.5, 1, 2, 3, 4, 6, 8]

    def __init__(self, *args):
        super(CompMARXSAOTraceenergies, self).__init__(*args)
        self.marxnames = ['marx{0}'.format(e) for e in self.parameter]
        self.saonames = ['sao{0}'.format(e) for e in self.parameter]

    @base.Ciao
    def step_1(self):
        '''Set up default marx.par file'''
        com = ['cp {0} marx.par'.format(self.conf.get('marx', 'marxparfile'))]
        for s in ['GratingType=NONE', 'DetectorType=HRC-I',
                  'DitherModel=INTERNAL',
                  'ExposureTime=10000', 'SourceFlux=0.2',
                  'SourceRA=0.', 'SourceDEC=0.',
                  'RA_Nom=0.', 'Dec_Nom=0.', 'Roll_Nom=0.']:
            com.append('pset marx.par {0}'.format(s))
        return com

    @base.Marx
    def step_2(self):
        '''run marx for different energies'''
        marxpars = []
        for e, n in zip(self.parameter, self.marxnames):
            marxpars.append({'MinEnergy': e, 'MaxEnergy': e, 'OutputDir': n})
        return marxpars

    @base.Marxasp
    def step_3(self):
        '''Generate asol file

        Since all simulations use the same pointing and exposure time, it is
        enough to run :marxtool:`marxasp` once.
        '''
        return {'MarxDir': self.marxnames[0]}

    @base.Marx2fits
    def step_4(self):
        '''Fits files from |marx| runs'''
        options = ['--pixadj=NONE'] * len(self.parameter)
        fitsnames = [n + '.fits' for n in self.marxnames]
        return options, self.marxnames, fitsnames

    @base.SAOTrace
    def step_10(self):
        '''Now use `SAOTrace`_'''
        asolh = fits.getheader(os.path.join(self.basepath,
                                            'marx_asol1.fits'), 1)

        limit = asolh['TSTOP'] - asolh['TSTART']
        limit = limit - 0.02
        com = []
        for t, e in zip(self.saonames, self.parameter):
            # The following line has a ' in a " delimited string which itself
            # is placed inside a ' delimited string like this: '"\'aaa\'"'
            com.append('trace-nest tag={tag} srcpars="point{{ position = {{ ra = 0., dec = 0., ra_aimpt=0., dec_aimpt=0. }}, spectrum = {{ {{ {energy}, 0.2 }} }} }} roll(0) dither_asol_marx{{ file = \'marx_asol1.fits\', ra = 0., dec = 0., roll = 0. }}" tstart={tstart}  limit={limit} limit_type=sec'.format(tag=t, tstart=asolh['TSTART'] + 0.01, energy=e, limit=limit))

        return com

    @base.Marx
    def step_11(self):
        '''run marx for all SAOTrace runs'''
        marxpars = []
        for s, n in zip(self.saonames, self.marxnames):
            marxpars.append({'SourceType': 'SAOSAC', 'SAOSACFile': s + '.fits',
                             'OutputDir': 'sao' + n,
                             'DitherModel': 'FILE',
                             'DitherFile': 'marx_asol1.fits'})
        return marxpars

    @base.Marx2fits
    def step_12(self):
        '''Fits files from |marx| + `SAOTrace`_ runs'''
        options = ['--pixadj=NONE'] * len(self.parameter)
        dirnames = ['sao' + n for n in self.marxnames]
        fitsnames = ['sao' + n + '.fits' for n in self.marxnames]
        return options, dirnames, fitsnames

    @base.Python
    def step_20(self):
        '''Image gallery of PSFs'''
        import numpy as np
        from matplotlib import pyplot as plt
        from mpl_toolkits.axes_grid1 import AxesGrid
        from matplotlib.ticker import StrMethodFormatter
        from astropy.table import Table

        fig = plt.figure(figsize=(10, 4))

        grid = AxesGrid(fig, 111,  # similar to subplot(142)
                        nrows_ncols=(2, len(self.parameter)),
                        axes_pad=0.0,
                        share_all=True,
                        label_mode="L",
                        cbar_location="top",
                        cbar_mode="single")

        for i, e in enumerate(self.parameter):

            for prog in ['marx', 'saomarx']:
                tab = Table.read(os.path.join(self.basepath,
                                              '{0}{1}.fits'.format(prog, e)),
                                 hdu=1)
                # The old question: Is an index the center of a pixel of the corner?
                # Differs by 0.5...
                d_ra = (tab['X'] - tab.meta['TCRPX9'] - 0.5) * tab.meta['TCDLT9'] * 3600.
                # And do I count from the left or right (RA is reversed on the sky)?
                d_dec = (tab['Y'] - tab.meta['TCRPX10'] + 0.5) * tab.meta['TCDLT10'] * 3600.

                if prog == 'marx':
                    offset = 0
                else:
                    offset = len(self.parameter)
                im = grid[i + offset].hist2d(d_ra, d_dec,
                                             range=[[-1, 1], [-1, 1]],
                                             bins=[20, 20],
                                             cmap=plt.get_cmap('OrRd'))
                levels = np.max(im[0]) * np.array([0.3, 0.6, 0.9])
                grid[i + offset].contour(0.5 * (im[1][1:] + im[1][:-1]),
                                         0.5 * (im[2][1:] + im[2][:-1]),
                                         # for some reason, the data is
                                         # ordered for other orientation
                                         im[0].T, levels=levels, colors='k',
                                         origin=im[3].origin)
                grid[i + offset].grid()
                if prog == 'marx':
                    grid[i].text(-.5, .5, '{0} keV'.format(e))
                else:
                    grid[i].xaxis.set_major_formatter(StrMethodFormatter('{x:2.1g}'))

        grid[0].yaxis.set_major_formatter(StrMethodFormatter('{x:2.1g}'))
        grid[len(self.parameter)].yaxis.set_major_formatter(StrMethodFormatter('{x:2.1g}'))
        grid.cbar_axes[0].colorbar(im[3])
        for cax in grid.cbar_axes:
            cax.toggle_label(False)

        # This affects all axes as share_all = True.
        grid.axes_llc.set_xticks([-.5, 0, .5])
        grid.axes_llc.set_yticks([-.5, 0, .5])

        grid[0].set_ylabel('Marx')
        grid[offset].set_ylabel('Marx +\nSAOTrace')
        grid[int(1.5 * offset)].set_xlabel('arcsec (measured from nominal source position)')

        fig.subplots_adjust(top=1, bottom=0, left=0.1, right=0.99)
        fig.savefig(self.figpath(self.figures.keys()[0]))

    @base.Python
    def step_21(self):
        '''Plots of radial distribution'''
        import os
        import numpy as np
        from matplotlib import pyplot as plt
        from astropy.table import Table

        fig = plt.figure()
        axecf = fig.add_subplot(111)
        color = plt.cm.jet_r(np.linspace(0, 1, len(self.parameter)))
        for i, e in enumerate(self.parameter):

            for prog in ['marx', 'saomarx']:
                tab = Table.read(os.path.join(self.basepath,
                                              '{0}{1}.fits'.format(prog, e)),
                                 hdu=1)
                # The old question: Is an index the center of a pixel of the corner?
                # Differs by 0.5...
                d_ra = (tab['X'] - tab.meta['TCRPX9'] - 0.5) * tab.meta['TCDLT9'] * 3600.
                # Count from the left or right (RA is reversed on the sky)?
                d_dec = (tab['Y'] - tab.meta['TCRPX10'] + 0.5) * tab.meta['TCDLT10'] * 3600.

                # The simulation is set up to make this simple: RA=DEC=0
                # so cos(DEC) = 1 and we can approximate with Euklidian distance
                r = np.linalg.norm(np.vstack([d_ra, d_dec]), axis=0)
                val, edges = np.histogram(r, range=[0, 5], bins=50)
                bin_mid_marx = 0.5 * (edges[:-1] + edges[1:])
                ecf_marx = 1.0 * val.cumsum() / val.sum()
                if prog == 'marx':
                    line, = axecf.plot(bin_mid_marx, ecf_marx, color=color[i],
                                       label='{0} keV'.format(e))
                else:
                    axecf.plot(bin_mid_marx, ecf_marx, color=line.get_color(), ls=':')
        axecf.set_xscale('power', power=0.5)
        axecf.legend(loc='lower right')
        axecf.set_ylabel('encircled count fraction')
        axecf.set_xlabel('radius [arcsec]')
        axecf.set_xticks([0, .1, .2, .4, .6, .8, 1, 2, 3, 4, 5])
        fig.savefig(self.figpath('ECF'))


class CompMARXSAOTraceoffaxis(CompMARXSAOTraceenergies):
    '''Compare |marx| and `SAOTrace`_ + |marx| simulations for different
    off-axis angles. For simplicity, this is done here for a single energy.
    We pick 4 keV because it is roughly the center of the Chandra band.
    '''

    title = 'Simulated off-axis PSF'

    summary = "|marx| and |marx| + `SAOTrace`_ simulations show essentially identical off-axis behavior in this test."

    figures = OrderedDict([('ECF', {'alternative': 'TBD',
                                    'caption': 'Radius of encircled count fraction for different off-axis angles. **solid line**: |marx| only, **dotted lines**: `SAOTrace`_ + |marx|. The values shown include the effect of the HRC-I positional uncertainty and the uncertainty in the aspect solution.'})])

    parameter = [0, .5, 1, 2, 3, 5, 7, 10, 15]

    def __init__(self, *args):
        super(CompMARXSAOTraceenergies, self).__init__(*args)
        self.marxnames = ['marx{0}'.format(e) for e in self.parameter]
        self.saonames = ['sao{0}'.format(e) for e in self.parameter]

    @base.Ciao
    def step_1(self):
        '''Set up default marx.par file'''
        com = ['cp {0} marx.par'.format(self.conf.get('marx', 'marxparfile'))]
        for s in ['GratingType=NONE', 'DetectorType=HRC-I',
                  'DitherModel=INTERNAL',
                  'ExposureTime=10000', 'SourceFlux=0.2',
                  'SourceRA=0.', 'SourceDEC=0.',
                  'RA_Nom=0.', 'Dec_Nom=0.', 'Roll_Nom=0.',
                  'MinEnergy=4', 'MaxEnergy=4']:
            com.append('pset marx.par {0}'.format(s))
        return com

    @base.Marx
    def step_2(self):
        '''run marx for different energies'''
        marxpars = []
        for p, n in zip(self.parameter, self.marxnames):
            marxpars.append({'SourceRA': p / 60., 'OutputDir': n})
        return marxpars

    @base.SAOTrace
    def step_10(self):
        '''Now use `SAOTrace`_'''
        asolh = fits.getheader(os.path.join(self.basepath,
                                            'marx_asol1.fits'), 1)

        limit = asolh['TSTOP'] - asolh['TSTART']
        limit = limit - 0.02
        com = []
        for t, e in zip(self.saonames, self.parameter):
            # The following line has a ' in a " delimited string which itself
            # is placed inside a ' delimited string like this: '"\'aaa\'"'
            com.append('trace-nest tag={tag} srcpars="point{{ position = {{ ra = {ra}, dec = 0., ra_aimpt=0., dec_aimpt=0. }}, spectrum = {{ {{ 4., 0.2 }} }} }} roll(0) dither_asol_marx{{ file = \'marx_asol1.fits\', ra = 0., dec = 0., roll = 0. }}" tstart={tstart}  limit={limit} limit_type=sec'.format(tag=t, tstart=asolh['TSTART'] + 0.01, ra=e / 60., limit=limit))

        return com

    @base.Python
    def step_20(self):
        '''No image gallery needed'''
        pass

    @base.Python
    def step_21(self):
        '''Plots of radial distribution'''
        import os
        import numpy as np
        from matplotlib import pyplot as plt
        from astropy.table import Table

        ecf = [0.5, 0.7, 0.9]
        fig = plt.figure()
        axecf = fig.add_subplot(111)
        out = np.zeros((len(self.parameter), 2, len(ecf)))

        for i, e in enumerate(self.parameter):

            for prog in ['marx', 'saomarx']:
                tab = Table.read(os.path.join(self.basepath,
                                              '{0}{1}.fits'.format(prog, e)),
                                 hdu=1)
                # The old question: Is an index the center of a pixel of the corner?
                # Differs by 0.5...
                d_ra = (tab['X'] - tab.meta['TCRPX9'] - 0.5) * tab.meta['TCDLT9'] * 3600.
                # Count from the left or right (RA is reversed on the sky)?
                d_dec = (tab['Y'] - tab.meta['TCRPX10'] + 0.5) * tab.meta['TCDLT10'] * 3600.

                # The simulation is set up to make this simple: RA=DEC=0
                # so cos(DEC) = 1 and we can approximate with Euklidian distance
                r = radial_distribution(d_ra, d_dec)
                if prog == 'marx':
                    j = 0
                else:
                    j = 1
                out[i, j, :] = np.percentile(r, np.array(ecf) * 100.)

        for i, e in enumerate(ecf):
            line, = axecf.plot(self.parameter, out[:, 0, i],
                               label='{0} %'.format(e*100))
            axecf.plot(self.parameter, out[:, 1, i], color=line.get_color(),
                       ls=':')
        axecf.legend(loc='lower right')
        axecf.set_ylabel('radius [arcsec]')
        axecf.set_xlabel('off-axis angle [arcmin]')
        axecf.grid()
        fig.savefig(self.figpath('ECF'))
