'''
The `point-spread function (PSF) for Chandra <http://cxc.harvard.edu/ciao/PSFs/psf_central.html>`_ describes how the light from a point source is spread over a larger area on the detector. Several effects contribute to this, e.g. the uncertainty in the pointing, imperfections in the mirror (specifically for large off-axis angles) and the pixalization of data on detector read-out.

The following tests compare |marx| simulations, `SAOTrace`_ simulations, and data to look at different aspects of the Chandra PSF.
'''
import os
import shutil
import numpy as np
from collections import OrderedDict
from matplotlib import pyplot as plt
from astropy.table import Table
from astropy.io import fits

from .. import base
from ..utils import colname_case as cc

from ..process_utils import (marxpars_from_asol, spectrum_from_fluxcorrection,
                             target_coos_from_asol)

tests = ['ACISSPSF', 'ACISIPSF', 'HRCIPSF', 'OffAxisPSF']

title = 'Point Spread Function (PSF)'


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
    centx = x - np.median(x)
    centy = y - np.median(y)

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

    summary = '''TBD'''

    source = {'x': 16405,
              'y': 16500,
              'r': 15}
    '''Source position and extraction radius.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''

    @property
    def source_reg(self):
        return "circle({0}, {1}, {2})".format(self.source['x'],
                                              self.source['y'],
                                              self.source['r'])

    @base.Python
    def step_0(self):

        shutil.copy(os.path.join(self.pkg_data, 'ARLac_input_spec_marx.tbl'),
                    os.path.join(self.basepath, 'input_spec_marx.tbl'))
        shutil.copy(os.path.join(self.pkg_data, 'ARLac_input_spec_saotrace.rdb'),
                    os.path.join(self.basepath, 'input_spec_saotrace.rdb'))

    @base.Marx
    def step_1(self):
        asol = self.get_data_file('asol')
        evt = self.get_data_file('evt2')
        pars = marxpars_from_asol(asol, evt)
        pars['SpectrumType'] = 'FILE'
        pars['SpectrumFile'] = 'input_spec_marx.tbl'
        pars['OutputDir'] = 'marx_only'
        pars['SourceFlux'] = -1
        return pars

    @base.Marx2fits
    def step_2(self):
        return ('--pixadj=NONE', 'marx_only', 'marx_only.fits')

    @base.SAOTraceLua
    def step_3(self):
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
        asol = self.get_data_file('asol')
        asolh = fits.getheader(asol, 1)

        limit = asolh['TSTOP'] - asolh['TSTART']
        # CXO time numbers are large and round-off error can appear
        # which make SAOTrace fail.
        limit = limit - 0.02
        return ['trace-nest tag=saotrace srcpars=saotrace_source.lua tstart={tstart}  limit={limit} limit_type=sec'.format(tstart=asolh['TSTART'] + 0.01, limit=limit)]

    @base.Marx
    def step_5(self):
        asol = self.get_data_file('asol')
        evt = self.get_data_file('evt2')
        pars = marxpars_from_asol(asol, evt)
        pars['SpectrumType'] = 'FILE'
        pars['SpectrumFile'] = 'input_spec_marx.tbl'
        pars['OutputDir'] = 'marx_saotrace'
        pars['SourceType'] = 'SAOSAC'
        pars['SAOSACFile'] = 'saotrace.fits'
        pars['SourceFlux'] = -1
        return pars

    @base.Marx2fits
    def step_6(self):
        return ('--pixadj=NONE', 'marx_saotrace', 'marx_saotrace.fits')

    @base.Python
    def step_20(self):

        filterargs = {'circle': (self.source['x'], self.source['y'],
                                 5 * self.source['r'])
                      }

        evt2file = self.get_data_file('evt2')
        files = [evt2file, 'marx_only.fits', 'marx_saotrace.fits']

        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(111)
        title = [evt2file, 'marx_only.fits', 'marx_saotrace.fits']

        # In this energy range we have the most counts and we limited
        # the simulations to the same range.
        obs = filter_events(Table.read(files[0]), **filterargs)
        r = radial_distribution(obs['x'], obs['y'])
        val, edges = np.histogram(r, range=[0, 5], bins=25)
        ax.semilogx(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'k', lw=3, label='Observation')

        simmarx = filter_events(Table.read(files[1]), **filterargs)
        r = radial_distribution(simmarx['X'], simmarx['Y'])
        val, edges = np.histogram(r, range=[0, 5], bins=25)
        ax.semilogx(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'r', lw=3, label='MARX')

        simmarx = filter_events(Table.read(files[2]), **filterargs)
        r = radial_distribution(simmarx['X'], simmarx['Y'])
        val, edges = np.histogram(r, range=[0, 10], bins=50)
        ax.semilogx(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'b', lw=3, label='SAOTrace + MARX')

        ax.set_xlabel('radius [pixel]')
        ax.set_ylabel('enclosed count fraction')
        ax.legend(loc='upper left')
        ax.set_xlim([0.1, 5])
        ax.grid()

        fig.savefig(self.figpath(self.figures.keys()[0]))


class ACISSPSF(HRCIPSF):
    '''
    The PSF depends on many things, some of which are common to all observations
    like the shape of the mirror, and some are due to detector effects.
    For ACIS detectors, the sub-pixel event repositioning (EDSER) can improve
    the quality of an image, by repositioning events based on the event grade.
    This correction depends on the type pf chip (FI or BI). This test compares
    the simulation of a point source on a BI ACIS-S chip to an observation.
    The observed object is TYC 8241 2652 1, a young star, and was observed in
    1/8 sub-array mode to reduce pile-up. The pile-up fraction in the data is about
    5% in the brightest pixel.
    '''

    title = 'On-axis PSF on an ACIS-BI chip'
    obsid = 15713
    download_all = True

    figures = OrderedDict([('ECF', {'alternative': 'TBD',
                                    'caption': 'Enclosed count fraction for observation and simulations.'})
                       ])

    summary = '''TBD'''

    source = {'x': 4096,
              'y': 4074,
              'r': 4}
    '''Source position and extraction radius.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''
    @base.Ciao
    def step_0(self):
        asolfile = self.get_data_file('asol')
        evtfile = self.get_data_file('evt2')
        commands = spectrum_from_fluxcorrection(asolfile, evtfile,
                                                self.source['x'],
                                                self.source['y'],
                                                self.source_reg)
        return commands

    @base.Marx2fits
    def step_2(self):
        return ('--pixadj=EDSER', 'marx_only', 'marx_only.fits')

    @base.Marx2fits
    def step_6(self):
        return ('--pixadj=EDSER', 'marx_saotrace', 'marx_saotrace.fits')

    @base.Marx2fits
    def step_10(self):
        return ('--pixadj=RANDOMIZE', 'marx_only', 'marx_only_rand.fits')

    @base.Marx2fits
    def step_11(self):
        return ('--pixadj=RANDOMIZE', 'marx_saotrace', 'marx_saotrace_rand.fits')

    @base.Ciao
    def step_12(self):
        evt2 = self.get_data_file('evt2')
        asol = self.get_data_file('asol')
        return ['acis_process_events infile={evt} outfile=obs_rand.fits acaofffile={asol} doevtgrade=no calculate_pi=no pix_adj=RANDOMIZE clobber=yes'.format(evt=evt2, asol=asol)]

    @base.Python
    def step_20(self):

        filterargs = {'energy': [300, 3000],
                      'circle': (self.source['x'], self.source['y'],
                                 5 * self.source['r'])
                      }

        evt2file = self.get_data_file('evt2')

        fig = plt.figure(figsize=(10,5))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        for ax, files, title in zip([ax1, ax2],
                                    [[evt2file, 'marx_only.fits', 'marx_saotrace.fits'],
                                     ['obs_rand.fits', 'marx_only_rand.fits', 'marx_saotrace_rand.fits']],
                                    ['EDSER', 'RANDOMIZE']):

            # In this energy range we have the most counts and we limited
            # the simulations to the same range.
            obs = filter_events(Table.read(files[0]), **filterargs)
            r = radial_distribution(obs['x'], obs['y'])
            val, edges = np.histogram(r, range=[0, 5], bins=25)
            ax.semilogx(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'k', lw=3, label='Observation')

            simmarx = filter_events(Table.read(files[1]), **filterargs)
            r = radial_distribution(simmarx['X'], simmarx['Y'])
            val, edges = np.histogram(r, range=[0, 5], bins=25)
            ax.semilogx(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'r', lw=3, label='MARX')

            simmarx = filter_events(Table.read(files[2]), **filterargs)
            r = radial_distribution(simmarx['X'], simmarx['Y'])
            val, edges = np.histogram(r, range=[0, 5], bins=25)
            ax.semilogx(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'b', lw=3, label='SAOTrace + MARX')

            ax.set_title(title)
            ax.set_xlabel('radius [pixel]')
            ax.set_ylabel('enclosed count fraction')
            ax.legend(loc='lower right')
            ax.set_xlim([0.1, 5])
            ax.grid()

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

    summary = '''TBD'''

    source = {'x': 4140,
              'y': 4102,
              'r': 4}


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
    obsid = 15713
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
        # Use one of the asol files only.
        # Good enough for this example and simpler than stack syntax
        asol2 = self.get_data_file('N002_asol')
        os.remove(asol2)
        asolfile = self.get_data_file('asol')
        evt2 = self.get_data_file('N002_evt2')
        os.remove(evt2)
        evtfile = self.get_data_file('evt2')
        commands = spectrum_from_fluxcorrection(asolfile, evtfile,
                                                self.source['x'],
                                                self.source['y'],
                                                self.source_reg)
        return commands

    @base.Ciao
    def step_20(self):
        return ['''ds9 -width 800 -height 500 -log -cmap heat {0} marx_only.fits marx_saotrace.fits -pan to 5256 6890 physical -bin about 5256 6890 -match frame wcs -match bin -frame 1 -regions command 'text 5:39:27.987 -69:43:52.31 # text=Observation font="helvetica 24"' -frame 2 -regions command 'text 5:39:27.987 -69:43:52.31 # text="only MARX" font="helvetica 24"' -frame 3 -regions command 'text 5:39:27.987 -69:43:52.31 # text=SAOTrace font="helvetica 24"' -saveimage {1} -exit'''.format(self.get_data_file('evt2'), self.figpath(self.figures.keys()[0]))]
