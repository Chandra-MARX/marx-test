'''
The `point-spread function (PSF) for Chandra <http://cxc.harvard.edu/ciao/PSFs/psf_central.html>`_ describes how the light from a point source is spread over a larger area on the detector. Several effects contribute to this, e.g. the uncertainty in the pointing, imperfections in the mirror (specifically for large off-axis angles) and the pixalization of data on detector read-out.

The following tests compare |marx| simulations, `SAOTrace`_ simulations, and data to look at different aspects of the Chandra PSF.
'''

import numpy as np
from collections import OrderedDict
from matplotlib import pyplot as plt
from astropy.table import Table
from astropy.io import fits

from .. import base
from ..utils import colname_case as cc

from ..process_utils import (marxpars_from_asol, spectrum_from_fluxcorrection,
                             target_coos_from_asol)

tests = ['ACISSPSF']

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


class ACISSPSF(base.MarxTest):
    '''The PSF depends on many things, some of which are common to all observations
    like the shape of the mirror, and some are due to detector effects.
    For ACIS detectors, the sub-pixel event repositioning (EDSER) can improve
    the quality of an image, by repositioning events based on the event grade.
    This correction depends on the type pf chip (FI or BI). This test compares
    the simulation of a point source on a BI ACIS-S chip to an observation.
    '''

    obsid = 15713
    download_all = True

    source = {'x': 4096,
              'y': 4074,
              'r': 4}
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

    figures = OrderedDict([('ECF', {'alternative': 'TBD',
                                      'caption': 'Enclosed count fraction for observation and simulations.'})
                       ])

    summary = '''TBD'''

    @base.Ciao
    def step_0(self):
        asolfile = self.get_data_file('asol')
        evtfile = self.get_data_file('evt2')
        commands = spectrum_from_fluxcorrection(asolfile, evtfile,
                                                self.source['x'],
                                                self.source['y'],
                                                self.source_reg)
        return commands

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
        return ('--pixadj=EDSER', 'marx_only', 'marx_only.fits')

    @base.SAOTraceLua
    def step_3(self):
        asol = self.get_data_file('asol')
        asolh = fits.getheader(asol, 1)
        ra, dec = target_coos_from_asol(asol)
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
           ra=ra, dec=dec)

    @base.SAOTrace
    def step_4(self):
        asol = self.get_data_file('asol')
        asolh = fits.getheader(asol, 1)

        return ['trace-nest tag=saotrace srcpars=saotrace_source.lua tstart={tstart}  limit={limit} limit_type=sec'.format(tstart=asolh['TSTART'], limit=asolh['TSTOP'] - asolh['TSTART'])]

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
        return ('--pixadj=EDSER', 'marx_saotrace', 'marx_saotrace.fits')

    @base.Python
    def step_7(self):

        filterargs = {'energy': [300, 3000],
                      'circle': (self.source['x'], self.source['y'],
                                 5 * self.source['r'])
                      }

        evt2file = self.get_data_file('evt2')
        obs = Table.read(evt2file)
        # In this energy range we have the most counts and we limited
        # the simulations to the same range.
        obs = filter_events(obs, **filterargs)
        r = radial_distribution(obs['x'], obs['y'])
        val, edges = np.histogram(r, range=[0, 5], bins=25)
        plt.plot(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'k', lw=3, label='Obs')

        simmarx = filter_events(Table.read('marx_only.fits'), **filterargs)
        r = radial_distribution(simmarx['X'], simmarx['Y'])
        val, edges = np.histogram(r, range=[0, 5], bins=25)
        plt.plot(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'r', lw=3, label='MARX')

        simmarx = filter_events(Table.read('marx_saotrace.fits'), **filterargs)
        r = radial_distribution(simmarx['X'], simmarx['Y'])
        val, edges = np.histogram(r, range=[0, 5], bins=25)
        plt.plot(0.1 + edges[:-1], 1.0 * val.cumsum() / val.sum(), 'r', lw=3, label='MARX + SAOTrace')

        plt.title('EDSER')
        plt.xlabel('radius [pixel]')
        plt.ylabel('enclosed count fraction')
        plt.savefig(self.figpath('ECF'))