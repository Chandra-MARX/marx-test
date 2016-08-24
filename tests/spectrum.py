'''
|marx| allows users to specify either a flat input spectrum or pass in a file.
Here, we generate input spectra with a spectral modeling program.
If everything works well and all detector effects are treated consistently,
we can extract the spectrum from the simulated |marx| data, fit it, and recover
the spectral parameters.
'''
from textwrap import dedent
import os
import numpy as np
from collections import OrderedDict
from matplotlib import pyplot as plt

from .. import base
from ..utils import colname_case
from ..process_utils import marxpars_from_asol, spectrum_from_fluxcorrection

tests = ['SpectrumAbsPow']

title = 'Reproducing an input spectrum'


class SpectrumAbsPow(base.MarxTest):
    '''
    '''

    figures = OrderedDict([('spec', {'alternative': 'TBD',
                                    'caption': 'figure caption here.'})
                           ])


    summary = 'TBD'

    source = {'x': 4096.5,
              'y': 4096.5,
              'r': 20}
    '''Source position and extraction radius.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''

    @property
    def src_reg(self):
        return "circle({0},{1},{2})".format(self.source['x'],
                                            self.source['y'],
                                            self.source['r'])

    @base.Sherpa
    def step_1(self):
        sherpa = '''
        # set source properties
        set_source(xsphabs.a * xspowerlaw.p)
        a.nH = 1.
        p.PhoIndex = 1.8
        p.norm = 0.001
        # get source
        my_src = get_source()

        # set energy grid
        bin_width = 0.01
        energies = np.arange(0.15, 12., bin_width)

        # evaluate source on energy grid
        flux = my_src(energies)

        save_arrays("marx_input.tbl", [energies[1:], flux[:-1] / bin_width], ["keV","photons/s/cm**2/keV"], ascii=True, clobber=True)
        '''
        return dedent(sherpa)

    @base.Marx
    def step_2(self):
        pars = {'SourceFlux': -1,
                'SpectrumType': "FILE",
                'SpectrumFile': "marx_input.tbl",
                'ExposureTime': 30000,
                'TStart': 2015.5,
                'OutputDir': 'plaw',
                'GratingType': "NONE",
                'DetectorType': "ACIS-S",
                'DitherModel': "INTERNAL",
                'RA_Nom': 30,
                'Dec_Nom': 40,
                'SourceRA': 30,
                'SourceDEC': 40,
                }
        return pars

    @base.Marx2fits
    def step_3(self):
        return ('--pixadj=NONE', 'plaw', 'plaw_evt2.fits')

    @base.Marxasp
    def step_4(self):
        pars = {'MarxDir': "plaw",
                'OutputFile': "plaw_asol1.fits"
                }
        return pars

    @base.Ciao
    def step_5(self):
        '''
        Use special settings to account for steps that |marx| does not include:
        ACIS QE maps and current (non FEF) rmfs.
        '''
        commands = '''
        phagrid="pi=1:1024:1"

        evtfile="plaw_evt2.fits"
        asolfile="plaw_asol1.fits"
        phafile="plaw_pha.fits"
        rmffile="plaw_rmf.fits"
        arffile="plaw_arf.fits"
        asphistfile="plaw_asp.fits"

        asphist infile="$asolfile" outfile="$asphistfile" evtfile="$evtfile" clobber=yes

        dmextract infile="$evtfile[sky={src}][bin $phagrid]" outfile="$phafile" clobber=yes

        ccdid=7;

        detname="ACIS-$ccdid;UNIFORM;bpmask=0"
        grating="NONE"
        # For ACIS-I, use engrid="0.3:11.0:0.003". This reflects a limitation of mkrmf.
        engrid="0.3:12.0:0.003"

        mkarf mirror="hrma" detsubsys="$detname" grating="$grating" \
        outfile="$arffile" obsfile="$evtfile" engrid="$engrid" asphistfile="$asphistfile" \
        sourcepixelx={x} sourcepixely={y} maskfile=NONE pbkfile=NONE dafile=NONE clobber=yes verbose=0

        # Mean chip position
        chipx=512;
        chipy=512;

        fef="$CALDB/data/chandra/acis/fef_pha/acisD2000-01-29fef_phaN0005.fits"

        cxfilter="chipx_hi>=$chipx,chipx_lo<=$chipx"
        cyfilter="chipy_hi>=$chipy,chipy_lo<=$chipy"
        mkrmf infile="$fef[ccd_id=$ccdid,$cxfilter,$cyfilter]" outfile="$rmffile" axis1="energy=$engrid" axis2="$phagrid" thresh=1e-8 clobber=yes verbose=0
        '''.format(src=self.src_reg, x=self.source['x'], y=self.source['y'])
        commands = dedent(commands)
        return commands.split('\n')

    @base.Ciao
    def step_6(self):
        '''Use default CIAO tools to extract a spectrum'''
        command = 'specextract "plaw_evt2.fits[sky={src}]" spec asp=plaw_asp.fits weight=no clobber=yes mskfile=NONE badpixfile=NONE'.format(src=self.src_reg)
        return [command]

    @base.Sherpa
    def step_7(self):
        sherpa='''
        load_pha(1, 'plaw_pha.fits')
        load_rmf(1, 'plaw_rmf.fits')
        load_arf(1, 'plaw_arf.fits')
        load_data(2, 'spec_grp.pi')
        group_counts(1, 25)
        group_counts(2, 25)
        ignore(None, 0.3)
        ignore(7, None)
        # set source properties
        set_source(xsphabs.a * xspowerlaw.p)
        a.nH = 1.
        p.PhoIndex = 1.8
        p.norm = 0.001

        set_source(2, a * p)
        plot_fit(1)
        fit(1)
        plot_model(1, overplot=True)
        plot_data(2, overplot=True)
        fit(2)
        plot_model(2, overplot=True)
        log_scale(X_AXIS)
        log_scale(Y_AXIS)
        print_window('{out}')
        '''.format(out=self.figpath(self.figures.keys()[0]))
        return dedent(sherpa)
