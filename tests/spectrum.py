'''
|marx| allows users to specify either a flat input spectrum or to pass in a file.
Here, we generate input spectra with the spectral modeling program `Sherpa`_.
We run |marx| to
simulate a point source and extract the data from that simulation using `CIAO`_,
read it back into `Sherpa`_ and fit the simulated data.
If everything works well and all detector effects are treated consistently,
we should recover the same spectral parameters (up to Poisson noise).
On the other hand, if |marx| and
`CIAO`_ are inconsistent the fit parameters will deviate from the input
parameters.

In the past, this has happened after the release of |marx| 5.0, which contains
some files to describe the ACIS contamination. This contamination changes with
time and several `CalDB`_ releases had happened before we released |marx| 5.1.
At that time, |marx| always predicted too many counts at low energies.

The following tests are designed to test **consistency** with `CIAO`_.
Since there is always some uncertainty about the intrinsic spectrum of
astrophysical sources, this test is best done with simulated input spectra.
'''

from textwrap import dedent
import json
import os
import numpy as np
from collections import OrderedDict

from marxtest import base

tests = ['SpectrumAbsPowACISS', 'SpectrumAbsPowACISS_offaxis'
         'SpectrumAPECACISI']

title = 'Reproducing an input spectrum'


class ModelComparisonError(Exception):
    pass


class SpectrumAbsPowACISS(base.MarxTest):
    '''This test checks the internal consistency of a |marx| spectral simulation
    by simulating a source placed on an back-illuminated chip of ACIS-S.
    The input spectrum is an absorbed powerlaw.
    '''

    title = 'Absorbed powerlaw on ACIS-S'

    figures = OrderedDict([('spec', {'alternative': 'Spectra and model shown in top of each other. Differences are described in the text.',
                                     'caption': 'Data reduced as required to reproduce the |marx| setup (see :ref:`caveats`) is shown in green, data reduced with standard `CIAO`_ setting in orange. The orange error bars are omitted for clarity. In this display, green and orange points should be on top of each other, but the next figure shows that real differences exist when ARF is needed to calculate the flux. Lines show the input model (red, using the ARF of the green data), the best-fit model to the green data (green line) and to the orange data (orange line).'}),
                           ('arf', {'alternative': 'ARFs differ in the range 2-5 keV',
                                    'caption': 'ARFs for both extraction methods (some color as above).'}),
                           ])

    expresults = [{'name': 'chi2marx', 'title': 'chi^2 of input model',
                   'description': 'chi^2 of the data assuming the input model. The data extraction takes into account MARX li/melkor/d1/guenther/marx/doc/source/tests/spectrum.rstmitations',
                   'value': 1},
                  {'name': 'chi2default', 'title': 'chi^2 of input model',
                   'description': 'chi^2 of the data assuming the input model for a naive extraction, **not** taking into account MARX limitations',
                   'value': 1},
                  {'name': 'arfdiff', 'title': 'ARFs: rel diff', 'description': 'Relative difference between a naive arf calcualted with CIAO standard settings and an ARF taking into account MARX limitations. The number is averaged over the entire energy range.',
                   'value': 0}]


    _summary = '''
    The table shows how the best fit of the simulated and extracted data compares to the
    parameter values that were used to generate the input spectrum. The two columns compare
    the fit results using an extraction that takes |marx|'s :ref:`caveats` into account
    and an extraction that just applies `CIAO`_ default values.

    In a Monte-Carlo simulation we can never expect the fit to hit the input values exactly.
    In fact, a perfect agreement between input and output numbers would indicate
    that the random number generator is not random enough. When running this example
    many times, the average distance should be roughly 1 sigma.
    '''
    source = {'x': 4096.5,
              'y': 4096.5,
              'r': 20,
              'ccd_id': 7,
              'detector_type': "ACIS-S"}
    '''Source position and extraction radius.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''

    input_model = {'set_source': 'set_source(xsphabs.a * xspowerlaw.p)',
                   'set_source2': 'set_source(2, a * p)',
                   'parnames': ('a.nH', 'p.PhoIndex', 'p.norm'),
                   'parvals': (1., 1.8, 0.001)}
    '''Define the spectral model that generates the input spectrum and is
    fit later.
    ``set_source2``: Set the same model to data id=2. Yes, this could be made more
    general and even be extracted automatically, but it's easier to just list
    the string here.
    '''

    def input_model2Sherpa(self):
        if isinstance(self.input_model['set_source'], list):
            out = self.input_model['set_source']
        else:
            out = [self.input_model['set_source']]
        for n, v in zip(self.input_model['parnames'],
                        self.input_model['parvals']):
            out.append('{0} = {1}'.format(n, v))
        # Spaces are necessary to match the indentation in the long string
        # in step_7
        return '\n        '.join(out)

    @property
    def src_reg(self):
        return "circle({0},{1},{2})".format(self.source['x'],
                                            self.source['y'],
                                            self.source['r'])

    @property
    def summary(self):
        with open(os.path.join(self.basepath, 'sherpaout.json')) as f:
            sher = json.load(f)

        tablines = []
        for n, v in zip(self.input_model['parnames'], self.input_model['parvals']):
            i = np.nonzero(np.array(sher['fitmarx']['parnames']) == n)[0][0]
            tablines.append('{0:12}  {1:6.2}  {2:8.2}  {3:8.2}  {4:8.2}  {5:8.2}  {6:8.2}  {7:8.2}'.format(n, v,
                                      sher['fitmarx']['parvals'][i],
                                      sher['fitmarx']['parmins'][i],
                                      sher['fitmarx']['parmaxes'][i],
                                      sher['fitdefault']['parvals'][i],
                                      sher['fitdefault']['parmins'][i],
                                      sher['fitdefault']['parmaxes'][i]))
        tablines = '\n        '.join(tablines)

        # spaces in next line needed to make dedent work
        summary = '''
        Fit results:

        ============  ======  ========  ========  ========  ========  ========  ========
        Parameter     Input   MARX special extraction       CIAO default extraction
        ------------  ------  ----------------------------  ----------------------------
        name          value   value     err_down  err_up    value     err_down  err_up
        ============  ======  ========  ========  ========  ========  ========  ========
        {tablines}
        ============  ======  ========  ========  ========  ========  ========  ========
        '''.format(tablines=tablines)

        return dedent(summary) + '\n' + dedent(self._summary)

    @base.Sherpa
    def step_1(self):
        '''Generate input spectrum'''
        sherpa = '''
        # set source properties
        {srcstring}
        # get source
        my_src = get_source()

        # set energy grid
        bin_width = 0.01
        energies = np.arange(0.15, 12., bin_width)

        # evaluate source on energy grid
        flux = my_src(energies)

        save_arrays("marx_input.tbl", [energies[1:], flux[:-1] / bin_width], ["keV","photons/s/cm**2/keV"], ascii=True, clobber=True)
        '''.format(srcstring=self.input_model2Sherpa())
        return dedent(sherpa)

    @base.Marx
    def step_2(self):
        '''Point source'''
        pars = {'SourceFlux': -1,
                'SpectrumType': "FILE",
                'SpectrumFile': "marx_input.tbl",
                'ExposureTime': 30000,
                'TStart': 2015.5,
                'OutputDir': 'marx',
                'GratingType': "NONE",
                'DetectorType': self.source['detector_type'],
                'DitherModel': "INTERNAL",
                'RA_Nom': 30,
                'Dec_Nom': 40,
                'SourceRA': 30,
                'SourceDEC': 40,
                }
        return pars

    @base.Marx2fits
    def step_3(self):
        '''turn into fits file'''
        return ('--pixadj=NONE', 'marx', 'marx_evt2.fits')

    @base.Marxasp
    def step_4(self):
        '''Make asol file for CIAO'''
        pars = {'MarxDir': "marx",
                'OutputFile': "marx_asol1.fits"
                }
        return pars

    @base.Ciao
    def step_5(self):
        '''Extract spectra

        Use special settings to account for steps that |marx| does not include:
        ACIS QE maps and current (non FEF) rmfs.
        '''
        commands = '''
        phagrid="pi=1:1024:1"

        evtfile="marx_evt2.fits"
        asolfile="marx_asol1.fits"
        phafile="marx_pha.fits"
        rmffile="marx_rmf.fits"
        arffile="marx_arf.fits"
        asphistfile="marx_asp.fits"

        asphist infile="$asolfile" outfile="$asphistfile" evtfile="$evtfile" clobber=yes

        dmextract infile="$evtfile[sky={src}][bin $phagrid]" outfile="$phafile" clobber=yes

        dmcoords $evtfile asol=$asolfile option=sky x={x} y={y}
        ccdid=$(pget dmcoords chip_id)
        chipx=$(pget dmcoords chipx)
        chipy=$(pget dmcoords chipy)

        detname="ACIS-$ccdid;UNIFORM;bpmask=0"
        grating="NONE"
        # For ACIS-I, use engrid="0.3:11.0:0.003". This reflects a limitation of mkrmf.
        engrid="0.3:12.0:0.003"

        mkarf mirror="hrma" detsubsys="$detname" grating="$grating" \
        outfile="$arffile" obsfile="$evtfile" engrid="$engrid" asphistfile="$asphistfile" \
        sourcepixelx={x} sourcepixely={y} maskfile=NONE pbkfile=NONE dafile=NONE clobber=yes verbose=0

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
        command = 'specextract "marx_evt2.fits[sky={src}]" spec asp=marx_asp.fits weight=no clobber=yes mskfile=NONE badpixfile=NONE'.format(src=self.src_reg)
        return [command]

    @base.Sherpa
    def step_7(self):
        '''Compare spectra to input model'''
        sherpa = '''
        load_pha(1, 'marx_pha.fits')
        load_rmf(1, 'marx_rmf.fits')
        load_arf(1, 'marx_arf.fits')
        load_data(2, 'spec_grp.pi')
        group_counts(1, 25)
        group_counts(2, 25)
        ignore(None, 0.3)
        ignore(7, None)
        # set source properties
        {srcstring}

        {srcstring2}

        # Check how good the input parameters work
        statres = get_stat_info()
        out = {{'chi2marx': statres[0].rstat,
                'chi2default': statres[1].rstat}}

        # Make the plots
        plot_data(1)
        plot_model(1, overplot=True)
        set_curve({{'*.color': 'forest'}})
        set_histogram({{'*.color': 'forest'}})
        fit(1)
        conf(1)
        c1 = get_conf_results()
        plot_model(1, overplot=True)
        plot_data(2, overplot=True)
        fit(2)
        conf(2)
        c2 = get_conf_results()
        plot_model(2, overplot=True)
        log_scale(X_AXIS)
        log_scale(Y_AXIS)
        set_curve({{'*.color': 'orange', 'err.*': 'false', 'symbol.style': 'uptriangle'}})
        set_histogram(['*.color', 'orange'])
        print_window('{out1}', ['export.clobber', 'True'])

        plot_arf(1)
        set_histogram({{'*.color': 'forest'}})
        plot_arf(2, overplot=True)
        set_histogram(['*.color', 'orange'])
        print_window('{out2}', ['export.clobber', 'True'])

        # compile conf outputs for saving
        out['fitmarx'] = {{'parnames': c1.parnames,
                           'parmins': c1.parmins,
                           'parvals': c1.parvals,
                           'parmaxes': c1.parmaxes}}
        out['fitdefault'] = {{'parnames': c2.parnames,
                              'parmins': c2.parmins,
                              'parvals': c2.parvals,
                              'parmaxes': c2.parmaxes}}

        # Calculate maximum relative difference in arf
        a1 = get_arf(1)
        a2 = get_arf(2)
        # bring on same energy grid
        engrid = np.arange(0.3, 10., 0.05)
        a1 = np.interp(engrid, a1.get_x(), a1.get_y())
        a2 = np.interp(engrid, a2.get_x(), a2.get_y())

        out['arfdiff'] = np.mean(np.abs(a2/a1 -1))

        # write numbers to file for later
        import json
        with open('{sherpaout}', 'w') as f:
            json.dump(out, f)

        '''.format(srcstring=self.input_model2Sherpa(),
                   srcstring2=self.input_model['set_source2'],
                   out1=self.figpath(self.figures.keys()[0]),
                   out2=self.figpath(self.figures.keys()[1]),
                   sherpaout=os.path.join(self.basepath, 'sherpaout.json'))
        return dedent(sherpa)

    @base.Python
    def step_8(self):
        '''Compare input with fit model

        The input model has several parameters. The calculation here summarizes
        all that into a single number: The mean deviation between input and fit
        result as measured in terms of Gaussian sigmas.
        That is a massive simplification since, the distribution is typically
        not Gaussian shaped and parameters might not be independent (e.g. there
        is ambiguity between the absorbing column density and the amount of
        cool plasma).
        '''
        with open(os.path.join(self.basepath, 'sherpaout.json')) as f:
            sherpaout = json.load(f)

        self.save_test_result('chi2marx', sherpaout['chi2marx'])
        self.save_test_result('chi2default', sherpaout['chi2default'])
        self.save_test_result('arfdiff', sherpaout['arfdiff'])


class SpectrumAbsPowACISS_offaxis(SpectrumAbsPowACISS):
    '''Same as above, but for an off-axis source.'''

    title = 'Powerlaw spectrum of an off-axis source'

    _summary = '''
    The table shows how the best fit of the simulated and extracted data
    compares to the parameter values that were used to generate the input
    spectrum.

    See above for a detailed description.
    '''
    source = {'x': 4200,
              'y': 5200,
              'r': 200,
              'ccd_id': 1,
              'detector_type': "ACIS-I"}
    '''Source position and extraction radius.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''

    @base.Marx
    def step_2(self):
        '''Point source'''
        pars = {'SourceFlux': -1,
                'SpectrumType': "FILE",
                'SpectrumFile': "marx_input.tbl",
                'ExposureTime': 30000,
                'TStart': 2015.5,
                'OutputDir': 'marx',
                'GratingType': "NONE",
                'DetectorType': self.source['detector_type'],
                'DitherModel': "INTERNAL",
                'RA_Nom': 30,
                'Dec_Nom': 40,
                'SourceRA': 29.9814917,  # hardcode pos here because much
                                         # easier than calling ciao dmcoords
                'SourceDEC': 40.1508083,
                }
        return pars



class SpectrumAPECACISI(SpectrumAbsPowACISS):
    '''Same as above, but for a front-illuminated ACIS-I chip and using
    an input spectrum with two thermal components similar to a stellar corona.
    '''

    title = 'Two thermal components on ACIS-I'

    _summary = '''
    The table shows how the best fit of the simulated and extracted data
    compares to the parameter values that were used to generate the input
    spectrum.

    See above for a detailed description.
    '''
    source = {'x': 4096.5,
              'y': 4096.5,
              'r': 20,
              'ccd_id': 3,
              'detector_type': "ACIS-I"}
    '''Source position and extraction radius.

    This is used to automatically extract a source spectrum that can be used as
    input for |marx| and `SAOTrace`_ runs.
    x,y,r are given in "physical" pixel units on the detector.
    '''

    input_model = {'set_source': ['set_source(xsvapec.a1 + xsvapec.a2)',
                                  'a1.Ne.frozen = False',
                                  'a2.Ne = a1.Ne'],
                   'set_source2': 'set_source(2, a1 + a2)',
                   'parnames': ('a1.kT', 'a1.Ne', 'a1.norm', 'a2.kT', 'a2.norm'),
                   'parvals': (0.7, 2.0, 0.0002, 2.0, 0.0004)}
