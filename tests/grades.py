'''
When an X-ray photon (or a highly energetic particle) hits the ACIS detector, it produces a charge cloud in the silicon of the chip. The size and the shape of the this charge cloud depends strongly on the energy of the event, the type of the event, and on the characteristics of the detector.

In many cases, the charge cloud is large enough to produce detectable signal in more than one pixel. This information is encoded in the "flight grade" (see `the observatory guide <http://cxc.harvard.edu/proposer/POG/html/chap6.html#tth_sEc6.14>`_ for details). Flight grades are used to filter out background events by removing those grades where a large fraction of all detected events is due to particles and also in the sub-pixel event repositioning (`EDSER <http://cxc.harvard.edu/ciao/why/acissubpix.html>`_). For example, grade "0" events have charge in one pixel only, and grade "2" indicates an event that has charge in two neighboring pixels. Most likely, the grade "0" event occurred close to the center of the pixel and the grade "2" event close to the boundary of the two pixels. Thus, simulating the right event grade distribution is a prerequisite to using any sub-pixel event repositioning algorithm on simulated data.
'''
import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table\

from .. import base
from ..utils import colname_case

tests = ['ACIS_BI_low_energy', 'ACIS_FI']

title = 'Flight grade distribution'


def plotpie(ax, data, col, colors=['r', 'c', 'm', 'orange', 'g'],
            **plotargs):
    '''Bin up events data and display in pie chart'''
    counts = np.bincount(colname_case(data, col))
    labels = np.arange(len(counts))
    ind = counts > 0
    ax.pie(counts[ind], labels=labels[ind], colors=colors,
           autopct='%2i%%', **plotargs)


def plot22(obs, sim, colname, energies=[[300, 1000], [1000, 2000]]):
    '''Place 2*2 piecharts on a figure.
    Parameters
    ----------
    obs : `astropy.table.Table`
        Observed events
    sim : `astropy.table.Table`
        Events simulated with MARX
    energies : list
        Two energy bands in eV defined as a list of lists.

    Returns
    -------
    fig : `matplotlib.figure.Figure`
        Figure with plot
    '''
    fig = plt.figure(figsize=(10, 10))
    for i, data in enumerate([obs, sim]):
        for j, en in enumerate(energies):
            ax = fig.add_subplot(2, 2, i * 2 + j + 1)
            obssimlab = 'Observation' if i == 0 else 'MARX Simulation'
            energylab = '{0} - {1} eV'.format(en[0], en[1])
            ax.set_title('{0}/n{1}'.format(obssimlab, energylab))
            energy = colname_case(data, 'energy')
            plotpie(ax, data[(energy > en[0]) & (energy < en[1])], colname)
    return fig


class ACIS_BI_low_energy(base.MarxTest):
    '''This test compares the distribution of event grades between an observation of a young star observed on the back-illuminated S3 chip and the simulation. To avoid pile-up the observation used a sub-array read-out, but this does not influence the grade distribution.

The grade distribution depends on the energy of the incoming photons. To get things exactly right, marx thus has to run with the source spectrum as an input spectrum. However, the grade distribution changes slowly with energy, so this only matters when grade distribution over a large range of energy are compared. Here, marx is run with a constant input spectrum and the comparison between simulation and observation is done in narrow energy bands only.

In this particular case, there is very little source signal above 2 keV, so the comparison is limited to low energies.
'''

    obsid = 15713

    download_all = False

    summary = '''Both plots clearly show that MARX does not do a good job either at reproducing the distribution of grades in the first place in general, nor in reproducing the energy dependence. In the current implementation the MARX splits each pixel into a grid of 3*3 subpixels and the grade distribution for each of those sub-pixels is hardcoded.'''

    @base.Ciao
    def step_0(self):
        evt2 = self.get_data_file('evt2')
        return ["punlearn dmcopy",
                'dmcopy "{0}[sky=circle(4096,4073,6)]" obs.fits clobber=yes'.format(evt2)]

    @base.Marx
    def step_1(self):
        asol = self.get_data_file('asol')
        pars = self.marxpars_from_asol(asol)
        pars['MinEnergy'] = 0.3
        pars['MaxEnergy'] = 2.0
        return pars

    @base.Marx2fits
    def step_2(self):
        return ('--pixadj=EDSER', 'point', 'marxsim.fits')

    @base.Python
    def step_3(self):
        obs = Table.read(os.path.join(self.basepath, 'obs.fits'), hdu=1)
        sim = Table.read(os.path.join(self.basepath, 'marxsim.fits'), hdu=1)

        fig1 = plot22(obs, sim, 'grade')
        fig1.savefig(self.figpath('grades'))
        self.figures['grades'] = {'title': 'ASCA grades',
                                  'alternative': 'Four pie charts that display the distribution of grades for observation and simulation for two different energy bands',
                                  'caption': 'These pie charts show the distribution of grades in ASCA nomenclature (0-6). In a real observation the distribution of grades depends on the energy of the incoming photons.'
                                  }

        fig2 = plot22(obs, sim, 'fltgrade')
        fig2.savefig(self.figpath('fltgrades'))
        self.figures['grades'] = {'title': 'Flight grades',
                                  'alternative': 'Four pie charts that display the distribution of grades for observation and simulation for two different energy bands',
                                  'caption': 'Same as the above, but using the more detailed Chandra flight grade numbers.'
                                  }


class ACIS_FI(ACIS_BI_low_energy):
    '''This compares the grade distribution for a front-illuminated (FI) chip.
    '''

    obsid = 4496