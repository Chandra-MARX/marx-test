'''
When an X-ray photon (or a highly energetic particle) hits the ACIS detector, it produced a charge cloud in the silicon of the chip. The size and the shape of the this charge cloud depends strongly on the energy of the event, the type of the event, and on the characteristics of the detector.

In many cases, the charge cloud is large enough to produce detectable signal in more than one pixel. This information is encoded in the "flight grade" (see <a href="http://cxc.harvard.edu/proposer/POG/html/chap6.html#tth_sEc6.14" the observatory guide for details</a>). Flight grades are used to filter out background events by removing those grades where a large fraction of all detected events is due to particles and also in the sub-pixel event repositioning (<a href="http://cxc.harvard.edu/ciao/why/acissubpix.html" EDSER</a>). For example, grade "0" events have charge in one pixel only, and grade "2" indicates an event that has charge in two neighboring pixels. Most likely, the grade "0" event occurred close to the center of the pixel and the grade "2" event close to the boundary of the two pixels. Thus, simulating the right event grade distribution is a prerequisite to using any sub-pixel event repositioning algorithm on simulated data.
'''
import os
from matplotlib import pyplot as plt


from .. import base
from ..run_external import external_settings

__all__ = ['ACIS_BI']


class ACIS_BI(base.MarxTest):
    '''This test compares the distribution of event grades between an observation of a young star observed on the back-illuminated S3 chip and the simulation. To avoid pile-up the observation used a sub-array read-out, but this does not influence the grade distribution.

The grade distribution depends on the energy of the incoming photons. To get things exactly right, marx thus has to run with the source spectrum as an input spectrum. However, the grade distribution changes slowly with energy, so this only matters when grade distribution over a large range of energy are compared. Here, marx is run with a constant input spectrum and the comparison between simulation and observation is done in narrow energy bands only.
'''

    obsid = 15713

    download_all = False

    @base.Ciao
    def step_0(self):
        evt2 = self.get_data_file('evt2')
        return ["punlearn dmcopy",
                'dmcopy "{0}[sky=circle(4096,4073,6)]" obs.fits'.format(evt2)]

    @base.Marx
    def step_1(self):
        asol = self.get_data_file('asol')
        pars = self.marxpars_from_asol(asol)
        pars['MinEnergy'] = 0.3
        pars['MaxEnergy'] = 2.0
        pars['parfile'] = external_settings['marxparfile']
        return pars

    @base.Marx2fits
    def step_2(self):
        return ('--pixadj=EDSER', 'point', 'marxsim.fits')

    @base.Python
    def step_3(self):
        obs = Table.read(os.path.join(self.basepath, 'obs.fits'))
        sim = Table.read(os.path.join(self.basepath, 'marxsim.fits'))

        colors = ['r', 'c', 'm', 'orange', 'g']

        def plotpie(ax, data, col, **plotargs):
            counts = np.bincount(data[col])
            labels = np.arange(len(counts))
            ind = counts > 0
            ax.pie(counts[ind], labels=labels[ind], colors = colors, autopct='%2i%%', **plotargs)

        def plot22(obs, sim, colname):
            fig = plt.figure(figsize=(10, 10))
            for i, data in enumerate([obs, sim]):
                for j, en in enumerate([[0.3, 1.], [1, 2.]]):
                    ax = fig.add_subplot(i, j, i * 2 + j + 1)
                    if j == 0:
                        ax.text(0.1, 0.5, 'energy: {0} - {1}'.format(en[0], en[1]),
                                horizontalalignment='center',
                                verticalalignment='center',
                                transform=ax.transAxes)
                        if i == 0:
                            ax.set_title('Observation')
                        else:
                            ax.set_title('MARX simulation')
                    plotpie(ax, data[(data['energy'] > en[0]) & (data['energy'] < en[1])], 'grade')
            return fig

        fig1 = plot22(obs, sim, 'grade')
        fig1.savefig(os.path.join(self.env['outpath'],
                                  self.name + '_grades.png'))

        fig2 = plot22(obs, sim, 'fltgrade')
        fig2.savefig(os.path.join(self.env['outpath'],
                                  self.name + '_fltgrades.png'))
