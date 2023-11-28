"""
Generate histograms of CO2 leakage for the six test cases and SGR models described in [ref to paper].
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from matplotlib.ticker import LogLocator, AutoLocator
from matplotlib import ticker


class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format
   
   
CaseOrder = [1,2,3,4,5,6]
labels_case = ['I','II','III','IV','V','VI']

k_model_order = [4,3,1]
k_c = ['1e-4mD','1e-3mD', '1mD']

for k_model_ind in range(len(k_model_order)):
    k_model = k_model_order[k_model_ind]

    fig, ax = plt.subplots(1,len(CaseOrder), sharex=True, sharey=True, constrained_layout=True)
    fig.set_size_inches(15.0, 3.0)
    fig.set_size_inches(2*len(CaseOrder), 2.5)

    for CaseIndex in range(len(CaseOrder)):

        CaseId = CaseOrder[CaseIndex]
        path_QoI_MC2000 = 'results/case_' + str(CaseId) + '/MC_Nsamples_2000_kmod_' + str(k_model) + '.txt'
        path_QoI_MC3000 = 'results/case_' + str(CaseId) + '/MC_Nsamples_3000_kmod_' + str(k_model) + '.txt'
            
        QoI_MC_2000 = np.loadtxt(path_QoI_MC2000)
        QoI_MC_3000 = np.loadtxt(path_QoI_MC3000)
        QoI_MC_hist = np.concatenate((QoI_MC_2000, QoI_MC_3000))
    
        ax[CaseIndex].hist(QoI_MC_hist, bins=30, density=True, color="lightsteelblue")

        #print("Mean(QoI): ", np.mean(QoI_MC_hist))
    
        ax[CaseIndex].set(xlabel="Leaked CO$_2$ (kSm$^3$)")
        ax[CaseIndex].set_title("Case " + labels_case[CaseIndex])
        ax[CaseIndex].margins(0.05)
        ax[CaseIndex].set_ylim(bottom=0)
        if k_model == 3:
            ax[CaseIndex].set_xlim(-10,750)
            ax[CaseIndex].set_ylim(0,4.5e-3)
        if k_model == 1:
            ax[CaseIndex].set_xlim(500,1400)
            ax[CaseIndex].set_ylim(0,1e-2)
    
        P10 = np.quantile(QoI_MC_hist, 0.1 , axis = None)
        P90 = np.quantile(QoI_MC_hist, 0.9 , axis = None)
        P50 = np.quantile(QoI_MC_hist, 0.5 , axis = None)
        bottom, top = plt.ylim()
        #ax[CaseIndex].vlines(x = [P10, np.mean(QoI_MC_hist), P90], ymin = bottom, ymax = top,
        #       colors = 'red', ls='--')
        ax[CaseIndex].vlines(x = [P10, P50, P90], ymin = bottom, ymax = top,
           colors = 'red', ls='--')
    
        matplotlib.rcParams['ytick.major.pad'] = 15
    
        plt.locator_params(axis='y', nbins=4)
    
        ax[CaseIndex].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1,1))
        ax[CaseIndex].yaxis.set_major_formatter(formatter)
    
        #print("Prop <1e-6 leak", np.size( np.where( QoI_MC_hist < 1e-6 ),axis=1)/5000)
    
        ax[CaseIndex].annotate("P10 = " + str(round(P10)), xy=(0.5, 0.8), xycoords='axes fraction')
        ax[CaseIndex].annotate("P50 = " + str(round(P50)), xy=(0.5, 0.7), xycoords='axes fraction')
        ax[CaseIndex].annotate("P90 = " + str(round(P90)), xy=(0.5, 0.6), xycoords='axes fraction')

    for axe in ax:
        axe.yaxis.set_major_formatter(OOMFormatter(-3, "%1.1f"))
        axe.ticklabel_format(axis='y', style='sci', scilimits=(-3,-3))

    
    figname = 'figures/QoI_MC_hist_kc_' + str(k_c[k_model_ind]) + '.png'
    plt.savefig(figname)
plt.show()
