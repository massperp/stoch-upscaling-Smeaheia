import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter

from matplotlib.ticker import LogLocator, AutoLocator, MultipleLocator, AutoMinorLocator
from matplotlib import ticker

   
CaseOrder = [1,2,3,4,5,6]
k_mods = [4,3,1]
kcs = [1e-4, 1e-3, 1]

N_maxs = [1000, 2000, 5000, 10000, 20000]
x_labels = ['$10^3$', '$5 \cdot 10^3$', '$10^4$', '$2 \cdot 10^4$']
x_labels = ['$10^3$', '$5 \cdot 10^3$']
curve_labels = ['Case I', 'Case II', 'Case III', 'Case IV', 'Case V', 'Case VI']

pow_of_10s = [0, 1, 2]
y_labels = ['$10^0$', '$10^{1}$', '$10^2$']

x_labels = ['$10^3$', '$2 \cdot 10^3$', '$5 \cdot 10^3$', '$10^4$', '$2 \cdot 10^{4}$']

speedup = np.ones((len(CaseOrder)*len(k_mods), len(N_maxs)))

for CaseIndex in range(len(CaseOrder)):
    CaseId = CaseOrder[CaseIndex]
    for k in range(len(k_mods)):
        k_model = k_mods[k]
        for N in range(len(N_maxs)):
            N_max = N_maxs[N]

            summary_path = "stratifications/strat_case_"+str(CaseId) + "_k" + str(k_model) + "_N" + str(N_max) + "/summary_ADSS.txt"
            
            text_file = open(summary_path, "r")

            # List of strings for each line
            lines = text_file.readlines()
            route = None
            for i in range(len(lines)):
                # Gets rid of brackets and newline and converts to numpy array
                arr = np.fromstring(lines[i].strip('[]\n'), dtype=float, sep= ',')
                # This just instantiates the array with the correct shape as soon as its known
                if i == 0:
                    route = np.zeros((len(lines), len(arr)))
                route[i] = arr

            QoI_ADSS_var = route.flatten()[1]
            
            path_QoI_MC2000 = 'results/case_' + str(CaseId) + '/MC_Nsamples_2000_kmod_' + str(k_model) + '.txt'
            path_QoI_MC3000 = 'results/case_' + str(CaseId) + '/MC_Nsamples_3000_kmod_' + str(k_model) + '.txt'
            
            QoI_MC_2000 = np.loadtxt(path_QoI_MC2000)
            QoI_MC_3000 = np.loadtxt(path_QoI_MC3000)
            QoI_MC_5000 = np.concatenate((QoI_MC_2000, QoI_MC_3000))
                        
            speedup[(CaseIndex)*len(k_mods)+k,N] = (np.var(QoI_MC_5000)/N_max)/QoI_ADSS_var
#print(speedup)


plt.rcParams['text.usetex'] = True
marker_and_line = ['--s', '--P', '--o', '--X','--v','--^']
fig, ax = plt.subplots(1,3, sharex=True, sharey=True, constrained_layout=True)
fig.set_size_inches(9.0, 3.0)
for mod_ind in range(len(k_mods)):
    for CaseIndex in range(len(CaseOrder)):
        ax[mod_ind].loglog(N_maxs, speedup[  (CaseIndex)*len(k_mods)+mod_ind,:], marker_and_line[CaseIndex], markersize=6, label = curve_labels[CaseIndex])
    plt.rcParams['axes.formatter.min_exponent'] = 2
    
    plt.ylim(np.power(10,pow_of_10s[0]), 7*np.power(10,pow_of_10s[-1]))
    plt.yticks(np.power(10,pow_of_10s), y_labels)
    
    plt.xticks(N_maxs, x_labels)
    plt.xscale('log', subs=[2, 3, 4, 5, 6, 7, 8, 9])
    
    ax[mod_ind].set_title(r'$k_{\rm c} = $' + str(kcs[mod_ind]) + ' mD')
    ax[mod_ind].set(xlabel = r'$N_{\rm s}$', ylabel = 'Speedup')
    
    ax[mod_ind].legend()

figname = 'figures/Speedups_Case_I_VI.png'
plt.savefig(figname)
plt.show()
