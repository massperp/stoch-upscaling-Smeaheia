"""
Script to generate Figure 2: depth-dependent SGR and histograms for upscaled permeability in "Copula modeling and uncertainty propagation in field-scale simulation of CO2 fault leakage":
"""

import numpy as np
import matplotlib.pyplot as plt
#from src.upscaling_models import fault_facies

number_of_facies = 20 #
number_of_samples = 50000
init_loc_Vette = 700

k_model_order = [4,3,1]
k_c = ['1e-4mD','1e-3mD', '1mD']
k_max_values = [0.5, 1.3, 45]
font_size = 14.6
plt.rcParams.update({'font.size': font_size})


for k_model_ind in range(len(k_model_order)):
    k_model = k_model_order[k_model_ind]
    k_max = k_max_values[k_model_ind]

    pathName = 'copula_data/k_model_' + str(k_model) + '/'
    K_log_all = np.loadtxt(pathName + "K_col_log_samples.txt")
    Perm_samples_depth_dep = np.exp(K_log_all)
    
    fig, ax = plt.subplots(figsize= (3.2, 3.0))
    plt.hist(Perm_samples_depth_dep, density=True, bins=40)
    ax.set_xlabel('K (mD)')
    ax.set_ylabel('Norm. frequency')
    #plt.title('Upscaled perm., depth-dep. SGR')
    plt.tight_layout()
    ax.set_xlim(0,k_max)

    figname = 'figures/upscaled_perm_hist_kc_' + str(k_c[k_model_ind]) + '.eps'
    plt.savefig(figname, dpi=300)


SGR_data_Bjornaraa = np.loadtxt('Vette_data/SGR_synth_Bjornaraa.txt')
depth_data_Bjornaraa = np.loadtxt('Vette_data/SGR_synth_depth_Bjornaraa.txt')
# Correct reference level in data
depth_data_Bjornaraa = depth_data_Bjornaraa - 500

fig, ax = plt.subplots(figsize= (3.2, 3.0))
plt.plot(depth_data_Bjornaraa, SGR_data_Bjornaraa, 'ko-',markersize=1)
ax.set_xlim(np.max(depth_data_Bjornaraa), np.min(depth_data_Bjornaraa))
plt.vlines(x = [-700, -1200], ymin = 29, ymax = 64.5,
           colors = 'purple', ls='--')
ax.set_xlabel('Depth (m)')
ax.set_ylabel('SGR (%)')
plt.tight_layout()

figname = 'figures/Bjornaraa_SGR.eps'
plt.savefig(figname, dpi=300)

plt.show()

