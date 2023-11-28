"""
Script to generate Figure 2 in [ref to paper]: depth-dependent SGR and histograms for upscaled permeability
"""

import numpy as np
import matplotlib.pyplot as plt
from src.upscaling_models import fault_facies

number_of_facies = 20 #
number_of_samples = 50000
init_loc_Vette = 700

k_model_order = [4,3,1]
k_c = ['1e-4mD','1e-3mD', '1mD']

for k_model_ind in range(len(k_model_order)):
    k_model = k_model_order[k_model_ind]

    if k_model == 4:
        k_max = 0.5
    
        data = {
        "mu_SGR": 0.50,
        "sigma_SGR": 0.14,
        "perm_clay": 1e-4,#1, #1e-3,
        "perm_sandstone": 1000,
        "fault_height_Vette": 500,
        }
    
    if k_model == 3:
        k_max = 1.3
    
        data = {
        "mu_SGR": 0.50,
        "sigma_SGR": 0.14,
        "perm_clay": 1e-3,#1, #1e-3,
        "perm_sandstone": 1000,
        "fault_height_Vette": 500,
        }
    if k_model == 1:
        k_max = 45
    
        data = {
        "mu_SGR": 0.50,
        "sigma_SGR": 0.14,
        "perm_clay": 1,#1, #1e-3,
        "perm_sandstone": 1000,
        "fault_height_Vette": 500,
        }

    Perm_samples_depth_dep = fault_facies(number_of_facies, number_of_samples, data, data_interpr='depth_dep',s_disc=1e-6, k_model=k_model,sect='Vette', init_loc=init_loc_Vette,  two_phase=False)

    fig, ax = plt.subplots(figsize= (3.2, 3.0))
    plt.hist(Perm_samples_depth_dep, density=True, bins=40)
    ax.set_xlabel('k (mD)')
    ax.set_ylabel('Probability density')
    #plt.title('Upscaled perm., depth-dep. SGR')
    plt.tight_layout()
    ax.set_xlim(0,k_max)

    figname = 'figures/upscaled_perm_hist_kc_' + str(k_c[k_model_ind]) + '.png'
    plt.savefig(figname)


SGR_data_Bjornaraa = np.loadtxt("SGR_synth_Bjornaraa.txt")
depth_data_Bjornaraa = np.loadtxt("SGR_synth_depth_Bjornaraa.txt")
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
#plt.title('Bjørnaraa et al SGR - depth')

figname = 'figures/Bjornaraa_SGR.png'
plt.savefig(figname)

plt.show()

