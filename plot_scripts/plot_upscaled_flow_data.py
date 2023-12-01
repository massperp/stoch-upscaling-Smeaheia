import numpy as np
from matplotlib import pyplot as plt

N_samp_ref = 10000
num_pts = 21
data_samples_physical_all = np.zeros((5, num_pts, N_samp_ref))
flow_funs = ["Pc_log_samples.txt", "S_log_samples.txt", "Kw_log_samples.txt", "Knw_samples.txt"]

labels_param = ['$\log(P_c)$ (Pa)','$S$','$K_{r}^{w}$','$K_{r}^{nw}$']

k_model_order = [4,3,1]
k_c = ['1e-4mD','1e-3mD', '1mD']

for k_model_ind in range(len(k_model_order)):
    k_model = 'k_model_' + str(k_model_order[k_model_ind])
    
    paths = ["copula_data/" + k_model + "/", "copula_data_lin/" + k_model + "/", "copula_data_lin/" + k_model + "/", "copula_data/" + k_model + "/"]

    xlabs = [r'$\log(s_{\rm d})$',r'$s_{\rm d}$',r'$s_{\rm d}$',r'$s_{\rm d}$']
    fig, ax = plt.subplots(1,4, sharex=False, sharey=False, constrained_layout=True)
    fig.set_size_inches(12.0, 3.0)
    
    for plot_ind in range(4):
    
        xlab = xlabs[plot_ind]
        s_disc = np.loadtxt(paths[plot_ind] + "s_disc.txt")

        data_samples_physical_all[plot_ind,:,:] = np.loadtxt(paths[plot_ind] + flow_funs[plot_ind]).T
    
        if plot_ind == 0:
            data_samples_physical_all[plot_ind,:,:] = np.log10(np.exp(data_samples_physical_all[plot_ind,:,:]))
        #if plot_ind == 1 or plot_ind == 2:
        #    data_samples_physical_all[plot_ind,:,:] = np.exp(data_samples_physical_all[plot_ind,:,:])

        if plot_ind == 1 or plot_ind == 2:
            xlab = r'$s_{\rm d}$'
            ax[plot_ind].plot(s_disc,np.exp(data_samples_physical_all[plot_ind,:,:]),'ko',markersize=1)
            ax[plot_ind].plot(s_disc, np.mean(np.exp(data_samples_physical_all[plot_ind,:,:]),axis=1),'-r',markersize=1)
        else:
            xlab = r'$\log(s_{\rm d})$'
            ax[plot_ind].plot(np.log10(s_disc),data_samples_physical_all[plot_ind,:,:],'ko',markersize=1)
            ax[plot_ind].plot(np.log10(s_disc),np.mean(data_samples_physical_all[plot_ind,:,:], axis=1),'-r',markersize=1)

        ax[plot_ind].set(xlabel = xlab, ylabel=labels_param[plot_ind])

    figname = 'figures/flow_funs_vs_sdisc_' + str(k_c[k_model_ind]) +'.png'
    plt.savefig(figname)
plt.show()

