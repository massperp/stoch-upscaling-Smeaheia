
"""
Plot 2D projections of upscaled data, copula models, and local PDF estimates.
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import scipy.stats as stats
from src.upscaled_flow_functions import upscaled_Vette_flow_functions as upscaled_Vette


#k_models = ['k_model_4','k_model_3','k_model_1']

labels_param = ['$\log(P_c)$ (Pa)','$S$','$K_{r}^{w}$']
labels_param_fname = ['logPc','S','Kw']
var_order = np.array([[1,0],[1,2]])

s_disc_mods = ['log', 'lin']
k_model_order = [4,3,1]
k_c = ['1e-4mD','1e-3mD', '1mD']

for k_model_ind in range(len(k_model_order)):
    k_model = k_model_order[k_model_ind]

    for row_ind in range(2):
        s_disc_mod = s_disc_mods[row_ind]
        if s_disc_mod == 'log':
            pathName = "copula_data/k_model_" + str(k_model) + "/"
        elif s_disc_mod == 'lin':
            pathName = "copula_data_lin/k_model_" + str(k_model) + "/"

        s_disc = np.loadtxt(pathName + "s_disc.txt")

        # Independent uniforms
        N_samp, num_dim = 10000, 5
        u_samples_ind = np.random.uniform(low=0, high=1, size=(N_samp, num_dim))
        model_logK_cop_samples, model_logPc_cop_samples, model_logS_cop_samples, model_logKw_cop_samples, model_Knw_cop_samples = upscaled_Vette(pathName, u_samples_ind, k_model, s_disc_mod)
      
      
        logPc_all = np.loadtxt(pathName + "Pc_log_samples.txt").T
        logS_all = np.loadtxt(pathName + "S_log_samples.txt").T
        logKw_all = np.loadtxt(pathName + "Kw_log_samples.txt").T

        n_x, N_samp_ref = logKw_all.shape


        data_samples_physical_all = np.zeros((3,n_x, N_samp_ref))
        data_samples_physical_all[:,:,:] = [np.log10(np.exp(logPc_all)), np.exp(logS_all), np.exp(logKw_all)]

        copula_test_samples_all = np.zeros((3,n_x,N_samp))
        copula_test_samples_all[:,:,:] = [np.log10(np.exp(model_logPc_cop_samples)), np.exp(model_logS_cop_samples), np.exp(model_logKw_cop_samples)]

    
        ind = var_order[row_ind,:]
    
        xlab = labels_param[ind[0]]
        ylab = labels_param[ind[1]]
    
        fig, ax = plt.subplots(1,2, sharex=True, sharey=True, constrained_layout=True)
        fig.set_size_inches(6.0, 3.0)
    
        cmap1 = mpl.cm.Blues(np.linspace(0,1,n_x))
        cmap2 = mpl.cm.OrRd(np.logspace(1e-6,1,n_x))
        cmap2 = mpl.cm.OrRd(np.linspace(1e-6,1,n_x))
        
        if s_disc_mod == 'lin':
            cmap1 = mpl.colors.ListedColormap(cmap1[14:,:-1])
            cmap2 = mpl.colors.ListedColormap(cmap2[14:,:-1])
            ax[0].scatter(data_samples_physical_all[ind[0],:,:].T, data_samples_physical_all[ind[1],:,:].T, s=1, c = np.tile(np.arange(0,n_x),(N_samp_ref,1)), cmap=cmap1)
        if s_disc_mod == 'log':
            cmap1 = mpl.colors.ListedColormap(cmap1[10:,:-1])
            cmap2 = mpl.colors.ListedColormap(cmap2[10:,:-1])
            ax[0].scatter(data_samples_physical_all[ind[0],:,:].T, data_samples_physical_all[ind[1],:,:].T, s=1, c = n_x*np.tile(np.logspace(1e-6,1,n_x ),(N_samp_ref,1)), cmap=cmap1)

        ax[0].set_title("Upscaled data")
        if s_disc_mod == 'lin':
            ax[1].scatter(copula_test_samples_all[ind[0],:].T, copula_test_samples_all[ind[1],:,:].T, s=1, c = np.tile(np.arange(0,n_x),(N_samp,1)), cmap=cmap2)
        if s_disc_mod == 'log':
            ax[1].scatter(copula_test_samples_all[ind[0],:].T, copula_test_samples_all[ind[1],:,:].T, s=1, c = n_x*np.tile(np.logspace(1e-6,1,n_x),(N_samp,1)), cmap=cmap2)
    
        ax[1].set_title("Copula samples")
        
        for ax_ind in ax.flat:
            ax_ind.set(xlabel = xlab, ylabel=ylab)
        
        figname = 'figures/' + labels_param_fname[ind[1]] + '_vs_'+labels_param_fname[ind[0]]+ '_kc_' + str(k_c[k_model_ind]) + '.png'
        
        plt.savefig(figname)


        # PDF plots
        if s_disc_mod == 'lin':
            s_ind = 10
        elif s_disc_mod == 'log':
            s_ind = 19 #10 #19
            
        fig, ax = plt.subplots(1,2, sharex=True, sharey=True, constrained_layout=True)
        #fig.set_size_inches(12.0, 3.0)
        fig.set_size_inches(6.0, 3.0)
    
    
        x = data_samples_physical_all[ind[0],s_ind,:]
        y = data_samples_physical_all[ind[1],s_ind,:]
        
        xmin, xmax = np.amin(x), np.amax(x)
        ymin, ymax = np.amin(y), np.amax(y)
    
        # Peform the kernel density estimate
        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x, y])
        kernel = stats.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)

        ax[0].set_xlim(xmin, xmax)
        ax[0].set_ylim(ymin, ymax)
        # Contourf plot
        cfset = ax[0].contourf(xx, yy, f, cmap='Blues')
        # Contour plot
        cset = ax[0].contour(xx, yy, f, colors='k')
        ax[0].set_title("PDF upscaled data")
        ax[0].plot(np.mean(x), np.mean(y ), 'rx', markersize=10)
        ax[0].axis('tight')
    
        #### Subplot 2
        x = copula_test_samples_all[ind[0],s_ind,:]
        y = copula_test_samples_all[ind[1],s_ind,:]
    
        values = np.vstack([x, y])
        kernel = stats.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)

        ax[1].set_xlim(xmin, xmax)
        ax[1].set_ylim(ymin, ymax)
        # Contourf plot
        cfset = ax[1].contourf(xx, yy, f, cmap='Blues')
  
        # Contour plot
        cset = ax[1].contour(xx, yy, f, colors='k')
        ax[1].set_title("PDF copula")
        ax[1].plot(np.mean(x), np.mean(y ), 'rx', markersize=10)
        ax[1].axis('tight')

        for ax_ind in ax.flat:
            ax_ind.set(xlabel = xlab, ylabel=ylab)

        figname = 'figures/' + labels_param_fname[ind[1]] + '_vs_'+labels_param_fname[ind[0]] + '_local_pdf_kc_' + str(k_c[k_model_ind]) +'_'+ '.png'

        plt.savefig(figname)
        
plt.show()

