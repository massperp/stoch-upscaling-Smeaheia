
import numpy as np
import pyvinecopulib as pv
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.interpolate import interp1d

    

def upscaled_Vette_flow_functions(pathName, u_samples_ind, k_model, s_disc_mod):

    #print("k_model: ", k_model)
    N_samp, num_dim = np.shape(u_samples_ind)

    s_disc = np.loadtxt(pathName + "s_disc.txt")
    n_x = len(s_disc)
    logK_ref_data = np.loadtxt(pathName + "K_col_log_samples.txt")
    logPc_ref_data = np.loadtxt(pathName + "Pc_col_log_samples.txt")
    logS_ref_data = np.loadtxt(pathName + "S_col_log_samples.txt")
    logKw_ref_data = np.loadtxt(pathName + "Kw_col_log_samples.txt")
    Knw_ref_data = np.loadtxt(pathName + "Knw_col_samples.txt")
    
    logS_all = np.loadtxt(pathName + "S_log_samples.txt").T
    logKw_all = np.loadtxt(pathName + "Kw_log_samples.txt").T
    
    N_samp_ref = len(logK_ref_data)
    #n_x, N_samp_ref = logKw_all.shape

    # For copula evaluations:
    data_samples = np.array([logK_ref_data - logK_ref_data.mean(), logPc_ref_data - logPc_ref_data.mean(), logS_ref_data - logS_ref_data.mean(), logKw_ref_data - logKw_ref_data.mean(), (Knw_ref_data-Knw_ref_data.mean())/np.std(Knw_ref_data)])

    if k_model == 1:
        mean_logK = 3.1251
        a_Pc = 8.44440405
        b_Pc = -1.492537
        a_Knw = -5.6038
        b_Knw = 1.6690047193430988
        # Index for size 21 s_disc
        s_nw_ind = 16 #0.0550
        if s_disc_mod == 'lin':
            s_nw_ind = 1
    
    if k_model == 3:
        mean_logK = -1.5215
        a_Pc = 9.05254
        b_Pc = -1.492537
        a_Knw = -28.1569
        b_Knw = 20.81520434713147
        # Index for size 21 s_disc
        s_nw_ind = 12 #11
        if s_disc_mod == 'lin':
            s_nw_ind = 0
    
    if k_model == 4:
        mean_logK = -3.1816
        a_Pc = 9.265245165846324
        b_Pc = -1.4925373134327877
    
        a_Knw = -50.948265146681564
        b_Knw = 50.948265146681564
        s_nw_ind = 11
        if s_disc_mod == 'lin':
            s_nw_ind = 0

    # Compute emp cumulative dist functions for the reference data for S and Kw
    ref_ECDF_logS = ECDF(logS_ref_data)
    ref_ECDF_logKw = ECDF(logKw_ref_data)

    ind_Kw_ref = np.argsort(logKw_ref_data)
    ind_S_ref = np.argsort(logS_ref_data)

    ind_S_loc = np.zeros((N_samp_ref, n_x))
    ind_S_loc_inv = np.zeros((N_samp_ref, n_x))
    ind_Kw_loc = np.zeros((N_samp_ref, n_x))
    ind_Kw_loc_inv = np.zeros((N_samp_ref, n_x))

    logS_all_ref_sorted = np.zeros((n_x,N_samp_ref))
    logKw_all_ref_sorted = np.zeros((n_x,N_samp_ref))
    for i in range(n_x):

        ind_S_loc[:,i] = np.argsort(logS_all[i,ind_S_ref])
        ind_S_loc_inv[:,i] = np.argsort(ind_S_loc[:,i])
 
        ind_Kw_loc[:,i] = np.argsort(logKw_all[i,ind_Kw_ref])
        ind_Kw_loc_inv[:,i] = np.argsort(ind_Kw_loc[:,i])
        
        logS_all_ref_sorted[i,:] = logS_all[i,np.int_(ind_S_loc[:,i])]
        logKw_all_ref_sorted[i,:] = logKw_all[i,np.int_(ind_Kw_loc[:,i])]

            
        
    # Instantiate copula object from precomputed model saved to file
    copula = pv.Vinecop(filename = pathName + 'copula_tree.txt', check=True)



    # Sample copula via inverse Rosenblatt of independent uniform samples
    cop_samp_from_ind = copula.inverse_rosenblatt(u_samples_ind)

    # Transform back simulations to the original scale
    copula_samples_from_unif = np.asarray([np.quantile(data_samples[i,:], cop_samp_from_ind[:, i]) for i in range(0, num_dim)])

    logS_cop_samp_unsorted = np.tile(np.mean(logS_ref_data),(1,N_samp)) + copula_samples_from_unif[2,:]
    logKw_cop_samp_unsorted = np.tile(np.mean(logKw_ref_data),(1,N_samp)) + copula_samples_from_unif[3,:]

    # EmpCDF of ref dist applied to copula samples
    cop_ECDF_logS_U  = np.zeros((1,N_samp))
    cop_ECDF_logKw_U  = np.zeros((1,N_samp))

    cop_ECDF_logS_U = np.reshape(ref_ECDF_logS(logS_cop_samp_unsorted),(1,N_samp))
    cop_ECDF_logKw_U = np.reshape(ref_ECDF_logKw( logKw_cop_samp_unsorted),(1,N_samp))

    model_logK_cop_samples = np.zeros((n_x,N_samp))
    model_logPc_cop_samples = np.zeros((n_x,N_samp))
    model_logS_cop_samples = np.zeros((n_x,N_samp))
    model_logKw_cop_samples = np.zeros((n_x,N_samp))
    model_Knw_cop_samples = np.zeros((n_x,N_samp))
    model_Knw_cop_samples[0:s_nw_ind+1,:] = 1 + a_Knw*np.tile(np.c_[s_disc[0:s_nw_ind+1]], (1,N_samp)) + b_Knw*np.tile(np.c_[s_disc[0:s_nw_ind+1]],(1,N_samp))*np.maximum(np.tile(copula_samples_from_unif[4,:],(s_nw_ind+1,1)), np.tile(np.c_[-np.divide(1+a_Knw*s_disc[0:s_nw_ind+1], b_Knw*s_disc[0:s_nw_ind+1])], (1,N_samp)))
    

 
    model_logK_cop_samples = np.tile(mean_logK,(n_x, N_samp)) + np.tile(copula_samples_from_unif[0,:],(n_x,1))
    
    model_logPc_cop_samples = a_Pc + b_Pc*np.tile(np.c_[np.log(s_disc)],(1,N_samp)) + np.tile(copula_samples_from_unif[1,:],(n_x,1))
    
    for i in range(n_x):
    
        loc_ECDF_logS = ECDF(logS_all_ref_sorted[i,:])
        slope_changes_logS = sorted(set(logS_all_ref_sorted[i,:]))
    
        sample_edf_values_at_slope_changes_logS = [ loc_ECDF_logS(item) for item in slope_changes_logS]
    
        loc_ECDF_logKw = ECDF(logKw_all_ref_sorted[i,:])
        slope_changes_logKw = sorted(set(logKw_all_ref_sorted[i,:]))
    
        sample_edf_values_at_slope_changes_logKw = [ loc_ECDF_logKw(item) for item in slope_changes_logKw]
        if i<n_x-1:
            inv_loc_ECDF_logS = interp1d(sample_edf_values_at_slope_changes_logS, slope_changes_logS, fill_value="extrapolate")
            
            model_logS_cop_samples[i:i+1,:] = inv_loc_ECDF_logS(ind_S_loc_inv[np.int_(np.rint(cop_ECDF_logS_U*N_samp_ref)),i]*1/N_samp_ref)#.T
        
            inv_loc_ECDF_logKw = interp1d(sample_edf_values_at_slope_changes_logKw, slope_changes_logKw, fill_value="extrapolate")
                
            model_logKw_cop_samples[i:i+1,:] = inv_loc_ECDF_logKw(ind_Kw_loc_inv[np.int_(np.rint(cop_ECDF_logKw_U*N_samp_ref)),i]*1/N_samp_ref)#.T

    return model_logK_cop_samples, model_logPc_cop_samples, model_logS_cop_samples, model_logKw_cop_samples, model_Knw_cop_samples
