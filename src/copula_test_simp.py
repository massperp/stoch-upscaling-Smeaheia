import sys
import time
import subprocess  # shlex
import numpy as np
from scipy.stats import uniform
import itertools
from typing import Tuple
import pyvinecopulib as pv
import os
from adaptive_stratification.stratification import AdaptiveStratification
import functools
from matplotlib import pyplot as plt

from statsmodels.distributions.empirical_distribution import ECDF
from scipy.interpolate import interp1d



def repeat_product(x, y):
    #return np.transpose([np.tile(x, len(y)),
    #                        np.repeat(y, len(x))])
    return np.transpose([np.repeat(x, len(y)),
                            np.tile(y, len(x))])


def copula_test_functions_2(data_samples: np.ndarray, copula, s_disc, xi: np.ndarray):
  #print("xi: ",xi)
  n_samp, num_dim = np.shape(xi)
  num_dim_copula = 5
  
  N_disc = len(s_disc)

  # Sample copula via inverse Rosenblatt
  cop_samp_from_ind = copula.inverse_rosenblatt(xi[:,0:5])
  # Transform back simulations to the original scale
  Y_sim_from_unif = np.asarray([np.quantile(data_samples[i,:], cop_samp_from_ind[:, i]) for i in range(0, num_dim_copula)])
  
  
  K_copula, Pc_copula, S_copula, Kw_copula, Knw_copula = upscaled_param_eval(s_disc, k_model_nr, ind_S_loc_inv, ind_Kw_loc_inv, logS_ref_data, logKw_ref_data, logS_all, logKw_all, Y_sim_from_unif)
  copula_samples_physical = np.array([K_copula, Pc_copula, S_copula, Kw_copula, Knw_copula])
      
  return copula_samples_physical


def upscaled_param_eval(s_disc, k_model_nr, ind_S_loc_inv, ind_Kw_loc_inv, logS_ref_data, logKw_ref_data, logS_all, logKw_all, copula_samples: np.ndarray):
    """
    Evaluation of the five stochastic physical parameters (permeability, capillary pressure, saturation, wetting relative permeability, nonwetting relative permeability)
    """
    N_disc = len(s_disc)
    N_dim, N_samples = np.shape(copula_samples)
    print("np.shape(copula_samples): ", np.shape(copula_samples))
    _, N_samp_ref = logKw_all.shape
    
   
    ref_ECDF_logS = ECDF(logS_ref_data)
    ref_ECDF_logKw = ECDF(logKw_ref_data)

    ind_Kw_ref = np.argsort(logKw_ref_data)
    ind_S_ref = np.argsort(logS_ref_data)

        

    logS_cop_samp_unsorted = np.tile(np.mean(logS_ref_data),(1,N_samples)) + copula_samples[2,:]
    logKw_cop_samp_unsorted = np.tile(np.mean(logKw_ref_data),(1,N_samples)) + copula_samples[3,:]


    cop_ECDF_logS_U  = np.zeros((1,N_samples))
    cop_ECDF_logKw_U  = np.zeros((1,N_samples))

    cop_ECDF_logS_U = np.reshape(ref_ECDF_logS(logS_cop_samp_unsorted),(1,N_samples))
    cop_ECDF_logKw_U = np.reshape(ref_ECDF_logKw( logKw_cop_samp_unsorted),(1,N_samples))

    model_logK_cop_samples = np.zeros((N_disc,N_samples))
    model_logPc_cop_samples = np.zeros((N_disc,N_samples))
    model_logS_cop_samples = np.zeros((N_disc,N_samples))
    model_logKw_cop_samples = np.zeros((N_disc,N_samples))
    model_Knw_cop_samples = np.zeros((N_disc,N_samples))

    
    # This could possibly be read in via json, but the values as well as the functions above are quite model specific, so it might only give a false sense of generality
    
    if k_model_nr == 1:

        mean_logK = 3.1251
        a_Pc = 8.44440405
        b_Pc = -1.492537
        a_Knw = -5.6038
        b_Knw = 1.6690047193430988
        # Index for size 21 s_disc
        s_nw_ind = 16 #0.0550
    

    if k_model_nr == 3:
        mean_logK = -1.5215
        a_Pc = 9.05254
        b_Pc = -1.492537
        a_Knw = -28.1569
        b_Knw = 20.81520434713147
        # Index for size 21 s_disc
        s_nw_ind = 12 #11
    
    if k_model_nr == 4:
        mean_logK = -3.1816
        a_Pc = 9.265245165846324
        b_Pc = -1.4925373134327877
    
        a_Knw = -50.948265146681564
        b_Knw = 50.948265146681564
        s_nw_ind = 11
    
    if k_model_nr == 5:
        mean_logK = -4.861098435299905
        a_Pc = 9.469112241474773
        b_Pc = -1.492537313432776
    
        a_Knw = -88.16181836208175
        b_Knw = 88.16181836208175
        s_nw_ind = 10
    
    K_samples = np.zeros((N_samples, N_disc))
    Pc_samples = np.zeros_like(K_samples)
    S_samples = np.zeros_like(K_samples)
    Kw_samples = np.zeros_like(K_samples)
    Knw_samples = np.zeros_like(K_samples)
    
    model_logK_cop_samples = np.tile(mean_logK,(N_disc, N_samples)) + np.tile(copula_samples[0,:],(N_disc,1))
    #model_logPc_cop_samples = a_Pc + b_Pc*np.tile(np.c_[s_disc],(1,N_samples)) + np.tile(copula_samples[1,:],(N_disc,1))
    model_logPc_cop_samples = a_Pc + b_Pc*np.tile(np.c_[np.log(s_disc)],(1,N_samples)) + np.tile(copula_samples[1,:],(N_disc,1))

    model_Knw_cop_samples[0:s_nw_ind+1,:] = 1 + a_Knw*np.tile(np.c_[s_disc[0:s_nw_ind+1]], (1,N_samples)) + b_Knw*np.tile(np.c_[s_disc[0:s_nw_ind+1]],(1,N_samples))*np.maximum(np.tile(copula_samples[4,:],(s_nw_ind+1,1)), np.tile(np.c_[-np.divide(1+a_Knw*s_disc[0:s_nw_ind+1], b_Knw*s_disc[0:s_nw_ind+1])], (1,N_samples)))
    
    for i in range(N_disc):
    
        loc_ECDF_logS = ECDF(logS_all_ref_sorted[i,:])
        slope_changes_logS = sorted(set(logS_all_ref_sorted[i,:]))
    
        sample_edf_va/Users/pepe/Documents/FRISK/Python_scripts/smeaheia-test/opm_run_ADSS.pylues_at_slope_changes_logS = [ loc_ECDF_logS(item) for item in slope_changes_logS]
    
        loc_ECDF_logKw = ECDF(logKw_all_ref_sorted[i,:])
        slope_changes_logKw = sorted(set(logKw_all_ref_sorted[i,:]))
    
        sample_edf_values_at_slope_changes_logKw = [ loc_ECDF_logKw(item) for item in slope_changes_logKw]
        if i<N_disc-1:
            inv_loc_ECDF_logS = interp1d(sample_edf_values_at_slope_changes_logS, slope_changes_logS, fill_value="extrapolate")
        
            model_logS_cop_samples[i:i+1,:] = inv_loc_ECDF_logS(ind_S_loc_inv[np.int_(np.rint(cop_ECDF_logS_U*N_samp_ref)),i]*1/N_samp_ref)#.T
        
            inv_loc_ECDF_logKw = interp1d(sample_edf_values_at_slope_changes_logKw, slope_changes_logKw, fill_value="extrapolate")
                
            model_logKw_cop_samples[i:i+1,:] = inv_loc_ECDF_logKw(ind_Kw_loc_inv[np.int_(np.rint(cop_ECDF_logKw_U*N_samp_ref)),i]*1/N_samp_ref)#.T

    
    K_samples = np.exp(model_logK_cop_samples).T
    Pc_samples = np.exp(model_logPc_cop_samples).T
    S_samples = np.exp(model_logS_cop_samples).T
    Kw_samples = np.exp(model_logKw_cop_samples).T
    Knw_samples = model_Knw_cop_samples.T

    return K_samples, Pc_samples, S_samples, Kw_samples, Knw_samples
    


        
    
if __name__ == "__main__":

    t = time.time()

        
    k_model_nr = 4 #[3,1] #[1, 3]
  
        
    #pathName = "/home/AD.NORCERESEARCH.NO/pepe/adaptive-stratification-python/adaptive_stratification/examples/opm/smeaheia-2/copula_data/k_model_" + str(k_model_nr) + "/"
    pathName = "/Users/pepe/Documents/FRISK/Python_scripts/smeaheia-test/copula_data/k_model_" + str(k_model_nr) + "/"
    logK_ref_data = np.loadtxt(pathName + "K_col_log_samples.txt")
    logPc_ref_data = np.loadtxt(pathName + "Pc_col_log_samples.txt")
    logS_ref_data = np.loadtxt(pathName + "S_col_log_samples.txt")
    logKw_ref_data = np.loadtxt(pathName + "Kw_col_log_samples.txt")
    Knw_ref_data = np.loadtxt(pathName + "Knw_col_samples.txt")

    logPc_all = np.loadtxt(pathName + "Pc_log_samples.txt").T
    logS_all = np.loadtxt(pathName + "S_log_samples.txt").T
    logKw_all = np.loadtxt(pathName + "Kw_log_samples.txt").T
    Knw_all = np.loadtxt(pathName + "Knw_samples.txt").T

    ind_Kw_loc = np.loadtxt(pathName + "ind_Kw_loc.txt")
    ind_S_loc = np.loadtxt(pathName + "ind_S_loc.txt")
    
    ind_Kw_loc_inv = np.loadtxt(pathName + "ind_Kw_loc_inv.txt")
    ind_S_loc_inv = np.loadtxt(pathName + "ind_S_loc_inv.txt")
    
    data_samples = np.array([logK_ref_data - logK_ref_data.mean(), logPc_ref_data - logPc_ref_data.mean(), logS_ref_data - logS_ref_data.mean(), logKw_ref_data - logKw_ref_data.mean(), (Knw_ref_data-Knw_ref_data.mean())/np.std(Knw_ref_data)])
    
    N_samp_ref = 10000
    N_disc = 21
    s_disc = np.logspace(-6,0,N_disc)
    logS_all_ref_sorted = np.zeros((N_disc,N_samp_ref))
    logKw_all_ref_sorted = np.zeros((N_disc,N_samp_ref))
    for i in range(N_disc):
        
        logS_all_ref_sorted[i,:] = logS_all[i,np.int_(ind_S_loc[:,i])]
        logKw_all_ref_sorted[i,:] = logKw_all[i,np.int_(ind_Kw_loc[:,i])]
    
    
    copula = pv.Vinecop(filename='copula_data/k_model_'+str(k_model_nr)+'/copula_tree.txt', check=True)
  

    
    tol = 1e-5
    u_disc = np.linspace(0.0+tol, 1.0-tol, num=51)
    u_disc_2 = np.linspace(0.0+tol, 1.0-tol, num=51)
    u_evals_2D = repeat_product(u_disc, u_disc_2)
    num_evals = 51*51
    
    #print("u_evals_2D ", u_evals_2D)

    test_fun_part_copula = functools.partial(copula_test_functions_2, data_samples, copula, s_disc)
        
    
    
    if (not os.path.exists("strat_test_3")):
        os.mkdir("strat_test_3")
    for ind in itertools.combinations(range(5), 2):
        target_dims = ind
        
        
        if (not os.path.exists('strat_test_3/variables_'+str(ind[0]+1)+str(ind[1]+1)+'_k'+str(k_model_nr))):
            os.mkdir('strat_test_3/variables_'+str(ind[0]+1)+str(ind[1]+1)+'_k'+str(k_model_nr))

        print("Var comb: ", ind)

        u_evals = 0.5*np.ones((num_evals, 5))
        u_evals[:,ind] = u_evals_2D
        
        copula_evals = test_fun_part_copula(u_evals)
        #print("np.shape(copula_evals): ", np.shape(copula_evals))
        
        for m in range(5):
            fname = ('strat_test_3/variables_'+str(ind[0]+1)+str(ind[1]+1)+'_k'+str(k_model_nr)+'/copula_evals_Y_' + str(m+1) + '.txt')
            np.savetxt(fname, copula_evals[m,:], fmt="%s")
        
    
    print("s_disc: ", s_disc)


