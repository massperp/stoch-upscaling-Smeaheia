
"""
Evaluate upscaled flow functions on uniform grid for pairs of input parameters U_i, U_j on [0, 1]^2 while keeping the other three U's at 0.5.
"""
import sys
import numpy as np
from scipy.stats import uniform
import itertools
import os
import functools
from upscaled_flow_functions import upscaled_Vette_flow_functions as upscaled_Vette

def repeat_product(x, y):
    return np.transpose([np.repeat(x, len(y)),
                            np.tile(y, len(x))])

def copula_test_functions(pathName, k_model, s_disc_mod, u_samples: np.ndarray):
  
  logK_copula, logPc_copula, logS_copula, logKw_copula, Knw_copula = upscaled_Vette(pathName, u_samples, k_model, s_disc_mod)
  
  copula_samples_physical = np.array([np.exp(logK_copula).T, np.exp(logPc_copula).T, np.exp(logS_copula).T, np.exp(logKw_copula).T, Knw_copula.T])
      
  return copula_samples_physical

        
if __name__ == "__main__":
        
    k_model_nr = 4 #[3,1] #[1, 3]
    pathName = "copula_data/k_model_" + str(k_model_nr) + "/"
    
    tol = 1e-5
    u_disc = np.linspace(0.0+tol, 1.0-tol, num=51)
    u_disc_2 = np.linspace(0.0+tol, 1.0-tol, num=51)
    u_evals_2D = repeat_product(u_disc, u_disc_2)
    num_evals = 51*51
    
    s_disc_mod = 'log'
    test_fun_part_copula = functools.partial(copula_test_functions, pathName, k_model_nr, s_disc_mod)
        
     
    if (not os.path.exists("flow_fun_evals")):
        os.mkdir("flow_fun_evals")
    for ind in itertools.combinations(range(5), 2):
        target_dims = ind
        
        if (not os.path.exists('flow_fun_evals/variables_'+str(ind[0]+1)+str(ind[1]+1)+'_k'+str(k_model_nr))):
            os.mkdir('flow_fun_evals/variables_'+str(ind[0]+1)+str(ind[1]+1)+'_k'+str(k_model_nr))

        print("Running variable combination: (", ind[0]+1, ",", ind[1]+1, ")")

        u_evals = 0.5*np.ones((num_evals, 5))
        u_evals[:,ind] = u_evals_2D
        
        copula_evals = test_fun_part_copula(u_evals)
        
        for m in range(5):
            fname = ('flow_fun_evals/variables_'+str(ind[0]+1)+str(ind[1]+1)+'_k'+str(k_model_nr)+'/copula_evals_Y_' + str(m+1) + '.txt')
            np.savetxt(fname, copula_evals[m,:], fmt="%s")
