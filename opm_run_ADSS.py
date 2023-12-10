"""
Script for running adaptive stratified sampling and standard Monte Carlo with OPM Flow.
The six test cases presented in "Copula modeling and uncertainty propagation in field-scale simulation of CO2 fault leakage" are defined in opm_wrapper.py. The user parameters alpha, N_max, and SR_const introduced below in main can be varied for, e.g., convergence studies
"""

import sys
import time
import subprocess  # shlex
import numpy as np
from typing import Tuple
import os
from adaptive_stratification.stratification import AdaptiveStratification
from opm_wrapper import opm_test_function as test_fun
import functools

def repeated_runs(test_fun, N_dim, N_max, SR_const, alpha, stype, N_rep) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Repeat stratification estimation N_rep times for reliability.

    Args:
        test_fun: function handle to quantity of interest.
        N_dim: number of stochastic variables, transformed to unit
               uniforms.
        N_max: total number of samples to be distributed.
        SR_const: average number of samples to be added to each
                           stratum in every iteration.
        alpha: proportion of samples to be distributed optimally.
        stype: decides the shape of the stratum to use.
        N_rep: Number of repetitions

    Returns:
        A tuple containing N_rep estimates and variances using a simple
        Monte Carlo, the N_rep estimates and variances of the adaptive
        stratification algorithm and the amounts of strata.
    """

    QoI_vec = np.zeros(N_rep)
    QoI_var_vec = np.zeros(N_rep)
    N_strat = np.zeros(N_rep, dtype='uint32')

    for j_rep in range(N_rep):
        
        try:
            t_calculation = time.time()
            estimator = AdaptiveStratification(test_fun, N_dim, N_max,
                                               SR_const, alpha, dynamic=False,
                                               type_strat=stype, rand_gen=None)
            QoI_vec[j_rep], all_strata, N_strat[j_rep], QoI_var_vec[j_rep] = (
                estimator.solve())  # solve_vis_steps(), solve_vis_result()
            elapsed_calculation = time.time() - t_calculation
            print(f'Elapsed time is {elapsed_calculation} seconds for this run.')
        except (Exception):
            (QoI_vec[j_rep], N_strat[j_rep],
             QoI_var_vec[j_rep]) = np.nan, 0, np.nan
            raise

    return (QoI_vec, QoI_var_vec, N_strat, all_strata)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        stype = sys.argv[1]
    elif len(sys.argv) == 1:
        stype = input("Choose Stratum Type: hyperrect or simplex? ")
    else:
        raise Exception('Programs expect zero or one command-line argument')


    if stype == 'hyperrect' or stype == 'simplex':
        fbase = f'./results/{stype}_'  # rot_p1_
    else:
        raise Exception('Please use hyperrect or simplex')

    #t = time.time()

    # The user parameters alpha, N_max, and SR_const can be varied for, e.g., convergence studies
    alpha = 0.5  # alpha: proportion of samples allocated optimally
    N_max = int(500)  #int(1e3)  # numbers of max samples
    SR_const = 50  # increase per adaptation iteration
    
    CaseIds = [1,2,3,4,5,6]
    
    # The speedup estimates can be replaced by empirical estimates based on repeated sampling by setting number of repetitions N_rep > 1.
    N_rep = 1
    
    k_model_order = [4,3,1]
    k_c = ['1e-4mD','1e-3mD', '1mD']

    for k_model_ind in range(len(k_model_order)):
        k_model_nr = k_model_order[k_model_ind]
        print("SGR-perm model with k_c = ", k_c[k_model_ind])
        pathName = "copula_data/k_model_" + str(k_model_nr) + "/"

        for CaseId in CaseIds:
            # N_dim: number of uniform input random variables. Required for ADSS.
            if CaseId==1:
                N_dim=1
            if CaseId==2:
                N_dim=7
            if CaseId==3:
                N_dim=5
            if CaseId==4:
                N_dim=6
            if CaseId==5:
                N_dim=11
            if CaseId==6:
                N_dim=12

           
            test_fun_part = functools.partial(test_fun, CaseId, k_model_nr, pathName)
            
            # Standard Monte Carlo simulation
            MC_samp = test_fun_part(np.random.default_rng().uniform(0, 1, (N_max,N_dim)))
    
            fname = ('results/case_' + str(CaseId) + '/MC_Nsamples_' + str(N_max) +'_kmod_' + str(k_model_nr) + '.txt')
            if (not os.path.exists('results/case_' + str(CaseId) )):
                os.makedirs('results/case_' + str(CaseId) )
            np.savetxt(fname, MC_samp, fmt="%s")
            print('Finished SMC for case ', str(CaseId))
    
            # Run adaptive stratified sampling
            (QoI_vec, QoI_var_vec, N_strat, all_strata) = repeated_runs(test_fun_part, N_dim, N_max, SR_const, alpha, stype, N_rep)
            
            fname = ('results/case_' + str(CaseId) + '/QoI_Nmax_' + str(N_max) +'_alpha_'+str(alpha)+'_kmod_' + str(k_model_nr) +'_rep_' + str(N_rep) + '.txt')
            if (not os.path.exists('results/case_' + str(CaseId) )):
                os.makedirs('results/case_' + str(CaseId) )
            np.savetxt(fname, QoI_vec, fmt="%s")
        
            #elapsed = time.time() - t
            #print(f'Total elapsed time is {elapsed} seconds.')
            print('Finished ADSS for case ', str(CaseId))
            
            # Print stratification results to file.
            p_all_strata = []
            N_all_strata = []
            pathN = "stratifications/strat_case_"+str(CaseId) + "_k" + str(k_model_nr) + "_N" + str(N_max)
            if (not os.path.exists(pathN)):
                os.makedirs(pathN)
  
            for i in range(N_strat[0]):
        
                fname = (pathN + '/Strat_num_f_samples' + str(i) + '.txt')
                np.savetxt(fname, all_strata[i].f_samples, fmt="%s")
                fname = (pathN + '/Strat_num_u_samples' + str(i) + '.txt')
                np.savetxt(fname, all_strata[i].samples, fmt="%s")
        
                p_all_strata.append(all_strata[i].p)
                N_all_strata.append(all_strata[i].N_samples)
        
            with open(pathN+"/summary_ADSS.txt", "w") as f:
            #f.write("QoI: " + str(QoI) + " QoI_var: " + str(QoI_var))
                f.write(str(QoI_vec) + "\n")
                f.write(str(QoI_var_vec) + "\n")
                np.savetxt(f, p_all_strata, fmt='%.5f')
                np.savetxt(f, N_all_strata, fmt='%.5f')
