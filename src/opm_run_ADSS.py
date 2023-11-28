
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
    #all_strata = np.zeros(N_rep)
    # main loop to repeat adaptive stratification Nrep times
    # rg = np.random.Generator(np.random.MT19937(42))
    for j_rep in range(N_rep):
        # We use a rng with a seeded MT to enable reproducibility
        # We need to put this in the loop, to get the same numbers.
        # rg = np.random.Generator(np.random.MT19937(42))
        
        try:
            t_calculation = time.time()
            estimator = AdaptiveStratification(test_fun, N_dim, N_max,
                                               SR_const, alpha, dynamic=False,
                                               type_strat=stype, rand_gen=None)
            QoI_vec[j_rep], all_strata, N_strat[j_rep], QoI_var_vec[j_rep] = (
                estimator.solve())  # solve_vis_steps(), solve_vis_result()
            print("type(all_strata): ", type(all_strata))
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

    t = time.time()

    alpha = 0.5  # alpha: proportion of samples allocated optimally
    N_max = int(300)  #int(1e3)  # numbers of max samples
    SR_const = 50  # increase per adaptation iteration
    
    CaseIds = [1,2,3,4,5,6]
    k_models = [3,4] #[1,3] #, 3]
    
    N_rep = 1
    for k_model_nr in k_models:
        
    
        pathName = "/home/AD.NORCERESEARCH.NO/pepe/adaptive-stratification-python/adaptive_stratification/examples/opm/stoch-upscaling-Smeaheia/copula_data/k_model_" + str(k_model_nr) + "/"
        pathName = "copula_data/k_model_" + str(k_model_nr) + "/"

    
        for CaseId in CaseIds:
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
            
            MC_samp = test_fun_part(np.random.default_rng().uniform(0, 1, (N_max,N_dim)))
    
            fname = ('results/case_' + str(CaseId) + '/MC_Nsamples_' + str(N_max) +'_kmod_' + str(k_model_nr) + '.txt')
            if (not os.path.exists('results/case_' + str(CaseId) )):
                os.makedirs('results/case_' + str(CaseId) )
            np.savetxt(fname, MC_samp, fmt="%s")
            print('Done with case ', str(CaseId))
    

            (QoI_vec, QoI_var_vec, N_strat, all_strata) = repeated_runs(test_fun_part, N_dim, N_max, SR_const, alpha, stype, N_rep)
            
            fname = ('results/case_' + str(CaseId) + '/QoI_Nmax_' + str(N_max) +'_alpha_'+str(alpha)+'_kmod_' + str(k_model_nr) +'_rep_' + str(N_rep) + '.txt')
            if (not os.path.exists('results/case_' + str(CaseId) )):
                os.makedirs('results/case_' + str(CaseId) )
            np.savetxt(fname, QoI_vec, fmt="%s")
        
            elapsed = time.time() - t
            print(f'Total elapsed time is {elapsed} seconds.')
            print('Done with case ', str(CaseId))

    
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
