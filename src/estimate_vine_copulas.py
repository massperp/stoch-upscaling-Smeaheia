import pyvinecopulib as pv
import numpy as np

k_model = 'k_model_5'

pathName = "/Users/pepe/Documents/FRISK/Python_scripts/smeaheia-test/copula_data/" + k_model + "/"


logK = np.loadtxt(pathName + "K_col_log_samples.txt")
logPc = np.loadtxt(pathName + "Pc_col_log_samples.txt")
logS = np.loadtxt(pathName + "S_col_log_samples.txt")
logKw = np.loadtxt(pathName + "Kw_col_log_samples.txt")
Knw = np.loadtxt(pathName + "Knw_col_samples.txt")


data_samples = np.array([logK - logK.mean(), logPc - logPc.mean(), logS - logS.mean(), logKw - logKw.mean(), (Knw-Knw.mean())/np.std(Knw)])
# Below is the actual copula fit. This takes time
# A possibility is to not use all 10000 data points, but something like the first 1000 samples of data_samples
u = pv.to_pseudo_obs(data_samples.T)
cop = pv.Vinecop(data=u)
print(cop)

# Save copula tree structure to file
pv.Vinecop.to_json(cop, pathName + "copula_tree.txt")
