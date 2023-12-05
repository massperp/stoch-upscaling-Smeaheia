"""
Wrapper script to run Cases I-VI with OPM
"""

from mako.template import Template
from joblib import Parallel, delayed
import time
import os
import sys
import numpy as np
from opm.io.ecl import ESmry
import relperm as rp
from scipy.stats import lognorm
from upscaled_flow_functions import upscaled_Vette_flow_functions as upscaled_Vette


def writePermeability(numXY, perms, samp_ind):
  layers = [7, 1, 5, 1, 5, 8]
  permAllL = []
  for layer, perm in zip(layers, perms):
   for i in range(0,layer):
     permAllL.append(perm)

  fn3 = "input/smeaheia_" + str(samp_ind) +".perm"
  g = open(fn3, "w")
  g.write("PERMX\n")
  for perm in permAllL:
    g.write( str(numXY) + "*" + str(perm) + "\n")
  g.write("/\n")
  g.write("PERMY\n")
  for perm in permAllL:
    g.write( str(numXY) + "*" + str(perm) + "\n")
  g.write("/\n")
  g.write("PERMZ\n")
  for perm in permAllL:
    g.write( str(numXY) + "*" + str(perm*0.1) + "\n")
  g.write("/\n")
  g.close()
  
def flow(samp_ind, K1, K2, Kb_aKm, perms, level):

  os.system("python3 vette.py " + str(samp_ind) + " " + str(Kb_aKm) + " " + str(level))
  
  mytemplate = Template(filename = 'SMEAHEIA_BOX.mako')
  
  
  well_I = [7,14,21]
  #well_I = [12,24,36]
  fip_num = [53460, 105300, 157140]
  poro_I = [33, 65, 97]
  numXY =[1980, 3900, 5820]
  var = {'K1' : K1, 'K2' : K2, 'ind' : samp_ind, 'WELL_I' : well_I[level], 'FIPNUM_NUM' : fip_num[level], 'PORO_I': poro_I[level], 'LEVEL' : level }
  FilledTemplate = mytemplate.render(**var)
  casename = "SMEAHEIA_BOX" + "_" + str(samp_ind)
  with open("input/" + casename + ".DATA", 'w') as f:
    f.write(FilledTemplate)

  writePermeability(numXY[level], perms, samp_ind)
  
  os.system("flow input/" + casename + ".DATA --linear-solver-reduction=1e-3 --tolerance-cnv-relaxed=1 --output-dir=output > " + "output/" + casename + ".out")

def anqr(samp_ind):
    casename = "output/SMEAHEIA_BOX_" + str(samp_ind)
    smry = ESmry(casename + ".SMSPEC")
    anqr1 = smry["ANQR:1"]
    anqr2 = smry["ANQR:2"]
    return [anqr1[-1],anqr2[-1]]

def anqp(samp_ind):
    casename = "output/SMEAHEIA_BOX_" + str(samp_ind)
    smry = ESmry(casename + ".SMSPEC")
    anqr1 = smry["ANQP:1"]
    anqr2 = smry["ANQP:2"]
    return [anqr1[-1], anqr2[-1]]

def anqt(samp_ind):
    casename = "output/SMEAHEIA_BOX_" + str(samp_ind)
    smry = ESmry(casename + ".SMSPEC")
    anqt1 = smry["ANQT:1"]
    anqt2 = smry["ANQT:2"]
    return [anqt1[-1], anqt2[-1]]
    
def rgip(samp_ind):
    casename = "output/SMEAHEIA_BOX_" + str(samp_ind)
    smry = ESmry(casename + ".SMSPEC")
    rgip1 = smry["RGIP:2"]
    rgip2 = smry["RGIP:3"]
    return [rgip1[-1], rgip2[-1]]

def writeRelativePermeability(data):

  [nsamples, ndisc] = np.shape(data[0])
  pc_entry = 3.5e-2
  exp = 0.67
  sor = 0.11
  maxkrg = 0.38
  
  for i in range(0,nsamples):
    fn = "input/vette_" + str(i) + ".satfun"
    f = open(fn, "w")
    f.write("SOF2\n")
    rp.writeSOF2(f, sor, exp)
    rp.writeSOF2_F(f, data[2][i], data[3][i], sor)

    f.write("\n")
    f.write("SGFN\n")
    rp.writeSGFN(f, sor, maxkrg, pc_entry, exp)
    rp.writeSGFN_F(f, data[2][i], data[4][i], data[1][i], sor, maxkrg)
    f.write("\n")
    f.close()

def writeRelativePermeabilityDefault(nsamples):
  pc_entry = 3.5e-2
  exp = 0.67
  sor = 0.11
  maxkrg = 0.38
  for i in range(0,nsamples):
    fn = "input/vette_" + str(i) + ".satfun"
    f = open(fn, "w")
    f.write("SOF2\n")
    rp.writeSOF2(f, sor, exp)
    rp.writeSOF2(f, sor, exp)
    f.write("\n")
    f.write("SGFN\n")
    rp.writeSGFN(f, sor, maxkrg, pc_entry, exp)
    rp.writeSGFN(f, sor, maxkrg, pc_entry, exp)
    f.write("\n")
    f.close()
    
  
def opm_test_function(CaseId, k_model_nr, pathName, xi):
  N_disc = 21
  set_n_jobs = 20

  n_samp, num_dim = np.shape(xi)
  s_disc_mod = 'log'

  # Sample copula via inverse Rosenblatt
  if n_samp > 0:
    
    # Parameters for Vette permeability
    if k_model_nr == 1:
        
        # Parameters for Troll permeability
        ln_shape, ln_loc, ln_scale = 0.5941719184341341, 1.744368631343102, 28.619257226121213
        K_Troll_mean = 35.724683018852446
        
        
    if k_model_nr == 3:
        # Parameters for Troll permeability
        ln_shape, ln_loc, ln_scale = 1.2334753451655716, 0.006524801945375478, 1.392124514043006
        K_Troll_mean = 3.025167
    
    if k_model_nr == 4:
        # Parameters for Troll permeability
        ln_shape, ln_loc, ln_scale = 1.4665525482014543, 0.0015096565275845822, 0.520655702642135
        K_Troll_mean = 1.867198573120268
    
        
    if (not os.path.exists("input")):
      os.mkdir("input")
    if (not os.path.exists("output")):
      os.mkdir("output")
      
    if CaseId > 2:

      logK_copula, logPc_copula, logS_copula, logKw_copula, Knw_copula = upscaled_Vette(pathName, xi[:,0:5], k_model_nr, s_disc_mod)
            
      copula_samples_physical = np.array([np.exp(logK_copula).T, np.exp(logPc_copula).T, np.exp(logS_copula).T, np.exp(logKw_copula).T, Knw_copula.T])
      
    # Write tables and adjust data to OPM flow
    ids = []
    cases = []
    ind = 0
    level = 0 # grid level (0 is without refinements)
    
    Kb_aKm = -1.0  # no deformation bands
       
    if (CaseId == 1 or CaseId == 2):
      # Parameters for Vette permeability
      if k_model_nr == 1:
        ln_shape_Vette, ln_loc_Vette, ln_scale_Vette = 0.3091053758612104, 0.0, 22.66301369533692
      
      if k_model_nr == 3:
        ln_shape_Vette, ln_loc_Vette, ln_scale_Vette = 0.8594025002008879, 0.0, 0.2170569658583355
      
      if k_model_nr == 4:
        ln_shape_Vette, ln_loc_Vette, ln_scale_Vette = 1.0661701511114192, 0.0, 0.04152786599598963
        
      K = lognorm.isf(xi[:,0], s=ln_shape_Vette, loc=ln_loc_Vette, scale=ln_scale_Vette)
          
      K = np.tile(np.c_[K],(1,N_disc))
      
      K_troll = K_Troll_mean
      
      
      logPc_all = np.loadtxt(pathName + "Pc_log_samples.txt").T
      logS_all = np.loadtxt(pathName + "S_log_samples.txt").T
      logKw_all = np.loadtxt(pathName + "Kw_log_samples.txt").T
      Knw_all = np.loadtxt(pathName + "Knw_samples.txt").T
        
      Pc_mean = np.mean(np.exp(logPc_all),axis=1)
      S_mean = np.mean(np.exp(logS_all),axis=1)
      Kw_mean = np.mean(np.exp(logKw_all),axis=1)
      Knw_mean = np.mean(Knw_all,axis=1)
        
     
      writeRelativePermeability(np.array([K, np.tile(Pc_mean,(n_samp,1)), np.tile(S_mean,(n_samp,1)), np.tile(Kw_mean,(n_samp,1)), np.tile(Knw_mean,(n_samp,1))]))
    # Case I
    if CaseId == 1:
      perms =[1000, 50, 1000, 50, 850, 25]
      for k1 in K:
        cases.append([ind, k1[0], K_troll, Kb_aKm, perms, level])
        ids.append(ind)
        ind = ind + 1
 
    # Case II
    elif (CaseId == 2):

      # Means and std's of the six layer perms
      perm_means = np.array([1000, 50, 1000, 50, 850, 25])
      #perm_std = 0.1*perm_means
      perm_std = np.array([100, 100, 100, 100, 100, 100])
      #Compute distribution parameters corresponding to the means and stds
      mu_layers = np.log(perm_means**2/(perm_means**2 + perm_std**2)**0.5)
      sigma_layers = (np.log(1+perm_std**2/perm_means**2))**0.5

      K_layers = lognorm.isf(xi[:,1:7], s=sigma_layers, scale=np.exp(mu_layers))

      for k1, perm in zip(K, K_layers):
        cases.append([ind, k1[0], K_troll, Kb_aKm, perm, level])
        ids.append(ind)
        ind = ind + 1
  
    # Case III
    elif (CaseId == 3):
      K = copula_samples_physical[0]
      K_troll = K_Troll_mean
      perms =[1000, 50, 1000, 50, 850, 25]
      writeRelativePermeability(copula_samples_physical)
      for k1 in K:
        cases.append([ind, k1[0], K_troll, Kb_aKm, perms, level])
        ids.append(ind)
        ind = ind + 1
 
    # Case IV
    elif (CaseId == 4):
      K = copula_samples_physical[0]
      K_Troll_fitted_ln = lognorm.isf(xi[:,5], ln_shape, ln_loc, ln_scale)
      
      perms =[1000, 50, 1000, 50, 850, 25]
      writeRelativePermeability(copula_samples_physical)

      for k1, k_troll in zip(K, K_Troll_fitted_ln):
        cases.append([ind, k1[0], k_troll, Kb_aKm, perms, level])
        ids.append(ind)
        ind = ind + 1

    # Case V
    elif (CaseId == 5):
      K = copula_samples_physical[0]
      K_troll = K_Troll_mean

      # Means and std's of the six layer perms
      perm_means = np.array([1000, 50, 1000, 50, 850, 25])
      perm_std = np.array([100, 100, 100, 100, 100, 100])
      #Compute distribution parameters corresponding to the means and stds
      mu_layers = np.log(perm_means**2/(perm_means**2 + perm_std**2)**0.5)
      sigma_layers = (np.log(1+perm_std**2/perm_means**2))**0.5
      K_layers = lognorm.isf(xi[:,5:11], s=sigma_layers, scale=np.exp(mu_layers))
      
      writeRelativePermeability(copula_samples_physical)

      for k1, perm in zip(K, K_layers):
        cases.append([ind, k1[0], K_troll, Kb_aKm, perm, level])
        ids.append(ind)
        ind = ind + 1
  
    # Case VI
    elif (CaseId == 6):
      K = copula_samples_physical[0]
      # Means and std's of the six layer perms
      perm_means = np.array([1000, 50, 1000, 50, 850, 25])
      #perm_std = 0.1*perm_means
      perm_std = np.array([100, 100, 100, 100, 100, 100])
      #Compute distribution parameters corresponding to the means and stds
      mu_layers = np.log(perm_means**2/(perm_means**2 + perm_std**2)**0.5)
      sigma_layers = (np.log(1+perm_std**2/perm_means**2))**0.5
      K_layers = lognorm.isf(xi[:,5:11], s=sigma_layers, scale=np.exp(mu_layers))
      K_Troll_fitted_ln = lognorm.isf(xi[:,5], ln_shape, ln_loc, ln_scale)
    
      writeRelativePermeability(copula_samples_physical)

      for k1, perm, k_troll in zip(K, K_layers, K_Troll_fitted_ln):
        cases.append([ind, k1[0], k_troll, Kb_aKm, perm, level])
        ids.append(ind)
        ind = ind + 1

        

    Parallel(n_jobs=set_n_jobs)(delayed(flow)(ind, k1, k2, Kb_aKm, perms, level) for ind, k1, k2, Kb_aKm, perms, level in cases)

    anqrs = []
    anqps = []
    anqts = []
    rgips = []
    rgips2 = []
    for ind in ids:
      anqrs.append(anqr(int(ind)))
      anqps.append(anqp(int(ind)))
      rgips.append(rgip(int(ind))[0]/1000)
      rgips2.append(rgip(int(ind))[1]/1000)
      anqts.append(anqt(int(ind))[1])

    return np.array(rgips)

  else:
    rgips = []
    return np.array(rgips)
