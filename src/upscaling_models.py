
from typing import Dict, Tuple, Union
import numpy as np
from scipy import interpolate
import itertools


def fault_perm(SGR: np.ndarray, data: Dict, k_model: str) -> np.ndarray:
    # Compute permeability from SGR
    num_samples, num_facies = SGR.shape

    # Clay-rich shale and sandstone permeabilities (in mD)
    k_c = data["perm_clay"]
    k_s = data["perm_sandstone"]
    # Compute permeability, Eq. 9 in Bjørnarå et al

    log_perm = SGR*np.log(k_c/k_s) + np.log(k_s)
    perm_fault = np.exp(log_perm)
    
    return perm_fault
    

def harmonic_mean(perm: np.ndarray, heights: np.ndarray) -> np.ndarray:
    # Calculate harmonic mean of permeability.
    
    if np.all(perm):
        harmonic_mean_perm = np.true_divide(
            np.sum(heights, axis=1), np.sum(np.true_divide(heights, perm), axis=1)
        )
    else:
        harmonic_mean_perm = 0

    return harmonic_mean_perm

def arithmetic_mean(perm: np.ndarray, heights: np.ndarray) -> np.ndarray:
    # Calculate arithmetic mean of permeability.
    return np.sum(heights * perm, axis = 1) / np.sum(heights, axis = 1)

    
def fault_heights_SGR(num_facies: int, num_samples: int, data: Dict, data_interpr, init_loc, sect
) -> Tuple[np.ndarray, np.ndarray]:
    # Geometry of facies distribution
    # Return transmisibility from upscaled permeability

    # Calculate mean height
    fault_height = data["fault_height_" + sect]
    mean_facies_height = fault_height / num_facies

    # Create uniformly distributed interior boundaries of fault facies
    # The perturbations are scaled with the mean facies height.
    facies_int_bds = np.tile(
        np.cumsum(np.repeat(mean_facies_height, num_facies - 1)), (num_samples, 1)
    ) + 0.8 * mean_facies_height * np.random.uniform(
        low=-0.5, high=0.5, size=(num_samples, num_facies - 1)
    )
    # heights = np.random.uniform(low=10, high=100, size=(num_facies,1))

    heights = np.diff(
        np.c_[
            np.tile(0, (num_samples, 1)),
            facies_int_bds,
            np.tile(fault_height, (num_samples, 1)),
        ]
    )
    # Sanity check
    assert np.all(heights > 0)
    
    if data_interpr == "depth_dep":
       
        mid_cells = init_loc + np.cumsum(heights, axis=1)-heights/2
     
        SGR_data_Bjornaraa = np.loadtxt("SGR_synth_Bjornaraa.txt")/100

        neg_depths = np.arange(500,2010,10)
        f = interpolate.interp1d(neg_depths, SGR_data_Bjornaraa)
    
        mu_SGR = f(mid_cells)
        sigma_SGR = data["sigma_SGR"]
        
    elif data_interpr == "no_data":
            # Sample SGR for the facies (independently)
            mu_SGR, sigma_SGR = data["mu_SGR"], data["sigma_SGR"]  # mean and standard deviation
        
    SGR = np.random.normal(mu_SGR, sigma_SGR, size=(num_samples, num_facies))
    
    return heights, SGR


##############################

    
def entry_pressure_from_perm(perm: np.ndarray, data: Dict) -> np.ndarray:
    # Get the capillary entry pressure calculated from the perm.
    entry_pressure_sand = data["entry_pressure_sand"]
    mDARCY = 9.869233e-13 * 1e-3
    k_s = data["perm_sandstone"] * mDARCY
    #leveret J
    return entry_pressure_sand * np.sqrt(k_s / perm)

def rel_perm_cap_pressure(
    perm: np.ndarray, heights: np.ndarray, data: Dict, s_disc
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Calculate capillary pressure (measured in Pascal) and relative permeabilities.
    # Use a Brooks-Corey correlation.

    num_samples, num_facies = perm.shape

    # Conversion from mD to SI
    DARCY = 9.869233e-13 * 1e-3
    perm *= DARCY

    # Brooks-Corey
    bc_exponent = data["brooks_corey_exponent"]

    # Maximum relative permeabilities
    kw_max = data.get("kw_max", 1)
    kwn_max = data.get("kwn_max", 1)

    pc_entry = entry_pressure_from_perm(perm, data)

    
    # Aim at an approximately uniform distribution of sampling points in the fine-scale
    # saturation (it will be uniform if the entry pressures are equal for all
    # facies in a sample).
    # To that end, calculate fine-scale capillary pressure with a uniform saturation.
    # NOTE: The cut-off value for saturation in effect will determine a maximum value for
    # capillary pressure.
    # NOTE: Taking the sample-wise max of entry-pressure is necessary to avoid negative
    # capillary pressures in some facies.

    
    Pc = np.outer(pc_entry.min(1), np.power(s_disc, -1 / bc_exponent)).T

    # Storage vectors for saturations and relative permeabilities
    S = np.zeros((num_samples, num_s_points))
    Kw = np.zeros_like(S)
    Knw = np.zeros_like(S)
    #indices = np.zeros_like(S)
    # Upscaled permeability, needed for scaling relative permeabilities. This may be
    # computed elsewhere, but the calculation should be fast.
    upscaled_perm = harmonic_mean(perm, heights)
    #s_temp = np.zeros((num_samples,num_s_points))
    for i, pc in enumerate(Pc):
        # Loop over the preset capillary pressures and:
        # 1) Calculate fine-scale saturation corresponding to a capillary equilibrium
        # 2) Calculate the upscaled saturation by averaging
        # 3) Compute fine-scale relative permeability for the calculated fine-scale saturation
        # 4) Compute upscaled relative permeability by harmonic averages

        # Brooks-Corey relation for permeability
        
        s_fine_raw = np.power(pc_entry / np.atleast_2d(pc.T).T, bc_exponent)
        s_fine = np.clip(s_fine_raw, a_min=0, a_max=1)


        # Arithmetic averaging of the saturation.
        S[:, i] = np.sum(s_fine * heights, axis=1) / heights.sum(axis=1)

        # Fine-scale effective permeability for the water phase
        Kw_fine = kw_max * perm * np.power(s_fine, (2 + 3 * bc_exponent) / bc_exponent)
        Knw_fine = (
            np.power(1 - s_fine, 2)
            * kwn_max
            * perm
            * (1 - np.power(s_fine, (2 + bc_exponent) / bc_exponent))
        )
        # Upscale phase-wise effective permeabilities by harmonic means, rescale to get
        # relative permeabilities
        Kw[:, i] = harmonic_mean(Kw_fine, heights) / upscaled_perm
        Knw[:, i] = harmonic_mean(Knw_fine, heights) / upscaled_perm
        # Sanity checks
        assert np.all(S[:, i] >= 0) and np.all(S[:, i] <= 1)
        
    return S, Kw, Knw, Pc.T  # Revert the transpose


def fault_facies(
    num_facies: int, num_samples: int, data: Dict, data_interpr, s_disc, k_model, sect, init_loc, two_phase: bool = False, harmonic: bool = True
) -> Union[
    np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
]:
    # Upscaling of permeability and two-phase flow functions

    # Facies distribution
    heights, SGR = fault_heights_SGR(num_facies, num_samples, data, data_interpr, init_loc, sect)
   
    # Fine-scale permeability
    fine_scale_perm = fault_perm(SGR, data, k_model)
    # Upscaled permeability
    if (harmonic):
        upscaled_permeability = harmonic_mean(fine_scale_perm, heights)
    else:
        upscaled_permeability = arithmetic_mean(fine_scale_perm, heights)
        
    if two_phase:
        # Two-phase flow functions consistent with the fine-scale permeability
        S, Kw, Knw, Pc = rel_perm_cap_pressure(fine_scale_perm, heights, data, s_disc)
        return upscaled_permeability, S, Kw, Knw, Pc
    else:
        return upscaled_permeability

def Troll_permeability(
    num_facies: int, num_samples: int, data: Dict, data_interpr, s_disc, init_loc_Troll, k_model, harmonic: bool = True
) -> np.ndarray:

    # Facies distribution
    heights, SGR = fault_heights_SGR(num_facies, num_samples, data, data_interpr, init_loc=init_loc_Troll, sect='Troll')
   
    # Fine-scale permeability
    fine_scale_perm = fault_perm(SGR, data, k_model)
    # Upscaled permeability
    if (harmonic):
        upscaled_permeability = harmonic_mean(fine_scale_perm, heights)
    else:
        upscaled_permeability = arithmetic_mean(fine_scale_perm, heights)
        
    return upscaled_permeability


data = {
    "mu_SGR": 0.50,
    "sigma_SGR": 0.14,
    #"perm_clay": 1, #1e-3,
    "perm_sandstone": 1000,
    "fault_height_Vette": 500,
    "fault_height_Troll": 100, #Changed from 10 to 100 by Per, Oct 3 -22
    # Brooks-Corey exponent. TH's backcalculation of CO2DataShare information (0.8)
    "brooks_corey_exponent": 0.67,
    # use 2/3 from "Relative permeability and trapping of CO2 and water in sandstone rocks at reservoir conditions"
    # gives n_brine = 6 and n_co2 = 4 which is ok.
    "entry_pressure_sand": 2.5e3
}


def estimate_coef(x, y):
    # number of observations/points
    n = np.size(x)
 
    # mean of x and y vector
    m_x = np.mean(x)
    m_y = np.mean(y)
 
    # calculating cross-deviation and deviation about x
    SS_xy = np.sum(y*x) - n*m_y*m_x
    SS_xx = np.sum(x*x) - n*m_x*m_x
 
    # calculating regression coefficients
    b_1 = SS_xy / SS_xx
    b_0 = m_y - b_1*m_x
 
    return (b_0, b_1)
    


# k_model_1: perm_clay: 1 mD, and the SGR-k model in Bjørnaraa et al, and the depth data inspired by them
# k_model_2: SGR-k model from Manzocchi et al
# k_model_3: perm_clay: 1e-3 mD, otherwise the SGR-k model in Bjørnaraa et al, and the depth data inspired by them

k_model = 'k_model_3'
s_disc_mod = 'lin'

if k_model == 'k_model_1':
    data["perm_clay"] = 1
if k_model == 'k_model_2':
    data["perm_clay"] = 1
if k_model == 'k_model_3':
    data["perm_clay"] = 1e-3
if k_model == 'k_model_4':
    data["perm_clay"] = 1e-4
if k_model == 'k_model_5':
    data["perm_clay"] = 1e-5
    
num_dim = 5
number_of_facies_Vette = 20
number_of_facies_Troll = int(number_of_facies_Vette/5) #Ensure equal expected geocell heights
number_of_samples = 10000
num_s_points = 21
#s_disc = np.logspace(-3, 0, num_s_points)
if s_disc_mod == 'log':
    s_disc = np.logspace(-6, 0, num_s_points)
    
if s_disc_mod == 'lin':
    s_disc = np.linspace(1e-6, 1, num_s_points)
init_loc_Vette = 700
init_loc_Troll = 1500




K, S, Kw, Knw, Pc = fault_facies(
    number_of_facies_Vette, number_of_samples, data, data_interpr='depth_dep',s_disc=s_disc, k_model=k_model,sect='Vette', init_loc=init_loc_Vette,  two_phase=True
)

#Uncomment when we want to wtite to file for copula fitting
"""
if s_disc_mod == 'log':
    pathName = "/Users/pepe/Documents/FRISK/Python_scripts/smeaheia-test/copula_data/" + k_model + "/"
if s_disc_mod == 'lin':
    pathName = "/Users/pepe/Documents/FRISK/Python_scripts/smeaheia-test/copula_data_lin/" + k_model + "/"
np.savetxt(pathName + "K_col_log_samples.txt", np.log(K), fmt="%s")
np.savetxt(pathName + "Pc_col_log_samples.txt", np.log(Pc[:,0]), fmt="%s") #10
np.savetxt(pathName + "S_col_log_samples.txt", np.log(S[:,0]), fmt="%s") #10
np.savetxt(pathName + "Kw_col_log_samples.txt", np.log(Kw[:,0]), fmt="%s") #10
np.savetxt(pathName + "Knw_col_samples.txt", Knw[:,0], fmt="%s") #10

np.savetxt(pathName + "Pc_log_samples.txt", np.log(Pc), fmt="%s")
np.savetxt(pathName + "S_log_samples.txt", np.log(S), fmt="%s")
np.savetxt(pathName + "Kw_log_samples.txt", np.log(Kw), fmt="%s")
np.savetxt(pathName + "Knw_samples.txt", Knw, fmt="%s")

np.savetxt(pathName + "s_disc.txt", s_disc, fmt="%s")
#"""
