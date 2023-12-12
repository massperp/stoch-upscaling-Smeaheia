
#  Stochastic upscaling for UP with OPM and ADSS
Python framework for stochastic upscaling, copula modelling of upscaled flow functions, combined with forward uncertainty propagation (UP) using adaptive stratified sampling (ADSS) through the OPM reservoir simulator.
This package contains scripts needed to reproduce the numerical results and relevant figures in **Copula modeling and uncertainty propagation in field-scale simulation of CO**$_2$ **fault leakage** (https://arxiv.org/abs/2312.05851).



After installing dependencies (see below), to perform numerical simulations, run
```
python3 run_opm_ADSS.py hyperrect from the main folder.
```

To reproduce the figures, run (without .py) from the main folder:
```
python3 -m plot_scripts.plot_script
```
where ```plot_script.py``` is to be replaced by the name of the script.


## The Smeaheia Dataset
The Smeaheia Dataset has been published by Gassnova SF and Equinor ASA
and is available from

   https://co2datashare.org/dataset/smeaheia-dataset

We recommend visiting that site for full information about the data set,
as well as more data than that provided here, which is only a simulation
model.

The data in this directory are owned by Equinor ASA and Gassnova SF

The license of the data is the SMEAHEIA DATASETS LICENSE
which can be found in the 'SMEAHEIA%20DATASET%20LICENSE_Gassnova%20and%20Equinor.pdf' file, as well as at:

   https://co2datashare.org/smeaheia-dataset/static/SMEAHEIA%20DATASET%20LICENSE_Gassnova%20and%20Equinor.pdf



## Installation
The following python packages are needed:
* joblib
* time
* sys
* matplotlib
* numpy
* scipy
* opm
* scipy
* os
* pyvinecopulib
* adaptive_stratification
* opm
* typing
* subprocess
* itertools
* functools
* matplotlib
* mako

