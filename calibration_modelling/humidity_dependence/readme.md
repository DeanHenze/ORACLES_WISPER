# Humidity dependence correction directory contents:

mako_humidity_dependence_cals.csv
----------------------------------
Data for all Mako isotope ratio humidity dependence calibrations for the ORACLES field 
project. The variables dD*_permil and d18O*_permil are shifted relative to the dD and d18O 
measurements so that dD*=0 and d18O*=0 for the measurement taken at the highest humidity 
for the respective calibration run. 

gulper_humidity_dependence_cals.csv
----------------------------------
Same as mako_humidity_dependence_cals.csv but for Gulper. This file also includes standard deviation columns.

fit_qdependence.py
------------------
Python script to fit model parameters to calibration data. The script fits models for both 
the Mako and Gulper Picarros.

qdependence_fit_results.cvs
---------------------------
Parameter fit results and uncertainties.
