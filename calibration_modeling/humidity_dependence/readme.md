This directory contains calibration data and processing code to tune model parameters for 
the humidity-dependence correction of Picarro isotope ratio measurements. There is data and 
processing code for two Picarros, "Mako" and "Gulper".


The python script ```fit_qdependence.py``` will produce parameter fits to the function 
described by Equation B2 of [Henze et al., 2022](https://doi.org/10.5194/essd-14-1811-2022). 
```fit_qdependence.py``` produces fits for both Mako and Gulper, stored in the output 
```qdependence_fit_results.csv```.


Calibration data to tune the model (from several lab tests) can be found in 
```mako_humidity_dependence_cals.csv``` and ```gulper_humidity_dependence_cals.csv```. The 
variables dD*_permil and d18O*_permil in those files are the dD and d18O measurements with 
a constant offset, such   that dD*=0 and d18O*=0 for the measurement taken at the highest 
humidity for that particular calibration run.
