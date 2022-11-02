## Time series QC and calibration

### Calibration procedure
The WISPER processing and calibration procedure is broadly:
* Preprocessing: QC, data-restructuring, removing outliers.
* Time synchronization to one of the wing-mounted cloud probes.
* Add isotope ratio precision columns as a function of humidity.
* Calibration of Picarro-1, which includes correction for humidity-dependence
  of the isotope ratios and then calibration to an absolute scale.
* Calibration of Picarro-2, which for 2016 is analogous to Picarro-1 calibration, 
  and for 2017, 2018 is a cross-calibration to Picarro-1.
  
### !! Additional Folders Required !!
This repo includes all necessary processing scripts. However, several folders containing the data files needed to reproduce the calibrated 
time series are not in this GitHub repo due to larger storage requirements, contact deanchenze@gmail.com to obtain them:
* ```./apply_cal+QC/WISPER_raw_data/```
* ```./apply_cal+QC/outlier_time_intervals/```
* ```./apply_cal+QC/P3_merge_data/```
