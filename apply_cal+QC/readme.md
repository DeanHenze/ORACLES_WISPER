## Time series QC and calibration

### !!Heads up!!
Several data file folders required to generate the calibrated time series data are not in this GitHub repo due to larger storage needs, 
contact deanchenze@gmail.com to obtain them:
* ```./WISPER_raw_data/```
* ```./outlier_time_intervals/```
* ```./P3_merge_data/```

### Calibration procedure
The WISPER processing and calibration procedure is broadly:
* Preprocessing: QC, data-restructuring, removing outliers.
* Time synchronization to one of the wing-mounted cloud probes.
* Add isotope ratio precision columns as a function of humidity.
* Calibration of Picarro-1, which includes correction for humidity-dependence
  of the isotope ratios and then calibration to an absolute scale.
* Calibration of Picarro-2, which for 2016 is analogous to Picarro-1 calibration, 
  and for 2017, 2018 is a cross-calibration to Picarro-1.
