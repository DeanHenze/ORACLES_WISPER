WISPER ORACLES data processing and calibration
===============================================

The WISPER processing and calibration procedure is broadly:
	(1) Preprocessing: QC, data-restructuring, removing outliers.
	(2) time synchronization to one of the wing-mounted cloud probes.
	(3) calibration of Picarro-1, which includes correction for humidity-dependence
	    of the isotope ratios and then calibration to an absolute scale.
	(4) calibration of Picarro-2, which for 2016 is analogous to Picarro-1 calibration, 
	    and for 2017, 2018 is a cross-calibration to Picarro-1.
	

Files/Folders on GitHub
-------------------------

outlier_time_intervals

P3_merge_data

WISPER_calibrated_data

WISPER_raw_data

mc_sampler

ORACLES_WISPER_file_header_*.txt

pic1_cal.py

pic2_cal.py

precisions.py

preprocess.py

run_fullcal.py

time_sync.py

uncertainty_estimation.py

	
Files/Folders not on GitHub
-----------------------------
