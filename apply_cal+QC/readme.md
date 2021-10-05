WISPER ORACLES data processing and calibration
===============================================

Quick guide
-----------

Run the script 'run_fullcal.py' with Python to get fully calibrated files. The calibrated files will 
be in the 'WISPER_calibrated_data' folder.

Run the 'uncertainty_estimation.py' script to get isotope ratio uncertainties on the WISPER 
measurements after calibration. A function is fit to the uncertainties (as described elsewhere) 
and the fit parameters are saved in './uncertainty_params.csv'.  


Calibration overview and .py files
----------------------------------

The WISPER processing and calibration procedure is broadly:
	(1) Preprocessing: QC, data-restructuring, removing outliers (preprocess.py).
	(2) time synchronization to one of the wing-mounted cloud probes (time_sync.py).
	(3) Add isotope ratio precision columns as a function of humidity (precisions.py).
	(4) calibration of Picarro-1, which includes correction for humidity-dependence
	    of the isotope ratios and then calibration to an absolute scale (pic1_cal.py).
	(5) calibration of Picarro-2, which for 2016 is analogous to Picarro-1 calibration, 
	    and for 2017, 2018 is a cross-calibration to Picarro-1 (pic2_cal.py).
	

Other files/folders on GitHub
-------------------------

mc_sampler.py:
	Monte-Carlo function used for uncertainty estimation.

	
Other files/folders not on GitHub
-----------------------------
P3_merge_data/:
	Includes .nc files of all P3 variables merged for each flight, taken from ESPO. Needed for 
	time_synchronization. 

outlier_time_intervals/:
	Includes .txt files of time intervals manually identified as having bad WISPER data.

WISPER_raw_data/:
	Raw WISPER files. One file per flight.

ORACLES_WISPER_file_header_*.txt:
	Contains header script to add to the beginning of each calibrated WISPER file. Format is 
	compatible with NASA ESPO standards.