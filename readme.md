WISPER ORACLES data processing and calibration
===============================================

ORACLES (ObseRvations of Aerosols above CLouds and their intEractionS) is a NASA earth science field experiment with three Intensive Observation Periods (IOP) designed to study key processes that determine the climate impacts of African BB aerosols. More information on the ORACLES field experiment: https://espo.nasa.gov/oracles/content/ORACLES.

The ORACLES experiment included aircraft measurements of the heavy water isotope ratios D/H and 18O/16O of atmospheric water vapor and cloud water. Isotope ratio measurements were taken with the Water Isotope System for Precipitation and Entrainment Research (WISPER). 

Folders on GitHub
-------------------------

apply_cal+QC: Python scripts for quality control and calibration of all ORACLES WIPSER files. Supporting files. A copy of the WISPER raw data files (not on GitHub).

calibration_modelling: Python scripts and data files for model selection and parameter tuning of all calibration formulas used in 'apply_cal+QC'. For each calibration step, data from all laboratory and field calibrations have been collected into a single file here. Includes plots of the results.
	
Folders not on GitHub
-----------------------------

field_data+notes: Flight notes on WISPER instrument performance. Raw data/log files from the Picarros. Raw data files from the WISPER control box. Raw calibration data from the field.

lab_data+notes: Raw calibration data from the lab in 2017, between the first and second field deployments.
