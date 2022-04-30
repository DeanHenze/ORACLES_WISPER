# WISPER ORACLES data products 
## data products, processing code, calibration information

Data products from the Water Isotope System for Precipitation and Entrainment Research 
(WISPER) during the ORACLES campaign. WISPER is designed to provide in-situ aircraft 
meausrements of atmospheric water concentration and its heavy isotope ratios D/H and 
18O/16O for both total water and cloud water concentrations. A detailed review of the 
instrument, measurements, data products, and calibration procedure can be found in 
[Henze et al., 2022](https://doi.org/10.5194/essd-14-1811-2022).

ORACLES (ObseRvations of Aerosols above CLouds and their intEractionS) is a NASA earth 
science field experiment with three Intensive Observation Periods (IOPs) designed to study 
key processes that determine the climate impacts of African biomass burning aerosols. 
More information on the ORACLES field experiment can be found in [Redemann et al., 2021](https://espo.nasa.gov/oracles/content/ORACLES).

The 3 WISPER data products provided are:
* 1Hz time series.
* Individual vertical profiles collocated with latitude, longitude, altitude, temperature, 
and pressure.
* Mean latitude-altitude curtains providing climatology for each IOP.  


The ORACLES experiment included aircraft measurements of the heavy water isotope ratios D/H 
and 18O/16O of atmospheric water vapor and cloud water. Isotope ratio measurements were 
taken with the Water Isotope System for Precipitation and Entrainment Research (WISPER). 

Folders on GitHub
-------------------------

apply_cal+QC: Python scripts for quality control and calibration of all ORACLES WIPSER files. Supporting files. A copy of the WISPER raw data files (not on GitHub).

calibration_modelling: Python scripts and data files for model selection and parameter tuning of all calibration formulas used in 'apply_cal+QC'. For each calibration step, data from all laboratory and field calibrations have been collected into a single file here. Includes plots of the results.
	
Folders not on GitHub
-----------------------------

field_data+notes: Flight notes on WISPER instrument performance. Raw data/log files from the Picarros. Raw data files from the WISPER control box. Raw calibration data from the field.

lab_data+notes: Raw calibration data from the lab in 2017, between the first and second field deployments.
