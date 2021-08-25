# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:46:16 2021

@author: Dean

Contains functions to apply full calibration to Picarro 1 data in a single 
WISPER file after preprocessing and time sync. This applies to ORACLES 2017 
and 2018 sampling periods (no Pic1 for 2016). 

Two steps:
    (1) Humidity dependence correction for dD and d18O, following the formula:
            dD_correction = a*(np.log(50000)-logq)**b, 
        where a and b are fit parameters determined elsewhere.
        
    (2) Absolute calibration of humidity, cloud liquid water content, and 
        isotope ratios. All absolute calibrations follow a line:
            corrected_value = m*uncorrected_value + k.
        
List of Functions
-----------------
apply_cal:
    Function to call to apply calibrations to a single WISPER file.

q_dep_cal, abs_cal:
    Used by 'apply_cal'. No need to call separately.

"""


# Third party:
import numpy as np # 1.19.2
  

## Parameters for steps (1) and (2) equations in the file header:
## -------------------------------
# Step (1), parameters a, b:
    # dD:
aD = [-0.438,-0.331] # Values for ORACLES [2017, 2018]
bD = [2.184,2.856]
    # d18O:
a18O = [-0.01327,-0.00646]
b18O = [3.71,4.577]
      
# Step (2), abs cal of iso ratios (same params for 2017 and 2018):
mD = 1.056412478; kD = -5.957469671 # dD slope and offset.
m18O = 1.051851852; k18O = -1.041851852 # d18O slope and offset.
    # Fudge terms to add to k18O:
ff = [1.25,0]

# Step (2), abs cal of humidity (same params for 2017 and 2018):
mq = 0.8512; kq = 0
## -------------------------------


def q_dep_cal(deltavals, qvals, a, b):
    """
    Apply humidity dependence correction for either dD or d18O. Returns 
    corrected values.
    
    deltavals: float or np.array (length N):
        Either uncalibrated dD or d18O [permil].

    qvals: float or np.array (length N):
         Uncalibrated humidity [ppmv].
         
    a, b: floats:
        Parameters in the humidity-dependence correction formula.
    """
    def qdep_correction(logq, a, b): # Humidity-dependence correction formula.
        return a*(np.log(50000)-logq)**b

    # Apply correction:
    correction = qdep_correction(np.log(qvals), a, b)
    return deltavals - correction
    
  
def abs_cal(x, m, k):
    """
    Formula for absolute calibration of water concentration or iso ratios. 
    Just a line. Returns corrected water conc., dD, d18O, respectively.
    
    x: float or np.array:
        Either humidity [ppmv], dD [permil], or d18O [permil].
        
    m, k: floats:
        Slope and offset.
    """
    return m*x + k

    
def apply_cal(data, date):
    """
    Apply calibrations to a single WISPER file.
    
    data: pandas DataFrame.
        Data for a single P3 flight, -9999.0 flags replaced with np.nan.
        
    date: str.
        P3 flight date, 'yyyymmdd'.
    """
    print("Calibrating Pic1 for P3 flight %s" % date)
    
    if date[:4] == '2017': i = 0
    if date[:4] == '2018': i = 1
    
    # Humidity-dependence corrections:    
        # dD:
    data['dD_tot1'] = q_dep_cal(data['dD_tot1'], data['h2o_tot1'], aD[i], bD[i])
    data['dD_cld'] = q_dep_cal(data['dD_cld'], data['h2o_cld'], aD[i], bD[i])
        # d18O:
    data['d18O_tot1'] = q_dep_cal(data['d18O_tot1'], data['h2o_tot1'], 
                                  a18O[i], b18O[i])
    data['d18O_cld'] = q_dep_cal(data['d18O_cld'], data['h2o_cld'], 
                                 a18O[i], b18O[i])
    
    # Absolute calibration of isotope ratios:
        # dD:
    data[['dD_tot1','dD_cld']] = abs_cal(data[['dD_tot1','dD_cld']], mD, kD)
        # d18O:
    data[['d18O_tot1','d18O_cld']] = abs_cal(data[['d18O_tot1','d18O_cld']], 
                                             m18O, k18O + ff[i]
                                             ) 
    
    # Absoluate calibration for humidity and cloud lwc:
    data[['h2o_tot1','h2o_cld']] = abs_cal(data[['h2o_tot1','h2o_cld']], mq, kq)
        # Recalc cloud lwc from corrected cloud h2o:
    data['cvi_lwc'] = 0.622*data['h2o_cld']/1000/data['cvi_enhance']
    
    return data