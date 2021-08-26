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
import matplotlib.pyplot as plt # 3.3.2
  

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

    
def apply_cal(data, date, testplots=False):
    """
    Apply calibrations to a single WISPER file.
    
    data: pandas DataFrame.
        Data for a single P3 flight, -9999.0 flags replaced with np.nan.
        
    testplots: bool, default=False.
        Set to True for test plots of the calibration.
        
    date: str.
        P3 flight date, 'yyyymmdd'.
    """
    print("Calibrating Pic1 for P3 flight %s" % date)
    
    if testplots: data_precal = data.copy()
    
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
    
    # Optional test plots:
    if testplots: test_plots(data_precal, data, date)
    
    return data


def test_plots(precal, postcal, date):
    """
    Test plots. 'precal'', 'postcal' (pandas.DataFrames) are the 
    WISPER data pre/post calibration for the ORACLES flight on 'date' 
    (str, 'yyyymmdd').
    """

    # Time series before and after cross-cal:
        # Light running mean:
    vars2smooth = ['dD_tot1','d18O_tot1','dD_cld','d18O_cld']
    precal[vars2smooth] = precal[vars2smooth].rolling(window=10).mean()    
    postcal[vars2smooth] = postcal[vars2smooth].rolling(window=10).mean()    
        # Plot total water and cloud water:
    plt.figure()
    ax11 = plt.subplot(2,1,1)
    ax11.plot(precal['Start_UTC'], precal['h2o_tot1'], 'b', label='pic2 pre')
    ax11.plot(postcal['Start_UTC'], postcal['h2o_tot1'], 'r', label='pic2 post')
    ax12 = plt.subplot(2,1,2)
    ax12.plot(precal['Start_UTC'], precal['h2o_cld'], 'b', label='pic2 pre')
    ax12.plot(postcal['Start_UTC'], postcal['h2o_cld'], 'r', label='pic2 post')
        # Plot total water dD and cloud water dD:
    plt.figure()
    ax21 = plt.subplot(2,1,1)
    ax21.plot(precal['Start_UTC'], precal['dD_tot1'], 'b', label='pic2 pre')
    ax21.plot(postcal['Start_UTC'], postcal['dD_tot1'], 'r', label='pic2 post')
    ax22 = plt.subplot(2,1,2)
    ax22.plot(precal['Start_UTC'], precal['dD_cld'], 'b', label='pic2 pre')
    ax22.plot(postcal['Start_UTC'], postcal['dD_cld'], 'r', label='pic2 post')
        # Plot total water d18O and cloud water d18O:
    plt.figure()
    ax31 = plt.subplot(2,1,1)
    ax31.plot(precal['Start_UTC'], precal['d18O_tot1'], 'b', label='pic2 pre')
    ax31.plot(postcal['Start_UTC'], postcal['d18O_tot1'], 'r', label='pic2 post')       
    ax32 = plt.subplot(2,1,2)
    ax32.plot(precal['Start_UTC'], precal['d18O_cld'], 'b', label='pic2 pre')
    ax32.plot(postcal['Start_UTC'], postcal['d18O_cld'], 'r', label='pic2 post')

    
    # Figure labels, legends:
    ax11.set_title('Test plot for Pic1 H2O calibration, '
                  'flight on %s' % date, fontsize=14) 
    ax11.set_ylabel(r'$H2O_{total}$ (ppmv)', fontsize=14)
    ax12.set_xlabel('Time (secs since UTC midnight)', fontsize=14)
    ax12.set_ylabel(r'$H2O_{cloud}$ (ppmv)', fontsize=14)
    
    ax21.set_title(r'Test plots for Pic1 $\delta$D calibration, '
                   'flight on %s' % date, fontsize=14)
    ax21.set_ylabel(r'$\delta D_{total}$ (permil)', fontsize=14)
    ax22.set_xlabel('Time (secs since UTC midnight)', fontsize=14)
    ax22.set_ylabel(r'$\delta D_{cloud}$ (permil)', fontsize=14)
        
    ax31.set_title(r'Test plots for Pic1 $\delta^{18} O$ calibration, '
                   'flight on %s' % date, fontsize=14)
    ax31.set_ylabel(r'$\delta^18 O_{total}$ (permil)', fontsize=14)
    ax32.set_xlabel('Time (secs since UTC midnight)', fontsize=14)
    ax32.set_ylabel(r'$\delta^18 O_{cloud}$ (permil)', fontsize=14)
    
    ax11.legend(fontsize=12)
    ax21.legend(fontsize=12)
    ax31.legend(fontsize=12)


















