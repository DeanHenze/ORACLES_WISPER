# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 11:50:41 2021

@author: Dean

Contains functions to apply full calibration to Picarro 2 data in a single 
WISPER file after preprocessing and time sync.

For 2016, the calibration is similar to calibration of Picarro 1. For 2017 and 
2018, Picarro 2 is cross-calibrated to Picarro 1 (both humidity and isotope 
ratios).


Main functions to run
---------------------
apply_cal: Main function to run. Calibrates Pic2 data for a single 
    WISPER file.


List of other Functions
-----------------------
apply_cal_2016file, apply_cal_20172018file: Called by 'apply_cal()'.

crosscalibrate_h2o, crosscalibrate_dD_d18O: Cross calibration functions 
    called by 'apply_cal_20172018file()'.
"""


# Third party:
import numpy as np # 1.19.2
import pandas as pd # 1.1.3
import matplotlib.pyplot as plt # 3.3.2

# Pic2 calibration fxns for ORACLES 2016 are taken from the pic1_cal script:
from pic1_cal import q_dep_cal, abs_cal


def apply_cal(wisperdata, date, testplots=False):
    """
    Calibration of Pic2 water concentration and isotope ratios for a single file.

    wisperdata: pandas.DataFrame. 
        WISPER data with Pic2 measurements. Should have all -9999.0 flags 
        replaced with np.nan.
        
    date: str. 
        ORACLES flight data 'yyyymmdd' corresponding to wisperdata.
        
    testplots: bool, default=False.
        Set to True for test plots of the calibration.
    
    Returns:
        wisperdata with the Pic2 measurements modified in place.
    """    
    print("Calibrating Pic1 for P3 flight %s" % date)

    if testplots: wisper_precal = wisperdata.copy()
    
    if date[:4] == '2016': 
        if date in ['20160830','20160831','20160902','20160904']:
            pic = 'Mako'
        else:
            pic = 'Gulper'
        wisperdata = applycal_2016file(wisperdata, pic)
        
    elif date[:4] in ['2017','2018']:
        wisperdata = applycal_20172018file(wisperdata, int(date[:4]))
        
    # Optional test plots:
    if testplots:
        test_plots(wisper_precal, wisperdata, date)
        
    return wisperdata


def applycal_2016file(wisperdata, pic):
    """
    Calibration of Pic2 water concentration and isotope ratios for a single 
    2016 file.
    
    wisperdata: See header for 'fullcal_singlefile()'.
        
    pic: str, either 'Mako' or 'Gulper' for the instrument used.
    
    Returns:
        wisperdata with the Pic2 water concentration modified in place.
    """
    # Model parameters for the humidity-dependence and absolute calibrations
    # are just hard-coded here. Two sets, since 2016 measurements are split 
    # between Mako and Gulper instruments:
    if pic=='Mako':
        # Parameters a, b for qdep_fitcurve:
        aD = -0.365; bD = 3.031; a18O = -0.00581; b18O = 4.961
        # Parameters slope, offset for abs cal of iso ratios:
        mD = 1.056412; kD = -5.957469; m18O = 1.051851; k18O = -1.041851
        # Fudge factor to add to k18O:
        ff=3.5; k18O = k18O + ff
        # Parameters slope, offset for abs cal of H2O:
        mq = 0.8512; kq = 0
        
    if pic=='Gulper':
        aD = 0.035; bD = 4.456; a18O = 0.06707; b18O = 1.889
        # Slopes for abs cal of iso ratios:        
        mD = 1.094037184; # kD = 2.540714192
        m18O = 1.06831472; # k18O = -7.5959267
        # Offsets are derived as outlined in the data paper appendix:
            # Histogram peaks of MBL isotope ratios during routine flights:
        pD_M = -75 # Mako dD peak, +/- 3 permil.
        pD_G = -94 # Gulper dD peak, +/- 3 permil.
        p18O_M = -11.5 # Mako d18O peak, +/- 0.5 permil.
        p18O_G = -16.7 # Gulper d18O peak, +/- 0.5 permil.
            # Offsets derived from peaks and cal slopes:
        kD = pD_M - mD*pD_G
        k18O = p18O_M - m18O*p18O_G
        # Parameters for abs cal of H2O:
        mq = 0.9085; kq = 0
    
    
    # Apply calibrations (Cal fxns are taken from the pic1_cal script):
        # Humidity dependence corrections:
    wisperdata['dD_tot2'] = q_dep_cal(wisperdata['dD_tot2'], 
                                      wisperdata['h2o_tot2'], aD, bD)
    wisperdata['d18O_tot2'] = q_dep_cal(wisperdata['d18O_tot2'], 
                                        wisperdata['h2o_tot2'], a18O, b18O)
        # Isotope ratio absolute calibration:
    wisperdata['dD_tot2'] = abs_cal(wisperdata['dD_tot2'], mD, kD)
    wisperdata['d18O_tot2'] = abs_cal(wisperdata['d18O_tot2'], m18O, k18O)    
        # Humidity absolute calibration:
    wisperdata['h2o_tot2'] = abs_cal(wisperdata['h2o_tot2'], mq, kq)

    return wisperdata


def applycal_20172018file(wisperdata, year):
    """
    Calibration of Pic2 water concentration and isotope ratios for a 2017 or 
    2018 file.

    wisperdata: See header for 'fullcal_singlefile()'.
        
    year: int. ORACLES year, either 2017 or 2018.
    
    Returns:
        wisperdata with the Pic2 measurements modified in place.
    """
    wisperdata = crosscalibrate_dD_d18O(wisperdata, year)
    wisperdata = crosscalibrate_h2o(wisperdata, year)
    return wisperdata


def crosscalibrate_h2o(wisperdata, year):
    """
    Cross-calibrates Pic2 water concentration for ORACLES 2017 and 2018.
    
    wisperdata, year: See header for 'applycal_20172018file()'.

    Returns:
        wisperdata with the Pic2 water concentration modified in place.
    """
    # Load slope for cross-calibration line:
    dir_pars = r"../calibration_modelling/pic2-pic1_cross-cal/"
    q_xcalslopes = pd.read_csv(dir_pars + "h2o_xcal_results.csv")
    slope = q_xcalslopes.loc[q_xcalslopes['year']==year, 'slope'].values
    
    # Return calibrated h2o:
    wisperdata['h2o_tot2'] = slope*wisperdata['h2o_tot2']
    return wisperdata


def crosscalibrate_dD_d18O(wisperdata, year):
    """    
    Cross-calibrates Pic2 dD and d18O for ORACLES 2017 and 2018. Uses the 
    results of the polynomial fit derived elsewhere. 
    
    wisperdata, year: See header for 'applycal_20172018file()' below.
    
    Returns:
        wisperdata with the Pic2 isotope ratios modified in place.
    """
    # Make a separate dataframe of variables for both the dD and d18O 
    # cross-cal models:
    logq = np.log(wisperdata['h2o_tot2'])
    df_xcalvars = pd.DataFrame({'logq': logq,
                                'dD': wisperdata['dD_tot2'],
                                'd18O': wisperdata['d18O_tot2'],
                                'logq*dD': logq*wisperdata['dD_tot2'],
                                'logq*d18O': logq*wisperdata['d18O_tot2'],
                                })
    
    
    # Load cross-calibration model parameters as pandas.DataFrame's and 
    # recast as dictionaries where keys are the predictor_var column:
    dir_pars = r"../calibration_modelling/pic2-pic1_cross-cal/"
    p_d18O_df = pd.read_csv(dir_pars + ("d18O_xcal_results_%i.csv" % year))
    p_dD_df = pd.read_csv(dir_pars + ("dD_xcal_results_%i.csv" % year))
    p_dD_dict = dict(zip(p_dD_df['predictor_var'], p_dD_df['coeff']))
    p_d18O_dict = dict(zip(p_d18O_df['predictor_var'], p_d18O_df['coeff']))


    def iso_crosscal(df_xcalvars, pars_dict):
        """
        Computes and returns the cross-calibrated data. A linear combo of the 
        xcalvars raised to the appropriate power and mulitplied by the 
        corresponding coefficient in 'pars_dict' (dictionary).
        """
        terms = [] # Will store each term (var^pow*coeff) in the linear combo.
    
        for k in pars_dict.keys():
            if k=='const':
                terms.append(pars_dict[k]*np.ones(np.shape(df_xcalvars)[0]))
            else:
                # Get the xcal-var name and power it needs to be raised to:
                varname, power = k.split('^')
                # Compute term:
                terms.append(pars_dict[k]*df_xcalvars[varname]**int(power))
        
        return np.sum(terms, axis=0)
    
    
    # Apply cross-calibration and return:
    wisperdata['dD_tot2'] = iso_crosscal(df_xcalvars, p_dD_dict)
    wisperdata['d18O_tot2'] = iso_crosscal(df_xcalvars, p_d18O_dict)
    return wisperdata

    
def test_plots(precal, postcal, date):
    """
    Test plots. 'precal'', 'postcal' (pandas.DataFrames) are the 
    WISPER data pre/post calibration for the ORACLES flight on 'date' 
    (str, 'yyyymmdd').
    """

    # Time series before and after cross-cal:
        # Light running mean:
    if date[:4] == '2016':
        vars2smooth = ['dD_tot2','d18O_tot2']
    elif date[:4] in ['2017','2018']:
        vars2smooth = ['dD_tot1','dD_tot2','d18O_tot1','d18O_tot2']
    precal[vars2smooth] = precal[vars2smooth].rolling(window=10).mean()    
    postcal[vars2smooth] = postcal[vars2smooth].rolling(window=10).mean()    
        # Plot H2O:
    plt.figure()
    ax1 = plt.axes()
    ax1.plot(precal['Start_UTC'], precal['h2o_tot2'], 'b', label='pic2 pre')
    ax1.plot(postcal['Start_UTC'], postcal['h2o_tot2'], 'r', label='pic2 post')
        # Plot dD, d18O:
    plt.figure()
    ax2 = plt.subplot(2,1,1)
    ax2.plot(precal['Start_UTC'], precal['dD_tot2'], 'b', label='pic2 pre')
    ax2.plot(postcal['Start_UTC'], postcal['dD_tot2'], 'r', label='pic2 post')
    ax3 = plt.subplot(2,1,2)
    ax3.plot(precal['Start_UTC'], precal['d18O_tot2'], 'b', label='pic2 pre')
    ax3.plot(postcal['Start_UTC'], postcal['d18O_tot2'], 'r', label='pic2 post')
        # If a 2017 and 2018 flight, also plot Pic1 vars:
    if date[:4] in ['2017','2018']:
        ax1.plot(precal['Start_UTC'], precal['h2o_tot1'], 'k', label='pic1')
        ax2.plot(precal['Start_UTC'], precal['dD_tot1'], 'k', label='pic1')
        ax3.plot(precal['Start_UTC'], precal['d18O_tot1'], 'k', label='pic1')
       
        
    # Pic1 H2O vs Pic2 H2O before and after cross-cal for 2017 or 2018:
    if date[:4] in ['2017','2018']:
        plt.figure()
        ax4 = plt.axes()
        ax4.plot(precal['h2o_tot2'], precal['h2o_tot1'], 'bo', label='pre cal')
        ax4.plot(postcal['h2o_tot2'], postcal['h2o_tot1'], 'ro', label='post cal')
            # 1-1 line:
        ax4.plot(np.linspace(200,20000,100), np.linspace(200,20000,100), 
                 'k-', label='1-1 line.')

    
    # Figure labels, legends:
    ax1.set_title('Test plot for Pic2 H2O calibration, '
                  'flight on %s' % date, fontsize=14) 
    ax1.set_xlabel('Time (secs since UTC midnight)', fontsize=14)
    ax1.set_ylabel('H2O (ppmv)', fontsize=14)
    ax2.set_title(r'Test plots for Pic2 $\delta$D, $\delta^18$O calibration, '
                  'flight on %s' % date, fontsize=14)
    ax2.set_ylabel(r'$\delta$D (permil)', fontsize=14)
    ax3.set_xlabel('Time (secs since UTC midnight)', fontsize=14)
    ax3.set_ylabel(r'$\delta^18$O (permil)', fontsize=14)
    
    ax1.legend(fontsize=12)
    ax2.legend(fontsize=12)

    if date[:4] in ['2017','2018']:
        ax4.set_title('Test plot for Pic2 H2O calibration, '
                      'flight on %s' % date, fontsize=14)
        ax4.set_xlabel('Pic2 H2O (ppmv)', fontsize=14)
        ax4.set_ylabel('Pic1 H2O (ppmv)', fontsize=14)
        ax4.legend(fontsize=12)
