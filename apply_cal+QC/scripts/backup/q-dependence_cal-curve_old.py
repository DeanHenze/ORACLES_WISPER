# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 13:09:33 2020

@author: DeanH808

Old calibration of humidity dependence of isotope ratio measurements for Mako. 
As of Jan 2021, this calibration is the one used for the WISPER data on the 
NASA ESPO website. 
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.odr import Model, RealData, ODR


def mako_q_dependence_cal():

# Load calibration data:
    relpath_data = ("..\\Calibration_Data\\"
                    "Mako_humidity_dependence_cal_all-dates.xlsx" )
    caldata = pd.read_excel(relpath_data, sheet_name=0, header=3, 
                            usecols=np.arange(0,10))
    
# Clean up rows/columns in the data:
    # Remove data for David's 2016 field calibrations:
    caldata.loc[caldata['date']==np.datetime64("2016-09-09")] = np.nan
    caldata.loc[caldata['date']==np.datetime64("2016-09-21")] = np.nan
    # Add extra columns:
    caldata['h2o_gkg'] = caldata['h2o_ppmv']*0.622/1000 # q in units g/kg.
    caldata['log_h2o'] = np.log(caldata['h2o_ppmv']) # log(q).
    # Add guess estimates of standard errors:
    caldata['SE_dD'] = 5000/caldata['h2o_ppmv']
    caldata['SE_d18O'] = 500/caldata['h2o_ppmv']
    # Drop any rows where the isotope measurements are NAN
    caldata.dropna(subset=['d18O*_permil','dD*_permil'], inplace=True)
 
    
# Fxn to get calibration data for all dates within and including the start and 
# end dates provided:
    def get_caldata_dates(startdate, enddate):
    # Convert dates to np.datetime64 instances:
        start_dt64 = np.datetime64(startdate)
        end_dt64 = np.datetime64(enddate)
    # Caldata rows for the input dates:
        date_mask = (caldata.date >= start_dt64) & (caldata.date <= end_dt64)
        return caldata.loc[date_mask]
    
    
# Function to fit for humidity dependence curves:
    def qdep_fitcurve(beta,x):
        return beta[0] + beta[1]/x**beta[2] + beta[3]/x**beta[4] + beta[5]/x**beta[6]


# Do an ODR fit of the fxn 'fit_curve' to the passed calibration data 
# (log_q and delt); Uses the scipy.odr module:
    def fit_odr(fit_curve, log_q, delt, beta0):
        log_q=np.array(log_q); delt=np.array(delt)
        # Prepare curve to be used with odr object:
        curve = Model(fit_curve)        
        # Prepare data to be used with odr object; remove NAN's:
        i_nan = np.isnan(log_q*delt)
        log_q=log_q[~i_nan]; delt=delt[~i_nan]
        obs = RealData(log_q, delt)
        # Run the ODR model and return fit parameters:
        odr = ODR(obs, curve, beta0=beta0)
        results = odr.run()
        return results.beta
        
    
# Combines fxns 'get_caldata_dates' and 'fit_odr' to fit calibration data for 
# all dates within and including the start and end dates provided. 
# Inputs: startdate and enddate are str 'yyyy-mm-dd'. The beta0 inputs 
# are the initial guesses for the fit curve parameters.
    def odr_fit_to_dates(startdate, enddate, beta0_dD, beta0_d18O): 
        caldata_dates = get_caldata_dates(startdate, enddate)
    # Fit results
        beta_dD = fit_odr(qdep_fitcurve, caldata_dates['log_h2o'], 
                          caldata_dates['dD*_permil'], beta0_dD)
        beta_d18O = fit_odr(qdep_fitcurve, caldata_dates['log_h2o'], 
                            caldata_dates['d18O*_permil'], beta0_d18O)
        return beta_dD, beta_d18O
 
       
# Plot calibration data for each ORACLES year:
    #--------------------------------------------------------------------------
    # The 3 dateranges of calibration data for the 3 ORACLES years:
    dateranges = [('2017-05-26', '2017-06-01'), ('2017-08-16', '2017-08-27'), 
                  ('2018-10-19', '2018-10-19')]

    # Plot data for all calibrations in each daterange:
    fig = plt.figure()
    axlist = []
    for i in range(6): axlist.append(plt.subplot(3,2,i+1))    

    for i in range(3):
    # Get calibration data for all runs in the current daterange:
        dr = dateranges[i]
        caldata_dates = get_caldata_dates(dr[0], dr[1])   

    # Scatter plots of data for each calibration run:
        for datetime_obj, data in caldata_dates.groupby('date'):
            
            axlist[2*i].scatter(data['h2o_ppmv'], data['dD*_permil'], s=10,  
                                label=str(datetime_obj.date()))
            axlist[2*i+1].scatter(data['h2o_ppmv'], data['d18O*_permil'], s=10, 
                              label=str(datetime_obj.date()))
    #--------------------------------------------------------------------------
               

# Fit qdep_fitcurve to data and plot results:        
     #--------------------------------------------------------------------------
   # Fit parameters:
    b0_D=[7,-500000,5,-5000,2, -500, 1] # Initial parameter guesses.
    b0_18O=[0.7,-50000,5,-50,2, -5, 1] # Initial parameter guesses.

    b_dD_original = []; b_d18O_original = [] # Store results here.
    for dr in dateranges:
        fit_results = odr_fit_to_dates(dr[0], dr[1], b0_D, b0_18O)
        b_dD_original.append(fit_results[0])
        b_d18O_original.append(fit_results[1])

    # Plot fit curves:
    q = np.linspace(500, 18000, 1000) # Humidity, units of ppmv
    q_gkg = q*0.622/1000 # Units of g/kg.
    logq = np.log(q) # Calc fit curves at these log(humidity)'s.
   
    for i in range(3):
        axlist[2*i].plot(q, qdep_fitcurve(b_dD_original[i], logq), 'r-')
        axlist[2*i+1].plot(q, qdep_fitcurve(b_d18O_original[i], logq), 'r-')
    #--------------------------------------------------------------------------


# Set plot axes labels, mods, limits:
    #--------------------------------------------------------------------------
    xlim = (350, 23000) # x limits for all axes.
    for i in range(3):
    # Axes limits and scale:
        axlist[2*i].set_xscale("log"); axlist[2*i+1].set_xscale("log")
        axlist[2*i].set_xlim(xlim); axlist[2*i].set_ylim(-60, 3)        
        axlist[2*i+1].set_xlim(xlim); axlist[2*i+1].set_ylim(-16, 1)
            
    # y-labels/ticks and positions:
        axlist[2*i].set_ylabel(r'$\delta$D'+u' (\u2030)', fontsize=12)
        axlist[2*i+1].set_ylabel(r'$\delta^{18}$O'+u' (\u2030)', fontsize=12)
        axlist[2*i+1].yaxis.set_label_position("right")
        axlist[2*i+1].yaxis.tick_right()
        axlist[2*i+1].set_yticks([0,-5,-10])

        axlist[2*i].legend()

    # x-labels:   
    for i in range(4): axlist[i].set_xticks([])
    axlist[4].set_xlabel('q (ppmv)', fontsize=12) 
    axlist[5].set_xlabel('q (ppmv)', fontsize=12)
    

    # Titles:
    axlist[0].set_title('Mako dD(q) calibrations', fontsize=12)
    axlist[1].set_title('Mako d18O(q) calibrations', fontsize=12)
    #--------------------------------------------------------------------------


if __name__ == '__main__':
    mako_q_dependence_cal()