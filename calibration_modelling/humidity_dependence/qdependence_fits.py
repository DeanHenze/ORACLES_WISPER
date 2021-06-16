# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 13:09:33 2020

@author: DeanH808

Calibration of humidity dependence of isotope ratio measurements for Mako. 

Fit the following function to calibration data:
    a*(np.log(50000)-q)**b
where q is humidity in ppmv, and a and b are fit parameters.

There are three sets of calibration data, one for each ORACLES year. Fit the 
above function to each of these sets.
"""



# Third party:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



def get_caldata_dates(caldata, startdate, enddate):
    """
    Return subset of calibration data for dates within and including the start 
    and end dates provided. Returns as a pd.DataFrame.
    
    caldata: pd.DataFrame, calibration data.
    
    startdate, enddate: both str's of the form 'yyyymmdd'.
    """
    # Convert dates to np.datetime64 instances:
    start_dt64 = np.datetime64(startdate)
    end_dt64 = np.datetime64(enddate)
    
    # Caldata rows for the input dates:
    date_mask = (caldata.date >= start_dt64) & (caldata.date <= end_dt64)
    return caldata.loc[date_mask]


   
def qdep_model(logq, a, b):
    """
    Model for dD and d18O humidity dependence curves. Returns correction(s) to 
    be applied for input log humidity value(s).
    
    logq: float or np.array. Natural log of humidity in ppmv.
    
    a, b: floats. Model parameters for either dD or d18O.
    """
    return a*(np.log(50000)-logq)**b



def fit_qdep(caldata, p0_D, bnds_D, p0_18O, bnds_18O):
    """
    Fit the humidity dependence model (qdep_model) to calibration data for 
    dD and d18O, using scipy's curve_fit (non-linear least squares). Returns 
    parameter fits and standard errors for dD and d18O in a single dictionary.
    
    caldata: pd.DataFrame. Contains all calibration data: model input 
        var log(humidity), isotope ratio delta measurements to fit the model 
        with, and standard errors on the measurements.
    
    p0_D, bnds_D: Initial parameter values and paramter-space bounds for the 
        dD fit (see scipy.optimize curve_fit docs).
        
    p0_18O, bnds_18O: Same as above but for d18O fit.
    """
    
    # Get parameter fits and covariance matrix using scipy's curve_fit:        
    pfit_dD, pcov_dD = \
        curve_fit(qdep_model, caldata['log_h2o'], 
                  caldata['dD*_permil'], p0=p0_D,
                  method='trf', sigma=caldata['dD_stderror'], bounds=bnds_D
                  )
    pfit_d18O, pcov_d18O = \
        curve_fit(qdep_model, caldata['log_h2o'], 
                  caldata['d18O*_permil'], p0=p0_18O,
                  method='trf', sigma=caldata['d18O_stderror'], bounds=bnds_18O
                  )

    # Parameter standard errors computed from diagnal of covar matrix:
    sig_pars_dD = np.sqrt(np.diag(pcov_dD))
    sig_pars_d18O = np.sqrt(np.diag(pcov_d18O))
    

    # Return results in a dictionary:
    keys = 'aD','bD','sig_aD','sig_bD','a18O','b18O','sig_a18O','sig_b18O'
    results = np.append(pfit_dD, [sig_pars_dD, pfit_d18O, sig_pars_d18O])
    return dict(zip(keys, results))
             
    
    
def run_calibrations():
    """
    Run all dD and d18O calibrations for Mako and Gulper. Mako has 
    calibrations for all years. Gulper has for only 2016. Plot calibration 
    data and fit curves for Mako. 
    """

    ## Mako calibration data:
    #--------------------------------------------------------------------------
    # Load calibration data:
    caldata_M = pd.read_csv("mako_humidity_dependence_cals.csv")
    
    # Row/column modifications, additions, QC:
        # Date column as datetime64 type:
    caldata_M['date'] = caldata_M['date'].astype('datetime64')
        # Remove bad data for David's 2016 field calibrations:
    caldata_M.loc[caldata_M['date']==np.datetime64("2016-09-09")] = np.nan
    caldata_M.loc[caldata_M['date']==np.datetime64("2016-09-21")] = np.nan
        # Add extra columns:
    caldata_M['h2o_gkg'] = caldata_M['h2o_ppmv']*0.622/1000 # q in units g/kg.
    caldata_M['log_h2o'] = np.log(caldata_M['h2o_ppmv']) # log(q).
        # Add estimates of standard errors:
    caldata_M['dD_stderror'] = 5000/caldata_M['h2o_ppmv']
    caldata_M['d18O_stderror'] = 500/caldata_M['h2o_ppmv']
        # Drop any rows where the isotope measurements are NAN
    caldata_M.dropna(subset=['d18O*_permil','dD*_permil'], inplace=True)
    #--------------------------------------------------------------------------

    
    ## Gulper calibration data:
    #--------------------------------------------------------------------------
    # Load calibration data:
    caldata_G = pd.read_csv("gulper_humidity_dependence_cals.csv")

    # Row/column mods, additions, QC:
    caldata_G['h2o_gkg'] = caldata_G['h2o_ppmv']*0.622/1000 # q in units g/kg.
    caldata_G['log_h2o'] = np.log(caldata_G['h2o_ppmv'])
            # Date column as datetime64 type:
    caldata_G['date'] = caldata_G['date'].astype('datetime64')
        # Remove rows for the 5/29 cal (looks like condensation in the tubes):
    caldata_G.loc[caldata_G['date']==np.datetime64('2017-05-29')] = np.nan
    caldata_G.dropna(subset=['d18O*_permil','dD*_permil'], inplace=True)
    #--------------------------------------------------------------------------
        
    
    ## Fit model for all calibrations (Picarro, year) and store in a pd.df:        
    #--------------------------------------------------------------------------
    # Store parameter fits and standard errors here:
    columns = ['picarro','year','aD','bD','sig_aD','sig_bD','a18O','b18O',
               'sig_a18O','sig_b18O']
    results_df = pd.DataFrame(np.zeros([4,10]), columns=columns)
    results_df['picarro'] = ['Mako','Mako','Mako','Gulper']
    results_df['year'] = ['2016','2017','2018','2016']
    
    
    # Fits for Mako for each ORACLES years:
        # Dateranges for Mako calibration data collection for the 3 years:
    dateranges = [('2017-05-26', '2017-06-01'), ('2017-08-16', '2017-08-27'), 
                  ('2018-10-19', '2018-10-19')]
    
    parameter_colkeys = ['aD','bD','sig_aD','sig_bD','a18O',
                         'b18O','sig_a18O','sig_b18O']
    
    for i in range(3):
        dr = dateranges[i]
        caldata_dates = get_caldata_dates(caldata_M, dr[0], dr[1])   
                
        fitresults = fit_qdep(caldata_dates, # This does the model fit. 
                              p0_D = [-10, 1], bnds_D = ([-30,0], [0,10]), 
                              p0_18O = [-1, 1], bnds_18O = ([-30,0], [0,10]))
        
        results_df.loc[i, parameter_colkeys] = list(fitresults.values())
        
        
    # The one Gulper fit:     
    fitresults_G = fit_qdep(caldata_G,
                            p0_D = [1, 1], bnds_D = ([0,-10], [30,10]), 
                            p0_18O = [0.1, 1], bnds_18O = ([0,-10], [30,10]))
    results_df.loc[3, parameter_colkeys] = list(fitresults_G.values())        
    #--------------------------------------------------------------------------


    ## Plot Mako calibration data and fit curves:
    #--------------------------------------------------------------------------
    plt.figure()
    axlist = []
        # Separate axes for each daterange and for each isotopologue:
    for i in range(6): axlist.append(plt.subplot(3,2,i+1))    


    for i in range(3):
        # Calibration data for all runs in the daterange:
        dr = dateranges[i]
        caldata_dates = get_caldata_dates(caldata_M, dr[0], dr[1])   

        # Scatter plots of data for each calibration run:
        for datetime_obj, data in caldata_dates.groupby('date'):
            
            axlist[2*i].scatter(data['h2o_ppmv'], data['dD*_permil'], s=10,  
                                label=str(datetime_obj.date()))
            axlist[2*i+1].scatter(data['h2o_ppmv'], data['d18O*_permil'], s=10, 
                                  label=str(datetime_obj.date()))
            
    
    # Calibration curves:
    q = np.linspace(500, 18000, 100) # Get model output for these humidities, ppmv.
    logq = np.log(q)

    for i in range(3):
        axlist[2*i].plot(q, qdep_model(logq, *results_df.loc[i, ['aD','bD']]), 'k-')
        axlist[2*i+1].plot(q, qdep_model(logq, *results_df.loc[i, ['a18O','b18O']]), 'k-')
    #--------------------------------------------------------------------------    


    ## Tidy up Mako plot for publication:
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
    axlist[0].set_title(r'Mako $\delta$D(q) calibrations', fontsize=12)
    axlist[1].set_title(r'Mako $\delta^{18}$O(q) calibrations', fontsize=12)
    #--------------------------------------------------------------------------
  

    ## Plot Gulper calibration data and fit curves:
    #--------------------------------------------------------------------------
    plt.figure()
    ax1_G = plt.subplot(1,2,1); ax2_G = plt.subplot(1,2,2)

    # Scatter plot of cal data. Plot separate cal runs in different colors:
    for datetime_obj, data in caldata_G.groupby('date'):
        
        ax1_G.scatter(data['h2o_ppmv'], data['dD*_permil'], s=10,  
                    label=str(datetime_obj.date()))
        ax2_G.scatter(data['h2o_ppmv'], data['d18O*_permil'], s=10, 
                    label=str(datetime_obj.date()))
        
        
    # Plot calibration curves:
    ax1_G.plot(q, qdep_model(logq, *results_df.loc[3, ['aD','bD']]), 'k-')
    ax2_G.plot(q, qdep_model(logq, *results_df.loc[3, ['a18O','b18O']]), 'k-')
    #--------------------------------------------------------------------------
    

    ## Save fit results to csv file:
    results_df.to_csv("qdependence_fit_results.csv", index=False)

        

if __name__=='__main__':
    run_calibrations()
