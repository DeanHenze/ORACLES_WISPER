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


# Built in:
import os

# Third party:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from openpyxl import load_workbook


"""
Save pandas dataframe (data) as a new excel sheet (with name 'newsheetname'), 
to a workbook (at path='path_workbook'). If workbook doesn't exist, create it.
"""
def add_excel_sheet(path_workbook, data, newsheetname):
    
    writer = pd.ExcelWriter(path_workbook, engine='openpyxl')

    if os.path.isfile(path_workbook): # Workbook exists. Add sheets.
        book = load_workbook(path_workbook)
        writer.book = book
        
        # ExcelWriter uses writer.sheets to access the sheet. If you leave it 
        # empty it will not overwrite an existing sheet with the same name.
        writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        
        data.to_excel(writer, newsheetname, index=False)
        writer.save()

    else: # Workbook doesn't exist. Create it and save dataframe.     
        data.to_excel(writer, newsheetname, index=False)
        writer.save()


"""
Plot humidity dependence of isotope ratios for Mako. 
Fit the data and plot the fit curve.

overwrite set to 'y' will re-write the appropriate sheet in the 
calibration_results.xlsx table.
"""
def mako_q_dependence_cal(overwrite='n'):

# Load calibration data:
    relpath_data = ("..\\Calibration_Data\\"
                    "Mako_humidity_dependence_cals.xlsx" )
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
    

# Function to fit for humidity dependence curves. Input log of humidity q in 
# ppmv:    
    def qdep_fitcurve(logq, a, b):
        return a*(np.log(50000)-logq)**b
 
    
# Plot calibration data for each ORACLES year:
    #--------------------------------------------------------------------------
    # The 3 dateranges of calibration data for the 3 ORACLES years:
    dateranges = [('2017-05-26', '2017-06-01'), ('2017-08-16', '2017-08-27'), 
                  ('2018-10-19', '2018-10-19')]

    fig = plt.figure()
    axlist = []
    # Separate axes for each daterange and for each isotopologue:
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
    
    # Will store all parameter fits and standard errors in a pandas df:
    columns = ['year','aD','sig_aD','bD','sig_bD','a18O','sig_a18O','b18O',
               'sig_b18O']
    results = pd.DataFrame(np.zeros([3,9]), columns=columns)
    results['year'] = ['2016','2017','2018']
    
    
    for i in range(3): # Loop through calibration dateranges for each year.
    # Get calibration data for all runs in the current daterange:
        dr = dateranges[i]
        caldata_dates = get_caldata_dates(dr[0], dr[1])   

    # Get fit parameters and parameter standard errors::
        bnds = ([-30,0], [-0,10]) # Bounds on parameter values to try.      
        popt_dD, pcov_dD = \
            curve_fit(qdep_fitcurve, caldata_dates['log_h2o'], 
                      caldata_dates['dD*_permil'], p0=[-10, 1],
                      method='trf', sigma=caldata_dates['SE_dD'], bounds=bnds
                      )
        popt_d18O, pcov_d18O = \
            curve_fit(qdep_fitcurve, caldata_dates['log_h2o'], 
                      caldata_dates['d18O*_permil'], p0=[-1, 1],
                      method='trf', sigma=caldata_dates['SE_d18O'], bounds=bnds
                      )

    # Add parameter fits to the df:
            # stand. errors computed from diagnal elements of covar matrix:
        se_b_dD = np.sqrt(np.diag(pcov_dD))
        se_b_d18O = np.sqrt(np.diag(pcov_d18O))  
        
        results.loc[i, ['aD','bD']] = popt_dD
        results.loc[i, ['sig_aD','sig_bD']] = se_b_dD
        results.loc[i, ['a18O','b18O']] = popt_d18O
        results.loc[i, ['sig_a18O','sig_b18O']] = se_b_d18O


    # Plot calibration curves with the scatter data:
    q = np.linspace(500, 18000, 1000) # Units of ppmv
    logq = np.log(q)

    for i in range(3):
        axlist[2*i].plot(q, qdep_fitcurve(logq, *results.loc[i, ['aD','bD']]), 'k-')
        axlist[2*i+1].plot(q, qdep_fitcurve(logq, *results.loc[i, ['a18O','b18O']]), 'k-')
    #--------------------------------------------------------------------------

    
# Plot axes labels, mods, limits:
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
    
    # Second x-axis for q units of g/kg:
    #cf = (0.622/1000) # Conversion factor
    #ax5twin = axlist[4].twiny(); ax6twin = axlist[5].twiny()
    #ax5twin.set_xlim(xlim)
    #xticks = axlist[4].get_xticks()
    #xticklabels = axlist[4].get_xticklabels()
    #ax5twin.set_xticks(xticks)
    #ax5twin.set_xticklabels(cf*xticks)
    #ax5twin.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
    #ax5twin.spines['bottom'].set_position(('outward', 20))
    
    # Titles:
    axlist[0].set_title(r'Mako $\delta$D(q) calibrations', fontsize=12)
    axlist[1].set_title(r'Mako $\delta^{18}$O(q) calibrations', fontsize=12)
    #--------------------------------------------------------------------------
  
        
# Print fit results and save fits to a xlsx file:
    for col in results.columns:
        if col=='year': continue
        if col in ['a18O','sig_a18O']:
            results.loc[:,col] = results.loc[:,col].round(decimals=5)
        else:
            results.loc[:,col] = results.loc[:,col].round(decimals=3)
    
    print('===========\nFit Results\n===========')
    print(results)
    
    relpath_caltable = r'../Calibration_Data/calibration_fits.xlsx' 

    if overwrite=='y':
        add_excel_sheet(relpath_caltable, results, 'Mako_delta(q)')
        print('Results written to excel file.')
    else:
        print('Results not written to excel file.')
        

if __name__=='__main__':
    overwrite = input('Overwrite calibration_table sheet with these '
                      'results?(y/n):')
    mako_q_dependence_cal(overwrite=overwrite)
