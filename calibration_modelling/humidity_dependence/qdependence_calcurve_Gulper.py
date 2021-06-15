# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 13:09:33 2020

@author: DeanH808

Calibration of humidity dependence of isotope ratio measurements for Gulper. 

Fit the following function to calibration data:
    a*(np.log(50000)-q)**b
where q is humidity in ppmv, and a and b are fit parameters.
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
def gulper_q_dependence_cal(overwrite='n'):

# Load calibration data as pandas df:
    relpath_data = ("..\\Calibration_Data\\"
                    "Gulper_humidity-dependence_cals.xlsx" )
    caldata = pd.read_excel(relpath_data, sheet_name=0, header=0, 
                            usecols=np.arange(0,10))
    
    
# Drop some rows, add some columns:
    # Add columns:
    caldata['h2o_gkg'] = caldata['h2o_ppmv']*0.622/1000 # q in units g/kg.
    caldata['log_h2o'] = np.log(caldata['h2o_ppmv'])
    # Remove rows for the 5/29 cal (looks like condensation in the tubes):
    caldata.loc[caldata['date']==np.datetime64('2017-05-29')] = np.nan
    caldata.dropna(subset=['d18O*_permil','dD*_permil'], inplace=True)
    

# Function to fit for humidity dependence curves. Input log of humidity q in 
# ppmv:    
    def qdep_fitcurve(logq, a, b):
        return a*(np.log(50000)-logq)**b
 
    
# Plot calibration data:
    #--------------------------------------------------------------------------
    fig = plt.figure()
    ax1 = plt.subplot(1,2,1); ax2 = plt.subplot(1,2,2)

    # Plot separately the cal run for each day:
    for datetime_obj, data in caldata.groupby('date'):
        
        ax1.scatter(data['h2o_ppmv'], data['dD*_permil'], s=10,  
                    label=str(datetime_obj.date()))
        ax2.scatter(data['h2o_ppmv'], data['d18O*_permil'], s=10, 
                    label=str(datetime_obj.date()))
    #--------------------------------------------------------------------------
    
    
# Fit qdep_fitcurve to data and plot results:        
    #--------------------------------------------------------------------------
    
    # Will store all parameter fits and standard errors in a 1-row pandas df:
    columns = ['year','aD','sig_aD','bD','sig_bD','a18O','sig_a18O','b18O',
               'sig_b18O']
    results = pd.DataFrame(np.zeros([1,9]), columns=columns)
    results['year'] = ['2016']
    
    
    # Get fit parameters and parameter standard errors::
    bnds = ([0,-10], [30,10]) # Bounds on parameter values to try.      
    popt_dD, pcov_dD = \
        curve_fit(qdep_fitcurve, caldata['log_h2o'], 
                  caldata['dD*_permil'], p0=[1, 1],
                  method='trf', sigma=caldata['dD_stderror'], bounds=bnds
                  )
    popt_d18O, pcov_d18O = \
        curve_fit(qdep_fitcurve, caldata['log_h2o'], 
                  caldata['d18O*_permil'], p0=[0.1, 1],
                  method='trf', sigma=caldata['d18O_stderror'], bounds=bnds
                  )

    # Add parameter fits to the df:
        # stand. errors computed from diagnal elements of covar matrix:
    se_b_dD = np.sqrt(np.diag(pcov_dD))
    se_b_d18O = np.sqrt(np.diag(pcov_d18O))  
    
    results.loc[0, ['aD','bD']] = popt_dD
    results.loc[0, ['sig_aD','sig_bD']] = se_b_dD
    results.loc[0, ['a18O','b18O']] = popt_d18O
    results.loc[0, ['sig_a18O','sig_b18O']] = se_b_d18O


    # Plot calibration curves with the scatter data:
    q = np.linspace(500, 18000, 1000) # Units of ppmv
    logq = np.log(q)

    ax1.plot(q, qdep_fitcurve(logq, *results.loc[0, ['aD','bD']]), 'k-')
    ax2.plot(q, qdep_fitcurve(logq, *results.loc[0, ['a18O','b18O']]), 'k-')
    #--------------------------------------------------------------------------

    
# Plot axes labels, mods, limits:
    #--------------------------------------------------------------------------
    xlim = (350, 23000) # x limits for axes.
    
    # Axes limits and scale:
    ax1.set_xscale("log"); ax2.set_xscale("log")
    ax1.set_xlim(xlim)      
    ax2.set_xlim(xlim)
        
    # y-labels/ticks and positions:
    ax1.set_ylabel(r'$\delta$D'+u' (\u2030)', fontsize=12)
    ax2.set_ylabel(r'$\delta^{18}$O'+u' (\u2030)', fontsize=12)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()

    ax1.legend()

    # x-labels:   
    ax1.set_xlabel('q (ppmv)', fontsize=12) 
    ax2.set_xlabel('q (ppmv)', fontsize=12)
    
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
    ax1.set_title(r'Gulper $\delta$D(q) calibrations', fontsize=12)
    ax2.set_title(r'Gulper $\delta^{18}$O(q) calibrations', fontsize=12)
    #--------------------------------------------------------------------------
  
        
# Print fit results and save fits to xlsx file:
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
        add_excel_sheet(relpath_caltable, results, 'Gulper_delta(q)')
        print('Results written to excel file.')
    else:
        print('Results not written to excel file.')
        

if __name__=='__main__':
    overwrite = input('Overwrite calibration_table sheet with these '
                      'results?(y/n):')
    gulper_q_dependence_cal(overwrite=overwrite)
