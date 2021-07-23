# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 11:50:41 2021

@author: Dean

Full calibration of WISPER Picarro-2 humidity and isotope ratios. For 2016, 
the calibration is similar to calibration of Picarro-1. For 2017 and 2018, 
Picarro-2 is cross-calibrated to Picarro-1 (both humidity and isotope ratios).
"""


# Built in:
import os

# Third party:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## Data paths:
path_pic1cal_dir = r"../WISPER_data/pic1_cal/"
path_pic2cal_dir = r"../WISPER_data/pic2_cal/"
if not os.path.isdir(path_pic2cal_dir): os.mkdir(path_pic2cal_dir)


## Load cross-cal fit parameters:
relpath_caltable = r'../Calibration_Data/calibration_fits.xlsx'
params_dD = pd.read_excel(relpath_caltable, sheet_name='pic2pic1_xcal_dD')
params_d18O = pd.read_excel(relpath_caltable, sheet_name='pic2pic1_xcal_d18O')
    # Convert to dictionaries and drop the df index:
params_dD = dict(zip(params_dD.columns, params_dD.values.squeeze()))
params_d18O = dict(zip(params_d18O.columns, params_d18O.values.squeeze()))


"""
Takes the WISPER (logq, dD) data from a flight and the cross-cal parameter fits 
(above) and returns the corrections to apply to Pic2 dD. 
"""
def xcal_correct_dD(logq, dD, p):
    return (p['logq^-3']*logq**(-3) + p['logq^-2']*logq**(-2)
            + p['logq^-1']*logq**(-1) + p['dD^1']*dD**(1) 
            + p['dD^2']*dD**(2) + p['dD/logq^1']*(dD/logq)**(1) 
            + p['dD/logq^2']*(dD/logq)**(2) + p['const']
            )


"""
Takes the WISPER (logq, d18O) data from a flight and the cross-cal parameter fits 
(above) and returns the corrections to apply to Pic2 d18O. 
"""
def xcal_correct_d18O(logq, d18O, p):
    return (p['logq^-3']*logq**(-3) + p['logq^-2']*logq**(-2)
            + p['logq^-1']*logq**(-1) + p['d18O^1']*d18O**(1) 
            + p['d18O^2']*d18O**(2) + p['d18O/logq^1']*(d18O/logq)**(1) 
            + p['d18O/logq^2']*(d18O/logq)**(2) + p['const']
            )


"""
Takes the WISPER (logq, dD) data from a flight and the cross-cal parameter fits 
(above) and returns the corrected Pic2 dD. 
"""
def xcal_dD(logq, dD, p):
    return (p['logq^3']*logq**(3) + p['logq^2']*logq**(2)
            + p['logq^1']*logq**(1) + p['dD^1']*dD**(1) 
            + p['dD^2']*dD**(2) + p['dD^3']*dD**(3)
            + p['(dD*logq)^1']*(dD*logq)**(1) 
            + p['(dD*logq)^2']*(dD*logq)**(2)
            + p['(dD*logq)^3']*(dD*logq)**(3)
            + p['const']
            )


"""
Takes the WISPER (logq, d18O) data from a flight and the cross-cal parameter fits 
(above) and returns the corrected Pic2 d18O. 
"""
def xcal_d18O(logq, d18O, p):
    return (p['logq^3']*logq**(3) + p['logq^2']*logq**(2)
            + p['logq^1']*logq**(1) + p['d18O^1']*d18O**(1) 
            + p['d18O^2']*d18O**(2) + + p['d18O^3']*d18O**(3)
            + p['(d18O*logq)^1']*(d18O*logq)**(1) 
            + p['(d18O*logq)^2']*(d18O*logq)**(2)
            + p['(d18O*logq)^3']*(d18O*logq)**(3)
            + p['const']
            )


## Apply calibration for each flight:
dates2017_good = ['20170815','20170817','20170818','20170821','20170824',
                  '20170826','20170828','20170830','20170831','20170902']

for date in dates2017_good:
    # Load WISPER data:
    fname = "WISPER_%s_pic1-cal.ict" % date
    data = pd.read_csv(path_pic1cal_dir + fname)
    data.replace(-9999, np.nan, inplace=True) # Fill missing flag with NAN.
    """
    # Compute dD corrections:
    dD_corrections = xcal_correct_dD(np.log(data['h2o_tot2']), 
                                     data['dD_tot2'].rolling(window=10).mean(), 
                                     params_dD
                                     )
    # Apply corrections:
    dD_corrected = data['dD_tot2'] - dD_corrections
    
    # Compute d18O corrections:
    d18O_corrections = xcal_correct_d18O(np.log(data['h2o_tot2']), 
                                         data['d18O_tot2'].rolling(window=10).mean(), 
                                         params_d18O
                                         )
    # Apply corrections:
    d18O_corrected = data['d18O_tot2'] - d18O_corrections
    """
    
    # Corrected dD and d18O:
    dD_corrected = xcal_dD(np.log(data['h2o_tot2']), 
                           data['dD_tot2'].rolling(window=10).mean(), 
                           params_dD
                           )
    d18O_corrected = xcal_d18O(np.log(data['h2o_tot2']), 
                               data['d18O_tot2'].rolling(window=10).mean(), 
                               params_d18O
                               )
    
    # Test plot:
    vars2smooth = ['dD_tot1','dD_tot2','d18O_tot1','d18O_tot2']
    data[vars2smooth] = data[vars2smooth].rolling(window=10).mean()    
        
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(data['Start_UTC'], data['dD_tot1'], 'k')
    plt.plot(data['Start_UTC'], data['dD_tot2'], 'b')
    plt.plot(data['Start_UTC'], dD_corrected, 'r')
    plt.subplot(2,1,2)
    plt.plot(data['Start_UTC'], data['d18O_tot1'], 'k')
    plt.plot(data['Start_UTC'], data['d18O_tot2'], 'b')
    plt.plot(data['Start_UTC'], d18O_corrected, 'r')

