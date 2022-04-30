# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 17:20:24 2021

@author: DeanH808

Single function in this script is called to add precision columns to a single 
WISPER data file.
"""


import numpy as np # 1.19.2
import matplotlib.pyplot as plt # 3.3.2


def add_precision_cols(data, date, test_plot=False):
    """
    data: pandas.DataFrame:
        Dataframe to add precision columns to.
        
    date: str.
        Date of ORACLES flight, 'yyyymmdd'.
        
    test_plot: bool.
        Set to True to make a test plot of computed precisions.
        
    Returns:
        data with precision columns.
    """
    
    # Calculate Pi2 precisions as functions of log humidity:
    data['std_dD_tot2'] = 6*10**6*np.log(data['h2o_tot2'])**-6.69
    data['std_d18O_tot2'] = 6170*np.log(data['h2o_tot2'])**-4.72
        
        
    # For 2017 and 2018, also make precision columns for Pic1:  
    if date[:4] in ['2017','2018']:
        pic1_h2okeys = ['h2o_tot1','h2o_cld']
        std_D_pic1 = 9*10**6*np.log(data[pic1_h2okeys])**-6.66 
        std_18O_pic1 = 451000*np.log(data[pic1_h2okeys])**-6.56 
    # Assign precision data to columns in 'data':
        data['std_dD_tot1']=std_D_pic1['h2o_tot1']
        data['std_dD_cld']=std_D_pic1['h2o_cld']
        data['std_d18O_tot1']=std_18O_pic1['h2o_tot1']
        data['std_d18O_cld']=std_18O_pic1['h2o_cld']


    # Optional test plot; plot data std's vs humidity:
    if test_plot:
        plt.figure('Precisions test')
        ax1=plt.subplot(1,2,1)   
        ax1.plot(np.log(data['h2o_tot2']), data['std_dD_tot2'], 'ro')    
        ax2=plt.subplot(1,2,2)   
        ax2.plot(np.log(data['h2o_tot2']), data['std_d18O_tot2'], 'ro')  
        if date[:4] in ['2017','2018']:
            ax1.plot(np.log(data['h2o_tot1']), data['std_dD_tot1'], 'bo')    
            ax1.plot(np.log(data['h2o_cld']), data['std_dD_cld'], 'bo') 
            ax2.plot(np.log(data['h2o_tot1']), data['std_d18O_tot1'], 'bo')    
            ax2.plot(np.log(data['h2o_cld']), data['std_d18O_cld'], 'bo')                 
                
        # Plot labels:
    if test_plot:
        plt.figure('Precisions test', fontsize=14)
        ax1.set_title('ORACLES WISPER precision test plot', fontsize=14)
        ax1.set_xlabel('log humidity, log(ppmv)', fontsize=14)
        ax1.set_ylabel('std dD, permil', fontsize=14)
        ax2.set_ylabel('std d18O, permil', fontsize=14)
        
        
    return data