# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 13:00:10 2021

@author: Dean

Collection of functions to apply pressure-dependent time shifts to WISPER 
data and create verification plots.

Functions
---------
wisper_tsync: 
    Function to call to apply time shifts to a single WISPER file.

test_plots:
    Plot WISPER data before and after shift.

data_with_pressure, time_shift: 
    Used by 'wisper_tsync'. No need to call separately.
"""

# Third party:
import numpy as np # 1.19.2
import pandas as pd # 1.1.3
import netCDF4 as nc # 1.5.3
import matplotlib.pyplot as plt # 3.3.2


def data_with_pressure(df_wisper, date):
    """
    Combine WISPER data with pressure data from the P3 merge file for 
    a single flight. Return combined data as a pandas df.
    
    Inputs
    ------
    df_wisper: pandas df.
        WISPER data for a single flight.
        
    date: str.
        Flight date 'yyyymmdd'.
    """    
    # Make sure missing value flags in WISPER data are replaced with NANs:
    df_wisper.replace(-9999.0, np.nan, inplace=True)

    # Load P3 merge data as nc dataset (has the pressure data):
    path_p3merge = r"./P3_merge_data/"
    year = date[0:4]
    if year=='2016': revnum = 'R25'
    if year=='2017': revnum = 'R18'
    if year=='2018': revnum = 'R8'
    fname_merge = "mrg1_P3_" + date + "_" + revnum + ".nc"
    merge = nc.Dataset(path_p3merge + fname_merge)
        # Time and pressure as a separate df, NAN-masked:
    p = pd.DataFrame({'Static_Pressure':merge.variables['Static_Pressure'][:], 
                      'Start_UTC':merge.variables['Start_UTC'][:]})
    p.replace(-9999.0, np.nan, inplace=True)
        
    # Return wisper merged with pressure data:
    return pd.merge(df_wisper, p, how='inner', on='Start_UTC', sort=True)


def time_shift(data, P_var, t_var, shift_vars, params):
    """  
    Time-shift variables (columns) in a pandas df by a linear pressure- 
    dependent function. Return shifted data as a pandas df.
    
    Inputs
    ------
    data: pandas df.
        Includes columns for the variables to shift as well as time and 
        pressure columns (other columns will be ignored).
    
    t_var, P_var: str's. 
        Column headers of the time and pressure variables in 'data'.
    
    shift_vars: list of str's.
        The column headers of the vars to shift.
    
    params: 2-element list/tuple.
        Slope (first element) and offset (second) for linear function of 
        pressure.
    """
    # Create a new df with the time, pressure, and shift_vars:
    data_shift = data.loc[:,[t_var,P_var]+shift_vars].copy()
    
    
    # Apply pressure dependent correction and round to the nearest integer:
        # First remove any unphyically high pressures:
    data_shift.loc[data_shift[P_var]>1100, P_var]=np.nan
        # Apply shift, to nearest second:
    data_shift.loc[:,t_var] =                                                 \
        data_shift[t_var] + params[0]*data_shift[P_var] + params[1]
    data_shift.loc[:,t_var] = data_shift[t_var].round(decimals=0)
    
    
    # Average over any duplicate timestamps as a result of the correction:
    data_shift = data_shift.groupby(by=t_var, as_index=False).mean()
    
    
    # Interpolate any time gaps resulting from the correction, if needed:
        # Check for gaps greater than 1s. If none, good to go:
    dt = np.array(data_shift[t_var].iloc[1:])-np.array(data_shift[t_var].iloc[:-1])
    
    if len(np.where(dt>1)[0])!=0:
        # If there are gaps > 1s...
        print('\t time gap greater than 1 s after the shift, fixing...')
        # Create reference time var with no gaps:
        t0 = data_shift[t_var].iloc[0]; tf = data_shift[t_var].iloc[-1]
        t_ref = pd.DataFrame({'t_ref':np.arange(t0,tf+1,1)})
        # Merge data with ref time to get Nan rows in the df where needed:
        data_shift = data_shift.merge(t_ref, how='outer', left_on=t_var, 
                                      right_on='t_ref',sort=True)   
        # Interpolate and verify everything went well:
        data_shift.interpolate(method='linear', axis=0, limit=3, inplace=True)
        if np.sum( pd.notnull(data_shift[t_var]) ) != tf-t0+1:
            print('\t WARNING: Unsuccessful averaging/interpolation.')
        else:
            print('\t Fixed!')
    
    
    # Merge shifted data into original df, replacing the original vars:
    data_dropped = data.drop(shift_vars, axis=1)
    return data_dropped.merge(data_shift[[t_var]+shift_vars], on=t_var, 
                              how='outer',sort=True)


def wisper_tsync(df_wisper, date):        
    """
    Apply pressure dependent time shifts to WISPER data from a single flight. 
    Return shifted data as a pandas df.
    
    If both Pic1 and Pic2 data present, first shift Pic2 total water 
    quantities to Pic1. Next, shift Pic1 and/or Pic2 vars to match COMA 
    humidity. Both shifts are pressure dependent. Finally, shift Pic1 cloud 
    water vars by a constant offset of -20.5 secs to sync them with cloud 
    probes' King-probe lwc. 
    
    The flight on 9/27/2018 has an extra 30s constant shift fudge factor added, 
    which was determined by inspection.
    
    Inputs
    -----
    df_wisper: pandas df.
        WISPER data for a single flight.
        
    date: str.
        Flight date 'yyyymmdd'.
    """
    print("Time synchronization for ORACLES flight date %s" % date)
    
    # Specify shift slope(s) and offset(s) depending on the year.
    #   m2C=slope for Pic2 shifting to COMA, m1C=slope for Pic1&2 shifting to COMA
    #   m21=slope for Pic2 shifting to Pic1,
    #   mcld=slope for Pic1 shifting to cloud probes:
    if date[0:4]=='2016':        
        if date[4:8] in ['0830','0831','0902','0904']: 
            m2C=-0.025; b2C=-2.82
        if date[4:8] in ['0912','0914','0918','0920','0925']: 
            m2C=0.006; b2C=-24.83
        if date[4:8] in ['0910']: 
            m2C=-0.0316; b2C=-92.8
        if date[4:8] in ['0924']: 
            m2C=-0.0234; b2C=-36.57
        
    if date[0:4]=='2017':
        m21=0.0159; b21=-19.8; m1C=-0.0088; b1C=-11.6; mcld=0; bcld=-20.5 
    
    if date[0:4]=='2018':
        m1C=-0.0024; b1C=-14.52; mcld=0; bcld=-20.5 
        if date[4:8] in ['0927']: m21=-0.0304; b21=-117.56 - 30 # The 30s is an extra fudge term.
        if date[4:8] in ['0930','1002','1003','1005','1007','1010']: 
            m21=0.0116; b21=-21.86
        # The extra 10s for b21 is a fudge term:
        if date[4:8] in ['1012','1015','1017','1019','1021','1023']: 
            m21=0.0057; b21=-187.11 - 10
    
    
    # Combine WISPER with pressure data:
    data = data_with_pressure(df_wisper, date)    
    
    
    # Apply shifts:    
    if date[0:4]=='2016':
        # Apply time shift for Pic2 relative to COMA:
        data = time_shift(data, 'Static_Pressure', 'Start_UTC', 
                          ['h2o_tot2','dD_tot2','d18O_tot2'], [m2C, b2C]
                          )  
            
    if date[0:4] in ['2017','2018']:
        # Shift for Pic2 relative to Pic1:
        data = time_shift(data, 'Static_Pressure', 'Start_UTC', 
                          ['h2o_tot2','dD_tot2','d18O_tot2'], [m21, b21]
                          )    
        # Shift for Pic1 and Pic2 relative to COMA:
        data = time_shift(data, 'Static_Pressure', 'Start_UTC', 
                          ['h2o_tot1','dD_tot1','d18O_tot1',
                           'h2o_tot2','dD_tot2','d18O_tot2'], 
                          [m1C, b1C]
                          ) 
        # Shift for Pic1 cloud vars relative to cloud probes:
        data = time_shift(data, 'Static_Pressure', 'Start_UTC', 
                          ['h2o_cld','dD_cld','d18O_cld'], [mcld, bcld]
                          )
        # Recalculate lwc using new h2o_cld shifted relative to the enhancement factor:
        data.loc[:,'cvi_lwc'] = 0.622*data['h2o_cld']/1000/data['cvi_enhance']
  
    
    # Return data after dropping pressure column:
    data.drop(['Static_Pressure'], axis=1, inplace=True)
    return data
    

def test_plot(df1, df2, date):
    """
    Test plot to check the time shift.
    
    df1, df2: pandas df's.
        WISPER data before (df1) and after (df2) the shift. 
        
    date: str.
        Flight date 'yyyymmdd'.
    """
    # 0.2 Hz versions of the timeseries as separate df's:
    df1_0p2Hz = df1.copy(); df2_0p2Hz = df2.copy()
    df1_0p2Hz.replace(-9999.0, np.nan, inplace=True)    
    df2_0p2Hz.replace(-9999.0, np.nan, inplace=True)    
    df1_0p2Hz = df1_0p2Hz.groupby(np.round(df1_0p2Hz['Start_UTC']/5)*5).mean()
    df2_0p2Hz = df2_0p2Hz.groupby(np.round(df2_0p2Hz['Start_UTC']/5)*5).mean()


    # ORACLES 2016 only: Pic2 before and after shift:
    if date[:4] == '2016':
        plt.figure()
        ax1 = plt.subplot(3,1,1)
        ax2 = plt.subplot(3,1,2)
        ax3 = plt.subplot(3,1,3)
        
        varsprefix = ['h2o','dD','d18O']
        for ax, v in zip([ax1,ax2,ax3], varsprefix):
            ax.plot(df1_0p2Hz['Start_UTC'], df1_0p2Hz[v+'_tot2'], 
                    'k-', label='Pic2, before')
            ax.plot(df2_0p2Hz['Start_UTC'], df2_0p2Hz[v+'_tot2'], 
                    'k--', label='Pic2, after')
            ax.legend(fontsize=12)
            
        ax1.set_ylabel("H2O (ppmv)", fontsize=12)
        ax2.set_ylabel("dD (permil)", fontsize=12)
        ax3.set_ylabel("d18O (permil)", fontsize=12)
        ax1.set_title("Pic2 shifted relative to Pic1", fontsize=12)        
    

    # Pic1 and Pic2 total water quantities before and after the total shifts:
    if date[:4] in ['2017','2018']:
        plt.figure()
        ax1 = plt.subplot(3,1,1)
        ax2 = plt.subplot(3,1,2)
        ax3 = plt.subplot(3,1,3)
        
        varsprefix = ['h2o','dD','d18O']
        for ax, v in zip([ax1,ax2,ax3], varsprefix):
            ax.plot(df1_0p2Hz['Start_UTC'], df1_0p2Hz[v+'_tot1'], 
                    'r-', label='Pic1, before')
            ax.plot(df1_0p2Hz['Start_UTC'], df1_0p2Hz[v+'_tot2'], 
                    'b-', label='Pic2, before')
            ax.plot(df2_0p2Hz['Start_UTC'], df2_0p2Hz[v+'_tot1'], 
                    'r--', label='Pic1, after')
            ax.plot(df2_0p2Hz['Start_UTC'], df2_0p2Hz[v+'_tot2'], 
                    'b--', label='Pic2, after')
            ax.legend(fontsize=12)
            
        ax1.set_ylabel("H2O (ppmv)", fontsize=12)
        ax2.set_ylabel("dD (permil)", fontsize=12)
        ax3.set_ylabel("d18O (permil)", fontsize=12)
        ax1.set_title("Pic1 and Pic2 after total shifts", fontsize=12)
    
    
    del df1_0p2Hz, df2_0p2Hz 
