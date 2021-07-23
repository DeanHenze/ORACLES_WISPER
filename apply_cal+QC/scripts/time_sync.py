# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 13:00:10 2021

@author: Dean
"""


# Built in:
import os

# Third party:
import numpy as np
import pandas as pd
import netCDF4 as nc


path_basicprep_dir = r"../WISPER_data/basic_prep/"
path_timesync_dir = r"../WISPER_data/time_sync/"
path_p3merge = r"../Ancillary/P3_Merge_Data/"


if not os.path.isdir(path_timesync_dir): os.mkdir(path_timesync_dir)


##_____________________________________________________________________________
"""
Returns the WISPER basic-prep data as a pandas df with air pressure added. 
Input date is str 'yyyymmdd'.
"""
##_____________________________________________________________________________
def data_with_pressure(date):
    
    global path_basicprep_dir, path_p3merge
    
    # Load WISPER data as pandas df, NAN-masked:
    fname_wisper = 'WISPER_'+date+'_basic-prep.ict'
    wisper = pd.read_csv(path_basicprep_dir + fname_wisper, header=0)
    wisper.replace(-9999.0, np.nan, inplace=True)

    # Load P3 merge data as nc dataset (has the pressure data):
    year = date[0:4]
    if year=='2016': revnum = 'R25'
    if year=='2017': revnum = 'R18'
    if year=='2018': revnum = 'R8'
    fname_merge = "mrg1_P3_" + date + "_" + revnum + ".nc"
    merge = nc.Dataset(path_p3merge + fname_merge)
    
    # Get time and pressure from the merge data into a pandas df, NAN-masked:
    p = pd.DataFrame({'Static_Pressure':merge.variables['Static_Pressure'][:], 
                      'Start_UTC':merge.variables['Start_UTC'][:]})
    p.replace(-9999.0, np.nan, inplace=True)
    
    # Return wisper merged with pressure data:
    return pd.merge(wisper, p, how='inner', on='Start_UTC', sort=True)
##_____________________________________________________________________________


##_____________________________________________________________________________
"""
Time-shift Pic2 total water vars (humidity and iso ratios) to synchronize 
them with Pic1 total water vars using a pressure dependent function. Then, 
shift both Pic1 and Pic2 vars to time-synchronize them with COMA humdity, 
again with a pressure dependence. 

Next, shift Pic1 cloud water vars by a constant offset of -20.5 secs to sync them
with cloud probes' King-probe lwc. 

Note that 9/27/2018 has an extra 30s constant shift fudge factor added, which 
was determined by inspection.

Input: dates (list of strings) of the form 'yyyymmdd' for flight dates. 
"""
##_____________________________________________________________________________
def wisper_tsync(dates):        

    global path_basicprep_dir, path_timesync_dir
    

## Time shift fxn:    
    """  
    Time-shift variables (columns) in a pandas df by a linear pressure 
    dependent function.
    
    Inputs:
        data: pandas df containing the time and pressure variables as well
              as the variable to shift ('data' can contain other vars).
        t_var, P_var: Strings, the column headers of the time and
            pressure variables.
        shift_vars: list of strings, the column headers of the vars to shift.
        params: 2-element list; first element is the slope of the linear
            function, and second element is the offset.
    """
    #--------------------------------------------------------------------------
    def time_shift(data, P_var, t_var, shift_vars, params):
    # Isolate the time, pressure, and shift variables in a new df:
        data_shift = data[[t_var,P_var]+shift_vars]
    # Mask any unreasonably high pressure values (these should be obviously unphysical):
        data_shift.loc[data_shift[P_var]>1100]=np.nan
    # Apply pressure dependent correction and round to the nearest integer:
        data_shift.loc[:,t_var] =                                                 \
            data_shift[t_var] + params[0]*data_shift[P_var] + params[1]
        #data_shift.loc[:,t_var] = np.round(data_shift[t_var].astype(float), 
        #                                   decimals=0)
        data_shift.loc[:,t_var] = data_shift[t_var].round(decimals=0)
    # average over any duplicate time stamps created from the correction:
        data_shift = data_shift.groupby(by=t_var, as_index=False).mean()
    # interpolate any time gaps created from the correction:
        # First check to make sure there are no gaps greater than 1s:
        dt = np.array(data_shift[t_var].iloc[1:])-np.array(data_shift[t_var].iloc[:-1])
        if len(np.where(dt>1)[0])!=0: print('time gap greater than 1s, fixing...')
        # create reference time var with no gaps:
        t0 = data_shift[t_var].iloc[0]; tf = data_shift[t_var].iloc[-1]
        t_ref = pd.DataFrame({'t_ref':np.arange(t0,tf+1,1)})
        # merge data with reference time to get Nan holes in time series
        # where needed:
        data_shift = data_shift.merge(t_ref, how='outer', left_on=t_var, 
                                      right_on='t_ref',sort=True)   
        # interpolate:
        data_shift.interpolate(method='linear', axis=0, limit=3, inplace=True)
    # verify that everything went well with averaging and interpolation:
        if np.sum( pd.notnull(data_shift[t_var]) ) != tf-t0+1:
            print('Warning: Unsuccessful averaging/interpolation.')
        else:
            print('Fixed!')
    # Merge shifted data with original data (without the original data shift
    # variables) and return:
        data_dropped = data.drop(shift_vars, axis=1)
        return data_dropped.merge(data_shift[[t_var]+shift_vars], on=t_var, how='outer',sort=True)
    #--------------------------------------------------------------------------


## Apply shifts and save as new files:    
    for date in dates:
        print('Time-sync data for %s' % date)

    # Specify shift slope(s) and offset(s) depending on the year.
    #   m2C=slope for Pic2 shifting to COMA, m1C=slope for Pic1&2 shifting to COMA
    #   m21=slope for Pic2 shifting to Pic1,
    #   mcld=slope for Pic1 shifting to cloud probes:
        if date[0:4]=='2016':        
            if date[4:8] in ['0830','0831','0902','0904']: m2C=-0.025; b2C=-2.82
            if date[4:8] in ['0912','0914','0918','0920','0925']: m2C=0.006; b2C=-24.83
            if date[4:8] in ['0910']: m2C=-0.0316; b2C=-92.8
            if date[4:8] in ['0924']: m2C=-0.0234; b2C=-36.57
        if date[0:4]=='2017':
            m21=0.0159; b21=-19.8; m1C=-0.0088; b1C=-11.6; mcld=0; bcld=-20.5 
        if date[0:4]=='2018':
            m1C=-0.0024; b1C=-14.52; mcld=0; bcld=-20.5 
            if date[4:8] in ['0927']: m21=-0.0304; b21=-117.56 - 30 # The 30s is an extra fudge factor.
            if date[4:8] in ['0930','1002','1003','1005','1007','1010']: m21=0.0116; b21=-21.86
            # The extra 10s for b21 is a fudge factor:
            if date[4:8] in ['1012','1015','1017','1019','1021','1023']: m21=0.0057; b21=-187.11 - 10
        
    # Load WISPER data as a pandas df, with the merge file pressure data as
    # an added column:
        data = data_with_pressure(date)    
    
    # Apply shifts:
        if date[0:4]=='2016':
        # Apply time shift for Pic2 relative to COMA:
            data = time_shift(data, 'Static_Pressure', 'Start_UTC', 
                   ['h2o_tot2','dD_tot2','d18O_tot2'], [m2C, b2C])  
        if date[0:4] in ['2017','2018']:
        # Apply time shift for Pic2 relative to Pic1:
            data = time_shift(data, 'Static_Pressure', 'Start_UTC', 
                   ['h2o_tot2','dD_tot2','d18O_tot2'], [m21, b21])    
        # Apply time shift for Pic1 and Pic2 relative to COMA:
            data = time_shift(data, 'Static_Pressure', 'Start_UTC', 
                   ['h2o_tot1','dD_tot1','d18O_tot1','h2o_tot2','dD_tot2','d18O_tot2'], 
                   [m1C, b1C]) 
        # Apply time shift for Pic1 cloud vars relative to cloud probes:
            data = time_shift(data, 'Static_Pressure', 'Start_UTC', 
                   ['h2o_cld','dD_cld','d18O_cld'], [mcld, bcld])
        # Recalculate lwc using new h2o_cld shifted relative to the enhancement factor:
            data.loc[:,'cvi_lwc'] = 0.622*data['h2o_cld']/1000/data['cvi_enhance']
  
    # Save after dropped pressure column and replacing nan values with the 
    # flag value -9999:
        data.drop(['Static_Pressure'], axis=1, inplace=True)
        data.fillna(-9999, inplace=True)    
        fname_tsync = 'WISPER_' + date + '_time-sync.ict'
        data.to_csv(path_timesync_dir + fname_tsync, index=False)
##_____________________________________________________________________________


if __name__=='__main__':
##_____________________________________________________________________________
## ORACLES research flight dates for days where WISPER took good data:
    dates2016_good = ['20160830','20160831','20160902','20160904','20160910',
                      '20160912','20160914','20160918','20160920','20160924',
                      '20160925']
    dates2017_good = ['20170812','20170813','20170815','20170817','20170818',
                      '20170821','20170824','20170826','20170828','20170830',
                      '20170831','20170902']
    dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
                      '20181012','20181015','20181017','20181019','20181021',
                      '20181023']


## Run prep functions for all good days:
    print("=============================\n"
          "Starting time synchronization\n"
          "=============================")
        
    wisper_tsync(dates2016_good)
    wisper_tsync(dates2017_good)
    wisper_tsync(dates2018_good)
    
    print("=============================\n"
          "Time synchronization complete\n"
          "=============================")
##_____________________________________________________________________________
