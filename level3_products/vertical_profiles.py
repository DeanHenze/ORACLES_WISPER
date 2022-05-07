# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 13:03:16 2022

@author: Dean

Create the vertical profiles dataset.


To do:
    - Remember to convert temperature from deg C to deg K.
"""



# Built in:
import sys

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code:
import getdata



def vertical_profiles(year):
    """
    Collects vertical profiles from all flights during a specific year 
    (int, yyyy) into an xarray dataset.
    
    Returns: xarray.dataset.
    """

    ## Table of vertical profiles times for the current year:
    profiletable = pd.read_csv("./vertical_profile_times.csv").copy() # all dates.
    profiletable['year'] = [d//10000 for d in profiletable['date']]
    profiletable = profiletable.loc[profiletable['year']==year]



    ## Keys in the P-3 time series datafiles for variables that will be 
    ## included in the final profiles dataset:
    keys_p3df = ['h2o_gkg', 'dD_permil', 'd18O_permil', 'dxs', 
                 'P_hPa', 'T_K', 'lat', 'lon']
    if year in [2017, 2018]: # Add cloud keys for these sampling periods.
        cldkeys = ['cwc', 'dD_cld', 'd18O_cld', 'dxs_cld']
        keys_p3df = keys_p3df + cldkeys
    
    
    
    ## Objects to collect profile data from each flight into:
        # (1) Dictionary to collect all 2D variables in (dimensions of profile 
        # number, alt). Initialize each entry with an empty pd.DataFrame, 
        # with index = 50 m altitude grid levels:
    alt_grid = np.arange(50, 7050, 50)
    vars2D = {}
    for k in keys_p3df: vars2D[k] = pd.DataFrame({}, index = alt_grid)
        # (2) Lists for 1D variables (dimension profile number):
    flightdate = []
    start_UTCsec = []
    end_UTCsec = []
    ascdesc_flag = [] # Flag for whether P-3 was ascending or descending
    


    ## Load and append data for all profiles in all flights:    
    for date, rows_date in profiletable.groupby(by='date'):
       
        print(date) 
       
        # Load WISPER data with added vars from the merge file:
        data_singledate = get_p3data(str(date))

        for i, row in rows_date.iterrows(): # 1 row per profile.
            
            # Isolate data for a single profile:
            data_prf = data_singledate.loc[
                (data_singledate['Start_UTC'] > row['start_utcsecs']) & 
                (data_singledate['Start_UTC'] < row['end_utcsecs'])
                ]
        
            # Append vars from profile to collection:
            append_gridded_prf(data_prf[keys_p3df + ['height_m']], vars2D)  

            # Append date and start/end times:
            flightdate.append(int(date))
            start_UTCsec.append(data_prf['Start_UTC'].iloc[0])
            end_UTCsec.append(data_prf['Start_UTC'].iloc[-1])
            
            # Append ascending vs. descending flag:
            dz = data_prf['height_m'].iloc[-1] - data_prf['height_m'].iloc[0]
            if dz > 0: ascdesc_flag.append(1)
            if dz < 0: ascdesc_flag.append(0)
    


    ## Create xarray dataset:
        # Empty dataset with altitude and profile-number as the coords:
    n_prf = np.arange(1, len(flightdate)+1, 1) # Profile number coord.
    prfs_xrds = xr.Dataset(coords = dict(profile = n_prf, alt = alt_grid))
        
        # Add 2D vars:
    keys_xr2dvars = ['q', 'dD', 'd18O', 'dxs', # The new 2D variable keys for  
                     'P', 'T', 'lat', 'lon']   # final dataset.
    if year in [2017, 2018]: # Add cloud keys.
        keys_xr2dvars = keys_xr2dvars + ['cwc', 'dD_cld', 'd18O_cld', 'dxs_cld']    
    
    for k1, k2 in zip(keys_xr2dvars, keys_p3df):
        prfs_xrds = prfs_xrds.assign(
            {k1:(["alt", "profile"], vars2D[k2].values)}
            )
        
        # Add 1D vars:
    prfs_xrds = prfs_xrds.assign({'date':(["profile"], flightdate)})    
    prfs_xrds = prfs_xrds.assign({'t_start':(["profile"], start_UTCsec)})    
    prfs_xrds = prfs_xrds.assign({'t_end':(["profile"], end_UTCsec)})   
    prfs_xrds = prfs_xrds.assign({'asc_desc_flag':(["profile"], ascdesc_flag)})   
        
        # Fill missing values with -9999:
    fillval = -9999
    prfs_xrds.fillna(fillval)
              
    
    
    ## Add xarray dataset attributes:
        # global attributes:
    title = ("Water and isotope ratio observations from the ORACLES MISSION,"
             " individual vertical profiles.")
    source = ("Aircraft in-situ data. Time intervals for individual vertical "
              " profiles are isolated and averaged only 50 m vertical levels."
              )
    institution = "Oregon State University"
    references = ("Henze, D., Noone, D., and Toohey, D.: Aircraft"
                  " measurements of water vapor heavy isotope ratios in"
                  " the marine boundary layer and lower troposphere"
                  " during ORACLES, Earth Syst. Sci. Data Discuss."
                  " [preprint], https://doi.org/10.5194/essd-2021-238,"
                  " in review, 2021")
    contact = ("Dean Henze <henzede@oregonstate.edu>"
               " and David Noone <david.noone@auckland.ac.nz>")
    dataset_doi = "10.5281/zenodo.5748368"
    creation_date = "2022-02-01 00:00:00"
    version = "1.0"
    if year==2016:
        campaign = "ORACLES-1"
        period = "2016-08-30 to 2016-09-25"
    if year==2017:
        campaign = "ORACLES-2"
        period = "2017-08-12 to 2017-09-02"
    if year==2018:
        campaign = "ORACLES-3"
        period = "2018-09-27 to 2018-10-23"
    region = "Southeast Atlantic"
    prfs_xrds.attrs = dict(
        title=title, source=source, period=period, region=region, 
        campaign=campaign, insitution=institution, 
        references=references, contact=contact, dataset_doi=dataset_doi, 
        creation_date=creation_date, version=version
                             )
 
        # coordinate attributes:
    prfs_xrds.lat.attrs = dict(
        long_name = "Vertical profile number.") 
    prfs_xrds.alt.attrs = dict(
        long_name = "GPS altitude above reference geoid", 
        units = "meters")

        # variable attributes:
    prfs_xrds.q.attrs = dict(
        long_name = "Total water mixing ratio", 
        units = "g/kg", 
        _FillValue = str(fillval))
    prfs_xrds.dD.attrs = dict(
        long_name = "Total water D/H ratio, expressed in delta-notation", 
        units = "permil", 
        _FillValue = str(fillval))
    prfs_xrds.d18O.attrs = dict(
        long_name = "Total water 18O/16O ratio, expressed in delta-notation", 
        units = "permil", 
        _FillValue = str(fillval))
    prfs_xrds.dxs.attrs = dict(
        long_name = "Total water deuterium excess, expressed in delta-notation", 
        units = "permil", 
        _FillValue = str(fillval))
    prfs_xrds.P.attrs = dict(
        long_name = "Static air pressure", 
        units = "hPa", 
        _FillValue = str(fillval))
    prfs_xrds.T.attrs = dict(
        long_name = "Static air temperature", 
        units = "K", 
        _FillValue = str(fillval))
    prfs_xrds.lat.attrs = dict(
        long_name = "Latitude", 
        units = "degrees north", 
        _FillValue = str(fillval))
    prfs_xrds.lon.attrs = dict(
        long_name = "Longitude", 
        units = "degrees east", 
        _FillValue = str(fillval))
    prfs_xrds.date.attrs = dict(
        long_name = "Sampling UTC date", 
        units = "yyyymmdd", 
        _FillValue = str(fillval))
    prfs_xrds.t_start.attrs = dict(
        long_name = "Start time of profile", 
        units = "seconds since UTC midnight", 
        _FillValue = str(fillval))
    prfs_xrds.t_end.attrs = dict(
        long_name = "End time of profile", 
        units = "seconds since UTC midnight", 
        _FillValue = str(fillval))
    prfs_xrds.asc_desc_flag.attrs = dict(
        long_name = "Flag value for whether the aircraft was ascending vs. "
                    "descending during the profile.", 
        units = "1 = ascending, 0 = descending", 
        _FillValue = str(fillval))
    
    if year in [2017, 2018]:
        prfs_xrds.cwc.attrs = dict(
            long_name = "Cloud water content (liquid + ice)", 
            units = "g/m3", 
            _FillValue = str(fillval)) 
        prfs_xrds.dD_cld.attrs = dict(
            long_name = "Cloud water D/H ratio, expressed in delta-notation", 
            units = "permil", 
            _FillValue = str(fillval))
        prfs_xrds.d18O_cld.attrs = dict(
            long_name = "Cloud water 18O/16O ratio, expressed in delta-notation", 
            units = "permil", 
            _FillValue = str(fillval))
        prfs_xrds.dxs_cld.attrs = dict(
            long_name = "Cloud water deuterium excess, expressed in delta-notation", 
            units = "permil", 
            _FillValue = str(fillval)) 
          
                
    return prfs_xrds
        


#def get_p3data(date):
    """
    Load and return P-3 data, as a pd.DataFrame, for the input flight 
    date. Compute and append dxs. For 2017 and 2018 sampling periods, also 
    compute and append cloud water content.
    """
#    p3data = odl.p3_with_mixsegs(date=date) # Load P-3 data.
#    p3data['dxs'] = p3data['dD_permil'] - 8*p3data['d18O_permil']
#    if date[0:4] in ['2017','2018']:
#        p3data['cwc_gm3'] = cvi_cwc(p3data['h2o_cld_gkg'].values, 
#                                    p3data['T_K'].values, 
#                                    p3data['P_hPa'].values*100, 
#                                    p3data['cvi_enhance'].values)
#        p3data['dxs_cld'] = p3data['dD_cld_permil'] - 8*p3data['d18O_cld_permil']
#    return p3data



def get_p3data(date):
    """
    Load WISPER data with added vars from the merge file. Add dxs for 
    total and cloud water.
    """
    data = getdata.wisperaddvars(date)
    data['dxs'] = data['dD_permil'] - 8*data['d18O_permil']
    if date[0:4] in ['2017','2018']:
        data['dxs_cld'] = data['dD_cld'] - 8*data['d18O_cld']
    return data
    
    

def append_gridded_prf(single_prf_df, all_prfs_dict):
    """
    Merge profile data for each variable in single_prf_df (pandas.DataFrame 
    for a single profile) to their corresponding keys in all_prfs_dict 
    (dictionary of pandas.DataFrame's). Average single profile vars onto 50 m 
    grid levels before merge.
    
    This fxn assumes the df key for height is 'height_m'.
    """
    # Average onto gridded levs:
    hght_gridded = (single_prf_df['height_m']/50).round()*50
    prf_50mlevs = single_prf_df.groupby(hght_gridded).mean()
    prf_50mlevs = prf_50mlevs.round(decimals=2)
    
    # Merge each var:
    for key in all_prfs_dict.keys():
        all_prfs_dict[key] = pd.merge(
            all_prfs_dict[key], prf_50mlevs[key], 
            how='left', left_index=True, right_index=True
            )
        
    return None
                        

    
if __name__ == "__main__":
    """
    Create and save profile datasets for each ORACLES year.
    """
    for year in [2016, 2017, 2018]: # ORACLES sampling years.
        print("=====================================================\n"
              "Creating ORACLES %i vertical profiles data product\n"
              "=====================================================" % year)
        profiles_xrds = vertical_profiles(year) # Profiles as xr.dataset.
        fname = r"./wisper_oracles_verticalprofiles_%i.nc" % year
        profiles_xrds.to_netcdf(fname, format="NETCDF4")
        print("Data product created and saved")
        
        
        
def test_output():
    """
    Rough code tests the saved data products from a call to main.
    """

    year = 2018
    path_prf = r"../output/wisper_oracles_verticalprofiles_%i.nc" % year
    prfs = xr.load_dataset(path_prf, decode_times=False)
    
    prfs = prfs.assign(mean_lat = prfs['lat'].mean(dim='alt'))
    low_lats = prfs.where(prfs['mean_lat']<-10, drop=True)
    mid_lats = prfs.where((prfs['mean_lat']>-10) & (prfs['mean_lat']<-7), drop=True)
    hi_lats = prfs.where(prfs['mean_lat']>-7, drop=True)
    plt.figure()
    varplot = 'dxs'
    for i in low_lats['profile']:
        prfsingle = low_lats.sel(profile=i)
        plt.plot(prfsingle[varplot], prfsingle['alt'], c='b')
    for i in mid_lats['profile']:
        prfsingle = mid_lats.sel(profile=i)
        plt.plot(prfsingle[varplot], prfsingle['alt'], c='r')
    for i in hi_lats['profile']:
        prfsingle = hi_lats.sel(profile=i)
        plt.plot(prfsingle[varplot], prfsingle['alt'], c='k')
        
        
    for i in [6, 36, 52]:
        plt.figure()
        prfsingle = prfs.sel(profile=i)
        ax1 = plt.axes()
        ax1.plot(prfsingle['T'], prfsingle['alt'], 'k')
        ax2 = ax1.twiny()
        ax2.plot(prfsingle['d18O'], prfsingle['alt'], 'b')
        ax3 = ax1.twiny()
        ax3.plot(prfsingle['q'], prfsingle['alt'], 'g') 
        
        
    if year != 2016:
        plt.figure()
        plt.hist(prfs.dD_cld.values.flatten(), histtype='step')
        plt.hist(prfs.dD.values.flatten(), histtype='step')

















