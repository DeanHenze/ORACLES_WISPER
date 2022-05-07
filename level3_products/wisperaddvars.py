# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:57:20 2022

@author: Dean

function 'data_singledate'
    Return WISPER data with a set of other variable from the merge file. 
    Returns data for a single flight.
"""



# Third party
import pandas as pd
import netCDF4 as nc # 1.3.1

# Local code
import convertq
import cvi_cwc



def wisperaddvars(date):
    """
    Return WISPER data with a set of other variable from the merge file. 
    Returns data for a single flight. Other variables included are:
        longitude, latitude, altitude
        temperature, pressure
    If a 2017 or 2018 sampling period flight, includes WISPER cwc and isotope 
    ratios. cwc is corrected for enhancement factor and density.
    
    Inputs
    ------
    date: str.
        Flight date, 'yyyymmdd'.
    """
    
    year = date[0:4]
    
    
    # WISPER vars:
    wisper = wisperdata(date)
        

    # Additional variables merged into wisper:
        # Path and filename head info for merge data:
    relpath_merged = r"../apply_cal+QC/P3_merge_data/"
    merged_revnum = {'2016':'R25', '2017':'R18', '2018':'R8'}[year]

        # Load merged files as nc.Dataset object and place a subset of the 
        # vars in a pandas df:
    merged_nc = nc.Dataset(
        relpath_merged + "mrg1_P3_%s_%s.nc" % tuple([date, merged_revnum])
        )
    if year in ['2016','2017']: altitude_key='MSL_GPS_Altitude'
    if year == '2018': altitude_key='GPS_Altitude'
    addvarkeys_nc = [altitude_key, 'Latitude', 'Longitude', 
                     'Static_Air_Temp', 'Static_Pressure'
                     ]
    varkeys_assign = ['alt', 'lat', 'lon', 'T_C', 'P_hPa']
    merged_pd = pd.DataFrame({})
    merged_pd['Start_UTC'] = merged_nc.variables['Start_UTC'][:]
    for knc, knew in zip(addvarkeys_nc, varkeys_assign):
        merged_pd[knew] = merged_nc.variables[knc][:]
        
        # Convert temperature to units of degK:
    merged_pd['T_K'] = merged_pd['T_C'] + 273
    merged_pd.drop(labels='T_C', axis=1, inplace=True)
        
        # Combine with WISPER:
    wisper_addvars = wisper.merge(merged_pd, on='Start_UTC', how='inner')


    # Replace cloud h2o var with enhancement-corrected cloud water 
    # content in units of g/m3:
    wisper_addvars['cwc'] = cvi_cwc.cvi_cwc(
        wisper_addvars['h2o_cld_gkg'].values, 
        wisper_addvars['T_K'].values, 
        wisper_addvars['P_hPa'].values*100, 
        wisper_addvars['cvi_enhance'].values)
    wisper_addvars.drop(labels='h2o_cld_gkg', axis=1, inplace=True)
        
    
    return wisper_addvars



def wisperdata(date):
    """
    Returns wisper data for a single flight. Single columns for each vapor 
    var (q, dD, d18O) where Pic1 measurements are used where available and 
    Pic2 is used otherwise. 
        
    Water vars are converted to g/kg units.
    """
        
    # Path and file headerline info:    
    year = date[0:4]
    relpath_wisper = r"../apply_cal+QC/WISPER_calibrated_data/"
    wisper_headerline = {'2016':70, '2017':85, '2018':85}[year]


    # Vapor vars
    # Get a single column for each vapor variable, filled with Pic1  
    # data where available and Pic2 otherwise:
    wisper = pd.read_csv(
        relpath_wisper + "WISPER_P3_%s_R2.ict" % date, 
        header=wisper_headerline
        )

    if year == '2016': # Only Pic2 data available.
        wisper_new = wisper[['Start_UTC','h2o_tot2','dD_tot2','d18O_tot2']]

    elif year in ['2017','2018']: # Pic1 and Pic2 data available.
        wisper_new = wisper[['Start_UTC','h2o_tot1','dD_tot1','d18O_tot1']] # Pic1 values.
        for k in wisper_new.columns[1:]:
            k2 = k[:-1]+'2'
            inan = wisper_new[k].isnull() # Where Pic1 has NAN
            wisper_new.loc[inan, k] = wisper.loc[inan, k2].copy() # Replace with Pic2.

    wisper_new.columns = ['Start_UTC','h2o_gkg','dD_permil','d18O_permil']


    # Convert water vapor from ppmv units to g/kg.
    wisper_new['h2o_gkg'] = convertq.ppmv_to_gkg(wisper_new['h2o_gkg'])
    
    
    # Add cloud vars:
    if year in ['2017','2018']:
        for cloudkey in ['h2o_cld','dD_cld','d18O_cld','cvi_enhance']:
            wisper_new[cloudkey] = wisper[cloudkey]
            
            
    # Convert cloud water from ppmv units to g/kg.
    wisper_new['h2o_cld_gkg'] = convertq.ppmv_to_gkg(wisper_new['h2o_cld'])
    wisper_new.drop(labels='h2o_cld', axis=1, inplace=True)
    

    return wisper_new




