# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:57:20 2022

@author: Dean

Return WISPER data with a set of other variable from the merge file. Returns 
data for a single flight.
"""



# Third party
import pandas as pd
import netCDF4 as nc # 1.3.1

# Local code
import convertq



def data_singledate(date, mergevarkeys_nc, mergevarkeys_return):
    """
    Inputs
    ------
    date: str.
        'yyyymmdd'.
    
    mergevarkeys_nc: list/tuple of str's.
        Keys for variables in the merge file to include.
    
    mergevarkeys_return: list/tuple of str's.
        New keys to assign the merge file vars. Same length as mergevarkeys_nc.
    """
    
    year = date[0:4]
    
    
    # Path and file info to load the datasets:
    relpath_wisper = r"../apply_cal+QC/WISPER_calibrated_data/"  
    wisper_headerline = {'2016':70, '2017':85, '2018':85}[year]
    relpath_merged = r"../apply_cal+QC/P3_merge_data/"
    merged_revnum = {'2016':'R25', '2017':'R18', '2018':'R8'}[year]
    
    
    # WISPER data
    #   1) Get a single column for each variable filled with Pic1 instrument 
    #      data where available and Pic2 otherwise:
    #   2) Convert water concentration from ppmv units to g/kg.
    wisper = pd.read_csv(
        relpath_wisper + "WISPER_P3_%s_R2.ict" % date, 
        header=wisper_headerline
        )

    if year == '2016': # Only Pic2 data available.
        wisper_new = wisper[['Start_UTC','h2o_tot2','dD_tot2','d18O_tot2']]
        #wisper_new['h2o_tot2'] = convertq.ppmv_to_gkg(wisper_new['h2o_tot2'])
        #wisper_new.columns = ['Start_UTC','h2o_gkg','dD_permil','d18O_permil']
 
    elif year in ['2017','2018']: # Pic1 and Pic2 data available.
        wisper_new = wisper[['Start_UTC','h2o_tot1','dD_tot1','d18O_tot1']] # Pic1 values.
        for k in wisper_new.columns[1:]:
            k2 = k[:-1]+'2'
            inan = wisper_new[k].isnull() # Where Pic1 has NAN
            wisper_new.loc[inan, k] = wisper.loc[inan, k2].copy() # Replace with Pic2.
        
        #wisper_new['h2o_tot1'] = convertq.ppmv_to_gkg(wisper_new['h2o_tot1'])
        #wisper_new.columns = ['Start_UTC','h2o_gkg','dD_permil','d18O_permil']
        
    wisper_new.columns = ['Start_UTC','h2o_gkg','dD_permil','d18O_permil']
    wisper_new['h2o_gkg'] = convertq.ppmv_to_gkg(wisper_new['h2o_gkg'])
    
       
    
    # Additional variables
    # Load merged files as nc.Dataset object and place a subset of the 
    # vars in a pandas df:
    merged_nc = nc.Dataset(
        relpath_merged + "mrg1_P3_%s_%s.nc" % tuple([date, merged_revnum])
        )
    merged_pd = pd.DataFrame({})
    for knc, knew in zip(mergevarkeys_nc, mergevarkeys_return):
        merged_pd[knew] = merged_nc.variables[knc][:]
    
    
    return wisper_new.merge(merged_pd, on='Start_UTC', how='inner')
    
    
    
    
    
    
    
    
    
    
    
    