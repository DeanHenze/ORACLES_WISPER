# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 16:42:03 2021

@author: Dean

Gather some basic statistics on the ORACLES isotope ratio measurements.
"""


# Built in:
import sys

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as linsegcmap

# My modules:
if r'../..' not in sys.path: sys.path.insert(0,r'../..')
from my_python_modules import oracles_data_loader as odl
from my_python_modules import math_henze as hmath


"""
Returns all data for a subset of P3 variables for the passed ORACLES year (str). 
Adds dxs.

Inputs:
    year (str): 'yyyymmdd'.
    block_secs (int): An optional number of seconds to group/average the data 
        by before oversampling - this reduces computational req's. 
    morethanx_gkg: float, default = 0.1. Keep only data where humidity is 
        greater than the passed value in g/kg.
"""
def get_data(year, block_secs=None, morethanx_gkg=0.5):
    
    ## Collect data from each flight into a single pandas df:
    vars2get = ['time_secUTC','lat','lon','height_m','h2o_gkg','dD_permil',
                'd18O_permil']
    
    data = pd.DataFrame({}, columns=vars2get) # Will hold all data.
    dates = odl.p3_flightdates(year)
    for date in dates:
        data_temp = odl.oracles_p3_merge_subset(date, with_updated_cal=True)
        data = data.append(data_temp[vars2get], ignore_index=True, sort=False)
        
    data.replace(-9999, np.nan, inplace=True) # Change missing value flag.    


    ## Optional average into blocks:
    if block_secs is not None:
        data['block'] = np.round(data.index.values/block_secs)*block_secs
        data = data.groupby('block').mean()
        

    ## Only keep data where humidity is greater than passed value:
    data = data.loc[data['h2o_gkg']>morethanx_gkg]
        

    ## Add dxs and height in km as new columns:
    data['dxs_permil'] = data['dD_permil'] - 8*data['d18O_permil']
    data['height_km'] = data['height_m']/1000
    
    
    return data


def mixed_layer_stats():
    
    ## WISPER ORACLES data for each year, averaged into 10s blocks:
    data16 = get_data('2016', block_secs=10, morethanx_gkg=0.25)
    data17 = get_data('2017', block_secs=10, morethanx_gkg=0.25)
    data18 = get_data('2018', block_secs=10, morethanx_gkg=0.25)
    
    
    ## Get 95% frequency intervals for isotope ratios:
    """
    Get bounds of inner 95% of passed data array. Return as a 2-tuple.
    """
    def freqintvl_95p(data):
        # Cut into quantiles each containing 2.5% of the data:
        dxs_qcut = pd.qcut(data, 20, retbins=True)
        # Range will be the two outermost quantiles:
        return np.round(dxs_qcut[1][[1,-2]], decimals=1)

    """
    Get and print 95% frequency intervals for all isotope ratio variables. 
    Data should be dict-like with the appropriate keys.
    """
    def print_95p_intvls(data):
        for v in ['dD_permil','d18O_permil','dxs_permil']:
            intvl_95p = freqintvl_95p(data[v])
            print("%s: %s" % tuple([v, str(intvl_95p)]))
    
        # mixed layer intervals for each year
    years = ('2016','2017','2018')
    for year, data in zip(years, (data16,data17,data18)):
        print('\n========================================\n'
              + year + ' 95% frequency intervals'
              '\n========================================')
        
        print('height < 500m'+'\n---------')
        ml = data.loc[data['height_m']<500]
        print_95p_intvls(ml)
        
        print('height > 2500m'+'\n---------')
        lt = data.loc[data['height_m']>2500]
        print_95p_intvls(lt)
        









