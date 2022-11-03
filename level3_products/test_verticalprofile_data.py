# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 15:51:53 2022

@author: Dean

Tests if the ```wisper_oracles_verticalprofiles_*.nc``` files were succesfully 
created. Plot all vertical profiles for several variables, for each sampling period.
"""



import xarray as xr
import matplotlib.pyplot as plt



## Load data
##-----------------------------------------------------------------------------
prfs2016 = xr.load_dataset("wisper_oracles_verticalprofiles_2016.nc", decode_times=False)
prfs2017 = xr.load_dataset("wisper_oracles_verticalprofiles_2017.nc", decode_times=False)
prfs2018 = xr.load_dataset("wisper_oracles_verticalprofiles_2018.nc", decode_times=False)
##-----------------------------------------------------------------------------



## All profiles for temperature, pressure, water, cwc, total water isotope ratios
##-----------------------------------------------------------------------------
fig, axset = plt.subplots(3, 6, figsize=(14,10))
varkeys_plot = ['T', 'P', 'q', 'cwc', 'dD', 'd18O']

# 2016 sampling period (no cwc data for this year)
for pnum in prfs2016['profile']:
    
    prfsingle = prfs2016.sel(profile=pnum)
    
    for vk, ax in zip(varkeys_plot, axset[0,:]):
        if vk=='cwc': continue
        ax.plot(prfsingle[vk], prfsingle['alt']/1000) # /1000 for units of km.
    

# 2017 sampling period
for pnum in prfs2017['profile']:

    prfsingle = prfs2017.sel(profile=pnum)
    
    for vk, ax in zip(varkeys_plot, axset[1,:]):
        ax.plot(prfsingle[vk], prfsingle['alt']/1000)
    

# 2018 sampling period
for pnum in prfs2018['profile']:

    prfsingle = prfs2018.sel(profile=pnum)
    
    for vk, ax in zip(varkeys_plot, axset[2,:]):
        ax.plot(prfsingle[vk], prfsingle['alt']/1000)
##-----------------------------------------------------------------------------



## Limits, labels, text:        
##-----------------------------------------------------------------------------
axset[2,0].set_xlim(260, 305)


axset[1,0].set_ylabel('altitude (km)', fontsize=14)
plt_xlabels = ['T (K)', 'P (hPa)', 'q (q/kg)', 
               'cwc (g/m3)', 'dD'+u' (\u2030)', 'd18O'+u' (\u2030)']
for ax, lab in zip(axset[2,:], plt_xlabels):
    ax.set_xlabel(lab, fontsize=14)
    

axset[0,3].text(0.5, 0.5, 'No cwc data \nfor 2016 IOP', 
                fontsize=12, va='center', ha='center', 
                transform=axset[0,3].transAxes)


axset[0,-1].text(
    1.05, 0.5, '2016 IOP', ha='left', va='center', 
    rotation=270, fontsize=16, transform=axset[0,-1].transAxes
    )
axset[1,-1].text(
    1.05, 0.5, '2017 IOP', ha='left', va='center', 
    rotation=270, fontsize=16, transform=axset[1,-1].transAxes
    )
axset[2,-1].text(
    1.05, 0.5, '2018 IOP', ha='left', va='center', 
    rotation=270, fontsize=16, transform=axset[2,-1].transAxes
    )
##-----------------------------------------------------------------------------



fig.savefig('verification_verticalprofiles_test.png')