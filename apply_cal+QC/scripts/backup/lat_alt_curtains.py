# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 13:48:57 2021

@author: Dean


Make latitude-altitude curtains of humidity, dD, and dxs for each ORACLES 
year.
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
    morethanx_gkg: float, default = 0.2. Keep only data where humidity is 
        greater than the passed value in g/kg.
"""
def get_data(year, block_secs=None, morethanx_gkg=0.2):
    
    vars2get = ['time_secUTC','lat','lon','height_m','h2o_gkg','dD_permil',
                'd18O_permil','CO_ppbv','T_K','RH','King_LWC_ad','P_hPa']
    
    data = pd.DataFrame({}, columns=vars2get) # Will hold all data.
    dates = odl.p3_flightdates(year)
    for date in dates:
        #data_temp = odl.oracles_p3_merge_subset(date)
        data_temp = odl.oracles_p3_merge_subset(date, with_updated_cal=True)
        data = data.append(data_temp[vars2get], ignore_index=True, sort=False)
        
    data.replace(-9999, np.nan, inplace=True) # Change missing value flag.    

    # Optional average into blocks:
    if block_secs is not None:
        data['block'] = np.round(data.index.values/block_secs)*block_secs
        data = data.groupby('block').mean()

    # Only keep data where humidity is greater than passed value:
    data = data.loc[data['h2o_gkg']>morethanx_gkg]

    # Add dxs and height in km as new columns:
    data['dxs_permil'] = data['dD_permil'] - 8*data['d18O_permil']
    data['height_km'] = data['height_m']/1000
    
    return data


"""
Returns latitude-altitude curtains for H2O, dD, d18O, and CO. Curtains are 
computed by oversampling. Data is passed as a pandas df and oversampled map 
is computed at the passed grid_lat and height (1d arrays)

Inputs:
    data: pandas df
    grid_lat, grid_height: Grid latitudes and heights (1d arrays). 
    h2o_weights (default=False): If set to True, use H2O data as additional 
        weighting factors when computing the oversampled dD and d18O maps.
"""
def get_curtains(data, grid_lat, grid_height, 
                 curtvars = ['h2o_gkg','dD_permil','d18O_permil',
                             'dxs_permil','CO_ppbv'], 
                 h2o_weights=False):

    curtains = {} # Will hold curtains:
    for v in curtvars:
        if (v in ['dD_permil','d18O_permil','dxs']) and h2o_weights:
            curtains[v] = hmath.oversampler(data[v], data['lat'], data['height_km'], 
                                            grid_lat, grid_height, w_add=data['h2o_gkg'])
        else:
            curtains[v] = hmath.oversampler(data[v], data['lat'], data['height_km'], 
                                            grid_lat, grid_height)
        
    return curtains


"""
Get vertical profiles of two sources of variance in WISPER measurements for an 
ORACLES sampling period:
    (1) typical standard deviations of h2o, dD, and dxs measurements - their 
        environmental variability.
    (2) Typical instrument precisions. Since they are fxn's of 
        humidity, humidity profiles for the sampling period are used to get 
        the precision profiles.
    
For (1) get standard deviations in 2D lat-height bins of width x 
height of 1 degree x 300m. Then get the mean standard deviation at each 
height bin (average over the latitudes). 

Results are returned as a dictionary of tuples containing pandas Series.
"""
def var_profiles(year): # str, 'yyyymmdd'
    
    ## Get data averaged by 10s blocks and remove data above 6km:
    data = get_data(year, block_secs=10)
    data.loc[data['height_m']>6000] = np.nan
    
    
    ## Data columns for latitude and height bin:
    data['lat_bin'] = np.round(data['lat'])
    data['hght_bin'] = np.round(data['height_m']/300)*300
    
    
    ## Standard deviation vertical profiles:
        # Get variances in each bin using pivot tables:
    var_h2o = pd.pivot_table(data, values='h2o_gkg', index=['hght_bin'], 
                             columns='lat_bin', aggfunc=np.var)
    var_dD = pd.pivot_table(data, values='dD_permil', index=['hght_bin'], 
                            columns='lat_bin', aggfunc=np.var)
    var_dxs = pd.pivot_table(data, values='dxs_permil', index=['hght_bin'], 
                             columns='lat_bin', aggfunc=np.var)
        # standard deviations taken as the square-root of the latitude-averaged 
        # variances:
    std_h2o_prf = (var_h2o.mean(axis=1, skipna=True))**0.5
    std_dD_prf = (var_dD.mean(axis=1, skipna=True))**0.5
    std_d18O_prf = (var_dxs.mean(axis=1, skipna=True))**0.5        
    
    
    ## Get representative vertical profiles of isotope ratio instrument 
    ## precisions, based off the humidity profiles of each latitude bin:
    """
        # First need the mean humidity profile:
    mean_h2o_prf = data[['hght_bin','h2o_gkg']].groupby('hght_bin').mean()
        # dD, d18O precisions are fxns of humidity in units of ppmv:
    logh2o = np.log(mean_h2o_prf['h2o_gkg']*1000/0.622)
    prec_D = (6*10**6)*(logh2o**-6.69)
    prec_18O = 15425*(logh2o**-4.72)
        # dxs precision using the partial derivative method:
    prec_dxs = (prec_D**2 + (8*prec_18O)**2)**0.5
    """ 
        # First get mean humidity profile in each latitude bin from pivot tabs: 
    mean_h2o_prfs = pd.pivot_table(data, values='h2o_gkg', index=['hght_bin'], 
                                   columns='lat_bin', aggfunc=np.mean)
        # Compute precision profiles for each latitude bin:
    precs_dD_prfs = pd.DataFrame({}, index=mean_h2o_prfs.index, 
                                 columns = mean_h2o_prfs.columns)
    precs_dxs_prfs = pd.DataFrame({}, index=mean_h2o_prfs.index, 
                                 columns = mean_h2o_prfs.columns)
    for col in mean_h2o_prfs.columns.values:
        # dD, d18O precisions are fxns of humidity in units of ppmv:
        logh2o = np.log(mean_h2o_prfs[col]*1000/0.622)
        precs_dD_prfs[col] = (6*10**6)*(logh2o**-6.69)
        prec_18O = 15425*(logh2o**-4.72)
        # dxs precision using the partial derivative method:
        precs_dxs_prfs[col] = (precs_dD_prfs[col]**2 + (8*prec_18O)**2)**0.5
        
        # Average precision profiles:
    prec_dD_meanprf = precs_dD_prfs.mean(axis=1)
    prec_dxs_meanprf = precs_dxs_prfs.mean(axis=1)

    
    return {'h2o':(std_h2o_prf,), 'dD':(std_dD_prf, prec_dD_meanprf), 
            'dxs':(std_d18O_prf, prec_dxs_meanprf)} 
    

"""
Returns contour levels and colormappings to use for the H2O, dD, and d18O 
curtain plots.
"""
def get_clevels_cmapping():
    
    ## H2O contour levels
    levs_h2o = np.arange(0,17,1.)
    
    ## Iso ratio contour levels:
    dD_min = -350; dD_max = -64 # dD_min=-260
    d18O_min = -45; d18O_max = -10.5 # d18O_min=-32
        
    a_dD = 1 # Desired spacing in permil of between the first 2 contours.
    n = np.log2((dD_max-dD_min)/a_dD) # Power of 2 that spans the range (dD_max-dD_min)/a_dD.
    levspacing_dD = a_dD*(2**np.linspace(1,n,15))
    levs_dD = np.flip(np.append([dD_max], dD_max-levspacing_dD))
    levs_dD = np.round(levs_dD)
    
    a_d18O = 0.1
    n = np.log2((d18O_max-d18O_min)/a_d18O)
    levspacing_d18O = a_d18O*(2**np.linspace(1,n,15))
    levs_d18O = np.flip(np.append([d18O_max], d18O_max-levspacing_d18O))
    levs_d18O = np.round(levs_d18O, decimals=1)
    
    levs_dxs = np.arange(-10,30,3)
    
    #levs_stddD = [5,10,30,80,100]
    
    ## Color mapping:
        # H2O is straightforward:
    cmap_name = 'gist_ncar' 
    cmap_h2o = plt.get_cmap(cmap_name)
  
        # Iso ratios are less straightforward:
    """
    Map uneven level spacings to even colormap spacings.
    levs: contour levels (array) in increasing order.
    cmap (str): matplotlib colormap to use.
    """
    def cmapping_to_even(levs, cmap_name):
        cmap = plt.get_cmap(cmap_name)
        normedlevs = (levs-levs[0])/(levs[-1]-levs[0]) # Norm to 0-1 scale.
        colors = cmap(np.linspace(0,1,len(levs))) # Cmap colors at even spacings.
        # Map:
        return linsegcmap.from_list(cmap_name, list(zip(normedlevs,colors)))
    
    cmap_dD = cmapping_to_even(levs_dD, cmap_name)
    cmap_d18O = cmapping_to_even(levs_d18O, cmap_name)
        
    
    return dict(h2o_gkg=(levs_h2o,cmap_h2o), dD_permil=(levs_dD, cmap_dD), 
                d18O_permil=(levs_d18O, cmap_d18O), 
                dxs_permil=(levs_dxs, cmap_h2o))


"""
Get MERRA monthly mean planetary boundary layer top interpolated along the 
ORACLES routine flight track lat/lon for the input year (str, 'yyyymmdd'). 
"""
def get_merra_pbltop(year):

    # Load MERRA mean data for the ORACLES year,month:
    if year=='2016': month='09'
    if year=='2017': month='08'
    if year=='2018': month='10'
    merra = odl.merra_monthly_met(month, year)
    
    if year=='2016':        
        # Routine flight track lat/lon was a line between these two points: 
        lon1=-1;lat1=-8.8; lon2=13; lat2=-21.8
        
        slope = (lat2-lat1)/(lon2-lon1)
        offset = lat1 - slope*lon1
        
        # Interpolate to this rout using xarray built-in functions:
        lons_rout = xr.DataArray(np.arange(lon1,lon2,0.5), dims='dummy')
        lats_rout = xr.DataArray(slope*lons_rout + offset, dims='dummy')
        pbl_rout = merra['PBLTOP'].interp(lon=lons_rout, lat=lats_rout)
        
    if year in ['2017','2018']:
        # Routine flight track was a north-south trajectory between latitudes 
        # 0 and 15S, and at a constant longitude of 5E:
        lats_rout = xr.DataArray(np.arange(-12,0,0.5), dims='dummy')
        lons_rout= xr.DataArray(np.ones(len(lats_rout))*5, dims='dummy')
        pbl_rout = merra['PBLTOP'].interp(lon=lons_rout, lat=lats_rout)
        
    # MERRA PBL height is in units of Pa, do a rough conversion to km, 
    # assuming surface pressure of 1013hPa and 100hPa/km:
    return (1013-pbl_rout/100)*(1/100)


"""
Curtain plots for the input year. 

Inputs:
    year (str): 'yyyymmdd' ORACLES year.
    curtains: the output from my 'get_curtains' fxn (above) for when it is fed 
        data for the desired year.  
    grid_lat and grid_hght (1d numpy arrays): the latitude, altitude grid 
        coordinates for the curtains.
"""
def curtain_plots(year, curtains, grid_lat, grid_hght):
    
    fig = plt.figure(figsize=(6.5,7.5))
    ax1 = fig.add_axes([0.25,0.74,0.45,0.25]) # H2O
    ax2 = fig.add_axes([0.25,0.44,0.45,0.25]) # dD
    ax3 = fig.add_axes([0.25,0.14,0.45,0.25]) # d18O
    
    cbax_1 = fig.add_axes([0.05,0.74,0.035,0.25]) # H2O colorbar
    cbax_2 = fig.add_axes([0.05,0.44,0.035,0.25]) # dD colorbar
    cbax_3 = fig.add_axes([0.05,0.14,0.035,0.25]) # d18O colorbar
    
    stdax1 = fig.add_axes([0.725,0.74,0.2,0.25]) # H2O std profile
    stdax2 = fig.add_axes([0.725,0.44,0.2,0.25]) # dD std profile
    stdax3 = fig.add_axes([0.725,0.14,0.2,0.25]) # d18O std profile
    

    # Contour levels:
    clevs_maps = get_clevels_cmapping()


    # Contour maps:
    c1 = ax1.contourf(grid_lat, grid_hght, curtains['h2o_gkg'], 
                      levels=clevs_maps['h2o_gkg'][0], 
                      cmap=clevs_maps['h2o_gkg'][1])
    c2 = ax2.contourf(grid_lat, grid_hght, curtains['dD_permil'], 
                      levels=clevs_maps['dD_permil'][0], 
                      cmap=clevs_maps['dD_permil'][1])
    #c3 = ax3.contourf(grid_lat, grid_hght, curtains['d18O_permil'], 
    #                  levels=clevs_maps['d18O_permil'][0], 
    #                 cmap=clevs_maps['d18O_permil'][1])
    c3 = ax3.contourf(grid_lat, grid_hght, curtains['dxs_permil'], 
                  levels=clevs_maps['dxs_permil'][0], 
                  cmap=clevs_maps['dxs_permil'][1])


    # Get MERRA month mean PBL as a fxn of latitude, plot on each map: 
    hpbl = get_merra_pbltop(year)
    for ax in [ax1,ax2,ax3]: ax.plot(hpbl['lat'].values, hpbl.values, 'k--')
    
    
    # CO=200 ppbv contour on each map:
    for ax in [ax1,ax2,ax3]:
        ax.contour(grid_lat, grid_hght, curtains['CO_ppbv'], 
                   levels=[180], colors='k', linewidths=4)
        
        
    # Mean vertical profiles of (1) h2o, dD, dxs standard deviations and for the 
    # curtains, (2) Instrument precision at a given height's mean humidity. 
    # Plotted to the right of their respective curtains:
    var_prfs = var_profiles(year)
    for key, ax in zip(var_prfs.keys(), (stdax1,stdax2,stdax3)):
        ax.plot(var_prfs[key][0], var_prfs[key][0].index, 'k-')
        if key != 'h2o': 
            ax.plot(var_prfs[key][1], var_prfs[key][1].index, 'k--')


    # Colorbars:
    fig.colorbar(c1, cax=cbax_1, orientation='vertical')
    fig.colorbar(c2, cax=cbax_2, orientation='vertical')
    fig.colorbar(c3, cax=cbax_3, orientation='vertical')
    
    
    # Plot labels, limits, scales:        
    ax3.set_xlabel('Latitude [deg]', fontsize=14)
    for ax in [ax1,ax2,ax3]: ax.set_ylabel('Altitude [km]', fontsize=14)
    
    #plt_sublabels = ['(a)','(b)','(c)']
    #for ax, l in zip([ax1,ax2,ax3], plt_sublabels):
    #    ax.text(-21.5, 5.5, l, fontsize=16, 
    #            ha='left', va='top', backgroundcolor='w', zorder=10)
        
    if year=='2016': xlims = (-22,-10.3)
    if year in ['2017','2018']: xlims = (-12,0)
    for ax in [ax1,ax2,ax3]: ax.set_xlim(*xlims)
    
    stdax2.set_xscale('log'); stdax3.set_xscale('log')
    for ax in [stdax1,stdax2,stdax3]: 
        ax.set_ylim(0,6000); ax.set_yticks([])
    
    cax_ylabels = ('q [g/kg]', r'$\delta$D [permil]', 'dxs [permil]')
    for ax, lab in zip((cbax_1,cbax_2,cbax_3), cax_ylabels):
        ax.set_ylabel(lab, fontsize=12, labelpad=0, rotation=90)
        ax.yaxis.set_label_position("left")
 
    
"""
Vertical profiles for all years, taken from the curtain data. Inputs 
'curtains16', 'curtains17', 'curtains18' are the output from my 'get_curtains' 
fxn (above) when it is fed data from the respective year. 'grid_lat' and 
'grid_hght' (1d numpy arrays) are the latitude, altitude grid coordinates for 
the curtains.
"""
def vertical_profiles(curtains16, curtains17, curtains18, grid_lat16, 
                      grid_lat1718, grid_hght):

    # Initialize figure and profile axes:
        #----------------------------------------------------------------------
    fig = plt.figure(figsize=(6.5,8))
        # Horizontal locations on fig for the left side of each axes in a row:
    hl_prfax = [0.1,0.4,0.55,0.7,0.85]
        # Vertical locs of bottoms:
    vl_prfax = [0.74,0.44,0.14]
        # Axes heights:
    h_prfax = 0.25
        # Axes widths:
    w_prfax = [0.29,0.14,0.14,0.14,0.14]
    
    axset16 = [fig.add_axes([hl_prfax[i], vl_prfax[0], w_prfax[i], h_prfax])
                   for i in range(len(hl_prfax))] # Row of axes for 2016.
    axset17 = [fig.add_axes([hl_prfax[i], vl_prfax[1], w_prfax[i], h_prfax])
                   for i in range(len(hl_prfax))] # Row of axes for 2017.
    axset18 = [fig.add_axes([hl_prfax[i], vl_prfax[2], w_prfax[i], h_prfax])
                   for i in range(len(hl_prfax))] # Row of axes for 2018.
        #----------------------------------------------------------------------
    
    
    # 2016 profiles at set latitudes:
    #--------------------------------------------------------------------------
    prflats16 = np.arange(-18,-8,2) # Latitudes to get profiles at.
    
    for lat in prflats16:
        # index of grid_lat value closest to lat.
        i = np.nanargmin(abs(grid_lat16-lat)) 
        
        axset16[0].plot(curtains16['h2o_gkg'][:,i], grid_hght, 
                 label=str(int(grid_lat16[i])))
        axset16[1].plot(curtains16['dD_permil'][:,i], grid_hght)
        axset16[2].plot(curtains16['dD_permil'][:,i], grid_hght)
        axset16[3].plot(curtains16['d18O_permil'][:,i], grid_hght)    
        axset16[4].plot(curtains16['d18O_permil'][:,i], grid_hght)    
    #--------------------------------------------------------------------------
        
        
    # 2017,2018 profiles at set latitudes:
    #--------------------------------------------------------------------------
    prflats1718 = np.arange(-12,-2,2) # Latitudes to get profiles at.
    
    for lat in prflats1718:
        # index of grid_lat value closest to lat.
        i = np.nanargmin(abs(grid_lat1718-lat)) 
        
        axset17[0].plot(curtains17['h2o_gkg'][:,i], grid_hght, 
                 label=str(int(grid_lat1718[i])))
        axset17[1].plot(curtains17['dD_permil'][:,i], grid_hght)
        axset17[2].plot(curtains17['dD_permil'][:,i], grid_hght)
        axset17[3].plot(curtains17['d18O_permil'][:,i], grid_hght)    
        axset17[4].plot(curtains17['d18O_permil'][:,i], grid_hght)   
        
        axset18[0].plot(curtains18['h2o_gkg'][:,i], grid_hght, 
                 label=str(int(grid_lat1718[i])))
        axset18[1].plot(curtains18['dD_permil'][:,i], grid_hght)
        axset18[2].plot(curtains18['dD_permil'][:,i], grid_hght)
        axset18[3].plot(curtains18['d18O_permil'][:,i], grid_hght)    
        axset18[4].plot(curtains18['d18O_permil'][:,i], grid_hght)   
    #--------------------------------------------------------------------------
    
    # Profile axes labels, ticks, legend:
    #--------------------------------------------------------------------------
        # Set x-limits:
            # 2016:
    axset16[0].set_xlim(0,16) # H2O axes
    axset16[1].set_xlim(-300,-250); axset16[2].set_xlim(-180,-50) # dD axes
    axset16[3].set_xlim(-35,-30); axset16[4].set_xlim(-25,-8) # d18O axes
            # 2017 and 2018:
    for axset in [axset17,axset18]: 
        axset[0].set_xlim(0,16) # H2O axes
        axset[1].set_xlim(-440,-275); axset[2].set_xlim(-180,-50) # dD axes
        axset[3].set_xlim(-47,-30); axset[4].set_xlim(-25,-8) # d18O axes
        
        # Set x-ticks for dD and d18O axes:
            # 2016:
    axset16[1].set_xticks([-260,-290]); axset16[2].set_xticks([-150,-75]) # dD
    axset16[3].set_xticks([-34,-32]); axset16[4].set_xticks([-20,-10]) # d18O
            # 2017 and 2018:
    for axset in [axset17,axset18]: 
        axset[1].set_xticks([-400,-300]); axset[2].set_xticks([-150,-75]) # dD
        axset[3].set_xticks([-45,-35]); axset[4].set_xticks([-20,-10]) # d18O
    
        # Set y-ticks:
    for axset in [axset16,axset17,axset18]:
        for i in range(1,len(axset16)):
            axset[i].set_yticks([])
            
    for axset in [axset16,axset17,axset18]:
        axset[1].spines['right'].set_visible(False)
        axset[2].spines['left'].set_visible(False)
        axset[3].spines['right'].set_visible(False)
        axset[4].spines['left'].set_visible(False)
        
        # Legend:
    for axset in [axset16,axset17,axset18]: axset[0].legend()
    
        # x,y labels:
    #axset18.set_xlabel('q [g/kg]', fontsize=14)
    #ax6.set_xlabel('dD [permil]', fontsize=14)
    #ax9.set_xlabel('d18O [permil]', fontsize=14)
    #--------------------------------------------------------------------------


    # Insets of isotope profiles zoomed-in to the lower 1km: 
    #--------------------------------------------------------------------------
    axset_insets17 = [fig.add_axes([hl_prfax[i]+0.05, vl_prfax[1]+0.03, 0.1, 0.1])
                          for i in (1,3)]
    axset_insets18 = [fig.add_axes([hl_prfax[i]+0.05, vl_prfax[2]+0.03, 0.1, 0.1])
                      for i in (1,3)]
    
    for lat in prflats1718:
        # index of grid_lat value closest to lat.
        i = np.nanargmin(abs(grid_lat1718-lat)) 

        axset_insets17[0].plot(curtains17['dD_permil'][:,i], grid_hght)
        axset_insets17[1].plot(curtains17['d18O_permil'][:,i], grid_hght)
        axset_insets18[0].plot(curtains18['dD_permil'][:,i], grid_hght)
        axset_insets18[1].plot(curtains18['d18O_permil'][:,i], grid_hght)


    for axset in [axset_insets17,axset_insets18]:
        axset[0].set_xlim(-85,-68); axset[0].set_ylim(0,1)
        axset[1].set_xlim(-12.5,-10.5); axset[1].set_ylim(0,1)
    #--------------------------------------------------------------------------


"""
Makes some plots that I use to estimate the boundary layer top for ORACLES 
2017 and 2018, where the boundary layers were decoupled.

Draws latitude-height curtains of (1) the vertical derivative of  
temperature, (2) a binary map which is 1 if the P3 sampled a liquid water 
content >0.1 at least once in that area, and 0 otherwise.

From the plots, it looks like there is fair-good coincidence of LWC presence 
and sharp potential temperature gradients for each year. I do not have code 
to draw this region on the humidity and isotope curtains, I just visually get 
the region from the figures in this fxn, then draw on the curtains after 
saving them.
"""
def estimate_BLtop(year):
        
    # Vertical derivative of temperature:
    grid_lat1718 = np.arange(-12,0,0.2)
    grid_hght = np.arange(0,6,0.1) # height in km.
    curtains = get_curtains(get_data(year,block_secs=10), 
                            grid_lat1718, grid_hght, curtvars=['T_K'])

    grid2D_lat1718, grid2D_hght = np.meshgrid(grid_lat1718, grid_hght)
    dT_dz = np.diff(curtains['T_K'], axis=0)/np.diff(grid2D_hght, axis=0)    
    
    plt.figure()
    plt.contourf(grid_lat1718, grid_hght[:-1]*1000, dT_dz, 
                 levels=np.arange(-10,10), cmap='seismic')
                 #levels=np.arange(0,16))
    
    plt.colorbar() 
    plt.ylim(0,5000)
    plt.title('ORACLES '+year+', dT/dz', fontsize=14)   
    plt.xlabel('Latitude', fontsize=14)
    plt.ylabel('Height [m]', fontsize=14)
    

    # Presence of cloud at least once during sampling period - determined 
    # from LWC:
    lwcdata = get_data(year, block_secs=10)[['lat','height_m','King_LWC_ad']]
    lwcdata['lat'] = np.round(lwcdata['lat'])
    lwcdata['height_m'] = np.round(lwcdata['height_m']/250)*250
    lwcdata['incloud'] = lwcdata['King_LWC_ad']>0.1
    lwccount = pd.pivot_table(lwcdata, values='incloud', index=['height_m'], columns=['lat'], 
                              aggfunc=np.sum)
    
    plt.figure()
    plt.pcolor(lwccount.columns.values, lwccount.index.values, lwccount>0)
    
    plt.ylim(0,5000)    
    plt.title('ORACLES '+year+r', lwc>0.1 at least once', fontsize=14)    
    plt.xlabel('Latitude', fontsize=14)
    plt.ylabel('Height [m]', fontsize=14)

    
"""
Main function to make all plots.
"""
def make_all_plots():
    
    # Plots that I use to estimate BL top for ORACLES 2017 and 2018:
    #print('Making plots to estimate MBL top for 2017')
    #estimate_BLtop('2017')
    #print('Making plots to estimate MBL top for 2018')
    #estimate_BLtop('2018')        

    # Weight isotope ratios by water concentration when making curtains?
    h2o_weights=True
    
    # Curtains for each year:
    for year in ['2016','2017','2018']:
        
        print(year+'\n'+'=======')
        
        grid_hght = np.arange(0,6,0.1)
        if year=='2016': grid_lat = np.arange(-22,-10,0.2)
        if year in ['2017','2018']: grid_lat = np.arange(-12,0,0.2)
        
        print('Computing curtains via oversampling')
        curtains = get_curtains(get_data(year, block_secs=10),
                                grid_lat, grid_hght, h2o_weights=h2o_weights)        
        print('Plotting curtains')
        curtain_plots(year, curtains, grid_lat, grid_hght)
   
    
if __name__=='__main__': make_all_plots()






















