# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 13:48:57 2021

@author: Dean


Make latitude-altitude curtains of humidity, dD, and dxs for each ORACLES 
year. Curtains are made using kernel density estimation with a gaussian 
weighting.
"""


# Built in:
import sys

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc # 1.3.1
from warnings import filterwarnings

# Local code:
import oversampler



## Suppress a harmless warning triggered in this script:
filterwarnings(
    action='ignore', category=DeprecationWarning, 
    message='`np.bool` is a deprecated alias'
    )



def lev3product_allyears():
    """ Get level 3 curtain products for each ORACLES year and save as 
    .nc files.
    
    Returns: None
    """
    for year in ['2016','2017','2018']:
        lev3curtains = lev3_product(year)
        path_save = r"./wisper_oracles_curtains_%s.nc" % year
        lev3curtains.to_netcdf(path_save, format="NETCDF4")
        
    print("Complete.")

    

def lev3_product(year):
    """ Creates standardized xarray dataset product (ready for writing to a 
    .nc file) for water mixing ratio 
    and isotope ratio curtains.
    
    Parameters
    ----------
    year : str
        ORACLES sampling year to create curtain product for.

    Returns
    -------
    xarray dataset.
    """

    print('Creating Level 3 data product for ORACLES %s.' % year)


    ## Gridpoint latitudes, altitudes to include in dataset:        
    dalt_grid = 100 # grid altitude spacing, meters.
    dlat_grid = 0.2 # grid latitude spacing, degrees.
    alt_grid = np.arange(0, 7000 + dalt_grid, dalt_grid)
    lat_grid = np.arange(-24, 1 + dlat_grid, dlat_grid)
    
    
    ## Create xarray Dataset object with empty entries over all lats & alts 
    ## for each variable:
    coords_ds = ['alt', 'lat']
    varnames_ds = ['q', 'dD', 'd18O', 'dxs', 'sig_q', 
                   'sig_dD', 'sig_d18O', 'sig_dxs', 'n_obs'
                   ]
    
    empty_grid = -9999.*np.ones([len(alt_grid), len(lat_grid)])
    empty_dsvars = [(coords_ds, empty_grid.copy()) for v in varnames_ds]
            
    curtains_ds = xr.Dataset(
        data_vars = dict(zip(varnames_ds, empty_dsvars)),
        coords = dict(lat = lat_grid, alt = alt_grid), 
        )
    

    ## Compute curtains for all variables and assign to xr dataset:
        # Gridpoint latitudes, altitudes to compute curtains at depend on 
        # the sampling year:        
    alt_curtain = alt_grid
    if year=='2016': 
        lat_curtain = np.arange(-24, -10 + dlat_grid, dlat_grid)
    elif year in ['2017','2018']: 
        lat_curtain = np.arange(-14, 1 + dlat_grid, dlat_grid)
        # Compute:
    curtains = compute_curtains(get_data(year, block_secs=30),
                                lat_curtain, alt_curtain, h2o_weights=True
                                ) 
        # Assign:
    varnames_curt = ['h2o_gkg', 'dD_permil', 'd18O_permil', 'dxs_permil']
    latcurtainslice = slice(    # The +/-10**-5 terms are to deal with  
        lat_curtain[0]-10**-5,  # floating-point errors for xarray. This 
        lat_curtain[-1]+10**-5  # ensures that the endpoints of the slice 
        )                       # get included.
    altcurtainslice = slice(alt_curtain[0]-10**-5, alt_curtain[-1]+10**-5)
    slicedict = {'lat':latcurtainslice, 'alt':altcurtainslice}
    
    for v_ds, v_curt in zip(varnames_ds[0:4], varnames_curt):
        mean_round = np.round(curtains[v_curt]['mean'], decimals=2)
        curtains_ds[v_ds].loc[slicedict] = mean_round
        std_round = np.round(curtains[v_curt]['stdev'], decimals=2)
        curtains_ds['sig_'+v_ds].loc[slicedict] = std_round
    
    n_obs_round = np.round(curtains['h2o_gkg']['n_obs_weighted'], decimals=2)
    curtains_ds['n_obs'].loc[slicedict] = n_obs_round
    
    
    ## Add xarray dataset attributes:
        # global attributes:
    title = ("Water and isotope ratio observations from the ORACLES MISSION,"
             " averaged onto latitude-altitude grids.")
    source = ("Aircraft in-situ data processed using"
              " kernel density estimation, with Gaussian kernel.")
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
    creation_date = "2021-12-03 00:00:00"
    version = "1.0"
    if year=='2016':
        campaign = "ORACLES-1"
        period = "2016-08-30 to 2016-09-25"
    if year=='2017':
        campaign = "ORACLES-2"
        period = "2017-08-12 to 2017-09-02"
    if year=='2018':
        campaign = "ORACLES-3"
        period = "2018-09-27 to 2018-10-23"
    region = "Southeast Atlantic"
    curtains_ds.attrs = dict(
        title=title, source=source, period=period, region=region, 
        campaign=campaign, insitution=institution, 
        references=references, contact=contact, dataset_doi=dataset_doi, 
        creation_date=creation_date, version=version
                             )
 
        # coordinate attributes:
    curtains_ds.lat.attrs = dict(
        long_name = "Latitude", 
        units = "degrees north")
    curtains_ds.alt.attrs = dict(
        long_name = "GPS altitude above reference geoid", 
        units = "meters")

        # variable attributes:
    fillval = "-9999"
    curtains_ds.q.attrs = dict(
        long_name = "Total water mixing ratio", 
        units = "g/kg", 
        _FillValue = fillval)
    curtains_ds.dD.attrs = dict(
        long_name = "Total water D/H ratio, expressed in delta-notation", 
        units = "permil", 
        _FillValue = fillval)
    curtains_ds.d18O.attrs = dict(
        long_name = "Total water 18O/16O ratio, expressed in delta-notation", 
        units = "permil", 
        _FillValue = fillval)
    curtains_ds.dxs.attrs = dict(
        long_name = "Total water deuterium excess, expressed in delta-notation", 
        units = "permil", 
        _FillValue = fillval)
    curtains_ds.sig_q.attrs = dict(
        long_name = "Standard deviation in total water mixing ratio", 
        units = "g/kg", 
        _FillValue = fillval)
    curtains_ds.sig_dD.attrs = dict(
        long_name = "Standard deviation in dD", 
        units = "permil", 
        _FillValue = fillval)
    curtains_ds.sig_d18O.attrs = dict(
        long_name = "Standard deviation in d18O", 
        units = "permil", 
        _FillValue = fillval)
    curtains_ds.sig_dxs.attrs = dict(
        long_name = "Standard deviation in deuterium excess", 
        units = "permil", 
        _FillValue = fillval)
    curtains_ds.n_obs.attrs = dict(
        long_name = "Weighted number of observations used for kernel density "
                    "estimation", 
        _FillValue = fillval)
    
    
    return curtains_ds



def compute_curtains(
        data, grid_lat, grid_alt, 
        curtvars = ['h2o_gkg','dD_permil','d18O_permil', 'dxs_permil'], 
        h2o_weights=False
        ):
    """
    Returns latitude-altitude curtains. Curtains are computed by oversampling. 
    Data is passed as a pandas df and oversampled map is computed at the 
    passed grid_lat and height (1d arrays)
    
    Inputs:
        data: pandas df
        grid_lat, grid_alt: Grid latitudes and heights (1d arrays). 
        h2o_weights (default=False): If set to True, use H2O data as additional 
            weighting factors when computing the oversampled dD and d18O maps.
    """

    print("Computing curtains via oversampling.")

    curtains = {} # Put curtain means here.
    for v in curtvars:

        if (v in ['dD_permil','d18O_permil','dxs']) and h2o_weights:
            curtains[v] = oversampler.oversampler(
                data[v], data['lat'], data['height_m'], 
                grid_lat, grid_alt, w_add=data['h2o_gkg'], 
                return_stdev='yes', ffact = 0.8
                )

        else:
            curtains[v] = oversampler.oversampler(
                data[v], data['lat'], data['height_m'], 
                grid_lat, grid_alt, return_stdev='yes', 
                ffact = 0.8
                )
            
    return curtains



def get_data(year, block_secs=30, morethanx_gkg=0.2):
    """
    Returns all P-3 data for the passed ORACLES year (str). Returns the 
    following variables:
        time, latitude, longitude, altitude, h2o (g/kg), dD (permil), 
        d18O (permil), dxs (permil).
    
    Inputs:
    -------
    year: str
        'yyyy'.
    block_secs: int or None
        Default=30. Time in seconds to group/average the data. Pass None to 
        keep original 1Hz data.
    morethanx_gkg: float
        Default = 0.2. Keep only data where humidity is greater than the 
        passed value in g/kg.
    """
    
    pd.options.mode.chained_assignment = None  # Suppress false warning.
    
    
    print("Retrieving data.")
    
    
    # Some info necessary to load the datasets:
    relpath_wisper = r"../apply_cal+QC/WISPER_calibrated_data/"  
    wisper_headerline = {'2016':70, '2017':85, '2018':85}[year]
    relpath_merged = r"../apply_cal+QC/P3_merge_data/"
    merged_revnum = {'2016':'R25', '2017':'R18', '2018':'R8'}[year]
    
    
    # Load and append data for each flight:
    data = pd.DataFrame({}) # Will hold all data.
    dates = p3_flightdates(year)
    
    for date in dates:
        
        # 1) Load wisper data. 
        # 2) Get a single column for each variable filled with Pic1 instrument 
        #    data where available and Pic2 otherwise:
        # 3) Convert water concentration from ppmv units to g/kg.
        wisper = pd.read_csv(
            relpath_wisper + "WISPER_P3_%s_R2.ict" % date, 
            header=wisper_headerline
            )

        if year == '2016': # Only Pic2 data available.
            wisper_updated = wisper[['Start_UTC','h2o_tot2','dD_tot2','d18O_tot2']]
            wisper_updated['h2o_tot2'] = convert_q(wisper_updated['h2o_tot2'])
            wisper_updated.columns = ['Start_UTC','h2o_gkg','dD_permil','d18O_permil']
 
        elif year in ['2017','2018']: # Pic1 and Pic2 data available.
            wisper_updated = wisper[['Start_UTC','h2o_tot1','dD_tot1','d18O_tot1']] # Pic1 values.
            for k in wisper_updated.columns[1:]:
                k2 = k[:-1]+'2'
                inan = wisper_updated[k].isnull() # Where Pic1 has NAN
                wisper_updated.loc[inan, k] = wisper.loc[inan, k2].copy() # Replace with Pic2.
            
            wisper_updated['h2o_tot1'] = convert_q(wisper_updated['h2o_tot1'])
            wisper_updated.columns = ['Start_UTC','h2o_gkg','dD_permil','d18O_permil']
           
        
        # Load merged files as nc.Dataset object and place a subset of the 
        # vars in a pandas df:
        merged_nc = nc.Dataset(
            relpath_merged + "mrg1_P3_%s_%s.nc" % tuple([date, merged_revnum])
            )
        if year in ['2016','2017']: altitude_key='MSL_GPS_Altitude'
        if year == '2018': altitude_key='GPS_Altitude'
        merged_pd = pd.DataFrame({ # Data to use in curtains.
            'Start_UTC':merged_nc.variables['Start_UTC'][:],
            'height_m':merged_nc.variables[altitude_key][:],
            'lat':merged_nc.variables['Latitude'][:],
            'lon':merged_nc.variables['Longitude'][:],
            })
        
        
        # Average both datasets into time blocks:
        if block_secs is not None:
            wisper_t30 = np.round(wisper_updated['Start_UTC']/30)*30
            wisper_updated = wisper_updated.groupby(wisper_t30, as_index=False).mean()
            merge_t30 = np.round(merged_pd['Start_UTC']/30)*30
            merged_pd = merged_pd.groupby(merge_t30, as_index=False).mean()      
        
        
        # Combine wisper and merge data, then append:
        data = data.append(
            wisper_updated.merge(merged_pd, on='Start_UTC', how='inner'), 
            ignore_index=True, sort=False
            )
        
    
    data.replace(-9999, np.nan, inplace=True) # Change missing value flag.    


    # Only keep data where humidity is greater than passed value:
    data = data.loc[data['h2o_gkg']>morethanx_gkg]

    # Add dxs as new column:
    data['dxs_permil'] = data['dD_permil'] - 8*data['d18O_permil']
    
    
    pd.options.mode.chained_assignment = 'warn'  # Restore warning.

    
    return data



def p3_flightdates(year):
    """
    Returns P-3 flight dates to use for curtains. Returns as list of str's. 
    """
    if year=='2016':
        return ['20160831','20160902','20160904',
                '20160910','20160912','20160914',
                '20160918','20160920','20160925'
                ]
    
    elif year=='2017':
        return ['20170815','20170817','20170818',
                '20170821','20170824','20170826',
                '20170828','20170830','20170831',
                '20170902']

    elif year=='2018':
        return ['20180927','20180930','20181003','20181007',
                '20181010','20181012','20181015','20181017',
                '20181019','20181021','20181023']
    
    
    
def convert_q(q):
    """
    Converts specific humidity from units of ppmv to g/kg.
    """
    return 0.622*q/1000

        

if __name__=='__main__': 
    lev3product_allyears()






















