# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 09:50:03 2021

@author: Dean
"""


# Third party:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc


## Data paths:
path_timesync_dir = r"../WISPER_data/time_sync/"
path_pic1cal_dir = r"../WISPER_data/pic1_cal/"
path_pic2cal_dir = r"../WISPER_data/pic2_cal/"
path_p3merge = r"../Ancillary/P3_Merge_Data/"


##_____________________________________________________________________________
"""
Returns the WISPER basic-prep data as a pandas df with air pressure added. 
Input date is str 'yyyymmdd'. rev (str) is the revision suffix on the data 
files.
"""
##_____________________________________________________________________________
def data_with_pressure(date, rev):
        
    # Load WISPER data as pandas df, NAN-masked:
    if rev=='time-sync': path_datadir=path_timesync_dir    
    if rev=='pic1-cal': path_datadir=path_pic1cal_dir    
    if rev=='pic2-cal': path_datadir=path_pic2cal_dir    
        
    fname_wisper = "WISPER_%s_%s.ict" % tuple([date, rev])
    wisper = pd.read_csv(path_datadir + fname_wisper, header=0)
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


## ORACLES flight dates:
dates2016_good_mako = ['20160830','20160831','20160902','20160904']
dates2016_good_gulper = ['20160910','20160912','20160914','20160918',
                         '20160920','20160924','20160925']
dates2017_good = ['20170815','20170817','20170818','20170821','20170824',
                  '20170826','20170828','20170830','20170831','20170902']

dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
              '20181012','20181015','20181017','20181019','20181021',
              '20181023']


## Histograms of dD, d18O, and DXS in the subcloud layer for 2016 flights - 
## comparison of Mako and Gulper measurements:
##-----------------------------------------------------------------------------
"""
Returns isotope data for all passed dates. Also includes calculated DXS. 
rev=data revision suffix (str).
"""
def subcloud_iso(dates, rev, pic='2'):
    
    # Get dD, and d18O from all flights:
    subcloud_iso = pd.DataFrame({}) # Put all measurements here.
    for date in dates:
        # Load data:
        data = data_with_pressure(date, rev)      
        # Data in subcloud layer estimated as data below 950hPa height:
        data_subcloud = data.loc[data['Static_Pressure']>950]
        # Append isotope ratio measurements:
        subcloud_iso = subcloud_iso.append(data_subcloud[['dD_tot'+pic,'d18O_tot'+pic]])

    # Calc dxs:
    subcloud_iso['dxs_tot'+pic] = subcloud_iso['dD_tot'+pic] - 8*subcloud_iso['d18O_tot'+pic]

    return subcloud_iso


# Subcloud isotope measurements separately for Mako and Gulper:
#mako = subcloud_iso(dates2016_good_mako, 'pic2-cal')
#gulper_beforecal = subcloud_iso(dates2016_good_gulper, 'time-sync')
#gulper_aftercal = subcloud_iso(dates2016_good_gulper, 'pic2-cal')
def mako_gulper_subcloud_compare():
    mako = subcloud_iso(['20160831','20160904'], 'pic2-cal')
    gulper_beforecal = subcloud_iso(['20160910','20160912','20160925'], 'time-sync')
    gulper_aftercal = subcloud_iso(['20160910','20160912','20160925'], 'pic2-cal')
    
    
    # Histograms:
    plt.figure()
    ax1 = plt.subplot(1,3,1); ax2 = plt.subplot(1,3,2); ax3 = plt.subplot(1,3,3)
    
    mako['dD_tot2'].hist(ax=ax1, bins=20, histtype='step')
    gulper_beforecal['dD_tot2'].hist(ax=ax1, bins=20, histtype='step')
    gulper_aftercal['dD_tot2'].hist(ax=ax1, bins=20, histtype='step')
    
    mako['d18O_tot2'].hist(ax=ax2, bins=20, histtype='step')
    gulper_beforecal['d18O_tot2'].hist(ax=ax2, bins=20, histtype='step')
    gulper_aftercal['d18O_tot2'].hist(ax=ax2, bins=20, histtype='step')
    
    mako.loc[(mako['dxs_tot2']>-20) & (mako['dxs_tot2']<40), 'dxs_tot2'].hist(ax=ax3, bins=20, histtype='step')
    gulper_beforecal['dxs_tot2'].hist(ax=ax3, bins=20, histtype='step')
    gulper_aftercal['dxs_tot2'].hist(ax=ax3, bins=20, histtype='step')
    
    del mako; del gulper_beforecal; del gulper_aftercal        
##-----------------------------------------------------------------------------


## Histograms of Mako dD, d18O, and DXS in the subcloud layer for each 
## sampling period:
##-----------------------------------------------------------------------------
def mako_subcloud_drift():

    mako16_before = subcloud_iso(dates2016_good_mako, 'time-sync', pic='2')
    mako17_before = subcloud_iso(dates2017_good, 'time-sync', pic='1')
    mako18_before = subcloud_iso(dates2018_good, 'time-sync', pic='1')
    
    mako16_after = subcloud_iso(dates2016_good_mako, 'pic2-cal', pic='2')
    mako17_after = subcloud_iso(dates2017_good, 'pic1-cal', pic='1')
    mako18_after = subcloud_iso(dates2018_good, 'pic1-cal', pic='1')
    
    plt.figure()
    ax1 = plt.subplot(2,3,1); ax2 = plt.subplot(2,3,2);
    ax3 = plt.subplot(2,3,3); ax4 = plt.subplot(2,3,4)
    ax5 = plt.subplot(2,3,5); ax6 = plt.subplot(2,3,6)
    
    varlabels_plot = ['dD_tot','d18O_tot','dxs_tot']
    for v, axset in zip(varlabels_plot, [(ax1,ax4),(ax2,ax5),(ax3,ax6)]):
        mako16_before[v+'2'].hist(ax=axset[0], histtype='step', bins='fd')
        mako17_before[v+'1'].hist(ax=axset[0], histtype='step', bins='fd')
        mako18_before[v+'1'].hist(ax=axset[0], histtype='step', bins='fd')
    
        mako16_after[v+'2'].hist(ax=axset[1], histtype='step', bins='fd')
        mako17_after[v+'1'].hist(ax=axset[1], histtype='step', bins='fd')
        mako18_after[v+'1'].hist(ax=axset[1], histtype='step', bins='fd') 
        
    for ax in (ax1,ax4): ax.set_xlim(-90,-45)
    for ax in (ax2,ax5): ax.set_xlim(-15,-7)
    for ax in (ax3,ax6): ax.set_xlim(0,50)
##-----------------------------------------------------------------------------











































