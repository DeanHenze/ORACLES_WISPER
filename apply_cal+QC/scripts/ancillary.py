# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 09:21:14 2021

@author: Dean

Contains some ancillary information useful for the other scripts.
    Lists of ORACLES flight dates.
    Paths to data files.
"""


## Paths to data files:
path_basicprep_dir = r"../WISPER_data/basic_prep/"
path_timesync_dir = r"../WISPER_data/time_sync/"
path_p3merge = r"../Ancillary/P3_Merge_Data/"


## ORACLES research flight dates for days where WISPER took good data:
dates2016_good_mako = ['20160830','20160831','20160902','20160904']
dates2016_good_gulper = ['20160910','20160912','20160914','20160918',
                         '20160920','20160924','20160925']
dates2017_good = ['20170815','20170817','20170818','20170821','20170824',
                  '20170826','20170828','20170830','20170831','20170902']
dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
                  '20181012','20181015','20181017','20181019','20181021',
                  '20181023']


## ORACLES research flight dates where WIPSER took some lower quality data:
dates2017_loQC = ['20170812','20170813']


## Get paths:
#def get_datapaths(towhat):
    
    #datapaths = 