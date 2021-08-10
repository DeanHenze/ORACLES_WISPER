# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 15:28:54 2021

@author: Dean

Once complete, will run the full preprocessing and calibration of all WISPER 
data files.
"""


import preprocess 
import time_sync


# Paths to raw and (to-be) calibrated file directories:
path_rawdir = r"./WISPER_raw_data/"
path_caldir = r"./WISPER_calibrated_data/"


# ORACLES research flight dates for days where WISPER took good data:
dates2016_good = ['20160830','20160831','20160902','20160904','20160910',
                  '20160912','20160914','20160918','20160920','20160924',
                  '20160925']
dates2017_good = ['20170812','20170813','20170815','20170817','20170818',
                  '20170821','20170824','20170826','20170828','20170830',
                  '20170831','20170902']
dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
                  '20181012','20181015','20181017','20181019','20181021',
                  '20181023']


# Calibration of a single file:
def calibrate_file(date):
    """
    Calibrate a WISPER file for the P3 flight on the passed date ('yyyymmdd').
    """
    
    # Preprocessing procedure:
    pre = preprocess.Preprocessor(date)
    pre.preprocess_file()
    
    # Time synchronization:
    data_syncd = time_sync.wisper_tsync(pre.preprodata, date)
    #time_sync.test_plot(pre.preprodata, data_syncd, date)
    
    # Save calibrated data:
    fname = "WISPER_calibrated_%s.ict" % date
    data_syncd.to_csv(path_caldir + fname, index=False)
    

for date in dates2016_good: calibrate_file(date)
for date in dates2017_good: calibrate_file(date)
for date in dates2018_good: calibrate_file(date)






