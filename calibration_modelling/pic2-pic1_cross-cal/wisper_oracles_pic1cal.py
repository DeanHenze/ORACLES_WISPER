# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 15:28:54 2021

@author: Dean

Performs preprocessing and full Pic1 calibration for all 2017 and 2018 WISPER 
files to use for the Pic2-Pic1 cross-calibration modelling. Calibrated files 
are placed in a folder within the 'pic2-pic1_cross-cal' directory.

Run fxn "calibrate_20172018_allfiles()" to get all 2017 and 2018 files 
calibrated for Pic1. Saves results in directory "path_pic1caldir".

Functions list
--------------
calibrate_20172018_file: 
    Run this to get calibrated data for a single ORACLES date.
    
calibrate_20172018_allfiles:
    Run this to get calibrated data for all ORACLES dates.
"""

# Built in:
import sys 
import os

# my calibration scripts:
if r'../../apply_cal+QC/' not in sys.path: 
    sys.path.insert(0, r'../../apply_cal+QC/')
import preprocess 
import time_sync
import pic1_cal


# Paths to raw and (to-be) calibrated file directories:
path_rawdir = r"../../apply_cal+QC/WISPER_raw_data/"
path_pic1caldir = r"./WISPER_pic1cal/"
if not os.path.isdir(path_pic1caldir): os.mkdir(path_pic1caldir)


# ORACLES flight dates where WISPER took good data:
dates2017_good = ['20170812','20170813','20170815','20170817','20170818',
                  '20170821','20170824','20170826','20170828','20170830',
                  '20170831','20170902']
dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
                  '20181012','20181015','20181017','20181019','20181021',
                  '20181023']


# Calibration of a single file:
def calibrate_20172018_file(date):
    """
    Calibrate a WISPER file for the P3 flight on the passed date ('yyyymmdd').
    """
    # Preprocessing procedure:
    pre = preprocess.Preprocessor(date)
    pre.preprocess_file()
    
    # Time synchronization:
    data_syncd = time_sync.wisper_tsync(pre.preprodata, date)
    #time_sync.test_plot(pre.preprodata, data_syncd, date)
    
    # Pic1 calibration:
    data_pic1cal = pic1_cal.apply_cal(data_syncd, date)
        
    # Save calibrated data:
    fname = "WISPER_pic1cal_%s.ict" % date
    data_pic1cal.to_csv(path_pic1caldir + fname, index=False)

    del pre, data_syncd, data_pic1cal
    

def calibrate_20172018_allfiles():
    for date in dates2017_good: calibrate_20172018_file(date)
    for date in dates2018_good: calibrate_20172018_file(date)
    