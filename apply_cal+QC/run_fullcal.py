# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 15:28:54 2021

@author: Dean

Once complete, will run the full preprocessing and calibration of all WISPER 
data files.
"""


# Built in:
import os

# Third party:
import numpy as np # 1.19.2

# My modules:
import preprocess 
import time_sync
import precision_columns as stdcols
import pic1_cal
import pic2_cal


# Paths to raw and (to-be) calibrated file directories:
path_rawdir = r"./WISPER_raw_data/"
path_caldir = r"./WISPER_calibrated_data/"
if not os.path.exists(path_caldir): os.makedirs(path_caldir)



# Calibration of a single file:
def calibrate_file(date):
    """
    Calibrate a WISPER file for the P3 flight on the input date ('yyyymmdd').
    """
    year = date[:4]
    
    # Preprocessing:
    pre = preprocess.Preprocessor(date) # Class instance loads the data.
    pre.preprocess_file() # Apply preprocessing.
    data_prepro = pre.preprodata.copy() # Copy to pass on for calibration.
    data_prepro.replace(-9999, np.nan, inplace=True)

    # Time synchronization:
    data_cal = time_sync.wisper_tsync(data_prepro, date)
    
    # Add precision columns:
    data_cal = stdcols.add_precision_cols(data_cal, date, test_plot=False)
    
    # Pic1 calibration (only relevant for 2017 and 2018 sampling periods):
    if year in ['2017','2018']:
        data_cal = pic1_cal.apply_cal(data_cal, date, testplots=False)
    else: 
        data_cal = data_cal
        
    # Pic2 calibration:
    data_cal = pic2_cal.apply_cal(data_cal, date, testplots=False)
    
    # Decimal rounding:
    if year == '2016':
        data_cal = data_cal.round({'Start_UTC':0, 'h2o_tot2':0, 
                                   'dD_tot2':1, 'std_dD_tot2':1, 
                                   'd18O_tot2':1, 'std_d18O_tot2':1})
    if year in ['2017','2018']:
        data_cal = data_cal.round({'Start_UTC':0, 'wisper_valve_state':0, 
                                   'h2o_tot1':0, 'h2o_tot2':0, 
                                   'h2o_cld':0, 'dD_tot1':1, 
                                   'std_dD_tot1':1, 
                                   'dD_tot2':1, 'std_dD_tot2':1, 
                                   'dD_cld':1, 'std_dD_cld':1, 
                                   'd18O_tot1':1, 'std_d18O_tot1':1,
                                   'd18O_tot2':1, 'std_d18O_tot2':1, 
                                   'd18O_cld':1, 'std_d18O_cld':1,
                                   'cvi_lwc':2, 'cvi_enhance':2, 
                                   'cvi_dcut50':2, 'cvi_inFlow':2,
                                   'cvi_xsFlow':2, 'cvi_userFlow':2})
        
    # Rearrange columns:
    if year=='2016':
        data_cal = data_cal[['Start_UTC','h2o_tot2', 'dD_tot2', 
                             'std_dD_tot2','d18O_tot2', 'std_d18O_tot2']]
    if year in ['2017','2018']:
        data_cal = data_cal[['Start_UTC', 'wisper_valve_state', 'h2o_tot1', 
                             'h2o_tot2', 'h2o_cld', 'dD_tot1', 
                             'std_dD_tot1', 'dD_tot2', 'std_dD_tot2', 
                             'dD_cld', 'std_dD_cld', 'd18O_tot1', 
                             'std_d18O_tot1', 'd18O_tot2', 'std_d18O_tot2', 
                             'd18O_cld', 'std_d18O_cld', 'cvi_lwc', 
                             'cvi_enhance', 'cvi_dcut50', 'cvi_inFlow', 
                             'cvi_xsFlow', 'cvi_userFlow']]
    
    # Save calibrated data in a temporary file:
    data_cal.fillna(-9999, inplace=True)  
    fname = "WISPER_calibrated_%s.ict" % date
    data_cal.to_csv(path_caldir + fname, index=False)
    
    # Add file header and save with final filename:
    add_file_header(date)
    
    del pre, data_cal
    
    
    
# Add file header:
def add_file_header(date):
    """
    Add a file header to a calibrated WISPER file.
    """
    
    # Get header text as list of strings from file:
    year = date[:4]
    fname_header = r'./ORACLES_WISPER_file_header_template_'+year+'.txt'  
    with open(fname_header, mode='r') as file_header:
        header = file_header.readlines() # read header text
    
    # Modify the "date of collection/revision" line in the header text:
        revDate = '20190901' # Latest revision date.
        header[6] = date[0:4]+', '+date[4:6]+', '+date[6:8]+', '                  \
                    +revDate[0:4]+', '+revDate[4:6]+', '+revDate[6:8]+'\n'
    
    # Open WISPER file and append header text:
        fname_calfile = path_caldir + "WISPER_calibrated_%s.ict" % date
        with open(fname_calfile, mode='r') as wisper:
            text_wisper = wisper.read()
        writetext = header+[text_wisper] # text to write to new file.
    
    # Write new file:
        revNum = 'R2' # Revision number on ESPO site.
        fname_final = path_caldir + 'WISPER_P3_'+date+'_'+revNum+'.ict'
        with open(fname_final, mode='w') as file_w:
            file_w.writelines(writetext)
    
    # Delete old file:
        os.unlink(fname_calfile)
    


# ORACLES research flight dates where WISPER took good data:
dates2016_good = ['20160830','20160831','20160902','20160904','20160910',
                  '20160912','20160914','20160918','20160920','20160924',
                  '20160925']
dates2017_good = ['20170812','20170813','20170815','20170817','20170818',
                  '20170821','20170824','20170826','20170828','20170830',
                  '20170831','20170902']
dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
                  '20181012','20181015','20181017','20181019','20181021',
                  '20181023']

# Apply calibration to all data:    
for date in dates2016_good: calibrate_file(date)
for date in dates2017_good: calibrate_file(date)
for date in dates2018_good: calibrate_file(date)






