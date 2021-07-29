# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 12:37:49 2021

@author: Dean

Some basic restructuring and QC of the WISPER data to get it ready for 
calibration.

For 2016 ORACLES flights:
------------------------
    - Take the subset of variables from the raw data to be used in the 
      final ESPO files.
    - Process to 1 Hz frequency. 
    - Check that there is a full set of timestamps and at 1 Hz. 
    - Mask outlying data with flagged value -9999.
    
For 2017 and 2018 ORACLES flights:
(Data have already been processed to 1 Hz frequency)
----------------------------------
    - Check that there is a full set of timestamps and at 1 Hz. 
    - Take the subset of variables from the raw data to be used in the 
      final ESPO files. 
    - separate Pic1 data into separate columns for CVI measurements and 
      SDI measurements.
    - Calculate additional CVI/WISPER system flow quantities.
    - Convert timestamps from date-time to UTC-seconds.
    - Mask outlying data with flagged value -9999.

Note that for 2017 and 2018, what I am calling Pic1(Pic2) vars in this script 
are labelled as pic0(pic1) vars in the raw files.

The data are masked in two ways. First, there are several automated flagging 
criteria hardcoded in this script. Second, there are .txt files where I 
have manually looked through the data and recorded intervals of clearly bad 
data. This script reads the intervals in from those files and masks them.
"""


# Built in:
import os

# Third party:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


path_wisper_raw_dir = r"./WISPER_raw_data/"
path_wisper_prep_dir = r"./WISPER_processed_data/"
path_outlier_tintvls_dir = r"./outlier_time_intervals/"


if not os.path.isdir(path_wisper_prep_dir): os.mkdir(path_wisper_prep_dir)


class Preprocessor(object):
    """
    Preprocesses a raw data file as outlined in this script header.
    
    Parameters
    ----------
    date: str.
        ORACLES flight date "yyyymmdd" for the datafile.
        
    rawdata: optional, pandas.core.frame.DataFrame.
        The raw dataset for the appropriate flight date. If passed, rawdata is 
        assigned directly. Otherwise, the raw datafile will be loaded into a 
        pandas df and assigned during __init__. 
        
    writeloc: str
        Path (relative or absolute) to save the preprocessed data to.
    """
    
    
    def __init__(self, date, rawdata=None, writeloc=""):
        
        # If a pandas df was already passed:
        if type(rawdata)==pd.core.frame.DataFrame:
            self.rawdata = rawdata
            
        # If no df passed, load the datafile as pandas df using 'date':
        elif rawdata is None:
            if date[:4]=='2016': 
                header=37
                fname_raw = 'WISPER_'+date+'_rawField.ict'
            if date[:4] in ['2017','2018']: 
                header=0
                fname_raw = 'WISPER_'+date+'_rawField.dat'
            self.rawdata = pd.read_csv(r"./WISPER_raw_data/" + fname_raw, 
                                       header=header)            
            
        self.date = date
        self.writeloc = writeloc
        
        
    def preprocess_2016file(self, save_output=False):
        """
        Full preprocessing of a 2016 file, as outlined in this script header. 
        Processed data is added as new attribute "preprodata". If "save_output" 
        is set to True, the preprocessed data is saved to path "writeloc".
        """
        
        print("Preprocessing data for ORACLES flight date %s" % self.date)

        preprodata = self.rawdata.copy()
        self.flag_na_switch(preprodata, flag=-9999.0, flag2nan=True)

        ## Keep only a subset of the vars:
        ## -------------------------------
            
        preprodata.drop([' pressure_cavity',' temp_cavity',' temp_das'],
                        axis=1, inplace=True)
            # Rename some vars to be consistent with the other two years:
        keys_old = [' H2O_ppmv',' dD_permil',' d18O_permil']
        keys_new = ['h2o_tot2','dD_tot2','d18O_tot2']
        preprodata.rename(columns=dict(zip(keys_old, keys_new)), inplace=True)


        ## Get data into 1 Hz frequency:
        ## -----------------------------
            
        # a) group/average timestamps within the same second:
        preprodata.loc[:,'Start_UTC'] = preprodata['Start_UTC'].round(0)    
        preprodata = preprodata.groupby(by='Start_UTC', as_index=False).mean()
    
        # b) Interpolate any 1 Hz gaps: 
            # Full set of reference timestamps at 1Hz:
        t0 = preprodata.iloc[0]['Start_UTC']
        tf = preprodata.iloc[-1]['Start_UTC']
        t_ref = pd.DataFrame({'t_ref':np.arange(t0,tf+1,1)}) 
            # merge data with reference to get rows of NANs where needed:
        preprodata = preprodata.merge(t_ref, how='outer', left_on='Start_UTC', 
                                      right_on='t_ref',sort=True)   
            # Use pandas to interpolate on NAN rows:
        preprodata.interpolate(method='linear', axis=0, limit=3, inplace=True)

        # c) Check that timestamps are complete and at 1 Hz, then drop 't_ref':
        self.timestamp_checker(preprodata)
        preprodata.drop(['t_ref'], axis=1, inplace=True)        
        
        
        ## Assign NAN to bad data intervals recorded in the .txt files:
        ## -----------------------------
        self.applyna_from_txtfiles(preprodata)
        
        
        ## Assign attribute and optional save:
        ## ----------------------------------
        self.preprodata = preprodata
        if save_output:
            self.flag_na_switch(preprodata, flag=-9999.0, nan2flag=True)
            self.preprodata.to_csv(self.writeloc, index=False)
        
    
    def timestamp_checker(self, data):
        """
        Check that the timestamps column is 1 Hz with no missing times. 
        Perform the following checks:
           1) The number of rows in data1 is equal to tf-t0+1, 
           2) The number of non-nan values in 'Start_UTC' is equal to tf-t0+1, 
           3) The time difference between adjacent rows is 1s for every row.
    
        For future readers: this function was not written very well, so if 
        you are confused it is my fault not yours. It does the job though.
        """
        
        utc = 'Start_UTC'
    
        t0 = data.iloc[0][utc]
        tf = data.iloc[-1][utc]
        t_ref = np.arange(t0,tf+1,1)

        # 1):
        dt = np.array(data.iloc[1:][utc]) - np.array(data.iloc[:-1][utc])
        n_full = np.size(data,0)
        bool_full = ( n_full == (tf-t0+1) )
        # 2):
        n_nonan = np.sum(pd.notnull(data[utc]))
        bool_nonan = ( n_nonan == (tf-t0+1) )
        # 3):
        n_dt_ls1 = np.size(np.where(dt<1),1)
        n_dt_gt1 = np.size(np.where(dt>1),1)
        bool_dt = ( (n_dt_ls1+n_dt_gt1) == 0 ) 
        
        if bool_dt and bool_full and bool_nonan:
            print("Completed time averaging/interpolation to 1 Hz "
                  " without needing modifications.")     
        
        else:
        # If any of the tests fail:
            print('Modification to time variable needed.')
            print('\t tf-t0+1 = '+str(tf-t0+1))
            print('\t Number of rows in 1s-averaged data after NAN filling = '
                  +str(n_full))
            print('\t Number of non-nan values in \'Start_UTC\' = '
                  +str(n_nonan))
            print('\t Number adjacent rows where dt is not 1s AND neither of '
                  'the rows has a NAN timestamp= '+str(n_dt_ls1+n_dt_gt1))
            
            if bool_full and bool_dt:
            # In this case the averaging/interpolation worked and there were 
            # just some intervals > 3 s with missing data. Replace the time 
            # column with 't_ref':   
                print('\t Filling in missing timestamps but keeping missing value'
                      ' flags for other WISPER vars. Good to go!')
                data.loc[:,'Start_UTC'] = t_ref


    def applyna_from_txtfiles(self, data):
        """
        There are .txt files where I have manually recorded intervals of 
        clearly bad data. This fxn sets the data in those intervals to NAN for 
        the appropriate variable columns.
        
        data: handle to the pandas df to mask.
        """

        d = self.date
        p_txtfiledir = r"./outlier_time_intervals/%s/" % d[:4]

        def apply_nan(t_badintvls, varkeys):
            """
            Assigns NAN to intervals in "data" (for the columns in 
            "varkeys") contained in "t_badintvls" (pandas df). 
            """
            for i,row in t_badintvls.iterrows():
            # Each row contains start/end times for a bad data interval:
                if row['Start']==0: 
                    data.loc[:row['End'], varkeys] = np.nan
                else:
                    data.loc[row['Start']:row['End'], varkeys] = np.nan


        if d[0:4] == '2016':
        # For 2016, only one .txt file to load for each flight:  
            path_txt = p_txtfiledir + "Pic2_Outlier_Times_" + d + ".txt"
            t_badintvls = pd.read_csv(path_txt, header=0)
            apply_nan(t_badintvls, ['h2o_tot2','dD_tot2','d18O_tot2'])
        
        
        if d[0:4] in  ['2017','2018']:
        # For other years, for each flight there are separate .txt files for 
        # _tot1, _tot2, and _cld vars:
            varkeys_list = [['h2o_tot1','dD_tot1','d18O_tot1'],
                            ['h2o_cld','dD_cld','d18O_cld','cvi_lwc'],
                            ['h2o_tot2','dD_tot2','d18O_tot2']]
            fname_list = ["Pic1_Tot_Outlier_Times_" + d + ".txt",
                          "Pic1_Cld_Outlier_Times_" + d + ".txt",
                          "Pic2_Outlier_Times_" + d + ".txt"]
            for varkeys, f in zip(varkeys_list, fname_list):
                t_badintvls = pd.read_csv(p_txtfiledir + f, header=0)
                apply_nan(t_badintvls, varkeys)
        

    def flag_na_switch(self, data, flag=None, flag2nan=False, nan2flag=False):
        """
        Switch from a flag value to a NAN or vise versa for "data". Set "flag" 
        as the flag value. Set one of "nan2flag" or "flag2nan" to True.
        """

        if flag2nan: 
            data.replace(flag, np.nan, inplace=True)
        if nan2flag: 
            data.fillna(flag, inplace=True)    
            data.replace(np.inf,flag, inplace=True)
            data.replace(-np.inf,flag, inplace=True)
            
            
    def test_plots_2016(self):
        
        rawdata = self.rawdata
        preprodata = self.preprodata
        
        plt.figure()
        ax1 = plt.subplot(3,1,1)
        ax2 = plt.subplot(3,1,2)
        ax3 = plt.subplot(3,1,3)
        
        ax1.plot(rawdata['Start_UTC'], rawdata[' H2O_ppmv'], 'k-')
        ax1.plot(preprodata['Start_UTC'], preprodata['h2o_tot2'], 'ro')

        ax2.plot(rawdata['Start_UTC'], rawdata[' dD_permil'], 'k-')
        ax2.plot(preprodata['Start_UTC'], preprodata['dD_tot2'], 'ro')
        
        ax3.plot(rawdata['Start_UTC'], rawdata[' d18O_permil'], 'k-')
        ax3.plot(preprodata['Start_UTC'], preprodata['d18O_tot2'], 'ro')


    def preprocess_20172018file(self, save_output=False):
        """
        Full preprocessing of a 2017 or 2018 file, as outlined in this script 
        header. Processed data is added as new attribute "preprodata". If 
        "save_output" is set to True, the preprocessed data is saved to path 
        "writeloc".
        """
        
        print("Preprocessing data for ORACLES flight date %s" % self.date)
        
        preprodata = self.rawdata.copy()


        ## Keep only a subset of the vars:
        ## -------------------------------
            
        varskeep = ['timestamp','pic0_qh2o','pic0_deld','pic0_delo',
                    'pic1_qh2o','pic1_deld','pic1_delo','state_valve',
                    'cvi_enhancement','cvi_dcut50','cvi_lwc','f_user_slpm',]
        preprodata = preprodata[varskeep]

        self.flag_na_switch(preprodata, flag=-9999.0, flag2nan=True)        
        self.flag_na_switch(preprodata, flag=-99.0, flag2nan=True)        
        

        ## In a new dataframe, restructure pic0 data so that SDI and CVI inlet 
        ## measurements are separate columns:
        ## -------------------------------

            # split original dataframe Pic1 variables:
        pic1_varkeys_raw = ['pic0_qh2o','pic0_deld','pic0_delo']
        pic1_SDI = data0.loc[data0['state_valve']==0, pic1_varkeys_raw]    
        pic1_CVI = data0.loc[data0['state_valve']==1, pic1_varkeys_raw]
            # Rename columns:
        pic1_varkeys_tot = ['h2o_tot1','dD_tot1','d18O_tot1']
        pic1_SDI.rename(columns = dict(zip(pic1_varkeys_raw, pic1_varkeys_tot)), 
                        inplace=True)
        pic1_varkeys_cld = ['h2o_cld','dD_cld','d18O_cld']
        pic1_CVI.rename(columns = dict(zip(pic1_varkeys_raw, pic1_varkeys_cld)), 
                        inplace=True)
            # Recombine to get restructured frame:
        data1 = pd.merge(pic1_SDI, pic1_CVI, how='outer', left_index=True, 
                     right_index=True)
        
        
        
        data1['cvi_enhance'] = data0['cvi_enhancement']
    data1['cvi_dcut50'] = data0['cvi_dcut50']
    data1['cvi_lwc'] = data0['cvi_lwc']
    data1['cvi_userFlow'] = data0['f_user_slpm'] 
    data1['wisper_valve_state'] = data0['state_valve']
    
    
        


               
        
        
        ## Assign NAN to bad data intervals recorded in the .txt files:
        ## -----------------------------
        self.applyna_from_txtfiles(preprodata)
        
        
        ## Assign attribute and optional save:
        ## ----------------------------------
        self.preprodata = preprodata
        if save_output:
            self.flag_na_switch(preprodata, flag_9999=True)
            self.preprodata.to_csv(self.writeloc, index=False)
            

def dataprep_17_18(date):    
    """
    Process a raw data file for 2017 or 2018.
        
    date (str) 'yyyymmdd'. flight date.
    """   
    
    print('Preparing data for %s' % date)

    global path_wisper_raw_dir, path_wisper_prep_dir
    global path_outlier_tintvls_dir


    # Load WISPER data and change missing value flags to np.nan:
    fname_raw = 'WISPER_'+date+'_rawField.dat'
    data0 = pd.read_csv(path_wisper_raw_dir + fname_raw)
    data0.replace(-9999.0, np.nan, inplace=True)
    data0.replace(-99.0, np.nan, inplace=True)
                    

    # In a new dataframe, restructure Pic1 data so that SDI and CVI inlet 
    # measurements are separate columns:
        # split original dataframe Pic1 variables:
    pic1_varkeys_raw = ['pic0_qh2o','pic0_deld','pic0_delo']
    pic1_SDI = data0.loc[data0['state_valve']==0, pic1_varkeys_raw]    
    pic1_CVI = data0.loc[data0['state_valve']==1, pic1_varkeys_raw]
        # Rename columns:
    pic1_varkeys_tot = ['h2o_tot1','dD_tot1','d18O_tot1']
    pic1_SDI.rename(columns = dict(zip(pic1_varkeys_raw, pic1_varkeys_tot)), 
                    inplace=True)
    pic1_varkeys_cld = ['h2o_cld','dD_cld','d18O_cld']
    pic1_CVI.rename(columns = dict(zip(pic1_varkeys_raw, pic1_varkeys_cld)), 
                    inplace=True)
        # Recombine to get restructured frame:
    data1 = pd.merge(pic1_SDI, pic1_CVI, how='outer', left_index=True, 
                 right_index=True)
        

    # Add additional columns to the new dataframe from above:    
        # Get seconds since midnight UTC from the raw data timestamps:
    start_utc = np.array([])    
    for i in data0['timestamp']:
        clock = i[len(i)-8:len(i)]
        time_secs = int(clock[0:2])*3600 + int(clock[3:5])*60 + int(clock[6:8])
        start_utc = np.append(start_utc,time_secs)
    data1['Start_UTC'] = start_utc
        
        # CVI variables taken directly from raw data:
    data1['cvi_enhance'] = data0['cvi_enhancement']
    data1['cvi_dcut50'] = data0['cvi_dcut50']
    data1['cvi_lwc'] = data0['cvi_lwc']
    data1['cvi_userFlow'] = data0['f_user_slpm'] 
    data1['wisper_valve_state'] = data0['state_valve']
        
        # CVI excess flow, in slpm:
    f_xs = data0['f_dry_slpm'] - data0['state_valve']*(350/1000)              \
        - data0['f_user_slpm'] - data0['f_byp_slpm'] 
    data1['cvi_xsFlow'] = f_xs.round(decimals=2)   

        # CVI inlet flow, in slmp:
    f_CVI_inlet = data0['state_valve']*(350/1000) + data0['f_user_slpm']      \
        + data0['f_byp_slpm']
    data1['cvi_inFlow'] = f_CVI_inlet.round(decimals=2)
    
        # Pic2 measurements:
    data1['h2o_tot2'] = data0['pic1_qh2o']
    data1['dD_tot2'] = data0['pic1_deld']
    data1['d18O_tot2'] = data0['pic1_delo']
    

    # Mask bad data:
        # Any rows where Pic1 SDI measurements are all 0 (ie this would be the 
        # case if Mako were shut off during flight):
    pic1_off = ( (data1['h2o_tot1']==0) & (data1['dD_tot1']==0) 
                 & (data1['d18O_tot1']==0) )
    data1.loc[pic1_off, ['h2o_tot1','dD_tot1','d18O_tot1']] = np.nan

        # All CVI variables and measurements on the CVI for the 10/10/2018 flight. 
        # (the CVI was not operated correctly on that day):
    if date=='20181010':
        data1.loc[:, ['h2o_cld','dD_cld','d18O_cld','cvi_lwc', 
                      'cvi_enhance', 'cvi_dcut50','cvi_inFlow', 
                      'cvi_xsFlow', 'cvi_userFlow']
                  ] = np.nan

        # CVI measurements where q<500ppmv (with the enhancement 
        # factor of ~30, 500ppmv corresponds to a very low amount of liquid):
    data1.loc[data1.h2o_cld<500, ['h2o_cld','dD_cld','d18O_cld']] = np.nan

        # Pic2 (Spiny) measurements in 2018 where pressure deviations from 
        # 35 torr are greater than 0.2 torr:
    if date[0:4]=='2018':
        data1.loc[ abs(data0['pic1_pcav']-35)>0.2, 
                  ['h2o_tot2', 'dD_tot2', 'd18O_tot2'] ] = np.nan
 
        # Some weird Pic1 behavior (mostly in 2017 due to the bad 
        # thermistor attachment), leading to some clearly outlying values:
    data1.loc[data1['dD_tot1']>100, ['h2o_tot1','dD_tot1','d18O_tot1']] = np.nan   
    data1.loc[data1['h2o_tot1']<0, ['h2o_tot1','dD_tot1','d18O_tot1']] = np.nan
        
        # Additional time intervals of outlying data identified in the 
        # separate .txt files:
    data1.set_index('Start_UTC', drop=False, inplace=True)    
    
    
    def ol_txt_mask(data, path_outliers_txt, varkeys_mask):
        t_outlier = pd.read_csv(path_outliers_txt, header=0) # Load .txt file.
        for i,row in t_outlier.iterrows():
            if row['Start']==0: 
                # Only mask variables in varkeys_mask:
                data.loc[:row['End'], varkeys_mask] = np.nan
            else:
                data.loc[row['Start']:row['End'], varkeys_mask] = np.nan
        return data
        
    data1 = ol_txt_mask(data1, # Pic1 total water mask from .txt file.
                        path_outlier_tintvls_dir + 
                        r'/'+date[0:4]+r'/Pic1_Tot_Outlier_Times_'+date+'.txt', 
                        ['h2o_tot1','dD_tot1','d18O_tot1'])
        
    data1 = ol_txt_mask(data1, # Pic1 cloud water mask from .txt file.
                        path_outlier_tintvls_dir + 
                        r'/'+date[0:4]+r'/Pic1_Cld_Outlier_Times_'+date+'.txt', 
                        ['h2o_cld','dD_cld','d18O_cld','cvi_lwc']) 
        
    data1 = ol_txt_mask(data1, # Pic2 total water mask from .txt file.
                        path_outlier_tintvls_dir + 
                        r'/'+date[0:4]+r'/Pic2_Outlier_Times_'+date+'.txt', 
                        ['h2o_tot2','dD_tot2','d18O_tot2'])                 

        
    # re-order columns:
    columns = ['Start_UTC','wisper_valve_state','h2o_tot1','h2o_tot2','h2o_cld',
               'dD_tot1','dD_tot2','dD_cld','d18O_tot1','d18O_tot2',
               'd18O_cld','cvi_lwc', 'cvi_enhance', 'cvi_dcut50',
               'cvi_inFlow', 'cvi_xsFlow', 'cvi_userFlow']
    data1 = data1[columns]
            
        
    # Save after replacing nan, inf, and -inf values with the flag value -9999:
    data1.fillna(-9999, inplace=True)    
    data1 = data1.replace(np.inf,-9999)
    data1 = data1.replace(-np.inf,-9999)

    fname_prep = 'WISPER_' + date + '_basic-prep.ict'
    data1.to_csv(path_wisper_prep_dir + fname_prep, index=False)


"""
if __name__=='__main__':
##_____________________________________________________________________________
## ORACLES research flight dates for days where WISPER took good data:
    dates2016_good = ['20160830','20160831','20160902','20160904','20160910',
                      '20160912','20160914','20160918','20160920','20160924',
                      '20160925']
    dates2017_good = ['20170812','20170813','20170815','20170817','20170818',
                      '20170821','20170824','20170826','20170828','20170830',
                      '20170831','20170902']
    dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
                      '20181012','20181015','20181017','20181019','20181021',
                      '20181023']


## Run prep functions for all good days:
    print("=============================\n"
          "Starting raw data preparation\n"
          "=============================")
        
    for date in dates2016_good: dataprep_16(date)
    for date in dates2017_good: dataprep_17_18(date)
    for date in dates2018_good: dataprep_17_18(date)
    
    print("=============================\n"
          "Raw data preparation complete\n"
          "=============================")
##_____________________________________________________________________________
"""   
#dates2016_good = ['20160830','20160831','20160902','20160904','20160910',
#                  '20160912','20160914','20160918','20160920','20160924',
#                  '20160925']

#for date in dates2016_good:
#    p = Preprocessor(date)
#    p.preprocess_2016file()
#    p.test_plots_2016()
    
    
    
    
    
    
    
    
    