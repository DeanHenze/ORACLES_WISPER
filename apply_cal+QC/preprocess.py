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

Example usage:
    date = '20170815'
    preprocessor = Preprocessor(date)
    preprocessor.preprocess_file(save=True)
    preprocessor.test_plots()
"""


# Built in:
import inspect
import os

# Third party:
import numpy as np # 1.19.2
import pandas as pd # 1.1.3
import matplotlib.pyplot as plt # 3.3.2


# Get the path of the directory containing this script (used in combination 
# with relative paths to locate datafiles):
filename = inspect.getframeinfo(inspect.currentframe()).filename
scriptpath = os.path.dirname(os.path.abspath(filename))


class Preprocessor(object):
    """
    Preprocesses a raw data file as outlined in this script header.
    
    Parameters
    ----------
    date: str.
        ORACLES flight date "yyyymmdd" for the datafile.
        
    rawdata: pandas.core.frame.DataFrame, default=None.
        The raw dataset for the appropriate flight date. If left as None, the 
        datafile will automatically be loaded and assigned to 'rawdata'.
            
    writeloc: str.
        Path (relative or absolute) to save the preprocessed data to.
    """
    def __init__(self, date, rawdata=None, writeloc=""):
        
        # If a pandas df was already passed:
        if type(rawdata)==pd.core.frame.DataFrame:
            self.rawdata = rawdata
            
        # If no df passed, load the datafile:
        elif rawdata is None:
            if date[:4]=='2016': 
                header=37
                fname_raw = 'WISPER_'+date+'_rawField.ict'
            if date[:4] in ['2017','2018']: 
                header=0
                fname_raw = 'WISPER_'+date+'_rawField.dat'
            rawdatadir = scriptpath + "\\WISPER_raw_data\\"
            self.rawdata = pd.read_csv(rawdatadir + fname_raw, header=header)            
            
        self.date = date
        self.writeloc = writeloc
    
        
    def preprocess_file(self, save=False):
        """
        Full preprocessing of a file, as outlined in this script header. 
        Processed data is added as new attribute "preprodata". If "save" 
        is set to True, the preprocessed data is saved to path "writeloc".
        """        
        if self.date[:4] == '2016': self.preprocess_2016file(save=save)
        if self.date[:4] in ['2017','2018']: self.preprocess_20172018file(save=save)
        
    
    def test_plots(self):
        """
        Plot WISPER quantities before and after preprocessing of a file. 
        """
        if self.date[:4] == '2016': self.test_plots_2016()
        if self.date[:4] in ['2017','2018']: self.test_plots_20172018()
        
 
    def preprocess_2016file(self, save=False):
        """
        Full preprocessing of a 2016 file. If "save" is set to True, the 
        preprocessed data is saved to path "writeloc".
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
        if save:
            self.flag_na_switch(preprodata, flag=-9999.0, nan2flag=True)
            self.preprodata.to_csv(self.writeloc, index=False)
        

    def preprocess_20172018file(self, save=False):
        """
        Full preprocessing of a 2017 or 2018 file. If "save" is set to True, 
        the preprocessed data is saved to path "writeloc".
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
        

        ## restructure pic0 data so that SDI and CVI inlet 
        ## measurements are separate columns:
        ## -------------------------------
            # split original dataframe Pic1 variables:
        pic1_keys_raw = ['pic0_qh2o','pic0_deld','pic0_delo']
        pic1_SDI = preprodata.loc[preprodata['state_valve']==0, pic1_keys_raw].copy()  
        pic1_CVI = preprodata.loc[preprodata['state_valve']==1, pic1_keys_raw].copy()
            # Rename columns:
        pic1_keys_tot = ['h2o_tot1','dD_tot1','d18O_tot1']
        pic1_SDI.rename(columns = dict(zip(pic1_keys_raw, pic1_keys_tot)), 
                        inplace=True)
        pic1_keys_cld = ['h2o_cld','dD_cld','d18O_cld']
        pic1_CVI.rename(columns = dict(zip(pic1_keys_raw, pic1_keys_cld)), 
                        inplace=True)
            # Recombine to get NANs where appropriate and add as new columns:
        preprodata = preprodata.join(
                                     pd.merge(pic1_SDI, pic1_CVI, how='outer', 
                                              left_index=True, right_index=True
                                              ),
                                     how='left'
                                     )

    
        ## Convert the raw data timestamps to seconds since midnight UTC:
        ## -------------------------------
        start_utc = np.array([])    
        for i in preprodata['timestamp']:
            clock = i[len(i)-8:len(i)]
            time_secs = int(clock[0:2])*3600 + int(clock[3:5])*60 + int(clock[6:8])
            start_utc = np.append(start_utc,time_secs)
        preprodata['Start_UTC'] = start_utc
        preprodata.drop(columns='timestamp', inplace=True)
        
        self.timestamp_checker(preprodata)


        ## Compute CVI inlet and excess flows:
        ## -------------------------------
            # excess flow, slpm:
        f_xs = self.rawdata['f_dry_slpm']  \
               - self.rawdata['state_valve']*(350/1000)  \
               - self.rawdata['f_user_slpm'] - self.rawdata['f_byp_slpm'] 
        preprodata['cvi_xsFlow'] = f_xs.round(decimals=2)   
    
            # inlet flow, slmp:
        f_CVI_inlet = self.rawdata['state_valve']*(350/1000)  \
                      + self.rawdata['f_user_slpm']  \
                      + self.rawdata['f_byp_slpm']
        preprodata['cvi_inFlow'] = f_CVI_inlet.round(decimals=2)
        
        
        ## Rename some of the columns:
        ## -------------------------------
        keys_old = ['pic1_qh2o','pic1_deld','pic1_delo',
                    'cvi_enhancement','f_user_slpm','state_valve']
        keys_new = ['h2o_tot2','dD_tot2','d18O_tot2','cvi_enhance',
                    'cvi_userFlow','wisper_valve_state']
        preprodata.rename(columns=dict(zip(keys_old, keys_new)), inplace=True)      


        ## Locate bad data and assign NANs:
        ## -----------------------------
            # Any rows where Pic1 SDI measurements are all 0 (i.e. for the few  
            # cases where Mako was shut off during flight):
        pic1_off = ( (preprodata['h2o_tot1']==0) & (preprodata['dD_tot1']==0) 
                     & (preprodata['d18O_tot1']==0) )
        preprodata.loc[pic1_off, ['h2o_tot1','dD_tot1','d18O_tot1']] = np.nan

            # All CVI variables and measurements on the CVI for the 10/10/2018 
            # flight (the CVI was not operated correctly on that day):
        if self.date=='20181010':
            preprodata.loc[:, ['h2o_cld','dD_cld','d18O_cld','cvi_lwc', 
                               'cvi_enhance', 'cvi_dcut50','cvi_inFlow', 
                               'cvi_xsFlow', 'cvi_userFlow']
                           ] = np.nan
        
            # CVI measurements where q<500ppmv (with the enhancement factor 
            # of ~30, 500ppmv corresponds to a very low amount of liquid):
        preprodata.loc[preprodata.h2o_cld<500, ['h2o_cld','dD_cld','d18O_cld']
                       ] = np.nan

            # Pic2 measurements where pressure deviations from 
            # 35 torr are greater than 0.2 torr (applicable mostly to 2018):
        preprodata.loc[ abs(self.rawdata['pic1_pcav']-35)>0.2, 
                       ['h2o_tot2', 'dD_tot2', 'd18O_tot2'] ] = np.nan
 
            # Some weird Pic1 behavior (mostly in 2017 due to the bad 
            # thermistor attachment), leading to some clearly outlying values:
        preprodata.loc[preprodata['dD_tot1']>100, pic1_keys_tot] = np.nan 
        preprodata.loc[preprodata['h2o_tot1']<0, pic1_keys_tot] = np.nan
        preprodata.loc[preprodata['dD_cld']>100, pic1_keys_cld] = np.nan 
        preprodata.loc[preprodata['h2o_cld']<0, pic1_keys_cld] = np.nan
        
            # Intervals of bad data recorded in the .txt files:
        self.applyna_from_txtfiles(preprodata)
        
        
        ## Assign attribute and optional save:
        ## ----------------------------------
        self.preprodata = preprodata
        if save:
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
            print("\t Time stamps are complete and at 1 Hz.")     
        
        else:
        # If any of the tests fail:
            print('\t Modification to time variable needed.')
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
        bad_intvls_dir = (scriptpath + "\\outlier_time_intervals\\" 
                          + d[:4] + "\\") # Path to directory of bad data intervals.
        data.set_index('Start_UTC', drop=False, inplace=True)    


        def apply_nan(data, t_badintvls, varkeys):
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
            path_txt = bad_intvls_dir + "Pic2_Outlier_Times_" + d + ".txt"
            t_badintvls = pd.read_csv(path_txt, header=0)
            apply_nan(data, t_badintvls, ['h2o_tot2','dD_tot2','d18O_tot2'])
        
        
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
                t_badintvls = pd.read_csv(bad_intvls_dir + f, header=0)
                apply_nan(data, t_badintvls, varkeys)
                
        
        data.reset_index(drop=True, inplace=True)
        

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
        """
        Plot WISPER quantities before and after preprocessing of a 2016 file. 
        """
        
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
        
        
    def test_plots_20172018(self):
        """
        Plot WISPER quantities before and after preprocessing of a 2017/18 file. 
        """
        
        raw = self.rawdata
        prepro = self.preprodata
        
        ## Time series of water concentration and isotope ratios:
        ## --------------------------
        plt.figure()
        ax1 = plt.subplot(3,1,1)
        ax2 = plt.subplot(3,1,2)
        ax3 = plt.subplot(3,1,3)
        
        # Use time variable from the preprocessed dataframe as x-axis for 
        # both raw and preprocessed plotting.
        ax1.plot(prepro['Start_UTC'], raw['pic0_qh2o'], 'k-', label='raw')
        ax1.plot(prepro['Start_UTC'], prepro['h2o_tot1'], 'ro', label='CVI')
        ax1.plot(prepro['Start_UTC'], prepro['h2o_cld'], 'bx', label='SDI')

        ax2.plot(prepro['Start_UTC'], raw['pic0_deld'], 'k-', label='raw')
        ax2.plot(prepro['Start_UTC'], prepro['dD_tot1'], 'ro', label='CVI')
        ax2.plot(prepro['Start_UTC'], prepro['dD_cld'], 'bx', label='SDI')
        
        ax3.plot(prepro['Start_UTC'], raw['pic0_delo'], 'k-', label='raw')
        ax3.plot(prepro['Start_UTC'], prepro['d18O_tot1'], 'ro', label='CVI')
        ax3.plot(prepro['Start_UTC'], prepro['d18O_cld'], 'bx', label='SDI')

        ax3.set_xlabel('Time (secs UTC midnight)', fontsize=14)
        ax1.set_ylabel('H2O, Pic1 (ppmv)', fontsize=14)
        ax2.set_ylabel('dD, Pic1 (permil)', fontsize=14)
        ax3.set_ylabel('d18O, Pic1 (permil)', fontsize=14)
        
        for ax in [ax1,ax2,ax3]: ax.legend(fontsize=12)
    

        ## Time series of CVI system flows:
        ## --------------------------
        plt.figure()
        ax4 = plt.subplot(3,1,1)
        ax4.plot(prepro['Start_UTC'], prepro['cvi_enhance'], 
                 label='cvi enhancement factor')
        ax5 = ax4.twinx()
        ax5.plot(prepro['Start_UTC'], prepro['wisper_valve_state'], 
                 'k--', label='valve state')
        ax4.legend(loc='upper right', fontsize=12); ax5.legend(loc='upper left')
        ax4.set_ylim(-5, 50)
        ax4.set_title('CVI system', fontsize=14)
        
        
        ax6 = plt.subplot(3,1,2)
        ax6.plot(prepro['Start_UTC'], prepro['cvi_inFlow'], label='cvi inlet')
        ax6.plot(prepro['Start_UTC'], prepro['cvi_userFlow'],label='user')
        ax7 = ax6.twinx()
        ax7.plot(prepro['Start_UTC'], prepro['wisper_valve_state'], 
                 'k--', label='valve state')
        ax6.legend(loc='upper right'); ax7.legend(loc='upper left')
        
        
        ax8 = plt.subplot(3,1,3)
        ax8.plot(prepro['Start_UTC'], prepro['cvi_xsFlow'], 
                 label='cvi excess')        
        ax9 = ax8.twinx()
        ax9.plot(prepro['Start_UTC'], prepro['wisper_valve_state'], 
                 'k--', label='valve state')
        ax8.set_ylim(-1, 1)
        ax8.set_xlabel('Time (secs UTC midnight)', fontsize=14)
        ax8.legend(loc='upper right'); ax9.legend(loc='upper left')        
        

## Some tester code:
#dates2016_good = ['20160830','20160831','20160902','20160904','20160910',
#                  '20160912','20160914','20160918','20160920','20160924',
#                  '20160925']

#for date in dates2016_good:
#    p = Preprocessor(date)
#    p.preprocess_2016file()
#    p.test_plots_2016()
