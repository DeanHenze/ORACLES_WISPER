# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:46:16 2021

@author: Dean

Full calibration of Picarro 1 data for all ORACLES flights where Picarro 1 
was taking good quality data.
"""


# Built in:
import os

# Third party:
import numpy as np
import pandas as pd
   

"""
Apply humidity dependence calibration.
"""
def q_dep_cal(deltavals, qvals, a, b):

# Formula for humidity dependence correction, for either dD or d18O. Input 
# log of humidity q in ppmv:    
    def qdep_correction(logq, a, b):
        return a*(np.log(50000)-logq)**b

# Apply correction:
    correction = qdep_correction(np.log(qvals), a, b)
    return deltavals - correction
    
  
"""
Formula for absolute calibration of humidity or iso ratios. Just a line. x is 
either humidity, dD, or d18O. Output is corrected humidity, dD, d18O resp.:
"""
def abscal_line(x, m, k):
    return m*x + k

    
## Parameters a, b for qdep_fitcurve for ORACLES 2017 and 2018:
aD = [-0.438,-0.331] # Values for 17, 18
bD = [2.184,2.856]
a18O = [-0.01327,-0.00646]
b18O = [3.71,4.577]
    
    
## Parameters for abs cal of iso ratios (same values for 2017 and 2018):
mD = 1.056412478; kD = -5.957469671
m18O = 1.051851852; k18O = -1.041851852
    # Fudge factors to add to k18O, for ORACLES 2017 and 2018:
ff = [1.25,0]


## Parameters for abs cal of humidity (same values for 2017 and 2018):
mq = 0.8512; kq = 0
    

## Path to time-synchronized data, ready for Pic1 cal:    
path_timesync_dir = r"../WISPER_data/time_sync/"
## Path to place Pic1-calibrated data in:
path_pic1cal_dir = r"../WISPER_data/pic1_cal/"

if not os.path.isdir(path_pic1cal_dir): os.mkdir(path_pic1cal_dir)


## ORACLES dates where Pic1 took good quality data:
dates2017_good = ['20170815','20170817','20170818','20170821','20170824',
                  '20170826','20170828','20170830','20170831','20170902']
dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
                  '20181012','20181015','20181017','20181019','20181021',
                  '20181023']


print("=============================\n"
      "Starting Picarro 1 calibration\n"
      "=============================")


## Calibrate 2017 flights:
##-----------------------------------------------------------------------------
for date in dates2017_good:
    
    print(date)
    
# Load data:
    fname_tsync = "WISPER_%s_time-sync.ict" % date
    data = pd.read_csv(path_timesync_dir + fname_tsync)
    data.replace(-9999, np.nan, inplace=True)

    
# Humidity-dependence correction for dD:
    data['dD_tot1'] = q_dep_cal(data['dD_tot1'], data['h2o_tot1'], aD[0], bD[0])
    data['dD_cld'] = q_dep_cal(data['dD_cld'], data['h2o_cld'], aD[0], bD[0])

# Humidity-dependence correction for d18O:
    data['d18O_tot1'] = q_dep_cal(data['d18O_tot1'], data['h2o_tot1'], 
                                  a18O[0], b18O[0])
    data['d18O_cld'] = q_dep_cal(data['d18O_cld'], data['h2o_cld'], 
                                 a18O[0], b18O[0])
    
# Abs cal for dD and d18O:
    data[['dD_tot1','dD_cld']] = abscal_line(data[['dD_tot1','dD_cld']], 
                                             mD, kD)
    data[['d18O_tot1','d18O_cld']] = abscal_line(data[['d18O_tot1','d18O_cld']], 
                                                 m18O, k18O + ff[0]) 
    
# Abs cal for humidity and cloud lwc:
    data[['h2o_tot1','h2o_cld']] = abscal_line(data[['h2o_tot1','h2o_cld']], 
                                               mq, kq)
    data['cvi_lwc'] = 0.622*data['h2o_cld']/1000/data['cvi_enhance']
    
# Save after replacing nan with flag value -9999:
    data.fillna(-9999, inplace=True)    
    fname_save = 'WISPER_' + date + '_pic1-cal.ict'
    data.to_csv(path_pic1cal_dir + fname_save, index=False)
##-----------------------------------------------------------------------------
    
    
## Calibrate 2018 flights:
##-----------------------------------------------------------------------------
for date in dates2018_good:
 
    print(date)

# Load data:
    fname_tsync = "WISPER_%s_time-sync.ict" % date
    data = pd.read_csv(path_timesync_dir + fname_tsync)
    data.replace(-9999, np.nan, inplace=True)

    
# Humidity-dependence correction for dD:
    data['dD_tot1'] = q_dep_cal(data['dD_tot1'], data['h2o_tot1'], aD[1], bD[1])
    data['dD_cld'] = q_dep_cal(data['dD_cld'], data['h2o_cld'], aD[1], bD[1])

# Humidity-dependence correction for d18O:
    data['d18O_tot1'] = q_dep_cal(data['d18O_tot1'], data['h2o_tot1'], 
                                  a18O[1], b18O[1])
    data['d18O_cld'] = q_dep_cal(data['d18O_cld'], data['h2o_cld'], 
                                 a18O[1], b18O[1])
    
# Abs cal for dD and d18O:
    data[['dD_tot1','dD_cld']] = abscal_line(data[['dD_tot1','dD_cld']], 
                                             mD, kD)
    data[['d18O_tot1','d18O_cld']] = abscal_line(data[['d18O_tot1','d18O_cld']], 
                                                 m18O, k18O + ff[1]) 
    
# Abs cal for humidity and cloud lwc:
    data[['h2o_tot1','h2o_cld']] = abscal_line(data[['h2o_tot1','h2o_cld']], 
                                               mq, kq)
    data['cvi_lwc'] = 0.622*data['h2o_cld']/1000/data['cvi_enhance']
    
# Save after replacing nan with flag value -9999:
    data.fillna(-9999, inplace=True)    
    fname_save = 'WISPER_' + date + '_pic1-cal.ict'
    data.to_csv(path_pic1cal_dir + fname_save, index=False)
##-----------------------------------------------------------------------------


print("=============================\n"
      "Picarro 1 calibration complete\n"
      "=============================")