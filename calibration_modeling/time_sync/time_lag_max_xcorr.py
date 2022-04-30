# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 13:01:07 2018

@author: Dean Henze, henzede@oregonstate.edu

Functions to find maximum time lagged cross correlation as a function of 
pressure between two ORACLES variables.
"""

import pandas as pd
import numpy as np
import os, inspect
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.odr import Model, RealData, ODR
from matplotlib.lines import Line2D


###____________________________________________________________________________
### Get the root path outside of directory 'subpath'
### Input: 
###     subpath: string, directory path.  
###____________________________________________________________________________
def get_root(subpath):
        # full path of this script (1st line identifies this script file):
    fname_script = inspect.getframeinfo(inspect.currentframe()).filename
    path_script = os.path.dirname(os.path.abspath(fname_script))
        # subtract script subpath to get external directory structure:
    root_ = path_script[:-len(subpath)]
    return root_
###____________________________________________________________________________


###____________________________________________________________________________
### Load wisper data and merge it with a specified subset of the merge file
### data.
### Inputs:
###         date: string 'yyyymmdd' of flight date to get data for.
###         fname_suffix: string. Data revision extension to open. I.e.
###                   'first_step', 'masked_outliers', ...
###         var_names: list of strings. Variable names in the merge nc file
###         to return.
###____________________________________________________________________________
def wisper_with_merge(date, fname_suffix, var_names):    
### Get root path to the ORACLES WISPER calibration directory: 
    path_ext = get_root(r'/Scripts')


### Depending on ORACLES year, use a different revision number filename suffix
### for the merge file:
    year = date[0:4]
    if year=='2016': revNum='R25'
    if year=='2017': revNum='R18'
    if year=='2018': revNum='R8'
    
    
### Load wisper data
    path_wisper = path_ext+r'/WISPER_Calibrated_Data/'+year+r'/' # WISPER path.
    fname_wisper = 'WISPER_'+date+'_'+fname_suffix+'.ict' # WISPER filename.
    iso = pd.read_csv(path_wisper+fname_wisper, header=0)
    iso.replace(-9999.0, np.nan, inplace=True)
    
    
### Load merge file data:
    path_merge = path_ext+r'/P3_Merge_Data/'
    fname_merge = 'mrg1_P3_'+date+'_'+revNum+'.nc'
    merge_df = pd.DataFrame({})
    merge_nc = nc.Dataset(path_merge+fname_merge)
    for name in var_names:
        merge_df[name] = merge_nc.variables[name][:]
    merge_df.replace(-9999.0, np.nan, inplace=True)

    
### Merge and return:
    return pd.merge(merge_df, iso, how='inner', on='Start_UTC')
###____________________________________________________________________________


###____________________________________________________________________________
### Method: Time-lagged cross correlation calc. ds1 shifted wrf ds2.
### Inputs: data sets 1 and 2, and the maximum shift (integer number of array elements),
### between the two to condiser for the calc.
###____________________________________________________________________________

def correlation(ds1,ds2,c_max): 
    
    import numpy as np

    ### Prelim stuff:
        # Assign data sets to numpy arrays:
    a1 = np.array(ds1); a2 = np.array(ds2)
        # Assign lengths of the arrays to l1 and l2:
    l1 = len(a1); l2 = len(a2)
        # Cut this method short if the arrays do not have the same length:
    if l1 != l2:
        print('Data sets do not have the same number of elements')
        return [],[],[]
        # cut program short if the array lengths are less than the user-specified
        # maximum offset:
    if l1 < c_max:
        print('Arrays are too short for correlation calc at given maximum '
              'offset.')
        return [],[],[]
    
    
    ### Calculate correlations for a1 shifted right wrt a2, from 0 to c_max:
    ###
        # corrR will hold the correlation calculations for right-shifts. The 
        # index of corrR corresponds to the amount by which a1 has been shifted
        # (include the 0 shift with the right shifts):
    corrR = np.zeros(c_max+1) # the +1 is to include the case of 0 shift
        # initialize array N_R, which will hold the number of non-nan samples  
        # used to calc the correlation for each right-shift:
    N_R = np.zeros(c_max+1) # the +1 is to include the case of 0 shift    
        # loop through correlation calc for each shift:
    for j in range(c_max+1):
        # subsections of a1 and a2 which overlap after shift:
        s1 = a1[0:l1-j]; s2 = a2[j:l2]
        # only use subset of s1 and s2 where both of them have non-Nan and 
        # non-INF values:
        ss1 = s1[np.isfinite(s1*s2)]
        ss2 = s2[np.isfinite(s1*s2)]
        # define n as the length of s1 and s2 after all reductions in size:
        n = len(ss1)
        # calc means and standard deviations:
        mu1 = np.nanmean(ss1); std1 = np.nanstd(ss1)
        mu2 = np.nanmean(ss2); std2 = np.nanstd(ss2)
        # calculate the correlation, and assign it to the appropriate index
        # of corrR:
        corrR[j] = np.nansum((ss1-mu1)*(ss2-mu2))/(n*std1*std2)
        # assign n to N_R:
        N_R[j] = n
        
    ### Calculate correlations for a1 shifted left wrt a2, from -c_max to -1.
    ### This is accomplished by using the same algorithm as above but now 
    ### shifting a2 right with respect to a1:
    ###
        # corrL will hold the correlation calculations for left-shifts. The 
        # 1st element of corrL contains the calc for -c_max, and the last 
        # element contains the calc for -1:
    corrL = np.zeros(c_max)
        # initialize array N_L, which will hold the number of non-nan samples  
        # used to calc the correlation for each left-shift:
    N_L = np.zeros(c_max) # the +1 is to include the case of 0 shift 
        # Loop through correlation calc for each shift:
    for j in range(0,c_max):
        # Subsection of a1 and a2. j runs from 0 to c_max-1 rather than 1 to 
        # c_max, so add 1 to j:
        s1 = a1[j+1:l1]; s2 = a2[0:l2-(j+1)] 
        # Only use subset of s1 and s2 where both of them have non-Nan and 
        # non_INF values:
        ss1 = s1[np.isfinite(s1*s2)]
        ss2 = s2[np.isfinite(s1*s2)]
        # Define n as the length of s1 and s2 after all reductions in size:
        n = len(ss1)
        # Calc means and standard deviations:
        mu1 = np.nanmean(ss1); std1 = np.nanstd(ss1)
        mu2 = np.nanmean(ss2); std2 = np.nanstd(ss2)
        # Calculate the correlation, and assign it to the appropriate index
        # of corrL (the extra -1 in the indexing brackets is to make python's
        # 0-indexing scheme work out):
        corrL[(c_max-1)-j] = np.nansum((ss1-mu1)*(ss2-mu2))/(n*std1*std2) 
        # assign n to N_L:
        N_L[(c_max-1)-j] = n
    
    # Assign to a var, c, all the shifts from -c_max to c_max performed above:
    c = np.arange(-c_max,c_max+1,1)
    # Return the shifts, correlation calc for each shift, and number of non-
    # nan samples used in the corr calc for each shift:
    return c, np.append(corrL,corrR), np.append(N_L,N_R)
###____________________________________________________________________________
  

###____________________________________________________________________________
### Function: Partition rows of a pandas df by air pressure, and find maximum 
### correlation between var1 and var2 for each pressure bin.
### Inputs:
###         data: pandas df containing the time and pressure variables as well
###               as the two variables to correlate ('data' can contain other vars).
###         P_range: List/array of pressure bin values.
###         t_var, P_var: Strings, the column header names of the time and
###             pressure variables.
###         var1, var2: Strings, the column header names of the two variables
###             to cross correlate.
###____________________________________________________________________________
def max_corr_P(data, P_range, t_var, P_var, var1, var2):         
    ###------------------------------------------------------------------------
    ### Partition data by static pressure bins and do correlation calcs.
    ###------------------------------------------------------------------------
        # static pressure var binned to nearest 50 mbar:
    data = data.assign(P_bin = np.round(data[P_var]/50)*50)
    #DataFrame.groupby(by=None, axis=0, level=None, as_index=True, sort=True, group_keys=True, squeeze=False, observed=False, **kwargs)
    
    
    ### Correlation calc on subsets of the data in 50 mbar increments:
        # store offsets, c, and correlation calcs here for each incrememnt here:
    c_list = []; corr_list = [] 
        # need to preserve the total # of rows, so save it using the time var:
    dat_i = data.index.values
    ref = pd.DataFrame({'Start_UTC':data[t_var]}, index=dat_i)
    for P in P_range:
            # subset of data in current pressure bin:
        data_P = data.loc[data['P_bin']==P]
            # fill the rest of the rows with NANs:
        data_P = ref.merge(data_P, left_on='Start_UTC', right_on=t_var, 
                           how='outer', sort=True)
        # apply correlation calc:
        c, corr, N = correlation(data_P[var1], data_P[var2], 
                                 c_max=60)
        # add offsets and calcs to lists:
        c_list = c_list+[c]; corr_list = corr_list+[corr]
    ###------------------------------------------------------------------------


    ### Returm maximum correlations for each pressure bin:
        # get indices of max correlations for each bin:
    corr_stars = np.apply_along_axis(np.nanargmax, 1, corr_list)
        # get offsets at the above indices:
    c_stars = []
    for i in range(len(c_list)):
        c_stars = c_stars+[c_list[i][corr_stars[i]]]
    return c_stars
###____________________________________________________________________________
 

###____________________________________________________________________________
### Apply 2-sigma filter: remove all data outside 2-sigma of the dataset.
###____________________________________________________________________________
def filter_2sig(a):
    a=np.array(a).astype(float)
    # Calculate mean and standard deviation
    mu = np.nanmean(a); sig = np.nanstd(a)
    # Replace data points which are 3-sig outlying with NAN:
    a[np.where(np.abs(a-mu)>2*sig)] = np.nan
    
    return a
###____________________________________________________________________________


###____________________________________________________________________________
### Find best fit line:
###____________________________________________________________________________
def fit_line(P_range, c_stars):
    # Straight line to fit:
    def line(beta,t):
        return beta[0]*t + beta[1]
    # Prepare curve to be used with odr object:
    curv = Model(line)        
    # Prepare data to be used with odr object; flatten lists and remove NAN's:
    P_range_forFit = [P_range for row in c_stars]
    P_range_forFit = np.array([item for sublist in P_range_forFit for item in sublist])
    c_stars_forFit = np.array([item for sublist in c_stars for item in sublist])
    P_range_forFit = P_range_forFit[~np.isnan(c_stars_forFit)]
    c_stars_forFit = c_stars_forFit[~np.isnan(c_stars_forFit)]
    obs = RealData(P_range_forFit, c_stars_forFit)
    # Run the ODR model and generate a fit line:
    odr = ODR(obs, curv, beta0=[1,1])
    results = odr.run()
    
    return results.beta[0], results.beta[1]
###____________________________________________________________________________
    
    
###____________________________________________________________________________
### Find fit line to pressure dependent time lag between two specified vars. 
### Inputs:
###     dates: list of strings of the form 'yyyymmdd' for the flight dates to
###         use in the calculations.
###     var1 and var2: strings, column header names in a pandas df for the 
###         variables to cross-correlate. Shifts are defined as var1 shifted
###         wrt to var2. These are either wisper headers or P3 merge file headers.
###____________________________________________________________________________
def t_shift_fit(dates, var1, var2):
### Get max correlations between Pic1 and COMA total water mixing ratio for 
### each date:
    # Pressure levels to partition data into:
    P_range = np.arange(550,1050,50)
    # Store max correlations for each date here:
    c_stars = [[] for date in dates]
    for i in range(len(dates)):
    # Load WISPER data merged with the time, pressure, COMA humidity, and cloud
    # probs lwc vars:
        data = wisper_with_merge(date=dates[i], fname_suffix='mask_outliers',
                                 var_names=['Start_UTC','Static_Pressure',
                                            'COMA_H2O_ppmv','King_LWC_ad'])
    # Correlation calc:
        c_stars[i] = max_corr_P(data, P_range, t_var='Start_UTC', 
               P_var='Static_Pressure', var1=var1, var2=var2)
    
    
    ### Plot of maximum correlations vs pressure with fit line:
    # Set font and size for plot labels:
    font = {'family' : 'serif',
            'weight' : 'normal',
            'size'   : 20}
    plt.rc('font', **font)
    # Plot and label:
    plt.figure()
    for row in c_stars:
        plt.plot(P_range, row, 'bo-')
    plt.title('Offsets w/maximum correlations by pressure bin'
              '\nORACLES 2017 - '+var1+' vs '+var2+' total water correlation analysis')
    plt.xlabel('Pressure bin [mbar]')
    plt.ylabel('Offset w/maximum correlation [s]')
    #plt.legend()
    
    
### Apply three passes of a 2-sig filter to the offset data:
    i,j = np.shape(c_stars)
    c_flat = [item for sublist in c_stars for item in sublist]
    c_filt = filter_2sig(c_flat)
    c_filt = filter_2sig(c_filt)
    c_filt = filter_2sig(c_filt).reshape((i,j))
    for row in c_filt:
        plt.plot(P_range, row, 'ko-')    
        
        
### Best fit line slope and intercept to the filtered data:
    m,b=fit_line(P_range, c_filt)
    print(m); print(b)
    plt.plot(P_range, m*P_range+b, 'r--')


### Plot legend
    custom_lines = [Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='k', lw=4),
                    Line2D([0], [0], color='r', lw=4)]
    # Round slope and intercept to two decimals:
    ms=str(np.round(m,decimals=4)); bs=str(np.round(b,decimals=2))
    plt.legend(custom_lines, ['Unfiltered','2-sig filtered','fit, y = '+ms+'x + '+bs])

###____________________________________________________________________________


###____________________________________________________________________________
### Find time shift for pic1 on the CVI via time lagged cross correlation of
### lwc measurements between pic1 and cloud probes' king probe.
###____________________________________________________________________________
def pic1_cloud_shift_fit(dates):
    # Store max correlations and average pressure for each segment of each date here:
    c_star = [[] for date in dates]
    P_mu = [[] for date in dates]
    
    
    for i in range(len(dates)):
        # Load WISPER data merged with the time, pressure, and cloud probs lwc vars:
        data = wisper_with_merge(date=dates[i], fname_suffix='mask_outliers',
                                 var_names=['Start_UTC','Static_Pressure',
                                            'King_LWC_ad'])
        data.set_index('Start_UTC', drop=False, inplace=True) 
        
        
        # If either the wisper or cloud probes lwc data is all NAN's skip this date:
        if (sum(np.isfinite(data['King_LWC_ad']))==0) or (sum(np.isfinite(data['cvi_lwc']))==0):
            c_star[i]=np.nan; P_mu[i]=np.nan
        else:
        # Isolate data where wisper was on the cvi:
            data_cvi = data.loc[data['wisper_valve_state']==1]
        # dt = difference in adjacent index (time) values, and dt>1 denotes 
        # different plumes:
            dt = data_cvi.index[1:] - data_cvi.index[:-1]
            i_bounds = np.where(dt>1)[0]
        # i_bounds is the end indicies for all plumes, so add beginning 
        # indicies (i+1 is the beginning index for the plume following the one 
        # with end index i):
            i_bounds_plus1 = i_bounds+1
            i_bounds = list(zip(i_bounds, i_bounds_plus1))
            i_bounds = [item for sublist in i_bounds for item in sublist]
        # Include indicies for very begining and very end times:
            i_bounds = np.insert(i_bounds, 0, 0)
            i_bounds = np.append(i_bounds, len(data_cvi.index)-1)
            i_bounds = [int(i) for i in list(i_bounds)]
        # Get bound times:
            t_bounds = [ 
                        [ data_cvi.index[i_bounds[i]], data_cvi.index[i_bounds[i+1]] ]
                        for i in np.arange(0,len(i_bounds),2) 
                        ]
            t_bounds = np.array(t_bounds)
        # Keep only time invervals where there is more than 1 minute of
        # data where lwc>0.02 g/kg:
            i_keep = []
            for t in t_bounds:
                test = data.loc[t[0]:t[1], 'cvi_lwc']
                test = test[np.isfinite(test)]
                if len(np.where(test>0.02)[0])>=60: i_keep = i_keep + [True]
                if len(np.where(test>0.02)[0])<60: i_keep = i_keep + [False]
            t_bounds = t_bounds[i_keep] 

            
        # For each time interval identified, find shift with the max correlation,
        # and compute the mean pressure:
            for t in t_bounds:
                c, corr, N = correlation(data.loc[t[0]:t[1], 'cvi_lwc'], 
                                         data.loc[t[0]:t[1], 'King_LWC_ad'], 
                                         c_max=35)
                i_star = np.nanargmax(corr)
                c_star[i] = c_star[i] + [c[i_star]]
                P_mu[i] = P_mu[i] + [np.nanmean(data.loc[t[0]:t[1], 'Static_Pressure'])]
                
    
### Plot results:
    c_flat = [item for sublist in c_star for item in sublist]            
    P_flat = [item for sublist in P_mu for item in sublist]  
    plt.figure()
    plt.plot(P_flat, c_flat, 'bo-')          
    
    
### Apply two passes of a 2-sig filter to the offset data and re-plot:
    c_filt = filter_2sig(c_flat)
    c_filt = filter_2sig(c_filt)
    plt.plot(P_flat, c_filt, 'ko-')          
        
    return c_star, P_mu
###____________________________________________________________________________


###____________________________________________________________________________
### Find fit line to pressure dependent time lag, c(P), for 
###     (1) Pic1 humidity wrt COMA humditiy,
###     (2) Pic2 humidity wrt COMA humditiy,
###     (3) Pic2 humidity wrt Pic1 humditiy,
###     (4) Pic1 cloud liquid water content wrt cloud probes lwc.
###____________________________________________________________________________

i_year = input('Input 1, 2, or 3 to find time shifts for 2016, 2017, or 2018 respectively: ')

### 2016 ----------------------------------------------------------------------
if i_year=='1':
### All ORACLES 2016 research flight dates. There are 4 sets: flights with
### Picarro instrument 'Mako', and 3 sets of flights for instrument 'Gulper'
### where the time lags are clearly different:
    # Dates Mako
    dates_16_m = ['0831','0902','0904']
    dates_16_m = ['2016'+date for date in dates_16_m]
    # Dates Gulper, first set
    dates_16_g1 = ['0912','0914','0918','0920','0925']
    dates_16_g1 = ['2016'+date for date in dates_16_g1]
    # Dates Gulper, second set
    dates_16_g2 = ['0910']
    dates_16_g2 = ['2016'+date for date in dates_16_g2]
    # Dates Gulper, third set
    dates_16_g3 = ['0924']
    dates_16_g3 = ['2016'+date for date in dates_16_g3]
    

### Fit lines:
    # (2), for Mako:
    #t_shift_fit(dates_16_m, var1='h2o_tot2', var2='COMA_H2O_ppmv')
    # (2), for Gulper, first set:
    #t_shift_fit(dates_16_g1, var1='h2o_tot2', var2='COMA_H2O_ppmv')
    # (2), for Gulper, second set:
    t_shift_fit(dates_16_g2, var1='h2o_tot2', var2='COMA_H2O_ppmv')
    # (2), for Gulper, third set:
    t_shift_fit(dates_16_g3, var1='h2o_tot2', var2='COMA_H2O_ppmv')
### ---------------------------------------------------------------------------


### 2017 ----------------------------------------------------------------------
if i_year=='2':
### All ORACLES 2017 research flight dates:
    dates_17 = ['0812','0813','0815','0817','0818','0821',
             '0824','0826','0828','0830','0831','0902']
    dates_17 = ['2017'+s for s in dates_17]


### Fit lines:
    # (1)  
    t_shift_fit(dates_17, var1='h2o_tot1', var2='COMA_H2O_ppmv')
    # (2)  
    #t_shift_fit(dates_17, var1='h2o_tot2', var2='COMA_H2O_ppmv')
    # (3)  
    #t_shift_fit(dates_17, var1='h2o_tot2', var2='h2o_tot1')


    c_star, P_mu = pic1_cloud_shift_fit(dates_17)
### ---------------------------------------------------------------------------
    
    
### 2018 ----------------------------------------------------------------------
if i_year=='3':
### All ORACLES 2017 research flight dates:
    dates_18 = ['0927','0930','1002','1003','1005','1007','1010','1012','1015',
                 '1017','1019','1021','1023']
    dates_18 = ['2018'+date for date in dates_18]


### Fit lines:
    # (1)  
    t_shift_fit(dates_18, var1='h2o_tot1', var2='COMA_H2O_ppmv')
    # (2)  
    t_shift_fit(dates_18, var1='h2o_tot2', var2='COMA_H2O_ppmv')
    # (3)  
    t_shift_fit(dates_18, var1='h2o_tot2', var2='h2o_tot1')


    c_star, P_mu = pic1_cloud_shift_fit(dates_18)
### ---------------------------------------------------------------------------

###____________________________________________________________________________

