# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 11:50:41 2021

@author: Dean

Contains functions to apply full calibration to Picarro 2 data in a single 
WISPER file after preprocessing and time sync.

For 2016, the calibration is similar to calibration of Picarro 1. For 2017 and 
2018, Picarro 2 is cross-calibrated to Picarro 1 (both humidity and isotope 
ratios).


Main functions to run:
---------------------
fullcal_singlefile: Main function to run. Calibrates Pic2 data for a single 
    WISPER file.


List of other Functions
-----------------
fullcal_2016file, fullcal_20172018file: Called by 'fullcal_singlefile()'.

crosscalibrate_h2o, crosscalibrate_dD_d18O: Cross calibration functions 
    called by 'fullcal_20172018file()'.
"""


# Third party:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Pic2 calibration fxns for ORACLES 2016 are taken from the pic1_cal script:
from pic1_cal import q_dep_cal, abscal_line


def fullcal_2016file(wisperdata, pic):
    """
    Calibration of Pic2 water concentration and isotope ratios for a single 
    2016 file.
    
    wisperdata: See header for 'fullcal_singlefile()' below.
        
    pic: str, either 'Mako' or 'Gulper' for the instrument used.
    
    Returns:
        wisperdata with the Pic2 water concentration modified in place.
    """
    # Model parameters for the humidity-dependence and absolute calibrations
    # are just hard-coded here. Two sets, since 2016 measurements are split 
    # between Mako and Gulper instruments:
    if pic=='Mako':
        # Parameters a, b for qdep_fitcurve:
        aD = -0.365; bD = 3.031; a18O = -0.00581; b18O = 4.961
        # Parameters slope, offset for abs cal of iso ratios:
        mD = 1.056412; kD = -5.957469; m18O = 1.051851; k18O = -1.041851
        # Fudge factor to add to k18O:
        ff=3.5; k18O = k18O + ff
        # Parameters slope, offset for abs cal of H2O:
        mq = 0.8512; kq = 0
        
    if pic=='Gulper':
        aD = 0.035; bD = 4.456; a18O = 0.06707; b18O = 1.889
        # Slopes for abs cal of iso ratios:        
        mD = 1.094037184; # kD = 2.540714192
        m18O = 1.06831472; # k18O = -7.5959267
        # Offsets are derived as outlined in the data paper appendix:
            # Histogram peaks of MBL isotope ratios during routine flights:
        pD_M = -75 # Mako dD peak, +/- 3 permil.
        pD_G = -94 # Gulper dD peak, +/- 3 permil.
        p18O_M = -11.5 # Mako d18O peak, +/- 0.5 permil.
        p18O_G = -16.7 # Gulper d18O peak, +/- 0.5 permil.
            # Offsets derived from peaks and cal slopes:
        kD = pD_M - mD*pD_G
        k18O = p18O_M - m18O*p18O_G
        # Parameters for abs cal of H2O:
        mq = 0.9085; kq = 0
    
    
    # Apply calibrations (Cal fxns are taken from the pic1_cal script):
        # Humidity dependence corrections:
    wisperdata['dD_tot2'] = q_dep_cal(wisperdata['dD_tot2'], 
                                      wisperdata['h2o_tot2'], aD, bD)
    wisperdata['d18O_tot2'] = q_dep_cal(wisperdata['d18O_tot2'], 
                                        wisperdata['h2o_tot2'], a18O, b18O)
        # Isotope ratio absolute calibration:
    wisperdata['dD_tot2'] = abscal_line(wisperdata['dD_tot2'], mD, kD)
    wisperdata['d18O_tot2'] = abscal_line(wisperdata['d18O_tot2'], m18O, k18O)    
        # Humidity absolute calibration:
    wisperdata['h2o_tot2'] = abscal_line(wisperdata['h2o_tot2'], mq, kq)

    return wisperdata


def crosscalibrate_h2o(wisperdata, year):
    """
    Cross-calibrates Pic2 water concentration for ORACLES 2017 and 2018.
    
    wisperdata, year: See header for 'fullcal_20172018file()' below.

    Returns:
        wisperdata with the Pic2 water concentration modified in place.
    """
    # Load slope for cross-calibration line:
    dir_pars = r"../calibration_modelling/pic2-pic1_cross-cal/"
    q_xcalslopes = pd.read_csv(dir_pars + "h2o_xcal_results.csv")
    slope = q_xcalslopes.loc[q_xcalslopes['year']==year, 'slope'].values
    
    # Return calibrated h2o:
    wisperdata['h2o_tot2'] = slope*wisperdata['h2o_tot2']


def crosscalibrate_dD_d18O(wisperdata, year):
    """    
    Cross-calibrates Pic2 dD and d18O for ORACLES 2017 and 2018. Uses the 
    results of the polynomial fit derived elsewhere. 
    
    wisperdata, year: See header for 'fullcal_20172018file()' below.
    
    Returns:
        wisperdata with the Pic2 isotope ratios modified in place.
    """
    # Make a separate dataframe of variables for both the dD and d18O 
    # cross-cal models:
    logq = np.log(wisperdata['h2o_tot2'])
    df_xcalvars = pd.DataFrame({'logq': logq,
                                'dD': wisperdata['dD_tot2'],
                                'd18O': wisperdata['d18O_tot2'],
                                'logq*dD': logq*wisperdata['dD_tot2'],
                                'logq*d18O': logq*wisperdata['d18O_tot2'],
                                })
    
    
    # Load cross-calibration model parameters and pandas.DataFrame's:
    dir_pars = r"../calibration_modelling/pic2-pic1_cross-cal/"
    pars_d18O = pd.read_csv(dir_pars + ("d18O_xcal_results_%i.csv" % year))
    pars_dD = pd.read_csv(dir_pars + ("dD_xcal_results_%i.csv" % year))


    def iso_crosscal(df_xcalvars, modelpars):
        """
        Computes and returns the cross-calibrated data. A linear combo of the 
        xcalvars raised to the appropriate power and mulitplied by the 
        corresponding coefficient in 'modelpars'.
        """
        terms = [] # Will store each term (var^pow*coeff) in the linear combo.
    
        for k in modelpars.keys():
            if k=='const':
                terms.append(modelpars[k]*np.ones(np.shape(df_xcalvars)[0]))
            else:
                # Get the xcal-var name and power it needs to be raised to:
                varname, power = k.split('^')
                # Compute term:
                terms.append(modelpars[k]*df_xcalvars[varname]**int(power))
        
        return np.sum(terms, axis=0)
    
    
    # Apply cross-calibration and return:
    wisperdata['dD_tot2'] = iso_crosscal(df_xcalvars, pars_dD)
    wisperdata['d18O_tot2'] = iso_crosscal(df_xcalvars, pars_d18O)
    return wisperdata


def fullcal_20172018file(wisperdata, year):
    """
    Calibration of Pic2 water concentration and isotope ratios for a 2017 or 
    2018 file.

    wisperdata: See header for 'fullcal_singlefile()' below.
        
    year: int. ORACLES year, either 2017 or 2018.
    
    Returns:
        wisperdata with the Pic2 measurements modified in place.
    """
    wisperdata = crosscalibrate_dD_d18O(wisperdata, year)
    wisperdata = crosscalibrate_h2o(wisperdata, year)
    return wisperdata


def fullcal_singlefile(wisperdata, date):
    """
    Calibration of Pic2 water concentration and isotope ratios for a single file.

    wisperdata: pandas.DataFrame. 
        WISPER data with Pic2 measurements. Should have all -9999.0 flags 
        replaced with np.nan.
        
    date: str. 
        ORACLES flight data 'yyyymmdd' corresponding to wisperdata.
    
    Returns:
        wisperdata with the Pic2 measurements modified in place.
    """    
    if date[:4] == '2016': 
        if date in ['20160830','20160831','20160902','20160904']:
            pic = 'Mako'
        else:
            pic = 'Gulper'
        wisperdata = fullcal_2016file(wisperdata, pic)
        
    elif date[:4] in ['2017','2018']:
        wisperdata = fullcal_20172018file(wisperdata, int(date[:4]))
        
    return wisperdata
    
    


"""
Cross-calibration of H2O, dD, d18O for 2017 or 2018.
"""
##_____________________________________________________________________________
def pic2_cal1718(year, dates, testplots=False):
    ## Load cross-cal fit parameters:
    relpath_caltable = r'../Calibration_Data/calibration_fits.xlsx'
    params_h2o = pd.read_excel(relpath_caltable, sheet_name='pic2pic1_xcal_h2o')
    params_dD = pd.read_excel(relpath_caltable, sheet_name='pic2pic1_xcal_dD')
    params_d18O = pd.read_excel(relpath_caltable, sheet_name='pic2pic1_xcal_d18O')
        # Set df indexes to year:
    params_h2o.set_index('year', inplace=True)
    params_dD.set_index('year', inplace=True)
    params_d18O.set_index('year', inplace=True)


    ## Cross-cal formulas:
    ##-------------------------------------------------------------------------
    """
    Takes in:
        (1) the Pic2 H2O data from a flight, 
        (2) the H2O cross-cal parameter fits (loaded above) 
    and returns the corrected Pic2 H2O. 
    """
    def xcal_h2o(q, p):
        return p['h2o_tot2']*q # Just a line passing through origin.
        
    """
    Takes in:
        (1) the Pic2 (logq, dD) data from a flight, 
        (2) the dD cross-cal parameter fits (loaded above) 
    and returns the corrected Pic2 dD. 
    """
    def xcal_dD(logq, dD, p):
        return (p['logq^3']*logq**(3) + p['logq^2']*logq**(2)
                + p['logq^1']*logq**(1) + p['dD^1']*dD**(1) 
                + p['dD^2']*dD**(2) + p['dD^3']*dD**(3)
                + p['(dD*logq)^1']*(dD*logq)**(1) 
                + p['(dD*logq)^2']*(dD*logq)**(2)
                + p['(dD*logq)^3']*(dD*logq)**(3)
                + p['const']
                )
    
    """
    Takes in:
        (1) the Pic2 (logq, d18O) data from a flight, 
        (2) the d18O cross-cal parameter fits (loaded above) 
    and returns the corrected Pic2 d18O. 
    """
    def xcal_d18O(logq, d18O, p):
        return (p['logq^3']*logq**(3) + p['logq^2']*logq**(2)
                + p['logq^1']*logq**(1) + p['d18O^1']*d18O**(1) 
                + p['d18O^2']*d18O**(2) + + p['d18O^3']*d18O**(3)
                + p['(d18O*logq)^1']*(d18O*logq)**(1) 
                + p['(d18O*logq)^2']*(d18O*logq)**(2)
                + p['(d18O*logq)^3']*(d18O*logq)**(3)
                + p['const']
                )
    ##-------------------------------------------------------------------------


    ## Apply calibrations:
    ##-------------------------------------------------------------------------
    for date in dates:
        print(date)
        
        # Load WISPER data:
        fname = "WISPER_%s_pic1-cal.ict" % date
        data = pd.read_csv(path_pic1cal_dir + fname)
        data.replace(-9999, np.nan, inplace=True) # Fill missing flag with NAN.
        
        # Corrected H2O:
        h2o_corrected = xcal_h2o(data['h2o_tot2'], dict(params_h2o.loc[year]))
        
        # Calibrations dD and d18O:
        dD_corrected = xcal_dD(np.log(data['h2o_tot2']), 
                               data['dD_tot2'].rolling(window=10).mean(), 
                               dict(params_dD.loc[year])
                               )
        d18O_corrected = xcal_d18O(np.log(data['h2o_tot2']), 
                                   data['d18O_tot2'].rolling(window=10).mean(), 
                                   dict(params_d18O.loc[year])
                                   )
        
        # Optional test plots:
        if testplots:
            # Pic1 H2O vs Pic2 H2O before and after cross-cal:
            plt.figure()
            plt.plot(data['h2o_tot2'], data['h2o_tot1'], 'bo')
            plt.plot(h2o_corrected, data['h2o_tot1'], 'ro')
            plt.plot(np.linspace(200,20000,100), np.linspace(200,20000,100), 'k-') # 1-1 line
        
            # Isotope ratio time series before and after cross-cal:
            vars2smooth = ['dD_tot1','dD_tot2','d18O_tot1','d18O_tot2']
            data[vars2smooth] = data[vars2smooth].rolling(window=10).mean()    
                
            plt.figure()
            plt.subplot(2,1,1)
            plt.plot(data['Start_UTC'], data['dD_tot1'], 'k')
            plt.plot(data['Start_UTC'], data['dD_tot2'], 'b')
            plt.plot(data['Start_UTC'], dD_corrected, 'r')
            plt.subplot(2,1,2)
            plt.plot(data['Start_UTC'], data['d18O_tot1'], 'k')
            plt.plot(data['Start_UTC'], data['d18O_tot2'], 'b')
            plt.plot(data['Start_UTC'], d18O_corrected, 'r')
        
    
        # Save calibrated data:
        data['h2o_tot2'] = h2o_corrected
        data['dD_tot2'] = dD_corrected
        data['d18O_tot2'] = d18O_corrected
        fname_save = "WISPER_%s_pic2-cal.ict" % date
        data.to_csv(path_pic2cal_dir+fname_save)
    ##-------------------------------------------------------------------------
##_____________________________________________________________________________
            

"""
Inputs:
    dates: list of strings 'yyyymmdd'.
    pic: Either 'Mako', or 'Gulper'
"""
##_____________________________________________________________________________
def pic2_cal16(dates, pic):
        
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
    
    
    ## Get parameters for calibration formulas:
    if pic=='Mako':
        # Parameters a, b for qdep_fitcurve:
        aD = -0.365; bD = 3.031; a18O = -0.00581; b18O = 4.961
        # Parameters for abs cal of iso ratios:
        mD = 1.056412; kD = -5.957469; m18O = 1.051851; k18O = -1.041851
        # Fudge factor to add to k18O:
        ff=3.5; k18O = k18O + ff
        # Parameters for abs cal of H2O:
        mq = 0.8512; kq = 0
        
    if pic=='Gulper':
        aD = 0.035; bD = 4.456; a18O = 0.06707; b18O = 1.889
        # Slopes for abs cal of iso ratios:        
        mD = 1.094037184; # kD = 2.540714192
        m18O = 1.06831472; # k18O = -7.5959267
        # Offsets are derived as outlined in the data paper appendix:
            # Histogram peaks of MBL isotope ratios during routine flights:
        pD_M = -75 # Mako dD peak, +/- 3 permil.
        pD_G = -94 # Gulper dD peak, +/- 3 permil.
        p18O_M = -11.5 # Mako d18O peak, +/- 0.5 permil.
        p18O_G = -16.7 # Gulper d18O peak, +/- 0.5 permil.
            # Offsets derived from peaks and cal slopes:
        kD = pD_M - mD*pD_G
        k18O = p18O_M - m18O*p18O_G
        # Parameters for abs cal of H2O:
        mq = 0.9085; kq = 0
        
        
    ## Apply calibrations and save as new datafiles:
    for date in dates:
        print(date)

        # Load WISPER data:
        fname = "WISPER_%s_time-sync.ict" % date
        data = pd.read_csv(path_timesync_dir + fname)
        data.replace(-9999, np.nan, inplace=True) # Fill missing flag with NAN.
            
        # Humidity dependence corrections:
        data['dD_tot2'] = q_dep_cal(data['dD_tot2'], data['h2o_tot2'], aD, bD)
        data['d18O_tot2'] = q_dep_cal(data['d18O_tot2'], data['h2o_tot2'], a18O, b18O)
        
        # Isotope ratio absolute calibration:
        data['dD_tot2'] = abscal_line(data['dD_tot2'], mD, kD)
        data['d18O_tot2'] = abscal_line(data['d18O_tot2'], m18O, k18O)
            
        # Humidity absolute calibration:
        data['h2o_tot2'] = abscal_line(data['h2o_tot2'], mq, kq)
        
        # Save:
        fname_save = "WISPER_%s_pic2-cal.ict" % date
        data.to_csv(path_pic2cal_dir+fname_save)

##_____________________________________________________________________________


dates2016_good_mako = ['20160830','20160831','20160902','20160904']
dates2016_good_gulper = ['20160910','20160912','20160914','20160918',
                         '20160920','20160924','20160925']
dates2017_good = ['20170815','20170817','20170818','20170821','20170824',
                  '20170826','20170828','20170830','20170831','20170902']

dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
              '20181012','20181015','20181017','20181019','20181021',
              '20181023']


print("========================\n"
      "Pic2 calibration for 2016\n"
      "========================\n")
#pic2_cal16(dates2016_good_mako, 'Mako')
#pic2_cal16(dates2016_good_gulper, 'Gulper')

print("========================\n"
      "Pic2 calibration for 2017\n"
      "========================\n")
#pic2_cal1718(2017, dates2017_good, testplots=False)
print("========================\n"
      "Pic2 calibration for 2018\n"
      "========================\n")
#pic2_cal1718(2018, dates2018_good, testplots=False)

















