# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 10:47:32 2021

@author: Dean

Cross calibration of Pic2 (Gulper) to Pic1 (Mako). This is done for all 2017 
and 2018 flights - where there is data - other than Aug. 12th and 13th, 2017. 
For those two flights, Mako's calibration was clearly off.

Notes:
    I get the statsmodels message:
        "The condition number is large, 8.17e+03. This might indicate that there are
        strong multicollinearity or other numerical problems." 
    when I have terms logq, logq**2 and a constant. While I don't get the message 
    when I remove the constant, the R2 goes down by 0.5!
"""


# Built in:
import os

# Third party:
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
from openpyxl import load_workbook

# My modules:
import ancillary as anc

"""
path_data_dir = r"../WISPER_data/pic1_cal/"


dates2017_good = ['20170815','20170817','20170818','20170821','20170824',
                  '20170826','20170828','20170830','20170831','20170902']
dates2018_good = ['20180927','20180930','20181003','20181007','20181010',
                  '20181012','20181015','20181017','20181019','20181021',
                  '20181023']



## Load all WISPER (q,dD,d18O) data for ORACLES dates into pandas df:
##-----------------------------------------------------------------------------
#dates = dates2017_good
dates = dates2018_good
fnames = ['WISPER_'+d+'_pic1-cal.ict' for d in dates]
paths_dates = [path_data_dir+f for f in fnames]

columns = ['h2o_tot1','h2o_tot2','dD_tot1','dD_tot2','d18O_tot1','d18O_tot2']
wisper = pd.DataFrame({}, columns=columns) # Append all data here.
for p in paths_dates:
    data_temp = pd.read_csv(p) # Load.
    data_temp.replace(-9999, np.nan, inplace=True)
    data_temp = data_temp.rolling(window=10).mean() # Running mean.
    wisper = wisper.append(data_temp[columns], ignore_index=True)

wisper = wisper.dropna(how='any') # Drop missing values
##-----------------------------------------------------------------------------
"""

"""
Load all WISPER (q,dD,d18O) data for either 2017 or 2018 with good data. 
Return as a pandas df.
"""
def get_wisperdata(year):

    path_data_dir = r"../WISPER_data/pic1_cal/"

    if year=='2017':
        dates_good = ['20170815','20170817','20170818','20170821','20170824',
                      '20170826','20170828','20170830','20170831','20170902']
    elif year=='2018':
        dates_good = ['20180927','20180930','20181003','20181007','20181010',
                      '20181012','20181015','20181017','20181019','20181021',
                      '20181023']
            
    fnames = ['WISPER_'+d+'_pic1-cal.ict' for d in dates_good]
    paths_dates = [path_data_dir+f for f in fnames] # Paths to each data file.
    
## Loop through each date and append data to a single pandas df:
    columns = ['h2o_tot1','h2o_tot2','dD_tot1',
               'dD_tot2','d18O_tot1','d18O_tot2'] # Keep these data columns.
    wisper = pd.DataFrame({}, columns=columns) # Append all data here.
    for p in paths_dates:
        data_temp = pd.read_csv(p) # Load.
        data_temp.replace(-9999, np.nan, inplace=True)
        data_temp = data_temp.rolling(window=10).mean() # Running mean.
        wisper = wisper.append(data_temp[columns], ignore_index=True)
    
    return wisper.dropna(how='any') # Drop missing values


"""
Fit a line to Pic1 vs Pic2 vapor concentration. Line is constrained to pass
through the origin. 'df' is the pandas dataframe of wisper measurements.
"""
def q_xcal(df):
    return sm.OLS(df['h2o_tot1'], df['h2o_tot2'], missing='drop').fit()
    
    
"""
Produces powers of predictor vars for a candidate model to use with 
statsmodels. 

predictvars (pandas dataframe): Contains predictor vars for model.
nord (dict): keys are a subset of the columns in predictvars. Elements are 
    2-tuples for min, max powers for that key. Powers can be negative.
"""
def get_xcal_model(predictvars, nord):

    modelvars = pd.DataFrame({}) # Will hold all vars to use with model.
    
    for key in nord.keys():
        # Powers (excluding 0) of predictor var to try:
        powers = list(range(nord[key][0],nord[key][1]+1))
        if 0 in powers: del powers[powers.index(0)]
        
        for i in powers: # Compute and append powers:
            modelvars[key+'^%i' % i] = predictvars[key]**i
            
    # Add constant offset var:
    modelvars['const'] = np.ones(len(predictvars))
            
    return modelvars
    

"""
Get cross-calibration formula for either dD or d18O.

Fit Pic1 isotope measurements to a polynomial of Pic2 log(q) and isotope 
measurements, including an interaction term. Use a weighted linear regression 
from the statsmodels package (weights are the water concentrations). 

'df' is the pandas dataframe of wisper measurements. 'iso' is the isotope to 
get cross-cal formula for (can be either 'dD' or 'd18O').
"""
def iso_xcal(df, iso):

    # Pandas df of predictor vars:
    interterm = np.log(df['h2o_tot2'])*df[iso+'_tot2'] # Interaction term.
    predictvars = pd.DataFrame({'logq':np.log(df['h2o_tot2']),
                                iso:df[iso+'_tot2'],
                                '('+iso+'*logq)':interterm})
    
    # Get all desired powers of each predictor var; add as new df columns:
    nord = {'logq':(0,3), iso:(0,3), '('+iso+'*logq)':(0,3)}
    predictvars_poly = get_xcal_model(predictvars, nord)
    
    # Return model fit:
    return sm.WLS(df[iso+'_tot1'], predictvars_poly, missing='drop', 
                      weights=df['h2o_tot2']).fit()

"""
## Pic1 vs Pic2 humidity scatter plot:
##-----------------------------------------------------------------------------
fig_q = plt.figure()
ax_q = plt.axes()
ax_q.scatter(wisper['h2o_tot2'], wisper['h2o_tot1'], s=1)
##-----------------------------------------------------------------------------


## Pic2 (logq, dD)-space scatter plot, colored by Pic2 dD. Likewise for d18O:
##-----------------------------------------------------------------------------
fig_dD = plt.figure()
ax1_dD = plt.subplot(1,2,1)
ax2_dD = plt.subplot(1,2,2)
s = ax1_dD.scatter(np.log(wisper['h2o_tot2']), wisper['dD_tot2'], 
                c=wisper['dD_tot1'], vmin=-650, vmax=20)
fig_dD.colorbar(s, ax=ax1_dD)    

fig_d18O = plt.figure()
ax1_d18O = plt.subplot(1,2,1)
ax2_d18O = plt.subplot(1,2,2)
s = ax1_d18O.scatter(np.log(wisper['h2o_tot2']), wisper['d18O_tot2'], 
                     c=wisper['d18O_tot1'], vmin=-80, vmax=20)
fig_d18O.colorbar(s, ax=ax1_d18O) 
##-----------------------------------------------------------------------------
"""

"""
Get cross-calibration formula fit parameters for water concentration and both 
isotopologues for an ORACLES year (either '2017' or '2018').
"""
def get_fits(year):
        
    wisperdata = get_wisperdata(year)
    
    # Results are returned as statsmodels results objects:
    model_q = q_xcal(wisperdata)
    model_dD = iso_xcal(wisperdata, 'dD')
    model_d18O = iso_xcal(wisperdata, 'd18O')
    
    # Print results:
    print("****************************************************\n"
          "****************************************************\n"
          "Cross-calibration fit parameters for ORACLES "+year+"\n"
          "****************************************************\n"
          "****************************************************\n")
    print(model_q.summary())
    print(model_dD.summary())
    print(model_d18O.summary())
    
    # Return parameter fits:
    return {'q':model_q.params, 'dD':model_dD.params, 'd18O':model_q.params}
    


## Look at fit parameter covariances and correlations:
#cov = cov = model.cov_params().values # Covariance matrix.
#Dinv = np.diag(1 / np.sqrt(np.diag(cov))) 
#corr = (Dinv @ cov) @ Dinv
#print(corr)


## Publication-ready plot of 2D cross-cal maps:
##-----------------------------------------------------------------------------    
"""
    # Cross calibration formulas: 
def xcal_dD(logq, dD, p):
    return (p['logq^3']*logq**(3) + p['logq^2']*logq**(2)
            + p['logq^1']*logq**(1) + p['dD^1']*dD**(1) 
            + p['dD^2']*dD**(2) + p['dD^3']*dD**(3)
            + p['(dD*logq)^1']*(dD*logq)**(1) 
            + p['(dD*logq)^2']*(dD*logq)**(2)
            + p['(dD*logq)^3']*(dD*logq)**(3)
            + p['const']
            )

def xcal_d18O(logq, d18O, p):
    return (p['logq^3']*logq**(3) + p['logq^2']*logq**(2)
            + p['logq^1']*logq**(1) + p['d18O^1']*d18O**(1) 
            + p['d18O^2']*d18O**(2) + + p['d18O^3']*d18O**(3)
            + p['(d18O*logq)^1']*(d18O*logq)**(1) 
            + p['(d18O*logq)^2']*(d18O*logq)**(2)
            + p['(d18O*logq)^3']*(d18O*logq)**(3)
            + p['const']
            )

    # Compute cross-calibration values on (logq, delta) grids:
logq1=np.log(100); logq2=np.log(22000); dD1=-450; dD2=0; d18O1=-55; d18O2=0
logq_grid, dD_grid = np.meshgrid((np.linspace(logq1, logq2, 200)), 
                                 np.linspace(dD1, dD2, 100))
logq_grid, d18O_grid = np.meshgrid((np.linspace(logq1, logq2, 200)), 
                                 np.linspace(d18O1, d18O2, 100))
dD_cal = xcal_dD(logq_grid, dD_grid, dict(model_dD.params))
d18O_cal = xcal_d18O(logq_grid, d18O_grid, dict(model_d18O.params))

    # Plot:
fig_pub = plt.figure()
ax1_pub = plt.subplot(1,2,1)    
ax2_pub = plt.subplot(1,2,2)

pc_dD = ax1_pub.scatter(logq_grid.flatten(), dD_grid.flatten(), 
                        c=dD_cal.flatten())
pc_d18O = ax2_pub.scatter(logq_grid.flatten(), d18O_grid.flatten(), 
                          c=d18O_cal.flatten())

fig_pub.colorbar(pc_dD, ax=ax1_pub)
fig_pub.colorbar(pc_d18O, ax=ax2_pub)

ax1_pub.set_xlim(logq1, logq2); ax1_pub.set_ylim(dD1, dD2)
ax2_pub.set_xlim(logq1, logq2); ax2_pub.set_ylim(d18O1, d18O2)

ax1_pub.set_xlabel('Pic2 log(q[ppmv])', fontsize=14)
ax1_pub.set_ylabel(r'Pic2 $\delta$D [permil]', fontsize=14)
ax2_pub.set_xlabel('Pic2 log(q[ppmv])', fontsize=14)
ax2_pub.set_ylabel(r'Pic2 $\delta^{18}$O [permil]', fontsize=14)    
"""
##-----------------------------------------------------------------------------    


## Save fit results to my calibration fits excel file:
##-----------------------------------------------------------------------------    
"""
Save pandas dataframe (data) as a new excel sheet (with name 'newsheetname'), 
to a workbook (at path='path_workbook'). If workbook doesn't exist, create it.
"""
def add_excel_sheet(path_workbook, data, newsheetname):
    
    writer = pd.ExcelWriter(path_workbook, engine='openpyxl')

    if os.path.isfile(path_workbook): # Workbook exists. Add sheets.
        book = load_workbook(path_workbook)
        writer.book = book
        
        # ExcelWriter uses writer.sheets to access the sheet. If you leave it 
        # empty it will not overwrite an existing sheet with the same name.
        writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        
        data.to_excel(writer, newsheetname, index=False)
        writer.save()

    else: # Workbook doesn't exist. Create it and save dataframe.     
        data.to_excel(writer, newsheetname, index=False)
        writer.save()


#relpath_caltable = r'../Calibration_Data/calibration_fits.xlsx' 
#add_excel_sheet(relpath_caltable, 
#                pd.DataFrame(dict(model_dD.params), index=[1]), 
#                'pic2pic1_xcal_dD')
#add_excel_sheet(relpath_caltable, 
#                pd.DataFrame(dict(model_d18O.params), index=[1]), 
#                'pic2pic1_xcal_d18O')
#add_excel_sheet(relpath_caltable, 
#                pd.DataFrame(dict(model_q.params), index=[1]), 
#                'pic2pic1_xcal_q')

get_fits('2017')
#get_fits('2018')
##-----------------------------------------------------------------------------






## q cross-cal
##-----------------------------------------------------------------------------
#modelvars = pd.DataFrame({'q':wisper['h2o_tot2']})
#model_q = sm.WLS(wisper['h2o_tot1'], modelvars, missing='drop', 
#                  weights=wisper['h2o_tot2']).fit()
#ypred = model_q.predict()
#print(model_q.summary())
#ax_q.scatter(wisper['h2o_tot2'], ypred, c='k', s=1)
##-----------------------------------------------------------------------------


   
    
"""
## dD cross-cal:
##-----------------------------------------------------------------------------
interterm = np.log(wisper['h2o_tot2'])*wisper['dD_tot2']
predictvars = pd.DataFrame({'logq':np.log(wisper['h2o_tot2']),
                            'dD':wisper['dD_tot2'],
                            '(dD*logq)':interterm})
nord = {'logq':(0,3), 'dD':(0,3), '(dD*logq)':(0,3)}
modelvars = get_xcal_model(predictvars, nord)

model_dD = sm.WLS(wisper['dD_tot1'], modelvars, missing='drop', 
                  weights=wisper['h2o_tot2']).fit()
ypred = model_dD.predict()
print(model_dD.summary())

s2 = ax2_dD.scatter(np.log(wisper['h2o_tot2']), wisper['dD_tot2'], 
                 c=ypred, vmin=-650, vmax=20)    
fig_dD.colorbar(s2, ax=ax2_dD)   
##-----------------------------------------------------------------------------


## d18O cross-cal:
##-----------------------------------------------------------------------------
interterm = np.log(wisper['h2o_tot2'])*wisper['d18O_tot2']
predictvars = pd.DataFrame({'logq':np.log(wisper['h2o_tot2']),
                            'd18O':wisper['d18O_tot2'],
                            '(d18O*logq)':interterm})
nord = {'logq':(0,3), 'd18O':(0,3), '(d18O*logq)':(0,3)}
modelvars = get_xcal_model(predictvars, nord)

model_d18O = sm.WLS(wisper['d18O_tot1'], modelvars, missing='drop', 
                   weights=wisper['h2o_tot2']).fit()
ypred = model_d18O.predict()
print(model_d18O.summary())

s2 = ax2_d18O.scatter(np.log(wisper['h2o_tot2']), wisper['d18O_tot2'], 
                      c=ypred, vmin=-80, vmax=20)    
fig_d18O.colorbar(s2, ax=ax2_d18O)   
##-----------------------------------------------------------------------------
"""




