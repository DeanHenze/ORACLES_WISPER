# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 17:11:38 2021

@author: Dean

Uncertainty estimation for the WISPER ORACLES measurements.
"""

# Built in:
from itertools import product as product
    
# Third party:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import statsmodels.api as sm

# My modules:
import cal_formulas


"""
Use monte carlo methods to generate a spread in predictions for a model.
    
Parameter space (dim N) is explored by sampling from N normal distributions 
centered on the respective parameter means and with standard deviations equal 
to the respective parameter standard errors. Random-normal errors may also be 
added to the input variables on each iteration (e.g. measurement errors).

Inputs:
    model (function, callable): The model to run monte carlo with. The model 
        should take as arguments all input variables first, followed by 
        model parameters. If the model takes each parameter as a separate 
        argument, set 'unpack_params' to True. Otherwise it is assumed that 
        the model takes all params in a single list/array.
        
    inputvars (array-like): First dimension is the number of model input 
        variables. Subsequent dimensions are the values of each variable. 
        E.g. for 3 input variables that each have 100 values, inputvars is 
        3x100. For 3 inputs vars that are each shape 50x100, inputvars is 
        3x50x100. In the case of the latter where each variable is 
        multidimensional, make sure the passed model is capable of handling 
        multidimensional inputs.
        
    param_means (array-like): Expected values for the parameters.
    
    sig_params (array-like): Standard uncertainties in the parameters. Same 
        length as param_means. Parameters can have uncertainties of 0.
        
    iterations (int): Number of model realizations.
    
    unpack_params: Default=False. Keep as False if the model takes all 
        parameters in a single list/array. Set to True if the models takes 
        parameters as separate arguments.
        
    sig_inputvars (optional array-like): Standard uncertainties in input vars 
        (e.g. measurement uncertainty). If passed, should be same shape as 
        inputvars.
"""
def mc_normalsampler(model, inputvars, param_means, sig_params, 
                     iterations, unpack_params=False, sig_inputvars=None):
    
    # Make sure some things are definitey np.arrays:
    inputvars = np.array(inputvars)
    if sig_inputvars is not None: sig_inputvars = np.array(sig_inputvars)
    
    
    # Model output will be stored in a 2d array; the rows are separate model 
    # realizations:
    if len(np.shape(inputvars))==1: inputvars = [inputvars]
    modout = np.zeros( np.append(iterations, np.shape(inputvars[0])) )

    
    # If no uncertainties for the input variables were opted for, make a 
    # bunch of 0's. I.e. 'errors' of 0 will be added for each mc iteration:
    if sig_inputvars is None:
        sig_inputvars = np.zeros(np.shape(inputvars))

    
    # Monte carlo iteration. How I delt with the unpack_params split is 
    # clunky but reduces the number of if-else statements during runtime:
    if unpack_params:
        for i in range(iterations):
            # Random-normal sample of the model parameters:
            psample = np.random.normal(loc=param_means, scale=sig_params)
            # Add random-normal errors to the input vars:
            eps = np.random.normal(loc=0, scale=sig_inputvars)
            inputvars_witherror = inputvars + eps
            # Run the model:
            modout[i,:] = model(*inputvars_witherror, *psample)
    
    else: # Same code block as above except model takes psample with no *.
        for i in range(iterations):
            psample = np.random.normal(loc=param_means, scale=sig_params)
            eps = np.random.normal(loc=0, scale=sig_inputvars)
            inputvars_witherror = inputvars + eps
            modout[i,:] = model(*inputvars_witherror, psample)
        
        
    # Mean and standard dev of model results at each set of values for 
    # the input vars:
    modelmean = np.mean(modout, axis=0)
    modelstd = np.std(modout, axis=0)
    
    return np.array([modelmean, modelstd, modout], dtype=object)


"""
2016 calibration parameters will be used to get uncertainty estimates. However, 
the uncertainties derived using the 2016 parameters are assumed to be 
representative for all ORACLES years.
"""
def pic1_uncertainties():
        
    """
    Runs full calibration of either dD or d18O for Pic1:
        d: dD or d18O values (np.array like).
        q: humidity values (np.array like, same shape as dD).
        calparams: 4-element list/array-like. First 2 elements are the 'a' and 
            'b' parameters in the humidity dependence formula. Second 2 elements 
            are the slope and intercept of the absolute calibration.
    """
    def pic1_isoratio_fullcal(d, q, calparams):
        
        d_qdep_correct = cal_formulas.q_dep_cal(d, q, 
                                                calparams[0], calparams[1])
        d_abscal = cal_formulas.abscal_line(d_qdep_correct, 
                                            calparams[2], calparams[3])
        return d_abscal
            
     
    """
    Following two functions compute uncertainties for dD and d18O 
    respectively. The results are plotted. Then, a polynomial to the 
    uncertainties as a fxn of q and the respective isotope ratio.
    
    params: 4-element list/array-like. First 2 elements are the 'a' and 
        'b' parameters in the humidity dependence formula. Second 2 elements 
        are the slope and intercept of the absolute calibration. 
    sig_params: 4-element list/array-like; uncertainties for params, as stdevs.
    sig_inputvars: default=False. If True, include measurement 
        uncertainties (i.e. instrument precision) in the total 
        computation of WISPER uncertainties.
    ax_cax: 2-tuple of matplotlib.pyplot axes. 1st is for plotting the monte 
        carlo derived standard deviations, 2nd is for the color scale.
        
    Returns the polynomial fit parameters for the monte carlo derived 
    uncertainties.
    """    
    def sigma_with_fit_dD(params, sig_params, sig_inputvars=False, 
                          ax_cax=None):
    ##-------------------------------------------------------------------------
        ## Monte Carlo computation of uncertainties over a regularly-spaced 
        ## humidity (ppmv) vs. dD (permil) map:
        q, d = np.meshgrid(np.linspace(1000,15000,100), np.linspace(-300,-60,150))
            
        if sig_inputvars: # Include instrument precisions?
            sig_q = np.zeros(np.shape(q)) # q is precise enough to ignore.
            sig_dD = cal_formulas.dD_precision(q)
            results = mc_normalsampler(pic1_isoratio_fullcal, [d,q], 
                           params, sig_params, 
                           8000,
                           sig_inputvars=[sig_dD, sig_q]
                           )
        else:
            results = mc_normalsampler(pic1_isoratio_fullcal, [d,q], 
                           params, sig_params, 
                           8000,
                           )
         
            
        ## Optional plot of Monte Carlo standard devs:  
        if ax_cax is not None:
            p = ax_cax[0].scatter(q, d, c=results[1], cmap='gist_ncar', 
                                  vmin=2, vmax= 10)
            plt.colorbar(p, cax=ax_cax[1], orientation='horizontal')
            
            
        ## Fit standard deviation map to polynomial fxn of q and isoratio:
            # Put vars in a pandas df to use with the statsmodel package:
        df_forfit = pd.DataFrame({'logq':np.log(q.flatten()),
                                  'dD':d.flatten(),
                                  'sig_dD':results[1].flatten()})
            
            # Add columns for powers of q, dD. Add column for constant offset:
        for p in (2,3):
            df_forfit['logq^'+str(p)] = df_forfit['logq']**p
            df_forfit['dD^'+str(p)] = df_forfit['dD']**p
        df_forfit['const'] = np.ones(len(df_forfit))

            # Run linear regression with statsmodels
        #predictorvars = ['const','logq','logq^2','dD']
        predictorvars = ['const','logq','logq^2','logq^3','dD']
        fit = sm.OLS(df_forfit['sig_dD'], df_forfit[predictorvars], missing='drop').fit()            
        #print(fit.summary())     
        #print(fit.params)
        print(np.round(fit.rsquared, decimals=2))


        return fit.params 
    ##-------------------------------------------------------------------------
                
            
    def sigma_with_fit_d18O(params, sig_params, sig_inputvars=False, 
                            ax_cax=None):
    ##-------------------------------------------------------------------------
        ## Monte Carlo computation of uncertainties over a regularly-spaced 
        ## humidity (ppmv) vs. d18O (permil) map:
        q, d = np.meshgrid(np.linspace(1000,15000,100), np.linspace(-80,-9,150))
            
        if sig_inputvars: # Include instrument precisions?
            sig_q = np.zeros(np.shape(q)) # q is precise enough to ignore.
            sig_d18O = cal_formulas.d18O_precision(q)
            results = mc_normalsampler(pic1_isoratio_fullcal, [d,q], 
                           params, sig_params, 
                           8000,
                           sig_inputvars=[sig_d18O, sig_q]
                           )
        else:
            results = mc_normalsampler(pic1_isoratio_fullcal, [d,q], 
                           params, sig_params, 
                           8000,
                           )
            

        ## Optional plot of Monte Carlo standard devs:  
        if ax_cax is not None:
            p = ax_cax[0].scatter(q, d, c=results[1], cmap='gist_ncar', 
                                  vmin=0.2, vmax=8)
            plt.colorbar(p, cax=ax_cax[1], orientation='horizontal')
            
            
        ## Fit standard deviation map to polynomial fxn of q and isoratio:
            # Put vars in a pandas df to use with the statsmodel package:
        df_forfit = pd.DataFrame({'logq':np.log(q.flatten()),
                                  'd18O':d.flatten(),
                                  'sig_d18O':results[1].flatten()})
            # Add columns for powers of q, dD. Add column for constant offset:
        for p in (2,3):
            df_forfit['logq^'+str(p)] = df_forfit['logq']**p
            df_forfit['d18O^'+str(p)] = df_forfit['d18O']**p 
        df_forfit['const'] = np.ones(len(df_forfit))

            # Run linear regression with statsmodels
        #predictorvars = ['const','logq','logq^2','d18O']
        predictorvars = ['const','logq','logq^2','logq^3','d18O']
        fit = sm.OLS(df_forfit['sig_d18O'], df_forfit[predictorvars], missing='drop').fit()            
        #print(fit.summary())  
        #print(fit.params)
        print(np.round(fit.rsquared, decimals=2))
        
        
        return fit.params 
    ##-------------------------------------------------------------------------
    

    fig = plt.figure(figsize=(6.5,2.75))
    ax_sigD = fig.add_axes([0.1,0.175,0.38,0.6])
    cax_sigD = fig.add_axes([0.1,0.9,0.38,0.05])
    ax_sig18O = fig.add_axes([0.6,0.175,0.38,0.6])
    cax_sig18O = fig.add_axes([0.6,0.9,0.38,0.05])
   
    
    ## dD uncertainties:
    ##-------------------------------------------------------------------------
    ## Collect expected values for calibration parameters into a list.
        # Load humidity-dependence calibration parameter fits for 2016:
    relpath_caltable = r'../Calibration_Data/calibration_fits.xlsx' 
    pars_qdep = pd.read_excel(relpath_caltable, sheet_name='Mako_delta(q)')
    pars_qdep_16 = pars_qdep.loc[pars_qdep['year']==2016]
    pars_qdep_D = pars_qdep_16[['aD','bD']].values[0] 
        # Hard code absolute calibration parameters:
    m_D = 1.0564; k_D = -5.957469671
        # Collect in list to pass to one of my above fxns:
    pars_all_D = np.append(pars_qdep_D, [m_D,k_D]) 


    ## Run monte carlo simulation for 3 sets of parameter + measurement 
    ## uncertainty combos:
    
        # (1) Uncertainties when using the 1Hz WISPER measurements for 
        # relative comparisons:    
    sigpars_qdep_D = pars_qdep_16[['sig_aD','sig_bD']].values[0] 
    sig_mD = 0.0564/2; sig_kD = 0 # No offset needed for relative comparisons.
    sigpars_all_D = np.append(sigpars_qdep_D, [sig_mD,sig_kD])
    sigWISP1_D = sigma_with_fit_dD(pars_all_D, sigpars_all_D, sig_inputvars=True, 
                                   ax_cax=(ax_sigD, cax_sigD))
    
        # (2) Uncertainties when averaging data to 0.1Hz or lower, or 
        # comparing PDFs:
    sigWISP2_D = sigma_with_fit_dD(pars_all_D, sigpars_all_D)
    
        # (3) Uncertainties when comparing WISPER 0.1Hz data or PDFs to other 
        # datasets, or to absolute theoretical values:
    sig_kD = 5. # Now we care about absolute offset.
    sigpars_all_D = np.append(sigpars_qdep_D, [sig_mD,sig_kD])
    sigWISP3_D = sigma_with_fit_dD(pars_all_D, sigpars_all_D)
    ##-------------------------------------------------------------------------
    
    
    ## d18O uncertainties:
    ##-------------------------------------------------------------------------
    ## Collect expected values for calibration parameters into a list.
        # Humidity-dependence calibration parameter fits for 2016:
    pars_qdep_18O = pars_qdep_16[['a18O','b18O']].values[0] 
        # Hard code absolute calibration params:
    m_18O = 1.051851852; k_18O = -1.041851852
        # Collect in list to pass to one of my above fxns:
    pars_all_18O = np.append(pars_qdep_18O, [m_18O,k_18O])             


    ## Run monte carlo simulation for 3 sets of parameter + measurement 
    ## uncertainty combos:
    
        # (1) Uncertainties when using the 1Hz WISPER measurements for 
        # relative comparisons:
    sigpars_qdep_18O = pars_qdep_16[['sig_a18O','sig_b18O']].values[0] 
    sig_m18O = 0.05185/2; sig_k18O = 0 # No offset needed for rel comparisons.
    sigpars_all_18O = np.append(sigpars_qdep_18O, [sig_m18O,sig_k18O])
    sigWISP1_18O = sigma_with_fit_d18O(pars_all_18O, sigpars_all_18O, 
                                       sig_inputvars=True, 
                                       ax_cax=(ax_sig18O, cax_sig18O))
    
        # (2) Uncertainties when averaging data to 0.1Hz or lower, or 
        # comparing PDFs:
    sigWISP2_18O = sigma_with_fit_d18O(pars_all_18O, sigpars_all_18O)

        # (3) Uncertainties when comparing WISPER 0.1Hz data or PDFs to other 
        # datasets, or to absolute theoretical values:
    sig_kD = 1/2. # Now we care about absolute offset.
    sigpars_all_18O = np.append(sigpars_qdep_18O, [sig_m18O,sig_k18O])
    sigWISP3_18O = sigma_with_fit_d18O(pars_all_18O, sigpars_all_18O)
    ##-------------------------------------------------------------------------
    
    
    ## Collect WISPER uncertainty map parameter fits into a pandas df:
    idx_levnames = ('use case','isotope')
    idx_labs = (['1','2','3'],['dD','d18O'])
    multi_idx = pd.MultiIndex.from_product(idx_labs, names=idx_levnames)
    sigWISP_df = pd.DataFrame(np.zeros([6,5]),
                              index=multi_idx, 
                              columns=('alph1','alph2','alph3','alph4','alph5')
                              )
    
    for case, pD, p18O in zip(['1','2','3'],
                        [sigWISP1_D, sigWISP2_D, sigWISP3_D], 
                        [sigWISP1_18O, sigWISP2_18O, sigWISP3_18O]):
        sigWISP_df.loc[(case,'dD')] = pD.values
        sigWISP_df.loc[(case,'d18O')] = p18O.values
        
    sigWISP_df = sigWISP_df.round(dict(zip(sigWISP_df.columns, [1,1,2,4])))
    print(sigWISP_df)       
 
    ## Compute WISPER uncertainties for some typical q, dD, d18O values in 
    ## the MBL, cloud layer, BB-loaded free tropo, and clean free tropo:
        # Fit fxn for either dD or d18O. 
    def wisper_sigfit(q, d, pars):
        return pars[0] + pars[1]*np.log(q) + pars[2]*np.log(q)**2 + \
               pars[3]*np.log(q)**3 + pars[4]*d
 
    print('MBL dD values\n=========')
    q_mbl = 15000; dD_mbl = -70; d18O_mbl = -10
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('1','dD')].values))
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('2','dD')].values))
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('3','dD')].values))
   
    print('BB-plume dD values\n=========')
    q_mbl = 6000; dD_mbl = -70; d18O_mbl = -10
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('1','dD')].values))
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('2','dD')].values))
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('3','dD')].values))
    
    print('Clean FT dD values\n=========')
    q_mbl = 3000; dD_mbl = -70; d18O_mbl = -10
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('1','dD')].values))
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('2','dD')].values))
    print(wisper_sigfit(q_mbl, dD_mbl, sigWISP_df.loc[('3','dD')].values))
    
    print('MBL d18O values\n=========')
    q_mbl = 15000; dD_mbl = -70; d18O_mbl = -10
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('1','d18O')].values))
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('2','d18O')].values))
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('3','d18O')].values))
   
    print('BB-plume dD values\n=========')
    q_mbl = 6000; dD_mbl = -100; d18O_mbl = -10
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('1','d18O')].values))
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('2','d18O')].values))
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('3','d18O')].values))
    
    print('Clean FT dD values\n=========')
    q_mbl = 1000; dD_mbl = -70; d18O_mbl = -10
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('1','d18O')].values))
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('2','d18O')].values))
    print(wisper_sigfit(q_mbl, d18O_mbl, sigWISP_df.loc[('3','d18O')].values))
    
    
    
    q, d18O = np.meshgrid(np.linspace(1000,15000,100), np.linspace(-80,-9,150))
    q, dD = np.meshgrid(np.linspace(1000,15000,100), np.linspace(-300,-60,150))
    mod_sigdD = wisper_sigfit(q, dD, sigWISP_df.loc[('1','dD')].values)
    mod_sigd18O = wisper_sigfit(q, d18O, sigWISP_df.loc[('1','d18O')].values)

    fig = plt.figure(figsize=(6.5,2.75))
    ax1 = fig.add_axes([0.1,0.175,0.38,0.6])
    cax1 = fig.add_axes([0.1,0.9,0.38,0.05])
    ax2 = fig.add_axes([0.6,0.175,0.38,0.6])
    cax2 = fig.add_axes([0.6,0.9,0.38,0.05])
    

    p1 = ax1.scatter(q, dD, c=mod_sigdD, cmap='gist_ncar', 
                          vmin=2, vmax= 10)
    plt.colorbar(p1, cax=cax1, orientation='horizontal')

    p2 = ax2.scatter(q, d18O, c=mod_sigd18O, cmap='gist_ncar', 
                     vmin=0.2, vmax=8)
    plt.colorbar(p2, cax=cax2, orientation='horizontal')
            
            
    """
        ## Collect expected values for calibration parameters into a list.
            # Humidity dependence correction params:
        year = 2016
        qdep_caldata = caltable.loc[caltable['year']==year]       

            # Collect into list:
        params = np.append(qdep_caldata[['a18O','b18O']].values[0], [m,k]) 
        
        
    sig_m18O = 0.05185/2
    sig_k18O = 0 
    sig_params_d18O = np.append(qdep_caldata[['sig_a18O','sig_b18O']].values[0], 
                                [sig_m18O,sig_k18O])
    get_uncertainties_with_fit_d18O(sig_params_d18O, sig_inputvars=False, 
                                    ax_cax=(ax_sig18O, cax_sig18O))
    get_uncertainties_with_fit_d18O(sig_params_d18O, sig_inputvars=True)
    
    sig_k18O = 1/2. 
    sig_params_d18O = np.append(qdep_caldata[['sig_a18O','sig_b18O']].values[0], 
                                [sig_m18O,sig_k18O])
    get_uncertainties_with_fit_d18O(sig_params_d18O, sig_inputvars=False, 
                                    ax_cax=(ax_sig18O, cax_sig18O))    
    """
    ##-------------------------------------------------------------------------


    """
    Test 1. 
    Run my monte carlo sampler through only the humidity-dependence 
    correction.
    """          
    # Test 1:
    """
    caldata2017 = caltable.loc[caltable['year']==2018]
    #m_D = 1.0; sig_mD = 0.01
    #kD = 0; sig_kD = 1
    m_D = 1.0564; sig_mD = 0.0564/2
    kD = -5.957469671; sig_kD = 5/2.
    params = np.append(caldata2017[['aD','bD']].values[0], [m_D,kD])
    sig_params = np.append(caldata2017[['sig_aD','sig_bD']].values[0], 
                           [sig_mD,sig_kD])

    dD = np.ones(100)*-70
    q = np.linspace(1000,15000,100)
    results1 = mc_normalsampler(cal_formulas.q_dep_cal, [dD,q], 
                               caldata2017[['aD','bD']].values[0], 
                               caldata2017[['sig_aD','sig_bD']].values[0], 
                               1000, 
                               unpack_params=True)
    
    plt.figure()
    plt.hist(results1[2][:,-1], histtype='step')
    plt.hist(results1[2][:,-10], histtype='step')
    plt.hist(results1[2][:,-30], histtype='step')
    plt.hist(results1[2][:,-60], histtype='step')
    plt.hist(results1[2][:,-100], histtype='step')
    
    
    # Test2
    dD = np.ones(100)*-70
    q = np.linspace(1000,15000,100)
    results2 = mc_normalsampler(pic1_isoratio_fullcal, [dD,q], 
                               params, sig_params, 
                               1000)    

    plt.figure()
    plt.hist(results2[2][:,-1], histtype='step')
    plt.hist(results2[2][:,-10], histtype='step')
    plt.hist(results2[2][:,-30], histtype='step')
    plt.hist(results2[2][:,-60], histtype='step')
    plt.hist(results2[2][:,-100], histtype='step')
    
    
    # Test3
    q, dD = np.meshgrid(np.linspace(1000,15000,100), np.ones(100)*-70)
    results3 = mc_normalsampler(pic1_isoratio_fullcal, [dD,q], 
                               params, sig_params, 
                               1000)
    plt.figure()
    plt.hist(results3[2][:,0,-1], histtype='step', label=q[0,-1])
    plt.hist(results3[2][:,0,-10], histtype='step', label=q[0,-10])
    plt.hist(results3[2][:,0,-30], histtype='step', label=q[0,-30])
    plt.hist(results3[2][:,0,-60], histtype='step', label=q[0,-60])
    plt.hist(results3[2][:,0,-100], histtype='step', label=q[0,-100])
    plt.legend()
    """
    
    # Test4
    """
    q, dD = np.meshgrid(np.linspace(1000,15000,100), np.linspace(-300,-70,100))
    results4 = mc_normalsampler(pic1_isoratio_fullcal, [dD,q], 
                               params, sig_params, 
                               1000)
    plt.figure()
    plt.hist(results4[2][:,-1,-1], histtype='step', label=q[0,-1])
    plt.hist(results4[2][:,-1,-10], histtype='step', label=q[0,-10])
    plt.hist(results4[2][:,-1,-30], histtype='step', label=q[0,-30])
    plt.hist(results4[2][:,-1,-60], histtype='step', label=q[0,-60])
    plt.hist(results4[2][:,-1,-100], histtype='step', label=q[0,-100])
    plt.legend()
    """