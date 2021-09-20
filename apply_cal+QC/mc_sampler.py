# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 15:57:12 2021

@author: Dean

Contains a monte carlo function to generate a spread in predictions for a model.
    
Parameter space (dim N) is explored by sampling from N normal distributions 
centered on the respective parameter means and with standard deviations equal 
to the respective parameter standard errors. Random-normal errors may also be 
added to the input variables on each iteration (e.g. measurement errors).
"""


# Third party:
import numpy as np


def mc_normalsampler(model, inputvars, param_means, sig_params, 
                     iterations, unpack_params=False, sig_inputvars=None):
    """
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