# -*- coding: utf-8 -*-
"""
Created on Sat May  7 10:04:34 2022

@author: Dean

Function: Compute cloud water content from the WISPER CVI.
"""


def cvi_cwc(q_cld, T, P, cvi_enhance):
    """
    Compute CVI-measured cloud water content given cloud water mixing ratio 
    (g/kg), ambient air temperature [K], pressure [Pa], and the CVI 
    enhancement factor. Approximate air density assuming air is completely dry.
    
    Inputs are either scalar or numpy.array-like. For those input which are 
    arrays, they must be same size.

    Returns
    -------
    Cloud water content as a numpy array, units of g/m3.    
    """
    Rd = 287.1 # specific gas constant for dry air, [J/K/kg].
    rho = P/(Rd*T) # density.
    return q_cld*rho/cvi_enhance