# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 15:08:11 2021

@author: Dean

Collection of calibration formulas for both Pic1 and Pic2: 
    - humidity dependence correction
    - absolute calibration
    - Pic2 -> Pic1 cross calibration
"""


import numpy as np


"""
Humidity dependence correction for either Pic1 or Pic2. Takes in values for 
uncorrected isotope ratios (delta-notation) and humidities, and returns 
corrected values.

Inputs:
    deltavals, qvals (each np.array-like, same shape): uncorrected Picarro 
        isotope ratio (delta-notation) and humidity.
    a,b (floats): fit parameters. They differ for Pic1 and Pic2.  
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
Formula for absolute calibration of humidity or iso ratios for any Picarro. 
Just a line. x is either humidity, dD, or d18O. Output is corrected humidity, 
dD, d18O resp.:
"""
def abscal_line(x, m, k):
    return m*x + k


"""
Returns precision in dD as a function of humidity q (np.array-like).
"""
def dD_precision(q):
    return (6*10**6)*(np.log(q)**-6.69)


"""
Returns precision in d18O as a function of humidity q (np.array-like).
"""
def d18O_precision(q):
    return 15425*(np.log(q)**-4.72)