# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:21:47 2020

@author: DeanH808

Script to estimate Picarro accuracy uncertainties for heavy water isotope 
measurements during the ORACLES campaign.

Estimation method:
=================

(1) Uncertainty due to absolute calibration uncertainty:
=================================================================
The absolute calibration of Picarro heavy water measurements scales linearly 
with the true delta value being measured. Based off near surface dD values 
from other studies in the Atlantic, I am assuming that the true value for the 
ORACLES marine boundary layer (MBL) should lie in the range [-60 permil, 
-80 permil]. Therefore it is fair to say that the uncertainty in 
our marine boundary layer (MBL) measurements is around 10 permil, by inspection 
of our data. To get uncertainties at lower dD values, I assume that the 
uncertainty scales as the deviation in the Picarro absolute calibration slope 
from 1. Call this deviation 'c'. 

For our 5Hz instrument before the ORACLES campaign, the abs cal slope was 
1.0564 and therefore c =  1 - 1.0564 = -0.0564. From this together with the 
assumed uncertainty for the MBL value, I am assuming that our uncertainty at 
dD = -70 permil is 10 permil, and that the uncertainty increases with 
decreasing dD by a scale factor of c=-0.0564. In my mind, this is in effect 
accounting for the possibility that the value of c during ORACLES is anywhere 
in the range [0,2*-0.0564] - i.e. I want to account for a possible doubling of 
c. 
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.odr import Model, RealData, ODR


# (1) Uncertainty due to absolute calibration uncertainty:
#------------------------------------------------------------------------------
# Range of dD values we're interested in getting uncertainties for:
valsdD = np.arange(-30,-250,-10)

# Absolute calibration slopes and offsets for Mako, given to me by David. 
# I think they are from calibrations done prior to ORACLES. Slopes for dD and 
# d18O are placed here in case d18O uncertainties are needed in the future:
m_D = 1.0564; m_18 = 1.05185
b_D = 6; b_18 = 1.04

# Deviation of above slopes from 1:
c_D = 1 - m_D; c_18 = 1 - m_18

# Assumed uncertainty in accuracy (permil) for dD marine boundary layer (MBL) 
# value, say at -70 permil: 
sigdD_mbl = 10

# Derive an offset b for a y=mx+b relationship based off the above MBL 
# uncertainty at -70 permil:
bprime_D = sigdD_mbl - (c_D)*(-70)

# Deviation of Mako's measurements from actual values, based off of cal slope 
# and derived b offset:
devsdD = abs((valsdD*(c_D) + bprime_D))

# Plot deviations as estimated accuracies as a function of dD; also plot the 
# abs cal slope:
plt.figure()
plt.subplot(1,2,1)
plt.plot(valsdD, valsdD*m_D + b_D, 'k-', linewidth=5) # Abs cal line.
plt.plot(valsdD, valsdD, 'k--', linewidth=2) # 1-1 line.
plt.fill_between(valsdD, valsdD + b_D, valsdD*m_D + b_D, color='grey')
plt.fill_between(valsdD, valsdD*m_D + b_D, valsdD*(m_D-c_D) + b_D, color='grey')
plt.title('Absolute calibration Pic1. \nSlope = '+str(m_D), fontsize=14)
plt.xlabel('Picarro dD', fontsize=14)
plt.ylabel('True dD', fontsize=14)
plt.subplot(1,2,2)
plt.plot(valsdD, devsdD, 'b-', linewidth=5)
plt.xlabel('Picarro measured dD (permil)', fontsize=14)
plt.ylabel('Estimated uncertainty in accuracy (permil)', fontsize=14)
plt.title('Estimated uncertainty in accuracy (permil)'
           '\n slope= 1 - '+str(m_D))
#------------------------------------------------------------------------------




    













