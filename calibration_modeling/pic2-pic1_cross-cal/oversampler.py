# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 14:26:30 2019

@author: Dean

Function to oversample a 2D map to get a smoothed map.
"""


import numpy as np # 1.13.3
from scipy.stats import iqr # 1.1.0


def oversampler(mapdata, xcoords, ycoords, gridx, gridy, w_add=None, ffact=1., 
                return_bandwidths=False):
    """
    Oversample a 2d map to get a smoothed map. Smoothed map is computed at set 
    grid locations by taking the weighted mean of data points near that location. 
    Weights are the gaussian of the distance of the data points from the grid loc.
    
    Estimate bandwidth for guassian weighting function based on Silverman's rule of 
    thumb for kernel density estimations (seems reasonable).
    
    Inputs:
        mapdata (np array, length n): 1d array of map data values.
        
        xcoords, ycoords (np arrays, each size n): 1d arrays of x,y-coords 
            for each data point in mapdata. 
            
        gridx, gridy (1d np arrays): x,y-locations on the grid for the 
            smoothed map. Do not need to be the same length.
            
        w_add (1d np array, length n; default=None): Optional weighting factors 
            for each data point in mapdata, which are applied in addition to the 
            gaussian distance weighting. 
        
        ffact (float, default=1.): Fudge factor added to the bandwidth 
            calc, i.e. to make it smaller if desired.
            
        return_bandwidths (default=False): If True, return x and y bandwiths as 
            a 2-tuple.
            
    Returns: oversampled map as a 2D numpy array, and optionally the x and y 
        bandwiths as a 2-tuple. 
    """
    # Optional additional weighting factors:
    if w_add is None: # No additional weights.
        w_add = np.ones(len(mapdata)) # Array of 1's won't change result
    else:
        w_add=w_add

    # Make sure inputs are numpy arrays and remove any NAN's:
    mapdata = np.asarray(mapdata)
    xcoords = np.asarray(xcoords); ycoords = np.asarray(ycoords) 
    
    all_finite = np.isfinite(mapdata*xcoords*ycoords*w_add)
    mapdata = mapdata[all_finite]
    xcoords = xcoords[all_finite]; ycoords = ycoords[all_finite]
    w_add = w_add[all_finite]
    
    # Estimate x,y bandwidth with Silverman's rule of thumb:
        # Interquartile ranges for x,y:
    rx = iqr(xcoords, nan_policy='omit'); ry = iqr(ycoords, nan_policy='omit')
        # bandwidths:
    hx = ffact*0.9*(rx/1.34)*len(xcoords)**(-1/5)
    hy = ffact*0.9*(ry/1.34)*len(ycoords)**(-1/5)  

    # Weighting function as gaussian of distance from grid x,y point. xdist 
    # and ydist, are the abs values of the distances in the x,y directions:
    def weights_gauss(xdist, ydist, sigx, sigy):
        A_norm = 1/(2*np.pi*sigx*sigy) # Normalization factor
        return A_norm*np.exp(-0.5*(xdist**2/sigx**2 + ydist**2/sigy**2))
        
    # Get weighted means at each gridpoint as the oversampled map:
    m = len(gridy); k = len(gridx)
    map_ovrsmpled = np.zeros([m, k])
    
    for i in range(m):
        for j in range(k):
        
        # x,y distances of each data point from current gridpoint: 
            xdist = abs(xcoords - gridx[j]); ydist = abs(ycoords - gridy[i])
            
        # Gaussian-distance weights:
            w_gauss = weights_gauss(xdist, ydist, hx, hy)
        
        # Weighted mean:
            S = np.sum(w_gauss*w_add)
            map_ovrsmpled[i,j] = np.sum(mapdata*w_gauss*w_add)/S
            

    if return_bandwidths: return map_ovrsmpled, (hx,hy)
    return map_ovrsmpled




    





