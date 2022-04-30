# -*- coding: utf-8 -*-
"""
@author: Dean


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
        
    return_stdev (default='no'): 
        Options of 'no', 'yes', 'yes_addweights'. If 'yes', return weighted 
        standard deviations at each gridpoint, computed using the same general 
        oversampling technique as for the means. If 'yes_addweights', compute 
        it with the additional weights from the input 'w_add'. If 'no', returns 
        NANs. In all cases, results are returned as an array with the same 
        shape as the oversampled map.
        
Returns: Dictionary with keys:
    'mean': oversampled map as a 2D numpy array, 
    'stdev': standard deviation at each gridpoint, same shape as 'mean' (optional),
    'bandwidths': x and y bandwiths as a 2-tuple. 
    'n_obs_weighted': Weighted number of observations at each gridoint, same 
        shape as 'mean'.
"""



import numpy as np # 1.13.3
from scipy.stats import iqr # 1.1.0



def oversampler(mapdata, xcoords, ycoords, gridx, gridy, 
                w_add=None, ffact=1., return_stdev='no'):

    ## Optional additional weighting factors:
    if w_add is None: # No additional weights.
        w_add = np.ones(len(mapdata)) # Array of 1's won't change result
    else:
        w_add=w_add


    ## Make sure inputs are numpy arrays and remove any NAN's:
    mapdata = np.asarray(mapdata)
    xcoords = np.asarray(xcoords); ycoords = np.asarray(ycoords) 
    
    all_finite = np.isfinite(mapdata*xcoords*ycoords*w_add)
    mapdata = mapdata[all_finite]
    xcoords = xcoords[all_finite]; ycoords = ycoords[all_finite]
    w_add = w_add[all_finite]
    
    
    ## Estimate x,y bandwidth with Silverman's rule of thumb:
        # Interquartile ranges for x,y:
    rx = iqr(xcoords, nan_policy='omit'); ry = iqr(ycoords, nan_policy='omit')
        # bandwidths:
    hx = ffact*0.9*(rx/1.34)*len(xcoords)**(-1/5)
    hy = ffact*0.9*(ry/1.34)*len(ycoords)**(-1/5)  


    ## Get weighted means, number of weighted obs, and optionally weighted 
    ## standard deviations at each gridpoint:
        
        # Weighting function as gaussian of distance from grid x,y point. xdist 
        # and ydist, are the abs values of the distances in the x,y directions:
    def weights_gauss(xdist, ydist, sigx, sigy):
        A_norm = 1/(2*np.pi*sigx*sigy) # Normalization factor
        return A_norm*np.exp(-0.5*(xdist**2/sigx**2 + ydist**2/sigy**2))    
        
    m = len(gridy); k = len(gridx)
    map_ovrsmpled = np.ones([m, k])*np.nan
    std_map = np.ones([m, k])*np.nan
    n_obs_weighted = np.ones([m, k])*np.nan
    
    for i in range(m):
        for j in range(k):
        
        # Gaussian-distance weights:
            xdist = abs(xcoords - gridx[j])
            ydist = abs(ycoords - gridy[i])
            w_gauss = weights_gauss(xdist, ydist, hx, hy)
        
        # Weighted total number of obs:
            n_obs_weighted[i,j] = np.sum(w_gauss)
        
        # Weighted mean:
            S = np.sum(w_gauss*w_add)
            map_ovrsmpled[i,j] = np.sum(mapdata*w_gauss*w_add)/S
            
        # Optional weighted standard deviation. It's a biased estimate but 
        # won't matter for a large number of data points:
            if return_stdev == 'no':
                continue
            elif return_stdev == 'yes': 
                S = np.sum(w_gauss) # Recompute sum of weights w/o additional w's.
                #wmean_ij = map_ovrsmpled[i,j]
                wmean_ij = np.sum(mapdata*w_gauss)/S # Recompute mean w/o additional w's.
                std_map[i,j] = (np.sum(w_gauss*(wmean_ij - mapdata)**2)/S)**0.5

            elif return_stdev == 'yes_addweights': 
                wmean_ij = map_ovrsmpled[i,j]
                std_map[i,j] = (np.sum(w_gauss*w_add*(wmean_ij - mapdata)**2)/S)**0.5                
         

    ## Collect results in dict and return:
    results = {'mean':map_ovrsmpled, 'stdev':std_map, 'bandwidths':(hx, hy), 
               'n_obs_weighted':n_obs_weighted}
    return results

