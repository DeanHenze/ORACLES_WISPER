# -*- coding: utf-8 -*-
"""
Created on Sun May 16 12:30:38 2021

@author: Dean
"""


# Built in:
import sys

# Third party:
import matplotlib.pyplot as plt # 2.0.2
#from mpl_toolkits.basemap import Basemap # 1.1.0
import matplotlib.colors as colors # 2.0.2
import numpy as np # 1.13.3

# My modules:
if '../..' not in sys.path: sys.path.insert(0,'../..')
from my_python_modules import oracles_data_loader as odl
from my_python_modules import iso_fxns as iso


# Set font for plot labels/annotations/captions:
font = {'family' : 'serif', 'weight' : 'normal'}; plt.rc('font', **font)


"""
Plot potential temperature, humidity, carbon monoxide, dD, and dxs for P3 
data from the passed flight date (str 'yyyymmdd') and vertical profile number 
(int). Input axset is a 5-element list-like of matplotlib axes objects to plot 
each variable on, in the order of the variables listed above.
"""
#______________________________________________________________________________
def plot_1prf(date, profile, axset):
    
# Load datafile and isolate data for the passed profile num. Load data using 
# one of my other modules:
    p3data = odl.p3_with_mixsegs(return_what='data', date=date) # Load data.
    prfdata = p3data.loc[p3data.profile==profile] # Data for specific profile.


# Bin/average the data by 25m altitude bins:
    # z bins as a DataFrame column:
    dz=50
    prfdata['binz'] = np.round(prfdata['height_m']/dz, decimals=0)*dz
    # average by bins with pandas magic:
    prfdata = prfdata.groupby(by='binz', axis=0, sort=True).mean()
    
    
# Plot vertical profiles of potential temperature, humidity, CO, and dD:
    hght_km = prfdata['height_m']/1000 # height in km
    axset[0].plot(prfdata['theta_K'], hght_km, c='k')
    axset[1].plot(prfdata['h2o_gkg'], hght_km, c='k')
    axset[2].plot(prfdata['CO_ppbv'], hght_km, c='k')
    axset[3].plot(prfdata['dD_permil'], hght_km, c='k')
    

# Compute and plot dxs:
    axset[4].plot(prfdata['dD_permil'] - 8*prfdata['d18O_permil'], hght_km, 
                  c='k', linestyle='dashed')
#______________________________________________________________________________


#______________________________________________________________________________
"""
Makes a matplot figure with a set of 5 axes to plot the profiles in the 
plot_1prf fxn above. Axes are positioned as a row of 4 axes with the 5th being 
a twiny axes added to the 4th. Also takes care of the axes spacing, aesthetics, 
labels. Returns the axes as a 5-element list.
"""
#______________________________________________________________________________
def make_fig_with_axes():

# Figure and axes:
    fig = plt.figure(figsize=(6.5,2.25))

    marg_l = 0.1 # Left margin
    marg_b = 0.2 # bottom margin
    w_ax = 0.215; h_ax = 0.6 # axes width and height
    d_interax = 0.005 # horizontal distance between adjacent axes

    ax1 = fig.add_axes([marg_l, marg_b, w_ax, h_ax])
    ax2 = fig.add_axes([marg_l+w_ax+d_interax, marg_b, w_ax, h_ax])
    ax3 = fig.add_axes([marg_l+2*(w_ax+d_interax), marg_b, w_ax, h_ax])
    ax4 = fig.add_axes([marg_l+3*(w_ax+d_interax), marg_b, w_ax, h_ax])
    ax5 = ax4.twiny()
    
    axset = [ax1,ax2,ax3,ax4,ax5]


# Remove some of the axes frame sides
    ax1.spines["right"].set_visible(False)
    for ax in [ax1,ax2,ax3,ax4]: ax.spines["top"].set_visible(False)
    ax4.spines["top"].set_visible(False)
    for ax in [ax2,ax3,ax4,ax5]:
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
    
    ax2.set_yticks([]); ax3.set_yticks([]); ax4.set_yticks([])
    

# tick marks, labels, limits:
    for ax in axset: ax.set_ylim(0,6)
    ax4.set_xlim(-250,-50)
    ax5.set_xticks([0,10,20]); ax5.set_xlim(0,30)
    
    
# Axes labels:
    ax1.set_ylabel('height (km)', fontsize=14)
    ax1.set_xlabel(r'$\theta$ (K)', fontsize=14, labelpad=0)
    ax2.set_xlabel('q (g/kg)', fontsize=14, labelpad=0)
    ax3.set_xlabel(r'c$_{CO}$ (ppbv)', fontsize=14, labelpad=0)
    ax4.set_xlabel(r'$\delta$D ' + u'(\u2030)', fontsize=14, labelpad=0)
    ax5.set_xlabel(r'dxs ' + u'(\u2030)', fontsize=14)
    
    
    return fig, axset
#______________________________________________________________________________
    

"""
Get vertical profiles for 6 cases, 2 from each year. Plot figures and save.
"""
# The 6 cases. Each 2-tuple is the flight date and the vertical profile number:
#cases = [('20160914',2),('20160831',6),('20170815',6),
#         ('20170826',3),('20181003',3),('20181019',2)]
cases = [('20160914',2),('20160831',6),('20170815',5),
         ('20170828',3),('20181003',3),('20181019',2)]

for case in cases:
    fig, axset = make_fig_with_axes()
    plot_1prf(case[0], case[1], axset)
    fig.savefig("../scratch/single-prf_"+case[0]+"_prf"+str(case[1])+".png")


















