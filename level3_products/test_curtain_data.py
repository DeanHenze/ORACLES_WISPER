

import numpy as np
import xarray as xr
import netCDF4 as nc # 1.3.1
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as linsegcmap



##-----------------------------------------------------------------------------
## Load curtain data for each sampling year:
curtains2016 = xr.load_dataset("wisper_oracles_curtains_2016.nc")
curtains2017 = xr.load_dataset("wisper_oracles_curtains_2017.nc")
curtains2018 = xr.load_dataset("wisper_oracles_curtains_2018.nc")
##-----------------------------------------------------------------------------



##-----------------------------------------------------------------------------
## Contour levels and color mapping for mean curtain maps

def get_clevels_cmapping():
    """
    Returns contour levels and colormappings to use for the H2O, dD, and d18O 
    curtain plots.
    """
    
    ## H2O contour levels
    levs_h2o = np.arange(0,17,1.)
    
    ## Iso ratio contour levels:
    dD_min = -350; dD_max = -64 # dD_min=-260
    d18O_min = -45; d18O_max = -10.5 # d18O_min=-32
        
    a_dD = 1 # Desired spacing in permil between the first 2 contours.
    n = np.log2((dD_max-dD_min)/a_dD) # Power of 2 that spans the range (dD_max-dD_min)/a_dD.
    levspacing_dD = a_dD*(2**np.linspace(1,n,15))
    levs_dD = np.flip(np.append([dD_max], dD_max-levspacing_dD))
    levs_dD = np.round(levs_dD)
    
    a_d18O = 0.1
    n = np.log2((d18O_max-d18O_min)/a_d18O)
    levspacing_d18O = a_d18O*(2**np.linspace(1,n,15))
    levs_d18O = np.flip(np.append([d18O_max], d18O_max-levspacing_d18O))
    levs_d18O = np.round(levs_d18O, decimals=1)
    
    levs_dxs = np.arange(-7,27,2)
        
    ## Color mapping:
        # H2O is straightforward:
    cmap_name = 'gist_ncar' 
    cmap_h2o = plt.get_cmap(cmap_name)
  
        # Iso ratios are less straightforward:
    def cmapping_to_even(levs, cmap_name):
        """
        Map uneven level spacings to even colormap spacings.
        levs: contour levels (array) in increasing order.
        cmap (str): matplotlib colormap to use.
        """
        cmap = plt.get_cmap(cmap_name)
        normedlevs = (levs-levs[0])/(levs[-1]-levs[0]) # Norm to 0-1 scale.
        colors = cmap(np.linspace(0,1,len(levs))) # Cmap colors at even spacings.
        # Map:
        return linsegcmap.from_list(cmap_name, list(zip(normedlevs,colors)))
    
    cmap_dD = cmapping_to_even(levs_dD, cmap_name)
    cmap_d18O = cmapping_to_even(levs_d18O, cmap_name)
        
    
    return dict(q=(levs_h2o,cmap_h2o), dD=(levs_dD, cmap_dD), 
                d18O=(levs_d18O, cmap_d18O), 
                dxs=(levs_dxs, cmap_h2o))


clevs = get_clevels_cmapping()
##-----------------------------------------------------------------------------



##-----------------------------------------------------------------------------
## Contour levels for standard deviations
clevs_std_2016 = dict(q=[1, 2, 3], dD=[10, 40, 70, 100], dxs=[5, 10, 15, 20])
clevs_std_20172018 = dict(q=[1, 2, 3], dD=[5, 20, 35, 50], dxs=[2, 6, 10, 14])
##-----------------------------------------------------------------------------



## Function to plot curtains with standard deviations
##-----------------------------------------------------------------------------
def plotcurtain(data, varkey, clevs, clevs_std, fig, ax, colorbar=True):
    
    # Filled contours for curtains (mean values):
    csf = ax.contourf(
        data['lat'], data['alt']/1000, data[varkey], 
        levels=clevs[varkey][0], cmap=clevs[varkey][1]
        )
    if colorbar: fig.colorbar(csf, orientation='vertical', ax=ax)
    
    # Black lines for standard deviations:
    cs_sig = ax.contour(
        data['lat'], data['alt']/1000, data['sig_'+varkey], 
        colors='black', levels=clevs_std[varkey]
        )
    ax.clabel(cs_sig, cs_sig.levels, inline=True, fmt='%i', fontsize=10)
##-----------------------------------------------------------------------------



##-----------------------------------------------------------------------------
## 2016 plots
fig, axes = plt.subplots(3, 3, figsize=(14, 8), tight_layout=True)

varkeys = ['q', 'dD', 'dxs']
pltlabel_keys = dict(q='h2o (g/kg)', dD=r'$\delta$D'+u' (\u2030)', dxs='dxs'+u' (\u2030)')
for k, ax in zip(varkeys, axes[:,0]):
    plotcurtain(curtains2016, k, clevs, clevs_std_2016, fig, ax)
    ax.text(
        1.03, 0.5, pltlabel_keys[k], ha='center', va='center', 
        rotation=90, fontsize=12, transform=ax.transAxes
        )
    
axes[2,0].set_xlabel(r'Latitude  ($^\circ$N)', fontsize=12)
axes[1,0].set_ylabel('Altitude (km)', fontsize=12)
axes[0,0].set_title('2016 IOP', fontsize=12)

for ax in axes[:,0]: ax.set_xlim(-24, -10)
##-----------------------------------------------------------------------------



##-----------------------------------------------------------------------------
## 2017 plots
for k, ax in zip(varkeys, axes[:,1]):
    plotcurtain(curtains2017, k, clevs, clevs_std_20172018, fig, ax)
    ax.text(
        1.03, 0.5, pltlabel_keys[k], ha='center', va='center', 
        rotation=90, fontsize=12, transform=ax.transAxes
        )
    
axes[2,1].set_xlabel(r'Latitude  ($^\circ$N)', fontsize=12)
axes[1,1].set_ylabel('Altitude (km)', fontsize=12)
axes[0,1].set_title('2017 IOP', fontsize=12)

for ax in axes[:,1]: ax.set_xlim(-14, 1)
##-----------------------------------------------------------------------------



##-----------------------------------------------------------------------------
## 2018 Results
for k, ax in zip(varkeys, axes[:,2]):
    plotcurtain(curtains2018, k, clevs, clevs_std_20172018, fig, ax)
    ax.text(
        1.03, 0.5, pltlabel_keys[k], ha='center', va='center', 
        rotation=90, fontsize=12, transform=ax.transAxes
        )
    
axes[2,2].set_xlabel(r'Latitude  ($^\circ$N)', fontsize=12)
axes[1,2].set_ylabel('Altitude (km)', fontsize=12)
axes[0,2].set_title('2018 IOP', fontsize=12)

for ax in axes[:,2]: ax.set_xlim(-14, 1)
##-----------------------------------------------------------------------------



fig.savefig('curtains_verification_test.png')