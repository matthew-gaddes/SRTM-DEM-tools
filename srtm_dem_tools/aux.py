#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:16:29 2022

@author: matthew
"""
#%%
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib.colors as colors
    import numpy as np
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap 


#%%

def open_crop_gebco_data(gebco_file, crop_area, padding = 0):
    """  Open the Gebco DEM (which includes bathymtery), and crop to an area based on lon and lat.  
    Also an option to pad it slightly (units are degrees).  
    """
    import numpy as np
    from netCDF4 import Dataset
    
    # load the bathymetry data
    print("Starting to open the gebco .nc which can be slow (as it's ~8GB)...", end = '')
    f = Dataset(gebco_file, 'a')
    bath_data = np.asarray(f.variables['elevation'][:,:])
    bath_lons = np.asarray(f.variables['lon'][:])
    bath_lats = np.asarray(f.variables['lat'][:])
    print("Done.")
    
    
    lon_args = np.ravel(np.argwhere(np.logical_and(crop_area['west']-padding < bath_lons, crop_area['east']+padding > bath_lons)))
    lat_args = np.ravel(np.argwhere(np.logical_and(crop_area['south']-padding < bath_lats, crop_area['north']+padding > bath_lats)))
    lon_args[-1]
    
    
    bath_data_crop = np.flipud(bath_data[lat_args[0] : lat_args[-1], lon_args[0] : lon_args[-1]])
    bath_lons_crop = bath_lons[lon_args[0] : lon_args[-1]]
    bath_lats_crop = bath_lats[lat_args[0] : lat_args[-1]]
    
    bath_lons_mg, bath_lats_mg = np.meshgrid(bath_lons_crop, bath_lats_crop)
    bath_lats_mg = np.flipud(bath_lats_mg)
    
    return bath_data_crop, bath_lons_mg, bath_lats_mg