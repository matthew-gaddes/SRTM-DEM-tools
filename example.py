# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:12:23 2020

@author: User
"""


import numpy.ma as ma

from lib import SRTM_dem_make, dem_show, water_pixel_masker


#%% Make and show an SRTM3 DEM

dem, lons, lats =  SRTM_dem_make(-4, -1, 53, 55, SRTM1_or3 = 'SRTM3', SRTM1_tiles_folder = './SRTM3/')                  # make the dem
dem_mask = water_pixel_masker(dem, lons, lats, 'i', verbose = True)                                                    # make a mask of water pixels
dem_ma = ma.array(dem, mask = dem_mask)                                                                                 # combine the dem and mask
dem_show(dem_ma,lons,lats,srtm = 3, units_deg = True)                                                                   # plot the DEM

#%% Make and show an SRTM1 DEM


ed_username = input(f'Please enter your USGS Earthdata username:  ')
ed_password = input(f'Please enter your USGS Earthdata password (NB characters will be visible!   ):  ')


dem, lons, lats =  SRTM_dem_make(-3, -1, 53, 55, SRTM1_or3 = 'SRTM1', SRTM1_tiles_folder = './SRTM1/',
                                 ed_username = ed_username, ed_password = ed_password)

dem_mask = water_pixel_masker(dem, lons, lats, 'i', verbose = True)
dem_ma = ma.array(dem, mask = dem_mask)
dem_show(dem_ma,lons,lats,srtm = 1, units_deg = True)                                                                   # plot the DEM

