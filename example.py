# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:12:23 2020

@author: User
"""



from dem_tools_lib import SRTM_dem_make, dem_show, water_pixel_masker
import numpy.ma as ma                                                                 # used for DEMS and water masks 


#%% Make and show a SRTM3 DEM

dem, lons, lats =  SRTM_dem_make(-4, -1, 53, 55, SRTM1_or3 = 'SRTM3', SRTM1_tiles_folder = './SRTM3/',
                                  water_mask_resolution = 'i', void_fill = False)                                    # make the dem

dem_show(dem,lons,lats,srtm = 3, units_deg = True, title = 'SRTM3 DEM')                                              # plot the DEM

#%% Make and show a SRTM1 DEM

ed_username = input(f'Please enter your USGS Earthdata username:  ')
ed_password = input(f'Please enter your USGS Earthdata password (NB characters will be visible!   ):  ')


dem, lons, lats =  SRTM_dem_make(-3, -1, 53, 55, SRTM1_or3 = 'SRTM1', SRTM1_tiles_folder = './SRTM1/',
                                    water_mask_resolution = 'i', ed_username = ed_username, ed_password = ed_password)


dem_show(dem,lons,lats,srtm = 1, units_deg = True, title = 'SRTM1 DEM')                                                                   # plot the DEM


#%% Or make a DEM, and mask the water in a separate step.   

dem, lons, lats =  SRTM_dem_make(11, 13, 41, 43, SRTM1_or3 = 'SRTM3', SRTM3_tiles_folder = './SRTM3/',
                                  water_mask_resolution = None, download = True)                                    # make the dem

standalone_mask =  water_pixel_masker(dem, lons, lats, 'i', verbose = True)

dem_show(ma.array(dem, mask = standalone_mask),lons,lats,srtm = 3, units_deg = True, title = 'SRTM3 DEM - standalone mask')                                  # plot the DEM

#%%
