# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:12:23 2020

@author: User
"""


import sys
import srtm_dem_tools
from srtm_dem_tools.constructing import SRTM_dem_make, SRTM_dem_make_batch
from srtm_dem_tools.plotting import dem_show
from srtm_dem_tools.water_bodies import water_pixel_masker

import numpy.ma as ma                                                                 # used for DEMS and water masks 
from pathlib import Path

gshhs_dir = Path("./GSHHG_coastlines_2.3.7/GSHHS_shp/")

ed_username = input(f'Please enter your USGS Earthdata username:  ')
ed_password = input(f'Please enter your USGS Earthdata password (NB characters will be visible!   ):  ')

#%% Make and show a SRTM3 DEM

print("\n\n\n###############################\nExample1\n###############################")

dem, lons, lats =  SRTM_dem_make({'west':-4, 'east':-1, 'south':53, 'north':55},
                                  SRTM1_or3 = 'SRTM3', SRTM1_tiles_folder = Path('./SRTM3/'),
                                  water_mask_resolution = 'i', void_fill = False, gshhs_dir = gshhs_dir,                                    # make the dem
                                  ed_username = ed_username, ed_password = ed_password)

dem_show(dem,lons,lats,title = '1: SRTM3 DEM')                                              # plot the DEM




#%% Make and show a SRTM1 DEM
print("\n\n\n###############################\nExample2\n###############################")

dem, lons, lats =  SRTM_dem_make({'west':-3, 'east':-1, 'south':53, 'north':55},
                                  SRTM1_or3 = 'SRTM1', SRTM1_tiles_folder = Path('./SRTM1/'), gshhs_dir = gshhs_dir,
                                  water_mask_resolution = 'i', ed_username = ed_username, ed_password = ed_password)


dem_show(dem,lons,lats, title = '2: SRTM1 DEM')                                                                   # plot the DEM


#%% Or make a DEM, and mask the water in a separate step.   

print("\n\n\n###############################\nExample3\n###############################")

dem2, lons2, lats2 =  SRTM_dem_make({'west':11, 'east':13, 'south':41, 'north':43},
                                    SRTM1_or3 = 'SRTM3', SRTM3_tiles_folder = Path('./SRTM3/'),
                                    water_mask_resolution = None, download = True,                                                      # make the dem
                                    ed_username = ed_username, ed_password = ed_password)


standalone_mask =  water_pixel_masker(dem2, (lons2[0,0], lats2[-1,0]), (lons2[-1,-1], lats2[0, -1]), 'h', gshhs_dir, verbose = True)               # make the water mask

dem_show(ma.array(dem2, mask = standalone_mask), lons2, lats2, title = '3: SRTM3 DEM - standalone mask')                                # plot the DEM



#%% The limits need not be integers, and it will automatically make the water mask in a faster way for small dems:

print("\n\n\n###############################\nExample4\n###############################")
    
dem, lons, lats =  SRTM_dem_make({'west':-5.1, 'east':-4.4, 'south':55.5, 'north':56.1},
                                  SRTM1_or3 = 'SRTM3', SRTM1_tiles_folder = Path('./SRTM3/'),
                                  water_mask_resolution = 'f', void_fill = False, gshhs_dir = gshhs_dir,                                   # make the dem
                                  ed_username = ed_username, ed_password = ed_password)

dem_show(dem, lons, lats, title = '4: SRTM3 DEM')                                              # plot the DEM


#%% Taal Volcano is also an intersting example as there's land a lake, in which there's an island, in which there's a lake, in which there's an island

print("\n\n\n###############################\nExample5\n###############################")

dem, lons, lats =  SRTM_dem_make({'west':120.7, 'east':121.2, 'south':13.8, 'north':14.1},
                                  SRTM1_or3 = 'SRTM3', SRTM1_tiles_folder = Path('./SRTM3/'),
                                  water_mask_resolution = 'f', void_fill = False, gshhs_dir = gshhs_dir,                                   # make the dem
                                  ed_username = ed_username, ed_password = ed_password)

dem_show(dem, lons, lats, title = '5: Taal (a good test for water body masking)')                                              # plot the DEM


#%% Instead of using the edges of the DEM (west, east etc.), a centre and side length (in m) can be supplied.  

print("\n\n\n###############################\nExample6\n###############################")

dem, lons, lats =  SRTM_dem_make({'centre': (-3.396, 37.041), 'side_length':(110e3, 100e3)},
                                  SRTM1_or3 = 'SRTM3', SRTM1_tiles_folder = Path('./SRTM3/'),
                                  water_mask_resolution = 'f', void_fill = False,  gshhs_dir = gshhs_dir,                                  # make the dem
                                  ed_username = ed_username, ed_password = ed_password)

dem_show(dem, lons, lats, title = '5: SRTM3 DEM, centre and side_length style')                                              # plot the DEM


#%% DEMs can also be made in batches by creating a list of dictionaries.  

print("\n\n\n###############################\nExample7\n###############################")

volcano_dems = [{'name' : 'Vulsini',       'centre' : (11.93, 42.6),    'side_length' : (20e3,20e3)},
                {'name' : 'Campi Flegrei', 'centre' : (14.139, 40.827), 'side_length' : (30e3,20e3)},
                {'name' : 'Etna', 'west' : 14.75, 'east' : 15.5, 'south'  : 37.4, 'north' : 38.0}]

volcano_dems2 = SRTM_dem_make_batch(volcano_dems, water_mask_resolution = 'f', ed_username = ed_username, ed_password = ed_password,        # make the DEMS
                                    gshhs_dir = gshhs_dir)

for volc_n, volcano in enumerate(volcano_dems2):                                                                                   # loop through to display them
    dem_show(volcano['dem'], volcano['lons_mg'], volcano['lats_mg'], title = f"{volc_n+5}: DEM # {volc_n}: {volcano['name']}")

