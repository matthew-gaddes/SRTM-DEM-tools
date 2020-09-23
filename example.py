# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:12:23 2020

@author: User
"""

import sys
print('Temporary import of matrix_show for debugging.  ')
sys.path.append('/home/matthew/university_work/python_stuff/python_scripts/')
from small_plot_functions import matrix_show


from dem_tools_lib import SRTM_dem_make, SRTM_dem_make_batch, dem_show, water_pixel_masker
import numpy.ma as ma                                                                 # used for DEMS and water masks 




#%% Make and show a SRTM3 DEM

dem, lons, lats =  SRTM_dem_make(-4, -1, 53, 55, SRTM1_or3 = 'SRTM3', SRTM1_tiles_folder = './SRTM3/',
                                  water_mask_resolution = 'i', void_fill = False)                                    # make the dem

dem_show(dem,lons,lats,title = 'SRTM3 DEM')                                              # plot the DEM

#%% Make and show a SRTM1 DEM

ed_username = input(f'Please enter your USGS Earthdata username:  ')
ed_password = input(f'Please enter your USGS Earthdata password (NB characters will be visible!   ):  ')


dem, lons, lats =  SRTM_dem_make(-3, -1, 53, 55, SRTM1_or3 = 'SRTM1', SRTM1_tiles_folder = './SRTM1/',
                                    water_mask_resolution = 'i', ed_username = ed_username, ed_password = ed_password)


dem_show(dem,lons,lats, title = 'SRTM1 DEM')                                                                   # plot the DEM


#%% Or make a DEM, and mask the water in a separate step.   

dem2, lons2, lats2 =  SRTM_dem_make(11, 13, 41, 43, SRTM1_or3 = 'SRTM3', SRTM3_tiles_folder = './SRTM3/',
                                  water_mask_resolution = None, download = True)                                    # make the dem

standalone_mask =  water_pixel_masker(dem2, (lons2[0,0], lats2[0,0]), (lons2[-1,-1], lats2[-1, -1]), 'h', verbose = True)

dem_show(ma.array(dem2, mask = standalone_mask), lons, lats, title = 'SRTM3 DEM - standalone mask')                                  # plot the DEM


#%% The limits need not be integers, and it will automatically make the water mask in a faster way for small dems:


dem, lons, lats =  SRTM_dem_make(-5.1, -4.4, 55.5, 56.1, SRTM1_or3 = 'SRTM3', SRTM1_tiles_folder = './SRTM3/',
                                  water_mask_resolution = 'f', void_fill = False)                                    # make the dem

dem_show(dem, lons, lats, title = 'SRTM3 DEM')                                              # plot the DEM



#%% DEMs can also be made in batches by creating a list of dictionaries.  

volcano_dems = [{'name' : 'Vulsini',       'centre' : (11.93, 42.6),    'side_length' : (20,20)},
                {'name' : 'Campi Flegrei', 'centre' : (14.139, 40.827), 'side_length' : (30,20)},
                {'name' : 'Etna', 'west' : 14.75, 'east' : 15.5, 'south'  : 37.4, 'north' : 38.0}]

volcano_dems2 = SRTM_dem_make_batch(volcano_dems, water_mask_resolution = 'f')                                  # make the DEMS

for volc_n, volcano in enumerate(volcano_dems2):                                                                                   # loop through to display them
    dem_show(volcano['dem'],volcano['lons'], volcano['lats'], title = f"DEM # {volc_n}: {volcano['name']}")

#%%





# Lolobau,-4.92,151.158
# Kuwae,-16.829,168.536
# Krakatau,-6.102,105.423
# Batur,-8.242,115.375
# Banda Api,-4.523,129.881
# Taal,14.002,120.993
# Kikai,30.789,130.308
# Aira,31.593,130.657
# Asosan,32.884,131.104
# Naruko,38.729,140.734
# Towada,40.51,140.88
# Nishinoshima,27.247,140.874
# Ioto,24.751,141.289
# Shikotsu,42.688,141.38