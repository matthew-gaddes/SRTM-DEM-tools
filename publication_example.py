#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 11:37:37 2022

@author: matthew
"""
import numpy as np
import numpy.ma as ma                                                                 # used for DEMS and water masks 
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes            
import cartopy
import cartopy.crs as ccrs

import srtm_dem_tools
from srtm_dem_tools.constructing import SRTM_dem_make, SRTM_dem_make_batch
from srtm_dem_tools.plotting import dem_show
from srtm_dem_tools.water_bodies import load_GSHHS_coastlines, water_pixel_masker
from srtm_dem_tools.aux import truncate_colormap, open_crop_gebco_data



#%% Things to set (default of SRTM1 can be changed in "Make a DEM")


crop_area = {'west'  : 14.0,                                    # in degrees
             'east'  : 14.75,
             "south" : 40.60,
             "north" : 41.10}

cbar_position = {'x' : 14.01,                           # longitude of bottom left of  cbar
                 'y' : 40.95,                           # latitude of bottom left of cbar
                 'width' : 0.01,                       # width, in degrees
                 'height' : 0.1}                         # height, in degrees

water_mask_resolution = 'f'                             # Same naming as GSHSS (shoreline) data:
                                                        #      f     full resolution: Original (full) data resolution.
                                                        #      h     high resolution: About 80 % reduction in size and quality.
                                                        #      i     intermediate resolution: Another ~80 % reduction.
                                                        #      l     low resolution: Another ~80 % reduction.
                                                        #      c     crude resolution: Another ~80 % reduction.
                                                        
figsize=(7.84, 6)                                      #, 7.84 standard width for JGR:SE figsize in inches

gshhs_dir = Path("./GSHHG_coastlines_2.3.7/GSHHS_shp/")
gebco_file = Path("/home/matthew/university_work/data/GEBCO_bathymetry/gebco_2021/GEBCO_2021_sub_ice_topo.nc")



#%% Get login details.  Not necessary (e.g. just enter 'a' and 'b') to run the example as the tiles are provided here.  

ed_username = input(f'Please enter your USGS Earthdata username:  ')
ed_password = input(f'Please enter your USGS Earthdata password (NB characters will be visible!   ):  ')

#%% Make a DEM (SRTM 1)

dem, lons, lats =  SRTM_dem_make(crop_area, SRTM1_or3 = 'SRTM1', SRTM1_tiles_folder = Path('./SRTM1/'), gshhs_dir = gshhs_dir,
                                  water_mask_resolution = water_mask_resolution, ed_username = ed_username, ed_password = ed_password)


#dem_show(dem,lons,lats, title = '2: SRTM1 DEM')                                                                   # plot the DEM


#%% Make a bathymetery DEM (GEBCO data)

bath_data, bath_lons, bath_lats, = open_crop_gebco_data(gebco_file, crop_area, padding = 0.1)                           # lons and lats are meshgrids (i.e. rank 2)

# f, ax = plt.subplots(1)
# ax.imshow(bath_data_crop)
# from srtm_dem_tools.water_bodies import water_pixel_masker

import srtm_dem_tools
from srtm_dem_tools.water_bodies import water_pixel_masker
bath_water_mask = water_pixel_masker(bath_data, (bath_lons[-1,0], bath_lats[-1,0]), (bath_lons[0,-1], bath_lats[0,-1]), 
                                      coast_resol = 'f', gshhs_dir = gshhs_dir, verbose = True, pixs_per_deg = 240,
                                      debug_mode = False)
                       
bath_data_ma = ma.array(bath_data, mask = np.invert(bath_water_mask))


#%% load vector shoreline data

_, intersect_coastline_polygons, intersect_cropped_coastline_polygons = load_GSHHS_coastlines(gshhs_dir, resolution = water_mask_resolution, bbox = crop_area, level = '1')

#%% Make the figure

# Setup for parts of figure
fig = plt.figure(figsize = figsize)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([crop_area['west'], crop_area['east'], crop_area['south'], crop_area['north']])
cmap_terrain = plt.get_cmap('terrain')
cmap_land = truncate_colormap(cmap_terrain, 0.2, 1)
cmap_water = truncate_colormap(cmap_terrain, minval=0.0, maxval=0.17, n=100)

# bathymetery 
### needs masking for vmin v max
ax.imshow(bath_data, extent=(np.min(bath_lons), np.max(bath_lons), np.min(bath_lats), np.max(bath_lats)),                       # plotting the non masked data ensures we don't have any white bits visible due to low resolution of this data
          cmap = cmap_water, vmin = ma.min(bath_data_ma), vmax = ma.max(bath_data_ma), transform=ccrs.PlateCarree())            # but getting v min nad vmas from the masked data ensures the colorbscale fits the data we see.      

# DEM part
#ax.imshow(dem, extent=(np.min(lons), np.max(lons), np.min(lats), np.max(lats)), transform=ccrs.PlateCarree(), cmap = cmap_terrain)

ls = LightSource(azdeg=315, altdeg=45)              # some illumination settings
dem_rgba = ls.shade(dem, cmap=cmap_land, vert_exag=1, blend_mode='hsv')               # 4 channel RGBA data (4th channel is opacity)
ax.imshow(dem_rgba, extent=(np.min(lons), np.max(lons), np.min(lats), np.max(lats)), transform=ccrs.PlateCarree())


# coastlines
for intersect_coastline_polygon in intersect_coastline_polygons:
    ax.plot(intersect_coastline_polygon[:,0], intersect_coastline_polygon[:,1], c = 'k', linewidth = 1.8)

# labels/grid lines etc.  
gl = ax.gridlines(draw_labels=True, linestyle = '--', color = 'k', alpha = 0.5) #, dms=True, x_inline=False, y_inline=False)
gl.bottom_labels = False
gl.right_labels = False

#colorbar

ifg_cb_axes = ax.inset_axes([cbar_position['x'], cbar_position['y'],
                             cbar_position['width'], cbar_position['height']], transform=ax.transData)      # [xpos, ypos, xwidth, ywidth], note the x_pos_scaler that reduces the size of the inset axes to make sure it remains in tehe right place

norm = mpl.colors.Normalize(vmin = (1/1000)*np.min(dem), vmax = (1/1000)*np.max(dem))
cb1 = mpl.colorbar.ColorbarBase(ifg_cb_axes, cmap = cmap_land, norm = norm,  orientation = 'vertical')
tick_locator = ticker.MaxNLocator(nbins=4)
cb1.locator = tick_locator
cb1.update_ticks()
cb1.set_label('Height (Km)', ha = 'center', fontsize = 8)