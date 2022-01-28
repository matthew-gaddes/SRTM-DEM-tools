#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 09:55:56 2022

@author: matthew
"""



def water_pixel_masker(dem, lon_lat_ll, lon_lat_ur, coast_resol, gshhs_dir, verbose = False, pixs_per_deg = 1201,
                       debug_mode = False):
    """
       A function to creat a mask of pixels over water. This can be very slow for big DEMs (best called on each tile,
       as per SRTM_dem_make).  
       Note the too avoid problems with edges, the fucntion works by expanding many regions to be one tile larger in all directions.  
       Note that basemap is rather an old package, and can be problematic on Windows machines.  
       
    Inputs:
        dem | rank 2 array | gridded data (e.g. a dem)
        lon_lat_ll | tuple| lon lat of lower left corner of dem.  Must be floats! 
        lon_lat_ur | tuple| lon lat of upper right corner of dem
        coast_resol | str | resolution of vector coastlines: c l i h f 
        gshhs_dir | Path | Path to directory of the GSHHS shorelines.  Available from: http://www.soest.hawaii.edu/pwessel/gshhg/
        pixs_per_deg | int | number of pixels in 1 degree.  1201 for 3 arc second SRTM3 data.  
        debug_mode | Boolean | If true, a simple figure of the data and the coastlines is produced.  
       
    Output:
        result_ary | rank 2 array | array to be used as pixel mask
    
    Hitory:
        2017/03/01 | MEG  | adapterd from dem_show_oceans_detailed_2
        2020/05/18 | MEG | Update to ensure that pixels at the edge that fall on the edge of a land polygon are included
                            (by expanding the plot by 1 deg in each direction)
        2020/06/02 | MEG | Fix bug when masking large DEMS (ie not just a single tile)
        2020/07/29 | MEG | Fix bug for masking of lakes.  
        2020/08/05 | MEG | Reduce extra area added around dem to speed up masking.  
        2020/08/13 | MEG | Fix bug in water masking for tiles with no water bodies.  
    """
         
    from matplotlib.path import Path
    import numpy as np
    import matplotlib.pyplot as plt
    import platform
    import os

    edge_fraction = 0.1
    
    if verbose:
        print('Creating a mask of pixels that lie in water (sea and lakes)... ', end = '')
       
    # 0: check inputs
    if lon_lat_ll[0] > lon_lat_ur[0]:
        raise Exception(f"The western (left) edge ({lon_lat_ll[0]}) appears to be east of the "
                        f"eastern (right) edge ({lon_lat_ur[0]}).  Exiting.  ")
        
    if lon_lat_ll[1] > lon_lat_ur[1]:
        raise Exception(f"The southern (lower) edge ({lon_lat_ll[1]}) appears to be north of the "
                        f"northern (upper) edge ({lon_lat_ur[1]}).  Exiting.  ")
            
    # 1: deal with sizes of various things and make a grid of points
    ny = dem.shape[0]                                                                                           # number of pixels vetically
    nx = dem.shape[1]                                                                                           # and horizontally            
    x = np.linspace(lon_lat_ll[0], lon_lat_ur[0], nx)
    y = np.linspace(lon_lat_ll[1], lon_lat_ur[1], ny)
    xx, yy = np.meshgrid(x,y)
    locations = np.hstack((np.ravel(xx)[:,np.newaxis], np.ravel(yy)[:,np.newaxis]))                             # x and y (lon and lat) for all pixels in the dem, size (n_pixels x 2)
    
    
    bbox = {"west" : lon_lat_ll[0] - edge_fraction, "east" : lon_lat_ur[0] + edge_fraction,
            "south": lon_lat_ll[1] - edge_fraction, "north" : lon_lat_ur[1] + edge_fraction}

    #plt.switch_backend('Qt5Agg')                                  # debugging, plot the coastlines     
    l1_mask = np.zeros(len(locations), dtype=bool)                                                                               # initialise as false, will be True (1) where land
    l2_mask = np.zeros(len(locations), dtype=bool)                                                                               # initialise as false, will be True (1) where land
    l3_mask = np.zeros(len(locations), dtype=bool)                                                                               # initialise as false, will be True (1) where land
    l4_mask = np.zeros(len(locations), dtype=bool)                                                                               # initialise as false, will be True (1) where land
    for product in ['L1','L2','L3','L4']:                                                   # the GSHHS products
        _, _, coastlines = load_GSHHS_coastlines(gshhs_dir, resolution = coast_resol, bbox = bbox, level = product[-1])            # open the GSHHS product (list of numpy arrays of lon lats of water/land boundary)

        if debug_mode:
            title = f'{product} coastlines'
            f1, ax = plt.subplots(1,1)                                                                 # debug plot
            f1.suptitle(title)
            f1.canvas.manager.set_window_title(title)
            for p in coastlines:
                ax.plot(p[:,0], p[:,1])
        
        for n, coastline in enumerate(coastlines):
            landmass = Path(coastline)  
            if product == 'L1':
                l1_mask += np.array(landmass.contains_points(locations))                           # True if pixel is in land
            elif product == 'L2':
                l2_mask += np.array(landmass.contains_points(locations))                           # True if pixel is in lake
            elif product == 'L3':
                l3_mask += np.array(landmass.contains_points(locations))                           # True if pixel is in island
            elif product == 'L4':
                l4_mask += np.array(landmass.contains_points(locations))                           # True if pixel is in pond (on island)
            if verbose:
                print(f"Completed product {product} coastline {n}")
                                                                    
    
    l12_mask = np.logical_and(l1_mask, np.invert(l2_mask))                        # combine first two
    l13_mask = np.logical_or(l12_mask, (l3_mask))                                 # and with 3
    l14_mask = np.logical_and(l13_mask, np.invert(l4_mask))                        # and with 4
    water_mask = np.invert(l14_mask)                                                  # water is where not land, so invert (now True for land)
    water_mask = np.flipud(np.reshape(water_mask, (ny, nx)))                                 # reshape, not sure why flipud as geographic but numpy is referring to top left pixel
    
    if debug_mode:    
        f, ax = plt.subplots(1,5)                                                     #debug plot
        ax[0].imshow(np.flipud(np.reshape(l1_mask, (ny, nx))))
        ax[1].imshow(np.flipud(np.reshape(l2_mask, (ny, nx))))
        ax[2].imshow(np.flipud(np.reshape(l3_mask, (ny, nx))))
        ax[3].imshow(np.flipud(np.reshape(l4_mask, (ny, nx))))
        ax[4].imshow(water_mask)
        
        # f, ax = plt.subplots(2)
        # ax.imshow(np.flipud(np.reshape(l12_mask, (ny, nx))))
        # f, ax = plt.subplots(1)
        # ax.imshow(np.flipud(np.reshape(l13_mask, (ny, nx))))
        # f, ax = plt.subplots(1)
        # ax.imshow(np.flipud(np.reshape(l14_mask, (ny, nx))))
    
    

    if verbose:
        print('Done!')
    
    return water_mask

#%%


def load_GSHHS_coastlines(gshhs_dir, resolution = 'i', bbox = None, level = '1'):
    """ Load the GSHHS (global self-consistent heirarchical high-resolution shorelines), available from:
        http://www.soest.hawaii.edu/pwessel/gshhg/
        Return as list of numpy arrays for shorelines.   Note, as this only works with shorelines, these are the L1 products 
        in the GSHHS docs, so it won't include lakes and rivers.  
        
    Inputs:
        gshhs_dir | Path | location of directory of shape files
        resoltuion | string | From GSHHS docs:
                                The geography data come in five resolutions:
                                    full resolution: Original (full) data resolution.
                                    high resolution: About 80 % reduction in size and quality.
                                    intermediate resolution: Another ~80 % reduction.
                                    low resolution: Another ~80 % reduction.
                                    crude resolution: Another ~80 % reduction.
        bbox | None or dict |  contains west east south north in degress, and ensures only an area of interest is returned.  
        level | int | From GSHHS docs:
                            L1: boundary between land and ocean, except Antarctica.
                            L2: boundary between lake and land.
                            L3: boundary between island-in-lake and lake.
                            L4: boundary between pond-in-island and island.
                            L5: boundary between Antarctica ice and ocean.
                            L6: boundary between Antarctica grounding-line and ocean.
            
            
        
    Returns:
        all_coastline_polygons | list of arrays | all coastlines in the file, as a list of numpy arrays
        intersect_coastline_polygons | list of arrays |  all coastlines that cross the bounding box, as list of arrays
        intersect_cropped_coastline_polygons | list of arrays |  all coastlines inside the bounding box (and use the edge of the bounding box for bits that go outside it), as list of arrays
        
    History:
        2022_01_26 | MEG | Written
        
    
    """
    import numpy as np
    import shapefile                                                                  # from pyshp, to load shape file
    from shapely.geometry import Polygon, Point                                       # to work with coastlines.  
    import matplotlib.pyplot as plt
    
    def polygon_to_nparray(polygon):
        """ Convert a shapely Polygon to a numpy array (lons in column 0, lats in column 1)
        """
        x, y = polygon.exterior.coords.xy                                 # get the coords of that
        x = np.array(x)                                                                     # make into a numpy array, rank 1
        y = np.array(y)                                                                     # make into a numpy array, rank 1
        xy = np.concatenate((x[:, np.newaxis], y[:, np.newaxis]), axis = 1)
        return xy
    

    # build the bounding box (bbox), if required.  
    if bbox is not None:
        bbox = Polygon([(bbox['west'], bbox['south']), 
                        (bbox['east'], bbox['south']),
                        (bbox['east'], bbox['north']),
                        (bbox['west'], bbox['north'])])
        
    shape_file = str(gshhs_dir / f"{resolution}" / f"GSHHS_{resolution}_L{level}.shp")
    sf = shapefile.Reader(shape_file)
    n_landmasses = len(sf.shapes())


    all_coastline_polygons = []                                                                     # all coastlines in the file, as a list of numpy arrays
    if bbox is not None:
        intersect_coastline_polygons = []                                                           # all coastlines that cross the bounding box, as list of arrays
        intersect_cropped_coastline_polygons = []                                                   # all coastlines inside the bounding box (and use the edge of the bounding box for bits that go outside it), as list of arrays

    for n_landmass in range(n_landmasses):                                                          # loop through each landmass (could be an island, could be all of eurasia)
        one_landmass = sf.shape(n_landmass).points                                                  # get the lon and lats of the points defining its coastline
        one_landmass_polygon = Polygon(one_landmass)                                                # convert to shapely polygon
        all_coastline_polygons.append(np.array(one_landmass))                                       # convert to numpy array and add to list
        
        # plt.switch_backend("Qt5Agg")                                                              # useful to debug - visualise the landmass and the bounding box.  
        # f, ax = plt.subplots()
        # ax.plot(np.array(one_landmass)[:,0], np.array(one_landmass)[:,1])
        # bbox_as_nparray = polygon_to_nparray(bbox)
        # ax.plot(bbox_as_nparray[:,0], bbox_as_nparray[:,1])

        if bbox is not None:
            if one_landmass_polygon.intersects(bbox):                                               # if any part of the coastline goes into the bounding box 
                intersect_coastline_polygons.append(np.array(one_landmass))                         # convert it to a numpy array and store
                
                intersect_shapely_polygon = one_landmass_polygon.intersection(bbox)                 # find which bits are in the bbox, this can be a Polygon or a MultiPolygon (if the coastline has now broken into two parts (e.g. a horshoe shaped island could become two islands if the bbox only capture the bottom part of it))
                if intersect_shapely_polygon.type == 'Polygon':
                    intersect_cropped_coastline_polygons.append(polygon_to_nparray(intersect_shapely_polygon))       # conver the Polygon to a numpy array and store
                elif intersect_shapely_polygon.type == 'MultiPolygon':
                    # f, ax = plt.subplots()                                                    # useful to debug and visualise why it's a MultiPolygon
                    # for p in intersect_shapely_polygon:
                    #     p_as_nparray = polygon_to_nparray(p)
                    #     ax.plot(p_as_nparray[:,0], p_as_nparray[:,1])
                    for p in intersect_shapely_polygon.geoms:                                             # loop through each polygon in the MultiPolygon (which is a separate coastline)
                        intersect_cropped_coastline_polygons.append(polygon_to_nparray(p))         # and append
              

    if bbox is None:
        return all_coastline_polygons
    else:
        return all_coastline_polygons, intersect_coastline_polygons, intersect_cropped_coastline_polygons

#%%



# def water_pixel_masker(dem, lon_lat_ll, lon_lat_ur, coast_resol, verbose = False, pixs_per_deg = 1201,
#                        debug_mode = False):
#     """
#        A function to creat a mask of pixels over water. This can be very slow for big DEMs (best called on each tile,
#        as per SRTM_dem_make).  
#        Note the too avoid problems with edges, the fucntion works by expanding many regions to be one tile larger in all directions.  
#        Note that basemap is rather an old package, and can be problematic on Windows machines.  
       
#     Inputs:
#         dem | rank 2 array | gridded data (e.g. a dem)
#         lon_lat_ll | tuple| lon lat of lower left corner of dem.  Must be floats! 
#         lon_lat_ur | tuple| lon lat of upper right corner of dem
#         coast_resol | str | resolution of vector coastlines: c l i h f 
#         pixs_per_deg | int | number of pixels in 1 degree.  1201 for 3 arc second SRTM3 data.  
#         debug_mode | Boolean | If true, a simple figure of the data and the coastlines is produced.  
       
#     Output:
#         result_ary | rank 2 array | array to be used as pixel mask
    
#     Hitory:
#         2017/03/01 | MEG  | adapterd from dem_show_oceans_detailed_2
#         2020/05/18 | MEG | Update to ensure that pixels at the edge that fall on the edge of a land polygon are included
#                             (by expanding the plot by 1 deg in each direction)
#         2020/06/02 | MEG | Fix bug when masking large DEMS (ie not just a single tile)
#         2020/07/29 | MEG | Fix bug for masking of lakes.  
#         2020/08/05 | MEG | Reduce extra area added around dem to speed up masking.  
#         2020/08/13 | MEG | Fix bug in water masking for tiles with no water bodies.  
#     """
         
#     from matplotlib.path import Path
#     import numpy as np
#     import matplotlib.pyplot as plt
#     import platform
#     import os

    
#     pdb.set_trace()
    
    
#     edge_fraction = 0.1
    
#     if platform.system() == 'Windows':                                                                              # check if windows as can be tricky with basemap
#         try:
#             print("Basemap is coming to the end of its life and support is now difficult with windows.  "
#                   "The os.environ argument PROJ_LIB often has to be set manually, but this is usually straightforward"
#                   " as most of it is just the path to your anaconda environment.  Trying to set it now, but if "
#                   " if this fails, it is most likely that you'll have to edit /user/ to your username in the water_pixel_masker fucntion... " , end = '')
#             os.environ['PROJ_LIB'] = r'C:\Users\User\anaconda3\envs\SRTM-DEM-tools\Library\share'                  # https://stackoverflow.com/questions/52911232/basemap-library-using-anaconda-jupyter-notebooks-keyerror-proj-lib/54087410#54087410
#             print(" Success!  ")
#         except:
#             print(f'Continuining anyway - look out for the PROJ_LIB error')
#     from mpl_toolkits.basemap import Basemap
    
#     if verbose:
#         print('Creating a mask of pixels that lie in water (sea and lakes)... ', end = '')
       
#     # 0: check inputs
#     if lon_lat_ll[0] > lon_lat_ur[0]:
#         raise Exception(f"The western (left) edge ({lon_lat_ll[0]}) appears to be east of the "
#                         f"eastern (right) edge ({lon_lat_ur[0]}).  Exiting.  ")
        
#     if lon_lat_ll[1] > lon_lat_ur[1]:
#         raise Exception(f"The southern (lower) edge ({lon_lat_ll[1]}) appears to be north of the "
#                         f"northern (upper) edge ({lon_lat_ur[1]}).  Exiting.  ")
            
#     # 1: deal with sizes of various things and make a grid of points
#     ny = dem.shape[0]                                                                                           # number of pixels vetically
#     nx = dem.shape[1]                                                                                           # and horizontally            
#     x = np.linspace(lon_lat_ll[0], lon_lat_ur[0], nx)
#     y = np.linspace(lon_lat_ll[1], lon_lat_ur[1], ny)
#     xx, yy = np.meshgrid(x,y)
#     locations = np.hstack((np.ravel(xx)[:,np.newaxis], np.ravel(yy)[:,np.newaxis]))                             # x and y (lon and lat) for all pixels in the dem, size (n_pixels x 2)
     
#     # 2 Make the basemap figure
#     ll_extent = [lon_lat_ll[0]-edge_fraction, lon_lat_ur[0]+edge_fraction, 
#                  lon_lat_ll[1]-edge_fraction, lon_lat_ur[1]+edge_fraction]                    # get the west east south north extent of the region needed by basemap.  Expand by 0.1 of a degree as edges can be difficult
#     plt.figure()                                                                                                # need a figure instance for basemap
#     map = Basemap(projection='cyl', llcrnrlat=ll_extent[2],urcrnrlat=ll_extent[3],
#                                     llcrnrlon=ll_extent[0],urcrnrlon=ll_extent[1], resolution=coast_resol)      # make the figure with coastlines, edges have already been expanded by ll extent
#     #import sys; sys.exit()
#     try:
#         map.drawcoastlines()
#         no_water_tile = False
#     except:
#         no_water_tile = True

    
#     if no_water_tile:
#         water_mask = np.zeros(dem.shape)
#     else:
#         land_mask = np.zeros(len(locations), dtype=bool)                                        # initialise as false
#         land_polygons = [Path(p.boundary) for p in map.landpolygons]                            # get a list of the land polygons.  Each "island" of land is own item
#         for polygon in land_polygons:                                                           # loop through each of these
#             land_mask += np.array(polygon.contains_points(locations))                           # True if pixel is in land
#         land_mask =  np.flipud(np.reshape(land_mask, (ny, nx)))                                 # reshape, not sure why flipud?
#         lake_mask = np.zeros(len(locations), dtype=bool)                                        # initialise as false
#         lake_polygons = [Path(p.boundary) for p in map.lakepolygons]                            # get a list of lakes
#         for polygon in lake_polygons:                                                           # loop through each of these     
#             lake_mask += np.array(polygon.contains_points(locations))                           # True if in lake   
#         lake_mask = np.flipud(np.reshape(lake_mask, (ny, nx)))                                  # reshape, again not sure why flipud
#         land_lake_mask = np.logical_and(land_mask, np.invert(lake_mask))                        # land is where land is true and lake is not true
#         water_mask = np.invert(land_lake_mask)                                                  # water is where not land

#     if debug_mode == False:
#         plt.close()
#     else:
#         #im1 = map.pcolormesh(x,y[::-1],dem,shading='flat',cmap=plt.cm.terrain,latlon=True, vmin = 0)
#         im1 = map.pcolormesh(x,y[::-1],land_lake_mask,cmap=plt.cm.terrain,latlon=True,)
    
#     if verbose:
#         print('Done!')
    
#     return water_mask