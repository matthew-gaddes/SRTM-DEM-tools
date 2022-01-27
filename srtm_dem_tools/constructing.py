# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:08:03 2020

@author: User
"""


from pathlib import Path
import pdb

def SRTM_dem_make(dem_loc_size, SRTM1_or3 = 'SRTM3', water_mask_resolution = None,
                  SRTM1_tiles_folder = Path("./SRTM1/"), SRTM3_tiles_folder = Path('./SRTM3/'),
                  ed_username = None, ed_password = None, download = True, void_fill = False,
                  gshhs_dir = None):
    """
    Given lons and lats (integers), make a multi tile DEM from either SRTM1 or 3 data.   Note that there are two approached used for masking water bodies in the DEM.  For large DEMS, the mask must be calculated on a tile by tile basis 
    (or it becomes painfully slow, or exceeds RAM requirements), but for small DEMS it can be much faster to create the DEM, crop it (e.g. to a 20x20km size), and then calculate the mask.  
    This decision is made by determine_masking_stratergy
    
    Inputs:
        dem_loc_size | dict | This can be one of two styles:
                              1) 'centre', a lon lat tuple and 'side_length', a tuple of the x and y side lengths in m.  
                              2) 'west', 'east', 'south' 'north' bounding edges, in degrees.  
                                  west | -179 -> 180 | west of GMT is negative
                                  east | -179 -> 180 | west of GMT is negative
                                  south | -90 -> 90  | northern hemishpere is positive
                                  north | -90 -> 90  | northern hemishpere is positive
        SRTM1_or3 | string | either use SRTM3 (~90m pixels) or SRTM1 (~30m pixels)
        water_mask_resolution | None or string | If not none, the DEM will be returned as a masked array, with water masked.  
                                        Resolution of vector coastlines: c l i h f   (ie coarse down to fine)
        SRTM1_tiles_folder | Path | folder where SRTM1 .hgt file are kept
        SRTM3_tiles_folder | Path| folder where SRTM3 .hgt file are kept
        
        ed_username | string | Earthdata username, to apply: https://earthdata.nasa.gov/eosdis/science-system-description/eosdis-components/earthdata-login
        ed_password | string | Earthdata password
        download | boolean | if False, function will not try to download DEM tiles.  Good if have a large number of tile already downloaded.  
        void_fill | boolean | If true, will try to linearly interpolate across voids in data.  
        gshhs_dir | Path or None | Path to directory of the GSHHS shorelines.  Available from: http://www.soest.hawaii.edu/pwessel/gshhg/
    
    Output:
        dem | rank 2 array | the DEM, usually as a masked array with water bodies masked, but if resolution is "none" then this is just a numpy array.  
        lons_mg | rank 2 array | The lons of each pixel in the DEM.  Product of np.meshgrid
        lats_mg | rank 2 array | The lats of each pixel in the DEM.  Products of np.meshgrid.  
    History:
        2020/04/?? | MEG | Much of this was written during Covid-19 lockdown for a project work with Strava file.  
        2020/09/21 | MEG | Update to handle non-integer extents (ie. teh western edge need no longer be an integer).
        2020/09/30 | MEG | Change format so that DEMs can be specified as either a centre and side length, or as before as bounds (west, east, south, north)
        2020/10/01 | MEG | Fix a bug in the lats_mg that was causing the lowest lats to be on the top row.  
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import numpy.ma as ma
    import geopy
    from geopy.distance import geodesic
    
    from srtm_dem_tools.water_bodies import water_pixel_masker
    #import os 
    #plt.switch_backend('agg')                                                                       # with this backend, no windows are created during figure creation.                 
    
    
    def determine_masking_stratergy(west, east, south, north, area_threshold = 1.0):
        """ Determine the area of a DEM (in number of 1'x1' squares) and return if this is bigger than 
        a certain threshold.  
        Inputs:
            west | float or int | western edge of dem.  Same for other three.  
            area_threshold | float | if the area of the DEM is larger than this, return that masking should
                                      be done before cropping
        Returns:
            mask_before_crop | boolean | True is area of DEM is larger than the area_threshold
        History:
            2020/09/22 | MEG | Written.  
            2021/02/24 | MEG | Update to use a single function to download either SRTM1 or 3 tiles.  
            
        """
        del_x = east - west
        del_y = north - south
        if (del_x * del_y) > area_threshold:
            mask_before_crop = True
        else:
            mask_before_crop = False
        return mask_before_crop
    

    ########### Begin
    if gshhs_dir is None:
        water_mask_resolution = None
        print(f"As no directory was provided for the GSHHS coastline shapefiles, water masking can't be completed and 'water_mask_resolution' has been set to None.  ")
    null = -32768                                                               # from SRTM documentation   
    # -2: The DEM location and size can be given in one of two formats.  Unpack each of these to be in the same form
    dkeys = dem_loc_size.keys()                                                                               # get the keys that describe the DEM we are currently working on.  
    if ('centre' in dkeys) and ('side_length' in dkeys):                                                      # DEMs are described in one of two ways - either a centre and side length
        origin = geopy.Point(dem_loc_size['centre'][1], dem_loc_size['centre'][0])                            # geopy uses lat then lon notation, here for the centre of the DEM
        west = geodesic(meters=(dem_loc_size['side_length'][0])/2).destination(origin, 270)[1]                  # 
        east = geodesic(meters=(dem_loc_size['side_length'][0])/2).destination(origin, 90)[1]                  # 
        south = geodesic(meters=(dem_loc_size['side_length'][1])/2).destination(origin, 180)[0]                 # 
        north = geodesic(meters=(dem_loc_size['side_length'][1])/2).destination(origin, 000)[0]                 # 
        
    elif ('west' in dkeys) and ('east' in dkeys) and ('south' in dkeys) and ('north' in dkeys):                 # or in terms of bounds (edges)
        west = dem_loc_size['west']                                                                             # get edges of DEM
        east = dem_loc_size['east']
        south = dem_loc_size['south']
        north = dem_loc_size['north']
    else:
        raise Exception(f"'dem_loc_size' must be a dictionary containing either 'centre' (a lon lat tuple) and 'side_length' (x direciton y direction tuple), or  "
                        f"'west', 'east', 'south', 'north' (lons and lats as floats).  Exiting...")
    

    # -1 Check various input arguments
    if west > east:
        raise Exception(f"'west' ({west}) must always be smaller than 'east' ({east}).  If west of Grenwich, values should be negative.  Exiting. ")
    if south > north:
        raise Exception(f"'south' ({south}) must always be smaller than 'north' ({north}).  If south of the equator, values should be negative.  Exiting. ")

    print(f"Making a dem between {west} and {east} longitude, and {south} and {north} latitude.  ")
    if water_mask_resolution is None:
        print(f"{SRTM1_or3} tiles will be used, and water bodies won't be masked.  ")
    else:
        print(f"{SRTM1_or3} tiles will be used, and water bodies will be masked.  ")
        
    if (ed_username is None) or (ed_password is None):
        print('Both an Earthdata (USGS) username and password must be supplied to download data.  Continuing as the tiles may be stored locally.  ')
   

    # 0: determine resolution and check inputs
    if SRTM1_or3 == 'SRTM3':
        pixs2deg = 1201
        SRTM3 = True
        tiles_folder = SRTM3_tiles_folder
        SRTM_resolution = '3'
    elif SRTM1_or3 == 'SRTM1':
        pixs2deg = 3601
        SRTM3 = False
        tiles_folder = SRTM1_tiles_folder
        SRTM_resolution = '1'
    else:
        raise Exception(f"'SRTM1_or3' must be eitehr SRTM1 or SRTM3.  Exiting...")
    
    if download is True or download is False:
        pass
    else:
        raise Exception(f"'download' must be eitehr 'True' or 'False'.  Exiting...")
        
    if water_mask_resolution == 'None':
        raise Exception(f"'water_mask_resolution' can be 'None' if you don't want to create one, but that is None "
                        f", and not the string 'None'.  Exiting...")
        
    
    west_i, east_i, south_i, north_i = get_tile_edges(west, east, south, north)                                                      # edges can be floats.  Get the edges in whole numbers of tiles
    mask_before_crop = determine_masking_stratergy(west, east, south, north, area_threshold = 1.0)                                  # determine when water masking should occur.  If True, making done to each tile.  
    
    if not mask_before_crop:
        print(f"The area of the DEM is sufficently low that water masking and void filling (if requested) will be done for the whole DEM, and not for each tile.  ")
    
    # 1: Initiliase the big DEM:
    lats = np.arange(south_i, north_i, 1)
    lons = np.arange(west_i, east_i, 1)
    num_x_pixs = lons.size * pixs2deg
    num_y_pixs = lats.size * pixs2deg    

    dem = null * np.ones((num_y_pixs, num_x_pixs))                              # make the blank array of null values
    if water_mask_resolution is not None:
        water_mask = np.zeros((num_y_pixs, num_x_pixs))                                # make the blank array of 0s (ie nothing is masked)
    
    # 2: Work through each tile
    for lon in lons:                                                                                  # one column first, make the name for the tile to try and download
        for lat in lats:                                                                              # and then rows for that column
            replaced_with_null = False                                                                  # reset for each tile
            download_success = False                                                                # reset to haven't been able to download        
            tile_name = dem_tile_namer(lon, lat)                                                    # get name of tile in format used by USGS
            print(f"{tile_name} : Trying to open locally...", end = "")
            try:
                tile_elev = open_hgt_file(tile_name, pixs2deg, pixs2deg, tiles_folder)
                print(' Done!')
            except:
                print(' Failed.  ')
                if download == True:
                    print(f"{tile_name} : Trying to download it...  ", end = "" )
                    try:
                        srtm1or3_tile_downloader(tile_name, ed_username, ed_password, SRTM_resolution, hgt_path = tiles_folder)
                        download_success = True
                    except:
                        download_success = False

                    if download_success:
                        print( 'Done!')
                        print(f"{tile_name} : Opening the .hgt file", end = "" )
                        tile_elev = open_hgt_file(tile_name, pixs2deg, pixs2deg, tiles_folder)
                        print( ' Done!')
                    else:
                        print(' Failed.  ')
                        
                else:
                    pass
                    
                if (download == False) or (download_success == False):
                    print(f"{tile_name} : Replacing with nulls (probably a water only tile).  ", end = "" )
                    tile_elev = null * np.ones((pixs2deg,pixs2deg))
                    replaced_with_null = True   
                    print( 'Done!')
                else:
                    pass
                        
            # 2: if required, fill voids in the tile
            if void_fill and (np.min(tile_elev) == (-32768)) and (replaced_with_null is False) and mask_before_crop:                             # if there is a void in the tile, and we've said that we want to fill voids, and it's not a tile we can't find and have filled with nulls
                print(f"{tile_name} : Filling voids in the tile... ", end = "")
                tile_elev = fill_gridata_voids(tile_elev)                
                print(' Done!')
    
            # 3: Make the water mask for that tile
            if  (water_mask_resolution != None) and (mask_before_crop):                                                     # if we're masking water, and doing it before cropping.  
                if replaced_with_null:                                                                                      # if it was replaced by null, should be water so don't need to make a mask
                    print(f"{tile_name} : Assuming water tile and masking it all...", end = '')    
                    tile_mask = np.ones((pixs2deg,pixs2deg))                                                                 # make a mask in which all pixels are masked (as it's assumed to be a water tile)
                else:
                    print(f"{tile_name} : Creating a mask of water areas...", end = '')    
                    tile_mask = water_pixel_masker(tile_elev, (lon, lat), (lon+1, lat+1), water_mask_resolution, gshhs_dir, 
                                                   pixs_per_deg = pixs2deg, verbose = False)                                 # or make a water mask for that tile, note that bounds are just set to be for that tile.  
                print(" Done!")
            else:
                pass
                
            # 4: stitch the current tile and water mask into their respective full arrays 
            dem[num_y_pixs-((lat+1-lats[0])*pixs2deg) :num_y_pixs-((lat-lats[0])*pixs2deg), (lon-lons[0])*pixs2deg:(lon+1-lons[0])*pixs2deg] = tile_elev
            if  water_mask_resolution is not None and mask_before_crop:
                water_mask[num_y_pixs-((lat+1-lats[0])*pixs2deg) :num_y_pixs-((lat-lats[0])*pixs2deg), (lon-lons[0])*pixs2deg:(lon+1-lons[0])*pixs2deg] = tile_mask
            print('')                                                                                                        # next tile output will be on a new line, which makes it easier to read.  

    # 5: crop the dem                                                                                               # 
    ll_pixels = ll2xy(((lons[0], lats[0])), pixs2deg, np.array([[west, south]]))                                    # Conver the lon and lat limits of the DEM to pixels coordinates in the dem we've made, ll = lower left corner
    ur_pixels = ll2xy(((lons[0], lats[0])), pixs2deg, np.array([[east, north]]))                                    # ur = upper right corner
    
    if mask_before_crop:                                                                                                        # if we have already done maskign (ie before_crop is true)
        if  water_mask_resolution is not None:
            dem = ma.array(dem, mask = water_mask)                                                                                # possbily conver the DEM from an array to a masked array.
        dem = dem[dem.shape[0]-ur_pixels[0][1] : dem.shape[0]-ll_pixels[0][1] , ll_pixels[0][0] : ur_pixels[0][0] ]               # crop the DEM, which may be a masked array (see above), or it may just be an array
                                                                                                                                  # Nb matrix notation starts from the top left, so y has conversions from xy coordinate of pixel (which is from lower left)
    else:                                                                                                                         # if mask before crop is false, we still need to mask and void fill
        if  water_mask_resolution is not None:                                                                                    # 1: the water masking (determine if needed)
            print("Masking the water bodies in the entire DEM (as its area was determined to be below 1 tile)... ", end = '')
            dem = dem[dem.shape[0]-ur_pixels[0][1] : dem.shape[0]-ll_pixels[0][1] , ll_pixels[0][0] : ur_pixels[0][0] ]           # crop the DEM, which is still an array.  
            dem_mask = water_pixel_masker(dem, (west, south), (east, north), water_mask_resolution, gshhs_dir,
                                                       pixs_per_deg = pixs2deg, verbose = False)                                  # make the water mask for the entire DEM.         
            dem = ma.array(dem, mask = dem_mask)                                                                                  # convert the DEM from an array to a masked array.
            print("Done!")
        
        if void_fill and (np.min(ma.compressed(dem)) == (-32768)):                                                                  # 2: The void fill (determine if needed)
            print("Filling voids in the entire DEM (as its area was determined to be below 1 tile)... ", end = '')
            dem_filled = fill_gridata_voids(ma.getdata(dem))                                                                        # getdata turns the ma back to a numpy array
            dem = ma.array(dem_filled, mask = dem_mask)                                                                             # conver the void filled numpy array we just made back to the masked array.  
            print("Done!")

    # 6: make long and lats for each pixel in the DEM
    lons_mg, lats_mg = np.meshgrid(np.linspace(west, east-(1/pixs2deg), dem.shape[1]), np.linspace(south, north-(1/pixs2deg), dem.shape[0]))        # lons and lats are for the lower left corner of eahc pixel, so we stop (east and north)
    lats_mg = np.flipud(lats_mg)                                                                                                                    # flip so that lowest lats are at bottom of array                                            
    #plt.switch_backend('Qt5Agg')                                                                                                            # return to the standard interactive backend 
    
    return dem, lons_mg, lats_mg
  
#%%

def SRTM_dem_make_batch(list_dems,  SRTM1_or3 = 'SRTM3', water_mask_resolution = 'None',
                        SRTM1_tiles_folder = Path('./SRTM1/'), SRTM3_tiles_folder = Path('./SRTM3/'),
                        ed_username = None, ed_password = None, download = True, void_fill = False, gshhs_dir = None):
    """ Create multiple DEMS in one go.  
    Inputs:
        list_dems | list of dicts | each dem is an item in the list and can be described as either:
                                    1) 'centre', a lon lat tuple and 'side_length', a tuple of the x and y side lengths in m.  
                                    2) 'west', 'east', 'south' 'north' bounding edges, in degrees.  
        For all other args, see the doc for SRTM_dem_make
    Returns:
        list_dems | lis of dicts | Now also contains 'dem', 'lats' and 'lons' for each DEM.  Lats and lons are coordinates of the lower left corner of every pixel.  
    History:
        2020/09/23 | MEG | Written
    """
    import numpy as np
    import copy
    
    n_dems = len(list_dems)                                                                                             # get the total number of DEMs to make
    list_dems_with_dem = copy.deepcopy(list_dems)
    
    for n_dem, list_dem in enumerate(list_dems):                                                                                      # loop through each DEM to make it.  
        print(f"Starting DEM {n_dem} of {n_dems}.  ")       
        try:
            dem_loc_size = {}                                                                       # initiate a dictionary that will describe the DEM location and size.  
            try:
                dem_loc_size['centre'] = list_dem['centre']                                         # either as a centre and side length
                dem_loc_size['side_length'] = list_dem['side_length']
            except:
                dem_loc_size['west'] = list_dem['west']                                             # or just as bounds
                dem_loc_size['east'] = list_dem['east']
                dem_loc_size['south'] = list_dem['south']
                dem_loc_size['north'] = list_dem['north']
            
            #import pdb; pdb.set_trace()
                
            dem, lons, lats =  SRTM_dem_make(dem_loc_size,
                                             SRTM1_or3, water_mask_resolution, SRTM1_tiles_folder, SRTM3_tiles_folder, 
                                             ed_username, ed_password, download, void_fill, gshhs_dir)                                         # make the DEM, using either of the two ways the DEM can be described.  
            
            list_dems_with_dem[n_dem]['dem'] = dem                                                                                       # and write the products to the dictionary.  
            list_dems_with_dem[n_dem]['lons_mg'] = lons                                                                                     # cont'd
            list_dems_with_dem[n_dem]['lats_mg'] = lats                                                                                     # cont'd
                    
        except:
            print(f'Failed to create DEM {n_dem} but continuing to the next DEM.  ')                                                                                    # print warning message if the try statement fails.  
                
    return list_dems_with_dem

  
    
#%% should this be rewritten 

def get_tile_edges(west, east, south, north):
    """ Given DEM limits as floats (ie not integers), determine the limits of the DEM in integers.  
    e.g. if western limit is 3.5, this will return 3 as the western extent of the westernmost tile required.  
    Inputs:
        west | float or int | western edge of dem.  
        as above for other edges.  
    Returns:
        west_i | int | western edge of tile required to span the DEM
        as above for the other edges.  
    History:
        2020/09/21 | MEG | Written """
    import numpy as np
    
    west_i = int(np.floor(west))
    east_i = int(np.ceil(east))
    south_i = int(np.floor(south))
    north_i = int(np.ceil(north))
    return west_i, east_i, south_i, north_i

#%%

def srtm1or3_tile_downloader(tile_name, ed_username, ed_password, SRTM_resolution = '3', hgt_path = None, verbose = False):
    """Given the name of an SRTM1 tile, download it, unzip it.    
    An Earthdata account is needed to succesfully download tiles.  
    Inputs:
        tile_name | string | e.g. N51W004
        ed_username | string | Earthdata username.  To apply: https://earthdata.nasa.gov/eosdis/science-system-description/eosdis-components/earthdata-login
        ed_password | string | Earthdata password
        hgt_path | string | path to folder of SRTM1 tiles
    Returns:
        .hgt file
    History:
        2020/05/10 | MEG | Written
        2021/02/24 | MEG | Update to handle both SRTM1 and 3 tiles.  
        
    """
    import requests
    import zipfile
    import os
    
    SRTM1_tile_location = 'http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/'               # Where USGS keeps SRTM1 tiles
    SRTM3_tile_location = 'https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL3.003/2000.02.11/'                # ditto
    
    # 0:  Set strings according to if we're using SRTM 1 or 3 data
    if SRTM_resolution == '3':
        tile_location = SRTM3_tile_location
        tile_extension = 'SRTMGL3.hgt.zip'
        if hgt_path is None:
            hgt_path = './SRTM3_tiles/'
    elif SRTM_resolution == '1':
        tile_location = SRTM1_tile_location
        tile_extension = 'SRTMGL1.hgt.zip'
        if hgt_path is None:
            hgt_path = './SRTM1_tiles/'
    else:
        raise Exception("'SRTM_resolution' must be itehr '3' or '1' (i.e. strings, and not numbers).  Exiting. ")
    
    # 1: download the zip file
    if verbose:
        print(f"{tile_name}: Downloading the zip file...", end = '')
    zip_filename = f"./{tile_name}.hgt.zip"                                                     # name to save downloaded file to
    with requests.Session() as session:
        session.auth = (ed_username, ed_password)                                               # login steps
        r1 = session.request('get', f"{tile_location}{tile_name}.{tile_extension}")  
        r = session.get(r1.url, auth=(ed_username, ed_password))                                      # download file
        if r.ok:
            with open(zip_filename, 'wb') as fd:                                                # write the .hgt.zip file
                for chunk in r.iter_content(chunk_size=1024*1024):
                    fd.write(chunk)
        else:
            raise Exception('Download failed.  ')
    if verbose:
        print('Done!')
                    
    # 2: unzip to get a hgt, and then delete zip file
    if verbose:
        print(f"{tile_name}: Unzipping...", end = '')                                               # unzip the file
    with zipfile.ZipFile(zip_filename ,"r") as zip_ref:                                
        zip_ref.extractall(hgt_path)
    os.remove(zip_filename)                                                                     # remove the redundant zip file
    if verbose:
        print("Done!")
  
   
#%%

def fill_gridata_voids(tile, void_value = -32768):
    """A function to fill voids in grided data using  scipy's interpolate function.
    Input:
        tile | rank 2 array | gridded data, containing voids
        void_value | ? | value that voids are filled with.  e.g. -32768 for SRTM data.  
    Returns:
        tile_filled  rank 2 array | gridded data, now without voids.  
    History:
        2020/08/13 | MEG | Modified from script. 
        2021/08/04 | MEG | Fix bug in grid_x and grid_y (x and y dimensions were transposed, which only caused an error when interpolating parts of a tile that weren't equal in size in the x and y dimension)
    """
    import numpy as np
    from scipy.interpolate import griddata             # for interpolating over any voids
    
    grid_x, grid_y = np.mgrid[0:tile.shape[0]:1, 0:tile.shape[1]:1 ]                                    # make a grid of the x and y coords for each point (units are number of pixels, note that it takes n_rows then n_cols (ie y before x)
    valid_points = np.argwhere(tile != void_value)                                                         # void are filled with -32768 (for srtm).  The length of this vector (n_pixels x 2) should be less than the number of pixels in tile, if some are voids.  If no voids, equal.  
    tile_at_valid_points = tile[valid_points[:,0], valid_points[:,1]]                                       # tile values at each of the valid locations.  
    tile_filled = griddata(valid_points, tile_at_valid_points, (grid_x, grid_y), method='linear')
    return tile_filled
     
#%%
def dem_tile_namer(lon, lat):    
    """ Given longitude and latitude in the form of - for west south, conver to 
    awlays positive format prefixed by NESW format used by USGS.  
    """
    if lat >= 0 and lon >= 0:                       # north east quadrant
        tile_name = 'N' + str(lat).zfill(2) + 'E' + str(lon).zfill(3)                                    # zfill pads to the left with zeros so always 2 or 3 digits long. 
    if lat >= 0 and lon < 0:                        # north west quadant
        tile_name = 'N' + str(lat).zfill(2) + 'W' + str(-lon).zfill(3)
    if lat < 0 and lon >= 0:                        # south east quadrant
        tile_name = 'S' + str(-lat).zfill(2) + 'E' + str(lon).zfill(3)
    if lat < 0 and lon < 0:                         # south east quadrant
        tile_name = 'S' + str(-lat).zfill(2) + 'W' + str(-lon).zfill(3)
    return tile_name
    

#%%

def open_hgt_file(tile_name, pixels_y, pixels_x, tile_folder = './SRTM1/', verbose = False):
    """ Does this work with SRTM3 tiles?
    """
    import numpy as np
    
    #import ipdb; ipdb.set_trace()
    if verbose:
        print(f"{tile_name}: Opening the hgt file...", end = '')
    elevations = np.fromfile(tile_folder / f"{tile_name}.hgt", np.dtype('>i2'))                    # get as a rank 1 array
    tile_array = elevations.reshape((pixels_y, pixels_x))                                        # convert to rank 2
    if verbose:
        print('Done!')
    return tile_array



#%%

def ll2xy(bottom_left_ll, pix2deg, points_ll):
    """    
    Input:
        bottom_left_ll | 1x2 np.array | lon lat of bottom left pixel of xy space
        deg2pix | int | number of pixels in 1 deg (e.g. 1201 for SRTM3)
        points_ll  | nx2   np.array | n >= 1 for it to work (ie no 1d arrays, must be at least 1x2).  lons in column 1, lats in column 2
    Output:
        points_xy | nx2 | (x y) in pixels from lower left corner as intergers 
                                Careful, as matrix indices are from top left forner 
        
    xy space has to be orientated so that north is vertical (ie angles are not supported)
    
    2016/12/14 | MEG | written
    2020/08/06 | MEG | Change so that ll is lonlat.  

    """
    import numpy as np
    
    n_data, dims = points_ll.shape
    points_diff = points_ll - bottom_left_ll              # difference in degrees from bottom left 
    points_xy = points_diff * pix2deg
    #points_xy = np.roll(points_diff_pix, 1, 1)          # lat lon is yx, switch to xy
    points_xy = points_xy.astype(int)                   # pixels must be integers 
    return points_xy                      

   