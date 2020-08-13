# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:08:03 2020

@author: User
"""




def SRTM_dem_make(west, east, south, north, SRTM1_or3 = 'SRTM3', water_mask_resolution = None,
                  SRTM1_tiles_folder = './SRTM1/', SRTM3_tiles_folder = './SRTM3/',
                  ed_username = None, ed_password = None, download = True, void_fill = False):
    """
    Given lons and lats (integers), make a multi tile DEM from either SRTM1 or 3 data.  
    Inputs:
        west | -179 -> 180 | west of GMT is negative
        east | -179 -> 180 | west of GMT is negative
        south | -90 -> 90  | northern hemishpere is positive
        north | -90 -> 90  | northern hemishpere is positive
        SRTM1_or3 | string | either use SRTM3 (~90m pixels) or SRTM1 (~30m pixels)
        water_mask_resolution | None or string | If not none, the DEM will be returned as a masked array, with water masked.  
                                        Resolution of vector coastlines: c l i h f   (ie coarse down to fine)
        SRTM1_tiles_folder | string | folder where SRTM1 .hgt file are kept
        SRTM3_tiles_folder | string | folder where SRTM3 .hgt file are kept
        
        ed_username | string | Earthdata username, needed for SRMT1 only.  To apply: https://earthdata.nasa.gov/eosdis/science-system-description/eosdis-components/earthdata-login
        ed_password | string | Earthdata password
        download | boolean | if False, function will not try to download DEM tiles.  Good if have a large number of tile already downloaded.  
        void_fill | boolean | If true, will try to linearly interpolate across voids in data.  
    
    Output:
        dem | rank 2 array | the dem
        lons | rank 1 array | longitude of bottom left of each 1' x1' grid cell
        lats | rank 1 array | latitude of bottom left of each 1' x1' grid cell
    """
    #import matplotlib.pyplot as plt
    import numpy as np
    import numpy.ma as ma
    #import os 

    if west > east:
        raise Exception(f"'west' ({west}) must always be smaller than 'east' ({east}).  If west of Grenwich, values should be negative.  Exiting. ")
    if south > north:
        raise Exception(f"'south' ({south}) must always be smaller than 'north' ({north}).  If south of the equator, values should be negative.  Exiting. ")

    print(f"Making a dem between {west} and {east} longitude, and {south} and {north} latitude.  ")
    if water_mask_resolution is None:
        print(f"{SRTM1_or3} tiles will be used, and water bodies won't be masked.  ")
    else:
        print(f"{SRTM1_or3} tiles will be used, and water bodies will be masked.  ")

    # Things to set
    null = -32768                                                               # from SRTM documentation   
    
    # 0: determine resolution and check inputs
    if SRTM1_or3 == 'SRTM3':
        samples = 1201
        SRTM3 = True
        tiles_folder = SRTM3_tiles_folder
    elif SRTM1_or3 == 'SRTM1':
        samples = 3601
        SRTM3 = False
        tiles_folder = SRTM1_tiles_folder
    else:
        raise Exception(f"'SRTM1_or3' must be eitehr SRTM1 or SRTM3.  Exiting...")
    
    if download is True or download is False:
        pass
    else:
        raise Exception(f"'download' must be eitehr 'True' or 'False'.  Exiting...")
        
    if water_mask_resolution == 'None':
        raise Exception(f"'water_mask_resolution' can be 'None' if you don't want to create one, but that is None "
                        f", and not the string 'None'.  Exiting...")
        
    
    
    # 1: Initiliase the big DEM:
    lats = np.arange(south, north, 1)
    lons = np.arange(west, east, 1)
    num_x_pixs = lons.size * samples
    num_y_pixs = lats.size * samples    

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
                tile_elev = open_hgt_file(tile_name, samples, samples, tiles_folder)
                print(' Done!')

            except:
                print(' Failed.  ')
                if download == True:
                    print(f"{tile_name} : Trying to download it...  ", end = "" )
                    if SRTM3:
                        try:
                            download_success = srtm3_tile_downloader(tile_name, tiles_folder)
                        except:
                            download_success = False
                    else:
                        try:
                            srtm1_tile_downloader(tile_name, ed_username, ed_password, tiles_folder)          # download tile
                            download_success = True
                        except:
                            download_success = False
                    if download_success:
                        print( 'Done!')
                        print(f"{tile_name} : Opening the .hgt file", end = "" )
                        tile_elev = open_hgt_file(tile_name, samples, samples, tiles_folder)
                        print( ' Done!')
                    else:
                        print(' Failed.  ')
                        
                else:
                    pass
                    
                if (download == False) or (download_success == False):
                    print(f"{tile_name} : Replacing with nulls (probably a water only tile).  ", end = "" )
                    tile_elev = null * np.ones((samples,samples))
                    replaced_with_null = True   
                    print( 'Done!')
                else:
                    pass
                        
            # 2: if required, fill voids in the tile
            if void_fill is True and np.min(tile_elev) == (-32768) and replaced_with_null is False:                             # if there is a void in the tile, and we've said that we want to fill voids, and it's not a tile we can't find and have filled with nulls
                print(f"{tile_name} : Filling voids in the tile... ", end = "")
                tile_elev = fill_gridata_voids(tile_elev)                
                print(' Done!')
    
            # 3: Make the water mask for that tile
            if  water_mask_resolution is not None:
                if replaced_with_null:                                                                                      # if it was replaced by null, should be water so don't need to make a mask
                    print(f"{tile_name} : Assuming water tile and masking it all...", end = '')    
                    tile_mask = np.ones((samples,samples))                                                                 # make a mask in which all pixels are masked (as it's assumed to be a water tile)
                else:
                    print(f"{tile_name} : Creating a mask of water areas...", end = '')    
                    tile_mask = water_pixel_masker(tile_elev, (lon, lat), (lon+1, lat+1), water_mask_resolution, 
                                                   pixs_per_deg = samples, verbose = False)                                 # or make a water mask for that tile
                print(" Done!")
            else:
                pass
                
            # 4: stitch the current tile and water mask into their respective full arrays 
            dem[num_y_pixs-((lat+1-lats[0])*samples) :num_y_pixs-((lat-lats[0])*samples), (lon-lons[0])*samples:(lon+1-lons[0])*samples] = tile_elev
            if  water_mask_resolution is not None:
                water_mask[num_y_pixs-((lat+1-lats[0])*samples) :num_y_pixs-((lat-lats[0])*samples), (lon-lons[0])*samples:(lon+1-lons[0])*samples] = tile_mask
    
            print('')                                                                                                        # next tile output will be on a new line, which makes it easier to read.  
            import matplotlib.pyplot as plt
            # f, ax = plt.subplots(1)
            # ax.imshow(water_mask)
            
            if lon == 0:
                f, ax = plt.subplots(1)
                ax.imshow(tile_mask)
            
    # combine the dem and the mask
    if  water_mask_resolution is not None:
        dem = ma.array(dem, mask = water_mask)
    
    
    return dem, lons, lats
    


#%%

def srtm1_tile_downloader(tile_name, ed_username, ed_password, hgt_path = './SRTM1_tiles', verbose = False):
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
            
        """
        import requests
        import zipfile
        import os
        
        tile_location = 'http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/'               # Where USGS keeps SRTM1 tiles
        
        # 0: download the zip file
        if verbose:
            print(f"{tile_name}: Downloading the zip file...", end = '')
        zip_filename = f"./{tile_name}.hgt.zip"                                                     # name to save downloaded file to
        with requests.Session() as session:
            session.auth = (ed_username, ed_password)                                               # login steps
            r1 = session.request('get', f"{tile_location}{tile_name}.SRTMGL1.hgt.zip")  
            r = session.get(r1.url, auth=(ed_username, ed_password))                                      # download file
            if r.ok:
                with open(zip_filename, 'wb') as fd:                                                # write the .hgt.zip file
                    for chunk in r.iter_content(chunk_size=1024*1024):
                        fd.write(chunk)
            else:
                raise Exception('Download failed.  ')
        if verbose:
            print('Done!')
                        
        # 1: unzip to get a hgt, and then delete zip file
        if verbose:
            print(f"{tile_name}: Unzipping...", end = '')                                               # unzip the file
        with zipfile.ZipFile(zip_filename ,"r") as zip_ref:                                
            zip_ref.extractall(hgt_path)
        os.remove(zip_filename)                                                                     # remove the redundant zip file
        if verbose:
            print("Done!")
  
#%%

def srtm3_tile_downloader(tile_name, hgt_path = './SRTM3_tiles/', verbose = False):
    """Given the name of an SRTM3 tile, download it, unzip it.        
    Inputs:
        tile_name | string | e.g. N51W004
        hgt_path | string | path to folder of SRTM1 tiles
    Returns:
        .hgt file
    History:
        2020/05/11 | MEG | Written
    """
    import wget
    import zipfile
    import os
    
    path_download = 'https://dds.cr.usgs.gov/srtm/version2_1/SRTM3/'            # other locations, but easiest to http from (except for the region thing)
    regions = ['Africa', 'Australia', 'Eurasia', 'Islands', 'North_America', 'South_America']
    
    
    # 0: Download the .hgt.zip
    if verbose:
        print(f"{tile_name}: Downloading the zip file...", end = '')
    download_success = False                                                                            # initialise
    for region in regions:                                                                             # Loop through all the SRTM regions (as from a lat and lon it's hard to know which folder they're in)
        if download_success is False:
            try:
                #import ipdb; ipdb.set_trace()
                filename = wget.download(f'{path_download}{region}/{tile_name}.hgt.zip', bar = False)                                                    # srtm data is in different folders for each region.  Try all alphabetically until find the ones it's in
                download_success = True
            except:
                download_success = False
    if verbose:
        print('Done!')
        
    # 1: unzip to get a hgt, and then delete zip file (if did download it)
    if download_success:
        if verbose:
            print(f"{tile_name}: Unzipping...", end = '')                                               # unzip the file
        with zipfile.ZipFile(f'./{tile_name}.hgt.zip',"r") as zip_ref:                                
            zip_ref.extractall(hgt_path)                                                            # unzip into folder of hgt files
        os.remove(f'./{tile_name}.hgt.zip')                                                 # remove the redundant zip file
        if verbose:
            print("Done!")
        
    return download_success    
  
   
#%%

def fill_gridata_voids(tile, void_value = -32786):
    """A function to fill voids in grided data using  scipy's interpolate function.
    Input:
        tile | rank 2 array | gridded data, containing voids
        void_value | ? | value that voids are filled with.  e.g. -32768 for SRTM data.  
    Returns:
        tile_filled  rank 2 array | gridded data, now without voids.  
    History:
        2020/08/13 | MEG | Modified from script.  
    """
    import numpy as np
    from scipy.interpolate import griddata             # for interpolating over any voids

    grid_x, grid_y = np.mgrid[0:tile.shape[1]:1, 0:tile.shape[0]:1 ]
    valid_points = np.argwhere(tile != void_value)                                                         # void are filled with -32768 perhaps?  less than -1000 should catch these
    tile_at_valid_points = tile[valid_points[:,0], valid_points[:,1]]
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
    elevations = np.fromfile(f'{tile_folder}{tile_name}.hgt', np.dtype('>i2'))                    # get as a rank 1 array
    tile_array = elevations.reshape((pixels_y, pixels_x))                                        # convert to rank 2
    if verbose:
        print('Done!')
    return tile_array


#%%


def dem_show(matrix,lons,lats,srtm, units_deg= True, title = None):
    """Visualise a DEM using lat lon as the axis
    Inputs:
        matrix | rank2 array | dem data to be plotted
        lons | list | lons of bottom left of each tile
        lats | list | lats of bottom left of each tile
        srtm1 | 1 or 3 | sets which resolution working with
         units_deg | boolean | IF ture, axes are in degrees and not pixels.  
         title | None | or string
    Returns:
        Figure
    History:
        2020/05/11 | MEG | added to project
        2020/05/11 | MEG | Add flag to swtich between degrees and pixels
        2020/06/04 | MEG | Add title option
    
    """
    
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap  
    
    cmap = plt.get_cmap('terrain')
    new_cmap = truncate_colormap(cmap, 0.2, 1)
    if srtm == 3:
        samples = 1200                                                          # pixels per deg
    else:
        samples = 3601                                                          # pixels per deg, from SRTM docs

    lats2 = np.arange(lats[0], lats[-1] +2 ,1)              # to help with plotting

    f, ax = plt.subplots()
    if title is not None:
        f.suptitle(title)
        f.canvas.set_window_title(title)
    plot_data = ax.imshow(matrix,interpolation='none', aspect=1, vmin = 0, vmax=np.max(matrix), cmap=new_cmap)
    f.colorbar(plot_data)
    if units_deg:
        ax.set_xticks([i for i in range(0,((lons.size)+1)*samples,samples)])                        # tick labels only for lats
        ax.set_xticklabels(lons) 
        ax.set_yticks([i for i in range(0,((lats2.size))*samples,samples)])                         # tick labels only for lats
        ax.set_yticklabels(lats2[::-1]) 

#%%


def water_pixel_masker(dem, lon_lat_ll, lon_lat_ur, coast_resol, verbose = False, pixs_per_deg = 1201,
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
    
    if platform.system() == 'Windows':                                                                              # check if windows as can be tricky with basemap
        try:
            print('Trying the Windows fix to get basemap to keep working...', end = '')
            os.environ['PROJ_LIB'] = r'C:\Users\User\anaconda3\envs\strava_sd_keras\Library\share'                  # https://stackoverflow.com/questions/52911232/basemap-library-using-anaconda-jupyter-notebooks-keyerror-proj-lib/54087410#54087410
            print(" Done.  ")
        except:
            print(f'Continuining anyway - look out for the PROJ_LIB error')
    from mpl_toolkits.basemap import Basemap
    
    if verbose:
        print('Creating a mask of pixels that lie in water (sea and lakes)... ', end = '')
                
    # 1: deal with sizes of various things and make a grid of points
    ny = dem.shape[0]                                                                                           # number of pixels vetically
    nx = dem.shape[1]                                                                                           # and horizontally            
    x = np.linspace(lon_lat_ll[0], lon_lat_ur[0], nx)
    y = np.linspace(lon_lat_ll[1], lon_lat_ur[1], ny)
    xx, yy = np.meshgrid(x,y)
    locations = np.hstack((np.ravel(xx)[:,np.newaxis], np.ravel(yy)[:,np.newaxis]))
     
    # 2 Make the basemap figure
    ll_extent = [lon_lat_ll[0]-edge_fraction, lon_lat_ur[0]+edge_fraction, 
                 lon_lat_ll[1]-edge_fraction, lon_lat_ur[1]+edge_fraction]                    # get the west east south north extent of the region needed by basemap.  Expand by 0.1 of a degree as edges can be difficult
    plt.figure()                                                                                                # need a figure instance for basemap
    map = Basemap(projection='cyl', llcrnrlat=ll_extent[2],urcrnrlat=ll_extent[3],
                                    llcrnrlon=ll_extent[0],urcrnrlon=ll_extent[1], resolution=coast_resol)      # make the figure with coastlines, edges have already been expanded by ll extent
    #import sys; sys.exit()
    try:
        map.drawcoastlines()
        no_water_tile = False
    except:
        no_water_tile = True

    
    if no_water_tile:
        water_mask = np.zeros(dem.shape)
    else:
        land_mask = np.zeros(len(locations), dtype=bool)                                        # initialise as false
        land_polygons = [Path(p.boundary) for p in map.landpolygons]                            # get a list of the land polygons.  Each "island" of land is own item
        for polygon in land_polygons:                                                           # loop through each of these
            land_mask += np.array(polygon.contains_points(locations))                           # True if pixel is in land
        land_mask =  np.flipud(np.reshape(land_mask, (ny, nx)))                                 # reshape, not sure why flipud?
        lake_mask = np.zeros(len(locations), dtype=bool)                                        # initialise as false
        lake_polygons = [Path(p.boundary) for p in map.lakepolygons]                            # get a list of lakes
        for polygon in lake_polygons:                                                           # loop through each of these     
            lake_mask += np.array(polygon.contains_points(locations))                           # True if in lake   
        lake_mask = np.flipud(np.reshape(lake_mask, (ny, nx)))                                  # reshape, again not sure why flipud
        land_lake_mask = np.logical_and(land_mask, np.invert(lake_mask))                        # land is where land is true and lake is not true
        water_mask = np.invert(land_lake_mask)                                                  # water is where not land

    if debug_mode == False:
        plt.close()
    else:
        #im1 = map.pcolormesh(x,y[::-1],dem,shading='flat',cmap=plt.cm.terrain,latlon=True, vmin = 0)
        im1 = map.pcolormesh(x,y[::-1],land_lake_mask,cmap=plt.cm.terrain,latlon=True,)
    
    if verbose:
        print('Done!')
    
    return water_mask



