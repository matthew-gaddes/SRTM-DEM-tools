# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:08:03 2020

@author: User
"""




def SRTM_dem_make(west, east, south, north, SRTM1_or3 = 'SRTM3',
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
    #import os 
    from scipy.interpolate import griddata             # for interpolating over any voids

    # Things to set
    null = -32768                                                               # from SRTM documentation   
    
    # 0: determine resolution
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
    
    # 1: Initiliase the big DEM:
    lats = np.arange(south, north, 1)
    lons = np.arange(west, east, 1)
    num_x_pixs = lons.size * samples
    num_y_pixs = lats.size * samples    

    dem = null * np.ones((num_y_pixs, num_x_pixs))             # make the blank array of null values
    
    # 2: Work through each tile
    for lon in lons:                                                                                  # one column first, make the name for the tile to try and download
        for lat in lats:                                                                              # and then rows for that column
            void_fill_skip = False                                                                  # reset for each tile
            download_success = False                                                                # reset to haven't been able to download        
            tile_name = dem_tile_namer(lon, lat)                                                    # get name of tile in format used by USGS
            try:
                print(f"{tile_name} : Trying to open locally...", end = "")
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
                    
                if (download == 0) or (download_success == False):
                    print(f"{tile_name} : Replacing with null values (the tile probably doesn't exist and covers only water.  ", end = "" )
                    tile_elev = null * np.ones((samples,samples))
                    void_fill_skip = True   
                    print( 'Done!')
                else:
                    pass
                        
                
            #2: if required, fill voids in the tile
            if void_fill is True and np.min(tile_elev) == (-32768) and void_fill_skip is False:                             # if there is a void in the tile, and we've said that we want to fill voids.  
                print(f"{tile} : Filling voids in the tile... ", end = "")
                grid_x, grid_y = np.mgrid[0:samples:1  ,0:samples:1 ]
                valid_points = np.argwhere(tile_elev > (-32768))                                               # void are filled with -32768 perhaps?  less than -1000 should catch these
                tile_elev_at_valid_points = tile_elev[valid_points[:,0], valid_points[:,1]]
                tile_elev = griddata(valid_points, tile_elev_at_valid_points, (grid_x, grid_y), method='linear')
                print(' Done!')
    
            # 3: stitch the current tile into the full DEM    
            dem[num_y_pixs-((lat+1-lats[0])*samples) :num_y_pixs-((lat-lats[0])*samples), (lon-lons[0])*samples:(lon+1-lons[0])*samples] = tile_elev
    
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
    
    if verbose:
        print(f"{tile_name}: Opening the hgt file...", end = '')
    elevations = np.fromfile(f'{tile_folder}/{tile_name}.hgt', np.dtype('>i2'))                    # get as a rank 1 array
    tile_array = elevations.reshape((pixels_y, pixels_x))                                        # convert to rank 2
    if verbose:
        print('Done!')
    return tile_array


#%%


def dem_show(matrix,lons,lats,srtm, units_deg= True):
    """Visualise a DEM using lat lon as the axis
    Inputs:
        matrix | rank2 array | dem data to be plotted
        lons | list | lons of bottom left of each tile
        lats | list | lats of bottom left of each tile
        srtm1 | 1 or 3 | sets which resolution working with
         units_deg | boolean | IF ture, axes are in degrees and not pixels.  
    Returns:
        Figure
    History:
        2020/05/11 | MEG | added to project
        2020/05/11 | MEG | Add flag to swtich between degrees and pixels
    
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

    plt.figure()
    plt.imshow(matrix,interpolation='none', aspect=1, vmin = 0, vmax=np.max(matrix), cmap=new_cmap)
    cbar = plt.colorbar()
    if units_deg:
        ax = plt.gca()
        ax.set_xticks([i for i in range(0,((lons.size)+1)*samples,samples)])                        # tick labels only for lats
        ax.set_xticklabels(lons) 
        ax.set_yticks([i for i in range(0,((lats2.size))*samples,samples)])                         # tick labels only for lats
        ax.set_yticklabels(lats2[::-1]) 

#%%



def water_pixel_masker(data, lons, lats, coast_resol, verbose = False):
    """
       A function to creat a mask of pixels over water. This can be very slow for big DEMs
    Inputs:
        data | rank 2 array | gridded data (e.g. a dem)
        lons | list | lons of bottom left of each tile
        lats | list | lats of bottom left of each tile
        coast_resol | str | resolution of vector coastlines: c l i h f 
       
    Output:
        result_ary | rank 2 array | array to be used as pixel mask
        
    2017/03/01 | adapterd from dem_show_oceans_detailed_2
    """
     
    
    from matplotlib.path import Path
    import numpy as np
    import numpy.ma as ma
    import matplotlib.pyplot as plt
    import platform
    import os
    
    if platform.system() == 'Windows':                                                                              # check if windows as can be tricky with basemap
        try:
            print('Trying the Windows fix to get basemap to keep working...', end = '')
            os.environ['PROJ_LIB'] = r'C:\Users\User\anaconda3\envs\strava_sd_keras\Library\share'                  # https://stackoverflow.com/questions/52911232/basemap-library-using-anaconda-jupyter-notebooks-keyerror-proj-lib/54087410#54087410
            print(" Done.  ")
        except:
            print(f'Continuining anyway - look out for the PROJ_LIB error')
    from mpl_toolkits.basemap import Basemap
    
    if verbose:
        print('Creating a mask of pixels that lie in water (sea and lakes).  This can be slow')
    
    
    
    ll_extent = [lons[0], (lons[-1]+1), lats[0], (lats[-1]+1)]
    ny = data.shape[0]; nx = data.shape[1]
    plt.figure()                                                                # need a figure instance for basemap
    map = Basemap(projection='cyl', llcrnrlat=ll_extent[2],urcrnrlat=ll_extent[3],llcrnrlon=ll_extent[0],urcrnrlon=ll_extent[1], resolution=coast_resol)
    map.drawcoastlines()
    
    
    mesh_lons, mesh_lats = map.makegrid(nx, ny)               # get lat/lons of in evenly spaced grid with increments to match data matrix
    lons = np.ravel(mesh_lons)
    lats = np.ravel(mesh_lats)
    x, y = map(lons, lats)
    locations = np.c_[x, y]                                                     # concatenate to one array
    
    result = np.zeros(len(locations), dtype=bool)                                  # initialise as false
    land_polygons = [Path(p.boundary) for p in map.landpolygons]                    # check if land
    for polygon in land_polygons:
        result += np.array(polygon.contains_points(locations))
    lake_polygons = [Path(p.boundary) for p in map.lakepolygons]                    # check if in lake
    for polygon in lake_polygons:
        result = np.invert(np.array(polygon.contains_points(locations)))              # pixels in lakes are 1s, so subtract these.  
    result = np.invert(result)                                                       # true if in see, so masked out
    result_ary = np.flipud(np.reshape(result, (ny, nx)))
    plt.close()                                                                 # close figure instance
    return result_ary




