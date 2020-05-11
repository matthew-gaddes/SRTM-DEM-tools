# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:08:03 2020

@author: User
"""




def SRTM1_dem_make(dem_limits, ed_username, ed_password, tiles_folder = './SRTM1/'):
    """
    Given lons and lats (integers), make a multi tile DEM.  
    Inputs:
        dem_limits | list | [west east south north]
        ed_username | string | Earthdata username.  To apply: https://earthdata.nasa.gov/eosdis/science-system-description/eosdis-components/earthdata-login
        ed_password | string | Earthdata password
        path_tile | string | folder where SRTM1 .hgt file are kept
    
    Output:
        dem | rank 2 array | the dem
        lons | rank 1 array | longitude of bottom left of each 1' x1' grid cell
        lats | rank 1 array | latitude of bottom left of each 1' x1' grid cell
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import os 
    from scipy.interpolate import griddata             # for interpolating over any voids

    def srtm1_tile_downloader(tile_name, ed_username, ed_password, hgt_path = './SRTM1_tiles'):
        """Given the name of an SRTM1 tile, download it, unzip it.    
        An Earthdata account is needed to succesfully download tiles.  
        Inputs:
            tile_name | string | e.g. N51W004
            ed_username | string | Earthdata username.  To apply: https://earthdata.nasa.gov/eosdis/science-system-description/eosdis-components/earthdata-login
            ed_password | string | Earthdata password
            hgt_path | string | path to folder of SRTM1 tiles
            pixels_x | int | number of pixels in x direction
            pixels_y | int | number of pixels in y direction
        Returns:
            dem_array | rank2 numpy array | tile as a numpy array
        History:
            2020/05/10 | MEG | Written
            
        """
        import requests
        import zipfile
        import os
        
        tile_location = 'http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/'               # Where USGS keeps SRTM1 tiles
        
        # 0: download the zip file
        print(f"{tile_name}: Downloading the zip file...", end = '')
        zip_filename = f"./{tile_name}.hgt.zip"                                                     # name to save downloaded file to
        with requests.Session() as session:
            session.auth = (ed_username, ed_password)                                               # login steps
            r1 = session.request('get', f"{tile_location}{tile_name}.SRTMGL1.hgt.zip")  
            r = session.get(r1.url, auth=(username, password))                                      # download file
            if r.ok:
                with open(zip_filename, 'wb') as fd:                                                # write the .hgt.zip file
                    for chunk in r.iter_content(chunk_size=1024*1024):
                        fd.write(chunk)
        print('Done!')
                        
        # 1: unzip to get a hgt, and then delete zip file
        print(f"{tile_name}: Unzipping...", end = '')                                               # unzip the file
        with zipfile.ZipFile(zip_filename ,"r") as zip_ref:                                
            zip_ref.extractall(hgt_path)
        os.remove(zip_filename)                                                                     # remove the redundant zip file
        print("Done!")

    
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
        return tile
    
    def open_hgt_file(tile_name, tile_folder = './SRTM1/', pixels_y = 3601, pixels_x = 3601):
        """ Does this work with SRTM3 tiles?
        """
        print(f"{tile_name}: Opening the hgt file...", end = '')
        elevations = np.fromfile(f'{tile_folder}/{tile_name}.hgt', np.dtype('>i2'))                    # get as a rank 1 array
        tile_array = elevations.reshape((pixels_y, pixels_x))                                        # convert to rank 2
        print('Done!')
        return tile_array
        
    
    # make big dem
    west = dem_limits[0]
    east = dem_limits[1]
    south = dem_limits[2]
    north = dem_limits[3]
    samples = 3601                                                                                  # pixels in 1' of latitude/longitude. 1200 for SRTM3, 3601 for SRTM1
    null = -32768
    
    lats = np.arange(south, north, 1)
    lons = np.arange(west, east, 1)
    num_x_pixs = lons.size * samples
    num_y_pixs = lats.size * samples    
    dem = null * np.ones((num_y_pixs, num_x_pixs))             # make the blank array of null values
    
    # get list of tiles to download
    for lon in lons:                                                                                  # one column first, make the name for the tile to try and download
        for lat in lats:                                                                              # and then rows for that column
            tile_name = dem_tile_namer(lon, lat)                                                    # get name of tile in format used by USGS
        
            try:
                tile_elev = open_hgt_file(tile_name, tiles_folder)                                   # try to open the .hgt file
            except:
                srtm1_tile_downloader(tile_name, ed_username, ed_password, hgt_path = './SRTM1_tiles')          # download tile
                tile_elev = open_hgt_file(tile_name, tiles_folder)                                          # try to open the .hgt file again, after downloading it
                
                
                
                
                
                
                if download == 1:
                    print(tile + " : Couldn't find in local folder so trying to download it...  ", end = "" )
                    for region in regions:
                        if download_success is False:
                            try:
                                tile_downloader(region , tile)                                                    # srtm data is in different folders for each region.  Try all alphabetically until find the ones it's in
                                download_success = True
                            except:
                                download_success = False
                    if download_success:
                        with zipfile.ZipFile(tile + '.hgt.zip' ,"r") as zip_ref:                                # if we downloaded a file, need to now delete it
                            zip_ref.extractall(path_tiles)
                        os.remove(tile + '.hgt.zip')                                                                # delete the downloaded zip file
                        tile_elev = read_hgt_file(path_tiles + '/' + tile + '.hgt', samples+1)                  # look for tile in tile folder
                        void_fill_skip = False
                    else:
                        print(tile + " : Couldn't download tile so replacing with null values.  " )
                        tile_elev = null * np.ones((samples,samples))   
                        void_fill_skip = True                               # it's all a void, so don't try and fill it
                else:                                                                                   # if download is disabled, just fill with null value _-32768
                    print(tile + " : Tile downloading disabled so replacing with null values.  " )
                    tile_elev = null * np.ones((samples,samples))   
                    void_fill_skip = True  

    
            if void_fill is True and np.min(tile_elev) == (-32768) and void_fill_skip is False:
                print('Filling voids in the DEM using linear interpolation... ', end = "")
                grid_x, grid_y = np.mgrid[0:1200:1  ,0:1200:1 ]
                valid_points = np.argwhere(tile_elev > (-32768))                                               # void are filled with -32768 perhaps?  less than -1000 should catch these
                tile_elev_at_valid_points = tile_elev[valid_points[:,0], valid_points[:,1]]
                tile_elev = griddata(valid_points, tile_elev_at_valid_points, (grid_x, grid_y), method='linear')
                print('Done!')
                
            dem[num_y_pixs-((k+1-lats[0])*samples) :num_y_pixs-((k-lats[0])*samples), (j-lons[0])*samples:(j+1-lons[0])*samples] = tile_elev[0:1200,0:1200]
    return dem, lons, lats
    

