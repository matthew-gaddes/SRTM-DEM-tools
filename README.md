# SRTM-DEM-tools
Tools for making and manipulating SRTM1 and SRTM3 DEMs.  Note that only a small amount of the GSHHS data is included here to give a working example.  The full dataset can be downloaded from http://www.soest.hawaii.edu/pwessel/gshhg/

# New in January 2022:
Publication style figures can now be created with very few settings.  See "publication_example.py"
![image](https://user-images.githubusercontent.com/10498635/151543682-7b063e31-9b1b-44cb-b191-2b17460a13c2.png)




# Installation
A suitable python environment can be created using the conda command (for version 3.X.X):<br>
<code>conda env create --file srtm_dem_tools_v3.yml</code>

<br>

# Usage
DEMs an be described either using using west/east/south/north limits: <br>
<code> SRTM_dem_make({'west':-4, 'east':-1, 'south':53, 'north':55})</code><br>

Or, using a centre location and side length (in metres): <br>
<code> SRTM_dem_make({'centre': (-3.396, 37.041), 'side_length':(110e3, 100e3)}, </code><br>


There are also other useful functions: <br>
- water_pixel_masker   | To determine if pixels lie within water bodes.  This is usually called from within SRTM_dem_make.  
- SRTM_dem_make_batch  | To create several DEMs in one step.  
- dem_show             | To create simple figures showing the DEMs that have been created.  

# Examples

"example.py" includes several examples, including of an SRTM3 DEM with water bodies masked:

![SRTM3_DEM](https://user-images.githubusercontent.com/10498635/83517618-12473c80-a4d1-11ea-9645-37e6b74ffa34.png)


and of an SRTM1 DEM:

![SRTM1_DEM](https://user-images.githubusercontent.com/10498635/83517667-2b4fed80-a4d1-11ea-877d-16fb8a0d8efe.png)

Note, however, that to download SRTM1 tiles, an account with USGS Earthdata will be required - see: https://urs.earthdata.nasa.gov/
