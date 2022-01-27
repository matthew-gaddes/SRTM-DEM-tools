#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 11:37:37 2022

@author: matthew
"""

import sys
import srtm_dem_tools
from srtm_dem_tools.constructing import SRTM_dem_make, SRTM_dem_make_batch
from srtm_dem_tools.plotting import dem_show
from srtm_dem_tools.water_bodies import water_pixel_masker

import numpy.ma as ma                                                                 # used for DEMS and water masks 
from pathlib import Path




#%% Things to set



#%%

