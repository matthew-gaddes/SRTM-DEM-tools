#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 09:55:05 2022

@author: matthew
"""

#%%

def dem_show(matrix, lons_mg, lats_mg, title = None):
    """Visualise a DEM using lat lon as the axis
    Inputs:
        matrix | rank2 array | dem data to be plotted
        lons_mg | rank 2 array | The lons of each pixel in the DEM.  Product of np.meshgrid
        lats_mg | rank 2 array | The lats of each pixel in the DEM.  Products of np.meshgrid.  
        title | None or string | The figure title.  
    Returns:
        Figure
    History:
        2020/05/11 | MEG | added to project
        2020/05/11 | MEG | Add flag to swtich between degrees and pixels
        2020/06/04 | MEG | Add title option
        2020/09/22 | MEG | Change lons and lats from list of integers to meshgrids
        2020/10/01 | MEG | Fix a bug in how lats on the tick labels were handled.  
    """
    
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.ticker as mticker
    import numpy as np
    
    def create_tick_labels_in_deg(tick_labels, degs):
        """tick_labels are in pixel number, but this function converts them to degrees. 
        Inputs:
            tick_labels | rank 1 array | tick lables, such as you would get from ax.get_xticks().  These would be pixel numbers.  e.g. [0, 500, 1000]
            degs | rank 1 array | lons or lats of each pixel, as you would get from running np.meshgrid.  Should be the same size as the dem, and the lon or
                                  lat values at whichever pixel number has a tick lable is returned for use as new tick labels.  
        Returns:
            tick_labels_degs | list | tick labels in degrees, rounded to 2 dp.  
        History:
            2020/09/22 | MEG | Written
        """
        import numpy as np
        
        tick_labels_degs = []
        for tick_label in tick_labels.astype(int):
            if tick_label < 0:
                tick_labels_degs.append(' ')
            else:
                try: 
                    tick_labels_degs.append(np.round(degs[tick_label], 2))
                except:
                    tick_labels_degs.append(' ')
        return tick_labels_degs
    
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap  

    ################## Begin
    cmap = plt.get_cmap('terrain')                                                              # get a colourmap
    new_cmap = truncate_colormap(cmap, 0.2, 1)                                                  # remove the blues at the lowest range, so that 0 plots as green.  

    f, ax = plt.subplots()                                                                      # initatite figure.  
    if title is not None:
        f.suptitle(title)                                                                       # possibly add a title
        f.canvas.set_window_title(title)
    plot_data = ax.imshow(matrix,interpolation='none', aspect=1, vmin = 0, vmax=np.max(matrix), cmap=new_cmap)
    f.colorbar(plot_data)                                                                        # add a colour bar.  
        
    old_ticks = ax.get_xticks()                                                                 # upate the x ticks, now suits version ? of matplotlib
    ax.xaxis.set_major_locator(mticker.FixedLocator(old_ticks))
    ax.set_xticklabels(create_tick_labels_in_deg(old_ticks, lons_mg[-1,:]))
    
    old_ticks = ax.get_yticks()                                                                 # update the y ticks, now suits version ? of matplotlib
    ax.yaxis.set_major_locator(mticker.FixedLocator(old_ticks))
    ax.set_yticklabels(create_tick_labels_in_deg(old_ticks, lats_mg[-1,:]))

