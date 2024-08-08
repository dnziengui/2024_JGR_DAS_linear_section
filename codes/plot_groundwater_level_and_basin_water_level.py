#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:14:18 2024

@author: dnziengui
"""

import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from MyToolbox.utils.utils import set_fonts
from MyToolbox.cc.hydro.functions import (
    _load_aggregated_data_piezo,
    _labels_for_piezo_and_basin,
    _dataframe_columns_and_time_selections,
    _plot_piezo_data3
)

# Set fonts for the plots
set_fonts(18)

#%%
#===============================
#         parameters  
#===============================
# Configuration
save_figure = False  # Flag to save figures
figure_directory = "../figures"



#%%
#===============================
#       plot figure 
#===============================
# Plotting parameters
colormap = 'tab20c'
label_size = 18

# Create the figure and axis for the plot
fig, ax = plt.subplots(figsize=(16, 5), constrained_layout=True)

# Plot piezometric and basin data
ax, ax2 = _plot_piezo_data3(
    data_piezo, data_basin, masked_datetime_index, offset_datetime, labels_piezo, labels_basin,
    fig=fig, ax=ax, cmap=colormap, TMIN=start_time, TMAX=end_time, lw_piezo=5, lw_basin=4,
    freq_ticks='1D', legend=False, zeros0=False, labelsize=label_size, color_piezo="dodgerblue",
    color_basin="cyan", ax1_ylim=None, ax2_ylim=(-0.1, 2.5)
)



plt.show()


#%%
#===============================
#      save figure 
#===============================

# Save the figure if save_figure is True
if save_figure:

    output_filename = os.path.join(figure_directory, "02_chronogram_piezometric_and_basin_water_level.pdf")
    fig.savefig(output_filename)

plt.show()