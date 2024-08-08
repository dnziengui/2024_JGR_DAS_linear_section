#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:42:20 2024

@author: destin nziengui Bâ (destin.nziengui@febus-optics.com / destin.nziengui-ba@univ-grenoble-alpes.fr)

Codes to plot raw and filtered virtual shot gathers 
"""

#%%
#===============================
#         modules 
#===============================
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from MyToolbox.Seismic_filtering.FK import (
    interstation_distances,
    interstation_distance_binning,
    super_gather
)
from MyToolbox.Seismic_filtering.FK_filters import (
    FK_filter_0_,
    FK_filter_1_,
    FK_filter_2_
)
from MyToolbox.Seismic_filtering.FK import (
    fk_filtering_plot,
    fk_filtering_plot_dispersion_curves,
    fk_filtering_1,
    fk_filtering_plot_article
)


from MyToolbox.utils.utils import set_fonts


# Set fonts for plotting
set_fonts(18)

#%% d
#===============================
#   DATA LOADING & PARAMETERS   
#===============================

#
# open .h5 file  
#

DATA_PATH  ="../data"
fname      = "Virtual-shot-gathers_source-CH1125_FrequencyBand-1-20Hz.h5"
fname      = os.path.join(DATA_PATH,fname)
f          = h5py.File(fname, 'r')

#
# extract data 
#

vsg_raw_reference= np.squeeze(f['/data/raw_vsg/reference'][:])
vsg_raw_stimulation = np.squeeze(f['/data/raw_vsg/stimulation'][:])
vsg_stacked_reference = np.squeeze(f['/data/stacked_vsg/reference'][:])
vsg_stacked_stimulation = np.squeeze(f['/data/stacked_vsg/stimulation'][:])


#
# extract metadata 
#
lags        = np.squeeze(f['/metadata/lags'][:])                            # time lags (in sec)
distances   = np.squeeze(f['/metadata/distance_from_virtual_source'][:])    # distance from the virtual source ( in meters)
receivers   =  np.squeeze(f['/metadata/receivers'][:])                      # virtual receivers 
fs          = np.squeeze(f['/metadata/frequency_sampling'][:])              # frequency sampling 
dx          = np.squeeze(f['/metadata/spatial_sampling'][:])                # spatial sampling 

#
# figure parameters 
#
save_figure = True  # Flag to save figures
figure_directory = "../figures"




#%%

#===============================
#   PLOT FIGURES 
#===============================

#
# Define FK filtering parameters
#
# 
npadd = [1024]
nquadrants = 2
fk_filter = np.array(FK_filter_2_)


#
# plot raw virtual shot gathers 
#

fig, _ = fk_filtering_plot_article(
    
vsg_raw_reference , vsg_raw_stimulation, lags, distances,
    option="raw", filter_=fk_filter, fs=None, dx=None,
    npadd=npadd, nquadrants=nquadrants, lags_params=None,
    distances_params=None, frequencies_params=None,
    wavenumbers_params=None, npts_fk_filter=15,
    amax_norm_xcorr=[5., 1.]
)


if save_figure:

    output_filename = os.path.join(figure_directory, "04_raw_virtual_shots_gathers.pdf")
    fig.savefig(output_filename)


#
# plot filtered virtual shot gathers 
#

fig, _ = fk_filtering_plot_article(
    vsg_raw_reference , vsg_raw_stimulation, lags, distances,
    option="filter", filter_=fk_filter, fs=None, dx=None,
    npadd=npadd, nquadrants=nquadrants, lags_params=None,
    distances_params=None, frequencies_params=None,
    wavenumbers_params=None, npts_fk_filter=15,
    amax_norm_xcorr=[1., 1.]
)


if save_figure:

    output_filename = os.path.join(figure_directory, "S_filtered_virtual_shot_gathers.pdf")
    fig.savefig(output_filename)


