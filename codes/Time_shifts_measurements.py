#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:27:48 2024

@author: dnziengui
"""

import os 
import numpy as np  
import matplotlib.pyplot as plt 
#from functions import load_FK_filtered_premanip_and_manip_data,_preprocessing
#from functions import compute_time_shift
from MyToolbox.Seismic_filtering.FK import plot_wf_and_vels
from MyToolbox.utils.utils import set_fonts
from MyToolbox.delays.xwt import xwt
from MyToolbox.delays.compute_dt import extract_and_window_data,compute_weights ,compute_weighted_mean_and_std,TimeShiftXWT
from MyToolbox.tomo.tomo_fns import (
    create_dvv_colormap,
)
import matplotlib.gridspec as gridspec

from MyToolbox.utils.utils import tic, toc 
import h5py 

# set fontsize 
set_fonts(14)

#
# --- general parameters 
#


#%%

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


xcorr_premanip = np.squeeze(f['/data/stacked_vsg/reference'][:])[:,::-1]
xcorr_manip = np.squeeze(f['/data/stacked_vsg/stimulation'][:])[:,::-1]


#
# extract metadata 
#
lags        = np.squeeze(f['/metadata/lags'][:])                            # time lags (in sec)
distances   = np.squeeze(f['/metadata/distance_from_virtual_source'][:])    # distance from the virtual source ( in meters)
receivers   =  np.squeeze(f['/metadata/receivers'][:])                      # virtual receivers 
fs          = np.squeeze(f['/metadata/frequency_sampling'][:])              # frequency sampling 
dx          = np.squeeze(f['/metadata/spatial_sampling'][:])                # spatial sampling 




#%%

#===========================================
#
#        COMPUTE TIME SHIFTS
#
#===========================================

#
# --- set the parameters for time shift computation 
#

fmin      = 2.9548229548229585
fmax      = 15.45787545787546
fs        = 200.
nfreqs    = 257
ns        = 10
nt        = .5
vpo       = 12
threshold = [0.01,0.]

#
# --- choose the distance 
#
idist     = 90


#
# --- set the parameters for time windowing 
#
vels      = [150,800]
d         = distances[idist]
tmin,tmax = d/vels[1],d/vels[0]

#
# --- set reference and the current waveforms
#


ref,curr=xcorr_premanip[idist,:], xcorr_manip[idist,:]

#
# --- compute ttime shifts 
#

tic()

DelaysXWT = TimeShiftXWT(ref, curr, lags, tmin, tmax, fmin, fmax, nfreqs,ns,nt,vpo,fs,threshold)
DelaysXWT.process()
toc()


#
# --- extract variables 
#



WXamp,  Wcoh, WXdt = DelaysXWT.xwt.amplitude ,DelaysXWT.xwt.coherence ,DelaysXWT.xwt.time_shift 
weights            = DelaysXWT.xwt.weights
dt_mean,dt_std     = DelaysXWT.ftime_shift.mean,DelaysXWT.ftime_shift.std

#%%

#===========================================
#
#            PLOT 
#
#===========================================


# parameters 
time     = lags[DelaysXWT.time_mask]
freqs    = DelaysXWT.freqs
dvv_cmap = create_dvv_colormap(n=100).reversed()


# configure figure 
nrow, ncol = 4, 3
figsize    = (8.27, 11.27)
figsize    = (9, 9)
fig        = plt.figure(figsize=figsize)
gs         = gridspec.GridSpec(nrow,ncol,figure=fig)


# configure axis  
ax_waveforms       = plt.subplot(gs[0,:]) 
ax_time_shift      = plt.subplot(gs[1,0])  
ax_amplitude       = plt.subplot(gs[1,1]) 
ax_weights         = plt.subplot(gs[1,2]) 
ax_time_shift_mean = plt.subplot(gs[2:,:]) 


#

tic()
DelaysXWT.plot_waveforms(ax_waveforms,xlim=(0,2))
DelaysXWT.plot_time_shift(fig,ax_time_shift,xlim=(tmin,tmax),cmap=dvv_cmap)
DelaysXWT.plot_amplitude(fig,ax_amplitude,xlim=(tmin,tmax))
DelaysXWT.plot_weights(fig,ax_weights,xlim=(tmin,tmax))
DelaysXWT.plot_time_shift_mean(ax_time_shift_mean, xlim=(fmin,fmax),ylim=(-0.01,0.06))

toc()

fig.tight_layout()

figpath_egu24="/home/dnziengui/Desktop/codes/general/conferences/EGU_2024/figures"
egu24=False 


