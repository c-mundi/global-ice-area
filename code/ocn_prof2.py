#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 2025
ocn_prof2.py

expanded analysis time (storm duration, 1 week, 2 weeks)

@author: mundi
"""
#%% import and fpaths​
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
import cmocean.tools as cmo_tools
import calendar
from datetime import datetime, timedelta
import time as timeIN
import string

import functions as fx

import pandas as pd
from scipy.interpolate import griddata
import warnings

ice_fname =  '/Users/mundi/Desktop/seaice/' #'south/'

root = '/Users/mundi/Desktop/month-hemi/'
root_paths = [root+'nh_data/', root+'sh_data/']

hemi_names = ['Arctic', 'Antarctic']

years = np.arange(2010,2020) 
yr_title = str(years[0])+'-'+str(years[-1])
months = np.arange(1,12+1)

mcolors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
          '#980043','#7a0177','#253494','#2c7fb8','#41b6c4','#a1dab4']

ice_lims=[20,80]

#%% data loop

ocn_data_miz = []

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    cm_path = path1+'ocn_prof/'

    ### LOAD DEPTH
    with xr.open_dataset(cm_path+'depth.nc') as dds:
        DEPTH = dds['depth'].values
        
    data_miz = {mm: [] for mm in np.arange(1,12+1)}
    
    ### STORM LOOP
    for year in years:
        # CENSUS: storm info
        census_file = path1+'census/'+'census_'+str(year)+'.csv'
        [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
        
        # open ice area
        ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
        ice_area80 = ds_area['ice_area80'].values
        ice_area15 = ds_area['ice_area15'].values
        box_area = ds_area['box_area'].values
        ice_sorter = ice_area80
        
        for sn, start in enumerate(startdate): 
            month = start.month
            
            ice_frac = ice_sorter[sn]*100/box_area[sn]             
            if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                continue
    
            ## LOAD FILES
            start_str = start.strftime('%Y-%m-%d')
            
            try:
                prof_miz = np.load(cm_path+start_str+'_'+str(sn)+'_miz_series.npy')
                data_miz[month].append(prof_miz)
                
            except FileNotFoundError:
                print('missing:', start_str)
                continue
            
    ocn_data_miz.append(data_miz)
    
#%% monthly mean profiles

for di, diff in enumerate(['1-Week', '2-Week','3-week']):

    #### set up plot
    fontsize=12
    fig, axes_all = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 8))
    fig.suptitle('Change in Temperature Profiles within the MIZ '+yr_title, fontsize=fontsize+2)
    fig.text(0.5, 0.925, diff+' Difference', fontsize=fontsize+1, fontstyle='italic',
             horizontalalignment='center', verticalalignment='center')
    alph = iter(list(string.ascii_lowercase))
    
    for ax in axes_all:
        ax.set_ylim([0.5,78])
        ax.axvline(0, lw=1, color='gray', ls=':')
        
        ax.invert_yaxis()
        ax.set_ylabel('Depth (m)', fontsize=fontsize)
        ax.set_xlabel('Change in Temperature'+r' ($^\circ$C)', 
                      fontsize=fontsize-1)
        
        ax.tick_params(axis='both', which='major', labelsize=fontsize)  
        ax.text(0.0225, 0.966, next(alph), transform=ax.transAxes, fontsize=fontsize, 
                bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                      'edgecolor':'k', 'lw':0.75},zorder=50)
        
        ax.set_xlim([-1.25, 1.25])
    
    #### data loop
    month_mean_diffs = {0:[],1:[]}
    for loc_ind, loc in enumerate(hemi_names):
        
        data_miz = ocn_data_miz[loc_ind]
    
        ax = axes_all[loc_ind]
        ax.set_title(loc)
        
        for month, mcolor in zip(months,mcolors):
            if di==0: line_list = [prof[14]-prof[7] for prof in data_miz[month]]
            elif di==1: line_list = [prof[-1]-prof[7] for prof in data_miz[month]]
            elif di==2: line_list = [prof[-1]-prof[0] for prof in data_miz[month]]
                    
            if len(line_list)<10: month_mean_diffs[loc_ind].append([np.nan]*len(DEPTH));continue
            mean_line = np.nanmean(line_list, axis=0)
            std_line = np.nanstd(line_list, axis=0)
            
            ax.plot(mean_line, DEPTH, lw=2, color=mcolor, ls='-',label=str(month))
            
            month_mean_diffs[loc_ind].append(mean_line)
            
    axes_all[0].legend(handlelength=1.25, handletextpad=0.33, loc='lower left',
              bbox_to_anchor=(-0.2,-0.2), ncol=12, title='Months') 
    

#%% end
    