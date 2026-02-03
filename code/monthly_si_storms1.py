#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 2025
monthly_si_storms1.py

-> do months with more storms have greater MIZ ice area
- 3 month window / 1 month lag

@author: mundi
"""
#%% imports and filepaths
import numpy as np
import matplotlib.pyplot as plt
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress

import functions as fx

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

census_name = 'census_'
contour_name = '_contours'
si_name = '_seaice'

ice_fname = '/Users/mundi/Desktop/seaice/'

decades = [np.arange(1982, 1992), np.arange(2010,2020)]
decade_names = ['Early Satellite Era ', 'Present Day ']

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]

hemi_names= ['Arctic', 'Antarctic']

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

alph = ['a','b','c','d','e','f']

#%% organize data and plot timeseries
'''seaice_impact1.py'''

savepath = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/monthly_miz/'
_, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
_, si_lon_sh, si_lat_sh = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)

area = False  #T = area containing MIZ ice, F = actual ice area (*concentration)

if area: fname = '_area'
else: fname = '_miz'

#### yearly timeseries plots
fig2, axes2 = plt.subplots(2,2, figsize=(20,10),sharey=True)
fig2.suptitle('\nMonthly MIZ Ice Area and Storm Counts', fontweight='bold', fontsize=16)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=14)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=14, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Year',fontsize=14)
    
fig2.subplots_adjust(wspace=0.1)


#### data
miz_all = []
counts_all = []
for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    
    miz_data = {}
    storm_counts = {}
    for yi, years in enumerate(decades):
        axes2[loc_ind][yi].set_title(loc+', '+str(years[0])+'-'+str(years[-1]))
        ax2 = axes2[loc_ind][yi].twinx()
        for year in years:
            # get miz ice area data
            miz_data[year] = {}
            if loc_ind == 0: daily_miz = np.load(savepath+str(year)+fname+'.npy') 
            elif loc_ind == 1: daily_miz = np.load(savepath+str(year)+fname+'_sh.npy')
            
            # get storm counts
            storm_counts[year] = []
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines([year], path1+'census/', path1+'area/', path1+'seaice/')
            
            # separate into months
            mind=0
            for month in months:
                num_days = calendar.monthrange(2010, month)[1]
                mslice = daily_miz[mind:mind+num_days]
                miz_data[year][month] = mslice
                mind+=num_days
                
                storm_counts[year].append( len(start_day[month]) )
            
            #### line plots: full miz ice area/storm count timeseries
            axes2[loc_ind][yi].plot(year+np.linspace(0,1,12), storm_counts[year])
            ax2.plot(year+np.linspace(0,1,12), 
                     [np.nanmean(miz_data[year][month]) for month in months], 
                     color='k')
            
    miz_all.append(miz_data)
    counts_all.append(storm_counts)
    
    
ax2.plot([],[], color='b', label='Yearly storm counts')
ax2.plot([],[], color='k', label='MIZ area')
ax2.legend(fontsize=16)
    
    
#%% grouped-month window...
mwindow = 5

#### set up plot
fig2, axes2 = plt.subplots(2,2, figsize=(20,10),sharey=True)
fig2.suptitle('\nMonthly MIZ Ice Area and Storm Counts, Window='+str(mwindow), 
              fontweight='bold', fontsize=16)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=14)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=14, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Year',fontsize=14)
fig2.subplots_adjust(wspace=0.1)

#### data loop
for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    
    miz_data = miz_all[loc_ind]
    storm_counts = counts_all[loc_ind]

    for yi, years in enumerate(decades):
        # big decade list, then group after loops
        miz_list = []
        count_list = []
        for year in years:
            for month in months:
                daily_miz = list(miz_data[year][month])
                while len(daily_miz)<31:
                    daily_miz.append(np.nan)
                month_counts = storm_counts[year][month-1]
                
                miz_list.append(np.array(daily_miz))
                count_list.append(month_counts)
                
        # group months
        miz_values = []
        count_values = []
        for idx in range(len(miz_list[:-3])):
            mean_miz = np.nanmean(miz_list[idx:idx+mwindow])
            total_count = np.nansum(count_list[idx:idx+mwindow])
            
            miz_values.append(mean_miz)
            count_values.append(total_count)
            
        # plot
        axes2[loc_ind][yi].plot(count_values, label='Storm counts')
        ax2 = axes2[loc_ind][yi].twinx()
        ax2.plot(miz_values, color='k', label='MIZ area')
  
axes2[loc_ind][yi].legend(fontsize=16)              
ax2.legend(fontsize=16)
        
#%%% seasonality

mgroups = [10,11]

#### set up plot
fig2, axes2 = plt.subplots(2,2, figsize=(20,10),sharey=True)
fig2.suptitle('\nMonthly MIZ Ice Area and Storm Counts, Window='+str(mwindow)+', Months='+str(mgroups), 
              fontweight='bold', fontsize=16)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=14)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=14, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Year',fontsize=16)
fig2.subplots_adjust(wspace=0.1)

#### data loop
for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    
    miz_data = miz_all[loc_ind]
    storm_counts = counts_all[loc_ind]

    for yi, years in enumerate(decades):
        # big decade list, then group after loops
        miz_list = []
        count_list = []
        for year in years:
            mean_miz = [np.nanmean(miz_data[year][month]) for month in mgroups]
            month_counts = [storm_counts[year][month-1] for month in mgroups]
                
            miz_list.append(np.array(mean_miz))
            count_list.append(month_counts)
                
        # group months
        miz_values = []
        count_values = []
        for idx in range(len(miz_list)):
            mean_miz = np.nanmean(miz_list[idx:idx+mwindow])
            total_count = np.nansum(count_list[idx:idx+mwindow])
            
            miz_values.append(mean_miz)
            count_values.append(total_count)
            
        # plot
        axes2[loc_ind][yi].plot(count_values, marker='o')
        ax2 = axes2[loc_ind][yi].twinx()
        ax2.plot(miz_values, color='k', marker='o')

        ax2.set_title(loc+', '+str(years[0])+'-'+str(years[-1]),fontsize=16)

#%% end
