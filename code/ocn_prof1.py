#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 2025
ocean_prof1.py

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
months = np.arange(1,12+1)

cm_path_nh = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/copernicus-data2/'
cm_path_sh = '/Users/mundi/Desktop/antarctica/copernicus-data/'

ice_lims=[20,80]

#%% data loop

ocn_data_all = []
ocn_data_miz = []

for loc_ind, loc in enumerate(hemi_names):
    cm_path = cm_path_nh if loc_ind==0 else cm_path_sh
    path1 = root_paths[loc_ind]

    ### LOAD DEPTH
    with xr.open_dataset(cm_path+'depth.nc') as dds:
        DEPTH = dds['depth'].values
        
    data_all = {mm: [] for mm in np.arange(1,12+1)}
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
                prof_all = np.load(cm_path+start_str+'_'+str(sn)+'_all.npy')
                data_all[month].append(np.nanmean(prof_all, axis=0))
                # data_all[month].append(prof_all[-1] - prof_all[0])
                
                prof_miz = np.load(cm_path+start_str+'_'+str(sn)+'_miz.npy')
                data_miz[month].append(np.nanmean(prof_miz, axis=0))
                # data_miz[month].append(prof_miz[-1] - prof_miz[0])
                
            except FileNotFoundError:
                print('missing:', start_str)
                continue
            
                
    ocn_data_all.append(data_all)
    ocn_data_miz.append(data_miz)
    
#%% plot all profiles
def make_grid_plot(nrow, ncol, title=''):
    fontsize=12
    fig, axes_all = plt.subplots(nrow, ncol, sharex=True, sharey=True, figsize=(3*nrow,5*ncol))
    fig.suptitle('\n\n'+title, fontsize=fontsize+2)
    alph = iter(list(string.ascii_lowercase))
    
    for ax in axes_all[:,0]:
        ax.set_ylabel('Depth (m)', fontsize=fontsize)
        ax.set_ylim([0.5,80])
        ax.invert_yaxis()
        # ax.yaxis.tick_right()
    for ax in axes_all[-1,:]:
        ax.set_xlabel('Temperature'+r' ($^\circ$C)', fontsize=fontsize-1)
    for axl in axes_all:
        for ax1 in axl:
            # ax1.set_xlim([-2,1])
            # ax1.set_xlim([-0.25, 0.25])
            # ax1.axvline(0, ls='-', color='k', lw=0.55)
            ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
            
            ax1.text(0.0225, 0.915, next(alph), transform=ax1.transAxes, fontsize=fontsize, 
                    bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                          'edgecolor':'k', 'lw':0.75},zorder=50)
    return fig, axes_all


for loc_ind, loc in enumerate(hemi_names):
    
    data_all = ocn_data_all[loc_ind]
    data_miz = ocn_data_miz[loc_ind]
    
    yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
    
    for data, title, color in zip([data_all, data_miz], [': Full Storm Area ',': MIZ '], ['#ca0020','#0571b0']):
    
        fig, axes_all = make_grid_plot(4,3, title=loc+title+yr_title)
        for mi, ax in enumerate(axes_all.flatten()): 
            ax.set_title(calendar.month_name[mi+1])
            
            data_mm = data[mi+1]
            
            for line in data_mm:
                ax.plot(line, DEPTH, lw=0.55, color=color)
                
            # if len(data[mi+1])>3:###!!!
            #     ax.plot(np.nanmean(data_mm, axis=0), DEPTH, lw=1.5, color='navy')
                
#%%% mean + spread
for loc_ind, loc in enumerate(hemi_names):
    
    data_all = ocn_data_all[loc_ind]
    data_miz = ocn_data_miz[loc_ind]

    yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
    
    fig, axes_all = make_grid_plot(4,3, title=loc+' Monthly Mean '+yr_title)
    
    for data, title, color in zip([data_all, data_miz], 
                                  ['Full Storm Area ','MIZ '], 
                                  ['#ca0020','#0571b0']):
        for mi, ax in enumerate(axes_all.flatten()): 
            ax.set_title(calendar.month_name[mi+1])
            
            data_mm = data[mi+1]
            
            if len(data_mm)==0: continue
            
            mean_line = np.nanmean(data_mm, axis=0)
            std_line = np.nanstd(data_mm, axis=0)
            
            ax.plot(mean_line, DEPTH, lw=2, 
                    color=color, ls='-', marker='o')
            
            ax.plot(mean_line+std_line, DEPTH, lw=0.5, color=color, ls='-')
            ax.plot(mean_line-std_line, DEPTH, lw=0.5, color=color, ls='-')
            
            ax.fill_betweenx(DEPTH, mean_line, mean_line+std_line, 
                             alpha=0.275, color=color, zorder=-5)
            ax.fill_betweenx(DEPTH, mean_line,mean_line-std_line, 
                            alpha=0.275, color=color, zorder=-5)
 

#%%% mean month grouping

for loc_ind, loc in enumerate(hemi_names[0:1]):
    
    data_all = ocn_data_all[loc_ind]
    data_miz = ocn_data_miz[loc_ind]

    month_colors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
              '#980043','#7a0177','#253494','#2c7fb8','#41b6c4','#a1dab4']
    
    # mgroups = [[12,1,2,3],[4,5,6,7,8],[9,10,11]]     
    # mcolors = [month_colors[0], month_colors[5], month_colors[9]]
  
    # mgroups = [[9,10,11,12,1,2,3],[4,5,6,7,8]]     
    # mcolors = [month_colors[9], month_colors[5], month_colors[0]]
    
    mgroups = [[12,1,2,3,4],[5,6,7], [8,9],[10,11]]     
    mcolors = [month_colors[0], month_colors[5], month_colors[7], month_colors[9]]
    
    fontsize=12
    fig, axes_all = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 8))
    fig.suptitle('\n'+loc+' '+yr_title, fontsize=fontsize+2)
    alph = iter(list(string.ascii_lowercase))
    for ax in axes_all:
        # ax.set_xlim([-2,2])
        ax.set_ylim([0.5,78])
        
        ax.invert_yaxis()
        ax.set_ylabel('Depth (m)', fontsize=fontsize)
        ax.set_xlabel('Temperature'+r' ($^\circ$C)', fontsize=fontsize-1)
        
        ax.tick_params(axis='both', which='major', labelsize=fontsize)  
        ax.text(0.0225, 0.966, next(alph), transform=ax.transAxes, fontsize=fontsize, 
                bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                      'edgecolor':'k', 'lw':0.75},zorder=50)
    
    
    for ax, data, title in zip(axes_all, [data_all, data_miz], ['Full Storm Area','MIZ']):
        ax.set_title(title)
        
        for month_group, mcolor in zip(mgroups,mcolors):
            line_list = []
            for mm in month_group:
                for line in data[mm]:
                    line_list.append(line)
                    
            if len(line_list)==0: continue
            mean_line = np.nanmean(line_list, axis=0)
            std_line = np.nanstd(line_list, axis=0)
            
            ax.plot(mean_line, DEPTH, lw=2, 
                    color=mcolor, ls='-', marker='o')
            
            ax.plot(mean_line+std_line, DEPTH, lw=0.5, color=mcolor, ls='-')
            ax.plot(mean_line-std_line, DEPTH, lw=0.5, color=mcolor, ls='-')
            
            ax.fill_betweenx(DEPTH, mean_line, mean_line+std_line, 
                             alpha=0.275, color=mcolor, zorder=-5)
            ax.fill_betweenx(DEPTH, mean_line,mean_line-std_line, 
                            alpha=0.275, color=mcolor, zorder=-5)
            
            ax.plot([],[], color=mcolor, label=str(month_group)+'  '+str(round(mean_line[-1]-mean_line[0], 3)))
            
        ax.legend(handlelength=1, handletextpad=0.33, loc='lower left',
                  bbox_to_anchor=(0,-0.2))
        
#%%%% individual month lines

for loc_ind, loc in enumerate(hemi_names):
    
    data_all = ocn_data_all[loc_ind]
    data_miz = ocn_data_miz[loc_ind]

    month_colors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
              '#980043','#7a0177','#253494','#2c7fb8','#41b6c4','#a1dab4']
    
    mgroups = [[mm] for mm in months]
    mcolors = month_colors
    
    fontsize=12
    fig, axes_all = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 8))
    fig.suptitle('\n'+loc+' '+yr_title, fontsize=fontsize+2)
    alph = iter(list(string.ascii_lowercase))
    for ax in axes_all:
        # ax.set_xlim([-2,2])
        ax.set_ylim([0.5,78])
        
        ax.invert_yaxis()
        ax.set_ylabel('Depth (m)', fontsize=fontsize)
        ax.set_xlabel('Temperature'+r' ($^\circ$C)', fontsize=fontsize-1)
        
        ax.tick_params(axis='both', which='major', labelsize=fontsize)  
        ax.text(0.0225, 0.966, next(alph), transform=ax.transAxes, fontsize=fontsize, 
                bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                      'edgecolor':'k', 'lw':0.75},zorder=50)
    
    
    for ax, data, title in zip(axes_all, [data_all, data_miz], ['Full Storm Area','MIZ']):
        ax.set_title(title)
        
        for month_group, mcolor in zip(mgroups,mcolors):
            line_list = []
            for mm in month_group:
                for line in data[mm]:
                    line_list.append(line)
                    
            if len(line_list)==0: continue
            mean_line = np.nanmean(line_list, axis=0)
            std_line = np.nanstd(line_list, axis=0)
            
            ax.plot(mean_line, DEPTH, lw=2, 
                    color=mcolor, ls='-')
            
            ax.plot([],[], color=mcolor, label=str(month_group)+'  '+str(round(mean_line[-1]-mean_line[0], 3)))
            
        ax.legend(handlelength=1, handletextpad=0.33, loc='lower left',
                  bbox_to_anchor=(0,-0.33), ncol=2)
                            

#%% scatter
from scipy.stats import linregress

# tbh neither of these are that good (and they are very similar)

#%%% avg "slope" vs ice impact

fig, axes = plt.subplots(2,2, figsize=(8,8))
for ax in axes.flatten():
    ax.axvline(0, lw=0.55, color='gray', ls=':')
    ax.axhline(0, lw=0.55, color='gray', ls=':')
mcolors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
          '#980043','#7a0177','#253494','#2c7fb8','#41b6c4','#a1dab4']

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    data_all = ocn_data_all[loc_ind]
    data_miz = ocn_data_miz[loc_ind]
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
    
    for idx, data, title, color in zip([0,1], [data_all, data_miz], 
                                       ['Full Storm Area ','MIZ '], ['#ca0020','#0571b0']):
        axes[loc_ind][idx].set_title(loc+': '+title, color=color)
        xs, ys = [],[]
        for mi, mm in enumerate(months): 
            data_mm = data[mm]
            si_change = si_changes[mm]
            clim_change = clim_changes[mm]
            
            for line, sia, clim in zip(data_mm, si_change, clim_change):
                del_clim = clim[14] - clim[7]
                del_sia = sia[14] - sia[7]
                del_sia = del_sia - del_clim
                del_line = np.nanmean(line[-5:]) - np.nanmean(line[:8])
                
                xs.append(del_line)
                ys.append(del_sia/1e6)
                axes[loc_ind][idx].plot(del_line, del_sia/1e6, marker='o', markersize=5,
                                      color = mcolors[mi])
                
        mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
        m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
        x = np.linspace(np.nanmin(xs), np.nanmax(xs), 50)
        axes[loc_ind][idx].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
        
        ax6a = axes[loc_ind][idx].twinx()
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');
                
for era in [0,1]: axes[-1][era].set_xlabel('Change in Ocean Temperature (C)')
for loc_ind in [0,1]: axes[loc_ind][0].set_ylabel('Relative '+r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)')

#%%% differnce in first day bottom, last day top vs. ice impact

fig, axes = plt.subplots(2,2, figsize=(8,8))
mcolors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
          '#980043','#7a0177','#253494','#2c7fb8','#41b6c4','#a1dab4']

#### get data

ocn_profs = []
for loc_ind, loc in enumerate(hemi_names):
    cm_path = cm_path_nh if loc_ind==0 else cm_path_sh
    path1 = root_paths[loc_ind]

    profs = {'all':{mm:[] for mm in months}, 
             'miz':{mm:[] for mm in months}}
    for year in years:
        census_file = path1+'census/'+'census_'+str(year)+'.csv'
        [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
        ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
        ice_area80 = ds_area['ice_area80'].values
        ice_area15 = ds_area['ice_area15'].values
        box_area = ds_area['box_area'].values
        ice_sorter = ice_area80
        
        for sn, start in enumerate(startdate): 
            month = start.month
            start_str = start.strftime('%Y-%m-%d')
            ice_frac = ice_sorter[sn]*100/box_area[sn]             
            if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                continue
            
            ## LOAD FILES
            try:
                prof_all = np.load(cm_path+start_str+'_'+str(sn)+'_all.npy')
                profs['all'][month].append([prof_all[0], prof_all[1]])
                
                prof_miz = np.load(cm_path+start_str+'_'+str(sn)+'_miz.npy')
                profs['miz'][month].append([prof_miz[0], prof_miz[1]])
            except: continue
    ocn_profs.append(profs)

#### plot data
for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    profs = ocn_profs[loc_ind]
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
    
    for idx, data, title, color in zip([0,1], [profs['all'], profs['miz']], 
                                       ['Full Storm Area ','MIZ '], ['#ca0020','#0571b0']):
        axes[loc_ind][idx].set_title(loc+': '+title, color=color)
        for mi, mm in enumerate(months): 
            data_mm = data[mm]
            si_change = si_changes[mm]
            clim_change = clim_changes[mm]
            
            for lines, sia, clim in zip(data_mm, si_change, clim_change):
                del_clim = clim[14] - clim[7]
                del_sia = (sia[14] - sia[7])
                del_sia = del_sia - del_clim
                del_line = np.nanmean(lines[-1][-5:]) - np.nanmean(lines[0][:8])
                
                axes[loc_ind][idx].plot(del_line, del_sia/1e6, marker='o', markersize=5,
                                      color = mcolors[mi])
                
for era in [0,1]: axes[-1][era].set_xlabel('Change in Ocean Temperature (C)')
for loc_ind in [0,1]: axes[loc_ind][0].set_ylabel('Relative '+r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)')

#%% change in profile over storm

ocn_data_miz = []

for loc_ind, loc in enumerate(hemi_names):
    cm_path = cm_path_nh if loc_ind==0 else cm_path_sh
    path1 = root_paths[loc_ind]

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
                
                prof_miz = np.load(cm_path+start_str+'_'+str(sn)+'_miz.npy')
                # data_miz[month].append(np.nanmean(prof_miz, axis=0))
                data_miz[month].append(prof_miz[-1] - prof_miz[0])
                
            except FileNotFoundError:
                print('missing:', start_str)
                continue
            
                
    ocn_data_all.append(data_all)
    ocn_data_miz.append(data_miz)


#%%% profile plots
for loc_ind, loc in enumerate(hemi_names):
    
    data_all = ocn_data_all[loc_ind]
    data_miz = ocn_data_miz[loc_ind]
    
    yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
    
    for data, title, color in zip([data_miz], [': MIZ '], ['#0571b0']):
    
        fig, axes_all = make_grid_plot(4,3, title=loc+title+yr_title+'\nChange Over Storm Duration')
        for mi, ax in enumerate(axes_all.flatten()): 
            ax.set_title(calendar.month_name[mi+1])
            ax.set_xlim([-0.75, 0.75])
            
            data_mm = data[mi+1]
            
            for line in data_mm:
                ax.plot(line, DEPTH, lw=0.55, color=color)
                
            # if len(data[mi+1])>3:###!!!
            #     ax.plot(np.nanmean(data_mm, axis=0), DEPTH, lw=1.5, color='navy')

#%%% mean + spread
for loc_ind, loc in enumerate(hemi_names):
    
    data_all = ocn_data_all[loc_ind]
    data_miz = ocn_data_miz[loc_ind]

    yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
    
    fig, axes_all = make_grid_plot(4,3, title=loc+' Change Over Storm Duration '+yr_title)
    
    for data, title, color in zip([data_miz], ['MIZ '], ['#0571b0']):
        for mi, ax in enumerate(axes_all.flatten()): 
            ax.set_title(calendar.month_name[mi+1])
            ax.set_xlim([-0.66, 0.66])
            ax.axvline(0, lw=0.55, color='k')
            
            data_mm = data[mi+1]
            
            if len(data_mm)<10: continue
            
            mean_line = np.nanmean(data_mm, axis=0)
            std_line = np.nanstd(data_mm, axis=0)
            
            ax.plot(mean_line, DEPTH, lw=2, 
                    color=color, ls='-', marker='o')
            
            ax.plot(mean_line+std_line, DEPTH, lw=0.5, color=color, ls='-')
            ax.plot(mean_line-std_line, DEPTH, lw=0.5, color=color, ls='-')
            
            ax.fill_betweenx(DEPTH, mean_line, mean_line+std_line, 
                             alpha=0.275, color=color, zorder=-5)
            ax.fill_betweenx(DEPTH, mean_line,mean_line-std_line, 
                            alpha=0.275, color=color, zorder=-5)

#%%% monthly means
mcolors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
          '#980043','#7a0177','#253494','#2c7fb8','#41b6c4','#a1dab4']
mgroups = [[mm] for mm in np.arange(1,12+1)]

#### set up plot
fontsize=12
fig, axes_all = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 8))
fig.suptitle('\nChange in Temperation Profiles within the MIZ '+yr_title, fontsize=fontsize+2)
alph = iter(list(string.ascii_lowercase))

for ax in axes_all:
    ax.set_ylim([0.5,78])
    ax.axvline(0, lw=2, color='k', ls='-')
    
    ax.invert_yaxis()
    ax.set_ylabel('Depth (m)', fontsize=fontsize)
    ax.set_xlabel('Change in Temperature'+'\nOver Storm Duration'+r' ($^\circ$C)', 
                  fontsize=fontsize-1)
    
    ax.tick_params(axis='both', which='major', labelsize=fontsize)  
    ax.text(0.0225, 0.966, next(alph), transform=ax.transAxes, fontsize=fontsize, 
            bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                  'edgecolor':'k', 'lw':0.75},zorder=50)

#### data loop
month_mean_diffs = {0:[],1:[]}
for loc_ind, loc in enumerate(hemi_names):
    
    data_miz = ocn_data_miz[loc_ind]

    ax = axes_all[loc_ind]
    ax.set_title(loc)
    
    for month_group, mcolor in zip(mgroups,mcolors):
        line_list = []
        for mm in month_group:
            for line in data_miz[mm]:
                line_list.append(line)
                
        if len(line_list)<10: month_mean_diffs[loc_ind].append([np.nan]*len(DEPTH));continue
        mean_line = np.nanmean(line_list, axis=0)
        std_line = np.nanstd(line_list, axis=0)
        
        ax.plot(mean_line, DEPTH, lw=2, color=mcolor, ls='-',label=str(month_group[0]))
        
        month_mean_diffs[loc_ind].append(mean_line)
        
    #     ax.plot([],[], color=mcolor, label=str(month_group)+'  '+str(round(mean_line[0]-mean_line[-1], 3)))
        
    # ax.legend(handlelength=1, handletextpad=0.33, loc='lower left',
    #           bbox_to_anchor=(0,-0.33), ncol=2)


axes_all[0].legend(handlelength=1.25, handletextpad=0.33, loc='lower left',
          bbox_to_anchor=(-0.2,-0.2), ncol=12, title='Months')

#%%% scatter: upper ocean mean change with ice changes

fig, axes = plt.subplots(1,2, figsize=(8,4))
for ax in axes.flatten():
    ax.axvline(0, lw=0.55, color='gray', ls=':')
    ax.axhline(0, lw=0.55, color='gray', ls=':')
mcolors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
          '#980043','#7a0177','#253494','#2c7fb8','#41b6c4','#a1dab4']

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
    
    axes[loc_ind].set_title(loc+': MIZ')
    xs, ys = [],[]
    for mi, mm in enumerate(months): 
        impact = mean_lines[mi][-1]
        ocn_prof = month_mean_diffs[loc_ind][mi]
        ocn_diff = np.nanmean(ocn_prof[0:8]) #8=upper 10m
        
        xs.append(ocn_diff)
        ys.append(impact)
        axes[loc_ind].plot(ocn_diff, impact, marker='o', markersize=8,
                              color = mcolors[mi])
                
    mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
    m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
    x = np.linspace(np.nanmin(xs), np.nanmax(xs), 50)
    axes[loc_ind].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
    
    ax6a = axes[loc_ind].twinx()
    ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
    ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
    ax6a.axis('off');
                
for ax in axes: ax.set_xlabel('Change in Ocean Temperature (C)')
axes[0].set_ylabel('Normalized Change in MIZ Ice Area')
#'Relative '+r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)'





                 
#%% end             
                