#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 2026
major_revisions1.py

response to reviewers

@author: mundi
"""

#%% file paths and imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress
import time as timeIN
from glob import glob
import warnings

from matplotlib.gridspec import GridSpec
from osgeo import gdal
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata

import functions as fx

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

census_name = 'census_'
contour_name = '_contours'
si_name = '_seaice'

ice_fname = '/Users/mundi/Desktop/seaice/south/'

decades = [np.arange(1982, 1992), np.arange(2010,2020)]
decade_names = ['Early Satellite Era ', 'Present Day ']

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]

hemi_names= ['Arctic', 'Antarctic']

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

fontsize=11
fs = 14

loc_ind=1
path1= root_paths[loc_ind]

__, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010,1,1)

alph = ['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k']
month_colors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
                '#980043','#7a0177', '#253494','#2c7fb8','#41b6c4','#a1dab4']


#%%% bar fxn

def seaice_lines(years, census_path, area_path, si_path, AREA_IND, sample_size=10, p_thresh=984, si_name='_seaice'):
    ice_lims = [20,80]
    
    # addon_num = 0 # storm area
    sia_var_addons = [['_1000', '_990', '_contour'],['2_1000', '2_990', '2_contour']]
    
    miz_ind = 1 #[0=daily, 1=total]
    
    # data dict
    lines = {}
    start_day = {}
    end_day = {}
    si_changes = {} 
    clim_changes = {}
    for idx in np.arange(0,12): # months
        lines[idx+1] = []
        start_day[idx+1] = []
        end_day[idx+1] = []
        si_changes[idx+1] = []
        clim_changes[idx+1] = []
    
    # data loop
    for year in years:
        
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
        pressure = fx.readCensus(census_file, convertDT=True)[-1]
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in zip(startdate, enddate):
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(fx.daterange(week_ago, two_week, dt=24))
            storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
        try:
            ds = xr.open_dataset(si_path + str(year) + si_name + '.nc')
        except:
            print('- skip: '+si_path + str(year) + si_name + '.nc')
            continue
        
        for storm_num, storm in enumerate(storm_ranges):
             month = int(storm[0].month)
             
             ### remove storms that don't interact with the ice
             with warnings.catch_warnings():
                 warnings.simplefilter('ignore')
                 ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
             if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                continue
             try:
                 sia = ds['sia_miz'+sia_var_addons[miz_ind][AREA_IND]].values[storm_num]
             except (KeyError, IndexError):
                 print()
                 print('***', year, len(storm_ranges), 
                       np.shape(ds['sia_miz'+sia_var_addons[miz_ind][AREA_IND]].values))
                 print()
                 break
             if pressure[storm_num] > p_thresh: continue
        
             # get time series   
             try:
                 sia_clim = ds['sia_clim_miz'+sia_var_addons[miz_ind][AREA_IND]].values[storm_num]
             except KeyError:
                 sia_clim = ds['si_clim'].values[storm_num]
             ss = sia-sia_clim
             with warnings.catch_warnings():
                 warnings.simplefilter('ignore')
                 standardized_area = (ss-ss[0])/(np.nanmax(ss)-np.nanmin(ss))
             pcd = (standardized_area)
        
             # plot
             lines[month].append(pcd)
             start_day[month].append(storm[0])
             end_day[month].append(storm[-1])
             
             si_changes[month].append(sia)
             clim_changes[month].append(sia_clim)
        
    mean_lines = []  
    for idx in np.arange(0,12):
        
        if len(lines[idx+1])> sample_size: 
            mean_line = np.nanmean(lines[idx+1], axis=0)
            mean_lines.append(mean_line)
        else:
            mean_lines.append(np.nan*np.ones(np.shape(analysis_ranges[0])))
            continue
        
    return mean_lines, lines, start_day, end_day, si_changes, clim_changes


#%% "MIZ ice area" definition
# are results sensitive to choice of MIZ area?

# storm_tseries2_miz.py -> miz sensitivity runs

#%%% different miz area durations

years = np.arange(2010,2020)
yr_title = str(years[0])+'-'+str(years[-1])

#### set up plot
fig2 = plt.figure(figsize=(27, 14))
gs = GridSpec(4,3, hspace=0.5)

axes2=[]
for li, loc in enumerate(['Arctic', 'Antarctic']):
    ax1 = fig2.add_subplot(gs[2*li:2*li+2,0])
    ax2 = fig2.add_subplot(gs[2*li:2*li+2,1], sharey=ax1)
    axes2.append([ax1,ax2])
axes2=np.array(axes2)

axes_add = []
for li, loc in enumerate(['Arctic', 'Antarctic']):
    ax1 = fig2.add_subplot(gs[2*li,-1])
    ax2 = fig2.add_subplot(gs[2*li+1,-1], sharey=ax1)
    axes_add.append([ax1,ax2])
axes_add=np.array(axes_add)

fig2.suptitle('\n'+yr_title,fontsize=fs+8, fontweight='bold')
for i, ax2 in enumerate(list(axes2.flatten())+list(axes_add.flatten())):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.axvline(0, ls=':', color='gray', lw=1)
    ax2.set_xlim(-7,14)
    ax2.set_xticks(xxx)
    ax2.set_xticklabels(xlabels, minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
for ax2 in axes2[-1,:]: ax2.set_xlabel('Days Since Storm Start',fontsize=fs)
axes_add[-1][-1].set_xlabel('Days Since Storm Start',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Normalized Relative Change in Ice Area',fontsize=fs+1)

    
#### data and plot mean lines

for li, loc in enumerate(['Arctic', 'Antarctic']):
    
    fig2.text(0.0725, li*-0.4+0.70, loc, rotation=90,
             va='center', ha='center', fontsize=fontsize+10, fontweight='bold')

    path1 = root_paths[li]
    
    si_titles = ['3-day Window Around Storm Duration', 'MIZ Area on the Start Day of the Storm']
    for i, si_name in enumerate(['_seaice_3', '_seaice_start']):
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice_miz/', si_name=si_name)
    
        mean_lines_original = fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')[0]
    
        lines_anoms = {}
        for mi, mm in enumerate(months):
            axes2[li][i].plot(xxx, mean_lines[mi], color = month_colors[mi], lw=2, label = calendar.month_name[mm])
    
            #### difference plots
            axes_add[li][i].plot(xxx, mean_lines[mi]-mean_lines_original[mi], color = month_colors[mi], lw=2, label = calendar.month_name[mm])
        
    
        axes2[li][i].set_title(si_titles[i],fontsize=fs+2 )
        axes_add[li][i].set_title(si_titles[i],fontsize=fs+2 )
        if li==0 and i==0:axes_add[li][i].set_title('Difference from Original'+'\n\n'+si_titles[i],fontsize=fs+2 )
        
        #### monthly lines
        # fig, axes = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
        # fig.suptitle(loc+': '+str(years[0])+'-'+str(years[-1])+'\n'+si_titles[i])
        
        # axes = axes.flatten()
        # for ai, ax1 in enumerate(axes): 
        #     ax1.set_title(calendar.month_name[ai+1])
        #     ax1.set_xlim([-7,14])
        #     # ax1.set_ylim([-1.05,1.05])
        #     ax1.axhline(0, lw=0.55, color='k')
        #     ax1.axvline(0, lw=0.55, color='k')
            
        # for mi, mm in enumerate(months):
        #     axes[mi].plot(xxx, np.array(lines[mm]).T, color = 'gray', lw=0.55)
        #     axes[mi].plot(xxx, mean_lines[mi], color = month_colors[mi], lw=2)
        #     axes[mi].plot(xxx, mean_lines_original[mi], ls='--', color = month_colors[mi], lw=2)
            
            
            
    axes_add[0][-1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                      edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                      bbox_to_anchor=(1.33,0.85));
    
    
#%%%% bar plots

#### settings

# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

XSPACING = [0,0.2,0.4]
INDS = [7+3, 7+7, 7+14]
WIDTH = XSPACING[1]-XSPACING[0]
HATCH = '/////'
alphbar = ['a','b','c','d']

# mean cyclone impact vs summed across all storms
plot_mean = False
# normalized impact vs total change in area
norm_lines = False
# months 1->12 vs shifted SH seasonal cycle
shift = False 
# use same trials for min/max for all months vs use each month's min/max
set_minmax = False
set_trials = [3,1]

# errorbar 
capsize=2
eb_color= hemi_colors #['k','k']


#### set up
fig2, axes2 = plt.subplots(2,1, figsize=(10,12), sharex=True, sharey=True)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
    if plot_mean: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    else: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)

axes2[0].set_title('MIZ Ice Area Changes in Each Hemisphere', fontsize=fs+2)
axes2[1].set_title('Global Sum', fontsize=fs+2)

# axes2[0].text(0.05, 41,
#               'Error bars indicate spread of MIZ ice area change when using a 3-day window and \n'+
#               'a fixed area based on the start day of the storm. Bars use the original methodology.', 
#               )

plt.subplots_adjust(hspace=0.275)

import matplotlib.patches as mpatches
circ1 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,
                       label=str(decades[1][0])+'-'+str(decades[1][-1]))
circ2 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,hatch=HATCH,
                       label=str(decades[0][0])+'-'+str(decades[0][-1]))

ax2.legend(handles = [circ2, circ1], loc='lower left',
           fontsize=fs-1, handletextpad=0.5, handlelength=1.5)
    
    
fig2.subplots_adjust(wspace=0.1)

########################
#### data and plot
########################

values = {}
valmin = {}; valmax = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    
    path1 = root_paths[li]
    path2 = root_paths[li]+'sensitivity/'
    
            #name   #data path #seaice path #area-ind #p-thresh # sea ice name
    trials = {'original':[path1, path1+'seaice/', 0, [984, 957], '_seaice'],
             'si3':[path1, path1+'seaice_miz/', 0, [984, 957], '_seaice_3'],
             'si_start':[path1, path1+'seaice_miz/', 0, [984, 957],'_seaice_start'],
             }
    
    axes2[0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    # data loop, bar plot
    era_values = {}
    min_values = {}; max_values = {}
    for yi, years in enumerate([np.arange(1982,1992),np.arange(2010,2020)]):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        if yi == 0: trial = {'original': trials['original']}
        elif yi == 1: trial = trials

        era_values[yi] = [[] for x in np.arange(0, len(trial.keys()))]
        min_values[yi] = []
        max_values[yi] = []
        
        for k, name in enumerate(list(trial.keys())):
            attr = trial[name]
            
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                seaice_lines(years, path1+'census/', path1+'area/', attr[1], 
                             AREA_IND=attr[2], p_thresh=attr[3][li], si_name=attr[4])
                
            for mi, mm in enumerate(months):
                rel_area = []
                for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                    ss = sia-sia_clim
                    rel_area.append(ss-ss[0])
                if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
                else: MEANLINE = np.nansum(rel_area, axis=0)/1e6
            
                if len(lines[mi+1])<10: 
                    era_values[yi][k].append([np.nan]*3)
                    continue 
            
                if not shift: xvals = mi+np.array(XSPACING)
                else:
                    xvals = mi+np.array(XSPACING) + 6
                    if xvals[0]>=12: xvals -= 12
                
                diffs = [MEANLINE[ii] for ii in INDS]
                era_values[yi][k].append(diffs)
                
                # plot bars/ comparison dots
                if name == 'original':
                    if yi==1:
                        axes2[0].bar(xvals, diffs, width=WIDTH,
                                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
                    else:
                        axes2[0].bar(xvals[-1]+WIDTH, diffs[-1], width=WIDTH,
                                facecolor=hemi_colors[li], alpha=0.5, 
                                edgecolor=hemi_colors[li], hatch=HATCH, lw=2)


        #### find min/max values after each trial is calculated
        if yi==1:
            for mi, mm in enumerate(months):
                if not shift: xvals = mi+np.array(XSPACING)
                else: 
                    xvals = mi+np.array(XSPACING) + 6
                    if xvals[0]>=12: xvals -= 12
                
                # plot spread bars
                diff_arr = np.array([era_values[yi][k][mi] for k in np.arange(0, len(trial))])
                if set_minmax:
                    min_arr = diff_arr[set_trials[0]]
                    max_arr = diff_arr[set_trials[1]]
                else:
                    min_arr = np.nanmin(diff_arr, axis=0)
                    max_arr = np.nanmax(diff_arr, axis=0)
                    
                axes2[0].errorbar(xvals, diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                       fmt='none', capsize=capsize, ecolor=eb_color[li])
                min_values[yi].append(min_arr)
                max_values[yi].append(max_arr)
                

        ### TOTAL
        sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(len(xvals))] )
        if yi==1:
            axes2[0].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        else:
            axes2[0].bar(xvals[-1]+WIDTH+1.33+li, sums[-1], width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors[li], alpha=0.5, edgecolor=hemi_colors[li], hatch=HATCH, lw=2)
            
        if yi==1:
            min_sums = np.array( [np.nansum(np.array(min_values[yi])[:,dr]) for dr in range(len(xvals))] )
            max_sums = np.array( [np.nansum(np.array(max_values[yi])[:,dr]) for dr in range(len(xvals))] )
            axes2[0].errorbar(xvals+1.33+li, sums, yerr=[np.abs(min_sums-sums), np.abs(max_sums-sums)],
                                   fmt='none', capsize=capsize, ecolor=eb_color[li])
       
        axes2[0].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
    
#### sum
for yi, years in enumerate(decades):
    sums = []
    for vi, VAL in enumerate([values, valmin, valmax]): # each hemi
        val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
        val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
        sums.append( val0 + val1 )
    
        if vi==0:
            heights = []
            for mi, sum1 in enumerate(sums[0][0]):
                xvals = mi+np.array(XSPACING)
                heights.append(sum1)
                
                    # original
                fcs = [hemi_colors2[0][i] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
                ecs = [hemi_colors[0] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors[1] for i in [0,1,2]]
                
                if yi==1:
                    axes2[1].bar(xvals, sum1, width=WIDTH,
                            facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)
                else:
                    axes2[1].bar(xvals[-1]+WIDTH, sum1[-1], width=WIDTH,
                            facecolor=ecs[-1], alpha=0.5, edgecolor=ecs[-1], hatch=HATCH,lw=2)

            # total sum
            if yi==1:
                axes2[1].bar(xvals+1.33+0.5, np.nansum(sums[0][0], axis=0), width=WIDTH,
                        facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
            else:
                axes2[1].bar(xvals[-1]+WIDTH+1.33+0.5, np.nansum(sums[0][0], axis=0)[-1], width=WIDTH,
                        facecolor='k', alpha=0.5, edgecolor='k', hatch=HATCH, lw=2)
            
            axes2[1].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
            
        # add spread bars!
        elif vi==1 and yi==1:
            val0_min = val0
            val1_min = val1
        elif vi==2 and yi==1:
            for mi in np.arange(0,12):
                xvals = mi+np.array(XSPACING)
                
                axes2[1].errorbar(xvals, heights[mi], 
                                yerr=[np.abs(sums[1][mi]-heights[mi]), np.abs(sums[2][mi]-heights[mi])],
                                fmt='none', capsize=capsize, ecolor='#262626')
                
            # total sum 
            ht1 =  np.nansum(sums[0][0], axis=0)
            axes2[1].errorbar(xvals+1.33+0.5,  ht1, 
                            yerr=[np.abs(np.nansum(sums[1], axis=0)-ht1), np.abs(np.nansum(sums[2], axis=0)-ht1)],
                            fmt='none', capsize=capsize+0.5, ecolor='k')
  

#### plot labels
axes2[0].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
axes2[1].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
for ax in axes2.flatten():
    ax.xaxis.set_tick_params(labelbottom=True)
    ax.yaxis.set_tick_params(labelleft=True)
    
#%%%% compare spreads....

#### settings
# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

XSPACING = [0,0.2,0.4]
INDS = [7+3, 7+7, 7+14]
WIDTH = XSPACING[1]-XSPACING[0]
HATCH = '/////'
alphbar = ['a','b','c','d']

# mean cyclone impact vs summed across all storms
plot_mean = False
# normalized impact vs total change in area
norm_lines = False
# months 1->12 vs shifted SH seasonal cycle
shift = False 
# use same trials for min/max for all months vs use each month's min/max
set_minmax = True
set_trials = [3,1]

# errorbar 
capsize=2
eb_color= hemi_colors #['k','k']
s_size=27


#### set up
fig2, axes2 = plt.subplots(2,1, figsize=(10,12), sharex=True, sharey=True)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
    if plot_mean: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    else: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)
# ax2.set_xlabel('Month',fontsize=fs)

axes2[0].set_title('MIZ Ice Area Changes in Each Hemisphere', fontsize=fs+2)
axes2[1].set_title('Global Sum', fontsize=fs+2)

plt.subplots_adjust(hspace=0.275)

import matplotlib.patches as mpatches
circ1 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,
                       label=str(decades[1][0])+'-'+str(decades[1][-1]))
circ2 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,hatch=HATCH,
                       label=str(decades[0][0])+'-'+str(decades[0][-1]))

ax2.legend(handles = [circ2,circ1], loc='lower left',
           fontsize=fs-1, handletextpad=0.5, handlelength=1.5)
    
    
fig2.subplots_adjust(wspace=0.1)

for li, loc in enumerate(['Arctic','Antarctic']): axes2[0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
axes2[0].scatter([],[], s=s_size, marker='d', color='gray', label='New MIZ Definition')
axes2[0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)

#### data and plot

values = {}
valmin = {}; valmax = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    
    path1 = root_paths[li]
    path2 = root_paths[li]+'sensitivity/'
    
            #name   #data path #seaice path #area-ind #p-thresh # colors
    trials = {'original':[path1, path1+'seaice/', 0, [984, 957], 
                         ['lightcoral','indianred','maroon']],
             'area_990':[path1, path1+'sensitivity_seaice/', 1, [984, 957],
                         ['royalblue','mediumblue','navy']],
             'area_cont':[path1, path1+'sensitivity_seaice/', 2, [984, 957],
                          ['mediumaquamarine','mediumseagreen','green']],
             'p_1000':[path2, path2+'seaice/', 0, [994,967],
                       ['gold','goldenrod','darkgoldenrod']],
             # 'p_cont':[path2, path2+'seaice/', 2, [994,967],
             #           ['violet','mediumorchid','darkviolet']],
             }
    
    # data loop, bar plot
    era_values = {}
    min_values = {}; max_values = {}
    for yi, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        if yi == 0: trial = {'original': trials['original']}
        elif yi == 1: trial = trials

        era_values[yi] = [[] for x in np.arange(0, len(trial.keys()))]
        min_values[yi] = []
        max_values[yi] = []
        
        for k, name in enumerate(list(trial.keys())):
            attr = trial[name]
            
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                seaice_lines(years, attr[0]+'census/', attr[0]+'area/', attr[1], 
                             AREA_IND=attr[2], p_thresh=attr[3][li])

            for mi, mm in enumerate(months):
                rel_area = []
                for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                    ss = sia-sia_clim
                    rel_area.append(ss-ss[0])
                if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
                else: MEANLINE = np.nansum(rel_area, axis=0)/1e6
            
                if len(lines[mi+1])<10: 
                    era_values[yi][k].append([np.nan]*3)
                    continue 
            
                if not shift: xvals = mi+np.array(XSPACING)
                else:
                    xvals = mi+np.array(XSPACING) + 6
                    if xvals[0]>=12: xvals -= 12
                
                diffs = [MEANLINE[ii] for ii in INDS]
                era_values[yi][k].append(diffs)
                
                # plot bars/ comparison dots
                if name == 'original':
                    if yi==1:
                        axes2[0].bar(xvals, diffs, width=WIDTH,
                                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
                    else:
                        axes2[0].bar(xvals[-1]+WIDTH, diffs[-1], width=WIDTH,
                                facecolor=hemi_colors[li], alpha=0.5, 
                                edgecolor=hemi_colors[li], hatch=HATCH, lw=2)


        #### find min/max values after each trial is calculated
        if yi==1:
            for mi, mm in enumerate(months):
                if not shift: xvals = mi+np.array(XSPACING)
                else: 
                    xvals = mi+np.array(XSPACING) + 6
                    if xvals[0]>=12: xvals -= 12
                
                # plot spread bars
                diff_arr = np.array([era_values[yi][k][mi] for k in np.arange(0, len(trial))])
                if set_minmax:
                    min_arr = diff_arr[set_trials[0]]
                    max_arr = diff_arr[set_trials[1]]
                else:
                    min_arr = np.nanmin(diff_arr, axis=0)
                    max_arr = np.nanmax(diff_arr, axis=0)
                    
                axes2[0].errorbar(xvals, diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                       fmt='none', capsize=capsize, ecolor=eb_color[li])
                min_values[yi].append(min_arr)
                max_values[yi].append(max_arr)
                

        ### TOTAL
        sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(len(xvals))] )
        if yi==1:
            axes2[0].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        else:
            axes2[0].bar(xvals[-1]+WIDTH+1.33+li, sums[-1], width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors[li], alpha=0.5, edgecolor=hemi_colors[li], hatch=HATCH, lw=2)
            
        if yi==1:
            min_sums = np.array( [np.nansum(np.array(min_values[yi])[:,dr]) for dr in range(len(xvals))] )
            max_sums = np.array( [np.nansum(np.array(max_values[yi])[:,dr]) for dr in range(len(xvals))] )
            axes2[0].errorbar(xvals+1.33+li, sums, yerr=[np.abs(min_sums-sums), np.abs(max_sums-sums)],
                                   fmt='none', capsize=capsize, ecolor=eb_color[li])
       
        axes2[0].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
    
#### sum
for yi, years in enumerate(decades):
    sums = []
    for vi, VAL in enumerate([values, valmin, valmax]): # each hemi
        val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
        val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
        sums.append( val0 + val1 )
    
        if vi==0:
            heights = []
            for mi, sum1 in enumerate(sums[0][0]):
                xvals = mi+np.array(XSPACING)
                heights.append(sum1)
                
                    # original
                fcs = [hemi_colors2[0][i] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
                ecs = [hemi_colors[0] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors[1] for i in [0,1,2]]
                
                if yi==1:
                    axes2[1].bar(xvals, sum1, width=WIDTH,
                            facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)
                else:
                    axes2[1].bar(xvals[-1]+WIDTH, sum1[-1], width=WIDTH,
                            facecolor=ecs[-1], alpha=0.5, edgecolor=ecs[-1], hatch=HATCH,lw=2)

            # total sum
            if yi==1:
                axes2[1].bar(xvals+1.33+0.5, np.nansum(sums[0][0], axis=0), width=WIDTH,
                        facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
            else:
                axes2[1].bar(xvals[-1]+WIDTH+1.33+0.5, np.nansum(sums[0][0], axis=0)[-1], width=WIDTH,
                        facecolor='k', alpha=0.5, edgecolor='k', hatch=HATCH, lw=2)
            
            axes2[1].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
            
        # add spread bars!
        elif vi==1 and yi==1:
            val0_min = val0
            val1_min = val1
        elif vi==2 and yi==1:
            for mi in np.arange(0,12):
                xvals = mi+np.array(XSPACING)
                
                axes2[1].errorbar(xvals, heights[mi], 
                                yerr=[np.abs(sums[1][mi]-heights[mi]), np.abs(sums[2][mi]-heights[mi])],
                                fmt='none', capsize=capsize, ecolor='#262626')
                
            # total sum 
            ht1 =  np.nansum(sums[0][0], axis=0)
            axes2[1].errorbar(xvals+1.33+0.5,  ht1, 
                            yerr=[np.abs(np.nansum(sums[1], axis=0)-ht1), np.abs(np.nansum(sums[2], axis=0)-ht1)],
                            fmt='none', capsize=capsize+0.5, ecolor='k')
 
# -----------------------------------------------------------------------------------
#### repeat for new miz definition markers
# -----------------------------------------------------------------------------------

values = {}
valmin = {}; valmax = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    path1 = root_paths[li]
    
    yi=0
        
    # trial = {'original': trials['original']}
        
            #name   #data path #seaice path #area-ind #p-thresh # sea ice name
            #'original':[path1, path1+'seaice/', 0, [984, 957], '_seaice'],
    trial = {'si3':[path1, path1+'seaice_miz/', 0, [984, 957], '_seaice_3'],
             'si_start':[path1, path1+'seaice_miz/', 0, [984, 957],'_seaice_start'],
             }
    
    era_values = {yi:[[] for x in np.arange(0, len(trial.keys()))]}
    min_values = {yi:[]}; max_values = {yi:[]}

    for k, name in enumerate(list(trial.keys())):
        attr = trial[name]
        
        # adjusted climatology
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', trial[name][1], si_name=trial[name][-1])

        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            else: MEANLINE = np.nansum(rel_area, axis=0)/1e6
        
            if len(lines[mi+1])<10: 
                era_values[yi][k].append([np.nan]*3)
                continue 
        
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
            
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi][k].append(diffs)
            
            # plot bars/ comparison dots
            axes2[0].scatter(xvals, diffs, s=s_size, marker='d',
                    facecolor=hemi_colors2[li], alpha=1, edgecolor=hemi_colors[li])
            

    ### TOTAL
    sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(len(xvals))] )
    if yi==0:
        axes2[0].scatter(xvals+1.33+li, sums, s=s_size, marker='d',
                facecolor=hemi_colors2[li], alpha=1, edgecolor=hemi_colors[li])
    
    axes2[0].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
    
    values[li] = era_values
    
#### sum
sums = []
for vi, VAL in enumerate([values]): # each hemi
    val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
    val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
    sums.append( val0 + val1 )

    if vi==0:
        heights = []
        for sumx in sums[0]:
            for mi, sum1 in enumerate(sumx):
                xvals = mi+np.array(XSPACING)
                heights.append(sum1)
                
                    # original
                fcs = [hemi_colors2[0][i] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
                ecs = [hemi_colors[0] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors[1] for i in [0,1,2]]
                
                axes2[1].scatter(xvals, sum1, s=s_size, marker='d',
                        facecolor=fcs, alpha=1, edgecolor=ecs)

            # total sum
            axes2[1].scatter(xvals+1.33+0.5, np.nansum(sumx, axis=0), s=s_size, marker='d',
                    facecolor=shades, alpha=1, edgecolor='k')
            axes2[1].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
            
#### axes
axes2[0].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
axes2[1].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
for ax in axes2.flatten():
    ax.xaxis.set_tick_params(labelbottom=True)
    ax.yaxis.set_tick_params(labelleft=True)
    

#%%% miz anomalies from the 3-week miz area
# Have you considered calculating the MIZ extent within the Storm Area 
# for each day of this period and deriving anomalies relative to the 
# mean MIZ over the full period? 

#### set up plot
fig2, axes2 = plt.subplots(2,2, figsize=(18,14), sharey='row')
fig2.suptitle('\nMIZ area anomalies relative to the 3-week mean MIZ',fontsize=fs+6, fontweight='bold')
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.axvline(0, ls=':', color='gray', lw=1)
    ax2.set_xlim(-7,14)
    ax2.set_xticks(xxx)
    ax2.set_xticklabels(xlabels, minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
for ax2 in axes2[1,:]: ax2.set_xlabel('Days Since Storm Start',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel(r'Daily MIZ anomalies from 3-week mean ($\times10^5$km$^2$)',fontsize=fs)
    
#### data and plot mean lines
for li, loc in enumerate(['Arctic', 'Antarctic']):
    for yi, years in enumerate(decades):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        
        if li==0:
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/', miz_ind=0) #[0=daily, 1=total]
        elif li==1:
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice_miz0/', si_name='_seaice0')

        lines_anoms = {}
        for mi, mm in enumerate(months):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                
                lines_anoms[mm] = [((si-clim)-np.nanmean((si-clim)))/1e5 for si, clim in zip(si_changes[mm], clim_changes[mm])]
                # lines_anoms[mm] = [((si)-np.nanmean((si)))/1e5 for si, clim in zip(si_changes[mm], clim_changes[mm])]
            
            if len(lines[mm])>10:
                ml = np.nanmean(lines_anoms[mm], axis=0)
                
                axes2[li][yi].plot(xxx, ml, color = month_colors[mi], lw=2, label = calendar.month_name[mm])

        axes2[li][yi].set_title(loc+': '+yr_title,fontsize=fs+2 )

        #### monthly lines
        fig, axes = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
        fig.suptitle(loc+': '+str(years[0])+'-'+str(years[-1])+'\n'+r'Daily MIZ anomalies from 3-week mean (km$^2$)')
        
        axes = axes.flatten()
        for ai, ax1 in enumerate(axes): 
            ax1.set_title(calendar.month_name[ai+1])
            ax1.set_xlim([-7,14])
            # ax1.set_ylim([-1.05,1.05])
            ax1.axhline(0, lw=0.55, color='k')
            ax1.axvline(0, lw=0.55, color='k')
            
        for mi, mm in enumerate(months):
            if len(lines_anoms[mm]) > 10:
                axes[mi].plot(xxx, np.array(lines_anoms[mm]).T, color = 'gray', lw=0.55)
                ml = np.nanmean(lines_anoms[mm], axis=0)
                axes[mi].plot(xxx, ml, color = month_colors[mi], lw=2)
        
        
        
axes2[0][-1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                  edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                  bbox_to_anchor=(1.33,0.85));
    

#%% Climatology definition
'''
I think that removing a 10-year climatological mean might bias some 
of the results for the 2010-2019 period, where the trend is negative 
in both hemispheres, and suggest some early 10’s storms are having 
more of an effect than they are and the later 10’s storms are having 
less of an effect than they are. Perhaps a 10-year climatological 
mean centered on the date of the storm could remove that? 
'''
years = np.arange(2010,2020)
yr_title = str(years[0])+'-'+str(years[-1])

#### set up plot
fig2, axes2 = plt.subplots(2,3, figsize=(27,14), sharey=True)
fig2.suptitle('\n'+yr_title, fontsize=fs+6, fontweight='bold')
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.axvline(0, ls=':', color='gray', lw=1)
    ax2.set_xlim(-7,14)
    ax2.set_xticks(xxx)
    ax2.set_xticklabels(xlabels, minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
for ax2 in axes2[-1,:]: ax2.set_xlabel('Days Since Storm Start',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Normalized Relative Change in Ice Area',fontsize=fs+1)
    
#### data and plot mean lines
    
for li, loc in enumerate(['Arctic', 'Antarctic']):
    # if li==1: break

    si_titles = ['10-year Centered Climatology', '2010-2019 Climatology', 'Difference']
    for i, si_title in enumerate(si_titles):
        axes2[li][i].set_title(loc+'\n'+si_title,fontsize=fs+2 )
    
    path1 = root_paths[li]
    
    mean_lines_original = fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')[0]
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice_clim/', si_name='_seaice_clim')

    lines_anoms = {}
    for mi, mm in enumerate(months):
        axes2[li][0].plot(xxx, mean_lines[mi], color = month_colors[mi], lw=2)
        axes2[li][1].plot(xxx, mean_lines_original[mi], color = month_colors[mi], lw=2)
        axes2[li][2].plot(xxx, mean_lines[mi] - mean_lines_original[mi], 
                      color = month_colors[mi], lw=2, label = calendar.month_name[mm])


    #### monthly lines
    if li==0:
        fig, axes = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
        fig.suptitle(si_titles[li]+': '+str(years[0])+'-'+str(years[-1]))
        
        axes = axes.flatten()
        for ai, ax1 in enumerate(axes): 
            ax1.set_title(calendar.month_name[ai+1])
            ax1.set_xlim([-7,14])
            # ax1.set_ylim([-1.05,1.05])
            ax1.axhline(0, lw=0.55, color='k')
            ax1.axvline(0, lw=0.55, color='k')
            
        for mi, mm in enumerate(months):
            axes[mi].plot(xxx, np.array(lines[mm]).T, color = 'gray', lw=0.55)
            axes[mi].plot(xxx, mean_lines[mi], color = month_colors[mi], lw=2)
            axes[mi].plot(xxx, mean_lines_original[mi], ls='--', color = month_colors[mi], lw=2)
        
        
        
axes2[0][-1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                  edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                  bbox_to_anchor=(1.33,0.85));

#%%% bars

#### settings

# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

s_size = 27

XSPACING = [0,0.2,0.4]
INDS = [7+3, 7+7, 7+14]
WIDTH = XSPACING[1]-XSPACING[0]
HATCH = '/////'
alphbar = ['a','b','c','d']

# mean cyclone impact vs summed across all storms
plot_mean = False
# normalized impact vs total change in area
norm_lines = False
# months 1->12 vs shifted SH seasonal cycle
shift = False 
# use same trials for min/max for all months vs use each month's min/max
set_minmax = True
set_trials = [3,1]

# errorbar 
capsize=2
eb_color= hemi_colors #['k','k']


#### set up
fig2, axes2 = plt.subplots(2,1, figsize=(10,12), sharex=True, sharey=True)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
    if plot_mean: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    else: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)
# ax2.set_xlabel('Month',fontsize=fs)

axes2[0].set_title('MIZ Ice Area Changes in Each Hemisphere', fontsize=fs+2)
axes2[1].set_title('Global Sum', fontsize=fs+2)

plt.subplots_adjust(hspace=0.275)

import matplotlib.patches as mpatches
circ1 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,
                       label=str(decades[1][0])+'-'+str(decades[1][-1]))
ax2.legend(handles = [circ1], loc='lower left',
           fontsize=fs-1, handletextpad=0.5, handlelength=1.5)
    
    
fig2.subplots_adjust(wspace=0.1)

for li, loc in enumerate(['Arctic','Antarctic']): axes2[0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
axes2[0].scatter([],[], s=s_size, marker='d', color='gray', label='Adjusted Climatology')
axes2[0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)


#### data and plot

for xi in ['original', 'clim']:

    values = {}
    valmin = {}; valmax = {}
    for li, loc in enumerate(['Arctic','Antarctic']):
        
        path1 = root_paths[li]
        path2 = root_paths[li]+'sensitivity/'
        
                #name   #data path #seaice path #area-ind #p-thresh # colors
        trials = {'original':[path1, path1+'seaice/', 0, [984, 957], 
                             ['lightcoral','indianred','maroon']],
                 'area_990':[path1, path1+'sensitivity_seaice/', 1, [984, 957],
                             ['royalblue','mediumblue','navy']],
                 'area_cont':[path1, path1+'sensitivity_seaice/', 2, [984, 957],
                              ['mediumaquamarine','mediumseagreen','green']],
                 'p_1000':[path2, path2+'seaice/', 0, [994,967],
                           ['gold','goldenrod','darkgoldenrod']],
                 # 'p_cont':[path2, path2+'seaice/', 2, [994,967],
                 #           ['violet','mediumorchid','darkviolet']],
                 }
        
        
        # data loop, bar plot
        era_values = {}
        min_values = {}; max_values = {}
        for yi, years in enumerate(decades[::-1]):
            if yi==1: continue
            
            yr_title = str(years[0])+'-'+str(years[-1])
            
            if yi == 0: trial = {'original': trials['original']}
            elif yi == 1: trial = trials
    
            era_values[yi] = [[] for x in np.arange(0, len(trial.keys()))]
            min_values[yi] = []
            max_values[yi] = []
            
            for k, name in enumerate(list(trial.keys())):
                attr = trial[name]
                
                if xi=='original':
                    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                        seaice_lines(years, attr[0]+'census/', attr[0]+'area/', attr[1], 
                                     AREA_IND=attr[2], p_thresh=attr[3][li])
                elif xi=='clim':
                    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice_clim/', si_name='_seaice_clim')

    
                for mi, mm in enumerate(months):
                    rel_area = []
                    for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                        ss = sia-sia_clim
                        rel_area.append(ss-ss[0])
                    if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
                    else: MEANLINE = np.nansum(rel_area, axis=0)/1e6
                
                    if len(lines[mi+1])<10: 
                        era_values[yi][k].append([np.nan]*3)
                        continue 
                
                    if not shift: xvals = mi+np.array(XSPACING)
                    else:
                        xvals = mi+np.array(XSPACING) + 6
                        if xvals[0]>=12: xvals -= 12
                    
                    diffs = [MEANLINE[ii] for ii in INDS]
                    era_values[yi][k].append(diffs)
                    
                    # plot bars/ comparison dots
                    if name == 'original':
                        if xi=='original':
                            if yi==0:
                                axes2[0].bar(xvals, diffs, width=WIDTH,
                                        facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
                            else:
                                axes2[0].bar(xvals[-1]+WIDTH, diffs[-1], width=WIDTH,
                                        facecolor=hemi_colors[li], alpha=0.5, 
                                        edgecolor=hemi_colors[li], hatch=HATCH, lw=2)
                        elif xi=='clim':
                            axes2[0].scatter(xvals, diffs, s=s_size, marker='d',
                                    facecolor=hemi_colors2[li], alpha=1, edgecolor=hemi_colors[li])
    
            #### find min/max values after each trial is calculated
            if yi==1:
                for mi, mm in enumerate(months):
                    if not shift: xvals = mi+np.array(XSPACING)
                    else: 
                        xvals = mi+np.array(XSPACING) + 6
                        if xvals[0]>=12: xvals -= 12
                    
                    # plot spread bars
                    diff_arr = np.array([era_values[yi][k][mi] for k in np.arange(0, len(trial))])
                    if set_minmax:
                        min_arr = diff_arr[set_trials[0]]
                        max_arr = diff_arr[set_trials[1]]
                    else:
                        min_arr = np.nanmin(diff_arr, axis=0)
                        max_arr = np.nanmax(diff_arr, axis=0)
                        
                    axes2[0].errorbar(xvals, diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                           fmt='none', capsize=capsize, ecolor=eb_color[li])
                    min_values[yi].append(min_arr)
                    max_values[yi].append(max_arr)
                    
    
            ### TOTAL
            sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(len(xvals))] )
            if yi==0:
                if xi=='original':
                    axes2[0].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                            facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
                elif xi=='clim':
                    axes2[0].scatter(xvals+1.33+li, sums, s=s_size, marker='d',
                            facecolor=hemi_colors2[li], alpha=1, edgecolor=hemi_colors[li])
            
            else:
                axes2[0].bar(xvals[-1]+WIDTH+1.33+li, sums[-1], width=XSPACING[1]-XSPACING[0],
                        facecolor=hemi_colors[li], alpha=0.5, edgecolor=hemi_colors[li], hatch=HATCH, lw=2)
                
            if yi==1:
                min_sums = np.array( [np.nansum(np.array(min_values[yi])[:,dr]) for dr in range(len(xvals))] )
                max_sums = np.array( [np.nansum(np.array(max_values[yi])[:,dr]) for dr in range(len(xvals))] )
                axes2[0].errorbar(xvals+1.33+li, sums, yerr=[np.abs(min_sums-sums), np.abs(max_sums-sums)],
                                       fmt='none', capsize=capsize, ecolor=eb_color[li])
           
            axes2[0].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
            
        values[li] = era_values
        valmin[li] = min_values
        valmax[li] = max_values
        
    #### sum
    for yi, years in enumerate(decades[::-1]):
        if yi==1: continue
        sums = []
        for vi, VAL in enumerate([values, valmin, valmax]): # each hemi
            val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
            val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
            sums.append( val0 + val1 )
        
            if vi==0:
                heights = []
                for mi, sum1 in enumerate(sums[0][0]):
                    xvals = mi+np.array(XSPACING)
                    heights.append(sum1)
                    
                        # original
                    fcs = [hemi_colors2[0][i] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
                    ecs = [hemi_colors[0] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors[1] for i in [0,1,2]]
                    
                    if yi==0:
                        if xi=='original':
                            axes2[1].bar(xvals, sum1, width=WIDTH,
                                    facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)
                        elif xi=='clim':
                            axes2[1].scatter(xvals, sum1, s=s_size, marker='d',
                                    facecolor=fcs, alpha=1, edgecolor=ecs)
                    else:
                        axes2[1].bar(xvals[-1]+WIDTH, sum1[-1], width=WIDTH,
                                facecolor=ecs[-1], alpha=0.5, edgecolor=ecs[-1], hatch=HATCH,lw=2)
    
                # total sum
                if yi==0:
                    
                    
                    if xi=='original':
                        axes2[1].bar(xvals+1.33+0.5, np.nansum(sums[0][0], axis=0), width=WIDTH,
                                facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
                    elif xi=='clim':
                        axes2[1].scatter(xvals+1.33+0.5, np.nansum(sums[0][0], axis=0), s=s_size, marker='d',
                                facecolor=shades, alpha=1, edgecolor='k')
                    
                else:
                    axes2[1].bar(xvals[-1]+WIDTH+1.33+0.5, np.nansum(sums[0][0], axis=0)[-1], width=WIDTH,
                            facecolor='k', alpha=0.5, edgecolor='k', hatch=HATCH, lw=2)
                
                axes2[1].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
                
            # add spread bars!
            elif vi==1 and yi==1:
                val0_min = val0
                val1_min = val1
            elif vi==2 and yi==1:
                for mi in np.arange(0,12):
                    xvals = mi+np.array(XSPACING)
                    
                    axes2[1].errorbar(xvals, heights[mi], 
                                    yerr=[np.abs(sums[1][mi]-heights[mi]), np.abs(sums[2][mi]-heights[mi])],
                                    fmt='none', capsize=capsize, ecolor='#262626')
                    
                # total sum 
                ht1 =  np.nansum(sums[0][0], axis=0)
                axes2[1].errorbar(xvals+1.33+0.5,  ht1, 
                                yerr=[np.abs(np.nansum(sums[1], axis=0)-ht1), np.abs(np.nansum(sums[2], axis=0)-ht1)],
                                fmt='none', capsize=capsize+0.5, ecolor='k')
                
                print('* total spread')
                print(ht1*100/np.nansum(sums[1], axis=0))
    
    
    # print('* fractional contributions')
    # nh_cont = np.array( [np.nansum(np.array(values[0][yi][0])[:,dr]) for dr in range(len(xvals))] )
    # sh_cont = np.array( [np.nansum(np.array(values[1][yi][0])[:,dr]) for dr in range(len(xvals))] )
    # total = nh_cont + sh_cont
    
    # print('NH: '+str([round(x,2) for x in nh_cont]))
    # print('SH: '+str([round(x,2) for x in sh_cont]))
    # print('% NH contribution: '+str([round(x/total[i],2) for i,x in enumerate(nh_cont)]))
    # print('% SH contribution: '+str([round(x/total[i],2) for i,x in enumerate(sh_cont)]))


axes2[0].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
axes2[1].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
for ax in axes2.flatten():
    ax.xaxis.set_tick_params(labelbottom=True)
    ax.yaxis.set_tick_params(labelleft=True)
    
    
    
    
#%%%% try again with spread

#### settings
# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

XSPACING = [0,0.2,0.4]
INDS = [7+3, 7+7, 7+14]
WIDTH = XSPACING[1]-XSPACING[0]
HATCH = '/////'
alphbar = ['a','b','c','d']

# mean cyclone impact vs summed across all storms
plot_mean = False
# normalized impact vs total change in area
norm_lines = False
# months 1->12 vs shifted SH seasonal cycle
shift = False 
# use same trials for min/max for all months vs use each month's min/max
set_minmax = True
set_trials = [3,1]

# errorbar 
capsize=2
eb_color= hemi_colors #['k','k']


#### set up
fig2, axes2 = plt.subplots(2,1, figsize=(10,12), sharex=True, sharey=True)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
    if plot_mean: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    else: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)
# ax2.set_xlabel('Month',fontsize=fs)

axes2[0].set_title('MIZ Ice Area Changes in Each Hemisphere', fontsize=fs+2)
axes2[1].set_title('Global Sum', fontsize=fs+2)

plt.subplots_adjust(hspace=0.275)

import matplotlib.patches as mpatches
circ1 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,
                       label=str(decades[1][0])+'-'+str(decades[1][-1]))
circ2 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,hatch=HATCH,
                       label=str(decades[0][0])+'-'+str(decades[0][-1]))

ax2.legend(handles = [circ2,circ1], loc='lower left',
           fontsize=fs-1, handletextpad=0.5, handlelength=1.5)
    
    
fig2.subplots_adjust(wspace=0.1)

for li, loc in enumerate(['Arctic','Antarctic']): axes2[0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
axes2[0].scatter([],[], s=s_size, marker='d', color='gray', label='Adjusted Climatology')
axes2[0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)

#### data and plot

values = {}
valmin = {}; valmax = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    
    path1 = root_paths[li]
    path2 = root_paths[li]+'sensitivity/'
    
            #name   #data path #seaice path #area-ind #p-thresh # colors
    trials = {'original':[path1, path1+'seaice/', 0, [984, 957], 
                         ['lightcoral','indianred','maroon']],
             'area_990':[path1, path1+'sensitivity_seaice/', 1, [984, 957],
                         ['royalblue','mediumblue','navy']],
             'area_cont':[path1, path1+'sensitivity_seaice/', 2, [984, 957],
                          ['mediumaquamarine','mediumseagreen','green']],
             'p_1000':[path2, path2+'seaice/', 0, [994,967],
                       ['gold','goldenrod','darkgoldenrod']],
             # 'p_cont':[path2, path2+'seaice/', 2, [994,967],
             #           ['violet','mediumorchid','darkviolet']],
             }
    
    # data loop, bar plot
    era_values = {}
    min_values = {}; max_values = {}
    for yi, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        if yi == 0: trial = {'original': trials['original']}
        elif yi == 1: trial = trials

        era_values[yi] = [[] for x in np.arange(0, len(trial.keys()))]
        min_values[yi] = []
        max_values[yi] = []
        
        for k, name in enumerate(list(trial.keys())):
            attr = trial[name]
            
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                seaice_lines(years, attr[0]+'census/', attr[0]+'area/', attr[1], 
                             AREA_IND=attr[2], p_thresh=attr[3][li])

            for mi, mm in enumerate(months):
                rel_area = []
                for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                    ss = sia-sia_clim
                    rel_area.append(ss-ss[0])
                if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
                else: MEANLINE = np.nansum(rel_area, axis=0)/1e6
            
                if len(lines[mi+1])<10: 
                    era_values[yi][k].append([np.nan]*3)
                    continue 
            
                if not shift: xvals = mi+np.array(XSPACING)
                else:
                    xvals = mi+np.array(XSPACING) + 6
                    if xvals[0]>=12: xvals -= 12
                
                diffs = [MEANLINE[ii] for ii in INDS]
                era_values[yi][k].append(diffs)
                
                # plot bars/ comparison dots
                if name == 'original':
                    if yi==1:
                        axes2[0].bar(xvals, diffs, width=WIDTH,
                                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
                    else:
                        axes2[0].bar(xvals[-1]+WIDTH, diffs[-1], width=WIDTH,
                                facecolor=hemi_colors[li], alpha=0.5, 
                                edgecolor=hemi_colors[li], hatch=HATCH, lw=2)


        #### find min/max values after each trial is calculated
        if yi==1:
            for mi, mm in enumerate(months):
                if not shift: xvals = mi+np.array(XSPACING)
                else: 
                    xvals = mi+np.array(XSPACING) + 6
                    if xvals[0]>=12: xvals -= 12
                
                # plot spread bars
                diff_arr = np.array([era_values[yi][k][mi] for k in np.arange(0, len(trial))])
                if set_minmax:
                    min_arr = diff_arr[set_trials[0]]
                    max_arr = diff_arr[set_trials[1]]
                else:
                    min_arr = np.nanmin(diff_arr, axis=0)
                    max_arr = np.nanmax(diff_arr, axis=0)
                    
                axes2[0].errorbar(xvals, diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                       fmt='none', capsize=capsize, ecolor=eb_color[li])
                min_values[yi].append(min_arr)
                max_values[yi].append(max_arr)
                

        ### TOTAL
        sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(len(xvals))] )
        if yi==1:
            axes2[0].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        else:
            axes2[0].bar(xvals[-1]+WIDTH+1.33+li, sums[-1], width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors[li], alpha=0.5, edgecolor=hemi_colors[li], hatch=HATCH, lw=2)
            
        if yi==1:
            min_sums = np.array( [np.nansum(np.array(min_values[yi])[:,dr]) for dr in range(len(xvals))] )
            max_sums = np.array( [np.nansum(np.array(max_values[yi])[:,dr]) for dr in range(len(xvals))] )
            axes2[0].errorbar(xvals+1.33+li, sums, yerr=[np.abs(min_sums-sums), np.abs(max_sums-sums)],
                                   fmt='none', capsize=capsize, ecolor=eb_color[li])
       
        axes2[0].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
    
#### sum
for yi, years in enumerate(decades):
    sums = []
    for vi, VAL in enumerate([values, valmin, valmax]): # each hemi
        val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
        val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
        sums.append( val0 + val1 )
    
        if vi==0:
            heights = []
            for mi, sum1 in enumerate(sums[0][0]):
                xvals = mi+np.array(XSPACING)
                heights.append(sum1)
                
                    # original
                fcs = [hemi_colors2[0][i] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
                ecs = [hemi_colors[0] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors[1] for i in [0,1,2]]
                
                if yi==1:
                    axes2[1].bar(xvals, sum1, width=WIDTH,
                            facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)
                else:
                    axes2[1].bar(xvals[-1]+WIDTH, sum1[-1], width=WIDTH,
                            facecolor=ecs[-1], alpha=0.5, edgecolor=ecs[-1], hatch=HATCH,lw=2)

            # total sum
            if yi==1:
                axes2[1].bar(xvals+1.33+0.5, np.nansum(sums[0][0], axis=0), width=WIDTH,
                        facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
            else:
                axes2[1].bar(xvals[-1]+WIDTH+1.33+0.5, np.nansum(sums[0][0], axis=0)[-1], width=WIDTH,
                        facecolor='k', alpha=0.5, edgecolor='k', hatch=HATCH, lw=2)
            
            axes2[1].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
            
        # add spread bars!
        elif vi==1 and yi==1:
            val0_min = val0
            val1_min = val1
        elif vi==2 and yi==1:
            for mi in np.arange(0,12):
                xvals = mi+np.array(XSPACING)
                
                axes2[1].errorbar(xvals, heights[mi], 
                                yerr=[np.abs(sums[1][mi]-heights[mi]), np.abs(sums[2][mi]-heights[mi])],
                                fmt='none', capsize=capsize, ecolor='#262626')
                
            # total sum 
            ht1 =  np.nansum(sums[0][0], axis=0)
            axes2[1].errorbar(xvals+1.33+0.5,  ht1, 
                            yerr=[np.abs(np.nansum(sums[1], axis=0)-ht1), np.abs(np.nansum(sums[2], axis=0)-ht1)],
                            fmt='none', capsize=capsize+0.5, ecolor='k')
 
# -----------------------------------------------------------------------------------
#### repeat for clim markers
# -----------------------------------------------------------------------------------

values = {}
valmin = {}; valmax = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    path1 = root_paths[li]
    
    yi=0
    ears = np.arange(2010,2020)
        
    trial = {'original': trials['original']}
    
    era_values = {yi:[[] for x in np.arange(0, len(trial.keys()))]}
    min_values = {yi:[]}; max_values = {yi:[]}

    for k, name in enumerate(list(trial.keys())):
        attr = trial[name]
        
        # adjusted climatology
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice_clim/', si_name='_seaice_clim')

        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            else: MEANLINE = np.nansum(rel_area, axis=0)/1e6
        
            if len(lines[mi+1])<10: 
                era_values[yi][k].append([np.nan]*3)
                continue 
        
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
            
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi][k].append(diffs)
            
            # plot bars/ comparison dots
            axes2[0].scatter(xvals, diffs, s=s_size, marker='d',
                    facecolor=hemi_colors2[li], alpha=1, edgecolor=hemi_colors[li])
            

    ### TOTAL
    sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(len(xvals))] )
    if yi==0:
        axes2[0].scatter(xvals+1.33+li, sums, s=s_size, marker='d',
                facecolor=hemi_colors2[li], alpha=1, edgecolor=hemi_colors[li])
    
    axes2[0].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
    
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
    
#### sum
sums = []
for vi, VAL in enumerate([values, valmin, valmax]): # each hemi
    val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
    val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
    sums.append( val0 + val1 )

    if vi==0:
        heights = []
        for mi, sum1 in enumerate(sums[0][0]):
            xvals = mi+np.array(XSPACING)
            heights.append(sum1)
            
                # original
            fcs = [hemi_colors2[0][i] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
            ecs = [hemi_colors[0] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors[1] for i in [0,1,2]]
            
            axes2[1].scatter(xvals, sum1, s=s_size, marker='d',
                    facecolor=fcs, alpha=1, edgecolor=ecs)

        # total sum
        axes2[1].scatter(xvals+1.33+0.5, np.nansum(sums[0][0], axis=0), s=s_size, marker='d',
                facecolor=shades, alpha=1, edgecolor='k')
        axes2[1].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
            
#### axes
axes2[0].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
axes2[1].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
for ax in axes2.flatten():
    ax.xaxis.set_tick_params(labelbottom=True)
    ax.yaxis.set_tick_params(labelleft=True)

#%% cyclone variability [2014]
# this study needs to address the roles of cyclones during the freezing 
# and melting seasons before 2014 and after 2014, in order to understand 
# if variability in cyclones plays a role in the puzzling variability in 
# Antarctic sea ice area before and after 2014

# mean lines every year? 5-year periods?

#%%% split fans

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(16,14), sharey=True)
fig2.suptitle('\nComparing Storm Impacts Before & After 2014', 
              fontweight='bold', fontsize=fs+6)

alph = ['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k']
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.axvline(0, ls=':', color='gray', lw=1)
    ax2.set_xlim(-7,14)
    ax2.set_xticks(xxx)
    ax2.set_xticklabels(xlabels, minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Days Since Storm Start',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Normalized Relative Change in Ice Area',fontsize=fs+1)
    
#### data and plot
print('- range of fan plots (7-, 14-day)')
for li, loc in enumerate(['arctic','antarctic']):
    for yi, years in enumerate([np.arange(2010,2015), np.arange(2015,2020)]):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/', sample_size=3)
        
        ### monthly mean
        ranges = {7:[], 14:[]}
        axes2[li][yi].set_title(hemi_names[li]+': '+yr_title, fontsize=fs+2)
        for idx, ml in enumerate(mean_lines):
            if li==1 and idx in [0,1,2]: continue
            axes2[li][yi].plot([],[], color = month_colors[idx], lw=2, 
                               ls = '-',label = calendar.month_name[idx+1])
            
            axes2[li][yi].plot(xxx, ml, color = month_colors[idx], lw=2, 
                               ls = '-' if len(lines[idx+1])>=10 else '--')
            
            
            if len(lines[idx+1])<10:
                yadd=0 if yi==0 else 0.05
                axes2[li][yi].text(xxx[-1]+0.1, ml[-1]+yadd, 'n='+str(len(lines[idx+1])), fontsize=fs-1) 
            
            # yadd = ytext_scale(idx, li, yi)
            # mlabel = calendar.month_abbr[idx+1] #+' ('+str(len(lines[idx+1]))+')'
            # axes2[li][yi].text(xxx[-1]+0.25, ml[-1]+yadd, mlabel, fontsize=fs-1) 
            ranges[7].append(ml[7+7])
            ranges[14].append(ml[7+14])
            
        range1 = [round(np.nanmin(ranges[7]),4), round(np.nanmax(ranges[7]),4)]    
        range2 = [round(np.nanmin(ranges[14]),4), round(np.nanmax(ranges[14]),4)]    
        print(loc, yr_title, range1, range2)
        
axes2[0][1].plot([],[], lw=0, label=' ')
axes2[0][1].plot([],[], color='gray', lw=2, ls='-', label=r'Sample Size ≥ 10 Storms')
axes2[0][1].plot([],[], color='gray', lw=2, ls='--', label='Sample Size < 10 Storms')

axes2[0][1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                  edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                  bbox_to_anchor=(1.66,0.33));



#%%% net bar plot

#### settings

# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

# XSPACING = [0,0.2,0.4]
XSPACING = np.array([0,0.125,0.25, 0.5,0.625,0.75])
INDS = [7+3, 7+7, 7+14]
WIDTH = XSPACING[1]-XSPACING[0]
HATCH = '/////'
alphbar = ['a','b','c','d']

# mean cyclone impact vs summed across all storms
plot_mean = False
# normalized impact vs total change in area
norm_lines = False
# months 1->12 vs shifted SH seasonal cycle
shift = False 
# use same trials for min/max for all months vs use each month's min/max
set_minmax = True
set_trials = [3,1]

# errorbar 
capsize=2
eb_color= hemi_colors #['k','k']

# new comparison years
newyrs = [np.arange(2010,2015), np.arange(2015,2020)]


#### set up
fig2, axes2 = plt.subplots(2,1, figsize=(16,10), sharex=True, sharey=False)
axes2 = axes2.flatten()
for i, ax2 in enumerate(axes2[0:2]):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1])+(WIDTH*2))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.01, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
    if plot_mean: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    else: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)
    
    for xi in np.arange(0,13): ax2.axvline(xi-WIDTH, lw=0.55, color='darkgray')
    
# ax2.set_xlabel('Month',fontsize=fs)

axes2[0].set_title('MIZ Ice Area Changes in Each Hemisphere', fontsize=fs+2)
axes2[1].set_title('Global Sum', fontsize=fs+2)

plt.subplots_adjust(hspace=0.25)

import matplotlib.patches as mpatches

circ1 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,
                       label=str(newyrs[0][0])+'-'+str(newyrs[0][-1]))
circ2 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,hatch=HATCH,
                       label=str(newyrs[1][0])+'-'+str(newyrs[1][-1]))

ax2.legend(handles = [circ1,circ2], loc='lower left',
           fontsize=fs-1, handletextpad=0.5, handlelength=1.5)
    
fig2.subplots_adjust(wspace=0.1)



#### data and plot

values = {}
valmin = {}; valmax = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    
    path1 = root_paths[li]
    path2 = root_paths[li]+'sensitivity/'
    
            #name   #data path #seaice path #area-ind #p-thresh # colors
    trials = {'original':[path1, path1+'seaice/', 0, [984, 957], 
                         ['lightcoral','indianred','maroon']],
             'area_990':[path1, path1+'sensitivity_seaice/', 1, [984, 957],
                         ['royalblue','mediumblue','navy']],
             'area_cont':[path1, path1+'sensitivity_seaice/', 2, [984, 957],
                          ['mediumaquamarine','mediumseagreen','green']],
             'p_1000':[path2, path2+'seaice/', 0, [994,967],
                       ['gold','goldenrod','darkgoldenrod']],
             # 'p_cont':[path2, path2+'seaice/', 2, [994,967],
             #           ['violet','mediumorchid','darkviolet']],
             }
    
    axes2[0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    
    # data loop, bar plot
    era_values = {}
    min_values = {}; max_values = {}
    for yi, years in enumerate(newyrs):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        if yi == 0: trial = trials
        elif yi == 1: trial = trials

        era_values[yi] = [[] for x in np.arange(0, len(trial.keys()))]
        min_values[yi] = []
        max_values[yi] = []
        
        for k, name in enumerate(list(trial.keys())):
            attr = trial[name]
            
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                seaice_lines(years, attr[0]+'census/', attr[0]+'area/', attr[1], 
                             AREA_IND=attr[2], p_thresh=attr[3][li], sample_size=0)

            for mi, mm in enumerate(months):
                rel_area = []
                for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                    ss = sia-sia_clim
                    rel_area.append(ss-ss[0])
                if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
                else: MEANLINE = np.nansum(rel_area, axis=0)/1e6
            
                if len(lines[mi+1])<10: 
                    era_values[yi][k].append([np.nan]*3)
                    continue 
            
                if not shift: xvals = mi+np.array(XSPACING)
                else:
                    xvals = mi+np.array(XSPACING) + 6
                    if xvals[0]>=12: xvals -= 12
                
                diffs = [MEANLINE[ii] for ii in INDS]
                era_values[yi][k].append(diffs)
                
                # plot bars/ comparison dots
                if name == 'original':
                    if yi==0:
                        axes2[0].bar(xvals[0:3], diffs, width=WIDTH,
                                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
                    else:
                        axes2[0].bar(xvals[3:], diffs, width=WIDTH,
                                facecolor=hemi_colors2[li], alpha=0.5, 
                                edgecolor=hemi_colors2[li], hatch=HATCH, lw=2)


        #### find min/max values after each trial is calculated
        for mi, mm in enumerate(months):
            if not shift: xvals = mi+np.array(XSPACING)
            else: 
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
            
            # plot spread bars
            diff_arr = np.array([era_values[yi][k][mi] for k in np.arange(0, len(trial))])
            diff_arr=np.squeeze(diff_arr)
            if set_minmax:
                min_arr = diff_arr[set_trials[0]]
                max_arr = diff_arr[set_trials[1]]
            else:
                min_arr = np.nanmin(diff_arr, axis=0)
                max_arr = np.nanmax(diff_arr, axis=0)
                
            if yi==0:
                axes2[0].errorbar(xvals[0:3], diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                   fmt='none', capsize=capsize, ecolor=eb_color[li])
            else:
                axes2[0].errorbar(xvals[3:], diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                   fmt='none', capsize=capsize, ecolor=eb_color[li])
                
            min_values[yi].append(min_arr)
            max_values[yi].append(max_arr)
                

        ### TOTAL
        sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(3)] )
        if yi==0:
            axes2[0].bar(xvals[0:3]+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
            xvo = xvals[0:3]
        else:
            axes2[0].bar(xvals[3:]+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], hatch=HATCH, lw=2)
            xvo = xvals[3:]
            
        min_sums = np.array( [np.nansum(np.array(min_values[yi])[:,dr]) for dr in range(3)] )
        max_sums = np.array( [np.nansum(np.array(max_values[yi])[:,dr]) for dr in range(3)] )
        axes2[0].errorbar(xvo+1.33+li, sums, yerr=[np.abs(min_sums-sums), np.abs(max_sums-sums)],
                               fmt='none', capsize=capsize, ecolor=eb_color[li])
       
        axes2[0].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
    
#### sum
for yi, years in enumerate(decades):
    sums = []
    for vi, VAL in enumerate([values, valmin, valmax]): # each hemi
        val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
        val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
        sums.append( val0 + val1 )
    
        if vi==0:
            heights = []
            for mi, sum1 in enumerate(sums[0][0]):
                xvals = mi+np.array(XSPACING)
                heights.append(sum1)
                
                    # original
                fcs = [hemi_colors2[0][i] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
                ecs = [hemi_colors[0] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors[1] for i in [0,1,2]]
                
                if yi==0:
                    axes2[1].bar(xvals[0:3], sum1, width=WIDTH,
                            facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)
                else:
                    axes2[1].bar(xvals[3:], sum1, width=WIDTH,
                            facecolor=fcs, alpha=0.5, edgecolor=ecs, hatch=HATCH,lw=2)

            # total sum
            if yi==0:
                axes2[1].bar(xvals[0:3]+1.33+0.5, np.nansum(sums[0][0], axis=0), width=WIDTH,
                        facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
                print('* total sum values')
                for i, ind in enumerate(INDS):
                    print('-','day',str(ind-7)+':', round(np.nansum(sums[0][0], axis=0)[i],3))
            else:
                axes2[1].bar(xvals[3:]+1.33+0.5, np.nansum(sums[0][0], axis=0), width=WIDTH,
                        facecolor=shades, alpha=0.5, edgecolor='k', hatch=HATCH, lw=2)
                print('-','80s:', round(np.nansum(sums[0][0], axis=0)[-1],3))
            
            axes2[1].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
            
        # add spread bars!
        elif vi==1 and yi==1:
            val0_min = val0
            val1_min = val1
        elif vi==2:
            for mi in np.arange(0,12):
                if yi==0: xvals = mi+np.array(XSPACING[0:3])
                elif yi==1: xvals = mi+np.array(XSPACING[3:])
                
                axes2[1].errorbar(xvals, heights[mi], 
                                yerr=[np.abs(sums[1][mi]-heights[mi]), np.abs(sums[2][mi]-heights[mi])],
                                fmt='none', capsize=capsize, ecolor='#262626')
                
            # total sum 
            ht1 =  np.nansum(sums[0][0], axis=0)
            axes2[1].errorbar(xvals+1.33+0.5,  ht1, 
                            yerr=[np.abs(np.nansum(sums[1], axis=0)-ht1), np.abs(np.nansum(sums[2], axis=0)-ht1)],
                            fmt='none', capsize=capsize+0.5, ecolor='k')
            
            print('* total spread')
            print(ht1*100/np.nansum(sums[1], axis=0))

# x-axis
axes2[0].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
axes2[1].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
for ax in axes2.flatten():
    ax.xaxis.set_tick_params(labelbottom=True)
    ax.yaxis.set_tick_params(labelleft=True)


print('* fractional contributions')
nh_cont = np.array( [np.nansum(np.array(values[0][yi][0])[:,dr]) for dr in range(len(xvals))] )
sh_cont = np.array( [np.nansum(np.array(values[1][yi][0])[:,dr]) for dr in range(len(xvals))] )
total = nh_cont + sh_cont

print('NH: '+str([round(x,2) for x in nh_cont]))
print('SH: '+str([round(x,2) for x in sh_cont]))
print('% NH contribution: '+str([round(x/total[i],2) for i,x in enumerate(nh_cont)]))
print('% SH contribution: '+str([round(x/total[i],2) for i,x in enumerate(sh_cont)]))

#%%%% storm counts
fig3, ax3 = plt.subplots(1,1, figsize=(16,3))

hemi_colors3 = [['#018571','#80cdc1'],['#a6611a','#dfc27d']]

for li, loc in enumerate(['Arctic','Antarctic']):
    for yi, years in enumerate([np.arange(2010,2015), np.arange(2015,2019)]):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/', sample_size=0)


        ax3.plot(months-1+(WIDTH*3), [len(lines[mm]) for mm in months], lw=1.75,
                   color=hemi_colors3[li][yi], marker = ['o','^'][yi], markersize=8,
                   label = loc+': '+yr_title)

ax3.set_title('Monthly Storm Counts', fontsize=fs+2)

ax3.set_xlim(-0.5,14.25)
ax3.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1])+(WIDTH*2)))
ax3.yaxis.set_tick_params(labelleft=True)
ax3.tick_params(axis='both', labelsize=fs)
ax3.text(0.01, 1.025, '('+alphbar[i]+')',transform=ax3.transAxes, 
          fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                              'edgecolor':'white', 'lw':0.75},zorder=50)

ax3.set_ylabel('Number of Storms', fontsize=fs)
ax3.set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)], fontsize=fs)
ax3.xaxis.set_tick_params(labelbottom=True)
ax3.yaxis.set_tick_params(labelleft=True)

ax3.legend(loc='upper right', fontsize=fs-1, handletextpad=0.33, handlelength=0.75)

'''divide bars by storm counts...'''

#%%% combine: bar+counts

#### settings

# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

# XSPACING = [0,0.2,0.4]
XSPACING = np.array([0,0.125,0.25, 0.5,0.625,0.75])
INDS = [7+3, 7+7, 7+14]
WIDTH = XSPACING[1]-XSPACING[0]
HATCH = '/////'
alphbar = ['a','b','c','d']

# mean cyclone impact vs summed across all storms
plot_mean = False
# normalized impact vs total change in area
norm_lines = False
# months 1->12 vs shifted SH seasonal cycle
shift = False 
# use same trials for min/max for all months vs use each month's min/max
set_minmax = True
set_trials = [3,1]

# errorbar 
capsize=2
eb_color= hemi_colors #['k','k']

# new comparison years
newyrs = [np.arange(2010,2015), np.arange(2015,2020)]


#### set up
fig2, axes_all = plt.subplots(3,1, figsize=(16,14), height_ratios=[1,1,0.8],
                              sharex=False, sharey=False)

axes2 = np.array([axes_all[0], axes_all[1]])
for i, ax2 in enumerate(axes2[0:2]):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1])+(WIDTH*2))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.01, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
    if plot_mean: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    else: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)
    
    for xi in np.arange(0,13): ax2.axvline(xi-WIDTH, lw=0.55, color='darkgray')
    
# ax2.set_xlabel('Month',fontsize=fs)

axes2[0].set_title('MIZ Ice Area Changes in Each Hemisphere', fontsize=fs+2)
axes2[1].set_title('Global Sum', fontsize=fs+2)

plt.subplots_adjust(hspace=0.25)

import matplotlib.patches as mpatches

circ1 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,
                       label=str(newyrs[0][0])+'-'+str(newyrs[0][-1]))
circ2 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,hatch=HATCH,
                       label=str(newyrs[1][0])+'-'+str(newyrs[1][-1]))

ax2.legend(handles = [circ1,circ2], loc='lower left',
           fontsize=fs-1, handletextpad=0.5, handlelength=1.5)
    
fig2.subplots_adjust(wspace=0.1)



#### data and plot

values = {}
valmin = {}; valmax = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    
    path1 = root_paths[li]
    path2 = root_paths[li]+'sensitivity/'
    
            #name   #data path #seaice path #area-ind #p-thresh # colors
    trials = {'original':[path1, path1+'seaice/', 0, [984, 957], 
                         ['lightcoral','indianred','maroon']],
             'area_990':[path1, path1+'sensitivity_seaice/', 1, [984, 957],
                         ['royalblue','mediumblue','navy']],
             'area_cont':[path1, path1+'sensitivity_seaice/', 2, [984, 957],
                          ['mediumaquamarine','mediumseagreen','green']],
             'p_1000':[path2, path2+'seaice/', 0, [994,967],
                       ['gold','goldenrod','darkgoldenrod']],
             # 'p_cont':[path2, path2+'seaice/', 2, [994,967],
             #           ['violet','mediumorchid','darkviolet']],
             }
    
    axes2[0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    
    # data loop, bar plot
    era_values = {}
    min_values = {}; max_values = {}
    for yi, years in enumerate(newyrs):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        if yi == 0: trial = trials
        elif yi == 1: trial = trials

        era_values[yi] = [[] for x in np.arange(0, len(trial.keys()))]
        min_values[yi] = []
        max_values[yi] = []
        
        for k, name in enumerate(list(trial.keys())):
            attr = trial[name]
            
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                seaice_lines(years, attr[0]+'census/', attr[0]+'area/', attr[1], 
                             AREA_IND=attr[2], p_thresh=attr[3][li], sample_size=0)

            for mi, mm in enumerate(months):
                rel_area = []
                for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                    ss = sia-sia_clim
                    rel_area.append(ss-ss[0])
                if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
                else: MEANLINE = np.nansum(rel_area, axis=0)/1e6
            
                if len(lines[mi+1])<10: 
                    era_values[yi][k].append([np.nan]*3)
                    continue 
            
                if not shift: xvals = mi+np.array(XSPACING)
                else:
                    xvals = mi+np.array(XSPACING) + 6
                    if xvals[0]>=12: xvals -= 12
                
                diffs = [MEANLINE[ii] for ii in INDS]
                era_values[yi][k].append(diffs)
                
                # plot bars/ comparison dots
                if name == 'original':
                    if yi==0:
                        axes2[0].bar(xvals[0:3], diffs, width=WIDTH,
                                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
                    else:
                        axes2[0].bar(xvals[3:], diffs, width=WIDTH,
                                facecolor=hemi_colors2[li], alpha=0.5, 
                                edgecolor=hemi_colors2[li], hatch=HATCH, lw=2)


        #### find min/max values after each trial is calculated
        for mi, mm in enumerate(months):
            if not shift: xvals = mi+np.array(XSPACING)
            else: 
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
            
            # plot spread bars
            diff_arr = np.array([era_values[yi][k][mi] for k in np.arange(0, len(trial))])
            diff_arr=np.squeeze(diff_arr)
            if set_minmax:
                min_arr = diff_arr[set_trials[0]]
                max_arr = diff_arr[set_trials[1]]
            else:
                min_arr = np.nanmin(diff_arr, axis=0)
                max_arr = np.nanmax(diff_arr, axis=0)
                
            if yi==0:
                axes2[0].errorbar(xvals[0:3], diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                   fmt='none', capsize=capsize, ecolor=eb_color[li])
            else:
                axes2[0].errorbar(xvals[3:], diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                   fmt='none', capsize=capsize, ecolor=eb_color[li])
                
            min_values[yi].append(min_arr)
            max_values[yi].append(max_arr)
                

        ### TOTAL
        sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(3)] )
        if yi==0:
            axes2[0].bar(xvals[0:3]+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
            xvo = xvals[0:3]
        else:
            axes2[0].bar(xvals[3:]+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], hatch=HATCH, lw=2)
            xvo = xvals[3:]
            
        min_sums = np.array( [np.nansum(np.array(min_values[yi])[:,dr]) for dr in range(3)] )
        max_sums = np.array( [np.nansum(np.array(max_values[yi])[:,dr]) for dr in range(3)] )
        axes2[0].errorbar(xvo+1.33+li, sums, yerr=[np.abs(min_sums-sums), np.abs(max_sums-sums)],
                               fmt='none', capsize=capsize, ecolor=eb_color[li])
       
        axes2[0].axvline(1.1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
    
#### sum
for yi, years in enumerate(decades):
    sums = []
    for vi, VAL in enumerate([values, valmin, valmax]): # each hemi
        val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
        val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
        sums.append( val0 + val1 )
    
        if vi==0:
            heights = []
            for mi, sum1 in enumerate(sums[0][0]):
                xvals = mi+np.array(XSPACING)
                heights.append(sum1)
                
                    # original
                fcs = [hemi_colors2[0][i] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
                ecs = [hemi_colors[0] if np.abs(val0[0][mi][i])>np.abs(val1[0][mi][i]) else hemi_colors[1] for i in [0,1,2]]
                
                if yi==0:
                    axes2[1].bar(xvals[0:3], sum1, width=WIDTH,
                            facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)
                else:
                    axes2[1].bar(xvals[3:], sum1, width=WIDTH,
                            facecolor=fcs, alpha=0.5, edgecolor=ecs, hatch=HATCH,lw=2)

            # total sum
            if yi==0:
                axes2[1].bar(xvals[0:3]+1.33+0.5, np.nansum(sums[0][0], axis=0), width=WIDTH,
                        facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
                print('* total sum values')
                for i, ind in enumerate(INDS):
                    print('-','day',str(ind-7)+':', round(np.nansum(sums[0][0], axis=0)[i],3))
            else:
                axes2[1].bar(xvals[3:]+1.33+0.5, np.nansum(sums[0][0], axis=0), width=WIDTH,
                        facecolor=shades, alpha=0.5, edgecolor='k', hatch=HATCH, lw=2)
                print('-','80s:', round(np.nansum(sums[0][0], axis=0)[-1],3))
            
            axes2[1].axvline(1.1+xvals[0], ls='--', color='k', lw=0.85)
            
        # add spread bars!
        elif vi==1 and yi==1:
            val0_min = val0
            val1_min = val1
        elif vi==2:
            for mi in np.arange(0,12):
                if yi==0: xvals = mi+np.array(XSPACING[0:3])
                elif yi==1: xvals = mi+np.array(XSPACING[3:])
                
                axes2[1].errorbar(xvals, heights[mi], 
                                yerr=[np.abs(sums[1][mi]-heights[mi]), np.abs(sums[2][mi]-heights[mi])],
                                fmt='none', capsize=capsize, ecolor='#262626')
                
            # total sum 
            ht1 =  np.nansum(sums[0][0], axis=0)
            axes2[1].errorbar(xvals+1.33+0.5,  ht1, 
                            yerr=[np.abs(np.nansum(sums[1], axis=0)-ht1), np.abs(np.nansum(sums[2], axis=0)-ht1)],
                            fmt='none', capsize=capsize+0.5, ecolor='k')
            
            print('* total spread')
            print(ht1*100/np.nansum(sums[1], axis=0))

# x-axis
axes2[0].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
axes2[1].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
for ax in axes2.flatten():
    ax.xaxis.set_tick_params(labelbottom=True)
    ax.yaxis.set_tick_params(labelleft=True)
    
    
#-------- storm counts
ax3 = axes_all[-1]
hemi_colors3 = [['#018571','#80cdc1'],['#a6611a','#dfc27d']]

for li, loc in enumerate(['Arctic','Antarctic']):
    for yi, years in enumerate([np.arange(2010,2015), np.arange(2015,2019)]):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/', sample_size=0)


        ax3.plot(months-1+(WIDTH*3), [len(lines[mm]) for mm in months], lw=1.75,
                   color=hemi_colors3[li][yi], marker = ['o','^'][yi], markersize=8,
                   label = loc+': '+yr_title)

ax3.set_title('Monthly Storm Counts', fontsize=fs+2)

ax3.set_xlim(-0.5,14.25)
ax3.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1])+(WIDTH*2)))
ax3.yaxis.set_tick_params(labelleft=True)
ax3.tick_params(axis='both', labelsize=fs)
ax3.text(0.01, 1.025, '('+alphbar[i]+')',transform=ax3.transAxes, 
          fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                              'edgecolor':'white', 'lw':0.75},zorder=50)

ax3.set_ylabel('Number of Storms', fontsize=fs)
ax3.set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)], fontsize=fs)
ax3.xaxis.set_tick_params(labelbottom=True)
ax3.yaxis.set_tick_params(labelleft=True)

ax3.legend(loc='upper right', fontsize=fs-1, handletextpad=0.33, handlelength=0.75)

'''divide bars by storm counts...'''


#%% daily era5 data? - skip


#%% miz vs total area for case study...
def seaice_series(storm_event, miz_points, in_bbox):
    dt1 = storm_event[0]-timedelta(days=7)
    dt2 = storm_event[0]+timedelta(days=14)
    analysis_range = fx.daterange(dt1, dt2, dt=24)
    
    si = []; clim = []
    for dt in analysis_range:
        sic = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day, latlon=False)
        sic_in = np.ma.masked_array(sic, mask=in_bbox).filled(np.nan)
        sic_miz = np.ma.masked_where(miz_points<1, sic_in).filled(np.nan)
        si.append(np.nansum(sic_miz*25*25))
        
        decade_si = []
        for yy in np.arange(2010, 2020):
            sic = fx.load_seaice_sh(ice_fname, yy, dt.month, dt.day, latlon=False)
            sic_in = np.ma.masked_array(sic, mask=in_bbox).filled(np.nan)
            sic_miz = np.ma.masked_where(miz_points<1, sic_in).filled(np.nan)
            decade_si.append(np.nansum(sic_miz*25*25))
        clim.append(np.nanmean(decade_si))
        
    return np.array(si), np.array(clim)

def full_series(storm_event, in_bbox):
    dt1 = storm_event[0]-timedelta(days=7)
    dt2 = storm_event[0]+timedelta(days=14)
    analysis_range = fx.daterange(dt1, dt2, dt=24)
    
    si = []; clim = []
    for dt in analysis_range:
        sic = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day, latlon=False)
        sic_in = np.ma.masked_array(sic, mask=in_bbox).filled(np.nan)
        si.append(np.nansum(sic_in*25*25))
        
        decade_si = []
        for yy in np.arange(2010, 2020):
            sic = fx.load_seaice_sh(ice_fname, yy, dt.month, dt.day, latlon=False)
            sic_in = np.ma.masked_array(sic, mask=in_bbox).filled(np.nan)
            decade_si.append(np.nansum(sic_in*25*25))
        clim.append(np.nanmean(decade_si))
        
    return np.array(si), np.array(clim)

def non_miz_series(storm_event, miz_points, in_bbox):
    dt1 = storm_event[0]-timedelta(days=7)
    dt2 = storm_event[0]+timedelta(days=14)
    analysis_range = fx.daterange(dt1, dt2, dt=24)
    
    si = []; clim = []
    for dt in analysis_range:
        sic = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day, latlon=False)
        sic_in = np.ma.masked_array(sic, mask=in_bbox).filled(np.nan)
        sic_miz = np.ma.masked_where(miz_points==1, sic_in).filled(np.nan)
        si.append(np.nansum(sic_miz*25*25))
        
        decade_si = []
        for yy in np.arange(2010, 2020):
            sic = fx.load_seaice_sh(ice_fname, yy, dt.month, dt.day, latlon=False)
            sic_in = np.ma.masked_array(sic, mask=in_bbox).filled(np.nan)
            sic_miz = np.ma.masked_where(miz_points==1, sic_in).filled(np.nan)
            decade_si.append(np.nansum(sic_miz*25*25))
        clim.append(np.nanmean(decade_si))
        
    return np.array(si), np.array(clim)

    
def get_miz_area(storm_event):
    miz=[0.15,0.80]
    t1 = storm_event[0] - timedelta(days=1)
    t2 = storm_event[-1] + timedelta(days=1)
    storm_range = fx.daterange(t1, t2, dt=24)
    
    miz_points = np.zeros(np.shape(si_lon))
    for date in storm_range:
        sic = fx.load_seaice_sh(ice_fname, date.year, date.month, date.day, latlon=False)
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
        
    return miz_points

def check_lims(arr, lon, lat, minlon, maxlon, minlat, maxlat):
    arr = np.where(lon<minlon, np.nan, arr)
    arr = np.where(lon>maxlon, np.nan, arr)
    arr = np.where(lat<minlat, np.nan, arr)
    arr = np.where(lat>maxlat, np.nan, arr)
    return arr

#%%% storm info
storm_event = fx.daterange(datetime(2011,4,15), datetime(2011,4,18), dt=24)

ncname = datetime(2011,4,15).strftime('%Y_%m%d')+'_contours.nc'
with xr.open_dataset(path1+'contours/'+ncname) as cs:
    all_coords = []
    for key in list(cs.keys()):
        coord = cs[key].values
        date = datetime.strptime(key.split('_')[1]+key.split('_')[2], '%Y%m%d')
        if key.split('_')[-2] == '980':
            all_coords.append(coord)
    
# plot bbox
minlat = np.nanmin([np.nanmin(cc[:,1]) for cc in all_coords])
maxlat = np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords])
minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
bbox1 = fx.get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False)
inside = fx.find_points_in_contour(bbox1, si_lon, si_lat) 

miz_points = get_miz_area(storm_event)
arr1 = np.ma.masked_array(np.ones(np.shape(si_lon)), mask=~inside).filled(np.nan)
arr1 = np.where(miz_points<1, np.nan, arr1)
miz_pts_in = check_lims(miz_points, si_lon, si_lat, minlon, maxlon, minlat, maxlat)

# change in sic
si_start, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[0].year, storm_event[0].month, storm_event[0].day, latlon=True)
si_end = fx.load_seaice_sh(ice_fname, storm_event[-1].year, storm_event[-1].month, storm_event[-1].day, latlon=False)
change = np.where(si_end-si_start==0, np.nan, si_end-si_start)

# storm series
si, clim = seaice_series(storm_event, miz_points, ~inside)

# full series
si_full, clim_full = full_series(storm_event, ~inside)

# non miz (niz)
si_niz, clim_niz = non_miz_series(storm_event, miz_points, ~inside)

#%%% data plot

fig = plt.figure(figsize=(6.5,7.5)) 
gs = GridSpec(nrows=2, ncols=1, height_ratios=[1,0.45])
# fig.suptitle('\n'+storm_event[-1].strftime('%Y %b %d')+' minus '+storm_event[0].strftime('%d'), fontsize=12)

#### map
ax1 = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
ax1.set_title('Change in Sea Ice Concentration\n'+storm_event[-1].strftime('%Y %B %d')+' minus '+storm_event[0].strftime('%d'), fontsize=12)
ax1.set_extent([-60,60,-59,-90], ccrs.PlateCarree())
ax1.gridlines(draw_labels=False)
ax1.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
ax1.add_feature(cfeature.LAKES, facecolor='0.85', zorder=1)
ax1.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=2)
ax1.plot(bbox1[:,0], bbox1[:,1], transform=ccrs.PlateCarree(), color='navy')
pcm = ax1.pcolormesh(si_lon, si_lat, change,
                    cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, zorder=-5,
                    transform=ccrs.PlateCarree())
ax1 = fx.plot_geocontour(ax1, si_lon, si_lat, miz_pts_in, levels=[0.99], color='k', lw=2, ls='solid')

# colorbar
cax3 = fig.add_axes([0.265,0.425,0.5,0.015]) 
cbar3 = fig.colorbar(pcm, cax=cax3, orientation='horizontal')
cax3.set_xlabel('Change in Sea Ice Concentration', fontsize=9)
cax3.tick_params(labelsize=9) #

#### timeseries
ax2 = fig.add_subplot(gs[1])
ax2.set_title('Sea Ice Area Time Series')
ax2.axvline(0, lw=1, color='gray', ls=':')
ax2.axvline(len(storm_event), lw=1, color='gray', ls=':')
ax2.set_xlim([-7,14])
ax2.set_xticks(np.arange(-7,14+1,1))
ax2.set_xticklabels(xlabels, minor=False, rotation=0)
ax2.set_xlabel('Days Since Storm Start')
ax2.set_ylabel(r'$\times 10^5$ km$^2$')

ax2.plot(xxx, si/1e5, color='navy', ls='-', lw=1.5, label='MIZ Area')
# ax2.plot(xxx, clim/1e5, color='navy', ls='--', lw=1.25, label='2010-2019')

ax2.plot(xxx, si_full/1e5, color='maroon', ls='-', lw=1.5, label='Full Storm Area')
# ax2.plot(xxx, clim_full/1e5, color='maroon', ls='--', lw=1.25, label='2010-2019')

ax2.plot(xxx, si_niz/1e5, color='darkgreen', ls='-', lw=1.5, label='Non-MIZ Storm Area')
# ax2.plot(xxx, clim_niz/1e5, color='darkgreen', ls='--', lw=1.25, label='2010-2019')

ax2.legend(loc='lower left', handletextpad=0.5, handlelength=1.25, ncol=3,
           bbox_to_anchor=(0.025, -0.55))

#%%%% statistics

start_frac = si[7]/si_full[7]
end_frac = si[11]/si_full[11]
final_frac = si[-1]/si_full[-1]
print('starting percent:', round(start_frac*100,1))
print('ending percent:', round(end_frac*100,1))
print('final percent:', round(final_frac*100,1))

print('original difference:', round((si-clim)[11] - (si-clim)[7], 2))

full_change = si_full[11]- si_full[7]
miz_change = si[11]-si[7]
niz_change = si_niz[11]-si_niz[7]
print('full difference:', round(full_change,1))
print('miz difference:', round(miz_change,1))
print('niz difference:', round(niz_change,1))
print('percent storm change:', round(miz_change*100/full_change, 5),
      round(niz_change*100/full_change, 5))


#%% decade storm differences
# see [storm_strength.py] for individual analysis
import pickle

def check_significance(data1, data2):
    from scipy.stats import ttest_ind, ks_2samp, mannwhitneyu

    t, pt = ttest_ind(data1, data2)
    ks, pks = ks_2samp(data1, data2)
    mw, pmw = mannwhitneyu(data1, data2)
     
    pvals = [pt, pks, pmw]
    conf_levels = [100*(1-p) for p in pvals]
    
    return np.array(conf_levels)

def set_box_color(bp, color, lw=1.5):
    plt.setp(bp['boxes'], color=color, lw=lw)
    plt.setp(bp['whiskers'], color=color, lw=lw)
    plt.setp(bp['caps'], color=color, lw=lw)
    plt.setp(bp['medians'], color=color, lw=lw)
    plt.setp(bp['fliers'], color=color)

#%%% combined plot...
decade_colors = ['#2166ac', '#b2182b']

o=1.25
fig = plt.figure(figsize=(16/o, 15/o))
subfigs = fig.subfigures(1, 2, wspace=-0.1)

axes_all = [subfigs[0].subplots(5, 1), subfigs[1].subplots(5, 1)]

for li, loc in enumerate(['Arctic','Antarctic']):
    path1 = root_paths[li]

    axes = axes_all[li]

    plt.subplots_adjust(hspace=0.4)
    subfigs[li].suptitle('\n\n\n'+loc, fontsize=14, fontweight='bold')
    
    if li==1:
        for yi, years in enumerate(decades):
            ystr = str(years[0])+'-'+str(years[-1])
            axes[0].plot([],[], lw=4, color=decade_colors[yi], label=ystr)
        axes[0].plot([],[], lw=4, color='gray', label='Significant Difference')
        leg = axes[0].legend(loc='upper left', handletextpad=0.5, handlelength=0.5, ncol = 1, fontsize=9)
        leg.get_frame().set_edgecolor('white')
        
    #### storm counts
    axes[0].set_ylabel('Number of Storms')
    axes[0].set_title('Annual Storm Counts', fontsize=11)
    compare = []
    for yi, years in enumerate(decades):
        ystr = str(years[0])+'-'+str(years[-1])
        flierprops = dict(marker='o', markerfacecolor=decade_colors[yi], markersize=4,
                  linestyle='none', markeredgecolor=decade_colors[yi], alpha=0.5)
        
        storm_counts = []
        for year in years:
            start, end, pressure = fx.storm_pressure([year], path1+'census/', path1+'area/')
            storm_counts.append([len(start[mm]) for mm in months])
        compare.append(np.array(storm_counts))

        # mean_counts = np.nanmean(storm_counts, axis=0)
        # std_counts = np.nanstd(storm_counts, axis=0)
        # axes[0].plot(months, mean_counts, 
        #         marker='o', color = decade_colors[yi], label = ystr)
        # axes[0].errorbar(months, mean_counts, yerr=std_counts, 
        #             color = decade_colors[yi], capsize=5)
        
        for month in months:
            bp = axes[0].boxplot(compare[yi][:,month-1], positions=[month+(yi*0.33)-0.15],
                                  flierprops=flierprops)
            set_box_color(bp, decade_colors[yi], lw=1)
        
    ymin, ymax = axes[0].get_ylim()
    step = 0.075*ymax
    for mm in months:
        cl = check_significance(compare[0][:,mm-1], compare[1][:,mm-1])
        if np.any(np.array(cl)>=95):
            axes[0].axvspan(mm-0.25,mm+0.25, color='gray', alpha=0.25, zorder=-1)
            for i, (c, t) in enumerate(zip(cl, ['t', 'KS', 'MW'])):
                if c>=95: axes[0].text(mm-0.5, ymax-(step*(i+1)), t+':'+str(round(c,1)), fontsize=7.5)

    
    
    #### SLP
    axes[1].set_ylabel('Pressure (hPa)')
    axes[1].set_title('Minimum Storm Sea Level Pressure', fontsize=11)
    compare = []
    for yi, years in enumerate(decades):
        ystr = str(years[0])+'-'+str(years[-1])
    
        start, end, pressure = fx.storm_pressure(years, path1+'census/', path1+'area/')
    
        pressure_values = [pressure[mm] if len(pressure[mm])>10 else [np.nan] for mm in months]
        compare.append(pressure_values)
    
        # pressure_mean = [np.nanmean(pressure[mm]) if len(pressure[mm])>10 else np.nan for mm in months]
        # pressure_std = [np.nanstd(pressure[mm]) if len(pressure[mm])>10 else np.nan for mm in months]
        # axes[1].plot(months, pressure_mean, 
        #         marker='o', color = decade_colors[yi], label = ystr)
        # axes[1].errorbar(months, pressure_mean, yerr=pressure_std, 
        #             color = decade_colors[yi], capsize=5)
        
        for month in months:
            if np.all(~np.isnan(pressure_values[month-1])):
                bp = axes[1].boxplot(pressure_values[month-1], positions=[month+(yi*0.33)-0.15],
                                      flierprops=flierprops)
                set_box_color(bp, decade_colors[yi], lw=1)
        
    ymin, ymax = axes[1].get_ylim()
    step = 0.05*ymax
    for mm in months:
        if len(compare[yi][mm-1])<10: continue
        cl = check_significance(compare[0][mm-1], compare[1][mm-1])
        if np.any(np.array(cl)>=95):
            axes[1].axvspan(mm-0.25,mm+0.25, color='gray', alpha=0.25, zorder=-1)
            for i, (c, t) in enumerate(zip(cl, ['t', 'KS', 'MW'])):
                if c>=95: axes[1].text(mm-0.5, ymax-(step*(i+1)), ' '+t+':'+str(round(c,1)), fontsize=7.5)

    
    
    #### ACE max
    # [ace_area, ace_max, mke1, mke2]
    titles = ['$ACE_{area}$','$ACE_{max}$','$MKE$','$MKE_{2}$']
    units = '$m^2 s^{-2}$'  
    full_names = ['Accumulated Cyclone Energy (area)',
                  'Accumulated Cyclone Energy (maximum)',
                  'Mean Kinetic Energy', 'Mean Kinetic Energy']
    
    for ax1, idx in zip([axes[2],axes[3],axes[4]],[1,0,2]):
        
        ax1.set_ylabel('{}'.format(titles[idx])+' [{}]'.format(units))
        ax1.set_title(full_names[idx], fontsize=11)
    
        compare = []
        for yi, years in enumerate(decades):
            flierprops = dict(marker='o', markerfacecolor=decade_colors[yi], markersize=4,
                      linestyle='none', markeredgecolor=decade_colors[yi], alpha=0.5)
            
            data = {mm:[] for mm in months}
            
            for year in years:
                out_data = pickle.load(open(path1+'wind/energy_'+str(year)+'.pkl', 'rb'))
                for mm in months:
                    data[mm] += [sublist[idx] for sublist in out_data[mm]]
                        
            ### plot
            for month in months:
                if li==1 and month in [1,2,3]: continue
                bp = ax1.boxplot(data[month], positions=[month+(yi*0.33)-0.15],
                                      flierprops=flierprops)
                set_box_color(bp, decade_colors[yi], lw=1)
                
            compare.append(data)
    
        ymin, ymax = ax1.get_ylim()
        step = 0.075*ymax
        for mm in months:
            if li==1 and mm in [1,2,3]: continue
            cl = check_significance(compare[0][mm], compare[1][mm])
            if np.any(cl >= 95): 
                ax1.axvspan(mm-0.25,mm+0.25, color='gray', alpha=0.25, zorder=-1)
                for i, (c, t) in enumerate(zip(cl, ['t', 'KS', 'MW'])):
                    # if mm==3 and c>=95: 
                    #     ax1.text(mm-0.3, ymax-(step*(i+1)), t+':'+str(round(c,1)), fontsize=7.5)
                    if mm==1 and c>=95: 
                        ax1.text(mm-0.5, ymax-(step*(i+1))+40, t+':'+str(round(c,1)), fontsize=7.5)
                    elif idx==1 and mm==6 and c>=95: 
                        ax1.text(mm-0.5, ymax-(step*(i+1))+50, t+':'+str(round(c,1)), fontsize=7.5)
                    elif c>=95:
                        ax1.text(mm-0.5, ymax-(step*(i+1)), t+':'+str(round(c,1)), fontsize=7.5)

    # axes, each loc column
    for ax in axes.flatten():
        ax.set_xlim([0,13])
        ax.set_xticks(months)
        ax.set_xticklabels(months) 
    axes[-1].set_xlabel('Month')
    
# alph = iter(list(string.ascii_lowercase))
# for i, ax in enumerate(axes_all[0]):
#     for ax1 in [ax, axes_all[1][i]]:
#         ax1.text(0.0, 1.05, '('+next(alph)+')', 
#                  transform=ax1.transAxes, fontsize=fontsize-2, 
#                  zorder=50)

alph = iter(list(string.ascii_lowercase))
for li in [0,1]:
    for ax1 in axes_all[li]:
        ax1.text(0.0, 1.05, '('+next(alph)+')', 
                 transform=ax1.transAxes, fontsize=fontsize-2, 
                 zorder=50)


#%% end

