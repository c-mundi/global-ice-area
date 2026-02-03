#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 2025
seaice_sens4.py

- test sensitivy to different storm areas
- 1000 hpa, 990 hpa - bounding boxes, 1000 hpa contour area

- more recent decade!
* storm area
* increased pressure threshold

- plot to look like original figure

@author: mundi
"""
#%% imports and files
import numpy as np
import matplotlib.pyplot as plt
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress

import xarray as xr
from datetime import timedelta
import warnings

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

#%%% appearance

fontsize=11

AREA_NAMES = ['1000 hPa Bounding Box','990 hPa Bounding Box','1000 hPa Contour']
AREA_COLORS = ['maroon', 'navy', 'green']

# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

fs=14
XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]
alphbar = ['a','b','c','d']

#%%% functions
def seaice_lines(years, census_path, area_path, si_path, AREA_IND, sample_size=10, p_thresh=984):
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
            ds = xr.open_dataset(si_path + str(year) + '_seaice' + '.nc')
        except:
            print('- skip: '+si_path + str(year) + '_seaice' + '.nc')
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

#%%% settings

# mean cyclone impact vs summed across all storms
plot_mean = False

# normalized impact vs total chaneg in area
norm_lines = False

# months 1->12 vs shifted SH seasonal cycle
shift = False 

# use same trials for min/max for all months vs use each month's min/max
set_minmax = True
set_trials = [3,1]


# errorbar 
capsize=2
eb_color= hemi_colors #['k','k']

#%% summed annual change!

XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True, sharey=True)
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: 
    if plot_mean: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    else: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^6$ km$^2$)',fontsize=fs+1)
    
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
    
    # if li==0:
    #     for key in trial.keys():
    #         axes2[0][0].plot([],[], lw=6,marker='o', color=trial[key][-1][-1], label=key)
    #     axes2[0][0].legend(loc='upper left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    axes2[1][0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[1][0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    
    #### data loop, bar plot
    era_values = {}
    min_values = {}; max_values = {}
    for yi, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        axes2[li][yi].set_title(yr_title, fontsize=fs+2)
        
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
                    axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                            facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
                else:
                    for xa, ya, iii in zip(xvals, diffs, np.arange(0, len(xvals))):
                        # axes2[li][yi].plot(xa, ya, color=attr[-1][iii], marker='o')
                        pass
            
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
                    
                axes2[0][yi].errorbar(xvals, diff_arr[0], yerr=[np.abs(min_arr-diff_arr[0]), np.abs(max_arr-diff_arr[0])], 
                                       fmt='none', capsize=capsize, ecolor=eb_color[li])
                min_values[yi].append(min_arr)
                max_values[yi].append(max_arr)
                
                # # plot which trial min/max
                # tk = list(trial.keys())
                
                # for arg, ht in zip([np.argmin(diff_arr, axis=0), np.argmax(diff_arr, axis=0)],[[-75,-200],[125,300]]):
                #     for index, xa, iii in zip(arg, xvals, [0,1,2]):
                #         if np.all(~np.isnan(diff_arr)):
                #             axes2[0][yi].plot(xa, ht[li], color=trial[tk[index]][-1][iii], marker='o')

        ### TOTAL
        sums = np.array( [np.nansum(np.array(era_values[yi][0])[:,dr]) for dr in range(len(xvals))] )
        axes2[0][yi].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        
        if yi==1:
            min_sums = np.array( [np.nansum(np.array(min_values[yi])[:,dr]) for dr in range(len(xvals))] )
            max_sums = np.array( [np.nansum(np.array(max_values[yi])[:,dr]) for dr in range(len(xvals))] )
            axes2[0][yi].errorbar(xvals+1.33+li, sums, yerr=[np.abs(min_sums-sums), np.abs(max_sums-sums)],
                                   fmt='none', capsize=capsize, ecolor=eb_color[li])
        
        if not shift: axes2[0][yi].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        axes2[0][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
    
#### sum
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    
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
                
                axes2[1][yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                        facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)

            # total sum
            axes2[1][yi].bar(xvals+1.33+0.5, np.nansum(sums[0][0], axis=0), width=XSPACING[1]-XSPACING[0],
                    facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
            axes2[1][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
            
        #### add spread bars!
        elif vi==1 and yi==1:
            val0_min = val0
            val1_min = val1
        elif vi==2 and yi==1:
            for mi in np.arange(0,12):
                xvals = mi+np.array(XSPACING)
                
                axes2[1][yi].errorbar(xvals, heights[mi], 
                                yerr=[np.abs(sums[1][mi]-heights[mi]), np.abs(sums[2][mi]-heights[mi])],
                                fmt='none', capsize=capsize, ecolor='#262626')
                
            # total sum 
            ht1 =  np.nansum(sums[0][0], axis=0)
            axes2[1][yi].errorbar(xvals+1.33+0.5,  ht1, 
                            yerr=[np.abs(np.nansum(sums[1], axis=0)-ht1), np.abs(np.nansum(sums[2], axis=0)-ht1)],
                            fmt='none', capsize=capsize+0.5, ecolor='k')
    


#%% combined decades

XSPACING = [0,0.2,0.4]
INDS = [7+3, 7+7, 7+14]
WIDTH = XSPACING[1]-XSPACING[0]
HATCH = '/////'

#### set up
fig2, axes2 = plt.subplots(2,1, figsize=(10,10), sharex=True, sharey=True)
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    
    if plot_mean: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    else: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^6$ km$^2$)',fontsize=fs+1)
ax2.set_xlabel('Month',fontsize=fs)

import matplotlib.patches as mpatches
circ1 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,
                       label=str(decades[1][0])+'-'+str(decades[1][-1]))
circ2 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,hatch=HATCH,
                       label=str(decades[0][0])+'-'+str(decades[0][-1]))

ax2.legend(handles = [circ2,circ1], loc='lower left',
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
    
    
    #### data loop, bar plot
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
        
        if not shift: axes2[0].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
        else: axes2[0].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        axes2[0].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
    
#### sum
for yi, years in enumerate(decades):
    axes2[1].set_title('Sum', fontsize=fs+2)
    
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
            
        #### add spread bars!
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
    



#%% end
