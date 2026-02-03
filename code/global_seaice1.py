#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 2025
global_seaice1.py

- anual sea ice cycle
- global cylone impacts

@author: mundi
"""

#%% imports and files
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import time as timeIN
from glob import glob

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

#%% data organization

fs = 14

alph = ['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k']

month_colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum


# sea ice decade shades

name1 = 'ice_area'
name2 = 'extent_' # 'miz_'
decades_all = [np.arange(1982,1992), np.arange(1990, 2000),
           np.arange(2000, 2010), np.arange(2010,2020)]

sh_colors = ['#dfc27d','#bf812d','#8c510a','#543005']
nh_colors = ['#80cdc1','#35978f','#01665e','#003c30']
hemi_shades = [nh_colors, sh_colors]
shades = ['#bdbdbd','#969696','#636363','#252525']


#%% global sea ice extent

#### set up plot
fig2, axes2 = plt.subplots(3,1, figsize=(10,12), sharex=True)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.set_xlim(0,365)
    ax2.set_xticks([datetime(2010,mm,1).timetuple().tm_yday-1 for mm in months])
    ax2.set_xticklabels([calendar.month_abbr[mm] for mm in months], minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2: ax2.set_ylabel('Ice Extent '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)

#### organize data
for era, years in enumerate(decades_all):
    yrstr = str(years[0])+'-'+str(years[-1])

    area_hemi = {}
    annual_means = []
    for loc_ind, loc in enumerate(hemi_names):
        area_hemi[loc_ind] = []
        for yi, year in enumerate(years):
            areas = np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy')
            if len(areas)>365:
                areas2 = []
                for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                    if dt!=datetime(year,2,29): 
                        areas2.append(areas[di])
                areas=areas2
                
            area_hemi[loc_ind].append(np.array(areas)/1e6)
        annual_means.append(np.nanmean(area_hemi[loc_ind], axis=0))
        
    # top plot
    axes2[0].set_title('Annual Mean Ice Extent', fontweight='bold', fontsize=fs)
    for loc_ind, loc in enumerate(hemi_names):
        axes2[0].plot(np.arange(0,365),annual_means[loc_ind], color=hemi_shades[loc_ind][era], 
                      lw=2.5) #label = loc+' '+yrstr
        
    # bottom plot
    axes2[1].set_title('Global Sum', fontweight='bold', fontsize=fs)
    annual_mean = annual_means[0] + annual_means[1]
    annual_std = np.nanstd(area_hemi[0]+area_hemi[1], axis=0)
    axes2[1].plot(np.arange(0,365), annual_mean, color=shades[era], lw=2.5, label=yrstr)
    leg1 = axes2[1].legend(fontsize=fs-1, handletextpad=0.5, handlelength=1.25)
    
# organize legends
for loc_ind, loc in enumerate(hemi_names):
    axes2[0].plot([],[], color=hemi_colors[loc_ind], label=loc, alpha=0.75, lw=5)
axes2[0].legend(fontsize=fs-1, handletextpad=0.5, handlelength=1.25, loc='upper left')

for line in leg1.get_lines():
    line.set_linewidth(5)
    
#-----------------------
#### add storm counts!
#-----------------------
ax1 = axes2[-1]
ax1.set_title('Storm Counts', fontweight='bold', fontsize=fs)
ax1.set_ylabel('Number of Storms',fontsize=fs+1)
xvals = [datetime(2010,mm,1).timetuple().tm_yday-1 for mm in months]

for li, loc in enumerate(['NH','SH']):
    for yi, years in enumerate([decades[0], decades[-1]]):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')

        counts = [len(start_day[mm]) if len(start_day[mm])>10 else np.nan for mm in months]
        
        ax1.plot(xvals, counts, hemi_colors[li], alpha=0.5+(0.5*yi), marker='o', lw=2.5,
                 label = loc+' '+yr_title)
    
ax1.legend(loc='upper left', ncol=1, fontsize=fs-2, 
           handletextpad=0.5, handlelength=1.25,)

#%%% full time series

fullx=np.arange(0,365*10)
#### set up plot
fig2, axes2 = plt.subplots(2,1, figsize=(12,6), sharex=True)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.set_xlim(0,fullx[-1])
    # ax2.set_xticks([datetime(2010,mm,1).timetuple().tm_yday-1 for mm in months])
    # ax2.set_xticklabels([calendar.month_abbr[mm] for mm in months], minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2: ax2.set_ylabel('Ice Extent '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)

#### organize data
for era, years in enumerate(decades_all):
    yrstr = str(years[0])+'-'+str(years[-1])

    area_hemi = {}
    for loc_ind, loc in enumerate(hemi_names):
        area_hemi[loc_ind] = []
        for yi, year in enumerate(years):
            areas = np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy')
            if len(areas)>365:
                areas2 = []
                for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                    if dt!=datetime(year,2,29): 
                        areas2.append(areas[di])
                areas=areas2
                
            area_hemi[loc_ind] += list(np.array(areas)/1e6)
        
    # top plot
    axes2[0].set_title('Annual Mean Ice Extent', fontweight='bold', fontsize=fs)
    for loc_ind, loc in enumerate(hemi_names):
        axes2[0].plot(fullx, area_hemi[loc_ind], color=hemi_shades[loc_ind][era], 
                      lw=2.5) #label = loc+' '+yrstr
        
    # bottom plot
    axes2[1].set_title('Global Sum', fontweight='bold', fontsize=fs)
    annual_sum = np.array(area_hemi[0]) + np.array(area_hemi[1])
    axes2[1].plot(fullx, annual_sum, color=shades[era], lw=2.5, label=yrstr)
    
slope, intercept, r, p, se = linregress(fullx, annual_sum)
axes2[1].plot(fullx, (slope*fullx)+intercept, color=shades[era], ls='--', lw=2.5, label='m='+str(slope*365))
leg1 = axes2[1].legend(loc='lower left', fontsize=fs-1, handletextpad=0.5, handlelength=1.25,
                       ncol=5, bbox_to_anchor=(-0.05,-0.45))



#%%% compute sea ice gradient
# dx = 0.1; y = [1, 2, 3, 4, 4, 5, 6] # dx constant
# np.gradient(y, dx) # dy/dx 2nd order accurate
# array([10., 10., 10.,  5.,  5., 10., 10.])

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
            

#### set up plot
fig2, axes2 = plt.subplots(2,1, figsize=(10,8), sharex=True)
fig2.suptitle('With derivatives plotted', fontstyle='italic')
for i, ax2 in enumerate(axes2.flatten()):
    ax2.set_xlim(0,365)
    ax2.set_xticks([datetime(2010,mm,1).timetuple().tm_yday-1 for mm in months])
    ax2.set_xticklabels([calendar.month_abbr[mm] for mm in months], minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2: ax2.set_ylabel('Ice Extent '+r'($\times 10^6$ km$^2$)',fontsize=fs+1)

ax0_d = axes2[0].twinx()
ax1_d = axes2[1].twinx()
for axd in [ax0_d, ax1_d]: axd.axhline(0, lw=0.5, color='k')

#### organize data
for era, years in enumerate(decades):
    
    yrstr = str(years[0])+'-'+str(years[-1])

    area_hemi = {}
    annual_means = []
    for loc_ind, loc in enumerate(hemi_names):
        area_hemi[loc_ind] = []
        for yi, year in enumerate(years):
            areas = np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy')
            if len(areas)>365:
                areas2 = []
                for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                    if dt!=datetime(year,2,29): 
                        areas2.append(areas[di])
                areas=areas2
                
            area_hemi[loc_ind].append(np.array(areas)/1e6)
        annual_means.append(np.nanmean(area_hemi[loc_ind], axis=0))
        
    # top plot
    axes2[0].set_title('Annual Mean Ice Extent', fontweight='bold', fontsize=fs)
    
    for loc_ind, loc in enumerate(hemi_names):
        axes2[0].plot(np.arange(0,365),annual_means[loc_ind], color=hemi_shades[loc_ind][era], 
                      lw=2.5, alpha=0.5) #label = loc+' '+yrstr
        
        # ax0_d.plot(np.arange(0,365), np.gradient(annual_means[loc_ind], 1),
        #            color=hemi_shades[loc_ind][era], lw=2.5, ls='--')
        
        ax0_d.plot(np.arange(0,365), smooth(np.gradient(annual_means[loc_ind], 1),30),
                   color=hemi_shades[loc_ind][era], lw=2.5, ls='--')
        
    # bottom plot
    axes2[1].set_title('Global Sum', fontweight='bold', fontsize=fs)
    global_mean = annual_means[0] + annual_means[1]
    axes2[1].plot(np.arange(0,365), global_mean, color=shades[era], lw=2.5, label=yrstr)
    leg1 = axes2[1].legend(fontsize=fs-1, handletextpad=0.5, handlelength=1.25)

    
    # ax1_d.plot(np.arange(0,365), np.gradient(global_mean, 5),
    #            color=shades[era], lw=2.5, ls='--')
    
    ax1_d.plot(np.arange(0,365), smooth( np.gradient(global_mean, 5), 30),
               color=shades[era], lw=2.5, ls='--')

#%% summed annual change + ice deriv

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
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    
fig2.subplots_adjust(wspace=0.1)

#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    axes2[0][0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[0][0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    era_values = {}
    for yi, years in enumerate(decades):
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[0][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            
            if len(lines[mi+1])<10: era_values[yi].append([np.nan]*3); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.66, edgecolor=hemi_colors[li], lw=2)

        ### TOTAL
        sums = [np.nansum(np.array(era_values[yi])[:,dr]) for dr in range(len(xvals))]
        axes2[0][yi].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                facecolor=hemi_colors2[li], alpha=0.66, edgecolor=hemi_colors[li], lw=2)
        
        if not shift: axes2[0][yi].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        axes2[0][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    
#### sum
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums = val0 + val1
    
    for mi, sum1 in enumerate(sums):
        xvals = mi+np.array(XSPACING)
        
        fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
        ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
        
        axes2[1][yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                facecolor=fcs, alpha=0.66, edgecolor=ecs, lw=2)

    axes2[1][yi].bar(xvals+1.33+0.5, np.nansum(sums, axis=0), width=XSPACING[1]-XSPACING[0],
            facecolor=shades, alpha=0.66, edgecolor='k', lw=2)
    axes2[1][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
    
#### add sea ice derivatives
for yi, years in enumerate(decades):

    # get sea ice data
    area_hemi = {}
    annual_means = []
    for loc_ind, loc in enumerate(hemi_names):
        area_hemi[loc_ind] = []
        for year in years:
            areas = np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy')
            if len(areas)>365:
                areas2 = []
                for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                    if dt!=datetime(year,2,29): 
                        areas2.append(areas[di])
                areas=areas2
            area_hemi[loc_ind].append(np.array(areas)/1e6)
        annual_means.append(np.nanmean(area_hemi[loc_ind], axis=0))
        
    # top plots
    
    ax0d = axes2[0][yi].twinx()
    ax1d = axes2[1][yi].twinx()
    for axd in [ax0d, ax1d]:
        axd.set_ylim(-0.15,0.15)
        if yi==1: 
            axd.set_ylabel('Sea Ice Extent Derivative '+r'km$^2$ day$^{-1}$',fontsize=fs+1)
            axd.tick_params(axis='both', labelsize=fs)
        else: axd.set_yticklabels([])
   
    for loc_ind, loc in enumerate(hemi_names):
        ax0d.plot(np.linspace(0,12,365), smooth(np.gradient(annual_means[loc_ind],1),30),
                   color=hemi_colors[loc_ind], lw=2, ls='-', alpha=0.5, zorder=-100)
        
    # bottom plot
    global_mean = annual_means[0] + annual_means[1]
    ax1d.plot(np.linspace(0,12,365), smooth( np.gradient(global_mean,1), 30),
               color='dimgray', lw=2, ls='-', alpha=0.5, zorder=-100,
               label = 'Change in Ice Extent')
    if yi==0:
        ax1d.legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    
#%%% summed annual change!

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
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    
fig2.subplots_adjust(wspace=0.1)

#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    axes2[0][0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[0][0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    era_values = {}
    for yi, years in enumerate(decades):
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[0][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            
            if len(lines[mi+1])<10: era_values[yi].append([np.nan]*3); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)

        ### TOTAL
        sums = [np.nansum(np.array(era_values[yi])[:,dr]) for dr in range(len(xvals))]
        axes2[0][yi].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        
        if not shift: axes2[0][yi].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        axes2[0][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    
#### sum
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums = val0 + val1
    
    for mi, sum1 in enumerate(sums):
        xvals = mi+np.array(XSPACING)
        
        fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
        ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
        
        axes2[1][yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)

    axes2[1][yi].bar(xvals+1.33+0.5, np.nansum(sums, axis=0), width=XSPACING[1]-XSPACING[0],
            facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
    axes2[1][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
    
#%% differences in ice latitude...
# ACC

from concurrent.futures import ThreadPoolExecutor
import time
from functools import wraps
import gc

def timethis(func):
    """ 
    Print the execution time for a function call
    """
    @wraps(func)
    def wrapped_method(*args, **kwargs):
        time_start = time.time()
        output = func(*args, **kwargs)
        time_end = time.time()
        if time_end-time_start < 120:
            print(f"{func.__name__}: {(time_end-time_start)} s")
        else:
            print(f"{func.__name__}: {(time_end-time_start)/60} min")

        return output

    return wrapped_method

def get_date_list(loc_ind, year):
    date_list = []
    for month in np.arange(1,12+1):
        for day in np.arange(1,calendar.monthrange(year, month)[-1]+1):
            if month==2 and day==29: continue
            dt = datetime(year, month, day)
            date_list.append([loc_ind, dt])
    return date_list

def get_latitude(loc_date):
    loc_ind = loc_date[0]
    date = loc_date[1]
    
    if loc_ind==0:
        si, si_lon, si_lat = fx.load_seaice(ice_fname, date.year, date.month, date.day)
        fig, ax = fx.background_plot_nh(returnfig=True)
    elif loc_ind==1:
        si, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', date.year, date.month, date.day)
        fig, ax = fx.background_plot_sh(returnfig=True)
    
    contours = ax.contour(si_lon, si_lat, si, levels=[0.15],
                          colors='k', transform=ccrs.PlateCarree())
    plt.close(fig)
    cont = contours.allsegs[0]
    
    lats = []
    for cc in cont:
        for val in cc[:,1]:
            lats.append(val)
        
    return np.nanmean(lats)

@timethis 
def get_latitude_threaded(date_list):
    with ThreadPoolExecutor(max_workers=10) as executor:
        return executor.map(get_latitude, date_list)

@timethis
def get_latitudes_loop(date_list):
    lat_list = []
    for loc_date in date_list:
        loc_ind = loc_date[0]
        date = loc_date[1]
        if loc_ind==0:
            si, si_lon, si_lat = fx.load_seaice(ice_fname, date.year, date.month, date.day)
            fig, ax = fx.background_plot_nh(returnfig=True)
        elif loc_ind==1:
            si, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', date.year, date.month, date.day)
            fig, ax = fx.background_plot_sh(returnfig=True)
        
        try:
            contours = ax.contour(si_lon, si_lat, si, levels=[0.15],
                                  colors='k', transform=ccrs.PlateCarree())
            plt.close(fig)
            cont = contours.allsegs[0]
        except ValueError:
            print('--> '+date.strftime('%Y-%m-%d'))
            lat_list.append(np.nan)
            continue
        
        lats = []
        for cc in cont:
            for val in cc[:,1]:
                lats.append(val)
        lat_list.append(np.nanmean(lats))
    return lat_list

#%%% compare decades... mean ice extent latitude (+spread error bars)

#### set up plot
fig2, axes2 = plt.subplots(2,1, figsize=(10,8), sharex=True)
fig2.suptitle('Annual Mean Ice Extent Latitude', fontweight='bold',fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.set_xlim(0,365)
    ax2.set_xticks([datetime(2010,mm,1).timetuple().tm_yday-1 for mm in months])
    ax2.set_xticklabels([calendar.month_abbr[mm] for mm in months], minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2, direction, hemi in zip(axes2.flatten(), ['N','S'], hemi_names): 
    ax2.set_ylabel('Latitude '+r'($^\circ$'+direction+')',fontsize=fs)
    ax2.set_title(hemi, fontsize=fs)

#### organize data
for era, years in enumerate(decades_all):
    yrstr = str(years[0])+'-'+str(years[-1])
    print(yrstr)

    area_hemi = {}
    annual_means = []
    for loc_ind, loc in enumerate(hemi_names):
        area_hemi[loc_ind] = []
        for yi, year in enumerate(years):
            lat_file = root_paths[loc_ind]+'seaice/'+'lat_'+str(year)+'.npy'
            try: lat_series = np.load(lat_file)
            except FileNotFoundError:
                print('*', loc, year)
                date_list = get_date_list(loc_ind, year)
                # lat_futures = get_latitude_threaded(date_list)
                # lat_series = [x for x in lat_futures]
                
                lat_series = get_latitudes_loop(date_list)
                
                np.save(lat_file, lat_series)
                
            out1 = gc.collect()
            area_hemi[loc_ind].append(lat_series)
        annual_means.append(np.nanmean(area_hemi[loc_ind], axis=0))
        
    # top plot
    for loc_ind, loc in enumerate(hemi_names):
        axes2[loc_ind].plot(np.arange(0,365),annual_means[loc_ind], color=hemi_shades[loc_ind][era], 
                      lw=2.5) #label = loc+' '+yrstr
            
#%% more timescales (3,7,14)

INDS = [7+x+1 for x in np.arange(0,14)]
XSPACING = np.linspace(0,0.95, len(INDS))

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(25,10), sharex=True, sharey=True)
fig2.suptitle('\nMonthly Mean [DAILY] Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    
fig2.subplots_adjust(wspace=0.1)

#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    axes2[0][0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[0][0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    era_values = {}
    for yi, years in enumerate(decades):
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[0][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            
            if len(lines[mi+1])<10: era_values[yi].append([np.nan]*len(INDS)); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)

        ### TOTAL
        sums = [np.nansum(np.array(era_values[yi])[:,dr]) for dr in range(len(xvals))]
        axes2[0][yi].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        
        if not shift: axes2[0][yi].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        axes2[0][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    
#### sum

# line plot
fig, ax_lp = plt.subplots(1,1, figsize=(10, 5))
ax_lp.set_title('Daily Changes in MIZ Ice Area')
ax_lp.set_xticks(np.arange(0, len(INDS))+1)

# data loop
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums = val0 + val1
    
    for mi, sum1 in enumerate(sums):
        xvals = mi+np.array(XSPACING)
        
        fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
        ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
        
        axes2[1][yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)

    axes2[1][yi].bar(xvals+1.33+0.5, np.nansum(sums, axis=0), width=XSPACING[1]-XSPACING[0],
            facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
    axes2[1][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)  
    
    ddays = np.arange(0, len(np.nansum(sums, axis=0)))+1
    cchanges = np.nansum(sums, axis=0)
    ax_lp.plot(ddays, cchanges,
               color=['mediumblue', 'firebrick'][yi], marker='o',
               label=str(years[0])+'-'+str(years[-1]))
    
    max_change = np.max(cchanges)
    ax_lp.axvline(ddays[np.where(cchanges==max_change)], ls=':', lw=0.55,
                  color = ['mediumblue', 'firebrick'][yi])
    
ax_lp.legend(loc='upper left')

#%% timeseries of mean sic

try:
    timeseries2 = np.load('/Users/mundi/Desktop/month-hemi/mean_ice_sic.npy')
except FileNotFoundError:
    from concurrent.futures import ThreadPoolExecutor
    import time
    from functools import wraps
    
    def timethis(func):
        """ 
        Print the execution time for a function call
        """
        @wraps(func)
        def wrapped_method(*args, **kwargs):
            time_start = time.time()
            output = func(*args, **kwargs)
            time_end = time.time()
            if time_end-time_start < 120:
                print(f"{func.__name__}: {(time_end-time_start)} s")
            else:
                print(f"{func.__name__}: {(time_end-time_start)/60} min")
    
            return output
    
        return wrapped_method
    
    @timethis
    def get_date_list(years):
        date_list = []
        for year in years:
            for month in np.arange(1,12+1):
                for day in np.arange(1,calendar.monthrange(year, month)[-1]+1):
                    date_list.append(datetime(year, month, day))
        return date_list
        
    date_list = get_date_list(np.arange(1982,2020))
    
    def get_daily_sic(date):
        
        if date.month==1 and date.day==1: print(date.year)
    
        sic_nh = fx.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
        sic_sh = fx.load_seaice_sh(ice_fname+'south/', date.year, date.month, date.day, latlon=False)
        
        sic_nh = np.where(sic_nh <=0, np.nan, sic_nh) # mask ice-free areas
        sic_sh = np.where(sic_sh <=0, np.nan, sic_sh)
    
        return [np.nanmean(sic_nh), np.nanmean(sic_sh)]
    
    @timethis 
    def get_timeseries_threaded(date_list):
        with ThreadPoolExecutor(max_workers=10) as executor:
            return executor.map(get_daily_sic, date_list)
        
        
    tseries = get_timeseries_threaded(date_list)
    timeseries2 = np.array([x for x in tseries])
    
    np.save('/Users/mundi/Desktop/month-hemi/mean_ice_sic.npy', np.array(timeseries2))
        
#%%%% plot
fig, ax = plt.subplots(1,1)
for idx, hemi in enumerate(['Arctic', 'Antarctic']):
    ax.plot(timeseries2[:,idx], color = hemi_colors[idx])
    
    x2 = np.arange(0, len(timeseries2[:,idx]))
    mask = ~np.isnan(np.array(x2)) & ~np.isnan(timeseries2[:,idx])
    # m, b, r, p, se = linregress(np.array(x6[i])[mask], np.array(y6[i])[mask])
    
    slope, intercept, r, p, se = linregress(x2[mask], timeseries2[:,idx][mask])
    ax.plot(x2, (slope*x2)+intercept, color=hemi_colors[idx], ls='--', lw=2.5, label='m='+str(slope*365))
leg1 = ax.legend(loc='lower left', fontsize=fs-1, handletextpad=0.5, handlelength=1.25,
                   ncol=5, bbox_to_anchor=(-0.05,-0.45))
    


#%%xxxxxxxxxxxxxxxxxxxxxxxx
import warnings
def seaice_lines(years, census_path, area_path, si_path, AREA_IND, 
                 sample_size=10, p_thresh=984):
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


def align_zeros(axes):

    ylims_current = {}   #  Current ylims
    ylims_mod     = {}   #  Modified ylims
    deltas        = {}   #  ymax - ymin for ylims_current
    ratios        = {}   #  ratio of the zero point within deltas

    for ax in axes:
        ylims_current[ax] = list(ax.get_ylim())
                        # Need to convert a tuple to a list to manipulate elements.
        deltas[ax]        = ylims_current[ax][1] - ylims_current[ax][0]
        ratios[ax]        = -ylims_current[ax][0]/deltas[ax]
    
    for ax in axes:      # Loop through all axes to ensure each ax fits in others.
        ylims_mod[ax]     = [np.nan,np.nan]   # Construct a blank list
        ylims_mod[ax][1]  = max(deltas[ax] * (1-np.array(list(ratios.values()))))
                        # Choose the max value among (delta for ax)*(1-ratios),
                        # and apply it to ymax for ax
        ylims_mod[ax][0]  = min(-deltas[ax] * np.array(list(ratios.values())))
                        # Do the same for ymin
        ax.set_ylim(tuple(ylims_mod[ax]))
        
#%% paper figs: global bar
print(); print('*** GLOBAL BARS ***')
# suptitle('\Monthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)

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

#### data and plot

values = {}
valmin = {}; valmax = {}
wk2_change = {}
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
    for yi, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        if yi == 0: trial = {'original': trials['original']}
        elif yi == 1: trial = trials

        wk2_change[str(li)+'_'+str(yi)] = []
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
                    wk2_change[str(li)+'_'+str(yi)].append(0)
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

                    wk2_change[str(li)+'_'+str(yi)].append(diffs[-1])

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
                print('* total sum values')
                for i, ind in enumerate(INDS):
                    print('-','day',str(ind-7)+':', round(np.nansum(sums[0][0], axis=0)[i],3))
            else:
                axes2[1].bar(xvals[-1]+WIDTH+1.33+0.5, np.nansum(sums[0][0], axis=0)[-1], width=WIDTH,
                        facecolor='k', alpha=0.5, edgecolor='k', hatch=HATCH, lw=2)
                print('-','80s:', round(np.nansum(sums[0][0], axis=0)[-1],3))
            
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


print('* fractional contributions')
nh_cont = np.array( [np.nansum(np.array(values[0][yi][0])[:,dr]) for dr in range(len(xvals))] )
sh_cont = np.array( [np.nansum(np.array(values[1][yi][0])[:,dr]) for dr in range(len(xvals))] )
total = nh_cont + sh_cont

print('NH: '+str([round(x,2) for x in nh_cont]))
print('SH: '+str([round(x,2) for x in sh_cont]))
print('% NH contribution: '+str([round(x/total[i],2) for i,x in enumerate(nh_cont)]))
print('% SH contribution: '+str([round(x/total[i],2) for i,x in enumerate(sh_cont)]))

# -------------------------------------
#### add insolation curve to top plot
# -------------------------------------
import pickle
savepath = '/Users/mundi/Desktop/month-hemi/'

daily_values = pickle.load(open(savepath+'insolation_daily.pkl', 'rb'))
monthly_values = pickle.load(open(savepath+'insolation_monthly.pkl', 'rb'))

monthly_x = []
for month in np.arange(1,12+1):
    daynum_range = [datetime(year,month,1).timetuple().tm_yday,
                    datetime(year,month,calendar.monthrange(year,month)[-1]).timetuple().tm_yday]
    monthly_x.append(np.nanmean(daynum_range))
    

axi = axes2[0].twinx()
axi.tick_params(axis='y', labelsize=fs)
axi.set_ylabel(r'MIZ-Area-Weighted Insolation (W m$^{-2}$)', fontsize=fs+1)

for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    axi.plot(np.arange(XSPACING[1], 12+XSPACING[1]), np.nanmean(monthly_values[loc_ind], axis=0), 
             lw=1.5, ls=['-','--'][loc_ind], color=['#01665e','#8c510a'][loc_ind], label=loc, zorder=-5000)

axes2[0].set_ylim([-40,40])
axi.set_ylim([-515,515])
align_zeros([axes2[0], axi])

axi.set_yticks([0,150,300,450]);

# xtick labels (monthly)
if not shift: 
    # axes2[0].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
    axes2[0].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
    axes2[1].set_xticklabels([calendar.month_abbr[m] for m in np.arange(1,12+1)]+['Total'], fontsize=fs)
    for ax in axes2.flatten():
        ax.xaxis.set_tick_params(labelbottom=True)
        ax.yaxis.set_tick_params(labelleft=True)
else: 
    # axes2[0].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
    axes2[0].set_xticklabels([calendar.month_abbr[m] for m in np.arange(7,12+1)] +
                             [calendar.month_abbr[m] for m in np.arange(1,7)]+['Total'])



#%%% percent change from annual cycle
print(); print('* Percent change from annual cycle')

NAME = 'extent_' # 'miz_' #

for yi, years in enumerate(decades):
    yrstr = str(years[0])+'-'+str(years[-1])
    print(yrstr)
    
    #### caluclate global sum
    area_hemi = {}
    annual_means = []
    for loc_ind, loc in enumerate(hemi_names):
        area_hemi[loc_ind] = []
        for yx, year in enumerate(years):
            areas = np.load(root_paths[loc_ind]+'seaice/'+NAME+name1+'_'+str(year)+'.npy')
            if len(areas)>365:
                areas2 = []
                for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                    if dt!=datetime(year,2,29): 
                        areas2.append(areas[di])
                areas=areas2
                
            area_hemi[loc_ind].append(np.array(areas)/1e6)
        annual_means.append(np.nansum(area_hemi[loc_ind], axis=0)) # *sum* all area
        
    annual_mean = annual_means[0] + annual_means[1]
    
    annual_range = np.nanmax(annual_mean)-np.nanmin(annual_mean)
    print('- annual range: '+str(round(annual_range,3)))
    
    #### calculate storm sum
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums =  val0 + val1 
    SUM = np.nansum(sums[0], axis=0)[-1]
    print('- storm impact:', round(SUM,3))
    
    #### fraction
    print('-- Fraction:', round(SUM*100/annual_range,3),'%')
          
#%%%% monthly deviations?
print(); print('* MONTHLY Percent change from annual cycle')

NAME = 'miz_' #'extent_' #

decade_means=[]
for yi, years in enumerate(decades):
    yrstr = str(years[0])+'-'+str(years[-1])
    print(yrstr)
    
    #### caluclate global sum
    area_hemi = {}
    annual_means = []
    for loc_ind, loc in enumerate(hemi_names):
        area_hemi[loc_ind] = []
        for yx, year in enumerate(years):
            areas = np.load(root_paths[loc_ind]+'seaice/'+NAME+name1+'_'+str(year)+'.npy')
            if len(areas)>365:
                areas2 = []
                for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                    if dt!=datetime(year,2,29): 
                        areas2.append(areas[di])
                areas=areas2
                
            area_hemi[loc_ind].append(np.array(areas)/1e6)
        annual_means.append(np.nansum(area_hemi[loc_ind], axis=0)) # *sum* all area
        
    annual_mean = annual_means[0] + annual_means[1]
    decade_means.append(annual_mean)
    change = np.array(wk2_change['0_'+str(yi)][0:12])+np.array(wk2_change['1_'+str(yi)][0:12])

    
    #### monthly split
    mcount = 0
    for mi, mm in enumerate(np.arange(1,12+1)):
        print('&', calendar.month_name[mm])
        mdays = calendar.monthrange(2010, mm)[1]
        mvalues = annual_mean[mcount:mcount+mdays]
        monthly_range = np.nanmean(mvalues) #np.nanmax(mvalues)-np.nanmin(mvalues)
        mcount+=mdays

        print('- '+NAME+'range: '+str(round(monthly_range,3)))
        print('- storm impact:', round(change[mi],3))
        print('-- Fraction:', round(change[mi]*100/monthly_range,3),'%')
    
#%% end

