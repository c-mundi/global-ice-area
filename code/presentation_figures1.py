#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 2025
presentation_figures1.py

agu 2025

@author: mundi
"""

#%% imports and paths
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

loc_ind=1
path1= root_paths[loc_ind]

#%% world view cleaning
storm_event = fx.daterange(datetime(2011,4,15), datetime(2011,4,18), dt=24)

# set up plot
fig=plt.figure(figsize=[12,10]) 
ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.set_extent([-45,30,-59,-87], ccrs.PlateCarree())
ax.gridlines(draw_labels=False)
ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=1)
ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=2)
ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=3)

# circle
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
# ax.set_boundary(circle, transform=ax.transAxes)

# plot image
img_file = '/Users/mundi/Downloads/snapshot-2011-04-16T00_00_00Z.tif'
with gdal.Open(img_file) as ds:
    data = ds.ReadAsArray()
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()   
    extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
              gt[3] + ds.RasterYSize * gt[5], gt[3])
    img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)), extent=extent, origin='upper', zorder=-1)

# add sea ice contours
# cc = ['dimgray', 'k']
cc = ['#67a9cf', '#2166ac']
si1, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[0].year, storm_event[0].month, storm_event[0].day)
fx.plot_geocontour(ax, si_lon, si_lat, si1, [0.15], color=cc[0], lw=3)
si2, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[-1].year, storm_event[-1].month, storm_event[-1].day)
fx.plot_geocontour(ax, si_lon, si_lat, si2, [0.15], color=cc[1], lw=3)

for c1, dt in zip(cc, [storm_event[0], storm_event[-1]]):
    ax.plot([],[], color=c1, lw=5, label=dt.strftime('%m-%d-%Y'))
ax.legend(loc='lower right',handletextpad=0.5, handlelength=1, labelspacing=0.15, 
          fontsize=18 )

#%%% shaded miz
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

def get_storm_bbox_mask():
    ncname = datetime(2011,4,15).strftime('%Y_%m%d')+'_contours.nc'
    with xr.open_dataset(path1+'contours/'+ncname) as cs:
        all_coords = []
        for key in list(cs.keys()):
            coord = cs[key].values
            if key.split('_')[-2] == '980':
                all_coords.append(coord)
        
    # plot bbox
    minlat = np.nanmin([np.nanmin(cc[:,1]) for cc in all_coords])
    maxlat = np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords])
    minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
    maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
    bbox1 = fx.get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False)
    return fx.find_points_in_contour(bbox1, si_lon, si_lat) 


# set up plot
fig=plt.figure(figsize=[12,10]) 
ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.set_extent([-45,30,-59,-87], ccrs.PlateCarree())
ax.gridlines(draw_labels=False)
ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=1)
ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=2)
ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=3)

# circle
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
# ax.set_boundary(circle, transform=ax.transAxes)

# plot image
img_file = '/Users/mundi/Downloads/snapshot-2011-04-16T00_00_00Z.tif'
with gdal.Open(img_file) as ds:
    data = ds.ReadAsArray()
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()   
    extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
              gt[3] + ds.RasterYSize * gt[5], gt[3])
    img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)), extent=extent, origin='upper', zorder=-1)

# add miz shading
cc = ['#67a9cf', '#2166ac']
si1, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[0].year, storm_event[0].month, storm_event[0].day)
fx.plot_geocontour(ax, si_lon, si_lat, si1, [0.15], color=cc[0], lw=3)
si2, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[-1].year, storm_event[-1].month, storm_event[-1].day)
fx.plot_geocontour(ax, si_lon, si_lat, si2, [0.15], color=cc[1], lw=3)

miz_points = get_miz_area(storm_event)
inside = get_storm_bbox_mask()
arr1 = np.ma.masked_array(np.ones(np.shape(si_lon)), mask=~inside).filled(np.nan)
arr1 = np.where(miz_points<1, np.nan, arr1)

ax.pcolormesh(si_lon, si_lat, arr1, transform=ccrs.PlateCarree(),
              cmap=cmo.balance_r, vmin=-1,vmax=1, alpha=0.5)

for c1, dt in zip(cc, [storm_event[0], storm_event[-1]]):
    ax.plot([],[], color=c1, lw=5, label=dt.strftime('%m-%d-%Y'))
ax.legend(loc='lower right',handletextpad=0.5, handlelength=1, labelspacing=0.15, 
          fontsize=18 )

#%% global plots
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


#%%% separate annual sea ice/count plots

#### set up plot
# fig2, axes2 = plt.subplots(3,1, figsize=(10,12), sharex=True)

axes2 = []
for idx in range(3):
    fig2, ax = plt.subplots(1,1, figsize=(10,4), sharex=True)
    axes2.append(ax)

for i, ax2 in enumerate(axes2):
    ax2.set_xlim(0,365)
    ax2.set_xticks([datetime(2010,mm,1).timetuple().tm_yday-1 for mm in months])
    ax2.set_xticklabels([calendar.month_abbr[mm] for mm in months], minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    # ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
    #           fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
    #                               'edgecolor':'white', 'lw':0.75},zorder=50)
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

#%%% global bar
'''seaice_sens#.py'''
from paper_figures1 import seaice_lines

print(); print('*** GLOBAL BARS ***')

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
fig2, axes2 = plt.subplots(2,1, figsize=(10,10), sharey=True) #, sharex=True)
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
    ax2.set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)

axes2[0].set_title('MIZ Ice Area Changes in Each Hemisphere', fontsize=fs+2)
axes2[1].set_title('Global Sum', fontsize=fs+2)

import matplotlib.patches as mpatches
circ1 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,
                       label=str(decades[1][0])+'-'+str(decades[1][-1]))
circ2 = mpatches.Patch( facecolor='gray', edgecolor='k', alpha=0.5,hatch=HATCH,
                       label=str(decades[0][0])+'-'+str(decades[0][-1]))

ax2.legend(handles = [circ2,circ1], loc='lower left',
           fontsize=fs-1, handletextpad=0.5, handlelength=1.5)
    
    
fig2.subplots_adjust(hspace=0.3)

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


#%% end
