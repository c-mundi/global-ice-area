#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 2025
locations1.py

% antarctica/start2.py
- analyze if variables are locally variable

@author: mundi
"""
#%% imports and filepaths
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import functions as fx
from datetime import datetime, timedelta
import string, calendar

import cmocean.cm as cmo
import cartopy.crs as ccrs
from scipy.stats import linregress
from matplotlib.gridspec import GridSpec

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

census_name = 'census_'
contour_name = '_contours'
si_name = '_seaice'

ice_fname = '/Users/mundi/Desktop/seaice/'
ice_fname = '/Users/mundi/Desktop/seaice/south/'

decades = [np.arange(2010,2020), np.arange(1982, 1992)]
decade_names = ['Early Satellite Era ', 'Present Day ']

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]


xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]


#%% constants

def lb_title(lon_bounds):
    title = ''
    for lb in lon_bounds[0:-1]:
        title += str(int(lb)) if lb<360 else str(int(lb-360))
        title += r'$^\circ$' +' '+ r'$\rightarrow$'+' '
    title += str(int(lon_bounds[-1])) if lon_bounds[-1]<360 else str(int(lon_bounds[-1]-360))
    title += r'$^\circ$' 
    return str(title)

import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

_, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010, 1, 1, latlon=True)
si_lon = np.where(si_lon<0, si_lon+360, si_lon)
ones = np.ones(np.shape(si_lon))

# lon_bounds = np.linspace(0,360,9)
lon_bounds = [60,150,210,300,420]
nbins = len(lon_bounds)-1

markersize=12
cmap = cmo.deep
norm = plt.Normalize(vmin=1, vmax=nbins)

hemi_names = ['Arctic', 'Antarctic']

alph = ['a','b','c','d','e','f']
colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

decade_colors = ['#b2182b', '#2166ac']
shade_colors = ['#ef8a62', '#67a9cf']

#%% location map
fig1, axm = fx.background_plot_sh(returnfig=True)
fig1.suptitle('\n\n'+lb_title(lon_bounds), fontsize=24)
axm.set_boundary(circle, transform=axm.transAxes)

for li, lb in enumerate(lon_bounds[0:-1]):
    lower = lb ;  upper = lon_bounds[li+1]
    if upper > 360:
        lower -= 360
        upper -= 360
        si_lon = np.where(si_lon>=180, si_lon-360, si_lon)
    else:
        si_lon = np.where(si_lon<0, si_lon+360, si_lon)
        
    one_sector = np.where(np.logical_and(si_lon>=lower, si_lon<upper), ones, np.nan)
    axm.pcolormesh(si_lon, si_lat, (li+1)*one_sector, cmap=cmap, vmin=1,vmax=nbins,
                  transform=ccrs.PlateCarree(), zorder=-700, alpha=0.25)

#%% functions

ice_lims = [20,80]
var_addons = [['_1000', '_990', '_contour'],['2_1000', '2_990', '2_contour']]
addon_num = 0 # storm area
sia_var_addons = var_addons
miz_ind = 1 #[0=daily, 1=total]
miz_names = ['daily_', 'total_']

def indiv_lines_lonsort(years, lon_bounds, LB, LI, census_path, area_path, si_path):
    storm_counts = {}
    lines = {}
    start_day = {}
    end_day = {}
    
    si_changes = {} 
    clim_changes = {}
    
    for idx in np.arange(0,12):
        lines[idx+1] = []
        start_day[idx+1] = []
        end_day[idx+1] = []
        si_changes[idx+1] = []
        clim_changes[idx+1] = []
    
    for year in years:
        storm_counts[year] = [0]*len(months)
        
        census_file = census_path+census_name+str(year)+'.csv'
        [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
        [[startlon, startlat],[endlon,endlat]] = fx.readCensus(census_file, convertDT=True)[1]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        for startdt, enddt in timing_grid:
            storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
        try:
            ds = xr.open_dataset(si_path + str(year) +'_seaice.nc')
        except:
            print('- skip: '+si_path + str(year) +'_seaice.nc')
            continue
        
        for storm_num, strm in enumerate(storm_ranges):
             month = int(strm[0].month)
             
             ### remove storms that don't interact with the ice
             ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
             if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                continue
             try:
                 sia = ds['sia_miz'+sia_var_addons[miz_ind][addon_num]].values[storm_num]
             except IndexError as EE:
                 print('-- sia error', EE)
                 continue
             
             x1 = startlon[storm_num] #!!! mean lon?
             if x1<0: x1+=360
             for li, lb in enumerate(lon_bounds):
                if x1 < lb: break
             if li == 0: li=len(lon_bounds)-1
             
             if lon_bounds[li-1] != LB: continue
    
             storm_counts[year][month-1] += 1 # mm -> index value
             
             # get time series                                              
             sia_clim = ds['sia_clim_miz'+sia_var_addons[miz_ind][addon_num]].values[storm_num]
             ss = sia-sia_clim
             standardized_area = (ss-ss[0])/(np.nanmax(ss)-np.nanmin(ss))
             pcd = (standardized_area)
        
             # plot
             lines[month].append(pcd)
             start_day[month].append(strm[0])
             end_day[month].append(strm[-1])
             
             si_changes[month].append(sia)
             clim_changes[month].append(sia_clim)
        
    mean_lines = []  
    for idx in np.arange(0,12):
        if len(lines[idx+1])!=0:
            mean_line = np.nanmean(lines[idx+1], axis=0)
            mean_lines.append(mean_line)
        else:
            mean_lines.append(np.nan*np.ones(np.shape(xxx)))
            continue
        
    return mean_lines, lines, start_day, end_day, si_changes, clim_changes


def ytext_scale(idx, startyear, endyear):
    yadd=0
    
    # if startyear==2010 and endyear==2019:
    #     if idx==7: #aug
    #         yadd = 0.025
    #     elif idx==8: # sep
    #         yadd = -0.025
    
    return yadd

#%% monthly lines (indiv plots for each loc)

path1 = root_paths[1]
for li, lb in enumerate(lon_bounds[0:-1]):

    #### monthly mean
    fig2, axes2 = plt.subplots(1,2, figsize=(16,6.5), sharey=True)
    for i, ax2 in enumerate(axes2):
        ax2.axhline(0, ls='-', color='k', lw=1)
        ax2.axvline(0, ls=':', color='gray', lw=1)
        ax2.set_xlim(-7,14)
        ax2.set_ylim([-0.95, 0.95])
        ax2.set_xticks(xxx)
        ax2.set_xticklabels(xlabels, minor=False, rotation=0)
        ax2.set_xlabel('Days Since Storm Start')
        ax2.set_ylabel('Normalized Relative Change in Ice Area')
        ax2.axvspan(-7,14, color=cmap(norm(li+1)), alpha=0.15, zorder=-500)
        ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
                  fontsize=12, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                      'edgecolor':'white', 'lw':0.75},zorder=50)
    axes2[1].yaxis.set_tick_params(labelleft=True)
    fig2.suptitle('Longitude: '+str(int(lb))+r'$^\circ$' +'-'+str(int(lon_bounds[li+1] if lon_bounds[li+1]<360 else lon_bounds[li+1]-360))+r'$^\circ$')
    
    
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            indiv_lines_lonsort(years, lon_bounds, lb, li, path1+'census/', path1+'area/', path1+'seaice/')
    
        ### monthly mean
        axes2[era].set_title('Change in MIZ Area For All Storms'+' '+yr_title)
        for idx, ml in enumerate(mean_lines):
            if len(lines[idx+1]) > 10: ###!!!
                axes2[era].plot(xxx, ml, color = colors[idx], label = calendar.month_name[idx+1])
            else: 
                axes2[era].plot([],[], color = colors[idx], label = calendar.month_name[idx+1])
                continue
            yadd = ytext_scale(idx, years[0], years[-1])
            axes2[era].text(xxx[-1]+0.25, ml[-1]+yadd, calendar.month_abbr[idx+1]) 
        
        axes2[0].legend(loc='lower left', ncol=2, handletextpad=0.5, handlelength=1,
                          edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%% do more side-by-side location comparisons
# bar plots (3,7,14-day change)

fig2, axes2 = plt.subplots(2, nbins, figsize=(16,12), sharey=True)

for loc_ind, ax_group2 in zip([0,1], axes2):
    path1 = root_paths[loc_ind]

    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        
        for li, lb in enumerate(lon_bounds[0:-1]):
            ax2 = ax_group2[li]
            ax2.axvspan(-1,12, color=cmap(norm(li+1)), alpha=0.15, zorder=-500)
        
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                indiv_lines_lonsort(years, lon_bounds, lb, li, path1+'census/', path1+'area/', path1+'seaice/')

            for mi, mean_line in enumerate(mean_lines):
                if len(lines[mi+1])<5: continue ###!!!
                    
                diffs = [mean_line[7+3], mean_line[7+7], mean_line[7+14]]
                ax2.bar(mi+np.array([0,0.25,0.5]), diffs, width=0.25,
                        color=decade_colors[era], alpha=0.5)

            ax2.set_xlim([-0.5,12])
            ax2.set_xticks(np.arange(0.5,12.5))
            ax2.set_xticklabels(labels = np.arange(1,12+1))

ycoords = [0.7, 0.3]
for hname, loc in zip(hemi_names,ycoords):
    fig2.text(0.0725, loc, hname, 
             va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')

for years, c in zip(decades, decade_colors):
    axes2[-1][0].plot([],[], color=c, label=str(years[0])+'-'+str(years[-1]), lw=5, alpha=0.5)
axes2[-1][0].legend(loc='lower left', handletextpad=0.5, handlelength=0.5,
                    edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%% end
