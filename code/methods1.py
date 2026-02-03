#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 2025
methods1.py

create a figure to describe paper methods

- cyclone min sea pressure comparison
- storm area, example storm
- timeseries/normalization procedure 

>> for each hemisphere

## add in storm/miz areas!!

@author: mundi
"""
#%% imports and files
import numpy as np
import xarray as xr
import calendar
from datetime import datetime, timedelta
import time as timeIN
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.gridspec import GridSpec
import functions as fx
import cmocean.cm as cmo


nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

ice_path = '/Users/mundi/Desktop/seaice/'
ice_paths = [ice_path, ice_path+'south/']

decades = [np.arange(1982, 1992), np.arange(2010,2020)]
decade_names = ['Early Satellite Era ', 'Present Day ']

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]
month_labels = [calendar.month_abbr[mm] for mm in months]+[calendar.month_abbr[months[0]]]

hemi_names= ['Arctic', 'Antarctic']

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

#%% background slp

colors= [['salmon','maroon'], ['lightsteelblue','steelblue']]

total_days = 0
month_ticks = [0]
for month in months:
    num_days = calendar.monthrange(2010, month)[1]
    total_days += num_days
    month_ticks.append(total_days)
    
xdays = np.arange(0, total_days, 1)
shifted_ticks = np.array(month_ticks[:-1])+15

#### plot

fig, axes = plt.subplots(2,1, figsize=(12,10))

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    ax=axes[loc_ind]
    ax.set_xlabel('Month')
    ax.set_ylabel('Pressure [hPa]')
    ax.set_title(loc+' Daily Mean Sea Level Pressure')
    
    
    for era, years in enumerate([np.arange(2010,2019+1), np.arange(1980, 1990)]):
        if loc_ind==0:
            if era==0: text_y = 999
            elif era==1: text_y = 997
        elif loc_ind==1:
            if era==0: text_y = 975
            elif era==1: text_y = 973
        
        all_years = []
        for yy, year in enumerate(years):
            try:
                annual_series = np.load(path1+'background_slp/'+str(year)+'_slp.npy')
            except FileNotFoundError:
                print(+'>> '+loc+' '+str(year)+'_slp.npy')
                continue
            
            if year%4==0: # leap years
                annual_series[31+28] = np.nanmean([annual_series[31+28],annual_series[31+29]])
                annual_series = list(annual_series)
                annual_series.pop(31+29)
            
            ax.plot(xdays, annual_series, color=colors[era][0])
            all_years.append(np.array(annual_series))
            
        ax.plot(xdays, np.nanmean(all_years, axis=0), color = colors[era][1], zorder=50, lw=2.5, 
                label=str(years[0])+'-'+str(years[-1])+'')
        
        # analyze mean values
        mean_series = np.nanmean(all_years, axis=0)
        chunks = []
        for sp, split in enumerate(month_ticks[:-1]):
            chunks.append(mean_series[split:month_ticks[sp+1]])
            
        for cc, chunk in enumerate(chunks):
            ax.text(shifted_ticks[cc]-5, text_y, round(np.nanmean(chunk),1), color = colors[era][1])
    
    ax.legend()
    for mt in month_ticks: ax.axvline(mt, color='gray', ls=':',zorder=-20)
    ax.set_xticks(shifted_ticks, month_labels[:-1])


#%% example storm
yy, mm = 2017, 10
dd = [23,22]

# yy, mm = 2017, 11
# dd = [7,8]

###

# yy, mm = 2018, 11
# dd = [12,11]

# yy, mm = 2017, 11
# dd = [7,8]

###

fig = plt.figure(figsize=(15,10)) 
gs = GridSpec(nrows=2, ncols=2, width_ratios=[1,1], height_ratios=[1,0.5])

# maps
axes1 = [fig.add_subplot(gs[0,0], projection=ccrs.NorthPolarStereo()),
        fig.add_subplot(gs[0,1], projection=ccrs.SouthPolarStereo())]

# timeseries
axes2 = [fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])]

for li, ax in enumerate(axes1):
    
    if li==0: ax = fx.setup_plot_nh(ax, labels=False)
    elif li==1: ax = fx.setup_plot_sh(ax, labels=False)
    
    ice_fname = ice_paths[li]
    path1 = root_paths[li]
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines([yy], path1+'census/', path1+'area/', path1+'seaice/')
        
    storm_areas = fx.get_storm_bbox([yy], path1+'census/', path1+'area/', path1+'contours/')

    #### storm area

    for start,end, bbox, si, clim, line in zip(start_day[mm],end_day[mm], storm_areas[mm],
                                         si_changes[mm], clim_changes[mm], lines[mm]):
        if start.day==dd[li]: break
    
    ax.plot(bbox[:,0], bbox[:,1], transform = ccrs.PlateCarree())
    ax.set_title(loc+': '+start.strftime('%Y %b %d')+'-'+end.strftime('%d'))
    
    # change in sea ice
    
    if li==0: 
        si_start, si_lon, si_lat = fx.load_seaice(ice_fname, yy, mm, dd[li], latlon=True)
        si_end = fx.load_seaice(ice_fname, yy, mm, end.day, latlon=False)
    elif li==1: 
        si_start, si_lon, si_lat = fx.load_seaice_sh(ice_fname, yy, mm, dd[li], latlon=True)
        si_end = fx.load_seaice_sh(ice_fname, yy, mm, end.day, latlon=False)
    
    pcm = ax.pcolormesh(si_lon, si_lat, np.where(si_end-si_start==0, np.nan, si_end-si_start),
                          cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, zorder=-5,
                          transform=ccrs.PlateCarree())
    
    # contours, location of min pressure
    stormstr1 = start.strftime('%Y_%m%d')
    ncname = stormstr1+'_contours.nc'
    with xr.open_dataset(path1+'contours/'+ncname) as cs:
        for key in list(cs.keys()):
            coord = cs[key].values
            ax.plot(coord[:,0], coord[:,1], color='k', lw=0.75,
                    transform=ccrs.PlateCarree())


    #### timeseries/normalization procedure 
    axes2[li].plot(xxx, si, color='k')
    axes2[li].plot(xxx, clim, color='b', ls='--')
    
    axes2[li].axhline(0, ls='-', color='k', lw=1)
    axes2[li].axvline(0, ls=':', color='gray', lw=1)
    axes2[li].set_xlim(-7,14)
    axes2[li].set_xticks(xxx)
    axes2[li].set_xticklabels(xlabels, minor=False)

    ax2 = axes2[li].twinx()
    ax2.plot(xxx, line, color='k', lw=3)

    fx.align_zeros([axes2[li], ax2])


# color bar
cax1 = fig.add_axes([0.495,0.5,0.025,0.33]) 
cbar1 = fig.colorbar(pcm, cax=cax1, orientation='vertical')
# cbar1.set_label(r'Change in Sea Ice Concentration', fontsize=8)
cax1.set_xlabel('\nChange in\nSea Ice\nConcentration')
cax1.tick_params(labelsize=10)




#%% FIND STORMS
# are there storms on the same day in each hemisphere? (or at lease close)

#%%%% storm lists
# may/june, oct/nov

sample_month = 11

for year in np.arange(2010,2020):
    print(); print(year); print()
    
    for li, loc in enumerate(['Arctic','Antarctic']):
        
        print(loc)
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines([year], path1+'census/', path1+'area/', path1+'seaice/')
            
        print(start_day[sample_month])
            
#%%%% good ones?

yy, mm = 2017, 10
dd = [23,24]

yy, mm = 2018, 11
dd = [12,11]

yy, mm = 2017, 11
dd = [7,8]

        
#%%%% sample maps

yy, mm = 2010, 11
dd = [26,26]

fig = plt.figure(figsize=(12,8)) 
gs = GridSpec(nrows=1, ncols=2, width_ratios=[1,1], height_ratios=[1])

axes = [fig.add_subplot(gs[0,0], projection=ccrs.NorthPolarStereo()),
        fig.add_subplot(gs[0,1], projection=ccrs.SouthPolarStereo())]

proj = [ccrs.NorthPolarStereo(),ccrs.SouthPolarStereo()]

for li, loc in enumerate(['Arctic','Antarctic']):
    
    ax = fig.add_subplot(gs[0,li], projection=proj[li])
    if li==0: ax = fx.setup_plot_nh(ax, labels=False)
    elif li==1: ax = fx.setup_plot_sh(ax, labels=False)
    
    path1 = root_paths[li]
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines([yy], path1+'census/', path1+'area/', path1+'seaice/')
        
    storm_areas = fx.get_storm_bbox([yy], path1+'census/', path1+'area/', path1+'contours/')

    for start,end, bbox in zip(start_day[mm],end_day[mm], storm_areas[mm]):
        if start.day!=dd[li]: continue
    
        ax.plot(bbox[:,0], bbox[:,1], transform = ccrs.PlateCarree())
        ax.set_title(loc+': '+start.strftime('%Y-%m %d')+'-'+end.strftime('%d'))
        
        break







#%% end
