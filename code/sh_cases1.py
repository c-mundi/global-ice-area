#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 2025
sh_cases1.py

# oct 11 [Updated AMS Polar talk?] email
--> find and prep case studies

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
from matplotlib.gridspec import GridSpec
import time as timeIN
from glob import glob

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

#%%% functions

def geoplot_bbox(ax, bbox, color='k', lw=2):   
    x = bbox[:,0]
    y = bbox[:,1]
    #do masked-array on the lon
    x_greater = np.ma.masked_greater(x, -0.01)
    x_lesser = np.ma.masked_less(x, 0)    
    # apply masks to other associate arrays: lat
    y_greater = np.ma.MaskedArray(y, mask=x_greater.mask)
    y_lesser = np.ma.MaskedArray(y, mask=x_lesser.mask)
    
    # plot
    for x1, y1 in zip([x_greater, x_lesser], [y_greater, y_lesser]):
        ax.plot(x1, y1, color=color, lw=lw, transform = ccrs.PlateCarree())
        
def get_miz_area(storm_event):
    miz = [0.15, 0.80]

    t1 = storm_event[0] - timedelta(days=1)
    t2 = storm_event[-1] + timedelta(days=1)
    storm_range = fx.daterange(t1, t2, dt=24)
    
    miz_points = np.zeros(np.shape(si_lon))
    for date in storm_range:
        sic = fx.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
   
    return miz_points

#%% spaghetti plot

month_groups = [[4,5], [11,12]]

#### set up plot
fig, axes_all = plt.subplots(len(month_groups), len(decades), 
                           figsize=(8,6.5), sharex=True, sharey=True)
fig.suptitle('Change in MIZ Area For All Storms', fontsize=fontsize+1)
for ax in axes_all[:,0]:
    ax.set_ylabel('Normalized Relative\nChange in Ice Area', fontsize=fontsize)
for ax in axes_all[1,:]:
    ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
    ax.set_xlabel('Days Since Storm Start', fontsize=fontsize-1)
alph = iter(list(string.ascii_lowercase))
for axl in axes_all:
    for ax1 in axl:
        ax1.set_xlim(-7,14)
        ax1.set_ylim(-1.05,1.05)
        ax1.axhline(0, ls='-', color='k', lw=1)
        ax1.axvline(0, ls='-', color='k', lw=0.75)
        ax1.set_xticks(xxx)
        ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
        ax1.text(0.0225, 0.915, next(alph), transform=ax1.transAxes, fontsize=fontsize, 
                bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                      'edgecolor':'k', 'lw':0.75},zorder=50)

#### data plotting
for era, years in enumerate(decades):
    ystr = str(years[0])+'-'+str(years[-1])
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
    for mi, month_group in enumerate(month_groups):
        
        all_lines = []
        for month in month_group:
            all_lines+=lines[month]

        axes_all[mi][era].plot(xxx, np.array(all_lines).T, lw=0.55, color='gray')
        
        axes_all[mi][era].plot(xxx, np.nanmean(all_lines, axis=0), lw=2, color='maroon')
        
        mnames = [calendar.month_abbr[mm]+', ' for mm in month_group]
        axes_all[mi][era].set_title(ystr+': '+mnames[0]+mnames[1][:-2])

#%%----
#%% case studies
clim_years=decades[1]

#### prints out list of possible storms

sample_month = 4

sample_year = 2011

print(); print(sample_year, 'possible storms'); print()

path1 = root_paths[1] # SH only
mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
    fx.indiv_lines([sample_year], path1+'census/', path1+'area/', path1+'seaice/')
    
print(start_day[sample_month])

#%%% possible cases

# 2011, 5 18-20
# 2014, 5 19-23
# 2018, 5 10-12
# 2019, 5 26-29

# 2010, 4 14-16 #https://worldview.earthdata.nasa.gov/?v=-6565710.223427985,-5449658.934634143,3662221.7004043013,1394524.1896547538&p=antarctic&l=Reference_Labels_15m(hidden),Reference_Features_15m(hidden),Coastlines_15m,AMSRUE_Sea_Ice_Concentration_12km(min=15,max=81),OCI_PACE_True_Color(hidden),VIIRS_NOAA21_CorrectedReflectance_TrueColor(hidden),VIIRS_NOAA20_CorrectedReflectance_TrueColor(hidden),VIIRS_SNPP_CorrectedReflectance_TrueColor(hidden),MODIS_Aqua_CorrectedReflectance_TrueColor(hidden),MODIS_Terra_CorrectedReflectance_TrueColor&lg=true&t=2010-04-15-T20%3A51%3A01Z
# 2011, 4 15-18 #https://worldview.earthdata.nasa.gov/?v=-5641594.715727818,-2595829.205198206,4079435.560021003,3909152.6288080346&p=antarctic&l=Reference_Labels_15m(hidden),Reference_Features_15m(hidden),Coastlines_15m,AMSRUE_Sea_Ice_Concentration_12km(min=15,max=81),OCI_PACE_True_Color(hidden),VIIRS_NOAA21_CorrectedReflectance_TrueColor(hidden),VIIRS_NOAA20_CorrectedReflectance_TrueColor(hidden),VIIRS_SNPP_CorrectedReflectance_TrueColor(hidden),MODIS_Aqua_CorrectedReflectance_TrueColor,MODIS_Terra_CorrectedReflectance_TrueColor&lg=true&t=2011-04-19-T20%3A51%3A01Z

#%%% example may
print();print('Plotting maps of possible storms')

for dt in start_day[sample_month]:
    yy = dt.year
    mm = dt.month
    dd = dt.day

    fig = plt.figure(figsize=(7.5,10)) 
    gs = GridSpec(nrows=2, ncols=1, height_ratios=[1,0.5])
    # maps
    ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
    # timeseries
    axes2 = fig.add_subplot(gs[1])
    
    #### map plotting    
    ax = fx.setup_plot_sh(ax, labels=False)
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines([yy], path1+'census/', path1+'area/', path1+'seaice/')
        
    storm_areas = fx.get_storm_bbox([yy], path1+'census/', path1+'area/', path1+'contours/')
    
    # storm area
    for start,end, bbox, si, clim, line in zip(start_day[mm],end_day[mm], storm_areas[mm],
                                         si_changes[mm], clim_changes[mm], lines[mm]):
        if start.day==dd: break
    ax.plot(bbox[:,0], bbox[:,1], transform = ccrs.PlateCarree())
    ax.set_title(start.strftime('%Y %b %d')+'-'+end.strftime('%d'))
    
    # change in sea ice
    si_start, si_lon, si_lat = fx.load_seaice_sh(ice_fname, yy, mm, dd, latlon=True)
    si_end = fx.load_seaice_sh(ice_fname, yy, mm, end.day, latlon=False)
    
    change = np.where(si_end-si_start==0, np.nan, si_end-si_start)
    
    pcm = ax.pcolormesh(si_lon, si_lat, change,
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
    axes2.plot(xxx, si, color='k')
    axes2.plot(xxx, clim, color='b', ls='--')
    
    axes2.axhline(0, ls='-', color='k', lw=1)
    axes2.axvline(0, ls=':', color='gray', lw=1)
    axes2.set_xlim(-7,14)
    axes2.set_xticks(xxx)
    axes2.set_xticklabels(xlabels, minor=False)
    
    ax2 = axes2.twinx()
    ax2.plot(xxx, line, color='k', lw=3)
    
    fx.align_zeros([axes2, ax2])
    
    
    # color bar
    cax1 = fig.add_axes([0.1,0.5,0.025,0.33]) 
    cbar1 = fig.colorbar(pcm, cax=cax1, orientation='vertical')
    # cbar1.set_label(r'Change in Sea Ice Concentration', fontsize=8)
    cax1.set_xlabel('\nChange in\nSea Ice\nConcentration')
    cax1.tick_params(labelsize=10)



#%%% simultaneous storms?
print();print('Finding simultaneous storms and plotting maps and time series')

color_prev = 'green'
color1 = 'blue'

sample_month = 10

#np.arange(2010,2020)
for year in [2012]: #np.arange(2010,2020):
    print(year)
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines([year], path1+'census/', path1+'area/', path1+'seaice/')
        
    storm_areas = fx.get_storm_bbox([year], path1+'census/', path1+'area/', path1+'contours/')
        
    for idx, start in enumerate(start_day[sample_month]):
        end = end_day[sample_month][idx]
        
        storm_range = fx.daterange(start,end, dt=24)
        
        if idx==0: prev_storm_range = storm_range; continue
        
        # check for overlpas
        for dt in storm_range:
            if dt in prev_storm_range:
                prev_line = lines[sample_month][idx-1]
                line1 = lines[sample_month][idx]
                
                # check signs:
                if (prev_line[14]*line1[14])<0 or (prev_line[-1]*line1[-1])<0:
                    overlap = [dt1.day for dt1 in storm_range if dt1 in prev_storm_range]
                    if len(overlap)<2: continue
                
                    #### PLOT
                    ### map
                    fig = plt.figure(figsize=(7.5,10)) 
                    gs = GridSpec(nrows=2, ncols=1, height_ratios=[1,0.5])
                    ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
                    ax = fx.setup_plot_sh(ax, labels=False)
                    ax.set_title(str(overlap[-1])+' - '+str(overlap[0])+' '+dt.strftime('%b, %Y'))
                    
                    # change in sea ice
                    si_start, si_lon, si_lat = fx.load_seaice_sh(ice_fname, year, sample_month, overlap[0], latlon=True)
                    si_end = fx.load_seaice_sh(ice_fname, year, sample_month, overlap[-1], latlon=False)
                    change = np.where(si_end-si_start==0, np.nan, si_end-si_start)
                    pcm = ax.pcolormesh(si_lon, si_lat, change,
                                        cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, zorder=-5,
                                        transform=ccrs.PlateCarree())
                    
                    # storm areas
                    geoplot_bbox(ax, storm_areas[sample_month][idx-1], color=color_prev)
                    geoplot_bbox(ax, storm_areas[sample_month][idx], color=color1)
                    
                    # contours, location of min pressure
                    for dt2, color in zip([storm_range[0], prev_storm_range[0]],[color1, color_prev]):
                        ncname = dt2.strftime('%Y_%m%d')+'_contours.nc'
                        with xr.open_dataset(path1+'contours/'+ncname) as cs:
                            for key in list(cs.keys()):
                                coord = cs[key].values
                                ax.plot(coord[:,0], coord[:,1], color=color, lw=0.75,
                                        transform=ccrs.PlateCarree())
                    
                    ### timeseries
                    ax4 = fig.add_subplot(gs[1])
                    ax4.axhline(0, lw=0.5, color='gray', ls=':')
                    ax4.axvline(0, lw=0.5, color='gray', ls=':')
                    ax4.plot(xxx, prev_line, color=color_prev,
                             label=prev_storm_range[0].strftime('%Y-%m-%d')+'->'+prev_storm_range[-1].strftime('%d'))
                    ax4.plot(xxx, line1, color=color1,
                             label=storm_range[0].strftime('%Y-%m-%d')+'->'+storm_range[-1].strftime('%d'))
                    ax4.legend(loc='lower left')
                    
                    break
                
        # next storm comparison
        prev_storm_range = storm_range

#### list of "good" id'd storms
# 2016, 5 06,08 
# 2016, 5 16,17
# 2016, 5 27,27
# 2018, 5 20,24

# 2013, 6 03,04
# 2018, 6 16,18
# 2011, 8 17,22
# 2019, 8 08,10
# 2011 10 08,09
# 2012 10 25,28

# 2017 10 04,06

#%%% overlapping (consecutive storms)

print();print('Finding Consecutive Storms and plotting maps and time series')

color_prev = 'green'
color1 = 'blue'

sample_month = 9

#np.arange(2010,2020)
for year in [2017]: #np.arange(2010,2020):
    print(year, end=' ')
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines([year], path1+'census/', path1+'area/', path1+'seaice/')
        
    storm_areas = fx.get_storm_bbox([year], path1+'census/', path1+'area/', path1+'contours/')
        
    nstorms = 0
    for idx, start in enumerate(start_day[sample_month]):
        end = end_day[sample_month][idx]
        
        storm_range = fx.daterange(start-timedelta(days=7), start+timedelta(days=14), dt=24)
        
        if idx==0: prev_storm_range = storm_range; continue
        
        # check for overlpas
        for dt in storm_range[7:14]:
            if dt in prev_storm_range[10:]:
                
                # check loc:
                mlon_prev = np.nanmean(storm_areas[sample_month][idx-1][:,0])
                mlon_now = np.nanmean(storm_areas[sample_month][idx][:,0])
                if np.abs(mlon_prev-mlon_now)<50:
                    prev_line = lines[sample_month][idx-1]
                    line1 = lines[sample_month][idx]
                    
                    overlap = [dt1.day for dt1 in storm_range if dt1 in prev_storm_range]
                    
                    #### PLOT
                    ### map
                    fig = plt.figure(figsize=(7.5,10)) 
                    gs = GridSpec(nrows=2, ncols=1, height_ratios=[1,0.5])
                    ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
                    ax = fx.setup_plot_sh(ax, labels=False)
                    ax.set_title(str(overlap[-1])+' - '+str(overlap[0])+' '+dt.strftime('%b, %Y'))
                    
                    # change in sea ice
                    try: si_start, si_lon, si_lat = fx.load_seaice_sh(ice_fname, year, sample_month, overlap[0], latlon=True)
                    except:
                        if overlap[0]>30: si_start, si_lon, si_lat = fx.load_seaice_sh(ice_fname, year, sample_month-1, overlap[0], latlon=True)
                    si_end = fx.load_seaice_sh(ice_fname, year, sample_month, overlap[-1], latlon=False)
                    change = np.where(si_end-si_start==0, np.nan, si_end-si_start)
                    pcm = ax.pcolormesh(si_lon, si_lat, change,
                                        cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, zorder=-5,
                                        transform=ccrs.PlateCarree())
                    
                    # storm areas
                    geoplot_bbox(ax, storm_areas[sample_month][idx-1], color=color_prev)
                    geoplot_bbox(ax, storm_areas[sample_month][idx], color=color1)
                    
                    # contours, location of min pressure
                    for dt2, color in zip([storm_range[7], prev_storm_range[7]],[color1, color_prev]):
                        ncname = dt2.strftime('%Y_%m%d')+'_contours.nc'
                        with xr.open_dataset(path1+'contours/'+ncname) as cs:
                            for key in list(cs.keys()):
                                coord = cs[key].values
                                ax.plot(coord[:,0], coord[:,1], color=color, lw=0.75,
                                        transform=ccrs.PlateCarree())
                    
                    ### timeseries
                    ax4 = fig.add_subplot(gs[1])
                    ax4.axhline(0, lw=0.5, color='gray', ls=':')
                    ax4.axvline(0, lw=0.5, color='gray', ls=':')
                    ax4.plot(xxx, prev_line, color=color_prev,
                             label=prev_storm_range[0].strftime('%Y-%m-%d')+'->'+prev_storm_range[-1].strftime('%d'))
                    ax4.plot(xxx, line1, color=color1,
                             label=storm_range[0].strftime('%Y-%m-%d')+'->'+storm_range[-1].strftime('%d'))
                    ax4.legend(loc='lower left')
                    
                    ax4.set_title(start.strftime('%b %d')+', '+start_day[sample_month][idx-1].strftime('%b %d'))
                    
                    nstorms+=1
                    break
                    
        # next storm comparison
        prev_storm_range = storm_range

    print(nstorms)

#### list of "good" id'd storms (week before)
# 2011 09-25, 10-01
# 2017 08-31, 09-09 

#%%%% plot overlapping case...
plevel = 970 

case1_dt = datetime(2017, 9, 7)
case2_dt = datetime(2017, 9, 16) 

era_data_file = '/Users/mundi/Desktop/era_data/antarctic_2017_09_10_msl_wind.nc'

with xr.open_dataset(era_data_file) as ds:
    slp = ds['msl'].sel(valid_time=slice(case1_dt-timedelta(days=7), case2_dt+timedelta(days=14)))/100 #hpa
    slp_daily = slp.resample(valid_time='1D').mean(dim='valid_time')
    lon, lat = np.meshgrid(ds['longitude'].values, ds['latitude'].values)


print(); print('Plotting Consecutive maps and time series'); print()

storm_range1 = fx.daterange(case1_dt, datetime(2017, 9, 9), dt=24)
storm_range2 = fx.daterange(case2_dt, datetime(2017, 9, 18), dt=24)

clim_years = np.arange(2010,2020)

storm_colors = ['#e9a3c9','#a1d76a']

#### set up plots
fig = plt.figure(figsize=(19, 10))
fig.suptitle('Consecutive Storms', size='xx-large', fontweight='bold')
subfigs = fig.subfigures(1, 3)
subfigs[0].suptitle('\n\nFirst Storm', size='x-large')
subfigs[1].suptitle('\n\nSecond Storm', size='x-large')
subfigs[2].suptitle('\n\nNet Changes', size='x-large')

bboxs = []
miz_areas = []
for case_ind, storm_range in enumerate([storm_range1, storm_range2]):
    gs = GridSpec(nrows=3, ncols=1, height_ratios=[0.95,0.33, 0.33])
    dt_s = storm_range[0]
    dt_e = storm_range[-1]
    
    ##########
    #### map
    ##########
    ax = subfigs[case_ind].add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
    ax = fx.setup_plot_sh(ax, labels=False)
    ax.set_title(dt_s.strftime('%b %d')+'-'+dt_e.strftime('%d %Y'))
    
    # plot daily contours
    ncname = storm_range[0].strftime('%Y_%m%d')+'_contours.nc'
    with xr.open_dataset(path1+'contours/'+ncname) as cs:
        all_coords = []
        for key in list(cs.keys()):
            coord = cs[key].values
            all_coords.append(coord)
            ax.plot(coord[:,0], coord[:,1], color=color, lw=0.75,
                    transform=ccrs.PlateCarree())
        
    # plot bbox
    bbox = fx.get_bbox_edges(all_coords)
    bboxs.append(bbox)
    inside = fx.find_points_in_contour(bbox, si_lon, si_lat) ### inside bbox!!(~)
    geoplot_bbox(ax, bbox, color='magenta' if case_ind==0 else 'green')

    # plot change in sea ice
    si_start = fx.load_seaice_sh(ice_fname, dt_s.year, dt_s.month, dt_s.day, latlon=False)
    si_end = fx.load_seaice_sh(ice_fname, dt_e.year, dt_e.month, dt_e.day, latlon=False)
    diff = si_end-si_start
    diff = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), diff)
    ax.pcolormesh(si_lon, si_lat, np.ma.masked_array(diff, mask=inside), 
                  cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, transform=ccrs.PlateCarree() )
    
    ################
    #### timeseries
    ################
    full_time = fx.daterange(dt_s-timedelta(days=7), dt_s+timedelta(days=14), dt=24)
    
    # set up plot
    ax_t = subfigs[case_ind].add_subplot(gs[1])
    ax_t.axvline(0, lw=0.55, color='gray', ls=':')
    ax_t.set_xlim(-7,14)
    ax_t.set_xticks(xxx)
    ax_t.tick_params(axis='both', which='major')  
    ax_t.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
    
    # storm ranges
    axt2 = ax_t.twiny()
    xrange = [ft.day if ft.month==full_time[0].month else ft.day+31 for ft in full_time]
    axt2.set_xlim([xrange[0], xrange[-1]])
    axt2.set_xticks(xrange)
    axt2.set_xticklabels([ft.day for ft in full_time])
    if case_ind==0:
        axt2.axvspan(storm_range1[0].day+31, storm_range1[-1].day+1+31, color=storm_colors[0], alpha=0.33)
        axt2.axvspan(storm_range2[0].day+31, storm_range2[-1].day+1+31, color=storm_colors[1], alpha=0.33)
    else:
        axt2.axvspan(storm_range1[0].day, storm_range1[-1].day+1, color=storm_colors[0], alpha=0.33)
        axt2.axvspan(storm_range2[0].day, storm_range2[-1].day+1, color=storm_colors[1], alpha=0.33)
    
    # calculate lines
    miz_points = get_miz_area(storm_range); miz_areas.append(miz_points)
    tseries = []
    clim = []
    for dt in full_time:
        si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day)
        si_bbox = np.ma.masked_array(si, mask=inside)
        si_miz = np.ma.masked_where(miz_points==0, si_bbox).filled(np.nan)
        tseries.append(np.nansum(si_miz*25*25))
        daily = []
        for yy in clim_years:
            si = fx.load_seaice_sh(ice_fname, yy, dt.month, dt.day, latlon=False)
            si = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), si)
            si_bbox = np.ma.masked_array(si, mask=inside)
            si_miz = np.ma.masked_where(miz_points==0, si_bbox).filled(np.nan)
            daily.append(np.nansum(si_miz*25*25))
        clim.append(np.nanmean(daily))
    tseries = np.array(tseries); clim=np.array(clim)
    
    # plot timeseries
    ax_t.plot(xxx, tseries/1e6, color='k')
    ax_t.plot(xxx, clim/1e6, color='dimgray', ls='--')
    if case_ind==0: ax_t.set_ylabel(r'$\Delta$ MIZ Ice Area'+'\n'+r'($\times10^6$ km$^2$)')
    

    # normalized
    ax_n = subfigs[case_ind].add_subplot(gs[2])
    ax_n.axvline(0, lw=0.55, color='gray', ls=':')
    ax_n.axhline(0, lw=0.55, color='gray', ls=':')
    ax_n.set_xlim(-7,14)    
    ax_n.set_xticks(xxx)
    ax_n.tick_params(axis='both', which='major')  
    ax_n.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
    ax_n.set_xlabel('Days Since Storm Start', fontsize=fontsize-1)
    if case_ind==0: ax_n.set_ylabel('Normalized\n'+r'$\Delta$ MIZ Ice Area')
    
    
    rel = tseries-clim
    normalized = (rel-rel[0])/(np.max(rel)-np.min(rel))
    ax_n.plot(xxx, normalized, color='k')
    
    ax2 = ax_n.twiny()
    ax2.set_xlim([xrange[0], xrange[-1]])
    
    if case_ind==0:
        ax2.axvline(storm_range1[0].day+31, color=storm_colors[0], alpha=0.75)
        ax2.axvline(storm_range1[-1].day+1+31, color=storm_colors[0], alpha=0.75)
        ax2.axvline(storm_range2[0].day+31, color=storm_colors[1], alpha=0.75)
        ax2.axvline(storm_range2[-1].day+1+31, color=storm_colors[1], alpha=0.75)
    else:
        ax2.axvline(storm_range1[0].day, color=storm_colors[0], alpha=0.75)
        ax2.axvline(storm_range1[-1].day+1, color=storm_colors[0], alpha=0.75)
        ax2.axvline(storm_range2[0].day, color=storm_colors[1], alpha=0.75)
        ax2.axvline(storm_range2[-1].day+1, color=storm_colors[1], alpha=0.75)

    ax2.get_xaxis().set_visible(False)
    ax2.spines[['right', 'top','bottom','left']].set_visible(False)

####################################################
#### --> net changes
####################################################
gs = GridSpec(nrows=3, ncols=1, height_ratios=[0.95,0.33, 0.33])

#########
#### map
#########
ax = subfigs[-1].add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
ax = fx.setup_plot_sh(ax, labels=False)
ax.set_title(storm_range1[0].strftime('%b %d')+'-'+storm_range2[-1].strftime('%d %Y'))

# plot change in sea ice
si_start = fx.load_seaice_sh(ice_fname, storm_range1[0].year, storm_range1[0].month, storm_range1[0].day, latlon=False)
si_end = fx.load_seaice_sh(ice_fname, storm_range2[-1].year, storm_range2[-1].month, storm_range2[-1].day, latlon=False)
diff = si_end-si_start
# diff = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), diff)
# diff = np.ma.masked_array(diff, mask=inside)
ax.pcolormesh(si_lon, si_lat, diff, 
              cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, transform=ccrs.PlateCarree() )

# plot overlapped bbox
blons, blats = [],[]
for bb in bboxs:
    blons.append(np.nanmin(bb[:,0]))
    blons.append(np.nanmax(bb[:,0]))
    blats.append(np.nanmin(bb[:,1]))
    blats.append(np.nanmax(bb[:,1]))
    
blons.sort(); blats.sort()

bbox_net = fx.get_edge_lines(blons[1], blons[2], blats[1], blats[2])

inside3 = fx.find_points_in_contour(bbox_net, si_lon, si_lat) ### inside bbox!!(~)
geoplot_bbox(ax, bbox_net, color='k')


################
#### timeseries
################

full_time = fx.daterange(storm_range1[0]-timedelta(days=7), storm_range2[-1]+timedelta(days=14), dt=24)
xlabs = np.array([dt.day for dt in full_time])
xlabels_f = [str(val) if idx%2==0 else '' for idx, val in enumerate(xlabs)]
xvals = np.arange(0, len(full_time))

tseries = []
clim = []
for dt in full_time:
    si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day)
    si_bbox = np.ma.masked_array(si, mask=inside)
    si_miz = np.ma.masked_where(miz_areas[0]==0, si_bbox).filled(np.nan)
    si_miz = np.ma.masked_where(miz_areas[1]==0, si_miz).filled(np.nan)
    tseries.append(np.nansum(si_miz*25*25))
    daily = []
    for yy in clim_years:
        si = fx.load_seaice_sh(ice_fname, yy, dt.month, dt.day, latlon=False)
        si = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), si)
        si_bbox = np.ma.masked_array(si, mask=inside3)
        si_miz = np.ma.masked_where(miz_areas[0]==0, si_bbox).filled(np.nan)
        si_miz = np.ma.masked_where(miz_areas[1]==0, si_miz).filled(np.nan)
        daily.append(np.nansum(si_miz*25*25))
    clim.append(np.nanmean(daily))
tseries = np.array(tseries); clim=np.array(clim)

# set up plot
ax_t = subfigs[-1].add_subplot(gs[1])
ax_t.axvline(0, lw=0.55, color='gray', ls=':')
ax_t.set_xlim(xvals[0],xvals[-1])
ax_t.set_xticks(xvals)
ax_t.set_xticklabels(xlabels_f, fontsize=fontsize)
ax_t.tick_params(axis='both', which='major')  

# plot timeseries
ax_t.plot(xvals, tseries/1e6, color='k')
ax_t.plot(xvals, clim/1e6, color='dimgray', ls='--')

ax_t.axvspan(np.where(xlabs==storm_range1[0].day)[0][0], 
             (np.where(xlabs==storm_range1[-1].day)[0]+1)[0],
             color=storm_colors[0], alpha=0.33)
ax_t.axvspan(np.where(xlabs==storm_range2[0].day)[0][0], 
             (np.where(xlabs==storm_range2[-1].day)[0]+1)[0],
             color=storm_colors[1], alpha=0.33)


# normalized
ax_n = subfigs[-1].add_subplot(gs[2])
ax_n.axvline(0, lw=0.55, color='gray', ls=':')
ax_n.axhline(0, lw=0.55, color='gray', ls=':')
ax_n.set_xlim(xvals[0],xvals[-1])    
ax_n.set_xticks(xvals)
ax_n.tick_params(axis='both', which='major')  
ax_n.set_xticklabels(xlabels_f, minor=False, rotation=0, fontsize=fontsize)
ax_n.set_xlabel('Days Since Storm Start', fontsize=fontsize-1)

rel = tseries-clim
normalized = (rel-rel[0])/(np.max(rel)-np.min(rel))
ax_n.plot(xvals, normalized, color='k')

ax_n.axvline(np.where(xlabs==storm_range1[0].day)[0], color=storm_colors[0], alpha=0.75)
ax_n.axvline(np.where(xlabs==storm_range1[-1].day)[0]+1, color=storm_colors[0], alpha=0.75)
ax_n.axvline(np.where(xlabs==storm_range2[0].day)[0], color=storm_colors[1], alpha=0.75)
ax_n.axvline(np.where(xlabs==storm_range2[-1].day)[0]+1, color=storm_colors[1], alpha=0.75)


#%% sample case-pair plotting
plevel = 970 

start_dt = datetime(2017, 10, 1)
case1_dt = datetime(2017, 10, 4) #4->8
case2_dt = datetime(2017, 10, 6) #6->8
end_dt = datetime(2017, 10, 8)

era_data_file = '/Users/mundi/Desktop/era_data/antarctic_2017_09_10_msl_wind.nc'

with xr.open_dataset(era_data_file) as ds:
    slp = ds['msl'].sel(valid_time=slice(start_dt, end_dt))/100 #hpa
    slp_daily = slp.resample(valid_time='1D').mean(dim='valid_time')
    lon, lat = np.meshgrid(ds['longitude'].values, ds['latitude'].values)
    

#### sample pressure plots

for p_daily in slp_daily:
    datestr = str(p_daily.valid_time.values).split('T')[0]

    fig = plt.figure(figsize=(7.5,7.5)) 
    gs = GridSpec(nrows=1, ncols=1)
    ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
    ax = fx.setup_plot_sh(ax, labels=False)
    ax.set_title(datestr)
   
    ax.pcolormesh(lon, lat, p_daily, transform=ccrs.PlateCarree(),
                  cmap=cmo.haline, vmin=950, vmax=1000)
    ax.contour(lon, lat, p_daily, levels=[plevel], colors='magenta',
              transform=ccrs.PlateCarree())
    ax.contour(lon, lat, p_daily, levels=[957], colors='lime',
              transform=ccrs.PlateCarree())
    
    dt = datestr.split('-')
    si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt[0], dt[1], dt[2])
    fx.plot_geocontour(ax, si_lon, si_lat, si, [0.15,0.80], color='k', lw=1)

#%%% storm1 contours
import matplotlib.path as mpltPath

print(); print('Analyzing storm 1!')

case1_dt = datetime(2017, 10, 4)
storm_range1 = fx.daterange(case1_dt, end_dt, dt=24)

with xr.open_dataset(era_data_file) as ds:
    slp = ds['msl'].sel(valid_time=slice(case1_dt, end_dt))/100 #hpa
    slp_lon = slp.sel(longitude=slice(-180,0))
    slp_c1 = slp_lon.resample(valid_time='1D').mean(dim='valid_time')
    lon1, lat1 = np.meshgrid(slp_c1['longitude'].values, slp_c1['latitude'].values)

all_coords1 = []
for p_daily in slp_c1:
    fig, ax = fx.background_plot_sh(returnfig=True)
    contours = ax.contour(lon1, lat1, p_daily, levels=[plevel],
                          colors='k', transform=ccrs.PlateCarree())
    plt.close(fig)
    
    # min pressure
    indp = np.where(p_daily == np.nanmin(p_daily))
    minlon = np.squeeze(lon1[indp])
    minlat = np.squeeze(lat1[indp])
    points = np.array([ [minlon, minlat] ])
    
    cont = contours.allsegs[0]
    for cc in cont:
       # if len(cc) > 2 : path = mpltPath.Path(cc)
       # else: continue
       # pts_inside = path.contains_points(points) # radius=30
       
       # if any(pts_inside):
    
           all_coords1.append( cc )
    
#### plot storm areas
fig = plt.figure(figsize=(7.5,10)) 
gs = GridSpec(nrows=2, ncols=1, height_ratios=[1,0.5])
ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
ax = fx.setup_plot_sh(ax, labels=False)
ax.set_title(case2_dt.strftime('%d %b %Y'))

for cont in all_coords1:
    ax.plot(cont[:,0], cont[:,1], transform = ccrs.PlateCarree())
with fx.HidePrint(): bbox1 = fx.get_bbox_edges(all_coords1)
inside1 = fx.find_points_in_contour(bbox1, si_lon, si_lat)
geoplot_bbox(ax, bbox1, color='magenta')


#### calculate storm timeseries and plot
miz_points1 = get_miz_area(storm_range1)
tseries1 = []
clim1 = []
for dt1 in fx.daterange(case1_dt-timedelta(days=7), case1_dt+timedelta(days=14), dt=24):
    si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt1.year, dt1.month, dt1.day)
    si_bbox = np.ma.masked_array(si, mask=inside1)
    si_miz = np.ma.masked_where(miz_points1==0, si_bbox).filled(np.nan)
    tseries1.append(np.nanmean(si_miz*25*25))
    daily1 = []
    for yy in clim_years:
        si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, yy, dt2.month, dt2.day)
        si_bbox = np.ma.masked_array(si, mask=inside1)
        si_miz = np.ma.masked_where(miz_points1==0, si_bbox).filled(np.nan)
        daily1.append(np.nanmean(si_miz*25*25))
    clim1.append(np.nanmean(daily1))
        
tseries1 = np.array(tseries1); clim1=np.array(clim1)

# plot
ax_t1 = fig.add_subplot(gs[1])
# f1, ax_t1 = plt.subplots(1,1, figsize=(10,5))
ax_t1.axvline(0, lw=0.55, color='gray', ls=':')
ax_t1.axhline(0, lw=0.55, color='gray', ls=':')
ax_t1.plot(xxx, tseries1-clim1)

ax1 = ax_t1.twinx()
ax1.axhline(0, lw=0.55, color='green', ls=':')
t = tseries1-clim1
ax1.plot(xxx, (t-t[0])/(np.max(t)-np.min(t)), color='green')

#%%% storm2 contours
print(); print('Analyzing storm 2!')
case2_dt = datetime(2017, 10, 4)

storm_range2 = fx.daterange(case2_dt, end_dt, dt=24)

with xr.open_dataset(era_data_file) as ds:
    slp = ds['msl'].sel(valid_time=slice(case2_dt, end_dt))/100 #hpa
    slp_lon = slp.sel(longitude=slice(0,160))
    slp_c2 = slp_lon.resample(valid_time='1D').mean(dim='valid_time')
    lon2, lat2 = np.meshgrid(slp_c2['longitude'].values, slp_c2['latitude'].values)

all_coords2 = []
for p_daily in slp_c2:
    fig, ax = fx.background_plot_sh(returnfig=True)
    contours = ax.contour(lon2, lat2, p_daily, levels=[plevel],
                          colors='k', transform=ccrs.PlateCarree())
    plt.close(fig)
    
    all_coords2 += contours.allsegs[0]
    
#### plot storm areas
fig = plt.figure(figsize=(7.5,10)) 
gs = GridSpec(nrows=2, ncols=1, height_ratios=[1,0.5])
ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
ax = fx.setup_plot_sh(ax, labels=False)
ax.set_title(case2_dt.strftime('%d %b %Y'))

for cont in all_coords2:
    ax.plot(cont[:,0], cont[:,1], transform = ccrs.PlateCarree())
with fx.HidePrint(): bbox2 = fx.get_bbox_edges(all_coords2)
inside2 = fx.find_points_in_contour(bbox2, si_lon, si_lat)
geoplot_bbox(ax, bbox2, color='magenta')


#### calculate storm timeseries and plot
miz_points2 = get_miz_area(storm_range2)
tseries2 = []
clim2 = []
for dt2 in fx.daterange(case2_dt-timedelta(days=7), case2_dt+timedelta(days=14), dt=24):
    si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt2.year, dt2.month, dt2.day)
    si_bbox = np.ma.masked_array(si, mask=inside2)
    si_miz = np.ma.masked_where(miz_points2==0, si_bbox).filled(np.nan)
    tseries2.append(np.nanmean(si_miz*25*25))
    daily2 = []
    for yy in clim_years:
        si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, yy, dt2.month, dt2.day)
        si_bbox = np.ma.masked_array(si, mask=inside2)
        si_miz = np.ma.masked_where(miz_points2==0, si_bbox).filled(np.nan)
        daily2.append(np.nanmean(si_miz*25*25))
    clim2.append(np.nanmean(daily2))
        
tseries2 = np.array(tseries2); clim2=np.array(clim2)

# plot
ax_t2 = fig.add_subplot(gs[1])
ax_t2.axvline(0, lw=0.55, color='gray', ls=':')
ax_t2.axhline(0, lw=0.55, color='gray', ls=':')
ax_t2.plot(xxx, tseries2-clim2)


#%% second case-pair option! plotting
plevel = 970 

start_dt = datetime(2012, 10, 24)
case1_dt = datetime(2012, 10, 25) #25->29
case2_dt = datetime(2012, 10, 28) #28->31
end_dt = datetime(2012, 11, 1)

era_data_file = '/Users/mundi/Desktop/era_data/antarctic_2012_10_11_msl_wind.nc'

with xr.open_dataset(era_data_file) as ds:
    slp = ds['msl'].sel(valid_time=slice(start_dt, end_dt))/100 #hpa
    slp_daily = slp.resample(valid_time='1D').mean(dim='valid_time')
    lon, lat = np.meshgrid(ds['longitude'].values, ds['latitude'].values)
    

#### sample pressure plots

for p_daily in slp_daily:
    datestr = str(p_daily.valid_time.values).split('T')[0]

    fig = plt.figure(figsize=(7.5,7.5)) 
    gs = GridSpec(nrows=1, ncols=1)
    ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
    ax = fx.setup_plot_sh(ax, labels=False)
    ax.set_title(datestr)
   
    ax.pcolormesh(lon, lat, p_daily, transform=ccrs.PlateCarree(),
                  cmap=cmo.haline, vmin=950, vmax=1000)
    ax.contour(lon, lat, p_daily, levels=[plevel], colors='magenta',
              transform=ccrs.PlateCarree())
    ax.contour(lon, lat, p_daily, levels=[957], colors='lime',
              transform=ccrs.PlateCarree())
    
    dt = datestr.split('-')
    si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt[0], dt[1], dt[2])
    fx.plot_geocontour(ax, si_lon, si_lat, si, [0.15,0.80], color='k', lw=1)

#%%% storm1 contours
import matplotlib.path as mpltPath

print(); print('Analyzing storm 1!')

case1_dt = datetime(2012, 10, 25)
storm_range1 = fx.daterange(case1_dt, datetime(2012, 10, 29), dt=24)

with xr.open_dataset(era_data_file) as ds:
    slp = ds['msl'].sel(valid_time=slice(case1_dt, storm_range1[-1]))/100 #hpa
    slp_lon = slp.sel(longitude=slice(-120,0))
    slp_c1 = slp_lon.resample(valid_time='1D').mean(dim='valid_time')
    lon1, lat1 = np.meshgrid(slp_c1['longitude'].values, slp_c1['latitude'].values)

all_coords1 = []
for p_daily in slp_c1:
    fig, ax = fx.background_plot_sh(returnfig=True)
    contours = ax.contour(lon1, lat1, p_daily, levels=[plevel],
                          colors='k', transform=ccrs.PlateCarree())
    plt.close(fig)
    
    # min pressure
    indp = np.where(p_daily == np.nanmin(p_daily))
    minlon = np.squeeze(lon1[indp])
    minlat = np.squeeze(lat1[indp])
    points = np.array([ [minlon, minlat] ])
    
    cont = contours.allsegs[0]
    for cc in cont:
       if len(cc) > 2 : path = mpltPath.Path(cc)
       else: continue
       pts_inside = path.contains_points(points) # radius=30
       
       if any(pts_inside):
    
           all_coords1.append( cc )
    
#### plot storm areas
fig = plt.figure(figsize=(7.5,10)) 
gs = GridSpec(nrows=2, ncols=1, height_ratios=[1,0.5])
ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
ax = fx.setup_plot_sh(ax, labels=False)
ax.set_title(case1_dt.strftime('%d %b %Y'))

for cont in all_coords1:
    ax.plot(cont[:,0], cont[:,1], transform = ccrs.PlateCarree())

# with fx.HidePrint(): bbox1 = fx.get_bbox_edges(all_coords1)

minlat = np.nanmin([np.nanmin(cc[:,1]) for cc in all_coords1])
maxlat = np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords1])
minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords1])
maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords1])
bbox1 = fx.get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False)

inside1 = ~fx.find_points_in_contour(bbox1, si_lon, si_lat) ### inside bbox!!(~)
geoplot_bbox(ax, bbox1, color='magenta')

ax.pcolormesh(si_lon, si_lat, np.ma.masked_array(si, mask=inside1), transform=ccrs.PlateCarree() )

#### calculate storm timeseries and plot
miz_points1 = get_miz_area(storm_range1)
tseries1 = []
clim1 = []
for dt1 in fx.daterange(case1_dt-timedelta(days=7), case1_dt+timedelta(days=14), dt=24):
    si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt1.year, dt1.month, dt1.day)
    si_bbox = np.ma.masked_array(si, mask=inside1)
    si_miz = np.ma.masked_where(miz_points1==0, si_bbox).filled(np.nan)
    tseries1.append(np.nanmean(si_miz*25*25))
    daily1 = []
    for yy in clim_years:
        si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, yy, dt1.month, dt1.day)
        si_bbox = np.ma.masked_array(si, mask=inside1)
        si_miz = np.ma.masked_where(miz_points1==0, si_bbox).filled(np.nan)
        daily1.append(np.nanmean(si_miz*25*25))
    clim1.append(np.nanmean(daily1))
        
tseries1 = np.array(tseries1); clim1=np.array(clim1)

# plot
ax_t1 = fig.add_subplot(gs[1])
# f1, ax_t1 = plt.subplots(1,1, figsize=(10,5))
ax_t1.axvline(0, lw=0.55, color='gray', ls=':')
ax_t1.axhline(0, lw=0.55, color='gray', ls=':')
ax_t1.plot(xxx, tseries1-clim1)

ax1 = ax_t1.twinx()
ax1.axhline(0, lw=0.55, color='green', ls=':')
t = tseries1-clim1
ax1.plot(xxx, (t-t[0])/(np.max(t)-np.min(t)), color='green')

#%%% storm2 contours
print(); print('Analyzing storm 2!')

case2_dt = datetime(2012, 10, 29)
storm_range2 = fx.daterange(case2_dt, datetime(2012, 10, 31), dt=24)

with xr.open_dataset(era_data_file) as ds:
    # ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 
    # ds = ds.sortby(ds.longitude)
    # ds.longitude
    lon_name = 'longitude' 
    # Adjust lon values to make sure they are within (-180, 180)
    ds['_longitude_adjusted'] = xr.where(
        ds[lon_name] < 0,
        ds[lon_name] + 360,
        ds[lon_name])
    # reassign the new coords to as the main lon coords
    # and sort DataArray using new coordinate values
    ds = (
        ds
        .swap_dims({lon_name: '_longitude_adjusted'})
        .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
        .drop_vars(lon_name))
    ds = ds.rename({'_longitude_adjusted': lon_name})
    
    slp = ds['msl'].sel(valid_time=slice(case2_dt, storm_range2[-1]))/100 #hpa
    slp_lon = slp.sel(longitude=slice(60,300))
    slp_c2 = slp_lon.resample(valid_time='1D').mean(dim='valid_time')
    lon2, lat2 = np.meshgrid(slp_c2['longitude'].values, slp_c2['latitude'].values)

all_coords2 = []
for p_daily in slp_c2:
    fig, ax = fx.background_plot_sh(returnfig=True)
    contours = ax.contour(lon2, lat2, p_daily, levels=[plevel],
                          colors='k', transform=ccrs.PlateCarree())
    plt.close(fig)
    
    # min pressure
    indp = np.where(p_daily == np.nanmin(p_daily))
    minlon = np.squeeze(lon2[indp])
    minlat = np.squeeze(lat2[indp])
    points = np.array([ [minlon, minlat] ])
    
    cont = contours.allsegs[0]
    for cc in cont:
       if len(cc) > 2 : path = mpltPath.Path(cc)
       else: continue
       pts_inside = path.contains_points(points) # radius=30
       
       if any(pts_inside):
    
           all_coords2.append( cc )
    
    # all_coords2 += contours.allsegs[0]
    
#### plot storm areas
fig = plt.figure(figsize=(7.5,10)) 
gs = GridSpec(nrows=2, ncols=1, height_ratios=[1,0.5])
ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
ax = fx.setup_plot_sh(ax, labels=False)
ax.set_title(case2_dt.strftime('%d %b %Y'))

for cont in all_coords2:
    ax.plot(cont[:,0], cont[:,1], transform = ccrs.PlateCarree())
with fx.HidePrint(): bbox2 = fx.get_bbox_edges(all_coords2)
inside2 = fx.find_points_in_contour(bbox2, si_lon, si_lat)
geoplot_bbox(ax, bbox2, color='magenta')


#### calculate storm timeseries and plot
miz_points2 = get_miz_area(storm_range2)
tseries2 = []
clim2 = []
for dt2 in fx.daterange(case2_dt-timedelta(days=7), case2_dt+timedelta(days=14), dt=24):
    si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt2.year, dt2.month, dt2.day)
    si_bbox = np.ma.masked_array(si, mask=inside2)
    si_miz = np.ma.masked_where(miz_points2==0, si_bbox).filled(np.nan)
    tseries2.append(np.nanmean(si_miz*25*25))
    daily2 = []
    for yy in clim_years:
        si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, yy, dt2.month, dt2.day)
        si_bbox = np.ma.masked_array(si, mask=inside2)
        si_miz = np.ma.masked_where(miz_points2==0, si_bbox).filled(np.nan)
        daily2.append(np.nanmean(si_miz*25*25))
    clim2.append(np.nanmean(daily2))
        
tseries2 = np.array(tseries2); clim2=np.array(clim2)

# plot
ax_t2 = fig.add_subplot(gs[1])
# f2, ax_t2 = plt.subplots(1,1, figsize=(10,5))
ax_t2.axvline(0, lw=0.55, color='gray', ls=':')
ax_t2.axhline(0, lw=0.55, color='gray', ls=':')
ax_t2.plot(xxx, tseries2-clim2)

ax2 = ax_t2.twinx()
ax2.axhline(0, lw=0.55, color='green', ls=':')
t = tseries2-clim2
ax2.plot(xxx, (t-t[0])/(np.max(t)-np.min(t)), color='green')

#%%--
#%%% functions
from scipy.interpolate import griddata

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


def series(var, storm_event, miz, in_bbox):
    dt1 = storm_event[0]-timedelta(days=7)
    dt2 = storm_event[0]+timedelta(days=14)
    analysis_range = fx.daterange(dt1, dt2, dt=24)
    
    out = []
    for dt in analysis_range:
        vx = var.sel(valid_time=dt).values
        vx_in = np.ma.masked_array(vx, mask=in_bbox).filled(np.nan)
        vx_miz = np.ma.masked_where(miz<1, vx_in).filled(np.nan)
        out.append(np.nanmean(vx_miz))
        
    return np.array(out)

#%% Worldview plotting
from osgeo import gdal
import cartopy.crs as ccrs
import cartopy.feature as cfeature

### set up plot
fig=plt.figure(figsize=[15,15]) 
ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.set_extent([-180,180,-53,-90], ccrs.PlateCarree())
ax.gridlines(draw_labels=False)
ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=100)
ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=101)
ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=102)
 

### load and plot worldview image
img_file = '/Users/mundi/Downloads/snapshot-2010-04-15T00_00_00Z.tif'
ds = gdal.Open(img_file)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
proj = ds.GetProjection()   

extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])

img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)), extent=extent,
                origin='upper')
    

### add sea ice contours
si1, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010,4,15)
fx.plot_geocontour(ax, si_lon, si_lat, si1, [0.15], color='b', lw=1)
si2, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010,4,18)
fx.plot_geocontour(ax, si_lon, si_lat, si2, [0.15], color='r', lw=1)

### add storm contours
ncname = datetime(2010,4,14).strftime('%Y_%m%d')+'_contours.nc'
with xr.open_dataset(path1+'contours/'+ncname) as cs:
    all_coords = []
    for key in list(cs.keys()):
        coord = cs[key].values
        all_coords.append(coord)
        ax.plot(coord[:,0], coord[:,1], color='magenta', lw=2,
                transform=ccrs.PlateCarree())
    
# plot bbox
bbox = fx.get_bbox_edges(all_coords)
bboxs.append(bbox)
inside = fx.find_points_in_contour(bbox, si_lon, si_lat) ### inside bbox!!(~)
geoplot_bbox(ax, bbox, color='navy')


#%%% another possible case

### set up plot
fig=plt.figure(figsize=[15,15]) 
ax = plt.axes(projection=ccrs.SouthPolarStereo())
# ax.set_extent([-180,180,-53,-90], ccrs.PlateCarree())
ax.set_extent([-60,60,-59,-90], ccrs.PlateCarree())
ax.gridlines(draw_labels=False)
ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=100)
ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=101)
ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=102)
 

### load and plot worldview image
img_file = '/Users/mundi/Downloads/snapshot-2011-04-16T00_00_00Z.tif'
ds = gdal.Open(img_file)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
proj = ds.GetProjection()   

extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])

img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)), extent=extent,
                origin='upper')
    

### add sea ice contours
si1, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2011,4,15)
fx.plot_geocontour(ax, si_lon, si_lat, si1, [0.15], color='b', lw=1)
si2, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2011,4,18)
fx.plot_geocontour(ax, si_lon, si_lat, si2, [0.15], color='r', lw=1)

### add storm contours
colors = iter(['#fbb4b9','#f768a1','#c51b8a','#7a0177'])
# colors = iter(['#d7b5d8','#df65b0','#dd1c77','#980043'])
ncname = datetime(2011,4,15).strftime('%Y_%m%d')+'_contours.nc'
with xr.open_dataset(path1+'contours/'+ncname) as cs:
    all_coords = []
    for key in list(cs.keys()):
        coord = cs[key].values
        all_coords.append(coord)
        if key.split('_')[-2] == '980':
            ax.plot(coord[:,0], coord[:,1], color='k', lw=3,
                    transform=ccrs.PlateCarree())
            ax.plot(coord[:,0], coord[:,1], color=next(colors), lw=2,
                    transform=ccrs.PlateCarree())
    
# plot bbox
# bbox = fx.get_bbox_edges(all_coords)
minlat = np.nanmin([np.nanmin(cc[:,1]) for cc in all_coords])
maxlat = np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords])
minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
bbox1 = fx.get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False)

inside = fx.find_points_in_contour(bbox1, si_lon, si_lat) ### inside bbox!!(~)
geoplot_bbox(ax, bbox1, color='navy', lw=4)

# sic
cmap = cmo.balance_r
cmap.set_bad(alpha=0)
diff= np.where(si2-si1 == 0, np.nan, si2-si1)
ax.pcolormesh(si_lon, si_lat, np.ma.masked_array(diff, mask=~inside).filled(np.nan), 
              cmap=cmap, vmin=-0.3, vmax=0.3, alpha=0.5,
              transform=ccrs.PlateCarree())

# winds?

#%%% subplots (!)
####!!!!
nlines=2
fig = plt.figure(figsize=(10,6.5))
gs = GridSpec(1+nlines, 2, height_ratios=[1]+[0.25]*nlines)

storm_event = fx.daterange(datetime(2011,4,15), datetime(2011,4,18), dt=24)
f1 = '/Users/mundi/Desktop/era_data/antarctic_2011_04_msl_wind_t2m.nc'
f2 = '/Users/mundi/Desktop/era_data/antarctic_2011_04_waves.nc'

fig.suptitle('\n\n\n'+storm_event[0].strftime('%B %d, %Y'))

####################
#### worldview map
####################
ax = fig.add_subplot(gs[0,0], projection=ccrs.SouthPolarStereo())
ax.set_extent([-60,60,-59,-90], ccrs.PlateCarree())
ax.gridlines(draw_labels=False)
ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=1)
ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=2)
ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=3)
ax.set_title('(a) Storm Area')

# plot image
img_file = '/Users/mundi/Downloads/snapshot-2011-04-16T00_00_00Z.tif'
ds = gdal.Open(img_file)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
proj = ds.GetProjection()   
extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])
img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)), extent=extent, origin='upper', zorder=-1)

# add sea ice contours
si1, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2011,4,15)
fx.plot_geocontour(ax, si_lon, si_lat, si1, [0.15], color='k', lw=1)
si2, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2011,4,18)
fx.plot_geocontour(ax, si_lon, si_lat, si2, [0.15], color='dimgray', lw=1)

# add storm contours
colors = iter(['#fbb4b9','#f768a1','#c51b8a','#7a0177'])
ncname = datetime(2011,4,15).strftime('%Y_%m%d')+'_contours.nc'
with xr.open_dataset(path1+'contours/'+ncname) as cs:
    all_coords = []
    for key in list(cs.keys()):
        coord = cs[key].values
        all_coords.append(coord)
        if key.split('_')[-2] == '980':
            ax.plot(coord[:,0], coord[:,1], color='k', lw=2,
                    transform=ccrs.PlateCarree())
            ax.plot(coord[:,0], coord[:,1], color=next(colors), lw=1.5,
                    transform=ccrs.PlateCarree())
    
# plot bbox
minlat = np.nanmin([np.nanmin(cc[:,1]) for cc in all_coords])
maxlat = np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords])
minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
bbox1 = fx.get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False)
inside = fx.find_points_in_contour(bbox1, si_lon, si_lat) 
geoplot_bbox(ax, bbox1, color='navy', lw=2)


##################
#### sea ice map
##################
ax = fig.add_subplot(gs[0,1], projection=ccrs.SouthPolarStereo())
ax.set_extent([-60,60,-59,-90], ccrs.PlateCarree())
ax.gridlines(draw_labels=False)
ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=1)
ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=2)
ax.set_title('(b) Change in Sea Ice Concentration')

cmap = cmo.balance_r
cmap.set_bad(alpha=0)
diff= np.where(si2-si1 == 0, np.nan, si2-si1)
psm = ax.pcolormesh(si_lon, si_lat, 
                    np.ma.masked_array(diff, mask=~inside).filled(np.nan), 
                    cmap=cmap, vmin=-0.3, vmax=0.3, alpha=1,
                    transform=ccrs.PlateCarree())
geoplot_bbox(ax, bbox1, color='navy', lw=2)

miz_points =  get_miz_area(storm_event)
miz_points = np.ma.masked_array(miz_points, mask=~inside).filled(np.nan)
ax = fx.plot_geocontour(ax, si_lon, si_lat, miz_points, levels=[0.99], color='k', lw=1, ls='solid')

cax1 = fig.add_axes([0.55,0.475,0.35,0.015]) 
cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
# cbar1.set_label(r'Change in Sea Ice Concentration', fontsize=10)
# cax1.tick_params(labelsize=10)


# plot winds
with xr.open_dataset(f1) as ds:
    ds = ds.resample(valid_time='1D').mean() # daily mean winds!
    ds = ds.sel(latitude=slice(-60,-80))
    u = ds['u10'].sel(valid_time=datetime(2011,4,16))
    v = ds['v10'].sel(valid_time=datetime(2011,4,16))
    x,y = np.meshgrid(ds['longitude'], ds['latitude'])
    
miz_grd = griddata((si_lon.flatten(), si_lat.flatten()), miz_points.flatten(),
                            (x.flatten(), y.flatten()), method='nearest')
miz_grd = miz_grd.reshape(np.shape(x))
    
inside2 = fx.find_points_in_contour(bbox1, x,y) 
ax, Q = fx.plot_winds(x,y,u,v, ax, mask=~inside2, sc=300, interval=12)

qk = ax.quiverkey(Q, 0.875, 0.95, 15, '\n'+r'$15 \frac{m}{s}$', 
                  labelpos='E', angle=90,labelsep=0.075, 
                   coordinates='axes', transform=ccrs.PlateCarree())



#############
#### lines
#############
for idx in range(nlines):
    ax = fig.add_subplot(gs[idx+1,1])
    ax.set_xticks(xxx)
    ax.set_xticklabels(xlabels, minor=False, rotation=0)
    ax.set_xlim(-7,14)
    ax.axvline(0, color='k', lw=0.55, ls=':')
    
    if idx==0: # sea ice
        ax.set_title('(d) Time Series')
        si, clim = seaice_series(storm_event, miz_points, ~inside)
        ax.plot(xxx, si/1e5, color='navy', ls='-', lw=1.25)
        ax.plot(xxx, clim/1e5, color='navy', ls='--', lw=1.25)
    elif idx==1: # winds, sst
        v_winds = series(ds['v10'], storm_event, miz_grd, ~inside2)
        ax.plot(xxx, v_winds, color='navy', ls='-', lw=1.25)
        ax.axhline(0, color='k', lw=0.55, ls=':')
        
        #sst?
        
        
ax.set_xlabel('Days Since Storm Start')


################
#### spaghetti
################
ax1 = fig.add_subplot(gs[1:,0])
path1= root_paths[1]
month_group = [4,5]
mnames = [calendar.month_abbr[mm]+', ' for mm in month_group]

ax1.set_title('(c) Change in MIZ Area For All Storms: '+str(month_group))
ax1.set_ylabel('Normalized Relative\nChange in Ice Area')
ax1.set_xlabel('Days Since Storm Start')
ax1.set_xlim(-7,14)
ax1.set_ylim(-1.05,1.05)
ax1.axhline(0, ls='-', color='k', lw=1)
ax1.axvline(0, ls='-', color='k', lw=0.75)
ax1.set_xticks(xxx)
ax1.set_xticklabels(xlabels, minor=False, rotation=0)
ax1.tick_params(axis='both', which='major')  

mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
    fx.indiv_lines(np.arange(2010,2020), path1+'census/', path1+'area/', path1+'seaice/')
    
all_lines = []
for month in month_group: all_lines+=lines[month]

ax1.plot(xxx, np.array(all_lines).T, lw=0.55, color='gray')
ax1.plot(xxx, np.nanmean(all_lines, axis=0), lw=2, color='maroon', label='2010-2019 Mean')

ss = si-clim
pcd = (ss-ss[0])/(np.nanmax(ss)-np.nanmin(ss))
ax1.plot(xxx, pcd, lw=2, color='navy', label='Case Study')

ax1.legend(loc='lower left', handletextpad=0.5, handlelength=1.25)

#%% end
