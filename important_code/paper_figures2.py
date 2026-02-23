#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 2025
paper_figures2.py
updated 26 jan 2026

1: case study, sample spaghetti
2: fan
3: taylor sample
4: global bar plots (global_seaice1)

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

#%%% bar fxn
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

#%% case study, sample spaghetti
'''sh_cases#.py (1)'''

print(); print('*** CASE STUDY PLOT ***')

nlines=2
fig = plt.figure(figsize=(12,8))
gs = GridSpec(1+nlines, 2, height_ratios=[1]+[0.25]*nlines)
# plt.subplots_adjust(hspace=0.3)

storm_event = fx.daterange(datetime(2011,4,15), datetime(2011,4,18), dt=24)
f1 = '/Users/mundi/Desktop/era_data/antarctic_2011_04_msl_wind_t2m.nc'
f2 = '/Users/mundi/Desktop/era_data/antarctic_2011_04_waves.nc'

# fig.suptitle('\n\n\n'+storm_event[0].strftime('%B %d, %Y'))

#########################
#### worldview map
print('* worldview map')
#########################
ax = fig.add_subplot(gs[0,0], projection=ccrs.SouthPolarStereo())
ax.set_extent([-60,60,-59,-90], ccrs.PlateCarree())
ax.gridlines(draw_labels=False)
ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=1)
ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=2)
ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=3)
ax.set_title('(a) Storm Area')

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
cc = ['#67a9cf', '#2166ac']
si1, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[0].year, storm_event[0].month, storm_event[0].day)
fx.plot_geocontour(ax, si_lon, si_lat, si1, [0.15], color=cc[0], lw=1)
si2, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[-1].year, storm_event[-1].month, storm_event[-1].day)
fx.plot_geocontour(ax, si_lon, si_lat, si2, [0.15], color='k', lw=1)

# add storm contours
colors = iter(['#fbb4b9','#f768a1','#c51b8a','#7a0177'])
ncname = datetime(2011,4,15).strftime('%Y_%m%d')+'_contours.nc'
with xr.open_dataset(path1+'contours/'+ncname) as cs:
    all_coords = []
    for key in list(cs.keys()):
        coord = cs[key].values
        date = datetime.strptime(key.split('_')[1]+key.split('_')[2], '%Y%m%d')
        if key.split('_')[-2] == '980':
            all_coords.append(coord)
            ax.plot(coord[:,0], coord[:,1], color='k', lw=2,
                    transform=ccrs.PlateCarree())
            ax.plot(coord[:,0], coord[:,1], color=next(colors), lw=1.5,
                    transform=ccrs.PlateCarree(),
                    label = date.strftime('%-d'))
ax.legend(loc='lower right', ncol = len(all_coords), handletextpad=0.5, handlelength=0.85, 
          title=date.strftime('%B %Y'), labelspacing=0.15, columnspacing=0.2)
    
# plot bbox
minlat = np.nanmin([np.nanmin(cc[:,1]) for cc in all_coords])
maxlat = np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords])
minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
bbox1 = fx.get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False)
inside = fx.find_points_in_contour(bbox1, si_lon, si_lat) 
geoplot_bbox(ax, bbox1, color='navy', lw=2)


miz_points = get_miz_area(storm_event)
arr1 = np.ma.masked_array(np.ones(np.shape(si_lon)), mask=~inside).filled(np.nan)
arr1 = np.where(miz_points<1, np.nan, arr1)

ax.pcolormesh(si_lon, si_lat, arr1, transform=ccrs.PlateCarree(),
              cmap=cmo.balance_r, vmin=-1,vmax=1, alpha=0.5)


#######################
#### sea ice map
print('* sea ice map')
#######################
ax = fig.add_subplot(gs[0,1], projection=ccrs.SouthPolarStereo())
ax.set_extent([-60,60,-59,-90], ccrs.PlateCarree())
ax.gridlines(draw_labels=False)
ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=1)
ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=2)
ax.set_title(r'(b) $\Delta$SIC: '+storm_event[-1].strftime('%b %d')+'-'+storm_event[0].strftime('%d'))

cmap = cmo.balance_r
cmap.set_bad(alpha=0)
diff= np.where(si2-si1 == 0, np.nan, si2-si1)
pcm = ax.pcolormesh(si_lon, si_lat, 
                    np.ma.masked_array(diff, mask=~inside).filled(np.nan), 
                    cmap=cmap, vmin=-0.3, vmax=0.3, alpha=1,
                    transform=ccrs.PlateCarree())
geoplot_bbox(ax, bbox1, color='navy', lw=2)

miz_points = np.ma.masked_array(miz_points, mask=~inside).filled(np.nan)
ax = fx.plot_geocontour(ax, si_lon, si_lat, miz_points, levels=[0.99], color='k', lw=1, ls='solid')

cax1 = fig.add_axes([0.55,0.475,0.35,0.0125]) 
cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
cbar1.set_label(r'Change in Sea Ice Concentration')#, fontsize=10)
# cax1.tick_params(labelsize=10)


# plot winds
wind_dt = datetime(2011,4,16)
with xr.open_dataset(f1) as ds:
    ds = ds.resample(valid_time='1D').mean() # daily mean winds!
    ds = ds.sel(latitude=slice(-60,-80))
    u = ds['u10'].sel(valid_time=wind_dt)
    v = ds['v10'].sel(valid_time=wind_dt)
    x,y = np.meshgrid(ds['longitude'], ds['latitude'])
    
miz_grd = griddata((si_lon.flatten(), si_lat.flatten()), miz_points.flatten(),
                            (x.flatten(), y.flatten()), method='nearest')
miz_grd = miz_grd.reshape(np.shape(x))
    
inside2 = fx.find_points_in_contour(bbox1, x,y) 
ax, Q = fx.plot_winds(x,y,u,v, ax, mask=~inside2, sc=300, interval=12)

qk = ax.quiverkey(Q, 0.875, 0.95, 15, '\n'+r'$15 \frac{m}{s}$'+'\n'+wind_dt.strftime('%-m-%-d'), 
                  labelpos='E', angle=90,labelsep=0.075, 
                   coordinates='axes', transform=ccrs.PlateCarree())

# add loc of min pressure
print('finding minimum pressure')
min_p = 9999
for dt1 in storm_event:
    msl = ds['msl'].sel(valid_time=dt1).values/100
    msl_in = np.ma.masked_array(msl, mask=~inside2)
    if np.nanmin(msl_in)<min_p:
        min_p = np.nanmin(msl_in)
        xmin = x[np.where(msl_in==min_p)]
        ymin = y[np.where(msl_in==min_p)]
        print('~', min_p, dt1, xmin, ymin)
        
ax.plot(xmin, ymin, marker='*', mfc='yellow', mec='k', ms=14,
        transform=ccrs.PlateCarree(), zorder=100,
        label='Min. Pressure: '+str(round(min_p, 1))+' hPa')
ax.legend(loc='lower right', handletextpad=0.5, handlelength=1, labelspacing=0.15)


#######################
#### timeseries
print('* time series')
#######################

for idx in range(nlines):
    ax = fig.add_subplot(gs[idx+1,0])
    ax.set_xticks(xxx)
    ax.set_xlim(-7,14)
    ax.axvline(0, color='k', lw=0.55, ls=':')
    ax.axvline((storm_event[-1]-storm_event[0]).days + 1, color='k', lw=0.55, ls=':')
    
    if idx==0: # sea ice
        ax.set_title('(c) MIZ Ice Area')
        ax.set_xticklabels(['']*len(xxx))
        si, clim = seaice_series(storm_event, miz_points, ~inside)
        ax.plot(xxx, si/1e5, color='navy', ls='-', lw=1.5)
        ax.plot(xxx, clim/1e5, color='navy', ls='--', lw=1.25, label='2010-2019')
        ax.set_ylabel(r'$\times 10^5$ km$^2$')
        ax.legend(loc='lower right', handletextpad=0.5, handlelength=1.25)
    elif idx==1: # winds, sst
        ax.set_title('(d) Mean Meridional Winds')
        ax.set_xticklabels(xlabels, minor=False, rotation=0)
        v_winds = series(ds['v10'], storm_event, miz_grd, ~inside2)
        ax.plot(xxx, v_winds, color='navy', ls='-', lw=1.5)
        ax.axhline(0, color='k', lw=0.55, ls=':')
        ax.set_ylabel(r'm s$^{-1}$')
    # elif idx==2:
    #     #sst?
    #     mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
    #         fx.indiv_lines(np.arange(2010,2020), path1+'census/', path1+'area/', path1+'seaice/')
    #     SSTs = np.load(path1+'sst/'+'tseries_'+'0'+'-'+'4'+'.npy')
    #     for sd, sst in zip(start_day[4], SSTs):
    #         if sd==storm_event[0]:
    #             ax.plot(xxx, sst, color='dimgray', ls='-', lw=1.25)
    #             break
        
ax.set_xlabel('Days Since Storm Start')

print('timeseries calculations')

rel = (si/1e5) - (clim/1e5)
diff = rel[7+4] - rel[7]
print('- difference:', round(diff, 3))


####################
#### spaghetti
print('* spahetti')
####################
ax1 = fig.add_subplot(gs[1:,1])
path1= root_paths[1]
month_group = [4,5]
mnames = [calendar.month_abbr[mm]+', ' for mm in month_group]

ax1.set_title('(e) Change in MIZ Area For All Storms: '+str(month_group))
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

# ---- print out spaghetti stats ---- #
print('spaghetti stats')
print('- number of storms: '+str(len(all_lines)))
print('- number that decrease: '+str(len([x for x in all_lines if x[-1]<0])))
print('- number that increase: '+str(len([x for x in all_lines if x[-1]>0])))
print('- percent that increase: '+str( round(len([x for x in all_lines if x[-1]>0])*100/len(all_lines),3) ))


#%% fan plots
'''seaice_impact#.py'''

print(); print('*** FAN PLOTS ***')


alph = ['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k']
month_colors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
                '#980043','#7a0177', '#253494','#2c7fb8','#41b6c4','#a1dab4']

def ytext_scale(idx, li, yi):
    yadd=0
    idx += 1 # actual month
    
    if li==0: # arctic
        if yi == 0: # left (past)
            if idx==11: yadd=0.01 # november
            elif idx==2: yadd=-0.02 # feb
            elif idx==12: yadd=0.01 # dec
        elif yi ==1: # right (recent)
            if idx==4: yadd=0.05 # apr
            elif idx==3:yadd=0.01 # march
            elif idx==1: yadd=-0.05 # jan
            elif idx == 5: yadd=-0.05 # may
            elif idx == 8: yadd=-0.05 # aug
            elif idx==9: yadd=-0.05 # sep
            elif idx==11:yadd=0.015 # nov
    elif li==1: #antarctic
        if yi == 0: # left (past)
            if idx==6: yadd=0.02
            elif idx==7: yadd=-0.02
        elif yi ==1: # right (recent)
            if idx==8: yadd=0.02 #aug
            elif idx==9: yadd=-0.02 # sep
    
    return yadd


#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(16,14), sharey=True)
# plt.subplots_adjust(wspace=0.3)
# fig2.suptitle('\n\nMonthly Mean Change in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
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
for li, loc in enumerate(['arctic','aa']):
    for yi, years in enumerate(decades):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        ### monthly mean
        ranges = {7:[], 14:[]}
        axes2[li][yi].set_title(hemi_names[li]+': '+yr_title, fontsize=fs+2)
        for idx, ml in enumerate(mean_lines):
            axes2[li][yi].plot(xxx, ml, color = month_colors[idx], lw=2, label = calendar.month_name[idx+1])
            yadd = ytext_scale(idx, li, yi)
            mlabel = calendar.month_abbr[idx+1] #+' ('+str(len(lines[idx+1]))+')'
            axes2[li][yi].text(xxx[-1]+0.25, ml[-1]+yadd, mlabel, fontsize=fs-1) 
            ranges[7].append(ml[7+7])
            ranges[14].append(ml[7+14])
            
        range1 = [round(np.nanmin(ranges[7]),4), round(np.nanmax(ranges[7]),4)]    
        range2 = [round(np.nanmin(ranges[14]),4), round(np.nanmax(ranges[14]),4)]    
        print(loc, yr_title, range1, range2)
        

axes2[0][1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                  edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                  bbox_to_anchor=(1.5,0.33));

#%% taylor plot - see other script (!)
'''taylor#.py'''

print(); print('*** TAYLOR PLOTS ***')
print('--> see taylor5.py')


#%% global bar
'''seaice_sens#.py | global_seaice1.py'''

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
            
            print('* total spread')
            print(ht1*100/np.nansum(sums[1], axis=0))


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
    daynum_range = [datetime(2010,month,1).timetuple().tm_yday,
                    datetime(2010,month,calendar.monthrange(2010,month)[-1]).timetuple().tm_yday]
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





#%% end
