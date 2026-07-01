#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 2026

demo.py

recreates main paper figures with source data

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


decades = [np.arange(1982, 1992), np.arange(2010,2020)]

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]

hemi_names= ['Arctic', 'Antarctic']

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

fontsize=11
fs = 14

# ------------------------ #
## source data directory ###
# ------------------------ #
MAIN = '/Users/mundi/Desktop/month-hemi/code/source_data/'
# ------------------------ #
ice_fname = MAIN+'seaice/'

#%%% functions

def csv_open(file):
    import csv
    results = []
    with open(file) as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)
    return np.array(results)

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

#%% FIGURE 1
print(); print('*** CASE STUDY PLOT ***')
DIR1 = MAIN+'figure_1/'

_, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2011, 4, 15)

nlines=2
fig = plt.figure(figsize=(12,8))
gs = GridSpec(1+nlines, 2, height_ratios=[1]+[0.25]*nlines)

storm_event = fx.daterange(datetime(2011,4,15), datetime(2011,4,18), dt=24)
miz_points = get_miz_area(storm_event)


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
img_file = DIR1+'snapshot-2011-04-16T00_00_00Z.tif' # nasa worldview file
with gdal.Open(img_file) as ds:
    data = ds.ReadAsArray()
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()   
    extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
              gt[3] + ds.RasterYSize * gt[5], gt[3])
    img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)), extent=extent, origin='upper', zorder=-1)

# add sea ice contours
si1, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[0].year, storm_event[0].month, storm_event[0].day)
si2, si_lon, si_lat = fx.load_seaice_sh(ice_fname, storm_event[-1].year, storm_event[-1].month, storm_event[-1].day)

# add storm contours
colors = iter(['#fbb4b9','#f768a1','#c51b8a','#7a0177'])

for date in storm_event:
    cfile = DIR1 + 'coords_'+date.strftime('%Y_%m%d')+'_980_1.csv'
    coord = csv_open(cfile)
    
    ax.plot(coord[:,0], coord[:,1], color='k', lw=2,
            transform=ccrs.PlateCarree())
    ax.plot(coord[:,0], coord[:,1], color=next(colors), lw=1.5,
            transform=ccrs.PlateCarree(),
            label = date.strftime('%-d'))
            
ax.legend(loc='lower right', ncol = len(storm_event), handletextpad=0.5, handlelength=0.85, 
          title=date.strftime('%B %Y'), labelspacing=0.15, columnspacing=0.2)


# plot bbox
bbox1 = csv_open(DIR1 + 'bounding_box.csv')
inside = fx.find_points_in_contour(bbox1, si_lon, si_lat) 
geoplot_bbox(ax, bbox1, color='navy', lw=2)

# miz area
arr2 = np.ma.masked_array(miz_points, mask=~inside)
ax = fx.plot_geocontour(ax, si_lon, si_lat, arr2, levels=[0.99], color='k', lw=1.25, ls='solid')

arr1 = np.ma.masked_array(np.ones(np.shape(si_lon)), mask=~inside).filled(np.nan)
arr1 = np.where(miz_points<1, np.nan, arr1)
ax.pcolormesh(si_lon, si_lat, arr1, transform=ccrs.PlateCarree(),
              cmap=cmo.balance_r, vmin=-1,vmax=1, alpha=0.33)


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
cbar1.set_label(r'Change in Sea Ice Concentration')


# plot winds
u = csv_open(DIR1+'u_winds.csv')
v = csv_open(DIR1+'v_winds.csv')
x = csv_open(DIR1+'x_winds.csv')
y = csv_open(DIR1+'y_winds.csv')
    
miz_grd = griddata((si_lon.flatten(), si_lat.flatten()), miz_points.flatten(),
                            (x.flatten(), y.flatten()), method='nearest')
miz_grd = miz_grd.reshape(np.shape(x))
    
inside2 = fx.find_points_in_contour(bbox1, x,y) 
ax, Q = fx.plot_winds(x,y,u,v, ax, mask=~inside2, sc=300, interval=12)

qk = ax.quiverkey(Q, 0.875, 0.95, 15, '\n'+r'$15 \frac{m}{s}$'+'\n'+'4-16', 
                  labelpos='E', angle=90,labelsep=0.075, 
                   coordinates='axes', transform=ccrs.PlateCarree())

# add loc of min pressure
loc = csv_open(DIR1+ 'location_min_pressure.csv')
xmin = loc[0] 
ymin = loc[1]
min_p = float(loc[2][0])       

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
        si_csv = csv_open(DIR1 + "miz_ice_area_tseries.csv")
        si = si_csv[0]
        clim = si_csv[1]
    
        ax.set_title('(c) MIZ Ice Area')
        ax.set_xticklabels(['']*len(xxx))
        ax.plot(xxx, si/1e5, color='navy', ls='-', lw=1.5)
        ax.plot(xxx, clim/1e5, color='navy', ls='--', lw=1.25, label='2010-2019')
        ax.set_ylabel(r'$\times 10^5$ km$^2$')
        ax.legend(loc='lower right', handletextpad=0.5, handlelength=1.25)
        
    elif idx==1: # winds, sst
        v_winds = csv_open(DIR1+ "wind_tseries.csv")
        ax.set_title('(d) Mean Meridional Winds')
        ax.set_xticklabels(xlabels, minor=False, rotation=0)
        ax.plot(xxx, v_winds.T, color='navy', ls='-', lw=1.5)
        ax.axhline(0, color='k', lw=0.55, ls=':')
        ax.set_ylabel(r'm s$^{-1}$')
        
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
month_group = [4,5]
mnames = [calendar.month_name[mm] for mm in month_group]
mtitle = mnames[0]+', '+mnames[1]

ax1.set_title('(e) Change in MIZ Area For All Storms: '+mtitle)
ax1.set_ylabel('Normalized Relative\nChange in Ice Area')
ax1.set_xlabel('Days Since Storm Start')
ax1.set_xlim(-7,14)
ax1.set_ylim(-1.05,1.05)
ax1.axhline(0, ls='-', color='k', lw=1)
ax1.axvline(0, ls='-', color='k', lw=0.75)
ax1.set_xticks(xxx)
ax1.set_xticklabels(xlabels, minor=False, rotation=0)
ax1.tick_params(axis='both', which='major')  

all_lines = csv_open(DIR1 + "normalized_tseries.csv")

ax1.plot(xxx, np.array(all_lines).T, lw=0.55, color='gray')
ax1.plot(xxx, np.nanmean(all_lines, axis=0), lw=2, color='maroon', label='2010-2019 Mean')

ss = si-clim
pcd = (ss-ss[0])/(np.nanmax(ss)-np.nanmin(ss))

pcd =  csv_open(DIR1 + "casestudy_tseries.csv")
ax1.plot(xxx, pcd, lw=2, color='navy', label='Case Study')

ax1.legend(loc='lower left', handletextpad=0.5, handlelength=1.25)

# ---- print out spaghetti stats ---- #
print('spaghetti stats')
print('- number of storms: '+str(len(all_lines)))
print('- number that decrease: '+str(len([x for x in all_lines if x[-1]<0])))
print('- number that increase: '+str(len([x for x in all_lines if x[-1]>0])))
print('- percent that increase: '+str( round(len([x for x in all_lines if x[-1]>0])*100/len(all_lines),3) ))

#%% fan plots
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
for li, loc in enumerate(['nh','sh']):
    for yi, years in enumerate(decades):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        mean_lines = csv_open(MAIN+'figure_2/'+loc+'_'+yr_title+'.csv')
        
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

#%% FIGURE 3

## POLAR COORDINATE SCATTER PLOT

#%% FIGURE 4

print(); print('*** GLOBAL BARS ***')

#%% end