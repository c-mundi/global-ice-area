#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 2025
airetemp1.py

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

import functions as fx

import pandas as pd
from scipy.interpolate import griddata
import warnings

ice_fname =  '/Users/mundi/Desktop/seaice/' #'south/'

root = '/Users/mundi/Desktop/month-hemi/'
root_paths = [root+'nh_data/', root+'sh_data/']

hemi_names = ['Arctic', 'Antarctic']

decades = [np.arange(2010,2020), np.arange(1982, 1992)]
months = np.arange(1,12+1)

#%%% get temp function
def get_era5_airtemp(year, loc_ind):
    '''
    Parameters
    ----------
    year : list
       ex. ['2010']


    Returns
    -------
    daily_time : pandas datetime
        daily dates
    daily_avg_slp : list
        lsit of daily averaged slp

    '''
    if type(year) != list:
        year = [str(year)]
        
    months = ["01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"]
    months = months[0:1]
    days = ["01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12",
            "13", "14", "15",
            "16", "17", "18",
            "19", "20", "21",
            "22", "23", "24",
            "25", "26", "27",
            "28", "29", "30",
            "31"]

    import cdsapi
    import xarray as xr
    import io
    from urllib.request import urlopen# start the client
    cds = cdsapi.Client()# dataset you want to read
    dataset = "reanalysis-era5-single-levels"# flag to download data
    params = {
        "data_format": "netcdf",
        "download_format": "unarchived",
        "product_type": ["reanalysis"],
        "variable": ["2m_temperature"],
        'year':year,
        'month':months,
        'day':days,
        "time": ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',
                 '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00',
                 '16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'],
        "grid":[0.25,0.25],
        "area":[90, -180, 60, 180] if loc_ind==0 else [-55, -180, -90, 180]
        }
    # retrieves the path to the file
    fl = cds.retrieve(dataset, params)
    # load into memory
    with urlopen(fl.location) as f:
        bytes_ = f.read()
        ds = xr.open_dataset(io.BytesIO(bytes_), decode_times=True)
        
    # time_in = ds['valid_time'].values
    # hours_in = [(tt-time_in[0])/60/60 for tt in time_in] # seconds -> hours
    # starting_day = datetime(int(year[0]), int(month[0]), int(days[0]))
    # time = np.array([starting_day + timedelta(hours=int(hr)) for hr in hours_in])
    
    # ds.coords['valid_time'] = time
    ds = ds.resample(valid_time='1D').mean()
    
    return ds


#%% sample temps

if False:
    loc_ind = 1
    
    print('Requesting data...')
    ds = get_era5_airtemp(2010, loc_ind)
    print('Retrieved!')

    #### sample plot
    ds_mon = ds.resample(valid_time="MS").mean()
    lon, lat = np.meshgrid(ds['longitude'].values, ds['latitude'].values)
    
    for mm in [1]:
        ds_plot = ds_mon['t2m'].isel(valid_time=mm-1)
    
        if loc_ind==0: 
            fig, ax = fx.background_plot_nh(extent=[-160,90,50,60], returnfig=True, 
                                            labels=False, central_lon=0)
        elif loc_ind==1:
            fig, ax = fx.background_plot_sh(extent=[-180,180, -53,-90], returnfig=True, labels=True)
        
        ax.set_title(calendar.month_name[mm], fontsize=20)
        
        pcm = ax.pcolormesh(lon, lat, ds_plot.values - 273.15, transform=ccrs.PlateCarree(), 
                            cmap=cmo.thermal, vmin=-35,vmax=10)
        
        
        cax1 = fig.add_axes([0.25,0.033,0.4,0.04]) 
        cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
        cbar1.set_label('2-m air temp (C)', fontsize=18)
        cax1.tick_params(labelsize=18)
        
#%% plot storm timeseries

#%%% shade by mean magnitude
cmap = cmo.thermal #cmo.curl

fs=14
xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

def make_grid_plot(nrow, ncol, title='', ylabel=''):
    fontsize=12
    fig, axes_all = plt.subplots(nrow, ncol, sharex=True, sharey=True, figsize=(4*nrow,3.25*ncol))
    fig.suptitle(title, fontsize=fontsize+1)
    alph = iter(['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'])
    
    for ax in axes_all[:,0]:
        ax.set_ylabel(ylabel, fontsize=fontsize)
    for ax in axes_all[-1,:]:
        ax.set_xticks(xxx)
        ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
        ax.set_xlabel('Days Since Storm Start', fontsize=fontsize-1)
    for axl in axes_all:
        for ax1 in axl:
            ax1.set_xlim(-7,14)
            ax1.axhline(0, ls='-', color='k', lw=1)
            ax1.axvline(0, ls='-', color='k', lw=1)
            ax1.set_xticks(xxx)
            ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
            
            ax1.text(0.0225, 0.915, next(alph), transform=ax1.transAxes, fontsize=fontsize, 
                    bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                          'edgecolor':'k', 'lw':0.75},zorder=50)

    return fig, axes_all

mean_swh_lines = []
for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    if loc_ind==0: vmin=-5; vmax=5
    elif loc_ind==1: vmin=-5; vmax=5
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    
    
    mean_lines_mm = {mm:[] for mm in months}
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        fig, axes_all = make_grid_plot(4,3, title='\n\n'+loc+' Air Temperature'+' '+yr_title,
                                       ylabel='[C]')
        
        try:
            t2m_series = fx.era_lines(years, path1+'census/', path1+'area/', path1+'t2m/', 't2m')[-1]
        except FileNotFoundError as fnfe:
            print(fnfe)
            
        axes_w = axes_all.flatten()
        for mi, month in enumerate(months):
            t2m_lines = [np.array(s)-273.15 for s in t2m_series[month] if len(s)==22]
            
            if len(t2m_lines)<10: continue
            
            for line in t2m_lines: axes_w[mi].plot(xxx, line, lw=0.55, color='gray')
                
            mean_line = np.nanmean(t2m_lines, axis=0)
            std_line = np.nanstd(t2m_lines, axis=0)
            
            axes_w[mi].plot(xxx, mean_line, lw=1.75, color='maroon')
            mean_lines_mm[month].append(mean_line)
            
            wave_mean = np.nanmean(mean_line[7:10]) 
            axes_w[mi].axvspan(-7,14, color=cmap(norm(wave_mean)), alpha=0.33, zorder=-500)
            
            axes_w[mi].set_title(calendar.month_name[month]+' (n='+str(len(t2m_lines))+')')
            
        mean_swh_lines.append(mean_lines_mm)
                
        cax1 = fig.add_axes([0.33,-0.033,0.4,0.04]) 
        pcm = axes_w[mi].pcolormesh(np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,2)), 
                                    cmap=cmap, vmin=vmin, vmax=vmax)
        cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
        cbar1.set_label('Mean Air Temp (C)', fontsize=12)
        cax1.tick_params(labelsize=12-1)
        
#%%% shade by mean deviation
cmap = cmo.balance #cmo.curl

fs=14
xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

def make_grid_plot(nrow, ncol, title='', ylabel=''):
    fontsize=12
    fig, axes_all = plt.subplots(nrow, ncol, sharex=True, sharey=True, figsize=(4*nrow,3.25*ncol))
    fig.suptitle(title, fontsize=fontsize+1)
    alph = iter(['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'])
    
    for ax in axes_all[:,0]:
        ax.set_ylabel(ylabel, fontsize=fontsize)
    for ax in axes_all[-1,:]:
        ax.set_xticks(xxx)
        ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
        ax.set_xlabel('Days Since Storm Start', fontsize=fontsize-1)
    for axl in axes_all:
        for ax1 in axl:
            ax1.set_xlim(-7,14)
            ax1.axhline(0, ls='-', color='k', lw=1)
            ax1.axvline(0, ls='-', color='k', lw=1)
            ax1.set_xticks(xxx)
            ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
            
            ax1.text(0.0225, 0.915, next(alph), transform=ax1.transAxes, fontsize=fontsize, 
                    bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                          'edgecolor':'k', 'lw':0.75},zorder=50)

    return fig, axes_all

mean_swh_lines = []
for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    if loc_ind==0: vmin=-0.5; vmax=0.5
    elif loc_ind==1: vmin=-0.5; vmax=0.5
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    
    
    mean_lines_mm = {mm:[] for mm in months}
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        fig, axes_all = make_grid_plot(4,3, title='\n\n'+loc+' Air Temperature'+' '+yr_title,
                                       ylabel='[C]')
        
        try:
            swh_series = fx.era_lines(years, path1+'census/', path1+'area/', path1+'t2m/', 't2m')[-1]
        except FileNotFoundError as fnfe:
            print(fnfe)
            
        axes_w = axes_all.flatten()
        for mi, month in enumerate(months):
            swh_lines = [np.array(s)-273.15 for s in swh_series[month] if len(s)==22]
            
            if len(swh_lines)<10: continue
            
            for line in swh_lines: axes_w[mi].plot(xxx, line, lw=0.55, color='gray')
                
            mean_line = np.nanmean(swh_lines, axis=0)
            std_line = np.nanstd(swh_lines, axis=0)
            
            axes_w[mi].plot(xxx, mean_line, lw=1.75, color='maroon')
            mean_lines_mm[month].append(mean_line)
            
            wave_mean = mean_line[10] - mean_line[7]
            print(loc, era, month, round(wave_mean,3))
            
            axes_w[mi].axvspan(-7,14, color=cmap(norm(wave_mean)), alpha=0.33, zorder=-500)
            
            axes_w[mi].set_title(calendar.month_name[month]+' (n='+str(len(swh_lines))+')')
            
        cax1 = fig.add_axes([0.33,-0.033,0.4,0.04]) 
        pcm = axes_w[mi].pcolormesh(np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,2)), 
                                    cmap=cmap, vmin=vmin, vmax=vmax)
        cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
        cbar1.set_label('Temperature Change (C)', fontsize=12)
        cax1.tick_params(labelsize=12-1)
        
    mean_swh_lines.append(mean_lines_mm)

#%% mean lines analysis? scatter?
from scipy.stats import linregress
mcolors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']


fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex=True)

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    mean_lines_mm = mean_swh_lines[loc_ind]
    
    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        xs, ys = [], []
        for month in months:
            if len(mean_lines_mm[month])==0:continue
            try:
                swh_ml = mean_lines_mm[month][era]
            except IndexError: #???
                continue
                
            swh_dev = swh_ml[10] - swh_ml[7]
            print(loc, era, month, round(swh_dev,3))
            
            ml = mean_lines[month-1][-1]

            axes[loc_ind][era].plot(swh_dev, ml, color=mcolors[month-1],
                                    marker='o', markersize=8)
            xs.append(swh_dev); ys.append(ml)
            
        mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
        m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        ax6a = axes[loc_ind][era].twinx()
        ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.axvline(0, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

for era in [0,1]: axes[-1][era].set_xlabel('Change in Temperature (C)')
for loc in [0,1]: axes[loc][0].set_ylabel('Normalized Change in MIZ Ice Area')

        
#%% end

