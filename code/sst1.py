#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 2025
sst1.py

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


SHIFT_MIZ = True
lat_shift = 2
only_lower_lat = True

STORM_AREA = True

ice_fname =  '/Users/mundi/Desktop/seaice/' #'south/'
sst_path = '/Users/mundi/Desktop/data/'
sst_names = ['sst.wkmean.1981-1989.nc','sst.wkmean.1990-present.nc']

root = '/Users/mundi/Desktop/month-hemi/'
root_paths = [root+'nh_data/', root+'sh_data/']

hemi_names = ['Arctic', 'Antarctic']

decades = [np.arange(2010,2020), np.arange(1982, 1992)]
months = np.arange(1,12+1)

#%%% plotting info
total_days = 0
month_ticks = [0]
for month in months:
    num_days = calendar.monthrange(2010, month)[1]
    total_days += num_days
    month_ticks.append(total_days)
    
xxx = np.arange(0, total_days, 1)
shifted_ticks = np.array(month_ticks[:-1])+15

month_labels = [calendar.month_abbr[mm] for mm in months]+[calendar.month_abbr[months[0]]]

decade_colors = ['#b2182b', '#2166ac']
shade_colors = ['#ef8a62', '#67a9cf']

#%%% functions​

def daterange(start_date, end_date, dt=24):
    alldates=[]
    delta = timedelta(hours=dt)
    while start_date <= end_date:
        alldates.append(start_date)
        start_date += delta
    return alldates


def get_shifted_miz_mask(miz, lat_shift=2):
    method = 'nearest'
    shift_mask = np.zeros(np.shape(si_lon))
    shift_mask = np.where(np.isnan(miz), 0, 1)
    
    if not only_lower_lat:
        lat_neg = si_lat - lat_shift
        si_points_neg = [(x,y) for x,y in zip(si_lon.flatten(), lat_neg.flatten())]
        si_3 = griddata(si_points_neg, miz.flatten(), si_points, method=method)
        si_3 = si_3.reshape(np.shape(si_lat))
        shift_mask = np.where(~np.isnan(si_3), 1, shift_mask)
    
    lat_pos = si_lat + lat_shift
    si_points_pos = [(x,y) for x,y in zip(si_lon.flatten(), lat_pos.flatten())]
    si_2 = griddata(si_points_pos, miz.flatten(), si_points, method=method)
    si_2 = si_2.reshape(np.shape(si_lat))
    shift_mask = np.where(~np.isnan(si_2), 1, shift_mask)
    
    return shift_mask


#%% data

hemi_map = []

for loc_ind, loc in enumerate(hemi_names):

    if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
    elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)
    si_points = [(x,y) for x,y in zip(si_lon.flatten(), si_lat.flatten())]
    
    savepath = root_paths[loc_ind]+'sst/'
    
    decade_map = []
    
    if not SHIFT_MIZ:
        fname =  '_map'
    else: 
        if not only_lower_lat:
            fname = '_shifted_map_'+str(lat_shift)
        else:
            fname = '_lowered_map_'+str(lat_shift)
        
        
    print('Loading Data: '+loc)
    for years in decades:
    
        year_map = []
    
        for year in years:
            try:
                daily_map = np.load(savepath+str(year)+fname+'.npy')  
                # raise FileNotFoundError
            except FileNotFoundError:
                print(year, end='')
                year_timer = timeIN.time()
                
                ### open file
                if year in np.arange(1981,1990): sst_file = sst_path+sst_names[0]
                elif year in np.arange(1990, 2024): sst_file = sst_path+sst_names[1]
                else: raise IndexError('check sst file')
                
                with xr.open_dataset(sst_file, decode_times=True) as ds:
                    lon, lat = np.meshgrid(ds['lon'].values, ds['lat'].values)
                    lon = np.where(lon>180, lon-360, lon)
                    time = ds['time'].values
                    sst = ds['sst']
                    time1 = [pd.to_datetime(str(tx)) for tx in time]
                sst_points = [(x,y) for x,y in zip(lon.flatten(), lat.flatten())]
                
                ### starting info
                start_time = datetime(year, months[0], 1)
                ti = 0
                while time1[ti] <= start_time:
                    ti+=1
                    if ti == len(time1)+5: print(start_time); raise IndexError('missing date in sst file')
                ti -= 1 # go back an index to be within the week
                # print('Starting time:', time[ti])
                sst1 = sst.isel(time = ti).values
                sst_wk = griddata(sst_points, sst1.flatten(), si_points)
                sst_wk = sst_wk.reshape(np.shape(si_lon))
        
                ### start date loop
                daily_map = []
                for month in months:
                    print('.', end='')
                    
                    for day in np.arange(1,calendar.monthrange(year, month)[-1]+1):
                        date = datetime(year, month, day)
                        try:
                            week_range = daterange(time1[ti], time1[ti+1])[:-1] # index behind
                        except IndexError:
                            if month==12 and year==1989: week_range = daterange(time1[ti], datetime(1990,1,1))
                            else: week_range = daterange(time1[ti], time1[ti]+timedelta(days=6))
                        # print(date.strftime('%Y-%m-%d')+': ' + str(date not in week_range) )
                        # print('-- '+ week_range[0].strftime('%Y-%m-%d')+' : '+week_range[-1].strftime('%Y-%m-%d'))
                        
                        if date not in week_range:
                            ### load new sea ice
                            ti += 1
                            sst1 = sst.isel(time = ti).values
                            # print('-- Loading new sst ('+ str(time[ti])+')')
                            ### interp grid
                            sst_wk = griddata(sst_points, sst1.flatten(), si_points)
                            sst_wk = sst_wk.reshape(np.shape(si_lon))
                           
                        ### get miz
                        if month==2 and day==29: continue
                        if loc_ind ==0: si = fx.load_seaice(ice_fname, year, month, day, latlon=False)
                        elif loc_ind==1: si = fx.load_seaice_sh(ice_fname, year, month, day, latlon=False)
                        
                        if np.nanmean(si)!=0:
                            si = np.where(si_lat>-55, np.nan, si)
                            sie = np.where(si<0.15, np.nan, si)
                            miz = np.where(sie>0.80, np.nan, sie) 
                            
                            ### append sst
                            if SHIFT_MIZ:
                                shift_mask = get_shifted_miz_mask(miz, lat_shift)
                                daily_map.append( np.where(shift_mask==1, sst_wk, np.nan) )
                            else:
                                daily_map.append( np.where(np.isnan(miz), np.nan, sst_wk) )
                            
                        else:
                            daily_map.append( np.nan*np.ones(np.shape(si_lon)) )
                   
                print(' '+str(round((timeIN.time()-year_timer)/60,1))+' min')
    
                np.save(savepath+str(year)+fname+'.npy', daily_map)  
                
            year_map.append(daily_map)
    
        decade_map.append(year_map)
    
    decade_map = np.array(decade_map)
    
    hemi_map.append(decade_map)

print('Done')

#%% timeseries
print('Plotting time series - ', end='')

for loc_ind, decade_map in enumerate(hemi_map):

    #### set up
    total_days = 0
    month_ticks = [0]
    for month in months:
        num_days = calendar.monthrange(2010, month)[1]
        total_days += num_days
        month_ticks.append(total_days)
    x_days = np.arange(0, total_days, 1)
    shifted_ticks = np.array(month_ticks[:-1])+15
    
    era_colors = ['#ca0020', '#0571b0']
    
    def plot_spread(ax, x, mean_line, std_line, color, zorder=0):
        ax.plot(x, mean_line, color=color, zorder=zorder)
        ax.plot(x, mean_line + std_line, color=color, ls='--', zorder=zorder)
        ax.plot(x, mean_line - std_line, color=color, ls='--', zorder=zorder)
        ax.fill_between(x, mean_line, mean_line+std_line, 
                          alpha=0.275, color=color, zorder=zorder)
        ax.fill_between(x, mean_line, mean_line-std_line, 
                            alpha=0.275, color=color, zorder=zorder)
    
    
    #### plot data
    fig, axt = plt.subplots(1,1, figsize=(8,4))
    axt.set_title(hemi_names[loc_ind]+': Daily Mean SST in the MIZ')
    axt.set_ylabel(r'SST ($^\circ$C)')
    for mt in month_ticks: axt.axvline(mt, color='gray', ls=':',zorder=-20)
    
    for era, years in enumerate(decades):
        ystr = str(years[0])+'-'+str(years[-1])
        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            monthly_data = np.nanmean(decade_map[era], axis=(-2,-1))
        
        mean_line = np.nanmean(monthly_data,axis=0)
        std_line = np.nanstd(monthly_data, axis=0)
        
        # axt.plot(x_days, mean_line, color=era_colors[era])
        plot_spread(axt, x_days, mean_line, std_line, era_colors[era])
        
        mean_path = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/z-monthly-aa-comparison/'
        np.save(mean_path+'aa_'+str(years[0])+'-'+str(years[-1])+'_mean_sst.npy', mean_line)
    
    axt.set_xticks(shifted_ticks, month_labels[:-1])


print('Done')

#%% sst histograms
print('Plotting histograms... ', end='')

for loc_ind, decade_map in enumerate(hemi_map):

    fig_h, axes_h = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
    axes_h = axes_h.flatten()
    for ai, axh in enumerate(axes_h): axh.set_title(calendar.month_name[ai+1])
    bins = np.arange(-2.5,2.5,0.2)
    
    if loc_ind==1:
        if not SHIFT_MIZ: txt_height = 30
        else: txt_height = 19
    elif loc_ind==0:
        if not SHIFT_MIZ: txt_height = 9 #11
        else: txt_height = 7
    
    for era, years in enumerate(decades):
        if era==0:
            xx = -2.5
        elif era==1: 
            xx = 1
        ystr = str(years[0])+'-'+str(years[-1])
        axes_h[-1].plot([],[], color = shade_colors[era], label=ystr, lw=3, alpha=0.66)
        
        monthly_data = decade_map[era]
        mind=0
        for mi, month in enumerate(months):
            num_days = calendar.monthrange(2010, month)[1]
            mslice = monthly_data[:,mind:mind+num_days, :,:].flatten()
            
            hist, bin_edges = np.histogram(mslice, bins=bins)
            axes_h[mi].bar(bin_edges[:-1], (hist/len(mslice[~np.isnan(mslice)]))*100, width=0.2, 
                           facecolor=shade_colors[era],edgecolor=decade_colors[era], alpha=0.66)
            
            axes_h[mi].text(xx, txt_height, str(round(np.nanmean(mslice),1))+r'$^\circ$', color=decade_colors[era])
            
            mind+=num_days
    
    fig_h.text(0.05, 0.5, 'Frequency (%)',  va='center', ha='center', rotation=90, weight='bold')
    fig_h.text(0.5, 0.05, r'SST ($^\circ$C)',  va='center', ha='center', weight='bold')
    
    axes_h[-1].legend(ncol=2, bbox_to_anchor=(0.95, -0.175), handletextpad=0.5, handlelength=1)
    
    if loc_ind==1: fig_h.suptitle(hemi_names[loc_ind]+'\n'+'SSTs in the MIZ have decreased in the recent decade')
    if loc_ind==0: fig_h.suptitle(hemi_names[loc_ind]+'\n'+'SSTs in the MIZ have decreased in the recent decade')

print('Done')

#%% timeseries???
from concurrent.futures import ThreadPoolExecutor

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]


def calc_sst_area(start):
    if start.month==2 and start.day==29:
        start = datetime(start.year, 3, 1)
    try:
        yr_ind1 = [i for i,dt in enumerate(years) if dt==start.year][0]
    except IndexError:
        print('nan', start.year, years[0], years[1])
        return np.nan
    
    daily_range = fx.daterange(datetime(2010, 1,1), datetime(2010,12,31), dt=24)
    dt_ind1 = [i for i,dt in enumerate(daily_range) if (dt.month == start.month and dt.day == start.day)][0]
    
    try: sst_start = monthly_data[yr_ind1, dt_ind1, :,:]
    except IndexError as ie: 
        print(start, yr_ind1, dt_ind1)
        print(ie)
        sst_start = monthly_data[yr_ind1, -1, :,:]
        
    if not STORM_AREA:
        return np.nanmean(sst_start) 
    else:
        si_lon2 = np.where(si_lon<0, si_lon+360, si_lon)
        inside_points = fx.find_points_in_contour(bbox, si_lon2, si_lat)
        start_bbox = np.where(inside_points, np.nan, sst_start)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return np.nanmean(start_bbox) 

def get_tseries(STORMINFO): # sst_scatter
    ''' computes sst difference'''
    start, end, years, bbox, si_lon, monthly_data = STORMINFO
    
    week_ago = start - timedelta(days=7)
    two_week = start + timedelta(days=14) # relative to start date, since different storm lengths
    analysis_range = daterange(week_ago, two_week, dt=24)
    
    def get_sst_threaded(analysis_range):
        with ThreadPoolExecutor(max_workers=4) as executor:
            return executor.map(calc_sst_area, analysis_range)

    results = get_sst_threaded(analysis_range)
    SERIES = [x for x in results]
            
    return SERIES


# def get_tseries(STORMINFO): # sst_scatter
#     ''' computes sst difference'''
#     start, end, years, bbox, si_lon = STORMINFO
    
#     week_ago = start - timedelta(days=7)
#     two_week = start + timedelta(days=14) # relative to start date, since different storm lengths
#     analysis_range = daterange(week_ago, two_week, dt=24)
    
#     SERIES = []
#     for start in analysis_range:
#         if start.month==2 and start.day==29:
#             start = datetime(start.year, 3, 1)
            
#         try:
#             yr_ind1 = [i for i,dt in enumerate(years) if dt==start.year][0]
#         except IndexError:
#             SERIES.append(np.nan)
#             continue
        
#         daily_range = fx.daterange(datetime(2010, 1,1), datetime(2010,12,31), dt=24)
#         dt_ind1 = [i for i,dt in enumerate(daily_range) if (dt.month == start.month and dt.day == start.day)][0]
        
#         try: sst_start = monthly_data[yr_ind1, dt_ind1, :,:]
#         except IndexError as ie: 
#             print(start, yr_ind1, dt_ind1)
#             print(ie)
#             sst_start = monthly_data[yr_ind1, -1, :,:]
            
#         if not STORM_AREA:
#             SERIES.append( np.nanmean(sst_start) )
#         else:
#             si_lon2 = np.where(si_lon<0, si_lon+360, si_lon)
#             inside_points = fx.find_points_in_contour(bbox, si_lon2, si_lat)
#             start_bbox = np.where(inside_points, np.nan, sst_start)
#             with warnings.catch_warnings():
#                 warnings.simplefilter('ignore')
#                 SERIES.append( np.nanmean(start_bbox) )
            
#     return SERIES

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

#%%% tseries: data and plot
from scipy.stats import linregress
# cmap = cmo.thermal #cmo.curl
cmap = cmo.balance

start_timer = timeIN.time()
print(); print('Starting timeseries plots (this may take a few mins!)') ###!!!! can this be streamlined??? parallel/save data??

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    decade_map = hemi_map[loc_ind]
    
    if loc_ind==0: 
        _, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
        # norm = plt.Normalize(vmin=0, vmax=2.5) # mean sst
        norm = plt.Normalize(vmin=-0.33, vmax=0.33) # change in sst
    elif loc_ind==1: 
        _, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)
        # norm = plt.Normalize(vmin=-1, vmax=-0.5)# mean sst
        norm = plt.Normalize(vmin=-0.1, vmax=0.1) # change in sst
    
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        fig, axes_all = make_grid_plot(4,3, title='\n\n'+loc+' SST'+' '+yr_title,
                                       ylabel=r'SST ($^\circ$C)')
        print(loc, yr_title)

        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        storm_areas = fx.get_storm_bbox(years, path1+'census/', path1+'area/', path1+'contours/')
        
        axes_s = axes_all.flatten()
        monthly_data = decade_map[era]
        for mi, mm in enumerate(list(lines.keys())):
            monthly_starts = start_day[mm]
            monthly_ends = end_day[mm]
            areas = storm_areas[mm]
            
            # organize SSTs
            try:
                SSTs = np.load(path1+'sst/'+'tseries_'+str(era)+'-'+str(mm)+'.npy')
                for sst_series in SSTs: axes_s[mi].plot(xxx, sst_series, lw=0.55, color='gray')
            except FileNotFoundError:
                print('-', mm); month_timer = timeIN.time()
                SSTs = []
                for start, end, bbox in zip(monthly_starts, monthly_ends, areas):
                    STORMINFO = [start, end, years, bbox, si_lon, monthly_data]
                    sst_series = get_tseries(STORMINFO)
                    SSTs.append(sst_series)
                   
                    axes_s[mi].plot(xxx, sst_series, lw=0.55, color='gray')
                    
                print('--', np.shape(np.array(SSTs)))
                print(len(np.unique(SSTs)), np.nanmean(SSTs))
                np.save(path1+'sst/'+'tseries_'+str(era)+'-'+str(mm)+'.npy', np.array(SSTs))
                print('...', str(round((timeIN.time()-month_timer)/60,1)), 'min')
                
            ### plot
            if len(SSTs)>10:
                mean_line = np.nanmean(SSTs, axis=0)
                std_line = np.nanstd(SSTs, axis=0)
                
                m, b, r, p, se = linregress(xxx, mean_line)
                
                axes_s[mi].plot(xxx, mean_line, lw=1.75, color='maroon', 
                                label=str(round(m*1e3,2))+'e-3')
                
                sst_mean = np.nanmean(mean_line[14] - mean_line[7]) 
                axes_s[mi].axvspan(-7,14, color=cmap(norm(sst_mean)), alpha=0.33, zorder=-500)
            
                axes_s[mi].legend(loc='upper right', fontsize=12)
            
            axes_s[mi].set_title(calendar.month_name[mm]+' (n='+str(len(SSTs))+')')
            
        cax1 = fig.add_axes([0.33,-0.033,0.4,0.04]) 
        pcm = axes_s[mi].pcolormesh(np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,2)), 
                                    cmap=cmap, vmin=norm.vmin, vmax=norm.vmax)
        cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
        cbar1.set_label('Change in SST (degC)', fontsize=12)
        cax1.tick_params(labelsize=12-1)

            
     
tdiff = (timeIN.time() - start_timer)/60
print('Finished time series plots:', round(tdiff, 2), 'minutes')

#%% simple scatter (sst trend, si change)
print(); print('Simple scatter: sst vs. si')
mcolors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex=True)

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    decade_map = hemi_map[loc_ind]

    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))

        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
         
        monthly_data = decade_map[era]
        xs, ys = [], []
        for mi, mm in enumerate(list(lines.keys())):
            ml = mean_lines[mi][-1]
            
            # organize SSTs
            SSTs = np.load(path1+'sst/'+'tseries_'+str(era)+'-'+str(mm)+'.npy')
                
            ### plot
            if len(SSTs)>10:
                mean_line = np.nanmean(SSTs, axis=0)
                std_line = np.nanstd(SSTs, axis=0)
                slope, b, r, p, se = linregress(xxx, mean_line)
                
                axes[loc_ind][era].plot(slope, ml, color=mcolors[mi],
                                        marker='o', markersize=8)
                xs.append(slope); ys.append(ml)
                
        m, b, r, p, se = linregress(xs, ys)
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        ax6a = axes[loc_ind][era].twinx()
        ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.axvline(0, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

for era in [0,1]: axes[-1][era].set_xlabel(r'Mean Change in SST ($^\circ$C day$^{-1}$)')
for loc in [0,1]: axes[loc][0].set_ylabel('Normalized Change in MIZ Ice Area')

#%% calculate climatology
import os

def get_clim_series(STORMINFO):
    [start, end, years, bbox, si_lon, monthly_data] = STORMINFO
    
    daily_range = fx.daterange(datetime(2010, 1,1), datetime(2010,12,31), dt=24) #reference
    
    # get storm area
    si_lon2 = np.where(si_lon<0, si_lon+360, si_lon)
    inside_points = fx.find_points_in_contour(bbox, si_lon2, si_lat)
    
    # daily loop 
    clim_series = []
    for dt in fx.daterange(start-timedelta(days=7), start+timedelta(days=14), dt=24):
        if dt.month==2 and dt.day==29:
            dt = datetime(dt.year, 3, 1)
        try:
            yr_ind1 = [i for i,dt1 in enumerate(years) if dt1==dt.year][0]
        except IndexError:
            print('nan', dt.year, years[0], years[1])
            return np.nan
        
        dt_ind1 = [i for i,dt1 in enumerate(daily_range) if (dt1.month == dt.month and dt1.day == dt.day)][0]
        
        try: sst_start = monthly_data[yr_ind1, dt_ind1, :,:]
        except IndexError as ie: 
            print(dt, yr_ind1, dt_ind1)
            print(ie)
            sst_start = monthly_data[yr_ind1, -1, :,:]
            
        ### storm area  # decade map (monthly_data) isolates miz (and extends latitude?)
        start_bbox = np.where(inside_points, np.nan, sst_start)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mean_sst = np.nanmean(start_bbox) 
        clim_series.append(mean_sst)
        
    return clim_series


#%%% data and plot
start_timer = timeIN.time()
print(); print('Timeseries plots - Minus Climatology ') 

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    decade_map = hemi_map[loc_ind]
    
    clim_path = path1+'sst/clim/'
    if not os.path.exists(clim_path):
        os.makedirs(clim_path)
    
    if loc_ind==0: 
        _, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
        # norm = plt.Normalize(vmin=0, vmax=2.5) # mean sst
        norm = plt.Normalize(vmin=-0.33, vmax=0.33) # change in sst
    elif loc_ind==1: 
        _, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)
        # norm = plt.Normalize(vmin=-1, vmax=-0.5)# mean sst
        norm = plt.Normalize(vmin=-0.33, vmax=0.33) # change in sst
    
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        fig, axes_all = make_grid_plot(4,3, title='\n\n'+loc+' SST - clim'+' '+yr_title,
                                       ylabel=r'SST ($^\circ$C)')
        print(loc, yr_title)

        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        storm_areas = fx.get_storm_bbox(years, path1+'census/', path1+'area/', path1+'contours/')
        
        axes_s = axes_all.flatten()
        monthly_data = decade_map[era]
        for mi, mm in enumerate(list(lines.keys())):
            monthly_starts = start_day[mm]
            monthly_ends = end_day[mm]
            areas = storm_areas[mm]
            
            # load sst
            SSTs = np.load(path1+'sst/'+'tseries_'+str(era)+'-'+str(mm)+'.npy')
            
            # get climatology
            try:
                clims = np.load(clim_path+'series_'+str(era)+'-'+str(mm)+'.npy')
            except FileNotFoundError:
                print('-', mm); month_timer = timeIN.time()
                clims = []
                for start, end, bbox in zip(monthly_starts, monthly_ends, areas):
                    STORMINFO = [start, end, years, bbox, si_lon, monthly_data]
                    clim_series = get_clim_series(STORMINFO)
                    clims.append(sst_series)
                   
                    axes_s[mi].plot(xxx, sst_series, lw=0.55, color='gray')
                    
                print('--', np.shape(np.array(clims)))
                print(len(np.unique(clims)), np.nanmean(clims))
                np.save(clim_path+'series_'+str(era)+'-'+str(mm)+'.npy', np.array(clims))
                print('...', str(round((timeIN.time()-month_timer)/60,1)), 'min')
                
                
                
            ### plot
            for sst_series, clim_series in zip(SSTs, clims): 
                axes_s[mi].plot(xxx, sst_series-clim_series, lw=0.55, color='gray')
            
            if len(SSTs)>10:
                mean_line = np.nanmean(SSTs-clims, axis=0)
                std_line = np.nanstd(SSTs-clims, axis=0)
                
                m, b, r, p, se = linregress(xxx, mean_line)
                
                axes_s[mi].plot(xxx, mean_line, lw=1.75, color='maroon', 
                                label=str(round(m*1e3,2))+'e-3')
                
                sst_mean = np.nanmean(mean_line[14] - mean_line[7]) 
                axes_s[mi].axvspan(-7,14, color=cmap(norm(sst_mean)), alpha=0.33, zorder=-500)
            
                axes_s[mi].legend(loc='upper right', fontsize=12)
            
            axes_s[mi].set_title(calendar.month_name[mm]+' (n='+str(len(SSTs))+')')
            
        cax1 = fig.add_axes([0.33,-0.033,0.4,0.04]) 
        pcm = axes_s[mi].pcolormesh(np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,2)), 
                                    cmap=cmap, vmin=norm.vmin, vmax=norm.vmax)
        cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
        cbar1.set_label('Change in SST (degC)', fontsize=12)
        cax1.tick_params(labelsize=12-1)

            
     
tdiff = (timeIN.time() - start_timer)/60
print('Finished time series plots:', round(tdiff, 2), 'minutes')


#%% end
