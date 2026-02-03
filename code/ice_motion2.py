#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 2025
ice_motion2.py

> apply ice motion to storm cases

@author: mundi
"""

#%% imports and file paths
import numpy as np
import xarray as xr
from glob import glob
from datetime import datetime, timedelta
import os
import functions as fx
import warnings

ice_fname = '/Users/mundi/Desktop/seaice/'
ice_path = '/Users/mundi/Desktop/seaice/ice_motion/'

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

savedir = 'icemotion/'

census_name = 'census_'
contour_name = '_contours.nc'
si_name = '_seaice'

decades = [np.arange(1982,1992), np.arange(2010,2020)]
months = np.arange(1,12+1)
hemi_names= ['Arctic', 'Antarctic']

from scipy.stats import linregress
month_colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

#%%% functions
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


def load_ice_motion(ice_path, loc_ind, year, month, day):
    
    fname = ['nh','sh'][loc_ind]+ '_25km_'+str(year)+'0101_'+str(year)+'1231_v4.1.nc'
    file = glob(ice_path+'*'+fname)[0]
    ds = xr.open_dataset(file)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        datetimeindex = ds.indexes['time'].to_datetimeindex(time_unit='s')
        ds['time'] = datetimeindex
    
    lon = ds['longitude'].values
    lat = ds['latitude'].values

    u = ds['u'].sel(time=slice(datetime(year, month, day), datetime(year, month, day))).values
    v = ds['v'].sel(time=slice(datetime(year, month, day), datetime(year, month, day))).values
    
    # transform to lat/lon
    # https://nsidc.org/data/user-resources/help-center/how-convert-horizontal-and-vertical-components-east-and-north
    if loc_ind==0:
        E = ( u*np.cos(np.deg2rad(lon)) ) + ( v*np.sin(np.deg2rad(lon)) )
        N = ( -1*u*np.sin(np.deg2rad(lon)) ) + ( v*np.cos(np.deg2rad(lon)) )
    elif loc_ind==1:
        E = ( u*np.cos(np.deg2rad(lon)) )  -  ( v*np.sin(np.deg2rad(lon)) )
        N = ( u*np.sin(np.deg2rad(lon)) )  +  ( v*np.cos(np.deg2rad(lon)) )
        
    return lon, lat, E, N

#%%% calculations

def get_miz_area(loc_ind, storm_event, lon, lat):
    from scipy.interpolate import griddata
    miz=(0.15,0.80)
    
    if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010,1,1, latlon=True)
    elif loc_ind == 1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010,1,1, latlon=True)
    
    t1 = storm_event[0] - timedelta(days=1)
    t2 = storm_event[-1] + timedelta(days=1)
    storm_range = fx.daterange(t1, t2, dt=24)
    
    miz_points = np.zeros(np.shape(si_lon))
    for date in storm_range:
        if loc_ind==0: sic = fx.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
        elif loc_ind == 1:sic = fx.load_seaice_sh(ice_fname, date.year, date.month, date.day, latlon=False)
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
        
    miz_grd = griddata((si_lon.flatten(), si_lat.flatten()), miz_points.flatten(),
                                (lon.flatten(), lat.flatten()), method='nearest')
    miz_grd = miz_grd.reshape(np.shape(lon))
    return miz_grd #miz_points

def get_storm_areas(year, root_path):
    ice_lims = [20,80]
    
    census_file = root_path+'census/'+census_name+str(year)+'.csv'
    [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
    
    timing_grid = []
    for xx in range(0,len(startdate)):
        timing_grid.append((startdate[xx], enddate[xx]))
    storm_ranges = []
    for startdt, enddt in timing_grid:
        storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
    
    # open ice area
    ds_area = xr.open_dataset(root_path+'area/' + str(year) +'_area.nc')
    ice_sorter = ds_area['ice_area80'].values
    box_area = ds_area['box_area'].values
    
    try:
        ds = xr.open_dataset(root_path+'seaice/' + str(year) + si_name + '.nc')
        ds.close()
    except:
        print('- skip: '+root_path+'seaice/' + str(year) + si_name + '.nc')
        return [], []
    
    storm_bbox = []
    stormstr_prev=''
    for storm_num, strm in enumerate(storm_ranges):
        # remove storms that don't interact with the ice
        ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
        if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
            storm_bbox.append([])
            continue
        
        # load storm contours
        stormstr1 = strm[0].strftime('%Y_%m%d')
        # duplicate storm start date?
        if stormstr1==stormstr_prev:
            stormstr = stormstr1 + '*'
            dupe=True
        else:
            dupe=False
            stormstr=stormstr1
        stormstr_prev = stormstr1
        
        #### storm area
        ncname = stormstr + contour_name +'.nc'
        try:
            if not dupe: ###!!!
                cs = xr.open_dataset(root_path+'contours/'+ncname)
            else:
                raise FileNotFoundError
        except FileNotFoundError:
            contour_files = glob(root_path+'contours/'+'*'+strm[0].strftime('%Y_%m%d')+'*.nc')
            cs = xr.open_dataset(contour_files[-1])
            
        all_contours = []
        for key in list(cs.keys()):
            coord = cs[key].values
            all_contours.append(coord)
        cs.close()
        
        if all_contours == []: 
            print('-', stormstr) 
            storm_bbox.append([])
            continue
        
        with fx.HidePrint(): bbox_edges = fx.get_bbox_edges(all_contours) 
        
        storm_bbox.append( bbox_edges )
        
    return storm_bbox, storm_ranges

def calc_motion(data):
    dt = data[0]
    miz_points = data[1]
    inside_bbox = data[2]
    
    # load values
    lon, lat, E, N = load_ice_motion(ice_path, loc_ind, dt.year, dt.month, dt.day)
    total = np.sqrt( (E**2) + (N**2) )
    
    # ice speed timeseries
    series_out = []
    for tw, this in enumerate([total, E, N]):
         in_tot = np.ma.masked_array(np.squeeze(this), mask=inside_bbox).filled(np.nan)
         total_miz = np.ma.masked_where(miz_points==0, in_tot).filled(np.nan)
         
         with warnings.catch_warnings():
             warnings.simplefilter('ignore')
             series_out.append(np.nanmean(total_miz))
    return series_out
    

#%% data loop

@timethis
def yearly_data(year, loc_ind, ice_path):
    print();print(year)
    name = 'icemotion_'+str(year)
    
    try:
        motion_series = np.load(savepath + name +'.npy')
    except FileNotFoundError:
        lon, lat, E, N = load_ice_motion(ice_path, loc_ind, 2010,1,1)
        
        storm_bbox, storm_ranges = get_storm_areas(year, path1)
        
        motion_series = []
        for bbox, storm_event in zip(storm_bbox, storm_ranges):
            if len(bbox)==0:
                motion_series.append([[np.nan]*3]*22)
                continue
            
            ### storm info
            week_ago = storm_event[0] - timedelta(days=7)
            two_week = storm_event[0] + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_range = fx.daterange(week_ago, two_week, dt=24)
            
            ### miz  area
            miz_points = get_miz_area(loc_ind, storm_event, lon, lat)
            inside_bbox = fx.find_points_in_contour(bbox, lon, lat)
           
            # make data list for function
            data = []
            for dt in analysis_range:
                data.append([dt, miz_points, inside_bbox])
       
            ### calculate daily values
            def run_threaded_timeseries(data):
                with ThreadPoolExecutor(max_workers=7) as executor:
                    return executor.map(calc_motion, data)
                
            tseries = run_threaded_timeseries(data)
            motion_series.append( [x for x in tseries] )
        
        np.save(savepath+name +'.npy', motion_series)
        print('SAVED: '+ name +'.npy')


for loc_ind, loc in enumerate(hemi_names):
    print(loc)
    path1 = root_paths[loc_ind]
    savepath = path1+savedir
    if not os.path.exists(savepath):
        os.makedirs(savepath)
        
    for years in decades:
        for year in years:
            yearly_data(year, loc_ind, ice_path)
            
#%%% processing
def storm_ice_motion(years, census_path, area_path, icemotion_path):
    ice_lims=[20,80]
    start_day = {}
    end_day = {}
    motion = {}
    
    months = np.arange(1,12+1)
    for mm in months:
        start_day[mm] = []
        end_day[mm] = []
        motion[mm] = []
        
    for year in years:
        # open yearly ice motion file
        motion_series = np.load(icemotion_path + 'icemotion_'+str(year) +'.npy')
        
        # storm info
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in timing_grid:
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(fx.daterange(week_ago, two_week, dt=24))
            storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
        for storm_num, strm in enumerate(storm_ranges):
             month = int(strm[0].month)
             
             ### remove storms that don't interact with the ice
             ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
             if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                continue
        
             start_day[month].append(strm[0])
             end_day[month].append(strm[-1])
             
             motion[month].append(motion_series[storm_num])
        
    return start_day, end_day, motion


#%% ---
#%% plots
import cmocean.cm as cmo
import matplotlib.pyplot as plt
import string
import calendar

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

wnames0 = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 
          'xlims':[[-4.5,3],[-1,3]], 'vline':0}
wnames1 = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 
          'xlims':[[-4.5,2],[-1,2]], 'vline':0}
wnames2 = {'idx':0, 'ylims':[-2,14], 'name':'Total', 
          'xlims':[[4,10],[6,10]], 'vline':None}

wnames = [wnames0, wnames1, wnames2][1]
          
#%%% monthly plot

cmap = cmo.balance #cmo.curl
norm = plt.Normalize(vmin=wnames['ylims'][0], vmax=wnames['ylims'][1])

def make_grid_plot(nrow, ncol, title='', ylabel=''):
    fontsize=12
    fig, axes_all = plt.subplots(nrow, ncol, sharex=True, sharey=True, figsize=(4*nrow,3.25*ncol))
    fig.suptitle(title, fontsize=fontsize+1)
    alph = iter(list(string.ascii_lowercase))
    
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

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        fig, axes_all = make_grid_plot(4,3, title='\n\n'+loc+' '+wnames['name']+' Ice Motion'+' '+yr_title,
                                       ylabel=r'[m s$^{-1}$]')
        
        ### get ice motion lines
        motion_series = storm_ice_motion(years, path1+'census/', path1+'area/', path1+'icemotion/')[-1]
            
        axes_w = axes_all.flatten()
        for mi, month in enumerate(months):
            
            moves = [arr for arr in motion_series[month] if len(arr)==22]  
            im_lines = np.array(moves)[:,:,wnames['idx']]
            
            if len(im_lines)<10: continue
            
            for line in im_lines: axes_w[mi].plot(xxx, line, lw=0.55, color='gray')
                
            mean_line = np.nanmean(im_lines, axis=0)
            std_line = np.nanstd(im_lines, axis=0)
            
            axes_w[mi].plot(xxx, mean_line, lw=1.75, color='maroon')
            
            im_mean = np.nanmean(mean_line[7:10]) 
            axes_w[mi].axvspan(-7,14, color=cmap(norm(im_mean)), alpha=0.33, zorder=-500)
            
            axes_w[mi].set_title(calendar.month_name[month]+' (n='+str(len(im_lines))+')')
            axes_w[mi].set_ylim(wnames['ylims'])          

#%%% simple scatter
figs, axs = plt.subplots(2,2, figsize=(8,8), sharey='row')
figs.suptitle('Monthly '+wnames['name']+' Ice Motion')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        axs[loc_ind][era].set_title(loc+': '+yr_title)
        # axs[loc_ind][era].set_xlim(wnames['xlims'][loc_ind])
        axs[loc_ind][era].axhline(0, lw=0.55, color='gray', ls=':')
        axs[loc_ind][era].axvline(0, lw=0.55, color='gray', ls=':')
        
        motion_series = storm_ice_motion(years, path1+'census/', path1+'area/', path1+'icemotion/')[-1]
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        xs, ys = [], []
        for mi, month in enumerate(months):
            marker='o'
            
            moves = [arr for arr in motion_series[month] if len(arr)==22]  
            im_lines = np.array(moves)[:,:,wnames['idx']]
            
            if len(im_lines)<10: 
                continue
            
            mean_im_line = np.nanmean(im_lines, axis=0)
            im_mean = np.nanmean(mean_im_line[7:14])
            
            mean_impact = mean_lines[mi][-1]
            
            axs[loc_ind][era].plot(im_mean, mean_impact, 
                                   marker=marker, markersize=10, color=month_colors[mi])
            xs.append(im_mean)
            ys.append(mean_impact)
            
        #### lines of best fit? correlations?
        m, b, r, p, se = linregress(xs,ys)
        ax6 =  axs[loc_ind][era]
        x = np.linspace(np.nanmin(xs), np.nanmax(xs), 50)
        ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
        ax6a = ax6.twinx()
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

# final plot
for era, years in enumerate(decades): axs[-1][era].set_xlabel('Ice Motion')
for loc_ind, loc in enumerate(hemi_names): axs[loc_ind][0].set_ylabel('Normalized Change in MIZ Ice Area')



#%% ---
#%% end
