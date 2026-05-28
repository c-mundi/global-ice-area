#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 2026
storm_strength.py

-  assess changes in storm strength between two decades

@author: mundi
"""
#%% file paths and imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import calendar
import time as timeIN
from glob import glob
import warnings

import functions as fx
from functools import wraps
import pickle

ROOT = '/Users/mundi/Desktop/'
# ROOT = '/home/mundi/'

nh_path = ROOT+'month-hemi/nh_data/'
sh_path = ROOT+'month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

census_name = 'census_'
contour_name = '_contours'
si_name = '_seaice'

decades = [np.arange(1982, 1992), np.arange(2010,2020)] 

# YEARS = [2010,2011,2012]
# YEARS = [2013,2014,2015]
# YEARS = [2016,2017,2018]
# YEARS = [2019,1982,1983]
# YEARS = [1985,1986,1987]
# YEARS = [1988,1989,1990,1991]

YEARS = np.arange(2010,2020) #np.arange(1982,1992)

months = np.arange(1,12+1)

#%%% functions
def get_era5_wind(year, month, days, LOC):
    '''
    Parameters
    ----------
    year : list
       ex. ['2010']
    month : list
        ex. ['08']
    day : lsit
        ex. ['15','16']

    Returns
    -------
    daily_time : pandas datetime
        daily dates
    daily_avg_slp : list
        lsit of daily averaged slp

    '''
    if type(year) != list:
        year = [str(year)]
    if type(month) != list:
        if month < 10: monthstr = '0'+str(int(month))
        else: monthstr = str(int(month))
        month = [monthstr]
    if type(days) != list and type(days[0])!=str:
        days = ['0'+str(int(day)) if day<10 else str(int(day)) for day in days]

    import cdsapi
    import xarray as xr
    import io
    from urllib.request import urlopen# start the client
    cds = cdsapi.Client()# dataset you want to read
    dataset = "reanalysis-era5-single-levels"# flag to download data
    download_flag = False# api parameters 
    params = {
        "data_format": "netcdf",
        "download_format": "unarchived",
        "product_type": ["reanalysis"],
        "variable": [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind"
                    ],
        'year':year,
        'month':month,
        'day':days,
        "time": ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',
                 '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00',
                 '16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'],
        "grid":[0.25,0.25],
        "area":[90, -180, 60, 180] if LOC==0 else [-60, -180, -90, 180]
        }
    # retrieves the path to the file
    fl = cds.retrieve(dataset, params)
    # download the file 
    if download_flag:
        fl.download("./output.nc")
    # load into memory
    with urlopen(fl.location) as f:
        bytes_ = f.read()
        ds = xr.open_dataset(io.BytesIO(bytes_), decode_times=False)
        
    time_in = ds['valid_time'].values
    hours_in = [(tt-time_in[0])/60/60 for tt in time_in] # seconds -> hours
    starting_day = datetime(int(year[0]), int(month[0]), int(days[0]))
    time = np.array([starting_day + timedelta(hours=int(hr)) for hr in hours_in])
    
    ds.coords['valid_time'] = time
    ds = ds.resample(valid_time='1D').mean()
    
    return ds

def LZ(day):
    if day>=10:
        return str(day)
    elif day<10:
        return '0'+str(day)
    
def timethis(func):
    """ 
    Print the execution time for a function call
    """
    @wraps(func)
    def wrapped_method(*args, **kwargs):
        time_start = timeIN.time()
        output = func(*args, **kwargs)
        time_end = timeIN.time()
        if time_end-time_start < 120:
            print(f"{func.__name__}: {(time_end-time_start)} s")
        else:
            print(f"{func.__name__}: {(time_end-time_start)/60} min")

        return output

    return wrapped_method

def check_significance(data1, data2):
    from scipy.stats import ttest_ind, ks_2samp, mannwhitneyu

    t, pt = ttest_ind(data1, data2)
    ks, pks = ks_2samp(data1, data2)
    mw, pmw = mannwhitneyu(data1, data2)
     
    pvals = [pt, pks, pmw]
    conf_levels = [100*(1-p) for p in pvals]
    
    return np.array(conf_levels)

def set_box_color(bp, color, lw=1.5):
    plt.setp(bp['boxes'], color=color, lw=lw)
    plt.setp(bp['whiskers'], color=color, lw=lw)
    plt.setp(bp['caps'], color=color, lw=lw)
    plt.setp(bp['medians'], color=color, lw=lw)
    plt.setp(bp['fliers'], color=color)

#%% data loop
lon=None

for li, loc in enumerate(['arctic','aa']):
    if li==0: continue
    path1 = root_paths[li]
    
    for year in YEARS:
        
        try: out_data = pickle.load(open(path1+'wind/energy_'+str(year)+'.pkl', 'rb'))
        except FileNotFoundError:
            try:
                out_data = {mm:[] for mm in months}
            
                census_file = path1+'census/'+census_name+str(year)+'.csv'
                [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
                
                for mi, month in enumerate(months):
                    # load monthly era5
                    try:
                        days_in = [LZ(dd+1) for dd in range(calendar.monthrange(year, month)[1])]
                        ds_w1 = get_era5_wind(year, month, days_in, li)
                        dt2 = datetime(year, month, int(days_in[-1]))+timedelta(days=1)
                        ds_w2 = get_era5_wind(dt2.year, dt2.month, [LZ(dd+1) for dd in np.arange(0,7)], li)
                        ds_w = xr.concat([ds_w1, ds_w2], dim="valid_time")
                    except:
                        print('trying era again')
                        days_in = [LZ(dd+1) for dd in range(calendar.monthrange(year, month)[1])]
                        ds_w1 = get_era5_wind(year, month, days_in, li)
                        dt2 = datetime(year, month, int(days_in[-1]))+timedelta(days=1)
                        ds_w2 = get_era5_wind(dt2.year, dt2.month, [LZ(dd+1) for dd in np.arange(0,7)], li)
                        ds_w = xr.concat([ds_w1, ds_w2], dim="valid_time")
                        
                    
                    if lon is None:
                       lon, lat = np.meshgrid(ds_w['longitude'].values, ds_w['latitude'].values)
                    
                    # monthly storm data loop
                    enddt = [ed for i, ed in enumerate(enddate) if startdate[i].month==month]
                    startdt = [sd for i, sd in enumerate(startdate) if sd.month==month]
                    prev_sd = ''
                    for idx, (sd, ed) in enumerate(zip(startdt, enddt)):
                        
                        # get storm area
                        ncname = sd.strftime('%Y_%m%d')+'_contours.nc'
                        try: cs = xr.open_dataset(path1+'contours/'+ncname)
                        except FileNotFoundError:
                            cs_files = glob(xr.open_dataset(path1+'contours/'+ sd.strftime('%Y_%m%d')+'*_contours.nc'))
                            cs_files.sort()
                            if len(cs_files)==1: cfile = cs_files[0]
                            elif len(cs_files)>1 and sd!=prev_sd:cfile = cs_files[0]
                            elif len(cs_files)>1 and sd==prev_sd:cfile = cs_files[-1]
                        all_contours = []
                        for key in list(cs.keys()):
                            coord = cs[key].values
                            all_contours.append(coord)
                        cs.close()
                        
                        with fx.HidePrint(): bbox_edges = fx.get_bbox_edges(all_contours) 
                        inside_points1000 = fx.find_points_in_contour(bbox_edges, lon, lat)
                        
                        # isolate era5 data over storm
                        ds1 = ds_w.sel(valid_time=slice(sd, ed))
                        
                        u_wind = ds1['u10'].values
                        v_wind = ds1['v10'].values
                        wind_tot = np.array( np.sqrt((u_wind**2) + (v_wind**2)) )
    
                        # calculate wind metrics
                        total_winds = []                    
                        for daily_w in wind_tot:
                            wind_in_tot = np.ma.masked_array(daily_w, mask=inside_points1000).filled(np.nan)
                            total_winds.append(wind_in_tot)
                            
                        ace_area = np.nanmean(total_winds)**2
                        ace_max = np.nanmax(total_winds)**2
                        
                        u_in = [np.ma.masked_array(u, mask=inside_points1000).filled(np.nan) for u in u_wind]
                        u_tmean = np.nanmean(u_in, axis=0)**2
                        u_smean = np.nanmean(u_in, axis=(1,2))**2
                        v_in = [np.ma.masked_array(v, mask=inside_points1000).filled(np.nan) for v in v_wind]
                        v_tmean = np.nanmean(v_wind, axis=0)**2
                        v_smean = np.nanmean(v_in, axis=(1,2))**2
                            
                        mke1 = (1/2)*(np.nanmean(u_tmean+v_tmean))
                        mke2 = (1/2)*(np.nanmean(u_smean+v_smean))
                        
                        # append and save?
                        out_data[month].append([ace_area, ace_max, mke1, mke2])
    
                        # end storm loop
                        prev_sd = sd
            except Exception as ee:
                print('error with '+str(year)+' data')
                print(ee)
                print('continuing with next year...')
                
                    
            # save year file
            pickle.dump(out_data, open(path1+'wind/energy_'+str(year)+'.pkl', 'wb'))

#%% end data loop
# raise ValueError('done')

#%% wind tseries analysis

### set up plot
fig, axes = plt.subplots(2,1, sharex=True, sharey=True)

### data loop
for loc_ind, loc in enumerate(['Arctic','Antarctic']):
    path1 = root_paths[loc_ind]
    
    for yi, years in enumerate(decades):

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if loc_ind ==0:
                wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
            elif loc_ind==1: 
                wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        dist = {mm:[] for mm in months}
        for mi, month in enumerate(months):
            wind_lines = {}
            try: 
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[month] if len(arr)==22]  
                    for i, k in enumerate(['u','v','tot']):
                        wind_lines[k] = np.array(windies)[:,:,i]
                elif loc_ind==1: 
                    for i, k in enumerate(['u','v','tot']):
                        windies = [arr[i] for arr in wind_series[month]]  
                        wind_lines[k] = np.array(windies)[:,]
                    
            except IndexError: print('-->', month); continue
            
            if len(wind_lines['tot'])<10: 
                continue
            
            squared = (wind_lines['u']**2) + (wind_lines['v']**2)
            dist[month] = [np.nanmean(L[7:10]) for L in squared]
        
            
        ### plot
        axes[loc_ind].set_title(loc)
        for month in months:
            axes[loc_ind].boxplot(dist[month], positions=[month+(yi*0.33)-0.15])

for ax in axes.flatten():
    ax.set_xlim([0,13])
    ax.set_xticks(months)
    ax.set_xticklabels(months)          
            
            
#%% energy analysis

# [ace_area, ace_max, mke1, mke2]
titles = ['$ACE_{area}$','$ACE_{max}$','$MKE_{1}$','$MKE_{2}$']
idx = 2

units = '$m^2 s^{-2}$'
decade_colors = ['#2166ac', '#b2182b']

### set up plot
fig, axes = plt.subplots(2,1, sharex=True, sharey=True, figsize=(8,6))

### data loop
for li, loc in enumerate(['Arctic','Antarctic']):
    path1 = root_paths[li]
    
    compare = []
    for yi, years in enumerate(decades):
        flierprops = dict(marker='o', markerfacecolor=decade_colors[yi], markersize=4,
                  linestyle='none', markeredgecolor=decade_colors[yi], alpha=0.5)
        
        data = {mm:[] for mm in months}
        
        for year in years:
            
            try:
                out_data = pickle.load(open(path1+'wind/energy_'+str(year)+'.pkl', 'rb'))
                
                for mm in months:
                    data[mm] += [sublist[idx] for sublist in out_data[mm]]
                    
            except: pass
            
        
        ### plot
        axes[li].set_title(loc)
        for month in months:
            bp = axes[li].boxplot(data[month], positions=[month+(yi*0.33)-0.15],
                                  flierprops=flierprops)
            set_box_color(bp, decade_colors[yi], lw=1)
            
        compare.append(data)

    ymin, ymax = axes[li].get_ylim()
    step = 0.05*ymax
    for mm in months:
        cl = check_significance(compare[0][mm], compare[1][mm])
        # print(mm, cl)
        if np.any(cl >= 95):
            axes[li].axvspan(mm-0.25,mm+0.25, color='gray', alpha=0.25, zorder=-1)
            for i, (c, t) in enumerate(zip(cl, ['t', 'KS', 'MW'])):
                if c>=95:
                    axes[li].text(mm-0.5, ymax-(step*(i+1)), t+':'+str(round(c,1)), fontsize=7.5) #, transform=axes[li].transAxes)

for ax in axes.flatten():
    ax.set_xlim([0,13])
    ax.set_xticks(months)
    ax.set_xticklabels(months) 
    ax.set_ylabel('{}'.format(titles[idx])+' [{}]'.format(units))
axes[1].set_xlabel('Month')

#%% xxxxxxxxxxxxxxxxxxx

#%% storm pressure analysis
# spaghetti_pressure1.py
from scipy.stats import ttest_ind, ks_2samp, mannwhitneyu
import string
fontsize=11
hemi_names = ['Arctic','Antarctic']

cl_data_p = {}
cl_data_c = {}
for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    decade_colors = ['#2166ac', '#b2182b']
    
    for mi, month in enumerate(months):
        
        era_data =[]
        count_data = []
        for era, years in enumerate(decades):
            ystr = str(years[0])+'-'+str(years[-1])
        
            # data
            start, end, pressure = fx.storm_pressure(years, path1+'census/', path1+'area/')
        
            p_list = pressure[month]
            if len(p_list)<10: continue
            era_data.append(p_list)
        
            yrly_data = []
            for y in years:
                yrly_data.append([s for s in start[month] if s.year==y])
            count_data.append([len(yx) for yx in yrly_data])
            
        # statistics 
        if len(era_data)==2:
            t, pt = ttest_ind(era_data[0], era_data[1])
            ks, pks = ks_2samp(era_data[0], era_data[1])
            mw, pmw = mannwhitneyu(era_data[0], era_data[1])
            
            pvals = [pt, pmw, pks]
            cl1 = [100*(1-p) for p in pvals]
            cl_data_p[str(loc_ind)+'_'+str(month)] = cl1
            
            
            t, pt = ttest_ind(count_data[0], count_data[1])
            ks, pks = ks_2samp(count_data[0], count_data[1])
            mw, pmw = mannwhitneyu(count_data[0], count_data[1])
            
            pvals = [pt, pmw, pks]
            cl2 = [100*(1-p) for p in pvals]
            cl_data_c[str(loc_ind)+'_'+str(month)] = cl2
            

#%%%% plot mean pressures, noting significance

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
# fig.suptitle('Mean Minimum Storm Pressure')

alph = iter(list(string.ascii_lowercase))
for ax1 in axes:
    ax1.text(0.0, 1.05, '('+next(alph)+')', 
             transform=ax1.transAxes, fontsize=fontsize-2, 
             zorder=50)

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    decade_colors = ['#2166ac', '#b2182b']
    
    ax = axes[loc_ind]
    ax.set_title(loc)
    ax.set_ylabel('Mean Pressure (hPa)')
    ax.set_xticks(np.arange(1,12+1))
    
    for era, years in enumerate(decades):
        ystr = str(years[0])+'-'+str(years[-1])
    
        start, end, pressure = fx.storm_pressure(years, path1+'census/', path1+'area/')
    
        pressure_values = [np.nanmean(pressure[mm]) if len(pressure[mm])>10 else np.nan for mm in months]
        pressure_std = [np.nanstd(pressure[mm]) if len(pressure[mm])>10 else np.nan for mm in months]
        
        ax.plot(months, pressure_values, 
                marker='o', color = decade_colors[era], label = ystr)
        
        ax.errorbar(months, pressure_values, yerr=pressure_std, 
                    color = decade_colors[era], capsize=5)
        
        for mm in months:
            if len(pressure[mm])<10: continue
            if np.any(np.array(cl_data_p[str(loc_ind)+'_'+str(mm)])>=95):
                ax.axvspan(mm-0.25,mm+0.25, color='gray', alpha=0.25, zorder=-1)
         
# legend on last plot
ax.plot([],[], lw=5, color='gray', label='Significant Difference')
ax.legend(loc='upper left', handletextpad=0.5, handlelength=1.25, ncol = 1)
ax.set_xlabel('Month');

#%%% storm counts (per year)
decade_colors = ['#2166ac', '#b2182b']

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
fig.suptitle('Storm Counts per Year')

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    
    ax = axes[loc_ind]
    ax.set_title(loc)
    ax.set_ylabel('Number of Storms')
    ax.set_xticks(np.arange(1,12+1))
    
    for era, years in enumerate(decades):
        ystr = str(years[0])+'-'+str(years[-1])
        
        storm_counts = []
        for year in years:
            start, end, pressure = fx.storm_pressure([year], path1+'census/', path1+'area/')
            storm_counts.append([len(start[mm]) for mm in months])

        mean_counts = np.nanmean(storm_counts, axis=0)
        std_counts = np.nanstd(storm_counts, axis=0)
        
        ax.plot(months, mean_counts, 
                marker='o', color = decade_colors[era], label = ystr)
        
        ax.errorbar(months, mean_counts, yerr=std_counts, 
                    color = decade_colors[era], capsize=5)
        
        for mm in months:
            try:
                if np.any(np.array(cl_data_c[str(loc_ind)+'_'+str(mm)])>=95):
                    ax.axvspan(mm-0.25,mm+0.25, color='gray', alpha=0.25, zorder=-1)
            except KeyError:
                pass

         
# legend on last plot
ax.plot([],[], lw=5, color='gray', label='Significant Difference')
ax.legend(loc='upper left', handletextpad=0.5, handlelength=1.25, ncol = 1)
ax.set_xlabel('Month');


#%% end
