#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  8 2025
functions.py

@author: mundi
"""
#%% imports
import numpy as np
import xarray as xr
import netCDF4
from datetime import datetime, timedelta
from glob import glob

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import sys, os
import warnings

#%% COLORS

def month_colors():
    colors = []
    return colors

def hemi_colors():
    colors = []
    return colors 

def decade_colors():
    colors = []
    return colors


#%% data

def daterange(start_date, end_date, dt=6):
    alldates=[]
    delta = timedelta(hours=dt)
    while start_date <= end_date:
        alldates.append(start_date)
        start_date += delta
    return alldates

def readCensus(file, convertDT=False):
    import csv
    
    # Load file
    csv_file = open(file,'r')
    startdate, enddate = [],[]
    startlon, endlon = [], []
    startlat, endlat = [],[]
    pressure = []
    
    # Read off and discard first line, to skip headers
    csv_file.readline()
    
    # Split columns while reading
    for a, b, c, d, e, f, g in csv.reader(csv_file, delimiter=','):
        # Append each variable to a separate list
        startdate.append(a) 
        startlat.append(float(b))
        startlon.append(float(c))
        pressure.append(float(d))
        enddate.append(e)
        endlat.append(float(f))
        endlon.append(float(g))
    csv_file.close()
    
    if convertDT:
        startDT, endDT = [],[]                            
        for i, pres in enumerate(pressure):
            startDT.append( datetime(int(startdate[i][:4]), int(startdate[i][5:7]), 
                                           int(startdate[i][8:10]), 0) )
                           # int(startdate[i][-2:]))
            endDT.append( datetime(int(enddate[i][:4]), int(enddate[i][5:7]), 
                                           int(enddate[i][8:10]), 0) )
                            # int(enddate[i][-2:])
                            
        # startdate, enddate = startDT, endDT
        
        ### sort times in ascending order
        startdate, enddate, startlon, startlat, endlon, endlat, pressure = \
            map(list, zip(*sorted(zip(startDT, endDT, startlon, startlat, endlon, endlat, pressure))))

    
    dates = [startdate, enddate]
    coords = [[startlon, startlat],[endlon,endlat]]
    return dates, coords, pressure 

#%% plotting

def legend_without_duplicate_labels(ax, loc='best',fontsize=20):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    
    return ax.legend(*zip(*unique), loc=loc, fontsize=fontsize)

#%%% maps
def background_plot_sh(extent=[-180,180, -53,-90], returnfig=False, title=[], labels=True):
    
    fig=plt.figure(figsize=[15,15]) 
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax.coastlines('50m',edgecolor='black',linewidth=0.75)
    ax.set_extent(extent, ccrs.PlateCarree())
    try:
        ax.gridlines(draw_labels=labels)
    except:
        print('Unable to create grid lines on map')
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    ax.add_feature(cfeature.LAKES, facecolor='0.85')
    
    if title:
        if type(title)!= str: title=str(title)
        ax.set_title(title)
    
    if not returnfig: return ax
    if returnfig: return fig, ax
    
def setup_plot_sh(ax, extent=[-180,180, -53,-90], title=[], labels=True):
    ax.coastlines('50m',edgecolor='black',linewidth=0.75)
    ax.set_extent(extent, ccrs.PlateCarree())
    try:
        ax.gridlines(draw_labels=labels)
    except:
        print('Unable to create grid lines on map')
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    ax.add_feature(cfeature.LAKES, facecolor='0.85')
    if title:
        if type(title)!= str: title=str(title)
        ax.set_title(title)
    return ax


def background_plot_nh(extent=[-160,90,50,60], returnfig=False, title=[], labels=True, central_lon=-45):
    # alt extent: [-50,90,60,85]
    
    fig=plt.figure(figsize=[15,15]) 
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=central_lon))
    ax.coastlines('50m',edgecolor='black',linewidth=0.75)
    ax.set_extent(extent, ccrs.PlateCarree())
    try:
        ax.gridlines(draw_labels=labels)
    except:
        print('Unable to create grid lines on map')
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    ax.add_feature(cfeature.LAKES, facecolor='0.85')
    
    if title:
        if type(title)!= str: title=str(title)
        ax.set_title(title, fontsize=22)
    
    if not returnfig: return ax
    if returnfig: return fig, ax


def setup_plot_nh(ax, extent=[-160,90,50,60], title=[], labels=True):
    ax.coastlines('50m',edgecolor='black',linewidth=0.75)
    ax.set_extent(extent, ccrs.PlateCarree())
    try:
        ax.gridlines(draw_labels=labels)
    except:
        print('Unable to create grid lines on map')
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    ax.add_feature(cfeature.LAKES, facecolor='0.85')
    if title:
        if type(title)!= str: title=str(title)
        ax.set_title(title)
    return ax

def geoplot_2d(x,y,z=None):    
    #do masked-array on the lon
    x_greater = np.ma.masked_greater(x, -0.01)
    x_lesser = np.ma.masked_less(x, 0)    
    # apply masks to other associate arrays: lat
    y_greater = np.ma.MaskedArray(y, mask=x_greater.mask)
    y_lesser = np.ma.MaskedArray(y, mask=x_lesser.mask)
    
    if z is None:
        return [x_greater, x_lesser], [y_greater, y_lesser]
    else:
        # apply masks to other associate arrays
        z_greater = np.ma.MaskedArray(z, mask=x_greater.mask)
        z_lesser = np.ma.MaskedArray(z, mask=x_lesser.mask)
        return [x_greater, x_lesser], [y_greater, y_lesser], [z_greater, z_lesser]
    
def plot_geocontour(ax, lon, lat, var, levels, color='k', lw=3, ls='solid'):
    #do masked-array on the lon
    lon_greater = np.ma.masked_greater(lon, -0.01)
    lon_lesser = np.ma.masked_less(lon, 0)    
    # apply masks to other associate arrays: lat
    lat_greater = np.ma.MaskedArray(lat, mask=lon_greater.mask)
    lat_lesser = np.ma.MaskedArray(lat, mask=lon_lesser.mask)
    # apply masks to other associate arrays: daily ice
    si_greater = np.ma.MaskedArray(var, mask=lon_greater.mask)
    si_lesser = np.ma.MaskedArray(var, mask=lon_lesser.mask)

    # contours
    ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
              linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
              linestyles=ls) 
    ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
              linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
              linestyles=ls)
    return ax

#%%% winds

def plot_winds(x,y,u,v, ax, mask=[], sc=300, interval=12):  

    u = np.squeeze(u)
    v = np.squeeze(v)    

    if len(mask)>0:
        u=np.ma.masked_array(u, mask=mask)
        v=np.ma.masked_array(v, mask=mask)
    
    #do masked-array on the lon
    x_greater = np.ma.masked_greater(x, -0.01)
    x_lesser = np.ma.masked_less(x, 0)    
    # apply masks to other associate arrays: lat
    y_greater = np.ma.MaskedArray(y, mask=x_greater.mask)
    y_lesser = np.ma.MaskedArray(y, mask=x_lesser.mask)
    # apply masks to other associate arrays: u,v
    u_greater = np.ma.MaskedArray(u, mask=x_greater.mask)
    u_lesser = np.ma.MaskedArray(u, mask=x_lesser.mask)
    v_greater = np.ma.MaskedArray(v, mask=x_greater.mask)
    v_lesser = np.ma.MaskedArray(v, mask=x_lesser.mask)
    
    
    ax.quiver(x_greater[::interval, ::interval], y_greater[::interval, ::interval],
            u_greater[::interval, ::interval], v_greater[::interval, ::interval],
            color='k', transform = ccrs.PlateCarree(), scale = sc, zorder = 6,
            headwidth=2.5, headlength=3)
    Q = ax.quiver(x_lesser[::interval, ::interval], y_lesser[::interval, ::interval],
            u_lesser[::interval, ::interval], v_lesser[::interval, ::interval],
            color='k', transform = ccrs.PlateCarree(), scale = sc, zorder = 6,
            headwidth=2.5, headlength=3)
    
    return ax, Q
    
#%%% lines

def plot_spread(ax, x, mean_line, std_line, color='C0', ls1='-', ls2='--'):
    
    ax.plot(x, mean_line, lw=2, 
            color=color, ls=ls1)
    
    ax.plot(x, mean_line+std_line, lw=0.5, color=color, ls=ls2)
    ax.plot(x, mean_line-std_line, lw=0.5, color=color, ls=ls2)
    
    ax.fill_between(x, mean_line, mean_line+std_line, 
                     alpha=0.275, color=color, zorder=-5)
    ax.fill_between(x, mean_line,mean_line-std_line, 
                    alpha=0.275, color=color, zorder=-5)
    
    return ax

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

#%% SEA ICE

def load_netcdf(filepath, in_vars):
    """open netcdf file, load variables from list in_vars and output dictionary of variables"""

    out_vars = {}

    open_netcdf = netCDF4.Dataset(filepath, mode = 'r')
    #print open_netcdf
    for var in in_vars:
        out_vars[var] = open_netcdf.variables[var][:]
    open_netcdf.close()

    return out_vars

def load_seaice_sh(root_dir, year, month, day, latlon=True):
    from pyproj import Transformer
    
    # convert date inputs to strings for filename if needed
    if not isinstance(year, str):
        year = str(year)
    if not isinstance(month, str):
        if month<10: month = '0'+str(int(month))
        else: month = str(int(month))
    if not isinstance(day, str):
        if day<10: day = '0'+str(int(day))
        else: day = str(int(day))
    
    # get file(s)
    all_files = glob(root_dir + '*' + year+month+day + '*.nc')
    all_files.sort()

    if not all_files:
        print('Error with S.Hemisphere filename: ' + root_dir + '*' + year+month+day + '*.nc')
        raise NameError(' bad filename in load_seaice_sh')
        
    cdr_dic = xr.open_dataset(all_files[0])
    seaice = np.squeeze( cdr_dic['nsidc_bt_seaice_conc'].values )
    
    for flag in [251,252,253,254,255]:
        seaice= np.where(seaice==flag/100, np.nan, seaice)
   
    if latlon:
        x,y = np.meshgrid( cdr_dic['xgrid'].values, cdr_dic['ygrid'].values )
        transformer = Transformer.from_crs("EPSG:3412", "EPSG:4326", always_xy=True)
        lon, lat = transformer.transform(x, y)
        return seaice, lon, lat
    else:
        return seaice

def load_seaice_v4(root_dir, year, month, day, latlon=True):
    seaice_daily, xgrid, ygrid = [],[],[]
    
    # convert date inputs to strings for filename if needed
    if not isinstance(year, str):
        year = str(year)
    if not isinstance(month, str):
        if month<10: month = '0'+str(int(month))
        else: month = str(int(month))
    if not isinstance(day, str):
        if day<10: day = '0'+str(int(day))
        else: day = str(int(day))
    
    # get file(s)
    all_files = glob(root_dir + '*' + year+month+day + '*.nc')
    all_files.sort()
    variable_names = ['nsidc_bt_seaice_conc','xgrid','ygrid']

    if not all_files:
        print('Error with V4 filename: ' + root_dir + '*' + year+month+day + '*.nc')
        raise NameError(' bad filename in load_seaice_v4')
        
    for n, filename in enumerate(all_files):
        try:
            cdr_dic = load_netcdf(filename, variable_names)
            seaice = cdr_dic['nsidc_bt_seaice_conc']
        except KeyError:
            print('V3 sea ice data used !')
            variable_names = ['goddard_bt_seaice_conc','xgrid','ygrid']
            cdr_dic = load_netcdf(filename, variable_names)
            seaice = cdr_dic['goddard_bt_seaice_conc']
        
        if latlon:
            ygrid    = cdr_dic['ygrid']
            xgrid    = cdr_dic['xgrid']
        else:
            ygrid=[];xgrid=[]
       
        seaice_daily = seaice.mean(axis=0)
    
    return seaice_daily, xgrid, ygrid


def load_seaice(root_dir, year, month, day, latlon=True):
    
    seaice_daily, x, y = load_seaice_v4(root_dir, year, month, day, latlon=True)
    
    if not latlon: return seaice_daily
    
    ds = xr.open_dataset(root_dir + 'seaice_lonlat_v03.nc', decode_times=False)
    lon = ds['longitude'].values
    lat = ds['latitude'].values
    
    return seaice_daily, lon, lat

#%% LINES
ice_lims = [20,80]
var_addons = [['_1000', '_990', '_contour'],['2_1000', '2_990', '2_contour']]
addon_num = 0 # storm area
sia_var_addons = var_addons
miz_ind = 1 #[0=daily, 1=total]
miz_names = ['daily_', 'total_']

def indiv_locations(years, census_path, area_path):
    
    start_day = {}
    end_day = {}
    latitudes = {}
    longitudes = {}
    for idx in np.arange(0,12):
        start_day[idx+1] = []
        end_day[idx+1] = []
        latitudes[idx+1] = []
        longitudes[idx+1] = []
        
    for year in years:
        
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate], [[startlon, startlat],[endlon,endlat]] = readCensus(census_file, convertDT=True)[0:2]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in timing_grid:
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(daterange(week_ago, two_week, dt=24))
            storm_ranges.append(daterange(startdt, enddt, dt=24))  
            
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
             latitudes[month].append(np.nanmean([startlat[storm_num], endlat[storm_num]]))
             longitudes[month].append(np.nanmean([startlon[storm_num], endlon[storm_num]]))
            
    return start_day, end_day, longitudes, latitudes
        

def indiv_lines(years, census_path, area_path, si_path, sample_size=10):
    storm_counts = {}
    lines = {}
    start_day = {}
    end_day = {}
    
    si_changes = {} 
    clim_changes = {}
    
    latitudes = {}
    
    months = np.arange(1,12+1)
    si_name = '_seaice'
    census_name = 'census_'
    
    for idx in np.arange(0,12):
        lines[idx+1] = []
        start_day[idx+1] = []
        end_day[idx+1] = []
        si_changes[idx+1] = []
        clim_changes[idx+1] = []
        latitudes[idx+1] = []
    
    for year in years:
        storm_counts[year] = [0]*len(months)
        
        census_file = census_path+census_name+str(year)+'.csv'
        [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in timing_grid:
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(daterange(week_ago, two_week, dt=24))
            storm_ranges.append(daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
        try:
            ds = xr.open_dataset(si_path + str(year) + si_name + '.nc')
        except:
            print('- skip: '+si_path + str(year) + si_name + '.nc')
            continue
        
        for storm_num, strm in enumerate(storm_ranges):
             month = int(strm[0].month)
             
             ### remove storms that don't interact with the ice
             with warnings.catch_warnings():
                 warnings.simplefilter('ignore')
                 ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
             if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                continue
             try:
                 sia = ds['sia_miz'+sia_var_addons[miz_ind][addon_num]].values[storm_num]
             except (KeyError, IndexError):
                 print()
                 print('***', year, len(storm_ranges), 
                       np.shape(ds['sia_miz'+sia_var_addons[miz_ind][addon_num]].values))
                 print()
                 break
        
             storm_counts[year][month-1] += 1 # mm -> index value
             
             # get time series   
             try:
                 sia_clim = ds['sia_clim_miz'+sia_var_addons[miz_ind][addon_num]].values[storm_num]
             except KeyError:
                 sia_clim = ds['si_clim'].values[storm_num]
             ss = sia-sia_clim
             with warnings.catch_warnings():
                 warnings.simplefilter('ignore')
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
        
        if len(lines[idx+1])> sample_size: 
            mean_line = np.nanmean(lines[idx+1], axis=0)
            mean_lines.append(mean_line)
        else:
            mean_lines.append(np.nan*np.ones(np.shape(analysis_ranges[0])))
            continue
        
    return mean_lines, lines, start_day, end_day, si_changes, clim_changes

def wind_lines(years, census_path, area_path, wind_path, NAME = 'winds_'):
    storm_counts = {}
    start_day = {}
    end_day = {}
    
    winds = {}
    
    months = np.arange(1,12+1)
    
    for idx in np.arange(0,12):
        start_day[idx+1] = []
        end_day[idx+1] = []
        winds[idx+1] = []
        
    for year in years:
        storm_counts[year] = [0]*len(months)
        
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in timing_grid:
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(daterange(week_ago, two_week, dt=24))
            storm_ranges.append(daterange(startdt, enddt, dt=24))  
        
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
        
             storm_counts[year][month-1] += 1 # mm -> index value
             
             # get time series   
             try:
                 wind_in = np.load(wind_path+NAME+str(year)+'-'+str(month)+'.npy', allow_pickle=True)
             except FileNotFoundError:
                 print('winds_'+str(year)+'-'+str(month)+'.npy')
                 continue
        
             # plot
             start_day[month].append(strm[0])
             end_day[month].append(strm[-1])
             
             winds[month].append(wind_in[storm_counts[year][month-1]-1])
        
    return start_day, end_day, winds

def sh_winds(years, census_path, area_path, wind_fpath):
    ice_lims = [20,80]
    dec_data = {mm: [] for mm in np.arange(1,12+1)}
    for year in years:
        # WIND FILE
        wind_ds = xr.open_dataset(wind_fpath+str(year)+'_wind.nc')
        
        # CENSUS: storm info
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in zip(startdate, enddate):
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(daterange(week_ago, two_week, dt=24))
            storm_ranges.append(daterange(startdt, enddt, dt=24))  
        
        # ICE AREA/SEA ICE
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
    
        # STORM LOOP
        for storm_num, strm in enumerate(storm_ranges):
            month = int(strm[0].month)
             
            # remove storms that don't interact with the ice
            ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
            if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
               continue
   
            # wind analyis
            u_wind = wind_ds['total_u_1000'][storm_num].values
            v_wind = wind_ds['total_v_1000'][storm_num].values
            t_wind = wind_ds['total_winds_1000'][storm_num].values
             
            dec_data[month].append([t_wind, u_wind, v_wind])

    return dec_data

def era_lines(years, census_path, area_path, era_path, var):
    storm_counts = {}
    start_day = {}
    end_day = {}
    
    swhs = {}
    
    months = np.arange(1,12+1)
    
    for idx in np.arange(0,12):
        start_day[idx+1] = []
        end_day[idx+1] = []
        swhs[idx+1] = []
        
    for year in years:
        storm_counts[year] = [0]*len(months)
        
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in timing_grid:
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(daterange(week_ago, two_week, dt=24))
            storm_ranges.append(daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
        for storm_num, strm in enumerate(storm_ranges):
             month = int(strm[0].month)
             
             ### remove storms that don't interact with the ice
             with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
             if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                continue
        
             storm_counts[year][month-1] += 1 # mm -> index value
             
             # get time series   
             try:
                 swh_in = np.load(era_path+var+'_'+str(year)+'-'+str(month)+'.npy', allow_pickle=True)
             except FileNotFoundError:
                 print(var+'_'+str(year)+'-'+str(month)+'.npy')
                 continue
        
             # plot
             start_day[month].append(strm[0])
             end_day[month].append(strm[-1])
             
             swhs[month].append(swh_in[storm_counts[year][month-1]-1])
        
    return start_day, end_day, swhs

def storm_pressure(years, census_path, area_path):
    storm_counts = {}
    start_day = {}
    end_day = {}
    
    pressure = {}
    
    months = np.arange(1,12+1)
    
    for idx in np.arange(0,12):
        start_day[idx+1] = []
        end_day[idx+1] = []
        pressure[idx+1] = []
        
    for year in years:
        storm_counts[year] = [0]*len(months)
        
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
        slp = readCensus(census_file, convertDT=True)[-1]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in timing_grid:
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(daterange(week_ago, two_week, dt=24))
            storm_ranges.append(daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
        for storm_num, strm in enumerate(storm_ranges):
             month = int(strm[0].month)
             
             ### remove storms that don't interact with the ice
             with warnings.catch_warnings():
                 warnings.simplefilter('ignore')
                 ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
                 if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                    continue
        
             storm_counts[year][month-1] += 1 # mm -> index value
             
             start_day[month].append(strm[0])
             end_day[month].append(strm[-1])
             
             pressure[month].append(slp[storm_num])
        
    return start_day, end_day, pressure

def ocn_profiles(years, path1, all_or_miz):
    if type(years)==int: years=[years]

    ### LOAD DEPTH
    with xr.open_dataset(path1+'ocn_prof/'+'depth.nc') as dds:
        DEPTH = dds['depth'].values
        
    data = {mm: [] for mm in np.arange(1,12+1)}
    
    ### STORM LOOP
    for year in years:
        # CENSUS: storm info
        census_file = path1+'census/'+'census_'+str(year)+'.csv'
        [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
        
        # open ice area
        ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
        ice_area80 = ds_area['ice_area80'].values
        # ice_area15 = ds_area['ice_area15'].values
        box_area = ds_area['box_area'].values
        ice_sorter = ice_area80
        
        for sn, start in enumerate(startdate): 
            month = start.month
            
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                ice_frac = ice_sorter[sn]*100/box_area[sn]             
            if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                continue
    
            ## LOAD FILES
            start_str = start.strftime('%Y-%m-%d')
            
            try:
                if all_or_miz=='all': fname = '_all.npy'
                elif all_or_miz=='miz': fname = '_miz.npy'
                elif all_or_miz=='series': fname = '_miz_series.npy'
                else:raise NameError('functions.ocn_profiles: incorrect data type (all or miz)')
                
                prof = np.load(path1+'ocn_prof/'+start_str+'_'+str(sn)+fname)
                data[month].append(prof)
                
            except FileNotFoundError:
                print('missing:', start_str)
                continue
            
    return data, DEPTH

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
        [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
        
        storm_ranges = []
        for startdt, enddt in zip(startdate, enddate):
            storm_ranges.append(daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
        for storm_num, strm in enumerate(storm_ranges):
             month = int(strm[0].month)
             
             ### remove storms that don't interact with the ice
             with warnings.catch_warnings():
                 warnings.simplefilter('ignore')
                 ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
             if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                continue
        
             start_day[month].append(strm[0])
             end_day[month].append(strm[-1])
             
             motion[month].append(motion_series[storm_num])
        
    return start_day, end_day, motion


def get_sst_lines(yi, sst_path):
    return {mm:np.load(sst_path+'tseries_'+str(yi)+'-'+str(mm)+'.npy') for mm in np.arange(1,12+1)}
        

#%%%% sea ice
def get_storm_areas(year, root_path, lon, lat):
    ice_lims = [20,80]
    
    census_file = root_path+'census/'+'census_'+str(year)+'.csv'
    [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
    
    timing_grid = []
    for xx in range(0,len(startdate)):
        timing_grid.append((startdate[xx], enddate[xx]))
    storm_ranges = []
    for startdt, enddt in timing_grid:
        storm_ranges.append(daterange(startdt, enddt, dt=24))  
    
    # open ice area
    ds_area = xr.open_dataset(root_path+'area/' + str(year) +'_area.nc')
    ice_sorter = ds_area['ice_area80'].values
    box_area = ds_area['box_area'].values
    
    try:
        ds = xr.open_dataset(root_path+'seaice/' + str(year) + '_seaice.nc')
        ds.close()
    except:
        print('- skip: '+root_path+'seaice/' + str(year) + '_seaice.nc')
        return [], []
    
    storm_areas = []
    stormstr_prev=''
    for storm_num, strm in enumerate(storm_ranges):
        # remove storms that don't interact with the ice
        ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
        if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
            storm_areas.append([])
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
        ncname = stormstr + '_contours.nc'
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
            storm_areas.append([])
            continue
        
        with HidePrint(): bbox_edges = get_bbox_edges(all_contours) 
        storm_areas.append( find_points_in_contour(bbox_edges, lon, lat) )
        
    return storm_areas, storm_ranges
        
def get_sic_grad(loc, years, grad_path, root_path):
    
    gradients = {mm: [] for mm in np.arange(1,12+1)}
    
    for month in np.arange(1,12+1):
        for year in years:
            savename = loc+'_sic_grad_all_'+str(year)+'_'+str(month)+'.npy'
            miz_series = np.load(grad_path+savename)
            miz_series = iter(miz_series)
            
            # load census and sea ice info
            census_name = 'census_'+ str(year) +'.csv'
            [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
                readCensus(root_path+'census/'+census_name, convertDT=True)
                
            # open ice area
            with xr.open_dataset(root_path+'area/'+str(year)+'_area.nc') as ds_area:
                ice_area80 = ds_area['ice_area80'].values
                # ice_area15 = ds_area['ice_area15'].values
                box_area = ds_area['box_area'].values
                ice_sorter = ice_area80
            
            #### loop thru storms
            for sn, start in enumerate(startdate):
                if start.month != month: continue
                val= next(miz_series)
                
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    ice_frac = ice_sorter[sn]*100/box_area[sn]             
                if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                    continue
                
                gradients[month].append(val)
            
    return gradients

def get_si_conc(loc_ind, years, grad_path, root_path, ice_fname):
    import pickle
    
    savepath = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/monthly_miz/'
    if loc_ind==0: _, si_lon, si_lat = load_seaice(ice_fname, 1, 1, 1, latlon=True)
    elif loc_ind==1: _, si_lon_sh, si_lat_sh = load_seaice_sh(ice_fname, 1, 1, 1, latlon=True)

    miz_month = {mm:[] for mm in np.arange(1,12+1)}
    for year in years:
        
        try:
            with open(root_path+'seaice/'+str(year)+'sic.pkl', 'rb') as f:
                sic_in = pickle.load(f)
                
        except FileNotFoundError:
            
            sic_in = {mm:[] for mm in np.arange(1,12+1)}
            
            ### load storm census
            if loc_ind==0: lon1, lat1 = si_lon, si_lat 
            elif loc_ind==1: lon1, lat1 = si_lon_sh, si_lat_sh
            storm_areas, storm_ranges = get_storm_areas(year, root_path, lon1, lat1)
            
            ### daily miz
            if loc_ind == 0:
                daily_map = np.load(savepath+str(year)+'_map.npy') 
            elif loc_ind == 1:
                daily_map = np.load(savepath+str(year)+'_map_sh.npy')   
                  
            for area, storm in zip(storm_areas, storm_ranges):
                if np.shape(area) !=  np.shape(lon1): continue
            
                miz_list = []
                for date in daterange(storm[0]-timedelta(days=7), storm[0]+timedelta(days=14), dt=24):
                    # get day number to index daily_map
                    dt = datetime(2010, date.month, date.day if (date.day!=29 and date.month!=2) else 28 ) # avoid leap stuff
                    daynum = dt.timetuple().tm_yday
                    sic_map = np.ma.masked_array( daily_map[daynum-1], mask=area)
                    miz_list.append( np.nanmean(sic_map) )
                sic_in[storm[0].month].append(miz_list)
                
            with open(root_path+'seaice/'+str(year)+'sic.pkl', 'wb') as f:
                pickle.dump(sic_in, f, pickle.HIGHEST_PROTOCOL)
            
        for month in np.arange(1,13):
            miz_month[month] += sic_in[month]
            
    return miz_month

#%%% storm areas
def get_storm_bbox(years, census_path, area_path, contour_path):
    
    storm_areas = {}
    for idx in np.arange(0,12):
        storm_areas[idx+1] = []
    
    for year in years:
        ### FULL YEAR
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate] = readCensus(census_file, convertDT=True)[0]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        for startdt, enddt in timing_grid:
            storm_ranges.append(daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        ds_area.close()
        
        ### FY
        stormstr_prev =  ''
        for storm_num, strm in enumerate(storm_ranges):
            month = int(strm[0].month)
             
            ### remove storms that don't interact with the ice
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
            if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
               continue

            stormstr1 = strm[0].strftime('%Y_%m%d')
            # duplicate storm start date?
            if stormstr1==stormstr_prev:
                # stormstr = stormstr1 + next(dupe)
                files = glob(contour_path+stormstr1+'*'+'_contours.nc')
                for f in files:
                    datestr = f.split('_')[-2]
                    if len(datestr)>4:
                        cs = xr.open_dataset(f)
                        break
            else:
                ncname = stormstr1+'_contours.nc'
                try:
                    cs = xr.open_dataset(contour_path+ncname)
                except FileNotFoundError:
                    print('-- '+ncname)
                    continue
            stormstr_prev = stormstr1
            
            all_contours = []
            for key in list(cs.keys()):
                coord = cs[key].values
                all_contours.append(coord)
            cs.close()
            
            ### get bbox
            with HidePrint(): bbox_edges = get_bbox_edges(all_contours) 
            storm_areas[month].append(bbox_edges)
            
    return storm_areas


#%% contours / bbox
class HidePrint:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def find_points_in_contour(coords, var_x, var_y, var=None):
    import matplotlib.path as mpltPath
    
    coord_x = np.where(coords[:,0] < 0, coords[:,0]+360, coords[:,0])
    coord_y = coords[:,1]
    
    ### turn contour into polygon
    polygon = np.vstack((coord_x, coord_y)).T
    path = mpltPath.Path(polygon)
    
    points = np.vstack((var_x.flatten(), var_y.flatten()))
    points=points.T
    
    ### get inside points and plot
    inside = path.contains_points(points)
    inside = np.reshape(inside, np.shape(var_x))
    
    if np.nanmin(var_x) < 0:
        var_x2 = np.where(var_x < 0, var_x+360, var_x)
    else:
        var_x2 = np.where(var_x > 180, var_x-360, var_x)
    points2 = np.vstack((var_x2.flatten(), var_y.flatten())).T
    inside2  = path.contains_points(points2)
    inside2 = np.reshape(inside2, np.shape(var_x))
    
    inside_points = np.invert(np.logical_or(inside, inside2))

    return inside_points

def get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False):
    ### create new bbox edges
    #
    edge1x = np.linspace(minlon, minlon, n)
    edge1y = np.linspace(minlat, maxlat, n)
    edge2x = np.linspace(minlon, maxlon, n)
    edge2y = np.linspace(maxlat, maxlat, n)
    edge3x = np.linspace(maxlon, maxlon, n)
    edge3y = np.linspace(maxlat, minlat, n)
    edge4x = np.linspace(maxlon, minlon, n)
    edge4y = np.linspace(minlat, minlat, n)
    if reverse:
        edge2x = np.concatenate((np.linspace(minlon,0,round(n/3)),np.linspace(0,180,round(n/3))))
        edge2x = np.concatenate( ( edge2x,np.linspace(-180,maxlon,round(n/3)) ) ) 
        #
        edge4x = np.concatenate((np.linspace(maxlon,-180,round(n/3)),np.linspace(180,0,round(n/3))))
        edge4x = np.concatenate( ( edge4x,np.linspace(0,minlon,round(n/3)) ) ) 
    #
    bbox_lon = np.append(edge1x, edge2x)
    bbox_lon = np.append(bbox_lon, edge3x)
    bbox_lon = np.append(bbox_lon, edge4x)
    #
    bbox_lat = np.append(edge1y, edge2y)
    bbox_lat = np.append(bbox_lat, edge3y)
    bbox_lat = np.append(bbox_lat, edge4y)
    #
    bbox_edges = np.squeeze(np.array([[bbox_lon],[bbox_lat]])).T
    
    return bbox_edges

def get_bbox_edges(all_contours):
    minlon, minlat, maxlon, maxlat = 999,999,-999,-999

    isDivided = False
    alerted=False
    for cidx, contour in enumerate(all_contours):
        lons = contour[:,0]
        lats = contour[:,1]
        
        ### get rid of boundaries too close to pole
        for li, lat in enumerate(lats):
            if lat>82.5: ###!!! new maxlat thresh -- 85/82.5
                lons[li]=np.nan
                lats[li]=np.nan
                if not alerted:
                    print('bbox too close too pole; artifical boundary applied (85N)')
                    alerted=True

        ### convert longitude to 0-360 system
        lons1 = lons.copy()
        lons1 = np.where(lons1<0, lons1+360, lons1)
        lons1.sort()
        lons1 = lons1[~np.isnan(lons1)]
        
        if len(lons1) == 0: 
            print('skip')
            continue
        
        ### find e/w lons
        if np.nanmax(lons1) - np.nanmin(lons1) > 180:
            for li, ll in enumerate(lons1):
                if li == len(lons1)-1: break
                if lons1[li+1] - ll > 20: ###!!! threshold here
                    if not isDivided:
                        eastlon = ll
                        westlon = lons1[li+1]
                        isDivided = True
                        break
                    else:
                        if ll > eastlon:
                            eastlon = ll
                        if lons1[li+1] < westlon:
                            westlon = lons1[li+1]
                        break
        else:
            if lons1[0] < minlon:
                minlon = lons1[0]
                print('minlon ', minlon)
            if lons1[-1] > maxlon:
                maxlon = lons1[-1]
                print('maxlon ', maxlon)
        
        ### get min/max lat
        if np.nanmin(lats) < minlat:
            minlat = np.nanmin(lats)
        if np.nanmax(lats) > maxlat:
            maxlat = np.nanmax(lats)
     
        # end contour loop
        
    if not isDivided:    
        print('done easy: ', [minlon, maxlon, minlat, maxlat])
        return get_edge_lines(minlon, maxlon, minlat, maxlat, reverse=False)
    else: # isDivided
        if minlon == 999:
            print('only divide: ', [westlon, eastlon, minlat, maxlat])
            bbox_edges = get_edge_lines(westlon, eastlon, minlat, maxlat, reverse=True)
            bbox_edges = np.where(bbox_edges>180, bbox_edges-360, bbox_edges)
            return bbox_edges
        else:
            print('combo!', [minlon, maxlon], [westlon, eastlon])
            bbox_edges = get_edge_lines(westlon, maxlon, minlat, maxlat, reverse=True)
            if (eastlon<maxlon) and (westlon<minlon):
                bbox_edges = get_edge_lines(westlon,eastlon, minlat, maxlat, reverse=True)
            bbox_edges = np.where(bbox_edges>180, bbox_edges-360, bbox_edges)
            return bbox_edges


#%% end
