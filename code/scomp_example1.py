#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 2025
scomp_example1.py

- intro to wrappers
- parellel processing

@author: mundi
"""
#%% imports and file paths
import numpy as np
import glob
import xarray as xr
import calendar
import matplotlib.pyplot as plt

from concurrent.futures import ThreadPoolExecutor
import time
from functools import wraps

REMOTE=False

if not REMOTE: ice_fname = '/Users/mundi/Desktop/seaice/'
elif REMOTE: ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/daily/'

#%% functions

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

#%% example wrappers

def some_function():
    print("Ran some_function")
    
def my_decorator(func_to_run):
    def wrapper():
        print("Ran wrapper")
        func_to_run()
        print("Finished wrapper")
    return wrapper

f = my_decorator(some_function)
f()

#### another
@my_decorator
def hello():
    print("Hello")
hello()

#%% data example

files = glob.glob(ice_fname+'*2010*nc')

ds_si = xr.open_dataset(files[0])

print('array shape:', ds_si.sizes)

var = 'cdr_seaice_conc'
seaice = ds_si[var].values

for flag in [251,252,253,254,255]:
    seaice= np.where(seaice==flag/100, np.nan, seaice)
    
si_extent = np.where(seaice<0.15, np.nan, seaice)
    
x, y = ds_si['xgrid'], ds_si['ygrid']

#### lon/lat
from pyproj import Transformer
x,y = np.meshgrid(ds_si['xgrid'].values, ds_si['ygrid'].values )
transformer = Transformer.from_crs("EPSG:3411", "EPSG:4326", always_xy=True)
lon, lat = transformer.transform(x, y)

#### map plot
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo

def setup_plot(extent = [-180,180,60,90]):
    
    fig=plt.figure(figsize=[8,8]) 
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    
    import matplotlib.path as mpath
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    # geography
    ax.coastlines('50m',edgecolor='black',linewidth=0.75)
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    ax.add_feature(cfeature.LAKES, facecolor='0.85')
    return fig, ax

fig, ax = setup_plot()
plotted = ax.pcolormesh(lon, lat, np.squeeze(si_extent), cmap=cmo.ice,
              transform=ccrs.PlateCarree())

cax1 = fig.add_axes([0.33,0.033,0.4,0.033]) 
cbar1 = fig.colorbar(plotted, cax=cax1, orientation='horizontal')
cbar1.set_label(r'Sea Ice Concentration', fontsize=18)
cax1.tick_params(labelsize=18)
    
#%%
# # @timethis 
# def run_threaded_timeseries(analysis_range):
#     with ThreadPoolExecutor(max_workers=10) as executor:
#         return executor.map(calc_winds, analysis_range)
    
# tseries = run_threaded_timeseries(analysis_range)
# return [x for x in tseries]


def get_daily_extent(date_str):

    file = glob.glob(ice_fname + '*' + date_str + '*.nc')[0]
    
    ds_si = xr.open_dataset(file)
    
    seaice = ds_si['cdr_seaice_conc'].values

    for flag in [251,252,253,254,255]:
        seaice= np.where(seaice==flag/100, np.nan, seaice)

    si_extent = np.where(seaice<0.15, np.nan, seaice)
    
    return np.nansum(si_extent*25*25) # multiply by grid area


#%% slow run

@timethis
def get_timeseries_slow(years):
    timeseries = []
    for year in years:
        for month in np.arange(1,12+1):
            for day in np.arange(1,calendar.monthrange(year, month)[-1]+1):
                
                date_str = str(year)
                date_str += '0'+str(int(month)) if month<10 else str(month)
                date_str += '0'+str(int(day)) if day<10 else str(day)
                
                daily_extent = get_daily_extent(date_str)
                   
                timeseries.append(daily_extent)
                
    return timeseries

timeseries1 = get_timeseries_slow([2010,2011])

plt.plot(timeseries1)

#%% fast run

@timethis 
def get_timeseries_threaded(date_list):
    with ThreadPoolExecutor(max_workers=10) as executor:
        return executor.map(get_daily_extent, date_list)
    
@timethis
def get_date_list(years):
    date_list = []
    for year in years:
        for month in np.arange(1,12+1):
            for day in np.arange(1,calendar.monthrange(year, month)[-1]+1):
                
                date_str = str(year)
                date_str += '0'+str(int(month)) if month<10 else str(month)
                date_str += '0'+str(int(day)) if day<10 else str(day)
                date_list.append(date_str)
    return date_list
    
date_list = get_date_list([2010,2011])
tseries = get_timeseries_threaded(date_list)
timeseries2 = [x for x in tseries]


#%% end
