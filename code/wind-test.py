#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 2025

@author: mundi
"""

#%% imports, file names
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import time as timeIN
import calendar
import functions as fx
from datetime import datetime, timedelta
from scipy.interpolate import griddata
import warnings

REMOTE=True

if not REMOTE: 
    ROOT = '/Users/mundi/Desktop/'
    ice_path = ROOT +'seaice/'
    ice_paths = [ice_path, ice_path+'south/']
else: 
    ROOT = '/home/mundi/'
    ice_root = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/'
    ice_paths = [ice_root+ice_add for ice_add in ['daily/','south/']]

nh_path = ROOT +'month-hemi/nh_data/'
sh_path = ROOT +'month-hemi/sh_data/'
root_paths = [nh_path, sh_path]


savepath = ROOT+'month-hemi/linear_winds/'

census_name = 'census_'
contour_name = '_contours.nc'
si_name = '_seaice'

#%%% functions

def get_era5_wind(year, months, days, LOC):
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
    if type(months) != list and type(months[0])!=str:
        months = ['0'+str(int(month)) if month<10 else str(int(month)) for month in months]
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
        'month':months,
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
    starting_day = datetime(int(year[0]), int(months[0]), int(days[0]))
    time = np.array([starting_day + timedelta(hours=int(hr)) for hr in hours_in])
    
    ds.coords['valid_time'] = time
    ds = ds.resample(valid_time='1D').mean()
    
    return ds

#%% data
loc_ind = 0

year = 2010

for month in [1,2]:

    ## get wind dataset
    ds = get_era5_wind(year, [month, month+1], 
                       np.arange(1,calendar.monthrange(year, month)[-1]+1), loc_ind)
    
    lon, lat = np.meshgrid(np.array(ds['longitude'].values), np.array(ds['latitude'].values))




#%% end
