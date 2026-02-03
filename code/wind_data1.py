#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 2025
wind_data2a.py

get arctic winds winds!
> fewer api calls
> PARALLEL
* monthly files

@author: mundi
"""
import os
import numpy as np
import glob
import xarray as xr
from datetime import datetime, timedelta
import traceback
from scipy.interpolate import griddata
import string
import functions as mf
import calendar

from concurrent.futures import ThreadPoolExecutor
import time
from functools import wraps

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
        time_start = time.time()
        output = func(*args, **kwargs)
        time_end = time.time()
        if time_end-time_start < 120:
            print(f"{func.__name__}: {(time_end-time_start)} s")
        else:
            print(f"{func.__name__}: {(time_end-time_start)/60} min")

        return output

    return wrapped_method

STARTTIME= time.time()

#%% STARTING INFO
REMOTE = True

LOC = 0 #0=arctic,1=antarctic

# years = np.arange(2010, 2019+1)
# years = np.arange(1982, 1991+1)
# years=[2016]
years = [1991]
# years = np.arange(2017,2020)
# years = np.arange(1986, 1992)


# MONTHS = np.arange(6,12+1)
# MONTHS = [10,11,12]
# MONTHS = np.arange(1,6)
# MONTHS = np.arange(8,12+1)
MONTHS = [4]

if not REMOTE:
    root_path = '/Users/mundi/Desktop/month-hemi/'
elif REMOTE:
    root_path = '/home/mundi/month-hemi/'
data_path = ['nh_data/', 'sh_data/']

savepath = root_path+ data_path[LOC]
savedir = 'wind/'
ncnameadd= '_wind'
if not os.path.exists(savepath+ savedir):
    os.makedirs(savepath+savedir)

census_path =  root_path+ data_path[LOC] + 'census/'

if not REMOTE:
    ice_add = 'south/' if LOC==1 else ''
    ice_fname = '/Users/mundi/Desktop/seaice/'+ice_add
elif REMOTE:
    ice_add = 'south/' if LOC==1 else 'daily/'
    ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/'+ice_add

contour_path =  root_path+ data_path[LOC] + 'contours/'
    
census_name = 'census_'
contour_name = '_contours.nc'
si_name = '_seaice'
    
### get starting sea ice
if LOC==0: _, si_lon, si_lat = mf.load_seaice(ice_fname, 2000, 10, 1, latlon=True)
else: _, si_lon, si_lat = mf.load_seaice_sh(ice_fname, 2000, 10, 1, latlon=True)
area = 25*25  
miz = [0.15,0.80] #[0,0.58] #

def rstr(val, n=3): return str(round(val, n))

#%% era winds function

def get_era5_wind(year, month, days):
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


#%% calculation function!

def get_miz_area(storm_event):
    t1 = storm_event[0] - timedelta(days=1)
    t2 = storm_event[-1] + timedelta(days=1)
    storm_range = mf.daterange(t1, t2, dt=24)
    
    miz_points = np.zeros(np.shape(si_lon))
    for date in storm_range:                    ### load sea ice (hemis!!!)
        sic = mf.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
        
    miz_grd = griddata((si_lon.flatten(), si_lat.flatten()), miz_points.flatten(),
                                (lon.flatten(), lat.flatten()), method='nearest')
    miz_grd = miz_grd.reshape(np.shape(lon))
    return miz_grd #miz_points
    

def calc_wind_timeseries(storm_event):
    from datetime import timedelta
    import functions as mf
    from glob import glob
    import xarray as xr
    from scipy.interpolate import griddata
    
    stormstr_prev=''
        
    # print('*', storm_event[0])
    
    #### storm info
    stormstr1 = storm_event[0].strftime('%Y_%m%d')
    week_ago = storm_event[0] - timedelta(days=7)
    two_week = storm_event[0] + timedelta(days=14) # relative to start date, since different storm lengths
    analysis_range = mf.daterange(week_ago, two_week, dt=24)
    
    # duplicate storm start date?
    if stormstr1==stormstr_prev:
        stormstr = stormstr1 + '*'
        dupe=True
    else:
        dupe=False
        stormstr=stormstr1
    stormstr_prev = stormstr1
    
    #### miz area
    miz_points = get_miz_area(storm_event)
        
    #### storm area
    ncname = stormstr + contour_name
    try:
        if not dupe: ###!!!
            cs = xr.open_dataset(contour_path+ncname)
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        contour_files = glob(contour_path+'*'+storm_event[0].strftime('%Y_%m%d')+'*.nc')
        cs = xr.open_dataset(contour_files[-1])
        
    all_contours = []
    for key in list(cs.keys()):
        coord = cs[key].values
        all_contours.append(coord)
    cs.close()
    
    if all_contours == []: print('empty contours'); return []
    
    with mf.HidePrint(): bbox_edges = mf.get_bbox_edges(all_contours) 
    # inside_points1000 = mf.find_points_in_contour(bbox_edges, si_lon, si_lat)
    inside_points1000 = mf.find_points_in_contour(bbox_edges, lon, lat)
    
    #### WINDS!
    def calc_winds(dt):
        # load winds!
        try:
            u_wind = ds_w['u10'].sel(valid_time=dt).values
            v_wind = ds_w['v10'].sel(valid_time=dt).values
        except:
            try: # previous file
                u_wind = prev_wind['u10'].sel(valid_time=dt).values
                v_wind = prev_wind['v10'].sel(valid_time=dt).values
            except: # next file
                try:
                    u_wind = next_wind['u10'].sel(valid_time=dt).values
                    v_wind = next_wind['v10'].sel(valid_time=dt).values
                except:
                    print(dt)
                    ds = get_era5_wind(dt.year, dt.month, [dt.day])
                    u_wind = ds['u10'].sel(valid_time=dt).values
                    v_wind = ds['v10'].sel(valid_time=dt).values
        wind_tot = np.array( np.sqrt((u_wind**2) + (v_wind**2)) )
        
        # wind timeseries
        wind_out = []
        for tw, thiswind in enumerate([wind_tot, u_wind, v_wind]):
             wind_in_tot = np.ma.masked_array(thiswind, mask=inside_points1000).filled(np.nan)
             total_wind_miz = np.ma.masked_where(miz_points==0, wind_in_tot).filled(np.nan)
             wind_out.append(np.nanmean(total_wind_miz))
        return wind_out
    
    # @timethis 
    def run_threaded_timeseries(analysis_range):
        with ThreadPoolExecutor(max_workers=10) as executor:
            return executor.map(calc_winds, analysis_range)
        
    tseries = run_threaded_timeseries(analysis_range)
    return [x for x in tseries]
                
        
def get_data(data_range):
    start = data_range[0] + data_range[1]
    
    month = start.month
    if data_range[1].days==0:
        days_in = [LZ(dd+1) for dd in range(calendar.monthrange(years[0]-1, month)[1])]
    elif data_range[1].days < 0:
        day_range = calendar.monthrange(start.year, start.month)[1]
        days_in = [LZ(dd) for dd in np.arange(day_range-10, day_range+1)]
    elif data_range[1].days > 0:
        days_in = [LZ(dd) for dd in np.arange(1, 15)]
        
        
    ds_w = get_era5_wind(start.year, start.month, days_in)
    return ds_w
    
        
#%% YEAR LOOP!

for year in years:
    
    for month in MONTHS:
        
        try:
            print();print('***',year,'-',month,'***')
            name = 'winds_'+str(year)+'-'+str(month)
            
            ### get start/end dates 
            [startdate, enddate] = mf.readCensus(census_path+census_name + str(year)+'.csv' , convertDT=True)[0]
            
            storm_ranges = []
            for startdt, enddt in zip(startdate, enddate):
                if startdt.month != month: continue
                storm_ranges.append(mf.daterange(startdt, enddt, dt=24))  
                
            if storm_ranges == []: 
                print('No storms!')
                np.save(savepath + savedir + name +'.npy', [])
                continue
                
            # load wind files
            @timethis 
            def get_winds_threaded(data_range):
                with ThreadPoolExecutor(max_workers=3) as executor:
                    return executor.map(get_data, data_range)
            
            results = get_winds_threaded([[storm_ranges[0][0], timedelta(days=0)],
                                          [storm_ranges[0][0], timedelta(days=-10), ],
                                          [storm_ranges[0][0], timedelta(days=32), ]])
            wind_ds = [x for x in results]
            
            ds_w = wind_ds[0]
            prev_wind = wind_ds[1]
            next_wind = wind_ds[2]
            lon, lat = np.meshgrid(np.array(ds_w['longitude'].values), np.array(ds_w['latitude'].values))
            
                
            @timethis 
            def run_threaded_winds(storm_ranges):
                with ThreadPoolExecutor(max_workers=10) as executor:
                    return executor.map(calc_wind_timeseries, storm_ranges)
                
            print('Getting results... ')
            results1 = run_threaded_winds(storm_ranges)
            print('Done getting results')
    
            print('Organizing futures... ')
            wind_series = [x for x in results1]
            print('Retreived futures, on to saving')
            
            try:
                np.save(savepath + savedir + name +'.npy', wind_series)
                print('SAVED: '+ savedir + name +'.npy', '(1)')
            except:
                np.save(savepath + savedir + name +'.npy', wind_series)
                print('SAVED: '+ savedir + name +'.npy', '(2)')
                
        except:
            print('---->> rerun month: '+str(year)+'-'+str(month))

    
    
    
#%% end
print(); print()
print('Runtime:')
tdiff = (time.time() - STARTTIME)/60
if tdiff>60:
    tdiff /= 60
    units='hours'
else:
    units = 'minutes'

print(round(tdiff, 2), units)

# def interp_winds(wind):
#     wind_grd = griddata((lon.flatten(), lat.flatten()), wind.flatten(),
#                             (si_lon.flatten(), si_lat.flatten()))
#     return wind_grd.reshape(np.shape(si_lon))

# @timethis 
# def run_threaded_interp(winds):
#      with ThreadPoolExecutor(max_workers=3) as executor:
#          return executor.map(interp_winds, winds)

