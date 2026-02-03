#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 2025

era_tseries1.py
> significant wave height
> air temperature

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
import gc

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

#%% VARIABLE

# var = 'swh'
var = 't2m'

if var=='swh': long_name = "significant_height_of_combined_wind_waves_and_swell"
elif var=='t2m': long_name = "2m_temperature"
else: raise NameError('Incorrect Variable Name')

#%% STARTING INFO
REMOTE = True

LOC = 0 #0=arctic,1=antarctic

# years = np.arange(2010, 2019+1)
# years = np.arange(1982, 1991+1)

# years = [2013]
years=[1985]
# years = np.arange(1982,1992)
MONTHS = [11]
# MONTHS = [12]
# MONTHS = [4,10]
# MONTHS = np.arange(1,12+1)

if not REMOTE:
    root_path = '/Users/mundi/Desktop/month-hemi/'
elif REMOTE:
    root_path = '/home/mundi/month-hemi/'
data_path = ['nh_data/', 'sh_data/']

savepath = root_path+ data_path[LOC]
savedir = var+'/'
ncnameadd= '_'+var
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

#%% era function
def get_era5_var(var, loc_ind, year, month, days):
    '''
    Parameters
    ----------
    year : list
       ex. ['2010']
    months, days: lists

    Returns
    -------
    ds: xarray dataset wit variable
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
    params = {
        "data_format": "netcdf",
        "download_format": "unarchived",
        "product_type": ["reanalysis"],
        "variable": [long_name],
        'year':year,
        'month':month,
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
        
    ds = ds.resample(valid_time='1D').mean() # daily
    
    return ds

#%% calculation functions

def get_miz_area(storm_event):
    t1 = storm_event[0] - timedelta(days=1)
    t2 = storm_event[-1] + timedelta(days=1)
    storm_range = mf.daterange(t1, t2, dt=24)
    
    miz_points = np.zeros(np.shape(si_lon))
    for date in storm_range:
        sic = mf.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
        
    miz_grd = griddata((si_lon.flatten(), si_lat.flatten()), miz_points.flatten(),
                                (lon.flatten(), lat.flatten()), method='nearest')
    miz_grd = miz_grd.reshape(np.shape(lon))
    return miz_grd #miz_points


def calc_timeseries(storm_event):
    from datetime import timedelta
    import functions as mf
    from glob import glob
    import xarray as xr
    from scipy.interpolate import griddata
    
    stormstr_prev=''
        
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
    inside_points1000 = mf.find_points_in_contour(bbox_edges, lon, lat)
    
    #### WINDS!
    def calc_values(dt):
        # load winds!
        try:
            values = var_ds[var].sel(valid_time=dt).values
        except:
            try: # previous file
                values = prev_ds[var].sel(valid_time=dt).values
            except: # next file
                try:
                    values = next_ds[var].sel(valid_time=dt).values
                except:
                    print(dt)
                    ds1 = get_era5_var(var, LOC, dt.year, dt.month, [dt.day])
                    values = ds1[var].sel(valid_time=dt).values
        
        # wind timeseries
        var_in_tot = np.ma.masked_array(values, mask=inside_points1000).filled(np.nan)
        var_miz = np.ma.masked_where(miz_points==0, var_in_tot).filled(np.nan)
        var_out = np.nanmean(var_miz)
        return var_out
    
    # def run_threaded_timeseries(analysis_range):
    #     with ThreadPoolExecutor(max_workers=5) as executor:
    #         return executor.map(calc_values, analysis_range)
        
    # tseries = run_threaded_timeseries(analysis_range)
    # return [x for x in tseries]
    
    return [calc_values(dt) for dt in analysis_range]
        
                
        
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
        
    ds = get_era5_var(var, LOC, start.year, start.month, days_in)
    return ds
    
#%% YEAR LOOP!

for year in years:
    try:
        for month in MONTHS:
            gc.collect()
            
            try:
                print();print('***',year,'-',month,'***')
                name = var+'_'+str(year)+'-'+str(month)
                
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
                def get_var_threaded(data_range):
                    with ThreadPoolExecutor(max_workers=3) as executor:
                        return executor.map(get_data, data_range)
                
                results = get_var_threaded([[storm_ranges[0][0], timedelta(days=0)],
                                              [storm_ranges[0][0], timedelta(days=-10), ],
                                              [storm_ranges[0][0], timedelta(days=32), ]])
                ds_ranges = [x for x in results]
                
                var_ds = ds_ranges[0]
                prev_ds = ds_ranges[1]
                next_ds = ds_ranges[2]
                lon, lat = np.meshgrid(np.array(var_ds['longitude'].values), np.array(var_ds['latitude'].values))
                
                    
                # @timethis 
                # def run_threaded_winds(storm_ranges):
                #     with ThreadPoolExecutor(max_workers=5) as executor:
                #         return executor.map(calc_timeseries, storm_ranges)
                    
                # print('Getting results... ')
                # results1 = run_threaded_winds(storm_ranges)
                # print('Done getting results')
        
                # print('Organizing futures... ')
                # wind_series = [x for x in results1]
                # print('Retreived futures, on to saving')
                
                # try:
                #     np.save(savepath + savedir + name +'.npy', wind_series)
                #     print('SAVED: '+ savedir + name +'.npy', '(1)')
                # except:
                #     np.save(savepath + savedir + name +'.npy', wind_series)
                #     print('SAVED: '+ savedir + name +'.npy', '(2)')
                
                wind_series = []
                for storm_range in enumerate(storm_ranges):
                    wind_series.append( calc_timeseries(storm_range) )
                    
                try:
                    np.save(savepath + savedir + name +'.npy', wind_series)
                    print('SAVED: '+ savedir + name +'.npy', '(1)')
                except:
                    np.save(savepath + savedir + name +'.npy', wind_series)
                    print('SAVED: '+ savedir + name +'.npy', '(2)')
                
                
                    
            except Exception as e1:
                print(e1)
                print('---->> trying again... ')
                gc.collect()
                
                try:
                    print();print('***',year,'-',month,'***')
                    name = var+'_'+str(year)+'-'+str(month)
                    
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
                    def get_var_threaded(data_range):
                        with ThreadPoolExecutor(max_workers=3) as executor:
                            return executor.map(get_data, data_range)
                    
                    results = get_var_threaded([[storm_ranges[0][0], timedelta(days=0)],
                                                  [storm_ranges[0][0], timedelta(days=-10), ],
                                                  [storm_ranges[0][0], timedelta(days=32), ]])
                    ds_ranges = [x for x in results]
                    
                    var_ds = ds_ranges[0]
                    prev_ds = ds_ranges[1]
                    next_ds = ds_ranges[2]
                    lon, lat = np.meshgrid(np.array(var_ds['longitude'].values), np.array(var_ds['latitude'].values))
                    
                        
                    @timethis 
                    def run_threaded_winds(storm_ranges):
                        with ThreadPoolExecutor(max_workers=10) as executor:
                            return executor.map(calc_timeseries, storm_ranges)
                        
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
                except Exception as e2:
                    print(e2)
                    print('---->> rerun month: '+str(year)+'-'+str(month))
                    continue

    except:
        continue
    
    
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


#%% check data
if False: # check missing data
    import numpy as np

    # var = 'swh'
    var = "t2m"
    
    # root_path = '/home/mundi/month-hemi/'
    root_path = '/Users/mundi/Desktop/month-hemi/'
    data_path = ['nh_data/', 'sh_data/']

    savedir = var+'/'
    print('-->', var)
    for loc_ind in [0,1]:
        path1 = root_path+ data_path[loc_ind]
        print('***', loc_ind, '***')
        for years in [np.arange(2010,2020), np.arange(1982,1992)]:
            for year in years:
                for month in np.arange(1,12+1):
                   try: np.load(path1+savedir+var+'_'+str(year)+'-'+str(month)+'.npy', allow_pickle=True)
                   except FileNotFoundError: print(str(year)+'-'+str(month))






#%% end
