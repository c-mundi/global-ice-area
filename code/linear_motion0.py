#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 2025

get data

@author: mundi
"""
#%% imports, file names
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import functions as fx
from datetime import datetime, timedelta
from scipy.interpolate import griddata
import warnings
from glob import glob
import gc


ROOT = '/Users/mundi/Desktop/'
ice_path = ROOT +'seaice/'
ice_paths = [ice_path, ice_path+'south/']

ice_motion_path = '/Users/mundi/Desktop/seaice/ice_motion/'

nh_path = ROOT +'month-hemi/nh_data/'
sh_path = ROOT +'month-hemi/sh_data/'
root_paths = [nh_path, sh_path]


savepath = ROOT+'month-hemi/linear_motion/'
savepath_miz = '/Users/mundi/Desktop/month-hemi/linear_miz/'

census_name = 'census_'
contour_name = '_contours.nc'
si_name = '_seaice'

#%%% controls

years = np.arange(2010,2020)
# years = np.arange(1982, 1992)

months = np.arange(1,12+1)

calc_time = 0 # 0 = storm duration, 1 = 7 days, 2 = 3 weeks

xxx = np.arange(-7,14+1,1)
ice_lims = [20,80]
miz = [0.15,0.80]

lon_num = 100
lon_x = np.linspace(0, 1, lon_num)

DATA = 'v'

#%%% functions
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
        
    ds.close()
        
    return lon, lat, E, N


def timediff(calc_time, start, enddate):
    if calc_time==0:
        dt1 = start
        dt2 = enddate
    elif calc_time==1:
        dt1 = start
        dt2 = start + timedelta(days=7)
    elif calc_time==2:
        dt1 = start - timedelta(days=7)
        dt2 = start + timedelta(days=14)
    else:
        raise NameError("Incorrect Time Range")
    return dt1, dt2

#%%% calculations

def get_miz_area(storm_range, si_lon, si_lat, lon, lat, loc_ind, ice_fname):
    miz_points = np.zeros(np.shape(si_lon))
    for date in storm_range:
        if loc_ind==0: sic = fx.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
        elif loc_ind==1: sic = fx.load_seaice_sh(ice_fname, date.year, date.month, date.day, latlon=False)
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
    miz_grd = griddata((si_lon.flatten(), si_lat.flatten()), miz_points.flatten(),
                                (lon.flatten(), lat.flatten()), method='nearest')
    miz_grd = miz_grd.reshape(np.shape(lon))
    
    return miz_grd

def get_motion_values(dt1, dt2, loc_ind, ice_motion_path):
    east, north = [],[]
    for date in fx.daterange(dt1, dt2, dt=24):
        lon, lat, E, N = load_ice_motion(ice_motion_path, loc_ind, date.year, date.month, date.day)
        east.append(E); north.append(N)
    
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore", category=RuntimeWarning)
        u_wind = np.squeeze( np.nanmean(east, axis=0) )
        v_wind = np.squeeze( np.nanmean(north, axis=0) )
        
    return u_wind, v_wind

#%% data loop

from concurrent.futures import ThreadPoolExecutor

@timethis
def get_yearly_data(IN):
    year = IN[0]
    loc_ind = IN[1]
    root_path = IN[2]
    si_lon, si_lat = IN[3]
    lon, lat = IN[4]
    ice_fname, ice_motion_path = IN[5]
    
    print(year)
    
    savename = str(loc_ind)+'_'+str(calc_time)+'_'+str(year)+'_'+DATA+'.npy'
    
    # load census and sea ice info
    census_name = 'census_'+ str(year) +'.csv'
    [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
        fx.readCensus(root_path+'census/'+census_name, convertDT=True)
    
    # storm loop
    lines = []
    for sn, start in enumerate(startdate):
        end = enddate[sn]
        gc.collect()
        
        if sn==np.floor(len(startdate)/2): print('half',year)
        
        # storm info
        stormstr = start.strftime('%Y_%m%d')
        storm_range = fx.daterange(start, end, dt=24)
        
        # contours and bbox
        cs = xr.open_dataset(root_path+'contours/' + stormstr + contour_name)
        all_contours = []
        for key in list(cs.keys()):
            coord = cs[key].values
            all_contours.append(coord)
        cs.close()
        del cs
        with fx.HidePrint(): bbox_edges = fx.get_bbox_edges(all_contours) 
        
        # compute miz
        miz_grd = get_miz_area(storm_range, si_lon, si_lat, lon, lat, loc_ind, ice_fname)
        
        # calculate mean winds
        dt1,dt2 = timediff(calc_time, start, end)
        u_wind, v_wind = get_motion_values(dt1, dt2, loc_ind, ice_motion_path)
        if DATA=='u': wind = u_wind
        elif DATA=='v': wind = v_wind
        else: raise NameError('wrong DATA value: u or v only')

        # isolate miz
        inside_bbox = fx.find_points_in_contour(bbox_edges, lon, lat)
        w_bbox = np.ma.masked_array(wind, mask=inside_bbox)
        w_miz = np.ma.masked_where(miz_grd==0, w_bbox).filled(np.nan)
        
        # set up grid
        lon_bnds = np.linspace( np.nanmin(bbox_edges[:,0]), np.nanmax(bbox_edges[:,0]), lon_num)
        lat_bnds = np.linspace( np.nanmin(bbox_edges[:,1]), np.nanmax(bbox_edges[:,1]), lon_num)
        new_x, new_y = np.meshgrid(lon_bnds, lat_bnds)   

        #### calculate final values
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore", category=RuntimeWarning)
            
            w_gridded = griddata((lon.flatten(), lat.flatten()), w_miz.flatten(),
                                   (new_x.flatten(), new_y.flatten()))
            w_gridded = np.reshape(w_gridded, np.shape(new_x))
            LINE = np.nanmean(w_gridded, axis=0)
        lines.append(LINE)
            
    # export data
    np.save(savepath+savename, np.array(lines))
    return lines
    
#### data loop

@timethis
def threaded_get_lines(years):
    with ThreadPoolExecutor(max_workers=5) as executor:
        return executor.map(get_yearly_data, years)

if __name__ == '__main__':
    for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
        print(); print('***', loc, '***')
        root_path = root_paths[loc_ind]
        ice_fname = ice_paths[loc_ind]
        
        if loc_ind == 0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010,1,1, latlon=True)
        elif loc_ind ==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010,1,1, latlon=True)
        
        ## get ice motion lon/lat
        lon, lat, _, _ = load_ice_motion(ice_motion_path, loc_ind, 2010, 1, 1)
        
        # for year in years:
        #     print(year)
        #     savename = loc+'_'+str(calc_time)+'_'+str(year)+'_'+DATA+'.npy'
            
        #     try:
        #         lines = np.load(savepath+savename)
        #     except FileNotFoundError:
        #         lines = get_yearly_data(year)
        
        process_years = []
        for year in years:
            try:
                savename = str(loc_ind)+'_'+str(calc_time)+'_'+str(year)+'_'+DATA+'.npy'
                lines = np.load(savepath+savename)
            except FileNotFoundError:
                file_paths = (ice_fname, ice_motion_path)
                process_years.append( (year,loc_ind, root_path,(si_lon, si_lat),(lon,lat), file_paths) )
                
        ## PARALLEL DATA
        tseries = threaded_get_lines(process_years)
        lines =  [x for x in tseries] 

#%% end
