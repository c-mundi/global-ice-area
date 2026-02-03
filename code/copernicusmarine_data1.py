#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 2025
copernicusmarine_data1.py

get ocean profiles (glorys)

@author: mundi
"""
#%% controls
REMOTE=False

loc_ind =1 #0=arctic,1=antarctic

years = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
# years = [2010, 2011, 2012, 2013, 2014]
# years = [2015, 2016, 2017, 2018, 2019]
# years=[2017, 2018, 2019]
# years = [2016]

#%% imports and data
import numpy as np
import xarray as xr
import os, glob
from datetime import timedelta

import copernicusmarine
copernicusmarine.login(username='cmundi', password='Oceans!52', force_overwrite=True)
# copernicusmarine.login(username='cmundi', password='Oceans!52', overwrite_configuration_file=True)

if not REMOTE:
    import functions as mf
    root = '/Users/mundi/Desktop/month-hemi/'+['nh_data/','sh_data/'][loc_ind]
    ice_fname =  '/Users/mundi/Desktop/seaice/'
    if loc_ind==1: ice_fname+='south/'
    cm_path = root+'ocn_prof/'
else:
    import functions as mf 
    root = '/home/mundi/month-hemi/'+['nh_data/','sh_data/'][loc_ind]
    ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/'
    ice_fname+='south/' if loc_ind==1 else 'daily/'   
    cm_path = cm_path = root+'ocn_prof/'

census_path = root+'census/'
contour_path = root +'contours/'
contour_name = '_contours.nc'
ncpath = root+'area/'


miz = [0.15,0.80]
ice_lims = [20,80]
    
cm_login = False
if not os.path.exists(cm_path):
    os.makedirs(cm_path)
    
    
#%%% fxns

from scipy.interpolate import griddata

if loc_ind==0: _, si_lon, si_lat = mf.load_seaice(ice_fname, 2010, 1, 1, latlon=True)
elif loc_ind==1: _, si_lon, si_lat = mf.load_seaice_sh(ice_fname, 2010, 1, 1, latlon=True)
si_points = [(x,y) for x,y in zip(si_lon.flatten(), si_lat.flatten())]

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
   
    
#%% DATA FUNCTION

def get_cm_data(sn, start, end, bbox_edges):
    # storm info
    starting_str = start.strftime('%Y-%m-%d')
    
    start1 = start - timedelta(days=7)
    end1 = start + timedelta(days=14)
    
    start_str = start1.strftime('%Y-%m-%d')
    end_str = end1.strftime('%Y-%m-%d')
    
    # if os.path.isfile(cm_path+starting_str+'_'+str(sn)+'_miz_series.npy'):
    #     print('*** '+starting_str)
    #     return cm_path+start_str+'_'+str(sn)+'.nc'
    # else:
    #     print('.......... missing: '+starting_str)
    
    try: 
        np.load(cm_path+starting_str+'_'+str(sn)+'_miz_series.npy')
        print('******* '+starting_str)
        return cm_path+start_str+'_'+str(sn)+'.nc'
    except FileNotFoundError:
        print('.......... missing: '+starting_str)
    
    
    ### get bbox limits
    # LONS
    prev_x = [0,0]
    out_lon = []
    for x in list(bbox_edges[:,0]):
        if x == prev_x[0] and x==prev_x[1]:
            out_lon.append(x)
        prev_x[0]=prev_x[1]
        prev_x[1] = x
    # LATS
    prev_x = [0,0]
    out_lat = []
    for x in list(bbox_edges[:,1]):
        if x == prev_x[0] and x==prev_x[1]:
            out_lat.append(x)
        prev_x[0]=prev_x[1]
        prev_x[1] = x
    
    #### GET DATA
    try:
        request_dataframe = copernicusmarine.subset(
            dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
            minimum_longitude = np.nanmin(out_lon),
            maximum_longitude = np.nanmax(out_lon),
            minimum_latitude = np.nanmin(out_lat),
            maximum_latitude = np.nanmax(out_lat),
            variables = ["thetao"],
            start_datetime = start_str,
            end_datetime = end_str,
            minimum_depth=0,
            maximum_depth=90,
            output_filename = start_str+'_'+str(sn)+'.nc',
            output_directory = cm_path
        )
    except Exception as EE:
        print();print();print('Error Loading '+start_str); print()
        print(EE)
        # try again
        print('Download 2...')
        request_dataframe = copernicusmarine.subset(
            dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
            minimum_longitude = np.nanmin(out_lon),
            maximum_longitude = np.nanmax(out_lon),
            minimum_latitude = np.nanmin(out_lat),
            maximum_latitude = np.nanmax(out_lat),
            variables = ["thetao"],
            start_datetime = start_str,
            end_datetime = end_str,
            minimum_depth=0,
            maximum_depth=90,
            output_filename = start_str+'_'+str(sn)+'.nc',
            output_directory = cm_path
        )
        
    #### GET DATA
    try:
        ds_ocn = xr.open_dataset(cm_path+start_str+'_'+str(sn)+'.nc', decode_times=True)
    except:
        files = glob.glob(cm_path+'*'+start_str+'*_*'+str(sn)+'*.nc')
        if len(files)!=0:
            ds_ocn = xr.open_dataset(files[-1], decode_times=True)
        else: 
            # try downloading again
            print('Download 3...')
            request_dataframe = copernicusmarine.subset(
                dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
                minimum_longitude = np.nanmin(out_lon),
                maximum_longitude = np.nanmax(out_lon),
                minimum_latitude = np.nanmin(out_lat),
                maximum_latitude = np.nanmax(out_lat),
                variables = ["thetao"],
                start_datetime = start_str,
                end_datetime = end_str,
                minimum_depth=0,
                maximum_depth=90,
                output_filename = start_str+'_'+str(sn)+'.nc',
                output_directory = cm_path
            )
            ds_ocn = xr.open_dataset(cm_path+start_str+'_'+str(sn)+'.nc', decode_times=True)
    lon = ds_ocn['longitude'].values
    lat = ds_ocn['latitude'].values
    # depth = ds_ocn['depth'].values
    xx,yy = np.meshgrid(lon, lat)
    sst_points = [(x,y) for x,y in zip(xx.flatten(), yy.flatten())]
    temp = ds_ocn['thetao']
    ds_ocn.close()
    
    # compute miz
    miz_points = np.zeros(np.shape(si_lon))
    storm_range = mf.daterange(start, end, dt=24)
    for date in storm_range:
        sic = mf.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
    
    ocn = griddata(si_points, miz_points.flatten(), sst_points)
    ocn = ocn.reshape(np.shape(xx))
    
    ### save data ###!!! miz vs all storm area
    # mean_profile = temp.mean(dim=('latitude','longitude')).values 
    # np.save(cm_path+start_str+'_'+str(sn)+'_all.npy', mean_profile)
    
    mean_miz = temp.where(ocn==1).mean(dim=('latitude','longitude')).values
    np.save(cm_path+starting_str+'_'+str(sn)+'_miz_series.npy', mean_miz)
    print('*saved!*', starting_str)
    
    # delete big file!
    os.remove(cm_path+start_str+'_'+str(sn)+'.nc')
    print('deleted big file and moving on')
    
    return cm_path+start_str+'_'+str(sn)+'.nc'
   

#%% data loop

def ocn_prof(storm_info):
    sn = storm_info[0]
    start = storm_info[1]
    
    ice_frac = ice_sorter[sn]*100/box_area[sn]             
    if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
        return ''
    
    # contours and bbox
    stormstr = start.strftime('%Y_%m%d')
    cs = xr.open_dataset(contour_path + stormstr + contour_name)
    all_contours = []
    for key in list(cs.keys()):
        coord = cs[key].values
        all_contours.append(coord)
    cs.close()
    del cs
    with mf.HidePrint(): bbox_edges = mf.get_bbox_edges(all_contours) 
    
    # DATA
    return get_cm_data(sn, start, enddate[sn], bbox_edges)


for year in years:
    print(); print(year)
    print('****')
    
    # CENSUS: storm info
    census_file = census_path+'census_'+str(year)+'.csv'
    [startdate, enddate] = mf.readCensus(census_file, convertDT=True)[0]
    
    # open ice area
    ds_area = xr.open_dataset(ncpath + str(year) +'_area.nc')
    ice_area80 = ds_area['ice_area80'].values
    ice_area15 = ds_area['ice_area15'].values
    box_area = ds_area['box_area'].values
    ice_sorter = ice_area80
    
    
    #### loop thru storms
    prev_month=0
    load_summer = True
    for sn, start in enumerate(startdate): 
        month = start.month
        if month!=prev_month:print(month, end=' ')
        prev_month=month
        
        ice_frac = ice_sorter[sn]*100/box_area[sn]             
        if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
            continue
        
        # contours and bbox
        stormstr = start.strftime('%Y_%m%d')
        cs = xr.open_dataset(contour_path + stormstr + contour_name)
        all_contours = []
        for key in list(cs.keys()):
            coord = cs[key].values
            all_contours.append(coord)
        cs.close()
        del cs
        with mf.HidePrint(): bbox_edges = mf.get_bbox_edges(all_contours) 
        
        # DATA
        get_cm_data(sn, start, enddate[sn], bbox_edges)
           
    #### PARALLEL
    # storm_info = [(sn, start) for sn, start in enumerate(startdate)]
    
    # # @timethis 
    # def run_threaded_ocn(storm_info):
    #     with ThreadPoolExecutor(max_workers=5) as executor:
    #         return executor.map(ocn_prof, storm_info)
        
    # tseries = run_threaded_ocn(storm_info)
    # out = [x for x in tseries]


  






#%% end
