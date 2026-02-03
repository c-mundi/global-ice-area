#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 2026
insolation.py

practice calculations for polar solar radiation

https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture11%20--%20Insolation.html

@author: mundi
"""

#%% imports and files

import numpy as np
import matplotlib.pyplot as plt
from climlab import constants as const
from climlab.solar.insolation import daily_insolation
from climlab.solar.orbital import OrbitalTable

import calendar
import functions as fx
from datetime import datetime, timedelta

ice_fname = '/Users/mundi/Desktop/seaice/'

import cartopy.crs as ccrs
import time

#%% daily - one year
print('daily calculation')
start_timer = time.time()

AREA = 25*25
year = 2010

insolation = []

monthly_data = {mm:[] for mm in np.arange(1,12+1)}
for dt in fx.daterange(datetime(year,1,1), datetime(year,12,31), dt=24):
    
    if dt.day==1: print(dt.month, end=' ')
    
    daynum = dt.timetuple().tm_yday
    si, si_lon, si_lat = fx.load_seaice(ice_fname, dt.year, dt.month, dt.day)
    
    # get daily miz
    si_extent = np.where(si<0.15, np.nan, si)
    si_miz = np.where(si>0.80, np.nan, si_extent)
    miz_mask = np.where(np.isnan(si_miz), True, False)
    
    # get miz latitudes
    lats_miz = np.ma.masked_array(np.where(si_lat<60, np.nan, si_lat), miz_mask).filled(np.nan).flatten()
    lats_miz = lats_miz[~np.isnan(lats_miz)]
    
    # for each lat
    watts = [daily_insolation(lat,daynum)*AREA for lat in lats_miz]
    total_area = np.nansum(np.ma.masked_array(np.ones(np.shape(si_lat))*AREA, miz_mask).filled(np.nan))
    area_weighted = np.nansum(watts)/total_area
    
    # print(area_weighted)
    # print(daily_insolation(np.nanmean(lats_miz),daynum))
    
    insolation.append(area_weighted)
    monthly_data[dt.month].append(area_weighted)    

monthly_mean = [np.nanmean(monthly_data[mm]) for mm in np.arange(1,12+1)]

print(); print('time: '+str(round(time.time()-start_timer, 2))+' seconds'); print()

#%% monthly mean(s)
print('monthly calculation')
start_timer = time.time()

AREA = 25*25
year = 2010

monthly_insolation = []
monthly_x = []

for month in np.arange(1,12+1):
    print(month, end=' ')
    
    end_d = calendar.monthrange(year,month)[-1]
    start_daynum = datetime(year,month,1).timetuple().tm_yday
    end_daynum = datetime(year,month,end_d).timetuple().tm_yday
    
    monthly_seaice = []
    for dt in fx.daterange(datetime(year,month,1), datetime(year, month, end_d)):
        si, si_lon, si_lat = fx.load_seaice(ice_fname, dt.year, dt.month, dt.day)
        monthly_seaice.append(si)
    si = np.mean(monthly_seaice, axis=0)
    
    # get daily miz
    si_extent = np.where(si<0.15, np.nan, si)
    si_miz = np.where(si>0.80, np.nan, si_extent)
    miz_mask = np.where(np.isnan(si_miz), True, False)
    
    # get miz latitudes
    lats_miz = np.ma.masked_array(np.where(si_lat<60, np.nan, si_lat), miz_mask).filled(np.nan).flatten()
    lats_miz = lats_miz[~np.isnan(lats_miz)]
    
    # for each lat
    monthly_weighted = []
    for daynum in np.arange(start_daynum, end_daynum+1):
        watts = [daily_insolation(lat,daynum)*AREA for lat in lats_miz]
        total_area = np.nansum(np.ma.masked_array(np.ones(np.shape(si_lat))*AREA, miz_mask).filled(np.nan))
        monthly_weighted.append( np.nansum(watts)/total_area )
    
    area_weighted = np.nanmean(monthly_weighted)
    monthly_insolation.append(area_weighted)
        
    monthly_x.append(np.nanmean([start_daynum, end_daynum]))

print(); print('time: '+str(round(time.time()-start_timer, 2))+' seconds'); print()

#%% sample plot

plt.figure()

plt.plot(np.arange(0,365), insolation, color='k', label='daily')
plt.plot(monthly_x, monthly_mean, marker='o', color='gray', label='monthly mean')

plt.plot(monthly_x, monthly_insolation, marker='X', color='maroon', label='monthly calculation')


#%% antarctica

print('daily calculation')
start_timer = time.time()

AREA = 25*25

insolation_sh = []

monthly_data_sh = {mm:[] for mm in np.arange(1,12+1)}
for dt in fx.daterange(datetime(year,1,1), datetime(year,12,31), dt=24):
    
    if dt.day==1: print(dt.month, end=' ')
    
    daynum = dt.timetuple().tm_yday
    si, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', dt.year, dt.month, dt.day)
    
    # get daily miz
    si_extent = np.where(si<0.15, np.nan, si)
    si_miz = np.where(si>0.80, np.nan, si_extent)
    miz_mask = np.where(np.isnan(si_miz), True, False)
    
    # get miz latitudes
    lats_miz = np.ma.masked_array(np.where(si_lat>-55, np.nan, si_lat), miz_mask).filled(np.nan).flatten()
    lats_miz = lats_miz[~np.isnan(lats_miz)]
    
    # for each lat
    watts = [daily_insolation(lat,daynum)*AREA for lat in lats_miz]
    total_area = np.nansum(np.ma.masked_array(np.ones(np.shape(si_lat))*AREA, miz_mask).filled(np.nan))
    area_weighted = np.nansum(watts)/total_area
    
    insolation_sh.append(area_weighted)
    monthly_data_sh[dt.month].append(area_weighted)    

monthly_mean_sh = [np.nanmean(monthly_data_sh[mm]) for mm in np.arange(1,12+1)]

print(); print('time: '+str(round(time.time()-start_timer, 2))+' seconds'); print()

#%%% plot

plt.figure()

plt.plot(np.arange(0,365), insolation_sh, color='k', label='daily')
plt.plot(monthly_x, monthly_mean_sh, marker='o', color='gray', label='monthly mean')

#%% xxxxxxxxxxxxxxxx

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

#%% convert to function

@timethis
def get_annual_insolation(loc_ind, year, ice_fname, AREA=25*25):
    daily_values = []
    monthly_data = {mm:[] for mm in np.arange(1,12+1)}
    
    # get sea ice lon/lat
    if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, datetime(year,1,1).year, datetime(year,1,1).month, datetime(year,1,1).day)
    elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, datetime(year,1,1).year, datetime(year,1,1).month, datetime(year,1,1).day)
    
    # daily loop
    for dt in fx.daterange(datetime(year,1,1), datetime(year,12,31), dt=24):
        # load daily ice
        daynum = dt.timetuple().tm_yday
        if loc_ind==0: si = fx.load_seaice(ice_fname, dt.year, dt.month, dt.day, latlon=False)
        elif loc_ind==1: si = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day, latlon=False)
        
        # get daily miz
        si_extent = np.where(si<0.15, np.nan, si)
        si_miz = np.where(si>0.80, np.nan, si_extent)
        miz_mask = np.where(np.isnan(si_miz), True, False)
        
        # get miz latitudes
        if loc_ind==0: lats_miz = np.ma.masked_array(np.where(si_lat<60, np.nan, si_lat), miz_mask).filled(np.nan).flatten()
        elif loc_ind==1: lats_miz = np.ma.masked_array(np.where(si_lat>-55, np.nan, si_lat), miz_mask).filled(np.nan).flatten()
        lats_miz = lats_miz[~np.isnan(lats_miz)]
        
        # for each lat, calculate solar energy
        watts = [daily_insolation(lat,daynum)*AREA for lat in lats_miz]
        total_area = np.nansum(np.ma.masked_array(np.ones(np.shape(si_lat))*AREA, miz_mask).filled(np.nan))
        area_weighted = np.nansum(watts)/total_area # divide by total MIZ area
        
        daily_values.append(area_weighted)
        monthly_data[dt.month].append(area_weighted)    
    
    monthly_mean = [np.nanmean(monthly_data[mm]) for mm in np.arange(1,12+1)]

    return daily_values, monthly_mean


daily_curves = {0:[], 1:[]}
monthly_curves = {0:[], 1:[]}
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    print('***', loc, '***')
    for year in np.arange(2010,2020):
        print(year, end=' ')
        daily_values, monthly_mean = get_annual_insolation(loc_ind, year, 
                                                               ice_fname if loc_ind==0 else ice_fname+'south/')  
        daily_curves[loc_ind].append(daily_values)
        monthly_curves[loc_ind].append(monthly_mean)

#### pickle data
import pickle
savepath = '/Users/mundi/Desktop/month-hemi/'
pickle.dump(daily_curves, open(savepath+'insolation_daily.pkl', 'wb'))
pickle.dump(monthly_curves, open(savepath+'insolation_monthly.pkl', 'wb'))

#%%% plot

daily_values = pickle.load(open(savepath+'insolation_daily.pkl', 'rb'))
monthly_values = pickle.load(open(savepath+'insolation_monthly.pkl', 'rb'))

monthly_x = []
for month in np.arange(1,12+1):
    daynum_range = [datetime(year,month,1).timetuple().tm_yday,
                    datetime(year,month,calendar.monthrange(year,month)[-1]).timetuple().tm_yday]
    monthly_x.append(np.nanmean(daynum_range))
    


plt.figure()

for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    
    for dcurve in daily_values[loc_ind]:
        plt.plot(np.arange(0,len(dcurve)), dcurve, color='gray', lw=0.55)
    
    plt.plot(monthly_x, np.nanmean(monthly_values[loc_ind], axis=0), lw=2, marker='o',
             color=['#01665e','#8c510a'][loc_ind], label=loc)


    plt.xlabel('Day')
    plt.ylabel(r'MIZ-Area-Weighted Insolation (W m$^{-2}$)')
    plt.legend(ncol=2)


#%% end