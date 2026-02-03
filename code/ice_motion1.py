#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 2025
ice_motion1.py

> explore ice motion data
> develop function for vector transformations

data: https://nsidc.org/data/nsidc-0116/versions/4

@author: mundi
"""
#%% imports and file paths
import numpy as np
import xarray as xr
from glob import glob
from datetime import datetime
import functions as fx

ice_fname = '/Users/mundi/Desktop/seaice/'
ice_path = '/Users/mundi/Desktop/seaice/ice_motion/'

decades = [np.arange(1982,1992), np.arange(2010,2020)]
hemi_names= ['Arctic', 'Antarctic']


#%% sample data
sample_year = 2010
fname = 'nh'+ '_25km_'+str(sample_year)+'0101_'+str(sample_year)+'1231_v4.1.nc'
file = glob(ice_path+'*'+fname)[0]

ds = xr.open_dataset(file)

datetimeindex = ds.indexes['time'].to_datetimeindex(time_unit='s')
ds['time'] = datetimeindex

lon = ds['longitude'].values
lat = ds['latitude'].values
u = ds['u'].sel(time=slice(datetime(2010,8,15), datetime(2010,8,17))).values
v = ds['v'].sel(time=slice(datetime(2010,8,15), datetime(2010,8,17))).values

#%%% plot
import cartopy.crs as ccrs

ax = fx.background_plot_nh(extent=(0,360,60,90))

interval=7

ax.quiver(lon[::interval, ::interval], lat[::interval, ::interval], 
          u[0][::interval, ::interval], v[0][::interval, ::interval], 
          transform=ccrs.PlateCarree())

# ax.quiver(x_greater[::interval, ::interval], y_greater[::interval, ::interval],
#         u_greater[::interval, ::interval], v_greater[::interval, ::interval],
#         color='k', transform = ccrs.PlateCarree(), scale = sc, zorder = 6,
#         headwidth=2.5, headlength=3)

si, si_lon, si_lat = fx.load_seaice(ice_fname, 2010,8,15)
[si_lon1, si_lon2], [si_lat1, si_lat2], [si1,si2] = fx.geoplot_2d(si_lon, si_lat, si)
ax.contour(si_lon1, si_lat1, si1, levels=[0.15], transform=ccrs.PlateCarree())
ax.contour(si_lon2, si_lat2, si2, levels=[0.15], transform=ccrs.PlateCarree())

#%%% plot 2
ax = fx.background_plot_nh(extent=(0,360,60,90))

U = u[0]; V = v[0]
# E:   u * cos L  +  v * sin L
# N:  -u * sin L  +  v * cos L
E = ( U*np.cos(np.deg2rad(lon)) ) + ( V*np.sin(np.deg2rad(lon)) )
N = ( -1*U*np.sin(np.deg2rad(lon)) ) + ( V*np.cos(np.deg2rad(lon)) )

interval=7
ax.quiver(lon[::interval, ::interval], lat[::interval, ::interval], 
          E[::interval, ::interval], N[::interval, ::interval], 
          transform=ccrs.PlateCarree())

# ax.quiver(x_greater[::interval, ::interval], y_greater[::interval, ::interval],
#         u_greater[::interval, ::interval], v_greater[::interval, ::interval],
#         color='k', transform = ccrs.PlateCarree(), scale = sc, zorder = 6,
#         headwidth=2.5, headlength=3)

si, si_lon, si_lat = fx.load_seaice(ice_fname, 2010,8,15)
[si_lon1, si_lon2], [si_lat1, si_lat2], [si1,si2] = fx.geoplot_2d(si_lon, si_lat, si)
ax.contour(si_lon1, si_lat1, si1, levels=[0.15], transform=ccrs.PlateCarree())
ax.contour(si_lon2, si_lat2, si2, levels=[0.15], transform=ccrs.PlateCarree())

#%% function

def load_ice_motion(ice_path, loc_ind, year, month, day):
    
    fname = ['nh','sh'][loc_ind]+ '_25km_'+str(year)+'0101_'+str(year)+'1231_v4.1.nc'
    file = glob(ice_path+'*'+fname)[0]

    ds = xr.open_dataset(file)

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


#%% end
# https://nsidc.org/data/user-resources/help-center/how-convert-horizontal-and-vertical-components-east-and-north
