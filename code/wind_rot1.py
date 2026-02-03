#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 2025
wind_rot1.py

find NW/SE components

@author: mundi
"""
#%% imports & files
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import functions as fx
import cartopy.crs as ccrs


file = '/Users/mundi/Desktop/era_data/2011_casestudy_winds.nc'

ref_vec = (-1,1)

# projection 
# np.dot(x, y) / np.linalg.norm(y)

test = (-1,-1)

value = np.dot(test, ref_vec) / np.linalg.norm(ref_vec)
# print(value)

#%%% functions
def plot_winds(ax, x,y,u,v, sc=300, interval=12,headwidth=2.5, headlength=3):  

    u = np.squeeze(u)
    v = np.squeeze(v)    
    
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
            headwidth=headwidth, headlength=headlength)
    ax.quiver(x_lesser[::interval, ::interval], y_lesser[::interval, ::interval],
            u_lesser[::interval, ::interval], v_lesser[::interval, ::interval],
            color='k', transform = ccrs.PlateCarree(), scale = sc, zorder = 6,
            headwidth=headwidth, headlength=headlength)
    
    return ax

#%% sample winds

sample_time = datetime(2011,9,1)

with xr.open_dataset(file) as ds:
    u = ds['u10']
    v = ds['v10']
    lon = ds['longitude']
    lat = ds['latitude']
    
extent = [-70,30,50,90] #[-160,90,50,60]
fig, ax = fx.background_plot_nh(extent=extent, returnfig=True, 
                                title='Original Winds', labels=False, central_lon=0)

u1 = u.sel(time=sample_time).sel(longitude=slice(-60,0)).sel(latitude=slice(65,55))
v1 = v.sel(time=sample_time).sel(longitude=slice(-60,0)).sel(latitude=slice(65,55))

lon1 = lon.sel(longitude=slice(-60,0))
lat1 = lat.sel(latitude=slice(65,55))

lon, lat = np.meshgrid(lon1.values, lat1.values)

ax = plot_winds(ax, lon, lat, u1.values, v1.values, 
                sc=300, interval=12, headwidth=2.5, headlength=3)

#%%
#### projection?

values = []
for ux, vx in zip(u1.values.flatten(), v1.values.flatten()):
    values.append( np.dot((ux,vx), ref_vec) / np.linalg.norm(ref_vec) )

print(np.nanmean(values))

#%% end
