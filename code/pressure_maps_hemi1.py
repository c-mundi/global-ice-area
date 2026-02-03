#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 2025
pressure_maps_hemi1.py

makes min. pressure frequency maps

@author: mundi
"""
#%% imports and file paths
import numpy as np
import calendar
import time


decades = [np.arange(1982, 1991+1), np.arange(2010,2019+1)]
months = np.arange(1,12+1)

var_name = 'mean_sea_level_pressure'

min_pressures = [984, 957]

output = '/Users/mundi/Desktop/'
# output = '/home/mundi/'
output += 'month-hemi/pressure_counts/'

ice_path = '/Users/mundi/Desktop/seaice/'
ice_paths = [ice_path, ice_path+'south/']

hemi_names = ['Arctic', 'Antarctic']

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

#%%% functions

def get_era5_var(var_name, loc_ind, year, month, days):
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
        "variable": [var_name],
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

#%% data loop
get_new_data = True

p_counts = []
for loc_ind, loc in enumerate(hemi_names):
    min_p = min_pressures[loc_ind]
    
    data_dict = {mm:[[],[]] for mm in months}
    for era, years in enumerate(decades):
        for year in years:
            ### load pressure data
            for month in months:
                fname = output+str(loc_ind)+'_'+str(year)+'-'+str(month)+'_pcounts.npy'
                
                try:
                    count_array = np.load(fname)
                except FileNotFoundError:
                    if get_new_data:
                        print('Calculating '+str(year)+'-'+str(month)+' pressure map')
                        start_time = time.time()
                        
                        calendar_days = list(np.arange(1,calendar.monthrange(year, month)[1]+1))
                        days = [str(d) if d>=10 else '0'+str(d) for d in calendar_days]
                        
                        ds = get_era5_var(var_name, loc_ind, year, month, days)
                        slp = ds['msl']/100 #hpa
                        slp_daily = slp.resample(valid_time='1D').mean(dim='valid_time')
                        
                        # set up counting grid
                        count_array = np.zeros(np.shape(slp_daily.isel(valid_time=0).values))
                        
                        # loop through days
                        for ti, slp1 in enumerate(slp_daily):
                            min_inds = np.where(slp1 < min_p)
                            count_array[min_inds] += 1
                        
                        # save monthly map
                        np.save(fname, count_array)
                        print('> '+str(round((time.time()-start_time)/60,2))+' min')
                    else:
                        count_array = np.load(output+str(loc_ind)+'_1982-1_pcounts.npy')
                        print('Missing '+str(loc_ind)+'-'+str(month)+'-'+str(year))
                    
                data_dict[month][era].append(count_array)
                
    p_counts.append(data_dict)
  

#%% plot maps!
print('Plotting pressure maps:')
start_time = time.time()

import functions as fx
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean.cm as cmo
cmap_f = cmo.amp

import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

MONTH_GROUPS = [ [[12,1,2],[3,4,5],[6,7,8],[9,10,11]],
                 [[6,7,8], [9,10,11], [12,1,2],[3,4,5]]
                ]
month_groups = MONTH_GROUPS[0]

era = 0
years = decades[era]

plot_ice = True

def regrid(loc_ind, lon, lat, total_counts):
    x1 = np.arange(-180,181, 1)
    if loc_ind==0: y1 = np.arange(60,91,1)
    elif loc_ind==1: y1=np.arange(-90,-54,1)
    new_lon, new_lat = np.meshgrid(x1, y1)
    
    new_grid = np.zeros(np.shape(new_lon))
    
    for x, y, v in zip(lon.flatten(), lat.flatten(), total_counts.flatten()):
        xi = np.where(x>= x1)[0].max()
        yi = np.where(y>= y1)[0].max()
        new_grid[yi, xi]+=v

    return new_lon, new_lat, new_grid


#### set up plots
fig = plt.figure( figsize=(3*len(month_groups), 4*len(hemi_names)) )
fig.suptitle(str(years[0])+'-'+str(years[-1]), size='xx-large')
subfigs = fig.subfigures(2,1)

for loc_ind, loc in enumerate(hemi_names):
    print(' -', loc, end=': axes ')
    proj = ccrs.NorthPolarStereo() if loc_ind==0 else ccrs.SouthPolarStereo()
    
    path1 = root_paths[loc_ind]
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
    
    month_groups = MONTH_GROUPS[loc_ind]
    min_p = min_pressures[loc_ind]
    subfigs[loc_ind].suptitle('\n\n'+loc+': '+str(min_p)+' hPa', size='x-large')
    axes = subfigs[loc_ind].subplots(1, len(month_groups), subplot_kw={'projection':proj})
    
    if plot_ice: # sea ice grid
        if loc_ind==0: si, si_lon, si_lat = fx.load_seaice(ice_paths[loc_ind], 2010,8,1, latlon=True)
        elif loc_ind==1:  si, si_lon, si_lat = fx.load_seaice_sh(ice_paths[loc_ind], 2010,8,1, latlon=True)
        
    # get slp info (lon/lat from era files)
    lon = np.load(output+str(loc_ind)+'_lon.npy')
    lat = np.load(output+str(loc_ind)+'_lat.npy')
    data_dict = p_counts[loc_ind]
    
    for idx, ax in enumerate(axes):
        print(str(idx), end='')
        if loc_ind==0: ax = fx.setup_plot_nh(ax, extent=[-180,180,60,90], labels=False)
        elif loc_ind==1: ax = fx.setup_plot_sh(ax, [-180,180, -55,-90], labels=False)
        ax.set_boundary(circle, transform=ax.transAxes)

        #### load data
        total_counts = np.zeros(np.shape(lon))
        ndays = 0
        storm_counts = 0
        monthly_ice = []
        for mm in month_groups[idx]:
            storm_counts+=len(lines[mm])
            for yi, year in enumerate(years):
                total_counts+=data_dict[mm][era][yi]
                ndays += calendar.monthrange(year, mm)[-1]
                
                if plot_ice:
                    days = np.arange(1,calendar.monthrange(year,mm)[-1]+1)
                    for day in days:
                        if loc_ind==0: si = fx.load_seaice(ice_paths[loc_ind], year, mm, day, latlon=False)
                        elif loc_ind==1: si = fx.load_seaice_sh(ice_paths[loc_ind], year, mm, day, latlon=False)
                        monthly_ice.append(si)
        
        ax.set_title(str(month_groups[idx])+' n='+str(storm_counts), size='large')
        
        #### regrid
        new_lon, new_lat, new_grid = regrid(loc_ind, lon, lat, total_counts)
        pcm = ax.contourf(new_lon, new_lat, new_grid*100/ndays, 
                            transform=ccrs.PlateCarree(),
                            cmap=cmap_f, levels=np.arange(0,60,5), extend='max')
                            # vmin=0, vmax=60)
        
        #### plot map
        # pcm = ax.pcolormesh(lon, lat, total_counts*100/ndays, 
        #                     transform=ccrs.PlateCarree(),
        #                     cmap=cmap_f, vmin=0, vmax=15)
            
        #### add ice contour!
        if plot_ice: 
            ices = np.nanmean(monthly_ice, axis=0)
            fx.plot_geocontour(ax, si_lon, si_lat, ices, levels=[0.15], color='k', lw=1, ls='solid')
            fx.plot_geocontour(ax, si_lon, si_lat, ices, levels=[0.80], color='k', lw=1, ls=':')
    
        
#### colorbar
cax = fig.add_axes([0.33,0.01,0.33,.04]) 
cbar = fig.colorbar(pcm, cax=cax, orientation='horizontal')
cbar.set_label('Frequency (%)')
            
print(' > '+str(round((time.time()-start_time)/60,2))+' min')



#%% end
