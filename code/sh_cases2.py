#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 2025
sh_cases2.py

clean finalized case study figures

@author: mundi
"""
#%% imports and files
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.gridspec import GridSpec
import time as timeIN
from glob import glob

import functions as fx

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

census_name = 'census_'
contour_name = '_contours'
si_name = '_seaice'

ice_fname = '/Users/mundi/Desktop/seaice/south/'

decades = [np.arange(1982, 1992), np.arange(2010,2020)]
decade_names = ['Early Satellite Era ', 'Present Day ']

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]

hemi_names= ['Arctic', 'Antarctic']

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

fontsize=11

loc_ind=1
path1= root_paths[loc_ind]

#%%% open data files
_, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010, 1, 1)

def load_era_data(case_ind):
    era_data_file = '/Users/mundi/Desktop/era_data/antarctic_2012_10_11_msl_wind.nc'
    if case_ind==1:
        with xr.open_dataset(era_data_file) as ds:
            lon_name = 'longitude' 
            ds['_longitude_adjusted'] = xr.where(ds[lon_name] < 0,ds[lon_name] + 360,ds[lon_name])
            ds = (ds.swap_dims({lon_name: '_longitude_adjusted'})
                    .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
                    .drop_vars(lon_name))
            ds = ds.rename({'_longitude_adjusted': lon_name})
            
            slp1 = ds['msl']/100 #hpa
            slp_lon = slp1.sel(longitude=slice(60,300))
            slp = slp_lon.resample(valid_time='1D').mean(dim='valid_time')
            lon, lat = np.meshgrid(slp['longitude'].values, slp['latitude'].values)
    elif case_ind==0:
        with xr.open_dataset(era_data_file) as ds:
            slp = ds['msl']/100 #hpa
            slp_lon = slp.sel(longitude=slice(-120,0))
            slp = slp_lon.resample(valid_time='1D').mean(dim='valid_time')
            lon, lat = np.meshgrid(slp['longitude'].values, slp['latitude'].values)
    return lon, lat, slp

def get_winds(storm_range):
    era_data_file = '/Users/mundi/Desktop/era_data/antarctic_2012_10_11_msl_wind.nc'
    with xr.open_dataset(era_data_file) as ds:
        ds = ds.resample(valid_time='1D').mean(dim='valid_time')
        u10 = ds['u10'].sel(valid_time=slice(storm_range[0], storm_range[-1]))
        v10 = ds['v10'].sel(valid_time=slice(storm_range[0], storm_range[-1]))
        x10, y10 = np.meshgrid(ds['longitude'].values, ds['latitude'].values)
    return x10, y10, u10, v10

#%%% functions

def get_all_coords(storm_range, slp, lon, lat, plevel=970):
    import matplotlib.path as mpltPath
    slp_slice = slp.sel(valid_time=slice(storm_range[0], storm_range[-1]))
    
    all_coords = []
    for p_daily in slp_slice:
        fig, ax = fx.background_plot_sh(returnfig=True)
        contours = ax.contour(lon, lat, p_daily, levels=[plevel],
                              colors='k', transform=ccrs.PlateCarree())
        plt.close(fig)
        
        # min pressure
        indp = np.where(p_daily == np.nanmin(p_daily))
        minlon = np.squeeze(lon[indp])
        minlat = np.squeeze(lat[indp])
        points = np.array([ [minlon, minlat] ])
        
        cont = contours.allsegs[0]
        for cc in cont:
           if len(cc) > 2 : path = mpltPath.Path(cc)
           else: continue
           pts_inside = path.contains_points(points) # radius=30
           
           if any(pts_inside):
               all_coords.append( cc )
               
    return all_coords

def get_bbox(all_coords):
    minlat = np.nanmin([np.nanmin(cc[:,1]) for cc in all_coords])
    maxlat = np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords])
    minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
    maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
    return fx.get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False)

def geoplot_bbox(ax, bbox, color='k'):   
    x = bbox[:,0]
    y = bbox[:,1]
    #do masked-array on the lon
    x_greater = np.ma.masked_greater(x, -0.01)
    x_lesser = np.ma.masked_less(x, 0)    
    # apply masks to other associate arrays: lat
    y_greater = np.ma.MaskedArray(y, mask=x_greater.mask)
    y_lesser = np.ma.MaskedArray(y, mask=x_lesser.mask)
    
    # plot
    for x1, y1 in zip([x_greater, x_lesser], [y_greater, y_lesser]):
        ax.plot(x1, y1, color=color, lw=2, transform = ccrs.PlateCarree())
        
def get_miz_area(storm_event):
    miz = [0.15, 0.80]

    t1 = storm_event[0] - timedelta(days=1)
    t2 = storm_event[-1] + timedelta(days=1)
    storm_range = fx.daterange(t1, t2, dt=24)
    
    miz_points = np.zeros(np.shape(si_lon))
    for date in storm_range:
        sic = fx.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
   
    return miz_points


#%% simultaneous storm plots
print(); print('Plotting simultaneous maps and time series'); print()

storm_range1 = fx.daterange(datetime(2012, 10, 25), datetime(2012, 10, 29), dt=24)
storm_range2 = fx.daterange(datetime(2012, 10, 29), datetime(2012, 10, 31), dt=24)

clim_years = np.arange(2010,2020)

storm_colors = ['#e9a3c9','#a1d76a']

#### set up plots
fig = plt.figure(figsize=(12, 10))
fig.suptitle('Simultaneous Storms', size='xx-large', fontweight='bold')
subfigs = fig.subfigures(1, 2)
subfigs[0].suptitle('\nIce-Decreasing Storm', size='x-large')
subfigs[1].suptitle('\nIce-Increasing Storm', size='x-large')

for case_ind, storm_range in enumerate([storm_range1, storm_range2]):
    gs = GridSpec(nrows=3, ncols=1, height_ratios=[0.95,0.33, 0.33])
    dt_s = storm_range[0]
    dt_e = storm_range[-1]
    
    # set up data
    lon, lat, slp = load_era_data(case_ind)
    
    ##########
    #### map
    ##########
    ax = subfigs[case_ind].add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
    ax = fx.setup_plot_sh(ax, labels=False)
    ax.set_title(dt_s.strftime('%b %d')+'-'+dt_e.strftime('%d %Y'))
    
    # plot daily contours
    all_coords = get_all_coords(storm_range, slp, lon, lat, plevel=970)
    for cont in all_coords:
        ax.plot(cont[:,0], cont[:,1], color='k',
                transform = ccrs.PlateCarree())
        
    # plot bbox
    bbox = get_bbox(all_coords)
    inside = fx.find_points_in_contour(bbox, si_lon, si_lat) ### inside bbox!!(~)
    if case_ind==0: inside = ~inside
    geoplot_bbox(ax, bbox, color='magenta' if case_ind==0 else 'green')

    # plot change in sea ice
    si_start = fx.load_seaice_sh(ice_fname, dt_s.year, dt_s.month, dt_s.day, latlon=False)
    si_end = fx.load_seaice_sh(ice_fname, dt_e.year, dt_e.month, dt_e.day, latlon=False)
    diff = si_end-si_start
    diff = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), diff)
    ax.pcolormesh(si_lon, si_lat, np.ma.masked_array(diff, mask=inside), 
                  cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, transform=ccrs.PlateCarree() )
    
    ################
    #### timeseries
    ################
    full_time = fx.daterange(dt_s-timedelta(days=7), dt_s+timedelta(days=14), dt=24)
    
    # set up plot
    ax_t = subfigs[case_ind].add_subplot(gs[1])
    ax_t.axvline(0, lw=0.55, color='gray', ls=':')
    ax_t.set_xlim(-7,14)
    ax_t.set_xticks(xxx)
    ax_t.tick_params(axis='both', which='major')  
    ax_t.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
    
    # storm ranges
    axt2 = ax_t.twiny()
    axt2.set_xlim([full_time[0].day, full_time[-1].day+31])
    axt2.set_xticks([ft.day if ft.month==full_time[0].month else ft.day+31 for ft in full_time])
    axt2.set_xticklabels([ft.day for ft in full_time])
    axt2.axvspan(storm_range1[0].day, storm_range1[-1].day+1, color=storm_colors[0], alpha=0.33)
    axt2.axvspan(storm_range2[0].day, storm_range2[-1].day+1, color=storm_colors[1], alpha=0.33)
    
    # calculate lines
    miz_points = get_miz_area(storm_range)
    tseries = []
    clim = []
    for dt in full_time:
        si, si_lon, si_lat = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day)
        si_bbox = np.ma.masked_array(si, mask=inside)
        si_miz = np.ma.masked_where(miz_points==0, si_bbox).filled(np.nan)
        tseries.append(np.nansum(si_miz*25*25))
        daily = []
        for yy in clim_years:
            si = fx.load_seaice_sh(ice_fname, yy, dt.month, dt.day, latlon=False)
            si = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), si)
            si_bbox = np.ma.masked_array(si, mask=inside)
            si_miz = np.ma.masked_where(miz_points==0, si_bbox).filled(np.nan)
            daily.append(np.nansum(si_miz*25*25))
        clim.append(np.nanmean(daily))
    tseries = np.array(tseries); clim=np.array(clim)
    
    # plot timeseries
    ax_t.plot(xxx, tseries/1e6, color='k')
    ax_t.plot(xxx, clim/1e6, color='dimgray', ls='--')
    if case_ind==0: ax_t.set_ylabel(r'$\Delta$ MIZ Ice Area'+'\n'+r'($\times10^6$ km$^2$)')
    

    # normalized
    ax_n = subfigs[case_ind].add_subplot(gs[2])
    ax_n.axvline(0, lw=0.55, color='gray', ls=':')
    ax_n.axhline(0, lw=0.55, color='gray', ls=':')
    ax_n.set_xlim(-7,14)    
    ax_n.set_xticks(xxx)
    ax_n.tick_params(axis='both', which='major')  
    ax_n.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
    ax_n.set_xlabel('Days Since Storm Start', fontsize=fontsize-1)
    if case_ind==0: ax_n.set_ylabel('Normalized\n'+r'$\Delta$ MIZ Ice Area')
    
    
    rel = tseries-clim
    normalized = (rel-rel[0])/(np.max(rel)-np.min(rel))
    ax_n.plot(xxx, normalized, color='k')
    
    ax2 = ax_n.twiny()
    ax2.set_xlim([full_time[0].day, full_time[-1].day+31])
    ax2.axvline(storm_range1[0].day, color=storm_colors[0], alpha=0.5)
    ax2.axvline(storm_range1[-1].day+1, color=storm_colors[0], alpha=0.5)
    ax2.axvline(storm_range2[0].day, color=storm_colors[1], alpha=0.5)
    ax2.axvline(storm_range2[-1].day+1, color=storm_colors[1], alpha=0.5)

    ax2.get_xaxis().set_visible(False)
    ax2.spines[['right', 'top','bottom','left']].set_visible(False)

#%%% simultaneous daily differences
start_plot_time = timeIN.time()

print(); print('Plotting daily sea ice differences')
fig = plt.figure(figsize=(14,9))
fig.suptitle('Simultaneous Storms', size='xx-large', fontweight='bold')
subfigs = fig.subfigures(2, 1)
subfigs[0].suptitle('\n\nIce-Decreasing Storm', size='x-large')
subfigs[1].suptitle('Ice-Increasing Storm', size='x-large')

for case_ind, storm_range in enumerate([storm_range1, storm_range2]):
    print(); print('- Case', case_ind, end=': ')
    lon, lat, slp = load_era_data(case_ind)
    all_coords = get_all_coords(storm_range, slp, lon, lat, plevel=970)
    minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
    maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
    bbox = get_bbox(all_coords)
    inside = fx.find_points_in_contour(bbox, si_lon, si_lat)
    if case_ind==0: inside = ~inside
   
    miz_points = get_miz_area(storm_range)
    
    if len(storm_range)>3: storm_range = storm_range[::2]
    gs = GridSpec(nrows=1, ncols=len(storm_range))
    
    for di, dt in enumerate(storm_range):
        print(di, end=' ')
        # setup plot
        ax = subfigs[case_ind].add_subplot(gs[di], projection=ccrs.SouthPolarStereo(central_longitude=np.mean([minlon,maxlon])))
        ax = fx.setup_plot_sh(ax, extent = [minlon-2, maxlon+2, -85,-55], labels=False)
        geoplot_bbox(ax, bbox, color='k')

        # plot change in sea ice
        si1 = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day, latlon=False)
        dt2 = dt+timedelta(days=1)
        si2 = fx.load_seaice_sh(ice_fname, dt2.year, dt2.month, dt2.day, latlon=False)
        diff = si2-si1
        
        diff = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), diff)
        diff_bbox = np.ma.masked_array(diff, mask=inside)
        diff_miz = np.ma.masked_where(miz_points==0, diff_bbox).filled(np.nan)
        
        ax.set_title(dt2.strftime('%d-')+dt.strftime('%d %b'))
        pcm = ax.pcolormesh(si_lon, si_lat, np.ma.masked_array(diff_miz, mask=inside), 
                            cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, transform=ccrs.PlateCarree() )

        # add winds on top
        x10, y10, u10, v10 = get_winds(storm_range)
        inside_w = fx.find_points_in_contour(bbox, x10, y10)
        if case_ind==0: inside_w = ~inside_w
        ax, Q = fx.plot_winds(x10, y10, u10.sel(valid_time=dt), v10.sel(valid_time=dt), 
                              ax, mask=inside_w, 
                              sc= 270 if case_ind==0 else 300, 
                              interval = 14 if case_ind==0 else 10)
        
    # colorbar / wind scale
    qk = ax.quiverkey(Q, 0., 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')
    
# color bar
cax1 = fig.add_axes([0.33,0.0,0.33,0.025]) 
cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
cbar1.set_label(r'Change in Sea Ice Concentration', fontsize=10)
cax1.tick_params(labelsize=10)

timed = timeIN.time() - start_plot_time
print(); print(str(round(timed/60, 1))+' min'); print()


#%% may case study
print(); print('Plotting May Case Study'); print()

#### storm info
yy = 2019
mm = 5
dd = 26
dd2 = 29

dt5 = datetime(yy,mm,dd)

mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
    fx.indiv_lines([yy], path1+'census/', path1+'area/', path1+'seaice/')
    
full_time = fx.daterange(dt5-timedelta(days=7), dt5+timedelta(days=14), dt=24)
storm_range = fx.daterange(dt5, datetime(yy,mm,dd2), dt=24)

##########
#### PLOT
##########

fig = plt.figure(figsize=(7.5,10)) 
gs = GridSpec(nrows=3, ncols=1, height_ratios=[0.95,0.33,0.33])
# maps
ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
# timeseries
axes2 = fig.add_subplot(gs[1])
axes3 = fig.add_subplot(gs[2])

#### map  
ax = fx.setup_plot_sh(ax, labels=False, extent=[-180,180,-90,-58])

# storm area
for start,end, si, clim, line in zip(start_day[mm], end_day[mm], 
                                     si_changes[mm], clim_changes[mm], lines[mm]):
    if start.day==dd: break
ax.set_title(start.strftime('%Y %b %d')+'-'+end.strftime('%d'))

# contours, location of min pressure
stormstr1 = start.strftime('%Y_%m%d')
ncname = stormstr1+'_contours.nc'
with xr.open_dataset(path1+'contours/'+ncname) as cs:
    all_coords=[]
    for key in list(cs.keys()):
        coord = cs[key].values
        all_coords.append(coord)
        ax.plot(coord[:,0], coord[:,1], color='k', lw=0.75,
                transform=ccrs.PlateCarree())
bbox = get_bbox(all_coords)
geoplot_bbox(ax, bbox, color='k')
inside = ~fx.find_points_in_contour(bbox, si_lon, si_lat)

# change in sea ice
si_start, si_lon, si_lat = fx.load_seaice_sh(ice_fname, yy, mm, dd, latlon=True)
si_end = fx.load_seaice_sh(ice_fname, yy, mm, end.day, latlon=False)

change = np.where(si_end-si_start==0, np.nan, si_end-si_start)
change = np.ma.masked_where(si_lat<np.nanmin([np.nanmin(cc[:,1]) for cc in all_coords]), change)
change = np.ma.masked_array(change, mask=inside)

pcm = ax.pcolormesh(si_lon, si_lat, change,
                    cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, zorder=-5,
                    transform=ccrs.PlateCarree())

#### timeseries/normalization procedure 
axes2.plot(xxx, si/1e6, color='k')
axes2.plot(xxx, clim/1e6, color='k', ls='--')

axes2.axhline(0, ls='-', color='k', lw=1)
axes2.axvline(0, ls=':', color='gray', lw=1)
axes2.set_xlim(-7,14)
axes2.set_xticks(xxx)
axes2.set_xticklabels(xlabels, minor=False)
axes2.set_ylabel(r'$\Delta$ MIZ Ice Area'+'\n'+r'($\times10^6$ km$^2$)')

axes3.plot(xxx, line, color='k', lw=3)

axes3.axhline(0, ls='-', color='k', lw=1)
axes3.axvline(0, ls=':', color='gray', lw=1)
axes3.set_xlim(-7,14)
axes3.set_xticks(xxx)
axes3.set_xticklabels(xlabels, minor=False)
axes3.set_ylabel('Normalized\n'+r'$\Delta$ MIZ Ice Area')


# color bar
cax1 = fig.add_axes([0.125,0.55,0.025,0.33]) 
cbar1 = fig.colorbar(pcm, cax=cax1, orientation='vertical')
cax1.set_xlabel('\nChange in\nSea Ice\nConcentration')
cax1.tick_params(labelsize=10)


#%%% daily sea ice changes
print(); print('Plotting daily MAY sea ice differences')

fig, axes = plt.subplots(2, int(np.ceil(len(storm_range)/2)), figsize=(7,5), 
                         subplot_kw={'projection':ccrs.SouthPolarStereo()})
fig.suptitle(dt5.strftime('%b %Y'), size='xx-large', fontweight='bold')
    
for di, dt in enumerate(storm_range):
    print(di, end=' ')
    # setup plot
    ax = axes.flatten()[di]
    minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
    maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
    ax = fx.setup_plot_sh(ax, extent = [minlon-2, maxlon+2, -85,-55], labels=False)
    geoplot_bbox(ax, bbox, color='k')

    # plot change in sea ice
    si1 = fx.load_seaice_sh(ice_fname, dt.year, dt.month, dt.day, latlon=False)
    dt2 = dt+timedelta(days=1)
    si2 = fx.load_seaice_sh(ice_fname, dt2.year, dt2.month, dt2.day, latlon=False)
    diff = si2-si1
    
    diff = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), diff)
    diff_bbox = np.ma.masked_array(diff, mask=inside)
    diff_miz = np.ma.masked_where(miz_points==0, diff_bbox).filled(np.nan)
    
    ax.set_title(dt2.strftime('%d-')+dt.strftime('%d %b'))
    pcm = ax.pcolormesh(si_lon, si_lat, np.ma.masked_array(diff_miz, mask=inside), 
                        cmap=cmo.balance_r, vmin=-0.3, vmax=0.3, transform=ccrs.PlateCarree() )

# color bar
cax1 = fig.add_axes([0.33,0.0,0.33,0.025]) 
cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
cbar1.set_label(r'Change in Sea Ice Concentration', fontsize=10)
cax1.tick_params(labelsize=10)

timed = timeIN.time() - start_plot_time


#%%% net sea ice change
print(); print('Plotting NET MAY sea ice differences')

fig, axes = plt.subplots(1, 2, figsize=(8,3.5), 
                         subplot_kw={'projection':ccrs.SouthPolarStereo()})
fig.suptitle(dt5.strftime('%b %Y'), size='xx-large', fontweight='bold')


dt1 = storm_range[0]
dt2 = storm_range[-1]

# dt1 = full_time[7]
# dt2 = full_time[-1]

# intital
ax = axes[0]
minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
ax = fx.setup_plot_sh(ax, extent = [minlon-2, maxlon+2, -85,-55], labels=False)
geoplot_bbox(ax, bbox, color='k')

si1 = fx.load_seaice_sh(ice_fname, dt1.year, dt1.month, dt1.day, latlon=False)
si1 = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), si1)
si1_bbox = np.ma.masked_array(si1, mask=inside)
si1_miz = np.ma.masked_where(miz_points==0, si1_bbox).filled(np.nan)

ax.set_title(dt1.strftime('%d'))
pcm = ax.pcolormesh(si_lon, si_lat, si1_miz, 
                    cmap=cmo.ice, vmin=0, vmax=1, 
                    transform=ccrs.PlateCarree() )


# end
ax = axes[1]
minlon = np.nanmin([np.nanmin(cc[:,0]) for cc in all_coords])
maxlon = np.nanmax([np.nanmax(cc[:,0]) for cc in all_coords])
ax = fx.setup_plot_sh(ax, extent = [minlon-2, maxlon+2, -85,-55], labels=False)
geoplot_bbox(ax, bbox, color='k')

si2 = fx.load_seaice_sh(ice_fname, dt2.year, dt2.month, dt2.day, latlon=False)
si2 = np.ma.masked_where(si_lat>np.nanmax([np.nanmax(cc[:,1]) for cc in all_coords]), si2)
si2_bbox = np.ma.masked_array(si2, mask=inside)
si2_miz = np.ma.masked_where(miz_points==0, si2_bbox).filled(np.nan)

ax.set_title(dt2.strftime('%d'))
pcm = ax.pcolormesh(si_lon, si_lat, si2_miz, 
                    cmap=cmo.ice, vmin=0, vmax=1, 
                    transform=ccrs.PlateCarree() )

# for ax in [axes[0], axes[1]]:
#     ax.contour(si_lon, si_lat, si2_miz, levels=[0.15],
#                transform=ccrs.PlateCarree(), colors='m')


#%% end
