#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 2025

-> background for understanding winter ice maximum

@author: mundi
"""
#%% imports and files
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import calendar
import functions as fx
import cmocean.cm as cmo
import warnings

monthly_plot = False
normalize = False
sst_level = 2 #-0.5 # -1.5
lat_width = 4 #2

ice_fname =  '/Users/mundi/Desktop/seaice/south/'

sst_path = '/Users/mundi/Desktop/data/'
sst_files = [sst_path+'sst.wkmean.1981-1989.nc',sst_path+'sst.wkmean.1990-present.nc']

root = '/Users/mundi/Desktop/month-hemi/'
root_paths = [root+'nh_data/', root+'sh_data/']

hemi_names = ['Arctic', 'Antarctic']

decades = [np.arange(2010,2020), np.arange(1982, 1992)]
months = np.arange(1,12+1)


#%%% plotting functions

def load_seaice_sh(root_dir, year, month, day, latlon=True):
    import glob
    from pyproj import Transformer
    
    if day>0: day_list = [day]
    elif day==0: day_list = np.arange(1, calendar.monthrange(year, month)[-1]+1)
    
    si_list = []
    for day in day_list:
        # convert date inputs to strings for filename if needed
        if not isinstance(year, str):
            year = str(year)
        if not isinstance(month, str):
            if month<10: month = '0'+str(int(month))
            else: month = str(int(month))
        if not isinstance(day, str):
            if day<10: day = '0'+str(int(day))
            else: day = str(int(day))
        
        # get file(s)
        all_files = glob.glob(root_dir + '*' + year+month+day + '*.nc')
        all_files.sort()
    
        if not all_files:
            print('Error with S.Hemisphere filename: ' + root_dir + '*' + year+month+day + '*.nc')
            raise NameError(' bad filename in sip.load_seaice')
            
        cdr_dic = xr.open_dataset(all_files[0])
        si = np.squeeze( cdr_dic['nsidc_bt_seaice_conc'].values )
        
        for flag in [251,252,253,254,255]:
            si= np.where(si==flag/100, np.nan, si)
        
        si_list.append( si )
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        seaice= np.nanmean(si_list, axis=0)
   
    if latlon:
        x,y = np.meshgrid( cdr_dic['xgrid'].values, cdr_dic['ygrid'].values )
        transformer = Transformer.from_crs("EPSG:3412", "EPSG:4326", always_xy=True)
        lon, lat = transformer.transform(x, y)
        return seaice, lon, lat
    else:
        return seaice

def plot_geocontour(ax, lon, lat, var, levels, color='k', lw=3, ls='solid'):
    #do masked-array on the lon
    lon_greater = np.ma.masked_greater(lon, -0.01)
    lon_lesser = np.ma.masked_less(lon, 0)    
    # apply masks to other associate arrays: lat
    lat_greater = np.ma.MaskedArray(lat, mask=lon_greater.mask)
    lat_lesser = np.ma.MaskedArray(lat, mask=lon_lesser.mask)
    # apply masks to other associate arrays: daily ice
    si_greater = np.ma.MaskedArray(var, mask=lon_greater.mask)
    si_lesser = np.ma.MaskedArray(var, mask=lon_lesser.mask)

    # contours
    c1 = ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
              linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
              linestyles=ls) 
    c2 = ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
              linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
              linestyles=ls)
    return ax, c1.allsegs + c2.allsegs
    

def background_plot(extent=[-180,180, -53,-90], returnfig=False, title=[], labels=True):
    
    fig=plt.figure(figsize=[8,8]) 
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax.coastlines('50m',edgecolor='black',linewidth=0.75, zorder=102)
    ax.set_extent(extent, ccrs.PlateCarree())
    try:
        ax.gridlines(draw_labels=labels)
    except:
        print('Unable to create grid lines on map')
    ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=100)
    ax.add_feature(cfeature.LAKES, facecolor='0.85', zorder=101)
    
    if title:
        if type(title)!= str: title=str(title)
        ax.set_title(title, fontsize=20)
    
    if not returnfig: return ax
    if returnfig: return fig, ax


#%% map
mymonths = [8,9,10]
myyears = np.arange(2010,2020)

line_out = {mm:[[],[]] for mm in months} 
for year in myyears:
    print('*', year)
    for month in mymonths:
        ds = xr.open_dataset(sst_files[1]).sel(lat=slice(-50,-90)).sel(time=slice(datetime(year, month,1),datetime(year, month+1,1)))
        lon, lat = np.meshgrid(ds['lon'].values, ds['lat'].values)
        mean_sst = ds['sst'].mean(dim='time').values
        
        fig1, ax = background_plot(extent=[-180,180, -53,-90], returnfig=True, title=str(year)+'-'+str(month), labels=True)
        ax.contourf(lon, lat, mean_sst, transform=ccrs.PlateCarree(), levels=np.arange(-2,10,1), alpha=0.75)
        ax, contours = plot_geocontour(ax, lon, lat, mean_sst, levels=[sst_level], color='k', lw=3, ls='solid')
        
        sic, si_lon, si_lat =  load_seaice_sh(ice_fname, year, month, 0)
        ax, contours_si = plot_geocontour(ax, si_lon, si_lat, sic, [0.15], color='navy', lw=3, ls='dashed')
        if not monthly_plot: plt.close(fig1)
        
        #### fraction of winter storms that occur within 2-ish degree contour +/- a few degrees lat...
        
        full_clist = []
        for clist in contours:
            for cont in clist:
                full_clist += list(cont)
        sorted_clist = sorted(full_clist, key=lambda x: (x[0],x[1]))
        
        cpoints1 = []
        for i, x in enumerate(sorted_clist[1:]):
            if x[0]==sorted_clist[i-1][0]: 
                dupes = np.array(sorted_clist)[np.where(np.array([x[0] for x in sorted_clist])==sorted_clist[i-1][0])]
                cpoints1.append(dupes[np.argmax(dupes[:,1])])
            else:
                cpoints1.append(sorted_clist[i-1])
        cpoints2 = []
        for i, x in enumerate(cpoints1[1:]):
            if x[0]==sorted_clist[i-1][0]: 
                dupes = np.array(cpoints1)[np.where(np.array([x[0] for x in cpoints1])==cpoints1[i-1][0])]
                cpoints2.append(dupes[np.argmax(dupes[:,1])])
            else:
                cpoints2.append(sorted_clist[i-1])
        cpoints = np.array(cpoints2)
            
        if monthly_plot:
            ax.scatter(cpoints[:,0], cpoints[:,1]+2, transform=ccrs.PlateCarree(), color='cyan', s=4)
            ax.scatter(cpoints[:,0], cpoints[:,1]-2, transform=ccrs.PlateCarree(), color='magenta', s=4)
        
        
        # count=0
        # for i, x in enumerate(cpoints[1:]):
        #     if x[0]==cpoints[i-1][0]: 
        #         if x[1] != cpoints[i-1][1]:
        #             count+=1
        #             print(x, cpoints[i-1])
        # print(count, '/', len(cpoints))
        
        #### storm lists
            
        start_day, end_day, longitudes, latitudes = fx.indiv_locations([year], root_paths[1]+'census/', root_paths[1]+'area/')
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines([year], root_paths[1]+'census/', root_paths[1]+'area/', root_paths[1]+'seaice/')
            
        for i, start_x, start_y, impact, si, clim in zip(np.arange(0,len(start_day[month])), 
                                                longitudes[month], latitudes[month], lines[month], 
                                                si_changes[month], clim_changes[month]):
        
            # find closest lat/lon
            closest_lon = min(cpoints[:,0], key=lambda x:abs(x-start_x))
            closest_lat = [cp[1] for cp in cpoints if cp[0]==closest_lon][0]
            
            cmap = cmo.balance_r 
            norm = plt.Normalize(vmin=-1, vmax=1)
            color = cmap(norm(impact[-1]))
            
            if (start_y > closest_lat-lat_width) and (start_y < closest_lat+lat_width):
                # print('Within bounds!')
                marker = '*'
                ec='y'
                if normalize: line_out[month][0].append(impact)
                else: 
                    ss = si-clim
                    line_out[month][0].append( ss-ss[0] )
            else:
                marker='.'
                ec='r'
                if normalize: line_out[month][1].append(impact)
                else:
                    ss = si-clim
                    line_out[month][1].append( ss-ss[0] )
              
            if monthly_plot:
                ax.scatter(start_x, start_y, transform=ccrs.PlateCarree(), 
                           color=color, s=500, marker=marker, zorder=9999) #edgecolors=ec, 
                

#%%% spagetti plot
xxx = np.arange(-7,14+1,1)

fig, axes = plt.subplots(2, len(mymonths), figsize=(5*len(mymonths), 10), sharey=True)
fig.suptitle(str(myyears[0])+'-'+str(myyears[-1]))

for mi, mm in enumerate(mymonths):
    
    for idx, line_type in enumerate(['ACC','Interior']):
        axes[idx][mi].set_title(calendar.month_name[mm]+': '+line_type+'  ('+str(len(line_out[mm][idx]))+')')
        axes[idx][mi].axhline(0, lw=0.75, color='k')
        axes[idx][mi].axvline(0, lw=0.75, color='k')
        
        axes[idx][mi].plot(xxx, np.array(line_out[mm][idx]).T, 
                           lw=0.55, color='gray')
        axes[idx][mi].plot(xxx, np.nanmean(line_out[mm][idx],axis=0), 
                           lw=1.25, color='maroon')
        
#### mean/std comparison   
fig, axes2 = plt.subplots(1, len(mymonths), figsize=(5*len(mymonths), 5), sharey=True)
fig.suptitle(str(myyears[0])+'-'+str(myyears[-1]))

line_colors = ['C0','C1']

for mi, mm in enumerate(mymonths):
    axes2[mi].set_title(calendar.month_name[mm])
    axes2[mi].axhline(0, lw=0.75, color='k')
    axes2[mi].axvline(0, lw=0.75, color='k')
    
    for idx, line_type in enumerate(['ACC','Interior']):
        
        mean_line = np.nanmean(line_out[mm][idx],axis=0)
        std_line = np.nanstd(line_out[mm][idx],axis=0)
                
        fx.plot_spread(axes2[mi], xxx, mean_line, std_line, ls1='-', ls2='--', 
                       color=line_colors[idx])
        
        axes2[mi].plot([],[], color=line_colors[idx], label=line_type+', n='+str(len(line_out[mm][idx])))
        
    axes2[mi].legend(loc='upper left')

#%% end
