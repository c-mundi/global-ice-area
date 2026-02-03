#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 2025
linear_winds1.py

- like linear_miz, but with winds
* analyze east-west changes in winds over the MIZ

@author: mundi
"""

#%% imports, file names
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import time as timeIN
import calendar
import functions as fx
from datetime import datetime, timedelta
from scipy.interpolate import griddata
import warnings

REMOTE=False

if not REMOTE: 
    ROOT = '/Users/mundi/Desktop/'
    ice_path = ROOT +'seaice/'
    ice_paths = [ice_path, ice_path+'south/']
else: 
    ROOT = '/home/mundi/'
    ice_root = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/'
    ice_paths = [ice_root+ice_add for ice_add in ['daily/','south/']]

nh_path = ROOT +'month-hemi/nh_data/'
sh_path = ROOT +'month-hemi/sh_data/'
root_paths = [nh_path, sh_path]


savepath = ROOT+'month-hemi/linear_winds/'
savepath_miz = '/Users/mundi/Desktop/month-hemi/linear_miz/'

census_name = 'census_'
contour_name = '_contours.nc'
si_name = '_seaice'

#%%% controls

# years = np.arange(2010,2020)
years = np.arange(1982, 1992)
# years = np.arange(2010,2015)

months = np.arange(1,12+1)

calc_time = 0 # 0 = storm duration, 1 = 7 days, 2 = 3 weeks

xxx = np.arange(-7,14+1,1)
ice_lims = [20,80]
miz = [0.15,0.80]

lon_num = 100
lon_x = np.linspace(0, 1, lon_num)

#%%% functions

def get_era5_wind(year, months, days, LOC):
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
    if type(months) != list and type(months[0])!=str:
        months = ['0'+str(int(month)) if month<10 else str(int(month)) for month in months]
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
        'month':months,
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
    starting_day = datetime(int(year[0]), int(months[0]), int(days[0]))
    time = np.array([starting_day + timedelta(hours=int(hr)) for hr in hours_in])
    
    ds.coords['valid_time'] = time
    ds = ds.resample(valid_time='1D').mean()
    
    return ds


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


#%% data loop
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    print(); print('***', loc, '***')
    root_path = root_paths[loc_ind]
    ice_fname = ice_paths[loc_ind]
    
    if loc_ind == 0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010,1,1, latlon=True)
    elif loc_ind ==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010,1,1, latlon=True)
    si_lon1 = np.where(si_lon<0, si_lon+360, si_lon)
    
    for month in months:
        print()
        print('- '+str(month)+' -')
        month_group_timer = timeIN.time()
    
        for year in years:
            print(year)
            savename = loc+'_'+str(calc_time)+'_'+str(year)+'_'+str(month)+'.npy'
            
            try:
                lines = np.load(savepath+savename)
            except FileNotFoundError:
                ## get wind dataset
                daystrs = ['0'+str(int(day)) if day<10 else str(int(day)) for day in np.arange(1,calendar.monthrange(year, month)[-1]+1)]
                monstrs = ['0'+str(int(month)) if month<10 else str(int(month)) for month in [month, month+1]]
                ds = get_era5_wind([str(year)], monstrs, daystrs, loc_ind)
                
                lon, lat = np.meshgrid(np.array(ds['longitude'].values), np.array(ds['latitude'].values))
                try:
                    time_start = timeIN.time()
                    
                    # load census and sea ice info
                    census_name = 'census_'+ str(year) +'.csv'
                    [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
                        fx.readCensus(root_path+'census/'+census_name, convertDT=True)
                        
                    
                    #### loop thru storms
                    lines = []
                    for sn, start in enumerate(startdate):
                        if start.month != month: continue
                        end = enddate[sn]
                        
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
                        miz_points = np.zeros(np.shape(si_lon))
                        for date in storm_range:
                            if loc_ind==0: sic = fx.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
                            elif loc_ind==1: sic = fx.load_seaice_sh(ice_fname, date.year, date.month, date.day, latlon=False)
                            miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
                        miz_grd = griddata((si_lon.flatten(), si_lat.flatten()), miz_points.flatten(),
                                                    (lon.flatten(), lat.flatten()), method='nearest')
                        miz_grd = miz_grd.reshape(np.shape(lon))
                    
                        # calculate mean winds
                        dt1,dt2 = timediff(calc_time, start, end)
                        u_wind = ds['u10'].sel(valid_time=slice(dt1,dt2)).mean(dim='valid_time').values
                        v_wind = ds['v10'].sel(valid_time=slice(dt1,dt2)).mean(dim='valid_time').values
                    
                        # isolate miz
                        inside_bbox = fx.find_points_in_contour(bbox_edges, lon, lat)
                        u_bbox = np.ma.masked_array(u_wind, mask=inside_bbox)
                        u_miz = np.ma.masked_where(miz_grd==0, u_bbox).filled(np.nan)
                        v_bbox = np.ma.masked_array(v_wind, mask=inside_bbox)
                        v_miz = np.ma.masked_where(miz_grd==0, v_bbox).filled(np.nan)
                        
                        # set up grid
                        lon_bnds = np.linspace( np.nanmin(bbox_edges[:,0]), np.nanmax(bbox_edges[:,0]), lon_num)
                        lat_bnds = np.linspace( np.nanmin(bbox_edges[:,1]), np.nanmax(bbox_edges[:,1]), lon_num)
                        new_x, new_y = np.meshgrid(lon_bnds, lat_bnds)   
                    
                    
                        #### calculate final values
                        with warnings.catch_warnings(): 
                            warnings.simplefilter("ignore", category=RuntimeWarning)
                            
                            u_gridded = griddata((lon.flatten(), lat.flatten()), u_miz.flatten(),
                                                   (new_x.flatten(), new_y.flatten()))
                            u_gridded = np.reshape(u_gridded, np.shape(new_x))
                            ULINE = np.nanmean(u_gridded, axis=0)
                            
                            v_gridded = griddata((lon.flatten(), lat.flatten()), v_miz.flatten(),
                                                   (new_x.flatten(), new_y.flatten()))
                            v_gridded = np.reshape(v_gridded, np.shape(new_x))
                            VLINE = np.nanmean(v_gridded, axis=0)
                    
                    
                        lines.append( [ULINE, VLINE] )
                        
                    # export data
                    np.save(savepath+savename, np.array(lines))
                    print('Loaded Data: '+str(round((timeIN.time()-time_start)/60,1))+' min')
        
                except Exception as ee:
                    print(); print(); print(ee); print(year, month); print();print()
                    continue
                        
            
        print('--> '+ str(round((timeIN.time()-month_group_timer)/60,1))+' min')



#%% data organization

def get_miz_lines(loc, year, month, savepath, root_path):
    savename = loc+'_'+str(calc_time)+'_'+str(year)+'_'+str(month)+'.npy'
    wind_series = np.load(savepath+savename)
    wind_series = iter(wind_series)
    
    # load census and sea ice info
    census_name = 'census_'+ str(year) +'.csv'
    [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
        fx.readCensus(root_path+'census/'+census_name, convertDT=True)
        
    # open ice area
    with xr.open_dataset(root_path+'area/'+str(year)+'_area.nc') as ds_area:
        ice_area80 = ds_area['ice_area80'].values
        # ice_area15 = ds_area['ice_area15'].values
        box_area = ds_area['box_area'].values
        ice_sorter = ice_area80
    
    #### loop thru storms
    lines = {'u':[],'v':[]}
    for sn, start in enumerate(startdate):
        if start.month != month: continue
        winds = next(wind_series)
        
        ice_frac = ice_sorter[sn]*100/box_area[sn]             
        if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
            continue
        
        lines['u'].append(winds[0])
        lines['v'].append(winds[0])
    
    return lines


def get_miz_lines_seaice(loc, year, month, savepath, root_path):
    add = '' # all_
    savename = loc+'_miz_line_'+add+str(year)+'_'+str(month)+'.npy'
    miz_series = np.load(savepath+savename)
    miz_series = iter(miz_series)
    
    # load census and sea ice info
    census_name = 'census_'+ str(year) +'.csv'
    [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
        fx.readCensus(root_path+'census/'+census_name, convertDT=True)
        
    # open ice area
    with xr.open_dataset(root_path+'area/'+str(year)+'_area.nc') as ds_area:
        ice_area80 = ds_area['ice_area80'].values
        # ice_area15 = ds_area['ice_area15'].values
        box_area = ds_area['box_area'].values
        ice_sorter = ice_area80
    
    #### loop thru storms
    lines = []
    for sn, start in enumerate(startdate):
        if start.month != month: continue
        line = next(miz_series)
        
        ice_frac = ice_sorter[sn]*100/box_area[sn]             
        if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
            continue
        
        lines.append(line)
    
    return lines
        

#%% PLOT

years = np.arange(2010,2020)
# years = np.arange(1982, 1992)
ystr = str(years[0]) +'-'+ str(years[-1])

grouped_months = [[2,3,4],[5,6,7],[8,9,10],[11,12,1]]

wind_idx = 1
uv = ['u','v'][wind_idx]
wind_name = ['Zonal', 'Meridional'][wind_idx]

#%%% plot lines
# set up plot
nrows = len(grouped_months)+1
fig, axes = plt.subplots(nrows, 2, figsize=(12,3.5*nrows), sharey=True,
                         gridspec_kw={"height_ratios":[0.001]+([1]*(nrows-1))})
for ax1 in axes.flatten():
    ax1.set_xlim([lon_x[0], lon_x[-1]])
    ax1.set_xticks(np.arange(lon_x[0], lon_x[-1]+0.25, 0.25))
    ax1.axhline(0, lw=0.55, color='gray', ls=':')
    ax1.axvline(0.5, lw=0.55, color='gray', ls=':')
for ax1 in axes[-1][:]: ax1.set_xlabel('Longitude Extent Fraction')
for ax1 in axes[:,0]: ax1.set_ylabel(wind_name + ' Wind')

# organize data
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    
    root_path = root_paths[loc_ind]
    axes[0][loc_ind].axis("off")
    axes[0][loc_ind].set_title('\n'+loc, fontweight='bold')
    
    
    for mx, months in enumerate(grouped_months):
        axes[mx+1][loc_ind].set_title(ystr+'  '+str(months))
        data = []
        for month in months:
            for year in years:
                try:
                    lines = get_miz_lines(loc, year, month, savepath, root_path)
                    data += lines[uv]
                except FileNotFoundError:
                    print(loc, year, month, 'file not found')
                    
        for lin in data:
            axes[mx+1][loc_ind].plot(lon_x, lin, lw=0.5, color='gray')
            
        if len(data)<10: continue
            
        # plot mean line
        mean_line = np.nanmean(data, axis=0)
        axes[mx+1][loc_ind].plot(lon_x, mean_line, lw=2, color='k', label='n='+str(len(data)))
        axes[mx+1][loc_ind].legend(loc='upper right')
        
        # add shading
        pos_change = np.ma.masked_where(mean_line<0, mean_line).filled(np.nan)
        neg_change = np.ma.masked_where(mean_line>0, mean_line).filled(np.nan)
        axes[mx+1][loc_ind].fill_between(lon_x, neg_change, y2=0, color='salmon', alpha=0.5)
        axes[mx+1][loc_ind].fill_between(lon_x, pos_change, y2=0, color='lightskyblue', alpha=0.5)

#%%% east-west sorting: storm impact

# set up plot
nrows = (len(grouped_months)*2)+1
fig, axes = plt.subplots(nrows, 2*2, figsize=(18,2.5*nrows), sharey=True,sharex=True,
                         gridspec_kw={"height_ratios":[0.001]+([1]*(nrows-1))})
fig.suptitle('\n\n\n\n'+'Normalized Storm Impact: '+ystr, fontweight='bold')
for ax1 in axes.flatten():
    ax1.set_xlim([xxx[0], xxx[-1]])
    # ax1.set_xticks()
    ax1.axhline(0, lw=0.55, color='gray', ls=':')
    
# organize data
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    
    for mx, months in enumerate(grouped_months):
        for ax2 in [axes[0][2*loc_ind], axes[0][2*loc_ind+1]]:
            ax2.axis("off")
            ax2.set_title('\n'+loc, fontweight='bold')
        
        miz_data = []
        impact_data = []
        for year in years:
            mean_lines, impacts, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines([year], root_path+'census/', root_path+'area/', root_path+'seaice/')
                
            for month in months:
                lines = get_miz_lines(loc, year, month, savepath, root_path)
                miz_data += lines[uv]
    
                impact_data+=impacts[month]
            
        # sort data
        ax_sort = [axes[2*mx+1][2*loc_ind], axes[2*mx+1][(2*loc_ind)+1],
                   axes[2*mx+2][2*loc_ind], axes[2*mx+2][(2*loc_ind)+1]
                   ]
        
        sorted_lines = {'W_E':[],'E_W':[], 'both':[], 'neither':[]} # increases (first)
        for miz1, impact in zip(miz_data, impact_data):
            west = np.nanmean(miz1[0:50])
            east = np.nanmean(miz1[50:])
            
            if west>0 and east<0: sorted_lines['W_E'].append(impact)
            elif west<0 and east>0: sorted_lines['E_W'].append(impact)
            elif west>0 and east>0: sorted_lines['both'].append(impact)
            elif west<0 and east<0: sorted_lines['neither'].append(impact)
            
        for key, ax in zip(list(sorted_lines.keys()), ax_sort):
            if len(sorted_lines[key])==0:ax.set_title('none');continue
            for line in sorted_lines[key]:
                ax.plot(xxx, line, lw=0.5, color='gray')
            ax.plot(xxx, np.nanmean(sorted_lines[key], axis=0), 
                    lw=2, color='k', label='n='+str(len(sorted_lines[key])))
            ax.set_title(key+' '+str(months))



#%% linear miz scaled by linear winds

# set up plot
nrows = len(grouped_months)+1
fig, axes = plt.subplots(nrows, 2, figsize=(12,3.5*nrows), sharey=True,
                         gridspec_kw={"height_ratios":[0.001]+([1]*(nrows-1))})
for ax1 in axes.flatten():
    ax1.set_xlim([lon_x[0], lon_x[-1]])
    ax1.set_xticks(np.arange(lon_x[0], lon_x[-1]+0.25, 0.25))
    ax1.axhline(0, lw=0.55, color='gray', ls=':')
    ax1.axvline(0.5, lw=0.55, color='gray', ls=':')
    
    ax1.set_ylim([-1.5e4, 1.5e4]) 
    
for ax1 in axes[-1][:]: ax1.set_xlabel('Longitude Extent Fraction')
for ax1 in axes[:,0]: ax1.set_ylabel('MIZ Area Change' +'\nscaled by wind speed')

# organize data
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    axes[0][loc_ind].axis("off")
    axes[0][loc_ind].set_title('\n'+loc, fontweight='bold')
    
    
    for mx, months in enumerate(grouped_months):
        axes[mx+1][loc_ind].set_title(ystr+'  '+str(months))
        
        # get linear lines and scale
        miz_data = []
        for month in months:
            for year in years:
                lines = get_miz_lines_seaice(loc, year, month, savepath_miz, root_path)
                miz_data += lines
        # get linear wind lines
        wind_data = []
        for month in months:
            for year in years:
                lines = get_miz_lines(loc, year, month, savepath, root_path)
                wind_data += lines[uv]
                
        # scale and plot    
        scaled = []
        for wind_line, miz_line in zip(wind_data, miz_data):
            wind_line = np.where(np.abs(wind_line)<1, np.nan, wind_line)
            
            scaled_line = np.array(miz_line)/np.array(wind_line)
            axes[mx+1][loc_ind].plot(lon_x, scaled_line, lw=0.5, color='gray')
            scaled.append(scaled_line)
            
        # plot mean line
        mean_line = np.nanmean(scaled, axis=0)
        axes[mx+1][loc_ind].plot(lon_x, mean_line, lw=2, color='k', label='n='+str(len(wind_data)))
        axes[mx+1][loc_ind].legend(loc='upper right')
        
        # add shading
        pos_change = np.ma.masked_where(mean_line<0, mean_line).filled(np.nan)
        neg_change = np.ma.masked_where(mean_line>0, mean_line).filled(np.nan)
        axes[mx+1][loc_ind].fill_between(lon_x, neg_change, y2=0, color='salmon', alpha=0.5)
        axes[mx+1][loc_ind].fill_between(lon_x, pos_change, y2=0, color='lightskyblue', alpha=0.5)


#%%%% line correlations

# set up plot
nrows = len(grouped_months)+1
fig, axes = plt.subplots(nrows, 2, figsize=(12,3.5*nrows), sharey=True,
                         gridspec_kw={"height_ratios":[0.001]+([1]*(nrows-1))})
fig.suptitle('Correlation between Linear MIZ and Wind')
for ax1 in axes.flatten():
    # ax1.set_xlim([lon_x[0], lon_x[-1]])
    # ax1.set_xticks(np.arange(lon_x[0], lon_x[-1]+0.25, 0.25))
    # ax1.axhline(0, lw=0.55, color='gray', ls=':')
    ax1.axvline(0, lw=0.55, color='gray', ls=':')
    
    # # ax1.set_ylim([-2e4, 2e4]) 
    # ax1.set_ylim([-2.5e3, 2.5e3])
    pass
    
for ax1 in axes[-1][:]: ax1.set_xlabel('Correlation')
for ax1 in axes[:,0]: ax1.set_ylabel('Number of Storms')

# organize data
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    axes[0][loc_ind].axis("off")
    axes[0][loc_ind].set_title('\n'+loc, fontweight='bold')
    
    
    for mx, months in enumerate(grouped_months):
        axes[mx+1][loc_ind].set_title(ystr+'  '+str(months))
        
        # get linear lines and scale
        miz_data = []
        for month in months:
            for year in years:
                lines = get_miz_lines_seaice(loc, year, month, savepath_miz, root_path)
                miz_data += lines
        # get linear wind lines
        wind_data = []
        for month in months:
            for year in years:
                lines = get_miz_lines(loc, year, month, savepath, root_path)
                wind_data += lines[uv]
                
        # scale and plot    
        corrs = []
        for wind_line, miz_line in zip(wind_data, miz_data):
            # wind_line = np.where(np.abs(wind_line)<1, np.nan, wind_line)
            
            mask = ~np.isnan(np.array(wind_line)) & ~np.isnan(np.array(miz_line))
            corr = np.corrcoef(np.array(wind_line)[mask], np.array(miz_line)[mask])[0][1]
            corrs.append(corr)
            
        # plot mean line
        axes[mx+1][loc_ind].hist(corrs, bins=np.arange(-1,1.1, 0.1),
                                 label='n='+str(len(wind_data)))
        
        axes[mx+1][loc_ind].axvline(np.nanmean(corrs), lw=1, color='k', ls='-')
        
        axes[mx+1][loc_ind].legend(loc='upper right')
        


#%% end
