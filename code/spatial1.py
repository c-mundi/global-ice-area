#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 2025
spatial1.py

spatial changes in miz ice cover
- north and south hemisphere
- seasons vs months?
- decades

@author: mundi
"""
#%% imports and files
if True:
    import numpy as np
    import xarray as xr
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean.cm as cmo
    import cmocean.tools as cmtools
    from glob import glob
    from datetime import datetime, timedelta
    import calendar
    import string
    import gc
    import sys, os, warnings
    import time as timeIN
    from matplotlib.gridspec import GridSpec
    import functions as fx


plot_total=True

# months = [3,5,8,10,12] #
# months = [1,4,6,9,11]
months = [2,7,8,10,12]

decades = [np.arange(1982,1992), np.arange(2010,2020)]

decade_names = ['Early Satellite\nEra\n', 'Present Day\n']
if len(decades)!=len(decade_names):
    decade_names = [str(yy[0])+'-'+str(yy[-1])+'\n' for yy in decades]


nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]
ice_path = '/Users/mundi/Desktop/seaice/'
ice_paths = [ice_path, ice_path+'south/']

census_name = 'census_'
contour_name = '_contours.nc'
si_name = '_seaice'

import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

miz = [0.15, 0.80]

CMAP = cmo.balance_r

VMAX = 500#00 #7e3
VMIN = -VMAX

plot_ice_contour = True
plot_cyc_density = False



#%%% functions
#%%%% plotting
def setup_plot(ax, loc_ind, title=[], labels=True):
    
    if loc_ind==0:
        extent = [-180,180,60,90]
        y_block = 85
    elif loc_ind==1:
        extent = [-180,180,-90,-60]
        y_block = -85
        
    ax.set_boundary(circle, transform=ax.transAxes)
    
    # block off central arctic based on northern boundary
    xx = np.arange(0,360,1)
    yy = y_block*np.ones(np.shape(xx))
    mycirc = np.vstack((xx,yy)).T
    plt.rcParams['hatch.linewidth'] = 0.5

    ax.add_patch(Polygon(mycirc, closed=False,
                          fill=False, hatch='xxx',lw=0.5,color='gray',
                          transform=ccrs.PlateCarree())
                 )
    
    # geography
    ax.coastlines('50m',edgecolor='black',linewidth=0.75)
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    ax.add_feature(cfeature.LAKES, facecolor='0.85')
    if title:
        if type(title)!= str: title=str(title)
        ax.set_title(title)
    return ax

def plot_seaicecontour(ax, seaice, si_lon, label=[], linestyle='-', \
                       color='k', linewidth=2, levels=[0.15], zorder=10):
    
    si_lon = np.where(si_lon>180, si_lon-360, si_lon)

    #do masked-array on the lon
    lon_greater = np.ma.masked_greater(si_lon, -0.01)
    lon_lesser = np.ma.masked_less(si_lon, 0)
    # apply masks to other associate arrays: lat
    lat_greater = np.ma.MaskedArray(si_lat, mask=lon_greater.mask).filled(np.nan)
    lat_lesser = np.ma.MaskedArray(si_lat, mask=lon_lesser.mask).filled(np.nan)
    # apply masks to other associate arrays: daily ice
    si_greater = np.ma.MaskedArray(seaice, mask=lon_greater.mask).filled(np.nan)
    si_lesser = np.ma.MaskedArray(seaice, mask=lon_lesser.mask).filled(np.nan)
    
    lon_greater=lon_greater.filled(np.nan)
    lon_lesser=lon_lesser.filled(np.nan)

    # contours
    levels = levels # 15% ice extent definition
    contours1 = ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
                  linewidths = linewidth, zorder=zorder,linestyles=linestyle, 
                  transform=ccrs.PlateCarree()) 
    contours2 = ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
                  linewidths = linewidth, zorder=zorder,linestyles=linestyle, 
                  transform=ccrs.PlateCarree())
    return ax, [contours1,contours2]

def plot_geocontour(ax, lon, lat, var, levels, color='k', lw=3, ls='solid',label=False):
    #do masked-array on the lon
    lon_greater = np.ma.masked_greater(lon, -0.01)
    lon_lesser = np.ma.masked_less(lon, 0)    
    # apply masks to other associate arrays: lat
    lat_greater = np.ma.MaskedArray(lat, mask=lon_greater.mask)
    lat_lesser = np.ma.MaskedArray(lat, mask=lon_lesser.mask)
    # apply masks to other associate arrays: daily ice
    si_greater = np.ma.MaskedArray(var, mask=lon_greater.mask)
    si_lesser = np.ma.MaskedArray(var, mask=lon_lesser.mask)
    
    # si_greater = np.ma.masked_where(lat<58, si_greater)
    # si_lesser = np.ma.masked_where(lat<58, si_lesser)

    # contours
    cs1 = ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
              linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
              linestyles=ls) 
    cs2 = ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
              linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
              linestyles=ls)
    
    if label:
        ax.clabel(cs1, inline=True, fontsize=10)
        ax.clabel(cs2, inline=True, fontsize=10)
    
    return ax
#%%%% data
def get_storm_dates(root_path, year, month, return_conts=False):
    ice_lims = [20,80]
    
    census_file = root_path+'census/'+census_name+str(year)+'.csv'
    [startdate_in, enddate_in] = fx.readCensus(census_file, convertDT=True)[0]
    
    timing_grid = []
    for xx in range(0,len(startdate_in)):
        timing_grid.append((startdate_in[xx], enddate_in[xx]))
    storm_ranges = []
    for startdt, enddt in timing_grid:
        storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
    
    # open ice area
    ds_area = xr.open_dataset(root_path+'area/' + str(year) +'_area.nc')
    ice_sorter = ds_area['ice_area80'].values
    box_area = ds_area['box_area'].values
    
    try:
        ds = xr.open_dataset(root_path+'seaice/' + str(year) + si_name + '.nc')
        ds.close()
    except:
        print('- skip: '+root_path+'seaice/' + str(year) + si_name + '.nc')
        return [], []
    
    if return_conts: out_conts = []
    storm_areas = []
    startdate, enddate = [],[]
    stormstr_prev=''
    for storm_num, strm in enumerate(storm_ranges):
        if strm[0].month != month: continue
        
        # remove storms that don't interact with the ice
        ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
        if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
            continue
        
        # load storm contours
        stormstr1 = strm[0].strftime('%Y_%m%d')
        # duplicate storm start date?
        if stormstr1==stormstr_prev:
            stormstr = stormstr1 + '*'
            dupe=True
        else:
            dupe=False
            stormstr=stormstr1
        stormstr_prev = stormstr1
        
        #### storm area
        ncname = stormstr + contour_name +'.nc'
        try:
            if not dupe:
                cs = xr.open_dataset(root_path+'contours/'+ncname)
            else:
                raise FileNotFoundError
        except FileNotFoundError:
            contour_files = glob(root_path+'contours/'+'*'+strm[0].strftime('%Y_%m%d')+'*.nc')
            cs = xr.open_dataset(contour_files[-1])
            
        all_contours = []
        for key in list(cs.keys()):
            coord = cs[key].values
            all_contours.append(coord)
        cs.close()
        
        if all_contours == []: 
            print('-', stormstr) 
            continue
        
        if return_conts:
            out_conts.append(all_contours)
        
        with fx.HidePrint(): bbox_edges = fx.get_bbox_edges(all_contours) 
        storm_areas.append( bbox_edges )
        
        startdate.append(strm[0])
        enddate.append(strm[-1])
      
    if not return_conts: return startdate, enddate, storm_areas
    else: return startdate, enddate, storm_areas, out_conts

#%%%% sea ice

def get_miz_area(si_path, storm_range):
    # compute miz area
    miz_points = np.zeros(np.shape(si_lon))
    for date in storm_range:
        if si_path.split('/')[-2]!= 'south':
            sic = fx.load_seaice(si_path, date.year, date.month, date.day, latlon=False)
        else:
            sic = fx.load_seaice_sh(si_path, date.year, date.month, date.day, latlon=False)
        
        miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
        
    return miz_points

def isolate_area(inside_points, miz_points, si):
    si_in = np.where(inside_points, np.nan, si)
    si_miz = np.where(miz_points<1, np.nan, si_in)
    return si_miz*25*25

def compute_climatology(si_path, plot_date, inside_points, miz_points):
    yearly_changes = []
    
    if plot_date.year in np.arange(2010, 2020):
        clim_years = np.arange(2010, 2020)
    elif plot_date.year in np.arange(1982,1992):
        clim_years = np.arange(1982,1992)
    else:
        print('missing climatology years: ', plot_date.year)
    
    for yy in clim_years:
        
        # get dates
        try:
            dt = datetime(yy, plot_date.month, plot_date.day)
        except ValueError:
            dt = datetime(yy, 3, 1)
        
        # load yearly sea ice
        if si_path.split('/')[-2]!= 'south':
            si = fx.load_seaice(si_path, dt.year, dt.month, dt.day, latlon=False)
        else:
            si = fx.load_seaice_sh(si_path, dt.year, dt.month, dt.day, latlon=False)
        
        # storm area and miz
        si_miz = isolate_area(inside_points, miz_points, si)
        
        # append
        yearly_changes.append( si_miz )
        
    with warnings.catch_warnings(action="ignore"):
        return np.nanmean(yearly_changes, axis=0)   

def ice_changes(si_path, plot_date, miz_points, bb):
    '''
    Parameters
    ----------
    plot_date : datetime, storm end date
    bb : 2-n aarray, storm area (bounding boxes)
    seaice : xarray dataset
            opened cesm sea ice file for searching
    Returns
    -------
    storm difference map

    '''
    # get starting and ending ice concentrations
    if si_path.split('/')[-2]!= 'south':
        si = fx.load_seaice(si_path, plot_date.year, plot_date.month, plot_date.day, latlon=False)
    else:
        si = fx.load_seaice_sh(si_path, plot_date.year, plot_date.month, plot_date.day, latlon=False)
    
    # isolate storm areas
    inside_points = fx.find_points_in_contour(bb, si_lon, si_lat)
    
    si_area = isolate_area(inside_points, miz_points, si)
    
    # get climatological change
    clim_change = compute_climatology(si_path, plot_date, inside_points, miz_points)
   
    # storm difference
    storm_difference = si_area - clim_change
    
    return storm_difference


#%% DATA AND PLOT
#### set up plot
from matplotlib.patches import Polygon
fontsize = 22
ncols = len(months)
nrows = len(decades)*2
if plot_total: ncols+=1

fig = plt.figure(figsize=(5*ncols,4*nrows)) # (*2)=hemispheres
gs = GridSpec(nrows=nrows, ncols=ncols+2, 
              width_ratios=(1,1,1,1,1,1,0.55,0.55))


axes = []
for row in np.arange(nrows):
    ax_row = []
    for col in np.arange(ncols):
        if row in np.arange(len(decades)): projection=ccrs.NorthPolarStereo()
        elif row>=len(decades): projection=ccrs.SouthPolarStereo()
        ax_row.append(fig.add_subplot(gs[row, col], projection=projection))
    axes.append(ax_row)
axes=np.array(axes)


if len(decades)==1 and len(months)==1:
    axes = np.array([axes])


yiter = 1/(nrows+1)
ycoords = []
ystart = 1
while ystart >= yiter+0.01:
    ystart -= yiter
    ycoords.append(ystart)
ycoords = iter(ycoords)
    
for ax, mm in zip(axes.flatten()[0:len(months)], months):
    ax.set_title(calendar.month_name[mm]+'\n', fontsize=fontsize+2)

for li, lo in enumerate(['NH','SH']):
    for di, dec in zip(np.arange(0,len(decades)), decades):
        fig.text(0.0725, next(ycoords), decade_names[di]+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=fontsize+2)
    
# colorbar
expo = int( np.floor(np.log10(np.abs(VMAX))) )
if expo>2: 
    vmin1 = VMIN/(10**expo)
    vmax1 = VMAX/(10**expo)
else: 
    vmin1 = VMIN
    vmax1 = VMAX

cb = ax.pcolormesh(np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,2)), 
                   cmap = CMAP,
                   vmin=vmin1, vmax=vmax1)
cax1 = fig.add_axes([0.25,0.033,0.4,0.033]) 
cbar1 = fig.colorbar(cb, cax=cax1, orientation='horizontal')
if expo>2:
    cbar1.set_label(r'Change in Sea Ice Area ($\times 10^{}$ km$^2$)'.format(expo), fontsize=fontsize)
else:
    cbar1.set_label(r'Change in Sea Ice Area (km$^2$)', fontsize=fontsize)
cax1.tick_params(labelsize=fontsize)


##################################
#### month, decade, ensemble loop
##################################
for loc_ind, loc in enumerate(['Arctic','Antarctic']):
    print('--------->', loc)
    root_path = root_paths[loc_ind]
    ice_fname = ice_paths[loc_ind]
    
    if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010, 1, 1, latlon=True)
    elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010, 1, 1, latlon=True)

    for mm, month in enumerate(months):
        print('***', month , '***')
        month_start_time = timeIN.time()
        gc.collect()
         
        total_storm_count = 0
        
        for yx, years in enumerate(decades):
            print(str(years[0])+'-'+str(years[-1]))
            
            savename = 'plot_difference_'+str(month)+'_'+str(years[0])+'_'+str(years[-1])+'.npy'
            try:
                plot_change = np.load(root_path+'spatial/' +savename )
                
            except FileNotFoundError:
            
                decade_changes = []
                storm_durations = []                
                for year in years:

                    ### get storm dates and areas
                    startdate, enddate, storm_areas = get_storm_dates(root_path, year, month)
                
                    ### calculate ice changes
                    for sd, ed, bb in zip(startdate, enddate, storm_areas):
                        if len(bb)==0: print('#', sd); continue 
                        total_storm_count += 1
                        
                        # set up storm range for miz area determination
                        t1 = sd - timedelta(days=1)
                        t2 = ed + timedelta(days=1)
                        storm_range = fx.daterange(t1, t2, dt=24)
                        storm_duration = ed-sd
                        
                        # compute miz area
                        miz_points = get_miz_area(ice_fname, storm_range)
                    
                        # compute values and append for avg'ing
                        change_start = ice_changes(ice_fname, sd, miz_points, bb)
                        change_end = ice_changes(ice_fname, ed, miz_points, bb)
                        
                        storm_difference = change_end - change_start 
                    
                        decade_changes.append( storm_difference )
                        storm_durations.append(storm_duration.days * np.ones(np.shape(storm_difference)))
                        
                        gc.collect()
    
                changes = np.array(decade_changes)
                plot_change = np.nansum(changes, axis=0)
                print('>', np.nansum(plot_change))
                
                np.save(root_path+'spatial/'+savename, plot_change)
                print('*SAVED*', np.shape(plot_change))
        
                
            if len(months)>1 and len(decades)>1:
                plotax = axes[yx+(loc_ind*len(decades))][mm]
            elif len(months)==1 and len(decades)>1:
                plotax = axes[yx+(loc_ind*len(decades))]
            elif len(months)>1 and len(decades)==1:
                plotax = axes[mm]
            elif len(months)==1 and len(decades)==1:
                plotax = axes.item()
                
            plotax = setup_plot(plotax, loc_ind, labels=False)
                
            plot_change = np.where(plot_change==0, np.nan, plot_change)  
            if expo>2: plot_var = plot_change/(10**expo)
            else: plot_var = plot_change
            
            if ~np.all(np.isnan(plot_var)):
                change =  plotax.pcolormesh(si_lon, si_lat, plot_var, 
                                    cmap = CMAP,
                                    vmin=VMIN, vmax=VMAX,
                                   transform = ccrs.PlateCarree())
                
            net_change = int(round(np.nansum(plot_var),0))
            plotax.text(0.05, -0.15, f"{net_change:,}"+r' km$^2$',
                        transform=plotax.transAxes, fontsize=fontsize-2)
                
            
            #### add starting sea ice contour
            if plot_ice_contour:
                try:
                    # load file
                    fname = 'ice_extent_'+str(month)+'_'+str(years[0])+'_'+str(years[-1])+'.npy'
                    starting_ice = np.load(root_path+'spatial/'+fname)
                    
                except FileNotFoundError: # calculate
                    # calculate mean starting extent
                    print('... calculating '+str(years[0])+'-'+str(month)+' sea ice extent')
                    
                    starting_ices = []
                    for year in years:
                        if loc_ind==0:
                            si_start = fx.load_seaice(ice_fname, year, month, 1, latlon=False)
                        elif loc_ind==1:
                            si_start = fx.load_seaice_sh(ice_fname, year, month, 1, latlon=False)
                        starting_ices.append(si_start)
        
                    # ens. mean starting sea ice extent
                    with warnings.catch_warnings(action="ignore"): # empty slices
                        starting_ice = np.nanmean(starting_ices, axis=0)
                        
                    # save si extent
                    savename = 'ice_extent_'+str(month)+'_'+str(years[0])+'_'+str(years[-1])+'.npy'
                    np.save(root_path+'spatial/'+savename, starting_ice)
                    
                with warnings.catch_warnings(action="ignore"):
                    plot_seaicecontour(plotax, starting_ice, si_lon, 
                                       levels=[0.15], linewidth=1.75, zorder=500)
    
        print('> '+str(round((timeIN.time()-month_start_time)/60,1))+' min')
        print()
        
#%% total column
if plot_total:
    print('*** Total ***')
    for loc_ind, loc in enumerate(['Arctic','Antarctic']):
        root_path = root_paths[loc_ind]
        ice_fname = ice_paths[loc_ind]
        if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010, 1, 1, latlon=True)
        elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010, 1, 1, latlon=True)
        
        for yx, years in enumerate(decades):
            total_sum1 = []
            for mm, month in enumerate(months):
                savename = 'plot_difference_'+str(month)+'_'+str(years[0])+'_'+str(years[-1])+'.npy'
                total_sum1.append( np.load(root_path+'spatial/'+savename) )
                
            total_sum = [x for x in total_sum1 if ~np.all(np.isnan(x))]
            
            plot_change = np.nansum(total_sum,axis=0)
            plot_change = np.where(plot_change==0, np.nan, plot_change)  
            if ~np.all(np.isnan(plot_change)):
                axes[yx+(loc_ind*len(decades))][-1] = setup_plot(axes[yx+(loc_ind*len(decades))][-1], 
                                                                 loc_ind, labels=False)
                
                axes[yx+(loc_ind*len(decades))][-1].pcolormesh(si_lon, si_lat, plot_change, 
                                                            cmap = CMAP, vmin=VMIN, vmax=VMAX,
                                                            transform = ccrs.PlateCarree())
            
            net_change = int(round(np.nansum(plot_change),0))
            print(str(years[0])+'-'+str(years[-1])+':', net_change)
            axes[yx+(loc_ind*len(decades))][-1].text(0.05, -0.15, f"{net_change:,}"+r' km$^2$',
                                                     transform=axes[yx+(loc_ind*len(decades))][-1].transAxes, 
                                                     fontsize=fontsize-2)
    
        
        axes[0][-1].set_title('Net Change\n', fontsize=fontsize+2)     
        line = plt.Line2D((.666,.666),(.1,.875), color="k", linewidth=1.5)
        fig.add_artist(line)    
        print('*************')

#%% cyclone locations: contours of storm density

def cyclone_density(root_path, month, years, fname, si_lon):
    
    density_array = np.zeros(np.shape(si_lon))
    for year in years:
        
        ### get storm dates and areas
        startdate, enddate, storm_areas, out_conts = get_storm_dates(root_path, year, month, return_conts=True)
        
        #### get storm areas
        for sx, sd in enumerate(startdate):
            storm_conts = out_conts[sx]
            
            all_contours = np.zeros(np.shape(si_lon))
            for coord in storm_conts:
                inside_pts = fx.find_points_in_contour(coord, si_lon, si_lat)
                all_contours[~inside_pts]+=1
                                                        # weight by storm duration
            density_array[np.where(all_contours>0)] += (enddate[sx]-startdate[sx]).days
            
    ## save si extent
    np.save(root_path+'spatial/'+fname, density_array)
    
    return np.nanmean(density_array,axis=0)
    

if plot_cyc_density:
    for loc_ind, loc in enumerate(['Arctic','Antarctic']):
        root_path = root_paths[loc_ind]
        ice_fname = ice_paths[loc_ind]
        if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010, 1, 1, latlon=True)
        elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010, 1, 1, latlon=True)
        
        for mm, month in enumerate(months):
            for yx, years in enumerate(decades):
                
                try: # loading cyc density file
                    fname = 'cyc_density1_'+str(month)+'_'+str(years[0])+'_'+str(years[-1])+'.npy'
                    density_array = np.load(root_path+'spatial/'+fname)
                except FileNotFoundError:
                    print('... calculating '+str(years[0])+'-'+str(month)+' cyclone density')
                    try:
                        density_array = cyclone_density(root_path, month, years, fname, si_lon)
                    except Exception as EE:
                        print(years, month, EE)
                        continue
                
                if len(months)>1 and len(decades)>1:
                    plotax = axes[yx+(loc_ind*len(decades))][mm]
                elif len(months)==1 and len(decades)>1:
                    plotax = axes[yx+(loc_ind*len(decades))]
                elif len(months)>1 and len(decades)==1:
                    plotax = axes[mm]
                elif len(months)==1 and len(decades)==1:
                    plotax = axes
        
                xx = np.where(si_lon>180, si_lon-360, si_lon)
                density_array = np.where(si_lat>85, np.nan, density_array)
                plotax = plot_geocontour(plotax, xx, si_lat, density_array, 
                                         levels=np.arange(0,125,25), color='k', 
                                         lw=0.5, ls='solid', label=True)
            


#%% --
#%% monthly time series

axes = [fig.add_subplot(gs[0:len(decades), -2:]), fig.add_subplot(gs[len(decades):, -2:])]

decade_colors = ['maroon','navy']
all_months = np.arange(1,12+1)
plot_months = all_months #np.array([3,5,8,10,12])

for loc_ind, loc in enumerate(['Arctic','Antarctic']):
    root_path = root_paths[loc_ind]
    
    ax = axes[loc_ind]

    ax.set_ylim([0, 13])
    ax.axvline(0, color='#888888', ls='--', lw=2)
    ax.yaxis.tick_right()
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
    
    ax.set_yticks(all_months, [calendar.month_abbr[mm][0] for mm in all_months], fontsize=fontsize-2)
    scale=1e6
    ax.set_xticks(np.arange(-2.5,3.5,1))
    for item in ax.get_xticklabels(): item.set_fontsize(fontsize-2)
    ax.xaxis.set_tick_params(width=3)
    ax.yaxis.set_tick_params(width=3)
    
    
    #### calculations

    
    for yx, years in enumerate(decades):
        month_series = []
        for mm, month in enumerate(plot_months): ###!!! all months!
            savename = 'plot_difference_'+str(month)+'_'+str(years[0])+'_'+str(years[-1])+'.npy'
            grd = np.load(root_path+'spatial/'+savename)
            if np.all(np.isnan(grd)): grd = np.zeros(np.shape(si_lon))
            month_series.append(np.nansum(grd)/1e6)
            
        #### plot
        ax.plot(month_series, plot_months, color=decade_colors[yx], lw=3)
        
    
axes[-1].set_xlabel('Change in Sea Ice Area\n'+r'($\times 10^6$ km$^2$)', fontsize=fontsize+2)




#%% ---
raise NameError('Done!')

#%% individual plotting
loc_ind = 1
month = 3

root_path = root_paths[loc_ind]
ice_fname = ice_paths[loc_ind]
projections = [ccrs.NorthPolarStereo(), ccrs.SouthPolarStereo()]

if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010, 1, 1, latlon=True)
elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010, 1, 1, latlon=True)

#%%% month impacts

# fig_test, ax = plt.subplots(1,1, subplot_kw={'projection':ccrs.SouthPolarStereo()})

# ax = setup_plot(ax, loc_ind, labels=False)
# ax.set_boundary(circle, transform=ax.transAxes)

# # block off central arctic based on northern boundary
# xx = np.arange(0,360,1)
# yy = 85*np.ones(np.shape(xx))
# mycirc = np.vstack((xx,yy)).T
# plt.rcParams['hatch.linewidth'] = 0.5
# ax.add_patch(Polygon(mycirc, closed=False,
#                       fill=False, hatch='xxx',lw=0.5,color='gray',
#                       transform=ccrs.PlateCarree())
#              )


#%%% cyc freq
if False:

    fig_test, axes = plt.subplots(1,2, subplot_kw={'projection':projections[loc_ind]})
    for ax in axes: ax = setup_plot(ax, loc_ind, labels=False)
    
    for yx, years in enumerate(decades):
        
        try: # loading cyc density file
            fname = 'cyc_density1_'+str(month)+'_'+str(years[0])+'_'+str(years[-1])+'.npy'
            density_array = np.load(root_path+'spatial/'+fname)
        except FileNotFoundError:
            print('... calculating '+str(years[0])+'-'+str(month)+' cyclone density')
            try:
                density_array = cyclone_density(root_path, month, years, fname, si_lon)
            except Exception as EE:
                print(years, month, EE)
                continue
      
        xx = np.where(si_lon>180, si_lon-360, si_lon)
        density_array = np.where(si_lat>85, np.nan, density_array)
        # axes[yx] = plot_geocontour(axes[yx], xx, si_lat, density_array, 
        #                          levels=np.arange(0,50,10), color='k', 
        #                          lw=0.5, ls='solid', label=True)
        
        cb = axes[yx].pcolormesh(si_lon, si_lat, density_array,
                            transform=ccrs.PlateCarree(),
                            vmin=0, vmax=150)
        
        cax1 = fig_test.add_axes([0.25,0.033,0.4,0.033]) 
        cbar1 = fig_test.colorbar(cb, cax=cax1, orientation='horizontal')
        cax1.tick_params(labelsize=fontsize)
        
        axes[yx] = plot_geocontour(axes[yx], xx, si_lat, density_array, 
                                 levels=np.arange(0,50,10), color='k', 
                                 lw=1.5, ls='solid', label=True)


#%%% sea ice contours


















#%% end
