#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 2025
linear_miz1.py
(code adapted from: antarctica/casestudies_aa3.py --> location separation, antarctica/casestudies_aa1.py)

* analyze east-west changes in MIZ

@author: mundi
"""
#%% imports, file names
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import time as timeIN

import functions as fx
from datetime import datetime, timedelta
from scipy.interpolate import griddata
import cmocean.cm as cmo

from concurrent.futures import ThreadPoolExecutor
from functools import wraps

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]
ice_path = '/Users/mundi/Desktop/seaice/'
ice_paths = [ice_path, ice_path+'south/']

savepath = '/Users/mundi/Desktop/month-hemi/linear_miz/'

census_name = 'census_'
contour_name = '_contours.nc'
si_name = '_seaice'

#%%% controls

years = np.arange(2010,2020)
# years = np.arange(1982, 1992)
ystr = str(years[0]) +'-'+ str(years[-1])

clim_years = np.arange(2010,2020)
# clim_years = np.arange(1982, 1992)

# month_options = [[2,3,4]]#,[5,6,7],[8,9,10],[11,12,1]]
month_options = [[2,3,4],[5,6,7],[8,9,10],[11,12,1]]

calc_time = 0 # 0 = storm duration, 1 = 7 days, 2 = 3 weeks

normed_lines = False
name_add = '' if normed_lines else 'all_'

xxx = np.arange(-7,14+1,1)
ice_lims = [20,80]
miz = [0.15,0.80]

lon_num = 100
lon_x = np.linspace(0, 1, lon_num)

#%%% functions

def timethis(func):
    """ 
    Print the execution time for a function call
    """
    @wraps(func)
    def wrapped_method(*args, **kwargs):
        time_start = timeIN.time()
        output = func(*args, **kwargs)
        time_end = timeIN.time()
        if time_end-time_start < 120:
            print(f"{func.__name__}: {(time_end-time_start)} s")
        else:
            print(f"{func.__name__}: {(time_end-time_start)/60} min")

        return output

    return wrapped_method

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

def get_si_diff(calc_time, start, enddate):
    
    dt1, dt2 = timediff(calc_time, start, enddate)
    
    if loc_ind==0: 
        sic1 = fx.load_seaice(ice_fname, dt1.year, dt1.month, dt1.day, latlon=False)
        sic2 = fx.load_seaice(ice_fname, dt2.year, dt2.month, dt2.day, latlon=False)
    elif loc_ind==1: 
        sic1 = fx.load_seaice_sh(ice_fname, dt1.year, dt1.month, dt1.day, latlon=False)
        sic2 = fx.load_seaice_sh(ice_fname, dt2.year, dt2.month, dt2.day, latlon=False)
    
    return sic2 - sic1

def get_clim_diff(calc_time, start, enddate, clim_years):
    dt1, dt2 = timediff(calc_time, start, enddate)
    
    if start.month==12 and enddate.month==1: yy_add = 1
    else: yy_add = 0
    
    if dt1.month==2 and dt1.day==29:
        dt1 = datetime(dt1.year, 3, 1)
    if dt2.month==2 and dt2.day==29:
        dt2 = datetime(dt2.year, 3, 1)
    
    differences = []
    for yy in clim_years:
        if loc_ind==0: 
            sic1 = fx.load_seaice(ice_fname, yy, dt1.month, dt1.day, latlon=False)
            sic2 = fx.load_seaice(ice_fname, yy+yy_add, dt2.month, dt2.day, latlon=False)
        elif loc_ind==1: 
            sic1 = fx.load_seaice_sh(ice_fname, yy, dt1.month, dt1.day, latlon=False)
            sic2 = fx.load_seaice_sh(ice_fname, yy+yy_add, dt2.month, dt2.day, latlon=False)
        
        differences.append(sic2 - sic1)

    return np.nanmean(differences, axis=0)


#%% data loop
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    print(); print('***', loc, '***')
    root_path = root_paths[loc_ind]
    ice_fname = ice_paths[loc_ind]
    
    if loc_ind == 0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010,1,1, latlon=True)
    elif loc_ind ==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010,1,1, latlon=True)
    si_lon1 = np.where(si_lon<0, si_lon+360, si_lon)
    
    for months in month_options:
        for month in months:
            print()
            print('- '+str(month)+' -')
            month_group_timer = timeIN.time()
        
            for year in years:
                print(year)
                savename = loc+'_miz_line_'+name_add+str(year)+'_'+str(month)+'.npy'
                try:
                    lines = np.load(savepath+savename)
                except FileNotFoundError:
                    try:
                        time_start = timeIN.time()
                        
                        # load census and sea ice info
                        census_name = 'census_'+ str(year) +'.csv'
                        [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
                            fx.readCensus(root_path+'census/'+census_name, convertDT=True)
                            
                        
                        #### loop thru storms
                        lines = []
                        for sn, start in enumerate(startdate):
                            if start.month not in months: continue
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
                        
                            # get sic difference
                            si_diff = get_si_diff(calc_time, start, end)
                            
                            ## get clim change
                            clim_diff = get_clim_diff(calc_time, start, end, clim_years)
                            
                            ## total diff
                            diff = si_diff - clim_diff
                           
                            # isolate miz
                            inside_bbox = fx.find_points_in_contour(bbox_edges, si_lon, si_lat)
                            diff_bbox = np.ma.masked_array(diff, mask=inside_bbox)
                            diff_miz = np.ma.masked_where(miz_points==0, diff_bbox).filled(np.nan)
                            
                            # compute areas for averaging
                            diff_miz_area = diff_miz*25*25
                            # miz_area = np.nansum( np.where(diff_miz>-2, 1, np.nan)*25*25 )
                            
                            # set up grid
                            lon_bnds = np.linspace( np.nanmin(bbox_edges[:,0]), np.nanmax(bbox_edges[:,0]), lon_num)
                            lat_bnds = np.linspace( np.nanmin(bbox_edges[:,1]), np.nanmax(bbox_edges[:,1]), lon_num)
                            new_x, new_y = np.meshgrid(lon_bnds, lat_bnds)   
                        
                        
                            #### calculate final values
                            diff_gridded = griddata((si_lon1.flatten(), si_lat.flatten()), diff_miz_area.flatten(),
                                                   (new_x.flatten(), new_y.flatten()))
                        
                            diff_gridded = np.where(diff_gridded==0, np.nan, diff_gridded)
                            diff_gridded = np.reshape(diff_gridded, np.shape(new_x))
                            
                            pixel_count = np.where(diff_gridded>-999, 1, np.nan)
                            pixel_sum = np.nansum(pixel_count, axis=0)
                            
                            if normed_lines:
                                LINE = np.nansum(diff_gridded, axis=0)*100/pixel_sum  #miz_area? # why x100?
                            else:
                                LINE = np.nansum(diff_gridded, axis=0)
                        
                            lines.append( LINE )
                    
                        # export data
                        np.save(savepath+savename, lines)
                        print('Loaded Data: '+str(round((timeIN.time()-time_start)/60,1))+' min')
            
                    except Exception as ee:
                        print(); print(); print(ee); print(year, month); print();print()
                        continue
                        
            
        print('--> '+ str(round((timeIN.time()-month_group_timer)/60,1))+' min')
        
#%% data organization

def get_miz_lines(loc, year, month, savepath, root_path):
    savename = loc+'_miz_line_'+str(year)+'_'+str(month)+'.npy'
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
        
#%% plot lines

# set up plot
nrows = len(month_options)+1
fig, axes = plt.subplots(nrows, 2, figsize=(12,3.5*nrows), sharey=True,
                         gridspec_kw={"height_ratios":[0.001]+([1]*(nrows-1))})
for ax1 in axes.flatten():
    ax1.set_xlim([lon_x[0], lon_x[-1]])
    ax1.set_xticks(np.arange(lon_x[0], lon_x[-1]+0.25, 0.25))
    ax1.axhline(0, lw=0.55, color='gray', ls=':')
    ax1.axvline(0.5, lw=0.55, color='gray', ls=':')
for ax1 in axes[-1][:]: ax1.set_xlabel('Longitude Extent Fraction')
for ax1 in axes[:,0]: ax1.set_ylabel('MIZ Area Change' +' per MIZ Width' if normed_lines else '')

# organize data
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    axes[0][loc_ind].axis("off")
    axes[0][loc_ind].set_title('\n'+loc, fontweight='bold')
    
    
    for mx, months in enumerate(month_options):
        axes[mx+1][loc_ind].set_title(ystr+'  '+str(months))
        data = []
        for month in months:
            for year in years:
                lines = get_miz_lines(loc, year, month, savepath, root_path)
                data += lines
        for lin in data:
            axes[mx+1][loc_ind].plot(lon_x, lin, lw=0.5, color='gray')
            
        # plot mean line
        mean_line = np.nanmean(data, axis=0)
        axes[mx+1][loc_ind].plot(lon_x, mean_line, lw=2, color='k', label='n='+str(len(data)))
        axes[mx+1][loc_ind].legend(loc='upper right')
        
        # add shading
        pos_change = np.ma.masked_where(mean_line<0, mean_line).filled(np.nan)
        neg_change = np.ma.masked_where(mean_line>0, mean_line).filled(np.nan)
        axes[mx+1][loc_ind].fill_between(lon_x, neg_change, y2=0, color='salmon', alpha=0.5)
        axes[mx+1][loc_ind].fill_between(lon_x, pos_change, y2=0, color='lightskyblue', alpha=0.5)


#%% east-west sorting: storm impact

# set up plot
nrows = (len(month_options)*2)+1
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
    
    for mx, months in enumerate(month_options):
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
                miz_data += lines
    
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
        
#%% monthly scatter
month_colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

from scipy.ndimage import uniform_filter1d
from scipy.stats import linregress

#%%% calculate monthly miz lines

# set up plot
fig, axes = plt.subplots(3, 2, figsize=(12, 8), sharey=True,
                         gridspec_kw={"height_ratios":[0.001,1,1]})
plt.subplots_adjust(hspace=0.3)
for ax1 in axes.flatten():
    ax1.set_xlim([lon_x[0], lon_x[-1]])
    ax1.set_xticks(np.arange(lon_x[0], lon_x[-1]+0.25, 0.25))
    ax1.axhline(0, lw=0.55, color='gray', ls=':')
    ax1.axvline(0.5, lw=0.55, color='gray', ls=':')
for ax1 in axes[-1][:]: ax1.set_xlabel('Longitude Extent Fraction')
for ax1 in axes[:,0]: ax1.set_ylabel('MIZ Area Change' +' per MIZ Width' if normed_lines else '')

# organize data
smoothed_data = []
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    axes[0][loc_ind].axis("off")
    axes[0][loc_ind].set_title('\n'+loc, fontweight='bold')
    
    N_smooth = 1
    
    axes[1][loc_ind].set_title('Monthly Mean')
    axes[2][loc_ind].set_title('Smoothed: N='+str(N_smooth))
    
    smoothed = {}
    for mi, month in enumerate(np.arange(1,12+1)):
        data = []
        for year in years:
            lines = get_miz_lines(loc, year, month, savepath, root_path)
            data += lines

        # plot mean line
        if len(data)>10:
            mean_line = np.nanmean(data, axis=0)
            mask = ~np.isnan(np.array(lon_x)) & ~np.isnan(np.array(mean_line))
            
            axes[1][loc_ind].plot(lon_x[mask], mean_line[mask], lw=2, color=month_colors[mi], 
                                     label=str(month)+', n='+str(len(data)))
            axes[1][loc_ind].legend(loc='upper right', bbox_to_anchor=(1.15,1))
            
            # window smoothing?
            smoothed_line = uniform_filter1d(mean_line[mask], size=N_smooth, mode='nearest')
            axes[2][loc_ind].plot(lon_x[mask], smoothed_line, lw=2, color=month_colors[mi], 
                                     label=str(month)+', n='+str(len(data)))
            smoothed[month] = smoothed_line
        else:
            smoothed[month] = []
            
    smoothed_data.append(smoothed)

#%%% scatter vs impact (monthly)

#%%%% fraction of miz line pos/neg

fig, axes = plt.subplots(1,2,figsize=(10,5), sharey=True, sharex=True)
fig.suptitle('MIZ Fractions vs. Storm Impact')

for loc_ind, loc in  enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    
    smoothed_lines = smoothed_data[loc_ind]
    
    axes[loc_ind].set_title(loc+': '+ystr)
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, root_path+'census/', root_path+'area/', root_path+'seaice/')
        
    xs, ys = [], []
    for month in np.arange(1,12+1):
        if len(smoothed_lines[month]) > 0:
            smoothed_line = np.array( smoothed_lines[month] )
            frac = len(np.where(smoothed_line>0)[0])*100/len(smoothed_line)
            
            ml = mean_lines[month-1][-1]
    
            axes[loc_ind].plot(frac, ml, color=month_colors[month-1],
                                    marker='o', markersize=10)
            xs.append(frac); ys.append(ml)
        
    mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
    m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
    x = np.linspace(np.min(xs), np.max(xs), 50)
    axes[loc_ind].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
     
    ax6a = axes[loc_ind].twinx()
    ax6a.sharey(axes[loc_ind])     
    ax6a.axhline(0, lw=0.55, color='gray', ls=':')
    ax6a.axvline(50, lw=0.55, color='gray', ls=':')
    ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
    ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
    ax6a.axis('off');

for ax in axes: 
    ax.set_xlabel('Fraction of MIZ with Increases in Area')
    ax.set_ylabel('Normalized Change in MIZ Ice Area')

#%%%% X - amount* pos/neg

if False:
    fig, axes = plt.subplots(1,2,figsize=(10,5), sharey=True, sharex=True)
    fig.suptitle('MIZ Fractions vs. Storm Impact')
    
    for loc_ind, loc in  enumerate(['Arctic', 'Antarctic']):
        root_path = root_paths[loc_ind]
        
        smoothed_lines = smoothed_data[loc_ind]
        
        axes[loc_ind].set_title(loc+': '+ystr)
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, root_path+'census/', root_path+'area/', root_path+'seaice/')
            
        xs, ys = [], []
        for month in np.arange(1,12+1):
            if len(smoothed_lines[month]) > 0:
                smoothed_line = np.array( smoothed_lines[month] )
                pos = np.nansum( smoothed_line[np.where(smoothed_line>0)[0]] )
                neg = np.nansum( smoothed_line[np.where(smoothed_line<0)[0]] )
                
                pos = np.nanmean(smoothed_line[0:50])
                neg = np.nanmean(smoothed_line[50:])
                
                ml = mean_lines[month-1][-1]
        
                axes[loc_ind].plot(pos, neg, color=month_colors[month-1],
                                   marker='o' if ml>0 else '>', markersize=20*abs(ml))
                xs.append(pos); ys.append(neg)
            
        mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
        m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        ax6a = axes[loc_ind].twinx()
        ax6a.sharey(axes[loc_ind])     
        ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.axvline(50, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');
    
    for ax in axes: 
        ax.set_xlabel('Positive Area Changes')
        ax.set_ylabel('Negative Area Changes')

#%%% scatter with pie

def draw_pie(dist, 
             xpos, 
             ypos, 
             colors=['#f7f7f7','#cccccc','#969696','#525252'], 
             dist_labels = ['','','',''],
             size=500,
             ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()
    
    color_iter = iter(colors)

    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])

        ax.scatter([xpos], [ypos], marker=xy, s=size, fc=next(color_iter), zorder=500)

    ax2 = ax.twinx()
    for col, lab in zip(colors, dist_labels):
        ax2.plot([],[], lw=4, color=col, label=lab)
    ax2.legend(loc='center left', handlelength=1)
    ax2.axis('off')

    return ax


#### figure
fig, axes = plt.subplots(1,2,figsize=(10,5), sharey=True, sharex=True)
fig.suptitle('MIZ Fractions vs. Storm Impact')

for loc_ind, loc in  enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    
    smoothed_lines = smoothed_data[loc_ind]
    
    axes[loc_ind].set_title(loc+': '+ystr)
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, root_path+'census/', root_path+'area/', root_path+'seaice/')
        
    xs, ys = [], []
    for month in np.arange(1,12+1):
        #### get sorted data
        miz_data = []
        for year in years:
            lines = get_miz_lines(loc, year, month, savepath, root_path)
            miz_data += lines
            
        sorted_lines = {'W_E':[],'E_W':[], 'both':[], 'neither':[]} # increases (first)
        sorted_lines = [0,0,0,0]
        for miz1 in miz_data:
            west = np.nanmean(miz1[0:50])
            east = np.nanmean(miz1[50:])
            
            if west>0 and east<0: sorted_lines[0]+=1
            elif west<0 and east>0: sorted_lines[1]+=1
            elif west>0 and east>0: sorted_lines[2]+=1
            elif west<0 and east<0: sorted_lines[3]+=1
        
        #### plot scatter
        if len(smoothed_lines[month]) > 0:
            smoothed_line = np.array( smoothed_lines[month] )
            frac = len(np.where(smoothed_line>0)[0])*100/len(smoothed_line)
            
            ml = mean_lines[month-1][-1]
    
            axes[loc_ind].plot(frac, ml, color=month_colors[month-1],
                                    marker='o', markersize=22)
            xs.append(frac); ys.append(ml)
            
            draw_pie(sorted_lines, frac, ml, ax = axes[loc_ind], size=350,
                     dist_labels=['WI_ED', 'WD_EI', 'Incr', 'Decr'])
        
    mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
    m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
    x = np.linspace(np.min(xs), np.max(xs), 50)
    axes[loc_ind].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
     
    ax6a = axes[loc_ind].twinx()
    ax6a.sharey(axes[loc_ind])     
    ax6a.axhline(0, lw=0.55, color='gray', ls=':')
    ax6a.axvline(50, lw=0.55, color='gray', ls=':')
    ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
    ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
    ax6a.axis('off');

for ax in axes: 
    ax.set_xlabel('Fraction of MIZ with Increases in Area')
    ax.set_ylabel('Normalized Change in MIZ Ice Area')
    
#%%% scatter with pie - E/W incr/decr

# or --> prpl
# ['#e66101','#fdb863','#f7f7f7','#b2abd2','#5e3c99']

def draw_pie(dist, 
             xpos, 
             ypos, 
             colors=['#e66101','#fdb863','#5e3c99','#b2abd2','#969696','#525252'], 
             dist_labels = ['','','',''],
             size=500,
             ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()
    
    color_iter = iter(colors)

    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])

        ax.scatter([xpos], [ypos], marker=xy, s=size, fc=next(color_iter), zorder=500)

    ax2 = ax.twinx()
    for col, lab in zip(colors, dist_labels):
        ax2.plot([],[], lw=4, color=col, label=lab)
    ax2.legend(loc='center left', handlelength=1)
    ax2.axis('off')

    return ax


#### figure
fig, axes = plt.subplots(1,2,figsize=(10,5), sharey=True, sharex=True)
fig.suptitle('MIZ Fractions vs. Storm Impact')

for loc_ind, loc in  enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    
    smoothed_lines = smoothed_data[loc_ind]
    
    axes[loc_ind].set_title(loc+': '+ystr)
    
    mean_lines, si_lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, root_path+'census/', root_path+'area/', root_path+'seaice/')
        
    xs, ys = [], []
    for month in np.arange(1,12+1):
        #### get sorted data
        miz_data = []
        for year in years:
            lines = get_miz_lines(loc, year, month, savepath, root_path)
            miz_data += lines
            
        si_series = si_lines[month]
            
        sorted_lines = {'W_E_decr':[], 'W_E_incr':[],
                        'E_W_decr':[], 'E_W_incr':[], 
                        'both':[], 'neither':[]} 
        sorted_lines = [0,0,0,0,0,0]
        for miz1, si in zip(miz_data, si_series):
            west = np.nanmean(miz1[0:50])
            east = np.nanmean(miz1[50:])
            
            
            if west>0 and east<0: 
                if si[-1]<0:
                    sorted_lines[0]+=1
                else:
                    sorted_lines[1]+=1
            elif west<0 and east>0: 
                if si[-1]<0:
                    sorted_lines[2]+=1
                else:
                    sorted_lines[3]+=1
            elif west>0 and east>0: sorted_lines[-2]+=1
            elif west<0 and east<0: sorted_lines[-1]+=1
        
        #### plot scatter
        if len(smoothed_lines[month]) > 0:
            smoothed_line = np.array( smoothed_lines[month] )
            frac = len(np.where(smoothed_line>0)[0])*100/len(smoothed_line)
            
            ml = mean_lines[month-1][-1]
    
            axes[loc_ind].plot(frac, ml, color=month_colors[month-1],
                                    marker='o', markersize=22)
            xs.append(frac); ys.append(ml)
            
            draw_pie(sorted_lines, frac, ml, ax = axes[loc_ind], size=350,
                     dist_labels=['WI_ED-Decr', 'WI_ED-Incr', 
                                  'WD_EI-Decr', 'WD_EI-Incr', 
                                  'Incr', 'Decr'])
        
    mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
    m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
    x = np.linspace(np.min(xs), np.max(xs), 50)
    axes[loc_ind].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
     
    ax6a = axes[loc_ind].twinx()
    ax6a.sharey(axes[loc_ind])     
    ax6a.axhline(0, lw=0.55, color='gray', ls=':')
    ax6a.axvline(50, lw=0.55, color='gray', ls=':')
    ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
    ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
    ax6a.axis('off');

for ax in axes: 
    ax.set_xlabel('Fraction of MIZ with Increases in Area')
    ax.set_ylabel('Normalized Change in MIZ Ice Area')

#%% ---
#%% si concentration gradient
import calendar
from scipy.stats import linregress

savepath1 = '/Users/mundi/Desktop/month-hemi/miz_gradient/'
name_add='all_'

def get_sic(start):
    
    if loc_ind==0: 
        return fx.load_seaice(ice_fname, start.year, start.month, start.day, latlon=False)
    elif loc_ind==1: 
        return fx.load_seaice_sh(ice_fname, start.year, start.month, start.day, latlon=False)
    
def get_grads(loc, year, month, savepath, root_path):
    savename = loc+'_sic_grad_all_'+str(year)+'_'+str(month)+'.npy'
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
    gradients = []
    for sn, start in enumerate(startdate):
        if start.month != month: continue
        val= next(miz_series)
        
        ice_frac = ice_sorter[sn]*100/box_area[sn]             
        if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
            continue
        
        gradients.append(val)
    
    return gradients


#%%% data loop

for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    print(); print('***', loc, '***')
    root_path = root_paths[loc_ind]
    ice_fname = ice_paths[loc_ind]
    
    if loc_ind == 0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010,1,1, latlon=True)
    elif loc_ind ==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2010,1,1, latlon=True)
    si_lon1 = np.where(si_lon<0, si_lon+360, si_lon)
    
    for months in month_options:
        for month in months:
            print()
            print('- '+str(month)+' -')
            month_group_timer = timeIN.time()
        
            for year in years:
                print(year)
                savename = loc+'_sic_grad_'+name_add+str(year)+'_'+str(month)+'.npy'
                try:
                    lines = np.load(savepath1+savename)
                except FileNotFoundError:
                    try:
                        time_start = timeIN.time()
                        
                        # load census and sea ice info
                        census_name = 'census_'+ str(year) +'.csv'
                        [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
                            fx.readCensus(root_path+'census/'+census_name, convertDT=True)
                            
                        #### loop thru storms
                        lines = []
                        for sn, start in enumerate(startdate):
                            if start.month not in months: continue
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
                            
                            # set up grid
                            lon_bnds = np.linspace( np.nanmin(bbox_edges[:,0]), np.nanmax(bbox_edges[:,0]), lon_num)
                            lat_bnds = np.linspace( np.nanmin(bbox_edges[:,1]), np.nanmax(bbox_edges[:,1]), lon_num)
                            new_x, new_y = np.meshgrid(lon_bnds, lat_bnds)
                            
                            # compute miz
                            miz_points = np.zeros(np.shape(si_lon))
                            for date in storm_range:
                                if loc_ind==0: sic = fx.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
                                elif loc_ind==1: sic = fx.load_seaice_sh(ice_fname, date.year, date.month, date.day, latlon=False)
                                miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)
                        
                            # get sic
                            SIC = get_sic(start)
                            
                            # isolate miz
                            inside_bbox = fx.find_points_in_contour(bbox_edges, si_lon, si_lat)
                            sic_bbox = np.ma.masked_array(SIC, mask=inside_bbox)
                            sic_miz = np.ma.masked_where(miz_points==0, sic_bbox).filled(np.nan)
                            
                            #### calculate final values
                            sic_gridded = griddata((si_lon1.flatten(), si_lat.flatten()), sic_miz.flatten(),
                                                   (new_x.flatten(), new_y.flatten()))
                        
                            sic_gridded = np.where(sic_gridded==0, np.nan, sic_gridded)
                            sic_gridded = np.reshape(sic_gridded, np.shape(new_x))
                            
                            grads = []
                            for xsec in sic_gridded.T:
                                x1 = np.arange(0, len(xsec))
                                mask = ~np.isnan(np.array(x1)) & ~np.isnan(np.array(xsec))
                                try: miz_slope = linregress(x1[mask], xsec[mask])[0]
                                except ValueError: continue
                                grads.append( miz_slope )
                            
                            lines.append( np.nanmean(grads) )
                    
                        # export data
                        np.save(savepath1+savename, lines)
                        print('Loaded Data: '+str(round((timeIN.time()-time_start)/60,1))+' min')
            
                    except Exception as ee:
                        print(); print(); print(ee); print(year, month); print();print()
                        continue
                        
            
        print('--> '+ str(round((timeIN.time()-month_group_timer)/60,1))+' min')
        
#%%% histogram plot

for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    
    root_path = root_paths[loc_ind]
    ice_fname = ice_paths[loc_ind]
    
    if years[0]>2000: txt_height = 10 if loc_ind==0 else 12
    else: txt_height = 8 if loc_ind==0 else 12

    fig_h, axes_h = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
    fig_h.suptitle('\n'+loc+': '+str(years[0])+'-'+str(years[-1]))
    axes_h = axes_h.flatten()
    for ai, axh in enumerate(axes_h): 
        axh.set_title(calendar.month_name[ai+1])
        axh.axvline(0, lw=0.55, color='gray', ls=':')
    bins = np.arange(-2, 16, 1)
    
    MONTHS = np.arange(1,12+1)
    for mi, month in enumerate(MONTHS):
        
        grads = []
        for year in years:
            try: grads += get_grads(loc, year, month, savepath1, root_path)
            except FileNotFoundError: break
        
        if loc_ind==0: grads = np.array(grads)
        elif loc_ind==1: grads = -np.array(grads)
            
        hist, bin_edges = np.histogram(np.array(grads)*100, bins=bins)
        axes_h[mi].bar(bin_edges[:-1], hist, width=1, 
                       facecolor='dodgerblue',edgecolor='navy', alpha=0.66)
        
        # mean gradient text !!!
        axes_h[mi].text(np.nanmean(bins)-1, txt_height, 
                        str(round(np.nanmean(grads)*100,2)), color='navy')
    
    fig_h.text(0.05, 0.5, 'Number of Storms',  va='center', ha='center', rotation=90, weight='bold')
    fig_h.text(0.5, 0.05, 'SIC Gradient (% per width)',  va='center', ha='center', weight='bold')
    
    # axes_h[-1].legend(ncol=2, bbox_to_anchor=(0.95, -0.175), handletextpad=0.5, handlelength=1)

#%%% are decade distributions statistically different?
from scipy.stats import ttest_ind
MONTHS = np.arange(1,12+1)

for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    
    root_path = root_paths[loc_ind]
    ice_fname = ice_paths[loc_ind]
    
    txt_height = 5 if loc_ind==0 else 8

    fig_h, axes_h = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
    fig_h.suptitle('\n'+loc+': '+str(years[0])+'-'+str(years[-1]))
    axes_h = axes_h.flatten()
    for ai, axh in enumerate(axes_h): 
        axh.set_title(calendar.month_name[ai+1])
        axh.axvline(0, lw=0.55, color='gray', ls=':')
    bins = np.arange(-2, 16, 1)
    
    for mi, month in enumerate(MONTHS):
        
        grads1 = []
        for year in np.arange(2010,2020):
            try: grads1 += get_grads(loc, year, month, savepath1, root_path)
            except FileNotFoundError: break
        mean1 = np.nanmean(grads1)*100
        
        grads2 = []
        for year in np.arange(1982,1992):
            try: grads2 += get_grads(loc, year, month, savepath1, root_path)
            except FileNotFoundError: break
        mean2= np.nanmean(grads2)*100
        
        if loc_ind==0: 
            grads1 = np.array(grads1)
            grads2 = np.array(grads2)
        elif loc_ind==1: 
            grads1 = -np.array(grads1)
            grads2 = -np.array(grads2)
            
        mask1 = ~np.isnan(grads1)
        mask2 = ~np.isnan(grads2)
        pval = ttest_ind(grads1[mask1], grads2[mask2])[1]
        if pval<0.05: tadd = '*'
        else: tadd=''
            
            
        hist1, bin_edges = np.histogram(np.array(grads1)*100, bins=bins)
        hist2, bin_edges = np.histogram(np.array(grads2)*100, bins=bins)
        axes_h[mi].bar(bin_edges[:-1], hist2-hist1, width=1, 
                       facecolor='lightcoral',edgecolor='maroon', alpha=0.66)
        
        # mean gradient text !!!
        axes_h[mi].text(np.nanmean(bins)+3, txt_height, 
                        str(round(mean2-mean1,2))+tadd, color='maroon')
    
    fig_h.text(0.05, 0.5, 'Number of Storms',  va='center', ha='center', rotation=90, weight='bold')
    fig_h.text(0.5, 0.05, 'SIC Gradient (% per width)',  va='center', ha='center', weight='bold')
    
#%%% simple scatter


fig, axes = plt.subplots(1,2,figsize=(10,5), sharey=True, sharex=True)
fig.suptitle('SIC Gradient vs. Storm Impact')

for loc_ind, loc in  enumerate(['Arctic', 'Antarctic']):
    root_path = root_paths[loc_ind]
    
    
    axes[loc_ind].set_title(loc+': '+ystr)
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, root_path+'census/', root_path+'area/', root_path+'seaice/')
        
    xs, ys = [], []
    for month in np.arange(1,12+1):
        ml = mean_lines[month-1][-1]
        
        grads1 = []
        for year in years:
            try: grads1 += get_grads(loc, year, month, savepath1, root_path)
            except FileNotFoundError: break
        mean1 = np.nanmean(grads1)*100
        
        if loc_ind==1: mean1 *= -1
    
        axes[loc_ind].plot(mean1, ml, color=month_colors[month-1],
                                marker='o', markersize=10)
        xs.append(mean1); ys.append(ml)
        
    mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
    m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
    x = np.linspace(np.min(np.array(xs)[mask]), np.max(np.array(xs)[mask]), 50)
    axes[loc_ind].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
     
    ax6a = axes[loc_ind].twinx()
    ax6a.sharey(axes[loc_ind])     
    ax6a.axhline(0, lw=0.55, color='gray', ls=':')
    # ax6a.axvline(50, lw=0.55, color='gray', ls=':')
    ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
    ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
    ax6a.axis('off');

for ax in axes: 
    ax.set_xlabel('MIZ Gradient')
    ax.set_ylabel('Normalized Change in MIZ Ice Area')

#%% end
