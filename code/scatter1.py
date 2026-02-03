#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 2025
scatter1.py

scatter plots

@author: mundi
"""
#%% imports and filepaths
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import functions as fx
from datetime import datetime, timedelta
import string, calendar

import cmocean.cm as cmo
import cartopy.crs as ccrs
from scipy.stats import linregress
import warnings


nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

decades = [np.arange(2010,2020),np.arange(1982, 1992)]
decade_names = ['Present Day ', 'Early Satellite Era ']

months = np.arange(1,12+1)

ice_lims = [20,80]
var_addons = [['_1000', '_990', '_contour'],['2_1000', '2_990', '2_contour']]
addon_num = 0 # storm area
sia_var_addons = var_addons
miz_ind = 1 #[0=daily, 1=total]
miz_names = ['daily_', 'total_']

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

hemi_names= ['Arctic', 'Antarctic']

#%% info

def scatter_update(axes, i, mi, x6, y6, del_x, del_y, scale=False):
    if scale:
        del_x/=1e5
        del_y/=1e5
    axes[i].plot(del_x, del_y, marker='.', markersize=10, color=colors[mi])
    if not (np.isnan(del_x) or np.isnan(del_y)):
        x6[i].append(del_x)
        y6[i].append(del_y)
    return x6, y6



colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']


#%% plot functions

            
def monthly_plot_siclim(ax, x6, y6, mi, monthly_starts, monthly_ends, VAR_X, VAR_Y, scale=True):
   
    for start, end, var_x, var_y in zip(monthly_starts, monthly_ends, VAR_X, VAR_Y):
        ### storm duration
        try:
            del_x = var_x[7+(end-start).days+1] - var_x[7]
            del_y = (var_y[7+(end-start).days+1] - var_y[7])
            if scale: del_y = del_y - del_x
        except IndexError:
            del_x = var_x[-1] - var_x[7]
            del_y = (var_y[-1] - var_y[7]) 
            if scale: del_y = del_y - del_x
        x6, y6  = scatter_update(ax, 0, mi, x6, y6, del_x, del_y, scale)
        
        ### one week
        del_x = var_x[14] - var_x[7]
        del_y = (var_y[14] - var_y[7])
        if scale:del_y = del_y - del_x
        x6, y6  = scatter_update(ax, 1, mi, x6, y6, del_x, del_y, scale)
        
        ### three weeks
        del_x = var_x[-1] - var_x[0]
        del_y = (var_y[-1] - var_y[0])
        if scale: del_y = del_y - del_x
        x6, y6  = scatter_update(ax, 2, mi, x6, y6, del_x, del_y, scale)
        
    return x6, y6

def monthly_plot_si(axes, x6, y6, mi, monthly_starts, monthly_ends, si_clim, si, VAR):
    FLOATS = [np.float64,np.float32,float]
    
    for start, end, clim, sia, var in zip(monthly_starts, monthly_ends, si_clim, si, VAR):
        ### storm duration
        try:
            del_clim = clim[7+(end-start).days+1] - clim[7]
            del_sia = (sia[7+(end-start).days+1] - sia[7])
            del_sia = del_sia - del_clim
            if type(var) not in FLOATS: 
                del_var = var[7+(end-start).days+1] - var[7]
            else: del_var = var
        except IndexError:
            del_clim = clim[-1] - clim[7]
            del_sia = (sia[-1] - sia[7]) 
            del_sia = del_sia - del_clim
            if type(var) not in FLOATS: 
                del_var = var[-1] - var[7]
            else: del_var = var
        x6, y6  = scatter_update(axes, 0, mi, x6, y6, del_var, del_sia/1e5)
        
        ### one week
        del_clim = clim[14] - clim[7]
        del_sia = (sia[14] - sia[7])
        del_sia = del_sia - del_clim
        if type(var) not in FLOATS:  del_var = var[14] - var[7]
        else: del_var = var
        x6, y6  = scatter_update(axes, 1, mi, x6, y6, del_var, del_sia/1e5)
        
        ### three weeks
        del_clim = clim[-1] - clim[0]
        del_sia = (sia[-1] - sia[0])
        del_sia = del_sia - del_clim
        if type(var) not in FLOATS: del_var = var[-1] - var[7]
        else: del_var = var
        x6, y6  = scatter_update(axes, 2, mi, x6, y6, del_var, del_sia/1e5)
        
    return x6, y6

def monthly_plot_mean(axes, x6, y6, mi, monthly_starts, monthly_ends, si_clim, si, VAR):
   
    for start, end, clim, sia, var in zip(monthly_starts, monthly_ends, si_clim, si, VAR):
        ### storm duration
        try:
            del_clim = clim[7+(end-start).days+1] - clim[7]
            del_sia = (sia[7+(end-start).days+1] - sia[7])
            del_sia = del_sia - del_clim
            del_var = np.nanmean(var[7: 7+(end-start).days+1])
        except IndexError:
            del_clim = clim[-1] - clim[7]
            del_sia = (sia[-1] - sia[7]) 
            del_sia = del_sia - del_clim
            del_var = np.nanmean(var[7:])
        x6, y6  = scatter_update(axes, 0, mi, x6, y6, del_var, del_sia/1e5)
        
        ### one week
        del_clim = clim[14] - clim[7]
        del_sia = (sia[14] - sia[7])
        del_sia = del_sia - del_clim
        del_var = np.nanmean(var[7:14])
        x6, y6  = scatter_update(axes, 1, mi, x6, y6, del_var, del_sia/1e5)
        
        ### three weeks
        del_clim = clim[-1] - clim[0]
        del_sia = (sia[-1] - sia[0])
        del_sia = del_sia - del_clim
        del_var = np.nanmean(var)
        
        x6, y6  = scatter_update(axes, 2, mi, x6, y6, del_var, del_sia/1e5)
        
    return x6, y6

def scatter_plot_setup(xlabel, ylabel, suptitle, hline=None, vline=None):
    alph = ['a','b','c','d','e','f']
    
    fig6, axes6 = plt.subplots(2,3, figsize=(18,10), sharey=True)
    for i, ax6 in enumerate(axes6.flatten()):
        ax6.set_xlabel(xlabel)
        ax6.set_ylabel(ylabel)
    
        ax6.xaxis.set_tick_params(labelbottom=True)
        ax6.yaxis.set_tick_params(labelleft=True)
        
        if hline or type(hline)==int: ax6.axhline(hline, lw=0.75, color='darkgray')
        if vline or type(vline)==int: ax6.axvline(vline, lw=0.75, color='darkgray')
        
        ax6.text(0.0225, 1.025, '('+alph[i]+')',transform=ax6.transAxes, 
                  fontsize=12, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                      'edgecolor':'white', 'lw':0.75},zorder=50)
    fig6.suptitle(suptitle, fontweight='bold')
    
    return fig6, axes6

def plot_final_scatter(ax_row, x6, y6, xlims=None):
    titles6 = ['Storm Duration', 'One Week Post-Storm',
                '3-Week Analysis']
    
    for i, ax6 in enumerate(ax_row):
        mask = ~np.isnan(np.array(x6[i])) & ~np.isnan(np.array(y6[i]))
        m, b, r, p, se = linregress(np.array(x6[i])[mask], np.array(y6[i])[mask])
        # x = np.linspace(np.min(x6[i]), np.max(x6[i]), 100)
        if xlims: x = np.linspace(xlims[i][0], xlims[i][-1], 100)
        else: x = np.linspace(np.nanmin(x6[i]), np.nanmax(x6[i]), 100)
        ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
        ax6.set_title(titles6[i])
        if xlims: ax6.set_xlim(xlims[i])
        
        ax6a = ax6.twinx()
        ax6a.sharey(ax6)  
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

    return ax_row

#%% ----

#%% sea ice area
for loc_ind, location in enumerate(hemi_names):

    fig6, axes6 = scatter_plot_setup(xlabel=r'$\Delta$ Sea Ice Climatology ($\times10^5$ km$^2$)',
                                     ylabel='Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location,
                                     hline=0, vline=0)
            
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        
        path1 = root_paths[loc_ind]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
    
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
        
        for mi, mm in enumerate(list(lines.keys())):
            monthly_starts = start_day[mm]
            monthly_ends = end_day[mm]
            si_change = si_changes[mm]
            clim_change = clim_changes[mm]
            
            x6, y6 = monthly_plot_siclim(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, clim_change, si_change, scale=True)
            
            axes6[era][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
        
       
        xlims6 = [[-4,4],[-3.5,4],[-9,10]]
        axes6[era] = plot_final_scatter( axes6[era], x6, y6, xlims6)
       
            
    fig6.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig6.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
    axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
               edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%% sic

def sic_scatter(start, end, years, bbox, si_lon,si_lat, daily_map, STORMAREA=True):
    time_range = [start, start+timedelta(days=3),start+timedelta(days=7), start+timedelta(days=14)]
    time_range=[start, end]
    
    
    tseries = []
    for DT in time_range:
        if DT.month==2 and DT.day==29:
            DT = datetime(start.year, 3, 1)
    
        daily_range = fx.daterange(datetime(2010, 1,1), datetime(2010,12,31), dt=24)
        dt_ind1 = [i for i,dt in enumerate(daily_range) if (dt.month == DT.month and dt.day == DT.day)][0]
        
        try: sic_dt = daily_map[dt_ind1, :,:]
        except IndexError as ie: 
            print(DT, dt_ind1)
            print(ie)
            raise IndexError

        if not STORMAREA:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tseries.append( np.nanmean(sic_dt) )
        else:
            si_lon2 = np.where(si_lon<0, si_lon+360, si_lon)
            inside_points = fx.find_points_in_contour(bbox, si_lon2, si_lat)
            bbox_dt = np.where(inside_points, np.nan, sic_dt)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tseries.append( np.nanmean(bbox_dt) )
                
    # return [np.nan]*7+[tseries[0]]+[np.nan]*2+[tseries[1]]+[np.nan]*3+[tseries[2]]+[np.nan]*6 +[tseries[3]]
    return [np.nanmean(tseries)]*22

#%%% scatter (sic)
print('Starting SIC Calculations (may take a few min...)')
savepath = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/monthly_miz/'
ice_fname = '/Users/mundi/Desktop/seaice/'

for loc_ind, loc in enumerate(hemi_names):
    print('-->',loc)
    if loc_ind==0: continue
    path1 = root_paths[loc_ind]
    fig6, axes6 = scatter_plot_setup(xlabel='Mean Sea Ice Concentration',
                                     ylabel='Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle='\n\n'+loc,
                                     hline=0, vline=None)
    if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
    elif loc_ind==1: _, si_lon, si_lat, = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)
    
    for era, years in enumerate(decades):
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
            
        for year in years:
            print(year)
            if loc_ind == 0:
                daily_map = np.load(savepath+str(year)+'_map.npy') 
            elif loc_ind == 1:
                daily_map = np.load(savepath+str(year)+'_map_sh.npy')   
                
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines([year], path1+'census/', path1+'area/', path1+'seaice/')
            storm_areas = fx.get_storm_bbox([year], path1+'census/', path1+'area/', path1+'contours/')
        
            for mi, mm in enumerate(list(lines.keys())):
                monthly_starts = start_day[mm]
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                areas = storm_areas[mm]
                
                # organize SICs
                SIC = []
                for start, end, bbox in zip(monthly_starts, monthly_ends, areas):
                    sic1 = sic_scatter(start, end, years, bbox, si_lon, si_lat, daily_map)
                    SIC.append(sic1)
                    # print(sic1)
                
                # scatter
                x6, y6 = monthly_plot_mean(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                         clim_change, si_change, SIC)
                
                axes6[era][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
                
        for i, ax6 in enumerate(axes6[era]):
            if len(x6[i])!=0:
                m, b, r, p, se = linregress(x6[i], y6[i])
                x = np.linspace(np.min(x6[i]), np.max(x6[i]), 100)
                ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
                
                ax6a = ax6.twinx()
                ax6a.sharey(ax6)  
                ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
                ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
                ax6a.axis('off');

#%% seasonal: paired scatter

def seasonal_si(season, mi, monthly_starts, monthly_ends, VAR_X, VAR_Y, scale=True):
   
    for start, end, var_x, var_y in zip(monthly_starts, monthly_ends, VAR_X, VAR_Y):
        ### storm duration
        try:
            del_x = var_x[7+(end-start).days+1] - var_x[7]
            del_y = (var_y[7+(end-start).days+1] - var_y[7])
            if scale: del_y = del_y - del_x
        except IndexError:
            del_x = var_x[-1] - var_x[7]
            del_y = (var_y[-1] - var_y[7]) 
            if scale: del_y = del_y - del_x
        if not (np.isnan(del_x) or np.isnan(del_y)): season[0][start.year].append(del_y/1e5)
        
        ### one week
        del_x = var_x[14] - var_x[7]
        del_y = (var_y[14] - var_y[7])
        if scale:del_y = del_y - del_x
        if not (np.isnan(del_x) or np.isnan(del_y)): season[1][start.year].append(del_y/1e5)
        
        ### three weeks
        del_x = var_x[-1] - var_x[0]
        del_y = (var_y[-1] - var_y[0])
        if scale: del_y = del_y - del_x
        if not (np.isnan(del_x) or np.isnan(del_y)): season[2][start.year].append(del_y/1e5)
        
    return season

#### data loop

for loc_ind, location in enumerate(hemi_names[0:1]):

    fig5, axes5 = scatter_plot_setup(xlabel='Season 1: '+'Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     ylabel='Season 2: '+'Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location+' Seasonal Comparison',
                                     hline=0, vline=0)
        
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        
        path1 = root_paths[loc_ind]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
    
        ### scatter
        summer = {}
        spring = {}
        winter = {}
        fall = {}
        for i in np.arange(0,len(axes5[era])):
            spring[i]={y:[] for y in years}
            summer[i]={y:[] for y in years}
            fall[i]={y:[] for y in years}
            winter[i]={y:[] for y in years}
        
        for mi, mm in enumerate(list(lines.keys())):
            monthly_starts = start_day[mm]
            monthly_ends = end_day[mm]
            si_change = si_changes[mm]
            clim_change = clim_changes[mm]
            
            if mm in [6,7,8]: summer = seasonal_si(summer, mi, monthly_starts, monthly_ends, clim_change, si_change, scale=True)
            elif mm in [3,4,5]: spring = seasonal_si(spring, mi, monthly_starts, monthly_ends, clim_change, si_change, scale=True)
            elif mm in [12,1,2]: winter = seasonal_si(winter, mi, monthly_starts, monthly_ends, clim_change, si_change, scale=True)
            elif mm in [9,10,11]: fall = seasonal_si(fall, mi, monthly_starts, monthly_ends, clim_change, si_change, scale=True)
    
        #### PLOT
        x5, y5 = [],[] ###!!! change here
        for i, ax5 in enumerate(axes5[era]):
            x5.append([np.nanmean(fall[i][yy]) for yy in years])
            y5.append([np.nanmean(winter[i][yy]) for yy in years])
            ax5.plot(x5[i], y5[i], marker='o', lw=0)
            
        axes5[era] = plot_final_scatter(axes5[era], x5, y5, xlims=None)
            
    fig5.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig5.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
    axes5[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
               edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))


#%% storm location: latitude vs ice impact (no)
if False:

    for loc_ind, location in enumerate(hemi_names):
    
        fig6, axes6 = scatter_plot_setup(xlabel=r'Storm Latitude',
                                         ylabel='Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                         suptitle=location+': Latitudes', 
                                         hline=0, vline=None
                                         )
            
        #### DATA 
        for era, years in enumerate(decades):
            yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
            
            path1 = root_paths[loc_ind]
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
            latitudes = fx.indiv_locations(years, path1+'census/', path1+'area/')[-1]
        
            ### scatter
            x6, y6 = {},{}
            for i in np.arange(0,len(axes6[era])):
                x6[i]=[]
                y6[i]=[]
            
            for mi, mm in enumerate(list(lines.keys())):
                monthly_starts = start_day[mm]
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                lats = latitudes[mm]
                
                x6, y6 = monthly_plot_si(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                         clim_change, si_change, lats)
                
                axes6[era][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
            
           
            # xlims6 = [[-4,4],[-3.5,4],[-9,10]]
            axes6[era] = plot_final_scatter( axes6[era], x6, y6, xlims=None)
           
                
        fig6.subplots_adjust(hspace=0.3)
        ycoords = [0.7, 0.3]
        for dname, dec, loc in zip(decade_names, decades, ycoords):
            fig6.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                     va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
        axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
                   edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%% *SST differences*

def sst_scatter(start, end, years, bbox, si_lon):
    ### computes sst difference ###
    
    if start.month==2 and start.day==29:
        start = datetime(start.year, 3, 1)
    if end.month==2 and end.day==29:
        end = datetime(end.year, 3, 1)
    
    yr_ind1 = [i for i,dt in enumerate(years) if dt==start.year][0]
    
    daily_range = fx.daterange(datetime(2010, 1,1), datetime(2010,12,31), dt=24)
    dt_ind1 = [i for i,dt in enumerate(daily_range) if (dt.month == start.month and dt.day == start.day)][0]
    
    try: sst_start = monthly_data[yr_ind1, dt_ind1, :,:]
    except IndexError as ie: 
        print(start, yr_ind1, dt_ind1)
        print(ie)
        raise IndexError
    
    try:
        yr_ind2 = [i for i,dt in enumerate(years) if dt==end.year][0]
        dt_ind2 = [i for i,dt in enumerate(daily_range) if (dt.month == end.month and dt.day == end.day)][0]
        sst_end = monthly_data[yr_ind2, dt_ind2, :,:]
    except IndexError:
        sst_end = monthly_data[yr_ind1, -1, :,:]
        
    if not STORM_AREA:
        return np.nanmean( np.nanmean(sst_end) - np.nanmean(sst_start) )
    else:
        si_lon2 = np.where(si_lon<0, si_lon+360, si_lon)
        inside_points = fx.find_points_in_contour(bbox, si_lon2, si_lat)
        start_bbox = np.where(inside_points, np.nan, sst_start)
        end_bbox = np.where(inside_points, np.nan, sst_end)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return np.nanmean( np.nanmean(end_bbox) - np.nanmean(start_bbox) )

#%%% sst data
SHIFT_MIZ = True
lat_shift = 2
only_lower_lat = True
STORM_AREA = True

ice_fname =  '/Users/mundi/Desktop/seaice/'

hemi_map = []
for loc_ind, loc in enumerate(hemi_names):

    savepath = root_paths[loc_ind]+'sst/'
    
    if not SHIFT_MIZ:
        fname =  '_map'
    else: 
        if not only_lower_lat:
            fname = '_shifted_map_'+str(lat_shift)
        else:
            fname = '_lowered_map_'+str(lat_shift)
        
    decade_map = []
    for years in decades:
        year_map = []
        for year in years:
            try:
                daily_map = np.load(savepath+str(year)+fname+'.npy')  
            except FileNotFoundError:
                print('Missing SST:', savepath+str(year)+fname+'.npy')
            year_map.append(daily_map)
        decade_map.append(year_map)
    hemi_map.append(np.array(decade_map))

        
#%%% sst scatter plot
for loc_ind, location in enumerate(hemi_names):

    fig6, axes6 = scatter_plot_setup(xlabel=r'SST ($^\circ$C)',
                                     ylabel='Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location+': SSTs', 
                                     hline=0, vline=0
                                     )
    decade_map = hemi_map[loc_ind]
    
    if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
    elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)
        
    #### DATA 
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        
        path1 = root_paths[loc_ind]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        storm_areas = fx.get_storm_bbox(years, path1+'census/', path1+'area/', path1+'contours/')
    
        monthly_data = decade_map[era]

    
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
        
        for mi, mm in enumerate(list(lines.keys())):
            monthly_starts = start_day[mm]
            monthly_ends = end_day[mm]
            si_change = si_changes[mm]
            clim_change = clim_changes[mm]
            areas = storm_areas[mm]
            
            # organize SSTs
            SST = []
            for start, end, bbox in zip(monthly_starts, monthly_ends, areas):
                diff = sst_scatter(start, end+timedelta(days=7), years, bbox, si_lon) #!!! storm difference
                SST.append(diff)
            
            # scatter
            x6, y6 = monthly_plot_si(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                     clim_change, si_change, SST)
            
            axes6[era][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
        
       
        # xlims6 = [[-4,4],[-3.5,4],[-9,10]]
        axes6[era] = plot_final_scatter( axes6[era], x6, y6, xlims=None)
       
            
    fig6.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig6.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
    axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
               edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%%% mean sst

def sst_start(start, years, bbox, si_lon):
    ### computes sst difference ###
    
    if start.month==2 and start.day==29:
        start = datetime(start.year, 3, 1)
   
    yr_ind1 = [i for i,dt in enumerate(years) if dt==start.year][0]
    
    daily_range = fx.daterange(datetime(2010, 1,1), datetime(2010,12,31), dt=24)
    dt_ind1 = [i for i,dt in enumerate(daily_range) if (dt.month == start.month and dt.day == start.day)][0]
    
    try: sst_start = monthly_data[yr_ind1, dt_ind1, :,:]
    except IndexError as ie: 
        print(start, yr_ind1, dt_ind1)
        print(ie)
        raise IndexError
    
    if not STORM_AREA:
        return np.nanmean(sst_start)
    else:
        si_lon2 = np.where(si_lon<0, si_lon+360, si_lon)
        inside_points = fx.find_points_in_contour(bbox, si_lon2, si_lat)
        start_bbox = np.where(inside_points, np.nan, sst_start)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return np.nanmean(start_bbox) 

for loc_ind, location in enumerate(hemi_names):

    fig6, axes6 = scatter_plot_setup(xlabel=r'Mean SST ($^\circ$C)',
                                     ylabel='Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location+': SSTs', 
                                     hline=0, vline=0
                                     )
    decade_map = hemi_map[loc_ind]
    
    if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
    elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)
        
    #### DATA 
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        
        path1 = root_paths[loc_ind]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        storm_areas = fx.get_storm_bbox(years, path1+'census/', path1+'area/', path1+'contours/')
    
        monthly_data = decade_map[era]

    
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
        
        for mi, mm in enumerate(list(lines.keys())):
            monthly_starts = start_day[mm]
            monthly_ends = end_day[mm]
            si_change = si_changes[mm]
            clim_change = clim_changes[mm]
            areas = storm_areas[mm]
            
            # organize SSTs
            SST = []
            for start, end, bbox in zip(monthly_starts, monthly_ends, areas):
                diff = sst_start(start, years, bbox, si_lon) 
                SST.append(diff)
            
            # scatter
            x6, y6 = monthly_plot_si(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                     clim_change, si_change, SST)
            
            axes6[era][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
        
       
        # xlims6 = [[-4,4],[-3.5,4],[-9,10]]
        axes6[era] = plot_final_scatter( axes6[era], x6, y6, xlims=None)
       
            
    fig6.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig6.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
    axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
               edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%% winds
# wnames = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 'vline':0}
# wnames = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 'vline':0}
wnames = {'idx':0, 'ylims':[-2,14], 'name':'Total', 'vline':None}

for loc_ind, location in enumerate(hemi_names):
    
    fig6, axes6 = scatter_plot_setup(xlabel='Mean '+wnames['name'] +r' Wind Speed (m s$^{-1}$)',
                                     ylabel=r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location+': '+wnames['name']+' Wind', 
                                     vline=wnames['vline'], hline=0
                                     )
    
    for era, years in enumerate(decades):
        path1 = root_paths[loc_ind]
        
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for year in years:
            ### WIND
            if loc_ind ==0: 
                wind_series = fx.wind_lines([year],path1+'census/', path1+'area/', path1+'wind/')[-1]
            elif loc_ind ==1:
                wind_series = fx.sh_winds([year], path1+'census/', path1+'area/', path1+'wind/')
                
            # get time series   
            for mi, mm in enumerate(list(lines.keys())):
                monthly_starts = start_day[mm]
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                # try:
                #     wind_s = np.array(wind_series[mm])[:,:,wnames['idx']]
                # except: continue
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[mm] if len(arr)==22]  
                    if len(windies)==0: continue
                    wind_s = np.array(windies)[:,:,wnames['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames['idx']] for arr in wind_series[mm]]  
                    if len(windies)==0: continue
                    wind_s = np.array(windies)[:,]
                
                x6, y6 = monthly_plot_mean(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                         clim_change, si_change, wind_s)
                
            
        titles6 = ['Storm Duration', 'One Week Post-Storm',
                   '3-Week Analysis']
        for i, ax6 in enumerate(axes6[era]):
            if len(x6[i])!=0:
                m, b, r, p, se = linregress(x6[i], y6[i])
                x = np.linspace(np.min(x6[i]), np.max(x6[i]), 100)
                ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
                ax6.set_title(titles6[i])
                
                ax6a = ax6.twinx()
                ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
                ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
                ax6a.axis('off');
        
    fig6.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig6.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
        
        
for mi, mm in enumerate(months): axes6[0][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
       edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))


#%% waves

for loc_ind, location in enumerate(hemi_names):
    
    fig7, axes6 = scatter_plot_setup(xlabel='Mean Wave Height (m)',
                                     ylabel=r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location, 
                                     vline=None, hline=0
                                     )
    
    for era, years in enumerate(decades):
        path1 = root_paths[loc_ind]
        
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for year in years:
            ### WIND
            swh_series = fx.era_lines([year],path1+'census/', path1+'area/', path1+'swh/', 'swh')[-1]
                
            # get time series   
            for mi, mm in enumerate(list(lines.keys())):
                monthly_starts = start_day[mm]
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                
                swh_lines = [np.array(s) for s in swh_series[mm] if len(s)==22]

                                    # mean, si(difference)
                x6, y6 = monthly_plot_mean(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                         clim_change, si_change, swh_lines)
                
            
        titles6 = ['Storm Duration', 'One Week Post-Storm',
                   '3-Week Analysis']
        for i, ax6 in enumerate(axes6[era]):
            if len(x6[i])!=0:
                m, b, r, p, se = linregress(x6[i], y6[i])
                x = np.linspace(np.min(x6[i]), np.max(x6[i]), 100)
                ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
                ax6.set_title(titles6[i])
                
                ax6a = ax6.twinx()
                ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
                ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
                ax6a.axis('off');
        
    fig7.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig7.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
        
        
for mi, mm in enumerate(months): axes6[0][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
       edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))


#%% air temperature 

for loc_ind, location in enumerate(hemi_names):
    
    fig7, axes6 = scatter_plot_setup(xlabel='Temperature (C)',
                                     ylabel=r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location, 
                                     vline=0, hline=0
                                     )
    
    for era, years in enumerate(decades):
        path1 = root_paths[loc_ind]
        
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for year in years:
            ### WIND
            swh_series = fx.era_lines([year],path1+'census/', path1+'area/', path1+'t2m/', 't2m')[-1]
                
            # get time series   
            for mi, mm in enumerate(list(lines.keys())):
                monthly_starts = start_day[mm]
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                
                swh_lines = [np.array(s)-273.15 for s in swh_series[mm] if len(s)==22]

                                    # mean, si(difference)
                x6, y6 = monthly_plot_si(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                         clim_change, si_change, swh_lines)
                
            
        titles6 = ['Storm Duration', 'One Week Post-Storm',
                   '3-Week Analysis']
        for i, ax6 in enumerate(axes6[era]):
            if len(x6[i])!=0:
                m, b, r, p, se = linregress(x6[i], y6[i])
                x = np.linspace(np.min(x6[i]), np.max(x6[i]), 100)
                ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
                ax6.set_title(titles6[i])
                
                ax6a = ax6.twinx()
                ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
                ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
                ax6a.axis('off');
        
    fig7.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig7.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
        
        
for mi, mm in enumerate(months): axes6[0][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
       edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%% ice motion

mnames0 = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 
          'xlims':[[-4.5,3],[-1,3]], 'vline':0}
mnames1 = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 
          'xlims':[[-4.5,2],[-1,2]], 'vline':0}
mnames2 = {'idx':0, 'ylims':[-2,14], 'name':'Total', 
          'xlims':[[4,10],[6,10]], 'vline':None}

mnames = [mnames0, mnames1, mnames2][1]

for loc_ind, location in enumerate(hemi_names):
    
    fig7, axes6 = scatter_plot_setup(xlabel=mnames['name']+' Ice Motion '+r'(cm s$^{-1}$)',
                                     ylabel=r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location+': Ice Motion', 
                                     vline=0, hline=0
                                     )
    
    for era, years in enumerate(decades):
        path1 = root_paths[loc_ind]
        
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for year in years:
            ### motions series
            motion_series = fx.storm_ice_motion(years, path1+'census/', path1+'area/', path1+'icemotion/')[-1]
                
                
            # get time series   
            for mi, mm in enumerate(list(lines.keys())):
                monthly_starts = start_day[mm]
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                
                moves = [arr for arr in motion_series[mm] if len(arr)==22]  
                im_lines = np.array(moves)[:,:,mnames['idx']]

                                    # mean, si(difference)
                x6, y6 = monthly_plot_si(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                         clim_change, si_change, im_lines)
                
            
        titles6 = ['Storm Duration', 'One Week Post-Storm',
                   '3-Week Analysis']
        for i, ax6 in enumerate(axes6[era]):
            if len(x6[i])!=0:
                m, b, r, p, se = linregress(x6[i], y6[i])
                x = np.linspace(np.min(x6[i]), np.max(x6[i]), 100)
                ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
                ax6.set_title(titles6[i])
                
                ax6a = ax6.twinx()
                ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
                ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
                ax6a.axis('off');
        
    fig7.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig7.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
        
        
for mi, mm in enumerate(months): axes6[0][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
       edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%% storm pressure

for loc_ind, location in enumerate(hemi_names):
    
    fig7, axes6 = scatter_plot_setup(xlabel='Storm Pressure (hPa)',
                                     ylabel=r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location, 
                                     vline=None, hline=0
                                     )
    
    for era, years in enumerate(decades):
        path1 = root_paths[loc_ind]
        
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        _,_, pressures = fx.storm_pressure(years, path1+'census/', path1+'area/')
            
        for year in years:
            # get time series   
            for mi, mm in enumerate(list(lines.keys())):
                monthly_starts = start_day[mm]
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                slp = pressures[mm]

                                    # mean, si(difference)
                x6, y6 = monthly_plot_si(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                         clim_change, si_change, slp)
                
            
        titles6 = ['Storm Duration', 'One Week Post-Storm',
                   '3-Week Analysis']
        for i, ax6 in enumerate(axes6[era]):
            if len(x6[i])!=0:
                m, b, r, p, se = linregress(x6[i], y6[i])
                x = np.linspace(np.min(x6[i]), np.max(x6[i]), 100)
                ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
                ax6.set_title(titles6[i])
                
                ax6a = ax6.twinx()
                ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
                ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
                ax6a.axis('off');
        
    fig7.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig7.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
        
        
for mi, mm in enumerate(months): axes6[0][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
       edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%% Sea Ice Concentration

#%%% sic data​ (change over storms)

def get_storm_areas(year, root_path, lon, lat):
    from glob import glob
    ice_lims = [20,80]
    
    census_file = root_path+'census/'+'census_'+str(year)+'.csv'
    [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
    
    timing_grid = []
    for xx in range(0,len(startdate)):
        timing_grid.append((startdate[xx], enddate[xx]))
    storm_ranges = []
    for startdt, enddt in timing_grid:
        storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
    
    # open ice area
    ds_area = xr.open_dataset(root_path+'area/' + str(year) +'_area.nc')
    ice_sorter = ds_area['ice_area80'].values
    box_area = ds_area['box_area'].values
    
    try:
        ds = xr.open_dataset(root_path+'seaice/' + str(year) + '_seaice' + '.nc')
        ds.close()
    except:
        print('- skip: '+root_path+'seaice/' + str(year) + '_seaice' + '.nc')
        return [], []
    
    storm_areas = []
    stormstr_prev=''
    for storm_num, strm in enumerate(storm_ranges):
        # remove storms that don't interact with the ice
        ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
        if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
            storm_areas.append([])
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
        ncname = stormstr + '_contours' +'.nc'
        try:
            if not dupe: ###!!!
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
            storm_areas.append([])
            continue
        
        with fx.HidePrint(): bbox_edges = fx.get_bbox_edges(all_contours) 
        
        storm_areas.append( fx.find_points_in_contour(bbox_edges, lon, lat) )
        
    return storm_areas, storm_ranges
         

#%%%% data
savepath = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/monthly_miz/'
ice_fname = '/Users/mundi/Desktop/seaice/'
_, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
_, si_lon_sh, si_lat_sh = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)

area = True

if area: fname = '_area'
else: fname = '_miz'

miz_loc = []
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    print(loc, 'storm sic')
    miz_decade = []
    for years in decades:
        print('*', str(years[0])+'-'+str(years[-1]))
        
        miz_month = {mm:[] for mm in np.arange(1,12+1)}
        for year in years:
            
            ### load storm census
            path1 = root_paths[loc_ind]
            if loc_ind==0: lon1, lat1 = si_lon, si_lat 
            elif loc_ind==1: lon1, lat1 = si_lon_sh, si_lat_sh
            storm_areas, storm_ranges = get_storm_areas(year, root_paths[loc_ind], lon1, lat1)
            
            ### daily miz
            if loc_ind == 0:
                daily_map = np.load(savepath+str(year)+'_map.npy') 
            elif loc_ind == 1:
                daily_map = np.load(savepath+str(year)+'_map_sh.npy')   
                  
            for area, storm in zip(storm_areas, storm_ranges):
                if np.shape(area) !=  np.shape(lon1): continue
                
                dt_start = datetime(2010, storm[0].month, storm[0].day if (storm[0].day!=29 and storm[0].month!=2) else 28 ) # avoid leap stuff
                start_ind = dt_start.timetuple().tm_yday
                dt_end = datetime(2010, storm[-1].month, storm[-1].day if (storm[-1].day!=29 and storm[-1].month!=2) else 28 )
                end_ind = dt_end.timetuple().tm_yday
                dt_wk = datetime(2010, (storm[0]+timedelta(days=14)).month, 
                                 (storm[0]+timedelta(days=14)).day if ((storm[0]+timedelta(days=14)).day!=29 and (storm[0]+timedelta(days=14)).month==2) else 28 )
                wk_ind = dt_wk.timetuple().tm_yday
                
                mizs = [np.ma.masked_array( daily_map[i-1], mask=area) for i in [start_ind, end_ind, wk_ind]]
                miz_month[storm[0].month].append(mizs)
                
        miz_decade.append(miz_month)
    miz_loc.append(miz_decade)
    
#%%% individual storm scatter

sic_value = 0
sic_titles = {0:'Starting SIC', 1:'Ending SIC', 2:'2-Weeks Post-Storm', 
              3:'Storm Difference', 4:'2-Week Difference'}

for loc_ind, location in enumerate(hemi_names):
    
    fig7, axes6 = scatter_plot_setup(xlabel='SIC',
                                     ylabel=r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location+': '+sic_titles[sic_value], 
                                     vline=None if sic_value not in [3,4] else 0, 
                                     hline=0
                                     )
    
    for era, years in enumerate(decades):
        path1 = root_paths[loc_ind]
        
        ### scatter
        x6, y6 = {},{}
        for i in np.arange(0,len(axes6[era])):
            x6[i]=[]
            y6[i]=[]
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for year in years:
            # get time series   
            for mi, mm in enumerate(list(lines.keys())):
                monthly_starts = start_day[mm]
                if len(monthly_starts)<5:continue
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                
                if sic_value in [0,1,2]:
                    monthly_miz = np.array(miz_loc[loc_ind][era][mm])[:,sic_value] 
                elif sic_value==3:
                    end_miz = np.array(miz_loc[loc_ind][era][mm])[:,1] 
                    start_miz = np.array(miz_loc[loc_ind][era][mm])[:,0] 
                    monthly_miz = end_miz - start_miz
                elif sic_value==4:
                    end_miz = np.array(miz_loc[loc_ind][era][mm])[:,2] 
                    start_miz = np.array(miz_loc[loc_ind][era][mm])[:,0] 
                    monthly_miz = end_miz - start_miz
                else: raise IndexError('check sic_value: no matching values')
                    
                miz_sic = np.nanmean(monthly_miz, axis=(-2,-1))

                                    # mean, si(difference)
                x6, y6 = monthly_plot_si(axes6[era], x6, y6, mi, monthly_starts, monthly_ends, 
                                         clim_change, si_change, miz_sic)
                
            
        titles6 = ['Storm Duration', 'One Week Post-Storm',
                   '3-Week Analysis']
        for i, ax6 in enumerate(axes6[era]):
            if len(x6[i])!=0:
                m, b, r, p, se = linregress(x6[i], y6[i])
                x = np.linspace(np.min(x6[i]), np.max(x6[i]), 100)
                ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
                ax6.set_title(titles6[i])
                
                ax6a = ax6.twinx()
                ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
                ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
                ax6a.axis('off');
        
    fig7.subplots_adjust(hspace=0.3)
    ycoords = [0.7, 0.3]
    for dname, dec, loc in zip(decade_names, decades, ycoords):
        fig7.text(0.0725, loc, dname+ '('+str(dec[0])+' - '+str(dec[-1])+')', 
                 va='center', ha='center', fontsize=12, rotation=90, fontweight='bold')
        
        
for mi, mm in enumerate(months): axes6[0][0].plot([],[], color = colors[mi], label=calendar.month_name[mm])
axes6[0][0].legend(loc='lower right', ncol=3, handletextpad=0.5, handlelength=1,
       edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0))

#%%% monthly mean

sic_value = 0
sic_titles = {0:'Starting SIC', 1:'Ending SIC', 2:'2-Weeks Post-Storm', 
              3:'Storm Difference', 4:'2-Week Difference'}

fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex='row')
fig.suptitle('Sea Ice Concentration'+': '+sic_titles[sic_value])
for era in [0,1]: axes[-1][era].set_xlabel('SIC')
for loc in [0,1]: axes[loc][0].set_ylabel('Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)')  


for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]

    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))

        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
         
        xs, ys = [], []
        for mi, mm in enumerate(list(lines.keys())):
            
           clim_change = np.nanmean(clim_changes[mm], axis=0)
           clim1 = (clim_change[-1] - clim_change[0])/1e5
            
           si_change = np.nanmean(si_changes[mm], axis=0)
           si1 = (si_change[-1] - si_change[0])/1e5
           
           if sic_value in [0,1,2]:
               monthly_miz = np.array(miz_loc[loc_ind][era][mm])[:,sic_value] 
           elif sic_value==3:
               end_miz = np.array(miz_loc[loc_ind][era][mm])[:,1] 
               start_miz = np.array(miz_loc[loc_ind][era][mm])[:,0] 
               monthly_miz = end_miz - start_miz
           elif sic_value==4:
               end_miz = np.array(miz_loc[loc_ind][era][mm])[:,2] 
               start_miz = np.array(miz_loc[loc_ind][era][mm])[:,0] 
               monthly_miz = end_miz - start_miz
           else: raise IndexError('check sic_value: no matching values')
           mean_sic = np.nanmean(monthly_miz)
           
           xs.append(mean_sic)
           ys.append(si1- clim1)
               
           axes[loc_ind][era].plot(mean_sic, si1-clim1, color=colors[mi],
                                   marker='o', markersize=10, zorder=-20)
          
        m, b, r, p, se = linregress(xs, ys)
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        axes[loc_ind][era].axhline(0, lw=0.55, color='gray', ls=':')
        if sic_value in [3,4]:
            axes[loc_ind][era].axvline(0, lw=0.55, color='gray', ls=':')
        
        ax6a = axes[loc_ind][era].twinx()
        ax6a.sharey(axes[loc_ind][era])  
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');
        


#%% end
