#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 14:00:40 2025
scatter_pie1.py

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

#%% functions

#%%% plotting
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
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

    return ax_row

#%%% data

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



#%% sea ice area

#%%% all cyclones

for loc_ind, location in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    fig6, axes6 = scatter_plot_setup(xlabel=r'$\Delta$ Sea Ice Climatology ($\times10^5$ km$^2$)',
                                     ylabel='Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)',
                                     suptitle=location,
                                     hline=0, vline=0)
            
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        
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
    
#%% simple scatter (sst trend, si change)

fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex=True)

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
            
           # si_change = np.nanmean(np.array(si_changes[mm])-np.array(clim_changes[mm]), axis=0)
           si_change = np.nanmean(si_changes[mm], axis=0)
           si1 = (si_change[-1] - si_change[0])/1e5
           
           xs.append(clim1)
           ys.append(si1)
           
           axes[loc_ind][era].plot(clim1, si1, color=colors[mi],
                                   marker='o', markersize=8)
                
        m, b, r, p, se = linregress(xs, ys)
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        ax6a = axes[loc_ind][era].twinx()
        ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.axvline(0, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

for era in [0,1]: axes[-1][era].set_xlabel(r'$\Delta$ Sea Ice Climatology ($\times10^5$ km$^2$)')
for loc in [0,1]: axes[loc][0].set_ylabel('Relative ' + r'$\Delta$ Sea Ice ($\times10^5$ km$^2$)')  

#%% ---    
#%% PIE
#%% ---
#%% functions
def draw_pie(dist, 
             xpos, 
             ypos, 
             colors, 
             size=500,
             ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()
    
    colors = iter(colors)

    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])
        
        nstorms = np.sum(dist)

        ax.scatter([xpos], [ypos], marker=xy, s=size+(nstorms), fc=next(colors), zorder=500)

    return ax

def label_quadrants(ax, xs, ys, colors):
    
    x_pos = np.nanmean([0, np.nanmax(xs)])
    if x_pos <= 0 : x_pos=100
    x_neg = np.nanmean([0, np.nanmin(xs)])
    if x_neg >= 0 : x_neg=-100
    y_pos = np.nanmean([0, np.nanmax(ys)])
    if y_pos <= 0 : y_pos=100
    y_neg = np.nanmean([0, np.nanmin(ys)])
    if y_neg >= 0 : y_neg=-100
    
    qlabels = ['I', 'II', 'III', 'IV']
    qlocs = [(x_pos, y_pos), (x_neg, y_pos), (x_neg, y_neg), (x_pos, y_neg)]
    qcolors = iter(colors)
    
    for QL, LOC in zip(qlabels, qlocs):
        ax.text(LOC[0], LOC[1], QL, 
                horizontalalignment='center',
                verticalalignment='center',
                fontsize = 20,
                zorder=-100, 
                color = next(qcolors)
                )
    
    

#%%% info
quad_colors = ['#f7f7f7','#cccccc','#969696','#525252']
# ['#d9d9d9','#bdbdbd','#969696','#636363']

def quadrants(x_list, y_list):
    quads = {q:0 for q in [0,1,2,3]}
    for x, y in zip(x_list, y_list):
        if x>0 and y>0: quads[0]+=1
        elif x<0 and y>0: quads[1]+=1
        elif x<0 and y<0: quads[2]+=1
        elif x>0 and y<0: quads[3]+=1
    return quads

#%% plot

fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex='row')
fig.suptitle('Sea Ice')
for era in [0,1]: axes[-1][era].set_xlabel(r'$\Delta$ Sea Ice Climatology ($\times10^5$ km$^2$)')
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
           
           xs.append(clim1)
           ys.append(si1)
           
           quads = {q:0 for q in [0,1,2,3]}
           for si, clim in zip(si_changes[mm], clim_changes[mm]):
               si0 = si[-1] - si[0]
               clim0 = clim[-1] - clim[0] ###!!! there was a typo here
               
               if clim0>=0 and si0>=0: quads[0]+=1
               elif clim0<0 and si0>=0: quads[1]+=1
               elif clim0<0 and si0<0: quads[2]+=1
               elif clim0>=0 and si0<0: quads[3]+=1
              
               
           axes[loc_ind][era].plot(clim1, si1, color=colors[mi],
                                   marker='o', markersize=30, zorder=-20)
           axes[loc_ind][era].plot(clim1, si1, color='k',
                                   marker='o', markersize=24, zorder=-19)
          
           dist = [quads[q] for q in range(4)]
           draw_pie(dist, clim1, si1, ax = axes[loc_ind][era],
                    colors = quad_colors)
                
        m, b, r, p, se = linregress(xs, ys)
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        axes[loc_ind][era].axhline(0, lw=0.55, color='gray', ls=':')
        axes[loc_ind][era].axvline(0, lw=0.55, color='gray', ls=':')
        
        ax6a = axes[loc_ind][era].twinx()
        ax6a.sharey(axes[loc_ind][era])  
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');
        
        label_quadrants(axes[loc_ind][era], xs, ys, colors=quad_colors)


#%% ----
#%% mean wind and mean ice impact

# wnames2 = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 'xlims':[[-4.5,3],[-1,3]]}
wnames2 = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 'xlims':[[-4.5,2],[-1,2]]}
# wnames2 = {'idx':0, 'ylims':[-2,14], 'name':'Total', 'xlims':[[4,10],[6,10]]}

figs, axs = plt.subplots(2,2, figsize=(8,8), sharey=True)
figs.suptitle('Monthly '+wnames2['name']+' Wind')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        axs[loc_ind][era].set_title(loc+': '+yr_title)
        axs[loc_ind][era].set_xlim(wnames2['xlims'][loc_ind])
        axs[loc_ind][era].axhline(0, lw=0.55, color='gray', ls=':')
        axs[loc_ind][era].axvline(0, lw=0.55, color='gray', ls=':')
        
        if loc_ind ==0:
            wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
        elif loc_ind==1: 
            wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        xs, ys = [], []
        for mi, month in enumerate(months):
            
            try: 
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[month] if len(arr)==22]  
                    wind_lines = np.array(windies)[:,:,wnames2['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames2['idx']] for arr in wind_series[month]]  
                    wind_lines = np.array(windies)[:,]
            except IndexError: print('-->', month); continue
            
            if len(wind_lines)<10: 
                continue
            
            mean_line = np.nanmean(wind_lines, axis=0)
            std_line = np.nanstd(wind_lines, axis=0)

            wind_mean = np.nanmean(mean_line[7:10])
    
            axs[loc_ind][era].plot(wind_mean, mean_lines[mi][-1], 
                                   marker='o', markersize=22, color=colors[mi])
            xs.append(wind_mean)
            ys.append(mean_lines[mi][-1])
            
            # pie
            quads = quadrants([np.nanmean(wl[7:10]) for wl in wind_lines],
                              [si_l[-1] for si_l in lines[month]])
            dist = [quads[q] for q in range(4)]
            draw_pie(dist, wind_mean, mean_lines[mi][-1], ax = axs[loc_ind][era],
                     colors = quad_colors, size=350)
        label_quadrants(axs[loc_ind][era], xs, ys, quad_colors)
            
        #### lines of best fit? correlations?
        m, b, r, p, se = linregress(xs,ys)
        ax6 =  axs[loc_ind][era]
        x = np.linspace(np.nanmin(xs), np.nanmax(xs), 50)
        ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
        ax6a = ax6.twinx()
        ax6a.sharey(ax6)  
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

# final plot
for era, years in enumerate(decades): axs[-1][era].set_xlabel('Wind Speed')
for loc_ind, loc in enumerate(hemi_names): axs[loc_ind][0].set_ylabel('Normalized Change in MIZ Ice Area')


#%% sst trend
ice_fname =  '/Users/mundi/Desktop/seaice/'
SHIFT_MIZ = True
lat_shift = 2
only_lower_lat = True

#%%% data

hemi_map = []

for loc_ind, loc in enumerate(hemi_names):

    if loc_ind==0: _, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
    elif loc_ind==1: _, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)
    si_points = [(x,y) for x,y in zip(si_lon.flatten(), si_lat.flatten())]
    
    savepath = root_paths[loc_ind]+'sst/'
    
    decade_map = []
    
    if not SHIFT_MIZ:
        fname =  '_map'
    else: 
        if not only_lower_lat:
            fname = '_shifted_map_'+str(lat_shift)
        else:
            fname = '_lowered_map_'+str(lat_shift)
        
        
    print('Loading SST Data: '+loc)
    for years in decades:
    
        year_map = []
    
        for year in years:
            daily_map = np.load(savepath+str(year)+fname+'.npy')  
            year_map.append(daily_map)
    
        decade_map.append(year_map)
    
    decade_map = np.array(decade_map)
    
    hemi_map.append(decade_map)

print('Done')


#%%% plot
fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex=True)
fig.suptitle('SST')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    decade_map = hemi_map[loc_ind]

    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))

        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
         
        monthly_data = decade_map[era]
        xs, ys = [], []
        for mi, mm in enumerate(list(lines.keys())):
            ml = mean_lines[mi][-1]
            
            # organize SSTs
            SSTs = np.load(path1+'sst/'+'tseries_'+str(era)+'-'+str(mm)+'.npy')
                
            ### plot
            if len(SSTs)>10:
                mean_line = np.nanmean(SSTs, axis=0)
                slope, b, r, p, se = linregress(xxx[7:14], mean_line[7:14])
                
                axes[loc_ind][era].plot(slope, ml, color=colors[mi],
                                        marker='o', markersize=22)
                xs.append(slope); ys.append(ml)
                
            # pie                     ###!!! wk slope? full tseries slope?
            quads = quadrants([linregress(xxx[7:14],sst_line[7:14])[0] for sst_line in SSTs],
                              [si_l[-1] for si_l in lines[mm]])
            dist = [quads[q] for q in range(4)]
            draw_pie(dist, slope, ml, ax = axes[loc_ind][era],
                     colors = quad_colors, size=350)
        label_quadrants(axes[loc_ind][era], xs, ys, quad_colors)
                
        m, b, r, p, se = linregress(xs, ys)
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        ax6a = axes[loc_ind][era].twinx()
        ax6a.sharey(axes[loc_ind][era])  
        ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.axvline(0, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

for era in [0,1]: axes[-1][era].set_xlabel(r'Mean Change in SST ($^\circ$C day$^{-1}$)')
for loc in [0,1]: axes[loc][0].set_ylabel('Normalized Change in MIZ Ice Area')

#%%% plot: non-norm
fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex=True)
fig.suptitle('SST')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    decade_map = hemi_map[loc_ind]

    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))

        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
         
        monthly_data = decade_map[era]
        xs, ys = [], []
        for mi, mm in enumerate(list(lines.keys())):
            ml = mean_lines[mi][-1]
            
            # organize SSTs
            SSTs = np.load(path1+'sst/'+'tseries_'+str(era)+'-'+str(mm)+'.npy')
                
            ### plot
            if len(SSTs)>10:
                mean_line = np.nanmean(SSTs, axis=0)
                slope, b, r, p, se = linregress(xxx[7:14], mean_line[7:14])
                
                axes[loc_ind][era].plot(slope, ml, color=colors[mi],
                                        marker='o', markersize=22)
                xs.append(slope); ys.append(ml)
                
            # pie                     ###!!! wk slope? full tseries slope?
            quads = quadrants([linregress(xxx[7:14],sst_line[7:14])[0] for sst_line in SSTs],
                              [si_l[-1] for si_l in lines[mm]])
            dist = [quads[q] for q in range(4)]
            draw_pie(dist, slope, ml, ax = axes[loc_ind][era],
                     colors = quad_colors, size=350)
        label_quadrants(axes[loc_ind][era], xs, ys, quad_colors)
                
        m, b, r, p, se = linregress(xs, ys)
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        ax6a = axes[loc_ind][era].twinx()
        ax6a.sharey(axes[loc_ind][era])  
        ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.axvline(0, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

for era in [0,1]: axes[-1][era].set_xlabel(r'Mean Change in SST ($^\circ$C day$^{-1}$)')
for loc in [0,1]: axes[loc][0].set_ylabel('Normalized Change in MIZ Ice Area')

#%% air temps

#### data
mean_swh_lines = {}
all_swh_lines = {}
for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    mean_swh_lines[loc_ind] = {}
    all_swh_lines[loc_ind] = {}
    for era, years in enumerate(decades):
        mean_swh_lines[loc_ind][era] = {mm:[] for mm in months}
        all_swh_lines[loc_ind][era] = {mm:[] for mm in months}
        try:
            swh_series = fx.era_lines(years, path1+'census/', path1+'area/', path1+'t2m/', 't2m')[-1]
        except FileNotFoundError as fnfe:
            print(fnfe)
            
        for mi, month in enumerate(months):
            swh_lines = [np.array(s)-273.15 for s in swh_series[month] if len(s)==22]
            
            if len(swh_lines)<10: continue
            
            mean_line = np.nanmean(swh_lines, axis=0)
            std_line = np.nanstd(swh_lines, axis=0)
            
            mean_swh_lines[loc_ind][era][month] = mean_line
            all_swh_lines[loc_ind][era][month] += swh_lines
            

#### plot
fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex=True)
fig.suptitle('Air Temperature')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    mean_lines_mm = mean_swh_lines[loc_ind][era]
    
    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        xs, ys = [], []
        for month in months:
            if len(mean_lines_mm[month])==0:continue
            swh_ml = mean_lines_mm[month]
                
            swh_dev = swh_ml[10] - swh_ml[7]
            
            ml = mean_lines[month-1][-1]

            axes[loc_ind][era].plot(swh_dev, ml, color=colors[month-1],
                                    marker='o', markersize=22)
            xs.append(swh_dev); ys.append(ml)
            
            # pie
            quads = quadrants([swh1[10] - swh1[7] for swh1 in all_swh_lines[loc_ind][era][month]],
                              [si_l[-1] for si_l in lines[month]])
            dist = [quads[q] for q in range(4)]
            draw_pie(dist, swh_dev, ml, ax = axes[loc_ind][era],
                     colors = quad_colors, size=350)
        label_quadrants(axes[loc_ind][era], xs, ys, quad_colors)
                
            
        mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
        m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        ax6a = axes[loc_ind][era].twinx()
        ax6a.sharey(axes[loc_ind][era])  
        ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.axvline(0, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

for era in [0,1]: axes[-1][era].set_xlabel('Change in Temperature (C)')
for loc in [0,1]: axes[loc][0].set_ylabel('Normalized Change in MIZ Ice Area')


#%% waves

#### data
mean_swh_lines = {}
all_swh_lines = {}
for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    mean_swh_lines[loc_ind] = {}
    all_swh_lines[loc_ind] = {}
    for era, years in enumerate(decades):
        mean_swh_lines[loc_ind][era] = {mm:[] for mm in months}
        all_swh_lines[loc_ind][era] = {mm:[] for mm in months}
        try:
            swh_series = fx.era_lines(years, path1+'census/', path1+'area/', path1+'swh/', 'swh')[-1]
        except FileNotFoundError as fnfe:
            print(fnfe)
            
        for mi, month in enumerate(months):
            swh_lines = [np.array(s)-273.15 for s in swh_series[month] if len(s)==22]
            
            if len(swh_lines)<10: continue
            
            mean_line = np.nanmean(swh_lines, axis=0)
            std_line = np.nanstd(swh_lines, axis=0)
            
            mean_swh_lines[loc_ind][era][month] = mean_line
            all_swh_lines[loc_ind][era][month] += swh_lines
            

#### plot
fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex=True)
fig.suptitle('Wave Heights')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    mean_lines_mm = mean_swh_lines[loc_ind][era]
    
    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        xs, ys = [], []
        for month in months:
            if len(mean_lines_mm[month])==0:continue
            swh_ml = mean_lines_mm[month]
                
            swh_dev = swh_ml[10] - swh_ml[7]
            
            ml = mean_lines[month-1][-1]

            axes[loc_ind][era].plot(swh_dev, ml, color=colors[month-1],
                                    marker='o', markersize=22)
            xs.append(swh_dev); ys.append(ml)
            
            # pie
            quads = quadrants([swh1[10] - swh1[7] for swh1 in all_swh_lines[loc_ind][era][month]],
                              [si_l[-1] for si_l in lines[month]])
            dist = [quads[q] for q in range(4)]
            draw_pie(dist, swh_dev, ml, ax = axes[loc_ind][era],
                     colors = quad_colors, size=350)
        label_quadrants(axes[loc_ind][era], xs, ys, quad_colors)
                
            
        mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
        m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
        x = np.linspace(np.min(xs), np.max(xs), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
         
        ax6a = axes[loc_ind][era].twinx()
        ax6a.sharey(axes[loc_ind][era])   
        ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.axvline(0, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

for era in [0,1]: axes[-1][era].set_xlabel('Change in Wave Height (m)')
for loc in [0,1]: axes[loc][0].set_ylabel('Normalized Change in MIZ Ice Area')

#%% end

