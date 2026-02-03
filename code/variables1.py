#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 2025
variables1.py

organize figure story for paper

@author: mundi
"""
#%% imports and filename
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import time as timeIN
from glob import glob

import functions as fx

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

census_name = 'census_'
contour_name = '_contours'
si_name = '_seaice'

ice_fname = '/Users/mundi/Desktop/seaice/'

decades = [np.arange(1982, 1992), np.arange(2010,2020)]
decade_names = ['Early Satellite Era ', 'Present Day ']

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

ice_lims = [20,80]

pie_scatter = False

#%% data organization
fs = 14

alph = ['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k']

month_colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

hemi_colors = ['#01665e','#8c510a']
hemi_names= ['Arctic', 'Antarctic']

decade_colors = ['#2166ac', '#b2182b']
shade_colors = ['#67a9cf', '#ef8a62']

quad_colors = ['#f7f7f7','#cccccc','#969696','#525252']

#%%% functions

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

        ax.scatter([xpos], [ypos], marker=xy, s=size+(nstorms), 
                   fc=next(colors), alpha=0.66, zorder=500)

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
    
def quadrants(x_list, y_list):
    quads = {q:0 for q in [0,1,2,3]}
    for x, y in zip(x_list, y_list):
        if x>0 and y>0: quads[0]+=1
        elif x<0 and y>0: quads[1]+=1
        elif x<0 and y<0: quads[2]+=1
        elif x>0 and y<0: quads[3]+=1
    return quads


#%% si impact: spaghetti fans

def ytext_scale(idx, li, yi):
    yadd=0
    idx += 1 # actual month
    
    if li==0: # arctic
        if yi == 0: # left (past)
            if idx==11: yadd=0.01 # november
            elif idx==2: yadd=-0.02 # feb
            elif idx==12: yadd=0.01 # dec
        elif yi ==1: # right (recent)
            if idx==4: yadd=0.05 # apr
            elif idx==3:yadd=0.01 # march
            elif idx==1: yadd=-0.05 # jan
            elif idx == 5: yadd=-0.05 # may
            elif idx == 8: yadd=-0.05 # aug
            elif idx==9: yadd=-0.05 # sep
            elif idx==11:yadd=0.015 # nov
    elif li==1: #antarctic
        if yi == 0: # left (past)
            if idx==6: yadd=0.02
            elif idx==7: yadd=-0.02
        elif yi ==1: # right (recent)
            if idx==8: yadd=0.02 #aug
            elif idx==9: yadd=-0.02 # sep
    
    return yadd

#%%% make plot

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(16,14), sharey=True)
plt.subplots_adjust(wspace=0.25)
fig2.suptitle('\n\nMonthly Mean Change in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.axvline(0, ls=':', color='gray', lw=1)
    ax2.set_xlim(-7,14)
    ax2.set_xticks(xxx)
    ax2.set_xticklabels(xlabels, minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Days Since Storm Start',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Normalized Relative Change in Ice Area',fontsize=fs+1)
    
#### data and plot
for li, loc in enumerate(hemi_names):
    for yi, years in enumerate(decades):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        ### monthly mean
        axes2[li][yi].set_title(hemi_names[li]+': '+yr_title, fontsize=fs+2)
        for idx, ml in enumerate(mean_lines):
            axes2[li][yi].plot(xxx, ml, color = month_colors[idx], lw=2, label = calendar.month_name[idx+1])
            yadd = ytext_scale(idx, li, yi)
            mlabel = calendar.month_abbr[idx+1] #+' ('+str(len(lines[idx+1]))+')'
            axes2[li][yi].text(xxx[-1]+0.25, ml[-1]+yadd, mlabel, fontsize=fs-1) 
        

axes2[0][1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                  edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                  bbox_to_anchor=(1.5,0.33))

#%% WINDS

wnames = [{'idx':2, 'ylims':[-6,6], 'name':'Meridional', 'xlims':[[-4.5,2],[-1,2]]},
          {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 'xlims':[[-4.5,3],[-1,3]]}]

fig = plt.figure(figsize=(16, 8))
subfigs = fig.subfigures(1, 2, wspace=-0.075)

for wi, wnames2 in enumerate(wnames):
    subfigs[wi].suptitle('\nMonthly '+wnames2['name']+' Wind')
    axs = subfigs[wi].subplots(2, 2, sharey=True)

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
                                       marker='o', color=month_colors[mi],
                                       markersize=22 if pie_scatter else 10)
                xs.append(wind_mean)
                ys.append(mean_lines[mi][-1])
                
                # pie
                if pie_scatter:
                    quads = quadrants([np.nanmean(wl[7:10]) for wl in wind_lines],
                                      [si_l[-1] for si_l in lines[month]])
                    dist = [quads[q] for q in range(4)]
                    draw_pie(dist, wind_mean, mean_lines[mi][-1], ax = axs[loc_ind][era],
                             colors = quad_colors, size=350)
            if pie_scatter: label_quadrants(axs[loc_ind][era], xs, ys, quad_colors)
                
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
fig.suptitle('\n'+'Air Temperature')

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

            axes[loc_ind][era].plot(swh_dev, ml, color=month_colors[month-1],
                                    marker='o', markersize=22 if pie_scatter else 10)
            xs.append(swh_dev); ys.append(ml)
            
            # pie
            if pie_scatter:
                quads = quadrants([swh1[10] - swh1[7] for swh1 in all_swh_lines[loc_ind][era][month]],
                                  [si_l[-1] for si_l in lines[month]])
                dist = [quads[q] for q in range(4)]
                draw_pie(dist, swh_dev, ml, ax = axes[loc_ind][era],
                         colors = quad_colors, size=350)
        if pie_scatter: label_quadrants(axes[loc_ind][era], xs, ys, quad_colors)
                
            
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

#%% sst trend

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
fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=True, sharex='row')
fig.suptitle('\n'+'SST')

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
                
                axes[loc_ind][era].plot(slope, ml, color=month_colors[mi],
                                        marker='o', markersize=22 if pie_scatter else 10)
                xs.append(slope); ys.append(ml)
                
            # pie        
            if pie_scatter:             ###!!! wk slope? full tseries slope?
                quads = quadrants([linregress(xxx[7:14],sst_line[7:14])[0] for sst_line in SSTs],
                                  [si_l[-1] for si_l in lines[mm]])
                dist = [quads[q] for q in range(4)]
                draw_pie(dist, slope, ml, ax = axes[loc_ind][era],
                         colors = quad_colors, size=350)
        if pie_scatter: label_quadrants(axes[loc_ind][era], xs, ys, quad_colors)
                
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

#%% --> ocean profiles

#%%% data loop

ocn_data_miz = []

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]

    ### LOAD DEPTH
    with xr.open_dataset(path1+'depth.nc') as dds:
        DEPTH = dds['depth'].values
        
    data_all = {mm: [] for mm in np.arange(1,12+1)}
    data_miz = {mm: [] for mm in np.arange(1,12+1)}
    
    ### STORM LOOP
    for year in years:
        # CENSUS: storm info
        census_file = path1+'census/'+'census_'+str(year)+'.csv'
        [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
        
        # open ice area
        ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
        ice_area80 = ds_area['ice_area80'].values
        ice_area15 = ds_area['ice_area15'].values
        box_area = ds_area['box_area'].values
        ice_sorter = ice_area80
        
        for sn, start in enumerate(startdate): 
            month = start.month
            
            ice_frac = ice_sorter[sn]*100/box_area[sn]             
            if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                continue
    
            ## LOAD FILES
            start_str = start.strftime('%Y-%m-%d')
            
            try:
                prof_miz = np.load(path1+'ocn_prof/'+start_str+'_'+str(sn)+'_miz.npy')
                data_miz[month].append(np.nanmean(prof_miz, axis=0))
                # data_miz[month].append(prof_miz[-1] - prof_miz[0])
            except FileNotFoundError:
                print('missing:', start_str)
                continue
            
    ocn_data_miz.append(data_miz)
    
#%% plot all profiles

# either same plots (skinny storm lines and thick mean)
# or side-by-side subfigures

def make_grid_plot(nrow, ncol, title=''):
    fontsize=12
    fig, axes_all = plt.subplots(nrow, ncol, sharex=True, sharey=True, figsize=(3*nrow,5*ncol))
    fig.suptitle('\n\n'+title, fontsize=fontsize+2)
    alph = iter(list(string.ascii_lowercase))
    
    for ax in axes_all[:,0]:
        ax.set_ylabel('Depth (m)', fontsize=fontsize)
        ax.set_ylim([0.5,80])
        ax.invert_yaxis()
        # ax.yaxis.tick_right()
    for ax in axes_all[-1,:]:
        ax.set_xlabel('Temperature'+r' ($^\circ$C)', fontsize=fontsize-1)
    for axl in axes_all:
        for ax1 in axl:
            ax1.set_xlim([-2.5,5])
            # ax1.axvline(0, ls='-', color='k', lw=0.55)
            ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
            
            ax1.text(0.0225, 0.915, next(alph), transform=ax1.transAxes, fontsize=fontsize, 
                    bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                          'edgecolor':'k', 'lw':0.75},zorder=50)
    return fig, axes_all


fig, axes_all = make_grid_plot(4,3, title='MIZ '+yr_title)

for loc_ind, loc in enumerate(hemi_names):
    
    data_miz = ocn_data_miz[loc_ind]
    
    yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
    
    
    for mi, ax in enumerate(axes_all.flatten()): 
        ax.set_title(calendar.month_name[mi+1])
        
        data_mm = data_miz[mi+1]
        
        if len(data_miz[mi+1])>10:
            for line in data_mm:
                ax.plot(line, DEPTH, lw=0.25, color=hemi_colors[loc_ind])
       
            ax.plot(np.nanmean(data_mm, axis=0), DEPTH, lw=2.5, color=hemi_colors[loc_ind])


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
fig, axes = plt.subplots(2,2,figsize=(10,10))

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))
        
        mean_lines_mm = mean_swh_lines[loc_ind][era]
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        xs, ys = [], []
        for month in months:
            if len(mean_lines_mm[month])==0:
                xs.append(np.nan)
                ys.append(np.nan)
                continue
            swh_ml = mean_lines_mm[month]
            swh_dev = np.nanmean(swh_ml[7:10]) - np.nanmean(swh_ml)
            print(loc, era, month, round(swh_dev,3))
            
            ml = mean_lines[month-1][-1]

            if loc_ind==1: axes[loc_ind][era].plot(swh_dev, np.abs(ml), 
                                        color=month_colors[month-1],
                                        marker='o' if ml>0 else 'v', markersize=10)
            else:
                axes[loc_ind][era].plot(swh_dev, ml, color=month_colors[month-1],
                                        marker='o', markersize=10)
                
            xs.append(swh_dev); ys.append(ml if loc_ind==0 else np.abs(ml))
            
            # pie
            if pie_scatter:
                quads = quadrants([swh1[10] - swh1[7] for swh1 in all_swh_lines[loc_ind][era][month]],
                                  [si_l[-1] for si_l in lines[month]])
                dist = [quads[q] for q in range(4)]
                draw_pie(dist, swh_dev, ml, ax = axes[loc_ind][era],
                         colors = quad_colors, size=350)
            if pie_scatter: label_quadrants(axes[loc_ind][era], xs, ys, quad_colors)
            
            
        ax6a = axes[loc_ind][era].twinx()
        mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
        m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
        x = np.linspace(np.min(np.array(xs)[mask]), np.max(np.array(xs)[mask]), 50)
        axes[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=3, zorder=-20, alpha=0.66)
         
        # ax6a.axhline(0, lw=0.55, color='gray', ls=':')
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

for era in [0,1]: axes[-1][era].set_xlabel('Change in Wave Height (m)')
for loc in [0,1]: axes[loc][0].set_ylabel('Normalized Change in MIZ Ice Area')


#%% end
