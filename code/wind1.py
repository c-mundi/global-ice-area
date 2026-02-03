#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 2025
wind1.py

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


#%% data organization

fs = 14

alph = ['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k']

month_colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

# hemi_colors = ['b', 'r']
hemi_colors = ['#01665e','#8c510a']

hemi_names= ['Arctic', 'Antarctic']

decade_colors = ['#2166ac', '#b2182b']
shade_colors = ['#67a9cf', '#ef8a62']

widx = 0

wnames0 = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 
          'xlims':[[-4.5,3],[-1,3]], 'vline':0}
wnames1 = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 
          'xlims':[[-4.5,2],[-1,2]], 'vline':0}
wnames2 = {'idx':0, 'ylims':[-2,14], 'name':'Total', 
          'xlims':[[4,10],[6,10]], 'vline':None}

wnames = [wnames0, wnames1, wnames2][widx]

#%% monthly plot

cmap = cmo.balance #cmo.curl
norm = plt.Normalize(vmin=wnames['ylims'][0], vmax=wnames['ylims'][1])

def make_grid_plot(nrow, ncol, title='', ylabel=''):
    fontsize=12
    fig, axes_all = plt.subplots(nrow, ncol, sharex=True, sharey=True, figsize=(4*nrow,3.25*ncol))
    fig.suptitle(title, fontsize=fontsize+1)
    alph = iter(list(string.ascii_lowercase))
    
    for ax in axes_all[:,0]:
        ax.set_ylabel(ylabel, fontsize=fontsize)
    for ax in axes_all[-1,:]:
        ax.set_xticks(xxx)
        ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
        ax.set_xlabel('Days Since Storm Start', fontsize=fontsize-1)
    for axl in axes_all:
        for ax1 in axl:
            ax1.set_xlim(-7,14)
            ax1.axhline(0, ls='-', color='k', lw=1)
            ax1.axvline(0, ls='-', color='k', lw=1)
            ax1.set_xticks(xxx)
            ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
            
            ax1.text(0.0225, 0.915, next(alph), transform=ax1.transAxes, fontsize=fontsize, 
                    bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                          'edgecolor':'k', 'lw':0.75},zorder=50)

    return fig, axes_all

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        fig, axes_all = make_grid_plot(4,3, title='\n\n'+loc+' '+wnames['name']+' Wind'+' '+yr_title,
                                       ylabel=r'[m s$^{-1}$]')
        
        if loc_ind ==0:
            try:
                wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
            except FileNotFoundError as fnfe:
                pass
                # print(fnfe)
                # continue
        else: 
            try:
                wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
            except FileNotFoundError as fnfe:
                print(fnfe)
                pass
            
        axes_w = axes_all.flatten()
        for mi, month in enumerate(months):
            
            try: 
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[month] if len(arr)==22]  
                    wind_lines = np.array(windies)[:,:,wnames['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames['idx']] for arr in wind_series[month]]  
                    wind_lines = np.array(windies)[:,]
            except IndexError: print('-->', month); continue
            
            if len(wind_lines)<10: continue
            
            for line in wind_lines: axes_w[mi].plot(xxx, line, lw=0.55, color='gray')
                
            mean_line = np.nanmean(wind_lines, axis=0)
            std_line = np.nanstd(wind_lines, axis=0)
            
            axes_w[mi].plot(xxx, mean_line, lw=1.75, color='maroon')
            
            wind_mean = np.nanmean(mean_line[7:10]) 
            axes_w[mi].axvspan(-7,14, color=cmap(norm(wind_mean)), alpha=0.33, zorder=-500)
            
            axes_w[mi].set_title(calendar.month_name[month]+' (n='+str(len(wind_lines))+')')
            axes_w[mi].set_ylim(wnames['ylims'])
            
#%%% scale: monthly ice loss per unit wind?

            # zonal, merid, total
wnames9 = [wnames0, wnames1, wnames2][1]

cmap1 = cmo.balance #cmo.curl
norm1 = plt.Normalize(vmin=-4, vmax=4)

loc_means = []
for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    era_means = []
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        fig, axes_all = make_grid_plot(4,3, title='\nScaled Sea Ice Response\n'+loc+' '+wnames9['name']+' Wind'+' '+yr_title,
                                       ylabel='sea ice change\n'+r'per m s$^{-1}$')
        
        # load winds
        if loc_ind==0: wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
        elif loc_ind==1: wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
            
        # load sea ice
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        axes_w = axes_all.flatten()
        monthly_means = []
        for mi, month in enumerate(months):
            try: 
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[month] if len(arr)==22]  
                    wind_lines = np.array(windies)[:,:,wnames9['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames9['idx']] for arr in wind_series[month]]  
                    wind_lines = np.array(windies)[:,]
            except IndexError: print('-->', month); continue
            
            if len(wind_lines)<10: monthly_means.append(np.nan); continue
        
            # seaice_lines = lines[month]
            seaice_lines = (np.array(si_changes[month])-np.array(clim_changes[month]))/1e5
            
            scaled = []
            for wind, si in zip(wind_lines, seaice_lines): 
                
                axes_w[mi].plot(xxx, si/wind, lw=0.55, color='gray')
                scaled.append(si/wind)
                
            mean_line = np.nanmean(scaled, axis=0)
            std_line = np.nanstd(scaled, axis=0)
            
            axes_w[mi].plot(xxx, mean_line, lw=1.75, color='maroon')
            
            wind_mean = np.nanmean(mean_line[7:14]) 
            axes_w[mi].axvspan(-7,14, color=cmap1(norm1(wind_mean)), alpha=0.33, zorder=-500)
            monthly_means.append(wind_mean)
            
            axes_w[mi].set_title(calendar.month_name[month]+' (n='+str(len(wind_lines))+')')
            axes_w[mi].set_ylim(wnames9['ylims'])
    
        era_means.append(monthly_means)
    loc_means.append(era_means)
    
#%%% scaled monthly scatter/correlation

figs, axs = plt.subplots(2,2, figsize=(8,8), sharey='row', sharex='row')
figs.suptitle('Monthly Means Scaled By '+wnames9['name']+' Wind')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    for era, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        axs[loc_ind][era].set_title(loc+': '+yr_title)
        # axs[loc_ind][era].set_xlim(wnames['xlims'][loc_ind])
        axs[loc_ind][era].axhline(0, lw=0.55, color='gray', ls=':')
        axs[loc_ind][era].axvline(0, lw=0.55, color='gray', ls=':')
        
        monthly_mean_data = loc_means[loc_ind][era]
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for mi, mm in enumerate(months):
            axs[loc_ind][era].plot(monthly_mean_data[mi], mean_lines[mi][-1], 
                                   marker='o', markersize=10, color=month_colors[mi])

    
#%%  simple scatter: mean wind and mean ice impact

figs, axs = plt.subplots(2,2, figsize=(8,8), sharey='row')
figs.suptitle('Monthly '+wnames['name']+' Wind')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        axs[loc_ind][era].set_title(loc+': '+yr_title)
        axs[loc_ind][era].set_xlim(wnames['xlims'][loc_ind])
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
            marker='o'
            
            try: 
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[month] if len(arr)==22]  
                    wind_lines = np.array(windies)[:,:,wnames['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames['idx']] for arr in wind_series[month]]  
                    wind_lines = np.array(windies)[:,]
            except IndexError: print('-->', month); continue
            
            if len(wind_lines)<10: 
                continue
            
            mean_wind_line = np.nanmean(wind_lines, axis=0)
            wind_mean = np.nanmean(mean_wind_line[7:10])
            
            mean_impact = mean_lines[mi][-1]
            
            if widx==0 and loc_ind==1 and mean_impact<0:
                mean_impact = -mean_impact
                marker='v'
            
            axs[loc_ind][era].plot(wind_mean, mean_impact, 
                                   marker=marker, markersize=10, color=month_colors[mi])
            xs.append(wind_mean)
            ys.append(mean_impact)
            
        #### lines of best fit? correlations?
        m, b, r, p, se = linregress(xs,ys)
        ax6 =  axs[loc_ind][era]
        x = np.linspace(np.nanmin(xs), np.nanmax(xs), 50)
        ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
        ax6a = ax6.twinx()
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

# final plot
for era, years in enumerate(decades): axs[-1][era].set_xlabel('Wind Speed')
for loc_ind, loc in enumerate(hemi_names): axs[loc_ind][0].set_ylabel('Normalized Change in MIZ Ice Area')

#%% monthly scatters
MYMONTH = 6

figs, axs = plt.subplots(2,2, figsize=(8,8), sharey='row')
figs.suptitle('Monthly '+wnames['name']+' Wind: '+str(MYMONTH))

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        axs[loc_ind][era].set_title(loc+': '+yr_title)
        # axs[loc_ind][era].set_xlim(wnames['xlims'][loc_ind]) 
        axs[loc_ind][era].axhline(0, lw=0.55, color='gray', ls=':')
        axs[loc_ind][era].axvline(0, lw=0.55, color='gray', ls=':')
        
        if loc_ind==0:
            wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
        elif loc_ind==1: 
            wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        xs, ys = [], []
            
        if loc_ind==0: 
            windies = [arr for arr in wind_series[MYMONTH] if len(arr)==22]  
            wind_lines = np.array(windies)[:,:,wnames['idx']]
        elif loc_ind==1: 
            windies = [arr[wnames['idx']] for arr in wind_series[MYMONTH]]  
            wind_lines = np.array(windies)[:,]
        
        if len(wind_lines)<10: 
            continue
        
        for w_line, si_line in zip(wind_lines, lines[MYMONTH]):
            wind_mean = np.nanmean(w_line[7:10])
            mean_impact = si_line[-1]
        
            axs[loc_ind][era].plot(wind_mean, mean_impact, 
                                   marker='o', markersize=10, 
                                   color=month_colors[MYMONTH-1])
            xs.append(wind_mean)
            ys.append(mean_impact)
            
        #### lines of best fit? correlations?
        mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
        m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
        ax6 =  axs[loc_ind][era]
        x = np.linspace(np.nanmin(xs), np.nanmax(xs), 50)
        ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
        ax6a = ax6.twinx()
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

# final plot
for era, years in enumerate(decades): axs[-1][era].set_xlabel('Wind Speed')
for loc_ind, loc in enumerate(hemi_names): axs[loc_ind][0].set_ylabel('Normalized Change in MIZ Ice Area')


#%% total wind timeseries
## annual cycle?

fig9, axes9 = plt.subplots(2, 1, figsize=(10,6), sharex=True, sharey=True)

for loc_ind, location in enumerate(hemi_names[0:1]):
    
    for era, years in enumerate(decades):
        path1 = root_paths[loc_ind]
        tseries = []
        mm_tseries = {mm:[] for mm in months}
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for year in years:
            ### WIND
            if loc_ind ==0: 
                wind_series = fx.wind_lines([year],path1+'census/', path1+'area/', path1+'wind/')[-1]
            elif loc_ind ==1:
                wind_series = fx.sh_winds([year], path1+'census/', path1+'area/', path1+'wind/')
                
            # get time series   
            for mi, mm in enumerate(months):
                monthly_starts = start_day[mm]
                monthly_ends = end_day[mm]
                si_change = si_changes[mm]
                clim_change = clim_changes[mm]
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[mm] if len(arr)==22]  
                    if len(windies)==0: 
                        tseries.append(np.nan)
                        continue
                    wind_s = np.array(windies)[:,:,wnames['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames['idx']] for arr in wind_series[mm]]  
                    if len(windies)==0: 
                        tseries.append(np.nan)
                        continue
                    wind_s = np.array(windies)[:,]
                    
                means = [np.nanmean(line[7:14]) for line in wind_s]
                tseries.append(np.nanmean(means))
                mm_tseries[mm].append(np.nanmean(means))
                
        axes9[era].plot(tseries)
        for xv in np.arange(0,120,12): axes9[era].axvline(xv, lw=0.55, color='gray')
        axes9[-1].set_xlabel('Month Count')
        axes9[era].set_ylabel(wnames['name']+' Winds')
        axes9[era].set_title(str(years[0])+'-'+str(years[-1]))
        
        axes9[era].plot(np.arange(-12,0), [np.nanmean(mm_tseries[m1]) for m1 in months])
        
#%% sort storm timeseries by mean wind (all storms)

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    fig, axes_w = make_grid_plot(2,2, title=loc+' '+wnames['name']+' Wind',
                                   ylabel='Normalized Change in\nMIZ Ice Area')
    
    for era, years in enumerate(decades):
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        
        if loc_ind ==0:
            try:
                wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
            except FileNotFoundError as fnfe:
                print(fnfe)
        else: 
            try:
                wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
            except FileNotFoundError as fnfe:
                print(fnfe)
                
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
             fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        nly, sly = 0,0
        northerlies, southerlies = [],[]
        for mi, month in enumerate(months):
            
            try: 
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[month] if len(arr)==22]  
                    wind_lines = np.array(windies)[:,:,wnames['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames['idx']] for arr in wind_series[month]]  
                    wind_lines = np.array(windies)[:,]
            except IndexError: print('-->', month); continue
            
            if len(wind_lines)<10: continue
            
           
            for wline, iline in zip(wind_lines, lines[month]):
                if np.nanmean(wline[7:10])<-0.01:
                    axes_w[0][era].plot(xxx, iline, lw=0.55, color='gray')
                    northerlies.append(iline)
                    nly+=1
                elif np.nanmean(wline[7:10])>0.01:
                    axes_w[1][era].plot(xxx, iline, lw=0.55, color='gray')
                    southerlies.append(iline)
                    sly+=1
            
            
        axes_w[0][era].plot(xxx, np.nanmean(northerlies, axis=0), lw=2, color='maroon')
        axes_w[0][era].set_title('Northerly Winds'+' (n='+str(nly)+') ' +yr_title)
        axes_w[0][era].set_ylim([-1,1])
        
        axes_w[1][era].plot(xxx, np.nanmean(southerlies, axis=0), lw=2, color='maroon')
        axes_w[1][era].set_title('Southerly Winds'+' (n='+str(sly)+') ' +yr_title)
        axes_w[1][era].set_ylim([-1,1])

#%% sort storm timeseries by mean wind (monthly groups)

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    if loc_ind==0: 
        mgroups = [[9,10,11,12],[1,2,3],[4,5,6,7,8]]
        mcolors = ['b','g','m']
    elif loc_ind==1: 
        mgroups = [[4,5,6,7],[8,9,10],[11,12]]
        mcolors = ['m','g','b']
    
    for mgroup, mcolor, mtitle in zip(mgroups, mcolors, ['Increasing','Neutral','Decreasing']):
        
        fig, axes_w = make_grid_plot(2,2, title=loc+' '+wnames['name']+' Wind ('+mtitle+')',
                                       ylabel='Normalized Change\nin MIZ Ice Area')
    
        for era, years in enumerate(decades):
            yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
            
            if loc_ind ==0:
                try:
                    wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
                except FileNotFoundError as fnfe:
                    print(fnfe)
            else: 
                try:
                    wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
                except FileNotFoundError as fnfe:
                    print(fnfe)
                    
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                 fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
                
            nly, sly = 0,0
            northerlies, southerlies = [],[]
            for mi, month in enumerate(mgroup):
                
                try: 
                    if loc_ind==0: 
                        windies = [arr for arr in wind_series[month] if len(arr)==22]  
                        wind_lines = np.array(windies)[:,:,wnames['idx']]
                    elif loc_ind==1: 
                        windies = [arr[wnames['idx']] for arr in wind_series[month]]  
                        wind_lines = np.array(windies)[:,]
                except IndexError: print('-->', month); continue
                
                if len(wind_lines)<10: continue
                
               
                for wline, iline in zip(wind_lines, lines[month]):
                    if np.nanmean(wline[7:10])<-0.01:
                        axes_w[0][era].plot(xxx, iline, lw=0.55, color='gray')
                        northerlies.append(iline)
                        nly+=1
                    elif np.nanmean(wline[7:10])>0.01:
                        axes_w[1][era].plot(xxx, iline, lw=0.55, color='gray')
                        southerlies.append(iline)
                        sly+=1
                
                
            axes_w[0][era].plot(xxx, np.nanmean(northerlies, axis=0), lw=2, color=mcolor)
            axes_w[0][era].set_title('Northerly Winds'+' (n='+str(nly)+') ' +yr_title)
            axes_w[0][era].set_ylim([-1,1])
            
            axes_w[1][era].plot(xxx, np.nanmean(southerlies, axis=0), lw=2, color=mcolor)
            axes_w[1][era].set_title('Southerly Winds'+' (n='+str(sly)+') ' +yr_title)
            axes_w[1][era].set_ylim([-1,1])


        
#%% end// data check

if False: # check missing data
    for loc_ind, loc in enumerate(hemi_names[0:1]):
        path1 = root_paths[loc_ind]
        for years in decades:
            for year in years:
                for month in months:
                   try: np.load(path1+'wind/'+'winds_'+str(year)+'-'+str(month)+'.npy', allow_pickle=True)
                   except FileNotFoundError: print(str(year)+'-'+str(month))
                   
#%% -------------
#%% ROTATED WINDS

#%%% monthly timeseries
cmap1 = cmo.balance #cmo.curl
norm1 = plt.Normalize(vmin=-5, vmax=5)

for loc_ind, loc in enumerate(hemi_names): 
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        if era==0: continue ##!!!
        
        yr_title = '('+str(years[0])+'-'+str(years[-1])+')'
        fig, axes_all = make_grid_plot(4,3, title='\n\n'+loc+' Rotated Winds'+' '+yr_title,
                                       ylabel=r'[m s$^{-1}$]')
        
        wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind_rot/', NAME='wind_rot')[-1]
            
        axes_w = axes_all.flatten()
        for mi, month in enumerate(months):
            try: 
                windies = [arr for arr in wind_series[month] if len(arr)==22]  
                wind_lines = -np.array(windies) ##!!! negative??
            except IndexError: print('-->', month); continue
            
            # if len(wind_lines)<10: continue ###!!!
            if len(wind_lines)==0: continue
            
            for line in wind_lines: axes_w[mi].plot(xxx, line, lw=0.55, color='gray')
                
            mean_line = np.nanmean(wind_lines, axis=0)
            std_line = np.nanstd(wind_lines, axis=0)
            
            axes_w[mi].plot(xxx, mean_line, lw=1.75, color='maroon')
            
            wind_mean = np.nanmean(mean_line[7:10]) 
            axes_w[mi].axvspan(-7,14, color=cmap(norm(wind_mean)), alpha=0.33, zorder=-500)
            
            axes_w[mi].set_title(calendar.month_name[month]+' (n='+str(len(wind_lines))+')')
            axes_w[mi].set_ylim([norm1.vmin, norm1.vmax])
     
#%%% simple scatter
figs, axs = plt.subplots(2,2, figsize=(8,8), sharey=True)
figs.suptitle('Monthly Rotated Winds')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        axs[loc_ind][era].set_title(loc+': '+yr_title)
        # axs[loc_ind][era].set_xlim(wnames2['xlims'][loc_ind])
        axs[loc_ind][era].axhline(0, lw=0.55, color='gray', ls=':')
        axs[loc_ind][era].axvline(0, lw=0.55, color='gray', ls=':')
        
        wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind_rot/', NAME='wind_rot')[-1]
            
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        xs, ys = [], []
        for mi, month in enumerate(months):
            
            windies = [arr for arr in wind_series[month] if len(arr)==22]  
            wind_lines = -np.array(windies) ###!!! neg sign make more sense? +/+, -/- (pos slope)
            
            # if len(wind_lines)<10: continue ###!!!!
            if len(wind_lines)==0: continue
            
            mean_line = np.nanmean(wind_lines, axis=0)
            std_line = np.nanstd(wind_lines, axis=0)

            wind_mean = np.nanmean(mean_line[7:10])
    
            axs[loc_ind][era].plot(wind_mean, mean_lines[mi][-1], 
                                   marker='o', markersize=10, color=month_colors[mi])
            xs.append(wind_mean)
            ys.append(mean_lines[mi][-1])
            
        if len(xs)==0: continue
            
        #### lines of best fit? correlations?
        mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
        m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
        ax6 =  axs[loc_ind][era]
        x = np.linspace(np.nanmin(xs), np.nanmax(xs), 50)
        ax6.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
        ax6a = ax6.twinx()
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

# final plot
for era, years in enumerate(decades): axs[-1][era].set_xlabel('Wind Speed')
for loc_ind, loc in enumerate(hemi_names): axs[loc_ind][0].set_ylabel('Normalized Change in MIZ Ice Area')


#%% end
