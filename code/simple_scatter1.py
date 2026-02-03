#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 2025
simple_scatter1.py

@author: mundi
"""
#%% import and fpaths​
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import cmocean.tools as cmo_tools
import calendar
from datetime import datetime, timedelta
import time as timeIN

import functions as fx
import warnings

from scipy.stats import linregress

ice_fname =  '/Users/mundi/Desktop/seaice/' #'south/'

root = '/Users/mundi/Desktop/month-hemi/'
root_paths = [root+'nh_data/', root+'sh_data/']

hemi_names = ['Arctic', 'Antarctic']

decades = [np.arange(2010,2020), np.arange(1982, 1992)]
months = np.arange(1,12+1)

#%%% plotting
xxx = np.arange(-7,14+1,1)

mcolors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

def set_up_plot(xlabel='', ylabel='', axhline=None, axvline=None):
    fig, axes = plt.subplots(2,2,figsize=(10,10), sharey=False, sharex=True)
    
    for era in [0,1]: axes[-1][era].set_xlabel(xlabel)
    for loc in [0,1]: axes[loc][0].set_ylabel(ylabel)
    
    for ax in axes.flatten():
        if axhline or type(axhline)==int: ax.axhline(axhline, lw=0.55, color='gray', ls=':')
        if axvline or type(axhline)==int: ax.axvline(axvline, lw=0.55, color='gray', ls=':')
    
    return fig, axes

def best_fit_line(ax, xs, ys):
    mask = ~np.isnan(np.array(xs)) & ~np.isnan(np.array(ys))
    m, b, r, p, se = linregress(np.array(xs)[mask], np.array(ys)[mask])
    x = np.linspace(np.min(np.array(xs)[mask]), np.max(np.array(xs)[mask]), 50)
    ax.plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)
     
    ax6a = ax.twinx()
    ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
    ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
    ax6a.axis('off');
    
    return ax

#%% SST


fig, axes = set_up_plot(xlabel = r'Mean Change in SST ($^\circ$C day$^{-1}$)',
                        ylabel = 'Normalized Change in MIZ Ice Area',
                        axhline=0, axvline=0)

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]

    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))

        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        xs, ys = [], []
        for mi, mm in enumerate(list(lines.keys())):
            ml = mean_lines[mi][-1]
            
            # organize SSTs
            SSTs = np.load(path1+'sst/'+'tseries_'+str(era)+'-'+str(mm)+'.npy')
                
            ### plot
            if len(SSTs)>10:
                mean_line = np.nanmean(SSTs, axis=0)
                std_line = np.nanstd(SSTs, axis=0)
                slope, b, r, p, se = linregress(xxx, mean_line)
                
                axes[loc_ind][era].plot(slope, ml, color=mcolors[mi],
                                        marker='o', markersize=8)
                xs.append(slope); ys.append(ml)
                
        axes[loc_ind][era] = best_fit_line(axes[loc_ind][era], xs, ys)


#%% waves

fig, axes = set_up_plot(xlabel='Change in Wave Height (m)',
                        ylabel='Normalized Change in MIZ Ice Area')

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        axes[loc_ind][era].set_title(loc+': '+str(years[0])+'-'+str(years[-1]))
        
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        try:
            swh_series = fx.era_lines(years, path1+'census/', path1+'area/', path1+'swh/', 'swh')[-1]
        except FileNotFoundError as fnfe:
            print(fnfe)
            
        xs, ys = [], []
        for month in months:
            swh_lines = [np.array(s) for s in swh_series[month] if len(s)==22]
            
            if len(swh_lines)<10: 
                xs.append(np.nan)
                ys.append(np.nan)
                continue
            
            swh_ml= np.nanmean(swh_lines, axis=0)
            swh_dev = np.nanmean(swh_ml[7:10]) - np.nanmean(swh_ml)
            
            std_lines = np.nanstd([np.nanmean(SL[7:10]) - np.nanmean(SL) for SL in swh_lines])
            swh_values = [np.nanmean(SL[7:10]) - np.nanmean(SL) for SL in swh_lines]
            p25 = np.percentile(swh_values, 25)
            p75 = np.percentile(swh_values, 75) 
            # xerr = [[np.abs(p25)],[np.abs(p75)]]
            
            ml = mean_lines[month-1][-1]
            sdl = np.nanstd([LL[-1] for LL in lines[mm]])

            if loc_ind==1: 
                axes[loc_ind][era].plot(swh_dev, np.abs(ml), color=mcolors[month-1],
                                        marker='o' if ml>0 else 'v', markersize=8)
                # axes[loc_ind][era].errorbar(swh_dev, np.abs(ml), 
                #                             xerr=std_lines, yerr=sdl, color = mcolors[month-1])
                # axes[loc_ind][era].errorbar(swh_dev, np.abs(ml), 
                #                             xerr=xerr, color = mcolors[month-1])
                # axes[loc_ind][era].boxplot([swh_values], positions=[np.abs(ml)], orientation='horizontal')
                # axes[loc_ind][era].plot(p25, np.abs(ml), color=mcolors[month-1],
                #                         marker='s' if ml>0 else 'v', markersize=8)
                # axes[loc_ind][era].plot(p75, np.abs(ml), color=mcolors[month-1],
                #                         marker='s' if ml>0 else 'v', markersize=8)
            else:
                axes[loc_ind][era].plot(swh_dev, ml, color=mcolors[month-1],
                                        marker='o', markersize=8)
                # axes[loc_ind][era].errorbar(swh_dev, ml, 
                #                             xerr=std_lines, yerr=sdl, color = mcolors[month-1]) 
                # axes[loc_ind][era].errorbar(swh_dev, np.abs(ml), 
                #                             xerr=xerr, color = mcolors[month-1])
                # axes[loc_ind][era].boxplot([swh_values], positions=[np.abs(ml)], orientation='horizontal')
                # axes[loc_ind][era].plot(p25, ml, color=mcolors[month-1],
                #                         marker='s' if ml>0 else 'v', markersize=8)
                # axes[loc_ind][era].plot(p75, ml, color=mcolors[month-1],
                #                         marker='s' if ml>0 else 'v', markersize=8)
                
            xs.append(swh_dev); ys.append(ml if loc_ind==0 else np.abs(ml))
            
            
        axes[loc_ind][era] = best_fit_line(axes[loc_ind][era], xs, ys)













#%% end
