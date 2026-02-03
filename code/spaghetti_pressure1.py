#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 2025
spaghetti_pressure1.py



@author: mundi
"""
#%% imports and files
import numpy as np
import matplotlib.pyplot as plt
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress

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


#%% spaghetti plot

loc_ind=1
path1= root_paths[loc_ind]

month_groups = [[4,5], [11,12]]

#### set up plot
fig, axes_all = plt.subplots(len(month_groups), len(decades), 
                           figsize=(8,6.5), sharex=True, sharey=True)
fig.suptitle('Change in MIZ Area For All Storms', fontsize=fontsize+1)
for ax in axes_all[:,0]:
    ax.set_ylabel('Normalized Relative\nChange in Ice Area', fontsize=fontsize)
for ax in axes_all[1,:]:
    ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
    ax.set_xlabel('Days Since Storm Start', fontsize=fontsize-1)
alph = iter(list(string.ascii_lowercase))
for axl in axes_all:
    for ax1 in axl:
        ax1.set_xlim(-7,14)
        ax1.set_ylim(-1.05,1.05)
        ax1.axhline(0, ls='-', color='k', lw=1)
        ax1.axvline(0, ls='-', color='k', lw=0.75)
        ax1.set_xticks(xxx)
        ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
        ax1.text(0.0225, 0.915, next(alph), transform=ax1.transAxes, fontsize=fontsize, 
                bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                      'edgecolor':'k', 'lw':0.75},zorder=50)

#### data plotting
for era, years in enumerate(decades):
    ystr = str(years[0])+'-'+str(years[-1])
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
    for mi, month_group in enumerate(month_groups):
        
        all_lines = []
        for month in month_group:
            all_lines+=lines[month]

        axes_all[mi][era].plot(xxx, np.array(all_lines).T, lw=0.55, color='gray')
        
        axes_all[mi][era].plot(xxx, np.nanmean(all_lines, axis=0), lw=2, color='maroon')
        
        mnames = [calendar.month_abbr[mm]+', ' for mm in month_group]
        axes_all[mi][era].set_title(ystr+': '+mnames[0]+mnames[1][:-2])

#%% all months spaghetti

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        ystr = ' ('+str(years[0])+'-'+str(years[-1])+')'
    
        #### set up plot
        fig, axes_all = plt.subplots(3, 4, figsize=(10,6.5), sharex=True, sharey=True)
        axes_flat = axes_all.flatten()
        
        fig.suptitle('Change in MIZ Area For All Storms: '+loc+ystr, fontsize=fontsize+1)
        
        fig.text(0.075, 0.33, 'Normalized Relative Change in Ice Area', 
                 fontsize=fontsize, rotation=90)
        fig.text(0.4, 0.05, 'Days Since Storm Start', fontsize=fontsize)
    
        alph = iter(list(string.ascii_lowercase))
        for ax1 in axes_flat:
            ax1.set_xlim(-7,14)
            ax1.set_ylim(-1.05,1.05)
            ax1.axhline(0, ls='-', color='k', lw=1)
            ax1.axvline(0, ls='-', color='k', lw=0.75)
            ax1.set_xticks(xxx)
            ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
            # ax1.text(0.0225, 0.915, next(alph), transform=ax1.transAxes, fontsize=fontsize-2, 
            #         bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
            #               'edgecolor':'k', 'lw':0.75},zorder=50)
            ax1.text(0.0225, 0.915, '('+next(alph)+')', 
                     transform=ax1.transAxes, fontsize=fontsize-2, 
                     zorder=50)
        for ax in axes_all[-1,:]:
            ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
            
    
        #### data plotting
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for mi, month in enumerate(months):
            
            axes_flat[mi].plot(xxx, np.array(lines[month]).T, lw=0.55, color='gray')
            
            axes_flat[mi].plot(xxx, np.nanmean(lines[month], axis=0), lw=2, color='maroon')
            
            axes_flat[mi].set_title(calendar.month_abbr[month]+', n='+str(len(lines[month])),
                                    fontsize=fontsize)
            
#%%% non-normalized

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        ystr = ' ('+str(years[0])+'-'+str(years[-1])+')'
    
        #### set up plot
        fig, axes_all = plt.subplots(3, 4, figsize=(10,6.5), sharex=True, sharey=True)
        axes_flat = axes_all.flatten()
        
        fig.suptitle('Change in MIZ Area For All Storms: '+loc+ystr, fontsize=fontsize+1)
        
        fig.text(0.075, 0.33, 'Relative Change in Ice Area '+r'$\times10^5$ km$^2$', 
                 fontsize=fontsize, rotation=90)
        fig.text(0.4, 0.05, 'Days Since Storm Start', fontsize=fontsize)
    
        alph = iter(list(string.ascii_lowercase))
        for ax1 in axes_flat:
            if loc_ind==0: ax1.set_ylim(-2, 2)
            elif loc_ind==1: ax1.set_ylim(-3,3)
            ax1.set_xlim(-7,14)
            ax1.axhline(0, ls='-', color='k', lw=1)
            ax1.axvline(0, ls='-', color='k', lw=0.75)
            ax1.set_xticks(xxx)
            ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
            ax1.text(0.0225, 0.915, '('+next(alph)+')', 
                     transform=ax1.transAxes, fontsize=fontsize-2, 
                     zorder=50)
        for ax in axes_all[-1,:]:
            ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
            
    
        #### data plotting
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        for mi, month in enumerate(months):
            
            rel_lines = np.array(si_changes[month]) - np.array(clim_changes[month])
            change_lines = [(line-line[0])/1e5 for line in rel_lines]
            
            axes_flat[mi].plot(xxx, np.array(change_lines).T, lw=0.55, color='gray')
            
            axes_flat[mi].plot(xxx, np.nanmean(change_lines, axis=0), lw=2, color='maroon')
            
            axes_flat[mi].set_title(calendar.month_abbr[month]+', n='+str(len(lines[month])),
                                    fontsize=fontsize)
#%%---------------------

#%% storm pressure analysis
from scipy.stats import ttest_ind, ks_2samp, mannwhitneyu

cl_data = {}
for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    decade_colors = ['#2166ac', '#b2182b']
    
    if loc_ind==0: bins = np.arange(960,984, 2)
    elif loc_ind==1: bins = np.arange(930, 957, 2)
    
    #### set up plot
    fig, axes_all = plt.subplots(3, 4, figsize=(10,6.5), sharex=True, sharey=True)
    axes_flat = axes_all.flatten()
    
    fig.suptitle('Mean Minimum Storm Pressure: '+loc, fontsize=fontsize+1)
    
    fig.text(0.075, 0.33, 'Storm Counts', fontsize=fontsize, rotation=90)
    fig.text(0.4, 0.025, 'Minimum Pressure (hPa)', fontsize=fontsize)
    
    alph = iter(list(string.ascii_lowercase))
    for ax1 in axes_flat:
        ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
        ax1.text(0.0, 1.05, '('+next(alph)+')', 
                 transform=ax1.transAxes, fontsize=fontsize-2, 
                 zorder=50)
    
    for mi, month in enumerate(months):
        
        era_data =[]
        for era, years in enumerate(decades):
            ystr = str(years[0])+'-'+str(years[-1])
        
            #### data
            start, end, pressure = fx.storm_pressure(years, path1+'census/', path1+'area/')
        
            p_list = pressure[month]
            if len(p_list)<10: continue
            era_data.append(p_list)
            
            hist, bin_edges = np.histogram(p_list, bins=bins, density=False)
            axes_flat[mi].bar(bin_edges[:-1], hist, width=1.75, alpha=0.5,
                              facecolor=decade_colors[era], edgecolor=decade_colors[era])
            
            axes_flat[mi].axvline(np.nanmean(p_list), color = decade_colors[era], label=ystr)
            
            
        #### statistics 
        if len(era_data)==2:
            t, pt = ttest_ind(era_data[0], era_data[1])
            ks, pks = ks_2samp(era_data[0], era_data[1])
            mw, pmw = mannwhitneyu(era_data[0], era_data[1])
            
            pvals = [pt, pmw, pks]
            cl = [100*(1-p) for p in pvals]
            cl_data[str(loc_ind)+'_'+str(month)] = cl
            stats_text = 't: '+ str(round(cl[0],1)) +'%\nMW: '+str(round(cl[1],1))+'%\nKS: '+str(round(cl[-1],1))+'%'
            
            axes_flat[mi].text(0.0225, 0.715, stats_text, 
                         transform=axes_flat[mi].transAxes, fontsize=fontsize-2, 
                         zorder=50)
            
            # if any are above 95% add asterisk to title...
            cadd='*' if np.any(np.array(cl)>=95) else ''
            
            axes_flat[mi].set_title('   '+calendar.month_abbr[month]+cadd+', n='+str(len(era_data[0]))+', n='+str(len(era_data[1])),
                                    fontsize=fontsize)
            
    axes_flat[mi].legend(loc='lower left', handletextpad=0.5, handlelength=1.25,
                         ncol = 2, bbox_to_anchor=(-0.45, -0.45))
            
#%%% plot mean pressures, noting significance

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
# fig.suptitle('Mean Minimum Storm Pressure')

alph = iter(list(string.ascii_lowercase))
for ax1 in axes:
    ax1.text(0.0, 1.05, '('+next(alph)+')', 
             transform=ax1.transAxes, fontsize=fontsize-2, 
             zorder=50)

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    decade_colors = ['#2166ac', '#b2182b']
    
    ax = axes[loc_ind]
    ax.set_title(loc)
    ax.set_ylabel('Mean Pressure (hPa)')
    ax.set_xticks(np.arange(1,12+1))
    
    for era, years in enumerate(decades):
        ystr = str(years[0])+'-'+str(years[-1])
    
        start, end, pressure = fx.storm_pressure(years, path1+'census/', path1+'area/')
    
        pressure_values = [np.nanmean(pressure[mm]) if len(pressure[mm])>10 else np.nan for mm in months]
        pressure_std = [np.nanstd(pressure[mm]) if len(pressure[mm])>10 else np.nan for mm in months]
        
        ax.plot(months, pressure_values, 
                marker='o', color = decade_colors[era], label = ystr)
        
        ax.errorbar(months, pressure_values, yerr=pressure_std, 
                    color = decade_colors[era], capsize=5)
        
        for mm in months:
            if len(pressure[mm])<10: continue
            if np.any(np.array(cl_data[str(loc_ind)+'_'+str(mm)])>=95):
                ax.axvspan(mm-0.25,mm+0.25, color='gray', alpha=0.25, zorder=-1)
         
# legend on last plot
ax.plot([],[], lw=5, color='gray', label='Significant Difference')
ax.legend(loc='upper left', handletextpad=0.5, handlelength=1.25, ncol = 1)
ax.set_xlabel('Month');

#%% change in storm counts

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
fig.suptitle('Total Storm Counts')

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    decade_colors = ['#2166ac', '#b2182b']
    
    ax = axes[loc_ind]
    ax.set_title(loc)
    ax.set_ylabel('Number of Storms')
    ax.set_xticks(np.arange(1,12+1))
    
    for era, years in enumerate(decades):
        ystr = str(years[0])+'-'+str(years[-1])
    
        start, end, pressure = fx.storm_pressure(years, path1+'census/', path1+'area/')
    
        storm_counts = [len(start[mm]) if len(start[mm])>10 else np.nan for mm in months]
        
        
        ax.plot(months, storm_counts, 
                marker='o', color = decade_colors[era], label = ystr)
        
         
# legend on last plot
ax.plot([],[], lw=5, color='gray', label='Significant Difference')
ax.legend(loc='upper left', handletextpad=0.5, handlelength=1.25, ncol = 1)
ax.set_xlabel('Month');

#%%% per year analysis

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
fig.suptitle('Storm Counts per Year')

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    decade_colors = ['#2166ac', '#b2182b']
    
    ax = axes[loc_ind]
    ax.set_title(loc)
    ax.set_ylabel('Number of Storms')
    ax.set_xticks(np.arange(1,12+1))
    
    for era, years in enumerate(decades):
        ystr = str(years[0])+'-'+str(years[-1])
        
        storm_counts = []
        for year in years:
            start, end, pressure = fx.storm_pressure([year], path1+'census/', path1+'area/')
            storm_counts.append([len(start[mm]) for mm in months])

        mean_counts = np.nanmean(storm_counts, axis=0)
        std_counts = np.nanstd(storm_counts, axis=0)
        
        ax.plot(months, mean_counts, 
                marker='o', color = decade_colors[era], label = ystr)
        
        ax.errorbar(months, mean_counts, yerr=std_counts, 
                    color = decade_colors[era], capsize=5)
        
        # for mm in months:
        #     if len(pressure[mm])<10: continue
        #     if np.any(np.array(cl_data[str(loc_ind)+'_'+str(mm)])>=95):
        #         ax.axvspan(mm-0.25,mm+0.25, color='gray', alpha=0.25, zorder=-1)

         
# legend on last plot
ax.plot([],[], lw=5, color='gray', label='Significant Difference')
ax.legend(loc='upper left', handletextpad=0.5, handlelength=1.25, ncol = 1)
ax.set_xlabel('Month');

#%%% storm count - sensitivity

if False:
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle('Total Storm Counts: Sensitivity')
    
    for loc_ind, loc in enumerate(hemi_names):
        path1= root_paths[loc_ind]
        decade_colors = ['#2166ac', '#b2182b']
        
        ax = axes[loc_ind]
        ax.set_title(loc)
        ax.set_ylabel('Number of Storms')
        ax.set_xticks(np.arange(1,12+1))
        
        for era, years in enumerate(decades):
            if era==0:continue
            ystr = str(years[0])+'-'+str(years[-1])
        
            start, end, pressure = fx.storm_pressure(years, path1+'census/', path1+'area/')
            storm_counts = [len(start[mm]) if len(start[mm])>10 else np.nan for mm in months]
            ax.plot(months, storm_counts, 
                    marker='o', color = decade_colors[era], label = ystr)
            
            start, end, pressure = fx.storm_pressure(years, path1+'sensitivity/census/', path1+'sensitivity/area/')
            storm_counts = [len(start[mm]) if len(start[mm])>10 else np.nan for mm in months]
            ax.plot(months, storm_counts, 
                    marker='o', color = 'k', label = 'Adjusted Pressure Threshold')
            
             
    # legend on last plot
    ax.plot([],[], lw=5, color='gray', label='Significant Difference')
    ax.legend(loc='upper left', handletextpad=0.5, handlelength=1.25, ncol = 1)
    ax.set_xlabel('Month');


#%% end
