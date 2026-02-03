#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 2025
seaice_impact1.py

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

# census_path_nh = nh_path + 'census/'
# contour_path_nh = nh_path + 'contours/'
# ncarea_path_nh = nh_path + 'area/'
# si_path_nh = nh_path+ 'seaice/'

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

#%% plots!
#%% -------

#%% spaghetti fans

#%%% functions
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
# plt.subplots_adjust(wspace=0.3)
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
for li, loc in enumerate(['arctic','aa']):
    for yi, years in enumerate(decades):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        # ml_path = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/z-monthly-aa-comparison/'
        # ml_file = ml_path+loc+'_mean_lines_'+years+'.npy'
        # mean_lines = np.load(ml_file)

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

#%%%% combined decades fan

#### set up
fig2, axes2 = plt.subplots(1,2, figsize=(16,12), sharey=True)
# plt.subplots_adjust(wspace=0.3)
fig2.suptitle('Monthly Mean Change in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
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
    ax2.set_xlabel('Days Since Storm Start',fontsize=fs)
axes2[0].set_ylabel('Normalized Relative Change in Ice Area',fontsize=fs+1)
    
#### data and plot
for li, loc in enumerate(['arctic','aa']):
    path1 = root_paths[li]
    
    for yi, years in enumerate(decades):
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        
        ### monthly mean
        axes2[li].set_title(hemi_names[li], fontsize=fs+2)
        for idx, ml in enumerate(mean_lines):
            axes2[li].plot(xxx, ml, color = month_colors[idx], lw=2, 
                           ls = ['--','-'][yi])
            yadd = ytext_scale(idx, li, yi)
            mlabel = calendar.month_abbr[idx+1] 
            axes2[li].text(xxx[-1]+0.25, ml[-1]+yadd, mlabel, fontsize=fs-1) 
       
for mi, mm in enumerate(months):
    axes2[1].plot([],[], color = month_colors[mi], lw=2, label = calendar.month_name[mm])
for yi, years in enumerate(decades):
    yr_title = str(years[0])+'-'+str(years[-1])
    axes2[1].plot([],[], color = 'gray', ls = ['--','-'][yi], lw=2, 
                  label = yr_title)

axes2[1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                  edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                  bbox_to_anchor=(1.5,0.9))


#%%%% compare with extent rate of change
if False:
    '''
    storm and extent scales dont match well
    linregress slope vs
    '''
    from scipy.stats import linregress
    name1 = 'ice_area'
    name2 = 'extent_' # 'miz_'
    
    #### set up original plot
    fig2, axes2 = plt.subplots(2,2, figsize=(16,14), sharey=True)
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
    ## data plotting ##
    for li, loc in enumerate(['arctic','aa']):
        for yi, years in enumerate(decades):
            yr_title = str(years[0])+'-'+str(years[-1])
            path1 = root_paths[li]
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            # monthly mean
            axes2[li][yi].set_title(hemi_names[li]+': '+yr_title, fontsize=fs+2)
            # for idx, ml in enumerate(mean_lines):
                # axes2[li][yi].plot(xxx, ml, color = month_colors[idx], lw=2, label = calendar.month_name[idx+1])
            for mm in months:
                si = np.array(si_changes[mm]) - np.array(clim_changes[mm])
                change = [ice-ice[0] for ice in si]
                if len(change)<10: continue
                axes2[li][yi].plot(xxx, np.gradient(np.nanmean(change,axis=0),1), color = month_colors[mm-1], lw=2, label = calendar.month_name[mm])
    axes2[0][1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                      edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                      bbox_to_anchor=(1.5,0.33))
    
    #### organize extent data
    for era, years in enumerate(decades):
        yrstr = str(years[0])+'-'+str(years[-1])
    
        area_hemi = {}
        annual_means = []
        for loc_ind, loc in enumerate(hemi_names):
            area_hemi[loc_ind] = []
            for yi, year in enumerate(years):
                areas = np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy')
                if len(areas)>365:
                    areas2 = []
                    for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                        if dt!=datetime(year,2,29): 
                            areas2.append(areas[di])
                    areas=areas2
                    
                area_hemi[loc_ind].append(np.array(areas)/1e6)
            annual_means.append(np.nanmean(area_hemi[loc_ind], axis=0))
            
            #### compute monthly rate of change
            midx2=0
            for mm in months:
                midx1 = midx2
                midx2+= calendar.monthrange(2010, mm)[-1]
                portion = annual_means[loc_ind][midx1:midx2]
                slope, b, r, p, se = linregress(np.arange(0, len(portion)), portion)
                
                axes2[loc_ind][era].plot(xxx, (xxx*slope), color = month_colors[mm-1],
                                         lw=1.5, ls='--')
        

#%%% bar plots: day changes/interannual variability?

XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(18,12), sharey=True)
fig2.suptitle('\n\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,12)
    ax2.set_xticks(np.arange(XSPACING[1], 12+XSPACING[1]))
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Normalized Relative Change in Ice Area',fontsize=fs+1)
    
#### data and plot

for li, loc in enumerate(['arctic','aa']):
    shift = False if li==0 else True
    
    for yi, years in enumerate(decades):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[li][yi].set_title(hemi_names[li]+': '+yr_title, fontsize=fs+2)
        for mi, mean_line in enumerate(mean_lines):
            if len(lines[mi+1])<5: continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [mean_line[ii] for ii in INDS]
            axes2[li][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=month_colors[mi], alpha=0.5, edgecolor=month_colors[mi])
            
            ##### INTERANNUAL
            mlines = lines[mi+1]
            ## annual spread
            mstd = np.nanstd(mlines, axis=0)
            stdevs = [mstd[ii] for ii in INDS]
            axes2[li][yi].errorbar(xvals, diffs, yerr=stdevs, color=month_colors[mi], fmt='none')
            ## annual min/max
            # mmax = np.max(mlines, axis=0) 
            # maxs = [mmax[ii] for ii in INDS]
            # axes2[li][yi].scatter(xvals, maxs, marker='o', edgecolors=month_colors[mi], facecolors='white')
            # mmin = np.min(mlines, axis=0) 
            # mins = [mmin[ii] for ii in INDS]
            # axes2[li][yi].scatter(xvals, mins, marker='o', edgecolors=month_colors[mi], facecolors='white')
            
        if not shift: axes2[li][yi].set_xticklabels(np.arange(1,12+1), fontsize=fs)
        else: axes2[li][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
         
#%%% month/year heat map

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(18,12))
fig2.suptitle('\nMean Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    ax2.set_yticks(months)
    ax2.tick_params(axis='both', labelsize=fs)
    
for ax2 in axes2[1,:]: ax2.set_xlabel('Year',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Month',fontsize=fs)

#### data and plot

for li, loc in enumerate(['arctic','aa']):
    shift = False if li==0 else True
    
    axes = axes2[li]
    
    for yi, years in enumerate(decades):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        axes2[li][yi].set_xticks(years)
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[li][yi].set_title(hemi_names[li]+': '+yr_title, fontsize=fs+2)
        
        monthly_mean = []
        for mm in months:
            mm_lines = lines[mm]
            mm_start = start_day[mm]
            
            yy_lines = {y:[] for y in years}
            for lin, sdt in zip(mm_lines, mm_start):
                yy_lines[sdt.year].append(lin)
                
            monthly_mean.append([np.nanmean(yy_lines[yy], axis=0) if len(yy_lines[yy])>0 else np.nan*np.ones((22,)) for yy in years])
            
        # collapse yy_lines into an array
        monthly_mean = np.array(monthly_mean)
        # select timing
        marr = monthly_mean[:,:,7+7]
        
        cmap = cmo.balance_r
        cmap.set_bad('gray',1.)
        cb = axes2[li][yi].pcolormesh(years, months, marr, vmin=-1, vmax=1, cmap=cmap)
        
cax1 = fig2.add_axes([0.33,0.0,0.4,0.04]) 
cbar1 = fig2.colorbar(cb, cax=cax1, orientation='horizontal')
cbar1.set_label(r'Normalized Change in Sea Ice Area', fontsize=fs)
cax1.tick_params(labelsize=fs)

#%% -------
#%% HISTOGRAMS: sic 
# FINAL/cyclones_allmonths/monthly_miz_extent.py

#%%% sic data​ (both hemispheres)
savepath = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/monthly_miz/'
_, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
_, si_lon_sh, si_lat_sh = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)

area = False

if area: fname = '_area'
else: fname = '_miz'

hemi_miz = {}
hemi_map = {}

for loc_ind, loc in enumerate(hemi_names):
    if loc_ind==0:continue
    decade_miz = []
    decade_map = []
    print(); print('Loading '+loc+' Data', end=' ')
    for years in decades:
        print('.', end='')
    
        year_miz = []
        year_map = []
    
        for year in years:
            
            try:
                if loc_ind == 0:
                    daily_miz = np.load(savepath+str(year)+fname+'.npy') 
                    daily_map = np.load(savepath+str(year)+'_map.npy') 
                elif loc_ind == 1:
                    daily_miz = np.load(savepath+str(year)+fname+'_sh.npy') 
                    daily_map = np.load(savepath+str(year)+'_map_sh.npy')   
                    
            except FileNotFoundError:
            
                print(year, end='')
                year_timer = timeIN.time()
        
                daily_miz = []
                daily_map = []
        
                for month in months:
                    print('.', end='')
                    
                    for day in np.arange(1,calendar.monthrange(year, month)[-1]+1):
                        if month==2 and day==29: continue
                        
                        if loc_ind==0: si = fx.load_seaice(ice_fname, year, month, day, latlon=False)
                        if loc_ind==1: si = fx.load_seaice_sh(ice_fname+'south/', year, month, day, latlon=False)
                        
                        if np.nanmean(si)!=0:
                            if loc_ind==0: si = np.where(si_lat<60, np.nan, si)
                            elif loc_ind==1: si = np.where(si_lat_sh > -60, np.nan, si)
                            sie = np.where(si<0.15, np.nan, si)
                            miz = np.where(sie>0.80, np.nan, sie) 
                            
                            if area: miz_area= np.where(np.isnan(miz), 0,1 )
                            else: miz_area= miz
            
                            daily_miz.append( np.nansum(miz_area*25*25) )
                            daily_map.append( miz )
                        else:
                            daily_miz.append( np.nan )
                            if loc_ind==0: daily_map.append( np.nan*np.ones(np.shape(si_lon)) )
                            elif loc_ind==1: daily_map.append( np.nan*np.ones(np.shape(si_lon_sh)) )
        
                print(' '+str(round((timeIN.time()-year_timer)/60,1))+' min')
    
                if loc_ind==0:
                    np.save(savepath+str(year)+fname+'.npy', daily_miz)   
                    np.save(savepath+str(year)+'_map.npy', daily_map)  
                elif loc_ind==1:
                    np.save(savepath+str(year)+fname+'_sh.npy', daily_miz)   
                    np.save(savepath+str(year)+'_map_sh.npy', daily_map)  
                
            year_miz.append(daily_miz)
            year_map.append(daily_map)
    
        decade_miz.append(year_miz)
        decade_map.append(year_map)
    
    hemi_miz[loc_ind] = np.array(decade_miz)
    hemi_map[loc_ind] = np.array(decade_map)
    
    
    print(' Done'); print()
    
#%%%% plot 

for loc_ind, data in zip([0,1], [hemi_map[0], hemi_map[1]]):

    fig_h, axes_h = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
    axes_h = axes_h.flatten()
    for ai, axh in enumerate(axes_h): axh.set_title(calendar.month_name[ai+1])
    bins = np.arange(15, 81, 5)
    
    for era, years in enumerate(decades):
        if era==0:
            xx = 54
        elif era==1: 
            xx = 20
        ystr = str(years[0])+'-'+str(years[-1])
        axes_h[-1].plot([],[], color = shade_colors[era], label=ystr, lw=3, alpha=0.66)
        
        monthly_data = data[era]
        mind=0
        for mi, month in enumerate(months):
            num_days = calendar.monthrange(2010, month)[1]
            mslice = monthly_data[:,mind:mind+num_days, :,:]*100
            
            mslice_miz = np.where(mslice<15, np.nan, mslice)
            mslice_miz = np.where(mslice_miz>80, np.nan, mslice_miz)
            mslice_miz = mslice_miz.flatten()
            
            zorderH = 10 if era==0 else -5
            
            hist, bin_edges = np.histogram(mslice_miz, bins=bins)
            axes_h[mi].bar(bin_edges[:-1], (hist/len(mslice_miz[~np.isnan(mslice_miz)]))*100, width=5, 
                           facecolor=shade_colors[era],edgecolor=decade_colors[era], alpha=0.66,
                           zorder=zorderH)
            
            axes_h[mi].text(xx, 15, str(round(np.nanmean(mslice_miz),1))+'%', color=decade_colors[era])
            
            mind+=num_days
    
    fig_h.text(0.05, 0.5, 'Frequency (%)',  va='center', ha='center', rotation=90, weight='bold')
    fig_h.text(0.5, 0.05, 'SIC Concentration (%)',  va='center', ha='center', weight='bold')
    
    axes_h[-1].legend(ncol=2, bbox_to_anchor=(0.95, -0.175), handletextpad=0.5, handlelength=1)
    
    if loc_ind==0: fig_h.suptitle('Mean SIC within the Arctic MIZ has increased in the recent decade')
    elif loc_ind==1: fig_h.suptitle('Mean SIC has not substantially changed within the Antarctic MIZ')


    
#%%% sic data​ (change over storms)

def get_storm_areas(year, root_path, lon, lat):
    ice_lims = [20,80]
    
    census_file = root_path+'census/'+census_name+str(year)+'.csv'
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
        ds = xr.open_dataset(root_path+'seaice/' + str(year) + si_name + '.nc')
        ds.close()
    except:
        print('- skip: '+root_path+'seaice/' + str(year) + si_name + '.nc')
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
        ncname = stormstr + contour_name +'.nc'
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
_, si_lon, si_lat = fx.load_seaice(ice_fname, 1, 1, 1, latlon=True)
_, si_lon_sh, si_lat_sh = fx.load_seaice_sh(ice_fname+'south/', 1, 1, 1, latlon=True)

area = True

if area: fname = '_area'
else: fname = '_miz'

miz_loc = []
for loc_ind, loc in enumerate(hemi_names):
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
    
#%%%% plot: change in sic for each month
    
for loc_ind, data in zip([0,1], miz_loc):

    fig_h, axes_h = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
    axes_h = axes_h.flatten()
    for ai, axh in enumerate(axes_h): 
        axh.set_title(calendar.month_name[ai+1])
        axh.axvline(0, lw=0.55, ls=':', color='gray')
    bins = np.arange(-1,1, .05)
    
    for era, years in enumerate(decades):
        if era==0:
            xx = 54
        elif era==1: 
            xx = 20
        ystr = str(years[0])+'-'+str(years[-1])
        axes_h[-1].plot([],[], color = shade_colors[era], label=ystr, lw=3, alpha=0.66)
        
        monthly_data = data[era]
        for mi, month in enumerate(months):
            
            miz_month = monthly_data[month]
            sic_diff = [storm_data[1] - storm_data[0] for storm_data in miz_month] ###!!! storm duration vs week
            
            zorderH = 10 if era==0 else -5
            
            hist, bin_edges = np.histogram(np.array(sic_diff).flatten(), bins=bins, density=True)
            axes_h[mi].bar(bin_edges[:-1], hist, width=.05, 
                           facecolor=shade_colors[era],edgecolor=decade_colors[era], alpha=0.66,
                           zorder=zorderH)
            
    fig_h.text(0.05, 0.5, 'Frequency (%)',  va='center', ha='center', rotation=90, weight='bold')
    fig_h.text(0.5, 0.05, 'Change in SIC Concentration (%)',  va='center', ha='center', weight='bold')
    
    axes_h[-1].legend(ncol=2, bbox_to_anchor=(0.95, -0.175), handletextpad=0.5, handlelength=1)
    
    if loc_ind==0: fig_h.suptitle('SIC Change Within the Arctic MIZ')
    elif loc_ind==1: fig_h.suptitle('SIC Change Within the Antarctic MIZ')
    
#%%% monthly sic histograms

for era, years in enumerate(decades):

    fig, axes = plt.subplots(2,1, figsize=(10,5), sharex=True, sharey=True)
    fig.suptitle(str(years[0])+'-'+str(years[-1]), fontweight='bold')
    
    sample_nh = {mm:np.random.rand(100) for mm in months}
    sample_sh = {mm:np.random.rand(100) for mm in months}
    
    bins = np.arange(0, 1, 0.05)
    
    antarctic_shift = True
    
    for ax, data, shift in zip(axes, [hemi_map[0],hemi_map[1]], [False, antarctic_shift]):
        monthly_data = data[era] 
        mind=0
        for mm in months:
            num_days = calendar.monthrange(2010, month)[1]
            mslice = monthly_data[:,mind:mind+num_days, :,:]
            
            mslice_miz = np.where(mslice<.15, np.nan, mslice)
            mslice_miz = np.where(mslice_miz>.80, np.nan, mslice_miz)
            mslice_miz = mslice_miz.flatten()
            
            mind+=num_days
            
            ax.axvline(mm, color='gray', alpha=0.75, lw=1, zorder=-100)
            hist, bin_edges = np.histogram(mslice_miz, bins=bins, density=True) ###!!!
            
            if not shift: xvals = bin_edges[:-1]+(mm-1)
            else:
                xvals = bin_edges[:-1]+(mm-1) + 6
                if xvals[0]>=12: xvals -= 12
            yvals = hist#/len(mslice_miz)
            ax.bar(xvals, yvals, width=bins[1]-bins[0], 
                   facecolor=month_colors[mm-1],edgecolor=month_colors[mm-1], alpha=0.66)
            
    XLIM = [months[0]-1,months[-1]]
    axes[-1].set_xlim(XLIM)
    axes[-1].set_xticks(np.arange(months[0]-1,months[-1]+1), (['0','1|0']+['1|0','1|0']*5)+['1'])
    axes[-1].set_xlabel('Sea Ice Concentration')
    
    for ax, hname, hcolor in zip(axes, hemi_names, hemi_colors): 
        ax.set_ylabel(hname, color=hcolor, fontweight='bold')
    
    top_ax = axes[0].twiny()
    top_ax.set_xlim(XLIM)
    top_ax.set_xticks(np.arange(months[0]-0.5, months[-1], 1), month_abbrs,
                      color = hemi_colors[0] if antarctic_shift else 'k');
    if antarctic_shift:
        top_ax1 = axes[1].twiny()
        top_ax1.set_xlim(XLIM)
        top_ax1.set_xticks(np.arange(months[0]-0.5, months[-1], 1), 
                           [calendar.month_abbr[mm if mm<=12 else mm-12] for mm in months+6],
                          color = hemi_colors[1]);
        

#%% -------
#%% % of storm area (ice coverage)

from scipy.stats import ttest_ind

def get_storm_fractions(year, root_path, lon, lat):
    ice_lims = [20,80]
    
    census_file = root_path+'census/'+census_name+str(year)+'.csv'
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
    
    storm_areas = []
    for storm_num, strm in enumerate(storm_ranges):
        # remove storms that don't interact with the ice
        ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
        if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
            storm_areas.append(np.nan)
            continue
        storm_areas.append(ice_frac) 
       
        
    return storm_areas, storm_ranges

#%%% data 

frac_data = []
for loc_ind, loc in enumerate(hemi_names):
    frac_decade = []
    for years in decades:
        print('*', str(years[0])+'-'+str(years[-1]))
        
        area_month = {mm:[] for mm in np.arange(1,12+1)}
        for year in years:
            ### load storm census
            path1 = root_paths[loc_ind]
            if loc_ind==0: lon1, lat1 = si_lon, si_lat 
            elif loc_ind==1: lon1, lat1 = si_lon_sh, si_lat_sh
            storm_areas, storm_ranges = get_storm_fractions(year, root_paths[loc_ind], lon1, lat1)
            
            for area, storm in zip(storm_areas, storm_ranges):
                area_month[storm[0].month].append(area)
                
        frac_decade.append(area_month)
    frac_data.append(frac_decade)
                
#%%% plot histogram

for loc_ind, data in zip([0,1], frac_data):

    fig_h, axes_h = plt.subplots(3,4, figsize = (10,8), sharex=True, sharey=True)
    axes_h = axes_h.flatten()
    bins = np.arange(20, 81, 5)
    
    for mi, month in enumerate(months):
        
        month_comparison = []
        for era, years in enumerate(decades):
            if era==0: xx = 54
            elif era==1: xx = 20
            zorderH = 10 if era==0 else -5
            ystr = str(years[0])+'-'+str(years[-1])
            
            if mi==0: axes_h[-1].plot([],[], color = shade_colors[era], label=ystr, lw=3, alpha=0.66)
            
            monthly_data = data[era]
            month_comparison.append(monthly_data[month])
            
            hist, bin_edges = np.histogram(monthly_data[month], bins=bins)
            axes_h[mi].bar(bin_edges[:-1], hist*100/len(monthly_data[month]), width=5, 
                           facecolor=shade_colors[era],edgecolor=decade_colors[era], 
                           alpha=0.66, zorder=zorderH)
            
            axes_h[mi].text(xx, 15, str(round(np.nanmean(monthly_data[month]),1))+'%', color=decade_colors[era])
            

        mask0 = ~np.isnan(np.array(month_comparison[0]))
        mask1 = ~np.isnan(np.array(month_comparison[1]))
        pval = ttest_ind(np.array(month_comparison[0])[mask0], np.array(month_comparison[1])[mask1])[1]
                               
        txt = '*' if pval<0.05 else ''                                                          
        axes_h[mi].set_title(calendar.month_name[month]+txt)
    
    fig_h.text(0.05, 0.5, 'Frequency (%)',  va='center', ha='center', rotation=90, weight='bold')
    fig_h.text(0.5, 0.05, 'Percent Pack Ice Area (%)',  va='center', ha='center', weight='bold')
    
    axes_h[-1].legend(ncol=2, bbox_to_anchor=(0.95, -0.175), handletextpad=0.5, handlelength=1)
    
    if loc_ind==0: fig_h.suptitle('Arctic Ice Area Interaction Fraction')
    elif loc_ind==1: fig_h.suptitle('Antarctic Ice Area Interaction Fraction')

#%%% simple scatter

fig_s, axes_s = plt.subplots(2,2, figsize = (10,10), sharex=True, sharey=True)
fig_s.suptitle('\nFraction of Pack Ice')
for ax in axes_s.flatten():
    ax.axhline(0, lw=0.55, ls=':', color='gray')
for ax in axes_s[-1][:]: ax.set_xlabel('Fraction of Pack Ice (%)')
for ax in axes_s[:,0]: ax.set_ylabel('Normalized Change in MIZ Ice Area')

for loc_ind, data in zip([0,1], frac_data):
    path1 = root_paths[li]
    
    for era, years in enumerate(decades):
        if era==0: xx = 54
        elif era==1: xx = 20
        zorderH = 10 if era==0 else -5
        ystr = str(years[0])+'-'+str(years[-1])
        
        axes_s[loc_ind][era].set_title(hemi_names[loc_ind]+': '+ystr)
        
        x6, y6 = [], []
        for mi, month in enumerate(months):
        
            monthly_data = data[era][month]
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
            
            x_plot = np.nanmean(monthly_data)
            y_plot = mean_lines[mi][-1]
            marker='o'
            if loc_ind==1 and y_plot<0: 
                y_plot *= -1
                marker = 'v'
            
            axes_s[loc_ind][era].plot(x_plot, y_plot,
                                      marker=marker, lw=0, markersize=10,
                                      color= month_colors[mi])
            x6.append(x_plot)
            y6.append(y_plot)
            
        mask = ~np.isnan(np.array(x6)) & ~np.isnan(np.array(y6))
        m, b, r, p, se = linregress(np.array(x6)[mask], np.array(y6)[mask])
        x = np.linspace(np.nanmin(x6), np.nanmax(x6), 100)
        axes_s[loc_ind][era].plot(x, (m*x)+b, color='gray', ls='--', lw=4, zorder=-20, alpha=0.66)

        ax6a = axes_s[loc_ind][era].twinx()
        ax6a.sharey(axes_s[loc_ind][era])  
        ax6a.plot([],[], color='gray', ls='--', lw=4, alpha=0.66, label = r'R$^2$ = '+str(round(r**2, 2)))
        ax6a.legend(loc='upper left', handletextpad=0.5, handlelength=1.5)
        ax6a.axis('off');

#%% -------
#%% GLOBAL IMPACTS
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import time
from functools import wraps
import matplotlib.pyplot as plt
import calendar

def timethis(func):
    """ 
    Print the execution time for a function call
    """
    @wraps(func)
    def wrapped_method(*args, **kwargs):
        time_start = time.time()
        output = func(*args, **kwargs)
        time_end = time.time()
        if time_end-time_start < 120:
            print(f"{func.__name__}: {(time_end-time_start)} s")
        else:
            print(f"{func.__name__}: {(time_end-time_start)/60} min")

        return output

    return wrapped_method

#%%% HEMISPHERE TRENDS

#%%%% info

# years = np.arange(2010,2020)
# years = np.arange(1982, 1992)

# years = np.arange(1990,2000)
years = np.arange(2000,2010)

months = np.arange(1,12+1)
dates = [(yy, mm) for yy in years for mm in months]

root = '/Users/mundi/Desktop/month-hemi/'

#%%%% calculate daily miz area
print('MIZ AREA')

def calc_miz_area(date):
    import numpy as np
    import functions as fx
    
    year = date[0]
    month = date[1]
    day = date[2]
    
    if loc_ind==0: si = fx.load_seaice(ice_fname, year, month, day, latlon=False)
    elif loc_ind==1: si = fx.load_seaice_sh(ice_fname, year, month, day, latlon=False)
    
    if np.nanmean(si)!=0:
        if loc_ind==0: si = np.where(si_lat<55, np.nan, si)
        elif loc_ind==1: si = np.where(si_lat>-55, np.nan, si)
        sie = np.where(si<0.15, np.nan, si)
        miz = np.where(sie>0.80, np.nan, sie) 
        
        miz_area = np.where(np.isnan(miz), 0, miz )
        name = 'ice_area'
        
        daily_miz = np.nansum(miz_area*25*25) 
        
    else:
        daily_miz = np.nan 
    return daily_miz

@timethis 
def run_threaded_areas(dates):
    with ThreadPoolExecutor(max_workers=15) as executor:
        return executor.map(calc_miz_area, dates)
    
    
for loc_ind, loc in enumerate(hemi_names):
    print(loc)
    if loc_ind==0: 
        ice_fname = '/Users/mundi/Desktop/seaice/'
        with xr.open_dataset(ice_fname + 'seaice_lonlat_v03.nc', decode_times=False) as ds:
            si_lat = ds['latitude'].values
    elif loc_ind==1: 
        ice_fname = '/Users/mundi/Desktop/seaice/'+'south/'
        _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2000,1,1, latlon=True)
        
    for year in years:
        print(year, end=' - ')
        
        try: 
            name = 'ice_area'
            miz_areas = np.load(root_paths[loc_ind]+'seaice/'+'miz_'+name+'_'+str(year)+'.npy')
            print('loaded')
        except:
        
            dates = [(year, mm, dd) for mm in months for dd in np.arange(1, calendar.monthrange(year, mm)[1]+1)]
            
            results = run_threaded_areas(dates)
            
            miz_areas = [x for x in results]
            
            name = 'ice_area'
            np.save(root_paths[loc_ind]+'seaice/'+'miz_'+name+'_'+str(year)+'.npy', miz_areas)
    
#%%%% calculate daily sea ice extent
print('ICE EXTENT')

def calc_extent(date):
    import numpy as np
    import functions as fx
    
    year = date[0]
    month = date[1]
    day = date[2]
    
    if loc_ind==0: si = fx.load_seaice(ice_fname, year, month, day, latlon=False)
    elif loc_ind==1: si = fx.load_seaice_sh(ice_fname, year, month, day, latlon=False)
    
    if np.nanmean(si)!=0:
        if loc_ind==0: si = np.where(si_lat<55, np.nan, si)
        elif loc_ind==1: si = np.where(si_lat>-55, np.nan, si)
        sie = np.where(si<0.15, np.nan, si)
        
        sie_area = np.where(np.isnan(sie), 0, sie )
        daily_si = np.nansum(sie_area*25*25) 
    else:
        daily_si = np.nan
        
    return daily_si


@timethis 
def run_threaded_areas(dates):
    with ThreadPoolExecutor(max_workers=15) as executor:
        return executor.map(calc_extent, dates)
    
    
for loc_ind, loc in enumerate(hemi_names):
    print(loc)
    if loc_ind==0: 
        ice_fname = '/Users/mundi/Desktop/seaice/'
        with xr.open_dataset(ice_fname + 'seaice_lonlat_v03.nc', decode_times=False) as ds:
            si_lat = ds['latitude'].values
    elif loc_ind==1: 
        ice_fname = '/Users/mundi/Desktop/seaice/'+'south/'
        _, si_lon, si_lat = fx.load_seaice_sh(ice_fname, 2000,1,1, latlon=True)
        
    for year in years:
        print(year, end=' - ')
        name = 'ice_area'
        try: 
            miz_areas = np.load(root_paths[loc_ind]+'seaice/'+'extent_'+name+'_'+str(year)+'.npy')
            if loc_ind==1 and year>=2014: raise FileNotFoundError
            print('loaded')
        except:
            dates = [(year, mm, dd) for mm in months for dd in np.arange(1, calendar.monthrange(year, mm)[1]+1)]
            
            results = run_threaded_areas(dates)
            
            miz_areas = [x for x in results] 
            
            np.save(root_paths[loc_ind]+'seaice/'+'extent_'+name+'_'+str(year)+'.npy', miz_areas)

#%%%% regroup data + plot full timeseries
name1 = 'ice_area'
# name2 = 'miz_' 
name2 = 'extent_'

#### set up plot
fig, [ax1,ax2] = plt.subplots(2,1, figsize=(10,6), sharex=True, height_ratios=(0.75,1))
fig.suptitle('Global ice area: '+str('MIZ' if name2=='miz_' else 'Ice Extent'))
for ax, title in zip([ax1,ax2], ['Daily Ice Area', 'Global Sum']):
    # ax.legend(loc='upper left', bbox_to_anchor=(1,1),
    #           handletextpad=0.5, handlelength=1)
    ax.set_ylabel(r'Area ($\times 10^6$ km$^2$)')
    ax.set_title(title)

####  plot annual
for loc_ind, loc in enumerate(hemi_names):
    x0 = 0
    for yi, year in enumerate(years):
        areas = np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy')
        ax1.plot(np.arange(x0,x0+len(areas)), areas/1e6, color=hemi_colors[loc_ind])
        print(x0,x0+len(areas))
        x0+=len(areas)
        
ax1.set_xlim(0,3652)
for loc_ind, loc in enumerate(hemi_names):ax1.plot([],[], color=hemi_colors[loc_ind], label=loc)
ax1.legend()
            

####  plot sum
area_hemi = {}
for loc_ind, loc in enumerate(hemi_names):
    area_hemi[loc_ind] = []
    for yi, year in enumerate(years):
        area_hemi[loc_ind]+=list(np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy'))
   
ax2.plot((np.array(area_hemi[0])+np.array(area_hemi[1]))/1e6, color='maroon')
ax2.axhline(0, lw=0.75, color='gray', ls=':')
    
#%%%% annual mean plot

#### set up plot
fig2, axes2 = plt.subplots(2,2, figsize=(15,14), sharex=True)
fig2.suptitle('\n\nAnnual Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(0,365)
    ax2.set_xticks([datetime(2010,mm,1).timetuple().tm_yday-1 for mm in months])
    ax2.set_xticklabels([calendar.month_abbr[mm] for mm in months], minor=False, rotation=0,fontsize=fs)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Days',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel(r'Area Relative to Mean\n($\times 10^6$ km$^2$)',fontsize=fs+1)

#### organize data
for era, years in enumerate([np.arange(1982,1992), np.arange(2010,2020)]):
    yrstr = str(years[0])+'-'+str(years[-1])

    area_hemi = {}
    annual_means = []
    for loc_ind, loc in enumerate(hemi_names):
        area_hemi[loc_ind] = []
        for yi, year in enumerate(years):
            areas = np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy')
            if len(areas)>365:
                areas2 = []
                for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                    if dt!=datetime(year,2,29): 
                        areas2.append(areas[di])
                areas=areas2
                
            area_hemi[loc_ind].append(np.array(areas)/1e6)
       
        annual_means.append(np.nanmean(area_hemi[loc_ind], axis=0))
        
    axes2[0][era].set_title('Annual Ice Area: '+yrstr)
    for loc_ind, loc in enumerate(hemi_names):
        axes2[0][era] = fx.plot_spread(axes2[0][era], np.arange(0,365), 
                                  annual_means[loc_ind]-np.nanmean(annual_means[loc_ind]), np.nanstd(area_hemi[loc_ind], axis=0), 
                                  color=hemi_colors[loc_ind], ls1='-', ls2='--')
        axes2[0][era].plot([],[], color=hemi_colors[loc_ind], label=loc)
    axes2[0][era].legend()
        
    axes2[1][era].set_title('Global Sum: '+yrstr)
    annual_mean = (annual_means[0]-np.nanmean(annual_means[0])) + (annual_means[1]-np.nanmean(annual_means[1]))
    annual_std = np.nanstd(area_hemi[0]+area_hemi[1], axis=0)
    axes2[1][era] = fx.plot_spread(axes2[1][era], np.arange(0,365), 
                              annual_mean, annual_std, 
                              color='maroon', ls1='-', ls2='--')
    
#%%% CYCLONE IMPACTS
                        # dkbrwn, brwn, ltbrwn, ltblue, blue, dkblue
# hemi_colors2 = ['#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]

shades = ['#f0f0f0','#bdbdbd','#636363']

#%%%% bar plot (normalized)

XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True)
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,12)
    ax2.set_xticks(np.arange(XSPACING[1], 12+XSPACING[1]))
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
    ax2.set_ylim([-0.75,0.75])
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Normalized Relative\nChange in Ice Area',fontsize=fs+1)
    
#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    axes2[0][0].plot([],[], lw=3, color=hemi_colors[li], label=loc)
    axes2[0][0].legend(loc='lower left', fontsize=fs)
    
    era_values = {}
    for yi, years in enumerate(decades):
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[0][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mean_line in enumerate(mean_lines):
            if len(lines[mi+1])<5: era_values[yi].append([np.nan]*3); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [mean_line[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)

        if not shift: axes2[0][yi].set_xticklabels(np.arange(1,12+1), fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        # axes2[0][yi].set_ylim([-0.65,0.65])
    values[li] = era_values

#### sum
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    # axes2[1][yi].set_ylim([-0.65,0.65])
    
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums = val0 + val1
    
    for mi, sum1 in enumerate(sums):
        xvals = mi+np.array(XSPACING)
        
        fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
        ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
        
        axes2[1][yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)
        

#%%%% bar plot (non - normalized)

XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True, sharey='row')
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,12)
    ax2.set_xticks(np.arange(XSPACING[1], 12+XSPACING[1]))
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    
#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    axes2[0][0].plot([],[], lw=3, color=hemi_colors[li], label=loc)
    axes2[0][0].legend(loc='lower left', fontsize=fs)
    
    era_values = {}
    for yi, years in enumerate(decades):
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[0][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            
            if len(lines[mi+1])<10: era_values[yi].append([np.nan]*3); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)

        if not shift: axes2[0][yi].set_xticklabels(np.arange(1,12+1), fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
    values[li] = era_values

#### sum
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums = val0 + val1
    
    for mi, sum1 in enumerate(sums):
        xvals = mi+np.array(XSPACING)
        
        fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
        ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
        
        axes2[1][yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)
        
#%%%%% per storm 

XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True, sharey=True)
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,12)
    ax2.set_xticks(np.arange(XSPACING[1], 12+XSPACING[1]))
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\n'+r'Per Storm ($\times 10^5$ km$^2$)',fontsize=fs+1)
    
#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    axes2[0][0].plot([],[], lw=3, color=hemi_colors[li], label=loc)
    axes2[0][0].legend(loc='lower left', fontsize=fs)
    
    era_values = {}
    for yi, years in enumerate(decades):
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[0][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            
            if len(lines[mi+1])<10: era_values[yi].append([np.nan]*3); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = np.array([MEANLINE[ii] for ii in INDS])/len(lines[mm])
            era_values[yi].append(diffs)
            axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)

        if not shift: axes2[0][yi].set_xticklabels(np.arange(1,12+1), fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
    values[li] = era_values

#### sum
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums = val0 + val1
    
    for mi, sum1 in enumerate(sums):
        xvals = mi+np.array(XSPACING)
        
        fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
        ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
        
        axes2[1][yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)        

#%%%% --> summed annual change!

XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True, sharey=True)
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    
fig2.subplots_adjust(wspace=0.1)

#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    axes2[0][0].plot([],[], lw=3, color=hemi_colors[li], label=loc)
    axes2[0][0].legend(loc='lower left', fontsize=fs)
    
    era_values = {}
    for yi, years in enumerate(decades):
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[0][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            
            if len(lines[mi+1])<10: era_values[yi].append([np.nan]*3); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)

        ### TOTAL
        sums = [np.nansum(np.array(era_values[yi])[:,dr]) for dr in range(len(xvals))]
        axes2[0][yi].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        
        if not shift: axes2[0][yi].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        axes2[0][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    
#### sum
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums = val0 + val1
    
    for mi, sum1 in enumerate(sums):
        xvals = mi+np.array(XSPACING)
        
        fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
        ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
        
        axes2[1][yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)

    axes2[1][yi].bar(xvals+1.33+0.5, np.nansum(sums, axis=0), width=XSPACING[1]-XSPACING[0],
            facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
    axes2[1][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)

#%%% storm changes as fraction of total area

name1 = 'ice_area'
name2 = 'miz_'
mdivs = [datetime(2010,mm,1).timetuple().tm_yday-1 for mm in months]+[365]

#%%%% organize total area data
area_hemi = {}
for loc_ind, loc in enumerate(hemi_names):
    area_hemi[loc_ind] = {}
    for era, years in enumerate([np.arange(1982,1992), np.arange(2010,2020)]):
        yrstr = str(years[0])+'-'+str(years[-1])
        area_hemi[loc_ind][era]={mm:[] for mm in months}
        for yi, year in enumerate(years):
            areas = np.load(root_paths[loc_ind]+'seaice/'+name2+name1+'_'+str(year)+'.npy')
            if len(areas)>365:
                areas2 = []
                for di, dt in enumerate(fx.daterange(datetime(year,1,1),datetime(year,12,31), dt=24)):
                    if dt!=datetime(year,2,29): areas2.append(areas[di])
            
            for mi, mm in enumerate(months):
                area_hemi[loc_ind][era][mm].append(areas2[mdivs[mi]:mdivs[mi+1]])
                
#### monthly means
for loc_ind, loc in enumerate(hemi_names):
    for era, years in enumerate([np.arange(1982,1992), np.arange(2010,2020)]):
        for mi, mm in enumerate(months):
            area_hemi[loc_ind][era][mm] = np.nanmean(area_hemi[loc_ind][era][mm])
    
            
#%%%% bar plot

XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True, sharey='row')
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,12)
    ax2.set_xticks(np.arange(XSPACING[1], 12+XSPACING[1]))
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\nPer Total MIZ Area (%)',fontsize=fs+1)
    
#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    axes2[0][0].plot([],[], lw=3, color=hemi_colors[li], label=loc)
    axes2[0][0].legend(loc='lower left', fontsize=fs)
    
    era_values = {}
    for yi, years in enumerate(decades):
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[0][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            MEANLINE = np.nanmean(rel_area, axis=0)
            
            if len(lines[mi+1])<10: era_values[yi].append([np.nan]*3); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[0][yi].bar(xvals, (diffs/area_hemi[li][yi][mm])*100, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)

        if not shift: axes2[0][yi].set_xticklabels(np.arange(1,12+1), fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
    values[li] = era_values

#### sum
for yi, years in enumerate(decades):
    axes2[1][yi].set_title('Sum', fontsize=fs+2)
    
    val0 = np.where(np.isnan(np.array(values[0][yi])), 0, np.array(values[0][yi]))
    val1 = np.where(np.isnan(np.array(values[1][yi])), 0, np.array(values[1][yi]))
    sums = val0 + val1
    
    for mi, sum1 in enumerate(sums):
        xvals = mi+np.array(XSPACING)
        
        fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
        ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
        
        axes2[1][yi].bar(xvals, sum1*100/(area_hemi[0][yi][mm]+area_hemi[1][yi][mm]), 
                         width=XSPACING[1]-XSPACING[0],
                facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)          

#%% -------
def align_zeros(axes):

    ylims_current = {}   #  Current ylims
    ylims_mod     = {}   #  Modified ylims
    deltas        = {}   #  ymax - ymin for ylims_current
    ratios        = {}   #  ratio of the zero point within deltas

    for ax in axes:
        ylims_current[ax] = list(ax.get_ylim())
                        # Need to convert a tuple to a list to manipulate elements.
        deltas[ax]        = ylims_current[ax][1] - ylims_current[ax][0]
        ratios[ax]        = -ylims_current[ax][0]/deltas[ax]
    
    for ax in axes:      # Loop through all axes to ensure each ax fits in others.
        ylims_mod[ax]     = [np.nan,np.nan]   # Construct a blank list
        ylims_mod[ax][1]  = max(deltas[ax] * (1-np.array(list(ratios.values()))))
                        # Choose the max value among (delta for ax)*(1-ratios),
                        # and apply it to ymax for ax
        ylims_mod[ax][0]  = min(-deltas[ax] * np.array(list(ratios.values())))
                        # Do the same for ymin
        ax.set_ylim(tuple(ylims_mod[ax]))
        
#%% Compare Positive/Negative Storm Impacts

XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True, sharey='row')
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,12)
    ax2.set_xticks(np.arange(XSPACING[1], 12+XSPACING[1]))
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alph[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area '+r'($\times 10^5$ km$^2$)',fontsize=fs+1)

for pn, label in enumerate(['Positive', 'Negative']):
    axes2[0][0].plot([],[], lw=3, color=hemi_colors[pn], label=label)
    axes2[0][0].legend(loc='lower left', fontsize=fs)
    
#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    
    for yi, years in enumerate(decades):
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        axes2[li][yi].set_title(loc+': '+yr_title, fontsize=fs+2)
        count_ax =  axes2[li][yi].twinx()
        if yi==1:  count_ax.set_ylabel('Storm Counts',fontsize=fs+1)
        counts = {'pos':[], 'neg':[]}
        for mi, mm in enumerate(months):
            rel_area = {'pos':[], 'neg':[]}
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                s1 = ss-ss[0]
                if s1[-1]>0: rel_area['pos'].append(s1)
                elif s1[-1]<0: rel_area['neg'].append(s1)
            MEANLINE_pos = np.nansum(rel_area['pos'], axis=0) ###!!!! mean vs sum?!
            MEANLINE_neg = np.nansum(rel_area['neg'], axis=0)
            
            for pn in ['pos', 'neg']: counts[pn].append(len(rel_area[pn]))
            
            if len(lines[mm])<10: continue 
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs_pos = np.array([MEANLINE_pos[ii] for ii in INDS])
            axes2[li][yi].bar(xvals, diffs_pos/1e5, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[0], alpha=0.5, edgecolor=hemi_colors[0], lw=2)
            
            diffs_neg = np.array([MEANLINE_neg[ii] for ii in INDS])
            axes2[li][yi].bar(xvals, diffs_neg/1e5, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[1], alpha=0.5, edgecolor=hemi_colors[1], lw=2)
            
        # storm counts
        for pn, m, ii in zip(['pos', 'neg'], [1,-1], [0,1]):
            count_ax.plot(np.arange(1,12+1), np.array(counts[pn])*m, marker='o',
                    color=hemi_colors[ii], lw=1.5)
        align_zeros([count_ax, axes2[li][yi]])

        if not shift: axes2[0][yi].set_xticklabels(np.arange(1,12+1), fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))


#%% -------
#%% storm overlaps


for loc_ind, loc in enumerate(hemi_names):
    for yi, years in enumerate(decades):
        
        #### storm pairs
        path1 = root_paths[loc_ind]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
        longitudes, latitudes = fx.indiv_locations(years, path1+'census/', path1+'area/')[-2:]

        date_ranges = []
        locations = []
        lines_all = []
        for mm in months:  
            for sd, ed, lons, Ls in zip(start_day[mm], end_day[mm], longitudes[mm], lines[mm]):
                date_ranges.append([sd,ed])
                locations.append(lons)
                lines_all.append(Ls)
        
        storm_pairs = []
        line_pairs = []
        
        OD = 4 # overlap range (days)
        OL = 17 # overlap longitude
        
        prev_dt = date_ranges[0]
        prev_loc = locations[0]
        prev_line = lines_all[0]
        for dt_range, lon, Ls in zip(date_ranges[1:], locations[1:], lines_all[1:]):
            # overlap_t = [True if dt in prev_dt else False for dt in dt_range]
            overlap_t = True if (dt_range[0]>prev_dt[0]-timedelta(days=OD) and dt_range[0]<prev_dt[-1]+timedelta(days=OD)) else False
            overlap_l = True if (lon>=prev_loc-OL and lon<=prev_loc+OL) else False
            
            if np.any(overlap_t) and overlap_l: 
                storm_pairs.append([prev_dt, dt_range])
                line_pairs.append([prev_line, Ls])
                
            prev_dt = dt_range
            prev_loc = lon
            prev_line = Ls

        #### plot lines
        fig_o, axes_o = plt.subplots(4,3, figsize = (14,8), sharex=True, sharey=True)
        fig_o.text(0.05, 0.5, '\n\nNormalized Relative Change in MIZ Ice Area',  va='center', ha='center', rotation=90)
        fig_o.text(0.5, 0.05, r'Days Since $Second$ Storm Start',  va='center', ha='center')
        fig_o.suptitle('\n'+'Consecutive Storms: '+loc+', '+str(years[0])+'-'+str(years[-1]))
        axes_o = axes_o.flatten()
        for ai, axo in enumerate(axes_o): 
            axo.set_title(calendar.month_name[ai+1])
            axo.plot(xxx, mean_lines[ai], color='b', lw=1.25, label='Monthly Mean')
            axo.axhline(0, lw=0.55, ls=':', color='gray')
            axo.axvline(0, lw=0.55, ls=':', color='gray')
           
        prev_mo = 0
        mo_counter=0
        lines0, lines1 = [],[]
        for tpair, lpair in zip(storm_pairs, line_pairs):
            mo = tpair[0][0].month
            if prev_mo == mo: 
                mo_counter+=1
                lines0.append( lpair[0] )
                lines1.append( lpair[1] )
            else: 
                mo_counter=0
                if lines0 and lines1:
                    axes_o[prev_mo - 1].plot(xxx, np.nanmean(lines0, axis=0), lw=2, color='k', ls='--')
                    axes_o[prev_mo - 1].plot(xxx, np.nanmean(lines1, axis=0), lw=2, color='k', ls='-')
                lines0, lines1 = [],[]
            
            offset = (tpair[1][0] - tpair[0][0]).days
            
            axes_o[mo - 1].plot(xxx-offset, lpair[0], lw=0.75, color='C'+str(mo_counter), ls='--')
            axes_o[mo - 1].plot(xxx, lpair[1], lw=0.75, color='C'+str(mo_counter), ls='-')
            prev_mo=mo
            
            
        
        # axo.plot([],[], color='white', lw=0, label=' ')
        axo.plot([],[], color='gray', ls='--', label='First Storm')
        axo.plot([],[], color='gray', ls='-', label='Second Storm')
        axo.legend(loc='lower left', bbox_to_anchor=(-0.35,-0.5), ncol=3)
        
#%% annual cycle
def setup_plot(suptitle, titles, ylims):
    fig, axes = plt.subplots(2,1, figsize = (8,6))
    fig.suptitle(suptitle)
    
    for era, ax in enumerate(axes):
        
        ax.set_title(titles[era]+': '+str(decades[era][0])+'-'+str(decades[era][-1]))
        ax.set_ylabel('Normalized Change\nin MIZ Ice Area')
        ax.set_ylim(ylims[era])
        # ax.set_xticks(era_ticks[era])
        # ax.set_xticklabels(era_labels[era])
        # for xx in [yi*365 for yi, yr in enumerate(decades[era])]:
        #     ax.axvline(xx, color='gray', zorder=-12, lw=0.55)
        ax.axhline(0, lw=0.55, ls='--', color='gray', zorder=-5)
                
    return fig, axes
            
def moving_average(a, n=3):
    ret = np.nancumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
    

fig, axes = setup_plot('Daily Time Series', 
                       titles=['Arctic', 'Antarctic'], 
                       ylims=[[-1.1,1.1],[-1.1,1.1]])

for loc_ind, loc in enumerate(hemi_names):

    for yi, years in enumerate(decades):
        tseries = []
        
        path1 = root_paths[loc_ind]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
      
        for mm in np.arange(1,12+1):
            DAYS = np.arange(1,calendar.monthrange(1990, mm)[-1]+1)
            month_series = {dd:[] for dd in DAYS}
            
            for storm, line in zip(start_day[mm], lines[mm]):
                if storm.month==2 and storm.day==29:
                    storm=datetime(storm.year, 3,1)
                month_series[storm.day].append(line[-1])
            
            for dd in DAYS:
                tseries.append(np.nanmean(month_series[dd]))
                
        #### plot
        axes[loc_ind].plot(moving_average(tseries, n=30))
         
#%% ---
#%% wind/sst - storm impact pcolor
# fig 14, steele "loitering"
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JC011182


#%%% PLOT SETUP
sst_diff = False

# wnames = {'idx':1, 'ylims':[-6,6], 'name':'Zonal'}
wnames = {'idx':2, 'ylims0':[-10,6], 'ylims1':[-4,4], 'name':'Meridional'}
# wnames = {'idx':0, 'ylims':[-2,14], 'name':'Total'}

# sst_lims = [[-5,5],[-2,2]]
if sst_diff: sst_lims = [[-2,2],[-0.75,0.75]]
else: sst_lims = [[-5,5], [-2,2]]

wbin_num = 20
sbin_num = 20

def find_nearest(array, value):     
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin() 
    return idx

#%%% data/plot

fig8, axes8 = plt.subplots(2,2, figsize=(10,10))
for ax8 in axes8.flatten():
    ax8.axhline(0, color='k', lw=0.55, ls=':')
    ax8.axvline(0, color='k', lw=0.55, ls=':')
for ax8 in axes8[-1,:]: ax8.set_xlabel('SST')
for ax8 in axes8[:,0]: ax8.set_ylabel(wnames['name']+' Winds')

for loc_ind, location in enumerate(hemi_names):
    
    wind_bins = np.linspace(wnames['ylims'+str(loc_ind)][0], wnames['ylims'+str(loc_ind)][1], wbin_num)
    sst_bins = np.linspace(sst_lims[loc_ind][0], sst_lims[loc_ind][-1], sbin_num)

    wind_edges = np.linspace(np.min(wind_bins),np.max(wind_bins), len(wind_bins)+1)
    sst_edges = np.linspace(np.min(sst_bins),np.max(sst_bins), len(sst_bins)+1)
        
    #### DATA 
    for era, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        fill_arr = [[[] for i1 in range(len(sst_bins))] for i2 in range(len(wind_bins))]
        
        path1 = root_paths[loc_ind]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
    
        ## WIND
        if loc_ind ==0: 
            wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
        elif loc_ind ==1:
            wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
        
        for mi, mm in enumerate(list(lines.keys())):
            monthly_starts = start_day[mm]
            monthly_ends = end_day[mm]
            si_change = si_changes[mm]
            clim_change = clim_changes[mm]
            
            ## organize SSTs
            SSTs = np.load(path1+'sst/'+'tseries_'+str(era)+'-'+str(mm)+'.npy')
            
            ## organize winds
            try: 
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[month] if len(arr)==22]  
                    wind_lines = np.array(windies)[:,:,wnames['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames['idx']] for arr in wind_series[month]]  
                    wind_lines = np.array(windies)[:,]
            except IndexError: print('-->', month); continue
            
            # if len(wind_lines)<10: 
            #     continue
            
            for storm, sia, clim, wind1, sst1 in zip(monthly_starts, si_change, clim_change,
                                                     wind_lines, SSTs):
                
                pcm = sia-clim
                # pcm = (pcm-pcm[0])/(np.nanmax(pcm)-np.nanmin(pcm))
                
                wind_idx = find_nearest(wind_bins, wind1[7])
                sst_val = sst1[14]-sst1[7] if sst_diff else sst1[7]
                sst_idx = find_nearest(sst_bins, sst_val)
                fill_arr[wind_idx][sst_idx].append(pcm[14]-pcm[7])
                # fill_arr[wind_idx][sst_idx].append(mm)

        #### fill_arr: clean and plot
        new_arr = []
        for row in fill_arr:
            new_row = []
            for vals in row:
                new_row.append(np.nanmean(vals))
            new_arr.append(new_row) 
                
        ax8 = axes8[loc_ind][era]
        ax8.set_title(location+': '+yr_title)
        
        pcm8 = ax8.pcolormesh(sst_edges, wind_edges,  new_arr,
                             cmap=cmo.balance_r, 
                             vmin=-1e5, vmax=1e5,
                             )
        
cax1 = fig8.add_axes([0.33,-0.033,0.4,0.04]) 
cbar1 = fig8.colorbar(pcm8, cax=cax1, orientation='horizontal')
cbar1.set_label('MIZ Ice Area', fontsize=12)
cax1.tick_params(labelsize=12-1)

#%% repeat for winds - storm impact pcolor
# fig 13, steele "loitering"
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JC011182


#%%% PLOT SETUP

wnames2 = {'idx':1, 'ylims0':[-10,10], 'ylims1':[-8,8], 'name':'Zonal'}
wnames1 = {'idx':2, 'ylims0':[-10,10], 'ylims1':[-5,5], 'name':'Meridional'}

wbin_num = 15
sbin_num = 15

def find_nearest(array, value):     
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin() 
    return idx

fig8, axes8 = plt.subplots(2,2, figsize=(10,10))
for ax8 in axes8.flatten():
    ax8.axhline(0, color='k', lw=0.55, ls=':')
    ax8.axvline(0, color='k', lw=0.55, ls=':')
for ax8 in axes8[-1,:]: ax8.set_xlabel(wnames2['name']+' Winds')
for ax8 in axes8[:,0]: ax8.set_ylabel(wnames1['name']+' Winds')

#### loc/era loop 

for loc_ind, location in enumerate(hemi_names):
    
    wind_bins = np.linspace(wnames1['ylims'+str(loc_ind)][0], wnames1['ylims'+str(loc_ind)][1], wbin_num)
    s_bins = np.linspace(wnames2['ylims'+str(loc_ind)][0], wnames2['ylims'+str(loc_ind)][1], sbin_num)

    wind_edges = np.linspace(np.min(wind_bins),np.max(wind_bins), len(wind_bins)+1)
    s_edges = np.linspace(np.min(s_bins),np.max(s_bins), len(s_bins)+1)
        
    #### DATA 
    for era, years in enumerate(decades):
        yr_title = str(years[0])+'-'+str(years[-1])
        
        fill_arr = [[[] for i1 in range(len(s_bins))] for i2 in range(len(wind_bins))]
        
        path1 = root_paths[loc_ind]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
    
        ## WIND
        if loc_ind ==0: 
            wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
        elif loc_ind ==1:
            wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
        
        for mi, mm in enumerate(list(lines.keys())):
            monthly_starts = start_day[mm]
            monthly_ends = end_day[mm]
            si_change = si_changes[mm]
            clim_change = clim_changes[mm]
            
            ## organize winds
            if loc_ind==0: 
                windies = [arr for arr in wind_series[month] if len(arr)==22]  
                wind_lines1 = np.array(windies)[:,:,wnames1['idx']]
            elif loc_ind==1: 
                windies = [arr[wnames1['idx']] for arr in wind_series[month]]  
                wind_lines1 = np.array(windies)[:,]
                
            if loc_ind==0: 
                windies = [arr for arr in wind_series[month] if len(arr)==22]  
                wind_lines2 = np.array(windies)[:,:,wnames2['idx']]
            elif loc_ind==1: 
                windies = [arr[wnames2['idx']] for arr in wind_series[month]]  
                wind_lines2 = np.array(windies)[:,]
            
            
            for storm, sia, clim, wind1, wind2 in zip(monthly_starts, si_change, clim_change,
                                                     wind_lines1, wind_lines2):
                
                pcm = sia-clim
                # pcm = (pcm-pcm[0])/(np.nanmax(pcm)-np.nanmin(pcm))
                
                wind_idx = find_nearest(wind_bins, wind1[7])
                s_idx = find_nearest(s_bins, wind2[7])
                
                fill_arr[wind_idx][s_idx].append(pcm[14]-pcm[7])
                # fill_arr[wind_idx][sst_idx].append(mm)

        #### fill_arr: clean and plot
        new_arr = []
        for row in fill_arr:
            new_row = []
            for vals in row:
                new_row.append(np.nanmean(vals))
            new_arr.append(new_row) 
                
        ax8 = axes8[loc_ind][era]
        ax8.set_title(location+': '+yr_title)
        
        pcm8 = ax8.pcolormesh(s_edges, wind_edges,  new_arr,
                             cmap=cmo.balance_r, 
                             vmin=-7.5e4, vmax=7.5e4,
                             )
        
cax1 = fig8.add_axes([0.33,-0.033,0.4,0.04]) 
cbar1 = fig8.colorbar(pcm8, cax=cax1, orientation='horizontal')
cbar1.set_label('MIZ Ice Area', fontsize=12)
cax1.tick_params(labelsize=12-1)


#%% ----
#%% daily time series
print('Starting Daily Timeseries (10-panel plots)')

from concurrent.futures import ThreadPoolExecutor
import time
from functools import wraps

month_colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

#%%% functions
def timethis(func):
    """ 
    Print the execution time for a function call
    """
    @wraps(func)
    def wrapped_method(*args, **kwargs):
        time_start = time.time()
        output = func(*args, **kwargs)
        time_end = time.time()
        if time_end-time_start < 120:
            print(f"{func.__name__}: {(time_end-time_start)} s")
        else:
            print(f"{func.__name__}: {(time_end-time_start)/60} min")

        return output

    return wrapped_method

def get_storm_dates(root_path, year):
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
    
    storm_areas = []
    startdate, enddate = [],[]
    stormstr_prev=''
    for storm_num, strm in enumerate(storm_ranges):
        
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
        
        # storm area
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
        
        with fx.HidePrint(): bbox_edges = fx.get_bbox_edges(all_contours) 
        storm_areas.append( bbox_edges )
        
        startdate.append(strm[0])
        enddate.append(strm[-1])
      
    return startdate, enddate, storm_areas

def calc_ice_area(date):
    import calendar
    import numpy as np
    
    ice_fname = '/Users/mundi/Desktop/seaice/'
    if loc_ind ==1: ice_fname+='south/'
    
    year = date[0]
    month = date[1]
    day = date[2]
    
    try:
        if loc_ind == 0: si = fx.load_seaice(ice_fname, year, month, day, latlon=False)
        elif loc_ind ==1: si = fx.load_seaice_sh(ice_fname, year, month, day, latlon=False)
    except:
        print('$$$', year, month, day)
        raise ValueError
    
    if np.nanmean(si)!=0:
        sie = np.where(si<0.15, np.nan, si)
        miz = np.where(sie>0.80, np.nan, sie) 
        
        miz_area = np.where(np.isnan(miz), 0, miz*25*25 )
        daily_miz =  np.nansum(miz_area)
        
        sie_area = np.where(np.isnan(sie), 0, sie*25*25)
        daily_sie = np.nansum(sie_area) 
        
    else:
        daily_miz = np.nan
        daily_sie = np.nan
            
    return [daily_sie, daily_miz]

@timethis 
def run_threaded_areas(dates):
    with ThreadPoolExecutor(max_workers=15) as executor:
        return executor.map(calc_ice_area, dates)

#%%% plot
decades = [np.arange(2010,2020), np.arange(1982,1992)]
months = np.arange(1,12+1)

days = np.arange(0, 365)+1

for loc_ind, loc in enumerate(['Arctic','Antarctic']):
    print('--->', loc)
    root_path = root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        #### set up plot
        figure, axes = plt.subplots(10, 1, figsize=(10,18), sharex=True, sharey=True)
        figure.suptitle('\n\n'+loc+': '+str(years[0])+'-'+str(years[-1]))
        axes2 = [ax.twinx() for ax in axes]
        for xi, axis in enumerate(axes[1:]): axis.get_shared_y_axes().joined(axis, axes2[xi])
        xticks = [0]
        xticklabels = []
        for mm in months:
            num_days = calendar.monthrange(2010, mm)[1]
            xticks.append(xticks[-1]+num_days)
            xticklabels.append(calendar.month_abbr[mm][0])
        axes[-1].set_xticks(xticks)
        axes[-1].set_xticklabels(xticklabels+['']);
        
        #### data loop
        decade_miz, decade_sie = [],[] # intitialize for climatologies
        for yi, year in enumerate(years):
            print(year)
            axes[yi].set_title(year)
            
            ## sea ice timeseries
            dates = []
            for mm in months:
                days_mm = np.arange(1, calendar.monthrange(year, mm)[1]+1)
                if mm==2 and days_mm[-1]==29: days_mm = days_mm[:-1]
                for dd in days_mm:
                    dates.append((year, mm , dd))
                    
            results = run_threaded_areas(dates)
            areas = np.array([x for x in results])
            
            axes[yi].plot(days, areas[:,0], label='SIE', color='b')
            axes2[yi].plot(days, areas[:,1], label='MIZ', color='r')
            
            decade_sie.append(areas[:,0])
            decade_miz.append(areas[:,1])
            
            ## storm times
            startdate, enddate, storm_areas = get_storm_dates(root_path, year)
            for sd, ed in zip(startdate, enddate):
                sdn = sd.timetuple().tm_yday
                edn = ed.timetuple().tm_yday
                axes[yi].axvspan(sdn, edn, alpha=0.5, color=month_colors[sd.month-1])
                
        ## climatologies
        mean_miz = np.nanmean(decade_miz, axis=0)
        mean_sie = np.nanmean(decade_sie, axis=0)
        for ax_yi in axes: ax_yi.plot(days, mean_sie, color='navy')
        for ax_yi2 in axes2: ax_yi2.plot(days, mean_miz, color='maroon')

#%% total number of storms
print('Calculating total number of storms')

storm_count = 0
for li, loc in enumerate(['arctic','aa']):
    path1 = root_paths[li]
    for yi, years in enumerate(decades):
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        print('-',loc, years[0], np.nansum([len(lines[mm]) for mm in np.arange(1,12+1)]))
            
        for mm in np.arange(1,12+1):
            storm_count += len(lines[mm])
            
print(); print('Total storm count: '+str(storm_count))        
            
            

#%% end
