#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 2025 | Nov 10

seaice_sens2.py

- test sensitivy to different storm areas
- 1000 hpa, 990 hpa - bounding boxes, 1000 hpa contour area

@author: mundi
"""

#%% imports and files
import numpy as np
import matplotlib.pyplot as plt
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress

import xarray as xr
from datetime import timedelta
import warnings

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

hemi_names= ['Arctic', 'Antarctic']

xxx = np.arange(-7,14+1,1)
xlabels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]

fontsize=11

AREA_NAMES = ['1000 hPa Bounding Box','990 hPa Bounding Box','1000 hPa Contour']
AREA_COLORS = ['maroon', 'navy', 'green']

#%%% open new sea ice

def seaice_lines(years, census_path, area_path, si_path, AREA_IND, sample_size=10):
    ice_lims = [20,80]
    
    # addon_num = 0 # storm area
    sia_var_addons = [['_1000', '_990', '_contour'],['2_1000', '2_990', '2_contour']]
    
    miz_ind = 1 #[0=daily, 1=total]
    
    # data dict
    lines = {}
    start_day = {}
    end_day = {}
    si_changes = {} 
    clim_changes = {}
    for idx in np.arange(0,12): # months
        lines[idx+1] = []
        start_day[idx+1] = []
        end_day[idx+1] = []
        si_changes[idx+1] = []
        clim_changes[idx+1] = []
    
    # data loop
    for year in years:
        
        census_file = census_path+'census_'+str(year)+'.csv'
        [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
        
        timing_grid = []
        for xx in range(0,len(startdate)):
            timing_grid.append((startdate[xx], enddate[xx]))
        
        storm_ranges = []
        analysis_ranges = []
        for startdt, enddt in timing_grid:
            week_ago = startdt - timedelta(days=7)
            two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
            analysis_ranges.append(fx.daterange(week_ago, two_week, dt=24))
            storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
        
        # open ice area
        ds_area = xr.open_dataset(area_path + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
        try:
            ds = xr.open_dataset(si_path + str(year) + '_seaice' + '.nc')
        except:
            print('- skip: '+si_path + str(year) + '_seaice' + '.nc')
            continue
        
        for storm_num, strm in enumerate(storm_ranges):
             month = int(strm[0].month)
             
             ### remove storms that don't interact with the ice
             with warnings.catch_warnings():
                 warnings.simplefilter('ignore')
                 ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
             if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                continue
             try:
                 sia = ds['sia_miz'+sia_var_addons[miz_ind][AREA_IND]].values[storm_num]
             except (KeyError, IndexError):
                 print()
                 print('***', year, len(storm_ranges), 
                       np.shape(ds['sia_miz'+sia_var_addons[miz_ind][AREA_IND]].values))
                 print()
                 break
        
             # get time series   
             try:
                 sia_clim = ds['sia_clim_miz'+sia_var_addons[miz_ind][AREA_IND]].values[storm_num]
             except KeyError:
                 sia_clim = ds['si_clim'].values[storm_num]
             ss = sia-sia_clim
             with warnings.catch_warnings():
                 warnings.simplefilter('ignore')
                 standardized_area = (ss-ss[0])/(np.nanmax(ss)-np.nanmin(ss))
             pcd = (standardized_area)
        
             # plot
             lines[month].append(pcd)
             start_day[month].append(strm[0])
             end_day[month].append(strm[-1])
             
             si_changes[month].append(sia)
             clim_changes[month].append(sia_clim)
        
    mean_lines = []  
    for idx in np.arange(0,12):
        
        if len(lines[idx+1])> sample_size: 
            mean_line = np.nanmean(lines[idx+1], axis=0)
            mean_lines.append(mean_line)
        else:
            mean_lines.append(np.nan*np.ones(np.shape(analysis_ranges[0])))
            continue
        
    return mean_lines, lines, start_day, end_day, si_changes, clim_changes


#%% spaghetti plot - all months

for AREA_IND in [1,2]:
    for loc_ind, loc in enumerate(hemi_names):
        path1= root_paths[loc_ind]
        
        for era, years in enumerate(decades):
            if era==0:continue
            ystr = ' ('+str(years[0])+'-'+str(years[-1])+')'
        
            #### set up plot
            fig, axes_all = plt.subplots(3, 4, figsize=(10,6.5), sharex=True, sharey=True)
            axes_flat = axes_all.flatten()
            
            fig.suptitle('Change in MIZ Area For All Storms: '+loc+ystr+'\n'+AREA_NAMES[AREA_IND], fontsize=fontsize+1)
            
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
                ax1.text(0.0225, 0.915, '('+next(alph)+')', 
                         transform=ax1.transAxes, fontsize=fontsize-2, 
                         zorder=50)
            for ax in axes_all[-1,:]:
                ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
                
        
            #### data plotting
            mean_line_og, lines_og = fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')[0:2]
                
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', AREA_IND)
            
            
            for mi, month in enumerate(months):
                if len(lines[month])==0: continue
                
                axes_flat[mi].plot(xxx, np.array(lines[month]).T, lw=0.55, color='gray')
                axes_flat[mi].plot(xxx, np.nanmean(lines[month], axis=0), lw=2, color=AREA_COLORS[AREA_IND])
                
                axes_flat[mi].plot(xxx, mean_line_og[mi], lw=2, color=AREA_COLORS[0])
                
                axes_flat[mi].set_title(calendar.month_abbr[month]+', n='+str(len(lines[month]))+','+str(len(lines_og[month])),
                                        fontsize=fontsize)
                
#%%%% monthly fan plot
fs=14
alph1 = ['a','b','c','d','e','f', 'g', 'h', 'i', 'j', 'k']
month_colors = ['#238443','#78c679','#c2e699','#d7b5d8','#df65b0','#dd1c77',
                '#980043','#7a0177','#253494','#2c7fb8','#41b6c4','#a1dab4']

for AREA_IND in [1,2]:

    ### set up
    fig2, axes2 = plt.subplots(2,2, figsize=(16,14), sharey=True)
    fig2.suptitle('\nMonthly Mean Change in MIZ Ice Area\n'+AREA_NAMES[AREA_IND], 
                  fontweight='bold', fontsize=fs+2)
    for i, ax2 in enumerate(axes2.flatten()):
        ax2.axhline(0, ls='-', color='k', lw=1)
        ax2.axvline(0, ls=':', color='gray', lw=1)
        ax2.set_xlim(-7,14)
        ax2.set_xticks(xxx)
        ax2.set_xticklabels(xlabels, minor=False, rotation=0,fontsize=fs)
        ax2.yaxis.set_tick_params(labelleft=True)
        ax2.tick_params(axis='both', labelsize=fs)
        ax2.text(0.0225, 1.025, '('+alph1[i]+')',transform=ax2.transAxes, 
                  fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                      'edgecolor':'white', 'lw':0.75},zorder=50)
    for ax2 in axes2[1,:]: ax2.set_xlabel('Days Since Storm Start',fontsize=fs)
    for ax2 in axes2[:,0]: ax2.set_ylabel('Normalized Relative Change in Ice Area',fontsize=fs+1)
        
    ### data and plot
    for li, loc in enumerate(['arctic','aa']):
        for yi, years in enumerate(decades):
            if yi==0:continue
            
            yr_title = str(years[0])+'-'+str(years[-1])
            
            path1 = root_paths[li]
            mean_lines = fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')[0]
                
            new_lines = seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', AREA_IND)[0]
             
            ### monthly mean
            axes2[li][yi].set_title(hemi_names[li]+': '+yr_title, fontsize=fs+2)
            for idx, ml in enumerate(mean_lines):
                # original lines
                axes2[li][yi].plot(xxx, ml, color = month_colors[idx], lw=1, ls='--')
                
                # new area lines
                axes2[li][yi].plot(xxx, new_lines[idx], color = month_colors[idx],
                                   lw=2, label = calendar.month_name[idx+1])
                
                # # month labels
                # yadd = ytext_scale(idx, li, yi)
                # mlabel = calendar.month_abbr[idx+1] #+' ('+str(len(lines[idx+1]))+')'
                # axes2[li][yi].text(xxx[-1]+0.25, ml[-1]+yadd, mlabel, fontsize=fs-1) 
            
    
    axes2[0][1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                      edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                      bbox_to_anchor=(1.425,0.33))
    
#%%%% difference?

for AREA_IND in [1,2]:

    ### set up
    fig2, axes2 = plt.subplots(2,2, figsize=(16,14), sharey=True)
    fig2.suptitle('\nMonthly Mean Difference MIZ Ice Area Change\n'+AREA_NAMES[AREA_IND], 
                  fontweight='bold', fontsize=fs+2)
    for i, ax2 in enumerate(axes2.flatten()):
        ax2.axhline(0, ls='-', color='k', lw=1)
        ax2.axvline(0, ls=':', color='gray', lw=1)
        ax2.set_xlim(-7,14)
        ax2.set_xticks(xxx)
        ax2.set_xticklabels(xlabels, minor=False, rotation=0,fontsize=fs)
        ax2.yaxis.set_tick_params(labelleft=True)
        ax2.tick_params(axis='both', labelsize=fs)
        ax2.text(0.0225, 1.025, '('+alph1[i]+')',transform=ax2.transAxes, 
                  fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                      'edgecolor':'white', 'lw':0.75},zorder=50)
    for ax2 in axes2[1,:]: ax2.set_xlabel('Days Since Storm Start',fontsize=fs)
    for ax2 in axes2[:,0]: ax2.set_ylabel('Normalized Relative Change in Ice Area',fontsize=fs+1)
        
    ### data and plot
    for li, loc in enumerate(['arctic','aa']):
        for yi, years in enumerate(decades):
            if yi==0:continue
            
            yr_title = str(years[0])+'-'+str(years[-1])
            
            path1 = root_paths[li]
            mean_lines = fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')[0]
                
            new_lines = seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', AREA_IND)[0]
             
            ### monthly mean
            axes2[li][yi].set_title(hemi_names[li]+': '+yr_title, fontsize=fs+2)
            for idx, ml in enumerate(mean_lines):
                
                # new area lines
                axes2[li][yi].plot(xxx, new_lines[idx]-ml, color = month_colors[idx],
                                   lw=2, label = calendar.month_name[idx+1])
                
                # # month labels
                # yadd = ytext_scale(idx, li, yi)
                # mlabel = calendar.month_abbr[idx+1] #+' ('+str(len(lines[idx+1]))+')'
                # axes2[li][yi].text(xxx[-1]+0.25, ml[-1]+yadd, mlabel, fontsize=fs-1) 
            
    
    axes2[0][1].legend(loc='upper right', ncol=1, handletextpad=0.5, handlelength=1,
                      edgecolor=(1, 1, 1, 0), facecolor=(1, 1, 1, 0),fontsize=fs+2,
                      bbox_to_anchor=(1.425,0.33))
    
#%%%% histogram of 3,7,14 values (statistics)
from scipy.stats import ttest_ind, ks_2samp, mannwhitneyu

tseries_ind = -1 #14 #
tseries_lab = ['Week Before']+['']*6+['Day 0','','','3-day','','','','7-day']+['']*6+['14-day']
bins = np.arange(-1, 1.1, 0.1)

for AREA_IND in [1,2]:
    for loc_ind, loc in enumerate(hemi_names):
        path1= root_paths[loc_ind]
        
        for era, years in enumerate(decades):
            if era==0:continue
            ystr = ' ('+str(years[0])+'-'+str(years[-1])+')'
        
            #### set up plot
            fig, axes_all = plt.subplots(3, 4, figsize=(10,6.5), sharex=True, sharey=True)
            axes_flat = axes_all.flatten()
            
            fig.suptitle('Change in MIZ Area For All Storms: '+loc+ystr+'\n'+
                         AREA_NAMES[AREA_IND]+': '+tseries_lab[tseries_ind], fontsize=fontsize+1)
            
            fig.text(0.075, 0.33, 'Storm Counts', 
                     fontsize=fontsize, rotation=90)
            fig.text(0.4, 0.05, 'Normalized Relative Change in Ice Area', fontsize=fontsize)
        
            alph = iter(list(string.ascii_lowercase))
            for ax1 in axes_flat:
                ax1.axvline(0, ls='-', color='k', lw=0.75)
                ax1.tick_params(axis='both', which='major', labelsize=fontsize)  
                ax1.text(0.0, 1.05, '('+next(alph)+')', 
                         transform=ax1.transAxes, fontsize=fontsize-2, 
                         zorder=50)
        
            #### data plotting
            lines_og = fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')[1]
                
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', AREA_IND)
            
            
            for mi, month in enumerate(months):
                
                values_og = np.array([ln[tseries_ind] for ln in lines_og[month]])
                values = np.array([ln[tseries_ind] for ln in lines[month]])
                if len(values_og)<10: continue
                
                hist, bin_edges = np.histogram(values_og, bins=bins)
                axes_flat[mi].bar(bin_edges[:-1], hist, width=0.1, 
                               facecolor=AREA_COLORS[0],edgecolor=AREA_COLORS[0], 
                               alpha=0.33, zorder=0)
                
                hist, bin_edges = np.histogram(values, bins=bins)
                axes_flat[mi].bar(bin_edges[:-1], hist, width=0.1, 
                               facecolor=AREA_COLORS[AREA_IND],edgecolor=AREA_COLORS[AREA_IND], 
                               alpha=0.33, zorder=1)
                
                #### statistics
                mask1 = ~np.isnan(values_og) ####!!!! split masks!
                mask2 = ~np.isnan(values)
                t, pt = ttest_ind(values_og[mask1], values[mask2])
                ks, pks = ks_2samp(values_og[mask1], values[mask2])
                mw, pmw = mannwhitneyu(values_og[mask1], values[mask2])
                
                pvals = [pt, pmw, pks]
                cl = [100*(1-p) for p in pvals]
                stats_text = 't: '+ str(round(cl[0],1)) +'%\nMW: '+str(round(cl[1],1))+'%\nKS: '+str(round(cl[-1],1))+'%'
                
                axes_flat[mi].text(0.0225, 0.715, stats_text, 
                             transform=axes_flat[mi].transAxes, fontsize=fontsize-2, 
                             zorder=50)
                
                # if any are above 95% add asterisk to title...
                cadd='*' if np.any(np.array(cl)>=95) else ''
                axes_flat[mi].set_title('   '+calendar.month_abbr[month]+cadd+', n='+str(len(values_og)),
                                        fontsize=fontsize)
                
            
            # ### save
            # savepath = '/Users/mundi/Desktop/out_figs/todel/'
            # savename = str(AREA_IND)+'_'+str(loc_ind)+'_'+str(era)+'_seaice'
            # fig.savefig(savepath+savename+'.png', bbox_inches = "tight")
                
                
#%%% non-normalized

for AREA_IND in [1,2]:
    for loc_ind, loc in enumerate(hemi_names):
        path1= root_paths[loc_ind]
        
        for era, years in enumerate(decades):
            if era==0:continue
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
                
            mean_lines1, lines1, start_day1, end_day1, si_changes1, clim_changes1 = \
                seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', AREA_IND)
             
                
            for mi, month in enumerate(months):
                if len(lines1[month])==0:continue
                
                rel_lines = np.array(si_changes[month]) - np.array(clim_changes[month])
                change_lines = [(line-line[0])/1e5 for line in rel_lines]
                
                axes_flat[mi].plot(xxx, np.nanmean(change_lines, axis=0), lw=2, color=AREA_COLORS[0], zorder=99)
                
                rel_lines1 = np.array(si_changes1[month]) - np.array(clim_changes1[month])
                change_lines1 = [(line-line[0])/1e5 for line in rel_lines1]
                
                axes_flat[mi].plot(xxx, np.array(change_lines1).T, lw=0.55, color='gray')
                axes_flat[mi].plot(xxx, np.nanmean(change_lines1, axis=0), lw=2, color=AREA_COLORS[AREA_IND], zorder=100)
                
                
                
                axes_flat[mi].set_title(calendar.month_abbr[month]+', n='+str(len(lines[month])),
                                        fontsize=fontsize)
            
            # ### save
            # savepath = '/Users/mundi/Desktop/out_figs/todel/'
            # savename = str(AREA_IND)+'_'+str(loc_ind)+'_'+str(era)+'_nonnorm'
            # fig.savefig(savepath+savename+'.png', bbox_inches = "tight")
raise ValueError('end')
#%% 500 km radius
from cartopy.geodesic import Geodesic
gd = Geodesic()

ice_lims = [20,80]

def accumulated_miz(startday, enddt, si_lon, ice_fname):
    t1 = startday - timedelta(days=1)
    t2 = enddt + timedelta(days=1)
    miz_range = fx.daterange(t1, t2, dt=24)

    miz_points = np.zeros(np.shape(si_lon))
    for miz_dt in miz_range:
        sicm = fx.load_seaice(ice_fname, miz_dt.year, miz_dt.month, miz_dt.day, latlon=False)
        miz_points = np.where(((sicm>=0.15) & (sicm<=0.80)), 1, miz_points)

    return miz_points

def get_total_area(masked_grid, cell_area=25*25):
    ''' does not multiply by cell element! '''
    area_array = np.where(~np.isnan(masked_grid), cell_area, np.nan)
    return np.nansum(np.nansum(area_array))

def calc_seaice_area(my_sic_grid, cell_area=25*25):
    my_area = np.where(~np.isnan(my_sic_grid), my_sic_grid*cell_area, np.nan)
    return np.nansum(np.nansum(my_area))

#%%% get data
radius = 800

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    out_path = path1+'sensitivity/seaice/'
    print(loc)
    
    # get sea ice grid
    if loc_ind==0: 
        _, si_lon, si_lat = fx.load_seaice(ice_fname, 2010,6,1)
        seaice_path = ice_fname
    elif loc_ind==1: 
        _, si_lon, si_lat = fx.load_seaice_sh(ice_fname+'south/', 2010,6,1)
        seaice_path = ice_fname+'south/'

    # data loop
    for era, years in enumerate(decades):
        if era==0:continue
        si_list = []
        clim_list = []
        norm_list = []
        flag_list = []
        
        for year in years:
           file_name = out_path + str(loc_ind)+'_'+str(year)+'_'+str(radius)+'km_radius.npy'
           file2 = out_path + str(loc_ind)+'_'+str(year)+'_'+str(radius)+'km_flag.npy'
            
           try: # open saved data
               si_array = np.load(file_name)
               flags = np.load(file2)
           except FileNotFoundError:
               print('- ' + str(year)) 
                
               census_path = path1+'sensitivity/census/'
               ncpath_area = path1+'sensitivity/area/'
               ncadd = '_areas'
               
               ncpath = path1+'sensitivity/seaice/'
               ncadd='_500kmradius' 
               
               #### get location of storm minimum pressure
               # import pickle 
               # file = open(path1+'sensitivity/all_storms/'+str(year) + '_all_storms.pkl', 'rb')
               # storm_data = pickle.load(file)
               
               # open census 
               [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
                   fx.readCensus(census_path+'census_'+str(year)+'.csv', convertDT=True)

               storm_ranges = []
               analysis_ranges = []
               for startdt, enddt in zip(startdate, enddate):
                   week_ago = startdt - timedelta(days=7)
                   two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
                   analysis_ranges.append(fx.daterange(week_ago, two_week, dt=24))
                   storm_ranges.append(fx.daterange(startdt, enddt, dt=24))   
               
                         
               for storm_num, strm in enumerate(storm_ranges):
                   
                    daymonth_grid = [(dt.month, dt.day) for dt in analysis_ranges[storm_num]]
                    
                    stormstr = storm_ranges[storm_num][0].strftime('%Y_%m%d')
                    month = int(strm[0].month)
                    flag = 0
                    
                    # event1 = np.array(storm_data[storm_num])
                    # data_sort = event1[event1[:, 2].argsort()]
                    # [x1,y1,p1,t1] = data_sort[0] # get location of min pressure
                    # loc_lon, loc_lat = x1, y1
                    
                    loc_lon = startlon[storm_num]
                    loc_lat = startlat[storm_num]

                    
                    # calculate new storm boundaries
                    cp = gd.circle(lon=loc_lon, lat=loc_lat, radius=radius*1000.)
                    in_circle = fx.find_points_in_contour(cp, si_lon,si_lat)
                    miz_points = accumulated_miz(strm[0], strm[-1], si_lon, seaice_path)
                    
                    #### get ice area fraction
                    
                    ## total bbox area
                    if loc_ind==0: si_grid, si_lon, si_lat = fx.load_seaice(seaice_path, year, strm[0].month, strm[0].day, latlon=True)
                    elif loc_ind==1: si_grid, si_lon, si_lat = fx.load_seaice_sh(seaice_path, year, strm[0].month, strm[0].day, latlon=True)
                    si_in = np.ma.masked_array(si_grid, mask=in_circle).filled(np.nan)
                    box_area = get_total_area(si_in)
                    
                    ## get ice areas (>15%, >80%)
                    ice_area15 = calc_seaice_area(np.where(si_in<0.15, 0, si_in)) 
                    ice_area80 = calc_seaice_area(np.where(si_in<0.80, 0, si_in))
                    ## calculate ice area fraction
                    ice_frac = ice_area80*100/box_area
                   
                    if np.isnan(ice_frac) or np.isinf(ice_frac):
                       ice_frac=0
                    if (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
                       print('>> ice area: '+strm[0].strftime('%Y-%m-%d')+' '+str(round(ice_frac))+'%')
                       flag = 1
                        
                    ### calculate area
                    sia, sia_clim = [],[]
                    for month, day in daymonth_grid: 
                        # get daily sea ice
                        if loc_ind==0: sic = fx.load_seaice(seaice_path, year, month, day, latlon=False)
                        elif loc_ind==1: sic = fx.load_seaice_sh(seaice_path, year, month, day, latlon=False)
                       
                        # all marginal points
                        miz_masked = np.ma.masked_where(miz_points==0, sic).filled(np.nan)
                        miz_masked_radi = np.nansum((np.ma.masked_array(miz_masked, mask=in_circle).filled(np.nan))*25*25)
                        sia.append(miz_masked_radi)
        
                        all_yr2 = []
                        for yr in years:
                            if month==2 and day==29: day=28
                            
                            if loc_ind ==0: sict = fx.load_seaice(seaice_path, yr, month, day, latlon=False)
                            if loc_ind ==1: sict = fx.load_seaice_sh(seaice_path, yr, month, day, latlon=False)
                            
                            miz_mask = np.ma.masked_where(miz_points==0, sict).filled(np.nan)
                            all_yr2.append( np.nansum(np.ma.masked_array(miz_mask, mask=in_circle).filled(np.nan)*25*25) )
                        sia_clim.append(np.nanmean(all_yr2))
        
                    if (len(np.unique(sia))==1 and np.unique(sia)[0] == 0) or np.isnan(np.mean(sia[0:10+1])):
                        print('>> no ice: '+strm[0].strftime('%Y-%m-%d')+' '+str(np.unique(sia)))
                        flag = 2
                    
                    ### relativize
                    ss = np.array(sia) - np.array(sia_clim)
                    standardized_area = (ss-ss[0])/(np.nanmax(ss)-np.nanmin(ss))
                    pcd = (standardized_area)
    
                    ### append to storm list
                    si_list.append(sia)
                    clim_list.append(sia_clim)
                    norm_list.append(pcd)
                    flag_list.append(flag)
                    
               #### save to array
               si_array = np.array([si_list, clim_list, norm_list])
               np.save(file_name, si_array)
               flags = np.array(flag_list)
               np.save(file2, flags)
           
           
           
           
#%%% comparison plots

def seaice_radius(loc_ind, years, radius=500, sample_size=10):
    
    path1= root_paths[loc_ind]
    out_path = path1+'sensitivity/seaice/'
    
    # flags = {mm:[] for mm in np.arange(1, 12+1)}
    lines = {mm:[] for mm in np.arange(1, 12+1)}
    si_changes = {mm:[] for mm in np.arange(1, 12+1)}
    clim_changes = {mm:[] for mm in np.arange(1, 12+1)}
    
    
    for year in years:
        census_file = path1+'sensitivity/census/'+'census_'+str(year)+'.csv'
        [startdate, enddate] = fx.readCensus(census_file, convertDT=True)[0]
        
        file2 = out_path + str(loc_ind)+'_'+str(year)+'_'+str(radius)+'km_flag.npy'
        flag_list = np.load(file2) 
        
        file_name = out_path + str(loc_ind)+'_'+str(year)+'_'+str(radius)+'km_radius.npy'
        si_array = np.load(file_name)
        # si_array = np.array([si_list, clim_list, norm_list])
        si_list = si_array[0]
        clim_list = si_array[1]
        norm_list = si_array[2]
        
        # open ice area
        ds_area = xr.open_dataset(path1+'sensitivity/area/' + str(year) +'_area.nc')
        ice_sorter = ds_area['ice_area80'].values
        box_area = ds_area['box_area'].values
        
       
        for storm_num, start, flag, si, clim, norm in zip(np.arange(0, len(startdate)), startdate, flag_list, si_list, clim_list, norm_list):
            if flag ==2: continue
        
            ### remove storms that don't interact with the ice
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                ice_frac = ice_sorter[storm_num]*100/box_area[storm_num]             
            if np.isnan(ice_frac) or np.isinf(ice_frac) or (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
               continue
        
            ss = np.array(si) - np.array(clim)
            norm = (ss-ss[0])/(np.nanmax(ss)-np.nanmin(ss))

            lines[start.month].append(norm)
            si_changes[start.month].append(si)
            clim_changes[start.month].append(clim)
    
        # for mm in np.arange(1, 12+1):
        #     print(year, mm, len(lines9[mm]), len(lines[mm]))
    
    mean_lines = []  
    for idx in np.arange(0,12):
        if len(lines[idx+1])> sample_size: 
            mean_line = np.nanmean(lines[idx+1], axis=0)
            mean_lines.append(mean_line)
        else:
            mean_lines.append([np.nan]*22)


    return mean_lines, lines, si_changes, clim_changes        



#%%%% spaghetti
run_rad = 500

for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    out_path = path1+'sensitivity/seaice/'
    print(loc)
    
    # data loop
    for era, years in enumerate(decades):
        if era==0:continue
        ystr = ' ('+str(years[0])+'-'+str(years[-1])+')'
    
        #### set up plot
        fig, axes_all = plt.subplots(3, 4, figsize=(10,6.5), sharex=True, sharey=True)
        axes_flat = axes_all.flatten()
        
        fig.suptitle('Change in MIZ Area For All Storms: '+loc+ystr+'\n'+str(run_rad)+' km radius', fontsize=fontsize+1)
        
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
            ax1.text(0.0225, 0.915, '('+next(alph)+')', 
                     transform=ax1.transAxes, fontsize=fontsize-2, 
                     zorder=50)
        for ax in axes_all[-1,:]:
            ax.set_xticklabels(xlabels, minor=False, rotation=0, fontsize=fontsize)
            
    
        #### data plotting
        mean_line_og, lines_og = fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')[0:2]
            
        mean_lines, lines, si_changes, clim_changes = seaice_radius(loc_ind, years, radius=run_rad)  
        
        for mi, month in enumerate(months):
            if len(lines[month])<1: continue
            
            axes_flat[mi].plot(xxx, np.array(lines[month]).T, lw=0.55, color='gray')
            axes_flat[mi].plot(xxx, np.nanmean(lines[month], axis=0), lw=2, color='goldenrod')
            
            axes_flat[mi].plot(xxx, mean_line_og[mi], lw=2, color='maroon')
            
            axes_flat[mi].set_title(calendar.month_abbr[month]+', n='+str(len(lines[month]))+', '+str(len(lines_og[month])),
                                    fontsize=fontsize)
            
        
    
#%%%%% non-normalized
for loc_ind, loc in enumerate(hemi_names):
    path1= root_paths[loc_ind]
    
    for era, years in enumerate(decades):
        if era==0:continue
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
            if loc_ind==0: ax1.set_ylim(-0.75, 0.75)
            elif loc_ind==1: ax1.set_ylim(-2.25,2.25)
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
            
        mean_lines1, lines1, si_changes1, clim_changes1 = seaice_radius(loc_ind, years, radius=run_rad)  
            
        for mi, month in enumerate(months):
            if len(lines1[month])<3: continue
            
            rel_lines = np.array(si_changes[month]) - np.array(clim_changes[month])
            change_lines = [(line-line[0])/1e5 for line in rel_lines]
            axes_flat[mi].plot(xxx, np.nanmean(change_lines, axis=0), lw=2, color='maroon', zorder=99)
            
            rel_lines1 = np.array(si_changes1[month]) - np.array(clim_changes1[month])
            change_lines1 = [(line-line[0])/1e5 for line in rel_lines1]
            axes_flat[mi].plot(xxx, np.array(change_lines1).T, lw=0.55, color='gray')
            axes_flat[mi].plot(xxx, np.nanmean(change_lines1, axis=0), lw=2, color='goldenrod', zorder=100)
            
            
            
            axes_flat[mi].set_title(calendar.month_abbr[month]+', n='+str(len(lines[month])),
                                    fontsize=fontsize)
            
        ### save
        # savepath = '/Users/mundi/Desktop/out_figs/todel/'
        # savename = str(AREA_IND)+'_'+str(loc_ind)+'_'+str(era)+'_seaice'
        # fig.savefig(savepath+savename+'.png', bbox_inches = "tight")


#%% summed annual change! bar plots, separate hemis

# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

area_shades = [['royalblue','mediumblue','navy'],
               ['mediumaquamarine','mediumseagreen','green'],
               ['gold','goldenrod','darkgoldenrod']
               ]
fs=14
XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]
alphbar = ['a','b','c','d']

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True, sharey='row')
fig2.suptitle('\nMonthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    
fig2.subplots_adjust(wspace=0.1)

#### data and plot

values = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    # if li==1:continue
    
    axes2[li][0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[li][0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    era_values = {}
    for yi, years in enumerate(decades):
        if yi==0:continue
        era_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        mean_lines1, lines1, start_day1, end_day1, si_changes1, clim_changes1 = \
            seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', 1)
        mean_lines2, lines2, start_day2, end_day2, si_changes2, clim_changes2 = \
            seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', 2)
        # mean_lines3, lines3, si_changes3, clim_changes3 = seaice_radius(li, years, radius=800)  
             
            
            
        axes2[li][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            
            areas_rel =[]
            for si, clim in zip([si_changes1, si_changes2],#, si_changes3],
                                [clim_changes1, clim_changes2]): #, clim_changes3]):
                rel_area = []
                for sia, sia_clim in zip(si[mm], clim[mm]):
                    ss = sia-sia_clim
                    rel_area.append(ss-ss[0])
                areas_rel.append(np.nanmean(rel_area, axis=0)/1e5)
            
            if len(lines[mi+1])<10: era_values[yi].append([np.nan]*3); continue ###!!!
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[li][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.33, edgecolor=hemi_colors[li], lw=2)
            
            for ai, AREA_IND in enumerate([1,2]):#,3]):
                meanline_a = areas_rel[ai]
                diffs = [meanline_a[ii] for ii in INDS]
                # era_values[yi].append(diffs)
                # axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0], zorder=-5,
                #         facecolor=area_shades[ai], alpha=0.33, edgecolor=area_shades[ai], lw=2)
                for xa, ya, iii in zip(xvals, diffs, np.arange(0, len(xvals))):
                    axes2[li][yi].plot(xa, ya, color=area_shades[ai][iii], marker='o')
                

        ### TOTAL
        sums = [np.nansum(np.array(era_values[yi])[:,dr]) for dr in range(len(xvals))]
        axes2[li][yi].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        
        if not shift: axes2[0][yi].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        axes2[li][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values

#%%% range min/max
plot_mean = True

# bar plot shades (3,7,14-days)
hemi_colors = ['#01665e','#8c510a']
hemi_colors2 = [['#c7eae5','#5ab4ac','#01665e'],['#f6e8c3','#d8b365','#8c510a']]
shades = ['#f0f0f0','#bdbdbd','#636363'] # grayscale sum

area_shades = [['royalblue','mediumblue','navy'],
               ['mediumaquamarine','mediumseagreen','green'],
               ['gold','goldenrod','darkgoldenrod']
               ]
fs=14
XSPACING = [0,0.3,0.6]
INDS = [7+3, 7+7, 7+14]
alphbar = ['a','b','c','d']

#### set up
fig2, axes2 = plt.subplots(2,2, figsize=(15,10), sharex=True, sharey='row')
fig2.suptitle('Monthly Mean 3-, 7-, and 14-Day Changes in MIZ Ice Area\nMin/Max', fontweight='bold', fontsize=fs+2)
for i, ax2 in enumerate(axes2.flatten()):
    ax2.axhline(0, ls='-', color='k', lw=1)
    ax2.set_xlim(-0.5,14.25)
    ax2.set_xticks(list(np.arange(XSPACING[1], 12+XSPACING[1]))+[13.175])
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.tick_params(axis='both', labelsize=fs)
    ax2.text(0.0225, 1.025, '('+alphbar[i]+')',transform=ax2.transAxes, 
              fontsize=fs, bbox={'facecolor': 'white', 'alpha': 0, 'pad':5, 
                                  'edgecolor':'white', 'lw':0.75},zorder=50)
for ax2 in axes2[1,:]: ax2.set_xlabel('Month',fontsize=fs)
for ax2 in axes2[:,0]: ax2.set_ylabel('Relative Change in Area\n'+r'($\times 10^5$ km$^2$)',fontsize=fs+1)
    
fig2.subplots_adjust(wspace=0.1)

#### data and plot

values = {}
valmin = {}; valmax = {}
for li, loc in enumerate(['Arctic','Antarctic']):
    shift = False 
    # if li==1:continue
    
    axes2[li][0].plot([],[], lw=6, color=hemi_colors[li], label=loc)
    axes2[li][0].legend(loc='lower left', fontsize=fs, handletextpad=0.5, handlelength=1.25)
    
    era_values = {}
    min_values = {}; max_values = {}
    for yi, years in enumerate(decades):
        if yi==0:continue
        era_values[yi] = []
        min_values[yi] = []; max_values[yi] = []
        
        yr_title = str(years[0])+'-'+str(years[-1])
        
        path1 = root_paths[li]
        mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
            fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
            
        mean_lines1, lines1, start_day1, end_day1, si_changes1, clim_changes1 = \
            seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', 1)
        mean_lines2, lines2, start_day2, end_day2, si_changes2, clim_changes2 = \
            seaice_lines(years, path1+'sensitivity/census/', path1+'sensitivity/area/', path1+'sensitivity/seaice/', 2)
        # mean_lines3, lines3, si_changes3, clim_changes3 = seaice_radius(li, years, radius=800)  
             
            
            
        axes2[li][yi].set_title(yr_title, fontsize=fs+2)
        for mi, mm in enumerate(months):
            rel_area = []
            for sia, sia_clim in zip(si_changes[mm], clim_changes[mm]):
                ss = sia-sia_clim
                rel_area.append(ss-ss[0])
            if plot_mean: MEANLINE = np.nanmean(rel_area, axis=0)/1e5
            else: MEANLINE = np.nansum(rel_area, axis=0)/1e5
            
            areas_rel =[]
            for si, clim in zip([si_changes1, si_changes2],#, si_changes3],
                                [clim_changes1, clim_changes2]):#, clim_changes3]):
                rel_area = []
                for sia, sia_clim in zip(si[mm], clim[mm]):
                    ss = sia-sia_clim
                    rel_area.append(ss-ss[0])
                if plot_mean: areas_rel.append(np.nanmean(rel_area, axis=0)/1e5)
                else: areas_rel.append(np.nansum(rel_area, axis=0)/1e5)
            
            if len(lines[mi+1])<10: ###!!!
                era_values[yi].append([np.nan]*3)
                min_values[yi].append([np.nan]*3)
                max_values[yi].append([np.nan]*3)
                continue 
            
            if not shift: xvals = mi+np.array(XSPACING)
            else:
                xvals = mi+np.array(XSPACING) + 6
                if xvals[0]>=12: xvals -= 12
                
            diffs = [MEANLINE[ii] for ii in INDS]
            era_values[yi].append(diffs)
            axes2[li][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0],
                    facecolor=hemi_colors2[li], alpha=0.33, edgecolor=hemi_colors[li], lw=2)
            
            diff_arr= [diffs]
            for ai, AREA_IND in enumerate([1,2]):#,3]):
                meanline_a = areas_rel[ai]
                diffs = [meanline_a[ii] for ii in INDS]
                diff_arr.append(diffs)
                # era_values[yi].append(diffs)
                # axes2[0][yi].bar(xvals, diffs, width=XSPACING[1]-XSPACING[0], zorder=-5,
                #         facecolor=area_shades[ai], alpha=0.33, edgecolor=area_shades[ai], lw=2)
                for xa, ya, iii in zip(xvals, diffs, np.arange(0, len(xvals))):
                    axes2[li][yi].plot(xa, ya, color=area_shades[ai][iii], marker='o')
                    
            diff_arr=np.array(diff_arr)
            min_arr = np.nanmin(diff_arr, axis=0)
            max_arr = np.nanmax(diff_arr, axis=0)
            axes2[li][yi].errorbar(xvals, diffs, yerr=[np.abs(min_arr-diffs), np.abs(max_arr-diffs)], 
                                   fmt='none', capsize=3, ecolor='k')
            min_values[yi].append(min_arr)
            max_values[yi].append(max_arr)

        ### TOTAL
        sums = np.array( [np.nansum(np.array(era_values[yi])[:,dr]) for dr in range(len(xvals))] )
        axes2[li][yi].bar(xvals+1.33+li, sums, width=XSPACING[1]-XSPACING[0],
                facecolor=hemi_colors2[li], alpha=0.5, edgecolor=hemi_colors[li], lw=2)
        
        min_sums = np.array( [np.nansum(np.array(min_values[yi])[:,dr]) for dr in range(len(xvals))] )
        max_sums = np.array( [np.nansum(np.array(max_values[yi])[:,dr]) for dr in range(len(xvals))] )
        axes2[li][yi].errorbar(xvals+1.33+li, sums, yerr=[np.abs(min_sums-sums), np.abs(max_sums-sums)],
                               fmt='none', capsize=3, ecolor='k')
        
        if not shift: axes2[0][yi].set_xticklabels(list(np.arange(1,12+1))+['Total'], fontsize=fs)
        else: axes2[0][yi].set_xticklabels(list(np.arange(7,12+1))+list(np.arange(1,7)))
        axes2[li][yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
        
    values[li] = era_values
    valmin[li] = min_values
    valmax[li] = max_values
 
#%%% min/max total sum
fig, axes = plt.subplots(1,2, figsize=(15,5))
fig.suptitle('Total Sum (Min/Max)')
for yi, years in enumerate(decades):
    if yi==0:continue
    axes[yi].set_title(str(years[0])+'-'+str(years[-1]))
    axes[yi].axhline(0, ls='-', color='k', lw=1)
    
    sums = []
    for vi, VAL in enumerate([values, valmin, valmax]):
        val0 = np.where(np.isnan(np.array(VAL[0][yi])), 0, np.array(VAL[0][yi]))
        val1 = np.where(np.isnan(np.array(VAL[1][yi])), 0, np.array(VAL[1][yi]))
        sums.append( val0 + val1 )
    
        if vi==0:
            heights = []
            for mi, sum1 in enumerate(sums[0]):
                xvals = mi+np.array(XSPACING)
                heights.append(sum1)
                
                fcs = [hemi_colors2[0][i] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors2[1][i] for i in [0,1,2]]
                ecs = [hemi_colors[0] if np.abs(val0[mi][i])>np.abs(val1[mi][i]) else hemi_colors[1] for i in [0,1,2]]
                
                axes[yi].bar(xvals, sum1, width=XSPACING[1]-XSPACING[0],
                        facecolor=fcs, alpha=0.5, edgecolor=ecs, lw=2)

            axes[yi].bar(xvals+1.33+0.5, np.nansum(sums[0], axis=0), width=XSPACING[1]-XSPACING[0],
                    facecolor=shades, alpha=0.5, edgecolor='k', lw=2)
            axes[yi].axvline(1+xvals[0], ls='--', color='k', lw=0.85)
            
        #### add spread bars!
        elif vi==1:
            val0_min = val0
            val1_min = val1
        elif vi==2:
            for mi in np.arange(0,12):
                xvals = mi+np.array(XSPACING)
                
                axes[yi].errorbar(xvals, heights[mi], 
                                yerr=[np.abs(sums[1][mi]-heights[mi]), np.abs(sums[2][mi]-heights[mi])],
                                fmt='none', capsize=3, ecolor='k')
                
            # # total sum 
            # axes[yi].errorbar(xvals+1.33+0.5,  np.nansum(sums[0], axis=0), 
            #                 yerr=[np.abs(np.nansum(sums[1], axis=0)-heights[mi]), np.abs(np.nansum(sums[2], axis=0)-heights[mi])],
            #                 fmt='none', capsize=3, ecolor='k')
                
                
            
            

#%% end