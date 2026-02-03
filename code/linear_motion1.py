#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 2025
linear_motion1.py

- plotting

@author: mundi
"""

#%% imports, file names
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import time as timeIN
import calendar
import functions as fx
from datetime import datetime, timedelta
from scipy.interpolate import griddata
import warnings
from glob import glob


ROOT = '/Users/mundi/Desktop/'
ice_path = ROOT +'seaice/'
ice_paths = [ice_path, ice_path+'south/']

ice_motion_path = '/Users/mundi/Desktop/seaice/ice_motion/'

nh_path = ROOT +'month-hemi/nh_data/'
sh_path = ROOT +'month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

savepath = ROOT+'month-hemi/linear_motion/'

census_name = 'census_'
contour_name = '_contours.nc'
si_name = '_seaice'

#%%% controls

# years = np.arange(2010,2020)
years = np.arange(1982, 1992)
# years = np.arange(2010,2015)

months = np.arange(1,12+1)

calc_time = 0 # 0 = storm duration, 1 = 7 days, 2 = 3 weeks

xxx = np.arange(-7,14+1,1)
ice_lims = [20,80]
miz = [0.15,0.80]

lon_num = 100
lon_x = np.linspace(0, 1, lon_num)

DATA = 'v'


#%% data organization

def get_im_lines(loc_ind, year, month, savepath, root_path, DATA):
    savename = str(loc_ind)+'_'+str(calc_time)+'_'+str(year)+'_'+DATA+'.npy'
    im_series = np.load(savepath+savename)
    
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
        
        ice_frac = ice_sorter[sn]*100/box_area[sn]             
        if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
            continue
        
        lines.append(im_series[sn])
    
    return lines

#%% PLOT

# years = np.arange(2010,2020)
years = np.arange(1982, 1992)
ystr = str(years[0]) +'-'+ str(years[-1])

grouped_months = [[2,3,4],[5,6,7],[8,9,10],[11,12,1]]

wind_name = {'u':'Zonal', 'v':'Meridional'}[DATA]

#%%% plot lines
# set up plot
nrows = len(grouped_months)+1
fig, axes = plt.subplots(nrows, 2, figsize=(12,3.5*nrows), sharey=True,
                         gridspec_kw={"height_ratios":[0.001]+([1]*(nrows-1))})
for ax1 in axes.flatten():
    ax1.set_xlim([lon_x[0], lon_x[-1]])
    ax1.set_xticks(np.arange(lon_x[0], lon_x[-1]+0.25, 0.25))
    ax1.axhline(0, lw=0.55, color='gray', ls=':')
    ax1.axvline(0.5, lw=0.55, color='gray', ls=':')
for ax1 in axes[-1][:]: ax1.set_xlabel('Longitude Extent Fraction')
for ax1 in axes[:,0]: ax1.set_ylabel(wind_name + ' Wind')

# organize data
for loc_ind, loc in enumerate(['Arctic', 'Antarctic']):
    
    root_path = root_paths[loc_ind]
    axes[0][loc_ind].axis("off")
    axes[0][loc_ind].set_title('\n'+loc, fontweight='bold')
    
    
    for mx, months in enumerate(grouped_months):
        axes[mx+1][loc_ind].set_title(ystr+'  '+str(months))
        data = []
        for month in months:
            for year in years:
                try:
                    lines = get_im_lines(loc_ind, year, month, savepath, root_path, DATA)
                    data += lines
                except FileNotFoundError:
                    print(loc, year, month, 'file not found')
                    
        for lin in data:
            axes[mx+1][loc_ind].plot(lon_x, lin, lw=0.5, color='gray')
            
        if len(data)<10: continue
            
        # plot mean line
        mean_line = np.nanmean(data, axis=0)
        axes[mx+1][loc_ind].plot(lon_x, mean_line, lw=2, color='k', label='n='+str(len(data)))
        axes[mx+1][loc_ind].legend(loc='upper right')
        
        # add shading
        pos_change = np.ma.masked_where(mean_line<0, mean_line).filled(np.nan)
        neg_change = np.ma.masked_where(mean_line>0, mean_line).filled(np.nan)
        axes[mx+1][loc_ind].fill_between(lon_x, neg_change, y2=0, color='salmon', alpha=0.5)
        axes[mx+1][loc_ind].fill_between(lon_x, pos_change, y2=0, color='lightskyblue', alpha=0.5)



#%% end
