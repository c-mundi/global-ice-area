#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 2025
variable_corr1.py

> correlate "driver" varibles
> maybe heatmap? (seaborn)


@author: mundi
"""

#%% imports and files
import numpy as np
import matplotlib.pyplot as plt
import string, calendar
import cmocean.cm as cmo
from scipy.stats import linregress

import functions as fx
import seaborn as sns

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

grad_path = '/Users/mundi/Desktop/month-hemi/miz_gradient/'
ice_path = '/Users/mundi/Desktop/seaice/'
ice_paths = [ice_path, ice_path+'south/']

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

#%%% functions
def get_wind_data(widx):        
    if loc_ind==0:
        wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
    elif loc_ind==1: 
        wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')

    data_winds = {}  
    for mm in months:
        if loc_ind==0: 
            windies = [arr for arr in wind_series[mm] if len(arr)==22]  
            wind_lines = np.array(windies)[:,:,widx]
        elif loc_ind==1: 
            windies = [arr[widx] for arr in wind_series[mm]]  
            wind_lines = np.array(windies)[:,]
        data_winds[mm] = wind_lines
    return data_winds


def setup_data(loc_ind, yi, years, var_group):
    ''' see taylor digram code '''
    
    path1 = root_paths[loc_ind]
    data = {}
    
    mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
        fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')    
    nstorms = len(start_day)
    
    # sea ice impact
    if 'sia' in var_group:
        data['sia'] = lines
        # data['sia'] = {mm:si_changes[mm]-clim_changes[mm] for mm in months}
        
    # winds
    if 'wind_0' in var_group:
        widx = 0 #ZONAL
        data['wind_0'] = get_wind_data(widx)
        
    if 'wind_1' in var_group:
        widx = 1 #MERIDIONAL
        data['wind_1'] = get_wind_data(widx)
    
    # waves
    if 'waves' in var_group:
        data['waves']=fx.era_lines(years, path1+'census/', path1+'area/', path1+'swh/', 'swh')[-1]
    
    # air temperature
    if 't2m' in var_group:
        data['t2m'] = fx.era_lines(years, path1+'census/', path1+'area/', path1+'t2m/', 't2m')[-1]
    
    # ocn prof
    if 'ocn_prof' in var_group:
        if years[0]>2000:
            ocean = fx.ocn_profiles(years, path1, all_or_miz='series')[0]
            data['ocn_prof'] = {mm:np.nanmean( np.array(ocean[mm])[:,:,0:8], axis=-1) for mm in np.arange(1,12+1)}
        else:
            data['ocn_prof'] = {mm:[[np.nan]*22]*nstorms for mm in np.arange(1,12+1)}
            
    # sic gradient
    if 'sic_grad' in var_group:
        data['sic_grad'] = fx.get_sic_grad(loc, years, grad_path, path1)
    
    # si concentration
    if 'si_conc' in var_group:
        data['si_conc'] = fx.get_si_conc(loc_ind, years, grad_path, path1, ice_paths[loc_ind])
    
    # ice motion
    if 'ice_motion' in var_group:
        im_array = fx.storm_ice_motion(years, path1+'census/', path1+'area/', path1+'icemotion/')[-1]
        data['ice_motion'] = {mm:np.array(im_array[mm])[:,:,1] for mm in np.arange(1,12+1)}
                                             ###??? meridional ice motion index?
    # sst                             
    if 'sst' in var_group:
        # SSTs = np.load(root_paths[loc_ind]+'sst/'+'tseries_'+str(yi)+'-'+str(mm)+'.npy')
        data['sst'] = fx.get_sst_lines(yi, path1+'sst/')
            
    return data


#%% plot

if False:

    var_list = ['sia', 't2m', 'waves', 'wind_1']
    
    timescale = 3
    if timescale==0: # full 3-week
        idx1, idx2 = 0, 22
    elif timescale==1: # 3 days
        idx1, idx2 = 7, 10
    elif timescale==2: # 1 week
        idx1, idx2 = 7, 14
    elif timescale==3: # 2 weeks
        idx1, idx2 = 7, 22
    else:
        raise ValueError('Wrong timescale value')
        
    LAG = 0
    
    corr_type = 0 #0=indiv corrs->mean, 1=mean lines->corr
    
    cmap = cmo.balance.with_extremes(bad='gray')
    
    # month_groups = [[[4,5,6,7,8],[9,10,11,12],[1,2,3]],
    #                 [[10,11,12],[4,5,6,7],[8,9]]
    #                 ]
    # month_groups = [[np.arange(1,12+1)],[np.arange(1,12+1)]]
    month_groups = [[[4,5,6,7,8],[9,10,11,12]],
                    [[10,11,12],[4,5,6,7]]
                    ]
    
    #### set up figure
    fig = plt.figure(figsize = (9*len(month_groups[0]),10))
    subfigs = fig.subfigures(2, 2) # decades, hemis
    
    for loc_ind, loc in enumerate(hemi_names): 
        path1 = root_paths[loc_ind]
        
        month_grouping = month_groups[loc_ind]
        
        for era, years in enumerate(decades):
            ystr = str(years[0])+'-'+str(years[-1])
            
            subfigs[loc_ind][era].suptitle(loc+', '+ystr, size='x-large')
            axes = subfigs[loc_ind][era].subplots(1, len(month_grouping)) # month groups?
            if len(month_grouping)==1:axes=[axes]
            
            data = setup_data(loc_ind, years, var_list)
            
            #### compute correlations
            for mi, mgroup in enumerate(month_grouping):
                axes[mi].set_title(mgroup)
                monthly_corrs = []
                for month in mgroup:
                    corr_array = []
                    for var1 in var_list:
                        this_var = data[var1][month]
                        this_corr = []
                        for var2 in var_list:
                            comp_var = data[var2][month]
                            # mean_corr = np.nanmean([np.corrcoef(line1, line2)[0,1] for line1, line2 in zip(this_var, comp_var)])
                                                ###!!!! check line1, line2 for nan
                                                        #### adjust timescales here
                            print(var1, np.shape(this_var), var2, np.shape(comp_var))   
                                                    
                            if corr_type == 0:
                                all_corr = []
                                for line1, line2 in zip(this_var, comp_var):
                                    if len(line1)!=22: line1 = np.array([np.nan]*22)
                                    if len(line2)!=22: line2 = np.array([np.nan]*22)
                                    
                                    line1 = [np.nan]*LAG + list(line1)
                                    
                                    line1a = line1[idx1:idx2]
                                    line2a = line2[idx1:idx2]
                                    
                                    mask = ~np.isnan(np.array(line1a)) & ~np.isnan(np.array(line2a))
                                    line1b =np.array(line1a)[mask]
                                    line2b =np.array(line2a)[mask]
        
                                    all_corr.append( np.corrcoef(line1b, line2b)[0,1] )
                                mean_corr = np.nanmean(all_corr)
                            elif corr_type == 1:
                                all_lines1, all_lines2 = [],[]
                                for line1, line2 in zip(this_var, comp_var):
                                    if len(line1)!=22: line1 = np.array([np.nan]*22)
                                    if len(line2)!=22: line2 = np.array([np.nan]*22)
                                    
                                    line1 = [np.nan]*LAG + list(line1)
                                    
                                    line1a = line1[idx1:idx2]
                                    line2a = line2[idx1:idx2]
                                    
                                    mask = ~np.isnan(np.array(line1a)) & ~np.isnan(np.array(line2a))
                                    line1b =np.array(line1a)[mask]
                                    line2b =np.array(line2a)[mask]
                                    
                                    if len(line1b)==(idx2-idx1) and len(line2b)==(idx2-idx1):
                                        all_lines1.append(line1b)
                                        all_lines2.append(line2b)
                                    
                                mean_line1 = np.nanmean(all_lines1, axis=0)
                                mean_line2 = np.nanmean(all_lines2, axis=0)
        
                                mean_corr = np.corrcoef(mean_line1, mean_line2)[0,1]
                                
                            this_corr.append(mean_corr)
                        corr_array.append(this_corr)
                    monthly_corrs.append(np.array(corr_array))
                    
                ### plot mean heatmap
                mean_corr = np.nanmean(monthly_corrs, axis=0)
                mask1 = np.where(mean_corr==1, True, False)
                
                sns.heatmap(mean_corr, ax = axes[mi],
                            vmin=-0.3, vmax=0.3, cmap=cmap, 
                            cbar_kws={'label':'Correlation'}, annot=True,
                            xticklabels=var_list, yticklabels=var_list,
                            mask=mask1)
            
#%% reorganized! plot


# var_list = ['sia', 't2m', 'waves', 'wind_1']
# var_list = ['ocn_prof', 'si_conc', 'ice_motion', 'sst'] #'sic_grad',

# var_list = ['sia', 't2m', 'waves', 'wind_1',
#             'ocn_prof', 'si_conc', 'ice_motion', 'sst']

var_list = ['sia','si_conc','ocn_prof','sst', 
            't2m', 'waves','wind_1','ice_motion']

var_names = {'sia':'MIZ Ice\nArea', 'si_conc':'Mean SIC', 'ocn_prof':'Upper\nOcean\nTemp.',
             'sst':'SST','t2m':'Air\nTemp.','wind_1':'Meridional\nWind',
             'ice_motion':'Ice\nMotion','waves':'Waves'}

timescale = 3

if timescale==0: # full 3-week
    idx1, idx2 = 0, 22
elif timescale==1: # 3 days
    idx1, idx2 = 7, 10
elif timescale==2: # 1 week
    idx1, idx2 = 7, 14
elif timescale==3: # 2 weeks
    idx1, idx2 = 7, 22
else:
    raise ValueError('Wrong timescale value')
    
LAG = 0

corr_type = 1 #0=indiv corrs->mean, 1=mean lines->corr

if corr_type == 0 : vmin=-0.3; vmax=0.3
elif corr_type==1: vmin=-0.75; vmax=0.75

cmap = cmo.balance.with_extremes(bad='gray')

# month_groups = [[[4,5,6,7,8],[9,10,11,12],[1,2,3]],
#                 [[10,11,12],[4,5,6,7],[8,9]]
#                 ]
# month_groups = [[np.arange(1,12+1)],[np.arange(1,12+1)]]
# month_groups = [[[4,5,6,7,8],[9,10,11,12]],
#                 [[10,11,12],[4,5,6,7]]
#                 ]

month_groups = [[[9,10,11,12], [4,5,6,7,8]],
                [[4,5,6,7], [10,11,12]]
                ]
mgroup_names = ['Ice-Increasing', 'Ice-Decreasing']

alph = iter(list(string.ascii_lowercase))

#### set up figure
fig = plt.figure(figsize = (10*len(month_groups[0]),12))
subfigs = fig.subfigures(2, len(month_groups)) # hemis, months

fig.subplots_adjust(wspace=0.25)

for i1, title1 in enumerate(mgroup_names):
    fig.text([0.175, 0.675][i1], 1.015, title1+' Months', size='xx-large', style='italic')

for loc_ind, loc in enumerate(hemi_names): 
    path1 = root_paths[loc_ind]
    
    month_grouping = month_groups[loc_ind]
    
    for mi, mgroup in enumerate(month_grouping):
   
        subfigs[loc_ind][mi].suptitle(loc+', '+str(mgroup), size='x-large')
        axes = subfigs[loc_ind][mi].subplots(1, len(decades)) # decade subplots
        if len(month_grouping)==1:axes=[axes]
        
        #### compute correlations
        for era, years in enumerate(decades):
            ystr = str(years[0])+'-'+str(years[-1])
            
            data = setup_data(loc_ind, era, years, var_list)
        
            axes[era].set_title(ystr)
            monthly_corrs = []
            for month in mgroup:
                corr_array = []
                for var1 in var_list:
                    this_var = data[var1][month]
                    this_corr = []
                    for var2 in var_list:
                        comp_var = data[var2][month]
                    
                        if corr_type == 0:
                            all_corr = []
                            for line1, line2 in zip(this_var, comp_var):
                                if len(line1)!=22: line1 = np.array([np.nan]*22)
                                if len(line2)!=22: line2 = np.array([np.nan]*22)
                                
                                line1 = [np.nan]*LAG + list(line1)
                                
                                line1a = line1[idx1:idx2]
                                line2a = line2[idx1:idx2]
                                
                                mask = ~np.isnan(np.array(line1a)) & ~np.isnan(np.array(line2a))
                                line1b =np.array(line1a)[mask]
                                line2b =np.array(line2a)[mask]
    
                                all_corr.append( np.corrcoef(line1b, line2b)[0,1] )
                            mean_corr = np.nanmean(all_corr)
                        elif corr_type == 1:
                            all_lines1, all_lines2 = [],[]
                            for line1, line2 in zip(this_var, comp_var):
                                if len(line1)!=22: line1 = np.array([np.nan]*22)
                                if len(line2)!=22: line2 = np.array([np.nan]*22)
                                
                                line1 = [np.nan]*LAG + list(line1)
                                
                                line1a = line1[idx1:idx2]
                                line2a = line2[idx1:idx2]
                                
                                mask = ~np.isnan(np.array(line1a)) & ~np.isnan(np.array(line2a))
                                line1b =np.array(line1a)[mask]
                                line2b =np.array(line2a)[mask]
                                
                                if len(line1b)==(idx2-idx1) and len(line2b)==(idx2-idx1):
                                    all_lines1.append(line1b)
                                    all_lines2.append(line2b)
                                
                            mean_line1 = np.nanmean(all_lines1, axis=0)
                            mean_line2 = np.nanmean(all_lines2, axis=0)
    
                            mean_corr = np.corrcoef(mean_line1, mean_line2)[0,1]
                            
                        this_corr.append(mean_corr)
                    corr_array.append(this_corr)
                monthly_corrs.append(np.array(corr_array))
                
            ### plot mean heatmap
            mean_corr = np.nanmean(monthly_corrs, axis=0)
            mask1 = np.where(mean_corr>0.99999, True, False)
            
            ax = sns.heatmap(mean_corr, ax = axes[era],
                            vmin=vmin, vmax=vmax, cmap=cmap, cbar = False,
                            cbar_kws={'label':'Correlation'}, annot=True,
                            xticklabels=[var_names[v] for v in var_list], #var_list, 
                            yticklabels=[var_names[v] for v in var_list], #var_list,
                            mask=mask1)
            ax.tick_params(axis='x', rotation=45)
            ax.tick_params(axis='y', rotation=0)
            
            ax.text(0.0, 1.015, '('+next(alph)+')', 
                    transform=ax.transAxes, 
                    zorder=50)
            
            
    
# add separate colorbar
pcm = axes[era].pcolormesh(np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,2)), 
                           cmap=cmap, vmin=vmin, vmax=vmax)
cax1 = fig.add_axes([0.33,-0.055,0.33,0.033]) 
cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
cax1.set_xlabel('Correlation');
# cax1.tick_params(labelsize=10);
            

#%% combine decades?


var_list = ['sia','si_conc','ocn_prof','sst', 
            't2m', 'waves','wind_1','ice_motion']

var_names = {'sia':'MIZ Ice\nArea', 'si_conc':'Mean SIC', 'ocn_prof':'Upper\nOcean\nTemperature',
             'sst':'SST','t2m':'Air\nTemperature','wind_1':'Meridional\nWind',
             'ice_motion':'Ice\nMotion','waves':'Waves'}

timescale = 3

if timescale==0: # full 3-week
    idx1, idx2 = 0, 22
elif timescale==1: # 3 days
    idx1, idx2 = 7, 10
elif timescale==2: # 1 week
    idx1, idx2 = 7, 14
elif timescale==3: # 2 weeks
    idx1, idx2 = 7, 22
else:
    raise ValueError('Wrong timescale value')
    
LAG = 0

corr_type = 1 #0=indiv corrs->mean, 1=mean lines->corr

if corr_type == 0 : vmin=-0.3; vmax=0.3
elif corr_type==1: vmin=-0.75; vmax=0.75

cmap = cmo.balance.with_extremes(bad='gray')

month_groups = [[[9,10,11,12], [4,5,6,7,8]],
                [[4,5,6,7], [10,11,12]]
                ]
mgroup_names = ['Ice-Increasing', 'Ice-Decreasing']


alph = iter(list(string.ascii_lowercase))

#### set up figure
fig = plt.figure(figsize = (10*len(month_groups[0]),12))
subfigs = fig.subfigures(2, len(month_groups)) # hemis, months

fig.subplots_adjust(wspace=0.25)

for i1, title1 in enumerate(mgroup_names):
    fig.text([0.175, 0.675][i1], 1.015, title1+' Months', size='xx-large', style='italic')

for loc_ind, loc in enumerate(hemi_names): 
    path1 = root_paths[loc_ind]
    
    month_grouping = month_groups[loc_ind]
    
    for mi, mgroup in enumerate(month_grouping):
   
        subfigs[loc_ind][mi].suptitle(loc+', '+str(mgroup), size='x-large')
        axes = subfigs[loc_ind][mi].subplots(1, 1) # combine decades
        axes.set_title('Both Decades Combined (1982-1991, 2010-2019)')
        if len(month_grouping)==1:axes=[axes]
        
        #### combine decade data ###!!!
        data_combined = {v:{mm:[] for mm in months} for v in var_list}
        for era, years in enumerate(decades):
            data = setup_data(loc_ind, era, years, var_list)
            for var0 in var_list:
                for mm in months:
                    for l in data[var0][mm]:
                        data_combined[var0][mm].append(l)
        data = data_combined
    
        #### compute correlations
        monthly_corrs = []
        for month in mgroup:
            corr_array = []
            for var1 in var_list:
                this_var = data[var1][month]
                this_corr = []
                for var2 in var_list:
                    comp_var = data[var2][month]
                
                    if corr_type == 0:
                        all_corr = []
                        for line1, line2 in zip(this_var, comp_var):
                            if len(line1)!=22: line1 = np.array([np.nan]*22)
                            if len(line2)!=22: line2 = np.array([np.nan]*22)
                            
                            line1 = [np.nan]*LAG + list(line1)
                            
                            line1a = line1[idx1:idx2]
                            line2a = line2[idx1:idx2]
                            
                            mask = ~np.isnan(np.array(line1a)) & ~np.isnan(np.array(line2a))
                            line1b =np.array(line1a)[mask]
                            line2b =np.array(line2a)[mask]

                            all_corr.append( np.corrcoef(line1b, line2b)[0,1] )
                        mean_corr = np.nanmean(all_corr)
                    elif corr_type == 1:
                        all_lines1, all_lines2 = [],[]
                        for line1, line2 in zip(this_var, comp_var):
                            if len(line1)!=22: line1 = np.array([np.nan]*22)
                            if len(line2)!=22: line2 = np.array([np.nan]*22)
                            
                            line1 = [np.nan]*LAG + list(line1)
                            
                            line1a = line1[idx1:idx2]
                            line2a = line2[idx1:idx2]
                            
                            mask = ~np.isnan(np.array(line1a)) & ~np.isnan(np.array(line2a))
                            line1b =np.array(line1a)[mask]
                            line2b =np.array(line2a)[mask]
                            
                            if len(line1b)==(idx2-idx1) and len(line2b)==(idx2-idx1):
                                all_lines1.append(line1b)
                                all_lines2.append(line2b)
                            
                        mean_line1 = np.nanmean(all_lines1, axis=0)
                        mean_line2 = np.nanmean(all_lines2, axis=0)

                        mean_corr = np.corrcoef(mean_line1, mean_line2)[0,1]
                        
                    this_corr.append(mean_corr)
                corr_array.append(this_corr)
            monthly_corrs.append(np.array(corr_array))
            
        ### plot mean heatmap
        mean_corr = np.nanmean(monthly_corrs, axis=0)
        mask1 = np.where(mean_corr>0.99999, True, False)
        
        ax = sns.heatmap(mean_corr, ax = axes,
                        vmin=vmin, vmax=vmax, cmap=cmap, cbar = False,
                        cbar_kws={'label':'Correlation'}, annot=True,
                        xticklabels=[var_names[v] for v in var_list], #var_list, 
                        yticklabels=[var_names[v] for v in var_list], #var_list,
                        mask=mask1)
        ax.tick_params(axis='x', rotation=0)
        ax.tick_params(axis='y', rotation=0)
        
        ax.text(0.0, 1.015, '('+next(alph)+')', 
                transform=ax.transAxes, 
                zorder=50)
            
            
    
# add separate colorbar
pcm = axes.pcolormesh(np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,2)), 
                           cmap=cmap, vmin=vmin, vmax=vmax)
cax1 = fig.add_axes([0.33,-0.055,0.33,0.033]) 
cbar1 = fig.colorbar(pcm, cax=cax1, orientation='horizontal')
cax1.set_xlabel('Correlation');
# cax1.tick_params(labelsize=10);

#%% end
