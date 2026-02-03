#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 2025
storm_locations1.py

--> adapted from teststuff/cyclone_analysis_v4.py (figure 2 doi:10.1175/JCLI-D-24-0026.1)

- maps of storm starts
- heat maps of min pressure (like cesm paper)
- plot storm tracks (more zonal in southern hemisphere?)

@author: mundi
"""
#%% imports and files
import numpy as np
import xarray as xr
import calendar
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import functions as fx
import cmocean.cm as cmo
from calendar import monthrange

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

ice_path = '/Users/mundi/Desktop/seaice/'
ice_paths = [ice_path, ice_path+'south/']

decades = [np.arange(2010,2020),np.arange(1982, 1992)]

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]
month_labels = [calendar.month_abbr[mm] for mm in months]+[calendar.month_abbr[months[0]]]

month_colors = ['#238443','#78c679','#c2e699',
                '#d7b5d8','#df65b0','#dd1c77','#980043','#7a0177',
                '#253494','#2c7fb8','#41b6c4','#a1dab4'
                ]

hemi_names= ['Arctic', 'Antarctic']

import pickle
import cartopy.feature as cfeature
from matplotlib.cm import ScalarMappable
import matplotlib.colors as colors
era_ystr = [str(y[0])+'-'+str(y[-1]) for y in decades]
ice_lims = [20,80]

loc_ind = 0
path1 = root_paths[loc_ind]

#%%% functions

def monthly_grouping(year, months, startdate, monthly_groupings, path1):
    # open ice area
    ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
    ice_area80 = ds_area['ice_area80'].values
    box_area = ds_area['box_area'].values
    ice_sorter = ice_area80
    
    month_dict = {}
    for month in months:
        month_dict[month] = []
    for di, dt in enumerate(startdate):
       
        ice_frac = ice_sorter[di]*100/box_area[di]
        if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
            continue
            
        month_dict[dt.month].append(pressure[di])
        
    # flatten to list
    month_list = []
    for month in months:
        month_list.append(month_dict[month])
    monthly_groupings.append(month_list)
    
    return monthly_groupings

#%% make figure : map of storm starts
letters1 = iter(['c','d','e','f'])
letters2 = iter(['a','b'])

cmap = cmo.dense_r #cmo.ice_r #plt.cm.Blues
norm = plt.Normalize(vmin=950, vmax=990)
markersize=12

## sea ice plotting
si, si_lon, si_lat = fx.load_seaice(ice_paths[0], 2010,8,1, latlon=True)

norm_si = colors.Normalize(0,1)
color_map = [[norm_si(0), [0,0,0,0]],
            [norm_si(1), [211/255,211/255,211/255,0.5]]] 
si_cmap = colors.LinearSegmentedColormap.from_list("", color_map,N=2)

fig = plt.figure(figsize=[25,15])
gs = fig.add_gridspec(2, 3)
gs.update(hspace=0.3)
extent = [-160,90,50,60]

axes_all = []
for row in [0,1]:
    axes_col = []
    for col in [1,2]:
        ax = fig.add_subplot(gs[row, col], projection=ccrs.NorthPolarStereo(central_longitude=-45))
        ax.coastlines('50m',edgecolor='black',linewidth=0.75)
        ax.set_extent(extent, ccrs.PlateCarree())
        gl = ax.gridlines(draw_labels=False)
        gl.xlabels_top = False
        gl.xlabels_bottom = False
        ax.add_feature(cfeature.LAND, facecolor='0.75')
        ax.add_feature(cfeature.LAKES, facecolor='0.85')
        
        ax.text(0.0075, 1.055, next(letters1), transform=ax.transAxes, fontsize=20, 
                bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
                      'edgecolor':'k', 'lw':0.75},zorder=50)
        
        axes_col.append(ax)
    axes_all.append(axes_col)
axes_all = np.array(axes_all)

def slp_dist_plot(monthly_data, month, ax, fontsize=10, plot_cases=False):
    from scipy.stats import ttest_ind, ks_2samp, mannwhitneyu
    colors = ['#ca0020','#0571b0'] #RdBl
    
    if type(month)==int:
        month_index=[month-6]
        mname = calendar.month_name[month]
    elif type(month)==list or type(month)==np.array:
        month_index = [m-6 for m in month]
        if len(month)==2:
            m1 = calendar.month_name[month[0]]
            m2 = calendar.month_name[month[1]]
            if len(m1) > 4:
                m1 = m1[0:3]+'.'
            if len(m2) > 4:
                m2 = m2[0:3]+'.'
            mname = m1 +' and '+ m2
        else:
            mname = calendar.month_name[month[0]]+' - '+ calendar.month_name[month[-1]]
    
    era_data = []
    for era in [0,1]:
        data = []
        for midx in month_index:
            data = data + list(monthly_data[era][midx])
        
        mean_slp = np.mean(data)
        ax.axvline(mean_slp, color=colors[era], lw=2)
        
        bin_num = np.arange(955,987,2)#15 #12
        
        print(era_ystr[era]+': '+str(np.nanmin(data))+','+str(np.nanmax(data)) )
        
        hist, bins = np.histogram(data, bins=bin_num, density=True)
        bin_centers = (bins[1:]+bins[:-1])*0.5
        ax.plot(bin_centers, hist, color = colors[era], lw=2.5,
                label=era_ystr[era] + '\nn=' + str(len(data))+r', $\mu=$'+
                str(round(mean_slp,1))
                )
        
        era_data.append(data)
            
    t, pt = ttest_ind(era_data[0], era_data[1])
    ks, pks = ks_2samp(era_data[0], era_data[1])
    mw, pmw = mannwhitneyu(era_data[0], era_data[1])
    
    pvals = [pt, pmw, pks]
    cl = [100*(1-p) for p in pvals]
    
    stats_text = 't: '+ str(round(cl[0],1)) +'%, MW: '+str(round(cl[1],1))+'%, KS: '+str(round(cl[-1],1))+'%'
    
    ax.legend(fontsize=fontsize, loc='upper left')
    ax.set_title(mname + ' Intensity' +'\n'+stats_text, ha='left', x=-0, fontsize=fontsize)
    ax.set_xlabel('Pressure [hPa]', fontsize=fontsize)
    ax.set_ylabel('Density', fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.set_ylim([0,0.325])
    ax.set_xlim([955,985])
    
    return ax, [pt, pks, pmw]

monthly_data = []
monthstrs = ['June', 'July', 'Aug', 'Sep', 'Oct']
monthstrs2 = ['June and July', 'Aug. and Sep.']
for era, years in enumerate(decades):
    print(' ')
    storm_counter = 0
    # ices = ice_extents[era] ###!!!!
    
    if len(years)>1:
        titlestr = ' '+ str(years[0])+'-'+str(years[-1])+''
    else: 
        titlestr = ' '+ str(years[0])+''
        
    ticks, monthly_groupings = [], []
    p67, p89 = [],[]
    c67, c89 = 0,0
    for year in years:
        print(year)
        ticks.append(str(year))
        # load census and sea ice info
        census_name = 'census_'+ str(year) +'.csv'
        [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
            fx.readCensus(path1+'census/'+census_name, convertDT=True)
        # reorder data into monthly gorupings
        monthly_groupings = monthly_grouping(year, months, startdate, 
                                            monthly_groupings,path1)
        
       
        # open ice area
        ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
        ice_area80 = ds_area['ice_area80'].values
        ice_area15 = ds_area['ice_area15'].values
        box_area = ds_area['box_area'].values
        ice_sorter = ice_area80
        
        ### load storm info
        picklename = str(year) + '_all_storms.pkl'
        file = open(path1+'contours/'+picklename, 'rb')
        storm_data = pickle.load(file)
        file.close()
        
        ### load more info
        storm_ranges = []
        for startdt, enddt in zip(startdate, enddate):
            storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
        
        for sn, event in enumerate(startdate):
            ice_frac = ice_sorter[sn]*100/box_area[sn]             
            if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                continue
            
            p1=pressure[sn]
            p = p1
            x1,y1=startlon[sn], startlat[sn]
            t1 = startdate[sn]
            t=t1
            
            color = cmap(norm(p1))
            
            # get storm index
            j=7
            storm_indices = []
            for k in storm_ranges[sn]: 
               if k != '-':
                   storm_indices.append(j)
                   j+=1
                
            if t.month in [6,7]:
                axes_all[0][era].plot(x1,y1, color='k', marker='o', 
                                      markersize=markersize+2, zorder=50,
                                      transform=ccrs.PlateCarree())
                axes_all[0][era].plot(x1,y1, color=color, marker='o', alpha=1,
                                      markersize=markersize, zorder=50,
                                      transform=ccrs.PlateCarree())
                p67.append(p)
                c67+=1
            elif t.month in [8,9]:
                axes_all[1][era].plot(x1,y1, color='k', marker='o', 
                                      markersize=markersize+2, zorder=50,
                                      transform=ccrs.PlateCarree())
                axes_all[1][era].plot(x1,y1, color=color, marker='o', alpha=1,
                                      markersize=markersize, zorder=50,
                                      transform=ccrs.PlateCarree())
                p89.append(p)
                c89+=1

    # storm by month
    mdata = []
    storm_count = 0
    for mi, mm in enumerate(monthstrs):
         mdata_list = [item[mi] for item in monthly_groupings]
         mdata_flat = [item for sublist in mdata_list for item in sublist]
         mdata.append(mdata_flat)
         storm_count += len(mdata_flat)
    monthly_data.append(mdata)
        
    
    axes_all[0][era].set_title(monthstrs2[0]+titlestr+', n='+str(c67), fontsize=20)
    axes_all[1][era].set_title(monthstrs2[1]+titlestr+', n='+str(c89), fontsize=20)
    
    '''
    axes_all[0][era].pcolormesh(si_lon, si_lat, np.where(ices[0]<0.15,0,1), cmap=si_cmap, transform=ccrs.PlateCarree())
    axes_all[1][era].pcolormesh(si_lon, si_lat, np.where(ices[1]<0.15,0,1), cmap=si_cmap, transform=ccrs.PlateCarree())
    fx.plot_geocontour(axes_all[0][era], si_lon, si_lat, ices[0], levels=[0.15], color='k', lw=2, ls='solid')
    fx.plot_geocontour(axes_all[1][era], si_lon, si_lat, ices[1], levels=[0.15], color='k', lw=2, ls='solid')
    fx.plot_geocontour(axes_all[0][era], si_lon, si_lat, ices[0], levels=[0.80], color='k', lw=1.75, ls='dotted')
    fx.plot_geocontour(axes_all[1][era], si_lon, si_lat, ices[1], levels=[0.80], color='k', lw=1.75, ls='dotted')
    '''
    
### Add in ColorBar! (new cmap, check pressure bounds)
s_m = ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])

# fig.subplots_adjust(right=0.95)
cbar_ax = fig.add_axes([0.425, 0.075, 0.45, 0.025])
cbar = fig.colorbar(s_m, cax=cbar_ax, orientation='horizontal')
cbar.ax.tick_params(labelsize=18)
cbar.set_label('Minimum Pressure (hPa)', fontsize=20)


f67 = fig.add_subplot(gs[0, 0])
f67.text(-0.1, 1.055, next(letters2), transform=f67.transAxes, fontsize=20, 
        bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
              'edgecolor':'k', 'lw':0.75},zorder=50)
f89 = fig.add_subplot(gs[1, 0])
f89.text(-0.1, 1.055, next(letters2), transform=f89.transAxes, fontsize=20, 
        bbox={'facecolor': 'white', 'alpha': 1, 'pad':5, 
              'edgecolor':'k', 'lw':0.75},zorder=50)

f67, p_vals67 = slp_dist_plot(monthly_data, [6,7], ax=f67, fontsize=20)
f89, p_vals89 = slp_dist_plot(monthly_data, [8,9], ax=f89, fontsize=20)

handles, labels = f89.get_legend_handles_labels()
if len(labels)==4:
    order = [0,3,1,2]
    f89.legend([handles[idx] for idx in order],[labels[idx] for idx in order],fontsize=20,loc='upper left')

    
#%%% update seasons

plot_ice = True

for loc_ind, loc in enumerate(hemi_names):
    print(loc)
    path1 = root_paths[loc_ind]
    
    markersize=12
    cmap = cmo.dense_r #cmo.ice_r #plt.cm.Blues
    if loc_ind==0: 
        norm = plt.Normalize(vmin=950, vmax=990)
        si, si_lon, si_lat = fx.load_seaice(ice_paths[loc_ind], 2010,8,1, latlon=True)

    elif loc_ind==1: 
        norm = plt.Normalize(vmin=930, vmax=960)
        si, si_lon, si_lat = fx.load_seaice_sh(ice_paths[loc_ind], 2010,8,1, latlon=True)
    
    #### set up plot
    fig = plt.figure(figsize=[15,20])
    gs = fig.add_gridspec(3, 2)
    # gs.update(hspace=0.3)
    fig.suptitle('\n'+loc, fontweight='bold', fontsize=22)
    
    axes_all = []
    for row in [0,1,2]:
        axes_col = []
        for col in [0,1]:
            if loc_ind==0:
                ax = fig.add_subplot(gs[row, col], projection=ccrs.NorthPolarStereo(central_longitude=0))
                ax.set_extent([-160,90,50,60], ccrs.PlateCarree())
            elif loc_ind==1:
                ax = fig.add_subplot(gs[row, col], projection=ccrs.SouthPolarStereo())
                ax.set_extent([0,360,-90,-55], ccrs.PlateCarree())
                
            ax.coastlines('50m',edgecolor='black',linewidth=0.75)
            gl = ax.gridlines(draw_labels=False)
            gl.xlabels_top = False
            gl.xlabels_bottom = False
            ax.add_feature(cfeature.LAND, facecolor='0.75')
            ax.add_feature(cfeature.LAKES, facecolor='0.85')
            axes_col.append(ax)
        axes_all.append(axes_col)
    axes_all = np.array(axes_all)
    
    #### data
    
    if loc_ind==0: month_groups = [[9,10,11,12],[4,5,6,7,8],[1,2,3]]
    elif loc_ind==1: month_groups = [[4,5,6],[10,11,12],[8,9]]
    month_group_names = ['Fall', 'Spring/Summer', 'Winter']
    
    for era, years in enumerate(decades):
        print(' ')
        storm_counter = 0
        # ices = ice_extents[era] ###!!!!
        
        if len(years)>1:
            titlestr = ' '+ str(years[0])+'-'+str(years[-1])+''
        else: 
            titlestr = ' '+ str(years[0])+''
            
        ticks, monthly_groupings = [], []
        monthly_pressures = {mg:[] for mg in month_group_names}
        monthly_counts = {mg:0 for mg in month_group_names}
        monthly_ice = {mg:[] for mg in month_group_names}
        for year in years:
            print(year)
            ticks.append(str(year))
            # load census and sea ice info
            census_name = 'census_'+ str(year) +'.csv'
            [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
                fx.readCensus(path1+'census/'+census_name, convertDT=True)
    
            # open ice area
            ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
            ice_area80 = ds_area['ice_area80'].values
            ice_area15 = ds_area['ice_area15'].values
            box_area = ds_area['box_area'].values
            ice_sorter = ice_area80
            
            if plot_ice:
                print(str(year)+' ice', end='')
                for month in np.arange(1,12+1):
                    print('.', end='')
                    days = np.arange(1,monthrange(year,month)[-1]+1)
                    for day in days:
                        if loc_ind==0: si = fx.load_seaice(ice_paths[loc_ind], year, month, day, latlon=False)
                        elif loc_ind==1: si = fx.load_seaice_sh(ice_paths[loc_ind], year, month, day, latlon=False)

                        for mi, month_group in enumerate(month_groups):
                            if month in month_group:
                                monthly_ice[month_group_names[mi]].append(si)
                                break
                print(' done!')
            
            
            ### load more info
            storm_ranges = []
            for startdt, enddt in zip(startdate, enddate):
                storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
            
            for sn, event in enumerate(startdate):
                ice_frac = ice_sorter[sn]*100/box_area[sn]             
                if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                    continue
                
                p = pressure[sn]
                # x1,y1=startlon[sn], startlat[sn]
                x1 = np.nanmean([startlon[sn], endlon[sn]])
                y1 = np.nanmean([startlat[sn], endlat[sn]])
                t = startdate[sn]
                
                color = cmap(norm(p))
                       
                # check what season it is in
                for mi, month_group in enumerate(month_groups):
                    if t.month in month_group:
                        axes_all[mi][era].plot(x1,y1, color='k', marker='o', 
                                              markersize=markersize+2, zorder=50,
                                              transform=ccrs.PlateCarree())
                        axes_all[mi][era].plot(x1,y1, color=color, marker='o', alpha=1,
                                              markersize=markersize, zorder=50,
                                              transform=ccrs.PlateCarree())
                       
                        monthly_pressures[month_group_names[mi]].append(p)
                        monthly_counts[month_group_names[mi]]+=1
                        break
            
        for mi, mg_name in enumerate(month_group_names):
            axes_all[mi][era].set_title(mg_name+titlestr+', n='+str(monthly_counts[mg_name]), fontsize=20)
    
            if plot_ice: #### ice contour/areas
                ices = np.nanmean(monthly_ice[mg_name], axis=0)
                axes_all[mi][era].pcolormesh(si_lon, si_lat, np.where(ices<0.15,0,1), cmap=si_cmap, transform=ccrs.PlateCarree())
                fx.plot_geocontour(axes_all[mi][era], si_lon, si_lat, ices, levels=[0.15], color='k', lw=2, ls='solid')
                fx.plot_geocontour(axes_all[mi][era], si_lon, si_lat, ices, levels=[0.80], color='k', lw=1.75, ls='dotted')
        
    
    ### Add in ColorBar! (new cmap, check pressure bounds)
    s_m = ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    cbar_ax = fig.add_axes([0.3, 0.075, 0.45, 0.025])
    cbar = fig.colorbar(s_m, cax=cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label('Minimum Pressure (hPa)', fontsize=20)

#%%% storm tracks 

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    markersize=12
    cmap = cmo.dense_r #cmo.ice_r #plt.cm.Blues
    if loc_ind==0: norm = plt.Normalize(vmin=950, vmax=990)
    elif loc_ind==1: norm = plt.Normalize(vmin=930, vmax=960)
    
    
    #### set up plot
    fig = plt.figure(figsize=[15,20])
    gs = fig.add_gridspec(3, 2)
    # gs.update(hspace=0.3)
    fig.suptitle('\n'+loc+' Storm Tracks', fontweight='bold', fontsize=22)
    
    axes_all = []
    for row in [0,1,2]:
        axes_col = []
        for col in [0,1]:
            if loc_ind==0:
                ax = fig.add_subplot(gs[row, col], projection=ccrs.NorthPolarStereo(central_longitude=0))
                ax.set_extent([-160,90,50,60], ccrs.PlateCarree())
            elif loc_ind==1:
                ax = fig.add_subplot(gs[row, col], projection=ccrs.SouthPolarStereo())
                ax.set_extent([0,360,-90,-55], ccrs.PlateCarree())
                
            ax.coastlines('50m',edgecolor='black',linewidth=0.75)
            gl = ax.gridlines(draw_labels=False)
            gl.xlabels_top = False
            gl.xlabels_bottom = False
            ax.add_feature(cfeature.LAND, facecolor='0.75')
            ax.add_feature(cfeature.LAKES, facecolor='0.85')
            axes_col.append(ax)
        axes_all.append(axes_col)
    axes_all = np.array(axes_all)
    
    #### data
    
    if loc_ind==0: month_groups = [[9,10,11,12],[4,5,6,7,8],[1,2,3]]
    elif loc_ind==1: month_groups = [[4,5,6],[10,11,12],[8,9]]
    month_group_names = ['Fall', 'Spring/Summer', 'Winter']
    
    for era, years in enumerate(decades):
        print(' ')
        storm_counter = 0
        # ices = ice_extents[era] ###!!!!
        
        if len(years)>1:
            titlestr = ' '+ str(years[0])+'-'+str(years[-1])+''
        else: 
            titlestr = ' '+ str(years[0])+''
            
        ticks, monthly_groupings = [], []
        monthly_pressures = {mg:[] for mg in month_group_names}
        monthly_counts = {mg:0 for mg in month_group_names}
        for year in years:
            print(year)
            ticks.append(str(year))
            # load census and sea ice info
            census_name = 'census_'+ str(year) +'.csv'
            [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
                fx.readCensus(path1+'census/'+census_name, convertDT=True)
    
            # open ice area
            ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
            ice_area80 = ds_area['ice_area80'].values
            ice_area15 = ds_area['ice_area15'].values
            box_area = ds_area['box_area'].values
            ice_sorter = ice_area80
            
            ### load storm info
            picklename = str(year) + '_all_storms.pkl'
            file = open(path1+'all_storms/'+picklename, 'rb')
            storm_data = pickle.load(file)
            file.close()
            if len(storm_data)<20:continue
            
            ### load more info
            storm_ranges = []
            for startdt, enddt in zip(startdate, enddate):
                storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
            
            for sn, event in enumerate(startdate):
                ice_frac = ice_sorter[sn]*100/box_area[sn]             
                if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                    continue
                
                t = startdate[sn]
                p = pressure[sn]
                
                try:
                    event_data = np.array(storm_data[sn])
                except IndexError: continue
            
                lons = event_data[:,0]
                lons = np.where(lons>180, lons-360, lons)
                lats = event_data[:,1]
                
                [lon1, lon2], [lat1, lat2] = fx.geoplot_2d(lons, lats)
                
                # check what season it is in
                for mi, month_group in enumerate(month_groups):
                    if t.month in month_group:
                        for lonx, latx in zip([lon1,lon2],[lat1,lat2]):
                            axes_all[mi][era].plot(lonx, latx, color=month_colors[t.month-1], 
                                                  lw=0.55, transform=ccrs.PlateCarree())
                       
                        monthly_pressures[month_group_names[mi]].append(p)
                        monthly_counts[month_group_names[mi]]+=1
                        break
            
        for mi, mg_name in enumerate(month_group_names):
            axes_all[mi][era].set_title(mg_name+titlestr+', n='+str(monthly_counts[mg_name]), fontsize=20)
    
for mi, mcolor in enumerate(month_colors):
    axes_all[-1][-1].plot([],[], lw=5, color=mcolor, label = month_abbrs[mi])
axes_all[-1][-1].legend(ncol=6, loc='upper left', fontsize=20,
                        bbox_to_anchor = (-1.25, -0.05),
                        handletextpad=0.25, handlelength=1,)
#%% heat maps of min pressure

for loc_ind, loc in enumerate(hemi_names):
    path1 = root_paths[loc_ind]
    
    # set up grid
    if loc_ind==0: 
        lon_values = np.arange(0,360,5)
        lat_values = np.arange(60,90,5)
    elif loc_ind==1: 
        lon_values = np.arange(0,360,5)
        lat_values = np.arange(-90,-55,5)
        
    xs = list(lon_values)+[360]
    ys = list(lat_values)+[lat_values[-1]+(lat_values[-1]-lat_values[-2])]
    
    
    #### set up plot
    fig = plt.figure(figsize=[15,20])
    gs = fig.add_gridspec(3, 2)
    # gs.update(hspace=0.3)
    fig.suptitle('\n'+loc+' Storm Tracks', fontweight='bold', fontsize=22)
    
    axes_all = []
    for row in [0,1,2]:
        axes_col = []
        for col in [0,1]:
            if loc_ind==0:
                ax = fig.add_subplot(gs[row, col], projection=ccrs.NorthPolarStereo(central_longitude=0))
                ax.set_extent([-160,90,50,60], ccrs.PlateCarree())
            elif loc_ind==1:
                ax = fig.add_subplot(gs[row, col], projection=ccrs.SouthPolarStereo())
                ax.set_extent([0,360,-90,-55], ccrs.PlateCarree())
                
            ax.coastlines('50m',edgecolor='black',linewidth=0.75)
            gl = ax.gridlines(draw_labels=False)
            gl.xlabels_top = False
            gl.xlabels_bottom = False
            ax.add_feature(cfeature.LAND, facecolor='0.75')
            ax.add_feature(cfeature.LAKES, facecolor='0.85')
            axes_col.append(ax)
        axes_all.append(axes_col)
    axes_all = np.array(axes_all)
    
    #### data
    
    if loc_ind==0: month_groups = [[9,10,11,12],[4,5,6,7,8],[1,2,3]]
    elif loc_ind==1: month_groups = [[4,5,6],[10,11,12],[8,9]]
    month_group_names = ['Fall', 'Spring/Summer', 'Winter']
    
    for era, years in enumerate(decades):
        print(' ')
        storm_counter = 0
        # ices = ice_extents[era] ###!!!!
        
        if len(years)>1:
            titlestr = ' '+ str(years[0])+'-'+str(years[-1])+''
        else: 
            titlestr = ' '+ str(years[0])+''
            
        ticks, monthly_groupings = [], []
        monthly_pressures = {mg:[] for mg in month_group_names}
        monthly_counts = {mg:0 for mg in month_group_names}
        
        grids = {mg:np.zeros((len(lon_values),len(lat_values))) for mg in month_group_names}
        
        for year in years:
            print(year)
            ticks.append(str(year))
            # load census and sea ice info
            census_name = 'census_'+ str(year) +'.csv'
            [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
                fx.readCensus(path1+'census/'+census_name, convertDT=True)
    
            # open ice area
            ds_area = xr.open_dataset(path1+'area/' + str(year) +'_area.nc')
            ice_area80 = ds_area['ice_area80'].values
            ice_area15 = ds_area['ice_area15'].values
            box_area = ds_area['box_area'].values
            ice_sorter = ice_area80
            
            ### load storm info
            picklename = str(year) + '_all_storms.pkl'
            file = open(path1+'all_storms/'+picklename, 'rb')
            storm_data = pickle.load(file)
            file.close()
            if len(storm_data)<20:continue
            
            ### load more info
            storm_ranges = []
            for startdt, enddt in zip(startdate, enddate):
                storm_ranges.append(fx.daterange(startdt, enddt, dt=24))  
            
            for sn, event in enumerate(startdate):
                ice_frac = ice_sorter[sn]*100/box_area[sn]             
                if ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims):
                    continue
                
                t = startdate[sn]
                p = pressure[sn]
                
                try:
                    event_data = np.array(storm_data[sn])
                except IndexError: continue
            
                lons = event_data[:,0]
                # lons = np.where(lons>180, lons-360, lons)
                lats = event_data[:,1]
                
                
                # check what season it is in
                for mi, month_group in enumerate(month_groups):
                    if t.month in month_group:
                        
                        for x, y in zip(lons, lats):
                            # xloc = np.where(np.floor(x)==lon_values)
                            # yloc = np.where(np.floor(y)==lat_values)
                            
                            xloc = np.where(lon_values == max([i for i in lon_values if x >= i]))
                            yloc = np.where(lat_values == max([i for i in lat_values if y >= i]))
                            
                            grids[month_group_names[mi]][xloc, yloc]+=1
                       
                        monthly_pressures[month_group_names[mi]].append(p)
                        monthly_counts[month_group_names[mi]]+=1
                        break
            
        for mi, mg_name in enumerate(month_group_names):
            axes_all[mi][era].set_title(mg_name+titlestr+', n='+str(monthly_counts[mg_name]), fontsize=20)
    
            X, Y= np.meshgrid(xs, ys)
    
            axes_all[mi][era].pcolormesh(X.T, Y.T, grids[mg_name],
                                         cmap=cmo.amp, vmin=0, vmax=400,
                                         transform=ccrs.PlateCarree())
    


#%%% maps (interpolated)
'''
from datetime import datetime

def local_era5(slp_file, year, month, days, idx=6):
    if type(days)==int: days=[days]
    startdate = str(year)+'-'+LZ(month)+'-'+days[0]
    enddate = str(year)+'-'+LZ(month)+'-'+days[-1]
    
    ds = xr.open_dataset(slp_file)
    
    msl = ds['msl'].sel(time=slice(startdate, enddate))/100
    msl, cyc_x = add_cyclic_point(msl, coord=ds['longitude'].values)
    
    lon, lat = np.meshgrid(cyc_x, ds['latitude'].values)
    time = ds['time']
    
    msl_out, msl6, time6 = [],[],[]
    count = 0
    for xi, msl1 in enumerate(msl):
        if count == 0: 
            timestr = str(time[xi].values).split('.')[0]
            time6.append( datetime.strptime(timestr,'%Y-%m-%dT%H:%M:%S') )
        msl6.append(msl1) #.values)
        count += 1
        if count ==idx:
            msl_out.append(np.nanmean(msl6,axis=0))
            count=0; msl6 = []
    time6 = np.array(time6)
    ds.close()
    return time6, lon, lat, msl_out


years = [2010]
# years = np.arange(2011,2019+1)
# years = np.arange(1982, 1991+1)


for year in years:

    time = fx.daterange(datetime(year,6,1), datetime(year,9,30,23), dt=24)

    era_msl_file = era_path + 'msl_'+ str(year)+'-' # +'-6.nc'

    print('Starting '+str(year)+' ERA interpolation: ', end='')
    interp_start = time.time()
    
    ### get time vector
    time = fx.daterange(datetime(year,6,1), datetime(year,9,30,23), dt=6)
    
    ### load pressure data
    interp_msl = []
    for month in months:
        print(month, end=' ')
        calendar_days = list(np.arange(1,calendar.monthrange(year, month)[1]+1))
        alldays = [str(d) if d>=10 else '0'+str(d) for d in calendar_days]
        
        TIME, lon, lat, MSL_ERA = local_era5(era_msl_file+str(month)+'.nc', 
                                             year, month, alldays, idx=6)
        
        
        for msl_era in MSL_ERA:
            new_grid = griddata(np.array([lon.flatten(), lat.flatten()]).T, msl_era.flatten(),
                                np.array([cesm_lon.flatten(), cesm_lat.flatten()]).T)
            interp_msl.append( np.reshape(new_grid, np.shape(cesm_lon)) )
        
    print();print('... Done: ' + str(round((time.time() - interp_start)/60, 2))+' min')
    
    ### make data array for export
    da = xr.DataArray(
        np.array(interp_msl),
        dims=("time", "y", "x"),
        coords=[
            ("time", time),
            ("y", clat),
            ("x", clon)
        ],
    )
    da.to_netcdf(output_path+str(year)+'_arctic.nc')
'''

#%% plot storm tracks 





#%% end
