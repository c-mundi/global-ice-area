#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 2025
seaice_combine.py

merge *_seaice.nc for summer months and new runs

@author: mundi
"""
#%% imports and filepaths
import numpy as np
import xarray as xr
from datetime import datetime
import os

run = {'seaice':True, 'area':False}

years = list(np.arange(1982,1992)) + list(np.arange(2010,2020))

census_path0 = '/Users/mundi/Desktop/month-hemi/nh_data/census0/'

root_path = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/'#+'p-adjust/'
si_path = root_path + 'seaice/'

var_addons = [['_1000', '_990', '_contour'],['2_1000', '2_990', '2_contour']]
addon_num = 0 # storm area
sia_var_addons = var_addons
miz_ind = 1 #[0=daily, 1=total]
miz_names = ['daily_', 'total_']

summer_months = [6,7,8,9]


#%%% functions
def readCensus(file, convertDT=False):
    import csv
    
    # Load file
    csv_file = open(file,'r')
    startdate, enddate = [],[]
    startlon, endlon = [], []
    startlat, endlat = [],[]
    pressure = []
    
    # Read off and discard first line, to skip headers
    csv_file.readline()
    
    # Split columns while reading
    for a, b, c, d, e, f, g in csv.reader(csv_file, delimiter=','):
        # Append each variable to a separate list
        startdate.append(a) 
        startlat.append(float(b))
        startlon.append(float(c))
        pressure.append(float(d))
        enddate.append(e)
        endlat.append(float(f))
        endlon.append(float(g))
    csv_file.close()
    
    if convertDT:
        startDT, endDT = [],[]                            
        for i, pres in enumerate(pressure):
            startDT.append( datetime(int(startdate[i][:4]), int(startdate[i][5:7]), 
                                           int(startdate[i][8:10]), 0) )
                           # int(startdate[i][-2:]))
            endDT.append( datetime(int(enddate[i][:4]), int(enddate[i][5:7]), 
                                           int(enddate[i][8:10]), 0) )
                            # int(enddate[i][-2:])
                            
        # startdate, enddate = startDT, endDT
        
        ### sort times in ascending order
        startdate, enddate, startlon, startlat, endlon, endlat, pressure = \
            map(list, zip(*sorted(zip(startDT, endDT, startlon, startlat, endlon, endlat, pressure))))

    
    dates = [startdate, enddate]
    coords = [[startlon, startlat],[endlon,endlat]]
    return dates, coords, pressure 

#%% data merge (seaice)
 # sia = ds['sia_miz'+sia_var_addons[miz_ind][addon_num]].values[storm_num]

if run['seaice']:
    
    savepath = '/Users/mundi/Desktop/month-hemi/nh_data/seaice/'
    if not os.path.exists(savepath):
        os.makedirs(savepath)
 
    for year in years:
        si_series = []
        si_clim_series = []
        storm_count = 0
        load_summer=True
        
        # open census and sea ice
        [startdate, enddate] = readCensus(census_path0+'census_'+str(year)+'.csv', convertDT=True)[0]
        ds = xr.open_dataset(si_path+str(year)+'_seaice.nc')
        sia = ds['sia_miz'+sia_var_addons[miz_ind][addon_num]].values
        sia_clim = ds['sia_clim_miz'+sia_var_addons[miz_ind][addon_num]].values
        ds.close()
        
        print(year, len(startdate), len(sia), len(sia_clim))
        
        # storm loop
        for storm_num, sd in enumerate(startdate):
            if sd.month not in summer_months:
                try:
                    si_series.append(sia[storm_num])
                    si_clim_series.append(sia_clim[storm_num])
                except IndexError:
                    continue
                storm_count += 1
            else:
                if load_summer:
                    
                    census_path1 = '/Users/mundi/Desktop/teststuff/original_census/'
                    census_file1 = census_path1+'census_'+str(year)+'.csv'
                    [startdate1, enddate1] = readCensus(census_file1, convertDT=True)[0]
                    
                    try:
                        ds_si1 = xr.open_dataset('/Users/mundi/Desktop/FINAL/' + str(year) +'_seaice.nc')
                        sia1 = ds_si1['sia_miz'+sia_var_addons[miz_ind][addon_num]].values
                        sclim1 = ds_si1['sia_clim_miz'+sia_var_addons[miz_ind][addon_num]].values
                        ds_si1.close()
                    except:
                        print('- skip: '+'/Users/mundi/Desktop/FINAL/'+ str(year) +'_seaice.nc')
                        continue
                    
                    for storm_num1, sd1 in enumerate(startdate1):
                        if sd.month in summer_months:
                            si_series.append(sia1[storm_num1])
                            si_clim_series.append(sclim1[storm_num1])
                            storm_count += 1
                    
                    load_summer=False
            
        #### save new file
        # names of the dimensions in the required order
        dims = ('nstorms', 'time')
        # create coordinates to use for indexing along each dimension
        coords = {'time' :  np.arange(-7,14+1),
                  'nstorms' : np.arange(0,storm_count)}
        # attributes (metadata) of the data array
        attrs = { 'title' : '3-week sea ice area timeseries',
                  'units' : 'km2'}
        
        # ds_out = xr.DataArray(data = np.array(si_series),
        #                 dims = dims,
        #                 coords = coords,
        #                 attrs = attrs)
        
        ds_out = xr.Dataset(
                {"seaice": (["nstorms", "time"], np.array(si_series)),
                 "si_clim": (["nstorms", "time"], np.array(si_clim_series))},
                coords={"nstorms": (["nstorms"], np.arange(0,storm_count)),
                        'time' :   (["time"], np.arange(-7,14+1)),
                        },
                attrs = attrs
            )
        
        ds_out.to_netcdf(savepath+str(year)+"_seaice.nc")
    
#%% repeat for *area*

if run['area']:
    
    ncarea_path = root_path + 'area/'
    
    savepath_area = '/Users/mundi/Desktop/month-hemi/nh_data/area/'
    if not os.path.exists(savepath_area):
        os.makedirs(savepath_area)
    
    for year in years:
        box_areas = []
        ice_areas_15 = []
        ice_areas_80 = []
        
        storm_count = 0
        load_summer=True
        
        # open census and sea ice
        [startdate, enddate] = readCensus(census_path0+'census_'+str(year)+'.csv', convertDT=True)[0]
        
        
        # open ice area
        ds_area = xr.open_dataset(ncarea_path + str(year) +'_area.nc')
        ice_sorter80 = ds_area['ice_area80'].values
        ice_sorter15 = ds_area['ice_area15'].values
        box_area = ds_area['box_area'].values
        ds_area.close()
        
        # storm loop
        for storm_num, sd in enumerate(startdate):
            if sd.month not in summer_months:
                box_areas.append(box_area[storm_num])
                ice_areas_15.append(ice_sorter15[storm_num])
                ice_areas_80.append(ice_sorter80[storm_num])
                storm_count += 1
            else:
                if load_summer:
                    
                    census_path1 = '/Users/mundi/Desktop/teststuff/original_census/'
                    census_file1 = census_path1+'census_'+str(year)+'.csv'
                    [startdate1, enddate1] = readCensus(census_file1, convertDT=True)[0]
                    
                    try:
                       ### SUMMER
                       if year >= 2000:
                           ncpath_area1 = '/Users/mundi/Desktop/teststuff/modern_era/'
                       elif year < 2000:
                           ncpath_area1 = '/Users/mundi/Desktop/teststuff/hist_era/'
                       
                       # open ice area
                       ncpath_area1 = ncpath_area1+str(year)+'/'
                       ds_area1 = xr.open_dataset(ncpath_area1 + str(year) +'_area.nc')
                       ice_sorter1 = ds_area1['ice_area80'].values
                       ice_sorter2 = ds_area1['ice_area15'].values
                       box_area1 = ds_area1['box_area'].values
                       ds_area1.close()
                       
                    except:
                        print('- skip: '+'/Users/mundi/Desktop/FINAL/'+ str(year) +'_seaice.nc')
                        continue
                    
                    for storm_num1, sd1 in enumerate(startdate1):
                        if sd.month in summer_months:
                            box_areas.append(box_area1[storm_num1])
                            ice_areas_15.append(ice_sorter2[storm_num1])
                            ice_areas_80.append(ice_sorter1[storm_num1])
                            storm_count += 1
                    
                    load_summer=False
            
        #### save new file
        # names of the dimensions in the required order
        dims = ('nstorms')
        # create coordinates to use for indexing along each dimension
        coords = {'nstorms' : np.arange(0,storm_count)}
        # attributes (metadata) of the data array
        attrs = { 'title' : 'area fractions for each storm',
                  'units' : ''}
        
        ds_out = xr.Dataset(
                {
                    "box_area": (["nstorms"], np.array(box_areas)),
                    "ice_area15":  (["nstorms"], np.array(ice_areas_15)),
                    "ice_area80":  (["nstorms"], np.array(ice_areas_80))
                },
                coords={"nstorms": (["nstorms"], np.arange(0,storm_count))},
                attrs = attrs
            )
        
        
        ds_out.to_netcdf(savepath_area+str(year)+"_area.nc")

#%% end
