#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on July 29 2022
Last run 9/15/2022
Update new years 6/2/2023
cyclone_tracker4sens.py

Oct 2025
--> sensitivity to pressure threshold...


@author: mundi
"""
import numpy as np
#%% inputs

root1 = '/home/mundi/month-hemi/nh_data/'
# root1 = '/Users/mundi/Desktop/month-hemi/nh_data/'

ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/daily/'    
# ice_fname = '/Users/mundi/Desktop/seaice/'

census_save_path = root1+'sensitivity/'
nc_save_path = census_save_path+'contours/'

census_name = 'census_' 
contour_name = '_contours.nc'

# years = np.arange(2010,2020)
# years = np.arange(2010,2015)
# years = np.arange(2015,2020)

# years = [2012, 2013, 2014, 2016]
# years = [2017, 2018, 2019]
years = [2016]

# months = np.arange(1,12+1)
# months = [1,2]
months = [1,9]

#%% imports
if True:
    from datetime import datetime, timedelta
    import calendar
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import gc
    
import myfuncs as mf

import time as timeIN
TIMESTART = timeIN.time()

#%% fxns

def export2nc(nc_in, all_contours, all_cont_dts):
    import netCDF4
    import traceback
    
    ncfile = netCDF4.Dataset(nc_in,mode='w',format='NETCDF4_CLASSIC')
    try:
        ncfile.title='990, 1000 mb pressure over storm days'
        
        cp = ['990','1000']
        ncfile.createDimension('num_coords', 2)
        prev_stormstr=''
        num=0
        
        for c, cont in enumerate(all_contours):
            stormstr = all_cont_dts[c].strftime('%Y_%m%d')
            
            if c == len(all_cont_dts)-1:
                if all_cont_dts[c] == all_cont_dts[c-1]:
                    pstr='1000'
                else: pstr='990'
            elif (all_cont_dts[c+1] != all_cont_dts[c]):
                pstr = '1000'
            else:
                pstr = '990'
                
            if stormstr+'_'+pstr == prev_stormstr+'_'+pstr:
                num+=1
            else:
                num=0
            
            ## create dimensions
            ncfile.createDimension('points'+str(c+1), len(cont))
            
            try:
                pres_edges = ncfile.createVariable('coords_' +stormstr+'_'+pstr+'_'+str(num),
                                               np.float64, ('points'+str(c+1),'num_coords'))
            except RuntimeError as rte:
                print(rte)
                print(traceback.format_exc())
                print('~~~  '+'coords_' +stormstr+'_'+pstr+'_'+str(num)+ '  ~~~')
            
            prev_stormstr = stormstr
                
            # print('coords_' +stormstr+'_'+cp[c%2])
            pres_edges.long_name = 'Coordinates for ' + cp[c%2] +'mb on ' + stormstr
            pres_edges[:,:] = cont # load data to nc
            
            
    except Exception as e:
        print(''); print('* ERROR: exporting to nc')
        print(e)
        print(traceback.format_exc())
    finally:
        ncfile.close()
        print("*exported*")


def final_storm_check(starting_list, all_storms):
    
    starting_array = np.array(starting_list)
    # check final pressure of previous grouping
    pressures = starting_array[:,2]
    # and that storm duration is long enough
    dates = starting_array[:,-1]
    duration = dates[-1] - dates[0]
    # and close enough to ice to be considered
    storm_lat = np.mean(starting_array[:,1]) # average storm latitude
    ice_distance = mf.get_storm_distance(ice_fname, dates[0], storm_lat, bbox_edges=[])
    
    if (duration < min_duration) or (np.min(pressures) > min_pres) or (ice_distance > icedist_thresh):
        starting_list = [stormtime] # throw out, reset list
        return [stormtime], all_storms #continue
    # append old list; make new
    all_storms.append(starting_list)
    starting_list = [stormtime]
    
    return starting_list, all_storms

#%%% era data

def LZ(day):
    if day>=10:
        return str(day)
    elif day<10:
        return '0'+str(day)

def era5_hourly(year, month, days, daily=False, return_ds=False):
    '''
    Parameters
    ----------
    year : list
       ex. ['2010']
    months, days: lists

    Returns
    -------
    ds: xarray dataset with variable
    '''
    if type(year) != list:
        year = [str(year)]
    if type(month) != list:
        if month < 10: monthstr = '0'+str(int(month))
        else: monthstr = str(int(month))
        month = [monthstr]
    if type(days) != list and type(days[0])!=str:
        days = ['0'+str(int(day)) if day<10 else str(int(day)) for day in days]

    import cdsapi
    import xarray as xr
    import io
    from urllib.request import urlopen# start the client
    cds = cdsapi.Client()# dataset you want to read
    dataset = "reanalysis-era5-single-levels"# flag to download data
    params = {
        "data_format": "netcdf",
        "download_format": "unarchived",
        "product_type": ["reanalysis"],
        "variable": ['mean_sea_level_pressure'],
        'year':year,
        'month':month,
        'day':days,
        "time": ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',
                 '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00',
                 '16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'],
        "grid":[0.25,0.25],
        "area":[90, -180, 60, 180]
        }
    # retrieves the path to the file
    fl = cds.retrieve(dataset, params)
    # load into memory
    with urlopen(fl.location) as f:
        bytes_ = f.read()
        ds = xr.open_dataset(io.BytesIO(bytes_), decode_times=True)
        
    if daily: ds = ds.resample(valid_time='1D').mean() # daily
    ds['msl'] = ds['msl']/100
    
    if return_ds: return ds
    
    time = [str(dt) for dt in ds['valid_time'].values]
    TIME = []
    for dt in time:
        yy = int(dt[0:4])
        mm = int(dt[5:7])
        dd = int(dt[8:10])
        hh = int(dt.split('T')[1][0:2])
        TIME.append(datetime(yy,mm,dd,hh))
    
    lon, lat = np.meshgrid(ds['longitude'].values, ds['latitude'].values)
    MSL = ds['msl'].values
    
    return np.array(TIME), lon, lat, MSL

#%%% contours
from cartopy.util import add_cyclic_point

def get_contour_points2(x,y,pressure, year, month, day, level, storm_info, lon_thresh=10, lat_thresh=5):
    if type(level) != list: level=[level]
    
    # print(''); print(year, month, day); print('---------')
    # if type(year) != list: year=[str(year)]
    # if type(month)!= list: month=[str(month)]
    # if type(day) != list: day = [str(day)]
    
    xv = x[0,:]
    yv = y[:,0]
    
    cyclic_data, cyclic_lons = add_cyclic_point(pressure, coord=xv)
    x, y = np.meshgrid(xv,yv)
    # x = np.where(x<0, x+360, x) # make all lons positive for comparison
    
    ### find minimum pressures and compare to storm_info ###!!!!
    storm_x = storm_info[0]
    storm_y = storm_info[1]
    nearby_min = False
    
    while_count= 0
    while nearby_min==False:
        while_count+=1
        # print('wc', while_count, np.nanmin(pressure))
        
        indp = np.where(pressure == np.nanmin(pressure))
        minlon = x[indp]
        minlat = y[indp]
        
        # if while_count%50==0:
        #     print('wc', while_count, np.nanmin(pressure), indp)
        #     print('>', storm_x, storm_y, minlon, minlat)
        
        for pt_id, lo in enumerate(minlon):
            
            if minlon[pt_id]<0: minlon[pt_id] =minlon[pt_id]+360
            if abs(minlon[pt_id] - storm_x) < lon_thresh and abs(minlat[pt_id] - storm_y) < lat_thresh:
                min_x = float(minlon[0])
                min_y = float(minlat[0])
                min_p = np.nanmin(pressure)
                nearby_min = True
                break
            else:
                pressure[indp] = np.nan #[pt_id]
                break
            
        if while_count>121*1440: raise ValueError("'while' loop, get_contours")
        
    fig, ax = mf.background_plot(returnfig=True)
    ax.set_title(str(year)+' '+str(month)+' '+str(day) + ': ' +str(round(min_p,1))+' hPa' +',  '+str(round(min_x,2)), 
                 fontsize=26)
    
    contours_all = ax.contour(x, y, pressure, levels=level,
                      linewidths = 1.5, zorder=10, colors='k',
                      transform=ccrs.PlateCarree())
    
    
    ax.plot(min_x, min_y,'y*',transform=ccrs.PlateCarree(), markersize=40)    
    
    coords=[]
    for lvl, lev in enumerate(level):
        try:
            coords = coords + contours_all.allsegs[lvl]
        except IndexError:
            print(str(lev) +': ' + str(month)+', '+ str(day))
        
    min_x = np.squeeze(min_x)
    if min_x > 180:
        min_x = min_x - 360
    
    plt.close(fig)
    return coords, min_x, np.squeeze(min_y), ax

def detect_contours2(ds, storm_daymonth_grid, year, storm_info, intervals=[990,1000]):
    import matplotlib.path as mpltPath
    
    # return lists of selcted contours and their respective date
    all_cont_dt1, all_conts1 = [], []
    
    snum=-1
    for month, day in storm_daymonth_grid:
        snum+=1
        ### get pressure contours      
        intervals = [990,1000] #[980,985,990,995,1000,1005,1010] #
        current_date = datetime(int(year), int(month), int(day))
        
        if snum >= len(storm_info):break
        
        try:
            ds = ds.sel(valid_time=datetime(year, month, day))
            x, y = np.meshgrid(ds['longitude'].values, ds['latitude'].values)
            pressure = np.squeeze(ds['msl'].values)
            time=ds.valid_time.values
        except:
            time, x,y, pressure = era5_hourly(year, month, [LZ(day)], daily=True)
            # pressure=pressure[0]
            pressure = np.squeeze(pressure)
            
        print(month, day, time)
        
        coords, min_x, min_y, ax = get_contour_points2(x,y,pressure, year, month, day, intervals, 
                                                       storm_info[snum])
        # current_title = ax.get_title()
        
        points = np.array([ [min_x, min_y] ])
        for cc1 in coords :
            
            cc = np.array(cc1)
            # clon = cc2[:,0]; clat=cc2[:,1]
            # clon= np.where(clon>180, clon-360, clon)
            # cc = np.vstack((clon, clat)).T
            
            if len(cc) > 2 : path = mpltPath.Path(cc)
            else: continue
            pts_inside = path.contains_points(points, radius=15)
            
            # ax.set_title(current_title +'....' + str(min_x), fontsize=28)
            
            if any(pts_inside):
                # if abs( abs(np.nanmean(cc[:,0])) - abs(min_x) ) < 90:
                all_conts1.append(cc)
                all_cont_dt1.append(current_date)
                ax.plot(cc[:,0], cc[:,1], 'r', linewidth=3, transform=ccrs.PlateCarree())
                continue
            
            if abs(abs(min_x) - 180) < 20:
                second_x = 175
                second_pts = path.contains_points(np.array([ [second_x, min_y] ]))
                ax.plot(second_x, min_y,'c*',transform=ccrs.PlateCarree(), markersize=30) 
                if any(second_pts):
                    all_conts1.append(cc)
                    all_cont_dt1.append(current_date)
                    [x1,x2], [y1,y2] = mf.geoplot_2d(cc[:,0], cc[:,1])
                    ax.plot(x1, y1, 'b', linewidth=5, transform=ccrs.PlateCarree())
                    ax.plot(x2, y2, 'b', linewidth=5, transform=ccrs.PlateCarree())
                continue
    
    print()
    return all_cont_dt1, all_conts1   
    
#%% pressure threshold for storm lifetime
# weaker min pressure at storm start and storm final
pres_thresh = 6

min_pres = 984 + 10

### how far apart (degrees) two points below pressure threshold need to be to count as different storms
## used for comparing many points below threshold for single slp profile
grouping_thresh = 15

### how far apart (degrees) two absolute minima need to be to count as different storms
## used for comparing single points from different slp profiles
sorting_thresh = 50

### time between storms
timediff_thresh = timedelta(hours=6)

### minimum storm duration
min_duration = timedelta(days=2)

### distance between storm min and ice edge (km)
icedist_thresh = 1250 

#%% find locations of minimum pressures
# returns all_storms: lon/lat/minpressure each day

for year in years:
    year_start_time = timeIN.time()

    storm_counter = 1
    storm_sorter = {}
    sort_count = 0
    
    for month in months: 
        plt.close('all')
        day_counter = -1
        my_storms_daily = []
        # [14,15,16,17,18,19,20,21,22,23,24,25]: #range(20): #range(calendar.monthrange(year, month)[1])
        
        # load monthly data from server
        calendar_days = list(np.arange(1,calendar.monthrange(year, month)[1]+1))
        alldays = [str(d) if d>=10 else '0'+str(d) for d in calendar_days]
        TIME, lon, lat, MSL = era5_hourly([str(year)], [str(month)], alldays)
        
        for day in np.arange(1,calendar.monthrange(year, month)[1]+1): 
            day_counter += 1    
            current_date = datetime(year, month, day)
            datestr = current_date.strftime('%Y%m%d')
    
            # get daily value
            t_indices = []
            for tind, dt in enumerate(TIME):
                if dt.day == day:
                    t_indices.append(tind)
            msl = np.squeeze( MSL[t_indices, :, :] )
            time = np.squeeze( TIME[t_indices] )

            ### get storm locations
            hour_counter = -1
            my_storms_hourly = {}
            for tind, slp in enumerate(msl): # loop through each hour
                hour_counter += 1 
                my_storms_hourly[hour_counter] = []
                # find indices of points of low pressure
                indp = np.where(slp < min_pres + pres_thresh) 
                if slp[indp].size == 0: 
                    continue

                # sort by longitudes to group different storm events
                sorted_lons, sorted_lats, sorted_slp = zip(*sorted(zip(lon[indp], lat[indp], slp[indp])))
                sorted_lons = np.array(sorted_lons)
                sorted_lons = np.where(sorted_lons<0, sorted_lons+360, sorted_lons)
                
                ### group detected low pressure points by location
                ### find abs minimum for each pressure grouping
                # --- check lon for comparison
                for lx, ll in enumerate(sorted_lons):
                    if ll> 180: 
                        l1 = ll-360
                    else: 
                        l1=ll
                    if sorted_lons[lx-1] > 180:
                        l2 = sorted_lons[lx-1] -30
                    else:
                        l2 = sorted_lons[lx-1]
                    # --------------------------
                    if lx==0: 
                        starting_list = [[ll, sorted_lats[lx], sorted_slp[lx]]]
                        continue
                    if abs(l1-l2) <= grouping_thresh: # same storm grouping if pressure mins are within this many degrees
                        starting_list.append([ll, sorted_lats[lx], sorted_slp[lx]])
 
                    if abs(l1-l2) > grouping_thresh or lx == len(sorted_lons)-1: 
                        # if greater than that difference, get absolute min of previous grouping then make new group
                        # OR get minimum of last grouping
                        starting_array = np.array(starting_list)
                        maxslp, maxlon, maxlat = zip(*sorted(zip(starting_array[:,2], starting_array[:,0], starting_array[:,1])))
                        my_storms_hourly[hour_counter].append([maxlon[0],maxlat[0],maxslp[0], time[tind]])
                        # start new grouping
                        starting_list = [[ll, sorted_lats[lx], sorted_slp[lx]]] 
                starting_list= None
            # append mins detected each hour
            my_storms_daily.append(my_storms_hourly)
           
        ### sort storm locations for the day
        sortnum=0
        for hourlyinfo in my_storms_daily:
            for hr in hourlyinfo:
                for si, storm_list in enumerate(hourlyinfo[hr]):
                    if sort_count != 0:
                        sortnum = 0
                        while sortnum < sort_count:
                            
                            if storm_list[0] > 180:
                                l1 = storm_list[0] - 360
                            else:
                                l1 = storm_list[0]
                            if storm_sorter[sortnum][-1][0] > 180:
                                l2 = storm_sorter[sortnum][-1][0] - 360
                            else:
                                l2 = storm_sorter[sortnum][-1][0]
                            
                            if abs(l1-l2) <= sorting_thresh:
                                storm_sorter[sortnum].append(storm_list)
                                added=True
                                break
                            else:
                                sortnum+=1
                                if sortnum >= sort_count:
                                    added=False
                                    break
                        if not added:
                            storm_sorter[sort_count] = [storm_list]
                            added = False
                            sort_count += 1
                    elif sort_count == 0:
                        storm_sorter[sort_count] = [storm_list]
                        sort_count += 1

### storm_sorter --> all_storms

#%% chunk storms and organize
    all_storms = []
    for stormkey in storm_sorter:
        starting_list = []
        # loop through each individual point identified in each preliminary grouping
        for sidx, stormtime in enumerate(storm_sorter[stormkey]):
            if sidx == 0: 
                starting_list = [stormtime]
                continue
            # compare timestamp with previous timestamp (identify same storm grouping based on time)
            # if abs(stormtime[-1] - storm_sorter[stormkey][sidx-1][-1])<timediff_thresh:
            if abs(stormtime[-1] - starting_list[-1][-1])<timediff_thresh:
                starting_list.append(stormtime)
                
                # if final
                if sidx==len(storm_sorter[stormkey])-1:
                   starting_list, all_storms = final_storm_check(starting_list, all_storms)
                   if len(starting_list) == 1: continue
    
            else: # split storm grouping
                if len(starting_list) == 0: 
                    continue
                starting_list, all_storms = final_storm_check(starting_list, all_storms)
                if len(starting_list) == 1: continue
  #%% check for duplicate times - remove adjacent mins within the same storm
    all_storms_backup = all_storms.copy()
    all_storms1 = []
    stormnum=0
    for storm in all_storms:
        stormnum+=1
        checked_storm=[]
        first_time=True
        still_duplicates=True
        while still_duplicates:
            still_duplicates=False
            if first_time:
                working_storm = storm
                first_time=False
            else:
                working_storm = checked_storm
                checked_storm=[]
            # psuedo for loop
            ip = 0 
            while ip < len(working_storm)-1:
                if ip==0:
                    prev_date = working_storm[ip][-1]
                    ip+=1
                if working_storm[ip][-1] == prev_date:
                    still_duplicates=True
                    dupes=[working_storm[ip-1]]
                    # print('dupe found!')
                    # get true min
                    while working_storm[ip][-1] == prev_date:
                        # print('dupppees')
                        dupes.append(working_storm[ip])
                        ip+=1
                        if ip == len(working_storm)-1: break
                    dupes = np.array(dupes)
                    # sort by pressure and select min
                    dupes = dupes[dupes[:, 2].argsort()]
                    checked_storm.append(list(dupes[0]))
                else:
                    checked_storm.append(working_storm[ip-1])
                prev_date = working_storm[ip][-1]
                ip+=1
           
        working_storm = checked_storm
        all_storms1.append(working_storm)
        
    all_storms = all_storms1.copy()
    all_storms = sorted(all_storms, key=lambda t: t[0][-1]) # sort by starting date


#%% get census info for export

    startdates, enddates = [], []
    censusinfo = []
    storm_info = []
    for storm_event in all_storms:
        event = np.array(storm_event)
        startdates.append(event[0][-1])
        enddates.append(event[-1][-1])
        
        # find minimum pressure
        pressures = event[:,2]
        minimum = np.min(pressures)
        
        #start_lat, start_lon, minimum, finish_lat, finish_lon
        censusinfo.append([ event[0][1], event[0][0], minimum, event[-1][1], event[-1][0] ])
        
        
        info1 = [] # all hours of the day, to be averaged
        info2 = [] # average of each day for every storm 
        for tidx, timing in enumerate(event):
            if tidx==0: 
                info1 = [timing[0:3]]
                continue
            if timing[-1].day == event[tidx-1][-1].day:
                info1.append(timing[0:3])
                if tidx == len(event)-1:
                    info2.append(np.mean(info1,axis=0))
            else:
                info2.append(np.mean(info1,axis=0))
                info1 = [timing[0:3]]
        storm_info.append(info2)
        info2 = []


    #%% begin storm_id
    timing_grid = []
    for xx in range(0,len(startdates)):
        timing_grid.append((startdates[xx], enddates[xx]))
    
    storm_ranges = []
    analysis_ranges = []
    for startdt, enddt in timing_grid:
        week_ago = startdt - timedelta(days=7)
        two_week = startdt + timedelta(days=14) ###!!! startdt vs end (same duration)
        analysis_ranges.append(mf.daterange(week_ago, two_week, dt=24))
        storm_ranges.append(mf.daterange(startdt, enddt, dt=24))
    

    #%%% * START STORM LOOP * get contours
    print(' '); print(' ')
    stormstr_prev = ''
    
    # starting datatset
    # ds = era5_hourly([str(year-1),str(year)], ['12','01'], 
    #                  [LZ(d+1) for d in np.arange(0,31)], daily=True, return_ds=True)
    
    
    for storm_num, storm_event in enumerate(storm_ranges):         ###!!!!!!!
        print('Storm ' + str(storm_num+1) + '/' + str(len(storm_ranges)))
        gc.collect()
        
        ### get needed contours and create box
        stormstr1 = storm_event[0].strftime('%Y_%m%d')
        # daymonth_grid = [(dt.month, dt.day) for dt in analysis_ranges[storm_num]]
            
        storm_daymonth_grid = []
        for storm_day in storm_event:
            storm_daymonth_grid.append((storm_day.month, storm_day.day))
            
        month_list = np.unique([LZ(dt.month) for dt in analysis_ranges[storm_num]])
        day_list = np.unique([LZ(dt.day) for dt in analysis_ranges[storm_num]])
        ds = era5_hourly(year, list(month_list), list(day_list), daily=True, return_ds=True)
        
        try:
            ds2 = ds.sel(valid_time=slice(analysis_ranges[storm_num][0], analysis_ranges[storm_num][-1]))
        except:
            month_list = np.unique([LZ(dt.month) for dt in analysis_ranges[storm_num]])
            ds = era5_hourly([str(year)], list(month_list), 
                             [LZ(d+1) for d in np.arange(0,31)], daily=True, return_ds=True)
            try:
                ds2 = ds.sel(valid_time=slice(analysis_ranges[storm_num][0], analysis_ranges[storm_num][-1]))
            except:
                month_list = np.unique([LZ(dt.month) for dt in analysis_ranges[storm_num]])
                day_list = np.unique([LZ(dt.day) for dt in analysis_ranges[storm_num]])
                ds = era5_hourly(year, list(month_list), list(day_list), daily=True, return_ds=True)
                
                ds2 = ds.sel(valid_time=slice(analysis_ranges[storm_num][0], analysis_ranges[storm_num][-1]))
            
    
        all_cont_dts, all_contours = detect_contours2(ds2, storm_daymonth_grid, year, storm_info=storm_info[storm_num]) 
    
        if stormstr1==stormstr_prev:
            stormstr = stormstr1 + str(storm_num)
        else:
            stormstr=stormstr1
        stormstr_prev = stormstr
    
        ### export *.nc!
        nc_name = stormstr+contour_name
        nc_in = nc_save_path+ nc_name
        
        export2nc(nc_in, all_contours, all_cont_dts)


    try: ds.close() 
    except: pass
    try: ds2.close()
    except: pass
    gc.collect()
    
#%% export storm_info to new census database
    import csv
    header = ['start', 'start_lat', 'start_lon', 'minimum', 'finish', 'finish_lat', 'finish_lon']
    
    data = []
    for storm_num, storm_row in enumerate(storm_info):
        
        storm_data = [storm_ranges[storm_num][0]] # start date
        
        storm_data.append(round(storm_row[0][1],4)) # start lat
        storm_data.append(round(storm_row[0][0],4)) # start lon
        
        # find storm minimum pressure
        min_pres = 9999
        for day_entry in storm_row:
            if day_entry[-1] < min_pres:
                min_pres = day_entry[-1]
        storm_data.append(round(min_pres,3)) 
        
        storm_data.append(storm_ranges[storm_num][-1]) # finish date
        storm_data.append(round(storm_row[-1][1],4)) # finish lat
        storm_data.append(round(storm_row[-1][0],4)) # finsih lon
        
        # add data to full list for export
        data.append(storm_data)
    
    
    ### create csv and export
    with open(census_save_path+census_name+str(year)+'.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        
        # write the header
        writer.writerow(header)
    
        # write multiple rows
        writer.writerows(data)
    
#%% pickle some variables for later
    print('')
    print('FINISHED STORM DATA COLLECTION: '+str(year))
    print((timeIN.time() - year_start_time)/60, ' minutes')
    print('')
    
    ### save important stuff for later
    import pickle
    
    output = open(census_save_path + str(year)+'_all_storms.pkl', 'wb')
    pickle.dump(all_storms, output)
    output.close()
    
#%% end
print(' '); print(' ')
print('------------------')
print('       done!      ')
print('------------------')
print('------------------')
print('elapsed time: ')
print((timeIN.time() - TIMESTART)/60)