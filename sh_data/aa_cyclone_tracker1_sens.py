#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nov 15 2023 | Oct 2025
aa_cyclone_tracker1_sens.py

--> sensitivity to pressure threshold...

@author: mundi
"""
#%% inputs
from datetime import datetime, timedelta
import numpy as np
import gc

# years = np.arange(2010,2019+1)
# years = np.arange(2010, 2015)
years = np.arange(2015,2020)

months = np.arange(1,12+1)


# nc_save_path = '/Users/mundi/Desktop/antarctica/test_957_all/' #'_any2/'
# nc_save_path = '/Users/mundi/Desktop/antarctica/allmonths-test/'
nc_save_path = '/home/mundi/month-hemi/sh_data/sensitivity/'

# ice_fname = '/Users/mundi/Desktop/seaice/south/'
ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/south/'

census_name = 'census_' 

#%% pressure threshold for storm lifetime
# weaker min pressure at storm start and storm final (984)
pres_thresh = 6

min_pres = 957 + 10

### how far apart (degrees) two points below pressure threshold need to be to count as different storms
## used for comparing many points below threshold for single slp profile (15)
grouping_thresh = 15 #10

### how far apart (degrees) two absolute minima need to be to count as different storms
## used for comparing single points from different slp profiles (50)
sorting_thresh = 50 #30 #40

### time between storms
timediff_thresh = timedelta(hours=6)

### minimum storm duration
min_duration = timedelta(days=2)

### distance between storm min and ice edge (km)
icedist_thresh = 1250 

### ice area percentages
ice_lims = [0,100] #[20,80] #[10,90] #[0,100] #[20,80] #[20,80]
ice_criteria = [15, 'cell'] #[80,'ice'] # #15/80, 'cell'/'ice'

miz= [0.15,0.80] #[0,0.58]

### 
contour_levels = [970, 980]# [990, 1000]

if (ice_criteria[0] not in [15,80]) or (ice_criteria[1] not in ['ice', 'cell']):
    raise ValueError('ice_criteria = '+str(ice_criteria))


text = []
text.append('\nSTART: ' + str(years) + ', ' + str(months))
text.append('pressure: ' + str(min_pres) +' +- '+ str(pres_thresh))
text.append('grouping = '+str(grouping_thresh)+' | ' +'sorting = '+ str(sorting_thresh))
text.append('timing: '+str(timediff_thresh)+' hrs')
text.append('duration: '+str(min_duration))
text.append('ice distance: '+ str(icedist_thresh)+ ' km')
text.append('ice lims: '+ str(ice_lims)+ ' - '+str(ice_criteria[0])+'% SIC '+str(ice_criteria[1])+' area')
text.append('miz: '+str(miz))
text.append('contour levels: '+str(contour_levels) + ' (mb, for determining storm area)')
text.append(' ')
text.append('out: '+ nc_save_path)
text.append('census: '+ census_name)
text.append('\n\n')
text.append('*\n')


with open(nc_save_path+'readme.txt', 'w') as f:
    for txt in text:
        print(txt)
        f.write(txt)
        f.write('\n')

#%% imports
import calendar
import matplotlib.pyplot as plt
import xarray as xr
from cartopy.util import add_cyclic_point
import cartopy.crs as ccrs
    
import myfuncs_aa as mf

import time as timeIN
TIMESTART = timeIN.time()

#%% fxns
import sys, os
class HidePrint:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def export2nc(nc_in, all_contours, all_cont_dts, intervals = [990, 1000]):
    import netCDF4
    import traceback
    
    ncfile = netCDF4.Dataset(nc_in,mode='w',format='NETCDF4_CLASSIC')
    try:
        titlestr=''
        for interval in intervals:
            titlestr += str(interval)
            titlestr += ', '
        ncfile.title= titlestr[:-2] + ' mb pressure over storm days'
        
        ncfile.createDimension('num_coords', 2)
        prev_stormstr=''
        num=0
        
        for c, cont in enumerate(all_contours):
            stormstr = all_cont_dts[c].strftime('%Y_%m%d')
            
            if c == len(all_cont_dts)-1:
                if all_cont_dts[c] == all_cont_dts[c-1]:
                    pstr=str(intervals[1])
                else: pstr=str(intervals[0])
            elif (all_cont_dts[c+1] != all_cont_dts[c]):
                pstr = str(intervals[1])
            else:
                pstr = str(intervals[0])
                
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
                
            print('coords_' +stormstr+'_'+str(intervals[c%2]) )
            pres_edges.long_name = 'Coordinates for ' + str(intervals[c%2]) +'mb on ' + stormstr
            pres_edges[:,:] = cont # load data to nc
            
            
    except Exception as e:
        print(''); print('* ERROR: exporting to nc')
        print(e)
        print(traceback.format_exc())
    finally:
        ncfile.close()
        print("*exported*")


def get_storm_area(all_contours, storm_event, ice_fname, ice_criteria, miz):
    ### get bbox_edges
    with HidePrint():
        bbox_edges = mf.get_bbox_edges(all_contours) 
        
    out = {}
        
    ### area of bbox
    # get si lon/lat grid
    si_grid, si_lon, si_lat = mf.load_seaice(ice_fname, year, storm_event[0].month, storm_event[0].day, latlon=True)
    # get box mask
    in_mask = mf.find_points_in_contour(bbox_edges, si_lon, si_lat)
    # calculate area
    si_in = np.ma.masked_array(si_grid, mask=in_mask).filled(np.nan)
    box_area = mf.get_total_area(si_in)
    
    
    ### get ice areas (>15%, >80%)
    if ice_criteria[1] == 'ice':
        ice_area15 = mf.calc_seaice_area(np.where(si_in<miz[0], 0, si_in)) 
        ice_area80 = mf.calc_seaice_area(np.where(si_in<miz[1], 0, si_in))
    elif ice_criteria[1] == 'cell':
        ice_area15 = mf.get_total_area(np.where(si_in<miz[0], np.nan, si_in)) 
        ice_area80 = mf.get_total_area(np.where(si_in<miz[1], np.nan, si_in))
    else: 
        print('> get_storm_area(): ice_criteria ['+str(ice_criteria[1])+']')
        raise(ValueError)
    
    ### export info
    out['box_area'] = box_area
    out['ice_area15'] = ice_area15
    out['ice_area80'] = ice_area80
    
    ### calculate icearea fraction
    ice_frac = ice_area80*100/box_area

    return ice_frac, out, bbox_edges

def diagnostic_print(dates, duration, pressures, ice_distance):
    print('Removing Storm: '+ str(dates[0])+'-'+str(dates[-1]))
    if (duration < min_duration): print('Duration: '+str(duration))
    if (np.min(pressures) > min_pres): print('Pressure: '+str(np.min(pressures)))
    if (ice_distance > icedist_thresh): print('Ice Distance: '+str(ice_distance))
    print('-')

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
        diagnostic_print(dates, duration, pressures, ice_distance)
        starting_list = [stormtime] # throw out, reset list
        return [stormtime], all_storms #continue
    # append old list; make new
    all_storms.append(starting_list)
    starting_list = [stormtime]
    
    return starting_list, all_storms

#%%% era data

def era5_hourly(year, month, days, daily = False):
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
        "area":[-55, -180, -90, 180]
        }
    # retrieves the path to the file
    fl = cds.retrieve(dataset, params)
    # load into memory
    with urlopen(fl.location) as f:
        bytes_ = f.read()
        ds = xr.open_dataset(io.BytesIO(bytes_), decode_times=True)
        
    if daily: ds = ds.resample(valid_time='1D').mean() 
    
    time = [str(dt) for dt in ds['valid_time'].values]
    TIME = []
    for dt in time:
        yy = int(dt[0:4])
        mm = int(dt[5:7])
        dd = int(dt[8:10])
        hh = int(dt.split('T')[1][0:2])
        TIME.append(datetime(yy,mm,dd,hh))
    
    lon, lat = np.meshgrid(ds['longitude'].values, ds['latitude'].values)
    MSL = ds['msl'].values/100
    
    return np.array(TIME), lon, lat, MSL

#%%% contours

def get_contour_points2(year, month, day, level, storm_info, lon_thresh=10, lat_thresh=5):
    if type(level) != list: level=[level]
    
    print(''); print(year, month, day); print('---------')
    
    if type(year) != list: year=[str(year)]
    if type(month)!= list: month=[str(month)]
    if type(day) != list: day = [str(day)]
    time, x,y, pressure = era5_hourly(year, month, day, daily=True)
    pressure = np.squeeze(pressure)
            
    xv = x[0,:]
    yv = y[:,0]
    
    cyclic_data, cyclic_lons = add_cyclic_point(pressure, coord=xv)
    x, y = np.meshgrid(xv,yv)
    
    ### find minimum pressures and compare to storm_info ###!!!!
    storm_x = storm_info[0]
    storm_y = storm_info[1]
    nearby_min = False
    
    while nearby_min==False:
        indp = np.where(pressure == np.nanmin(pressure))
        minlon = x[indp]
        minlat = y[indp]
        
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
    fig, ax =  mf.background_plot(returnfig=True)
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


def detect_contours2(storm_daymonth_grid, daymonth_grid, year, storm_info, intervals=[990,1000]):
    import matplotlib.path as mpltPath
    
    # return lists of selcted contours and their respective date
    all_cont_dt1, all_conts1 = [], []
    
    snum=-1
    for month, day in storm_daymonth_grid:
        snum+=1
        ### get pressure contours      
        #intervals #[980,985,990,995,1000,1005,1010] #
        current_date = datetime(int(year), int(month), int(day))
        
        if snum >= len(storm_info):break
        
        coords, min_x, min_y, ax = get_contour_points2(year, month, day, intervals, 
                                                       storm_info[snum])
        
        points = np.array([ [min_x, min_y] ])
        for cc1 in coords :
            cc = np.array(cc1)
            
            if len(cc) > 2 : path = mpltPath.Path(cc)
            else: continue
            pts_inside = path.contains_points(points, radius=15)
            
            if any(pts_inside):
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
    
#%% find locations of minimum pressures
# returns all_storms: lon/lat/minpressure each day
print('--> daily min pressure locations')

for year in years:
    
    storm_counter = 1

    storm_sorter = {}
    sort_count = 0
    
    for month in months: 
        gc.collect()
        plt.close('all')
        day_counter = -1
        my_storms_daily = []
        
        print('*** '+str(month)+'-'+str(year)+' ***')
        month_time_start = timeIN.time()
        
        
        # load monthly data from server 
        calendar_days = list(np.arange(1,calendar.monthrange(year, month)[1]+1))
        # calendar_days = np.arange(1, 29, 1)
        alldays = [str(d) if d>=10 else '0'+str(d) for d in calendar_days]
        
        #### data(!)
        TIME, lon, lat, MSL = era5_hourly([str(year)], [str(month)], alldays)
            
        
        for day in calendar_days:
            day_counter += 1    
            current_date = datetime(year, month, day)
            datestr = current_date.strftime('%Y%m%d')
    
            # get daily value
            t_indices = []
            for tind, dt in enumerate(TIME):
                if dt.day == day and dt.month==month:
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
                    if ll>= 180: 
                        l1 = ll-360
                    else: 
                        l1=ll
                    if sorted_lons[lx-1] >= 180:
                        l2 = sorted_lons[lx-1] - 360
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
                        while sortnum < sort_count: #"loops" thru identified storms and matches up lons
                            
                            if storm_list[0] >= 180:
                                l1 = storm_list[0] - 360
                            else:
                                l1 = storm_list[0]
                            if storm_sorter[sortnum][-1][0] >= 180:
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

        
        print(round((timeIN.time() - month_time_start)/60,2), 'minutes')

### storm_sorter --> all_storms

#%% chunk storms and organize
    print('--> chunk')
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
    print('--> check for duplicates')
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
    print('--> census info')

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
        two_week = startdt + timedelta(days=14) # startdt vs end (same duration for all storms, 3wks)
        analysis_ranges.append(mf.daterange(week_ago, two_week, dt=24))
        storm_ranges.append(mf.daterange(startdt, enddt, dt=24))
    

    #%%% * START STORM LOOP * get contours
    print('--> storm loop')
    ### set up ice area export
    myvars= ['box_area', 'ice_area15', 'ice_area80']
    datavars, outvars = {}, {}
    for var in myvars:
          outvars[var]=[]
    
    print(' '); print(' ')
    stormstr_prev = ''
    new_storm_info, new_storm_ranges = [],[]
    new_all_storms = []
    alph = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    dupe = iter(alph*5)
    for storm_num, storm_event in enumerate(storm_ranges):         ###!!!!!!!
        print();print('-');
        print('Storm ' + str(storm_num+1) + '/' + str(len(storm_ranges)))
        print('-')
        
        ### get needed contours and create box
        stormstr1 = storm_event[0].strftime('%Y_%m%d')
        daymonth_grid = [(dt.month, dt.day) for dt in analysis_ranges[storm_num]]
            
        storm_daymonth_grid = []
        for storm_day in storm_event:
            storm_daymonth_grid.append((storm_day.month, storm_day.day))
    
        try:
            all_cont_dts, all_contours = detect_contours2(storm_daymonth_grid, daymonth_grid, year,
                                                             storm_info=storm_info[storm_num], intervals=contour_levels) 
        except Exception as e:
            print(e)
            all_cont_dts, all_contours = [],[]
            print('$$$ ERROR = detect_contours, ', storm_daymonth_grid[0])
            print('- len storm_info: '+str(len(storm_info)))
            print('- storm_num: ' +str(storm_num))
            
        if all_contours == []: continue
    
        #### add in interaction with sea ice [20,80]
        ice_frac, out, bbox_edges = get_storm_area(all_contours, storm_event, 
                                                   ice_fname, ice_criteria, miz)
        
        ices, si_lon, si_lat = mf.load_seaice(ice_fname, year, storm_event[0].month, storm_event[0].day, latlon=True)
        import cartopy.crs as ccrs
        fig=plt.figure(figsize=[15,15]) 
        ax = plt.axes(projection=ccrs.SouthPolarStereo(central_longitude=0))
        ax.coastlines('50m',edgecolor='black',linewidth=0.75)
        ax.set_extent([-180,180, -53,-90], ccrs.PlateCarree())
        ax.set_title(stormstr1+': '+str(ice_frac), fontsize=36)
        ax.pcolormesh(si_lon, si_lat, ices, cmap='Blues_r', transform=ccrs.PlateCarree())
        ax.plot(bbox_edges[:,0], bbox_edges[:,1], lw=2, color='r',
                transform=ccrs.PlateCarree(), zorder=100)
        mf.plot_geocontour(ax, si_lon, si_lat, ices, levels=[0.15], color='k', lw=2, ls='solid')
        print('- plotted!')
        
        if (ice_frac<np.min(ice_lims) or ice_frac>np.max(ice_lims)):
            print('--> ice interaction: '+str(round(ice_frac,2))+'%'); print();print()
            continue
        else: 
            new_storm_info.append(storm_info[storm_num])
            new_storm_ranges.append(storm_event)
            new_all_storms.append(all_storms[storm_num])
            for var in myvars:
                outvars[var].append(out[var])
                
        #### export
        if stormstr1==stormstr_prev:
            stormstr = stormstr1 + next(dupe)
        else:
            stormstr=stormstr1
        stormstr_prev = stormstr1
    
        ### export *.nc!
        nc_name = stormstr+'_contours.nc'
        nc_in = nc_save_path+ nc_name
        export2nc(nc_in, all_contours, all_cont_dts, intervals=contour_levels)

    ### update storm list
    storm_info = new_storm_info 
    storm_ranges = new_storm_ranges
    all_storms = new_all_storms
    coords = {'nstorms':(['nstorms'], np.arange(1,len(storm_ranges)+1))
              }
    
    #### export ice frac info to nc
    for var in myvars:
        datavars[var] = (['nstorms'], outvars[var])
    ds = xr.Dataset(data_vars=datavars, 
                    coords=coords)
    
    ds.to_netcdf(nc_save_path+str(year)+'_area.nc')
    
    try: ds.close()
    except: pass
        
#%% export storm_info to new census database
    print('--> export to census')
    import csv
    header = ['start', 'start_lat', 'start_lon', 'minimum', 'finish', 'finish_lat', 'finish_lon']
    
    data = []
    for storm_num, storm_row in enumerate(storm_info):
        
        storm_data = [storm_ranges[storm_num][0]] # start date
        
        storm_data.append(round(storm_row[0][1],4)) # start lat
        storm_data.append(round(storm_row[0][0],4)) # start lon
        
        # find storm minimum pressure
        min_p = 9999
        for day_entry in storm_row:
            if day_entry[-1] < min_p:
                min_p = day_entry[-1]
        storm_data.append(round(min_p,3)) 
        
        storm_data.append(storm_ranges[storm_num][-1]) # finish date
        storm_data.append(round(storm_row[-1][1],4)) # finish lat
        storm_data.append(round(storm_row[-1][0],4)) # finsih lon
        
        # add data to full list for export
        data.append(storm_data)
    
    
    ###!!! create csv and export
    with open(nc_save_path+census_name +str(year)+'.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        
        # write the header
        writer.writerow(header)
    
        # write multiple rows
        writer.writerows(data)
    
#%% pickle some variables for later
    print(''); print('FINISHED STORM DATA COLLECTION: '+str(year)); print('')
    
    ### save important stuff for later
    import pickle
    
    output = open(nc_save_path + str(year)+'_all_storms.pkl', 'wb')
    pickle.dump(all_storms, output)
    output.close()
    
    # output2 = open(nc_save_path + str(year)+'_storm_sorter.pkl', 'wb')
    # pickle.dump(storm_sorter, output2)
    # output2.close()
    
    plt.close('all')
    
    outx=gc.collect()
    
#%% end
print(' '); print(' ')
print('------------------')
print('       done!      ')
print('------------------')
print('------------------')
print('elapsed time: ')
print((timeIN.time() - TIMESTART)/60)

