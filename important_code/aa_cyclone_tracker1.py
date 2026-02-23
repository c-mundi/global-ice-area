#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nov 15 2023
aa_cyclone_tracker1.py

> do ice criteria need to change? (some percent of 80% ice... "thinner")
> total ice fraction?

- check ><= with cesm (storm_area|miz)

@author: mundi
"""
#%% inputs
from datetime import datetime, timedelta
import numpy as np

years = np.arange(2010,2019+2) #np.arange(1982, 1991+2) #np.arange(1982, 1991+2) #[2010] 

months = [12,1,2,3] #[12] #[1,2]


nc_save_path = '/Users/mundi/Desktop/antarctica/test_957_all/' #'_any2/'
# nc_save_path = '/home/mundi/antarctica/testout2/'

ice_fname = '/Users/mundi/Desktop/seaice/south/'
# ice_fname = '/boltzmann/data5/arctic/NOAA_CDR_SIC_V4/south/'

# eraname = '2010-2020_ndjf_era5-msl-hourly' #'1982-1992_ndjf_era5-msl-hourly'
# LOCALFILE = '/Users/mundi/Desktop/antarctica/DATA/'+eraname+'.nc'

LOCALPATH = '/Users/mundi/Desktop/era_data/'

census_name = 'census_test' #'test_' + 

#%% pressure threshold for storm lifetime
# weaker min pressure at storm start and storm final (984)
pres_thresh = 6

min_pres = 957 #960 #984 #968-bsckground_slp

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
    
#%% find locations of minimum pressures
# returns all_storms: lon/lat/minpressure each day
print('--> daily min pressure locations')

for year in years:
    
    storm_counter = 1

    storm_sorter = {}
    sort_count = 0
    
    for month in months: 
        plt.close('all')
        day_counter = -1
        my_storms_daily = []
        
        print('*** '+str(month)+'-'+str(year)+' ***')
        month_time_start = timeIN.time()
        
        LOCALFILE = LOCALPATH+'msl_'+str(year)+'-'+str(month)+'.nc'
        
        # load monthly data from server ###!!!
        # calendar_days = list(np.arange(1,calendar.monthrange(year, month)[1]+1))
        calendar_days = np.arange(1, 29, 1)
        alldays = [str(d) if d>=10 else '0'+str(d) for d in calendar_days]
        
        if not LOCALFILE:
            TIME, lon, lat, MSL = mf.era5_hourly([str(year)], [month], alldays)
        else:
            import pandas as pd
            ds1 = xr.open_dataset(LOCALFILE)
            TIME = pd.to_datetime( ds1['time'].values )
            lon, lat = np.meshgrid(ds1['longitude'].values, ds1['latitude'].values)
            MSL = ds1['msl'].values/100
            ds1.close()
        
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
    dupe = iter(['a','b','c','d','e','f','g','h','i','j','k','l','m'])
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
            all_cont_dts, all_contours = mf.detect_contours2(storm_daymonth_grid, daymonth_grid, year,
                                                             storm_info=storm_info[storm_num], intervals=contour_levels,
                                                             local_path=LOCALPATH, local_file='msl_') 
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
    
    
    ### create csv and export
    with open(nc_save_path+census_name +'_'+str(year)+'.csv', 'w', encoding='UTF8', newline='') as f:
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
    
#%% end
print(' '); print(' ')
print('------------------')
print('       done!      ')
print('------------------')
print('------------------')
print('elapsed time: ')
print((timeIN.time() - TIMESTART)/60)

raise NotImplementedError("End of Code")

#%% -------------
#%% TEST CODE
import cartopy.crs as ccrs
import myfuncs_aa as mf
import matplotlib

TEST= True

year, month = 2011, 12

if TEST:
    
    test_census = nc_save_path+census_name + str(year)+'.csv'
    #'/Users/mundi/Desktop/antarctica/census_test_2014.csv'
    
    test_inds = np.arange(24*0,24*14,12)
    
    [startdate, enddate], [[startlon, startlat],[endlon,endlat]], pressure = \
        mf.readCensus(test_census, convertDT=True)
    
    contour_levels = [980, 970, 960, 950]
    ccolors = ['k', 'b','lime', 'yellow']
    
    cmap = matplotlib.colormaps.get_cmap('rainbow')
    norm = matplotlib.colors.Normalize(vmin=0, vmax=len(startdate)-1)
    colors = [cmap(norm(ee)) for ee in range(len(startdate))]
    
    
    LOCALFILE = LOCALPATH+'msl_'+str(year)+'-'+str(month)+'.nc'
    ds1 = xr.open_dataset(LOCALFILE)
    TIME = pd.to_datetime( ds1['time'].values )
    lon, lat = np.meshgrid(ds1['longitude'].values, ds1['latitude'].values)
    MSL = ds1['msl'].values/100
    ds1.close()
    ttt = TIME[test_inds]
    
    for slpi, slp11 in enumerate(MSL[test_inds,:,:]): 
        ax=mf.background_plot()
        ax.set_title(str(ttt[slpi]), fontsize=30)
        ax.pcolormesh(lon,lat,np.where(slp11>min_pres,np.nan, slp11),transform=ccrs.PlateCarree())
        for cl, cc in zip(contour_levels, ccolors):
            ax.contour(lon, lat, slp11, levels=[cl], colors=cc, transform=ccrs.PlateCarree())
        
        for hi, sdt in enumerate(startdate):
            if ttt[slpi].day in np.unique([xx.day for xx in mf.daterange(sdt, enddate[hi])]):
                ax.plot(startlon[hi], startlat[hi],'o',transform=ccrs.PlateCarree(), 
                        markersize=30, color=colors[hi])
                ax.plot(endlon[hi], endlat[hi],'^',transform=ccrs.PlateCarree(), 
                        markersize=30, color=colors[hi])
                ax.plot(np.linspace(startlon[hi], endlon[hi],30), 
                        np.linspace(startlat[hi], endlat[hi],30),
                        transform=ccrs.PlateCarree(), color=colors[hi])

TEST2 = False
if TEST2:
     test_stormstr = '2014_1226'
     tds = xr.open_dataset(nc_save_path+test_stormstr+'_contours.nc')
     tkeys = list(tds.keys())
     tds.close()

#%%% .
raise NotImplementedError("End of Code")
#%%% flip through identified contours

for year in [2010]:
    
    # get start/end dates + add'l time
    [startdate, enddate] = mf.readCensus(nc_save_path+census_name +'_'+ str(year)+'.csv' , convertDT=True)[0]
    
    timing_grid = []
    for xx in range(0,len(startdate)):
        timing_grid.append((startdate[xx], enddate[xx]))
    
    storm_ranges = []
    analysis_ranges = []
    for startdt, enddt in timing_grid:
        week_ago = startdt - timedelta(days=7)
        two_week = startdt + timedelta(days=14) # relative to start date, since different storm lengths
        analysis_ranges.append(mf.daterange(week_ago, two_week, dt=24))
        storm_ranges.append(mf.daterange(startdt, enddt, dt=24))         
    
    #### * START STORM LOOP *
    print(' ')
    stormstr_prev=''
    dupe = iter(['a','b','c','d','e','f','g','h','i','j','k','l','m'])
    for storm_num, storm_event in enumerate(storm_ranges): 
        print(str(year)+': Storm ' + str(storm_num+1) + '/' + str(len(storm_ranges)), flush = True)
        stormstr1 = storm_event[0].strftime('%Y_%m%d')
        lengthstr = str(len(storm_event))
        
        # duplicate storm start date?
        if stormstr1==stormstr_prev:
            stormstr = stormstr1 + next(dupe)
            print('duplicate storm start')
        else:
            stormstr=stormstr1
        stormstr_prev = stormstr1
        
        daymonth_grid = [(dt.month, dt.day) for dt in analysis_ranges[storm_num]]
        
        storm_daymonth_grid = []
        for storm_day in storm_event:
            storm_daymonth_grid.append((storm_day.month, storm_day.day))
            
        #### load storm areas
        ncname = stormstr + '_contours.nc'
        cs = xr.open_dataset(nc_save_path+ncname)
        all_contours = []
        for key in list(cs.keys()):
            coord = cs[key].values
            all_contours.append(coord)
        cs.close()
        
        if all_contours == []: print('empty contours'); continue
        
        ### get bboxes : 1000,990
        with mf.HidePrint(): 
            bbox_edges = mf.get_bbox_edges(all_contours) 
        inside_points1000 = mf.find_points_in_contour(bbox_edges, si_lon, si_lat)



        ces, si_lon, si_lat = mf.load_seaice(ice_fname, year, storm_event[0].month, storm_event[0].day, latlon=True)
        fig=plt.figure(figsize=[15,15]) 
        ax = plt.axes(projection=ccrs.SouthPolarStereo(central_longitude=0))
        ax.coastlines('50m',edgecolor='black',linewidth=0.75)
        ax.set_extent([-180,180, -53,-90], ccrs.PlateCarree())
        ax.set_title(stormstr1, fontsize=36)
        ax.pcolormesh(si_lon, si_lat, ices, cmap='Blues_r', transform=ccrs.PlateCarree())
        ax.plot(bbox_edges[:,0], bbox_edges[:,1], lw=2, color='r',
                transform=ccrs.PlateCarree(), zorder=100)
        for co in all_contours:
            ax.plot(co[:,0], co[:,1], lw=2, color='k',
                    transform=ccrs.PlateCarree(), zorder=100)
        mf.plot_geocontour(ax, si_lon, si_lat, ices, levels=[0.15], color='k', lw=2, ls='solid')
        print('- plotted!')






#%% end 2
