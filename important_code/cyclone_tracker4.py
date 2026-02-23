#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on July 29 2022
Last run 9/15/2022
Update new years 6/2/2023
cyclone_tracker4.py

near-final version for remote runs

@author: mundi
"""
import numpy as np
#%% inputs

# nc_save_path = 'C:/Users/mundi/Documents/RESEARCH/Summer2022/cyclone_tracking_test/exports/'
# nc_save_path = 'C:/Users/mundi/Documents/RESEARCH/Summer2022/deleteme2/'
# ice_fname = 'C:/Users/mundi/Documents/RESEARCH/drive/SeaIce/daily/'

# nc_save_path = '/home/mundi/Summer2022/cyclone_tracker_out2/'
# nc_save_path = '/home/mundi/Summer2022/todel/census_test/'
# ice_fname = '/home/mundi/NOAA_CDR_SIC_V4/daily/' #2010-2020

# nc_save_path = '/home/mundi/Summer2022/cyclone_tracker_out_may/'
census_save_path = '/home/mundi/cyclones_allmonths/'
nc_save_path = census_save_path+'contours/'

# ice_fname = '/boltzmann/data5/arctic/NOAA_CDR_SIC_V4/daily/'
ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/daily/'    


census_name = 'census2' #'test_' + 

# years = list(np.arange(1990,1999+1)) + list(np.arange(2000,2009+1))
# years = list(np.arange(2010, 2020)) + list(np.arange(1980, 1990))
# years = np.arange(2010,2023)
# years = [2020]
# years = [1990,1991]
years = [1989]
# years = [1980, 1981, 1982, 1983, 1984]
# years = [1985, 1986, 1987, 1988, 1989]

# months = [5] #[6,7,8,9,10]
months = np.arange(1,12+1)
# months = np.arange(1,5+1)


#%% imports
if True:
    from datetime import datetime, timedelta
    import calendar
    import matplotlib.pyplot as plt
    
# import sys
# sys.path.append('/Users/mundi/Desktop/python/summer2022/')
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
                
            print('coords_' +stormstr+'_'+cp[c%2])
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
    
#%% pressure threshold for storm lifetime
# weaker min pressure at storm start and storm final
pres_thresh = 6

min_pres = 984

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
        TIME, lon, lat, MSL = mf.era5_hourly([str(year)], [str(month)], alldays)
        
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
    for storm_num, storm_event in enumerate(storm_ranges):         ###!!!!!!!
        print('Storm ' + str(storm_num+1) + '/' + str(len(storm_ranges)))
        
        ### get needed contours and create box
        stormstr1 = storm_event[0].strftime('%Y_%m%d')
        daymonth_grid = [(dt.month, dt.day) for dt in analysis_ranges[storm_num]]
            
        storm_daymonth_grid = []
        for storm_day in storm_event:
            storm_daymonth_grid.append((storm_day.month, storm_day.day))
    
        all_cont_dts, all_contours = mf.detect_contours2(storm_daymonth_grid, daymonth_grid, year, storm_info=storm_info[storm_num], local_file='') 
    
        if stormstr1==stormstr_prev:
            stormstr = stormstr1 + str(storm_num)
        else:
            stormstr=stormstr1
        stormstr_prev = stormstr
    
        ### export *.nc!
        nc_name = stormstr+'_contours2.nc'
        nc_in = nc_save_path+ nc_name
        
        export2nc(nc_in, all_contours, all_cont_dts)


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
    with open(census_save_path+census_name +'_'+str(year)+'.csv', 'w', encoding='UTF8', newline='') as f:
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
    
    output = open(census_save_path + str(year)+'_all_storms2.pkl', 'wb')
    pickle.dump(all_storms, output)
    output.close()
    
    # output2 = open(nc_save_path + str(year)+'_storm_sorter.pkl', 'wb')
    # pickle.dump(storm_sorter, output2)
    # output2.close()
    

#%% end
print(' '); print(' ')
print('------------------')
print('       done!      ')
print('------------------')
print('------------------')
print('elapsed time: ')
print((timeIN.time() - TIMESTART)/60)
