#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 2024
aa_storm_tseries1_sens.py
-> antarctic storm sensitivities

- sea ice only
- different storm areas

@author: mundi
"""
import os
import numpy as np
import glob

import xarray as xr
from datetime import datetime, timedelta

import traceback
import time as timeIN

import myfuncs_aa as mf

#%% STARTING INFO
REMOTE=True

# years = np.arange(2010, 2019+1)
# years = np.arange(1982, 1991+1)

# years = np.arange(2010, 2015)
# years = np.arange(2015, 2020)

# years = np.arange(1982, 1987)
# years = np.arange(1987, 1992)

# years = [2013,2014]
years = [1990]

# clim_years = np.arange(2010,2019+1)
clim_years = np.arange(1982,1991+1)

ncnameadd= '_seaice'

### SELECT STORM AREAS
run1000 = True
run990 = True
run_conts = True

### SELECT VARIABLES
run = {}
run['seaice'] = True


if REMOTE:
    root_path = '/home/mundi/month-hemi/sh_data/'#+'sensitivity/'
    ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/south/'    
else:
    root_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
    ice_fname = '/Users/mundi/Desktop/seaice/south/'
    
nc_path = root_path+'contours/'
contour_name = '_contours.nc'
census_path = root_path+'census/'
census_name = 'census_'

savepath = root_path+'sensitivity_seaice/' # 'seaice/' #
if not os.path.exists(savepath):
    os.makedirs(savepath)
    
print('- filepaths')

### get starting sea ice
_, si_lon, si_lat = mf.load_seaice(ice_fname, 1980, 10, 1, latlon=True)
area = 25*25  

miz = [0.15,0.80] #[0,0.58] #
    
#%% * START YEAR LOOP *
print('')

for year in years:
    SAVENAME = savepath+ str(year) + ncnameadd+ '.nc'
    
    # try:
    #     xr.open_dataset(SAVENAME)
    #     continue
    # except FileNotFoundError:
    #     pass
    
    savenc = True
    TIMESTART = timeIN.time()
    
    original_runkeys = ['seaice', 'winds', 'sst','sst_daily', 'glorys', 'air_temp']
    for ogk in original_runkeys:
        try:
            x = run[ogk]
        except KeyError:
            run[ogk] = False
            
    print('********')
    print('* '+str(year)+' *') 
    print('********')   
    print(run); print('')
    
    
    #%% get storm info for this year
    """
    timing_grid = list of tuples: (start_date, end_date)
    storm_ranges = list of DT lists: [storm_day1, storm_day2, ... storm_day_f]
    analysis_ranges = list of DT lists: [storm_day1-one week, ... storm_day_f + two weeks]
    """
    
    ### get start/end dates + add'l time
    [startdate, enddate] = mf.readCensus(census_path+census_name + str(year)+'.csv' , convertDT=True)[0]
    
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
        
    print('- census acquired')
    
    #%%% set up export nc
    
    analyrange_save, stormrange_save = [],[]
    
    miz_mask_list = []
    miz_clim_list = []
    
    ### initial xarray
    # define coordinates
    time_list = np.arange(-7,14+1,1)
    coords = {'time': (['time'], time_list),
              'nstorms':(['nstorms'], np.arange(1,len(storm_ranges)+1))
              }
    
    # define data with variable attributes
    data_vars = {}
    attrs = {'creation_date':datetime.now().strftime('%Y-%m-%d'), 
             'title':'Exported timeseries and storm info for ' + str(year)}
    
    print('- data export set up')
    
    #%% * START STORM LOOP *
    print(' ')
    stormstr_prev=''
    alph = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    dupe = iter(alph*5)
    for storm_num, storm_event in enumerate(storm_ranges): 
        print(str(year)+': Storm ' + str(storm_num+1) + '/' + str(len(storm_ranges)), flush = True)
        
        #%% get starting storm info
        stormstr1 = storm_event[0].strftime('%Y_%m%d')
        lengthstr = str(len(storm_event))
        
        ### duplicate storm start date?
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
            
        ### Save storm times
        if savenc:
            analyrange_save.append([analysis_ranges[storm_num][0].day,
                                    analysis_ranges[storm_num][-1].day])
            stormrange_save.append([storm_event[0].day, storm_event[-1].day])
            
        #%% get all 15-80 points
        
        t1 = storm_event[0] - timedelta(days=1)
        t2 = storm_event[-1] + timedelta(days=1)
        storm_range = mf.daterange(t1, t2, dt=24)
        
        miz_points = np.zeros(np.shape(si_lon))
        for date in storm_range:
            sic = mf.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
            miz_points = np.where(((sic>miz[0]) & (sic<=miz[1])), 1, miz_points)

        #%% load storm areas
        ncname = stormstr + contour_name
        cs = xr.open_dataset(nc_path+ncname)
        all_contours = []
        for key in list(cs.keys()):
            coord = cs[key].values
            all_contours.append(coord)
        cs.close()
        
        if all_contours == []: print('empty contours'); continue
        
        ### get bboxes : 1000,990
        if run1000:
            with mf.HidePrint(): bbox_edges = mf.get_bbox_edges(all_contours) 
            inside_points1000 = mf.find_points_in_contour(bbox_edges, si_lon, si_lat)
        
        if run990:
            try:
                ds = xr.open_dataset(nc_path+ncname)
            except FileNotFoundError:
                files = glob.glob(nc_path+'*'+storm_event[0].strftime('%Y_%m%d')+'*'+contour_name)
                try:
                    ds = xr.open_dataset(files[0])
                except:
                    print();print('run970 contours...' +storm_event[0].strftime('%Y_%m%d'))
                    print(nc_path+ncname); print()
                
            contours990 = []
            try:
                for v in list(ds.keys()):
                    if v[-5:-2] == '970':
                        contours990.append(np.array(ds[v]))
                if len(contours990)==0: raise IndexError
            except:
                for v in list(ds.keys()):
                    if v[-3:] == '970':
                        contours990.append(np.array(ds[v]))
            ds.close()
            # bounding area
            with mf.HidePrint(): bbox_edges2 = mf.get_bbox_edges(contours990)
            inside_points990 = mf.find_points_in_contour(bbox_edges2, si_lon, si_lat)
        
        if run_conts:
            ### use 1000-contour area
            inside_points3 = np.full(np.shape(si_lon), False, dtype=bool)
            for contour in all_contours:
                mask1 = mf.find_points_in_contour(contour, si_lon, si_lat)
                # append mask1 to mask (and/or) --> combine to single mask
                inside_points3 = inside_points3 | ~mask1
            inside_points_contour =  ~inside_points3
        
        #%% * SELECT STORM AREAS
        storm_area_list = []
        keys = []
        short_keys = []
        nan_list = []
        
        if run1000: 
            storm_area_list.append(inside_points1000)
            keys.append('1000 hPa Box')
            short_keys.append('1000')
            nan_list.append(np.nan)
            
        if run990: 
            storm_area_list.append(inside_points990)
            keys.append('990 hPa Box')
            short_keys.append('990')
            nan_list.append(np.nan)
        if run_conts: 
            storm_area_list.append(inside_points_contour)
            keys.append('1000 hPa Contour Area')
            short_keys.append('contour')
            nan_list.append(np.nan)
            
        #%%--------------------------- get data
        
        #%% SEA ICE AREA
        if run['seaice']:
            miz_mask_series = {}
            miz_mask_clim = {}
            for key in keys:
                miz_mask_series[key] = []
                miz_mask_clim [key] = []
            
            try:
               for month, day in daymonth_grid: 
                   try:
                       # get daily sea ice
                       sic = mf.load_seaice(ice_fname, year, month, day, latlon=False)
                       
                       # all marginal points
                       miz_masked = np.ma.masked_where(miz_points==0, sic).filled(np.nan)
                       
                       ### area list
                       miz_areas = []
                       
                       # restrict area - 1000 box  & sum
                       if run1000:
                           miz_masked_1000 = np.nansum((np.ma.masked_array(miz_masked, mask=inside_points1000).filled(np.nan))*area)
                           miz_areas.append(miz_masked_1000)
                           
                       # restrict area - 990 box  & sum
                       if run990:
                            miz_masked_990 = np.nansum((np.ma.masked_array(miz_masked, mask=inside_points990).filled(np.nan))*area)
                            miz_areas.append(miz_masked_990)
                         
                       # restrict area - 1000 contours  & sum
                       if run_conts:
                            miz_masked_c = np.nansum((np.ma.masked_array(miz_masked, mask=inside_points_contour).filled(np.nan))*area)
                            miz_areas.append(miz_masked_c)

                       # append to timeseries
                       for cc, carea in enumerate(miz_areas):
                           miz_mask_series[keys[cc]].append(carea)
                           
                           
                   except Exception as e:
                        # append to timeseries
                        print(str(month)+'-'+str(day))
                        print(e)
                        for cc, carea in enumerate(nan_list):
                            miz_mask_series[keys[cc]].append(carea)
                   
                   
                   # calculate climatology
                   for cc, inside_points in enumerate(storm_area_list):
                       miz_yr, all_yr = [],[]
                       for yr in clim_years:
                           
                           if month==2 and day == 29:  # leap year
                               sic1 = mf.load_seaice(ice_fname, yr, 2, 28, latlon=False)
                               sic2 = mf.load_seaice(ice_fname, yr, 3, 1, latlon=False)
                               sict = np.nanmean([sic1,sic2],axis=0)

                           # all marginal points
                           try:
                               if month!=2 and day!=29:
                                   sict = mf.load_seaice(ice_fname, yr, month, day, latlon=False)
                               elif month==2 and day == 29:  # leap year
                                   sic1 = mf.load_seaice(ice_fname, yr, 2, 28, latlon=False)
                                   sic2 = mf.load_seaice(ice_fname, yr, 3, 1, latlon=False)
                                   sict = np.nanmean([sic1,sic2],axis=0)
                                   
                               miz_mask = np.ma.masked_where(miz_points==0, sict).filled(np.nan)
                               all_yr.append( np.nansum(np.ma.masked_array(miz_mask, mask=inside_points).filled(np.nan)*area) )
                           except:
                               all_yr.append(np.nan)
                       miz_mask_clim[keys[cc]].append(np.nanmean(all_yr))
                 
               if savenc:  
                    miz_mask_list.append(miz_mask_series)
                    miz_clim_list.append(miz_mask_clim)
               print('... sea ice area completed')
                
            except Exception as e:
                print('')
                print('*ERROR* ... SEAICEAREA ... ' + stormstr)
                print(e)
                print(traceback.format_exc())
                print('')

    #%% ------------------------- end indiv storm
    
    # print output for comparison
    print('')
    print('--------------------------')
    print('DONE CALCULATING:')
    print('nstorms, ', np.shape(np.arange(1,len(storm_ranges)+1)))
    print('analysis ranges, ',np.shape(analysis_ranges[storm_num]))
    if run['seaice']: print('sea ice, ', np.shape(miz_mask_series[keys[0]]))
    print('')
    print('--------------------------')
    print('')
        
    #%% export all storms to nc
    print('')
    print('--------------------------')
    print('      EXPORT TO NC        ')
    print('--------------------------') 
    print(''); print('')
    
    if savenc:
        try:
            # try:
            #     analysis_and_storm_ranges = None
            #     analysis_str = [[dt.strftime('%Y%m%d') for dt in range1] for range1 in analysis_ranges]
            #     stormranges_convert = [[dt.strftime('%Y%m%d') for dt in range1] for range1 in storm_ranges]
            #     stormrange_str = []
            #     for stormstrs in stormranges_convert:
            #         appendstr = list(stormstrs)
            #         while len(appendstr) < 22: 
            #             appendstr.append('-')
                        
            #         stormrange_str.append(appendstr)
    
            #     data_vars['analysis_ranges'] = (['nstorms','time'], analysis_str, 
            #                   {'long_name':'range of dates before/after storm that were analyzed'})
            #     # data_vars['storm_ranges'] = (['nstorms','time'], stormrange_str, 
            #     #               {'long_name':'storm dates'})
            # except Exception as ee:
            #     print('')
            #     print('*ERROR*  START... exporting variables to nc')
            #     print(ee)
            #     print(traceback.format_exc())
            #     print('')
            
            if run['seaice']:
                try:
                    for kk, key in enumerate(keys):
                        sic1, clim1 = [],[]
                        short = short_keys[kk]
                        for list1, list2 in zip(miz_mask_list, miz_clim_list): 
                                sic1.append(list1[key]) 
                                clim1.append(list2[key])
                            
                        if np.shape(sic1) == (len(storm_ranges), len(time_list)):
                            coord1, coord2 = 'nstorms','time'    
                        else:
                            print('* using seaice coords: ', np.shape(sic1), len(storm_ranges))
                            coords = {'si_x': (['si_x'], np.arange(1,np.shape(sic1)[0]+1)),
                                      'si_y':(['si_y'], np.arange(1,np.shape(sic1)[1]+1))
                                      }
                            coord1, coord2 = 'si_x','si_y'
                            
                        data_vars['sia_miz2_'+short] = ([coord1, coord2], sic1, 
                                 {'units': 'km^2', 
                                  'long_name':'sea ice area timeseries, total miz area, '+key})
                        
                        data_vars['sia_clim_miz2_'+short] = ([coord1, coord2], clim1, 
                                 {'units': 'km^2', 
                                  'long_name':'sea ice area climatology, total miz area, '+key})
                    
                except Exception as ee:
                    print('')
                    print('*ERROR* SEAICE AREA... exporting variables to nc')
                    print(ee)
                    print(traceback.format_exc())
                    print('') 
                 
        except Exception as ee:
            print('')
            print('*ERROR* ... exporting variables to nc')
            print(ee)
            print(traceback.format_exc())
            print('')
            
    #%% create dataset -> nc
    
    try:
        ds = xr.Dataset(data_vars=data_vars, 
                        coords=coords, 
                        attrs=attrs)
        
        ds.to_netcdf(SAVENAME)
        
        print('\n netcdf saved!')
    except Exception as e:
        print('')
        print('*ERROR* ... saving nc file')
        print('--------------------------')
        
        # print('analysis ranges, ',np.shape(analysis_str))
        if run['seaice']: print('sea ice, ', np.shape(sic1), np.shape(clim1))
        
        print(''); print('')
        print(e)
        print(traceback.format_exc())
        print('')
    
            
    #%% end
    print(' '); print(' ')
    print('------------------')
    print('       done!      ')
    print('------------------')
    print(str(year)+' elapsed time: ')
    print(mf.rstr((timeIN.time() - TIMESTART)/60,1), 'minutes')  
    print(); print()
    
