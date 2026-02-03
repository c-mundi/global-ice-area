"""
Created on Thu Aug 31 2023
storm_tseries2sens.py
-> running sensititvity studies, Oct 2025

updated storm timeseries calculation and export
SEA ICE ONLY

@author: mundi
"""
print('START', flush = True)
import os, sys
import numpy as np

if True: # imports
    import xarray as xr
    import netCDF4
    import pandas as pd
    from datetime import datetime, timedelta
    import glob
    
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from scipy.interpolate import griddata
    
    import myfuncs as mf
    import traceback
    
    import time as timeIN
    
print('- imports', flush = True)

#%% STARTING INFO
REMOTE=True

# years = [1988]
# years = np.arange(2010, 2015)
# years = np.arange(2015, 2020)
# years = np.arange(1982, 1986+1)
# years = np.arange(1987, 1991+1)

# years = np.arange(2010,2019+1)
# years = np.arange(1982,1991+1)

years = [2017]

clim_years = np.arange(2010,2019+1)
# clim_years = np.arange(1982,1991+1)


### SELECT STORM AREAS
run1000 = True
run990 = True
run_conts = True

### SELECT VARIABLES
run = {}
run['seaice'] = True

ncnameadd= '_seaice'


if REMOTE: root_path = '/home/mundi/'
else: root_path = '/Users/mundi/Desktop/'
root_path += 'month-hemi/nh_data/'+'sensitivity/'

savepath = root_path + 'seaice/' #'sensitivity_seaice/'

if not os.path.exists(savepath):
    os.makedirs(savepath)
    
print('- starting info', flush = True)


#%% functions
### loads pressure contours to form bbox
def open_cont_nc(ncfile):
    out = netCDF4.Dataset(ncfile)
    all_contours = []
    
    # loop thru all the variables
    for v in out.variables:
        all_contours.append(np.array(out.variables[v]))

    return all_contours

def open_cont_dt(ncfile):
    out = netCDF4.Dataset(ncfile)
    cont_dts = []
    
    # loop thru all the variables
    for v in out.variables:
        dt_str = v[7:16]
        cont_dts.append(datetime.strptime(dt_str,'%Y_%m%d'))
    return cont_dts     


def rstr(val,n=3):
    return str(round(val,n))

class HidePrint:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

print('- functions acquired', flush = True)


#%% * START YEAR LOOP *
print('', flush = True);print('', flush = True);print('', flush = True)

for year in years:
    tf = True
    savenc = tf
    TIMESTART = timeIN.time()
    
    original_runkeys = ['seaice', 'winds', 'sst','sst_daily', 'glorys', 'air_temp']
    for ogk in original_runkeys:
        try:
            x = run[ogk]
        except KeyError:
            run[ogk] = False
            
            
    print('********', flush = True)
    print('* '+str(year)+' *', flush = True) 
    print('********', flush = True)   
    print(run, flush = True); print('', flush = True)
    
    
    si_levels = [0.15, 0.80]
    
    if years[0]<1993:
        run['glorys'] = False #always

    #%% datapaths and functions
    if True: 
        print('- starting fpaths...', flush = True)
        # root_path += 'month-hemi/nh_data/'
        if REMOTE:
            nc_path =  root_path+'contours/'
            contour_name = '_contours.nc'
            print(' -a', flush = True)
            census_path = root_path+'census/'
            census_name = 'census_'+ str(year) +'.csv' 
           
            ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/daily/'  
        else:
            contour_name = '_contours.nc'
            census_name = 'census_'+ str(year) +'.csv'
            ice_fname = '/Users/mundi/Desktop/seaice/'
            census_path = '/Users/mundi/Desktop/month-hemi/nh_data/census/'
            nc_path = root_path +'contours/'
            
            print(' -local', flush = True)

        print('- filepaths established', flush = True)
            
        ### get starting sea ice
        _, si_lon, si_lat = mf.load_seaice(ice_fname, 1980, 10, 1, latlon=True)
        area = 25*25  
        
        print('- set up files', flush = True)
        

    #%% get storm info for this year
    """
    timing_grid = list of tuples: (start_date, end_date)
    storm_ranges = list of DT lists: [storm_day1, storm_day2, ... storm_day_f]
    analysis_ranges = list of DT lists: [storm_day1-one week, ... storm_day_f + two weeks]
    """
    
    ### get start/end dates + add'l time
    [startdate, enddate] = mf.readCensus(census_path+census_name, convertDT=True)[0]
    
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
    
        
    print('census acquired', flush = True)
        
    #%%% set up export nc
    
    analyrange_save, stormrange_save = [],[]
    
    sic_miz_list, miz_mask_list = [],[]
    sic_clim_list, miz_clim_list = [],[]
    
    total_wind_list, daily_wind_list = [],[]
    total_u_list, daily_u_list = [],[]
    total_v_list, daily_v_list = [],[]

    total_temp_list, daily_temp_list = [],[]
    
    sst_list = []
    
    daily_sst_list, total_sst_list = [],[]
    
    total_theta_list, daily_theta_list = [],[]
    
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
    
    print('- data export set up', flush = True)
    
    #%% * START STORM LOOP *
    print(' ', flush = True); print(' ', flush = True)
    stormstr_prev=''
    dupe_ind = 0
    for storm_num, storm_event in enumerate(storm_ranges): 
        # if storm_num==1: break ####!!!!
                                                                # len(new_storm_list)?
        print(str(year)+': Storm ' + str(storm_num+1) + '/' + str(len(storm_ranges)), flush = True)
        #%% get starting storm info
        stormstr1 = storm_event[0].strftime('%Y_%m%d')
        lengthstr = str(len(storm_event))
        
        ### duplicate storm start date?
        if stormstr1==stormstr_prev:
            stormstr = stormstr1 + str(storm_num)
            print('duplicate storm start: ')
            dupe_ind+=1
        else:
            stormstr=stormstr1
            dupe_ind=0
        stormstr_prev = stormstr1
        
            
        ### Save storm times
        if savenc:
            analyrange_save.append([analysis_ranges[storm_num][0].day,
                                    analysis_ranges[storm_num][-1].day])
            stormrange_save.append([storm_event[0].day, storm_event[-1].day])
            
        #%% * SELECT STORM AREAS
        keys = []
        short_keys = []
        nan_list = []
        
        if run1000: 
            keys.append('1000 hPa Box')
            short_keys.append('1000')
            nan_list.append(np.nan)
        if run990: 
            keys.append('990 hPa Box')
            short_keys.append('990')
            nan_list.append(np.nan)
        if run_conts: 
            keys.append('1000 hPa Contour Area')
            short_keys.append('contour')
            nan_list.append(np.nan)
            
            
        #%% get all 15-80 points
        
        t1 = storm_event[0] - timedelta(days=1)
        t2 = storm_event[-1] + timedelta(days=1)
        storm_range = mf.daterange(t1, t2, dt=24)

        miz_points = np.zeros(np.shape(si_lon))

        for date in storm_range:
            sic = mf.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
            miz_points = np.where(((sic>=0.15) & (sic<=0.80)), 1, miz_points)
            
        #%% load storm areas
        ncname = stormstr + contour_name
        try:
            all_contours = open_cont_nc(nc_path+ncname)
            all_cont_dts = open_cont_dt(nc_path+ncname)
        except FileNotFoundError:
            files = glob.glob(nc_path+'*'+storm_event[0].strftime('%Y_%m%d')+'*'+contour_name)
            try:
                all_contours = open_cont_nc(files[dupe_ind])
                all_cont_dts = open_cont_dt(files[dupe_ind])
            except (IndexError, FileNotFoundError):
                print()
                print('MISSING CONTOUR FILE: '+ncname)
                
                sic_miz_series = {}
                miz_mask_series = {}
                sic_miz_clim = {}
                miz_mask_clim = {}
                for key in keys:
                    sic_miz_series[key] = [np.nan]*22
                    miz_mask_series[key] = [np.nan]*22
                    sic_miz_clim[key] = [np.nan]*22
                    miz_mask_clim [key] = [np.nan]*22
                
                sic_miz_list.append(sic_miz_series)
                miz_mask_list.append(miz_mask_series)
                sic_clim_list.append(sic_miz_clim)
                miz_clim_list.append(miz_mask_clim)
                continue
        
        ### get bboxes : 1000,990
        if run1000:
            with HidePrint(): bbox_edges = mf.get_bbox_edges(all_contours) 
            inside_points1000 = mf.find_points_in_contour(bbox_edges, si_lon, si_lat)
            
        if run990:
            try:
                ds = xr.open_dataset(nc_path+ncname)
            except FileNotFoundError:
                files = glob.glob(nc_path+'*'+storm_event[0].strftime('%Y_%m%d')+'*'+contour_name)
                try:
                    ds = xr.open_dataset(files[dupe_ind])
                except:
                    print();print('run990 contours...' +storm_event[0].strftime('%Y_%m%d'))
                    print(nc_path+ncname); print()
                
            contours990 = []
            try:
                for v in list(ds.keys()):
                    if v[-5:-2] == '990':
                        contours990.append(np.array(ds[v]))
                if len(contours990)==0: raise IndexError
            except:
                for v in list(ds.keys()):
                    if v[-3:] == '990':
                        contours990.append(np.array(ds[v]))
            ds.close()
            # bounding area
            with HidePrint(): bbox_edges2 = mf.get_bbox_edges(contours990)
            inside_points990 = mf.find_points_in_contour(bbox_edges2, si_lon, si_lat)
        
        if run_conts:
            ### use 1000-contour area
            inside_points3 = np.full(np.shape(si_lon), False, dtype=bool)
            for contour in all_contours:
                mask1 = mf.find_points_in_contour(contour, si_lon, si_lat)
                # append mask1 to mask (and/or) --> combine to single mask
                inside_points3 = inside_points3 | ~mask1
            inside_points_contour =  ~inside_points3
        
        
        ### STORM AREA
        storm_area_list = []
        if run1000: 
            storm_area_list.append(inside_points1000)
        if run990: 
            storm_area_list.append(inside_points990)
        if run_conts: 
            storm_area_list.append(inside_points_contour)
        
            
        #%%--------------------------- get data
        
        #%% SEA ICE AREA
        if run['seaice']:
            sic_miz_series = {}
            miz_mask_series = {}
            sic_miz_clim = {}
            miz_mask_clim = {}
            for key in keys:
                sic_miz_series[key] = []
                miz_mask_series[key] = []
                sic_miz_clim[key] = []
                miz_mask_clim [key] = []
            
            try:
               for dt in analysis_ranges[storm_num]:
                   year1 = dt.year
                   month1= dt.month
                   day1 = dt.day
                   
                   try:
                       # get daily sea ice
                       sic = mf.load_seaice(ice_fname, year1, month1, day1, latlon=False)
                       
                       # miz ice
                       sic_miz = mf.getIce_TwoContours(ice_fname, 0.15, 0.80, year1, month1, day1)
                       
                       # all marginal points
                       miz_masked = np.ma.masked_where(miz_points==0, sic).filled(np.nan)
                       
                       ### area list
                       sic_areas = []
                       miz_areas = []
                       
                       # restrict area - 1000 box  & sum
                       if run1000:
                           sic_miz_1000 = np.nansum((np.ma.masked_array(sic_miz, mask=inside_points1000).filled(np.nan))*area)
                           miz_masked_1000 = np.nansum((np.ma.masked_array(miz_masked, mask=inside_points1000).filled(np.nan))*area)
                           sic_areas.append(sic_miz_1000)
                           miz_areas.append(miz_masked_1000)
                       # restrict area - 990 box  & sum
                       if run990:
                           sic_miz_990 = np.nansum((np.ma.masked_array(sic_miz, mask=inside_points990).filled(np.nan))*area)
                           miz_masked_990 = np.nansum((np.ma.masked_array(miz_masked, mask=inside_points990).filled(np.nan))*area)
                           sic_areas.append(sic_miz_990)
                           miz_areas.append(miz_masked_990)
                        
                       # restrict area - 1000 contours  & sum
                       if run_conts:
                           sic_miz_c = np.nansum((np.ma.masked_array(sic_miz, mask=inside_points_contour).filled(np.nan))*area)
                           miz_masked_c = np.nansum((np.ma.masked_array(miz_masked, mask=inside_points_contour).filled(np.nan))*area)
                           sic_areas.append(sic_miz_c)
                           miz_areas.append(miz_masked_c)
                        
                       # append to timeseries
                       for cc, carea in enumerate(sic_areas):
                           sic_miz_series[keys[cc]].append(carea)
                       for cc, carea in enumerate(miz_areas):
                           miz_mask_series[keys[cc]].append(carea)
                   except:
                        # append to timeseries
                        for cc, carea in enumerate(nan_list):
                            sic_miz_series[keys[cc]].append(np.nan)
                        for cc, carea in enumerate(nan_list):
                            miz_mask_series[keys[cc]].append(carea)
                   
                   
                   # calculate climatology
                   for cc, inside_points in enumerate(storm_area_list):
                       miz_yr, all_yr = [],[]
                       for yr in clim_years:
                           if month1==2 and day1==29:
                               sict1 = mf.load_seaice(ice_fname, yr, 2,28, latlon=False)
                               sict2 = mf.load_seaice(ice_fname, yr, 3,1, latlon=False)
                               sict = np.nanmean(np.array([sict1,sict2]), axis=0)
                           else:
                               sict = mf.load_seaice(ice_fname, yr, month1, day1, latlon=False)
                           # miz
                           try:
                               if month1==2 and day1==29: day2=28
                               else: day2=day1
                               si_yr = mf.getIce_TwoContours(ice_fname, 0.15, 0.80, yr, month1, day2)
                               si_yr_in = np.nansum(np.ma.masked_array(si_yr, mask=inside_points).filled(np.nan)*area)
                               miz_yr.append(si_yr_in)
                           except:
                               miz_yr.append(np.nan)
                           # all marginal points
                           try:
                               miz_mask = np.ma.masked_where(miz_points==0, sict).filled(np.nan)
                               all_yr.append( np.nansum(np.ma.masked_array(miz_mask, mask=inside_points).filled(np.nan)*area) )
                           except:
                               all_yr.append(np.nan)
                       sic_miz_clim[keys[cc]].append(np.nanmean(miz_yr))
                       miz_mask_clim[keys[cc]].append(np.nanmean(all_yr))
                 
               if savenc:  
                    sic_miz_list.append(sic_miz_series)
                    miz_mask_list.append(miz_mask_series)
                    sic_clim_list.append(sic_miz_clim)
                    miz_clim_list.append(miz_mask_clim)
               print('... sea ice area completed')
                
            except Exception as e:
                print('')
                print('*ERROR* ... SEAICEAREA ... ' + stormstr)
                print(e)
                print(traceback.format_exc())
                print('')
                 
        #%% ------------------------- end indiv storm
        
        ### print output for comparison
        # print('')
        # print('--------------------------')
        print('DONE CALCULATING:', end=' ')
        print('nstorms,', np.shape(np.arange(1,len(storm_ranges)+1)), end=' - ')
        print('analysis ranges,',np.shape(analysis_ranges[storm_num]))
        if run['seaice']: print('sea ice,', np.shape(sic_miz_series[keys[0]]), np.shape(sic_miz_list))
        
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
            if run['seaice']:
                try:
                    for kk, key in enumerate(keys):
                        sic1, sic2 = [],[]
                        clim1, clim2 = [],[]
                        short = short_keys[kk]
                        for list1, list2, list3, list4 in zip(sic_miz_list, miz_mask_list, sic_clim_list, miz_clim_list): 
                                sic1.append(list1[key]) 
                                sic2.append(list2[key])
                                clim1.append(list3[key])
                                clim2.append(list4[key])
                            
                    
                        if np.shape(sic1) == (len(storm_ranges), len(time_list)):
                            coord1, coord2 = 'nstorms','time'    
                        else:
                            print('* using seaice coords: ', np.shape(sic1), np.shape(sic2), len(storm_ranges))
                            coords = {'si_x': (['si_x'], np.arange(1,np.shape(sic1)[0]+1)),
                                      'si_y':(['si_y'], np.arange(1,np.shape(sic1)[1]+1))
                                      }
                            coord1, coord2 = 'si_x','si_y'
                            
                        data_vars['sia_miz_'+short] = ([coord1, coord2], sic1, 
                                 {'units': 'km^2', 
                                  'long_name':'sea ice area timeseries, original miz, '+key})
                        data_vars['sia_miz2_'+short] = ([coord1, coord2], sic2, 
                                 {'units': 'km^2', 
                                  'long_name':'sea ice area timeseries, total miz area, '+key})
                        
                        data_vars['sia_clim_miz_'+short] = ([coord1, coord2], clim1, 
                                 {'units': 'km^2', 
                                  'long_name':'sea ice area climatology, original miz, '+key})
                        data_vars['sia_clim_miz2_'+short] = ([coord1, coord2], clim2, 
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
        
        ds.to_netcdf(savepath+ str(year) + ncnameadd+ '.nc')
        
        print('\n netcdf saved!')
    except Exception as ee:
        print('')
        print('*ERROR* ... saving nc file')
        print('--------------------------')
        # print('analysis ranges, ',np.shape(analysis_str))
        if run['seaice']: print('sea ice, ', np.shape(sic1), np.shape(clim1))
        
        print(''); print('')
        print(ee)
        print(traceback.format_exc())
        print('')
    
            
    #%% end
    print(' '); print(' ')
    print('------------------')
    print('       done!      ')
    print('------------------')
    print('------------------')
    print(str(year)+' elapsed time: ')
    print(rstr((timeIN.time() - TIMESTART)/60,1), 'minutes')  
        
