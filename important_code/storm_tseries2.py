"""
Created on Thu Aug 31 2023
storm_tseries2.py

updated storm timeseries calculation and export
- new miz metrics
- zonal winds, wind speeds, 2m air temp
- volume-mean ocean temperature

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

years = [1989]
# years = np.arange(2010, 2019+1)
# years = np.arange(2022, 2024+1)
# years = np.arange(2015, 2019+1)
# years = np.arange(1982, 1986+1)
# years = np.arange(1987, 1991+1)

# clim_years = np.arange(2010,2019+1,1)
clim_years = np.arange(1982,1991+1)


# ncnameadd = '_sample'
# ncnameadd='_ocean'
# ncnameadd = '_areas'
# ncnameadd='_winds2' #_airtemp
ncnameadd= '_seaice2'
# ncnameadd = '_final'

### SELECT STORM AREAS
run1000 = True
run990 = False
run_conts = False

### SELECT VARIABLES
run = {}
run['seaice'] = True
# run['winds'] = True
# run['sst'] = True # check daymonthgrid
# run['glorys'] = True
# run['air_temp'] = True
# run['sst_daily'] = True # check daymonth grid

if REMOTE:
    # savepath = '/home/mundi/Summer2022/tseries_sens/'+'decades/' #'final/' ###!!!
    savepath = '/home/mundi/cyclones_allmonths/seaice2/'
else:
    savepath = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/seaice2/'

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
        if REMOTE:
            # nc_path = '/home/mundi/Summer2022/cyctrack_sensitivity/'
            # nc_path = '/home/mundi/Summer2022/cyclone_tracker_out/'
            nc_path =  '/home/mundi/cyclones_allmonths/contours/'
            contour_name = '_contours.nc'
            print(' -a', flush = True)
            # census_path = '/home/mundi/Summer2022/cyctrack_sensitivity/'
            # census_path = nc_path
            census_path = '/home/mundi/cyclones_allmonths/'
            census_name = 'census_'+ str(year) +'.csv' 
            # print(' -b', flush = True)
            # if year in [1990, 1991]:
            #     print('NEW CENSUS PATH!!!') # cyc_track_out2
            #     census_path = '/home/mundi/Summer2022/cyclone_tracker_out2/'
            #     nc_path = census_path
            #     print(' -'+str(year), flush = True)
        
            ice_fname = '/langley/data12/old_data5/arctic/NOAA_CDR_SIC_V4/daily/'  
            # ice_fname = '/boltzmann/data5/arctic/NOAA_CDR_SIC_V4/daily/'
            # ice_fname = '/home/mundi/SeaIce/data/'
        else:
            # nc_path = '/Users/mundi/Desktop/teststuff/original_census/cyclone_tracker_out/'
            contour_name = '_contours.nc'
            # census_path = '/Users/mundi/Desktop/teststuff/original_census/'
            census_name = 'census_'+ str(year) +'.csv'
            ice_fname = '/Users/mundi/Desktop/seaice/'
            census_path = '/Users/mundi/Desktop/FINAL/cyclones_allmonths/p-adjust/'
            nc_path = census_path +'contours/'
            
            print(' -local', flush = True)
        
        if run['sst']:
            sst_path = '/home/mundi/Fall2021/sst/'
            if int(year) >= 1990:
                sst_fname = sst_path + 'sst.wkmean.1990-present.nc'
            elif int(year)>=1982: 
                sst_fname = sst_path + 'sst.wkmean.1981-1989.nc'
            else:
                run['sst'] = False
                print(' '); print('*** BAD SST YEAR: '+str(year)+' ***'); print(' ')
            print(' -oldsst', flush = True)
                
        if run['sst_daily']:
            if REMOTE:
                daily_sst_path = '/home/mundi/sst/daily/'
            else:
                daily_sst_path = '/Users/mundi/Desktop/data/Daily_SST/'
            daily_sst_fname = 'sst.day.mean.'
            print(' -dailysst', flush = True)
                
        print('- filepaths established', flush = True)
            
        ### get starting sea ice
        _, si_lon, si_lat = mf.load_seaice(ice_fname, 1980, 10, 1, latlon=True)
        area = 25*25  
        
        print('- set up files', flush = True)
        
    #%% yearly ocean info
    #%%% ocean database (GLORYS)
    if run['glorys']:
        pass
        # def copernicusmarine_datastore(dataset, username, password):
        #     from pydap.client import open_url
        #     from pydap.cas.get_cookies import setup_session
        #     import collections
        #     import collections.abc
        #     cas_url = 'https://cmems-cas.cls.fr/cas/login'
        #     session = setup_session(cas_url, username, password)
        #     session.cookies.set("CASTGC", session.cookies.get_dict()['CASTGC'])
        #     database = ['my', 'nrt']
        #     url = f'https://{database[0]}.cmems-du.eu/thredds/dodsC/{dataset}'
        #     try:
        #         data_store = xr.backends.PydapDataStore(open_url(url, session=session, user_charset='utf-8')) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits 
        #     except:
        #         url = f'https://{database[1]}.cmems-du.eu/thredds/dodsC/{dataset}'
        #         data_store = xr.backends.PydapDataStore(open_url(url, session=session, user_charset='utf-8')) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits
        #     return data_store

        # USERNAME = 'cmundi'
        # PASSWORD = 'Oceans!52' 
        # DATASET_ID = 'cmems_mod_glo_phy_my_0.083_P1D-m'
        # get_store=True

        # print('... starting')
        # while get_store:
        #     try:
        #         data_store = copernicusmarine_datastore(DATASET_ID, USERNAME, PASSWORD)
        #         get_store=False
        #     except Exception as e:
        #         print(e)
        #         timeIN.sleep(1.5)
        #         get_store=True
                
        # print('... accessed ocean reanlysis store')
        # oceanDB_time = timeIN.time()
    
    #%%% load daily sst for this year
    # daily_sst_path = '/home/mundi/sst/daily/'
    # daily_sst_fname = 'sst.day.mean.'
    
    if run['sst_daily']:
        dsst = xr.open_dataset(daily_sst_path+daily_sst_fname+str(year)+'.nc')
        dsst_lon, dsst_lat = np.meshgrid(dsst['lon'].values, dsst['lat'].values )
        daily_sst = dsst['sst']
        dsst.close()
        print('- acquired sst data!', flush = True)

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
    
    # new_storm_list = []
    # for stormi in storm_ranges:
    #     if stormi[0].month<=10:
    #         new_storm_list.append(stormi)
    #     else:
    #         continue
        
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
                print(); print()
                print('MISSING CONTOUR FILE: '+ncname)
                continue
        
        ### get bboxes : 1000,990
        if run1000:
            with HidePrint(): bbox_edges = mf.get_bbox_edges(all_contours) 
            inside_points1000 = mf.find_points_in_contour(bbox_edges, si_lon, si_lat)
            
        if run990:
            ds = xr.open_dataset(nc_path+ncname)
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
                           sict = mf.load_seaice(ice_fname, yr, month1, day1, latlon=False)
                           # miz
                           try:
                               si_yr = mf.getIce_TwoContours(ice_fname, 0.15, 0.80, yr, month1, day1)
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
                 
        #%% WINDS 
        if run['winds']:
            try:
                total_wind_series = {}
                daily_wind_series = {}
                total_u_series = {}
                daily_u_series = {}
                total_v_series = {}
                daily_v_series = {}
                
                for key in keys:
                    total_wind_series[key] = []
                    daily_wind_series[key] = []
                    total_u_series[key] = []
                    daily_u_series[key] = []
                    total_v_series[key] = []
                    daily_v_series[key] = []
                
                total_series_dicts = [total_wind_series, total_u_series, total_v_series]
                daily_series_dicts = [daily_wind_series, daily_u_series, daily_v_series]
                
                for dt in analysis_ranges[storm_num]:
                   daily_time, lon, lat, u_wind = mf.era5_daily(dt.year, dt.month, dt.day, variable = '10m_u_component_of_wind')
                   daily_time, lon, lat, v_wind = mf.era5_daily(dt.year, dt.month, dt.day, variable = '10m_v_component_of_wind')
                   u_wind = np.squeeze(u_wind); v_wind=np.squeeze(v_wind)
                   
                   wind_tot = np.sqrt((u_wind**2) + (v_wind**2))
                   
                   # regrid to si grid
                   wind_grd = griddata((lon.flatten(), lat.flatten()), wind_tot.flatten(),
                                           (si_lon.flatten(), si_lat.flatten()))
                   wind_grd = wind_grd.reshape(np.shape(si_lon))
                   u_grd = griddata((lon.flatten(), lat.flatten()), u_wind.flatten(),
                                           (si_lon.flatten(), si_lat.flatten()))
                   u_grd = u_grd.reshape(np.shape(si_lon))
                   v_grd = griddata((lon.flatten(), lat.flatten()), v_wind.flatten(),
                                           (si_lon.flatten(), si_lat.flatten()))
                   v_grd = v_grd.reshape(np.shape(si_lon))
                   
                   # get daily miz
                   sic_miz = mf.getIce_TwoContours(ice_fname, 0.15, 0.80, dt.year, dt.month, dt.day)
                   
                   # restrict area    # miz points
                   for kk, si_in in enumerate(storm_area_list):
                       for tw, thiswind in enumerate([wind_grd, u_grd, v_grd]):
                           wind_in_out = np.ma.masked_array(thiswind, mask=si_in).filled(np.nan)
                           wind_in_tot = wind_in_out.copy()
                           wind_in_daily = wind_in_out.copy()
                           
                           total_wind_miz = np.ma.masked_where(miz_points==0, wind_in_tot).filled(np.nan)
                           daily_wind_miz = np.ma.masked_array(wind_in_daily, mask=sic_miz.mask).filled(np.nan)
        
                           # wind timeseries
                           total_series_dicts[tw][keys[kk]].append(np.nanmean(total_wind_miz))
                           daily_series_dicts[tw][keys[kk]].append(np.nanmean(daily_wind_miz))
                    
                if savenc:
                    try:
                        total_wind_list.append(total_series_dicts[0])
                        daily_wind_list.append(daily_series_dicts[0])
                        total_u_list.append(total_series_dicts[1])
                        daily_u_list.append(daily_series_dicts[1])
                        total_v_list.append(total_series_dicts[2])
                        daily_v_list.append(daily_series_dicts[2])
                        print('... winds completed')
                    except:
                        print(stormstr, '... issue with appending winds to nc')

            except Exception as e:
                print('')
                print('*ERROR* ... WINDS ... ' + stormstr)
                print(e)
                print(traceback.format_exc())
                print('')
                
        #%% 2m air temp
        if run['air_temp']:
            total_temp_series = {}
            daily_temp_series = {}
            for key in keys:
                total_temp_series[key] = []
                daily_temp_series[key] = []
            try:
                for dt in analysis_ranges[storm_num]:
                   daily_time, lon, lat, temp = mf.era5_daily(dt.year, dt.month, dt.day, variable = '2m_temperature')
                   temp = np.squeeze(temp)
                   
                   # regrid to si grid
                   temp_grd = griddata((lon.flatten(), lat.flatten()), temp.flatten(),
                                           (si_lon.flatten(), si_lat.flatten()))
                   temp_grd = temp_grd.reshape(np.shape(si_lon))
                   
                   # get daily miz
                   sic_miz = mf.getIce_TwoContours(ice_fname, 0.15, 0.80, dt.year, dt.month, dt.day)
                   
                   # restrict area    # miz points
                   for kk, si_in in enumerate(storm_area_list):
                        t_in = np.ma.masked_array(temp_grd, mask=si_in).filled(np.nan)
                        t_in1 = t_in.copy()
                        t_in2 = t_in.copy()
                        
                        total_t_miz = np.ma.masked_where(miz_points==0, t_in1).filled(np.nan)
                        daily_t_miz = np.ma.masked_array(t_in2, mask=sic_miz.mask).filled(np.nan)
     
                        # temp timeseries
                        total_temp_series[keys[kk]].append(np.nanmean(total_t_miz))
                        daily_temp_series[keys[kk]].append(np.nanmean(daily_t_miz))
                 
                if savenc:
                    try:
                        total_temp_list.append(total_temp_series)
                        daily_temp_list.append(daily_temp_series)
                        print('... air temp completed')
                    except:
                        print(stormstr, '... issue with appending temps to nc')
            except Exception as e:
                print('')
                print('*ERROR* ... 2m air temp ... ' + stormstr)
                print(e)
                print(traceback.format_exc())
                print('')
                
        #%% SST
        if run['sst']:
            daymonth_grid=[]
            
            sst_series = {}
            for key in keys:
                sst_series[key] = []
            
            try:
                ds = xr.open_dataset(sst_fname)
                lon_in = ds['lon'].values
                lon_fix = np.where(lon_in>180, lon_in-360, lon_in)

                lon, lat = np.meshgrid(lon_fix, ds['lat'].values)
                time_in = ds['time'].values # convert?
                time = pd.to_datetime(time_in)
                sst = ds['sst'].values
                ds.close()
                
                
                t_ind = 0 # noaa
                while time[t_ind] < storm_event[0]:
                    t_ind += 1
                t_ind -= 2 # use previous file
                print(' '); print('Loading SSTs...')
                
                for month, day in daymonth_grid:
                    
                    try:
                        if datetime(year, month, day) >= time[t_ind+1] or day==daymonth_grid[0][1]: 
                            t_ind += 1
                            load_new = True
                            print(datetime(year, month, day), time[t_ind])
                        else:
                            load_new = False
                        
                        # get sea ice miz
                        sic_miz = mf.getIce_TwoContours(ice_fname, 0.15, 0.80, year, month, day)
                        
                        # convert sst grid to sea ice
                        if load_new:
                            sst_si = griddata((lon.flatten(), lat.flatten()), sst[t_ind].flatten(),
                                                    (si_lon.flatten(), si_lat.flatten()))
                            sst_si = np.reshape(sst_si, np.shape(si_lon))
                            
                        # restrict area - miz
                        sst_miz = np.ma.masked_array(sst_si, mask=sic_miz.mask).filled(np.nan)
                        
                        # restrict area - bbox
                        for kk, si_in in enumerate(storm_area_list):
                            try:
                                sst_in = np.ma.masked_array(sst_miz, mask=si_in).filled(np.nan)
                                # average temps
                                sst_series[keys[kk]].append(np.nanmean(sst_in))
                            except:
                                sst_series[keys[kk]].append(np.nan)
                    except Exception as e:
                            print(e)
                    
                if savenc:  
                    sst_list.append(sst_series)
                
                print('... sst completed')
                
            except Exception as e:
                print('')
                print('*ERROR* ... SST ... ' + stormstr)
                print(e)
                print(traceback.format_exc())
                print('')
                
        #%% DAILY SST
        # daily_sst_path = '/home/mundi/sst/daily/'
        # daily_sst_fname = 'sst.day.mean.'
        # dsst = xr.open_dataset(daily_sst_path+daily_sst_fname+str(year)+'.nc')
        # dsst_lon = dsst['lon'].values
        # dsst_lat = dsst['lat'].values 
        # daily_sst = dsst['sst']
        
        if run['sst_daily']:
            
            daymonth_grid=[]
            
            sst_series_dailymiz = {}
            sst_series_totalmiz = {}
            for key in keys:
                sst_series_dailymiz[key] = []
                sst_series_totalmiz[key]=[]
            try:
                for month, day in daymonth_grid:
                    myday= datetime(year,month,day)
                    datestr = myday.strftime('%Y-%m-%d')
                    
                    # get daily sea ice
                    mysst = daily_sst.sel(time=datestr).values
                    
                    # get sea ice miz
                    sic_miz = mf.getIce_TwoContours(ice_fname, 0.15, 0.80, year, month, day)
                    
                    # convert sst grid to sea ice
                    sst_si = griddata((dsst_lon.flatten(), dsst_lat.flatten()), mysst.flatten(),
                                            (si_lon.flatten(), si_lat.flatten()))
                    sst_si = np.reshape(sst_si, np.shape(si_lon))
                        
                    # restrict area - miz
                    sst_miz = np.ma.masked_array(sst_si, mask=sic_miz.mask).filled(np.nan)
                    total_sst_miz = np.ma.masked_where(miz_points==0, sst_si)
                    
                    # restrict area - bbox
                    for kk, si_in in enumerate(storm_area_list):
                        try:
                            sst_in = np.ma.masked_array(sst_miz, mask=si_in).filled(np.nan)
                            sst_in2 = np.ma.masked_array(total_sst_miz, mask=si_in).filled(np.nan)
                            # average temps
                            sst_series_dailymiz[keys[kk]].append(np.nanmean(sst_in))
                            sst_series_totalmiz[keys[kk]].append(np.nanmean(sst_in2))
                        except:
                            sst_series_dailymiz[keys[kk]].append(np.nan)
                            sst_series_totalmiz[keys[kk]].append(np.nan)
                    
                if savenc:  
                    daily_sst_list.append(sst_series_dailymiz)
                    total_sst_list.append(sst_series_totalmiz)
                
                print('... daily sst completed')
                
            except Exception as e:
                print('')
                print('*ERROR* ... DAILY SST ... ' + stormstr)
                print(e)
                print(traceback.format_exc())
                print('')
            
        #%% GLORYS
        if run['glorys']:
            print('*** uncomment ocean store')
            # total_theta_series = {}
            # daily_theta_series = {}
            # for key in keys:
            #     total_theta_series[key] = []
            #     daily_theta_series[key]=[]
            
            # try:
            #     incomplete_read=True
            #     ### Subsetting parameters
            #     while incomplete_read:
            #         try:
            #             TIME = slice(analysis_ranges[storm_num][0].strftime('%Y-%m-%d'),analysis_ranges[storm_num][-1].strftime('%Y-%m-%d'))
            #             DEPTH = slice(0,20)
            #             LATITUDE = slice(50,90)
            #             DS = xr.open_dataset(data_store).sel(time=TIME, latitude=LATITUDE).sel(depth=DEPTH)
            #             # get coords
            #             long, lat = np.meshgrid(DS['longitude'].values, DS['latitude'].values)
            #             depth = DS['depth'].values
            #             time = DS['time'].values
            #             ts = pd.to_datetime(time)
            #             # get temp
            #             temp = DS['thetao'].values
            #             ### close
            #             DS.close()
            #             incomplete_read = False
            #         except:
            #             incomplete_read= True
            #             oceanDB = (timeIN.time() - oceanDB_time)/60
            #             oceanDB_time = timeIN.time()
            #             reload = True
            #             while reload:
            #                 try:
            #                     print('... reloading GLORYS datastore: ', rstr(oceanDB,1), ' minutes')
            #                     data_store = copernicusmarine_datastore(DATASET_ID, USERNAME, PASSWORD)
            #                     reload = False
            #                 except Exception as e:
            #                     print('error loading: retry'); print(e)
            #                     timeIN.sleep(1.5)
            #                     reload=True
            #     ### get new grid masks
            #     shape3d = np.shape(temp[0])
            #     # miz
            #     miz_points_glorys = np.zeros(np.shape(long))
            #     for date in storm_range:
            #         sic = mf.load_seaice(ice_fname, date.year, date.month, date.day, latlon=False)
            #         sic = griddata((si_lon.flatten(), si_lat.flatten()), sic.flatten(),
            #                                 (long.flatten(), lat.flatten()))
            #         sic = np.reshape(sic, np.shape(long))
            #         miz_points_glorys = np.where(((sic>=0.15) & (sic<=0.80)), 1, miz_points_glorys)
            #     miz_points_glorys = np.ma.masked_where(miz_points_glorys==0, miz_points_glorys)
            #     miz_points_glorys = np.broadcast_to(miz_points_glorys.mask, shape3d)

            #     ### get bboxes : 1000,990
            #     glorys_storm_areas = []
            #     if run1000:
            #         with HidePrint(): bbox_edges = mf.get_bbox_edges(all_contours) 
            #         glorys_inside1000 = mf.find_points_in_contour(bbox_edges, long, lat)
            #         glorys_inside1000 = np.broadcast_to(glorys_inside1000, shape3d)
            #         glorys_storm_areas.append(glorys_inside1000)
            #     if run990:
            #         ds = xr.open_dataset(nc_path+ncname)
            #         contours990 = []
            #         try:
            #             for v in list(ds.keys()):
            #                 if v[-5:-2] == '990':
            #                     contours990.append(np.array(ds[v]))
            #             if len(contours990)==0: raise IndexError
            #         except:
            #             for v in list(ds.keys()):
            #                 if v[-3:] == '990':
            #                     contours990.append(np.array(ds[v]))
            #         ds.close()
    
            #         # bounding area
            #         with HidePrint(): bbox_edges2 = mf.get_bbox_edges(contours990)
            #         glorys_inside990 = mf.find_points_in_contour(bbox_edges2, long, lat)
            #         glorys_inside990 = np.broadcast_to(glorys_inside990, shape3d)
            #         glorys_storm_areas.append(glorys_inside990)
            #     if run_conts:
            #         ### use 1000-contour area
            #         inside_points3 = np.full(np.shape(long), False, dtype=bool)
            #         for contour in all_contours:
            #             mask1 = mf.find_points_in_contour(contour, long, lat)
            #             # append mask1 to mask (and/or) --> combine to single mask
            #             inside_points3 = inside_points3 | ~mask1
            #         glorys_inside_contour =  ~inside_points3
            #         glorys_inside_contour = np.broadcast_to(glorys_inside_contour, shape3d)
            #         glorys_storm_areas.append(glorys_inside_contour)
            #     ### begin loop
            #     print('timing: ', len(ts))
            #     for tslice, timing in enumerate(ts):
            #             # get slice
            #             temp1 = temp[tslice]

            #             # sea ice
            #             si_in = mf.load_seaice(ice_fname, year, timing.month, timing.day, latlon=False)
                        
            #             si = griddata((si_lon.flatten(), si_lat.flatten()), si_in.flatten(),
            #                                     (long.flatten(), lat.flatten()))
            #             si = np.reshape(si, np.shape(long))
            #             si_miz = np.ma.masked_where(si>0.80, si)
            #             si_miz2 = np.ma.masked_where(si_miz<0.15, si_miz)
                        
            #             si_3d_mask = np.broadcast_to(si_miz2.mask, np.shape(temp1))
            #             daily_miz_temp = np.ma.masked_array(temp1, mask=si_3d_mask).copy()
            #             miz_glorys_3d = np.broadcast_to(miz_points_glorys, np.shape(temp1))
            #             total_miz_temp = np.ma.masked_array(temp1, mask=miz_glorys_3d).copy()
                        
            #             for kk, si_in in enumerate(glorys_storm_areas):
            #                 try:
            #                     total_theta_series[keys[kk]].append(np.nanmean(np.ma.masked_array(total_miz_temp, mask = si_in).copy()))
            #                 except:
            #                     total_theta_series[keys[kk]].append(np.nan)
            #                 try:
            #                     daily_theta_series[keys[kk]].append(np.nanmean(np.ma.masked_array(daily_miz_temp, mask = si_in).copy()))
            #                 except:
            #                     daily_theta_series[keys[kk]].append(np.nan)
                        
            #     if savenc:
            #         total_theta_list.append(total_theta_series) 
            #         daily_theta_list.append(daily_theta_series)
                    
            #         for kk, si_in in enumerate([glorys_inside1000, glorys_inside990, glorys_inside_contour]):
            #             print(np.shape(total_theta_series[keys[kk]]))
            #             print(np.shape(daily_theta_series[keys[kk]]))
                    
            #     ### clear memory?
            #     del long, lat, depth, time, ts, temp
                
            # except Exception as e:
            #      print('')
            #      print('*ERROR* ... oceans-GLORYS ... ' + stormstr)
            #      print(e)
            #      print(traceback.format_exc())
            #      print('')
                        
                        
                        
                        
        #%% ------------------------- end indiv storm
        
        # print output for comparison
        print('')
        print('--------------------------')
        print('DONE CALCULATING:')
        print('nstorms, ', np.shape(np.arange(1,len(storm_ranges)+1)))
        print('analysis ranges, ',np.shape(analysis_ranges[storm_num]))
        if run['seaice']: print('sea ice, ', np.shape(sic_miz_series[keys[0]]))
        if run['winds']: print('winds, ',np.shape(total_wind_series[keys[0]]))
        if run['air_temp']: print('air temp, ',np.shape(total_temp_series[keys[0]]))
        if run['sst']: print('sst, ',np.shape(sst_series[keys[0]]))
        if run['sst_daily']: print('daily sst, ',np.shape(sst_series_totalmiz[keys[0]]))
        # if run['glorys']: print('ocean, ',np.shape(total_theta_series[keys[0]]))
        
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
            try:
                analysis_and_storm_ranges = None
                # analysis_str = [[dt.strftime('%Y%m%d') for dt in range1] for range1 in analysis_ranges]
                # stormranges_convert = [[dt.strftime('%Y%m%d') for dt in range1] for range1 in storm_ranges]
                # stormrange_str = []
                # for stormstrs in stormranges_convert:
                #     # appendstr = ['-']*(len(analysis_str[0])-len(stormranges_convert))
                    
                #     appendstr = list(stormstrs)
                #     while len(appendstr) < 22: 
                #         appendstr.append('-')
                        
                #     stormrange_str.append(appendstr)
                             
    
                # data_vars['analysis_ranges'] = (['nstorms','time'], analysis_str, 
                #              {'long_name':'range of dates before/after storm that were analyzed'})
                # data_vars['storm_ranges'] = (['nstorms','time'], stormrange_str, 
                #              {'long_name':'storm dates'})
            except Exception as ee:
                print('')
                print('*ERROR*  START... exporting variables to nc')
                print(ee)
                print(traceback.format_exc())
                print('')
            
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
                 
            if run['winds']:
                try:
                    for kk, key in enumerate(keys):
                        short = short_keys[kk]
                        wind1, wind2 = [], []
                        u1, u2 = [],[]
                        v1, v2 = [],[]
                        for list1, list2 in zip(total_wind_list, daily_wind_list): 
                                wind1.append(list1[key])
                                wind2.append(list2[key])
                                
                        if np.shape(wind1) == (len(storm_ranges), len(time_list)):
                            coord1, coord2 = 'nstorms','time'
                        else:
                            print('* using wind coords: ', np.shape(wind1), np.shape(wind2), len(storm_ranges))
                            coords = {'wind_x': (['wind_x'], np.arange(1,np.shape(wind1)[0]+1)),
                                      'wind_y':(['wind_y'], np.arange(1,np.shape(wind1)[1]+1))}
                            coord1, coord2 = 'wind_x', 'wind_y'
                            
                        data_vars['total_winds_'+short] = ([coord1, coord2],  wind1, 
                                 {'units': 'm/s', 
                                  'long_name':'wind timeseries, total miz area, '+key})
                        data_vars['daily_winds_'+short] = ([coord1, coord2],  wind2, 
                                 {'units': 'm/s', 
                                  'long_name':'wind timeseries, original miz, '+key})
                        
                        for list1, list2 in zip(total_u_list, daily_u_list): 
                                u1.append(list1[key])
                                u2.append(list2[key])
                        
                        if np.shape(u1) == (len(storm_ranges), len(time_list)):
                            coord1, coord2 = 'nstorms','time'
                        else:
                            print('* using wind coords: ', np.shape(u1), np.shape(u2), len(storm_ranges))
                            coords = {'wind_x': (['wind_x'], np.arange(1,np.shape(u1)[0]+1)),
                                      'wind_y':(['wind_y'], np.arange(1,np.shape(u1)[1]+1))}
                            coord1, coord2 = 'wind_x', 'wind_y'
                            
                        data_vars['total_u_'+short] = ([coord1, coord2],  u1, 
                                 {'units': 'm/s', 
                                  'long_name':'u-wind timeseries, total miz area, '+key})
                        data_vars['daily_u_'+short] = ([coord1, coord2],  u2, 
                                 {'units': 'm/s', 
                                  'long_name':'u-wind timeseries, original miz, '+key})
                        
                        for list1, list2 in zip(total_v_list, daily_v_list): 
                                v1.append(list1[key])
                                v2.append(list2[key])
                                
                        if np.shape(v1) == (len(storm_ranges), len(time_list)):
                            coord1, coord2 = 'nstorms','time'
                        else:
                            print('* using wind coords: ', np.shape(v1), np.shape(v2), len(storm_ranges))
                            coords = {'wind_x': (['wind_x'], np.arange(1,np.shape(v1)[0]+1)),
                                      'wind_y':(['wind_y'], np.arange(1,np.shape(v1)[1]+1))}
                            coord1, coord2 = 'wind_x', 'wind_y'
                            
                        data_vars['total_v_'+short] = ([coord1, coord2],  v1, 
                                 {'units': 'm/s', 
                                  'long_name':'v-wind timeseries, total miz area, '+key})
                        data_vars['daily_v_'+short] = ([coord1, coord2],  v2, 
                                 {'units': 'm/s', 
                                  'long_name':'v-wind timeseries, original miz, '+key})
                        
                except Exception as ee:
                    print('')
                    print('*ERROR* WINDS... exporting variables to nc')
                    print(ee)
                    print(traceback.format_exc())
                    print('')
                    
            if run['air_temp']:
                try:
                    for kk, key in enumerate(keys):
                        short = short_keys[kk]
                        t1, t2 = [],[]
                        for list1, list2 in zip(daily_temp_list, total_temp_list): 
                                t1.append(list1[key]) 
                                t2.append(list2[key])
                            
                        if np.shape(t1) == (len(storm_ranges), len(time_list)):
                            data_vars['temp_miz_'+short] = (['nstorms','time'], t1, 
                                     {'long_name':'2m air temperature, original miz, '+key})
                            data_vars['temp_miz2_'+short] = (['nstorms','time'], t2, 
                                     {'long_name':'2m air temperature, total miz area, '+key})
                        else:
                            print('* using air coords: ', np.shape(t1), np.shape(t2), len(storm_ranges))
                            coords = {'air_x': (['air_x'], np.arange(1,np.shape(t1)[0]+1)),
                                      'air_y':(['air_y'], np.arange(1,np.shape(t1)[1]+1))
                                      }
                            data_vars['air_temp_miz_'+short] = (['air_x','air_y'], t1, 
                                     {'long_name':'2m air temperature, original miz, '+key})
                            data_vars['air_temp_miz2_'+short] = (['air_x','air_y'], t2, 
                                     {'long_name':'2m air temperature, total miz area, '+key})
                        
                        
                        
                except Exception as ee:
                    print('')
                    print('*ERROR* AIR TEMP... exporting variables to nc')
                    print(ee)
                    print(traceback.format_exc())
                    print('') 
            
            if run['sst']:
                try:
                    for kk, key in enumerate(keys):
                        short = short_keys[kk]
                        sst = []
                        
                        for list1 in sst_list:
                            sst.append(list1[key]) 
                        
                        if np.shape(sst) == (len(storm_ranges), len(time_list)):
                            data_vars['sst_'+short] = (['nstorms','time'], sst, 
                                     {'units': 'deg C', 
                                      'long_name':'sea surface temperature, original miz, '+key})
                        else:
                            print('* using sst coords: ', np.shape(sst), len(storm_ranges))
                            coords = {'sst_x': (['sst_x'], np.arange(1,np.shape(sst)[0]+1)),
                                      'sst_y':(['sst_y'], np.arange(1,np.shape(sst)[1]+1))
                                      }
                            data_vars['sst_'+short] = (['sst_x','sst_y'], sst, 
                                     {'long_name':'sea surface temperature, original miz, '+key})

                except Exception as ee:
                    print('')
                    print('*ERROR* SST... exporting variables to nc')
                    print(ee)
                    print(traceback.format_exc())
                    print('')
                
            if run['sst_daily']: 
                try:
                    for kk, key in enumerate(keys):
                        short = short_keys[kk]
                        sst1, sst2 = [],[]
                        for list1, list2 in zip(daily_sst_list, total_sst_list): 
                                sst1.append(list1[key]) 
                                sst2.append(list2[key])
                            
                        if np.shape(sst1) == (len(storm_ranges), len(time_list)):
                            data_vars['daily_sst_'+short] = (['nstorms','time'], sst1, 
                                     {'long_name':'daily observed sst, original miz, '+key})
                            data_vars['daily_sst2_'+short] = (['nstorms','time'], sst2, 
                                     {'long_name':'daily observed sst, total miz area, '+key})
                        else:
                            print('* using daily-sst coords: ', np.shape(sst1), np.shape(sst2), len(storm_ranges))
                            coords = {'dsst_x': (['dsst_x'], np.arange(1,np.shape(sst1)[0]+1)),
                                      'dsst_y':(['dsst_y'], np.arange(1,np.shape(sst1)[1]+1))
                                      }
                            data_vars['daily_sst_'+short] = (['dsst_x','dsst_y'], sst1, 
                                     {'long_name':'daily observed sst, original miz, '+key})
                            data_vars['daily_sst2_'+short] = (['dsst_x','dsst_y'], sst2, 
                                     {'long_name':'daily observed sst, total miz area, '+key})
                        
                        
                        
                except Exception as ee:
                    print('')
                    print('*ERROR* DAILY SST... exporting variables to nc')
                    print(ee)
                    print(traceback.format_exc())
                    print('') 
                
            if run['glorys']:
                import pickle
                try:
                    for kk, key in enumerate(keys):
                        short = short_keys[kk]
                        tt1, tt2 = [],[]
                        for list1, list2 in zip(total_theta_list, daily_theta_list): 
                                tt1.append(list1[key]) 
                                tt2.append(list2[key])
                            
                        try:
                            with open(savepath+str(year)+'_ocean_'+key+'.pkl', 'wb') as f:
                                pickle.dump([tt1,tt2],f)
                            print('pickled: '+key)
                        except:
                            print('pickle error')
                            
                        if np.shape(tt1) == (len(storm_ranges), len(time_list)):
                            data_vars['ocean_temp_miz_'+short] = (['nstorms','time'], tt1, 
                                     {'long_name':'volume ocean temperature, original miz, '+key})
                            data_vars['ocean_temp_miz2_'+short] = (['nstorms','time'], tt2, 
                                     {'long_name':'volume ocean temperature, total miz area, '+key})
                        else:
                            print('* using ocean coords: ', np.shape(tt1), np.shape(tt2), len(storm_ranges))
                            coords = {'ocean_x': (['ocean_x'], np.arange(1,np.shape(tt1)[0]+1)),
                                      'ocean_y':(['ocean_y'], np.arange(1,np.shape(tt1)[1]+1))
                                      }
                            data_vars['ocean_temp_miz_'+short] = (['ocean_x','ocean_y'], tt1, 
                                     {'long_name':'volume ocean temperature, original miz, '+key})
                            data_vars['ocean_temp_miz2_'+short] = (['ocean_x','ocean_y'], tt2, 
                                     {'long_name':'volume ocean temperature, total miz area, '+key})
                        
                except Exception as ee:
                    print('')
                    print('*ERROR* OCEAN TEMP... exporting variables to nc')
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
        if run['winds']: print('winds, ', np.shape(wind1), np.shape(u1), np.shape(v1))
        if run['air_temp']: print('air temp, ', np.shape(t1), np.shape(t2))
        if run['sst']: print('sst, ',np.shape(sst))
        if run['glorys']: print('ocean, ', np.shape(tt1), np.shape(tt2))
        
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
        
