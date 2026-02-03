#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Updated 29 July 2022
myfuncs.py

@author: mundi
"""

#%% imports

import glob
import traceback
import netCDF4

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.pyplot import cm

from datetime import datetime, timedelta, date
import calendar

import sys
import os


si_levels = [0.15,0.80]

if True: # scientific format axis
    import matplotlib.ticker
    class OOMFormatter(matplotlib.ticker.ScalarFormatter):
        def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
            self.oom = order
            self.fformat = fformat
            matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
        def _set_order_of_magnitude(self):
            self.orderOfMagnitude = self.oom
        def _set_format(self, vmin=None, vmax=None):
            self.format = self.fformat
            if self._useMathText:
                self.format = r'$\mathdefault{%s}$' % self.format
                
 
#%% run function withou print output
# with HidePrint():
#   out = fxn.function(input)
     
class HidePrint:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


#%% read data files
def load_netcdf(filepath, in_vars):
    """open netcdf file, load variables from list in_vars and output dictionary of variables"""

    out_vars = {}

    open_netcdf = netCDF4.Dataset(filepath, mode = 'r')
    #print open_netcdf
    for var in in_vars:
        out_vars[var] = open_netcdf.variables[var][:]
    open_netcdf.close()

    return out_vars

#%% date time

def leading_zeros(day):
    if day>=10:
        return str(day)
    elif day<10:
        return '0'+str(day)

def daterange(start_date, end_date, dt=6):
    alldates=[]
    delta = timedelta(hours=dt)
    while start_date <= end_date:
        alldates.append(start_date)
        start_date += delta
    return alldates

#%% simple plotting

def setup_plot(ax, extent=[-160,90,50,60], title=[], labels=True):
    ax.coastlines('50m',edgecolor='black',linewidth=0.75)
    ax.set_extent(extent, ccrs.PlateCarree())
    try:
        ax.gridlines(draw_labels=labels)
    except:
        print('Unable to create grid lines on map')
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    ax.add_feature(cfeature.LAKES, facecolor='0.85')
    if title:
        if type(title)!= str: title=str(title)
        ax.set_title(title)
    return ax
    

def background_plot(extent=[-160,90,50,60], returnfig=False, title=[], labels=True, central_lon=-45):
    # alt extent: [-50,90,60,85]
    
    fig=plt.figure(figsize=[15,15]) 
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=central_lon))
    ax.coastlines('50m',edgecolor='black',linewidth=0.75)
    ax.set_extent(extent, ccrs.PlateCarree())
    try:
        ax.gridlines(draw_labels=labels)
    except:
        print('Unable to create grid lines on map')
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    ax.add_feature(cfeature.LAKES, facecolor='0.85')
    
    if title:
        if type(title)!= str: title=str(title)
        ax.set_title(title, fontsize=22)
    
    if not returnfig: return ax
    if returnfig: return fig, ax
    
def map_plot(X,Y,Z, ax=None):
    if ax == None: ax = background_plot()
    ax.pcolormesh(X,Y,Z, transform=ccrs.PlateCarree())
    return ax
    
def legend_without_duplicate_labels(ax, loc='best',fontsize=20):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    
    return ax.legend(*zip(*unique), loc=loc, fontsize=fontsize)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib.colors as colors
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def geoplot_2d(x,y,z=None):    
    #do masked-array on the lon
    x_greater = np.ma.masked_greater(x, -0.01)
    x_lesser = np.ma.masked_less(x, 0)    
    # apply masks to other associate arrays: lat
    y_greater = np.ma.MaskedArray(y, mask=x_greater.mask)
    y_lesser = np.ma.MaskedArray(y, mask=x_lesser.mask)
    
    if z is None:
        return [x_greater, x_lesser], [y_greater, y_lesser]
    else:
        # apply masks to other associate arrays
        z_greater = np.ma.MaskedArray(z, mask=x_greater.mask)
        z_lesser = np.ma.MaskedArray(z, mask=x_lesser.mask)
        return [x_greater, x_lesser], [y_greater, y_lesser], [z_greater, z_lesser]
    
    
def plot_bbox(bbox, ax=[], lw=2.5, color='k'):
    
    if ax==[]:
        ax = background_plot()
        
    bb_lon = bbox[:,0]
    bb_lat = bbox[:,1]
    
    [x_greater, x_lesser], [y_greater, y_lesser] = geoplot_2d(bb_lon, bb_lat)
    
    ax.plot(x_greater, y_greater, lw=lw, color=color, transform=ccrs.PlateCarree())
    ax.plot(x_lesser, y_lesser, lw=lw, color=color, transform=ccrs.PlateCarree())
    
    return ax

#%% simple calculations

def find_points_in_contour(coords, var_x, var_y, var=None):
    import matplotlib.path as mpltPath
    
    coord_x = np.where(coords[:,0] < 0, coords[:,0]+360, coords[:,0])
    coord_y = coords[:,1]
    
    ### turn contour into polygon
    polygon = np.vstack((coord_x, coord_y)).T
    path = mpltPath.Path(polygon)
    
    points = np.vstack((var_x.flatten(), var_y.flatten()))
    points=points.T
    
    ### get inside points and plot
    inside = path.contains_points(points)
    inside = np.reshape(inside, np.shape(var_x))
    
    if np.nanmin(var_x) < 0:
        var_x2 = np.where(var_x < 0, var_x+360, var_x)
    else:
        var_x2 = np.where(var_x > 180, var_x-360, var_x)
    points2 = np.vstack((var_x2.flatten(), var_y.flatten())).T
    inside2  = path.contains_points(points2)
    inside2 = np.reshape(inside2, np.shape(var_x))
    
    inside_points = np.invert(np.logical_or(inside, inside2))

    return inside_points


def calc_deriv(series, dt=3):
    
    deriv_series = [] #np.nan, np.nan
    for idx, ts in enumerate(series[:-dt]):
        ds = series[idx+dt] - series[idx]
        deriv_series.append(ds/dt)
        
    deriv_series.append(np.nan)
    
    return deriv_series




#%% era 5 download
# '10m_u_component_of_wind'
#u10

def get_era5(dataset_name='reanalysis-era5-single-levels', 
             var=None, 
             dates=None,
             pressure_level=None,
             grid=[0.25, 0.25],
             area=[90, -180, -90, 180],
             download_flag = False,
             download_file='./output.nc'
            ):
    raise NameError('New CDS API - use daily/hourly versions instead')
    ''' Get ERA5 reanalysis output from the web
    this script grabs ERA5 variables from the web and stores them 
    in an xarray dataset. 
    
    the ERA5 CDS API must be installed on the local machine.
    See section 4 here: https://cds.climate.copernicus.eu/api-how-to
    
    Parameters
    ----------                  
    dataset_name: str, default 'reanalysis-era5-single-levels'
        name of dataset to use. Options include:
        * 'reanalysis-era5-single-levels'
        * 'reanalysis-era5-single-levels-monthly-means'
        * 'reanalysis-era5-pressure-levels'
        * 'reanalysis-era5-pressure-levels-monthly-means'
        * 'reanalysis-era5-land'
        * 'reanalysis-era5-land-monthly-means'
        
    dates: list of strings or datetime64, default None
        example ['1980-01-01', '2020-12-31']
    var: str, default None
        name of variable to download
        example '2m_temperature'
    pressure_level: str, default None
        pressure level to grab data on
    grid: list, deafult [1.0, 1.0]
        spatial lat, lon grid resolution in deg
    area: list, default [90,-180,-90, 180]
        area extent download [N, W, S, E]
    download_flag = True or False, default False
        flag to download data or not
    download_file= str, default './output.nc'
        path to where data should be downloaed to.
        data only downloaded if download_flag is True
    Returns
    -------
    ds: xarrayDataSet
        all the data will be in an xarray dataset
        
    Example
    -------
    ds = get_era5(dataset_name='reanalysis-era5-single-levels-monthly-means', 
                 var='2m_temperature', 
                 dates=['2021-02-01'],
                 grid=[0.25, 0.25])
        
    Notes
    -------    
    # cdsapi code is here
    https://github.com/ecmwf/cdsapi/tree/master/cdsapi
    # information on api is here
    https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+%28CDS%29+API+Keywords
    # era5 dataset information is here
    https://confluence.ecmwf.int/display/CKB/The+family+of+ERA5+datasets
    '''
    import cdsapi
    import xarray as xr
    import pandas as pd
    from urllib.request import urlopen
    
    # test if acceptable pressure level
    acceptable_pressures = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, range(100, 1000, 25)]
    if pressure_level not in [str(lev) for lev in acceptable_pressures]:
        print(f"!! Pressure level must be in this list: {acceptable_pressures}")
    
    # start the cdsapi client
    c = cdsapi.Client()

    # parameters
    params = dict(
        format = "netcdf",
        product_type = "reanalysis",
        variable = var,
        grid = grid,
        area = area,
        date = list(dates.strftime('%Y-%m-%d %H:%M')) \
               if isinstance(dates, pd.core.indexes.datetimes.DatetimeIndex)\
               else dates,
        )

    # what to do if asking for monthly means
    if dataset_name in ["reanalysis-era5-single-levels-monthly-means", 
                        "reanalysis-era5-pressure-levels-monthly-means",
                        "reanalysis-era5-land-monthly-means"]:
        params["product_type"] = "monthly_averaged_reanalysis"
        _ = params.pop("date")
        params["time"] = "00:00"
        
        # if time is in list of pandas format
        if isinstance(dates, list):
            dates_pd = pd.to_datetime(dates)
            params["year"] = sorted(list(set(dates_pd.strftime("%Y"))))
            params["month"] = sorted(list(set(dates_pd.strftime("%m"))))
        else:
            params["year"] = sorted(list(set(dates.strftime("%Y"))))
            params["month"] = sorted(list(set(dates.strftime("%m"))))
            
        
    # if pressure surface
    if dataset_name in ["reanalysis-era5-pressure-levels-monthly-means",
                        "reanalysis-era5-pressure-levels"]:
        params["pressure_level"] = pressure_level
    
    # product_type not needed for era5_land
    if dataset_name in ["reanalysis-era5-land"]:
        _ = params.pop("product_type")
        
    # file object
    fl=c.retrieve(dataset_name, params) 
    
    # download the file 
    if download_flag:
        fl.download(f"{download_file}")
    
    # load into memory and return xarray dataset
    with urlopen(fl.location) as f:
        return xr.open_dataset(f.read())


def era5_daily(year, month, day, variable = 'mean_sea_level_pressure'):
    '''
    Parameters
    ----------
    year : list
       ex. ['2010']
    month : list
        ex. ['08']
    day : lsit
        ex. ['15','16']

    Returns
    -------
    daily_time : pandas datetime
        daily dates
    daily_avg_slp : list
        lsit of daily averaged slp

    '''
    if type(year) != list:
        year = [str(year)]
    if type(month) != list:
        if month < 10: monthstr = '0'+str(int(month))
        else: monthstr = str(int(month))
        month = [monthstr]
    if type(day) != list:
        if day < 10: daystr = '0'+str(int(day))
        else: daystr = str(int(day))
        day = [daystr]

    import pandas as pd
    import cdsapi
    import xarray as xr
    import io
    from urllib.request import urlopen# start the client
    cds = cdsapi.Client()# dataset you want to read
    dataset = "reanalysis-era5-single-levels"# flag to download data
    download_flag = False# api parameters 
    params = {
        "format": "netcdf",
        "product_type": 'reanalysis',
        "variable": variable,
        'year':year,
        'month':month,
        'day':day,
        "time": ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',
                 '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00',
                 '16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'],
        "grid":[0.25,0.25],
        "area":[90, -180, 60, 180]
        }
    # retrieves the path to the file
    fl = cds.retrieve(dataset, params)
    # download the file 
    if download_flag:
        fl.download("./output.nc")
    # load into memory
    with urlopen(fl.location) as f:
        bytes_ = f.read()
        ds = xr.open_dataset(io.BytesIO(bytes_), decode_times=False)
        
    lon, lat = np.meshgrid(ds['longitude'], ds['latitude'])
    
    if variable == 'msl' or variable =='mean_sea_level_pressure':
        key = 'msl'
    elif variable =='10m_v_component_of_wind':
        key = 'v10'
    elif variable =='10m_u_component_of_wind':
        key = 'u10'
    elif variable == '2m_temperature':
        key='t2m'
    else:
        key=variable
        

    # time = ds['time'].values
    # # ts = pd.to_datetime(time)
    # ts = [tt-time[0] for tt in time]
    
    time_in = ds['valid_time'].values
    ts = [(tt-time_in[0])/60/60 for tt in time_in] # seconds -> hours
    # starting_day = datetime(int(year[0]), int(month[0]), int(day[0]))
    # time = np.array([starting_day + timedelta(hours=int(hr)) for hr in hours_in])
    
    try:
        var = ds[key].values
    except:
        print(key)
        raise ValueError('ERA5_daily: invalid variable key')
    
    if variable == 'msl' or variable=='mean_sea_level_pressure':
        var = var/100

    start_day = int(day[0])
    my_pres = []
    daily_avg = []
    daily_time = []
    
    for tx, tt in enumerate(ts):
        # if tt.day == start_day:
        if tt != 23:
            my_pres.append(var[tx,:,:])
        else: # save this day and move to next
            my_pres.append(var[tx,:,:]) # append last hour
            daily_avg.append( np.mean(my_pres, axis=0) )
            daily_time.append(datetime(int(year[0]), int(month[0]), start_day, int(ts[0])))
            
            my_pres=[]
            start_day += 1
            
        # if tt==ts[-1]:
        #     daily_avg.append( np.mean(my_pres, axis=0) )
            
    if len(day) == 1:
        # daily_avg.append( np.mean(my_pres, axis=0) )
        daily_avg=np.squeeze(daily_avg)        

    return daily_time, lon, lat, daily_avg


def era5_hourly(year, month, day, variable = 'mean_sea_level_pressure'):
    '''
    Parameters
    ----------
    year : list
       ex. ['2010']
    month : list
        ex. ['08']
    day : list
        ex. ['15','16']

    Returns
    -------
    daily_time : pandas datetime
        daily dates
    daily_avg_slp : list
        lsit of hourly slp

    '''
    if type(year) != list:
        year = [str(year)]
    if type(month) != list:
        if month < 10: monthstr = '0'+str(int(month))
        else: monthstr = str(int(month))
        month = [monthstr]
    if type(day) != list:
        if day < 10: daystr = '0'+str(int(day))
        else: daystr = str(int(day))
        day = [daystr]

    import pandas as pd
    import cdsapi
    import xarray as xr
    import io
    from urllib.request import urlopen# start the client
    cds = cdsapi.Client()# dataset you want to read
    dataset = "reanalysis-era5-single-levels"# flag to download data
    download_flag = False# api parameters 
    params = {
        "product_type": 'reanalysis',
        "variable": variable,
        'year':year,
        'month':month,
        'day':day,
        "time": ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',
                 '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00',
                 '16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'],
        # "grid":[0.25,0.25],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area":[90, -180, 60, 180]
        }
    # retrieves the path to the file# retrieves the path to the file
    fl = cds.retrieve(dataset, params)
    # download the file 
    if download_flag:
        fl.download("./output.nc")
    # load into memory
    with urlopen(fl.location) as f:
        bytes_ = f.read()
        ds = xr.open_dataset(io.BytesIO(bytes_), decode_times=False)
        
    lon, lat = np.meshgrid(ds['longitude'], ds['latitude'])
    
    if variable == 'mean_sea_level_pressure' or variable=='msl':
        key = 'msl'
    elif variable =='10m_v_component_of_wind':
        key = 'v10'
    elif variable =='10m_u_component_of_wind':
        key = 'u10'
    else:
        raise ValueError('ERA5_hourly: invalid variable key')

    time_in = ds['valid_time'].values
    hours_in = [(tt-time_in[0])/60/60 for tt in time_in] # seconds -> hours
    starting_day = datetime(int(year[0]), int(month[0]), int(day[0]))
    time = np.array([starting_day + timedelta(hours=int(hr)) for hr in hours_in])
    
    # time = pd.to_datetime(time_in)
    
    var = ds[key].values
    
    if variable == 'msl' or key=='msl':
        var = var/100

    return time, lon, lat, var

#%% Pressure Contours
# ------------------------
#%%% plot

def plotConstPressure(all_contours, dateDT, levels=[990,1000], 
                      loc='best', fig=[], ax=[], fontsize=20,
                      colorbar = cm.cool, lw=3):
    
    if ax==[] or fig==[]:
        fig, ax = background_plot(returnfig=True)
    
    ### double colormap for plotting 1000 and 990 mb
    # colors = colorbar( np.linspace(0, 1, int((len(all_contours)+1)/2) ))
    # doublecol=[]
    # for col in colors:
    #     doublecol.append(col); doublecol.append(col)
    # colors = iter(np.array(doublecol))
    
    dateDT_backup = dateDT
    if len(dateDT_backup) == len(all_contours)/2:
        dateDT = []
        for dt in dateDT_backup:
            dateDT.append(dt); dateDT.append(dt)
            
        print('backup: ' + str(len(dateDT)))
        print(dateDT)
        
    
    ### new colormap iter
    colorlist = colorbar( np.linspace(0, 1, int(np.ceil((len(all_contours)+1)/2)) ))
    colors = iter( colorlist )
    mycol = next(colors)

    ### plot contour 
    for dd, dxy in enumerate(all_contours):       
        if len(dxy)>3: 
            # set up params
            lon_greater = np.ma.masked_greater(dxy[:,0], -0.01)
            lon_lesser = np.ma.masked_less(dxy[:,0], 0)    
            lat_greater = np.ma.MaskedArray(dxy[:,1], mask=lon_greater.mask)
            lat_lesser = np.ma.MaskedArray(dxy[:,1], mask=lon_lesser.mask)
            
            # get daily color
            if dd !=0:
                try:
                    if dateDT[dd].day != dateDT[dd-1].day:
                        try:
                            mycol = next(colors)
                            if len(levels) == 2: levstr = str(levels[-1])
                                
                        except StopIteration:
                            colors = iter( colorbar( np.linspace(0, 1, int(np.ceil((len(all_contours)+1)/2)) )) )
                            mycol = next(colors)
                    else:
                        if len(levels) == 2: levstr = str(levels[0])
                except:
                    print('- error with colors in contour plotting -')
                    print('index: '+ str(dd) + ', dateDT length: ' + str(len(dateDT)) + ', all_contours length: ' + str(len(all_contours)) )
                    mycol = colorlist[-1]
                    
            else: #first timestep
                if len(levels) == 2: levstr = str(levels[-1])
                    
            # if len(levels) == 2:
            #     levstr = str(levels[dd%2])
            #     if dd == len(all_contours)-1:
            #         levstr = str(levels[-1])
                    
            if len(levels) != 2:
                levstr = str(levels[dd])
                
            
            
            # plot (pos/neg)
            try:
                lab = dateDT[dd].strftime("%b %d")
            except:
                lab = '*' + dateDT[0].strftime("%b %d")
            ax.plot(lon_greater, lat_greater, transform=ccrs.PlateCarree(), linewidth=lw, 
                label = lab +', '+ levstr, color=mycol)  
            ax.plot(lon_lesser, lat_lesser, transform=ccrs.PlateCarree(), linewidth=lw, 
                label = lab +', '+ levstr, color=mycol)  
                 
    ### finish plot
    thislegend = legend_without_duplicate_labels(ax, loc=loc, fontsize=fontsize)
    ax.add_artist(thislegend)
    
    return fig, ax

#%%% manual detection
def load_pressure(file_name, year, month, day, daily=True):
    from netCDF4 import Dataset
    year=int(year) ; month=int(month); day=int(day)
    
    wd = Dataset(file_name)
    time = wd['time'][:]
    lat = wd['latitude'][:]
    lon = wd['longitude'][:]
        
    x, y = np.meshgrid(lon, lat)    
    start = datetime(1900,1,1)
    
    ### Get correct dates
    start = datetime(1900,1,1)
    inds, mytime = [], []
    for i, t in enumerate(time):
        hrs = int(t)
        delta = timedelta(hours=hrs)
        offset = start + delta   
        
        if offset.year==year and offset.month == month and offset.day==day:
            inds.append(i)
            mytime.append(offset)  
     
    if inds != []:
        pressure = wd['msl'][inds[0]:inds[-1]+1]
        pressure = pressure/100 # convert to hPa, mb
        
        ### Average pressure
        if daily:
            mytime=mytime[0]
            pressure = pressure.mean(axis=0)        
    else:
        print('error with loading pressure')
        pressure, mytime = [], []
        
    wd.close()
    pressure=np.array(pressure)
    return x, y, pressure, mytime

def get_contour_points(year, month, day, level, storm_info, local_file ='', lon_thresh=10, lat_thresh=5):
    raise ValueError('get_contour_points: switch to new verison')
    
    if type(level) != list: level=[level]
    
    print(''); print(year, month, day); print('---------')
    
    if local_file == '':
        if type(year) != list: year=[str(year)]
        if type(month)!= list: month=[str(month)]
        if type(day) != list: day = [str(day)]
        time, x,y, pressure = era5_daily(year, month, day)
        pressure=pressure[0]
    else:
        x, y, pressure, time = load_pressure(local_file, year, month, day, daily=True)
        if pressure == []:
            raise NameError('get_contour_points: empty pressure ' + str((year, month, day)))
        
    x = np.where(x<0, x+360, x) # make all lons positive for comparison
    
    ### find minimum pressures and compare to storm_info ###!!!!
    storm_x = storm_info[0]
    storm_y = storm_info[1]
    nearby_min = False
    while nearby_min==False:
        # print('new search')
        indp = np.where(pressure == np.nanmin(pressure))
        minlon = x[indp]
        minlat = y[indp]
       
        for pt_id, lo in enumerate(minlon):
            
            # print('---------------------------')
            # print(storm_x, storm_y)
            # print(minlon[pt_id], minlat[pt_id])
            # print(abs(minlon[pt_id] - storm_x)<lon_thresh, abs(minlat[pt_id] - storm_y) < lat_thresh)
            
            if abs(minlon[pt_id] - storm_x) < lon_thresh and abs(minlat[pt_id] - storm_y) < lat_thresh:
                min_x = float(minlon[0])
                min_y = float(minlat[0])
                min_p = np.nanmin(pressure)
                nearby_min = True
                break
            else:
                # print('clear')
                pressure[indp] = np.nan #[pt_id]
                break
                
      
    x = np.where(x>180, x-360, x)   ###!!!
    x1 = x
    x1 = np.where(x1<-179, np.nan, x1)
    x1 = np.where(np.logical_and(x1>179, x1<181), np.nan, x1)
    x1 = np.where(np.logical_and(x1<=1, x1>-1), np.nan, x1)

    # x1 = np.where(x1<180, x1+360, x1)
    # x1 = np.where(x1<1, np.nan, x1)
    
      
    #do masked-array on the lon
    x_greater = np.ma.masked_greater(x1, -0.01)
    x_lesser = np.ma.masked_less(x1, 0)    
    # apply masks to other associate arrays: lat
    y_greater = np.ma.MaskedArray(y, mask=x_greater.mask)
    y_lesser = np.ma.MaskedArray(y, mask=x_lesser.mask)
    # apply masks to other associate arrays: daily ice
    p_greater = np.ma.MaskedArray(pressure, mask=x_greater.mask)
    p_lesser = np.ma.MaskedArray(pressure, mask=x_lesser.mask)

    fig, ax =  background_plot(returnfig=True)
    ax.set_title(str(year)+' '+str(month)+' '+str(day) + ': ' +str(round(min_p,1))+' hPa' +',  '+str(round(min_x,2)), 
                 fontsize=26)
    # contours1 = ax.contour(x_greater, y_greater, p_greater, levels=level,
    #                   linewidths = 1.5, zorder=10, colors='k',
    #                   transform=ccrs.PlateCarree())
    # contours2 = ax.contour(x_lesser, y_lesser, p_lesser, levels=level,
    #                   linewidths = 1.5, zorder=10, colors='k',
    #                   transform=ccrs.PlateCarree())
    
    contours_all = ax.contour(x1, y, pressure, levels=level,
                      linewidths = 1.5, zorder=10, colors='k',
                      transform=ccrs.PlateCarree())
    
    
    ax.plot(min_x, min_y,'y*',transform=ccrs.PlateCarree(), markersize=40)    
    # plt.close(fig)
    
    coords=[]
    for lvl, lev in enumerate(level):
        # try:
        #     coords = coords + contours1.allsegs[lvl]
        # except IndexError:
        #     print(str(lev) +': ' + str(month)+', '+ str(day))
        #     pass
        # try:
        #     coords = coords + contours2.allsegs[lvl]
        # except IndexError:
        #     print(str(lev) +': ' + str(month)+', '+ str(day))
        #     pass
        try:
            coords = coords + contours_all.allsegs[lvl]
        except IndexError:
            print(str(lev) +': ' + str(month)+', '+ str(day))
        
    min_x = np.squeeze(min_x)
    if min_x > 180:
        min_x = min_x - 360
    
    return coords, min_x, np.squeeze(min_y), ax

def get_contours(storm_daymonth_grid, daymonth_grid, year, intervals=[990,1000], storm_info=[]):
    ''' have user select correct contours '''
    
    
    # return lists of selcted contours and their respective date
    all_cont_dt1, all_conts1 = [], []
    

    for month, day in storm_daymonth_grid:
        ### get pressure contours      
        intervals = [990,1000] #[980,985,990,995,1000,1005,1010] #
        current_date = datetime(int(year), int(month), int(day))
        
        coords = get_contour_points(year, month, day, intervals) 
        
        ### plot
        fig, ax = background_plot(returnfig=True)
        for txt, coord in enumerate(coords):
            ax.plot(coord[:,0], coord[:,1], transform=ccrs.PlateCarree())
            transform = ccrs.PlateCarree()._as_mpl_transform(ax)
            ax.annotate(str(txt+1), xy=(coord[0,0], coord[0,1]), xycoords=transform, fontsize=26)
        ax.set_title(current_date.strftime('%m %d %Y'), fontsize=22)
        
        
        if storm_info != []:
            
            ax.plot([storm_info[0][0],storm_info[-1][1]],[storm_info[0][1],storm_info[-1][1]], color='g',
                        linestyle = 'dashed', linewidth=2.5,
                        transform=ccrs.Geodetic())
            # plot start/end markers
            ax.plot(storm_info[0][0],storm_info[0][1], 'o--',  color='g', 
                        transform=ccrs.PlateCarree(), markersize=12)
            ax.plot(storm_info[-1][0],storm_info[-1][1], 's--', color='g', 
                        transform=ccrs.PlateCarree(), markersize=12)
            
        plt.show()
        
        ### ask user to pick correct contour for each date
        for c_str in ['990', '1000']:
            valid_ans=True
            while valid_ans:
                print('Select ' + c_str + ' contour of this storm day (or type append):')
                storm_min = input()
                
                try:
                    storm_min = int(storm_min)
                except ValueError:
                    if storm_min=='append': 
                        
                        all_in = 'start'
                        this_cont = np.array([[np.nan,np.nan]])
                        while all_in != 'end':
                            print('Enter all contour pieces or end to stop')
                            all_in = input()
                            try:
                                int(all_in)
                            except ValueError:
                                if all_in == 'end':
                                    continue
                            this_cont = np.concatenate((this_cont, coords[int(all_in)-1]))
                        all_conts1.append(this_cont)
                            
                        break
                    else:
                        print('Not a valid integer')
                        continue
                
                if int(storm_min) <= len(coords):
                    valid_ans = False
                    # get points of chosen contour
                    # append these coords
                    all_conts1.append( np.array(coords[int(storm_min)-1]) )
                    # all_cont_dt1.append(current_date)
                else: 
                    valid_ans = True
                    print('Try again')
        
    return all_cont_dt1, all_conts1   

#%%% automatic detection

def find_contour_points(x, y, slp, level=[990, 1000]):
   
    fig, ax = background_plot(returnfig=True)
    contours = ax.contour(x, y, slp, levels=level,
                      linewidths = 1.5, zorder=10, colors='k',
                      transform=ccrs.PlateCarree())
    plt.close(fig)
    
    coords=[]
    for lvl, lev in enumerate(level):
        coords = coords + contours.allsegs[lvl]
        
    return coords

def detect_contours(storm_daymonth_grid, daymonth_grid, year, storm_info, intervals=[990,1000], local_file=''):
    raise ValueError('detect_contours: switch to new verison')
    
    import matplotlib.path as mpltPath
    # return lists of selcted contours and their respective date
    all_cont_dt1, all_conts1 = [], []
    
    snum=-1
    for month, day in storm_daymonth_grid:
        snum+=1
        ### get pressure contours      
        intervals = [990,1000] #[980,985,990,995,1000,1005,1010] #
        current_date = datetime(int(year), int(month), int(day))
        
        coords, min_x, min_y, ax = get_contour_points(year, month, day, intervals, storm_info[snum], local_file=local_file)
        current_title = ax.get_title()
        
        points = np.array([ [min_x, min_y] ])
        for cc1 in coords :
            conts_detected = False
            
            cc2 = np.array(cc1)
            clon = cc2[:,0]; clat=cc2[:,1]
            clon= np.where(clon>180, clon-360, clon)
            cc = np.vstack((clon, clat)).T
            
            if len(cc) > 2 : path = mpltPath.Path(cc)
            else: continue
            pts_inside = path.contains_points(points) # radius=30
            
            ax.set_title(current_title +'....' + str(min_x), fontsize=28)
            
            if any(pts_inside):
                if abs( abs(np.nanmean(cc[:,0])) - abs(min_x) ) < 90:
                    all_conts1.append(cc)
                    all_cont_dt1.append(current_date)
                    conts_detected = True
                    # ax.plot(cc[:,0], cc[:,1], 'r', linewidth=3, transform=ccrs.PlateCarree())
                    
                    [x1,x2], [y1,y2] = geoplot_2d(cc[:,0], cc[:,1])
                    ax.plot(x1, y1, 'r', linewidth=5, transform=ccrs.PlateCarree())
                    ax.plot(x2, y2, 'r', linewidth=5, transform=ccrs.PlateCarree())
                    continue
                    
            else:
                if False: #check if this contour contains other contours
                    for cc2 in all_conts1:
                        try:
                            path2 = mpltPath.Path(cc2)
                            pts_inside2 = path.contains_path(path2)
                            pts_inside3 = path2.contains_path(path)
                            if pts_inside2 or pts_inside3:
                                all_conts1.append(cc)
                                all_cont_dt1.append(current_date)
                                
                                [x1,x2], [y1,y2] = geoplot_2d(cc[:,0], cc[:,1]) ###!!!
                                ax.plot(x1, y1, 'b', linewidth=5, transform=ccrs.PlateCarree())
                                ax.plot(x2, y2, 'b', linewidth=5, transform=ccrs.PlateCarree())
                        except IndexError:
                            print('Index Error: ' + current_date.strftime('%Y-%m-%d'))
                            continue
    
            if abs(abs(min_x) - 180) < 20:
                second_x = 177
                second_pts = path.contains_points(np.array([ [second_x, min_y] ]))
                ax.plot(second_x, min_y,'c*',transform=ccrs.PlateCarree(), markersize=40) 
                if any(second_pts):
                    all_conts1.append(cc)
                    all_cont_dt1.append(current_date)
                    conts_detected = True
                    [x1,x2], [y1,y2] = geoplot_2d(cc[:,0], cc[:,1])
                    ax.plot(x1, y1, 'b', linewidth=5, transform=ccrs.PlateCarree())
                    ax.plot(x2, y2, 'b', linewidth=5, transform=ccrs.PlateCarree())
                continue
                    
            if (min_x >= 0 and min_x<10) or (min_x<-170) or (min_x>350) or (min_x<=0 and min_x>-10):
                lefty = -10
                righty = 10
                third_pts = path.contains_points(np.array([ [righty, min_y] ]))
                fourth_pts = path.contains_points(np.array([ [lefty, min_y] ]))
                ax.plot(righty, min_y,'g^',transform=ccrs.PlateCarree(), markersize=30)
                ax.plot(lefty, min_y,'m^',transform=ccrs.PlateCarree(), markersize=30)
                if any(third_pts) or any(fourth_pts):
                    all_conts1.append(cc)
                    all_cont_dt1.append(current_date)
                    conts_detected = True
                    [x1,x2], [y1,y2] = geoplot_2d(cc[:,0], cc[:,1])
                    ax.plot(x1, y1, 'y', linewidth=5, transform=ccrs.PlateCarree())
                    ax.plot(x2, y2, 'y', linewidth=5, transform=ccrs.PlateCarree())
                else:
                    third_pts = path.contains_points(np.array([ [righty, min_y+3.5] ]))
                    fourth_pts = path.contains_points(np.array([ [lefty, min_y+3.5] ]))
                    ax.plot(righty, min_y,'g^',transform=ccrs.PlateCarree(), markersize=30)
                    ax.plot(lefty, min_y,'m^',transform=ccrs.PlateCarree(), markersize=30)
                    if any(third_pts) or any(fourth_pts):
                        all_conts1.append(cc)
                        all_cont_dt1.append(current_date)
                        conts_detected = True
                        [x1,x2], [y1,y2] = geoplot_2d(cc[:,0], cc[:,1])
                        ax.plot(x1, y1, 'm', linewidth=5, transform=ccrs.PlateCarree())
                        ax.plot(x2, y2, 'm', linewidth=5, transform=ccrs.PlateCarree())
                
            if not conts_detected:
                
                final_pts = path.contains_points(np.array([ [min_x, min_y+5] ]))
                ax.plot(min_x, min_y+5,'b^',transform=ccrs.PlateCarree(), markersize=40) 
                if any(final_pts):
                    all_conts1.append(cc)
                    all_cont_dt1.append(current_date)
                    conts_detected = True
                    [x1,x2], [y1,y2] = geoplot_2d(cc[:,0], cc[:,1])
                    ax.plot(x1, y1, 'c', linewidth=5, transform=ccrs.PlateCarree())
                    ax.plot(x2, y2, 'c', linewidth=5, transform=ccrs.PlateCarree())
    
    return all_cont_dt1, all_conts1   


#%%% new try
from cartopy.util import add_cyclic_point

def get_contour_points2(year, month, day, level, storm_info, local_file ='', lon_thresh=10, lat_thresh=5):
    if type(level) != list: level=[level]
    
    print(''); print(year, month, day); print('---------')
    
    if local_file == '':
        if type(year) != list: year=[str(year)]
        if type(month)!= list: month=[str(month)]
        if type(day) != list: day = [str(day)]
        time, x,y, pressure = era5_daily(year, month, day)
        # pressure=pressure[0]
        pressure = np.squeeze(pressure)
    else:
        x, y, pressure, time = load_pressure(local_file, year, month, day, daily=True)
        if 0 in np.shape(pressure): 
            raise NameError('get_contour_points: empty pressure ' + str((year, month, day)))
        
    xv = x[0,:]
    yv = y[:,0]
    
    cyclic_data, cyclic_lons = add_cyclic_point(pressure, coord=xv)
    x, y = np.meshgrid(xv,yv)
    # x = np.where(x<0, x+360, x) # make all lons positive for comparison
    
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
    fig, ax =  background_plot(returnfig=True)
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

def detect_contours2(storm_daymonth_grid, daymonth_grid, year, storm_info, intervals=[990,1000], local_file=''):
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
        
        coords, min_x, min_y, ax = get_contour_points2(year, month, day, intervals, 
                                                       storm_info[snum], 
                                                       local_file=local_file)
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
                    [x1,x2], [y1,y2] = geoplot_2d(cc[:,0], cc[:,1])
                    ax.plot(x1, y1, 'b', linewidth=5, transform=ccrs.PlateCarree())
                    ax.plot(x2, y2, 'b', linewidth=5, transform=ccrs.PlateCarree())
                continue
    
    print()
    return all_cont_dt1, all_conts1   



#%% storm census
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

#%% bounding box - storm area

def open_cont_nc(ncfile):
    # loads pressure contours to form bbox
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

def open_cont_nc2(ncfile):
    ds = xr.open_dataset(ncfile)
    
    all_contours = []
    all_cont_dts = []
    for key in list(ds.keys()):
        all_contours.append(ds[key].values)
        all_cont_dts.append(datetime.strptime(key[7:16],'%Y_%m%d'))
        
    ds.close()   
    return all_contours, all_cont_dts
        
    

def get_edge_lines(minlon, maxlon, minlat, maxlat, n=90, reverse=False):
    ### create new bbox edges
    #
    edge1x = np.linspace(minlon, minlon, n)
    edge1y = np.linspace(minlat, maxlat, n)
    edge2x = np.linspace(minlon, maxlon, n)
    edge2y = np.linspace(maxlat, maxlat, n)
    edge3x = np.linspace(maxlon, maxlon, n)
    edge3y = np.linspace(maxlat, minlat, n)
    edge4x = np.linspace(maxlon, minlon, n)
    edge4y = np.linspace(minlat, minlat, n)
    if reverse:
        edge2x = np.concatenate((np.linspace(minlon,0,round(n/3)),np.linspace(0,180,round(n/3))))
        edge2x = np.concatenate( ( edge2x,np.linspace(-180,maxlon,round(n/3)) ) ) 
        #
        edge4x = np.concatenate((np.linspace(maxlon,-180,round(n/3)),np.linspace(180,0,round(n/3))))
        edge4x = np.concatenate( ( edge4x,np.linspace(0,minlon,round(n/3)) ) ) 
    #
    bbox_lon = np.append(edge1x, edge2x)
    bbox_lon = np.append(bbox_lon, edge3x)
    bbox_lon = np.append(bbox_lon, edge4x)
    #
    bbox_lat = np.append(edge1y, edge2y)
    bbox_lat = np.append(bbox_lat, edge3y)
    bbox_lat = np.append(bbox_lat, edge4y)
    #
    bbox_edges = np.squeeze(np.array([[bbox_lon],[bbox_lat]])).T
    
    return bbox_edges

def get_bbox_edges_old(all_contours):
    minlon, maxlon, minlat, maxlat = None, None, None, None
    
    ### get rid of boundaries too close to pole
    alerted=False
    for cidx, contour in enumerate(all_contours):
        lons = contour[:,0]
        lats = contour[:,1]
        
        for li, lat in enumerate(lats):
            if lat>85:
                lons[li]=np.nan
                lats[li]=np.nan
                if not alerted:
                    print('bbox too close too pole; artifical boundary applied (85N)')
                    alerted=True

        ### convert longitude to 0-360 system
        lons1 = lons #np.array(lons)
        lons360 = np.where(lons1<0, lons1+360, lons1)
        # print([np.nanmin(lons360), np.nanmax(lons360)])
        
        if cidx == 0:
            minlon = np.nanmin(lons360)
            maxlon = np.nanmax(lons360)
            minlat = np.nanmin(lats)
            maxlat = np.nanmax(lats)
        else:
            if np.nanmin(lons360) < minlon:
                minlon = np.nanmin(lons360)
            if np.nanmax(lons360) > maxlon:
                maxlon = np.nanmax(lons360)
            if np.nanmin(lats) < minlat:
                minlat = np.nanmin(lats)
            if np.nanmax(lats) > maxlat:
                maxlat = np.nanmax(lats)

    bbox_edges = get_edge_lines(minlon, maxlon, minlat, maxlat, reverse=False)
    print([minlon, maxlon, minlat, maxlat])
    
    if (round(minlon) < 3) and round(maxlon, -1) > 358:
        bbox_edges = get_bbox_edges_backup(all_contours)
    
    return bbox_edges

def get_bbox_edges_backup(all_contours):
    ### secondary bbox calculator for get_bbox_edges based on storm location
    
    import warnings
   
    pos_minlon = 9999
    neg_minlon = 9999
    pos_maxlon = -9999
    neg_maxlon = -9999
    
    minlat = 9999
    maxlat = -9999
    
    rev1=False
    
    alerted=False
    for contour in all_contours:
        lons = contour[:,0]
        lats = contour[:,1]
        
        for li, lat in enumerate(lats):
            if lat>85:
                lons[li]=np.nan
                lats[li]=np.nan
                if not alerted:
                    print('bbox too close too pole; artifical boundary applied (85N)')
                    alerted=True
        
        ### break lons up into pos and neg
        lons_neg = np.where(lons>=0, np.nan, lons)
        lons_pos = np.where(lons<0, np.nan, lons)
        
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore", category=RuntimeWarning) # mean of empty slice
            # print('neg: ' + str([np.nanmin(lons_neg),np.nanmax(lons_neg)]))
        
        if any(~np.isnan(lons_neg)):
            if np.nanmin(lons_neg) < neg_minlon:
                neg_minlon = np.nanmin(lons_neg)
            if np.nanmax(lons_neg) > neg_maxlon:
                neg_maxlon = np.nanmax(lons_neg)
                
        if any(~np.isnan(lons_pos)):
            if np.nanmin(lons_pos) < pos_minlon:
                pos_minlon = np.nanmin(lons_pos)
            if np.nanmax(lons_pos) > pos_maxlon:
                pos_maxlon = np.nanmax(lons_pos)
            
        ### latitudes
        if np.nanmin(lats) < minlat:
            minlat = np.nanmin(lats)
        if np.nanmax(lats) > maxlat:
            maxlat = np.nanmax(lats)
      
    # final latitude check to avoid pole
    if maxlat>85: 
        print('final check: bbox too close too pole; artifical boundary applied (85N)')
        maxlat=85  

    ## final check
    if pos_minlon == 9999: pos_minlon = np.nan
    if neg_minlon == 9999: neg_minlon = np.nan
    if pos_maxlon == -9999: pos_maxlon = np.nan
    if neg_maxlon == -9999: neg_maxlon = np.nan
    
    minlon = neg_minlon
    if np.isnan(minlon): minlon = neg_maxlon
    if np.isnan(minlon): minlon = pos_minlon
    maxlon = pos_maxlon
    if np.isnan(maxlon): maxlon = pos_minlon
    if np.isnan(maxlon): maxlon = neg_maxlon   
    
    if round(minlon) == -180 and round(maxlon, -1) == 180:
        print('1 almost there...')
        minlon = neg_maxlon
        maxlon = pos_minlon
        if round(minlon) == 0 and round(maxlon, -1) == 0:
            print('2 almost there...')
            
            all_cons = np.array([[np.nan, np.nan]])
            for cont in all_contours:
                all_cons=np.concatenate((all_cons, np.array(cont)))
            all_lons = all_cons[:,0]
            all_sort = np.sort(all_lons)
            all_sort = all_sort[~np.isnan(all_sort)]
            all_sort2 = all_sort[::-1]
            
            maxlon = all_sort2[np.nanargmax(all_sort2<-90)]
            minlon = all_sort[np.nanargmax(all_sort>-90)]
            
            rev1=True

    print('***')
    print('final bounds: ' + str([[minlat, maxlat], [minlon,maxlon]]))
    print('***')
    
    
    bbox_edges = get_edge_lines(minlon, maxlon, minlat, maxlat, reverse=rev1)
    
    return bbox_edges

#%%% new bbox

def get_bbox_edges(all_contours):
    return get_bbox_edges_new(all_contours)    

def get_bbox_edges_new(all_contours):
    minlon, minlat, maxlon, maxlat = 999,999,-999,-999

    isDivided = False
    alerted=False
    for cidx, contour in enumerate(all_contours):
        lons = contour[:,0]
        lats = contour[:,1]
        
        ### get rid of boundaries too close to pole
        for li, lat in enumerate(lats):
            if lat>82.5: ###!!! new maxlat thresh -- 85/82.5
                lons[li]=np.nan
                lats[li]=np.nan
                if not alerted:
                    print('bbox too close too pole; artifical boundary applied (85N)')
                    alerted=True

        ### convert longitude to 0-360 system
        lons1 = lons.copy()
        lons1 = np.where(lons1<0, lons1+360, lons1)
        lons1.sort()
        lons1 = lons1[~np.isnan(lons1)]
        
        if len(lons1) == 0: 
            print('skip')
            continue
        
        ### find e/w lons
        if np.nanmax(lons1) - np.nanmin(lons1) > 180:
            for li, ll in enumerate(lons1):
                if li == len(lons1)-1: break
                if lons1[li+1] - ll > 20: ###!!! threshold here
                    if not isDivided:
                        eastlon = ll
                        westlon = lons1[li+1]
                        isDivided = True
                        break
                    else:
                        if ll > eastlon:
                            eastlon = ll
                        if lons1[li+1] < westlon:
                            westlon = lons1[li+1]
                        break
        else:
            if lons1[0] < minlon:
                minlon = lons1[0]
                print('minlon ', minlon)
            if lons1[-1] > maxlon:
                maxlon = lons1[-1]
                print('maxlon ', maxlon)
        
        ### get min/max lat
        if np.nanmin(lats) < minlat:
            minlat = np.nanmin(lats)
        if np.nanmax(lats) > maxlat:
            maxlat = np.nanmax(lats)
     
        # end contour loop
        
    if not isDivided:    
        print('done easy: ', [minlon, maxlon, minlat, maxlat])
        return get_edge_lines(minlon, maxlon, minlat, maxlat, reverse=False)
    else: # isDivided
        if minlon == 999:
            print('only divide: ', [westlon, eastlon, minlat, maxlat])
            bbox_edges = get_edge_lines(westlon, eastlon, minlat, maxlat, reverse=True)
            bbox_edges = np.where(bbox_edges>180, bbox_edges-360, bbox_edges)
            return bbox_edges
        else:
            print('combo!', [minlon, maxlon], [westlon, eastlon])
            bbox_edges = get_edge_lines(westlon, maxlon, minlat, maxlat, reverse=True)
            if (eastlon<maxlon) and (westlon<minlon):
                bbox_edges = get_edge_lines(westlon,eastlon, minlat, maxlat, reverse=True)
            bbox_edges = np.where(bbox_edges>180, bbox_edges-360, bbox_edges)
            return bbox_edges

#%% grid interp

def grid_interp(old_x, old_y, old_z, new_x, new_y):
    import scipy.interpolate
    
    irregular=False
    
    if len(old_x) == 361: # ice motion
        x=new_x
        irregular=True
    
    if not irregular and len(np.shape(old_x))>1:
        old_x = old_x[0,:]
        old_y = old_y[:,0]

    if len(old_x)==360: #sst 
        x = np.where(new_x<0, new_x+360, new_x)
    elif len(old_x)==1440: #era5
        x = new_x
    elif len(old_x) == 4500: #gofs:
        x=new_x
    elif len(old_x) == 2161: #glorys: 
        x=new_x
    elif len(old_x) == 3600: #RARE
        x=new_x
        
    else:
        if not irregular: print('grid_interp: check lons')
        
    target_x = x.flatten()
    target_y = new_y.flatten()
    print('target: ', np.shape(target_x), np.shape(target_y))
    
    if not irregular:
        interp_rgi = scipy.interpolate.RegularGridInterpolator((old_x,np.flipud(old_y)), np.flipud(old_z).T, 
                                                        method='nearest',
                                                        bounds_error=False, 
                                                        fill_value = None
                                                        )
    else:
        xl = old_x.flatten()
        yl = old_y.flatten()
        grin = (xl,yl)
        interp_rgi = scipy.interpolate.LinearNDInterpolator(grin, old_z.flatten())
            
    interp_grid = interp_rgi((target_x, target_y))
    interp_grid = np.reshape(interp_grid, np.shape(new_x))
    
    if np.min(old_y) > np.min(new_y):
        min_lat = np.min(old_y)
        lat_mask = np.ma.masked_where(new_y<min_lat, new_y)
        interp_grid = np.ma.masked_array(interp_grid, mask=lat_mask.mask)
    
    return  interp_grid



#%% SEA ICE

#%%% load data

def load_seaice_old(root_dir, year, month, day, latlon=True):
    seaice_daily, lon, lat = [],[],[]
    print('load_seaice: check data version (default 3)')
    # convert date inputs to strings for filename if needed
    if not isinstance(year, str):
        year = str(year)
    if not isinstance(month, str):
        if month<10: month = '0'+str(int(month))
        else: month = str(int(month))
    if not isinstance(day, str):
        if day<10: day = '0'+str(int(day))
        else: day = str(int(day))
    
    # get file(s)
    all_files = glob.glob(root_dir + '*' + year+month+day + '*.nc')
    all_files.sort()
    variable_names = ['goddard_bt_seaice_conc','xgrid','ygrid', 'latitude','longitude']

    if not all_files:
        print('Error with filename: ' + root_dir + '*' + year+month+day + '*.nc')
        raise NameError(' bad filename in sip.load_seaice')
        
    for n, filename in enumerate(all_files):
        cdr_dic = load_netcdf(filename, variable_names)
    
        seaice = cdr_dic['goddard_bt_seaice_conc']
        
        if latlon:
            lat    = cdr_dic['latitude']
            lon    = cdr_dic['longitude']
        else:
            lat=[];lon=[]
       
        seaice_daily = seaice.mean(axis=0)
    
    return seaice_daily, lon, lat

def load_seaice_v3(root_dir, year, month, day, latlon=True):
    seaice_daily, lon, lat = [],[],[]
    
    # convert date inputs to strings for filename if needed
    if not isinstance(year, str):
        year = str(year)
    if not isinstance(month, str):
        if month<10: month = '0'+str(int(month))
        else: month = str(int(month))
    if not isinstance(day, str):
        if day<10: day = '0'+str(int(day))
        else: day = str(int(day))
    
    # get file(s)
    all_files = glob.glob(root_dir + '*' + year+month+day + '*.nc')
    all_files.sort()
    variable_names = ['goddard_bt_seaice_conc','xgrid','ygrid', 'latitude','longitude']

    if not all_files:
        print('Error with filename: ' + root_dir + '*' + year+month+day + '*.nc')
        raise NameError(' bad filename in sip.load_seaice')
        
    for n, filename in enumerate(all_files):
        cdr_dic = load_netcdf(filename, variable_names)
    
        seaice = cdr_dic['goddard_bt_seaice_conc']
        
        if latlon:
            lat    = cdr_dic['latitude']
            lon    = cdr_dic['longitude']
            
            # ygrid    = cdr_dic['ygrid']
            # xgrid    = cdr_dic['xgrid']
        else:
            lat=[];lon=[]
       
        seaice_daily = seaice.mean(axis=0)
    
    return seaice_daily, lon, lat

def load_seaice_v4(root_dir, year, month, day, latlon=True):
    seaice_daily, xgrid, ygrid = [],[],[]
    
    # convert date inputs to strings for filename if needed
    if not isinstance(year, str):
        year = str(year)
    if not isinstance(month, str):
        if month<10: month = '0'+str(int(month))
        else: month = str(int(month))
    if not isinstance(day, str):
        if day<10: day = '0'+str(int(day))
        else: day = str(int(day))
    
    # get file(s)
    all_files = glob.glob(root_dir + '*' + year+month+day + '*.nc')
    all_files.sort()
    variable_names = ['nsidc_bt_seaice_conc','xgrid','ygrid']

    if not all_files:
        print('Error with V4 filename: ' + root_dir + '*' + year+month+day + '*.nc')
        raise NameError(' bad filename in sip.load_seaice')
        
    for n, filename in enumerate(all_files):
        try:
            cdr_dic = load_netcdf(filename, variable_names)
            seaice = cdr_dic['nsidc_bt_seaice_conc']
        except KeyError:
            print('V3 sea ice data used !')
            variable_names = ['goddard_bt_seaice_conc','xgrid','ygrid']
            cdr_dic = load_netcdf(filename, variable_names)
            seaice = cdr_dic['goddard_bt_seaice_conc']
        
        if latlon:
            ygrid    = cdr_dic['ygrid']
            xgrid    = cdr_dic['xgrid']
        else:
            ygrid=[];xgrid=[]
       
        seaice_daily = seaice.mean(axis=0)
    
    return seaice_daily, xgrid, ygrid


def load_seaice(root_dir, year, month, day, latlon=True):
    
    seaice_daily, x, y = load_seaice_v4(root_dir, year, month, day, latlon=True)
    
    if not latlon: return seaice_daily
    
    ds = xr.open_dataset(root_dir + 'seaice_lonlat_v03.nc', decode_times=False)
    lon = ds['longitude'].values
    lat = ds['latitude'].values
    
    return seaice_daily, lon, lat


#%%% contours and area calcs

def plot_geocontour(ax, lon, lat, var, levels, color='k', lw=3, ls='solid'):
    #do masked-array on the lon
    lon_greater = np.ma.masked_greater(lon, -0.01)
    lon_lesser = np.ma.masked_less(lon, 0)    
    # apply masks to other associate arrays: lat
    lat_greater = np.ma.MaskedArray(lat, mask=lon_greater.mask)
    lat_lesser = np.ma.MaskedArray(lat, mask=lon_lesser.mask)
    # apply masks to other associate arrays: daily ice
    si_greater = np.ma.MaskedArray(var, mask=lon_greater.mask)
    si_lesser = np.ma.MaskedArray(var, mask=lon_lesser.mask)

    # contours
    ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
              linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
              linestyles=ls) 
    ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
              linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
              linestyles=ls)
    return ax

def plot_seaicecontour(fpath, year, month, day, ax=[], label=[], loc='lower right',
                       color='k', linewidth=4, levels=[0.15], ls='solid',legend=True, 
                       zorder=False):
    from matplotlib.lines import Line2D
    
    ''' day=[] for monthly average'''
    if type(year) == str: year = int(year)
    
    if not ax: ax = background_plot()
    lw=linewidth
    
    if day == []:  # monthly average
        all_si_month = []
        num_days = calendar.monthrange(int(year), month)[1]
        days = [datetime(int(year), month, day) for day in range(1, num_days+1)]
        for thisday in days:
            #latlon=False if day.day != 1 else True
            si_day, si_lon, si_lat = load_seaice(fpath, year, month, thisday.day, latlon=True)
            all_si_month.append(si_day)
        si_day = np.ma.mean(all_si_month, axis=0)
    else:    
        # load ice data
        si_day, si_lon, si_lat = load_seaice(fpath, year, month, day, latlon=True)

    #do masked-array on the lon
    lon_greater = np.ma.masked_greater(si_lon, -0.01)
    lon_lesser = np.ma.masked_less(si_lon, 0)    
    # apply masks to other associate arrays: lat
    lat_greater = np.ma.MaskedArray(si_lat, mask=lon_greater.mask)
    lat_lesser = np.ma.MaskedArray(si_lat, mask=lon_lesser.mask)
    # apply masks to other associate arrays: daily ice
    si_greater = np.ma.MaskedArray(si_day, mask=lon_greater.mask)
    si_lesser = np.ma.MaskedArray(si_day, mask=lon_lesser.mask)
    

    # contours
    levels = levels # 15% ice extent definition
    if not zorder:
        ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
                  linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
                  linestyles=ls) 
        ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
                  linewidths = lw, zorder=10, transform=ccrs.PlateCarree(),
                  linestyles=ls)
    else:
        ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
                  linewidths = lw, zorder=-10, transform=ccrs.PlateCarree(),
                  linestyles=ls) 
        ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
                  linewidths = lw, zorder=-10, transform=ccrs.PlateCarree(),
                  linestyles=ls)
    
    # legend
    if legend:
        if not label:
            monthstr = '0'+str(month) if month < 10 else str(month)
            if day != []:
                daystr = '-0'+str(day) if day < 10 else '-'+str(day)
            else: daystr = ' Average'
            mystr = str(year) +'-'+ monthstr + daystr
            label = 'Sea Ice Extent: ' + mystr
            
        legend_elements = [Line2D([0], [0], color=color,linewidth=lw,
                                  label=label)]
        if loc != None or loc!=[] or loc=='':
            thislegend = ax.legend(handles=legend_elements, loc=loc)
            ax.add_artist(thislegend)
    
    return ax

def getIce_TwoContours(ice_fname, level1, level2, year, month, day, plot=False):
    
    # make sure correct input type; level1 is outer counter (<lev2)
    if type(level1)==float: level1 = [level1]
    if type(level2)==float: level2 = [level2]
    if level2 < level1:
        max_lev = level1
        level1 = level2
        level2 = max_lev

    # load data and mask
    si, si_lon, si_lat = load_seaice(ice_fname, year, month, day, latlon=True)
    
    si_masked = np.ma.masked_where(si < level1, si)
    si_masked = np.ma.masked_where(si > level2, si_masked)
    
    if plot:    
        fig, ax = background_plot(extent=[-160,90,50,60], returnfig=True)
        ax = plot_seaicecontour(ice_fname, year, month, day, ax=ax,levels=level1,color='red')
        ax = plot_seaicecontour(ice_fname, year, month, day, ax=ax,levels=level1, color='red', linewidth=2.5)
        ax.pcolormesh(si_lon,si_lat, si_masked, transform=ccrs.PlateCarree())

    if plot: return fig, ax
    
    return si_masked

def calc_seaice_area(my_sic_grid, cell_area=25*25):
    # my_area = [ [ elem*cell_area for elem in L if elem] for L in my_sic_grid ]
    my_area = np.where(~np.isnan(my_sic_grid), my_sic_grid*cell_area, np.nan)
    area = np.nansum(np.nansum(my_area))
    return area

def get_total_area(masked_grid, cell_area=25*25):
    ''' does not multiply by cell element! '''
    # area_array = [ [ cell_area for elem in L if elem] for L in masked_grid ]
    area_array = np.where(~np.isnan(masked_grid), cell_area, np.nan)
    my_area = np.nansum(np.nansum(area_array))
    return my_area

#%%% storm distance from ice edge
'''storm_location3.py'''

def seaicecontour(fpath, year, month, day, ax=[], label=[], \
                       color='k', linewidth=4, levels=[0.15], mask = None):
    import calendar
    if type(year) == str: year = int(year)
    
    if not ax: ax = background_plot(returnfig=True)
    lw=linewidth
    
    if day == []:  # monthly average
        all_si_month = []
        num_days = calendar.monthrange(int(year), month)[1]
        days = [datetime(int(year), month, day) for day in range(1, num_days+1)]
        for thisday in days:
            #latlon=False if day.day != 1 else True
            si_day, si_lon, si_lat = load_seaice(fpath, year, month, thisday.day, latlon=True)
            all_si_month.append(si_day)
        si_day = np.ma.mean(all_si_month, axis=0)
    else:    
        # load ice data
        si_day, si_lon, si_lat = load_seaice(fpath, year, month, day, latlon=True)

    #do masked-array on the lon
    lon_greater = np.ma.masked_greater(si_lon, -0.01)
    lon_lesser = np.ma.masked_less(si_lon, 0)    
    # apply masks to other associate arrays: lat
    lat_greater = np.ma.MaskedArray(si_lat, mask=lon_greater.mask)
    lat_lesser = np.ma.MaskedArray(si_lat, mask=lon_lesser.mask)
    # apply masks to other associate arrays: daily ice
    si_greater = np.ma.MaskedArray(si_day, mask=lon_greater.mask)
    si_lesser = np.ma.MaskedArray(si_day, mask=lon_lesser.mask)
    
    if mask is not None:
        si_greater = np.ma.MaskedArray(si_day, mask=mask)
        si_lesser = np.ma.MaskedArray(si_day, mask=mask)
        

    # contours
    levels = levels # 15% ice extent definition
    contours1 = ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
                  linewidths = lw, zorder=10,
                  transform=ccrs.PlateCarree()) #,extent=(-20,20,50,90),
    contours2 = ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
                  linewidths = lw, zorder=10,
                  transform=ccrs.PlateCarree())
    return ax, [contours1,contours2]


def get_storm_distance(ice_fname, current_date, storm_lat, bbox_edges=[]):
    # load 15% sea ice contour
    fig, ax = background_plot(returnfig=True)
    ax, conts = seaicecontour(ice_fname, current_date.year, current_date.month, 
                            current_date.day,ax=ax,levels=[0.15],color='black', 
                            linewidth=2.5)
    contlon, contlat = [],[]
    for CS in conts:
       for lst in CS.allsegs: #get contour coord arrays
           for ix, pts in enumerate(lst): #
               if ix==0: 
                   all_pts = pts
               else:
                   all_pts = np.concatenate((all_pts, pts), axis=0)
                   
               pts=np.array(pts)
               for yi, y in enumerate(pts[:,0]):
                    if y<0:
                        pts=np.concatenate((pts, np.array([[y+360],[pts[yi,1]]]).T ))

               ### average latitudes of contour (entire globe or latitudes within bbox?])
               if bbox_edges != []:
                   # get contour within box
                   cont_in = find_points_in_contour(bbox_edges, pts[:,0], pts[:,1])
                   contlat1 = np.ma.masked_array(pts[:,1],mask=cont_in)
                   contlon1 = np.ma.masked_array(pts[:,0],mask=cont_in)
                   
                   if len(np.unique(contlat1.mask)) > 1:
                       for j, jj in enumerate(contlat1.data[~contlat1.mask]):
                           contlat.append([contlat1.data[~contlat1.mask][j]])
                           contlon.append([contlon1.data[~contlon1.mask][j]])
               elif bbox_edges == []:
                   contlat = pts[:,1]
                   
    ### average latitudes of contour (entire globe or latitudes within bbox?])
    avglat = np.mean(contlat)

    ### calculate distance
    # subtract latitudes * some distance (111 km or so)
    km = 111
    dist = (avglat-storm_lat)*km
   
    plt.close(fig)
    return dist

#%% WINDS

def plot_winds(x,y,u,v, ax = [], plot_title='', mask=[], legend=True, timeDT='', sc=300, interval=12):  

    u = np.squeeze(u)
    v = np.squeeze(v)    

    if mask!=[]:
        u=np.ma.masked_array(u, mask=mask)
        v=np.ma.masked_array(v, mask=mask)
    
    if not ax: ax=background_plot()
    
    #do masked-array on the lon
    x_greater = np.ma.masked_greater(x, -0.01)
    x_lesser = np.ma.masked_less(x, 0)    
    # apply masks to other associate arrays: lat
    y_greater = np.ma.MaskedArray(y, mask=x_greater.mask)
    y_lesser = np.ma.MaskedArray(y, mask=x_lesser.mask)
    # apply masks to other associate arrays: u,v
    u_greater = np.ma.MaskedArray(u, mask=x_greater.mask)
    u_lesser = np.ma.MaskedArray(u, mask=x_lesser.mask)
    v_greater = np.ma.MaskedArray(v, mask=x_greater.mask)
    v_lesser = np.ma.MaskedArray(v, mask=x_lesser.mask)
    
    
    ax.quiver(x_greater[::interval, ::interval], y_greater[::interval, ::interval],
            u_greater[::interval, ::interval], v_greater[::interval, ::interval],
            color='k', transform = ccrs.PlateCarree(), scale = sc, zorder = 6,
            headwidth=2.5, headlength=3)
    ax.quiver(x_lesser[::interval, ::interval], y_lesser[::interval, ::interval],
            u_lesser[::interval, ::interval], v_lesser[::interval, ::interval],
            color='k', transform = ccrs.PlateCarree(), scale = sc, zorder = 6,
            headwidth=2.5, headlength=3)
    
    if legend:
        import matplotlib.patches as mpatches
        try: black_patch = mpatches.Patch(color='black', label=timeDT.strftime("%Y-%m-%d"))
        except: black_patch = mpatches.Patch(color='black', label=timeDT)
        plt.legend(handles=[black_patch])
        
    return ax

def plot_storm_ice_maps(ice_fname, year, storm_range, bbox_edges):
    
    ### constants
    vmin = -60
    vmax = 60
    storm_len = len(storm_range)
    
    ### set up plot(s)
    fig, axes = plt.subplots(1,storm_len, figsize=(15,15*storm_len),
                             subplot_kw={'projection':ccrs.NorthPolarStereo(central_longitude=-45)})
    for ax in axes:
        ax.coastlines('50m',edgecolor='black',linewidth=0.75)
        # ax.set_extent((-extent,extent,-extent, extent),
        #               crs=ccrs.NorthPolarStereo())
        ax.set_extent((np.min(bbox_edges[:,0])-2.5,np.max(bbox_edges[:,0])+2.5,
                       np.min(bbox_edges[:,1]-5),np.max(bbox_edges[:,1])+5),
                      crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, facecolor='0.75')
        ax.add_feature(cfeature.LAKES, facecolor='0.85')
        
    ### load winds
    dummy, x, y, dummy = era5_daily(year, ['01'], ['01'], variable = '10m_v_component_of_wind')
    # find winds within bbox
    in_mask_w = find_points_in_contour(bbox_edges, x, y)

    ### loop thru storm days, plotting
    for sidx, storm in enumerate(storm_range):
        
        ### load daily (meridional, n/s) wind
        monthstr= str(storm.month) 
        daystr = str(storm.day)
        if storm.month<10: monthstr='0'+monthstr
        if storm.day<10: daystr='0'+daystr
       
        ### plot change in sea ice
        sic1, lon, lat = load_seaice(ice_fname, storm.year, storm.month, storm.day)
        storm2 = storm + timedelta(days=1)
        sic2, _,_ = load_seaice(ice_fname, storm2.year, storm2.month, storm2.day, latlon=False)
        sic_diff = (sic2 - sic1)*100
        
        axes[sidx].pcolormesh(lon,lat, sic_diff, cmap='RdBu', 
                  vmin = vmin, vmax = vmax, transform=ccrs.PlateCarree())
        axes[sidx].set_title(storm2.strftime('%b %d').replace(' 0', ' ') 
                             +' - ' + storm.strftime('%d').replace('0', ''))
        
        ### winds
        # plot bbox
        axes[sidx].plot(bbox_edges[:,0], bbox_edges[:,1], color='midnightblue',
                        transform=ccrs.PlateCarree())
        
        # plot winds
        try:
            timeDT, x, y, u = era5_daily(year, [monthstr], [daystr], variable = '10m_u_component_of_wind')
            timeDT, x, y, v = era5_daily(year, [monthstr], [daystr], variable = '10m_v_component_of_wind')
            ax=plot_winds(x,y,u,v, ax = axes[sidx], plot_title='', mask=in_mask_w, legend=False)
        except Exception as e:
            print('... unable to make wind/ice plots :/')
            print(e)
            print(traceback.format_exc())
        
    return fig, axes


#%% SST

def get_sst_grid_old(lon, lat, in_mask, sst_in, si_lon, si_lat, si_in):
    from scipy.interpolate import NearestNDInterpolator
    
    # get masked sst grid
    x = np.ma.masked_array(lon, mask=in_mask)
    y = np.ma.masked_array(lon, mask=in_mask)
    xy = np.vstack([x.ravel(), y.ravel()])
    
    ### Find sst that correspond to each sea ice cell
    vv = np.vstack([sst_in.ravel()])
    myInterpolator = NearestNDInterpolator(xy.T, vv.T)

    sst_ice_grid = []
    for ix,iy in np.ndindex(si_lon.shape):
        sst_ice_grid.append(myInterpolator(si_lon[ix,iy],si_lat[ix,iy]))    
    sst_ice_grid = np.reshape(sst_ice_grid, si_lon.shape)
    sst_ice_grid = np.ma.masked_array(sst_ice_grid, mask=si_in) 
    
    return sst_ice_grid


def get_sst_grid(lon, lat, sst_in, si_lon, si_lat, si_in):
    
    new_sst_grid = grid_interp(lon, lat, sst_in, si_lon, si_lat)
    sst_ice_grid = np.ma.masked_array(new_sst_grid, si_in)
    
    return sst_ice_grid

def get_extra_deriv_days_sst(sst_fname, mask_fname, bbox_edges, ice_fname, year, 
                         daymonth_grid, sst_miz):
    sst_series_deriv = []
    
    start_dm = daymonth_grid[0]
    start_dt = datetime(year, start_dm[0], start_dm[1])
    for dd in [1,2,3]:
        new_days = [start_dt - timedelta(days=dd)]
    
    end_dm = daymonth_grid[-1]
    end_dt = datetime(year, end_dm[0], end_dm[1])
    for dd in [1]:
        new_days.append(end_dt+timedelta(days=dd))
        
    new_grid = [(dt.month, dt.day) for dt in new_days]
    new_series = sst_timeseries(sst_fname, mask_fname, bbox_edges, ice_fname, year, new_grid)
    
    sst_series_deriv = new_series[0:3] + sst_miz + new_series[3:]
        
    return sst_series_deriv

#%% -----------------------
#%% T I M E  S E R I E S
#%% -----------------------

#%%% Plotting

def plot_timeseries(number_of_axes, list_of_tseries, legend_labels, ylabels,zero_lines, 
                    daymonth_grid, year, census_path, title = '', duration='error',
                    plot_storm_shading=True, figsize=(12,9), linewidth=2.5,
                    legend_loc=['best'], legend_cols = [1]):
    '''list of tseries: [ax1,ax2] = [[curve1,curve2], [curve1]]'''
    ''' plots based on days relative to start '''
    from matplotlib.ticker import MultipleLocator
    
    figure, axes = plt.subplots(number_of_axes ,1, sharex=True, figsize=figsize)
    if number_of_axes == 1:
        axes = [axes]
    
    # make sure all lists are the same size for plotting
    if len(list_of_tseries) != number_of_axes:
        print('ERROR sip.plot_timeseries: incorrect number of axes')
        return figure, axes
    if len(zero_lines) != number_of_axes or type(zero_lines) != list:
        print('plot_timeseries error: zero lines')
        zero_lines = [False]*number_of_axes
    if len(ylabels) != number_of_axes:
        print('plot_timeseries error: ylabels')
        while len(ylabels) < number_of_axes:
            ylabels.append(' ')
    if len(legend_loc) != number_of_axes:
        legend_loc = list(legend_loc)*number_of_axes
    if len(legend_cols) != number_of_axes:
        legend_cols = list(legend_cols)*number_of_axes
        
# if storm_shading, set up colormap, get pressure value
    if plot_storm_shading:
        import matplotlib
        new_cmap = truncate_colormap(plt.get_cmap('BuPu_r'),0,0.75)
        norm = matplotlib.colors.Normalize(vmin=960, vmax=984)
        [startDT, endDT], dummy, pressure = readCensus(census_path+'census_'+str(year)+'.csv', convertDT=True)
        for idx, start in enumerate(startDT):
            if start.day == daymonth_grid[7][1] and start.month == daymonth_grid[7][0]:
                plot_pressure = float(pressure[idx])
    
    # day1 = daymonth_grid[0]; day2=daymonth_grid[-1]
    # alldates = daterange(datetime(int(year),day1[0],day1[1]),datetime(int(year),day2[0],day2[-1]), dt=24)
    # xxx = np.arange(0,len(alldates))
    alldates = np.arange(-7,14+1,1)
    labels = [-7] +['']*6 + [0] + ['']*6 + [7] + ['']*6 +[14]
    
    # plot curves
    for i in range(number_of_axes):
        curves = list_of_tseries[i]
        legend_lab = legend_labels[i]
        # check that there are enough legend entries
        if len(legend_lab) != len(curves):
            print('plot_timeseries error: legend')
            while len(legend_lab) < len(curves):
                legend_lab.append('_')
        # plot
        for j, curve in enumerate(curves):
            axes[i].plot(alldates, curve, linewidth=linewidth, label=legend_lab[j])
        axes[i].legend(loc = legend_loc[i], ncol = legend_cols[i])
        # storm shading
        if plot_storm_shading and duration!='error':
            axes[i].axvspan(0, int(duration), alpha=0.25, color=new_cmap(norm(plot_pressure)))
        # reference zero lines
        if zero_lines[i]:
            axes[i].plot(alldates, [0]*len(alldates), 'k--', linewidth=1)
           
           
    # prep axes, add labels
    axes[0].set_title(title) 
    axes[-1].set_xlabel('Days from Storm Start')
    for j in range(number_of_axes):
      axes[j].set_ylabel(ylabels[j])
    axes[-1].xaxis.set_minor_locator(MultipleLocator(1))
    
    axes[-1].set_xticks(alldates)
    axes[-1].set_xticklabels(labels, minor=False, rotation=0, fontsize=16)
    
    txt = 'Storm Duration: ' + str(duration) + ' days'
    figure.text(0.75, 0.95, txt, fontsize=14)
            
    return figure,axes

#%%% Sea Ice Area

def relative_loss_anomaly(series):
    
    rel_series = [np.nan]
    for ix, val in enumerate(series[:-1]):
        
        if ix == 0: continue
        
        deriv = series[ix+1] - series[ix-1] # divide by time - units?
        
        rel_series.append(deriv)
        
    rel_series.append(np.nan)
    
    return np.array(rel_series)

def seaice_calculations(ice_fname, daymonth_grid, year, bbox_edges, 
                        levels, yearrange_clim=True, clim_years=[2000,2019]):

    # get si lon/lat grid
    dummy, si_lon, si_lat = load_seaice(ice_fname, 2010, 8, 1, latlon=True)
    # get box mask
    in_mask = find_points_in_contour(bbox_edges, si_lon, si_lat)
    
    # load data
    annual_si = []
    seaice_area, seaice_frac = [],[]
    for month, day in daymonth_grid:
        
        # get average over all years for this day
        if yearrange_clim:
            annual_si_yearly = []
            for yr in np.arange(clim_years[0], clim_years[1]+1, 1):
                si_yr = getIce_TwoContours(ice_fname, levels[0], levels[1], yr, month, day, plot=False)
                si_yr_in = np.ma.masked_array(si_yr, mask=in_mask)
                siarea_yr = calc_seaice_area(si_yr_in)
                annual_si_yearly.append(siarea_yr)
            annual_si.append(np.nanmean(annual_si_yearly))
    
        ### load daily sea ice w/in region of interest
        si_masked = getIce_TwoContours(ice_fname, levels[0], levels[1], year, month, day)    
        si_in = np.ma.masked_array(si_masked, mask=in_mask)
    
        ### do prev calcs
        # sea ice area within contour
        siarea = calc_seaice_area(si_in)       
        # percent area covered by sea ice        
        withice = len(np.ma.nonzero(si_in)[0])
        noice = np.count_nonzero(si_in==0)
        try:
            sifrac = withice/(withice+noice)
        except ZeroDivisionError:
            sifrac = np.nan    
        # append
        seaice_area.append(siarea)
        seaice_frac.append(sifrac)
        
    if not yearrange_clim:
        return seaice_area
    
    ### do calcs
    
    # anomaly for this year
    anomaly = seaice_area - np.mean(seaice_area)
    # annual anomaly (clim_years)
    annual_si_diff = np.array(seaice_area) - np.array(annual_si)
    yearly_anomaly = annual_si - np.mean(annual_si)
    # # relative difference
    # rel_diff = anomaly - yearly_anomaly
    
    ### new relative anomaly
    sia_rel = relative_loss_anomaly(seaice_area)
    clim_rel = relative_loss_anomaly(annual_si)
    rel_diff = sia_rel - clim_rel
    
    return [seaice_area, anomaly, annual_si, annual_si_diff, yearly_anomaly, rel_diff]

def sia_timeseries(ice_fname, census_path, year, daymonth_grid, bbox_edges, levels, title='Daily Sea Ice Area', storm_len='', clim_years=[2000,2019]):
    
    [seaice_area, anomaly, annual_si, annual_si_diff, yearly_anomaly, rel_diff] = \
    seaice_calculations(ice_fname, daymonth_grid, year, bbox_edges, levels,clim_years=clim_years)
        
    number_of_axes = 2
    list_of_tseries = [[seaice_area, annual_si],[rel_diff]] 
    legend_labels = [[str(year),"Daily average" +str(clim_years[0])+"-"+str(clim_years[1])],['Relative Anomaly']]
    ylabels = [r'Total Sea Ice Area [km$^2$]',r'[km$^2$]']
    zero_lines = [False,True]
    
    
    figure, axes = plot_timeseries(number_of_axes, list_of_tseries, legend_labels, 
                                   ylabels,zero_lines,  daymonth_grid, year, census_path, title = title, 
                                   duration=storm_len,plot_storm_shading=True, figsize=(12,9), linewidth=2.5)
    
    for ax in axes:
        ax.yaxis.set_major_formatter(OOMFormatter(3, "%1.0f"))
    
    output = [seaice_area, annual_si, rel_diff]
    return figure, output


#%%% Wind
#grid_interp(old_x, old_y, old_z, new_x, new_y)

def wind_timeseries(storm_event, year, daymonth_grid, bbox_edges,
                    ice_fname, census_path, 
                    si_levels=[0.15,0.80], title='', durationstr=''): 
    
    ### load starting wind grid
    dummy, x, y, dummy = era5_daily(year, ['01'], ['01'], 
                                   variable = '10m_v_component_of_wind')
    
    ### load sea ice
    # get si lon/lat grid
    dummy, si_lon, si_lat = load_seaice(ice_fname, 2010, 8, 1, latlon=True)
    # get box mask
    in_mask_si = find_points_in_contour(bbox_edges, si_lon, si_lat)

    daily_avg_wind = []
    storm_avg_wind = []
    
    ### daymonth loop
    count=1
    storm_day_index = 0
    for month, day in daymonth_grid:
        count+=1
        current_date = datetime(int(year), month, day)
        
        ### load daily (meridional, n/s) wind
        v_avg = era5_daily(year, month, day, variable = '10m_v_component_of_wind')[-1]


        v_grid = grid_interp(x,y, v_avg, si_lon, si_lat)
        
        si_masked= getIce_TwoContours(ice_fname, 
                                    [si_levels[0]], [si_levels[1]], 
                                    year, month, day, plot=False)
        
        v_in = np.ma.masked_array(v_grid, mask=in_mask_si)
        v_miz = np.ma.masked_array(v_in, mask=si_masked.mask)
        
        '''here''' ###!!!
        v_miz2= v_miz.copy().filled(np.nan)
        # v_miz2.flags.writable=True
        
        vmean = np.nanmean(v_miz2)
        daily_avg_wind.append(vmean) 
              
        if current_date == storm_event[storm_day_index]:
            storm_day_index += 1
            storm_avg_wind.append(vmean)
            if storm_day_index == len(storm_event): storm_day_index=0


    ### plot winds
    textstr = 'Average wind during storm: '+ str(round(np.mean(storm_avg_wind),3)) + r'm s$^{-1}$'
    fig_w, axes_w = plot_timeseries(1, [[daily_avg_wind]], [['_']], [r'Average v-wind [m s$^{-1}$]'],[True], 
                        daymonth_grid, year, census_path, title = title, duration=durationstr,
                        plot_storm_shading=True, figsize=(12,9), linewidth=2.5,
                        legend_loc=['best'], legend_cols = [1])
    
    axes_w[-1].text(0,-0.125, textstr,  fontsize=18,
        verticalalignment='top',transform = axes_w[-1].transAxes,
         bbox=dict(boxstyle="round",
                   ec='white',fc='white')
        );
    
    return daily_avg_wind, np.mean(storm_avg_wind), fig_w, axes_w


def wind_timeseries_old(storm_event, year, daymonth_grid, bbox_edges,
                    ice_fname, census_path, 
                    si_levels=[0.15,0.80], title='', durationstr=''): 
    
    from scipy.interpolate import NearestNDInterpolator

    ### load starting wind grid
    dummy, x, y, dummy = era5_daily(year, ['01'], ['01'], 
                                   variable = '10m_v_component_of_wind')
    
    ### load sea ice
    # get si lon/lat grid
    dummy, si_lon, si_lat = load_seaice(ice_fname, 2010, 8, 1, latlon=True)
    # get box mask
    in_mask_si = find_points_in_contour(bbox_edges, si_lon, si_lat)
    in_mask_winds = find_points_in_contour(bbox_edges, x,y)
    # # get masked wind grid
    x = np.ma.masked_array(x, mask=in_mask_winds)
    y = np.ma.masked_array(y, mask=in_mask_winds)
    xy = np.vstack([x.ravel(), y.ravel()])
    
    daily_avg_wind = []
    storm_avg_wind = []
    
    ### daymonth loop
    count=1
    storm_day_index = 0
    for month, day in daymonth_grid:
        count+=1
        current_date = datetime(int(year), month, day)
        
        ### load daily (meridional, n/s) wind
        v_avg = era5_daily(year, month, day, variable = '10m_v_component_of_wind')[-1]

        if False:
            ### get daily sic values between contours and mask
            si_masked= getIce_TwoContours(ice_fname, 
                                              [si_levels[0]], [si_levels[1]], 
                                              year, month, day, plot=False)
            si_masked=np.ma.masked_array(si_masked, mask=in_mask_si)
        
            si_lon = np.ma.masked_array(si_lon, mask = si_masked.mask)
            si_lat = np.ma.masked_array(si_lat, mask = si_masked.mask)
    
            ### Find wind values that correspond to each sea ice cell
            vv = np.vstack([v_avg.ravel()])
            myInterpolator = NearestNDInterpolator(xy.T, vv.T)
        
            wind_ice = []
            for ix,iy in np.ndindex(si_lon.shape):
                wind_ice.append(myInterpolator(si_lon[ix,iy],si_lat[ix,iy]))    
            wind_ice=np.reshape(wind_ice, si_lon.shape)
            wind_ice1 = np.ma.masked_array(wind_ice, mask=si_masked.mask)
            
            try:
                wind_ice1.flags.writeable = True
                
                wind_ice2 = np.array(wind_ice1).copy()
                wind_ice2.flags.writeable = True
                wind_ice2 = np.ma.masked_array(wind_ice2, mask=wind_ice1.mask)
                
                ### average north/south component
                wind_ice3 = wind_ice2.filled(np.nan)
                MEAN_ = np.nanmean(wind_ice3) 
                daily_avg_wind.append(MEAN_) 
               
            except ValueError as ve:
              print('-- wind_ice nanmean read only issue')
              print(ve)
              daily_avg_wind.append(np.nan) 
              MEAN_ = np.nan
              
        else:
            v_grid = grid_interp(x,y, v_avg, si_lon, si_lat)
            
            si_masked= getIce_TwoContours(ice_fname, 
                                        [si_levels[0]], [si_levels[1]], 
                                        year, month, day, plot=False)
            
            v_in = np.ma.masked_array(v_grid, mask=in_mask_si)
            v_miz = np.ma.masked_array(v_in, mask=si_masked.mask)
            
            vmean = np.nanmean(v_miz)
            daily_avg_wind.append(vmean) 
              
        if current_date == storm_event[storm_day_index]:
            storm_day_index += 1
            storm_avg_wind.append(vmean)
            if storm_day_index == len(storm_event): storm_day_index=0


    ### plot winds
    textstr = 'Average wind during storm: '+ str(round(np.mean(storm_avg_wind),3)) + r'm s$^{-1}$'
    fig_w, axes_w = plot_timeseries(1, [[daily_avg_wind]], [['_']], [r'Average v-wind [m s$^{-1}$]'],[True], 
                        daymonth_grid, year, census_path, title = title, duration=durationstr,
                        plot_storm_shading=True, figsize=(12,9), linewidth=2.5,
                        legend_loc=['best'], legend_cols = [1])
    
    axes_w[-1].text(0,-0.125, textstr,  fontsize=18,
        verticalalignment='top',transform = axes_w[-1].transAxes,
         bbox=dict(boxstyle="round",
                   ec='white',fc='white')
        );
    
    return daily_avg_wind, np.mean(storm_avg_wind), fig_w, axes_w


#%%% SST

def sst_timeseries_old(sst_fname, mask_fname, bbox_edges, ice_fname, year, daymonth_grid):
    ### Load sst data        
    data = load_netcdf(sst_fname, ['lon','lat','time','sst'])
    lon, lat = np.meshgrid(data['lon'],data['lat'])
    time = data['time']
    sst = data['sst']
    
    start = date(1800,1,1) 
    time_DT = [timedelta(t)+start for t in time]   
    
    ### load land-sea mask for plotting
    mask_data = load_netcdf(mask_fname, ['mask'])
    land_mask = mask_data['mask']
    land_mask = np.where(land_mask==0, True, False)
    
    ### start day loop
    
    ### find relevant time for (weekly) sst data
    start_date = datetime(int(year), daymonth_grid[0][0], daymonth_grid[0][1])
    for idx, dt in enumerate(time_DT):
        if start_date.date() < dt: break
    idx = idx-2 # use previous date    
    
    # get starting sea ice grid
    _, si_lon, si_lat = load_seaice(ice_fname, 2010, 8, 1)  
    si_in = find_points_in_contour(bbox_edges, si_lon, si_lat)  
        
    switch_dates=[]    
    sst_series_miz = []
    
    for month, day in daymonth_grid:
        current_date = datetime(int(year), month, day)
        
        if current_date.date() >= time_DT[idx+1]:
            idx += 1 # move to new ice age file
            switch_dates.append(day) 
            
        ### load data
        sst1 = np.ma.masked_array(sst[idx,:,:], mask=land_mask)
        
        ### get average value within bbox
        sst_in = sst1.copy()
        # sst_in = np.ma.masked_array(sst_in, mask=in_mask)
        # sst_in.flags.writeable = True
        # avg_temp_inbbox = np.nanmean(sst_in)
                
        ### interpolate this new grid to sic grid
        sst_sic_grid = get_sst_grid(lon, lat, sst_in, si_lon, si_lat, si_in)
            
        ### for each day, get mask for sic in miz
        si_masked = getIce_TwoContours(ice_fname, si_levels[0], si_levels[1], year, month, day)
        
        sst_sic_grid.flags.writeable = True
        si_masked.flags.writeable = True
        
        ### append average temp within miz for each day
        sst_miz = np.ma.masked_array(sst_sic_grid, mask=si_masked.mask)
        # sstanom_miz = np.ma.masked_array(sstanom_sic_grid, mask=si_masked.mask)
        
        try:
            sst_miz.flags.writeable = True
            sst_miz = sst_miz.filled(np.nan)
            sst_miz.flags.writeable = True
            sst_series_miz.append(np.nanmean(sst_miz))
            # anom_series_miz.append(np.nanmean(sstanom_miz))
        except ValueError as ve:
            print('-- sst_miz nanmean read only issue')
            print(ve)
            sst_series_miz.append(np.nan)
    
    return sst_series_miz, sst_miz

def sst_timeseries_old2(sst_fname, mask_fname, bbox_edges, ice_fname, year, daymonth_grid):
    from scipy.interpolate import RegularGridInterpolator
    ### get starting sea ice grid
    _, si_lon, si_lat = load_seaice(ice_fname, 2010, 8, 1) 
    
    ### get bbox mask
    in_mask = find_points_in_contour(bbox_edges, si_lon, si_lat)
    
    ### load initial sst data
    data = load_netcdf(sst_fname, ['lon','lat','time','sst'])
    lon = data['lon']
    lat = data['lat']
    time = data['time']
    sst = data['sst']
    start = date(1800,1,1) 
    time_DT = [timedelta(t)+start for t in time]   

    ### find relevant time for (weekly) sst data
    start_date = date(int(year), daymonth_grid[0][0], daymonth_grid[0][1])
    for idx, dt in enumerate(time_DT):
        if start_date < dt: break
    idx = idx-2 # use previous date 
        
        
    ### get starting grid
    idx2 = 0
    while time_DT[idx2] < start_date: 
        idx2 += 1
    idx2 = idx2 -1
    
    # print('starting grid: ', daymonth_grid[0], time_DT[idx2])
    
    sstg = sst[idx2,:,:]    
    interp_rgi = RegularGridInterpolator((lon,np.flipud(lat)), np.flipud(sstg).T, 
                                        method='linear', bounds_error=False, fill_value = None
                                        )
    
    new_x = si_lon
    x = np.where(new_x<0, new_x+360, new_x)
    
    sst_interp = interp_rgi((x, si_lat))
    sst_sic_grid = np.reshape(sst_interp, np.shape(si_lon))
    
    # print(idx, idx2)
        
    sst_series_miz = []
    for month, day in daymonth_grid:
        current_date = datetime(int(year), month, day)
        # print(current_date, time_DT[idx])
    
        if (current_date.date() >= time_DT[idx+1]):
            idx += 1 # move to new ice age file
            
            # print('*SWITCH*', current_date, time_DT[idx], time_DT[idx2])
            
            sstg = sst[idx,:,:]
            interp_rgi = RegularGridInterpolator((lon,np.flipud(lat)), np.flipud(sstg).T, 
                                                        method='linear',
                                                        bounds_error=False, 
                                                        fill_value = None
                                                        )
            sst_interp = interp_rgi((si_lon, si_lat))
            sst_sic_grid = np.reshape(sst_interp, np.shape(si_lon))
        
                
        ### for each day, get mask for sic in miz
        si_masked = getIce_TwoContours(ice_fname, si_levels[0], si_levels[1], year, month, day)
        
        ### append average temp within miz for each day
        sst_in = np.ma.masked_array(sst_sic_grid, mask=in_mask)
        sst_miz = np.ma.masked_array(sst_in, mask=si_masked.mask)
        
        '''here'''###!!!
        sst_miz2= sst_miz.copy()
        sst_miz2 = np.array(sst_miz2.filled(np.nan))
        # sst_miz2.flags.writable=True
        
        sst_series_miz.append(np.nanmean(sst_miz2))
        
    return sst_series_miz


def sst_timeseries(sst_fname, mask_fname, bbox_edges, ice_fname, year, daymonth_grid):
    from scipy.interpolate import RegularGridInterpolator
    
    ### get starting sea ice grid
    _, si_lon, si_lat = load_seaice(ice_fname, 2010, 8, 1) 
    
    ### get bbox mask
    in_mask = find_points_in_contour(bbox_edges, si_lon, si_lat)
    
    ### load initial sst data
    data = load_netcdf(sst_fname, ['lon','lat','time','sst'])
    lon = data['lon']
    lat = data['lat']
    time = data['time']
    sst = data['sst']
    start = date(1800,1,1) 
    time_DT = [timedelta(t)+start for t in time]   


    ### find relevant time for (weekly) sst data
    start_date = date(int(year), daymonth_grid[0][0], daymonth_grid[0][1])
        
    idx2 = 0
    while time_DT[idx2] <= start_date: ###!!!
        idx2 += 1
    idx2 = idx2 - 1
    
    ### get starting grid
    sstg = sst[idx2,:,:]    
    interp_rgi = RegularGridInterpolator((lon,np.flipud(lat)), np.flipud(sstg).T, 
                                        method='linear', bounds_error=False, fill_value = None
                                        )
    new_x = si_lon
    x = np.where(new_x<0, new_x+360, new_x)
    
    sst_interp = interp_rgi((x, si_lat))
    sst_sic_grid = np.reshape(sst_interp, np.shape(si_lon))
    
    print('START SST: ', [daymonth_grid[0],time_DT[idx2]])
    
    sst_series_miz = []
    for month, day in daymonth_grid:
        current_date = datetime(int(year), month, day)
    
        if (current_date.date() >= time_DT[idx2]):
            idx2 += 1 # move to new ice age file
            sstg = sst[idx2,:,:]    
            interp_rgi = RegularGridInterpolator((lon,np.flipud(lat)), np.flipud(sstg).T, 
                                                method='linear', bounds_error=False, fill_value = None
                                                )
            new_x = si_lon
            x = np.where(new_x<0, new_x+360, new_x)
            
            sst_interp = interp_rgi((x, si_lat))
            sst_sic_grid = np.reshape(sst_interp, np.shape(si_lon))
            
                
        ### for each day, get mask for sic in miz
        si_masked = getIce_TwoContours(ice_fname, si_levels[0], si_levels[1], year, month, day)
        
        ### append average temp within miz for each day
        sst_in = np.where(in_mask, np.nan, sst_sic_grid)
        
        sst_miz = np.ma.masked_array(sst_in, mask=si_masked.mask)
        
        sst_miz2= sst_miz.copy()
        sst_miz2 = np.array(sst_miz2.filled(np.nan))
        sst_series_miz.append(np.nanmean(sst_miz2))
        
    return sst_series_miz






















