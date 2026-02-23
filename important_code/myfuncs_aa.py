#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 2023
myfuncs_aa.py
> functions used for cyclone tracking 

@author: mundi
"""
#%% imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

sh_ice = '/Users/mundi/Desktop/seaice/south/'
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import cmocean.cm as cmo

import sys, os

#%% starting info, plotting
def rstr(val, n=3): return str(round(val, n))

def daterange(start_date, end_date, dt=6):
    alldates=[]
    delta = timedelta(hours=dt)
    while start_date <= end_date:
        alldates.append(start_date)
        start_date += delta
    return alldates


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
    

def background_plot(extent=[-180,180, -53,-90], returnfig=False, title=[], labels=True):
    
    fig=plt.figure(figsize=[15,15]) 
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
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
    
    if not returnfig: return ax
    if returnfig: return fig, ax
    
#%% processing
class HidePrint:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

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

def get_total_area(masked_grid, cell_area=25*25):
    ''' does not multiply by cell element! '''
    # area_array = [ [ cell_area for elem in L if elem] for L in masked_grid ]
    area_array = np.where(~np.isnan(masked_grid), cell_area, np.nan)
    my_area = np.nansum(np.nansum(area_array))
    return my_area

def calc_seaice_area(my_sic_grid, cell_area=25*25):
    # my_area = [ [ elem*cell_area for elem in L if elem] for L in my_sic_grid ]
    my_area = np.where(~np.isnan(my_sic_grid), my_sic_grid*cell_area, np.nan)
    area = np.nansum(np.nansum(my_area))
    return area

#%% loading sea ice info 

def load_seaice(root_dir, year, month, day, latlon=True):
    
    if latlon:
        seaice, lon, lat = load_seaice_sh(root_dir, year, month, day, latlon)
        return seaice, lon, lat
    else:
        seaice = load_seaice_sh(root_dir, year, month, day, latlon)
        return seaice


def load_seaice_sh(root_dir, year, month, day, latlon=True):
    import glob
    from pyproj import Transformer
    
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

    if not all_files:
        print('Error with S.Hemisphere filename: ' + root_dir + '*' + year+month+day + '*.nc')
        raise NameError(' bad filename in sip.load_seaice')
        
    cdr_dic = xr.open_dataset(all_files[0])
    seaice = np.squeeze( cdr_dic['nsidc_bt_seaice_conc'].values )
    
    for flag in [251,252,253,254,255]:
        seaice= np.where(seaice==flag/100, np.nan, seaice)
   
    if latlon:
        x,y = np.meshgrid( cdr_dic['xgrid'].values, cdr_dic['ygrid'].values )
        transformer = Transformer.from_crs("EPSG:3412", "EPSG:4326", always_xy=True)
        lon, lat = transformer.transform(x, y)
        return seaice, lon, lat
    else:
        return seaice


def seaicecontour(fpath, year, month, day, ax=[], label=[], \
                       color='k', linewidth=4, levels=[0.15]):
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
            si_day, si_lon, si_lat = load_seaice_sh(fpath, year, month, thisday.day, latlon=True)
            all_si_month.append(si_day)
        si_day = np.ma.mean(all_si_month, axis=0)
    else:    
        # load ice data
        si_day, si_lon, si_lat = load_seaice_sh(fpath, year, month, day, latlon=True)

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
    contours1 = ax.contour(lon_greater, lat_greater, si_greater, colors=color, levels=levels, 
                  linewidths = lw, zorder=10,
                  transform=ccrs.PlateCarree()) #,extent=(-20,20,50,90),
    contours2 = ax.contour(lon_lesser, lat_lesser, si_lesser, colors=color, levels=levels, 
                  linewidths = lw, zorder=10,
                  transform=ccrs.PlateCarree())
    return ax, [contours1,contours2]




# year, month, day = 2010, 8, 15
# si, lon, lat = load_seaice_sh(sh_ice, year, month, day, latlon=True)
# fig, ax = background_plot(returnfig=True)
# ax.pcolormesh(lon,lat, si, transform=ccrs.PlateCarree(), cmap=cmo.ice,vmin=0,vmax=1)
# ax, conts = seaicecontour(sh_ice, year, month, day,
#                          ax=ax,levels=[0.15],color='r', linewidth=2.5)

#%% cyclone information

def find_points_in_contour(coords, var_x, var_y, var=None):
    import matplotlib.path as mpltPath
    
    ### turn contour into polygon
    polygon = np.vstack((coords[:,0], coords[:,1])).T
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


#%% contours
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
    return x, y, pressure, mytime


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
        if 0 in list(np.shape(pressure)):
            raise NameError('get_contour_points: empty pressure ' + str((year, month, day)))
        
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


def detect_contours2(storm_daymonth_grid, daymonth_grid, year, storm_info, intervals=[990,1000], local_path='',local_file=''):
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
                                                       storm_info[snum], 
                                                       local_file=local_path+local_file+str(year)+'-'+str(month)+'.nc')
        
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
                    [x1,x2], [y1,y2] = geoplot_2d(cc[:,0], cc[:,1])
                    ax.plot(x1, y1, 'b', linewidth=5, transform=ccrs.PlateCarree())
                    ax.plot(x2, y2, 'b', linewidth=5, transform=ccrs.PlateCarree())
                continue
    
    print()
    return all_cont_dt1, all_conts1   

#%% bbox

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


def get_bbox_edges(all_contours):
    minlon, minlat, maxlon, maxlat = 999,999,-999,-999

    isDivided = False
    alerted=False
    for cidx, contour in enumerate(all_contours):
        lons = contour[:,0]
        lats = contour[:,1]
        
        ### get rid of boundaries too close to pole
        for li, lat in enumerate(lats):
            if lat < -85: ###!!!
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


#%% era5

def era5_hourly(year, month, day, variable = 'msl'):
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
        "area":[-90, -180, -50, 180]
        }# retrieves the path to the file
    fl = cds.retrieve(dataset, params)# download the file 
    if download_flag:
        fl.download("./output.nc")# load into memory
    with urlopen(fl.location) as f:
        ds = xr.open_dataset(f.read(), decode_times=False)
        
    lon, lat = np.meshgrid(ds['longitude'], ds['latitude'])
    
    if variable == 'msl':
        key = 'msl'
    elif variable =='10m_v_component_of_wind':
        key = 'v10'
    elif variable =='10m_u_component_of_wind':
        key = 'u10'
    else:
        raise ValueError('ERA5_hourly: invalid variable key')

    time_in = ds['time'].values
    hours_in = [tt-time_in[0] for tt in time_in]
    starting_day = datetime(int(year[0]), int(month[0]), int(day[0]))
    time = np.array([starting_day + timedelta(hours=int(hr)) for hr in hours_in])
    
    # time = pd.to_datetime(time_in)
    
    var = ds[key].values
    
    if variable == 'msl':
        var = var/100

    return time, lon, lat, var

def era5_daily(year, month, day, variable = 'msl'):
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
        "area":[-90, -180, -50, 180]
        }# retrieves the path to the file
    fl = cds.retrieve(dataset, params)# download the file 
    if download_flag:
        fl.download("./output.nc")# load into memory
    with urlopen(fl.location) as f:
        ds = xr.open_dataset(f.read(),decode_times=False)
        
    lon, lat = np.meshgrid(ds['longitude'], ds['latitude'])
    
    if variable == 'msl':
        key = 'msl'
    elif variable =='10m_v_component_of_wind':
        key = 'v10'
    elif variable =='10m_u_component_of_wind':
        key = 'u10'
    elif variable == '2m_temperature':
        key='t2m'
    else:
        key=variable

    time = ds['time'].values
    ts = [tt-time[0] for tt in time]
    try:
        var = ds[key].values
    except:
        print(key)
        raise ValueError('ERA5_daily: invalid variable key')
    
    if variable == 'msl':
        var = var/100

    start_day = int(day[0])
    my_pres = []
    daily_avg = []
    daily_time = []
    
    for tx, tt in enumerate(ts):
        if tt != 23:
            my_pres.append(var[tx,:,:])
        else: # save this day and move to next
            my_pres.append(var[tx,:,:]) # append last hour
            daily_avg.append( np.mean(my_pres, axis=0) )
            daily_time.append(datetime(int(year[0]), int(month[0]), start_day, ts[0]))
            
            my_pres=[]
            start_day += 1
            
    if len(day) == 1:
        daily_avg=np.squeeze(daily_avg)        

    return daily_time, lon, lat, daily_avg