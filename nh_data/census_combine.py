#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 12:04:48 2025

@author: mundi
"""

#%% imports and filepaths
import numpy as np
import csv
from datetime import datetime
from glob import glob
import os

cpath0 = '/Users/mundi/Desktop/month-hemi/nh_data/census0/'
cfiles = glob(cpath0+'*.csv')

savepath = '/Users/mundi/Desktop/month-hemi/nh_data/census/'
if not os.path.exists(savepath):
    os.makedirs(savepath)

header = ['start', 'start_lat', 'start_lon', 'minimum', 'finish', 'finish_lat', 'finish_lon']

summer_months = [6,7,8,9]

#%% summer

def get_summer_data(data, backup):
    census_path1 = '/Users/mundi/Desktop/teststuff/original_census/'
    census_file1 = census_path1+'census_'+str(year)+'.csv'
    
    try:
        csv_file1 = open(census_file1,'r')
        csv_file1.readline()
    except FileNotFoundError:
        data.append(backup)
        print("* no summer")
        return data
    
    for a, b, c, d, e, f, g in csv.reader(csv_file1, delimiter=','):
        dt1 = datetime(int(a[:4]), int(a[5:7]), int(a[8:10]), 0)
        
        if dt1.month in summer_months:
            data.append([a, b, c, d, e, f, g])
            
    return data

#%% data loop

for file in cfiles:
    year = file.split('.csv')[0][-4:]
    print(year)

    load_summer = True
    missing_summer = False
    
    data = []
    
    # Read off and discard first line, to skip headers
    csv_file = open(file,'r')
    csv_file.readline()
    
    # Split columns while reading
    for a, b, c, d, e, f, g in csv.reader(csv_file, delimiter=','):
        dt = datetime(int(a[:4]), int(a[5:7]), int(a[8:10]), 0)
        
        if dt.month not in summer_months:
            data.append([a, b, c, d, e, f, g])
        else:
            if load_summer:
                
                data = get_summer_data(data, backup = [a, b, c, d, e, f, g])
                        
                load_summer=False

    if load_summer: # skipped summer?
        data = get_summer_data(data, backup = [])
    
    ### save new csv
    with open(savepath+'census_'+str(year)+'.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        
        # write the header
        writer.writerow(header)
    
        # write multiple rows
        writer.writerows(data)
     


#%% end