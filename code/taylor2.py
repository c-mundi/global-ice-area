#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 2025
taylor2.py

taylor diagram code
* applies to actual storm data

https://gist.github.com/ycopin/3342888
https://zenodo.org/records/5548061

@author: mundi
"""
#%% start
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import calendar
from scipy.stats import linregress

import functions as fx

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

month_colors = ['#238443','#78c679','#c2e699',
                '#d7b5d8','#df65b0','#dd1c77','#980043','#7a0177',
                '#253494','#2c7fb8','#41b6c4','#a1dab4'
                ]

months = np.arange(1,12+1)
month_names = [calendar.month_name[mm] for mm in months]
month_abbrs = [calendar.month_abbr[mm] for mm in months]

decades = [np.arange(1982, 1992), np.arange(2010,2020)]
decade_names = ['Early Satellite Era ', 'Present Day ']
dnames = ['80s','10s']
dmarkers = ['s','o']

hemi_colors = ['#01665e','#8c510a']
hemi_names= ['Arctic', 'Antarctic']

decade_colors = ['#2166ac', '#b2182b']
shade_colors = ['#67a9cf', '#ef8a62']

quad_colors = ['#f7f7f7','#cccccc','#969696','#525252']

xxx = np.arange(-7,14+1,1)

#%%% class
class TaylorDiagram(object):
    """
    Taylor diagram.
    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd,
                 fig=None, rect=111, label='_', srange=(0, 1.5), extend=False):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.
        Parameters:
        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        * extend: extend diagram to negative correlations
        """

        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF

        self.refstd = refstd            # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1])
        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.pi
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi/2
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd

        ghelper = FA.GridHelperCurveLinear(
            tr,
            extremes=(0, self.tmax, self.smin, self.smax),
            grid_locator1=gl1, tick_formatter1=tf1)

        if fig is None:
            fig = plt.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Scaled Environmental Characterstic")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("Normalized Absolute Storm Impact")

        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction(
            "bottom" if extend else "left")

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)          # Unused

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        ### Add reference point and stddev contour
        # l, = self.ax.plot([0], self.refstd, 'k*',
        #                   ls='', ms=10, label=label)
        # t = np.linspace(0, self.tmax)
        # r = np.zeros_like(t) + self.refstd
        # self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [] #[l]

    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """

        l, = self.ax.plot(np.arccos(corrcoef), stddev,
                          *args, **kwargs)  # (theta, radius)
        self.samplePoints.append(l)

        return l

    def add_grid(self, *args, **kwargs):
        """Add a grid."""

        self._ax.grid(*args, **kwargs)

    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """

        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax),
                             np.linspace(0, self.tmax))
        # Compute centered RMS difference
        rms = np.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*np.cos(ts))

        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours
    
#%%% functions
def translate(value, leftMin, leftMax):
    rightMin, rightMax = -1,1
    
    # Figure out how 'wide' each range is
    leftSpan = leftMax - leftMin
    rightSpan = rightMax - rightMin

    # Convert the left range into a 0-1 range (float)
    valueScaled = float(value - leftMin) / float(leftSpan)

    # Convert the 0-1 range into a value in the right range.
    out = rightMin + (valueScaled * rightSpan)
    
    # min.max values
    if out>rightMax: out=rightMax
    elif out<rightMin: out=rightMin
    
    return out
    
#%% real data: winds

widx = 1
wnames0 = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 
          'xlims':[[-4.5,3],[-1,3]], 'vline':0}
wnames1 = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 
          'xlims':[[-4.5,2],[-1,2]], 'vline':0}
wnames2 = {'idx':0, 'ylims':[-2,14], 'name':'Total', 
          'xlims':[[4,10],[6,10]], 'vline':None}
wnames = [wnames0, wnames1, wnames2][widx]

# windNorm = plt.Normalize(vmin=-4, vmax=4)


for MONTHS in [[[9,10,11,12],[4,5,6,7]],
               ]: # [4,5,6,7,8]

    fig = plt.figure()
    fig.suptitle('Ice-Increasing Months', size='large')  # Figure title
    
    dia = TaylorDiagram(0.55, fig=fig, label='Reference', extend=True)
    dia.add_grid()                                  # Add grid
    dia._ax.axis[:].major_ticks.set_tick_out(True)  # Put ticks outward
    
    for loc_ind, loc in enumerate(hemi_names):
        months = MONTHS[loc_ind]
        for yi, years in enumerate(decades):
            yr_title = str(years[0])+'-'+str(years[-1])
            # print(yr_title)
            path1 = root_paths[loc_ind]
    
            # sea ice
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
                
            # winds
            if loc_ind==0:
                wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
            elif loc_ind==1: 
                wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
            
            for mm in months:
                
                # wind data
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[mm] if len(arr)==22]  
                    wind_lines = np.array(windies)[:,:,wnames['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames['idx']] for arr in wind_series[mm]]  
                    wind_lines = np.array(windies)[:,]
                mean_wind_line = np.nanmean(wind_lines, axis=0)
                
                r = np.corrcoef(mean_wind_line, mean_lines[mm-1])[0,1]
                
                # var1 = np.nanmean(mean_wind_line[7:10])
                var1 = translate(np.nanmean(mean_wind_line[7:14]), -4,4)
                dia.add_sample(np.abs(mean_lines[mm-1][-1]), var1,
                               marker=dmarkers[yi], ms=8, ls='',
                               mfc=month_colors[mm-1], mec='k',
                               label=dnames[yi]+' '+str(mm)+' Winds')
                
                var2 = translate(np.nanmean(mean_wind_line[7:10]), -4,4)
                dia.add_sample(np.abs(mean_lines[mm-1][7]), var2,
                               marker=dmarkers[yi], ms=8, ls='',
                               mfc=month_colors[mm-1], mec=month_colors[mm-1],
                               label=dnames[yi]+' '+str(mm)+' Winds*')
       
    # Legend
    # fig.legend(dia.samplePoints, 
    #            [ p.get_label() for p in dia.samplePoints ],
    #            numpoints=1, prop=dict(size='small'), loc='upper right', 
    #            ncol=1, bbox_to_anchor = (1.15,1))
    
    ax2 = dia.ax.twinx()
    for mgroup in MONTHS:
        for mm in mgroup:
            ax2.plot(np.NaN, np.NaN, ls='-',lw=3,
                     label=calendar.month_abbr[mm], c=month_colors[mm-1])
    for yi, dname in enumerate(dnames):
        ax2.plot(np.nan, np.nan, lw=0, marker=dmarkers[yi], color='gray', label=dname)
    for outline, lname in zip(['k','gray'],['Long-term','Short-Term']):
        ax2.plot(np.nan, np.nan, lw=0, marker='.', mfc='gray', mec=outline, label=lname, markersize=10)
    
    ax2.get_yaxis().set_visible(False)
    ax2.spines[['right', 'top','bottom','left']].set_visible(False)
    ax2.legend(loc='upper right', bbox_to_anchor = (1.33,1.05),
               handletextpad=0.5, handlelength=1)
    
#%% seasonal subplots
def new_legend(dia, entries, bbox_to_anchor, ncol=2):
    ax2 = dia.ax.twinx()
    
    for entry in entries:
        ax2.plot(np.nan, np.nan, 
                 ls=entry['ls'],
                 lw=entry['lw'], 
                 marker=entry['marker'], 
                 color=entry['color'],
                 mfc=entry['mfc'], 
                 mec=entry['mec'],
                 label=entry['label'])
        
    ax2.get_yaxis().set_visible(False)
    ax2.spines[['right', 'top','bottom','left']].set_visible(False)
    ax2.legend(loc='upper right', bbox_to_anchor = bbox_to_anchor,
               handletextpad=0.5, handlelength=1, ncol=ncol)
            

dmarkers2 = ['P','D']
dmarkers3 = ['v', '^']

widx = 1
wnames0 = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 
          'xlims':[[-4.5,3],[-1,3]], 'vline':0}
wnames1 = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 
          'xlims':[[-4.5,2],[-1,2]], 'vline':0}
wnames2 = {'idx':0, 'ylims':[-2,14], 'name':'Total', 
          'xlims':[[4,10],[6,10]], 'vline':None}
wnames = [wnames0, wnames1, wnames2][widx]


fig = plt.figure(figsize=(8,14))
fig.suptitle(' '*115+wnames['name']+' Winds', fontweight='bold')
subfigs = fig.subfigures(3, 1)

MGROUPS = [[[9,10,11,12],[4,5,6,7]],
            [[4,5,6,7,8],[10,11,12]],
            [[1,2,3],[8,9]]
            ]
MNAMES = ['Ice-Increasing', 'Ice-Decreasing', 'Minimal Change']

for mg, MONTHS in enumerate(MGROUPS): 

    subfigs[mg].suptitle(MNAMES[mg]+' Months', size='large')  # Figure title
    
    dia = TaylorDiagram(0.55, fig=subfigs[mg], label='Reference', extend=True)
    dia.add_grid()                                  # Add grid
    dia._ax.axis[:].major_ticks.set_tick_out(True)  # Put ticks outward
    
    for loc_ind, loc in enumerate(hemi_names):
        months = MONTHS[loc_ind]
        for yi, years in enumerate(decades):
            yr_title = str(years[0])+'-'+str(years[-1])
            path1 = root_paths[loc_ind]
    
            # sea ice
            mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
                
            # winds
            if loc_ind==0:
                wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
            elif loc_ind==1: 
                wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
            
            # ocn prof
            if yi==1: ocn_data, DEPTH = fx.ocn_profiles(years, path1, all_or_miz='miz')
            
            # waves
            swh_series = fx.era_lines(years, path1+'census/', path1+'area/', path1+'swh/', 'swh')[-1]
                
            
            #### month loop
            for mm in months:
                
                # wind data
                if loc_ind==0: 
                    windies = [arr for arr in wind_series[mm] if len(arr)==22]  
                    wind_lines = np.array(windies)[:,:,wnames['idx']]
                elif loc_ind==1: 
                    windies = [arr[wnames['idx']] for arr in wind_series[mm]]  
                    wind_lines = np.array(windies)[:,]
                mean_wind_line = np.nanmean(wind_lines, axis=0)
                
                # plot winds
                var1 = translate(np.nanmean(mean_wind_line[7:22]), -4,4)
                dia.add_sample(np.abs(mean_lines[mm-1][-1]), var1,
                               marker=dmarkers[yi], ms=8, ls='',
                               mfc=month_colors[mm-1], mec='k',
                               label=dnames[yi]+' '+str(mm)+' Winds')
                
                var2 = translate(np.nanmean(mean_wind_line[7:14]), -4,4)
                dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                               marker=dmarkers[yi], ms=8, ls='',
                               mfc=month_colors[mm-1], mec=month_colors[mm-1],
                               label=dnames[yi]+' '+str(mm)+' Winds*')
                
                # sst data
                SSTs = np.load(path1+'sst/'+'tseries_'+str(yi)+'-'+str(mm)+'.npy')
                mean_line = np.nanmean(SSTs, axis=0)
                slope1, b, r, p, se = linregress(xxx[7:22], mean_line[7:22])
                slope2, b, r, p, se = linregress(xxx[7:14], mean_line[7:14])
                
                # plot sst
                var1 = translate(slope1, -0.04,0.04)
                dia.add_sample(np.abs(mean_lines[mm-1][-1]), var1,
                               marker=dmarkers2[yi], ms=8, ls='',
                               mfc=month_colors[mm-1], mec='k',
                               label=dnames[yi]+' '+str(mm)+' SST')
                var2 = translate(slope2, -0.04,0.04)
                dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                               marker=dmarkers2[yi], ms=8, ls='',
                               mfc=month_colors[mm-1], mec=month_colors[mm-1],
                               label=dnames[yi]+' '+str(mm)+' SST*')

                if yi==1:
                    # upper ocean data   
                    ocn_profs = ocn_data[mm]
                    prof_diffs = [profs[-1]-profs[0] for profs in ocn_profs]
                    upper_diffs = [np.nanmean(prof[0:8]) for prof in prof_diffs] #0:8->upper10m
                    
                    # plot upper ocean
                    var2 = translate(np.nanmean(upper_diffs), -0.25,0.25)
                    dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                                   marker='*', ms=10, ls='',
                                   mfc=month_colors[mm-1], mec=month_colors[mm-1],
                                   label=dnames[yi]+' '+str(mm)+' SST*')
                
                # wave data
                swh_lines = [np.array(s) for s in swh_series[mm] if len(s)==22]
                mean_wave_line = np.nanmean(swh_lines, axis=0)
                
                # plot waves
                var1 = translate(np.nanmean(mean_wave_line[7:22]), 0,4)
                dia.add_sample(np.abs(mean_lines[mm-1][-1]), var1,
                               marker=dmarkers3[yi], ms=8, ls='',
                               mfc=month_colors[mm-1], mec='k',
                               label=dnames[yi]+' '+str(mm)+' Waves')
                
                var2 = translate(np.nanmean(mean_wave_line[7:14]), 0,4)
                dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                               marker=dmarkers3[yi], ms=8, ls='',
                               mfc=month_colors[mm-1], mec=month_colors[mm-1],
                               label=dnames[yi]+' '+str(mm)+' Waves*')
                    
    # month legend 
    month_entries = []
    for mgroup in MONTHS:
        for mm in mgroup:
            c=month_colors[mm-1]
            month_entries.append({'ls':'-', 'lw':3, 'color':c,'mfc':c,'mec':c,
                                  'label':calendar.month_abbr[mm], 'marker':None})
    new_legend(dia, month_entries, bbox_to_anchor=(-0.075,1.0))
    
    # variable legend
    var_entries = []
    for yi, dname in enumerate(dnames):
        c='gray'
        var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dname+' Winds', 'marker':dmarkers[yi]})
        var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dname+' SST', 'marker':dmarkers2[yi]})
        var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dname+' Waves', 'marker':dmarkers3[yi]})
    var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dname+' Upper Ocean', 'marker':'*'})
    new_legend(dia, var_entries, bbox_to_anchor=(1.33,1.0))
    
    # time legend
    time_entries = []
    for outline, lname in zip(['k','gray'],['Long-term','Short-Term']):
        time_entries.append({'ls':'-', 'lw':0, 'color':'gray','mfc':'gray','mec':outline,
                            'label':lname, 'marker':'p'})
    new_legend(dia, time_entries, bbox_to_anchor=(-0.075,0.5), ncol=1)
    
    
#%% notes
## x download full timeseries for upper ocean profiles
## adjust code (function?) to turn on/off variables
## check number calculations

#%% more flexible code

#%%% fxns and other
def new_legend(dia, entries, bbox_to_anchor, ncol=2):
    ax2 = dia.ax.twinx()
    
    for entry in entries:
        ax2.plot(np.nan, np.nan, 
                 ls=entry['ls'],
                 lw=entry['lw'], 
                 marker=entry['marker'], 
                 color=entry['color'],
                 mfc=entry['mfc'], 
                 mec=entry['mec'],
                 label=entry['label'])
        
    ax2.get_yaxis().set_visible(False)
    ax2.spines[['right', 'top','bottom','left']].set_visible(False)
    
    
    handles, labels = ax2.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
        
    ax2.legend(*zip(*unique), loc='upper right', 
               bbox_to_anchor = bbox_to_anchor,
               handletextpad=0.5, handlelength=1, ncol=ncol)
            

dmarkers2 = ['P','D']
dmarkers3 = ['v', '^']

widx = 1
wnames0 = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 
          'xlims':[[-4.5,3],[-1,3]], 'vline':0}
wnames1 = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 
          'xlims':[[-4.5,2],[-1,2]], 'vline':0}
wnames2 = {'idx':0, 'ylims':[-2,14], 'name':'Total', 
          'xlims':[[4,10],[6,10]], 'vline':None}
wnames = [wnames0, wnames1, wnames2][widx]

#%%% plotting (new)

# variable_groups = [['wind','t2m'], ['sst','ocn_prof', 'waves']]
# varnames = ['Atmospheric Vairables', 'Ocean Forcings']

variable_groups = [['wind', 'waves'], ['sst','ocn_prof','t2m']]
varnames = ['Mechanics', 'Thermodynamics']

for v, var_group in enumerate(variable_groups):

    fig = plt.figure(figsize=(8,14))
    fig.text(0.5,1.025, varnames[v], fontweight='bold', size='x-large',
             horizontalalignment='center', verticalalignment='center')
    
    subfigs = fig.subfigures(3, 1)
    
    MGROUPS = [[[9,10,11,12],[4,5,6,7]],
                [[4,5,6,7,8],[10,11,12]],
                [[1,2,3],[8,9]]
                ]
    MNAMES = ['Ice-Increasing', 'Ice-Decreasing', 'Ice Maximum']
    
    for mg, MONTHS in enumerate(MGROUPS): 
    
        subfigs[mg].suptitle(MNAMES[mg]+' Months', size='large')  # Figure title
        
        dia = TaylorDiagram(0.55, fig=subfigs[mg], label='Reference', extend=True)
        dia.add_grid()                                  # Add grid
        dia._ax.axis[:].major_ticks.set_tick_out(True)  # Put ticks outward
        
        
        var_entries = []; c='gray'
        for loc_ind, loc in enumerate(hemi_names):
            months = MONTHS[loc_ind]
            for yi, years in enumerate(decades):
                yr_title = str(years[0])+'-'+str(years[-1])
                path1 = root_paths[loc_ind]
        
                # sea ice
                mean_lines, lines, start_day, end_day, si_changes, clim_changes = \
                    fx.indiv_lines(years, path1+'census/', path1+'area/', path1+'seaice/')
                    
                # winds
                if 'wind' in var_group:
                    if loc_ind==0:
                        wind_series = fx.wind_lines(years, path1+'census/', path1+'area/', path1+'wind/')[-1]
                    elif loc_ind==1: 
                        wind_series = fx.sh_winds(years, path1+'census/', path1+'area/', path1+'wind/')
                    
                # ocn prof
                if 'ocn_prof' in var_group:
                    if yi==1: ocn_data, DEPTH = fx.ocn_profiles(years, path1, all_or_miz='miz')
                
                # waves
                if 'waves' in var_group:
                    swh_series = fx.era_lines(years, path1+'census/', path1+'area/', path1+'swh/', 'swh')[-1]
                    
                # air temperature
                if 't2m' in var_group:
                    t2m_series = fx.era_lines(years, path1+'census/', path1+'area/', path1+'t2m/', 't2m')[-1]
                      
                
                #### month loop
                for mm in months:
                    
                    # wind data
                    if 'wind' in var_group:
                        if loc_ind==0: 
                            windies = [arr for arr in wind_series[mm] if len(arr)==22]  
                            wind_lines = np.array(windies)[:,:,wnames['idx']]
                        elif loc_ind==1: 
                            windies = [arr[wnames['idx']] for arr in wind_series[mm]]  
                            wind_lines = np.array(windies)[:,]
                        mean_wind_line = np.nanmean(wind_lines, axis=0)
                        
                        # plot winds
                        var1 = translate(np.nanmean(mean_wind_line[7:22]), -4,4)
                        dia.add_sample(np.abs(mean_lines[mm-1][-1]), var1,
                                       marker=dmarkers[yi], ms=8, ls='',
                                       mfc=month_colors[mm-1], mec='k',
                                       label=dnames[yi]+' '+str(mm)+' Winds')
                        
                        var2 = translate(np.nanmean(mean_wind_line[7:14]), -4,4)
                        dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                                       marker=dmarkers[yi], ms=8, ls='',
                                       mfc=month_colors[mm-1], mec=month_colors[mm-1],
                                       label=dnames[yi]+' '+str(mm)+' Winds*')
                        
                        var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dnames[yi]+' Winds', 'marker':dmarkers[yi]})
                        
                    # sst data
                    if 'sst' in var_group:
                        SSTs = np.load(path1+'sst/'+'tseries_'+str(yi)+'-'+str(mm)+'.npy')
                        mean_line = np.nanmean(SSTs, axis=0)
                        slope1, b, r, p, se = linregress(xxx[7:22], mean_line[7:22])
                        slope2, b, r, p, se = linregress(xxx[7:14], mean_line[7:14])
                        
                        # plot sst
                        var1 = translate(slope1, -0.04,0.04)
                        dia.add_sample(np.abs(mean_lines[mm-1][-1]), var1,
                                       marker=['P','D'][yi], ms=8, ls='',
                                       mfc=month_colors[mm-1], mec='k',
                                       label=dnames[yi]+' '+str(mm)+' SST')
                        var2 = translate(slope2, -0.04,0.04)
                        dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                                       marker=['P','D'][yi], ms=8, ls='',
                                       mfc=month_colors[mm-1], mec=month_colors[mm-1],
                                       label=dnames[yi]+' '+str(mm)+' SST*')
                        
                        var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dnames[yi]+' SST', 'marker':['P','D'][yi]})
    
                    # air temperature data
                    if 't2m' in var_group:
                        t2m_lines = [np.array(s)-273.15 for s in t2m_series[mm] if len(s)==22]
                        mean_temp_line = np.nanmean(t2m_lines, axis=0)
                        slope1, b, r, p, se = linregress(xxx[7:22], mean_temp_line[7:22])
                        slope2, b, r, p, se = linregress(xxx[7:14], mean_temp_line[7:14])
                        print(slope1, slope2)
                        # plot temperature
                        var1 = translate(slope1, -0.4,0.4)
                        dia.add_sample(np.abs(mean_lines[mm-1][-1]), var1,
                                       marker=['o','v'][yi], ms=8, ls='',
                                       mfc=month_colors[mm-1], mec='k',
                                       label=dnames[yi]+' '+str(mm)+' SST')
                        var2 = translate(slope2, -0.4,0.4)
                        dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                                       marker=['o','v'][yi], ms=8, ls='',
                                       mfc=month_colors[mm-1], mec=month_colors[mm-1],
                                       label=dnames[yi]+' '+str(mm)+' SST*')
                    
                    var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dnames[yi]+' Air Temp', 'marker':['o','v'][yi]})

    
                    # upper ocean data   
                    if yi==1 and 'ocn_prof' in var_group:
                        ocn_profs = ocn_data[mm]
                        prof_diffs = [profs[-1]-profs[0] for profs in ocn_profs]
                        upper_diffs = [np.nanmean(prof[0:8]) for prof in prof_diffs] #0:8->upper10m
                        
                        # plot upper ocean
                        var2 = translate(np.nanmean(upper_diffs), -0.25,0.25)
                        dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                                       marker='*', ms=10, ls='',
                                       mfc=month_colors[mm-1], mec=month_colors[mm-1],
                                       label=dnames[yi]+' '+str(mm)+' SST*')
                        
                        var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dnames[yi]+' Upper Ocean', 'marker':'*'})
                    
                    # wave data
                    if 'waves' in var_group:
                        swh_lines = [np.array(s) for s in swh_series[mm] if len(s)==22]
                        mean_wave_line = np.nanmean(swh_lines, axis=0)
                        
                        # plot waves
                        var1 = translate(np.nanmean(mean_wave_line[7:22]), 0,4)
                        dia.add_sample(np.abs(mean_lines[mm-1][-1]), var1,
                                       marker=['v', '^'][yi], ms=8, ls='',
                                       mfc=month_colors[mm-1], mec='k',
                                       label=dnames[yi]+' '+str(mm)+' Waves')
                        
                        var2 = translate(np.nanmean(mean_wave_line[7:14]), 0,4)
                        dia.add_sample(np.abs(mean_lines[mm-1][14]), var2,
                                       marker=['v', '^'][yi], ms=8, ls='',
                                       mfc=month_colors[mm-1], mec=month_colors[mm-1],
                                       label=dnames[yi]+' '+str(mm)+' Waves*')
                        
                        var_entries.append({'ls':'-', 'lw':0,'color':c, 'mfc':c,'mec':c,'label':dnames[yi]+' Waves', 'marker':dmarkers3[yi]})
                            
        # month legend 
        month_entries = []
        for mgroup in MONTHS:
            for mm in mgroup:
                c1=month_colors[mm-1]
                month_entries.append({'ls':'-', 'lw':3, 'color':c1,'mfc':c1,'mec':c1,
                                      'label':calendar.month_abbr[mm], 'marker':None})
        new_legend(dia, month_entries, bbox_to_anchor=(-0.075,1.0))
        
        # variable legend
        new_legend(dia, var_entries, bbox_to_anchor=(1.33,1.0))
        
        # time legend
        time_entries = []
        for outline, lname in zip(['k','gray'],['Long-term','Short-Term']):
            time_entries.append({'ls':'-', 'lw':0, 'color':'gray','mfc':'gray','mec':outline,
                                'label':lname, 'marker':'p'})
        new_legend(dia, time_entries, bbox_to_anchor=(-0.075,0.5), ncol=1)




#%% end
