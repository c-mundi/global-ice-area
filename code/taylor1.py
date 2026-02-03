#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 2025
taylor1.py

--> example taylor diagram code

https://gist.github.com/ycopin/3342888
https://zenodo.org/records/5548061

@author: mundi
"""
#%% start
import numpy as np
import matplotlib.pyplot as plt
import easy_mpl 
from easy_mpl import taylor_plot

month_colors = ['#238443','#78c679','#c2e699',
          '#d7b5d8','#df65b0','#dd1c77','#980043',
          '#7a0177',
          '#253494','#2c7fb8','#41b6c4','#a1dab4']

#%% sample data

         # J   F   M    A     M    J     J     A    S     O    N    D
impacts = [0, 0.1, 0.1, -0.4,-0.5,-0.6, -0.7, -0.8, 0.4, 0.5, 0.6, 0.7]


var_d = {'winds':[1]*12,
         'sst':[-0.1]*6+[0.1]*6,
         'waves':np.arange(1,12+1)
    }

figure = taylor_plot(observations=np.array(impacts),
                    simulations=var_d,
                    plot_bias=True,
                    cont_kws={'colors': '#1A74A5', 'linewidths': 1.5},
                    grid_kws={'axis': 'x', 'color': '#5CB994', 'lw': 2.0, 'ls': 'dotted'},             
                    marker_kws={'markersize': 10, 'markeredgewidth': 1.0,
                                'markeredgecolor': 'black', 'lw': 0.0},
                    title="Sample Taylor Plot")

#%% editable code
# https://gist.github.com/ycopin/3342888

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
        ax.axis["top"].label.set_text("Change in Variable")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("Normalized Storm Impact")

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
    
    
#%%% tests

#%%%% test2
def test2():
    """
    Climatology-oriented example (after iteration w/ Michael A. Rawlins).
    """

    # Reference std
    stdref = 48.491

    # Samples std,rho,name
    samples = [[25.939, 0.385, "Model A"],
               [29.593, 0.509, "Model B"],
               [33.125, 0.585, "Model C"],
               [29.593, 0.509, "Model D"],
               [71.215, 0.473, "Model E"],
               [27.062, 0.360, "Model F"],
               [38.449, 0.342, "Model G"],
               [35.807, 0.609, "Model H"],
               [17.831, 0.360, "Model I"]]

    fig = plt.figure()

    dia = TaylorDiagram(stdref, fig=fig, label='Reference', extend=True)
    # dia.samplePoints[0].set_color('r')  # Mark reference point as a red star

    # Add models to Taylor diagram
    for i, (stddev, corrcoef, name) in enumerate(samples):
        dia.add_sample(stddev, corrcoef,
                       marker='$%d$' % (i+1), ms=10, ls='',
                       mfc='k', mec='k',
                       label=name)

    # Add RMS contours, and label them
    # contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    # plt.clabel(contours, inline=1, fontsize=10, fmt='%.0f')

    dia.add_grid()                                  # Add grid
    dia._ax.axis[:].major_ticks.set_tick_out(True)  # Put ticks outward

    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ],
               numpoints=1, prop=dict(size='small'), loc='upper right')
    fig.suptitle("Taylor diagram", size='x-large')  # Figure title

    return dia

dia = test2()

#%%%% my test


         # J   F   M    A     M    J     J     A    S     O    N    D
impacts = [0, 0.1, 0.2, -0.4,-0.5,-0.6, -0.7, -0.8, 0.4, 0.5, 0.6, 0.7]


var_d = {'winds':[1]*12,
         'sst':[-0.1]*6+[0.1]*6,
         'waves':np.arange(1,12+1)
    }

fig = plt.figure()

dia = TaylorDiagram(0.33, fig=fig, label='Reference', extend=True)

# Add models to Taylor diagram
for mm, color in zip([0,1,2],['#238443','#78c679','#c2e699']):
    for i, key in enumerate(var_d.keys()):
        dia.add_sample(impacts[mm],var_d[key][mm],
                       marker='$%d$' % (i+1), ms=10, ls='',
                       mfc=color, mec=color,
                       label=key)
       
dia.add_grid()                                  # Add grid
dia._ax.axis[:].major_ticks.set_tick_out(True)  # Put ticks outward

# Add a figure legend and title
fig.legend(dia.samplePoints,
           [ p.get_label() for p in dia.samplePoints ],
           numpoints=1, prop=dict(size='small'), loc='upper right')
fig.suptitle("Taylor diagram", size='x-large')  # Figure title


#%% why not just a regular plot?
import functions as fx

nh_path = '/Users/mundi/Desktop/month-hemi/nh_data/'
sh_path = '/Users/mundi/Desktop/month-hemi/sh_data/'
root_paths = [nh_path, sh_path]

month_colors = ['#238443','#78c679','#c2e699',
                '#d7b5d8','#df65b0','#dd1c77','#980043','#7a0177',
                '#253494','#2c7fb8','#41b6c4','#a1dab4'
                ]
months = np.arange(1,12+1)

decades = [np.arange(1982, 1992), np.arange(2010,2020)]
decade_names = ['Early Satellite Era ', 'Present Day ']
dnames = ['80s','10s']
dmarkers = ['s','o']

hemi_colors = ['#01665e','#8c510a']
hemi_names= ['Arctic', 'Antarctic']

decade_colors = ['#2166ac', '#b2182b']
shade_colors = ['#67a9cf', '#ef8a62']

quad_colors = ['#f7f7f7','#cccccc','#969696','#525252']

widx = 0
wnames0 = {'idx':1, 'ylims':[-6,6], 'name':'Zonal', 
          'xlims':[[-4.5,3],[-1,3]], 'vline':0}
wnames1 = {'idx':2, 'ylims':[-6,6], 'name':'Meridional', 
          'xlims':[[-4.5,2],[-1,2]], 'vline':0}
wnames2 = {'idx':0, 'ylims':[-2,14], 'name':'Total', 
          'xlims':[[4,10],[6,10]], 'vline':None}
wnames = [wnames0, wnames1, wnames2][widx]

for months in [[4,5,6,7,8]]:

    fig,ax = plt.subplots(1,1)
    fig.suptitle(str(months), size='large')  # Figure title
    
    for loc_ind, loc in enumerate(hemi_names[0:1]):
        for yi, years in enumerate(decades):
            yr_title = str(years[0])+'-'+str(years[-1])
            print(yr_title)
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
                
                print(mm, mean_lines[mm-1][-1], mean_wind_line[-1])
                ax.plot(np.abs(mean_lines[mm-1][-1]), mean_wind_line[-1],
                        marker=dmarkers[yi], ms=8, ls='',
                        mfc=month_colors[mm-1], mec='k',
                        label=dnames[yi]+' '+str(mm)+' Winds')
       
    # Legend
    ax.legend(prop=dict(size='small'), loc='upper right')


#%% end
