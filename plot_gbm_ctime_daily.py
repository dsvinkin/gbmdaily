import os
from io import StringIO
import re
import numpy as np

from astropy.table import Table

import path_utils

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True 
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True 
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['axes.facecolor'] = 'w'  
mpl.rcParams['axes.edgecolor'] = 'k' 

mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['axes.titlesize'] = 18

mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.minor.size'] = 3

mpl.rcParams['figure.figsize'] = (8.27, 11.69) #A4

from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator

import matplotlib.pyplot as pl

class Plot:

    def __init__(self, n_panels):

        self.fig = pl.figure()
        self.lst_axis = []

        self._set_layout(n_panels)

    def _set_layout(self, n_panels):

        # Set panel parameters
        top, left, width = 0.1, 0.15, 0.7
        
        heigt_tot = 0.8
        heigt_channels = heigt_tot/n_panels

        left_ch_names = 0.8
        
        bottom = 1.0 - top - n_panels * heigt_channels
        
        # rect [left, bottom, width, height] 
        lst_rect = []
        for i in range(n_panels):
            lst_rect.append([left, bottom + (n_panels-i-1) * heigt_channels, width, heigt_channels])

        for i in range(len(lst_rect)):
            if i > 0:
                self.lst_axis.append(self.fig.add_axes(lst_rect[i],  sharex=self.lst_axis[0]))
            else:
                self.lst_axis.append(self.fig.add_axes(lst_rect[i]))

    def get_delta_y(self, y_min, y_max):

        dy = y_max - y_min

        lst_num = [20, 50, 100]

        delta_y = lst_num[0]
        n_ticks = y_max / delta_y

        i = 0
        while(n_ticks > 4):
            for n in lst_num:
               delta_y = int(n * 10**i)
               n_ticks = dy // delta_y + 2
               if n_ticks <=4:
                   break
            i = i + 1

        return delta_y

    def set_x_minor_ticks(self):
        """
        ax1.xaxis.set_minor_locator(minorLocator_x)
        ax1.yaxis.set_minor_locator(minorLocator_y_sum)  # y minor ticks
        ax1.set_xlim(arr_begin_end[0], arr_begin_end[1])
        """
        pass

    def add_plot(self, 
        idx_panel,
        arr_ti, 
        arr_counts, 
        bg_level,  
        arr_begin_end, 
        arr_vlines,
        text
    ):
        ax = self.lst_axis[idx_panel]
        ax.plot(arr_ti, arr_counts, drawstyle='steps-post',color='k', linewidth=0.5)

        if bg_level is not None:
            ax.axhline(bg_level, color='k', linestyle ='--', linewidth=0.5)

        arr_bool = np.logical_and(arr_ti >= arr_begin_end[0], arr_ti <= arr_begin_end[1])

        cnt_max = int(max(arr_counts[arr_bool]))
        cnt_min = int(min(arr_counts[arr_bool]))
        
        delta_y = self.get_delta_y(cnt_min, cnt_max)

        if idx_panel > 0:
            y_max = cnt_max // delta_y * delta_y + 3/2 * delta_y
        else:
            y_max = cnt_max // delta_y * delta_y + delta_y

        y_min = cnt_min // delta_y * delta_y
    
        # рисуем временной интервал
        #ax.vlines(arr_vlines, [y_min,y_min], [y_max,y_max], linestyles='dashed', color='k', linewidth=0.5)
    
        ax.set_yticks(np.arange(y_min, max(arr_counts) + delta_y, delta_y))
        ax.set_ylim(y_min, y_max)
        ax.set_xlim(arr_begin_end[0], arr_begin_end[1])

        ax.text(0.9, 0.85, text, horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=10)
        
        minorLocator_y = MultipleLocator(delta_y/2.0)
        ax.yaxis.set_minor_locator(minorLocator_y)  # y minor ticks

        if idx_panel < len(self.lst_axis) -1:
            pl.setp(ax.get_xticklabels(), visible=False)

    def save(self, file_name):
        self._set_xlabel()
        self._set_ylabel()
        pl.savefig(file_name, format='pdf')

    def set_xticks(self, step):
        for ax in self.lst_axis:
            ax.xaxis.set_major_locator(MultipleLocator(step))
            ax.xaxis.set_minor_locator(MultipleLocator(step/2.0))

    def set_title(self, text):
        self.fig.suptitle(text, fontsize=mpl.rcParams['axes.titlesize'])

    def _set_xlabel(self):
        self.lst_axis[-1].set_xlabel(r'T-T$_{0}$ (s)')

    def _set_ylabel(self, text='Counts'):
        """
        See
        https://stackoverflow.com/questions/8334164/matplotlib-single-axis-label-for-multi-panel-plot
        """
        ax = self.fig.add_axes( [0., 0., 1, 1] )
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.text( 
            .05, 0.5, text, rotation='vertical',
            horizontalalignment='center', verticalalignment='center',
            fontsize=mpl.rcParams['axes.labelsize']
        )

def get_ebounds():

    e_bounds = "4  11  26  50  101 290  538   997".split()
    dic_e = {}
    for i, e in enumerate(e_bounds[1:], 1):    
       dic_e["Ch{:d}".format(i)] = "{:s}-{:s} keV".format(e_bounds[i-1], e_bounds[i])

    return dic_e

def plot_ctime(tab, date, time, det, Ti, Tf, plot_name):

    lst_chan = tab.colnames[4:]
    n_chan = len(lst_chan)
    x_range = (Ti, Tf)
    #x_ticks = int((Ti-Tf)/5)

    e_bounds = get_ebounds()

    plot = Plot(n_chan)
    for i, ch in enumerate(lst_chan):
        #print(i, ch)
        #continue
        plot.add_plot(i, tab['Ti'], tab[ch], None, x_range, None, e_bounds[ch])

    str_title = "Fermi-GBM detector {:s}\n{:s} T$_0$={:s}".format(det, date, time)
    plot.set_title(str_title)
    plot.set_xticks(250)
    plot.save(plot_name)


if __name__ == '__main__':

    # 20190902 11370.000 s UT (03:09:30.000)
  
    date = '20190902'
    time = '11370.000 s UT (03:09:30.000)'
    
    lst_det = "n0 n1 n2 n3 n4 n5 n6 n7 n8 n9 na nb".split()
    resolution = '30s'
    ver = '01'

    Ti, Tf = -500.0, 500.0

    path = './' + date

    for det in lst_det:
        file_name = "glg_ctime_{:s}_{:s}_v{:s}_{:s}.txt".format(det, date[2:], ver, resolution)
        plot_name = "glg_ctime_{:s}_{:s}_v{:s}_{:s}.pdf".format(det, date[2:], ver, resolution)

        tab = Table.read(os.path.join(path, file_name), format='ascii')
        print(tab.colnames)

        arr_bool = np.logical_and(tab['Ti'] >= Ti, tab['Tf'] <= Tf)

        plot_ctime(tab[arr_bool], date, time, det, Ti, Tf, os.path.join(path, plot_name))