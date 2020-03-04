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

mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.minor.size'] = 4

mpl.rcParams['figure.figsize'] = (8.27, 11.69) #A4

from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator

import matplotlib.pyplot as plt

def plot_angles(tab, date, time):

    lst_det = tab.colnames[1:]

    for det in lst_det:
        # remove detectors with icncident angle > 70 deg
        if np.min(tab[det]) >=70:
            continue
        plt.plot(tab['T'], tab[det], label=det)

    plt.title("Fermi-GBM detector incident angles\n{:s} T$_0$={:s}".format(date, time))
    plt.xlabel("T-T$_0$ (s)")
    plt.ylabel("Incident angle (deg)")
    plt.xlim(-100,100)
    plt.ylim(0,80)

    ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(10))

    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_minor_locator(MultipleLocator(2))

    plt.legend(title='Detectors:')
    plt.savefig('angles.pdf')

if __name__ == '__main__':

    date = '20190902'
    time = '11370.0 s UT (03:09:30.0)'
    path = './' + date

    ang_file = path_utils.get_files(path, pattern='poshist_angles', prefix=True, all=False)
    tab_ang = Table.read(os.path.join(path, ang_file), format='ascii')
    #print(tab_ang.colnames)

    arr_bool = np.logical_and(tab_ang['T'] >= -100.0, tab_ang['T'] <= 100.0)
    plot_angles(tab_ang[arr_bool], date, time)
    

