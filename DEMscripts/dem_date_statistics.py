#!/usr/bin/env python
# Filename: dem_date_statistics 
"""
introduction: check the dem acquistion date

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 27 April, 2021
"""

import os,sys

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.timeTools as timeTools
import time
import vector_gpd

from dem_common import arcticDEM_reg_tif_dir


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FuncFormatter

shp_dir = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Strip_Index_Rel7')

global_bin_size = 1   # remember to change this one
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.

    # print(global_bin_size)
    # s = str(100 * y*bin_size)
    s = "%.0f"%(100 * y*global_bin_size)

    # The percent symbol needs escaping in latex
    # if matplotlib.rcParams['text.usetex'] is True:
    #     return s + r'$\%$'
    # else:
    #     return s + '%'
    if matplotlib.rcParams['text.usetex'] is True:
        return s
    else:
        return s

def get_file_list(b_urls=False):
    # get file name from urls
    if b_urls:
        t0 = time.time()
        ArcticDEM_strip_index = os.path.join(shp_dir, 'ArcticDEM_Strip_Index_Rel7.shp')
        # ref_raster='grid_20km_bin.tif'
        file_list = vector_gpd.read_attribute_values_list(ArcticDEM_strip_index,'fileurl')
        file_list = [ item for item in file_list if item.startswith('/mnt') is False ]      # remove invalid urls
        print('read ArcticDEM strip polygons, cost: ', time.time() - t0, ' seconds')
        print('Get %d dem urls ' % len(file_list))
    else:
        file_list = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir,'*_dem_reg.tif')
        print('Get %d dem_reg.tif from %s'%(len(file_list), arcticDEM_reg_tif_dir))

    return file_list

def plot_hist(value_list,save_path,bins = None, b_xtick_rotate=False, x_tick=None,x_label=None):
    # plot a histogram
    # bin_count = 12
    if bins is None:
        bins = np.arange(0.5, 12.5+1, 1)
    print('bins:', bins)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    n, bins, patches = ax.hist(value_list, bins=bins, alpha=0.75, ec="black", linewidth='1.5',
                               color='grey', hatch='', rwidth=1,density = True)  # density = True, # label=labels,

    # ax.legend(prop={'size': 12})
    if x_tick is not None:
        # ticks = None, labels = None
        plt.xticks(ticks=x_tick,labels=x_tick)
    ax.tick_params(axis='both', which='both', direction='out', length=7, labelsize=20)  # ,width=50 #,
    # if xlabelrotation is not None:
    #     ax.tick_params(axis='x', labelrotation=90)
    if b_xtick_rotate:
        ax.tick_params(axis='x', labelrotation=90)

    # if ylim is not None:
    #     ax.set_ylim(ylim)

    formatter = FuncFormatter(to_percent)
    # Set the formatter
    plt.gca().yaxis.set_major_formatter(formatter)

    if x_label is not None:
        plt.xlabel(x_label,fontsize=20)
    plt.ylabel('Percentage (%)',fontsize=20)

    plt.gcf().subplots_adjust(bottom=0.15)
    # plt.grid(True)
    plt.savefig(save_path)  #
    print('save fig to ', os.path.abspath(save_path))

def main():
    file_list = get_file_list(b_urls=True)

    year_dates = [ timeTools.get_yeardate_yyyymmdd(os.path.basename(item),pattern='[0-9]{8}_') for item in file_list]
    month_list = [item.month for item in year_dates]
    year_list = [item.year for item in year_dates]

    month_unique = sorted(set(month_list))
    print('month_unique:',month_unique)

    year_unique = sorted(set(year_list))
    print('year_unique:',year_unique)

    # save unique date to txt file
    dates_unique = set(year_dates)
    dates_unique = sorted(dates_unique)
    dates_unique_str = [ timeTools.date2str(item,'%Y-%m-%d')  for item in dates_unique]
    io_function.save_list_to_txt('dates_unique.txt',dates_unique_str)

    plot_hist(month_list, 'ArcticDEM_strip_month_hist.jpg',x_tick=month_unique,x_label='Month')

    year_bins = np.arange(min(year_list)-0.5, max(year_list)+0.5+1, 1)
    plot_hist(year_list, 'ArcticDEM_strip_year_hist.jpg',bins=year_bins,b_xtick_rotate=True,x_tick=year_unique,x_label='Year')




if __name__ == '__main__':
    main()
    pass
