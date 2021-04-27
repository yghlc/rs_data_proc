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

from dem_common import arcticDEM_reg_tif_dir


import numpy as np
import matplotlib.pyplot as plt

def main():

    file_list = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir,'*_dem_reg.tif')
    print('Get %d dem_reg.tif from %s'%(len(file_list), arcticDEM_reg_tif_dir))

    year_dates = [ timeTools.get_yeardate_yyyymmdd(item,pattern='[0-9]{8}_') for item in file_list]
    month_list = [item.month for item in year_dates]
    value_list = month_list

    # plot a histogram
    # bin_count = 12
    bins = np.arange(0, 12, 1)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    n, bins, patches = ax.hist(value_list, bins=bins, alpha=0.75, ec="black", linewidth='1.5',
                               color='grey', hatch='',  rwidth=1)  # density = True, # label=labels,

    # ax.legend(prop={'size': 12})
    plt.xticks(bins)
    ax.tick_params(axis='both', which='both', direction='out', length=7, labelsize=20)  # ,width=50 #,
    # if xlabelrotation is not None:
    #     ax.tick_params(axis='x', labelrotation=90)

    # if ylim is not None:
    #     ax.set_ylim(ylim)

    plt.gcf().subplots_adjust(bottom=0.15)
    # plt.grid(True)
    plt.savefig('ArcticDEM_strip_date_hist.jpg')  #





if __name__ == '__main__':
    main()
    pass
