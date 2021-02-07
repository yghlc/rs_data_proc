#!/usr/bin/env python
# Filename: plot_Images_histogram_test.py 
"""
introduction:
# run "pytest plot_Images_histogram_test.py  " or "pytest " for test, add " -s for allowing print out"
# "pytest can automatically search *_test.py files "

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 07 February, 2021
"""

import os,sys
import numpy as np

# the path of Landuse_DL
# code_dir = os.path.expanduser('~/codes/PycharmProjects/Landuse_DL')
code_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')   # get path of
print(code_dir)
sys.path.insert(0, code_dir)

import plot_Images_histogram

def test_get_max_min_histogram_percent():
    extent_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/Willow_River/extent/WR_extent_latlon.shp')
    img_dir = os.path.expanduser('~/Data/Arctic/canada_arctic/rsImages/planet_sr_images/2020_July_August')

    work_dir = os.path.expanduser('~/Data/Arctic/canada_arctic/rsImages/histograms')

    range = (1, 6000)
    bin_count = 500
    hist_allImg_list = []
    bin_edges = None

    pre_name = os.path.basename(img_dir)
    if extent_shp is not None:
        pre_name +=  '_' + os.path.splitext(os.path.basename(extent_shp))[0]

    hist_info = 'bin%d_range%d_%d'%(bin_count,range[0],range[1])
    print(hist_info)

    for idx in np.arange(1000):  # assume we have 1000 bands
        hist_txt = pre_name + '_hist_b%d' % idx + '_' + hist_info + '.txt'
        # print(hist_txt)
        hist_txt = os.path.join(work_dir, hist_txt)
        if os.path.isfile(hist_txt):
            hist = np.loadtxt(hist_txt)
            hist_allImg_list.append(hist)
        else:
            break
        _, bin_edges = np.histogram([1,2,3], bins=bin_count,range=range)



    min, max = plot_Images_histogram.get_max_min_histogram_percent(bin_edges, hist_allImg_list[0],min_percent=0.01, max_percent=0.99)
    print(min, max)



if __name__ == '__main__':

    pass