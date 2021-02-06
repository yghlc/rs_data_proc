#!/usr/bin/env python
# Filename: plot_planetImg_histogram 
"""
introduction: plot histogram of all the Planet images

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 07 January, 2021
"""

import os, sys
from optparse import OptionParser

HOME = os.path.expanduser('~')
# path of DeeplabforRS
codes_dir2 =  HOME +'/codes/PycharmProjects/DeeplabforRS'
sys.path.insert(0, codes_dir2)

import basic_src.io_function as io_function
import basic_src.basic as basic
import raster_io


code_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',)
sys.path.insert(0, os.path.join(code_dir,'planetScripts'))

# sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/ChangeDet_DL/dataTools'))
from get_planet_image_list import  get_Planet_SR_image_list_overlap_a_polygon
from get_planet_image_list  import read_extent_shapefile_epgs4326

import numpy as np
import matplotlib.pyplot as plt



def np_histogram_one_image(tif_path, axis_range, bin_count):

    # read tif
    img_data_allBands, nodata = raster_io.read_raster_all_bands_np(tif_path)  # (bandcount, height, width)
    nband, width, height = img_data_allBands.shape
    print('nband, width, height', nband, width, height)
    data = img_data_allBands.reshape(nband,-1)
    print('data.shape', data.shape)


    hist_list = []
    bin_edges = None
    for index in range(nband):
        one_band_data = data[index,:]
        if nodata is not None:
            one_band_data = one_band_data[ one_band_data!= nodata]  # remove nodata values
        hist, bin_edges = np.histogram(one_band_data,bins=bin_count,density=False,range=axis_range)
        hist_list.append(hist)
        # print(bin_edges)

    return hist_list, bin_edges

def get_Planet_image_tif_list(image_dir, ext_shp):

    # try to get all geojson under image_dir
    geojson_list = io_function.get_file_list_by_ext('.geojson',image_dir,bsub_folder=True)
    # remove incomplete scenes
    geojson_list = [item for item in geojson_list if 'incomplete_scenes' not in item]
    if len(geojson_list) < 1:
        raise ValueError('There is no geojson files in %s'%image_dir)

    extent_polygon = read_extent_shapefile_epgs4326(ext_shp)
    if len(extent_polygon) !=1:
        raise ValueError('The exetent is not one')
    extent = extent_polygon[0]

    cloud_cover_thr = 100
    image_path_list, cloud_cover_list = get_Planet_SR_image_list_overlap_a_polygon(extent, geojson_list, cloud_cover_thr)

    return image_path_list


def main(options, args):

    img_dir = args[0]
    if options.planet_geojson is True:
        tif_list = get_Planet_image_tif_list(img_dir, options.extent_shp)
    else:
        tif_list = io_function.get_file_list_by_ext('.tif',img_dir,bsub_folder=True)
    if len(tif_list) < 1:
        if options.extent_shp is None:
            print('no image files in %s'%img_dir)
        else:
            print('no image files in %s within extent of %s'%(img_dir,options.extent_shp))
        return True

    # tif_list = tif_list[:1] # for test

    pre_name = os.path.basename(img_dir)
    if options.extent_shp is not None:
        pre_name +=  '_' + os.path.splitext(os.path.basename(options.extent_shp))[0]

    # save tif list for checking
    io_function.save_list_to_txt(pre_name+'_tif.txt',tif_list)

    range = (1, 6000)
    bin_count = 500
    hist_allImg_list = []
    bin_edges = None

    hist_info = 'bin%d_range%d_%d'%(bin_count,range[0],range[1])

    for idx in np.arange(1000):  # assume we have 1000 bands
        hist_txt = pre_name + '_hist_b%d' % idx + '_' + hist_info + '.txt'
        if os.path.isfile(hist_txt):
            hist = np.loadtxt(hist_txt)
            hist_allImg_list.append(hist)
        else:
            break
        _, bin_edges = np.histogram([1,2,3], bins=bin_count,range=range)

    if len(hist_allImg_list) < 1:
        # read image data and get histograms for each band
        for idx, tif in enumerate(tif_list):
            print(idx,len(tif_list), tif)
            hist_list, bin_edges = np_histogram_one_image(tif, range, bin_count)

            bin_edges = bin_edges # we have the fix range, so the bin_edges should be fix
            if len(hist_allImg_list) < 1:
                hist_allImg_list = hist_list
                # print(hist_allImg_list)
            else:
                # accumulate the hist
                if len(hist_allImg_list) != len(hist_list):
                    raise ValueError('image band count of %s is different from the previous ones'%tif)
                for hist_all, hist in zip(hist_allImg_list,hist_list):
                    hist_all += hist
                    # print(hist_all)


    #save hist_allImg_list and bin_edges to txt file. Calculaing hist from over 600 GB images is time consuming.
    for idx, hist in enumerate(hist_allImg_list):
        hist_txt = pre_name + '_hist_b%d'%idx + '_' + hist_info + '.txt'
        # data = np.stack(bin_edges, hist)
        # print('data shape', data.shape)
        np.savetxt(hist_txt, hist,fmt='%d')


    fig = plt.figure(figsize=(6,4)) #
    ax1 = fig.add_subplot(111)
    # ax2 = ax1.twiny()    #have another x-axis

    # ax1.yaxis.tick_right()

    # line_list = []
    line_style = ['b-','g-','r-','r-.']
    values_x = bin_edges[:-1]  # ignore the last one
    for idx, hist in enumerate(hist_allImg_list):
    # plot histogram
        values_per = 100.0 * hist / np.sum(hist) # draw the percentage
        line, = ax1.plot(values_x, values_per,line_style[idx], label="Band %d"%(idx+1), linewidth=0.8)

    # ax2.set_xlabel("Elevation (m)",color="red",fontsize=15)
    # ax2.spines['bottom'].set_color('blue')
    # ax1.spines['top'].set_color('red')
    # ax2.xaxis.label.set_color('blue')
    ax1.tick_params(axis='x')

    # ## marked the values
    # threshold = [9500]
    # for dem in threshold:
    #     ax1.axvline(x=dem, color='r', linewidth=0.8, linestyle='--')
    #     ax1.text(dem + 100, 0.5, str(dem), rotation=90, fontsize=10, color='r')


    # plt.gcf().subplots_adjust(bottom=0.15)  #add space for the buttom
    # plt.gcf().subplots_adjust(top=0.8)  # the value range from [0,1], 1 is toppest, 0 is bottom
    # plt.gcf().subplots_adjust(left=0.15)
    # plt.gcf().subplots_adjust(right=0.15)

    ax1.legend(fontsize=10, loc="best")  # loc="upper left"
    ax1.set_xlabel('grey values')
    ax1.set_ylabel('%')

    # fig.legend((line_slope,line_dem),('Slope Histogram', 'DEM Histogram'))
    # ax1.legend(line_dem, ('Gray value'), fontsize=16,loc='upper center')

    save_path = pre_name  +'_bin%d'%bin_count + '_range_%d_%d'%(range[0],range[1]) + '.jpg'
    plt.savefig(save_path, dpi=200)  # 300

    pass

if __name__ == '__main__':
    usage = "usage: %prog [options] image_folder "
    parser = OptionParser(usage=usage, version="1.0 2021-01-07")
    parser.description = 'Introduction: plot histogram of images '

    parser.add_option("-e", "--extent_shp",
                      action="store", dest="extent_shp",
                      help="the path for extent, shapefile")

    parser.add_option("-p", "--planet_geojson",
                      action="store_true", dest="planet_geojson",default=False,
                      help="indicate that Planet images with geojson in the input folder")



    (options, args) = parser.parse_args()
    # print(options.planet_geojson)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

