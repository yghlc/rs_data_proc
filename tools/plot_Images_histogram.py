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


def get_max_min_histogram_percent(bin_edges, hist, min_percent=0.01, max_percent=0.99):
    '''
    get the max and min when cut of % top and bottom pixel values
    :param bin_edges:
    :param hist: numpy 1-d array.
    :param min_percent: percent
    :param max_percent: percent
    :return:
    '''

    if hist.ndim != 1:
        raise ValueError('Only accept one dimension array')

    if min_percent is not None and max_percent is not None:
        if min_percent >= max_percent:
            raise ValueError('min_percent >= max_percent')

    found_min = 0
    found_max = 0

    count = hist.size
    sum = np.sum(hist)
    accumulate_sum = 0
    for ii in range(count):
        accumulate_sum += hist[ii]
        if accumulate_sum/sum >= min_percent:
            found_min = bin_edges[ii]
            break

    accumulate_sum = 0
    for ii in range(count-1,0,-1):
        # print(ii)
        accumulate_sum += hist[ii]
        if accumulate_sum / sum >= (1 - max_percent):
            found_max = bin_edges[ii]
            break

    return found_min, found_max

def get_max_min_histogram_percent_allBands(bin_edges, hist_allBands, min_percent=0.01, max_percent=0.9):
    if min_percent is None or max_percent is None:
        return [],[]
    min_list = []
    max_list = []
    for hist in hist_allBands:
        min, max = get_max_min_histogram_percent(bin_edges, hist, min_percent=min_percent, max_percent=max_percent)
        min_list.append(min)
        max_list.append(max)
    return min_list, max_list

def main(options, args):

    img_dir_or_path = args[0]
    if os.path.isdir(img_dir_or_path):
        pre_name = os.path.basename(img_dir_or_path)
    else:
        pre_name = os.path.splitext(os.path.basename(img_dir_or_path))[0]
    if options.extent_shp is not None:
        pre_name +=  '_' + os.path.splitext(os.path.basename(options.extent_shp))[0]

    # range = (1, 6000)
    range = (options.value_range_min, options.value_range_max)
    bin_count = options.bin_count
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
        if options.planet_geojson is True:
            tif_list = get_Planet_image_tif_list(img_dir_or_path, options.extent_shp)
        elif os.path.isdir(img_dir_or_path):
            tif_list = io_function.get_file_list_by_ext('.tif', img_dir_or_path, bsub_folder=True)
        else:
            tif_list = [img_dir_or_path]
        if len(tif_list) < 1:
            if options.extent_shp is None:
                print('no image files in %s' % img_dir_or_path)
            else:
                print('no image files in %s within extent of %s' % (img_dir_or_path, options.extent_shp))
            return True

        # save tif list for checking
        io_function.save_list_to_txt(pre_name + '_tif.txt', tif_list)

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


    # find min and max grey values
    min_percent = options.hist_min_percent
    max_percent = options.hist_max_percent
    if min_percent is None or max_percent is None:
        pass
    else:
        min_list, max_list = get_max_min_histogram_percent_allBands(bin_edges,hist_allImg_list,min_percent=min_percent, max_percent=max_percent)
        min_max_info = 'min_%.3f_max_%.3f'%(min_percent,max_percent)
        save_min_max_txt = pre_name + '_' + hist_info + '_' + min_max_info + '.txt'
        with open(save_min_max_txt, 'w') as f_obj:
            f_obj.writelines('%s \n'%save_min_max_txt)
            f_obj.writelines('min (%.3f) and max (%.3f) values based on percentage on histogram \n'%(min_percent,max_percent))
            for min, max in zip(min_list, max_list):
                f_obj.writelines('min, max: %.4f, %.4f \n'%(min,max))

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

    ## marked the min and max values
    if min_percent is not None and max_percent is not None:
        for idx, (min, max) in enumerate(zip(min_list, max_list)):
            # print(idx, min, max)
            ax1.axvline(x=min, color=line_style[idx][0], linestyle=line_style[idx][1:], linewidth=0.5)
            ax1.axvline(x=max, color=line_style[idx][0], linestyle=line_style[idx][1:], linewidth=0.5)
            # ax1.text(dem + 100, 0.5, str(dem), rotation=90, fontsize=10, color='r')


    # plt.gcf().subplots_adjust(bottom=0.15)  #add space for the buttom
    # plt.gcf().subplots_adjust(top=0.8)  # the value range from [0,1], 1 is toppest, 0 is bottom
    # plt.gcf().subplots_adjust(left=0.15)
    # plt.gcf().subplots_adjust(right=0.15)

    ax1.legend(fontsize=10, loc="best")  # loc="upper left"
    ax1.set_xlabel('grey values')
    ax1.set_ylabel('%')

    # fig.legend((line_slope,line_dem),('Slope Histogram', 'DEM Histogram'))
    # ax1.legend(line_dem, ('Gray value'), fontsize=16,loc='upper center')

    if min_percent is not None and max_percent is not None:
        save_path = pre_name  + '_' + hist_info + '_' + min_max_info + '.jpg'
    else:
        save_path = pre_name + '_' + hist_info  + '.jpg'
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

    parser.add_option("-u", "--hist_max_percent",
                      action="store", dest="hist_max_percent",type=float,
                      help="the upper percent for choosing the max pixel value")

    parser.add_option("-l", "--hist_min_percent",
                      action="store", dest="hist_min_percent",type=float,
                      help="the lower percent for choosing the max pixel value")

    parser.add_option("-b", "--bin_count",
                      action="store", dest="bin_count",type=int, default=500,
                      help="the bin count of histogram")

    parser.add_option("", "--value_range_min",
                      action="store", dest="value_range_min",type=float,default=1,
                      help="the value range (min) of histogram")

    parser.add_option("", "--value_range_max",
                      action="store", dest="value_range_max",type=float,default=6000,
                      help="the value range (max) of histogram")

    #  range = (1, 6000)



    (options, args) = parser.parse_args()
    # print(options.planet_geojson)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

