#!/usr/bin/env python
# Filename: grey_image_segment.py 
"""
introduction: segment a grey image (8 bit)

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 01 March, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.io_function as io_function
import basic_src.basic as basic
import split_image

from image_segment import quickshift_segmentaion
from image_segment import mean_shift_segmentation

from skimage import measure

import multiprocessing
from multiprocessing import Pool
import numpy as np

import cv2

def get_stastics_from_array(in_array, nodata,range=None):
    data_1d = in_array.flatten()
    data_1d = data_1d[ data_1d != nodata]
    data_1d = data_1d[~np.isnan(data_1d)]  # remove nan value

    # current, only calculate the mean
    values = []
    values.append(np.mean(data_1d).astype(np.float64)) # if data_1d is empty, then it returns None.     # save to float64 for json
    values.append(np.std(data_1d).astype(np.float64))  # std
    values.append(data_1d.size)
    return values

def segment_a_patch(idx, patch, patch_count,img_path, org_raster,b_save_patch_label):

    print('tile: %d / %d' % (idx + 1, patch_count))
    image_name_no_ext = io_function.get_name_no_ext(img_path)
    patch_dir = image_name_no_ext + '_patch%d'%idx
    patch_label_path = os.path.join(patch_dir, image_name_no_ext+'_patch%d_label.tif'%idx)
    if b_save_patch_label is True:
        if os.path.isdir(patch_dir) is False:
            io_function.mkdir(patch_dir)
        if os.path.isfile(patch_label_path):
            print('%s exists, skip'%patch_label_path)
            return patch, patch_label_path, None, None

    # read imag
    one_band_img, nodata = raster_io.read_raster_one_band_np(img_path, boundary=patch)

    # # apply median filter to image (remove some noise)
    one_band_img = cv2.medianBlur(one_band_img, 3)  # with kernal=3, cannot accept int32

    # segmentation algorithm (the output of these algorithms is not alway good, need to chose the parameters carafully)
    # out_labels = watershed_segmentation(one_band_img)
    # out_labels = k_mean_cluster_segmentation(one_band_img)
    out_labels = quickshift_segmentaion(one_band_img,ratio=0.3, kernel_size=5, max_dist=10,
                           sigma=1, convert2lab=False)
    #
    #
    # out_labels = mean_shift_segmentation(one_band_img)

    # print('min and max labels of out_labels', np.min(out_labels), np.max(out_labels))

    if b_save_patch_label is True:
        # save the label
        raster_io.save_numpy_array_to_rasterfile(out_labels, patch_label_path, img_path, boundary=patch) # it copy nodata, need to unset it later
        return patch, patch_label_path, None, None

    # calculate the attributes based on orginal data for original data
    object_attributes = {}  # object id (label) and attributes (list)
    if org_raster is not None:
        org_img_b1, org_nodata = raster_io.read_raster_one_band_np(org_raster, boundary=patch)

        # get regions (the labels output by segmentation is not unique for superpixels)
        # regions = measure.regionprops(out_labels, intensity_image=org_img_b1)     # regions is based on out_labels, so it has the same issue.
        # print('region count from sk-image measure:',len(regions))

        label_list = np.unique(out_labels)
        # get statistics for each segmented object (label)
        for label in label_list:
            in_array = org_img_b1[ out_labels == label ]
            object_attributes[label] = get_stastics_from_array(in_array, org_nodata)

        return patch, out_labels, nodata, object_attributes

    return patch, out_labels, nodata, None

def segment_a_grey_image(img_path, save_dir,process_num, org_raster=None,b_save_patch_label=False):

    out_pre = os.path.splitext(os.path.basename(img_path))[0]
    label_path = os.path.join(save_dir, out_pre + '_label.tif')
    if os.path.isfile(label_path):
        basic.outputlogMessage('%s exist, skip segmentation'%label_path)
        return label_path

    height, width, band_num, date_type = raster_io.get_height_width_bandnum_dtype(img_path)
    print('input image: height, width, band_num, date_type',height, width, band_num, date_type)

    # if the original data is available, then calculate the attributes based on that
    if org_raster is not None:
        org_height, org_width, org_band_num, org_date_type = raster_io.get_height_width_bandnum_dtype(org_raster)
        if org_height != height or org_width != width:
            raise ValueError('%s and %s do not have the same size'%(img_path,org_raster))

    save_labes = np.zeros((height,width),dtype=np.int32)
    # divide the image the many small patches, then calcuate one by one, solving memory issues.
    image_patches = split_image.sliding_window(width,height, 1024, 1024,adj_overlay_x=0,adj_overlay_y=0)
    patch_count = len(image_patches)

    # for idx, patch in enumerate(image_patches):
    #     out_patch,out_labels = segment_a_patch(idx, patch, patch_count,img_path)
    #     # copy to the entire image
    #     row_s = patch[1]
    #     row_e = patch[1] + patch[3]
    #     col_s = patch[0]
    #     col_e = patch[0] + patch[2]
    #     save_labes[row_s:row_e, col_s:col_e] = out_labels

    theadPool = Pool(process_num)
    parameters_list = [ (idx, patch, patch_count,img_path, org_raster,b_save_patch_label) for idx, patch in enumerate(image_patches)]
    results = theadPool.starmap(segment_a_patch, parameters_list)

    patch_label_path_list = []
    patch_label_id_range = []
    object_attributes = {}  # object id (label) and attributes (list)
    for res in results:
        patch, out_labels, nodata, attributes = res
        if isinstance(out_labels, str) and os.path.isfile(out_labels):  #if it's a label file
            patch_label_path_list.append(out_labels)
        else:
            # copy to the entire image
            row_s = patch[1]
            row_e = patch[1] + patch[3]
            col_s = patch[0]
            col_e = patch[0] + patch[2]
            current_min = np.max(save_labes)
            print('current_max',current_min)
            patch_label_id_range.append(current_min)
            save_labes[row_s:row_e, col_s:col_e] = out_labels + current_min + 1
            if attributes is not None:
                update_label_attr = {}
                for key in attributes:
                    update_label_attr[ key + current_min] = attributes[key]
                # add to the attributes
                object_attributes.update(update_label_attr)

    theadPool.close()
    # # apply median filter (remove some noise), #  we should not use median  filter, because it's labels, not images.
    # label_blurs = cv2.medianBlur(np.float32(save_labes), 3)  # with kernal=3, cannot accept int32
    # # print(label_blurs, label_blurs.dtype)
    # save_labes = label_blurs.astype(np.int32)

    # return a list of labels saved in current working folder.
    if b_save_patch_label:
        return patch_label_path_list

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    # save attributes (if not empty)
    if object_attributes:
        attribute_path = os.path.join(save_dir, out_pre + '_attributes.txt')
        io_function.save_dict_to_txt_json(attribute_path,object_attributes)

    # save the label
    raster_io.save_numpy_array_to_rasterfile(save_labes, label_path, img_path) # do not set nodata
    # save id ranges to txt
    label_id_range_txt = os.path.splitext(label_path)[0] + '_IDrange.txt'
    patch_label_id_range = [str(item) for item in patch_label_id_range]
    io_function.save_list_to_txt(label_id_range_txt,patch_label_id_range)


    return label_path


def test_segment_a_grey_image():
    data_dir = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff')
    img_path = os.path.join(data_dir,'WR_extent_grid_ids_DEM_diff_grid9274_8bit.tif')
    # org_raster = os.path.join(data_dir,'WR_extent_grid_ids_DEM_diff_grid9274.tif')
    save_dir = os.path.join(data_dir,'segment_parallel_9274')
    process_num = 8
    b_save_patch_label = True
    segment_a_grey_image(img_path, save_dir, process_num, org_raster=None, b_save_patch_label=b_save_patch_label)


def main(options, args):

    img_path = args[0]
    io_function.is_file_exist(img_path)
    save_dir = options.save_dir
    process_num = options.process_num
    org_elevation_diff = options.elevation_diff

    segment_a_grey_image(img_path,save_dir,process_num,org_raster=org_elevation_diff)


if __name__ == "__main__":
    usage = "usage: %prog [options] image_path "
    parser = OptionParser(usage=usage, version="1.0 2021-3-1")
    parser.description = 'Introduction: segment a grey image using watershed algorithm   '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=1,
                      help="number of processes to create the mosaic")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir", default='./',
                      help="the folder to save results")

    parser.add_option("-c", "--elevation_diff",
                      action="store", dest="elevation_diff",
                      help="the original elevation difference")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
