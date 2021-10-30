#!/usr/bin/env python
# Filename: seg_polygonize_cal_attributes 
"""
introduction: convert segmentation results (labels) to polygons (polygonize), also do some calculation.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 31 March, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.basic as basic
import basic_src.io_function as io_function
import raster_statistic

import multiprocessing
from multiprocessing import Pool

def polygonize_one_label(idx,label_path,org_raster, stats, prefix,b_remove_nodata,process_num=1):

    save_dir = os.path.dirname(label_path)
    out_pre = io_function.get_name_no_ext(label_path)
    label_shp_path = os.path.join(save_dir, out_pre + '.shp')
    if os.path.isfile(label_shp_path):
        print('%s exist, skip'%label_shp_path)
        return idx,label_shp_path

    if b_remove_nodata is True:
        # remove nodato (it was copy from the input image)
        command_str = 'gdal_edit.py -unsetnodata ' + label_path
        res = os.system(command_str)
        if res != 0:
            return None,None

    # convert the label to shapefile
    command_string = 'gdal_polygonize.py -8 %s -b 1 -f "ESRI Shapefile" %s' % (label_path, label_shp_path)
    res = os.system(command_string)
    if res != 0:
        return None,None

    if org_raster is not None and stats is not None and prefix is not None:
        # get dem elevation information for each polygon,
        raster_statistic.zonal_stats_multiRasters(label_shp_path, org_raster, stats=stats,
                                                  prefix=prefix, process_num=process_num)

    return idx, label_shp_path

def polygonize_label_images(label_paths, org_raster=None, stats=None,prefix=None, process_num = 1, b_remove_nodata=False):
    '''
    convert labels to shape file and calculate attributes
    :param label_paths: a label path or multiple labels
    :param org_raster: orignal raster for calculating the attributes
    :param stats: e.g., ['mean', 'std','count']
    :param prefix: e.g., 'demD'
    :return:
    '''
    if isinstance(label_paths, list):
        label_path_list = label_paths
    else:
        label_path_list = [label_paths] # if only one label


    label_shp_list = []
    if process_num == 1:
        for idx,label in enumerate(label_path_list):
            _, out_shp = polygonize_one_label(label, org_raster,stats,prefix,b_remove_nodata,1)
            if out_shp is None:
                raise ValueError('failed in polygonize %s'%label)
            label_shp_list.append(out_shp)

    if process_num > 1:
        theadPool = Pool(process_num)
        parameters_list = [(idx,label, org_raster,stats,prefix,b_remove_nodata,1) for idx,label in enumerate(label_path_list)]
        results = theadPool.starmap(polygonize_one_label, parameters_list)
        for idx, res in results:
            if res is not None:
                label_shp_list.append(res)
            else:
                raise ValueError('failed in polygonize %s'%label_path_list[idx])
        theadPool.close()
    else:
        raise ValueError('Wrong process_num: %s'%str(process_num))

    return label_shp_list


def main(options, args):

    process_num = options.process_num
    b_remove_nodata = options.b_remove_nodata
    org_raster = options.org_raster

    label_path_list = [item for item in args if io_function.is_file_exist(item)]

    polygonize_label_images(label_path_list, org_raster=org_raster, stats=['mean', 'std','count'], prefix='demD',
                            process_num=process_num, b_remove_nodata=b_remove_nodata)


def test_polygonize_label_images():
    data_dir = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff')
    # img_path = os.path.join(data_dir,'WR_extent_grid_ids_DEM_diff_grid9274_8bit.tif')
    org_raster = os.path.join(data_dir,'WR_extent_grid_ids_DEM_diff_grid9274.tif')
    # save_dir = os.path.join(data_dir,'segment_parallel_9274')
    process_num = 8
    label_path_list = io_function.get_file_list_by_pattern('./','*grid9274*/*.tif')

    polygonize_label_images(label_path_list, org_raster=org_raster, stats=['mean', 'std', 'count'], prefix='demD',
                            process_num=process_num, b_remove_nodata=True)


if __name__ == '__main__':
    usage = "usage: %prog [options] label_path1 label_path2 ...  "
    parser = OptionParser(usage=usage, version="1.0 2021-3-31")
    parser.description = 'Introduction: convert segmentation results (labels) to polygons (polygonize), also do some calculation.  '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=1,
                      help="number of processes to create the mosaic")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir", default='./',
                      help="the folder to save results")

    parser.add_option("-o", "--org_raster",
                      action="store", dest="org_raster",
                      help="the original raster for calculating the attributes")

    parser.add_option("-r", "--b_remove_nodata",
                      action="store_true", dest="b_remove_nodata", default=False,
                      help="if set, remove nodata in the input label images")


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)


