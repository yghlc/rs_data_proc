#!/usr/bin/env python
# Filename: check_vector_file.py 
"""
introduction: check vector files, some files are incompleted due to unforseen issues.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 09 September, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function
import basic_src.basic as basic
import vector_gpd

from dem_common import grid_dem_headwall_shp_dir

machine_name = os.uname()[1]

from multiprocessing import Pool

def check_one_vector_file(idx,total,file_path, good_file_list):
    if os.path.basename(file_path) in good_file_list:
        return True
    try:
        print('checking %d/%d' % (idx, total))
        # src = raster_io.open_raster_read(tif_path)
        geometries = vector_gpd.read_polygons_gpd(file_path,b_fix_invalid_polygon=False)
        if len(geometries) < 1:
            basic.outputlogMessage('No geometries in %s'%file_path)
            return False
        geometries_list = [ item for item in geometries ]
        if None in geometries_list:
            basic.outputlogMessage('None geometries in %s'%file_path)
            return False
    except:
        basic.outputlogMessage('incomplete vector: %s' % file_path)
        return False

    return True


def main(options, args):
    # process_num = multiprocessing.cpu_count()
    process_num = options.process_num
    data_folder = grid_dem_headwall_shp_dir


    vector_files = io_function.get_file_list_by_ext('.shp',data_folder ,bsub_folder=True)
    save_invalid_txt_path = os.path.basename(data_folder) + '_incomplete_list.txt'
    save_good_txt_path = os.path.basename(data_folder) + '_good_list.txt'
    file_count = len(vector_files)
    basic.outputlogMessage('The number of vector files: %d'%file_count)

    good_files = []
    if os.path.isfile(save_good_txt_path):
        good_files.extend(io_function.read_list_from_txt(save_good_txt_path))
    incomplete_files = []

    # remove good one for the list
    if len(good_files)>0:
        vector_files = [item for item in vector_files if os.path.basename(item) not in good_files]

    if process_num == 1:
        # tifs = io_function.get_file_list_by_ext('.tif',arcticDEM_reg_tif_dir, bsub_folder=False)
        for idx, tif in enumerate(vector_files):
            if check_one_vector_file(idx, file_count, tif, good_files):
                good_files.append(os.path.basename(tif))
            else:
                incomplete_files.append(os.path.basename(tif))
    else:
        theadPool = Pool(process_num)  # multi processes
        parameters_list = [(idx, file_count, tif, good_files) for idx, tif in enumerate(vector_files)]
        results = theadPool.starmap(check_one_vector_file, parameters_list)  # need python3
        for tif, res in zip(vector_files, results):
            if res:
                good_files.append(os.path.basename(tif))
            else:
                incomplete_files.append(os.path.basename(tif))
        theadPool.close()

    io_function.save_list_to_txt(save_invalid_txt_path, incomplete_files)
    io_function.save_list_to_txt(save_good_txt_path, good_files)

    # tifs = io_function.get_file_list_by_ext('.tif', grid_dem_diff_dir, bsub_folder=False)
    # invalid_tif = check_tif(tifs)
    # for tif in invalid_tif:
    #     print('removing %s'%tif)
    # io_function.delete_file_or_dir(tif)


if __name__ == '__main__':
    if __name__ == '__main__':
        usage = "usage: %prog [options]  "
        parser = OptionParser(usage=usage, version="1.0 2021-9-9")
        parser.description = 'Introduction: find out incomplete vector files  '

        parser.add_option("-p", "--process_num",
                          action="store", dest="process_num", type=int, default=8,
                          help="number of processes to checking vector files")

        (options, args) = parser.parse_args()
        # print(options.create_mosaic)

        # if len(sys.argv) < 2 or len(args) < 1:
        #     parser.print_help()
        #     sys.exit(2)

        main(options, args)
        pass

