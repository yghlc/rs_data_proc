#!/usr/bin/env python
# Filename: test_parameters.py 
"""
introduction: testing different parameters: resolution, cohWinAz, cohWinRg

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 02 June, 2023
"""

import os, sys
from optparse import OptionParser

from sklearn.model_selection import ParameterGrid
import time
from datetime  import datetime
import re

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import slurm_utility
import basic_src.basic as basic
import basic_src.io_function as io_function
import vector_gpd

machine_name = os.uname()[1]
user_name = 'lingcao'           # default user name on computeCanada

s1_py = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/snap_sar_process/process_multiple_SAR.py')
asar_py = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/snap_sar_process/process_multiple_ERS_Envisat_esa.py')

def get_para_list_from_grid_serach(para_config):
    para_combinations = list(ParameterGrid(para_config))
    return para_combinations

def calculation(idx, config, sar_type, sar_list_txt, ext_shp, dem_path):

    # Hyperparameters
    resolution, cohWinAz, cohWinRg = config["resolution"],config["cohWinAz"], config["cohWinRg"]
    para_str = '%s_%s_%s'%(str(int(resolution)).zfill(2), str(int(cohWinAz)).zfill(2), str(int(cohWinRg)).zfill(2))
    max_job = 100
    script_dir = os.path.expanduser('~/Data/sar_coherence_mapping/scripts_sar_coh')
    save_dir = './%s_coh_results_%s'%(str(idx).zfill(3), para_str)
    work_dir = './%s_work_dir_%s'%(str(idx).zfill(3), para_str)

    if sar_type.lower() == 'sentinel-1':
        cmd_str = s1_py + ' ' + sar_list_txt + ' -a %s '%ext_shp + ' -d %s '%save_dir  + ' -t /tmp ' + ' -e %s '%dem_path \
                  + ' -j %d '%max_job + ' --script_dir=%s '%script_dir + ' --working_dir=%s '%work_dir + ' --thread_num=4 ' \
                  + ' -r %s '%str(resolution) + ' --cohWinAz=%s '%str(cohWinAz) + ' --cohWinRg=%s '%str(cohWinRg)
    elif sar_type.lower() in ['ers', 'envisat'] :
        cmd_str =  asar_py + ' ' + sar_list_txt + ' -a %s '%ext_shp + '--sar_type=%s'%sar_type + ' -d %s'%save_dir  \
                   + ' -t /tmp ' + ' -e %s '%dem_path \
                   + ' -j %d '%max_job + ' --script_dir=%s '%script_dir + ' --working_dir=%s '%work_dir + ' --thread_num=4 ' \
                   + ' -r %s '%str(resolution) + ' --cohWinAz=%s '%str(cohWinAz) + ' --cohWinRg=%s '%str(cohWinRg)
    else:
        raise ValueError('unknown sar_type %s' % str(sar_type))

    basic.os_system_exit_code(cmd_str)

    return save_dir, work_dir

def subset_image(input,output,bound_box_str):
    if os.path.isfile(output):
        print('%s exists, skip cropping'%output)
        return
    cmd_str = 'gdalwarp -of GTiff -co tiled=yes -co bigtiff=if_safer -te %s ' % bound_box_str + input + ' ' + output
    basic.os_system_exit_code(cmd_str)

def organize_subset_coh(save_dir_list, ext_shp, out_dir=None):
    if out_dir is None:
        out_dir = 'all_subset'
    if os.path.isdir(out_dir) is False:
        io_function.mkdir(out_dir)

    bounding_boxes = vector_gpd.get_vector_file_bounding_box(ext_shp)  # minx, miny, maxx, maxy
    box_str = ' '.join([str(item) for item in bounding_boxes])

    for idx, res_dir in enumerate(save_dir_list):
        print('%d/%d: copy and subset result imagery' % (idx+1, len(save_dir_list)))
        tif_list = io_function.get_file_list_by_pattern(res_dir, '*/*.tif')
        para_strs = re.findall('[0-9]+', os.path.basename(res_dir))
        para_strs = '_'.join(para_strs)
        # subset images
        for tif in tif_list:
            output = os.path.join(out_dir, para_strs+'_'+os.path.basename(tif))
            subset_image(tif, output,box_str)


def test_organize_subset_coh():
    w_dir = os.path.expanduser('~/Data/sar_coherence_mapping/test1/test_parameters')
    save_dir_list = ['000_coh_results_10_03_05','001_coh_results_15_03_05']
    ext_shp = os.path.expanduser('~/Data/Arctic/pan_Arctic/extent/SAR_coh_test_region/SAR_coh_test.shp')
    save_dir_list = [os.path.join(w_dir,item) for item in save_dir_list]
    print(save_dir_list)
    organize_subset_coh(save_dir_list, ext_shp)

def main(options, args):
    # get combination of parameters
    para_config = {
            "resolution": [10, 15, 20, 30, 40 ],
            "cohWinAz": [3, 5, 7, 9, 11], #
            "cohWinRg": [5, 10, 15, 20, 25] #
    }
    para_com_list = get_para_list_from_grid_serach(para_config)
    para_com_list = [item for item in para_com_list if item['cohWinAz'] <= item['cohWinRg']]   # remove cohWinAz > cohWinRg
    # print(para_list)
    for idx, para in enumerate(para_com_list):
        # if para['cohWinAz'] > para['cohWinRg']:
        #     continue
        print(idx, para)

    sar_list_txt = args[0]      # 'sar_list.txt'
    sar_type = options.sar_type #'sentinel-1'
    ext_shp = options.aoi_shp   # os.path.expanduser('~/Data/Arctic/pan_Arctic/extent/SAR_coh_test_region/SAR_coh_test.shp')
    dem_path = options.elevation_file  # os.path.expanduser('~/Data/CopernicusDEM/SAR_coh_test_CoDEM.tif')

    save_dir_list = []
    work_dir_list = []

    for idx, config in enumerate(para_com_list):
        save_dir, work_dir = calculation(idx, config,sar_type, sar_list_txt, ext_shp, dem_path)
        save_dir_list.append(save_dir)
        work_dir_list.append(work_dir)
        time.sleep(2)       # need to wait 2 second on compute canada

    while True:
        submit_job_count = slurm_utility.get_submit_job_count(user_name,job_name_substr='coh')
        if submit_job_count > 0:
            print(machine_name, datetime.now(),'Waiting jobs to be finished, submitted job count: %d'%submit_job_count)
            time.sleep(300)
        else:
            break

    # organize results
    organize_subset_coh(save_dir_list, ext_shp, out_dir='all_results_crop')


if __name__ == '__main__':

    usage = "usage: %prog [options] sar_files.txt "
    parser = OptionParser(usage=usage, version="1.0 2023-6-4")
    parser.description = 'Introduction: testing the parameters for calculating SAR coherence, using SNAP  '

    # test_organize_subset_coh()
    # sys.exit(0)

    parser.add_option("-a", "--aoi_shp",
                      action="store", dest="aoi_shp",
                      help="a shapefile containing AOI")

    parser.add_option("-e", "--elevation_file",
                      action="store", dest="elevation_file",
                      help="DEM file used for terrain correction, if not set, will use SRTM 1 sec ")

    parser.add_option("", "--sar_type",
                      action="store", dest="sar_type", default='Sentinel-1',
                      help="the type of SAR data: Sentinel-1, ERS, Envisat")


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
