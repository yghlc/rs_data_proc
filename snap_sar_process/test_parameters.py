#!/usr/bin/env python
# Filename: test_parameters.py 
"""
introduction: testing different parameters: resolution, cohWinAz, cohWinRg

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 02 June, 2023
"""

import os, sys

from sklearn.model_selection import ParameterGrid
import time
from datetime  import datetime

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import slurm_utility
import basic_src.basic as basic

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


def main():
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

    sar_type = 'sentinel-1'
    sar_list_txt = 'sar_list.txt'
    ext_shp = os.path.expanduser('~/Data/Arctic/pan_Arctic/extent/SAR_coh_test_region/SAR_coh_test.shp')
    dem_path = os.path.expanduser('~/Data/CopernicusDEM/SAR_coh_test_CoDEM.tif')

    save_dir_list = []
    work_dir_list = []

    for idx, config in enumerate(para_com_list):
        save_dir, work_dir = calculation(idx, config,sar_type, sar_list_txt, ext_shp, dem_path)
        save_dir_list.append(save_dir)
        work_dir_list.append(work_dir)
        time.sleep(1)       # need to wait 1 second on compute canada

    while True:
        submit_job_count = slurm_utility.get_submit_job_count(user_name)
        if submit_job_count > 0:
            print(machine_name, datetime.now(),'Waiting jobs to be finished, submitted job count: %d'%submit_job_count)
            time.sleep(300)
        else:
            break

    # organize results

    # # read miou and remove some files
    # over_miou_list = []
    # for work_dir, para_file in zip(work_dir_list, para_file_list):
    #     overall_miou = get_overall_miou_after_training(work_dir, para_file)
    #     over_miou_list.append(overall_miou)
    #     # remove_files(work_dir)
    #
    #     print('overall miou',os.path.basename(work_dir), overall_miou)


if __name__ == '__main__':
    main()
