#!/usr/bin/env python
# Filename: dem_segment_sam.py 
"""
introduction: using Segment anything model to segment elevation reduction using

~/codes/PycharmProjects/BigImageMapper/sam_dir/sam_predict.py
Please following the files (*.sh, *.ini) in ~/codes/PycharmProjects/BigImageMapper/sam_dir to set input files

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 07 February, 2025
"""

import os,sys

from datetime import datetime

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
import parameters

# ArcticDEM results on ygAlpha
ArcticDEM_results_dir = os.path.expanduser('~/data1/ArcticDEM_results')
tmp_dir = os.path.expanduser('~/tmpData')
bash_ini_dir = os.path.expanduser('~/Data/dem_diff_segment/sam_seg_script_ini_files')


def copy_modify_script_inifile(ini_dir, work_dir,dem_diff_color_dir, dem_diff_prompt_dir):

    files_copied = ['exe_sam.sh', 'sam_model_exp3.ini','main_para_exp3.ini','area_huang2023_demdiff.ini']
    file_paths = [ os.path.join(ini_dir,item) for item in files_copied]
    for item in file_paths:
        io_function.copyfiletodir(item,work_dir,overwrite=False)

    # rename 'area_huang2023_demdiff.ini'
    area_str = os.path.basename(work_dir)[:5]  # get ext?? , such ext18
    new_area_ini_file = 'area_huang2023_demdiff.ini'.replace('huang2023',area_str)
    io_function.move_file_to_dst('area_huang2023_demdiff.ini',new_area_ini_file)

    # modifiy parameters
    parameters.write_Parameters_file(new_area_ini_file,'area_name',area_str)
    parameters.write_Parameters_file(new_area_ini_file,'inf_image_dir',dem_diff_color_dir)
    parameters.write_Parameters_file(new_area_ini_file,'dem_diff_prompt_dir',dem_diff_prompt_dir)

    parameters.write_Parameters_file('main_para_exp3.ini','inference_regions',new_area_ini_file)

    return True


def create_colorRelief_DEM_diff(bash_ini_dir,dem_diff_dir,save_dir):

    color_txt = os.path.join(bash_ini_dir,'dem_diff_color_5to5m.txt')
    org_dir = os.getcwd()
    os.chdir(save_dir)
    basic.outputlogMessage(f'change current directory to {save_dir}')
    io_function.copyfiletodir(color_txt,save_dir,overwrite=False)

    py_script = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/DEMscripts/dem_diff_to_colorRelief.py')
    cmd_str = f'{py_script}  dem_diff_dir'
    basic.os_system_exit_code(cmd_str)

    # change directory back
    os.chdir(org_dir)
    basic.outputlogMessage(f'change current directory to {org_dir}')


def sam_segment_a_big_region(work_dir, dem_diff_dir, save_dir, tmp_output_dir):
    # create a working folder, then switch to it
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
    org_dir = os.getcwd()
    basic.outputlogMessage(f'current directory to {org_dir}')
    os.chdir(work_dir)
    basic.outputlogMessage(f'change current directory to {work_dir}')

    # check if it's done
    done_indicator = f'{os.path.basename(work_dir)}.done'
    if os.path.isfile(done_indicator):
        basic.outputlogMessage(f'this region: {work_dir} has been segmented, skip')
        return

    dem_diff_color_dir = os.path.join(tmp_output_dir,'demDiff_color')
    if os.path.isdir(dem_diff_color_dir) is False:
        io_function.mkdir(dem_diff_color_dir)

    # copy the script and ini file
    copy_modify_script_inifile(bash_ini_dir,work_dir,dem_diff_color_dir,dem_diff_dir)

    # create colorRelif DEM diff
    create_colorRelief_DEM_diff(bash_ini_dir, dem_diff_dir,dem_diff_color_dir)



    # run the script for segment
    cmd_str = './exe_sam.sh'
    basic.os_system_exit_code(cmd_str)


    io_function.is_folder_exist(dem_diff_dir)
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    if os.path.isdir(tmp_output_dir):
        io_function.mkdir(tmp_output_dir)

    # clean: (1) remove DEM diff colorRelif


    # post-processing and copy files for each grid


    # output a done indicator
    io_function.save_dict_to_txt_json(done_indicator,{'done_time':str(datetime.now())})



def test_sam_segment_a_big_region():
    work_folder =  'ext09_for_ArcticDEM_proc'
    dem_diff_dir = os.path.join(ArcticDEM_results_dir,work_folder,'grid_dem_diffs')
    save_dir = os.path.join(ArcticDEM_results_dir,work_folder,'grid_dem_diffs_sam_results')
    tmp_save_dir = os.path.join(tmp_dir,work_folder)

    work_dir = os.path.abspath(work_folder)
    sam_segment_a_big_region(work_dir,dem_diff_dir,save_dir,tmp_save_dir)


def main():
    test_sam_segment_a_big_region()

    # dem_diff_dir_list = io_function.get_file_list_by_pattern(ArcticDEM_results_dir,'ext??_*/*diffs*')
    # print(dem_diff_dir_list)


    pass

if __name__ == '__main__':
    main()