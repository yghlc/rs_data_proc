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
import time

from datetime import datetime

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
import parameters

from dem_common import get_grid_id_from_path

# ArcticDEM results on ygAlpha
ArcticDEM_results_dir = os.path.expanduser('~/data1/ArcticDEM_results')
# tmp_dir = os.path.expanduser('~/tmpData')
tmp_dir = os.path.expanduser('~/data3/tmp_results')
bash_ini_dir = os.path.expanduser('~/Data/dem_diff_segment/sam_seg_script_ini_files')


def copy_modify_script_inifile(ini_dir, work_dir,dem_diff_color_dir, dem_diff_prompt_dir):

    files_copied = ['exe_sam_get_prompt.sh','exe_sam_seg_postproc.sh', 'sam_model_exp3.ini',
                    'main_para_exp3.ini','area_huang2023_demdiff.ini']
    file_paths = [ os.path.join(ini_dir,item) for item in files_copied]
    for item in file_paths:
        io_function.copyfiletodir(item,work_dir,overwrite=True)

    # rename 'area_huang2023_demdiff.ini'
    area_str = os.path.basename(work_dir)[:5]  # get ext?? , such ext18
    new_area_ini_file = 'area_huang2023_demdiff.ini'.replace('huang2023',area_str)
    io_function.move_file_to_dst('area_huang2023_demdiff.ini',new_area_ini_file)

    # modifiy parameters
    parameters.write_Parameters_file(new_area_ini_file,'area_name',area_str)
    parameters.write_Parameters_file(new_area_ini_file,'inf_image_dir',dem_diff_color_dir)
    parameters.write_Parameters_file(new_area_ini_file,'dem_diff_prompt_dir',dem_diff_prompt_dir)

    parameters.write_Parameters_file('main_para_exp3.ini','inference_regions',new_area_ini_file)

    return new_area_ini_file, 'main_para_exp3.ini'


def create_colorRelief_DEM_diff(bash_ini_dir,dem_diff_dir,save_dir):

    color_txt = os.path.join(bash_ini_dir,'dem_diff_color_5to5m.txt')
    org_dir = os.getcwd()
    os.chdir(save_dir)
    basic.outputlogMessage(f'change current directory to {save_dir}')
    io_function.copyfiletodir(color_txt,save_dir,overwrite=True)

    py_script = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/DEMscripts/dem_diff_to_colorRelief.py')
    cmd_str = f'{py_script}  {dem_diff_dir}'
    basic.os_system_exit_code(cmd_str)

    # change directory back
    os.chdir(org_dir)
    basic.outputlogMessage(f'change current directory to {org_dir}')


def set_each_grid_as_a_region(area_ini, main_para_ini,dem_diff_color_dir,area_ini_dir = 'area_grid'):

    dem_diff_file_dir = parameters.get_directory(area_ini, 'dem_diff_prompt_dir')
    dem_diff_file_or_pattern = parameters.get_string_parameters(area_ini, 'dem_diff_prompt_or_pattern')
    dem_diff_file_list = io_function.get_file_list_by_pattern(dem_diff_file_dir, dem_diff_file_or_pattern)
    if len(dem_diff_file_list) < 1:
        raise IOError('No DEM Diff file found by \n dem_diff_file_dir: %s \n dem_diff_file_or_pattern: %s' % (
        dem_diff_file_dir, dem_diff_file_or_pattern))
    else:
        basic.outputlogMessage('find %d DEM diff files' % len(dem_diff_file_list))

    if os.path.isdir(area_ini_dir) is False:
        io_function.mkdir(area_ini_dir)

    area_grid_ini_list = []

    for idx, dem_diff_file in enumerate(dem_diff_file_list):
        grid_id = get_grid_id_from_path(dem_diff_file)
        grid_str = f'grid{grid_id}'

        area_grid_ini = os.path.join(area_ini_dir,f'area_{grid_str}_demdiff.ini')
        io_function.copy_file_to_dst(area_ini, area_grid_ini, overwrite=True)

        # dir_name = os.path.dirname(dem_diff_file)
        fle_name = os.path.basename(dem_diff_file)

        tif_color = io_function.get_name_by_adding_tail(fle_name, 'color')

        parameters.write_Parameters_file(area_grid_ini,'area_name',grid_str)
        # change dem_diff_color_dir? not necessary
        # change dem_diff_prompt_dir? no need
        parameters.write_Parameters_file(area_grid_ini,'inf_image_or_pattern',tif_color)
        parameters.write_Parameters_file(area_grid_ini,'dem_diff_prompt_or_pattern',fle_name)

        # set prompt path
        prompt_txt_list = io_function.get_file_list_by_pattern('prompts',f'*{grid_str}*.txt')
        if len(prompt_txt_list) == 1:
            parameters.write_Parameters_file(area_grid_ini,'prompt_path',prompt_txt_list[0])
        else:
            raise ValueError(f'the number of prompt txt is not 1: {str(prompt_txt_list)}')

        area_grid_ini_list.append(area_grid_ini)


    save_list_txt = 'area_grid_ini_list.txt'
    io_function.save_list_to_txt(save_list_txt,area_grid_ini_list)

    parameters.write_Parameters_file(main_para_ini, 'inference_regions', save_list_txt)
    return area_grid_ini_list


def copy_organize_seg_results(area_grid_ini_list,sam_seg_result,save_dir):

    # for region don't have results, only copy it's ini file
    for grid_area_ini in area_grid_ini_list:
        # copy ini file
        io_function.copyfiletodir(grid_area_ini,save_dir,overwrite=True)

        grid_str = parameters.read_Parameters_file(grid_area_ini,'area_name')
        res_folder_list = io_function.get_file_list_by_pattern(sam_seg_result, f'{grid_str}*')
        if len(res_folder_list) == 1:
            # copy shp file (post)
            res_shp_dir = res_folder_list[0]
            shp_list = io_function.get_file_list_by_pattern(res_shp_dir,"*post*.shp")
            if len(shp_list) != 1:
                raise IOError('the number of post shapefile is not one')
            save_shp_file = os.path.join(save_dir, f'DEM_diff_{grid_str}_sam_seg.shp')
            io_function.copy_shape_file(shp_list[0],save_shp_file)

        elif len(res_folder_list) == 0:
            basic.outputlogMessage(f'warning, No SAM segment results for {grid_area_ini}')
        else:
            raise ValueError(f'the number of res_folder_list is not 1, grid_str: {grid_str}')



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


    io_function.is_folder_exist(dem_diff_dir)
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    if os.path.isdir(tmp_output_dir) is False:
        io_function.mkdir(tmp_output_dir)

    dem_diff_color_dir = os.path.join(tmp_output_dir,'demDiff_color')
    if os.path.isdir(dem_diff_color_dir) is False:
        io_function.mkdir(dem_diff_color_dir)

    # copy the script and ini file
    area_ini, main_para_ini = copy_modify_script_inifile(bash_ini_dir,work_dir,dem_diff_color_dir,dem_diff_dir)

    # create colorRelif DEM diff
    create_colorRelief_DEM_diff(bash_ini_dir, dem_diff_dir,dem_diff_color_dir)


    # run the script for segment
    cmd_str = './exe_sam_get_prompt.sh'
    # basic.os_system_exit_code(cmd_str)
    res = os.system(cmd_str)
    # basic.os_system_exit_code(command_str)
    if res != 0:
        basic.outputlogMessage('Warning, a step (could meering all prompts into a file) may be failed, please check')
        time.sleep(10)

    # set each grid as a region for segmentation after "get_prompt.sh", allow parallel computing when gettting prompt
    area_grid_ini_list = set_each_grid_as_a_region(area_ini, main_para_ini,dem_diff_color_dir)

    # run the script for segment
    cmd_str = './exe_sam_seg_postproc.sh'
    basic.os_system_exit_code(cmd_str)

    # copy files for each grid (results are in result_backup)
    copy_organize_seg_results(area_grid_ini_list, 'result_backup',  save_dir)


    # clean: (1) remove DEM diff colorRelif

    # output a done indicator
    io_function.save_dict_to_txt_json(done_indicator,{'done_time':str(datetime.now())})

    # change director back
    os.chdir(org_dir)


def test_sam_segment_a_big_region():
    # work_folder =  'ext09_for_ArcticDEM_proc'
    work_folder =  'ext07_for_ArcticDEM_proc'
    dem_diff_dir = os.path.join(ArcticDEM_results_dir,work_folder,'grid_dem_diffs')
    save_dir = os.path.join(ArcticDEM_results_dir,work_folder,'grid_dem_diffs_sam_results')
    tmp_save_dir = os.path.join(tmp_dir,work_folder)

    work_dir = os.path.abspath(work_folder)
    sam_segment_a_big_region(work_dir,dem_diff_dir,save_dir,tmp_save_dir)


def main():
    # test_sam_segment_a_big_region()

    org_dir = os.getcwd()
    basic.outputlogMessage(f'current directory to {org_dir}')

    # select_ext_list = ['ext00','ext01','ext02','ext03','ext04','ext05','ext06','ext07','ext08','ext09']
    select_ext_list = ['ext00','ext01','ext02','ext03','ext04']

    dem_diff_dir_list = io_function.get_file_list_by_pattern(ArcticDEM_results_dir,'ext??_*/*diffs*')
    for dem_diff_dir in dem_diff_dir_list:
        os.chdir(org_dir)   # change to original dir, because inside sam_segment_a_big_region, they change to other folder
        basic.outputlogMessage(f'processing {dem_diff_dir}')
        work_folder = os.path.basename(os.path.dirname(dem_diff_dir))
        # select region
        work_folder_begin_str = work_folder[:5]
        if work_folder_begin_str not in select_ext_list:
            basic.outputlogMessage(f'{work_folder} is not in the selected list, skip')
            continue

        save_dir = os.path.join(ArcticDEM_results_dir, work_folder, 'grid_dem_diffs_sam_results')
        tmp_save_dir = os.path.join(tmp_dir, work_folder)

        work_dir = os.path.abspath(work_folder)
        sam_segment_a_big_region(work_dir, dem_diff_dir, save_dir, tmp_save_dir)



    pass

if __name__ == '__main__':
    main()