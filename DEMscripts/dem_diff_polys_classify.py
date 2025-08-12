#!/usr/bin/env python
# Filename: dem_diff_polys_classify.py 
"""
introduction: using deep learning (e.g., CLIP) to classify DEM diff polygons into different classes
for the entire ArcticDEM domain

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 10 August, 2025
"""

import os,sys
from optparse import OptionParser
import time

from datetime import datetime

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
import parameters
import raster_io

from dem_common import get_grid_id_from_path

# ArcticDEM results on ygAlpha
ArcticDEM_results_dir = os.path.expanduser('~/data1/ArcticDEM_results')
tmp_dir = os.path.expanduser('~/data3/tmp_results')
bash_ini_dir = os.path.expanduser('~/Data/slump_demdiff_classify/demDiff_classify/template_script_ini_files/')

from dem_segment_sam import create_colorRelief_DEM_diff
from dem_common import grid_20_id_raster, find_neighbours_grids

import torch

def copy_modify_script_inifile(ini_dir, work_dir, dem_diff_color_dir):

    files_copied = ['run_classify.sh', 'model_clip.ini','main_para_exp3.ini',
                    'area_ext00_alaska_extend_simple_DEMdiff.ini']
    file_paths = [ os.path.join(ini_dir,item) for item in files_copied]
    for item in file_paths:
        io_function.copyfiletodir(item,work_dir,overwrite=True)

    # rename area ini
    area_str = os.path.basename(work_dir)[:5]  # get ext?? , such ext18
    new_area_ini_file = f'area_{area_str}_demdiff.ini'
    io_function.move_file_to_dst('area_ext00_alaska_extend_simple_DEMdiff.ini',new_area_ini_file)

    # modifiy parameters
    parameters.write_Parameters_file(new_area_ini_file,'area_name',area_str)
    parameters.write_Parameters_file(new_area_ini_file,'inf_image_dir',dem_diff_color_dir)

    parameters.write_Parameters_file('main_para_exp3.ini','inference_regions',new_area_ini_file)

    return new_area_ini_file, 'main_para_exp3.ini'

def find_grid_dem_diff_poly_shp(grid_id, dem_diff_poly_shps, shp_grid_id_list):
    for shp, shp_grid_id in zip(dem_diff_poly_shps, shp_grid_id_list):
        # this is a potential bug, for example: grid10, is in grid101, grid102, .....
        # if grid_str in basename:
        #     return shp
        if grid_id == shp_grid_id:
            return shp

    return None

def find_dem_diff_include_neighbour(grid_ids_2d, grid_id, dem_diff_file_list, grid_id_list):

    neighbours_grids = find_neighbours_grids(grid_ids_2d,grid_id, connect=8)
    # print(f'grid {grid_id} has neighbours: {neighbours_grids}')
    # sys.exit(0)
    neighbours_grids.append(grid_id)

    select_idx_list = [grid_id_list.index(item) for item in neighbours_grids if item in grid_id_list]
    # only keep the existing files
    select_dem_diff_list = [dem_diff_file_list[item] for item in select_idx_list if os.path.isfile(dem_diff_file_list[item])]

    return select_dem_diff_list

def find_dem_diff_color(sel_dem_file_list,dem_diff_color_dir, copy_dir):
    if os.path.isdir(copy_dir) is False:
        io_function.mkdir(copy_dir)

    tif_color_list = []
    for dem_diff_file in sel_dem_file_list:
        file_name = os.path.basename(dem_diff_file)
        tif_color_filename = io_function.get_name_by_adding_tail(file_name, 'color')
        tif_color = os.path.join(dem_diff_color_dir, tif_color_filename)
        if os.path.isfile(tif_color) is False:
            raise IOError(f'{tif_color} does not exist')

        tif_color_list.append(tif_color)
        target = os.path.join(copy_dir,tif_color_filename)
        if os.path.isfile(target):
            print(f'{target} exist, skip creating the link')
            continue
        cmd_str = 'ln -s %s %s' % (tif_color, target)
        basic.os_system_exit_code(cmd_str)

    return tif_color_list



def set_each_grid_as_a_region_for_classify(area_ini, main_para_ini, dem_diff_file_dir, dem_diff_color_dir, dem_diff_poly_dir,
                                           area_ini_dir = 'area_grid_ini_files'):

    dem_diff_file_or_pattern = parameters.get_string_parameters(area_ini, 'inf_image_or_pattern')
    dem_diff_file_list = io_function.get_file_list_by_pattern(dem_diff_file_dir, dem_diff_file_or_pattern)
    grid_id_list = [get_grid_id_from_path(item) for item in dem_diff_file_list]
    dem_diff_poly_shps = io_function.get_file_list_by_pattern(dem_diff_poly_dir, '*.shp')
    shp_grid_id_list = [get_grid_id_from_path(item) for item in dem_diff_poly_shps]

    if len(dem_diff_file_list) < 1:
        raise IOError('No DEM Diff file found by \n dem_diff_file_dir: %s \n dem_diff_file_or_pattern: %s' % (
        dem_diff_file_dir, dem_diff_file_or_pattern))
    else:
        basic.outputlogMessage('find %d DEM diff files' % len(dem_diff_file_list))

    if os.path.isdir(area_ini_dir) is False:
        io_function.mkdir(area_ini_dir)

    grid_ids_2d, nodata = raster_io.read_raster_one_band_np(grid_20_id_raster)

    area_grid_ini_list = []

    for grid_id, dem_diff_file in zip(grid_id_list, dem_diff_file_list):
        # grid_id = get_grid_id_from_path(dem_diff_file)
        grid_str = f'grid{grid_id}'

        # find polygons shp
        dem_diff_poly_shp = find_grid_dem_diff_poly_shp(grid_id,dem_diff_poly_shps, shp_grid_id_list)
        if dem_diff_poly_shp is None:
            basic.outputlogMessage(f'Warning, grid {grid_id} does not contain DEM diff polygons')
            continue

        area_grid_ini = os.path.join(area_ini_dir,f'area_{grid_str}_demdiff.ini')

        if os.path.isfile(area_grid_ini):
            print(f'{area_grid_ini} exists, skip, remove it to generate a new one')
            area_grid_ini_list.append(area_grid_ini)
            continue

        io_function.copy_file_to_dst(area_ini, area_grid_ini, overwrite=True)

        sel_dem_file_list = find_dem_diff_include_neighbour(grid_ids_2d,grid_id,dem_diff_file_list,grid_id_list)

        grid_dem_diff_color_dir = os.path.join( 'soft_link_to_raster',f'{grid_str}_DEM_diff_color_link')

        sel_dem_file_list = find_dem_diff_color(sel_dem_file_list,dem_diff_color_dir,grid_dem_diff_color_dir)

        parameters.write_Parameters_file(area_grid_ini,'area_name',grid_str)

        parameters.write_Parameters_file(area_grid_ini, 'inf_image_dir', grid_dem_diff_color_dir)

        parameters.write_Parameters_file(area_grid_ini, 'all_polygons_labels', dem_diff_poly_shp)


        area_grid_ini_list.append(area_grid_ini)


    save_list_txt = 'area_grid_ini_list.txt'
    io_function.save_list_to_txt(save_list_txt,area_grid_ini_list)

    parameters.write_Parameters_file(main_para_ini, 'inference_regions', save_list_txt)
    return area_grid_ini_list



    pass

def copy_organize_classify_results(area_grid_ini_list,result_folder,save_dir):

    # for region don't have results, only copy it's ini file
    for grid_area_ini in area_grid_ini_list:
        # copy ini file
        io_function.copyfiletodir(grid_area_ini,save_dir,overwrite=True)

        grid_str = parameters.read_Parameters_file(grid_area_ini,'area_name')
        area_folder = parameters.get_area_name_remark_time(grid_area_ini)

        res_shp = os.path.join(result_folder,area_folder,area_folder + '-predicted_classID.shp')
        if os.path.isfile(res_shp):
            output_shp = os.path.join(save_dir, os.path.basename(res_shp))
            io_function.copy_shape_file(res_shp, output_shp)
        else:
            raise IOError(f'result shp: {res_shp} does not exist, grid_str: {grid_str}')

def clip_classify_a_big_region(para_file, work_dir, dem_diff_dir, dem_diff_poly_dir, save_dir, tmp_output_dir):

    # create a working folder, then switch to it
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
    org_dir = os.getcwd()
    basic.outputlogMessage(f'current directory to {org_dir}')
    os.chdir(work_dir)
    basic.outputlogMessage(f'change current directory to {work_dir}')

    # check if it's done
    classify_done_indicator = f'{os.path.basename(work_dir)}.done'
    copy_done_indicator = f'{os.path.basename(work_dir)}.copied_done'
    b_classified = False
    # b_copied = False
    if os.path.isfile(classify_done_indicator):
        basic.outputlogMessage(f'this region: {work_dir} has been classified')
        b_classified = True
        # return
    if os.path.isfile(copy_done_indicator):
        basic.outputlogMessage(f'this region: {work_dir} has been processed and is completed')
        return

    io_function.is_folder_exist(dem_diff_dir)
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    if os.path.isdir(tmp_output_dir) is False:
        io_function.mkdir(tmp_output_dir)

    dem_diff_color_dir = os.path.join(tmp_output_dir, 'demDiff_color')
    if os.path.isdir(dem_diff_color_dir) is False:
        io_function.mkdir(dem_diff_color_dir)

    # copy the script and ini file
    area_ini, main_para_ini = copy_modify_script_inifile(bash_ini_dir, work_dir, dem_diff_color_dir)

    # create colorRelif DEM diff
    if b_classified is False:
        create_colorRelief_DEM_diff(bash_ini_dir, dem_diff_dir, dem_diff_color_dir)

    # set each grid as a region for segmentation after "get_prompt.sh", allow parallel computing when gettting prompt
    area_grid_ini_list = set_each_grid_as_a_region_for_classify(area_ini, main_para_ini, dem_diff_dir, dem_diff_color_dir, dem_diff_poly_dir)

    num_gpus = torch.cuda.device_count()
    if num_gpus > 1:
        basic.outputlogMessage(f"Multiple ({num_gpus}) GPU exists, please manually change CUDA_VISIBLE_DEVICES to utilize these GPUs ")
    else:
        # run prediction
        cmd_str = './run_classify.sh'
        if b_classified is False:
            basic.os_system_exit_code(cmd_str)

    # run post-processing? (in ./run_classify.sh, it will modify original shapefile)


    # clean: (1) remove DEM diff colorRelif

    # output a done indicator (put this into the run_classify.sh, because we need to manually run this if multiple GPU exists)
    # io_function.save_dict_to_txt_json(done_indicator, {'done_time': str(datetime.now())})


    # organize the prediction results into ArcticDEM results
    expr_name = parameters.get_string_parameters(para_file, 'expr_name')
    outdir = os.path.join(parameters.get_directory(para_file, 'inf_output_dir'), expr_name)
    copy_organize_classify_results(area_grid_ini_list,outdir,save_dir)

    io_function.save_dict_to_txt_json(copy_done_indicator, {'copy_done_time': str(datetime.now())})

    # change director back
    os.chdir(org_dir)


def main(options, args):
    # test_sam_segment_a_big_region()

    task_name = options.task_name
    para_file = options.para_file

    org_dir = os.getcwd()
    basic.outputlogMessage(f'current directory to {org_dir}')

    ext_to_proc = 'ext_to_proc_list.txt'
    if len(args) > 0:
        select_ext_list = [item for item in args]
    elif os.path.isfile(ext_to_proc):
        select_ext_list = io_function.read_list_from_txt(ext_to_proc)
    else:
        select_ext_list = ['ext09'] # for testing
        # select_ext_list = ['ext00','ext01','ext02','ext03','ext04','ext05','ext06','ext07','ext08','ext09','ext13']
        #select_ext_list = ['ext00','ext01','ext02','ext03','ext04']

    dem_diff_dir_list = io_function.get_file_list_by_pattern(ArcticDEM_results_dir,'ext??_*/grid_dem_diffs')
    for dem_diff_dir in dem_diff_dir_list:
        os.chdir(org_dir)   # change to original dir, because inside sam_segment_a_big_region, they change to other folder
        basic.outputlogMessage(f'processing {dem_diff_dir}')
        work_folder = os.path.basename(os.path.dirname(dem_diff_dir))
        dem_diff_poly_dir = os.path.join(os.path.dirname(dem_diff_dir), 'grid_dem_diffs_sam_results')
        # select region
        work_folder_begin_str = work_folder[:5]
        if work_folder_begin_str not in select_ext_list:
            basic.outputlogMessage(f'{work_folder} is not in the selected list, skip')
            continue

        save_dir = os.path.join(ArcticDEM_results_dir, work_folder, f'grid_dem_diffs_{task_name}_res')
        tmp_save_dir = os.path.join(tmp_dir, work_folder)

        work_dir = os.path.abspath(work_folder)
        clip_classify_a_big_region(para_file, work_dir, dem_diff_dir, dem_diff_poly_dir, save_dir, tmp_save_dir)

if __name__ == '__main__':
    usage = "usage: %prog [options] ext00 ext01 ext02 ... "
    parser = OptionParser(usage=usage, version="1.0 2025-08-10")
    parser.description = 'Introduction: run image classification for all grids in the ArcticDEM domain '


    parser.add_option("-d", "--arcticDEM_res_dir",
                      action="store", dest="arcticDEM_res_dir",
                      help="the folder that contains ArcticDEM results")

    parser.add_option("-n", "--task_name",
                      action="store", dest="task_name", default='reduction_poly_classify',
                      help="the name of the task, using for saving folder name")

    parser.add_option("-p", "--para_file",
                      action="store", dest="para_file", default='main_para_exp3.ini',
                      help="the main parameter file")


    (options, args) = parser.parse_args()
    # if len(sys.argv) < 2 or len(args) < 1:
    #     parser.print_help()
    #     sys.exit(2)

    main(options, args)