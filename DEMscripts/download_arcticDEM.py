#!/usr/bin/env python
# Filename: download_arcticDEM 
"""
introduction: download arcticDEM (strip or mosaic) for a specific region

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 25 December, 2020
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import basic_src.basic as basic

from urllib.parse import urlparse

import time

from dem_common import tarball_dir,arcticDEM_reg_tif_dir,arcticDEM_tile_tarball_dir,arcticDEM_tile_reg_tif_dir

from dem_common import grid_no_dem_txt, process_log_dir
from ArcticDEM_unpack_registration import process_dem_tarball, is_ArcticDEM_tiles

from datetime import datetime
from multiprocessing import Process

machine_name = os.uname()[1]
# the maximum number of processes for downloading in parallel
max_task_count = 2
download_tasks = []
b_unpack_after_downloading = True

def get_total_size(url_list):
    total_size = 0
    for url in url_list:
        size = io_function.get_url_file_size(url)       # bytes
        if size is not False:
            total_size += size
    return total_size/(1024.0*1024.0*1024.0)    # GB


def save_id_grid_no_dem(grid_id):
    # grid_dem_diff_less2dem_txt
    if os.path.isdir(process_log_dir) is False:
        io_function.mkdir(process_log_dir)
    # update grid_dem_diff_less2dem_txt file
    id_list = []
    if os.path.isfile(grid_no_dem_txt):
        id_list = io_function.read_list_from_txt(grid_no_dem_txt)    # no need covert to int
    id_str = str(grid_id)
    if id_str in id_list:
        return True
    else:
        # save
        id_list.append(str(grid_id))
        io_function.save_list_to_txt(grid_no_dem_txt,id_list)
        basic.outputlogMessage('Save gird id (%d) to %s' % (grid_id,grid_no_dem_txt))
        return True

def wget_file_url(url,save_path):
    cmd_str = 'wget --no-check-certificate --output-document=%s  %s' % (save_path,url)
    status, result = basic.exec_command_string(cmd_str)
    return status, result

def run_a_process_download(url, tar_path, save_tif_dir, process_num=1, b_unpack=False):
    status, result = wget_file_url(url,tar_path)
    # try new url if it exists
    if status != 0 and "302 Moved Temporarily" in result:
        str_list = result.split()
        loc_idx = str_list.index('Location:')
        new_url = str_list[loc_idx + 1]  # find the new URL
        print('try to download the file using new url: %s' % new_url)
        status, result = wget_file_url(new_url,tar_path)

    if status != 0:
        if process_num == 1:
            print(result)
            basic.outputlogMessage('failed to download DEM: %s' % url)
        else:
            # when run in paralle, we may not able to save the log to disk
            print('\n failed info\n',result)
            print('failed to download DEM: %s' % url)
            print('\n\n')
        sys.exit(status)

    # unpack after downloading
    if b_unpack:
        tar_list = [tar_path]
        work_dir = './'
        b_rm_inter = True
        b_rm_tarball = True
        if is_ArcticDEM_tiles(tar_list):
            apply_registration = False
        else:
            apply_registration = True
        process_dem_tarball(tar_list, work_dir, save_tif_dir, remove_inter_data=b_rm_inter, rm_tarball=b_rm_tarball, apply_registration=apply_registration)


def download_dem_tarball(dem_index_shp, extent_polys, save_folder, pre_name, reg_tif_dir=None, poly_ids=None,b_arcticDEM_tile=False):
    # read dem polygons and url
    dem_polygons, dem_urls = vector_gpd.read_polygons_attributes_list(dem_index_shp, 'fileurl',b_fix_invalid_polygon=False)

    basic.outputlogMessage('%d dem polygons in %s' % (len(dem_polygons), dem_index_shp))

    dem_tar_ball_list = []
    reg_tifs_list = []
    curr_dir = os.getcwd()
    b_save_grid_id_noDEM = True
    if poly_ids is None:
        poly_ids = [idx for idx in range(len(extent_polys)) ]
        b_save_grid_id_noDEM = False    # if poly_ids is not the global unique id, then don't save it.

    if os.path.isfile('no_registration_strips.txt'):
        no_registration_strips = io_function.read_list_from_txt('no_registration_strips.txt')
    else:
        no_registration_strips = []

    # tarballs is being downloaded
    downloading_tarballs = []

    for count, (idx, ext_poly) in enumerate(zip(poly_ids, extent_polys)):
        basic.outputlogMessage('get data for the %d th extent (%d/%d)' % (idx, count, len(extent_polys)))

        save_txt_path = pre_name + '_dem_urls_poly_%d.txt' % idx
        if os.path.isfile(save_txt_path):
            urls = io_function.read_list_from_txt(save_txt_path)
            basic.outputlogMessage('read %d dem urls from %s' % (len(urls),save_txt_path))
        else:
            # get fileurl
            dem_poly_ids = vector_gpd.get_poly_index_within_extent(dem_polygons,ext_poly)
            basic.outputlogMessage('find %d DEM within %d th extent' % (len(dem_poly_ids), (idx )))
            urls = [dem_urls[id] for id in dem_poly_ids]

            # save to txt
            io_function.save_list_to_txt(save_txt_path, urls)
            basic.outputlogMessage('save dem urls to %s' % save_txt_path)

        if len(urls) > 0:

            # total_size_GB = get_total_size(urls)  # internet access, parallel running may cause problem. The info is not important
            # basic.outputlogMessage('the size of files will be downloaded is %.4lf GB for the %d th extent '%(total_size_GB,(idx+1)))
            # time.sleep(5)   # wait 5 seconds

            # download them using wget one by one
            for ii, url in enumerate(urls):
                tmp = urlparse(url)

                # in the Strip DEM, there are around 700 url are point to tif files, failed to download them
                # e.g. /mnt/pgc/data/elev/dem/setsm/ArcticDEM/geocell/v3.0/2m_temp/n59w137/SETSM_WV03_20150518_104001000B703200_104001000C715B00_seg8_2m_v3.0_dem.tif
                if url.startswith('/mnt') and url.endswith('.tif'):
                    basic.outputlogMessage("error: not a valid url: %s"%url)
                    continue

                filename = os.path.basename(tmp.path)
                save_dem_path = os.path.join(save_folder,filename)
                if reg_tif_dir is not None:
                    tar_base = os.path.basename(filename)[:-7]
                    # file_pattern = ['*dem_reg.tif', '*reg_dem.tif'] # Arctic strip and tile (mosaic) version
                    if b_arcticDEM_tile:
                        reg_tifs = io_function.get_file_list_by_pattern(reg_tif_dir, tar_base + '*reg_dem.tif')
                    else:
                        reg_tifs = io_function.get_file_list_by_pattern(reg_tif_dir,tar_base+'*dem_reg.tif')
                    if len(reg_tifs) > 0:
                        basic.outputlogMessage('warning, unpack and registrated tif for %s already exists, skip downloading' % filename)
                        reg_tifs_list.append(reg_tifs[0])
                        continue

                    if './'+tar_base in no_registration_strips:
                        basic.outputlogMessage(
                            'warning, %s is in no_registration_strips list, skip downloading' % filename)
                        continue

                if filename in downloading_tarballs:
                    basic.outputlogMessage('warning, %s is being downloaded by other processes'%filename)
                    continue

                if os.path.isfile(save_dem_path):
                    basic.outputlogMessage('warning, %s already exists, skip downloading'%filename)
                else:
                    # download the dem
                    basic.outputlogMessage('starting downloading %d th DEM (%d in total)'%((ii+1),len(urls)))
                    downloading_tarballs.append(filename)

                    # os.chdir(save_folder)

                    # run_a_process_download(url)  # download

                    ##################################################
                    # download in parallel
                    basic.check_exitcode_of_process(download_tasks)  # if there is one former job failed, then quit

                    while True:
                        job_count = basic.alive_process_count(download_tasks)
                        if job_count >= max_task_count:
                            print(machine_name, datetime.now(),'You are running %d or more tasks in parallel, wait ' % max_task_count)
                            time.sleep(60)  #
                            continue
                        break

                    # start the processing
                    sub_process = Process(target=run_a_process_download, args=(url, save_dem_path, reg_tif_dir, max_task_count,b_unpack_after_downloading))  # start a process, don't wait
                    sub_process.start()
                    download_tasks.append(sub_process)

                    basic.close_remove_completed_process(download_tasks)


                    # os.chdir(curr_dir)

                dem_tar_ball_list.append(save_dem_path)

        else:
            basic.outputlogMessage('Warning, can not find DEMs within %d th extent'%(idx))
            if b_save_grid_id_noDEM:
                save_id_grid_no_dem(idx)

    # wait until all task complete
    while True:
        job_count = basic.alive_process_count(download_tasks)
        if job_count > 0:
            print(machine_name, datetime.now(), 'wait until all task are completed, alive task account: %d ' % job_count)
            time.sleep(60)  #
        else:
            break


    return dem_tar_ball_list, reg_tifs_list

def main(options, args):

    extent_shp = args[0]
    dem_index_shp = args[1]
    b_arcticDEM_tile = False

    global max_task_count
    max_task_count = options.max_process_num

    if 'Tile' in os.path.basename(dem_index_shp):
        save_folder =arcticDEM_tile_tarball_dir
        reg_tif_dir = arcticDEM_tile_reg_tif_dir
        b_arcticDEM_tile = True
    else:
        save_folder = tarball_dir
        reg_tif_dir = arcticDEM_reg_tif_dir

    # use the user specific save_dir for saving downloaded tarballs
    if options.save_dir is not None:
        save_folder = options.save_dir
    if os.path.isdir(save_folder) is False:
        io_function.mkdir(save_folder)
    save_folder = os.path.abspath(save_folder)  # change to absolute path

    pre_name = os.path.splitext(os.path.basename(extent_shp))[0]
    pre_name += '_Tile' if 'Tile' in os.path.basename(dem_index_shp) else '_Strip'

    # extent polygons and projection (proj4)
    extent_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(extent_shp)
    dem_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(dem_index_shp)

    if extent_shp_prj != dem_shp_prj:
        basic.outputlogMessage('%s and %s do not have the same projection, will reproject %s'
                               %(extent_shp,dem_index_shp,os.path.basename(extent_shp)))
        epsg = map_projection.get_raster_or_vector_srs_info_epsg(dem_index_shp)
        # print(epsg)
        # extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,dem_shp_prj.strip())
        extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,epsg)
    else:
        extent_polys = vector_gpd.read_polygons_gpd(extent_shp)

    # read 'grid_id' if the extent shp is from grid shp file, if not, grid_id_list will be None
    grid_id_list = vector_gpd.read_attribute_values_list(extent_shp,'grid_id')

    if len(extent_polys) < 1:
        raise ValueError('No polygons in %s'%extent_shp)
    else:
        basic.outputlogMessage('%d extent polygons in %s'%(len(extent_polys),extent_shp))

    download_dem_tarball(dem_index_shp,extent_polys,save_folder,pre_name,reg_tif_dir=reg_tif_dir,
                         poly_ids=grid_id_list,b_arcticDEM_tile=b_arcticDEM_tile)


if __name__ == "__main__":

    usage = "usage: %prog [options] extent_shp dem_indexes_shp"
    parser = OptionParser(usage=usage, version="1.0 2020-12-25")
    parser.description = 'Introduction: download ArcticDEM within an extent  '
    # parser.add_option("-x", "--save_xlsx_path",
    #                   action="store", dest="save_xlsx_path",
    #                   help="save the sence lists to xlsx file")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",
                      help="the folder to save DEMs")

    parser.add_option("-m", "--max_process_num",
                      action="store", dest="max_process_num", type=int,default=16,
                      help="the maximum number of processes for downloading in parallel")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
