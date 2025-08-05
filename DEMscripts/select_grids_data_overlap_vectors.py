#!/usr/bin/env python
# Filename: select_grids_data_overlap_vectors.py 
"""
introduction: select grid data (DEM difference raster, polygons, and others) if they overlap input polygons or points

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 18 December, 2023
"""

import os, sys
from optparse import OptionParser
import time
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import pandas as pd
import geopandas as gpd
import vector_gpd

from datetime import datetime

arcticDEM_res_dir ="/Volumes/Seagate8T/ArcticDEM_results"
arcticdata_root="https://arcticdata.io/data/10.18739/A2ZS2KF4B/Lingcao-output"

def find_grids_overlap_vector_shp(grid_indexes_shp,vector_shp_list):

    if isinstance(vector_shp_list,list) is False:
        vector_shp_list = [vector_shp_list]

    overlap_touch_list = []
    for vector_shp in vector_shp_list:
        overlap_touch = vector_gpd.geometries_overlap_another_group(grid_indexes_shp,vector_shp)
        overlap_touch_list.append(overlap_touch)

    overlap_touch_all = pd.concat(overlap_touch_list)

    # remove duplicated geometries in overlap_touch
    overlap_touch_all = overlap_touch_all.drop_duplicates(subset=['geometry'])    # only check geometry
    overlap_touch_all.to_file('overlap_touch_all.shp')

    return overlap_touch_all

def save_grids_ids_to_txt(sel_grids_gpd, save_txt_path):
    grid_id_list = sel_grids_gpd['cell_id'].to_list()
    grid_id_list = [str(item) for item in grid_id_list]
    io_function.save_list_to_txt(save_txt_path,grid_id_list)

def save_grids_fileurl_to_txt(sel_grids_gpd, save_txt_path):
    fileurl_list = sel_grids_gpd['fileurl'].to_list()
    io_function.save_list_to_txt(save_txt_path, fileurl_list)

    ext_num_list = [ item.split('/')[1] for item in fileurl_list]
    ext_num_list = list(set(ext_num_list))
    io_function.save_list_to_txt('ext_num_list.txt', ext_num_list)

def save_arcticdata_url_dem_diff_raster(sel_grids_gpd, save_path):
    fileurl_list = sel_grids_gpd['fileurl'].to_list()
    if 'elevation-differences' in fileurl_list[0]:
        pass
    else:
        print('Warning, no elevation difference in the fireurl, skip saving the arcticdata url')
        return

    tmp_url_strs = [ item[:-7] for item in fileurl_list] # remove .tar.gz
    grid_strs = [ os.path.basename(item).split('_')[-1] for item in tmp_url_strs]

    all_urls = []
    # grid_ids_DEM_diff_grid9205.tif
    # grid_ids_date_diff_grid9205.tif
    # grid_ids_date_diff_grid9205.txt
    # grid_ids_date_diff_grid9205_newIndex.tif
    # grid_ids_date_diff_grid9205_oldIndex.tif
    # readme.txt
    for grid_str, m_url in zip(grid_strs, tmp_url_strs):
        file1 = 'grid_ids_DEM_diff_%s.tif'%grid_str
        file2 = 'grid_ids_date_diff_%s.tif'%grid_str
        file3 = 'grid_ids_date_diff_%s.txt'%grid_str
        file4 = 'grid_ids_date_diff_%s_newIndex.tif'%grid_str
        file5 = 'grid_ids_date_diff_%s_oldIndex.tif'%grid_str
        file6 = 'readme.txt'
        for filename in [file1, file2, file3, file4, file5, file6]:
            tmp_url = os.path.join(arcticdata_root,m_url,filename)
            all_urls.append(tmp_url)
    io_function.save_list_to_txt(save_path, all_urls)


def save_arcticdata_url_composited_images(sel_grids_gpd, save_path):
    fileurl_list = sel_grids_gpd['fileurl'].to_list()
    if 'composited-images' in fileurl_list[0]:
        pass
    else:
        print('Warning, no composited-images in the fireurl, skip saving the arcticdata url')
        return

    tmp_url_strs = [ item[:-7] for item in fileurl_list] # remove .tar.gz
    grid_strs = [ os.path.basename(item).split('_')[-1] for item in tmp_url_strs]

    all_urls = []
    # hillshade_HDLine_grid4574.tif
    # hillshade_HDLine_grid4574_count.tif
    # readme.txt
    for grid_str, m_url in zip(grid_strs, tmp_url_strs):
        file1 = 'hillshade_HDLine_%s.tif'%grid_str
        file2 = 'hillshade_HDLine_%s_count.tif'%grid_str
        file3 = 'readme.txt'
        for filename in [file1, file2, file3]:
            tmp_url = os.path.join(arcticdata_root,m_url,filename)
            all_urls.append(tmp_url)
    io_function.save_list_to_txt(save_path, all_urls)


def find_dem_difference_raster(sel_grids_gpd, save_txt_path):
    grid_id_list = sel_grids_gpd['cell_id'].to_list()
    fileurl_list = sel_grids_gpd['fileurl'].to_list()
    ext_num_list = [item.split('/')[1] for item in fileurl_list]

    dem_difference_file_list = []
    for grid_id, ext_num in zip(grid_id_list,ext_num_list):
        pattern = '%s*/grid_dem_diffs/*grid%d*.*'%(ext_num,grid_id)
        file_list = io_function.get_file_list_by_pattern(arcticDEM_res_dir,pattern)
        pattern2 = '%s*/grid_dem_diffs/*grid%d*_*' % (ext_num, grid_id)
        file_list2 = io_function.get_file_list_by_pattern(arcticDEM_res_dir, pattern2)
        file_list.extend(file_list2)
        if len(file_list) > 0:
            dem_difference_file_list.extend(file_list)
        else:
            print('warning: there is no files for grid %d in ext: %s'%(grid_id, ext_num))
            print(arcticDEM_res_dir, pattern)

    io_function.save_list_to_txt(save_txt_path, dem_difference_file_list)


def find_copy_dem_diff_polygons(sel_grids_gpd, save_folder):
    # find and crop polygons of elevation reductions from tarballs
    # to save disk space, I merged and compress polygons in each ext, then save to tarball

    grid_id_list = sel_grids_gpd['cell_id'].to_list()
    grid_polys = sel_grids_gpd['geometry'].to_list()
    fileurl_list = sel_grids_gpd['fileurl'].to_list()
    ext_num_list = [item.split('/')[1] for item in fileurl_list]

    ele_diff_poly_dir = os.path.join(arcticDEM_res_dir,'all_grid_dem_diffs_segment_results')
    # group grids, based on ext num
    ext_grid_ids = {}
    ext_grid_polys = {}
    for grid_id, ext_num, poly in zip(grid_id_list, ext_num_list, grid_polys):
        ext_grid_ids.setdefault(ext_num, []).append(grid_id)
        ext_grid_polys.setdefault(ext_num, []).append(poly)
    # print(ext_grids)

    # find and crop
    for key in ext_grid_ids.keys():
        print('finding and cropping polygons of elevations for ext: %s '%key)
        ext_grid_id_list = ext_grid_ids[key]
        ext_grid_list = ext_grid_polys[key]

        save_folder = '%s_grid_dem_diffs_segment_results' % key
        if os.path.isdir(save_folder):
            print('warning, folder for %s already exist, data may already there, skip'%key)
            continue

        # find the tarball
        gz_files = io_function.get_file_list_by_pattern(ele_diff_poly_dir,'%s_*dem_diff_segments*.gz'%key)
        print(gz_files)
        if len(gz_files) != 1:
            raise ValueError('Find %d gz files in %s, should be one '%(len(gz_files), ele_diff_poly_dir))
        # copy and unzip the tar ball
        gpkg_folder = io_function.unpack_tar_gz_file(gz_files[0], './')

        gpkg_path = io_function.get_file_list_by_pattern(gpkg_folder,'*.gpkg')[0]
        if os.path.isfile(gpkg_path) is False:
            raise IOError('%s does not exist'%gpkg_path)

        io_function.mkdir(save_folder)
        # crop and save polygons in each grid
        # gpd_dataframe = gpd.read_file(gpkg_path)    # this is a big file, only read once
        for g_id, g_poly in zip(ext_grid_id_list, ext_grid_list):
            print(datetime.now(), 'crop polygons for grid: %d'%g_id)
            # save_poly_path = os.path.join(save_folder,"dem_diffs_polygons_grid%d.gpkg"%g_id)
            # vector_gpd.clip_geometries(gpd_dataframe,save_poly_path,g_poly,format='GPKG')

            save_poly_path = os.path.join(save_folder, "dem_diffs_polygons_grid%d.gpkg" % g_id)
            vector_gpd.clip_geometries_ogr2ogr(gpkg_path,save_poly_path,g_poly.bounds,format='GPKG')

        # remove the gpkg
        io_function.delete_file_or_dir(gpkg_folder)

def find_copy_dem_diff_polygons_sam(sel_grids_gpd, save_folder):
    # find and crop polygons of elevation reductions (produced by Segment anything model)

    grid_id_list = sel_grids_gpd['cell_id'].to_list()
    grid_polys = sel_grids_gpd['geometry'].to_list()
    fileurl_list = sel_grids_gpd['fileurl'].to_list()
    ext_num_list = [item.split('/')[1] for item in fileurl_list]

    ele_diff_poly_dir = os.path.join(arcticDEM_res_dir,'all_grid_dem_diffs_sam_results')
    # group grids, based on ext num
    ext_grid_ids = {}
    ext_grid_polys = {}
    for grid_id, ext_num, poly in zip(grid_id_list, ext_num_list, grid_polys):
        ext_grid_ids.setdefault(ext_num, []).append(grid_id)
        ext_grid_polys.setdefault(ext_num, []).append(poly)
    # print(ext_grids)

    # find and crop
    for key in ext_grid_ids.keys():
        print('finding and cropping polygons of elevations for ext: %s '%key)
        ext_grid_id_list = ext_grid_ids[key]
        ext_grid_list = ext_grid_polys[key]

        save_folder = '%s_grid_dem_diffs_sam_results' % key
        if os.path.isdir(save_folder):
            print('warning, folder for %s already exist, data may already there, skip'%key)
            continue

        gpkg_path = io_function.get_file_list_by_pattern(ele_diff_poly_dir,'%s_*dem_diff_sam*.gpkg'%key)[0]
        if os.path.isfile(gpkg_path) is False:
            raise IOError('%s does not exist'%gpkg_path)

        io_function.mkdir(save_folder)
        for g_id, g_poly in zip(ext_grid_id_list, ext_grid_list):
            print(datetime.now(), 'crop polygons for grid: %d'%g_id)
            save_poly_path = os.path.join(save_folder, "dem_diffs_polygons_grid%d.gpkg" % g_id)
            vector_gpd.clip_geometries_ogr2ogr(gpkg_path,save_poly_path,g_poly.bounds,format='GPKG')



def test_find_grids_overlap_vector_shp():
    data_dir = os.path.expanduser('~/Data/dem_processing/products_derived_from_ArcticDEM/Index_shp')
    grid_indexes_shp = os.path.join(data_dir,'elevation-differences_Index/elevation-differences_Index.shp')
    # vector_shp = os.path.expanduser('~/codes/PycharmProjects/rts-site-data-review/rts_research_sites_point.shp')

    vector_shp = os.path.expanduser('~/Data/slump_demdiff_classify/rts_research_sites/rts_research_sites_point.shp')
    # vector_shp2 = os.path.expanduser('~/Data/slump_demdiff_classify/rts_research_sites/rts_research_sites_edit.shp')
    vector_shp2 = os.path.expanduser('~/Data/slump_demdiff_classify/rts_research_sites/rts_research_sites_arctic_small.shp')
    find_grids_overlap_vector_shp(grid_indexes_shp,[vector_shp,vector_shp2])
    # find_grids_overlap_vector_shp(grid_indexes_shp,vector_shp2)

def get_all_arcticdata_url_composited_images():
    composited_index_shp = \
        os.path.expanduser('~/Data/dem_processing/products_derived_from_ArcticDEM/Index_shp/composited-images_Index/composited-images_Index.shp')
    all_grids = gpd.read_file(composited_index_shp)
    save_arcticdata_url_composited_images(all_grids, 'all_grids_arcticdata_urls_compositedImages.txt')


def main(options, args):
    grid_indexes_shp= args[0]
    vector_shp_list = args[1:]
    if len(vector_shp_list) < 1:
        raise ValueError("No input vector files")

    global  arcticDEM_res_dir
    if options.arcticDEM_res_dir is not None:
        arcticDEM_res_dir = options.arcticDEM_res_dir

    #   cell_id (Integer64) = 58601
    #   tarball (String) = dem_diffs_2m_grid58601.tar.gz
    #   fileurl (String) = elevation-differences/ext11/dem_diffs_2m_grid58601.tar.gz
    #   POLYGON ((3320000 -120000,3340000 -120000,3340000 -140000,3320000 -140000,3320000 -120000))
    # read the grid index file
    # grid_polys, grid_attributes = vector_gpd.read_polygons_attributes_list(grid_indexes_shp,['cell_id','fileurl'],
    #                                                                        b_fix_invalid_polygon=False)

    sel_grids_gpd = find_grids_overlap_vector_shp(grid_indexes_shp, vector_shp_list)
    save_grids_ids_to_txt(sel_grids_gpd,'select_grids_ids.txt')
    save_grids_fileurl_to_txt(sel_grids_gpd,'select_grids_fileurls.txt')
    save_arcticdata_url_dem_diff_raster(sel_grids_gpd,'select_grids_arcticdata_urls_DEMdiff.txt')
    save_arcticdata_url_composited_images(sel_grids_gpd,'select_grids_arcticdata_urls_compositedImages.txt')

    find_dem_difference_raster(sel_grids_gpd, 'select_dem_diff_raster_list.txt')

    save_folder = 'ele_dem_polygons'
    # find_copy_dem_diff_polygons(sel_grids_gpd, save_folder)
    find_copy_dem_diff_polygons_sam(sel_grids_gpd, save_folder)

if __name__ == '__main__':
    usage = "usage: %prog [options] grid_indexes_shp vector_shp1 vector_shp2 ... "
    parser = OptionParser(usage=usage, version="1.0 2023-12-18")
    parser.description = 'Introduction: select grid data (DEM diff, polygons etc) if they overlap input polygons or points '

    parser.add_option("-s", "--save_dir",
                      action="store", dest="save_dir", default='select_grid_data',
                      help="the folder to save DEMs")

    parser.add_option("-d", "--arcticDEM_res_dir",
                      action="store", dest="arcticDEM_res_dir",
                      help="the folder that contains ArcticDEM results")

    # test_find_grids_overlap_vector_shp()
    # get_all_arcticdata_url_composited_images()
    # sys.exit(0)

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
