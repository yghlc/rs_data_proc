#!/usr/bin/env python
# Filename: download_32m_arcticDEM.py 
"""
introduction: download the mosaci version of ArcticDEM


authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 29 April, 2022
"""

import os, sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.RSImageProcess as RSImageProcess


from ArcticDEM_unpack_registration import process_dem_tarball

# for '2m', should use download_arcticDEM.py to download
mosaic_res_to_name={10:'10m',32:'32m',100:'100m',500:'500m', 1000:'1km'}

def wget_file_url(url,save_path):
    # --continue, continue if the previous attempt failed
    cmd_str = 'wget --continue --no-check-certificate --output-document=%s  %s' % (save_path,url)
    status, result = basic.exec_command_string(cmd_str)
    return status, result

def download_tarball_for_one_polygon(tarball_dir,tile_tif_dir,url_head,tile_list, mosaic_res=2):
    for idx, tile_num in enumerate(tile_list):
        dem_name = tile_num + '_' + mosaic_res_to_name[mosaic_res]+'_v3.0'
        tarball_name = dem_name  + '.tar.gz'
        tiff_name = dem_name  + '_reg_dem.tif'
        tiff_path = os.path.join(tile_tif_dir,tiff_name)
        tarball_path = os.path.join(tarball_dir, tarball_name)
        # print(tiff_path)
        if os.path.isfile(tiff_path) or os.path.isfile(tarball_path):
            print('geotiff or tarball for %s already exists, skip downloading'%dem_name)
        else:
            # download
            url = url_head + tile_num + '/'+tarball_name
            wget_file_url(url,tarball_path)

        # unpack
        tar_list = [tarball_path]
        work_dir = './'
        b_rm_inter = True
        b_rm_tarball = False
        process_dem_tarball(tar_list, work_dir, tile_tif_dir, remove_inter_data=b_rm_inter, rm_tarball=b_rm_tarball,
                            apply_registration=False)


def create_a_mosaic(pre_name,extent_id,save_dir,extent_poly,tile_list,tile_tif_dir,mosaic_res=2):

    # create mosaic
    tif_list = []
    for idx, tile_num in enumerate(tile_list):
        dem_name = tile_num + '_' + mosaic_res_to_name[mosaic_res] + '_v3.0'
        tiff_name = dem_name  + '_reg_dem.tif'
        tiff_path = os.path.join(tile_tif_dir,tiff_name)
        if os.path.isfile(tiff_path) is False:
            raise ValueError('%s not exists'%tiff_path)
        tif_list.append(tiff_path)

    thread_num = 8
    output_mosaic = os.path.join(save_dir,pre_name + '_ArcticDEM_mosaic_%d'% extent_id + '.tif')
    # create mosaic, can handle only input one file, but is slow
    if os.path.isfile(output_mosaic) is False:
        result = RSImageProcess.mosaic_crop_images_gdalwarp(tif_list, output_mosaic, resampling_method='average',
                                                            o_format='GTiff',
                                                            compress='lzw', tiled='yes', bigtiff='if_safer',
                                                            thread_num=thread_num)
    else:
        print('mosaic: %s exist, skip'%output_mosaic)
    # crop
    output_crop = os.path.join(save_dir,pre_name + '_ArcticDEM_mosaic_%d_crop'% extent_id + '.tif')
    if os.path.isfile(output_crop) is False:
        RSImageProcess.subset_image_by_polygon_box_image_min(output_crop, output_mosaic, extent_poly, resample_m='average',
                                                         o_format='GTiff',
                                                         xres=mosaic_res, yres=mosaic_res, compress='lzw', tiled='yes',
                                                         bigtiff='if_safer', thread_num=thread_num)
    else:
        print('Crop: %s exist, skip'%output_crop)


def main(options, args):
    extent_shp = args[0]
    dem_index_shp = args[1]
    mosaic_res = options.mosaic_resolution
    if mosaic_res not in mosaic_res_to_name.keys():
        raise ValueError('Only accept resolutions in %s'%str([item for item in mosaic_res_to_name.keys()]))

    save_folder = options.save_dir + '_' +mosaic_res_to_name[mosaic_res]
    if os.path.isdir(save_folder) is False:
        io_function.mkdir(save_folder)
    save_folder = os.path.abspath(save_folder)  # change to absolute path

    tarball_dir = os.path.expanduser('ArcticDEM_tile_tarball')
    dem_tif_dir = os.path.expanduser('ArcticDEM_tile_geotiff')
    if os.path.isdir(tarball_dir) is False:
        io_function.mkdir(tarball_dir)
    if os.path.isdir(dem_tif_dir) is False:
        io_function.mkdir(dem_tif_dir)

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

    poly_ids = [idx for idx in range(len(extent_polys))]
    extent_name = os.path.splitext(os.path.basename(extent_shp))[0]

    # read dem polygons and tile number
    dem_polygons, dem_tiles = vector_gpd.read_polygons_attributes_list(dem_index_shp, 'tile',b_fix_invalid_polygon=False)

    for count, (idx, ext_poly) in enumerate(zip(poly_ids, extent_polys)):
        basic.outputlogMessage('get data for the %d th extent (%d/%d)' % (idx, count, len(extent_polys)))

        save_txt_path = extent_name +'-'+ 'dem_tiles_poly_%d.txt' % idx
        if os.path.isfile(save_txt_path):
            tiles = io_function.read_list_from_txt(save_txt_path)
            basic.outputlogMessage('read %d dem tiles from %s' % (len(tiles),save_txt_path))
        else:
            # get fileurl
            dem_poly_ids = vector_gpd.get_poly_index_within_extent(dem_polygons,ext_poly)
            basic.outputlogMessage('find %d DEM within %d th extent' % (len(dem_poly_ids), (idx)))
            tiles = [dem_tiles[id] for id in dem_poly_ids]

            # save to txt
            io_function.save_list_to_txt(save_txt_path, tiles)
            basic.outputlogMessage('save dem urls to %s' % save_txt_path)

        # download and create a mosaic
        url_head = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/%s/'%mosaic_res_to_name[mosaic_res]
        download_tarball_for_one_polygon(tarball_dir, dem_tif_dir, url_head, tiles,mosaic_res=mosaic_res)

        # create a mosaic
        extent_name_idx = extent_name + '_%d'%idx
        create_a_mosaic(extent_name_idx, idx, save_folder, ext_poly,tiles,dem_tif_dir,mosaic_res=mosaic_res)



if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp dem_indexes_shp"
    parser = OptionParser(usage=usage, version="1.0 2020-12-25")
    parser.description = 'Introduction: download ArcticDEM mosaic within an extent at other resolution '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir", default='ArcticDEM_mosaic',
                      help="the folder to save DEMs")

    parser.add_option("-r", "--mosaic_resolution",
                      action="store", dest="mosaic_resolution", type=int, default=10,
                      help="the resolution of mosaic to download: 10, 32, 100, 500,1000")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
