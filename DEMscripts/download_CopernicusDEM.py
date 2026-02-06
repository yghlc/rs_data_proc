#!/usr/bin/env python
# Filename: download_CopernicusDEM.py 
"""
introduction: download Copernicus DEM: https://registry.opendata.aws/copernicus-dem/

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 15 May, 2023
"""
import os, sys
from optparse import OptionParser
from datetime import datetime

deeplabforRS = os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.io_function as io_function
import basic_src.basic as basic
import math

aws_s3 = 's3://copernicus-dem-30m/'

def extent_to_1degree_tiles(poly):
    minx, miny, maxx, maxy = poly.bounds  # (minx, miny, maxx, maxy)
    x_list = [ item for item in range(math.floor(minx), math.ceil(maxx)) ]
    y_list = [ item for item in range(math.floor(miny), math.ceil(maxy)) ]
    # Copernicus_DSM_COG_10_N59_00_E026_00_DEM
    tiles = []
    for y in y_list:
        lat_str = str(abs(y)).zfill(2)
        lat_str = 'N'+ lat_str if y > 0 else 'S'+lat_str
        for x in x_list:
            lon_str = str(abs(x)).zfill(3)
            lon_str = 'E' + lon_str if x > 0 else 'W' + lon_str
            # print(lon_str, lat_str)
            # tiles.append(lat_str+lon_str+'.SRTMGL1.hgt.zip')
            tiles.append('Copernicus_DSM_COG_10_%s_00_%s_00_DEM'%(lat_str, lon_str))

    return tiles

def aws_CopernicusDEM_download(tiles, cache_dir):
    if os.path.isdir(cache_dir) is False:
        io_function.mkdir(cache_dir)

    download_tiles = []

    for idx, tile in enumerate(tiles):
        print(datetime.now(),'(%d/%d): downloading %s'%(idx+1, len(tiles), tile))
        loc_str = aws_s3 + tile
        save_path = os.path.join(cache_dir,tile)
        if os.path.isfile(os.path.join(save_path,tile+'.tif')):
            print('%s exists, skip downloading'%tile)
            download_tiles.append(tile)
            continue
        cmd_str = 'aws s3 cp --no-sign-request --recursive %s %s'%(loc_str, save_path)
        print(cmd_str)
        basic.os_system_exit_code(cmd_str)
        download_tiles.append(tile)

    return download_tiles

def process_CopernicusDEM_tiles(cache_dir, download_tiles, save_path):
    if os.path.isfile(save_path):
        print('warning, %s exists, skip' % save_path)
        return False
    # unpack
    tile_tif_paths = [os.path.join(cache_dir, item, item + '.tif') for item in download_tiles]
    curr_dir = os.getcwd()

    # merge  # -n -32768 -a_nodata -32768
    cmd_str = 'gdal_merge.py -co COMPRESS=DEFLATE -co PREDICTOR=3 -o %s ' % save_path

    for tile in tile_tif_paths:
        if os.path.isfile(tile) is False:
            print('warning, %s does not exist' % tile)
            continue
        cmd_str += ' %s'%tile
    print('merging tiles')
    basic.os_system_exit_code(cmd_str)

    return save_path


def download_CopernicusDEM(extent_shp, save_path, cache_dir):
    '''
    Download Copernicus DEM 30m elevation tiles from AWS
    :param extent_shp: shapefile
    :param save_path: output path
    :param cache_dir:
    :return:
    '''
    # shapefile to 1 by 1 degrees.
    ext_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp, 'EPSG:4326')
    if len(ext_polys) == 1:
        # create file names
        tiles = extent_to_1degree_tiles(ext_polys[0])
        download_tiles = aws_CopernicusDEM_download(tiles, cache_dir)
        if len(download_tiles) > 0:
            process_CopernicusDEM_tiles(cache_dir, download_tiles, save_path)
        else:
            print('error, NO downloaded CopernicusDEM tiles')
    elif len(ext_polys) > 1:
        raise ValueError('currently, only support one polygon')
        pass
    else:
        raise ValueError('No extent polygons in %s' % extent_shp)


def main(options, args):
    extent_shp = args[0]
    assert os.path.isfile(extent_shp)

    save_path = options.save_path
    cache_dir = options.cache_dir
    if save_path is None:
        file_name = os.path.splitext(os.path.basename(extent_shp))[0]
        save_path = os.path.join(os.getcwd(), file_name + '_CoDEM.tif')
    else:
        save_path = os.path.abspath(save_path)
    print(datetime.now(), 'download Copernicus DEM, \n will save to %s'%save_path)
    if os.path.isfile(save_path):
        print("%s already exists, skip"%save_path)
        return

    if cache_dir is None:
        cache_dir = os.path.expanduser('~/Data/CopernicusDEM')
    if os.path.isdir(cache_dir) is False:
        os.makedirs(cache_dir)

    download_CopernicusDEM(extent_shp, save_path, cache_dir)

if __name__ == "__main__":
    usage = "usage: %prog [options] extent_shp "
    parser = OptionParser(usage=usage, version="1.0 2023-05-15")
    parser.description = 'Introduction: download Copernicus DEM  '

    parser.add_option("-d", "--save_path",
                      action="store", dest="save_path",
                      help="the path for saving a DEM file")

    parser.add_option("-a", "--cache_dir",
                      action="store", dest="cache_dir",
                      help="the cache directory")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)


