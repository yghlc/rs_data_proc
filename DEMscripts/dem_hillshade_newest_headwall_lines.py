#!/usr/bin/env python
# Filename: dem_hillshade_newest_headwall_lines.py 
"""
introduction:
1.from the multi-temporal DEM, choose the newest one on top, then create hillshade.
2. rasterize headwall lines to the hillshade.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 02 September, 2021
"""
import os,sys
import psutil
from optparse import OptionParser
import time
import re
import numpy as np

deeplabforRS = os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.timeTools as timeTools
import raster_io

from produce_DEM_diff_ArcticDEM import get_grid_20
from dem_mosaic_crop import get_dem_tif_ext_polygons
from dem_mosaic_crop import mosaic_crop_dem
from dem_mosaic_crop import mosaic_dem_list_gdal_merge
from dem_to_hillshade_slope_8bit import dem_to_hillshade

from dem_common import grid_hillshade_newest_HDLine_dir,grid_20_shp,\
    arcticDEM_reg_tif_dir,grid_dem_headwall_shp_dir,grd_hillshade_newest_on_top_dir

# some parameters
b_mosaic_id = True          # mosaic dem with the same id
b_mosaic_date = True        # mosaic dem within one day
b_mosaic_year = True        # mosaic dem for the same year, choose DEM close to July 1 on top.
b_apply_matchtag = False     # don't apply matchtag, it seems that matchtag make local value DEM after merging worse
b_mask_stripDEM_outlier = True  # mask outliers in strip DEM using the ArcticDEM tiles
b_mask_surface_water = True     # mask pixel of surface water


def draw_headwallLine_on_hillshade(hillshade_tif, headwall_line_shp,save_path):
    #
    proc = psutil.Process(os.getpid())

    lines, years = vector_gpd.read_polygons_attributes_list(headwall_line_shp,'dem_year',b_fix_invalid_polygon=False)
    sorted_years = sorted(set(years))
    if np.sum(np.isnan(sorted_years)) > 0:
        raise ValueError('dem_year: nan values in %s '% headwall_line_shp )
    print('Years:', sorted(set(years)))
    print('before rasterizing lines, used memory:', proc.memory_info()[0]/(1024*1024*1024.0),'GB')
    relative_years = [ item - 1970 for item in years]  # in computer, the date calculate from 1970 Jan 1st.
    xres, yres = raster_io.get_xres_yres_file(hillshade_tif)
    # convert lines to polygons
    line_polys = [ poly.buffer(xres) for poly in lines ]
    
    line_raster_array = raster_io.burn_polygons_to_a_raster(hillshade_tif,line_polys,relative_years,None)
    print('complete rasterizing lines, used memory:', proc.memory_info()[0]/(1024*1024*1024.0),'GB')
    print('line_raster_array, used memory:', line_raster_array.size*line_raster_array.itemsize/(1024*1024*1024.0),'GB')
    line_polys = None   # set as None, to save memory
    # assigned color
    # we choose some color from: https://www.rapidtables.com/web/color/RGB_Color.html#color-table
    # because hillshade is white, so we avoid light color
    color_table = {
        'Red':(255,0,0),
        'Lime':(0,255,0),  # bright green
        'Blue':(0,0,255),
        'Yellow':(255,255,0),
        'Cyan':(0,255,255),
        'Magenta':(255,0,255),
        'Maroon':(128,0,0),
        'Olive':(128,128,0),
        'Green':(0,128,0),
        'Purple':(128,0,128),
        'Teal':(0,128,128),
        'Navy':(0,0,128)
    }

    #from dem_date_statistics.py, we can see DEM year range from 2008 to 2017
    year_color_table = {
        2008:color_table['Red'],    # that is (255,0,0)
        2009:color_table['Lime'],
        2010:color_table['Blue'],
        2011:color_table['Yellow'],
        2012:color_table['Cyan'],
        2013:color_table['Magenta'],
        2014:color_table['Maroon'],
        2015:color_table['Olive'],
        2016:color_table['Green'],
        2017:color_table['Purple']
    }

    hillshade_np,nodata = raster_io.read_raster_all_bands_np(hillshade_tif)
    height, width = line_raster_array.shape
    img_3band = np.zeros((3,height, width),dtype=np.uint8)
    # copy hillshade to thre band
    for band in range(3):
        img_3band[band,:,:] = hillshade_np[0,:,:]

    print('created an image with three bands, used memory:', proc.memory_info()[0]/(1024*1024*1024.0),'GB')
    print('img_3band, used memory:', img_3band.size*img_3band.itemsize/(1024*1024*1024.0),'GB')

    # copy lines
    years_of_lines = np.unique(line_raster_array)
    years_of_lines = [ item + 1970 for item in years_of_lines if item > 0]
    line_raster_array = line_raster_array.astype(np.int32) + 1970    # from uint8 to year
    for year in years_of_lines:
        if year not in year_color_table.keys():
            print(year_color_table)
            raise ValueError('Year %d is not in the above color table'%year)
        # print(year)
        # print(year_color_table[year])
        loc = np.where(line_raster_array == year)
        # print(loc[0].size,width*height,loc[0].size/width*height)
        img_3band[0,loc[0],loc[1]] = year_color_table[year][0]
        img_3band[1,loc[0],loc[1]] = year_color_table[year][1]
        img_3band[2,loc[0],loc[1]] = year_color_table[year][2]


    print('complete assigning line colors, used memory:', proc.memory_info()[0]/(1024*1024*1024.0),'GB')

    # save to file
    raster_io.save_numpy_array_to_rasterfile(img_3band,save_path,hillshade_tif,nodata=None,compress='lzw',tiled='yes',bigtiff='if_safer')
    return raster_io.remove_nodata_from_raster_metadata(save_path)


def test_draw_headwallLine_on_hillshade():
    hillshade = os.path.expanduser('~/Data/dem_processing/grid_6174_tmp_files/dem_mosaic_newest_on_top_grid6174_hillshade.tif')
    headwall_line_shp = os.path.expanduser('~/Data/dem_processing/headwall_shps_grid6174/headwall_shp_multiDates_6174.shp')

    save_path = os.path.expanduser('~/Data/dem_processing/hillshade_line_grid6174.tif')
    draw_headwallLine_on_hillshade(hillshade, headwall_line_shp, save_path)


def get_headwall_line_shp(headwall_shps_dir):
    headwall_shp = io_function.get_file_list_by_ext('.shp',headwall_shps_dir,bsub_folder=True)
    # group headwall shp by id
    headwall_shp_group = {}
    for shp in headwall_shp:
        # headwall_shp_multiDates_9274.prj
        grid_id = int(re.findall('\d+',os.path.basename(shp))[0])
        if grid_id not in headwall_shp_group.keys():
            headwall_shp_group[grid_id] = shp
        else:
            raise ValueError('grid: %d already have a headwall Line shp'%grid_id)

    return headwall_shp_group

def get_hillshade_newest_top(grid_id,grid_poly,dem_ext_polys,reg_tifs,process_num, keep_dem_percent, o_res, pre_name):

    key = 'dem_mosaic_newest_on_top_grid%d' % grid_id
    hillshade_tif = os.path.join(grd_hillshade_newest_on_top_dir, key+'_hillshade.tif')
    if os.path.isfile(hillshade_tif):
        print('Warning, %s exists, skip get hillshade'%hillshade_tif)
        return hillshade_tif

    save_dir = 'grid_%d_tmp_files' % grid_id

    # check free disk space
    work_dir = './'
    free_GB = io_function.get_free_disk_space_GB(work_dir)
    total_wait_time = 0
    while free_GB < 50 and total_wait_time < 60 * 60 * 12:
        print(' The free disk space (%.4f) is less than 50 GB, wait 60 seconds' % free_GB)
        time.sleep(60)
        total_wait_time += 60
        free_GB = io_function.get_free_disk_space_GB(work_dir)

    # get subset of tifs
    dem_poly_index = vector_gpd.get_poly_index_within_extent(dem_ext_polys, grid_poly)
    if len(dem_poly_index) < 1:
        basic.outputlogMessage('warning, no dem tifs within %d grid, skip' % grid_id)
        return None

    dem_list_sub = [reg_tifs[index] for index in dem_poly_index]

    mosaic_tif_list = mosaic_crop_dem(dem_list_sub, save_dir, grid_id, grid_poly, b_mosaic_id, b_mosaic_date,
                                      process_num, keep_dem_percent, o_res, pre_name, resample_method='average',
                                      b_mask_matchtag=b_apply_matchtag, b_mask_stripDEM_outlier=b_mask_stripDEM_outlier,
                                      b_mask_surface_water=b_mask_surface_water, b_mosaic_year=b_mosaic_year)

    if len(mosaic_tif_list)<1:
        basic.outputlogMessage('warning, failed to get DEM mosaic for grid %d'%grid_id)
        return None
    # dem co-registration (cancel, the result in not good with the default setting)

    # mosaic multi-year DEM, put the newest one on the top
    dem_mosaic_yeardate = [timeTools.get_yeardate_yyyymmdd(os.path.basename(tif)[:8]) for tif in mosaic_tif_list]
    # sort, put the oldest one at the beginning, newest one at the end
    new_dem_list = [tif for _, tif in sorted(zip(dem_mosaic_yeardate, mosaic_tif_list))]
    # mosaic them using gdal_merge.py

    dem_mosaic_newest_top = mosaic_dem_list_gdal_merge(key, new_dem_list, save_dir, True)

    # convert to hillshade
    if dem_to_hillshade(dem_mosaic_newest_top, hillshade_tif) is False:
        return None

    return hillshade_tif

def combine_hillshade_headwall_line(grid_polys, grid_ids, pre_name, reg_tifs, headwall_line_shps_group, b_mosaic_id,
                                    b_mosaic_date, keep_dem_percent, o_res, process_num=1):

    dem_ext_polys = get_dem_tif_ext_polygons(reg_tifs)
    hillshade_HDLine_tifs = []
    # mosaic and crop
    for grid_id, grid_poly in zip(grid_ids, grid_polys):

        save_hillshade_headwall_path = os.path.join(grid_hillshade_newest_HDLine_dir, 'hillshade_HDLine_grid%d.tif'%grid_id)
        if os.path.isfile(save_hillshade_headwall_path):
            print('warning, %s exist, skip grid %d'%(save_hillshade_headwall_path,grid_id))
            continue

        hillshade_tif = get_hillshade_newest_top(grid_id, grid_poly, dem_ext_polys, reg_tifs, process_num,
                                                 keep_dem_percent, o_res,pre_name)
        if hillshade_tif is None:
            continue
        if grid_id not in headwall_line_shps_group.keys():
            basic.outputlogMessage('warning, grid %d does not have shapefile of headwall lines, skip'%grid_id)
            continue
        headwall_line_shp = headwall_line_shps_group[grid_id]
        # merge the hillshade and Headwall Line
        if draw_headwallLine_on_hillshade(hillshade_tif,headwall_line_shp,save_hillshade_headwall_path) is False:
            continue

        hillshade_HDLine_tifs.append(save_hillshade_headwall_path)


    return hillshade_HDLine_tifs




def get_existing_hillshade_newest_HDLine(hillshade_newest_HDLine_dir, grid_base_name, grid_ids):

    existing_grid_hillshade_HDLine_tif = []
    grid_id_no_hillshade_HDLine_tif = []
    for id in grid_ids:
        hillshade_newest_HDLine_tifs = io_function.get_file_list_by_pattern(hillshade_newest_HDLine_dir, '*_grid%d.tif' % id)
        if len(hillshade_newest_HDLine_tifs) == 1:
            existing_grid_hillshade_HDLine_tif.append(hillshade_newest_HDLine_tifs[0])
            continue
        elif len(hillshade_newest_HDLine_tifs) > 1:
            existing_grid_hillshade_HDLine_tif.append(hillshade_newest_HDLine_tifs[0])
            basic.outputlogMessage('warning, There are multiple hillshade (newest) HDLine tif for grid: %d' % id)
            for item in hillshade_newest_HDLine_tifs: basic.outputlogMessage(item)
            continue
        else:
            pass

        grid_id_no_hillshade_HDLine_tif.append(id)
    if len(existing_grid_hillshade_HDLine_tif) > 0:
        basic.outputlogMessage(
            '%d existing grid hillshade (newest) HDLine tifs for the input grid_ids or extent' % len(existing_grid_hillshade_HDLine_tif))
    else:
        basic.outputlogMessage('no existing grid hillshade (newest) HDLine tif')
    return existing_grid_hillshade_HDLine_tif, grid_id_no_hillshade_HDLine_tif


def main(options, args):
    extent_shp_or_ids_txt = args[0]
    process_num = options.process_num
    keep_dem_percent = options.keep_dem_percent
    o_res = options.out_res

    basic.setlogfile('produce_newest_hillshade_headwallLine_shp_ArcticDEM_log_%s.txt'%timeTools.get_now_time_str())

    if os.path.isdir(grid_hillshade_newest_HDLine_dir) is False:
        io_function.mkdir(grid_hillshade_newest_HDLine_dir)
    if os.path.isdir(grd_hillshade_newest_on_top_dir) is False:
        io_function.mkdir(grd_hillshade_newest_on_top_dir)

    # read grids and ids
    time0 = time.time()
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')
    print('time cost of read polygons and attributes', time.time() - time0)

    # get grid ids based on input extent
    grid_base_name = os.path.splitext(os.path.basename(extent_shp_or_ids_txt))[0]
    grid_polys, grid_ids = get_grid_20(extent_shp_or_ids_txt, all_grid_polys, all_ids)

    # check dem difference existence
    grid_hillshade_HDLine_shps, grid_id_no_hillshade_HDLine_shp = get_existing_hillshade_newest_HDLine(grid_hillshade_newest_HDLine_dir, grid_base_name, grid_ids)
    if len(grid_id_no_hillshade_HDLine_shp) > 0:
        # refine grid_polys
        if len(grid_ids) > len(grid_id_no_hillshade_HDLine_shp):
            id_index = [grid_ids.index(id) for id in grid_id_no_hillshade_HDLine_shp]
            grid_polys = [grid_polys[idx] for idx in id_index]

        # get all strip version of ArcticDEM
        reg_tifs = io_function.get_file_list_by_ext('.tif', arcticDEM_reg_tif_dir, bsub_folder=False)
        reg_tifs = [tif for tif in reg_tifs if 'matchtag' not in tif]  # remove matchtag

        # get all headwall line group, the key is the corresponding grid id
        headwall_shp_group = get_headwall_line_shp(grid_dem_headwall_shp_dir)

        # produce hillshade (newest on top), combine hillshade and headwall lines
        hillshade_HDLine_tifs = combine_hillshade_headwall_line(grid_polys,grid_id_no_hillshade_HDLine_shp,grid_base_name,reg_tifs,headwall_shp_group,
                                        b_mosaic_id,b_mosaic_date,keep_dem_percent,o_res,process_num=process_num)


    pass

if __name__ == '__main__':
    # test_draw_headwallLine_on_hillshade()
    # sys.exit(0)
    usage = "usage: %prog [options] extent_shp or grid_id_list.txt "
    parser = OptionParser(usage=usage, version="1.0 2021-3-6")
    parser.description = 'Introduction: produce DEM difference from multiple temporal ArcticDEM  '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")

    parser.add_option("-p", "--keep_dem_percent",
                      action="store", dest="keep_dem_percent", type=float, default=10.0,
                      help="keep dem with valid percentage greater than this value")

    parser.add_option("-o", "--out_res",
                      action="store", dest="out_res", type=float, default=2.0,
                      help="the resolution for final output")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
