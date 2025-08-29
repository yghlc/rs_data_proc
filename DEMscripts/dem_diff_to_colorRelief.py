#!/usr/bin/env python
# Filename: ArcticDEM_proc_grid.py
"""
introduction: Convert DEM difference to color relief

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 21 March, 2024
"""

import os,sys
from optparse import OptionParser
import time
machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
import raster_io

import numpy as np
import rasterio
from matplotlib.colors import LinearSegmentedColormap
# some folder paths

from dem_common import grid_dem_diffs_dir, grid_dem_diffs_color_dir

def dem_tif_to_colorReleif(input,output,out_format='GTiff',tif_compression='lzw'):
    # change the file extension
    file_extension = raster_io.get_file_extension(out_format)
    file_path, ext1 = os.path.splitext(output)
    output = file_path + file_extension

    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return True

    color_text_file = 'dem_diff_color_5to5m.txt'

    if out_format=='GTiff':
        command_str = f'gdaldem color-relief -of {out_format} -co compress={tif_compression} -co tiled=yes -co bigtiff=if_safer '
    else:
        command_str = f'gdaldem color-relief -of {out_format} '

    command_str +=  input + ' ' + ' %s '%color_text_file  + output

    print(command_str)
    res = os.system(command_str)
    # basic.os_system_exit_code(command_str)
    if res == 0:
        return True
    else:
        return False

def parse_qgis_color_file(txt_path):
    values = []
    colors = []
    nodata_color = None
    with open(txt_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            if parts[0].lower() == 'nv':
                nodata_color = [int(x) for x in parts[1:5]]  # RGBA
                continue
            # QGIS lines: value,R,G,B,A,value
            value = float(parts[0])
            rgb = [int(parts[1]), int(parts[2]), int(parts[3])]
            values.append(value)
            colors.append(rgb)
    return np.array(values), np.array(colors) / 255.0, nodata_color

def make_linear_colormap(values, rgbs):
    vmin, vmax = values.min(), values.max()
    norm_values = (values - vmin) / (vmax - vmin)
    cdict = {'red': [], 'green': [], 'blue': []}
    for v, rgb in zip(norm_values, rgbs):
        cdict['red'].append((v, rgb[0], rgb[0]))
        cdict['green'].append((v, rgb[1], rgb[1]))
        cdict['blue'].append((v, rgb[2], rgb[2]))
    return LinearSegmentedColormap('custom_cmap', cdict)

def colormap_numpy(dem, values, colors, nodata_mask=None, nodata_color=None):
    # dem: (H, W) array
    # values: (N,) array of breakpoints
    # colors: (N, 3) array of RGB (0-1 float)
    dem_flat = dem.ravel()
    rgb_flat = np.empty((dem_flat.size, 3), dtype=np.float32)
    for i in range(3):  # R, G, B
        rgb_flat[:, i] = np.interp(dem_flat, values, colors[:, i])
    rgb = rgb_flat.reshape(*dem.shape, 3)
    rgb = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)
    if nodata_mask is not None and nodata_color is not None:
        rgb[nodata_mask] = nodata_color[:3]
    return rgb

def colormap_lut(dem, values, colors, nodata_mask=None, nodata_color=None, bins=65536):
    # Map the DEM range to bins
    vmin, vmax = values.min(), values.max()
    lut = np.empty((bins, 3), dtype=np.uint8)
    # For each channel, interpolate the color for each bin center
    bin_centers = np.linspace(vmin, vmax, bins)
    for i in range(3):  # R, G, B
        lut[:, i] = np.clip(np.interp(bin_centers, values, colors[:, i]) * 255, 0, 255).astype(np.uint8)
    # Map dem values to lut indices
    dem_clip = np.clip(dem, vmin, vmax)
    indices = ((dem_clip - vmin) / (vmax - vmin) * (bins - 1)).astype(np.int32)
    rgb = lut[indices]
    if nodata_mask is not None and nodata_color is not None:
        rgb[nodata_mask] = nodata_color[:3]
    return rgb

def colormap_small_lut(dem, values, colors, nodata_mask=None, nodata_color=None):
    vmin = int(np.floor(values.min()))
    vmax = int(np.ceil(values.max()))
    lut_size = vmax - vmin + 1
    lut = np.empty((lut_size, 3), dtype=np.uint8)
    bin_centers = np.arange(vmin, vmax + 1)
    for i in range(3):
        lut[:, i] = np.clip(np.interp(bin_centers, values, colors[:, i]) * 255, 0, 255).astype(np.uint8)
    indices = np.clip(dem - vmin, 0, lut_size - 1)
    rgb = lut[indices]
    if nodata_mask is not None and nodata_color is not None:
        rgb[nodata_mask] = nodata_color[:3]
    return rgb

def build_small_lut_colormap(values, colors):
    vmin = int(np.floor(values.min()))
    vmax = int(np.ceil(values.max()))
    lut_size = vmax - vmin + 1
    lut = np.empty((lut_size, 3), dtype=np.uint8)
    bin_centers = np.arange(vmin, vmax + 1)
    for i in range(3):
        lut[:, i] = np.clip(np.interp(bin_centers, values, colors[:, i]) * 255, 0, 255).astype(np.uint8)
    return lut, vmin, lut_size

def colormap_direct(dem, values, colors, nodata_mask=None, nodata_color=None):
    # dem: (H, W) int16, values: (56,), colors: (56, 3)
    h, w = dem.shape
    rgb = np.empty((h, w, 3), dtype=np.uint8)
    for i in range(3):
        channel = np.interp(dem, values, colors[:, i]) * 255
        rgb[:, :, i] = np.clip(channel, 0, 255).astype(np.uint8)
    if nodata_mask is not None and nodata_color is not None:
        rgb[nodata_mask] = nodata_color[:3]
    return rgb

def dem_tif_to_colorReleif_npArray(input, color_relief_txt):
    # t0=time.time()
    values, rgbs, nodata_color = parse_qgis_color_file(color_relief_txt)
    # t1 = time.time()
    with rasterio.open(input) as src:
        dem = src.read(1)
        nodata = src.nodata
    nodata_mask = None
    # t2 = time.time()
    if nodata is not None:
        nodata_mask = (dem == nodata)
    # rgb = colormap_numpy(dem, values, rgbs, nodata_mask, nodata_color)
    # rgb = colormap_lut(dem, values, rgbs, nodata_mask, nodata_color)
    # rgb = colormap_small_lut(dem, values, rgbs, nodata_mask, nodata_color)  # this is the most efficent.

    # built color map
    lut, vmin, lut_size = build_small_lut_colormap(values,rgbs)
    # t3 = time.time()

    # to rgb
    indices = np.clip(dem - vmin, 0, lut_size - 1)
    rgb = lut[indices]
    if nodata_mask is not None and nodata_color is not None:
        rgb[nodata_mask] = nodata_color[:3]

    # rgb = colormap_direct(dem, values, rgbs, mask, nodata_color)
    # t4 = time.time()
    # print('time cost in dem_tif_to_colorReleif_npArray', t1-t0, t2-t1, t3-t2, t4-t3)
    return rgb  # its shape is: height,width, band_count


def test_dem_tif_to_colorReleif_npArray():
    from datetime import datetime
    import time
    data_dir = os.path.expanduser('~/Downloads/tmp')
    dem_tif = os.path.join(data_dir,'grid_ids_DEM_diff_grid44722.tif')
    color_txt = os.path.join(data_dir,'dem_diff_color_5to5m.txt')
    print(datetime.now(), 'start processing ')
    t0 = time.time()
    rgb_array = dem_tif_to_colorReleif_npArray(dem_tif, color_txt)

    t1 = time.time()
    print(datetime.now(),'complete converting to RGB')

    # for this data, colormap_small_lut is the most efficient.

    rgb_array = rgb_array.transpose(2, 0, 1)    # to band_count,height,width, for saving in rasterio.
    save_rgb = os.path.join(data_dir,'grid_ids_DEM_diff_grid44722_rgb.tif')
    raster_io.save_numpy_array_to_rasterfile(rgb_array,save_rgb,dem_tif,nodata=255)
    t2 = time.time()

    print(datetime.now(), 'saved files ')
    print('to RGB time cost', t1-t0 )
    print('save file time cost', t2-t1 )



# def test_dem_tif_to_colorReleif():
#     dem_diff_list = io_function.get_file_list_by_pattern('./','*.tif')
#     count = len(dem_diff_list)
#     for idx, tif in enumerate(dem_diff_list):
#         print('%d/%d convert %s to color relief '%(idx+1, count, tif))
#         tif_color = io_function.get_name_by_adding_tail(tif, 'color')
#         output = os.path.join(grid_dem_diffs_color_dir, os.path.basename(tif_color))
#         dem_tif_to_colorReleif(tif,output)

def test_dem_tif_to_colorReleif_one():

    tif = os.path.expanduser('~/Data/dem_processing/grid_dem_diffs/grid_ids_DEM_diff_grid13965.tif')
    tif_color = io_function.get_name_by_adding_tail(tif, 'color')
    # output = os.path.join(grid_dem_diffs_color_dir, os.path.basename(tif_color))
    dem_tif_to_colorReleif(tif,tif_color)

def one_dem_diff_to_colorRelief(demDiff_tif):
    if os.path.isdir(grid_dem_diffs_color_dir) is False:
        io_function.mkdir(grid_dem_diffs_color_dir)
    tif_color = io_function.get_name_by_adding_tail(demDiff_tif, 'color')
    output = os.path.join(grid_dem_diffs_color_dir, os.path.basename(tif_color))
    if dem_tif_to_colorReleif(demDiff_tif, output) is False:
        basic.outputlogMessage('failed to generate color relief from DEM difference')
        return False
    return True

def main(options, args):
    basic.setlogfile('log_convert_dem_diff_to_colorRelief.txt')
    out_format = options.out_format
    tif_compression = options.tif_compression
    out_dir = options.out_dir
    file_pattern = options.file_pattern

    if len(args) < 1:
        dem_diff_list = io_function.get_file_list_by_pattern(grid_dem_diffs_dir, file_pattern)
        if os.path.isdir(grid_dem_diffs_color_dir) is False:
            io_function.mkdir(grid_dem_diffs_color_dir)
    else:
        dem_diff_file_or_dir = args[0]
        if os.path.isfile(dem_diff_file_or_dir):
            dem_diff_list = [dem_diff_file_or_dir]
        else:
            dem_diff_list = io_function.get_file_list_by_pattern(dem_diff_file_or_dir,file_pattern)
            if len(dem_diff_list) < 1:
                basic.outputlogMessage(f'No DEM diff files in {dem_diff_file_or_dir} for colorRelief')
    # if out_dir is not set, using the one in dem_common.py
    if out_dir is None:
        out_dir = grid_dem_diffs_color_dir if os.path.isdir(grid_dem_diffs_color_dir) else './'

    count = len(dem_diff_list)
    failed_tifs = []
    for idx, tif in enumerate(dem_diff_list):
        print('%d/%d convert %s to Color Relief'%(idx+1, count, tif))
        tif_color = io_function.get_name_by_adding_tail(tif, 'color')
        output = os.path.join(out_dir, os.path.basename(tif_color))
        if dem_tif_to_colorReleif(tif,output,out_format=out_format,tif_compression=tif_compression) is False:
            failed_tifs.append(tif)

    if len(failed_tifs)>0:
        io_function.save_list_to_txt('failed_dem_diff_to_color.txt',failed_tifs)


if __name__ == '__main__':
    usage = "usage: %prog [options] dem_diff or dem_diff_dir "
    parser = OptionParser(usage=usage, version="1.0 2024-3-21")
    parser.description = 'Introduction: producing DEM color relief '

    parser.add_option("-f", "--out_format",
                      action="store", dest="out_format",default='GTIFF',
                      help="the format of output images, GTIFF, PNG, JPEG, VRT, etc")

    parser.add_option("-o", "--out_dir",
                      action="store", dest="out_dir",
                      help="the save directory for files")

    parser.add_option("-c", "--tif_compression",
                      action="store", dest="tif_compression",default='lzw',
                      help="the compression for tif format, JPEG compression reduce file size, "
                           "but JPEG format reduce file size more but loss some info")

    parser.add_option("-p", "--file_pattern",
                      action="store", dest="file_pattern",default='*DEM_diff_grid*.tif',
                      help="the file name pattern for search rasters in a folder ")

    (options, args) = parser.parse_args()
    main(options, args)
    # test_dem_tif_to_colorReleif_one()
    # test_dem_tif_to_colorReleif_npArray()
