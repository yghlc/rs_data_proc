#!/usr/bin/env python
# Filename: cal_volume_change.py 
"""
introduction: calculate the volumetric change from elevation difference, or estimate the volume from single DEM

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 04 January, 2026
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.basic as basic
import basic_src.io_function as io_function

from datetime import datetime
import raster_statistic

import numpy as np
from multiprocessing import Pool

from skimage.transform import resize

def cal_volume(in_array, nodata,res_x, res_y,v_range=None, scale=1.0):
    # print(in_array)
    # print(in_array.shape)
    data_1d = in_array.flatten()
    # print(data_1d)
    total_size =  data_1d.size
    data_1d = data_1d[ data_1d != nodata]
    size_2 = data_1d.size
    data_1d = data_1d[~np.isnan(data_1d)]  # remove nan value
    size_3 = data_1d.size

    nodata_count = total_size - size_2
    nan_count = size_2 - size_3

    # apply the scale before checking range
    data_1d = data_1d * scale

    if v_range is not None:
        lower = v_range[0]
        upper = v_range[1]
        if lower is None:
            data_1d = data_1d[data_1d <= upper]
        elif upper is None:
            data_1d = data_1d[data_1d >= lower]
        else:
            data_1d = data_1d[np.logical_and( data_1d >= lower, data_1d <= upper ) ]
    size_4 = data_1d.size
    valid_count = total_size - nodata_count - nan_count
    not_in_range_c = size_3 - size_4

    if data_1d.size < 1:
        return 0, 0, 0, 0, 0

    # print(data_1d)
    volume = np.sum(np.abs(data_1d)*res_x*res_y)

    # print(nodata)
    # print(total_size,size_2,size_3,size_4)
    # print(valid_count, not_in_range_c)

    return volume, valid_count, not_in_range_c, np.min(data_1d), np.max(data_1d)


def resample_raster_to_match_another(src_array, target_array, method='nearest'):
    '''resample a raster array to match another raster array, the target array is used as reference using skimage '''
    # 0 = nearest, 1 = bilinear

    # If your raster is categorical, such as land-cover class IDs, then do not use anti-aliasing, and use nearest neighbor
    # Note that when down-sampling an image, resize and rescale should perform Gaussian smoothing to avoid aliasing artifacts 
    # so, aliasing artifacts should be True 

    if method == 'nearest':
        resampled = resize(src_array, target_array.shape, order=0, preserve_range=True, anti_aliasing=False )
    elif method == 'bilinear':
        resampled = resize(src_array, target_array.shape, order=1, preserve_range=True, anti_aliasing=True )
    else:
        raise ValueError('method should be nearest or bilinear') 

    return resampled


def extrapolation_carbon_deeper_layers_lookup(carbon_mobilization_layers, max_num_layers, extrapolation='linear'):
    '''
    Extrapolate carbon values for deeper layers beyond the first three.
    :param carbon_mobilization_layers: the 3D array of carbon mobilization layers
    :param max_num_layers: the total number of layers to extrapolate
    :param extrapolation: set "linear" or "exponential_decay", will use the function to estimate carbon values beyond 300 cm
    :return: the updated carbon_mobilization_layers with extrapolated values
    '''
    # because most of the carbon values for are the same, so using lookup table for extrapolation is efficient
    # calculate the look up table for extrapolation based on the chosen method
    # find unique values for the first three layers, exlcuding nan values
    unique_values_layer1 = np.unique(carbon_mobilization_layers[:1, :, :])
    unique_values_layer1 = unique_values_layer1[~np.isnan(unique_values_layer1)]
    unique_values_layer2 = np.unique(carbon_mobilization_layers[1:2, :, :])
    unique_values_layer2 = unique_values_layer2[~np.isnan(unique_values_layer2)]
    unique_values_layer3 = np.unique(carbon_mobilization_layers[2:3, :, :])
    unique_values_layer3 = unique_values_layer3[~np.isnan(unique_values_layer3)]
    # print(unique_values_layer1, unique_values_layer2, unique_values_layer3)
    # create a dictionary to store the lookup table for each unique combination of the first three layers
    lookup_dict = {} 
    for val1 in unique_values_layer1:
        for val2 in unique_values_layer2:
            for val3 in unique_values_layer3:
                x = np.array([1, 2, 3])
                y = np.array([val1, val2, val3])
                x_new = np.arange(4, max_num_layers + 1) # for layers beyond the first three
                if extrapolation == 'linear':
                    # y=mx+b
                    m, b = np.polyfit(x, y, 1)
                    y_linear = m * x_new + b
                    # if smaller than 0, then set a 0
                    y_linear[y_linear < 0] = 0
                    lookup_dict[(val1, val2, val3)] = y_linear
                elif extrapolation == 'exponential_decay':
                    # y = A * exp(-k*x)
                    if np.any(y <= 0):
                        raise ValueError("Exponential decay fit requires all y values to be positive.")
                    slope, intercept = np.polyfit(x, np.log(y), 1)
                    A = np.exp(intercept)
                    k = -slope
                    y_exp = A * np.exp(-k * x_new)
                    lookup_dict[(val1, val2, val3)] = y_exp
                else:
                    raise ValueError('extrapolation should be linear or exponential_decay')

    # print("lookup_dict:", lookup_dict)
    return lookup_dict


def cal_carbon_mobilization(dem_diff_array,dem_diff_nodata, dem_diff_res, carbon_array, carbon_res, dem_diff_range=None, 
                            dem_diff_scale=1.0, layer_depth=1.0, extrapolation='linear'):
    '''
    calculate the carbon mobilization based on the volume change and carbon data
    :param dem_diff_array: the DEM difference array
    :param dem_diff_nodata: the nodata value of the DEM difference array
    :param dem_diff_res: the resolution of the DEM difference array
    :param carbon_array: the carbon array
    :param carbon_res: the original resolution of the carbon array (e.g., 1000 m)
    :param dem_diff_range: the range of the DEM difference values to consider for volume change
    :param dem_diff_scale: the scale factor to apply to the DEM difference values
    :param layer_depth: the depth of each layer to divide the DEM difference values into (default is 1.0 m)
    :param extrapolation: set "linear" or "exponential_decay", will use the function to estimate carbon values beyond 300 cm 
    :return: the total carbon mobilization
    '''

    # not nondata, not nan, and in the range of dem_diff_range
    valid_mask = ~np.isnan(dem_diff_array)
    if dem_diff_nodata is not None:
        valid_mask &= dem_diff_array != dem_diff_nodata
    # apply the scale to the dem diff array before applying range, as the dem diff array may be in cm, and we need to convert it to m
    dem_diff_array = dem_diff_array * dem_diff_scale
    if dem_diff_range is not None:
        lower, upper = dem_diff_range
        lower = -np.inf if lower is None else lower
        upper = np.inf if upper is None else upper
        valid_mask &= dem_diff_array >= lower
        valid_mask &= dem_diff_array <= upper
    # don't select, just mask out, keep the original shape, making each to check and debug
    # valid_locs = np.where(valid_mask)

    # print("dem_diff_array shape", dem_diff_array.shape, "carbon_array shape", carbon_array.shape)
    # set the non valid pixel as 0, so it will not be used in the calculation
    dem_diff_array[~valid_mask] = 0
    # select the valid locations from carbon_array, assuming carbon_array is 3D (band, height, width)  
    carbon_array[:, ~valid_mask] = np.nan

    # adjust the carbon array to the same resolution as dem diff array
    carbon_array = carbon_array * ((1.0*dem_diff_res / carbon_res) ** 2)  

    # calculate the layer number for each pixel, for nan values, we can set the layer number to 0, as they will not be used in the calculation 
    num_layer_each_pixel = np.floor((np.abs(dem_diff_array)+ 0.5*layer_depth) / layer_depth).astype(int)
    # print("num_layer_each_pixel:", num_layer_each_pixel)
    # print("max number of layers:", np.max(num_layer_each_pixel), ", min number of layers:", np.min(num_layer_each_pixel))
    max_num_layers = np.max(num_layer_each_pixel) 

    # create a 3D mask array for carbon_mobilization_layers based on num_layer_each_pixel and valid_mask
    valid_mask_carbon = (np.arange(max_num_layers)[:, None, None] < num_layer_each_pixel[None, :, :]) & valid_mask[None, :, :]
    # print("valid_mask_carbon shape", valid_mask_carbon.shape)

    # create a 3D array to store the carbon mobilization for each layer, the shape
    carbon_mobilization_layers = np.zeros((max_num_layers, dem_diff_array.shape[0],dem_diff_array.shape[1]), dtype=np.float32)
    # print("carbon_mobilization_layers shape", carbon_mobilization_layers.shape)

    # copy the first three layers from carbon_array to carbon_mobilization_layers directly, 
    for layer in range(min(3, max_num_layers)):
        carbon_mobilization_layers[layer, :, :] = carbon_array[layer, :, :]
    if max_num_layers > 3:
        # print("Extrapolating carbon values for layers beyond the first three")
        # for layers beyond the first three, we need to use the extrapolation function to estimate the carbon values
        deep_carbon_lookup = extrapolation_carbon_deeper_layers_lookup(carbon_mobilization_layers, max_num_layers, extrapolation=extrapolation)
        # fill in the deeper layers using the lookup
        for layer in range(3, max_num_layers):
            for i in range(carbon_mobilization_layers.shape[1]):
                for j in range(carbon_mobilization_layers.shape[2]):
                    key = tuple(carbon_mobilization_layers[:3, i, j])
                    if key in deep_carbon_lookup:
                        carbon_mobilization_layers[layer, i, j] = deep_carbon_lookup[key][layer - 3]

    # mask the carbon mobilization layers using the valid_mask_carbon
    carbon_mobilization_layers[~valid_mask_carbon] = 0

    # sum the carbon mobilization across all layers
    total_carbon = np.sum(carbon_mobilization_layers)
    # print("total_carbon", total_carbon)

    # sys.exit(0)  # for testing, exit after processing the first polygon


    return total_carbon




def cal_volume_one_polygon(idx, polygon, res,image_tiles, img_tile_polygons, carbon_tiles,carbon_tile_polygons, nodata=None,v_range=None,
                           all_touched=True,tile_min_overlap=None,scale=1.0):
    ''' calculate volume for a polygon, allow multiple polygons
    '''
    # output progress
    if idx%10 == 0:
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Working on poly:", idx)
    # # for testing, only process polygon with idx 52
    # if idx != 52:
    #     return 0, 0, 0, 0, 0, 0
    # read the raster
    out_image, out_transform,_ = raster_statistic.read_raster_data_img_list_geo(idx, polygon, image_tiles, img_tile_polygons, nodata=nodata,
                                          band=1, all_touched=all_touched, tile_min_overlap=tile_min_overlap)
    out_volume, pixel_count, not_in_range_c, min_v, max_v = cal_volume(out_image, nodata,  res, res, v_range=v_range,scale=scale)
    
    total_carbon = None
    if carbon_tiles is None:
        # NO carbon raster, return None for carbon data
        return out_volume, pixel_count, not_in_range_c, min_v, max_v, total_carbon
    else:
        carbon_nodata = raster_io.get_nodata(carbon_tiles[0])
        carbon_res, _ = raster_io.get_xres_yres_file(carbon_tiles[0])  # assume square pixel, so use xres
        _, _, carbon_band_count, _ = raster_io.get_height_width_bandnum_dtype(carbon_tiles[0])
        # print('carbon raster band count:', carbon_band_count)
        # read all bands
        carbon_data_list = []
        for band in range(1, carbon_band_count + 1):
            carbon_data, carbon_transform, _ = raster_statistic.read_raster_data_img_list_geo(idx, polygon, carbon_tiles, carbon_tile_polygons, nodata=carbon_nodata,
                                          band=band, all_touched=all_touched, tile_min_overlap=tile_min_overlap)
            # resample carbon raster to the same shape as dem diff raster
            if carbon_data.shape != out_image.shape:
                carbon_data = resample_raster_to_match_another(carbon_data, out_image, method='nearest')
            carbon_data_list.append(carbon_data)
        
        # combine all bands to a 3D array, for rasterio, the shape should be (band, height, width), so we need to stack along axis=0
        if len(carbon_data_list) > 1:
            carbon_data = np.stack(carbon_data_list, axis=0)
        else:
            carbon_data = carbon_data_list[0]

        # # check the resmapled data by saving it to a file, and check if the resampled data is correct
        # # as they crop to the same polygon, the transform should be the same, so we can use the transform from dem diff raster
        # raster_io.save_numpy_array_to_rasterfile(out_image, f'out_image_{idx}.tif', image_tiles[0],  nodata=nodata, transform=out_transform)
        # raster_io.save_numpy_array_to_rasterfile(carbon_data, f'carbon_data_resampled_{idx}.tif', carbon_tiles[0],  nodata=carbon_nodata, transform=out_transform)

        # calculate the total carbon mobilization based on the volume change and carbon data
        total_carbon = cal_carbon_mobilization(out_image,nodata, res, carbon_data, carbon_res, dem_diff_range=v_range, dem_diff_scale=scale)
        

        # calculate the volume change
        return out_volume, pixel_count, not_in_range_c, min_v, max_v, total_carbon

def cal_volumetric_change_from_dem_diff(boundary_shp,raster_file_or_files,carbon_raster=None, save_path=None,v_range=None,
                                        tile_min_overlap=None, nodata=None,buffer=None, all_touched=True, process_num=1):
    '''
    calculate the volumetric changes of all polygons in a vector file, allow multiple rasters
    :param boundary_shp:
    :param raster_file_or_files:
    :param carbon_raster: if not None, it will be used to calculate the carbon mobilization based on the volume change, 
                          if change carbon to soil, ice or others, the results should be related to soil or ice or others,
    :param save_path:
    :param v_range:
    :param tile_min_overlap:
    :param nodata: incorrect setting of nodata will end in low values of "INrangePer".
    :param buffer:
    :param all_touched:
    :param process_num:
    :return:
    '''

    io_function.is_file_exist(boundary_shp)
    if isinstance(raster_file_or_files,str):
        io_function.is_file_exist(raster_file_or_files)
        image_tiles = [raster_file_or_files]
    elif isinstance(raster_file_or_files,list):
        image_tiles = raster_file_or_files
    else:
        raise ValueError('unsupport type for %s'%str(raster_file_or_files))
    raster_prj = raster_io.get_projection(image_tiles[0])
    # check projection (assume we have the same projection), check them outside this function

    # get image box
    img_tile_boxes = [raster_io.get_image_bound_box(tile) for tile in image_tiles]
    img_tile_polygons = [vector_gpd.convert_image_bound_to_shapely_polygon(box) for box in img_tile_boxes]
    polygons = vector_gpd.read_polygons_gpd(boundary_shp)
    if len(polygons) < 1:
        basic.outputlogMessage('No polygons in %s' % boundary_shp)
        return False
    # polygons_json = [mapping(item) for item in polygons]  # no need when use new verion of rasterio
    if buffer is not None:
        polygons = [poly.buffer(buffer) for poly in polygons]

    res_x, res_y = raster_io.get_xres_yres_file(image_tiles[0])
    heigth, width, bandcount, dtype = raster_io.get_height_width_bandnum_dtype(image_tiles[0])
    scale = 1.0
    if dtype == 'int16':
        scale = 0.01    # DEM diff saved in 16 bit, that is cm, multiply by 0.01 to m.

    volume_res_list = []
    in_range_per = []
    dem_diff_min = []
    dem_diff_max = []
    carbon_move_list = []

    # check carbon raster file
    if carbon_raster is not None:
        io_function.is_file_exist(carbon_raster)
        # check projection
        carbon_r_prj = raster_io.get_projection(carbon_raster)
        if raster_prj != carbon_r_prj:
            raise ValueError('The projection of carbon raster and dem diff raster are not the same')
        carbon_tiles = [carbon_raster]
        carbon_tile_boxes = [raster_io.get_image_bound_box(tile) for tile in carbon_tiles]
        carbon_tile_polygons = [vector_gpd.convert_image_bound_to_shapely_polygon(box) for box in carbon_tile_boxes]
    else:
        carbon_tiles = None
        carbon_tile_polygons = None
        

    # process polygons one by one and the corresponding image tiles (parallel and save memory)
    # also to avoid error: daemonic processes are not allowed to have children
    if process_num == 1:

        for idx, polygon in enumerate(polygons):
            out_volume, pixel_count, not_in_range_c, min_v, max_v, carbon_move = cal_volume_one_polygon(idx, polygon, res_x, image_tiles, img_tile_polygons,carbon_tiles,carbon_tile_polygons,
                                                nodata=nodata, v_range=v_range,all_touched=all_touched,
                                                tile_min_overlap=tile_min_overlap,scale=scale)
            volume_res_list.append(out_volume)
            if pixel_count > 0:
                in_range_per.append(100*(pixel_count-not_in_range_c)/pixel_count)
            else:
                in_range_per.append(0)
            dem_diff_min.append(min_v)
            dem_diff_max.append(max_v)
            carbon_move_list.append(carbon_move)

    elif process_num > 1:
        threadpool = Pool(process_num)
        para_list = [(idx, polygon, res_x, image_tiles, img_tile_polygons, carbon_tiles, carbon_tile_polygons, nodata, v_range, all_touched,
                      tile_min_overlap,scale)
                     for idx, polygon in enumerate(polygons)]
        results = threadpool.starmap(cal_volume_one_polygon, para_list)
        # save results
        for res in results:
            out_volume, pixel_count, not_in_range_c, min_v, max_v, carbon_move = res
            volume_res_list.append(out_volume)
            if pixel_count > 0:
                in_range_per.append(100 * (pixel_count - not_in_range_c) / pixel_count)
            else:
                in_range_per.append(0)
            dem_diff_min.append(min_v)
            dem_diff_max.append(max_v)
            carbon_move_list.append(carbon_move)
        threadpool.close()
        threadpool.join()
    else:
        raise ValueError('Wrong process number: %s ' % str(process_num))


    # if os.path.isfile(save_path) is False:
    #     io_function.copy_file_to_dst(boundary_shp, save_path)

    vol_att = {'vol_chg': volume_res_list, 'INrangePer': in_range_per, 'D_diff_min':dem_diff_min, 'D_diff_max':dem_diff_max, 'carbon_mo': carbon_move_list}
    save_format = vector_gpd.guess_file_format_extension(save_path)
    if save_path != boundary_shp:
        vector_gpd.add_attributes_to_shp(boundary_shp, vol_att, save_as=save_path, format=save_format)
    else:
        vector_gpd.add_attributes_to_shp(boundary_shp,vol_att,format=save_format)


def test_cal_volumetric_change_from_dem_diff():
    print(datetime.now(),'Testing:', 'test_cal_volumetric_change_from_dem_diff')
    data_dir = os.path.expanduser('~/Data/rts_ArcticDEM_mapping/cal_volume_from_demDiff/')
    # poly_shp = os.path.join(data_dir,'dem_diff_sam_results_select','dem_diff_sam_selected_by_quickshift_ARTS_v6.gpkg')
    poly_shp = os.path.join(data_dir,'input_test.shp')
    dem_diff = os.path.join(data_dir,'valid_json_files_20251222_TP_dem_diff_2023_mosaic.vrt')

    # nodata = 0
    nodata = 32767

    save_path = os.path.join(data_dir,'testing_output.shp')
    cal_volumetric_change_from_dem_diff(poly_shp,dem_diff,save_path=save_path, v_range=[None,-2],nodata=nodata, process_num=8)


def main(options, args):

    # test_cal_volumetric_change_from_dem_diff()
    # sys.exit(0)

    boundary_shp = args[0]
    dem_diff_path = args[1]
    save_path = boundary_shp if options.save_path is None else options.save_path
    process_num = options.process_num
    # nodata on DEM difference
    bit16_nodata = 32767
    carbon_path = options.carbon_raster

    #TODO: add a option to set the range, for now, we use the default value [-2.0, None] to calculate the volume change



    cal_volumetric_change_from_dem_diff(boundary_shp, dem_diff_path, carbon_raster=carbon_path, save_path=save_path, v_range=[None,-2.0],
                                        tile_min_overlap=None, nodata=bit16_nodata, all_touched=True,
                                        process_num=process_num)


    basic.outputlogMessage("Processing completed.")


if __name__ == "__main__":
    usage = "usage: %prog [options] boundary_shp DEM_diff_path_or_dir "
    parser = OptionParser(usage=usage, version="1.0 2026-01-04")
    parser.description = 'Introduction: calculate volume change   '

    parser.add_option("-s", "--save_path",
                      action="store", dest="save_path",
                      help="The path to save the output as a new file")

    parser.add_option("-c", "--carbon_raster",
                      action="store", dest="carbon_raster",
                      help="The path to the carbon raster file, if set, it will be used to  calculate the carbon mobilization based on the volume change," \
                      "if change carbon to soil, ice or others, the results should be related to soil or ice or others")

    
    parser.add_option("-p", "--process_num",
                      action="store", dest="process_num", type=int, default=8,
                      help="number of processes to create the mosaic")




    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
