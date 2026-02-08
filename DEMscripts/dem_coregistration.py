#!/usr/bin/env python
# Filename: dem_coregistration 
"""
introduction: dem co-registration

using demcoreg, install it and its dependencies:

    git clone https://github.com/dshean/imview.git
    pip install -e imview

   git clone https://github.com/dshean/pygeotools.git
   pip install -e pygeotools

   git clone https://github.com/dshean/demcoreg.git
   pip install -e demcoreg

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 01 March, 2021
"""

import os, sys
from optparse import OptionParser
deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
import raster_io
import vector_gpd

from multiprocessing import Process

dem_dem_align = os.path.expanduser('~/codes/github_public_repositories/demcoreg/demcoreg/dem_align.py')

# use psutil to check available memory.
import psutil
import time

# check the resource used by python
# import resource
# resource.getrusage(resource.RUSAGE_SELF)

from dem_mosaic_crop import check_dem_valid_per

def choose_reference_dem(dem_list, dem_valid_per_txt):
    if dem_valid_per_txt is None:
        raise ValueError('NO information of valid percentage of DEMs, cannot choose a reference DEM')

    with open(dem_valid_per_txt, 'r') as f_obj:
        tif_valid_per_list = [line.strip().split() for line in f_obj.readlines()]
        tif_higest_per = None
        per_max = 0
        for tif, per in tif_valid_per_list:
            if float(per) > per_max:
                per_max = float(per)
                tif_higest_per = tif

        for dem_tif in dem_list:
            if tif_higest_per in dem_tif:
                return dem_tif

    return None

def check_coreg_results(dem_tif, save_dir):
    dem_align = os.path.join(save_dir, 'dem_coreg', os.path.basename(io_function.get_name_by_adding_tail(dem_tif, 'coreg')))
    if os.path.isfile(dem_align):
        return True
    return False

def check_align_folder(dem_tif):
    # by default, dem_align.py save the results to where dem_tif is
    res_dir = os.path.dirname(dem_tif)
    align_folder = os.path.splitext(os.path.basename(dem_tif))[0] + '_dem_align'
    align_dir = os.path.join(res_dir,align_folder)
    # after dem_align.py usually have 9 files
    align_outputs = io_function.get_file_list_by_pattern(align_dir,'*')
    # print(align_outputs)
    return align_outputs


def move_align_results(ref_dem, dem_tif, save_dir):

    coreg_save_dir = os.path.join(save_dir, 'dem_coreg')
    if os.path.isdir(coreg_save_dir) is False:
        io_function.mkdir(coreg_save_dir)

    align_outputs = check_align_folder(dem_tif)
    if len(align_outputs) < 9:
        raise ValueError('the output of dem_align.py is less than 9 files')

    dem_align = os.path.join(coreg_save_dir, os.path.basename(io_function.get_name_by_adding_tail(dem_tif, 'coreg')))
    # align DEM and a filt version, which one should I use? what filter they apply?
    # visually check one results (Banks east) , a the same location, align DEM and a filt one have exact values,
    # but the filt version have more nodata.  Let's use the filt version.
    # the nodata pixels usually are water pixels, but also some inside the thaw slumps
    align_filt = [out for out in align_outputs if out.endswith('align_filt.tif')][0]
    io_function.move_file_to_dst(align_filt,dem_align, overwrite=True)

    # copy reference dem if necessary
    ref_dem_copy = os.path.join(coreg_save_dir, os.path.basename(ref_dem))
    if os.path.isfile(ref_dem_copy) is False:
        io_function.copy_file_to_dst(ref_dem,ref_dem_copy)

    # move the elevation difference?
    ele_diff_folder = os.path.join(save_dir,'dem_diff_from_demcoreg')
    if os.path.isdir(ele_diff_folder) is False:
        io_function.mkdir(ele_diff_folder)
    dem_diff_filt = [out for out in align_outputs if out.endswith('align_diff_filt.tif')][0]
    io_function.movefiletodir(dem_diff_filt,ele_diff_folder, overwrite=True)

    coreg_png_plot_folder = os.path.join(save_dir,'demcoreg_png_plot')
    if os.path.isdir(coreg_png_plot_folder):
        io_function.mkdir(coreg_png_plot_folder)
    coreg_pngs = [out for out in align_outputs if out.endswith('.png')]
    for png in coreg_pngs:
        io_function.movefiletodir(png, coreg_png_plot_folder, overwrite=True)

    return True


# def co_registration_parallel(ref_dem, dem_list, save_dir, process_num):
#     print('ref_dem', ref_dem)
#     print('source dem:')
#     for dem_tif in dem_list:
#         print(dem_tif)
#
#     # dem_align.py requires large memory, for example, a region of 50000 by 50000 pixels, may requires more than 110 GB memory.
#
#     # parallel --progress --delay 10 -j 14 dem_align.py ${ref} {} ::: $(ls *_dem.tif | grep -v 2012)
#     commond_str = 'parallel --progress --delay 5 -j %d %s %s'%(process_num, dem_dem_align ,ref_dem)
#     commond_str += ' {} ::: ' + ' '.join(dem_list)
#     print(commond_str)
#     os.system(commond_str)

    # move results to another folder

def co_registration_one_dem(ref_dem, dem_tif, save_dir):
    if check_coreg_results(dem_tif,save_dir):
        basic.outputlogMessage('co-registeration results for %s exists, skip'%dem_tif)
        return 0

    align_outputs = check_align_folder(dem_tif)
    if len(align_outputs) >= 9:
        basic.outputlogMessage('%s has been co-registered, skip'%dem_tif)
    else:
        commond_str = dem_dem_align + ' ' + ref_dem + ' ' + dem_tif
        basic.outputlogMessage(commond_str)
        res = os.system(commond_str)
        if res != 0:
            sys.exit(1)

    return move_align_results(ref_dem, dem_tif, save_dir)


def co_registration_multi_process(ref_dem, dem_list, save_dir, process_num):
    print('ref_dem', ref_dem)
    print('source dem:')
    for dem_tif in dem_list:
        print(dem_tif)

    sub_tasks = []
    for dem_tif in dem_list:
        height, width, band_num, daet_type = raster_io.get_height_width_bandnum_dtype(dem_tif)
        # estimate memory need
        need_memory = height*width*4*12 # usually, each pixel need 4 Bytes (float), dem_align.py need more than memory 12 times of file size
        avai_memory = psutil.virtual_memory().available

        while need_memory > avai_memory:
            print('waiting more available memory, need: %.4f GB, available: %.4f GB'%(need_memory/(1000*1000*1000), avai_memory/(1000*1000*1000)))
            time.sleep(10)
            avai_memory = psutil.virtual_memory().available

        while basic.alive_process_count(sub_tasks) >= process_num:
            time.sleep(10)

        sub_process = Process(target=co_registration_one_dem, args=(ref_dem, dem_tif,save_dir))
        sub_process.start()
        sub_tasks.append(sub_process)
        time.sleep(1)
        if sub_process.exitcode is None:    # even the function has return, sub_process.alive still true for now
            time.sleep(10)   # wait 10 seconds
        if sub_process.exitcode is not None and sub_process.exitcode != 0:
            sys.exit(1)

def calculate_dem_valid_per(dem_list, dem_dir,process_num):
    # get the area_pixel_count by using the maximum size
    # it will save "dem_valid_percent.txt" under dem_dir
    area_pixel_count = 0
    for dem_file in dem_list:
        height, width, _, _ = raster_io.get_height_width_bandnum_dtype(dem_file)
        if height * width > area_pixel_count:
            area_pixel_count = height * width
    basic.outputlogMessage('Area pixel count: %d' % area_pixel_count)

    # do not remove dem files (move_dem_threshold=None)
    dem_tif_list = check_dem_valid_per(dem_list, dem_dir, process_num=process_num,
                                       move_dem_threshold=None, area_pixel_num=area_pixel_count)
    return dem_tif_list

def download_ICESat2(input_date,days_r,params,epsg,min_is2_points):

    first_iteration = True
    gdf = None
    for r in days_r:
        if first_iteration:
            first_iteration = False
        else:
            print(f"Expanding temporal search.")
        date_min = to_datetime(input_date) - datetime.timedelta(days=r)
        date_max = to_datetime(input_date) + datetime.timedelta(days=r)
        date_min = date_min.strftime("%Y-%m-%d")
        date_max = date_max.strftime("%Y-%m-%d")

        params["t0"] = date_min
        params["t1"] = date_max

        print(f"Querying points within {r} days... ", end="")

        gdf = icesat2.atl06p(params).to_crs(epsg)

        n_points = len(gdf)

        # Print statement
        if n_points == 0:
            print(f"No points found intersecting.", end=" ")
        elif n_points < min_is2_points:
            print(f"{n_points} points found. Below minimum threshold ({min_is2_points}).")
        else:
            print(f"{n_points} points found.")
            break

    if len(gdf) < min_is2_points:
        warn("Returning GeoDataFrame with number of points below minimum threshold.")

    return gdf

def download_ICESat2_in_bounds(bounds, epsg, save_path=None,max_point_count_year=None):
    # try to download all icesat 2 data from 2018, summer data.
    # max_point_count_year each year

    # bounds: bounds in format [xmin, ymin, xmax, ymax], if the
    #         `epsg` parameter is provided.

    # param epsg: EPSG code for the target_rxd dataset. If `target_rxd` is a tuple,
    #         this must be provided.
    #     :type epsg: int

    if save_path is not None and os.path.isfile(save_path):
        basic.outputlogMessage(f'{save_path} already exists, read it directly')
        is2_gdf = gpd.read_file(save_path)
        return is2_gdf

    # connect to sliderule
    # Get search region as shapely geometry in epsg:4326
    gdf_4326 = gpd.GeoDataFrame(geometry=[box(*bounds)], crs=epsg).to_crs(4326)

    icesat2.init("slideruleearth.io")
    sr_region = sliderule.toregion(gdf_4326)

    # set search parameters (apart from dates)
    params = {
        "poly": sr_region["poly"],
        "srt": 3,  # Surface. 0-land, 1-ocean, 2-seaice, 3-landice (default), 4-inlandwater
        "cnf": 1,  # Confidence. Default 1 (within 10 m). 2: Low. 3: Medium. 4: High.
        "ats": 10,  # Mininum along track spread. SR Default: 20. pDEMtoold default: 10
        "cnt": 10,  # Minimum photon count in segment. Default 10.
        "len": 40,  # Extent length. ATL06 default is 40 metres.
        "res": 20,  # Step distance. ATL06 default is 20 metres.
        "track": 0,  # Integer: 0: all tracks, 1: gt1, 2: gt2, 3: gt3. Default 0.
        "sigma_r_max": 5,  # Max robust dispersion [m]. Default 5.
    }


    is2_gdf_list = []
    for year in range(2018,  date.today().year):
        mid_date = f'{year}-07-01'
        print('\n\nmid-date:',mid_date)
        # Sanity check date
        dem_date = to_datetime(mid_date)
        cutoff_date = to_datetime("2018-10-04")
        if dem_date < cutoff_date:
            warn(
                f"You have searched for data beginning {dem_date}. No ICESat-2 data "
                "exists prior to 2018-10-04. This function will return no usable data."
            )
        days_r = [7, 14, 28, 60]  # maximum, May to September
        min_is2_points = 2000
        is2_gdf_year = download_ICESat2(mid_date, days_r, params, epsg, min_is2_points)

        # Only append non-empty results
        if is2_gdf_year is not None and len(is2_gdf_year) > 0:
            is2_gdf_year['year'] = year # add the year
            is2_gdf_list.append(is2_gdf_year)

    # combine all yearly GeoDataFrames
    if len(is2_gdf_list) > 0:
        is2_gdf = gpd.GeoDataFrame(pd.concat(is2_gdf_list, ignore_index=True), crs=is2_gdf_list[0].crs)
    else:
        return None

    # save to file if save_path is provided
    if save_path is not None:
        # Extension-based save (e.g., .gpkg, .shp, .geojson)
        is2_gdf.to_file(save_path)

    return is2_gdf

def co_registration_one_dem_is2(dem_file, is2_gdf, save_path, stable_mask=None):
    if os.path.isfile(save_path):
        print(f'Warning, co-registration results already exist: {save_path}, skip')
        return

    # read dem
    dem = rioxarray.open_rasterio(dem_file).squeeze()  # needed rioxarray for pDEMtools

    gdf_flt = is2_gdf.copy()

    # mask points in NoData regions
    coords_ds = Dataset(
        {
            "x": (["points"], is2_gdf.geometry.x.values),
            "y": (["points"], is2_gdf.geometry.y.values),
        }
    )
    # Determine raster value and mask at each point
    # problem: a few points in no-data region, but close to edge, got valid values? Just an issue in QGIS visualization?
    gdf_flt["values"] = dem.interp(coords_ds, method="nearest").values


    # problem: a few points in no-data region, but close to edge, got valid values
    # a few points got -9999 values
    # with rasterio.open(dem_file) as src:
    #     points = list(zip(is2_gdf.geometry.x.values, is2_gdf.geometry.y.values))
    #     values = [v[0] for v in src.sample(points)]
    # gdf_flt["values"] = values

    # For DataArray loaded via rioxarray (still:  a few points in no-data region)
    # transform = dem.rio.transform()
    # rows, cols = rasterio.transform.rowcol(transform, is2_gdf.geometry.x.values, is2_gdf.geometry.y.values)
    # mask = ((rows >= 0) & (rows < dem.shape[-2]) & (cols >= 0) & (cols < dem.shape[-1]))
    # values = np.full(len(rows), np.nan, dtype=float)
    # values[mask] = dem.data[rows[mask], cols[mask]]
    # gdf_flt["values"] = values

    # Filter values at each point
    gdf_flt = gdf_flt[~gdf_flt["values"].isnull()]

    # save point for checking:
    gdf_ft_save = io_function.get_name_by_adding_tail(save_path,'is2').replace('.tif','.gpkg')
    gdf_flt.to_file(gdf_ft_save)

    # quick testing, read the gpkg after removing the points in the nodata regions
    # after quick testing, the rms still very large (>50)
    # gdf_flt = gpd.read_file(gdf_ft_save)


    basic.outputlogMessage(f'After filtering Nodata, select {len(gdf_flt)} points from {len(is2_gdf)} for co-registration')

    # use gdf_flt for co-registration
    if len(gdf_flt) >= 2:
            # coregister DEM, with return_stats=True.
            dem, metadata = dem.pdt.coregister_is2(
                gdf_flt,
                stable_mask=stable_mask,
                max_horiz_offset=50,
                return_stats=True,
            )
    else:
        metadata = {"coreg_status": "failed",
                    "remark": f'Not enough ICESat-2 data (only found {len(is2_gdf)} points) '
                    }

    # Export to geotiff
    # out_fpath = io_function.get_name_by_adding_tail(out_fpath, 'dem')  # add dem in the end, for existing working flow
    dem.rio.to_raster(save_path, compress="ZSTD", predictor=3, zlevel=1)
    del dem

    # Export metadata as json
    json_fpath = save_path.rsplit(".", 1)[0] + ".json"
    with open(json_fpath, "w") as f_obj:
        json.dump(metadata, f_obj, indent=4)




def co_registration_icesat2_pDEMtools(dem_list, save_dir,extent_shp=None, process_num=8, bounds=None,epsg=None):
    # co_registration by using ICESat-2 datasets

    # bounds: bounds in format [xmin, ymin, xmax, ymax], if the
    #         `epsg` parameter is provided.

    # param epsg: EPSG code for the target_rxd dataset. If `target_rxd` is a tuple,
    #         this must be provided.
    #     :type epsg: int

    # download all ICESat-2 within the bounds
    if extent_shp is not None:
        polys = vector_gpd.read_polygons_gpd(extent_shp)
        if len(polys)!=1:
            raise ValueError('Currently, only support one extent')
        bounds = vector_gpd.get_polygon_bounding_box(polys[0])
        # print('bounds:', bounds)
        epsg = vector_gpd.get_projection(extent_shp,format='epsg')

    else:
        # assume they have the same bounds for dem_list
        if bounds is None:
            bounds_box = raster_io.get_image_bound_box(dem_list[0])
            # print('bounds_box:',bounds_box)
            bounds = tuple(bounds_box)
            # bounds = ( bounds[1],bounds[0], bounds[2], bounds[1] )
            # print('bounds:',bounds)
        if epsg is None:
            epsg = raster_io.get_projection(dem_list[0],format='epsg')
            # print('epsg:', epsg)

        # checking bounds and epsg consistent
        for dem_file in dem_list:
            tmp_bounds = raster_io.get_image_bound_box(dem_file)
            tmp_epsg = raster_io.get_projection(dem_file, format='epsg')
            if tmp_epsg!=epsg:
                raise ValueError(f'Projection ({epsg} vs {tmp_epsg}) insistence between {dem_list[0]} and {dem_file}')
            if bounds != tmp_bounds:
                raise ValueError(f'Bounds ({bounds} vs {tmp_bounds}) insistence between {dem_list[0]} and {dem_file}')

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    # ICESat-2 data
    icesat2_save_path = os.path.join(save_dir,'icesat2_data.gpkg')
    # download ICESat-2
    is2_gdf = download_ICESat2_in_bounds(bounds, epsg, save_path=icesat2_save_path, max_point_count_year=None)

    # clear ICESat-2 data, only keep these that are on stable surface
    # No idea on how to clear it, because we don't see repeat observation in the download data.Feb 8, 2026.

    # conduct co-registration
    for idx, dem_file in enumerate(dem_list):
        print(f'co-registration: {idx+1}/{len(dem_list)}: {os.path.basename(dem_file)}')
        co_reg_save = os.path.join(save_dir,os.path.basename(dem_file))
        co_registration_one_dem_is2(dem_file, is2_gdf, co_reg_save)




def test_co_registration_icesat2_pDEMtools():
    data_dir = os.path.expanduser('~/Data/dem_processing/registration_tifs')
    dem_dir = os.path.join(data_dir,'dem_grid0016000342')
    ext_shp = os.path.expanduser('~/Data/dem_processing/grid_shp/valid_json_files_20251222_TP_grid_10km_subsets/sub299.shp')
    dem_list = io_function.get_file_list_by_ext('.tif',dem_dir,bsub_folder=False)
    print(f'Found {len(dem_list)} dem files')

    process_num = 8
    save_dir = dem_dir + '_coreg'
    co_registration_icesat2_pDEMtools(dem_list, save_dir, extent_shp=ext_shp, process_num=process_num, bounds=None, epsg=None)


def main(options, args):

    save_dir = options.save_dir
    dem_dir_or_txt = args[0]
    ref_dem = options.ref_dem
    dem_valid_per_txt = options.dem_valid_per_txt
    process_num = options.process_num

    if os.path.isfile(dem_dir_or_txt):
        dem_list = io_function.read_list_from_txt(dem_dir_or_txt)
    else:
        dem_list = io_function.get_file_list_by_ext('.tif', dem_dir_or_txt, bsub_folder=False)
        if dem_valid_per_txt is None:
            dem_valid_per_txt = os.path.join(dem_dir_or_txt,'dem_valid_percent.txt')

        # if the file not exist, calculate it
        if os.path.isfile(dem_valid_per_txt) is False and os.path.isdir(dem_dir_or_txt):
            # calculate the valid percent (don't remove files, move_dem_threshold=None)
            calculate_dem_valid_per(dem_list, dem_dir_or_txt, process_num)

    dem_count = len(dem_list)
    if dem_count < 1:
        raise ValueError('No input dem files in %s' % dem_dir_or_txt)

    if ref_dem is None:
        ref_dem = choose_reference_dem(dem_list, dem_valid_per_txt)
        if ref_dem is None:
            raise ValueError('Cannot find a reference DEM')

    if ref_dem in dem_list:
        dem_list.remove(ref_dem)
    # co_registration_parallel(ref_dem,dem_list,save_dir,process_num)
    co_registration_multi_process(ref_dem, dem_list, save_dir, process_num)


if __name__ == '__main__':

    import pdemtools as pdt
    import rioxarray
    from xarray import DataArray, Dataset
    from datetime import date
    import geopandas as gpd
    import pandas as pd
    import rasterio
    import numpy as np

    from sliderule import sliderule, icesat2
    from pandas import to_datetime
    from warnings import warn
    from shapely.geometry import box
    import datetime
    import json

    test_co_registration_icesat2_pDEMtools()
    sys.exit(0)

    usage = "usage: %prog [options] dem_tif_dir or dem_list_txt "
    parser = OptionParser(usage=usage, version="1.0 2020-3-1")
    parser.description = 'Introduction: co-registration for multi-temporal DEMs '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir", default='./',
                      help="the folder to save results")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes")

    parser.add_option("-r", "--ref_dem",
                      action="store", dest="ref_dem",
                      help="the reference DEM, if not set, it will select with the most coverage")

    parser.add_option("-p", "--dem_valid_per_txt",
                      action="store", dest="dem_valid_per_txt",
                      help="a txt file storing the valid percentage of all the DEM")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)