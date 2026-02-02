#!/usr/bin/env python
# Filename: download_arcticDEM_pDEMtools.py
"""
introduction: download ArcticDEM and apply the co-registration using pDEMtools (https://pdemtools.readthedocs.io)

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 01 February, 2026
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.io_function as io_function
import basic_src.basic as basic

from dateutil import parser as datetime_parser # chnange the name to avoid conflict with OptionParser parser

from dem_common import arcticDEM_reg_tif_dir,arcticDEM_tile_reg_tif_dir

from download_arcticDEM import save_id_grid_no_dem

machine_name = os.uname()[1]
# # the maximum number of processes for downloading in parallel
# max_task_count = 1
# download_tasks = []


def download_dem_within_bounds(bounds,outdir,ext_id, dataset='arcticdem',date_start='2008-01-01', date_end='2026-12-31',
                               baseline_max_hours = 24, min_aoi_frac = 0.1,search_save='tmp.gpkg',out_res=None, b_unique_grid=False,
                               mask_coreg=False,custom_mask_fpath = None, geoid_correct = False):
    # this function is modified from: https://github.com/trchudley/pdemtools/blob/bb1e9d29721caf460ba56c1ea9b6ba9e388e5bb4/batch/batch_download_and_coregister_is2.py

    # bounds = (xmin, ymin, xmax, ymax)
    # dataset = "arcticdem" or  "rema"
    # mask_coreg: option to coregister only against stable ground, defaults to False

    # custom_mask_fpath:
    # filepath of custom mask of stable ground, if operating outside of Antarctica/Greenland
    # # This should be provided as a shapefile or equivalent geopandas-readable vector.
    # # Default is `None`, which will extract from BedMachine Greenland/Antarctica instead.

    # Filepath to BedMachine v5 netcdf file, if operating in Greenland/Antarctica. This
    # will allow for automated extraction of a stable ground mask for coregistration,
    # as well as a geoid for geoid correction. Leave blank if operating outside of
    # Greenland/Antarctica and/or providing your own mask and geoid data.
    #  - Download for Greenland: https://nsidc.org/data/idbmg4/versions/5
    #  - Download for Antarctica: https://nsidc.org/data/nsidc-0756/versions/3
    # bedmachine_fpath = ".../data/bedmachine_5/BedMachineGreenland-v5.nc"
    bedmachine_fpath = None

    # geoid_correct: option to geoid-correct data using the bedmachine geoid, defaults to False

    # filepath to custom geoid geotiff. Default is `None`, which will extract from BedMachine
    # Greenland/Antarctica instead.
    custom_geoid_fpath = None

    # this script will check
    override_datecheck = True

    date_start = date_start.replace('-','')
    date_end = date_end.replace('-','')

    dates = f"{date_start}/{date_end}"

    # sanity check date string:
    first_date_str = dates.split("/")[0]
    first_date = to_datetime(first_date_str, format="%Y%m%d")
    cutoff_date = to_datetime("2018-10-04")
    if first_date < cutoff_date:
        if override_datecheck:
            basic.outputlogMessage("Warning, "
                f"You have searched for data beginning {first_date_str}. No ICESat-2 data "
                "exists prior to 2018-10-04. As `override_datecheck` is set to True, "
                "any data prior to this date will be downloaded (uncorrected) anyway."
            )
        else:
            raise ValueError(
                f"You have searched for data beginning {first_date_str}. No ICESat-2 data "
                "exists prior to 2018-10-04. To override this check and download the (uncorrected) "
                "data anyway, set `override_datecheck = True`"
            )

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    done_indicator = outdir + f'_download.done'
    if os.path.isfile(done_indicator):
        print(f'downloading data for extent: {ext_id} has completed')
        return True

    # Search for DEM strips
    basic.outputlogMessage("\nSearching for DEM strips...")
    gdf = pdt.search(
        dataset=dataset,
        bounds=bounds,
        dates=dates,
        # months=[6, 7, 8, 9],
        # years=[2019],
        baseline_max_hours=baseline_max_hours,
        # sensors=["WV03", "WV02", "WV01"],
        # accuracy=2,
        min_aoi_frac=min_aoi_frac,
    )
    search_str = f'dataset:{dataset}, bounds: {bounds}, dates:{dates},baseline_max_hours:{baseline_max_hours}, min_aoi_frac:{min_aoi_frac}'
    if gdf is None:
        basic.outputlogMessage(f'Warning, NO search results for: {search_str}' )
        return

    n_strips = len(gdf)
    basic.outputlogMessage(f"{n_strips} strips found, with query conditions: {search_str}")

    gdf.to_file(search_save)
    basic.outputlogMessage(f'saved search results to {search_save}')

    if len(gdf) > 300:
        basic.outputlogMessage(f'error, too many ({len(gdf)}) DEMs within {ext_id} th extent with: {dataset} from {date_start} to {date_end}')
        raise ValueError(f'error, too many ({len(gdf)}) DEMs within {ext_id} th extent with: {dataset} from {date_start} to {date_end}')

    i = 1

    basic.outputlogMessage("Downloading DEM strips...")

    # because the script need to run the first loop, so many be not good to run in parallel
    first_loop = True

    for _, row in gdf.iterrows():

        # check if enough disk space available
        io_function.wait_until_enough_disk_space(outdir, min_disk_GB=5)

        date = row.pdt_datetime1.date()
        date_str = date.strftime("%Y%m%d")
        dem_id = row.pdt_id

        if b_unique_grid:
            out_fname = os.path.join(outdir, f'{dem_id}_grid{ext_id}')
        else:
            out_fname = os.path.join(outdir, f'{dem_id}_sub{ext_id}')

        # out_fname = os.path.join(outdir, f"{date_str}_{dem_id}")

        # If the file doesn't yet exist, download.
        if len(glob(f"{out_fname}*")) == 0:

            basic.outputlogMessage(f"\nDownloading {i}/{n_strips} {os.path.basename(out_fname)}...")

            # download DEM
            dem = pdt.load.from_search(row, bounds=bounds, bitmask=True)
            dem.compute()  # rioxarray uses lazy evaluation, so we can force the download using the `.compute()` function.

            # pad to full size of AOI
            dem = dem.rio.pad_box(*bounds, constant_values=np.nan)

            # If in the first loop, download the optional bedrock mask and geoid if requested
            if first_loop:

                # get bedrock mask
                if mask_coreg:
                    if custom_mask_fpath is None:
                        bedrock_mask = pdt.data.bedrock_mask_from_bedmachine(
                            bedmachine_fpath, dem
                        )
                    else:
                        bedrock_mask = pdt.data.bedrock_mask_from_vector(
                            custom_mask_fpath, dem
                        )
                else:
                    bedrock_mask = None

                if geoid_correct:
                    # Download and save geoid
                    if custom_geoid_fpath is None:
                        geoid = pdt.data.geoid_from_bedmachine(bedmachine_fpath, dem)
                    else:
                        geoid = pdt.data.geoid_from_raster(custom_geoid_fpath, dem)
                else:
                    geoid = None

                first_loop = False

            if to_datetime(date) < cutoff_date:
                metadata = {"coreg_status":"failed",
                            "remark":f'the date is {date}. No ICESat-2 data exists prior to 2018-10-04'
                            }
                basic.outputlogMessage(f"Warning, For {dem_id}, the date is {date}. No ICESat-2 data exists prior to 2018-10-04.")
            else:
                # get the IS2 coregistration data
                basic.outputlogMessage("Getting IS2 coregistration data...")
                is2_gdf = pdt.data.icesat2_atl06(dem, date, stable_mask=bedrock_mask)
                if len(is2_gdf) >= 2:
                    # coregister DEM, with return_stats=True.
                    dem, metadata = dem.pdt.coregister_is2(
                        is2_gdf,
                        stable_mask=bedrock_mask,
                        max_horiz_offset=50,
                        return_stats=True,
                    )
                else:
                    metadata = {"coreg_status": "failed",
                                "remark": f'Not enough ICESat-2 data (only found {len(is2_gdf)} points) '
                                }

            # geoid correct if necessary
            if geoid_correct:
                basic.outputlogMessage("Geoid-correcting...")
                dem = dem.pdt.geoid_correct(geoid)

            # check whether coreg worked, and construct filename appropriately
            if metadata["coreg_status"] == "failed":
                out_fpath = out_fname + ".tif"
            elif metadata["coreg_status"] == "coregistered":
                out_fpath = out_fname + "_coreg.tif"
            elif metadata["coreg_status"] == "dz_only":
                out_fpath = out_fname + "_coreg_dz.tif"
            else:
                basic.outputlogMessage(f"Unknown coregistration status {metadata['coreg_status']}. Saving wihtout extension.")
                out_fpath = out_fname + ".tif"

            # Export to geotiff
            dem.rio.to_raster(out_fpath, compress="ZSTD", predictor=3, zlevel=1)
            del dem


            # Export metadata as json
            json_fpath = out_fpath.rsplit(".", 1)[0] + ".json"
            with open(json_fpath, "w") as json_file:
                json.dump(metadata, json_file, indent=4)

            # quick check the resolution
            with rasterio.open(out_fpath) as src:
                dem_res_x, dem_res_y = src.res
            if abs(out_res - dem_res_x) > 0.000001:
                basic.outputlogMessage(f'Error, the resolution {dem_res_x} in file: {out_fpath} is not consistent with the settting: {out_res}')
                raise ValueError(f'Error, the resolution {dem_res_x} in file: {out_fpath} is not consistent with the settting: {out_res}')
        else:
            basic.outputlogMessage(f"Files exists: {out_fname}*, skip")
        i += 1

    basic.outputlogMessage(f"Finished downloading and co-registration for extent: {ext_id}")

    saved_file_list = io_function.get_file_list_by_pattern(outdir,'*.tif')
    # checking file count
    if n_strips != len(saved_file_list):
        raise ValueError(f'downloading for polygon {ext_id} is not completed')

    io_function.save_text_to_file(done_indicator, f'completed downloading {len(saved_file_list)} files')



def download_dem_pDENtools(extent_bounds_list, output_dir, date_start,date_end,pre_name,
                      poly_ids=None,save_prj_code=3413,b_unique_grid=False, out_res=None,
                           baseline_max_hours=24, min_overlap_per=0.05):

    # download data through STAC, not need to unpack
    b_save_grid_id_noDEM = True
    if poly_ids is None:
        poly_ids = [idx for idx in range(len(extent_bounds_list)) ]
        b_save_grid_id_noDEM = False    # if poly_ids is not the global unique id, then don't save it.

    for count, (idx, ext_bound) in enumerate(zip(poly_ids, extent_bounds_list)):
        if b_unique_grid:
            search_save = os.path.join(output_dir, pre_name + f'_grid{idx}.gpkg')
            ext_tif_save_dir = os.path.join(output_dir, f'dem_grid{idx}')
        else:
            search_save = os.path.join(output_dir, pre_name + f'_poly_{idx}.gpkg')
            ext_tif_save_dir = os.path.join(output_dir, f'dem_poly_{idx}')

        res = download_dem_within_bounds(ext_bound,ext_tif_save_dir, idx, date_start=date_start,date_end=date_end,
                                     baseline_max_hours = baseline_max_hours, min_aoi_frac = min_overlap_per, search_save=search_save,out_res=out_res,
                                         b_unique_grid=b_unique_grid)

        # if res is False, mean no data found
        if res is False and b_save_grid_id_noDEM is True:
            save_id_grid_no_dem(idx)


def main(options, args):

    basic.setlogfile(f'pDEMtools_log_pID{os.getpid()}.txt')

    extent_shp = args[0]
    basic.outputlogMessage(f'input extent: {extent_shp}')
    b_arcticDEM_mosaic = options.b_arcticDEM_mosaic # default is False

    # global max_task_count
    # max_task_count = options.max_process_num
    min_overlap_per = options.min_overlap_per
    baseline_max_hours = options.baseline_max_hours

    # reg_tif_dir is the folder to save tifs
    if b_arcticDEM_mosaic:
        reg_tif_dir = arcticDEM_tile_reg_tif_dir
    else:
        reg_tif_dir = arcticDEM_reg_tif_dir

    # use the user specific save_dir if it's set
    if options.save_dir is not None:
        reg_tif_dir = options.save_dir
    if os.path.isdir(reg_tif_dir) is False:
        io_function.mkdir(reg_tif_dir)
    reg_tif_dir = os.path.abspath(reg_tif_dir)  # change to absolute path

    prj_code = vector_gpd.get_projection(extent_shp)
    if prj_code != 3413:
        basic.outputlogMessage(f'Map projection of {extent_shp} is not EPSG:3413, will convert to it')
        ext_polys_3413 = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp, 'EPSG:3413')
    else:
        ext_polys_3413 = vector_gpd.read_polygons_gpd(extent_shp)

    # pDEMtools need bounds of polygons? for grid cells (rectangles), they should be the same
    ext_bounds_3413 = [vector_gpd.get_polygon_bounding_box(item) for item in ext_polys_3413]

    pre_name = os.path.splitext(os.path.basename(extent_shp))[0]
    pre_name += '_Mosaic' if b_arcticDEM_mosaic else '_Strip'

    # read 'grid_id' if the extent shp is from grid shp file, if not, grid_id_list will be None
    if vector_gpd.is_field_name_in_shp(extent_shp,'grid_id'):
        grid_id_list = vector_gpd.read_attribute_values_list(extent_shp,'grid_id')
    elif vector_gpd.is_field_name_in_shp(extent_shp,'RowCol_id'):
        grid_id_list = vector_gpd.read_attribute_values_list(extent_shp, 'RowCol_id')
    else:
        grid_id_list = None
    b_unique_grid = False if grid_id_list is None else True
    if len(ext_bounds_3413) < 1:
        raise ValueError('No polygons in %s'%extent_shp)
    else:
        basic.outputlogMessage('read %d extent polygons in %s for downloading using pDEMtools'%(len(ext_bounds_3413),extent_shp))

    date_start = options.date_start
    date_end = options.date_end

    # projection for ArcticDEM products
    save_prj_code = 3413
    out_resolution = options.out_res


    download_dem_pDENtools(ext_bounds_3413, reg_tif_dir, date_start, date_end, pre_name,
                      poly_ids=grid_id_list, save_prj_code=save_prj_code,b_unique_grid=b_unique_grid,out_res=out_resolution,
                           baseline_max_hours=baseline_max_hours, min_overlap_per=min_overlap_per)




if __name__ == '__main__':
    # import need library
    import pdemtools as pdt
    from warnings import warn
    from glob import glob
    from math import isnan
    from pandas import to_datetime

    import numpy as np
    import json
    import rasterio

    # import geopandas as gpd
    # import xarray as xr
    # import rioxarray


    # test_select_search_results_each_month()
    # sys.exit(0)

    usage = "usage: %prog [options] extent_shp "
    parser = OptionParser(usage=usage, version="1.0 2026-02-01")
    parser.description = 'Introduction: download ArcticDEM within an extent  '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",
                      help="the folder to save DEMs")

    parser.add_option("-s", "--date_start",
                      action="store", dest="date_start", type=str, default='2008-01-01',
                      help="the start date for downloading data")

    parser.add_option("-e", "--date_end",
                      action="store", dest="date_end", type=str, default='2026-12-31',
                      help="the end date for downloading data")

    parser.add_option("-r", "--out_res",
                      action="store", dest="out_res", type=float, default=2.0,
                      help="the spatial resolution for output (2.0)")

    parser.add_option("-o", "--min_overlap_per",
                      action="store", dest="min_overlap_per", type=float, default=0.05,
                      help="the threshold of minimum overlap between bounds and DEM files")

    parser.add_option("-b", "--baseline_max_hours",
                      action="store", dest="baseline_max_hours", type=float, default=720,
                      help="the maximum time baseline for query, ideally, one day, but low this to 30 days for more data")


    # parser.add_option("-p", "--max_process_num",
    #                   action="store", dest="max_process_num", type=int, default=4,
    #                   help="the maximum number of processes for downloading in parallel")

    parser.add_option("-m", "--b_arcticDEM_mosaic",
                      action="store_true", dest="b_arcticDEM_mosaic",default=False,
                      help="if set, will download the mosaic of ArcticDEM, other wise, download the strips")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
