

import os, sys
import ArcticDEM_proc_grid

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import raster_io

# def test_dem_diff_newest_oldest():
#     dir=os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff')
#     os.chdir(dir)
#
#     extent_id = 1
#     mosaic_yeardate_dir = os.path.join('./', 'dem_date_mosaic_sub_%d' % extent_id)
#     with open(os.path.join(mosaic_yeardate_dir, 'dem_valid_percent.txt')) as f_job:
#         tif_names = [line.split()[0] for line in f_job.readlines()]
#         dem_tif_list = [os.path.join(mosaic_yeardate_dir, item) for item in tif_names]
#
#     save_dem_diff = 'test' + '_DEM_diff_sub_%d.tif' % extent_id
#     save_date_diff = 'test' + '_date_diff_sub_%d.tif' % extent_id
#     ArcticDEM_proc_grid.dem_diff_newest_oldest(dem_tif_list, save_dem_diff, save_date_diff)

def test_check_dem_valid_per():
    dir=os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff')
    os.chdir(dir)

    tif = 'dem_stripID_mosaic_sub_1/20170211_10400100285F1700_1040010028337F00.tif'
    per = raster_io.get_valid_pixel_percentage(tif)
    print(tif, 'per',per)