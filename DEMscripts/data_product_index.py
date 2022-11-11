#!/usr/bin/env python
# Filename: data_product_index.py 
"""
introduction: create index file for data products

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 10 November, 2022
"""
import os,sys

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import basic_src.basic as basic
import vector_gpd

from optparse import OptionParser
import pandas as pd
from datetime import datetime

org_grid_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km.shp')


def create_product_index_shp(product_dir,product_name,save_path,prj_wkt):
    print(datetime.now(),'creating index for %s'%product_dir)
    # read cells in org_grid_shp
    polygons, grid_ids = vector_gpd.read_polygons_attributes_list(org_grid_shp,'id',b_fix_invalid_polygon=False)
    count = len(polygons)
    save_polys = []
    save_ids = []   # cell_id,
    file_names = [] # cell_id
    fileurls = []  # put as relative path
    product_names = []

    for idx, (poly, id) in enumerate(zip(polygons,grid_ids)):
        print('checking grid %d / %d'%(idx+1, count))
        files = io_function.get_file_list_by_pattern(product_dir,'ext*/*grid%d.tar.gz'%id)
        if len(files) == 0:
            continue
        elif len(files) == 1:
            save_polys.append(poly)
            save_ids.append(id)
            file_names.append(os.path.basename(files[0]))
            fileurls.append(files[0])
            product_names.append(product_name)
        else:
            raise ValueError('Multiple files for grid %d in %s'%(id, product_dir))


    # save
    save_pd = pd.DataFrame({'cell_id':save_ids, 'tarball':file_names,'fileurl':fileurls, 'Polygon':save_polys})
    vector_gpd.save_polygons_to_files(save_pd,'Polygon',prj_wkt,save_path)
    basic.outputlogMessage('saved to %s'%os.path.abspath(save_path))



def main(options, args):

    product_dir = args[0]
    io_function.is_folder_exist(product_dir)
    product_name = os.path.basename(product_dir)
    save_path = product_name + '_Index.shp' if options.save_path is None else options.save_path
    prj_wkt = map_projection.get_raster_or_vector_srs_info_proj4(org_grid_shp)
    create_product_index_shp(product_dir, product_name,save_path,prj_wkt)


if __name__ == '__main__':
    usage = "usage: %prog [options] product_folder "
    parser = OptionParser(usage=usage, version="1.0 2022-11-10")
    parser.description = 'Introduction: create index file for data products  '

    parser.add_option("-s", "--save_path",
                      action="store", dest="save_path",
                      help="the save path of index file")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
