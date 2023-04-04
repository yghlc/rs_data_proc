#!/usr/bin/env python
# Filename: process_multiple_SAR.py 
"""
introduction: to process multiple SAR granules downloaded from ASF

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 25 March, 2023
"""

import os,sys
from optparse import OptionParser
from datetime import datetime
import time
machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.basic as basic
import basic_src.io_function as io_function
import basic_src.timeTools as timeTools

import pandas as pd

import snap_s1_coherence

def save_sar_meta_to_shape(sar_meta_list,save_shp_path):

    geometry_list = []
    beamModeType_list = []
    file_name_list = []
    flightDirection_list = []
    frameNumber_list = []
    pathNumber_list = []
    orbit_list = []
    platform_list = []
    polarization_list = []
    processingLevel_list = []
    sensor_list = []
    startTime_list = []
    stopTime_list = []
    url_list = []

    for sar_meta in sar_meta_list:
        meta = sar_meta['sar_meta']
        geometry_list.append(meta['geometry'])
        beamModeType_list.append(meta['properties']['beamModeType'])
        file_name_list.append(meta['properties']['fileName'])
        flightDirection_list.append(meta['properties']['flightDirection'])
        pathNumber_list.append(meta['properties']['pathNumber'])
        frameNumber_list.append(meta['properties']['frameNumber'])
        orbit_list.append(meta['properties']['orbit'])
        platform_list.append(meta['properties']['platform'])
        polarization_list.append(meta['properties']['polarization'])
        processingLevel_list.append(meta['properties']['processingLevel'])
        sensor_list.append(meta['properties']['sensor'])
        startTime_list.append(meta['properties']['startTime'])
        stopTime_list.append(meta['properties']['stopTime'])
        url_list.append(meta['properties']['url'])

    polygons = [vector_gpd.json_geometry_to_polygons(item) for item in geometry_list]

    pd_meta = pd.DataFrame({'outline':polygons, 'beamMode':beamModeType_list, 'file_name':file_name_list, 'flightDire':flightDirection_list,
                            'pathNum':pathNumber_list,'frameNum':frameNumber_list,'orbit':orbit_list, 'platform':platform_list,
                            'polarizati':polarization_list,'procLevel':processingLevel_list,'sensor':sensor_list,
                            'startTime':startTime_list, 'stopTime':stopTime_list,'url':url_list})

    vector_gpd.save_polygons_to_files(pd_meta,'outline','EPSG:4326',save_shp_path)
    print('saved %s'%save_shp_path)



def get_sar_file_list(file_or_dir):
    if os.path.isdir(file_or_dir):
        sar_files = io_function.get_file_list_by_ext('.zip',file_or_dir,bsub_folder=False)   # Process VV files
    else:
        with open(file_or_dir,'r') as f_obj:
            sar_files = [line.strip() for line in f_obj.readlines()]
            sar_files = [ os.path.expanduser(item) for item in sar_files]
    if len(sar_files) == 0:
        raise ValueError("No SAR granules in %s"%file_or_dir)
    return sar_files

def organize_sar_pairs(sar_image_list, meta_data_path=None):
    # read sar metadata
    if meta_data_path is None:
        sar_dir = os.path.dirname(sar_image_list[0])
        meta_data_path = os.path.join(sar_dir,'download_data.json') # the default meta when download using ASF API

    file_name_list = [ os.path.basename(item) for item in sar_image_list]
    sar_metas = io_function.read_dict_from_txt_json(meta_data_path)['features']
    # only select meta data for input list
    sar_meta_list = []
    for sar_path, file_name in zip(sar_image_list,file_name_list):
        b_found = False
        for meta in sar_metas:
            if meta['properties']['fileName'] == file_name:
                sar_meta_list.append(meta)
                b_found = True
                break
        if b_found is False:
            basic.outputlogMessage('warning, %s dont have meta data'%file_name)

        # sar_meta_list = [ item for item in sar_metas if item['properties']['fileName'] in file_name_list ]

    # group then to groups with the same path and frame
    group_path_frame = {}
    for s_img, s_meta in zip(sar_image_list,sar_meta_list):
        s_img_meta = {'sar_path':s_img,
                      'sar_meta':s_meta}
        path_frame_str = '%d_%d'%(s_meta['properties']['pathNumber'],s_meta['properties']['frameNumber'])
        group_path_frame.setdefault(path_frame_str, []).append(s_img_meta)

    # print(sar_metas[0])
    for key in group_path_frame.keys():
        print('path_frame:',key)
        for item in group_path_frame[key]:
            print(item['sar_path'])
        # print(group_path_frame[key][0]['sar_path'])

    return group_path_frame


def test_organize_sar_pairs():
    sar_files = get_sar_file_list(os.path.expanduser('~/Data/sar_coherence_mapping/test1/snap_coh_run/s1_data'))
    organize_sar_pairs(sar_files)

def SAR_coherence_samePathFrame(path_frame,sar_meta_list, save_dir,res_meter, tmp_dir=None, wktAoi=None, dem_path=None,thread_num=16,process_num=1):
    save_dir = os.path.join(save_dir,path_frame)
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    # save meta data to disk for checking
    save_sar_meta_to_shape(sar_meta_list, os.path.join(save_dir,'%s_images.shp'%path_frame))

    total_count = len(sar_meta_list)
    if total_count < 2:
        basic.outputlogMessage('Warning, %s have SAR imagery less than 2, cannot calculate coherence'%str(path_frame))
        return False
    # check platform and flightDirection be consistent
    for idx in range(1,total_count):
        if sar_meta_list[0]['sar_meta']['properties']['platform'] != sar_meta_list[idx]['sar_meta']['properties']['platform']:
            [print(item['sar_meta']['properties']['platform']) for item in sar_meta_list]
            raise ValueError('inconsistent platform in SAR images of %s'%path_frame)
        if sar_meta_list[0]['sar_meta']['properties']['flightDirection'] != sar_meta_list[idx]['sar_meta']['properties']['flightDirection']:
            [print(item['sar_meta']['properties']['flightDirection']) for item in sar_meta_list]
            raise ValueError('inconsistent flightDirection in SAR images of %s'%path_frame)
        if sar_meta_list[0]['sar_meta']['properties']['polarization'] != sar_meta_list[idx]['sar_meta']['properties']['polarization']:
            [print(item['sar_meta']['properties']['polarization']) for item in sar_meta_list]
            raise ValueError('inconsistent polarization in SAR images of %s'%path_frame)

    # add date
    for item in sar_meta_list:
        item['acquire_time'] = pd.to_datetime([item['sar_meta']['properties']['startTime']])[0]
    print('before sorting')
    [print(item['sar_path']) for item in sar_meta_list]

    # sort by acquisition time
    sar_meta_list_sorted = sorted(sar_meta_list, key=lambda i: i['acquire_time'])
    print('after sorting')
    [print(item['sar_path']) for item in sar_meta_list_sorted]

    # calculate coherence pair by pair
    for idx in range(1, total_count):
        polarization_list = sar_meta_list_sorted[idx]['sar_meta']['properties']['polarization'].split('+')
        for polari in polarization_list:
            snap_s1_coherence.cal_coherence_from_two_s1(sar_meta_list_sorted[idx-1]['sar_path'], sar_meta_list_sorted[idx]['sar_path'],
                                  res_meter,save_dir, polarisation=polari, tmp_dir=tmp_dir, wktAoi=wktAoi, dem_path=dem_path,
                                                        thread_num=thread_num)



def multiple_SAR_coherence(sar_image_list,save_dir,res_meter, tmp_dir=None, wktAoi=None, dem_path=None,thread_num=16,process_num=1):

    group_path_frame = organize_sar_pairs(sar_image_list, meta_data_path=None)
    # process group by group
    for key in group_path_frame.keys():
        # print('path-frame:',key)
        # for item in group_path_frame[key]:
        #     print(item['sar_path'])
        SAR_coherence_samePathFrame(key,group_path_frame[key],save_dir,res_meter, tmp_dir=tmp_dir, wktAoi=wktAoi,
                                    dem_path=dem_path,thread_num=thread_num,process_num=process_num)


def main(options, args):
    # test_organize_sar_pairs()
    # return

    sar_image_list = get_sar_file_list(args[0])
    sar_image_list = [os.path.abspath(item) for item in sar_image_list]
    # verbose = options.verbose

    ext_shp = options.aoi_shp
    save_dir = options.save_dir
    tmp_dir = options.temp_dir if options.temp_dir is not None else save_dir
    out_res = options.save_pixel_size
    dem_file = options.elevation_file
    setting_json = options.env_setting

    process_num = options.process_num
    thread_num = options.thread_num


    if os.path.isfile(setting_json):
        env_setting = io_function.read_dict_from_txt_json(setting_json)
        snap_s1_coherence.baseSNAP = env_setting['snap_bin_gpt']
        print(datetime.now(), 'setting SNAP gpt:', snap_s1_coherence.baseSNAP)
        snap_s1_coherence.gdal_translate = env_setting['gdal_translate_bin']
        print(datetime.now(), 'gdal_translate:', snap_s1_coherence.gdal_translate)
    else:
        snap_s1_coherence.baseSNAP = os.getenv('SNAP_BIN_GPT')
        if snap_s1_coherence.baseSNAP is None:
            raise ValueError('SNAP_BIN_GPT is not in Environment Variables')
        snap_s1_coherence.gdal_translate = os.getenv('GDAL_TRANSLATE_BIN')
        if snap_s1_coherence.gdal_translate is None:
            raise ValueError('GDAL_TRANSLATE_BIN is not in Environment Variables')


    if ext_shp is not None:
        wktAoi = vector_gpd.shapefile_to_ROIs_wkt(ext_shp)
        wktAoi = '\"' + wktAoi[0] + '\"'  # only use the first one
    else:
        wktAoi = None

    # Polarisations = ['VH', 'VV']
    # read metadata
    multiple_SAR_coherence(sar_image_list, save_dir, out_res, tmp_dir=tmp_dir, wktAoi=wktAoi, dem_path=dem_file,
                           thread_num=thread_num,process_num=process_num)



if __name__ == '__main__':
    usage = "usage: %prog [options] sar_files.txt or image_directory "
    parser = OptionParser(usage=usage, version="1.0 2023-3-25")
    parser.description = 'Introduction: process multiple SAR granules '

    parser.add_option("-a", "--aoi_shp",
                      action="store", dest="aoi_shp",
                      help="a shapefile containing AOI")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='sar_coh_results',
                      help="the folder to save pre-processed results")

    parser.add_option("-t", "--temp_dir",
                      action="store", dest="temp_dir",
                      help="the temporal folder for saving intermediate data ")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=1,
                      help="number of processes to run the process")

    parser.add_option("", "--thread_num",
                      action="store", dest="thread_num", type=int, default=16,
                      help="number of SNAP thread, default is 16, we may change to others"
                           " such as 4 on a supercomputer depending on how much resources we got ")

    parser.add_option("-r", "--save_pixel_size",
                      action="store", dest="save_pixel_size",type=float,default=10.0,
                      help="the resolution for final output")

    parser.add_option("-e", "--elevation_file",
                      action="store", dest="elevation_file",
                      help="DEM file used for terrain correction, if not set, will use SRTM 1 sec ")

    parser.add_option("-s", "--env_setting",
                      action="store", dest="env_setting", default='env_setting.json',
                      help=" the setting of the software environment  ")


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)