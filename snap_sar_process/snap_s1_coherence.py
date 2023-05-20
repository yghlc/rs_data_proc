#!/usr/bin/env python
# Filename: snap_s1_coherence.py 
"""
introduction:  produce coherence images from Sentinel-1 imagery using SNAP

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 19 March, 2023
"""

import os,sys
from optparse import OptionParser
from datetime import datetime
import time
machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.basic as basic
import basic_src.io_function as io_function
import basic_src.timeTools as timeTools
import basic_src.map_projection as map_projection

import basic_src.RSImageProcess as RSImageProcess

import cmd_snap

def cal_coherence_from_two_s1(ref_sar, second_sar, res_meter,save_dir, polarisation='VH', tmp_dir=None, wktAoi=None, dem_path=None,
                              thread_num=16):
    # the ref_sar and second_sar should be have the same path, and frame, but don't check here.
    t0 = time.time()

    out_name = cmd_snap.get_granule_name_substr(ref_sar) + '_' + cmd_snap.get_granule_name_substr(second_sar) + '_%s_Coh_S1' % (polarisation)
    save_path = os.path.join(save_dir,out_name + '.tif')
    save_meta = os.path.join(save_dir,out_name + '.json')
    if os.path.isfile(save_path):
        if os.path.isfile(save_path):
            print(save_path + ' already exists, skip')
            return save_path
    snap_intermediate_files = []

    subswath_list = ['IW1','IW2','IW3']
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    if tmp_dir is None:
        tmp_dir = save_dir
    if os.path.isdir(tmp_dir) is False:
        io_function.mkdir(tmp_dir)

    io_function.write_metadata(['reference-image','second-image'], [ref_sar,second_sar], filename=save_meta)
    io_function.write_metadata(['output_resolution_m','polarisation'], [res_meter,polarisation], filename=save_meta)
    io_function.write_metadata(['wktAOI'], [wktAoi], filename=save_meta)

    swath_split_orb_dict = {}
    for sar_zip in [ref_sar, second_sar]:
        granule_name = os.path.basename(sar_zip).split('.')[0]
        for swatch in subswath_list:
            # TOPSAR-Split
            t1 = time.time()
            split = cmd_snap.run_TOPSAR_Split(sar_zip, granule_name, polarisation, swatch, tmp_dir, wktAoi,thread_num=thread_num)
            if split is None:
                continue
            io_function.write_metadata([swatch + '-' + 'split', 'Split-cost-time'], [split, time.time() - t1], filename=save_meta)
            snap_intermediate_files.append(split)
            # Apply-Orbit-File
            t1 = time.time()
            split_orb = cmd_snap.run_Apply_Orbit_File(split, granule_name, tmp_dir, sar_type='sentinel-1', thread_num=thread_num)
            io_function.write_metadata([swatch + '-' + 'Orb', 'Apply_Orbit-cost-time'], [split_orb,time.time() - t1], filename=save_meta)
            snap_intermediate_files.append(split_orb)
            swath_split_orb_dict.setdefault(swatch, []).append(split_orb)

    # print(swath_split_orb_dict)

    coh_res_list = []
    for swatch in subswath_list:
        if swatch not in swath_split_orb_dict.keys():
            continue
        # Back-Geocoding
        t1 = time.time()
        out_stack = cmd_snap.run_Back_Geocoding(swath_split_orb_dict[swatch][0],swath_split_orb_dict[swatch][1],polarisation,swatch,tmp_dir,
                                       dem_path,thread_num=thread_num)
        io_function.write_metadata([swatch + '-' + 'co-registration', 'Back_Geocoding-cost-time' ], [out_stack, time.time()-t1], filename=save_meta)
        snap_intermediate_files.append(out_stack)

        # Coherence
        t1 = time.time()
        out_coherence = cmd_snap.run_Coherence(out_stack,tmp_dir,thread_num=thread_num)
        io_function.write_metadata([swatch + '-' + 'coherence', 'Coherence-cost-time'], [out_coherence, time.time() - t1], filename=save_meta)
        snap_intermediate_files.append(out_coherence)

        # TOPSAR-Deburst
        t1 = time.time()
        out_deb = cmd_snap.run_TOPSAR_Deburst(out_coherence, tmp_dir,thread_num=thread_num)
        io_function.write_metadata([swatch + '-' + 'deburst', 'TOPSAR-Deburst-cost-time'], [out_deb, time.time() - t1], filename=save_meta)
        snap_intermediate_files.append(out_deb)

        coh_res_list.append(out_deb)

    ########################################
    # TOPSAR-Merge
    if len(coh_res_list) > 1:
        t1 = time.time()
        out_merge = cmd_snap.run_TOPSAR_Merge(ref_sar, second_sar,polarisation,coh_res_list,tmp_dir,thread_num=thread_num)
        io_function.write_metadata(['TOPSAR-Merge','TOPSAR-Merge-cost-time'], [out_merge, time.time() - t1], filename=save_meta)
        snap_intermediate_files.append(out_merge)
    else:
        out_merge = coh_res_list[0]

    # Terrain-Correction
    t1 = time.time()
    out_tc = cmd_snap.run_Terrain_Correction(out_merge,tmp_dir,res_meter,dem_path,thread_num=thread_num)
    io_function.write_metadata(['Terrain-correction', 'Terrain_Correction-cost-time'], [out_tc, time.time() - t1], filename=save_meta)
    snap_intermediate_files.append(out_tc)

    # export to tif
    t1 = time.time()
    cmd_snap.export_to_tiff(out_tc,save_path)
    io_function.write_metadata(['Coherence-tiff', 'export_tif-cost-time'], [save_path, time.time() - t1], filename=save_meta)

    t1 = time.time()
    tif_8bit = cmd_snap.sar_coh_to_8bit(save_path)
    io_function.write_metadata(['Coherence-tiff-8bit', 'to_8bit-cost-time'], [tif_8bit, time.time() - t1], filename=save_meta)

    io_function.write_metadata(['Process-time'], [str(datetime.now())], filename=save_meta)

    # clean files
    for tmp in snap_intermediate_files:
        print('removing %s'%tmp)
        io_function.delete_file_or_dir(tmp)
        tmp = tmp.replace('.dim', '.data')
        io_function.delete_file_or_dir(tmp)

    print(datetime.now(), 'Complete, took %s seconds' % (time.time() - t0))
    io_function.write_metadata(['cost-time-second','thread-number'], [time.time() - t0, thread_num], filename=save_meta)
    io_function.write_metadata(['hostname'], [os.uname()[1]], filename=save_meta)


def test_cal_coherence_from_two_s1():
    print(datetime.now(), 'run testing of test_cal_coherence_from_two_s1')
    work_dir=os.path.expanduser('~/Data/tmp_data/InSAR_snap_tutorial')
    ref_sar = os.path.join(work_dir,'S1B_IW_SLC__1SDV_20191208T163422_20191208T163449_019276_024651_E2F6.zip')
    sec_sar = os.path.join(work_dir,'S1B_IW_SLC__1SDV_20191220T163421_20191220T163448_019451_024BE6_79C7.zip')
    ext_shp = os.path.join(work_dir,'extent_aoi','test_aoi.shp')
    tmp_dir = os.path.join(work_dir,'tmp')
    save_dir = os.path.join(work_dir,'save_results')

    wktAoi = vector_gpd.shapefile_to_ROIs_wkt(ext_shp)
    wktAoi = '\"' + wktAoi[0]  + '\"'  # only use the first one
    Polarisations = ['VH', 'VV']
    out_res = 10    # 10 meter
    # print(wktAoi)
    cal_coherence_from_two_s1(ref_sar, sec_sar, out_res, save_dir, polarisation='VH', tmp_dir=tmp_dir, wktAoi=wktAoi)

def main(options, args):
    # if the input file is in a json file
    common_polars = ['VH', 'VV']
    if args[0].endswith('.json'):
        input_dict = io_function.read_dict_from_txt_json(args[0])
        ref_sar = input_dict['reference_sar']
        sec_sar = input_dict['second_sar']
        ref_polarisations = input_dict['ref_polarisations']
        sec_polarisations = input_dict['sec_polarisations']

        ext_shp = input_dict['aoi_shp']
        save_dir = input_dict['save_dir']
        if 'temp_dir' in input_dict.keys():
            tmp_dir = input_dict['temp_dir']
        else:
            tmp_dir = save_dir
        out_res = input_dict['save_pixel_size']
        dem_file = input_dict['elevation_file']
        setting_json = input_dict['env_setting']
        # process_num = input_dict['process_num']
        thread_num = input_dict['thread_num']

        ref_polars = ref_polarisations.split('+')
        sec_polars = sec_polarisations.split('+')
        common_polars = [item for item in set(ref_polars).intersection(sec_polars)]

    else:
        ref_sar = args[0]
        sec_sar = args[1]
        ext_shp = options.aoi_shp
        save_dir = options.save_dir
        tmp_dir = options.temp_dir if options.temp_dir is not None else save_dir
        out_res = options.save_pixel_size
        dem_file = options.elevation_file
        setting_json = options.env_setting
        # process_num = options.process_num
        thread_num = options.thread_num

    # global  baseSNAP, gdal_translate
    if os.path.isfile(setting_json):
        env_setting = io_function.read_dict_from_txt_json(setting_json)
        cmd_snap.baseSNAP = env_setting['snap_bin_gpt']
        print(datetime.now(), 'setting  bSNAP gpt:', cmd_snap.baseSNAP)
        cmd_snap.gdal_translate = env_setting['gdal_translate_bin']
        print(datetime.now(), 'gdal_translate:', cmd_snap.gdal_translate)
    else:
        cmd_snap.baseSNAP = os.getenv('SNAP_BIN_GPT')
        if cmd_snap.baseSNAP is None:
            raise ValueError('SNAP_BIN_GPT is not in Environment Variables')
        cmd_snap.gdal_translate = os.getenv('GDAL_TRANSLATE_BIN')
        if cmd_snap.gdal_translate is None:
            raise ValueError('GDAL_TRANSLATE_BIN is not in Environment Variables')

    # SLURM_TMPDIR only exists, after submit a SLURM job on compute canada
    SLURM_TMPDIR = os.getenv('SLURM_TMPDIR')
    if SLURM_TMPDIR is not None:
        tmp_dir = os.path.join(SLURM_TMPDIR,  os.path.basename(tmp_dir))
        basic.outputlogMessage('update tmp_dir to %s'%tmp_dir)


    # test_cal_coherence_from_two_s1()
    if ext_shp is not None:
        wktAoi = vector_gpd.shapefile_to_ROIs_wkt(ext_shp)
        wktAoi = '\"' + wktAoi[0] + '\"'  # only use the first one
    else:
        wktAoi = None

    # Polarisations = ['VH', 'VV']
    for polar in common_polars:
        cal_coherence_from_two_s1(ref_sar, sec_sar, out_res, save_dir, polarisation=polar, tmp_dir=tmp_dir, wktAoi=wktAoi, dem_path=dem_file,
                              thread_num=thread_num)


if __name__ == '__main__':
    usage = "usage: %prog [options] input.json OR reference_sar second_sar "
    parser = OptionParser(usage=usage, version="1.0 2023-3-19")
    parser.description = 'Introduction: produce coherence from Sentinel-1 using SNAP  '

    parser.add_option("-a", "--aoi_shp",
                      action="store", dest="aoi_shp",
                      help="a shapefile containing AOI")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='asf_data',
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


