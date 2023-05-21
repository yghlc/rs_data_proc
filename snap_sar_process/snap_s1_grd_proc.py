#!/usr/bin/env python
# Filename: snap_s1_grd_proc.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 21 May, 2023
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

import cmd_snap

def pre_process_s1_grd_file(input_grd, save_dir, res_meter, tmp_dir=None, dem_path=None, thread_num=16):
    t0 = time.time()
    out_name = os.path.splitext(os.path.basename(input_grd))[0]
    save_path = os.path.join(save_dir, out_name + '.tif')
    save_meta = os.path.join(save_dir, out_name + '.json')
    # TODO: need to fix:  if save_meta, and there are mutiple geotiff files, then this doesn't work
    if os.path.isfile(save_path):
        if os.path.isfile(save_path):
            print(save_path + ' already exists, skip')
            return save_path, save_meta
    snap_intermediate_files = []

    #
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    if tmp_dir is None:
        tmp_dir = save_dir
    if os.path.isdir(tmp_dir) is False:
        io_function.mkdir(tmp_dir)

    io_function.write_metadata(['input-s1-grd'], [input_grd], filename=save_meta)
    io_function.write_metadata(['output_resolution_m'], [res_meter], filename=save_meta)

    # Apply-Orbit-File
    t1 = time.time()
    out_orb = cmd_snap.run_Apply_Orbit_File(input_grd, out_name, tmp_dir,sar_type='sentinel-1', thread_num=thread_num)
    io_function.write_metadata(['Apply_Orbit-cost-time'], [time.time() - t1],filename=save_meta)
    snap_intermediate_files.append(out_orb)

    # Remove-GRD-Border-Noise
    t1 = time.time()
    out_gbn = cmd_snap.run_remove_border_noise(out_orb,out_name,tmp_dir,thread_num=thread_num)
    io_function.write_metadata(['Remove-GRD-Border-Noise-cost-time'], [time.time() - t1], filename=save_meta)
    snap_intermediate_files.append(out_gbn)

    # Calibration to sigma nought
    t1 = time.time()
    out_cal = cmd_snap.run_Calibration(out_gbn,out_name,tmp_dir,thread_num=thread_num)
    io_function.write_metadata(['Calibration-cost-time'], [time.time() - t1], filename=save_meta)
    snap_intermediate_files.append(out_cal)

    # speckle filter
    t1 = time.time()
    out_sp = cmd_snap.run_Speckle_Filter(out_cal,out_name,tmp_dir,thread_num=thread_num)
    io_function.write_metadata(['Speckle-Filter-cost-time'], [time.time() - t1], filename=save_meta)
    snap_intermediate_files.append(out_sp)

    # terrain correction
    t1 = time.time()
    out_tc = cmd_snap.run_Terrain_Correction(out_sp,tmp_dir,res_meter,dem_path,thread_num=thread_num)
    io_function.write_metadata(['Terrain_Correction-cost-time'], [time.time() - t1], filename=save_meta)
    snap_intermediate_files.append(out_tc)

    # export to tif
    t1 = time.time()
    out_tiffs = cmd_snap.export_to_tiff(out_tc,save_path)
    io_function.write_metadata(['Out-Grd-tiff', 'export_tif-cost-time'], [out_tiffs, time.time() - t1], filename=save_meta)

    io_function.write_metadata(['Process-time'], [str(datetime.now())], filename=save_meta)

    # clean files
    # for tmp in snap_intermediate_files:
    #     print('removing %s'%tmp)
    #     io_function.delete_file_or_dir(tmp)
    #     tmp = tmp.replace('.dim', '.data')
    #     io_function.delete_file_or_dir(tmp)

    print(datetime.now(), 'Complete, took %s seconds' % (time.time() - t0))
    io_function.write_metadata(['cost-time-second','thread-number'], [time.time() - t0, thread_num], filename=save_meta)
    io_function.write_metadata(['hostname'], [os.uname()[1]], filename=save_meta)

    # return save_path,save_meta

def main(options, args):

    if args[0].endswith('.json'):
        input_dict = io_function.read_dict_from_txt_json(args[0])
        print('not support yet')
        input_sar_grd = None
        save_dir = None
        tmp_dir = options.temp_dir if options.temp_dir is not None else save_dir
        out_res = None      #options.save_pixel_size
        dem_file = None     #options.elevation_file
        setting_json = None #options.env_setting
        thread_num = None   #options.thread_num
    else:
        input_sar_grd = args[0]
        save_dir = options.save_dir
        tmp_dir = options.temp_dir if options.temp_dir is not None else save_dir
        out_res = options.save_pixel_size
        dem_file = options.elevation_file
        setting_json = options.env_setting
        thread_num = options.thread_num

    cmd_snap.update_env_setting(setting_json)

    pre_process_s1_grd_file(input_sar_grd, save_dir, out_res, tmp_dir=tmp_dir, dem_path=dem_file, thread_num=thread_num)



if __name__ == '__main__':
    usage = "usage: %prog [options] input.json OR sar.zip "
    parser = OptionParser(usage=usage, version="1.0 2023-5-21")
    parser.description = 'Introduction: process S1 GRD data  '

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