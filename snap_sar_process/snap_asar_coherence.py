#!/usr/bin/env python
# Filename: snap_s1_coherence.py 
"""
introduction:  produce coherence images from ERS 1&2, Envisat

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 6 May, 2023
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

def cal_coherence_from_two_ERS(ref_sar, second_sar, res_meter,save_dir, polarisation='VV', tmp_dir=None, ext_shp=None, dem_path=None,
                              thread_num=16,coregister_graph='CoregistrationGraph.xml',rm_tmp_files=True):
    # the ref_sar and second_sar should be have the same path, and frame, but don't check here.
    t0 = time.time()

    ref_short_name = cmd_snap.get_ESA_ERS_granule_name_substr(ref_sar)
    sec_short_name = cmd_snap.get_ESA_ERS_granule_name_substr(second_sar)

    extension = os.path.splitext(ref_sar)[1]        # extension: E1, E2, N1 for (ERS 1&2, Envisat)
    extension = extension[1:]                   # remove dot

    out_name = ref_short_name + '_' + sec_short_name + '_%s_Coh_%s' % (polarisation, extension)
    save_path = os.path.join(save_dir,out_name + '.tif')
    save_meta = os.path.join(save_dir,out_name + '.json')
    if os.path.isfile(save_path):
        if os.path.isfile(save_path):
            print(save_path + ' already exists, skip')
            return save_path
    snap_intermediate_files = []

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    if tmp_dir is None:
        tmp_dir = save_dir
    if os.path.isdir(tmp_dir) is False:
        io_function.mkdir(tmp_dir)

    io_function.write_metadata(['reference-image','second-image'], [ref_sar,second_sar], filename=save_meta)
    io_function.write_metadata(['output_resolution_m','polarisation'], [res_meter,polarisation], filename=save_meta)
    io_function.write_metadata(['extent_shp'], [ext_shp], filename=save_meta)

    # Apply-Orbit-File
    t1 = time.time()
    ref_orb = cmd_snap.run_Apply_Orbit_File(ref_sar, ref_short_name, tmp_dir, thread_num=thread_num)
    io_function.write_metadata([ref_short_name + '-' + 'Orb', 'Apply_Orbit-cost-time'], [ref_orb, time.time() - t1],
                               filename=save_meta)
    snap_intermediate_files.append(ref_orb)

    t1 = time.time()
    sec_orb = cmd_snap.run_Apply_Orbit_File(second_sar, sec_short_name, tmp_dir, thread_num=thread_num)
    io_function.write_metadata([sec_short_name + '-' + 'Orb', 'Apply_Orbit-cost-time'], [sec_orb, time.time() - t1],
                               filename=save_meta)
    snap_intermediate_files.append(sec_orb)

    # Co-registration
    t1 = time.time()
    ERS_stack, graph = cmd_snap.CoregistrationGraph_ERS(ref_orb,sec_orb,tmp_dir,org_graph=coregister_graph)
    io_function.write_metadata(['Co-registration', 'Co-registration-graph', 'Co-registration-cost-time'],
                               [ERS_stack, os.path.basename(graph), time.time() - t1], filename=save_meta)
    snap_intermediate_files.append(ERS_stack)
    io_function.movefiletodir(graph,save_dir)

    # Coherence
    t1 = time.time()
    out_coherence = cmd_snap.run_Coherence(ERS_stack,tmp_dir)
    io_function.write_metadata(['Coherence', 'Coherence-cost-time'], [out_coherence, time.time() - t1], filename=save_meta)
    snap_intermediate_files.append(out_coherence)

    # Terrain-Correction
    t1 = time.time()
    out_tc = cmd_snap.run_Terrain_Correction(out_coherence,tmp_dir,res_meter,dem_path,thread_num=thread_num)
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
    if rm_tmp_files:
        for tmp in snap_intermediate_files:
            print('removing %s'%tmp)
            io_function.delete_file_or_dir(tmp)
            tmp = tmp.replace('.dim', '.data')
            io_function.delete_file_or_dir(tmp)

    print(datetime.now(), 'Complete, took %s seconds' % (time.time() - t0))
    io_function.write_metadata(['cost-time-second','thread-number'], [time.time() - t0, thread_num], filename=save_meta)
    io_function.write_metadata(['hostname'], [os.uname()[1]], filename=save_meta)


def test_cal_coherence_from_two_ERS():
    print(datetime.now(), 'run testing of test_cal_coherence_from_two_s1')
    data_dir=os.path.expanduser('~/Library/CloudStorage/OneDrive-UniversityofVictoria/transfer/ERS1_SNAP_tutorial')
    ref_sar = os.path.join(data_dir,'SAR_IMS_1PNESA19950421_061615_00000017G145_00163_19690_0000.E1')
    sec_sar = os.path.join(data_dir,'SAR_IMS_1PNESA19950526_061619_00000017G146_00163_20191_0000.E1')
    work_dir=os.path.expanduser('~/Data/tmp_data/InSAR_snap_tutorial')
    ext_shp = os.path.expanduser('~/Data/Arctic/pan_Arctic/extent/SAR_coh_test_region/ALDs_Dawson_Yukon_Lipovsky_2004.shp')
    tmp_dir = os.path.join(work_dir,'tmp')
    save_dir = os.path.join(work_dir,'save_results_ERS1')
    co_register_graph = os.path.join(work_dir,'CoregistrationGraph_edit.xml')

    # wktAoi = vector_gpd.shapefile_to_ROIs_wkt(ext_shp)
    # wktAoi = '\"' + wktAoi[0]  + '\"'  # only use the first one
    Polarisations = ['VH', 'VV']
    out_res = 10    # 10 meter
    # print(wktAoi)
    cal_coherence_from_two_ERS(ref_sar, sec_sar, out_res, save_dir, polarisation='VV', tmp_dir=tmp_dir, ext_shp=ext_shp,
                               coregister_graph=co_register_graph)

def main(options, args):
    # if the input file is in a json file
    common_polars = ['VV']
    if args[0].endswith('.json'):
        input_dict = io_function.read_dict_from_txt_json(args[0])
        ref_sar = input_dict['reference_sar']
        sec_sar = input_dict['second_sar']
        ref_polarisations = input_dict['ref_polarisations']
        sec_polarisations = input_dict['sec_polarisations']

        ext_shp = input_dict['aoi_shp']
        save_dir = input_dict['save_dir']
        tmp_dir = input_dict['temp_dir'] if 'temp_dir' in input_dict.keys() else save_dir

        out_res = input_dict['save_pixel_size']
        dem_file = input_dict['elevation_file']
        setting_json = input_dict['env_setting']
        # process_num = input_dict['process_num']
        thread_num = input_dict['thread_num']

        ref_polars = ref_polarisations.split('+')
        sec_polars = sec_polarisations.split('+')
        common_polars = [item for item in set(ref_polars).intersection(sec_polars)]

        coregister_graph = input_dict['coregister_graph'] if 'coregister_graph' in input_dict.keys() else None

        b_rm_tmp_files = input_dict['b_dont_remove_tmp_files'] if 'b_dont_remove_tmp_files' in input_dict.keys() else True

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
        coregister_graph = options.coregister_graph
        b_rm_tmp_files = options.b_dont_remove_tmp_files

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


    # Polarisations = ['VH', 'VV']
    for polar in common_polars:
        cal_coherence_from_two_ERS(ref_sar, sec_sar, out_res, save_dir, polarisation=polar, tmp_dir=tmp_dir, ext_shp=ext_shp, dem_path=dem_file,
                              thread_num=thread_num,coregister_graph=coregister_graph,rm_tmp_files=b_rm_tmp_files)


if __name__ == '__main__':
    usage = "usage: %prog [options] input.json OR reference_sar second_sar "
    parser = OptionParser(usage=usage, version="1.0 2023-5-6")
    parser.description = 'Introduction: produce coherence from ERS 1&2 and EnviSat, using SNAP  '

    # test_cal_coherence_from_two_ERS()
    # sys.exit(0)

    parser.add_option("-a", "--aoi_shp",
                      action="store", dest="aoi_shp",
                      help="a shapefile containing AOI")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='asf_data',
                      help="the folder to save pre-processed results")

    parser.add_option("-t", "--temp_dir",
                      action="store", dest="temp_dir",
                      help="the temporal folder for saving intermediate data ")

    parser.add_option("", "--coregister_graph",
                      action="store", dest="coregister_graph", default='CoregistrationGraph.xml',
                      help="the template graph file for ERS co-registration ")

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

    parser.add_option("", "--b_dont_remove_tmp_files",
                      action="store_false", dest="b_dont_remove_tmp_files",default=True,
                      help="if set, then dont remove intermediate files created by SNAP")


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)


