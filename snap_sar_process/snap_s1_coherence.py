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


# ---------------------------------------------------------------------------
# Where the Sentinel 1 Toolbox graphing tool exe and GDAL is located
baseSNAP = '/home/rcassotto/snap/bin/gpt'
gdal_translate = '/usr/local/bin/gdal_translate'


def run_TOPSAR_Split(input_sar_zip, granule_name, Polarisations, subswath, save_dir, awk_aoi):
    output = os.path.join(save_dir,granule_name + '_%s_%s.dim'%(subswath,Polarisations))
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    cmd_str = baseSNAP + ' TOPSAR-Split ' + ' -Ssource=%s '%input_sar_zip + ' -PfirstBurstIndex=1 -PlastBurstIndex=9999 ' + \
            ' -PselectedPolarisations=%s -Psubswath=%s -t %s '%(Polarisations,subswath,output)
    if awk_aoi is not None:
        cmd_str +=  '  -PwktAoi=%s '%awk_aoi

    # print(cmd_str)
    # basic.os_system_exit_code(cmd_str)
    res, results = basic.exec_command_string(cmd_str)   # need to catch error message if no overlap
    if res != 0:
        if "wktAOI does not overlap any burst" in results:
            basic.outputlogMessage('Warning, wktAOI does not overlap any burst')
            return None
        else:
            print(results)
            sys.exit(res)
    return output

def run_Apply_Orbit_File(input, granule_name, save_dir):
    if input.endswith('.dim'):
        output = io_function.get_name_by_adding_tail(input,'Orb')
    else:
        output = os.path.join(save_dir, granule_name + '_Orb.dim')
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output

    cmd_str = baseSNAP + ' Apply-Orbit-File -PcontinueOnFail=false -PorbitType="Sentinel Precise (Auto Download)" ' + \
              '  -t %s %s'%(output, input)

    basic.os_system_exit_code(cmd_str)
    return output

def get_granule_name_substr(file_path):
    name_strs = os.path.splitext(os.path.basename(file_path))[0].split('_')
    date_str = name_strs[5][:8]
    out_str = date_str + '_' + name_strs[9]
    return out_str

def run_Back_Geocoding(input_ref, input_second, polarisation, subswath, save_dir, dem_path):
    out_name = get_granule_name_substr(input_ref) + '_' + get_granule_name_substr(input_second) + '_%s_%s_Orb_Stack.dim'%(subswath,polarisation)
    output = os.path.join(save_dir, out_name)
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output

    cmd_str = baseSNAP + ' Back-Geocoding -PdemResamplingMethod=BILINEAR_INTERPOLATION -PresamplingType=BILINEAR_INTERPOLATION '
    if dem_path is not None:
        cmd_str += ""   #TODO: add external DEM
    else:
        cmd_str += ' -PdemName="SRTM 1Sec HGT" '
    cmd_str += ' %s %s -t %s'%(input_ref, input_second, output)
    basic.os_system_exit_code(cmd_str)
    return output

def run_Coherence(input_stack,save_dir, cohWinAz=3,cohWinRg=10):
    output = io_function.get_name_by_adding_tail(input_stack, 'Coh')
    output = os.path.join(save_dir, os.path.basename(output))
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    cmd_str = baseSNAP + ' Coherence -SsourceProduct=%s '%input_stack + ' -PsubtractFlatEarthPhase=true -PsquarePixel=true ' + \
            ' -PcohWinAz=%d -PcohWinRg=%d '%(cohWinAz,cohWinRg) + ' -t %s '%output
    basic.os_system_exit_code(cmd_str)
    return output

def run_TOPSAR_Deburst(input, save_dir):
    output = io_function.get_name_by_adding_tail(input, 'Deb')
    output = os.path.join(save_dir, os.path.basename(output))
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    cmd_str = baseSNAP + ' TOPSAR-Deburst -Ssource=%s -t %s'%(input, output)
    basic.os_system_exit_code(cmd_str)
    return output

def run_TOPSAR_Merge(ref_sar, second_sar, polarisation, subswatch_list, save_dir):
    out_name = get_granule_name_substr(ref_sar) + '_' + get_granule_name_substr(second_sar) + '_%s_Orb_Stack_Coh_Deb_merge.dim'%(polarisation)
    output = os.path.join(save_dir, out_name)
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output

    # if len(subswatch_list) < 2:
    #     pass

    subswatch_strs = ' '.join(subswatch_list)
    cmd_str = baseSNAP + ' TOPSAR-Merge ' + subswatch_strs + ' -t %s'%output
    basic.os_system_exit_code(cmd_str)
    return output

def run_Terrain_Correction(input, save_dir, out_res_meter, dem_path):
    output = io_function.get_name_by_adding_tail(input, 'TC')
    output = os.path.join(save_dir, os.path.basename(output))
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    cmd_str = baseSNAP + ' Terrain-Correction -Ssource=%s '%input + ' -PpixelSpacingInMeter=%f'%out_res_meter
    if dem_path is not None:
        cmd_str += ""   #TODO: add external DEM
    else:
        cmd_str += ' -PdemName="SRTM 1Sec HGT" '
    cmd_str += ' -t %s'%output
    basic.os_system_exit_code(cmd_str)
    return output

def export_to_tiff(input, save_path):
    input_dir = input.replace('.dim', '.data')
    img_paths = io_function.get_file_list_by_ext('.img',input_dir,bsub_folder=False)
    if len(img_paths) == 1:
        img_path = img_paths[0]
    elif len(img_paths) < 1:
        raise ValueError('NO img file in %s'%input_dir)
    else:
        raise ValueError('multiple img file in %s' % input_dir)

    # cmd_str = gdal_translate + ' -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer ' + img_path + ' ' + save_path
    cmd_str = gdal_translate + ' -of GTiff ' + img_path + ' ' + save_path
    basic.os_system_exit_code(cmd_str)
    return save_path

def cal_coherence_from_two_s1(ref_sar, second_sar, res_meter,save_dir, polarisation='VH', tmp_dir=None, wktAoi=None, dem_path=None):

    out_name = get_granule_name_substr(ref_sar) + '_' + get_granule_name_substr(second_sar) + '_%s_Coh' % (polarisation)
    save_path = os.path.join(save_dir,out_name + '.tif')
    save_meta = os.path.join(save_dir,out_name + '.json')
    if os.path.isfile(save_path):
        if os.path.isfile(save_path):
            print(save_path + ' already exists, skip')
            return save_path
    snap_intermediate_files = []

    io_function.write_metadata(['reference-image','second-image'], [ref_sar,second_sar], filename=save_meta)
    io_function.write_metadata(['output_resolution_m','polarisation'], [res_meter,polarisation], filename=save_meta)

    subswath_list = ['IW1','IW2','IW3']
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    if tmp_dir is None:
        tmp_dir = save_dir
    if os.path.isdir(tmp_dir) is False:
        io_function.mkdir(tmp_dir)

    swath_split_orb_dict = {}
    for sar_zip in [ref_sar, second_sar]:
        granule_name = os.path.basename(sar_zip).split('.')[0]
        for swatch in subswath_list:
            # TOPSAR-Split
            split = run_TOPSAR_Split(sar_zip, granule_name, polarisation, swatch, tmp_dir, wktAoi)
            if split is None:
                continue
            io_function.write_metadata([swatch + '-' + 'split'], [split], filename=save_meta)
            snap_intermediate_files.append(split)
            # Apply-Orbit-File
            split_orb = run_Apply_Orbit_File(split, granule_name, tmp_dir)
            io_function.write_metadata([swatch + '-' + 'Orb'], [split_orb], filename=save_meta)
            snap_intermediate_files.append(split_orb)
            swath_split_orb_dict.setdefault(swatch, []).append(split_orb)

    # print(swath_split_orb_dict)

    coh_res_list = []
    for swatch in subswath_list:
        if swatch not in swath_split_orb_dict.keys():
            continue
        # Back-Geocoding
        out_stack = run_Back_Geocoding(swath_split_orb_dict[swatch][0],swath_split_orb_dict[swatch][1],polarisation,swatch,tmp_dir,dem_path)
        io_function.write_metadata([swatch + '-' + 'co-registration'], [out_stack], filename=save_meta)
        snap_intermediate_files.append(out_stack)

        # Coherence
        out_coherence = run_Coherence(out_stack,tmp_dir)
        io_function.write_metadata([swatch + '-' + 'coherence'], [out_coherence], filename=save_meta)
        snap_intermediate_files.append(out_coherence)

        # TOPSAR-Deburst
        out_deb = run_TOPSAR_Deburst(out_coherence, tmp_dir)
        io_function.write_metadata([swatch + '-' + 'deburst'], [out_deb], filename=save_meta)
        snap_intermediate_files.append(out_deb)

        coh_res_list.append(out_deb)

    ########################################
    # TOPSAR-Merge
    out_merge = run_TOPSAR_Merge(ref_sar, second_sar,polarisation,coh_res_list,tmp_dir)
    io_function.write_metadata(['TOPSAR-Merge'], [out_merge], filename=save_meta)
    snap_intermediate_files.append(out_merge)

    # Terrain-Correction
    out_tc = run_Terrain_Correction(out_merge,tmp_dir,res_meter,dem_path)
    io_function.write_metadata(['Terrain_Correction'], [out_tc], filename=save_meta)
    snap_intermediate_files.append(out_tc)

    # export to tif
    export_to_tiff(out_tc,save_path)
    io_function.write_metadata(['Coherence_tiff'], [save_path], filename=save_meta)

    io_function.write_metadata(['Process_time'], [str(datetime.now())], filename=save_meta)

    # clean files
    for tmp in snap_intermediate_files:
        print('removing %s'%tmp)
        io_function.delete_file_or_dir(tmp)
        tmp = tmp.replace('.dim', '.data')
        io_function.delete_file_or_dir(tmp)


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
    # grd_file_list = get_grd_file_list(args[0])
    # save_dir = options.save_dir
    # temp_dir = options.temp_dir if options.temp_dir is not None else save_dir
    # pixel_size = options.save_pixel_size
    # dem_file = options.elevation_file
    setting_json = options.env_setting

    global  baseSNAP, gdal_translate
    if os.path.isfile(setting_json):
        env_setting = io_function.read_dict_from_txt_json(setting_json)
        baseSNAP = env_setting['snap_bin_gpt']
        print(datetime.now(), 'setting SNAP gpt:', baseSNAP)
        gdal_translate = env_setting['gdal_translate_bin']
        print(datetime.now(), 'gdal_translate:', gdal_translate)
    else:
        baseSNAP = os.getenv('SNAP_BIN_GPT')
        if baseSNAP is None:
            raise ValueError('SNAP_BIN_GPT is not in Environment Variables')
        gdal_translate = os.getenv('GDAL_TRANSLATE_BIN')
        if gdal_translate is None:
            raise ValueError('GDAL_TRANSLATE_BIN is not in Environment Variables')

    test_cal_coherence_from_two_s1()

    pass

if __name__ == '__main__':
    usage = "usage: %prog [options] grd_files.txt or grd_directory "
    parser = OptionParser(usage=usage, version="1.0 2023-3-19")
    parser.description = 'Introduction: produce coherence from Sentinel-1 using SNAP  '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='asf_data',
                      help="the folder to save pre-processed results")

    parser.add_option("-t", "--temp_dir",
                      action="store", dest="temp_dir",
                      help="the temporal folder for saving intermediate data ")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to run the process")

    parser.add_option("-o", "--out_res",
                      action="store", dest="out_res",type=float,default=10.0,
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


