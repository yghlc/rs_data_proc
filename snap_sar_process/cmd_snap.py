#!/usr/bin/env python
# Filename: cmd_snap.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 27 April, 2023
"""

import os,sys

deeplabforRS = os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.basic as basic
import basic_src.io_function as io_function
import raster_io
import vector_gpd

import xml.etree.ElementTree as ET
from datetime import datetime

# ---------------------------------------------------------------------------
# Where the Sentinel 1 Toolbox graphing tool exe and GDAL is located
baseSNAP = '/Applications/snap/bin/gpt'
gdal_translate = '/usr/local/bin/gdal_translate'

def update_env_setting(setting_json):
    global baseSNAP, gdal_translate
    if os.path.isfile(setting_json):
        env_setting = io_function.read_dict_from_txt_json(setting_json)
        baseSNAP = env_setting['snap_bin_gpt']
        print(datetime.now(), 'setting  bSNAP gpt:', baseSNAP)
        gdal_translate = env_setting['gdal_translate_bin']
        print(datetime.now(), 'gdal_translate:', gdal_translate)
    else:
        baseSNAP = os.getenv('SNAP_BIN_GPT')
        if baseSNAP is None:
            raise ValueError('SNAP_BIN_GPT is not in Environment Variables')
        gdal_translate = os.getenv('GDAL_TRANSLATE_BIN')
        if gdal_translate is None:
            raise ValueError('GDAL_TRANSLATE_BIN is not in Environment Variables')

def run_TOPSAR_Split(input_sar_zip, granule_name, Polarisations, subswath, save_dir, awk_aoi,thread_num=16):
    output = os.path.join(save_dir,granule_name + '_%s_%s.dim'%(subswath,Polarisations))
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    cmd_str = baseSNAP + ' TOPSAR-Split -q %d '%thread_num + ' -Ssource=%s '%input_sar_zip + ' -PfirstBurstIndex=1 -PlastBurstIndex=9999 ' + \
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

def run_Apply_Orbit_File(input, granule_name, save_dir,sar_type=None, thread_num=16):
    if input.endswith('.dim'):
        output = io_function.get_name_by_adding_tail(input,'Orb')
    else:
        output = os.path.join(save_dir, granule_name + '_Orb.dim')
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    # orbitType: if we set the wrong OrbitType, SNAP still can download correct orbit file for the input SAR data?
    # 'Sentinel Precise (Auto Download)', 'PRARE Precise (ERS1&2) (Auto Download)','DORIS Precise VOR (ENVISAT) (Auto Download)'
    # ERS, Envisat, Sentinel-1
    if sar_type is not None:
        if sar_type.lower() == 'sentinel-1':
            orbitType = "Sentinel Precise (Auto Download)"
        elif sar_type.lower() == 'ers':
            orbitType = "PRARE Precise (ERS1&2) (Auto Download)"
        elif sar_type.lower() == 'envisat':
            orbitType = "DORIS Precise VOR (ENVISAT) (Auto Download)"
        else:
            raise ValueError('unknown sar_type %s' % str(sar_type))
    else:
        # guess the orbitType from filename
        extension = os.path.splitext(input)[1]        # extension: E1, E2, N1 for (ERS 1&2, Envisat)
        extension = extension[1:]                   # remove dot
        if extension.lower() in ['zip', 'dim']:  # Sentinel-1,  dim: after run_TOPSAR_Split
            orbitType = "Sentinel Precise (Auto Download)"
        elif extension.lower() in ['e1', 'e2']: # ERS 1 & 2
            orbitType = "PRARE Precise (ERS1&2) (Auto Download)"
        elif extension.lower() == 'n1':
            orbitType = "DORIS Precise VOR (ENVISAT) (Auto Download)"
        else:
            raise ValueError('unknown orbitType')

    cmd_str = baseSNAP + ' Apply-Orbit-File -q %d -PcontinueOnFail=false -PorbitType="%s" '%(thread_num,orbitType) + \
              '  -t %s %s'%(output, input)

    basic.os_system_exit_code(cmd_str)
    return output

def get_granule_name_substr(file_path):
    name_strs = os.path.splitext(os.path.basename(file_path))[0].split('_')
    date_str = name_strs[5][:8]
    out_str = date_str + '_' + name_strs[9]
    return out_str

def get_ESA_ERS_granule_name_substr(file_path):
    name_strs = os.path.splitext(os.path.basename(file_path))[0].split('_')
    date_str = name_strs[2][6:]
    out_str = date_str
    return out_str

def run_Back_Geocoding(input_ref, input_second, polarisation, subswath, save_dir, dem_path,thread_num=16):
    out_name = get_granule_name_substr(input_ref) + '_' + get_granule_name_substr(input_second) + '_%s_%s_Orb_Stack.dim'%(subswath,polarisation)
    output = os.path.join(save_dir, out_name)
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output

    cmd_str = baseSNAP + ' Back-Geocoding -q %d -PdemResamplingMethod=BILINEAR_INTERPOLATION -PresamplingType=BILINEAR_INTERPOLATION '%thread_num
    if dem_path is not None:
        nodata = raster_io.get_nodata(dem_path)
        if nodata is None:
            basic.outputlogMessage('Warning, nodata is not set in %s, use -9999'%dem_path)
            nodata = -9999
        cmd_str += "-PdemName=\"External DEM\"  -PexternalDEMFile=%s  -PexternalDEMNoDataValue=%s "%(dem_path, str(nodata))
    else:
        cmd_str += ' -PdemName="SRTM 1Sec HGT" '
    cmd_str += ' %s %s -t %s'%(input_ref, input_second, output)
    basic.os_system_exit_code(cmd_str)
    return output

def run_Coherence(input_stack,save_dir, cohWinAz=3,cohWinRg=10,thread_num=16):
    output = io_function.get_name_by_adding_tail(input_stack, 'Coh')
    output = os.path.join(save_dir, os.path.basename(output))
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    cmd_str = baseSNAP + ' Coherence -q %d -SsourceProduct=%s '%(thread_num,input_stack) + ' -PsubtractFlatEarthPhase=true -PsquarePixel=true ' + \
            ' -PcohWinAz=%d -PcohWinRg=%d '%(cohWinAz,cohWinRg) + ' -t %s '%output
    basic.os_system_exit_code(cmd_str)
    return output

def run_TOPSAR_Deburst(input, save_dir,thread_num=16):
    output = io_function.get_name_by_adding_tail(input, 'Deb')
    output = os.path.join(save_dir, os.path.basename(output))
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    cmd_str = baseSNAP + ' TOPSAR-Deburst -q %d -Ssource=%s -t %s'%(thread_num,input, output)
    basic.os_system_exit_code(cmd_str)
    return output

def run_TOPSAR_Merge(ref_sar, second_sar, polarisation, subswatch_list, save_dir,thread_num=16):
    out_name = get_granule_name_substr(ref_sar) + '_' + get_granule_name_substr(second_sar) + '_%s_Orb_Stack_Coh_Deb_merge.dim'%(polarisation)
    output = os.path.join(save_dir, out_name)
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output

    # if len(subswatch_list) < 2:
    #     pass

    subswatch_strs = ' '.join(subswatch_list)
    cmd_str = baseSNAP + ' TOPSAR-Merge -q %d '%thread_num + subswatch_strs + ' -t %s'%output
    basic.os_system_exit_code(cmd_str)
    return output

def run_Terrain_Correction(input, save_dir, out_res_meter, dem_path,thread_num=16):
    output = io_function.get_name_by_adding_tail(input, 'TC')
    output = os.path.join(save_dir, os.path.basename(output))
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output
    cmd_str = baseSNAP + ' Terrain-Correction -q %d -Ssource=%s '%(thread_num,input) + ' -PpixelSpacingInMeter=%f'%out_res_meter
    if dem_path is not None:
        nodata = raster_io.get_nodata(dem_path)
        if nodata is None:
            basic.outputlogMessage('Warning, nodata is not set in %s, use -9999'%dem_path)
            nodata = -9999
        cmd_str += " -PdemName=\"External DEM\" -PexternalDEMFile=%s  -PexternalDEMNoDataValue=%s "%(dem_path, str(nodata))
    else:
        # available DEMs: SRTM 3Sec, SRTM 1Sec HGT, GETASSE30, Copernicus 90m Global DEM, Copernicus 30m Global DEM,
        #                   CDEM  # Copernicus 30m Global DEM is global, better,
        cmd_str += ' -PdemName="Copernicus 30m Global DEM" '
    cmd_str += ' -t %s'%output
    basic.os_system_exit_code(cmd_str)
    return output

def CoregistrationGraph_ERS(input_ref, input_second, save_dir,org_graph='CoregistrationGraph.xml', thread_num=16):

    out_name = io_function.get_name_no_ext(input_ref) + '_' + io_function.get_name_no_ext(input_second) + '_Stack.dim'
    output = os.path.join(save_dir, out_name)

    # copy and modify the graphy file (xml) based on input and output
    copy_graph = os.path.join(save_dir, 'CoregistrationGraph_%s.xml'%io_function.get_name_no_ext(out_name))
    io_function.copy_file_to_dst(org_graph,copy_graph,overwrite=True)

    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output, copy_graph


    tree = ET.parse(org_graph)
    root = tree.getroot()
    # modify the input
    filelist_elem = root.find(".//node[@id='ProductSet-Reader']/parameters/fileList")
    filelist_elem.text = '%s,%s'%(input_ref,input_second)
    # modify the output
    file_elem = root.find(".//node[@id='Write']/parameters/file")
    file_elem.text = output

    # save to xml
    tree.write(copy_graph)

    # run SNAP
    # ${gpt} CoregistrationGraph_edit.xml
    cmd_str = baseSNAP + ' %s '%copy_graph
    basic.os_system_exit_code(cmd_str)
    return output, copy_graph


# remove border noise
def run_remove_border_noise(input,granule_name, save_dir,thread_num=16):
    if input.endswith('.dim'):
        output = io_function.get_name_by_adding_tail(input,'GBN')
    else:
        output = os.path.join(save_dir, granule_name + '_GBN.dim')
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output

    cmd_str = baseSNAP + ' Remove-GRD-Border-Noise -q %d -SsourceProduct=%s '%(thread_num,input)
    cmd_str += ' -t %s'%output
    basic.os_system_exit_code(cmd_str)
    return output

def run_Calibration(input,granule_name, save_dir,thread_num=16):
    if input.endswith('.dim'):
        output = io_function.get_name_by_adding_tail(input,'CAL')
    else:
        output = os.path.join(save_dir, granule_name + '_CAL.dim')
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output

    cmd_str = baseSNAP + ' Calibration -q %d -PoutputBetaBand=false -PoutputSigmaBand=true '%thread_num
    cmd_str += ' -Ssource=%s '%(input)
    cmd_str += ' -t %s'%output
    basic.os_system_exit_code(cmd_str)
    return output

def run_Speckle_Filter(input,granule_name, save_dir,thread_num=16):
    if input.endswith('.dim'):
        output = io_function.get_name_by_adding_tail(input,'SP')
    else:
        output = os.path.join(save_dir, granule_name + '_SP.dim')
    if os.path.isfile(output):
        print(output +' already exists, skip')
        return output

    cmd_str = baseSNAP + ' Speckle-Filter -q %d -Pfilter="Refined Lee" '%thread_num
    cmd_str += ' -Ssource=%s '%(input)
    cmd_str += ' -t %s'%output
    basic.os_system_exit_code(cmd_str)
    return output


def export_to_tiff(input, save_path):
    input_dir = input.replace('.dim', '.data')
    img_paths = io_function.get_file_list_by_ext('.img',input_dir,bsub_folder=False)
    input_list = []
    output_list = []
    if len(img_paths) == 1:
        # img_path = img_paths[0]
        input_list.append(img_paths[0])
        output_list.append(save_path)
    elif len(img_paths) < 1:
        raise ValueError('NO img file in %s'%input_dir)
    else:
        #
        basic.outputlogMessage('warning, multiple img file in %s, will convert all of them to geotiff' % input_dir)
        # raise ValueError('multiple img file in %s' % input_dir)
        for idx, img_path in enumerate(img_paths):
            input_list.append(img_path)
            out = io_function.get_name_by_adding_tail(save_path, io_function.get_name_no_ext(img_path))
            basic.outputlogMessage('%d: %s' % (idx+1, out))
            output_list.append(out)

    for img_path, save_path in zip(input_list, output_list):
        ## trim nodata region
        if os.path.isfile(save_path):
            print('%s already exists, skip exporting to geotiff'%save_path)
            continue
        # calculate the valid region
        mask_tif = os.path.join(input_dir, io_function.get_name_no_ext(img_path) +'_mask.tif')
        cmd_str = 'gdal_calc.py --type=Byte --NoDataValue=0  --calc="A!=0" -A %s --outfile=%s'%(img_path,mask_tif)
        basic.os_system_exit_code(cmd_str)
        assert os.path.isfile(mask_tif)

        # polygonize
        outline_shp = os.path.join(input_dir, io_function.get_name_no_ext(img_path) +'_outline.shp')
        cmd_str = 'gdal_polygonize.py -8  %s  -b 1 -f "ESRI Shapefile" %s'%(mask_tif,outline_shp)
        basic.os_system_exit_code(cmd_str)
        assert os.path.isfile(outline_shp)

        bounding_boxes = vector_gpd.get_vector_file_bounding_box(outline_shp) # minx, miny, maxx, maxy
        box_str = ' '.join([str(item) for item in bounding_boxes])

        # crop, and translate to Geotiff
        # cmd_str = gdal_translate + ' -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer ' + img_path + ' ' + save_path
        # cmd_str = gdal_translate + ' -of GTiff ' + img_path + ' ' + save_path
        # cmd_str = 'gdalwarp -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer -cutline %s -crop_to_cutline '%outline_shp + img_path + ' ' + save_path
        # cmd_str = 'gdalwarp -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer -te %s '%box_str + img_path + ' ' + save_path
        cmd_str = 'gdalwarp -of GTiff -co tiled=yes -co bigtiff=if_safer -te %s '%box_str + img_path + ' ' + save_path
        basic.os_system_exit_code(cmd_str)
    if len(output_list) > 1:
        return output_list
    else:
        return save_path

def sar_coh_to_8bit(input,save_path=None):
    if save_path is None:
        save_path = io_function.get_name_by_adding_tail(input,'8bit')
    # py_8bit = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py')
    # cmd_str = py_8bit + ' -u 0.99 -l 0.01 %s %s '%(input, save_path)
    # basic.os_system_exit_code(cmd_str)
    return save_path

def test_export_to_tiff():
    input = os.path.expanduser('~/Data/sar_coherence_mapping/ALDs_Dawson_Yukon_Lipovsky_2004/tmp_coh386_2295D0/20040613_Orb_20040822_Orb_Stack_Coh_TC.dim')
    save_path = os.path.expanduser('~/Data/sar_coherence_mapping/ALDs_Dawson_Yukon_Lipovsky_2004/crop.tif')
    export_to_tiff(input, save_path)

if __name__ == '__main__':
    pass