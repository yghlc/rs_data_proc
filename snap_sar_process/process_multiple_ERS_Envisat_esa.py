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
import slurm_utility
import parallel_run_slurm

import pandas as pd

import snap_asar_coherence

working_dir = os.path.abspath('./')
script_dir = os.path.expanduser('~/Data/sar_coherence_mapping/scripts_sar_coh')
group_account = 'def-tlantz'    # need to set account as this on computeCanada
user_name = 'lingcao'           # default user name on computeCanada
b_run_job_local=False
max_job_count=2

setting_json = None

from process_multiple_SAR import get_sar_file_list

def save_sar_meta_to_shape(sar_meta_list,save_shp_path):

    geometry_list = []
    sar_id_list = []
    dateTime_list = []
    platform_list = []
    swath_list = []
    flightDirection_list = []
    frameNumber_list = []
    pathNumber_list = []
    polarization_list = []
    url_list = []

    for sar_meta in sar_meta_list:
        meta = sar_meta['sar_meta']
        sar_id_list.append(meta['id'])
        geometry_list.append(meta['footprint'])
        flightDirection_list.append(meta['orbit'])  #  flightDirection
        pathNumber_list.append(meta['path'])
        frameNumber_list.append(meta['frame'])
        swath_list.append(meta['swath'])
        platform_list.append(meta['platform'])
        polarization_list.append(meta['polarisation'])
        dateTime_list.append(meta['date'])
        url_list.append(meta['url'])

    polygons = [vector_gpd.wkt_string_to_polygons(item) for item in geometry_list]

    pd_meta = pd.DataFrame({'outline':polygons,  'sar_id':sar_id_list, 'flightDire':flightDirection_list,
                            'pathNum':pathNumber_list,'frameNum':frameNumber_list, 'platform':platform_list,
                            'polarizati':polarization_list, 'swath':swath_list,
                            'dateTime':dateTime_list,'url':url_list})

    vector_gpd.save_polygons_to_files(pd_meta,'outline','EPSG:4326',save_shp_path)
    print('saved %s'%save_shp_path)

    sar_path_list = [item['sar_path'] for item in sar_meta_list ]
    save_path_txt = os.path.splitext(save_shp_path)[0] + '_list.txt'
    io_function.save_list_to_txt(save_path_txt,sar_path_list)


def process_one_pair(sar_meta_list_sorted, ref_idx, sec_idx, path_frame_str, res_meter, save_dir, tmp_dir, ext_shp, dem_path,
                     cohWinAz, cohWinRg, thread_num):

    parallel_run_slurm.b_run_job_local = b_run_job_local
    parallel_run_slurm.slurm_username = user_name
    parallel_run_slurm.wait_if_reach_max_jobs(max_job_count, 'coh')

    flight_direction = sar_meta_list_sorted[ref_idx]['sar_meta']['orbit'][:1] # A or D
    job_name = 'coh%s%s%d' % (path_frame_str,flight_direction,ref_idx)
    parallel_run_slurm.check_length_jobname(job_name,length=14)   # compute canada allow longer name (14 characters)

    # save input parameters to json

    work_dir = os.path.abspath(os.path.join(working_dir, job_name))      # work dir for a pair
    tmp_dir = os.path.join(tmp_dir, 'tmp_' + job_name)

    json_path = os.path.join(work_dir, 'sar_pair_for_coh.json')

    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
        os.chdir(work_dir)
        basic.outputlogMessage('tmp_dir is %s' % tmp_dir)

        input_para_dict = {
            'reference_sar': sar_meta_list_sorted[ref_idx]['sar_path'],
            'second_sar': sar_meta_list_sorted[sec_idx]['sar_path'],
            'ref_polarisations': sar_meta_list_sorted[ref_idx]['sar_meta']['polarisation'],
            'sec_polarisations': sar_meta_list_sorted[sec_idx]['sar_meta']['polarisation'],
            'aoi_shp': ext_shp,
            'save_dir': save_dir,
            'temp_dir': tmp_dir,
            'save_pixel_size': res_meter,
            'cohWinAz':cohWinAz,
            'cohWinRg':cohWinRg,
            'elevation_file': dem_path,
            'env_setting': setting_json,
            'thread_num': thread_num,
            'coregister_graph':'CoregistrationGraph.xml'
        }
        io_function.save_dict_to_txt_json(json_path,input_para_dict)

        # bash for run
        sh_list = ['sar_coh_pair_asar.sh', 'job_sar_coh_pair_thread_asar.sh','env_setting.json','CoregistrationGraph.xml']
        parallel_run_slurm.copy_curc_job_files(script_dir, work_dir, sh_list)
        slurm_utility.modify_slurm_job_sh('job_sar_coh_pair_thread_asar.sh', 'job-name', job_name)

    else:
        os.chdir(work_dir)
        basic.outputlogMessage('tmp_dir is %s' % tmp_dir)

        # SLURM_TMPDIR may changed, need to update it
        input_para_dict = io_function.read_dict_from_txt_json(json_path)
        input_para_dict['temp_dir'] = tmp_dir
        io_function.save_dict_to_txt_json(json_path,input_para_dict)

        if b_run_job_local is False:
            submit_job_names = slurm_utility.get_submited_job_names(user_name)
            if job_name in submit_job_names:
                print('The folder: %s already exist and the job has been submitted, skip submitting a new job' % work_dir)
                return

        # job is completed
        if os.path.isfile('done.txt'):
            print('The job in the folder: %s is Done' % work_dir)
            return

    # submit the job
    parallel_run_slurm.submit_job_curc_or_run_script_local('job_sar_coh_pair_thread_asar.sh', 'sar_coh_pair_asar.sh')

    os.chdir(curr_dir_before_start)


def organize_sar_pairs_ERS_Envisat_esa(sar_image_list, meta_data_path=None):
    # organize the sar data download from ESA
    if len(sar_image_list) < 1:
        return None
    sar_dir = os.path.dirname(sar_image_list[0])
    if meta_data_path is None:
        meta_data_path = os.path.join(sar_dir, 'download_data.json')  # the default filename of metadata when downloading
    else:
        if os.path.isfile(meta_data_path) is False:
            meta_data_path = os.path.join(sar_dir, meta_data_path)
    if os.path.isfile(meta_data_path) is False:
        basic.outputlogMessage('error: %s does not exists'%meta_data_path)
        return None
    file_name_list = [os.path.splitext(os.path.basename(item))[0] for item in sar_image_list]        # exclude extension
    sar_metas = io_function.read_dict_from_txt_json(meta_data_path)
    # only select meta data that have download data
    sel_sar_meta_list = []
    sel_sar_list = []
    for meta in sar_metas:
        if meta['id'] in file_name_list:
            idx = file_name_list.index(meta['id'])
            sel_sar_meta_list.append(meta)
            sel_sar_list.append(sar_image_list[idx])
        else:
            basic.outputlogMessage('Warning, %s is not downloaded or ignored in the input list(txt)' % meta['id'])

    # group them based on the same path and frame
    group_path_frame = {}
    for s_img, s_meta in zip(sel_sar_list, sel_sar_meta_list):
        s_img_meta = {'sar_path': s_img,
                      'sar_meta': s_meta}
        path_frame_str = '%d_%d' % (s_meta['path'], s_meta['frame'])
        group_path_frame.setdefault(path_frame_str, []).append(s_img_meta)

    # print(sar_metas[0])
    for key in group_path_frame.keys():
        print('path_frame:', key)
        for item in group_path_frame[key]:
            print(item['sar_path'])
        # print(group_path_frame[key][0]['sar_path'])

    return group_path_frame

def test_organize_sar_pairs_ERS_Envisat_esa():
    sar_type = 'ERS'
    sar_files = get_sar_file_list(os.path.expanduser('~/Data/asar_ERS_Envisat/ERS'), sar_type=sar_type)
    organize_sar_pairs_ERS_Envisat_esa(sar_files,meta_data_path='ALDs_Dawson_Yukon_Lipovsky_2004_meta.json')


def SAR_coherence_samePathFrame(path_frame,sar_meta_list, save_dir,res_meter, tmp_dir=None, ext_shp=None, dem_path=None,
                                cohWinAz=3,cohWinRg=10,thread_num=16,process_num=1):
    # add decending or ascending to path_frame
    path_frame_direction = path_frame+'_'+sar_meta_list[0]['sar_meta']['orbit'][:3]
    save_dir = os.path.join(save_dir,path_frame_direction)
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
        if sar_meta_list[0]['sar_meta']['orbit'] != sar_meta_list[idx]['sar_meta']['orbit']:
            [print(item['sar_meta']['orbit']) for item in sar_meta_list]
            raise ValueError('inconsistent flightDirection in SAR images of %s'%path_frame)


    # add date
    for item in sar_meta_list:
        item['acquire_time'] = pd.to_datetime([item['sar_meta']['date']])[0]
    print('before sorting')
    [print(item['sar_path']) for item in sar_meta_list]

    # sort by acquisition time
    sar_meta_list_sorted = sorted(sar_meta_list, key=lambda i: i['acquire_time'])
    print('after sorting')
    [print(item['sar_path']) for item in sar_meta_list_sorted]

    # calculate coherence pair by pair
    for idx in range(1, total_count):

        process_one_pair(sar_meta_list_sorted,idx-1,idx,path_frame,res_meter,save_dir,tmp_dir,ext_shp,dem_path,cohWinAz,cohWinRg,thread_num)
        time.sleep(1)   # wait one second, as suggested by computeCanada



def multiple_SAR_coherence(sar_type, sar_image_list,save_dir,res_meter, tmp_dir=None, ext_shp=None, dem_path=None,
                           cohWinAz=3,cohWinRg=10,thread_num=16,process_num=1):

    ext_base_name = io_function.get_name_no_ext(ext_shp)
    ROIs_wkt = vector_gpd.shapefile_to_ROIs_wkt(ext_shp)
    if len(ROIs_wkt) < 1:
        raise ValueError('There is zero AOI')
    for idx, roi_wkt in enumerate(ROIs_wkt):
        if len(ROIs_wkt) == 1:
            data_meta_file = '%s_meta.json'%ext_base_name
        else:
            data_meta_file = '%s_meta_%d.json' % (ext_base_name, idx)

        if sar_type.lower() == 'ers':
            group_path_frame = organize_sar_pairs_ERS_Envisat_esa(sar_image_list, meta_data_path=data_meta_file)
        elif sar_type.lower() == 'envisat':
            group_path_frame = organize_sar_pairs_ERS_Envisat_esa(sar_image_list, meta_data_path=data_meta_file)
        else:
            raise ValueError('unknown sar_type %s' % str(sar_type))

        if group_path_frame is None:
            continue

        # process group by group
        for key in group_path_frame.keys():
            # print('path-frame:',key)
            # for item in group_path_frame[key]:
            #     print(item['sar_path'])
            SAR_coherence_samePathFrame(key,group_path_frame[key],save_dir,res_meter, tmp_dir=tmp_dir, ext_shp=ext_shp,
                                        dem_path=dem_path,cohWinAz=cohWinAz,cohWinRg=cohWinRg,thread_num=thread_num,process_num=process_num)


def main(options, args):
    # test_organize_sar_pairs_ERS_Envisat_esa()
    # return

    sar_type = options.sar_type
    sar_image_list = get_sar_file_list(args[0],sar_type=sar_type)
    sar_image_list = [os.path.abspath(item) for item in sar_image_list]
    # verbose = options.verbose

    ext_shp = options.aoi_shp
    sar_type = options.sar_type
    save_dir = os.path.abspath(options.save_dir)
    tmp_dir = options.temp_dir if options.temp_dir is not None else save_dir
    tmp_dir = os.path.abspath(tmp_dir)
    out_res = options.save_pixel_size
    dem_file = options.elevation_file

    cohWinAz = options.cohWinAz
    cohWinRg = options.cohWinRg

    process_num = options.process_num
    thread_num = options.thread_num

    global setting_json
    setting_json = options.env_setting
    global working_dir
    if options.working_dir is not None:
        working_dir = os.path.abspath(options.working_dir)   # set working path as a absolute path, otherwise, can cause problem
    global script_dir
    if options.script_dir is not None:
        script_dir = os.path.abspath(options.script_dir)
    global b_run_job_local
    b_run_job_local =  options.b_run_job_local
    global max_job_count
    max_job_count = options.max_job_count

    print('setting_json:',setting_json)
    print('working_dir:',working_dir)
    print('script_dir:',script_dir)
    print('max_job_count:', max_job_count)
    print('b_run_job_local:',b_run_job_local)



    # Polarisations = ['VH', 'VV']
    # read metadata
    multiple_SAR_coherence(sar_type, sar_image_list, save_dir, out_res, tmp_dir=tmp_dir, ext_shp=ext_shp, dem_path=dem_file,
                           cohWinAz=cohWinAz,cohWinRg=cohWinRg,thread_num=thread_num,process_num=process_num)

    # wait all local task finished
    while basic.b_all_process_finish(parallel_run_slurm.local_tasks) is False:
        print(datetime.now(),'wait 5 minutes to let all local tasks to complete')
        time.sleep(60*5)



if __name__ == '__main__':
    usage = "usage: %prog [options] sar_files.txt or image_directory "
    parser = OptionParser(usage=usage, version="1.0 2023-5-13")
    parser.description = 'Introduction: process multiple SAR granules '

    parser.add_option("-a", "--aoi_shp",
                      action="store", dest="aoi_shp",
                      help="a shapefile containing AOI")

    parser.add_option("", "--sar_type",
                      action="store", dest="sar_type", default='ERS',
                      help="the type of SAR data: ERS, Envisat")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='sar_coh_results',
                      help="the folder to save pre-processed results")

    parser.add_option("-t", "--temp_dir",
                      action="store", dest="temp_dir",
                      help="the temporal folder for saving intermediate data ")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=1,
                      help="number of processes to run the process")

    parser.add_option("", "--cohWinAz",
                      action="store", dest="cohWinAz", type=int, default=3,
                      help="The window size (Azimuth) for calculating coherence ")

    parser.add_option("", "--cohWinRg",
                      action="store", dest="cohWinRg", type=int, default=10,
                      help="The window size (Range) for calculating coherence ")


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

    parser.add_option("-j", "--max_job_count",
                      action="store", dest="max_job_count", type=int, default=2,
                      help="number of jobs to submit at the same time")

    parser.add_option("", "--b_run_job_local",
                      action="store_true", dest="b_run_job_local",default=False,
                      help="if set (True), will run the job on the machine instead of submitting a slurm job")

    parser.add_option("-u", "--user_name",
                      action="store", dest="user_name",
                      help="the username of the server")

    parser.add_option("", "--working_dir",
                      action="store", dest="working_dir",
                      help="the working directory")
    parser.add_option("", "--script_dir",
                      action="store", dest="script_dir",
                      help="the directory that template scripts and setting file")


    curr_dir_before_start = os.getcwd()
    print('\ncurrent folder before running a job: ', curr_dir_before_start, '\n')


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)