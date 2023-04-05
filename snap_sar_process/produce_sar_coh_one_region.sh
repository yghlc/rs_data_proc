#!/bin/bash

# Introduction:  produce multiple SAR coherence using SNAP for one regions

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 5 April, 2023

# for test on my mac
cd ~/Data/sar_coherence_mapping/test1/snap_coh_run

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

py=~/codes/PycharmProjects/rs_data_proc/snap_sar_process/process_multiple_SAR.py

ext=~/Data/Arctic/pan_Arctic/extent/SAR_coh_test_region/SAR_coh_test.shp
dem=~/Data/sar_coherence_mapping/test1/arcticDEM_5m/ArcticDEM_SAR_coh_test_10m_epsg4326.tif
save_dir=~/Data/sar_coherence_mapping/results_sar_coh
tmp_dir=/tmp
thread_num=4
max_jobs=1
work_dir=~/Data/sar_coherence_mapping/working_dir
script_dir=~/Data/sar_coherence_mapping/scripts_sar_coh


${py} sar_list.txt -d ${save_dir} -t ${tmp_dir} -j ${max_jobs} -e ${dem} -a ${ext} --thread_num=${thread_num} \
  --working_dir=${work_dir} --script_dir=${script_dir} --b_run_job_local





