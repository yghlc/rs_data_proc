#!/usr/bin/env bash

## Introduction: downloaded Planet images

# run this script in ~/Data/Arctic/canada_arctic/rsImages

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 2 February, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

work_dir=~/Data/Arctic/canada_arctic/rsImages
code_dir=~/codes/PycharmProjects/rs_data_proc

cd ${work_dir}


shp_file=~/Data/Arctic/canada_arctic/Willow_River/extent/WR_extent.shp
#shp_file=$1

#cloud_cover_thr=0.1
#cloud_cover_thr=0.03
#cloud_cover_thr=0.05
#cloud_cover_thr=0.10
cloud_cover_thr=0.01 # try download cloudcover less than 0.01 first
item_type=PSScene4Band
process_num=9

account=lingcao.huang@colorado.edu

### download one coverage of images in July&August 2020.
#for year in 2020 ; do
#	echo $year
#    # download images of $year
#    start_date=${year}-07-01
#    end_date=${year}-08-31
#    save_folder=planet_sr_images/${year}_Jul_Aug
#
#    ${code_dir}/planetScripts/download_planet_img.py ${shp_file} ${save_folder} \
#    -s ${start_date} -e ${end_date} -c ${cloud_cover_thr} -i ${item_type} -a ${account}
#done

# download daily images (multiple coverages) in July&August 2020.
start_date=2020-07-01       # start time is 2020-07-01 0:0:0
end_date=2020-09-01         # end time is 2020-09-01 0:0:0
d=${start_date}
while [ "$d" != ${end_date} ]; do
#  echo $d
    if [ $(uname -s) = Darwin ]; then
        # mac option for d decl (the +1d is equivalent to + 1 day)
        n_d=$(date -j -v +1d -f "%Y-%m-%d" ${d} +%Y-%m-%d)
    elif [ $(uname -s) = Linux ]; then
        n_d=$(date -I -d "$d + 1 day")
    else
        echo "Error, supported system, only run on Mac or Linux"
        exit
    fi

    echo $d ${n_d}      # start date and end date
    # wait a few seconds
#    sleep 3

    save_folder=planet_sr_images/${d}

    ${code_dir}/planetScripts/download_planet_img.py ${shp_file} ${save_folder} \
    -s ${d} -e ${n_d} -c ${cloud_cover_thr} -i ${item_type} -a ${account} -p ${process_num}



    # update for next loop
    d=${n_d}
#    exit

done
