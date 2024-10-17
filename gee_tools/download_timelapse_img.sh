#!/bin/bash
set -eE -o functrace

# download a few images each year
py=~/codes/PycharmProjects/rs_data_proc/gee_tools/get_timelapse_img_gee.py
input_shp=../rgc2024-2-thawslumps_buff500.shp

#s_year=1980
#e_year=2024
#image_type='landsat7_pan'

s_year=2015
e_year=2024
image_type='sentinel2_rgb'

buffer_size=1000
max_count=3
save_dir='images_for_2slumps'

# Loop from start year to end year
for (( year=s_year; year<=e_year; year++ )); do
    echo "Year: $year"
    start_date=${year}-06-01
    end_date=${year}-9-30
    months=6,7,8
    could_cover=0.3

    ${py} ${input_shp} ${save_dir} -s ${start_date} -e ${end_date} -m ${months} -c ${could_cover} \
     -b ${buffer_size} -n ${max_count} -i ${image_type}
done



