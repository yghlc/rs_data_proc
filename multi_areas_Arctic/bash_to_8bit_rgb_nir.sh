#!/bin/bash

## Introduction:  bash convert Planet image from 4-band to 3-band (RGB and near-infrared).

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 6 Feb, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace


# set LD_LIBRARY_PATH for gdal_contrast_stretch on server
export LD_LIBRARY_PATH=~/programs/miniconda3/lib:$LD_LIBRARY_PAT
export PATH=~/programs/dans-gdal-scripts_install/bin:$PATH

py=~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py

b1_min=100
b1_max=2000
b2_min=300
b2_max=2200
b3_min=500
b3_max=2400
b4_min=1000
b4_max=4000

# 0 for nodata
dst_min=1
dst_max=255

function to8bit(){
    region=$1

    in_dir=daily_mosaic/${region}_daily_mosaic
    dir_8bit=${in_dir}_8bit
    rgb_out_dir=${in_dir}_8bit_rgb
    nirGB_out_dir=${in_dir}_8bit_nirGB

    for dd in $(ls -d ${in_dir}/*_???????? ); do
        echo $dd
        folder=$(basename ${dd})

        save_dir_8bit=${dir_8bit}/${folder}
        rgb_save_dir=${rgb_out_dir}/${folder}
        # near-infrared + green + blue
        nirGB_save_dir=${nirGB_out_dir}/${folder}

        mkdir -p ${save_dir_8bit}
        mkdir -p ${rgb_save_dir}
        mkdir -p ${nirGB_save_dir}

        for tif in $(ls ${dd}/*.tif);do
            echo $tif

            filename=$(basename "$tif")
		    extension="${filename##*.}"
		    filename_noext="${filename%.*}"

		    # to 8bit
            out8bit=${dir_8bit}/${filename_noext}_8bit.tif
#            gdal_translate -ot Byte -scale 0 2000 1 255 -of VRT ${tif} ${out8bit}

            # we need to set scale for each band based on histogram
            # gdal_transform: scale does not conform to prescribed range https://github.com/OSGeo/gdal/issues/1813
            # although we set dst_min as 1, but in the output, still some 0 pixel (not nodata in the input), the reason is above
#            gdal_translate -ot Byte -scale ${b1_min} ${b1_max} ${dst_min} ${dst_max} \
#            -scale ${b2_min} ${b2_max} ${dst_min} ${dst_max} \
#            -scale ${b3_min} ${b3_max} ${dst_min} ${dst_max} \
#            -scale ${b4_min} ${b4_max} ${dst_min} ${dst_max} \
#             -of VRT ${tif} ${out8bit}

#            ${py} -s ${b1_min} ${b1_max} ${dst_min} ${dst_max} \
#            -s ${b2_min} ${b2_max} ${dst_min} ${dst_max} \
#            -s ${b3_min} ${b3_max} ${dst_min} ${dst_max} \
#            -s ${b4_min} ${b4_max} ${dst_min} ${dst_max} \
#            ${tif} ${out8bit}

            # for a small region, we may use gdal_contrast_stretch, but if there is cloud, it makes the result bad
#            gdal_contrast_stretch -percentile-range 0.01 0.99 ${tif} ${out8bit}
            # use histogram normalization
            gdal_contrast_stretch -histeq 40 ${tif} ${out8bit}

            # get RGB
            outrgb=${rgb_save_dir}/${filename_noext}_8bit_rgb.tif
            gdal_translate -b 3 -b 2 -b 1 -of GTiff -a_nodata 0 -co compress=lzw -co tiled=yes -co bigtiff=if_safer \
            ${out8bit} ${outrgb}

            # get near-infrared + green + blue
            ournirGB=${nirGB_save_dir}/${filename_noext}_8bit_nirGB.tif
            gdal_translate -b 4 -b 2 -b 1 -of GTiff -a_nodata 0 -co compress=lzw -co tiled=yes -co bigtiff=if_safer \
            ${out8bit} ${ournirGB}

#            rm ${out8bit}

        done

    done
}

region=WR
to8bit ${region}

region=Banks_east
to8bit ${region}

region=Ellesmere_Island
to8bit ${region}
