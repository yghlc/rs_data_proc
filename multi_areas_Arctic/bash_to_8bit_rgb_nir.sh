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

function to8bit(){
    region=$1

    in_dir=${region}_daily_mosaic
    rgb_out_dir=${in_dir}_8bit_rgb
    nirGB_out_dir=${in_dir}_8bit_nirGB

    for dd in $(ls -d ${in_dir}/*_???????? ); do
        echo $dd
        folder=$(basename ${dd})

        rgb_save_dir=${rgb_out_dir}/${folder}
        # near-infrared + green + blue
        nirGB_save_dir=${nirGB_out_dir}/${folder}

        mkdir -p ${rgb_save_dir}
        mkdir -p ${nirGB_save_dir}

        for tif in $(ls ${dd}/*.tif);do
            echo $tif

            filename=$(basename "$tif")
		    extension="${filename##*.}"
		    filename_noext="${filename%.*}"

		    # to 8bit
            out8bit=${rgb_out_dir}/${filename_noext}_8bit.tif
#            gdal_translate -ot Byte -scale 0 2000 1 255 -of VRT ${tif} ${out8bit}  # we need to set scale for each band based on histogram
            # for a small region, we may use gdal_contrast_stretch, but if there is cloud, it makes the result bad
            gdal_contrast_stretch -percentile-range 0.01 0.99 ${tif} ${out8bit}

            # get RGB
            outrgb=${rgb_save_dir}/${filename_noext}_8bit_rgb.tif
            gdal_translate -b 3 -b 2 -b 1 -of GTiff -a_nodata 0 -co compress=lzw -co tiled=yes -co bigtiff=if_safer \
            ${out8bit} ${outrgb}

            # get near-infrared + green + blue
            ournirGB=${nirGB_save_dir}/${filename_noext}_8bit_nirGB.tif
            gdal_translate -b 4 -b 2 -b 1 -of GTiff -a_nodata 0 -co compress=lzw -co tiled=yes -co bigtiff=if_safer \
            ${out8bit} ${ournirGB}

            rm ${out8bit}

        done

    done
}

region=WR
to8bit ${region}

region=Banks_east
to8bit ${region}

region=Ellesmere_Island
to8bit ${region}
