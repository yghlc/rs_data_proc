#!/usr/bin/env bash

## Introduction: Apply histogram normalization

# run this script in ~/Data/Arctic/canada_arctic/rsImages

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 14 June, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

# set LD_LIBRARY_PATH for gdal_contrast_stretch on server
export LD_LIBRARY_PATH=~/programs/miniconda3/lib:$LD_LIBRARY_PATH
export PATH=~/programs/dans-gdal-scripts_install/bin:$PATH

#gdal_contrast_stretch -histeq <target_stddev>  src.tif dst.tif   Histogram normalize to a target bell curve

function histnormalize() {
    src=$1
    dst=$2
    # gdal_contrast_stretch -histeq <target_stddev>  src.tif dst.tif
    gdal_contrast_stretch -histeq 100 ${src} ${dst}
}

#histnormalize a bash
# cancel
# use the histnormalize during the preprocessing step of converting to 8bit.





