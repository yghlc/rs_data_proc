#!/bin/bash

## Introduction:  compose two sar images (8bit), one binary (0 or 255) to 3-band images

# run in ~/Bhaltos2/lingcaoHuang/flooding_area/Houston/compose_3bands_8bit

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

sar8bit_dir=~/Bhaltos2/lingcaoHuang/flooding_area/Houston/Houston_SAR_GRD_FLOAT_gee/S1_Houston_prj_8bit_select
binary_dir=~/Bhaltos2/lingcaoHuang/flooding_area/Houston/Houston_binary/Harvey_8_29_255

nodata=0

#up
b1=S1A_IW_GRDH_1SDV_20170805T002644_20170805T002710_017781_01DCB4_8FDF_prj_8bit.tif
b2=S1A_IW_GRDH_1SDV_20170829T002645_20170829T002710_018131_01E74D_3220_prj_8bit.tif
b3=S1A_IW_GRDH_1SDV_20170829T002645_20170829T002710_018131_01E74D_3220_Sigma0_VH_Ptf_binaryLM_prj_255.tif

out1=S1A_IW_GRDH_1SDV_20170829T002645_20170829T002710_018131_01E74D_3220_prj_3bands_8bit.tif


gdal_merge.py -separate -o ${out1} -n ${nodata} -a_nodata ${nodata} ${sar8bit_dir}/${b1} \
 ${sar8bit_dir}/${b2}  ${binary_dir}/${b3}

#down
b1=S1A_IW_GRDH_1SDV_20170805T002619_20170805T002644_017781_01DCB4_1716_prj_8bit.tif
b2=S1A_IW_GRDH_1SDV_20170829T002620_20170829T002645_018131_01E74D_D734_prj_8bit.tif
b3=S1A_IW_GRDH_1SDV_20170829T002620_20170829T002645_018131_01E74D_D734_Sigma0_VH_Ptf_binaryLM_prj_255.tif

out2=S1A_IW_GRDH_1SDV_20170829T002620_20170829T002645_018131_01E74D_D734_prj_3bands_8bit.tif
gdal_merge.py -separate -o ${out2} -n ${nodata} -a_nodata ${nodata} ${sar8bit_dir}/${b1} \
 ${sar8bit_dir}/${b2}  ${binary_dir}/${b3}








