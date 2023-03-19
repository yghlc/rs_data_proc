#!/bin/bash

# Introduction:  produce SAR coherence from Sentinel-1 using SNAP

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 18 March, 2023

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

gpt=/Applications/snap/bin/gpt

# on Ubuntu, there is a problem: org.jblas ERROR Couldn't load copied link file: java.lang.UnsatisfiedLinkError: libgfortran.so.5: cannot open shared object file: No such file or directory.
# see more in Error of Interferogram in SNAP
# fix that by: export LD_LIBRARY_PATH=/home/lihu9680/programs/miniconda3/envs/sar/lib:$LD_LIBRARY_PATH

ref_img=S1B_IW_SLC__1SDV_20191208T163422_20191208T163449_019276_024651_E2F6.zip
sec_img=S1B_IW_SLC__1SDV_20191220T163421_20191220T163448_019451_024BE6_79C7.zip

stack=Process/S1B_IW_SLC__1SDV_20191208T163422_20191208T163449_019276_024651_E2F6_Orb_Stack.dim

function subswatch_coh() {
    subswath=$1
    polar=$2
    # img1 and img2 are those after Apply-Orbit-File
    img1=$3
    img2=$4

    ${gpt} TOPSAR-Split -Ssource=${img1}.dim -PfirstBurstIndex=1 -PlastBurstIndex=9999 -PselectedPolarisations=${polar} -Psubswath=${subswath} -t ${img1}_OB_${subswath}_${polar}
    ${gpt} TOPSAR-Split -Ssource=${img2}.dim -PfirstBurstIndex=1 -PlastBurstIndex=9999 -PselectedPolarisations=${polar} -Psubswath=${subswath} -t ${img2}_OB_${subswath}_${polar}

    ${gpt} Back-Geocoding -PdemName="SRTM 1Sec HGT" -PdemResamplingMethod=BILINEAR_INTERPOLATION -PresamplingType=BILINEAR_INTERPOLATION \
            ${img1}_OB_${subswath}_${polar}.dim ${img2}_OB_${subswath}_${polar}.dim  -t ${subswath}_${polar}_stack

     ${gpt} Coherence -SsourceProduct=${subswath}_${polar}_stack.dim -PcohWinAz=3 -PcohWinRg=10 -PsubtractFlatEarthPhase=true -PsquarePixel=true -t ${subswath}_${polar}_stack_coh

     ${gpt} TOPSAR-Deburst -Ssource=${subswath}_${polar}_stack_coh.dim  -t ${subswath}_${polar}_stack_coh_Deb
}

# 1
${gpt} Apply-Orbit-File -PcontinueOnFail=false -PorbitType="Sentinel Precise (Auto Download)" -t 20191208_OB ${ref_img}
${gpt} Apply-Orbit-File -PcontinueOnFail=false -PorbitType="Sentinel Precise (Auto Download)" -t 20191220_OB ${sec_img}


# 2
#${gpt} TOPSAR-Split -Ssource=20191208_OB.dim -PfirstBurstIndex=1 -PlastBurstIndex=9999 -PselectedPolarisations=VH -Psubswath=IW1 -t 20191208_OB_IW1_VH
#${gpt} TOPSAR-Split -Ssource=20191220_OB.dim -PfirstBurstIndex=1 -PlastBurstIndex=9999 -PselectedPolarisations=VH -Psubswath=IW1 -t 20191220_OB_IW1_VH

# 3 co-registration 
#${gpt} Back-Geocoding -PdemName="SRTM 1Sec HGT" -PdemResamplingMethod=BILINEAR_INTERPOLATION -PresamplingType=BILINEAR_INTERPOLATION  20191208_OB_IW1_VH.dim 20191220_OB_IW1_VH.dim  -t IW1_VH_reg

# 4 Enhanced-Spectral-Diversity (optional, need to fine tune parameters)
#${gpt} Enhanced-Spectral-Diversity -Ssource=IW1_VH_reg.dim -t IW1_VH_reg_stack

# 5
#${gpt} Coherence -SsourceProduct=IW1_VH_reg_stack.dim -PcohWinAz=3 -PcohWinRg=10 -PsubtractFlatEarthPhase=true -PsquarePixel=true -t IW1_VH_reg_stack_coh

# 6
#${gpt} TOPSAR-Deburst -Ssource=IW1_VH_reg_stack_coh.dim  -t IW1_VH_reg_stack_coh_Deb

# repeat 2-6 for IW2, IW3, then run TOPSAR-Merge

# GoldsteinPhaseFiltering (optional)
#${gpt} GoldsteinPhaseFiltering -SsourceProduct=IW1_VH_reg_stack_coh_Deb.dim -t IW1_VH_reg_stack_coh_Deb_fil 

## calculate coherence for each swatch
subswatch_coh IW1 VH 20191208_OB 20191220_OB
subswatch_coh IW2 VH 20191208_OB 20191220_OB
subswatch_coh IW3 VH 20191208_OB 20191220_OB

# Merge
${gpt} TOPSAR-Merge IW1_VH_stack_coh_Deb.dim IW2_VH_stack_coh_Deb.dim  IW3_VH_stack_coh_Deb.dim  -t final_VH_coh_Deb

# Terrain Correction
#${gpt} Terrain-Correction -Ssource=IW1_VH_reg_stack_coh_Deb.dim -PpixelSpacingInMeter=10 -PdemName="SRTM 1Sec HGT" -t IW1_VH_reg_stack_coh_Deb_TC
#${gpt} Terrain-Correction -Ssource=IW1_VH_reg_stack_coh_Deb_fil.dim -PpixelSpacingInMeter=10 -PdemName="SRTM 1Sec HGT" -t IW1_VH_reg_stack_coh_Deb_fil_TC

${gpt} Terrain-Correction -Ssource=final_VH_coh_Deb.dim -PpixelSpacingInMeter=10 -PdemName="SRTM 1Sec HGT" -t final_VH_coh_Deb_TC