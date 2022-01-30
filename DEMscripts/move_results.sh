#!/bin/bash


#dstDir=/rc_scratch/lihu9680/ArcticDEM_results
dstDir=../ArcticDEM_results

function movefile(){
    id_txt=$1
    echo ${id_txt}
    savedir=$(basename ${id_txt})
    savedir=${savedir:0:-13}
    echo $savedir

    mkdir -p ${dstDir}/${savedir}/grid_dem_diffs
    mkdir -p ${dstDir}/${savedir}/grid_dem_diffs_segment_results
    mkdir -p ${dstDir}/${savedir}/dem_headwall_shp_grid
    mkdir -p ${dstDir}/${savedir}/dem_hillshade_newest_top_grid
    mkdir -p ${dstDir}/${savedir}/dem_hillshade_newest_HWLine_grid

    for id in $(cat ${id_txt}); do 
        echo $id
	mv grid_dem_diffs/*grid${id}.*  ${dstDir}/${savedir}/grid_dem_diffs/.
	mv grid_dem_diffs/*grid${id}_*Index.tif  ${dstDir}/${savedir}/grid_dem_diffs/.
	
	mv grid_dem_diffs_segment_results/*grid${id}  ${dstDir}/${savedir}/grid_dem_diffs_segment_results/.
	mv dem_headwall_shp_grid/*grid${id} ${dstDir}/${savedir}/dem_headwall_shp_grid/.

	mv dem_hillshade_newest_top_grid/*grid${id}_hillshade.tif ${dstDir}/${savedir}/dem_hillshade_newest_top_grid/.
	mv dem_hillshade_newest_HWLine_grid/*grid${id}.tif ${dstDir}/${savedir}/dem_hillshade_newest_HWLine_grid/.
	mv dem_hillshade_newest_HWLine_grid/*grid${id}_count.tif ${dstDir}/${savedir}/dem_hillshade_newest_HWLine_grid/.
    done

}

movefile  grid_ids_txt/alaska_extend_simple_grid_ids.txt

exit

for txt in $(ls grid_ids_txt/*_done); do 

    id_txt=${txt:0:-5}
    movefile ${id_txt}

done

