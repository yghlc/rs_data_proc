#!/bin/bash

dst=/home/lhuang/Data/tmp

for ff in $(ls *.tar); do

	filename_noext="${ff%.*}"
	echo $ff  $filename_noext
	mkdir -p ${dst}/${filename_noext}
	tar -xvf $ff -C ${dst}/${filename_noext}

done
