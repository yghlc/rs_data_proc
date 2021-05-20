#!/bin/bash


for d in $(ls -d headwall_shps\** ); do

	echo $d
	end=$(echo $d |cut -c15-)
	new=headwall_shps${end}
	echo $new

	mv $d $new

done
