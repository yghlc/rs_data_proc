#!/bin/bash

# rename folders from "headwall_shps*_grid00000"  to "headwall_shps*_grid00000"

for d in $(ls -d headwall_shps\** ); do

	echo $d
	end=$(echo $d |cut -c15-)
	new=headwall_shps${end}
	echo $new

	mv $d $new

done
