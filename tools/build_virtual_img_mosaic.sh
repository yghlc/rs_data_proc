#!/bin/bash
set -eE -o functrace

# creating a virtual mosaic of many grid images.


function v_mosaic() {
    local input_dir=$1       # Directory containing input raster files
    local output_vrt=$2      # Output VRT file name

    # Check if required arguments are provided
    if [[ -z "$input_dir" || -z "$output_vrt" ]]; then
        echo "Usage: build_mosaic <input_directory> <output_vrt>"
        return 1
    fi

    # Check if the input directory exists
    if [[ ! -d "$input_dir" ]]; then
        echo "Error: Input directory '$input_dir' does not exist."
        return 1
    fi

    # Find all raster files in the input directory
    raster_files=$(find "$input_dir" -type f \( -iname "*.tif" -o -iname "*.tiff" -o -iname "*.png" -o -iname "*.jpg" \))

    # Check if any raster files were found
    if [[ -z "$raster_files" ]]; then
        echo "Error: No raster files found in directory '$input_dir'."
        return 1
    fi

    # Build the mosaic using gdalbuildvrt
    echo "Building mosaic..."
    gdalbuildvrt -resolution average -r nearest "$output_vrt" $raster_files

    # Check if the VRT was created successfully
    if [[ $? -eq 0 ]]; then
        echo "Mosaic created successfully: $output_vrt"
    else
        echo "Error: Failed to create mosaic."
        return 1
    fi
}

img_dir=$1
save_vrt=$2

v_mosaic ${img_dir} ${save_vrt}
