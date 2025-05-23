#!/usr/bin/env python
# Filename: select_polygons_in_many_big_files.py 
"""
introduction: they are many polygons (> 100 million) saved in nearly 20 big GPKG files after segmentation and merging,
Now, randomly select a few of them  (~1000) for testing algorithm

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 11 May, 2025
"""

import os, sys
from optparse import OptionParser

import fiona
import geopandas as gpd
import random
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import gc  # Garbage collection module

import subprocess
from osgeo import ogr

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
# import basic_src.io_function as io_function
import vector_gpd

# Step 1: Load feature counts from the text file
def load_feature_counts(file_path):
    feature_counts = []
    file_names = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split(", Feature Count: ")
            file_name = parts[0].strip()
            count = int(parts[1].strip())
            file_names.append(file_name)
            feature_counts.append(count)
    return file_names, feature_counts

# Step 2: Generate random global indices
def generate_random_indices(total_count, num_samples):
    return sorted(random.sample(range(total_count), num_samples))

# Step 3: Map global indices to file and local indices
def map_indices_to_files(global_indices, cumulative_counts):
    file_indices = []
    for idx in global_indices:
        file_idx = np.searchsorted(cumulative_counts, idx, side="right") - 1
        local_idx = idx - cumulative_counts[file_idx]
        file_indices.append((file_idx, local_idx))
    return file_indices

def group_indices_by_file(file_indices):
    file_to_indices = {}
    # Group indices by file to minimize file reads
    for file_idx, local_idx in file_indices:
        if file_idx not in file_to_indices:
            file_to_indices[file_idx] = []
        file_to_indices[file_idx].append(local_idx)

    return file_to_indices

# # Step 4: Extract selected polygons using geopandas
def extract_selected_polygons(file_names, file_indices):
    selected_polygons = []

    dataset_crs = None  # Variable to store the CRS

    # Group indices by file to minimize file reads
    file_to_indices = group_indices_by_file(file_indices)

    # Process each file only once
    for file_idx, indices in file_to_indices.items():
        file = file_names[file_idx]
        print(f"Processing file: {file}, extracting {len(indices)} polygons")

        # Load the GeoDataFrame for the current file
        gdf = gpd.read_file(file)

        # Set the CRS if it's not already set
        if dataset_crs is None:
            dataset_crs = gdf.crs

        # Extract only the polygons at the required indices
        for idx in indices:
            selected_polygons.append(gdf.iloc[idx])

        # Release memory by deleting the GeoDataFrame and forcing garbage collection
        del gdf
        gc.collect()

    return selected_polygons, dataset_crs


# fiona, which allows selective reading of specific features (rows) from a file
# without loading the entire file into memory.
# You can use fiona directly to read only the rows you need, more efficient if the file is large
## not working well, although save memory. Need more effort in code engineering

#def extract_selected_polygons(file_names, file_indices):
#    selected_polygons = []
#    file_to_indices = {}
#
#    # Group indices by file for efficient processing
#    for file_idx, local_idx in file_indices:
#        if file_idx not in file_to_indices:
#            file_to_indices[file_idx] = []
#        file_to_indices[file_idx].append(local_idx)
#
#    # Process each file only once
#    for file_idx, indices in file_to_indices.items():
#        file = file_names[file_idx]
#        print(f"Processing file: {file}, extracting {len(indices)} polygons")
#
#        # Sort the indices for efficient sequential access
#        indices = sorted(indices)
#
#        # Open the file with fiona and read sequentially
#        with fiona.open(file) as src:
#            for i, feature in enumerate(src):  # Read features sequentially
#                if i in indices:  # Check if the current index is in the list
#                    if feature and feature.get("geometry"):  # Ensure the geometry exists
#                        selected_polygons.append(feature)
#                    else:
#                        print(f"Warning: Invalid or missing geometry at index {i} in file {file}")
#
#                # Stop early if all indices are processed
#                if len(indices) == 0:
#                    break
#
#    return selected_polygons


def random_select_from_gpkg_ogr2ogr(input_file, output_file, layer_name, num_samples):
    """
    Randomly selects polygons from a GeoPackage file using ogr2ogr and saves them to a new GeoPackage.

    Args:
        input_file (str): Path to the input GeoPackage file.
        output_file (str): Path to the output GeoPackage file.
        layer_name (str): The name of the layer to process.
        num_samples (int): The number of random polygons to select.

    Returns:
        None
    """
    # Construct the SQL query for random selection
    sql_query = f"SELECT * FROM {layer_name} ORDER BY RANDOM() LIMIT {num_samples}"

    # Build the ogr2ogr command
    ogr2ogr_cmd = [
        "ogr2ogr",
        "-f", "GPKG",             # Output format: GeoPackage
        output_file,              # Output file
        input_file,               # Input file
        "-sql", sql_query,        # SQL query for random selection
        "-nln", layer_name        # Name of the output layer
    ]

    # Execute the command
    try:
        subprocess.run(ogr2ogr_cmd, check=True)
        print(f"Randomly selected {num_samples} features from layer '{layer_name}' in '{input_file}' and saved to '{output_file}'.")
    except subprocess.CalledProcessError as e:
        print(f"Error running ogr2ogr: {e}")
        raise RuntimeError(f"ogr2ogr command failed: {e}")

# Step 5: Save selected polygons to a new GPKG file
def save_selected_polygons(selected_polygons, output_file, crs):
    if not selected_polygons:
        print("No valid polygons were selected. Skipping save.")
        return

    # Convert the list of selected polygons to a GeoDataFrame
    gdf_selected = gpd.GeoDataFrame(selected_polygons)
    # Set the CRS in the GeoDataFrame
    gdf_selected = gdf_selected.set_crs(crs)

    # Save to a Shapefile using GeoPandas
    gdf_selected.to_file(output_file, driver="ESRI Shapefile")
    print(f"Saved {len(selected_polygons)} polygons to {output_file}")


def random_select_polygons_in_multi_gpkg(file_names, samp_counts_in_files, output_file):
    """
    Randomly selects polygons from multiple GeoPackage files using ogr2ogr,
    merges the results into a single GeoPackage, and removes intermediate files.

    Args:
        file_names (list): List of input GeoPackage file paths.
        samp_counts_in_files (list): List of sample counts for each file.
        output_file (str): Path to the output GeoPackage file.

    Returns:
        None
    """

    if os.path.isfile(output_file):
        print(f'The final output: {output_file} already exist, cancel extracting polygons')
        return

    # Temporary files list to store intermediate small GeoPackages
    temp_files = []

    for file_idx, input_file in enumerate(file_names):
        sample_count = samp_counts_in_files[file_idx]

        # Open the input GeoPackage to get the layer name
        driver = ogr.GetDriverByName("GPKG")
        datasource = driver.Open(input_file, 0)  # Open in read-only mode
        if not datasource:
            raise RuntimeError(f"Could not open file: {input_file}")

        # Get the first layer name
        layer = datasource.GetLayer()
        layer_name = layer.GetName()
        datasource = None  # Close the datasource

        # Define a temporary file name for the small GeoPackage
        temp_file = f"temp_{file_idx}.gpkg"
        temp_files.append(temp_file)

        # Call random_select_from_gpkg_ogr2ogr to create the small GeoPackage
        random_select_from_gpkg_ogr2ogr(input_file, temp_file, layer_name, sample_count)


    vector_gpd.merge_vector_files(temp_files,output_file,format='ESRI Shapefile')

    # Remove temporary files
    for temp_file in temp_files:
        os.remove(temp_file)
        print(f"Removed temporary file: {temp_file}")

    print(f"All selected polygons have been merged into '{output_file}'.")



# Main Function
def main(options, args):
    # Input files
    feature_count_file = args[0]
    output_file = options.save_path
    num_samples = options.num_samples  # Number of polygons to select

    # Step 1: Load feature counts
    file_names, feature_counts = load_feature_counts(feature_count_file)
    total_count = sum(feature_counts)

    # Step 2: Generate random global indices
    random_indices = generate_random_indices(total_count, num_samples)

    # Step 3: Map global indices to file and local indices
    cumulative_counts = np.cumsum([0] + feature_counts)
    file_indices = map_indices_to_files(random_indices, cumulative_counts)

    # # Step 4: Extract selected polygons
    # selected_polygons, dataset_crs = extract_selected_polygons(file_names, file_indices)

    # # Step 5: Save selected polygons
    # save_selected_polygons(selected_polygons, output_file, dataset_crs)


    # Group indices by file to minimize file reads
    file_to_indices = group_indices_by_file(file_indices)
    # print(file_to_indices.keys())
    samp_counts_in_files = [len(file_to_indices[item]) for item in file_to_indices.keys()]
    random_select_polygons_in_multi_gpkg(file_names, samp_counts_in_files, output_file)




# Run the script
if __name__ == '__main__':
    usage = "usage: %prog [options] feature-count.txt "
    parser = OptionParser(usage=usage, version="1.0 2025-5-11")
    parser.description = 'Introduction: select random from many GPKG files '

    parser.add_option("-s", "--save_path",
                      action="store", dest="save_path", default='random_polygons.shp',
                      help="the folder to save the output")

    parser.add_option("-n", "--num_samples",
                      action="store", dest="num_samples", type=int, default=1000,
                      help="the number of sample for selection")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
