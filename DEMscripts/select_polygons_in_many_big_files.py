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
import gc  # Garbage collection module

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

# # Step 4: Extract selected polygons using geopandas
# def extract_selected_polygons(file_names, file_indices):
#     selected_polygons = []
#     file_to_indices = {}
#
#     # Group indices by file to minimize file reads
#     for file_idx, local_idx in file_indices:
#         if file_idx not in file_to_indices:
#             file_to_indices[file_idx] = []
#         file_to_indices[file_idx].append(local_idx)
#
#     # Process each file only once
#     for file_idx, indices in file_to_indices.items():
#         file = file_names[file_idx]
#         print(f"Processing file: {file}, extracting {len(indices)} polygons")
#
#         # Load the GeoDataFrame for the current file
#         gdf = gpd.read_file(file)
#
#         # Extract only the polygons at the required indices
#         for idx in indices:
#             selected_polygons.append(gdf.iloc[idx])
#
#         # Release memory by deleting the GeoDataFrame and forcing garbage collection
#         del gdf
#         gc.collect()
#
#     return selected_polygons

def extract_selected_polygons(file_names, file_indices):
    selected_polygons = []
    file_to_indices = {}

    # Group indices by file to minimize file reads
    for file_idx, local_idx in file_indices:
        if file_idx not in file_to_indices:
            file_to_indices[file_idx] = []
        file_to_indices[file_idx].append(local_idx)

    # Process each file only once
    for file_idx, indices in file_to_indices.items():
        file = file_names[file_idx]
        print(f"Processing file: {file}, extracting {len(indices)} polygons")

        # fiona, which allows selective reading of specific features (rows) from a file
        # without loading the entire file into memory.
        # You can use fiona directly to read only the rows you need

        # Open the file with Fiona for selective reading
        with fiona.open(file) as src:
            for idx in indices:
                feature = src[idx]  # Read only the specific feature
                selected_polygons.append(feature)

    # Convert the selected features to a GeoDataFrame
    gdf_selected = gpd.GeoDataFrame.from_features(selected_polygons)

    return gdf_selected

# Step 5: Save selected polygons to a new GPKG file
def save_selected_polygons(selected_polygons, output_file):
    gdf_selected = gpd.GeoDataFrame(selected_polygons)
    gdf_selected.to_file(output_file, driver="ESRI Shapefile")  # driver="GPKG"
    print(f"Saved {len(selected_polygons)} polygons to {output_file}")

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

    # Step 4: Extract selected polygons
    selected_polygons = extract_selected_polygons(file_names, file_indices)

    # Step 5: Save selected polygons
    save_selected_polygons(selected_polygons, output_file)

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
