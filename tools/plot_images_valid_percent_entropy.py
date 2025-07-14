#!/usr/bin/env python
# Filename: plot_images_valid_percent_entropy.py
"""
introduction: calculate the image shannon_entropy in do some basis statistics

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 14 July, 2025
"""

import os, sys
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))

import basic_src.io_function as io_function
import basic_src.basic as basic
import raster_io

def get_image_valid_percent(image_list, nodata_user=0):
    valid_percent_list = []
    for idx, img in enumerate(tqdm(image_list, desc="Calculating image valid percent")):
        valid_per = raster_io.get_valid_pixel_percentage(img,nodata=nodata_user)
        valid_percent_list.append(valid_per)
    return valid_percent_list

def get_image_entropy(image_list, nodata_user=0):
    valid_percent_list = []
    img_entropy_list = []
    for idx, img in enumerate(tqdm(image_list, desc="Calculating image valid and entropy")):
        valid_per, entropy = raster_io.get_valid_percent_shannon_entropy(img,nodata_input=nodata_user)
        valid_percent_list.append(valid_per)
        img_entropy_list.append(entropy)
    return valid_percent_list, img_entropy_list


def save_values_to_dict_file(image_list, values,save_path):
    save_dict = {}
    for img, v in zip(image_list,values):
        save_dict[os.path.basename(img)] = v
    io_function.save_dict_to_txt_json(save_path,save_dict)

def read_values_from_dict_file(file_path):
    data_dict = io_function.read_dict_from_txt_json(file_path)
    return data_dict


def basic_statistics_and_histogram(values, bins=10, stats_file="statistics.txt", histogram_file="histogram.png"):
    """
    Calculate basic statistics, plot a histogram, and save the results to files.

    Parameters:
        values (list or array-like): List of numerical values.
        bins (int): Number of bins for the histogram. Default is 10.
        stats_file (str): File path to save the statistics. Default is 'statistics.txt'.
        histogram_file (str): File path to save the histogram plot. Default is 'histogram.png'.

    Returns:
        dict: A dictionary containing basic statistics.
    """
    # Ensure the input is a NumPy array for easier computation
    values = np.array(values)

    # Calculate basic statistics
    stats = {
        "count": len(values),
        "mean": np.mean(values),
        "median": np.median(values),
        "std_dev": np.std(values),
        "min": np.min(values),
        "max": np.max(values),
        "25th_percentile": np.percentile(values, 25),
        "75th_percentile": np.percentile(values, 75)
    }

    # Save statistics to a file
    with open(stats_file, "w") as f:
        f.write("Basic Statistics:\n")
        for key, value in stats.items():
            f.write(f"{key.capitalize()}: {value:.2f}\n")
    print(f"Statistics saved to: {stats_file}")

    # Plot histogram
    plt.figure(figsize=(8, 6))
    sns.histplot(values, bins=bins, kde=True, color='skyblue', edgecolor='black', alpha=0.7)
    plt.title("Histogram of Values", fontsize=16)
    plt.xlabel("Values", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    # Save histogram plot to a file
    plt.savefig(histogram_file, dpi=300)
    print(f"Histogram saved to: {histogram_file}")

    # Show the plot (optional)
    # plt.show()

    return stats


def main(options, args):

    img_dir_or_path = args[0]
    file_pattern = options.file_pattern
    nodata_user = options.nodata
    if os.path.isdir(img_dir_or_path):
        image_list = io_function.get_file_list_by_pattern(img_dir_or_path,file_pattern)
    else:
        image_list = [img_dir_or_path]

    pre_name = os.path.splitext(os.path.basename(img_dir_or_path))[0]

    save_entropy = pre_name+'_entropy.json'
    if os.path.isfile(save_entropy):
        print(f'Image entropy has been calculate, read it from {save_entropy}')
        data_dict = read_values_from_dict_file(save_entropy)
        image_entropy_list = list(data_dict.values())
    else:
        valid_percent_list, image_entropy_list = get_image_entropy(image_list, nodata_user=nodata_user)
        save_values_to_dict_file(image_list,image_entropy_list,save_entropy)

    basic_statistics_and_histogram(image_entropy_list,500, stats_file='statistics_entropy.txt',
                                   histogram_file='histogram_entropy.jpg')


    save_valid_percent = pre_name+'_validPercent.json'
    if os.path.isfile(save_valid_percent):
        print(f'Image valid percent has been calculated, read it from {save_valid_percent}')
        data_dict = read_values_from_dict_file(save_valid_percent)
        valid_percent_list = list(data_dict.values())
    else:
        # valid_percent_list = get_image_valid_percent(image_list, nodata_user=nodata_user)
        save_values_to_dict_file(image_list,valid_percent_list,save_valid_percent)

    basic_statistics_and_histogram(valid_percent_list,500, stats_file='statistics_valid_percent.txt',
                                   histogram_file='histogram_valid_percent.jpg')





if __name__ == '__main__':
    usage = "usage: %prog [options] image_folder "
    parser = OptionParser(usage=usage, version="1.0 2021-01-07")
    parser.description = 'Introduction: plot histogram of images '


    parser.add_option("", "--file_pattern",
                      action="store", dest="file_pattern",default='*.tif',
                      help="the pattern to get raster list in a folder")

    parser.add_option("", "--nodata",
                      action="store", dest="nodata",type=int,default='0',
                      help="the user defined Nodata")


    (options, args) = parser.parse_args()
    # print(options.planet_geojson)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)


