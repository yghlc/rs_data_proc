#!/usr/bin/env python
# Filename: dem_headwall_medial_axis.py 
"""
introduction: try to extract headwall based on skimage.morphology.medial_axis

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 26 August, 2021
"""


import os,sys
from optparse import OptionParser
import time

from skimage import morphology
import numpy as np
import pandas as pd

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import vector_gpd
import raster_io
import rasterio
from rasterio.features import rasterize
import math
import cv2

import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import basic_src.basic as basic


def slope_tif_to_slope_bin(slope_tif,slope_bin_path,slope_threshold):
    if os.path.isfile(slope_bin_path):
        print('%s exist'%slope_bin_path)
    else:
        slope_data, nodata = raster_io.read_raster_one_band_np(slope_tif)
        bin_slope = np.zeros_like(slope_data,dtype=np.uint8)
        bin_slope[slope_data > slope_threshold] = 1
        bin_slope[slope_data > 88] = 0          # if slope is too large, it may caused by artifacts, so remove them

        # save
        slope_bin = bin_slope*255
        #  set nodata as 0
        if raster_io.save_numpy_array_to_rasterfile(slope_bin,slope_bin_path,slope_tif,nodata=0,compress='lzw',tiled='yes',bigtiff='if_safer') is not True:
            return None

    return slope_bin_path

def slope_bin_to_medial_axis_raster(in_image_path, out_image_path):

    if os.path.isfile(out_image_path):
        print('%s exists, skip slope_bin_to_medial_axis_raster'%out_image_path)
        return out_image_path


    image_np, nodata = raster_io.read_raster_one_band_np(in_image_path)

    out_np = morphology.medial_axis(image_np, mask=None, return_distance=False)
    # out_np,dist = morphology.medial_axis(image_np, mask=None, return_distance=True)
    # bool numpy to uint8
    out_np = out_np.astype(np.uint8)
    # dist = dist.astype(np.float32)

    # save to file
    raster_io.save_numpy_array_to_rasterfile(out_np, out_image_path, in_image_path, compress='lzw',
                                             tiled='yes', bigtiff='if_safer')
    # save distances to file (no need)
    # out_dist_path = io_function.get_name_by_adding_tail(out_image_path,'dist')
    # raster_io.save_numpy_array_to_rasterfile(dist, out_dist_path, in_image_path, compress='lzw',
    #                                          tiled='yes', bigtiff='if_safer')

    return out_image_path

def medial_axis_raster_to_vector(in_medial_axis_tif, out_vector_shp, raster_res=2):

    if vector_gpd.raster2shapefile(in_medial_axis_tif,out_shp=out_vector_shp,connect8=True) is None:
        return None

    vector_shp_buff = io_function.get_name_by_adding_tail(out_vector_shp,'buff')
    if os.path.isfile(vector_shp_buff):
        print('%s exists, skip buffering'%vector_shp_buff)
        return vector_shp_buff


    polys = vector_gpd.read_polygons_gpd(out_vector_shp,b_fix_invalid_polygon=False)

    # calculate the attributes before buffer
    # from the area, we can tell how many pixels in each line (line segment), each pixel have size of 2*2 m^2
    id_list = [item for item in range(len(polys))]
    area_list = [ item.area for item in polys]
    length_list = [ item.length for item in polys]

    pixel_size = raster_res*raster_res
    pixel_count_list = [ item/pixel_size for item in area_list]

    # buffer 0.1 meters, so the width of polygons is around 2.02 meters (2-m ArcticDEM), still have width of one pixel
    polys_buff = [item.buffer(0.01) for item in polys]

    # after buffer, for complex medial axis, it may have holes in the polygons, get hole count
    hole_count_list =[len(list(item.interiors)) for item in polys_buff]

    # save, overwrite out_vector_shp
    wkt = map_projection.get_raster_or_vector_srs_info_proj4(out_vector_shp)
    save_pd = pd.DataFrame({'id':id_list,'area':area_list,'length':length_list,'pixels':pixel_count_list,
                            'holes':hole_count_list,'Polygon':polys_buff})
    vector_gpd.save_polygons_to_files(save_pd, 'Polygon', wkt, vector_shp_buff)

    return vector_shp_buff


def remove_based_on_length_pixel(medial_axis_shp,min_length, max_length,wkt, rm_length_shp):
    polygons, lengths_pixel = vector_gpd.read_polygons_attributes_list(medial_axis_shp,'pixels',b_fix_invalid_polygon=False)

    remain_polygons_idx = []
    # remove relative large but narrow ones.
    remove_count = 0
    for idx, (poly,length) in enumerate(zip(polygons,lengths_pixel)):
        # remove too long or too short ones
        if length > max_length or length < min_length:
            remove_count += 1
            continue
        remain_polygons_idx.append(idx)

    basic.outputlogMessage('remove %d polygons based on length in pixel, remain %d ones saving to %s' %
                           (remove_count, len(remain_polygons_idx), rm_length_shp))

    if len(remain_polygons_idx) < 1:
        return False

    vector_gpd.save_shapefile_subset_as(remain_polygons_idx,medial_axis_shp,rm_length_shp)

    return rm_length_shp


def remove_based_on_hole(medial_axis_shp, max_hole,wkt, rm_hole_shp):

    polygons, holes_count = vector_gpd.read_polygons_attributes_list(medial_axis_shp,'holes',b_fix_invalid_polygon=False)

    remain_polygons_idx = []
    # remove relative large but narrow ones.
    remove_count = 0
    for idx, (poly,holes) in enumerate(zip(polygons,holes_count)):
        # remove too long or too short ones
        if holes > max_hole:
            remove_count += 1
            continue
        remain_polygons_idx.append(idx)

    basic.outputlogMessage('remove %d polygons based on holes, remain %d ones saving to %s' %
                           (remove_count, len(remain_polygons_idx), rm_hole_shp))

    if len(remain_polygons_idx) < 1:
        return False

    vector_gpd.save_shapefile_subset_as(remain_polygons_idx,medial_axis_shp,rm_hole_shp)

    return rm_hole_shp

# def calculate_remove_based_on_line_segments(medial_axis_shp, max_line_segments,wkt, rm_line_segment_shp):
#
#     # note: eventually, this is the same to hole count and does not provide new information.
#     # so it's not necessary to use it.
#
#     # calculate the number of line segments
#     polygons = vector_gpd.read_polygons_gpd(medial_axis_shp,b_fix_invalid_polygon=False)
#     line_segments_list = []
#     for idx, poly in enumerate(polygons):
#         # out_line = poly.exterior
#         in_lines = list(poly.interiors)
#         count = 1 + len(in_lines)
#         line_segments_list.append(count)
#     add_attributes = {'lines':line_segments_list}
#     vector_gpd.add_attributes_to_shp(medial_axis_shp,add_attributes)
#
#     # remove based on the number of line segments
#     remain_polygons_idx = []
#     # remove relative large but narrow ones.
#     remove_count = 0
#     for idx, lines in enumerate(line_segments_list):
#         # remove too long or too short ones
#         if lines > max_line_segments:
#             remove_count += 1
#             continue
#         remain_polygons_idx.append(idx)
#
#     basic.outputlogMessage('remove %d polygons based on the count of line segments, remain %d ones saving to %s' %
#                            (remove_count, len(remain_polygons_idx), rm_line_segment_shp))
#
#     if len(remain_polygons_idx) < 1:
#         return False
#
#     vector_gpd.save_shapefile_subset_as(remain_polygons_idx,medial_axis_shp,rm_line_segment_shp)
#
#     return rm_line_segment_shp


###################################################################################################
# copy the following codes from:
# https://github.com/gabyx/WormAnalysis/blob/master/SkeletonTest/Skeletonize.ipynb
import collections
import itertools
import networkx as nx

class Vertex:
    def __init__(self, point, degree=0, edges=None):
        self.point = np.asarray(point)
        self.degree = degree
        self.edges = []
        self.visited = False
        if edges is not None:
            self.edges = edges

    def __str__(self):
        return str(self.point)
class Edge:
    def __init__(self, start, end=None, pixels=None):
        self.start = start
        self.end = end
        self.pixels = []
        if pixels is not None:
            self.pixels = pixels
        self.visited = False

def buildTree(img, start=None):
    # copy image since we set visited pixels to black
    img = img.copy()
    shape = img.shape
    nWhitePixels = np.sum(img)

    # neighbor offsets (8 nbors)
    nbPxOff = np.array([[-1, -1], [-1, 0], [-1, 1],
                        [0, -1], [0, 1],
                        [1, -1], [1, 0], [1, 1]
                        ])

    queue = collections.deque()

    # a list of all graphs extracted from the skeleton
    graphs = []

    blackedPixels = 0
    # we build our graph as long as we have not blacked all white pixels!
    while nWhitePixels != blackedPixels:

        # if start not given: determine the first white pixel
        if start is None:
            it = np.nditer(img, flags=['multi_index'])
            while not it[0]:
                it.iternext()

            start = it.multi_index

        startV = Vertex(start)
        queue.append(startV)
        # print("Start vertex: ", startV)

        # set start pixel to False (visited)
        img[startV.point[0], startV.point[1]] = False
        blackedPixels += 1

        # create a new graph
        G = nx.Graph()
        G.add_node(startV)

        # build graph in a breath-first manner by adding
        # new nodes to the right and popping handled nodes to the left in queue
        while len(queue):
            currV = queue[0]  # get current vertex
            # print("Current vertex: ", currV)

            # check all neigboor pixels
            for nbOff in nbPxOff:

                # pixel index
                pxIdx = currV.point + nbOff

                if (pxIdx[0] < 0 or pxIdx[0] >= shape[0]) or (pxIdx[1] < 0 or pxIdx[1] >= shape[1]):
                    continue  # current neigbor pixel out of image

                if img[pxIdx[0], pxIdx[1]]:
                    # print( "nb: ", pxIdx, " white ")
                    # pixel is white
                    newV = Vertex([pxIdx[0], pxIdx[1]])

                    # add edge from currV <-> newV
                    G.add_edge(currV, newV, object=Edge(currV, newV))
                    # G.add_edge(newV,currV)

                    # add node newV
                    G.add_node(newV)

                    # push vertex to queue
                    queue.append(newV)

                    # set neighbor pixel to black
                    img[pxIdx[0], pxIdx[1]] = False
                    blackedPixels += 1

            # pop currV
            queue.popleft()
        # end while

        # empty queue
        # current graph is finished ->store it
        graphs.append(G)

        # reset start
        start = None

    # end while

    return graphs, img

def getEndNodes(g):
    # return [ n for n in nx.nodes_iter(g) if nx.degree(g,n) == 1 ]
    return [ n for n in nx.nodes(g) if nx.degree(g,n) == 1 ]


def getLongestPath(graph, endNodes):
    """
        graph is a fully reachable graph = every node can be reached from every node
    """

    if len(endNodes) < 2:
        raise ValueError("endNodes need to contain at least 2 nodes!")

    # get all shortest paths from each endpoint to another endpoint
    allEndPointsComb = itertools.combinations(endNodes, 2)

    maxLength = 0
    maxPath = None

    for ePoints in allEndPointsComb:

        # get shortest path for these end points pairs
        sL = nx.dijkstra_path_length(graph,
                                     source=ePoints[0],
                                     target=ePoints[1])

        # dijkstra can throw if now path, but we are sure we have a path

        # store maximum
        if (sL > maxLength):
            maxPath = ePoints
            maxLength = sL

    if maxPath is None:
        raise ValueError("No path found!")

    return nx.dijkstra_path(graph,
                            source=maxPath[0],
                            target=maxPath[1]), maxLength
##########################################################################################


def calculate_line_segment_polygon_pixels(polygon, raster_res=2):
    '''
    get number of line segments
    :param polygon:
    :param raster_res:
    :param min_line_pixel: line segment shorter than min_line_pixel would be ignored
    :return:
    '''
    # rasterize the polygon to pixels (original from medial axis, width is one pixel)
    # based on the pixels, find the longest line segment, and other segments.

    # raster_io.burn_polygons_to_a_raster(ref_raster, [polygon], 1, save_raster)

    minx, miny, maxx, maxy = vector_gpd.get_polygon_bounding_box(polygon)
    poly_transform = rasterio.transform.from_origin(minx,maxy,raster_res,raster_res)

    height = math.ceil((maxy - miny)/raster_res)
    width = math.ceil((maxx - minx)/raster_res)

    save_dtype = rasterio.uint8

    burn_out = np.zeros((height, width))
    # rasterize the shapes
    burn_shapes = [(item_shape, item_int) for (item_shape, item_int) in zip([polygon], [1])]

    out_label = rasterize(burn_shapes, out=burn_out, transform=poly_transform,
                          fill=0, all_touched=False, dtype=save_dtype)

    # print('pixel count',np.sum(out_label))

    graphs, imgB = buildTree(out_label)
    if len(graphs) != 1:
        raise ValueError('should only have one graph')

    line_graph = graphs[0]

    longest_path,maxLength = getLongestPath(line_graph,getEndNodes(line_graph))
    # print('note of longest_path and length',len(longest_path), maxLength)

    # remove longest path
    line_graph.remove_nodes_from(longest_path)
    # print('after removing longest path, end node count:',len(getEndNodes(graphs[0])))

    # print(len(line_graph.nodes))
    # print(len(line_graph.edges))

    # remove isolated nodes
    remove_list = []
    for idx,node in enumerate(line_graph.nodes):
        # print(line_graph.edges(node))
        # nbs = nx.neighbors(line_graph, node)
        # print(nbs[0], nbs[1])
        neighbours = [n for n in line_graph[node]]
        # print(idx,neighbours)
        if len(neighbours) < 1:
            remove_list.append(node)
    line_graph.remove_nodes_from(remove_list)


    # print(len(line_graph.nodes))
    # print(len(line_graph.edges))

    # # plot to check graph
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # import random
    # # draw all graphs
    # fig = plt.figure(figsize=(12, 12))
    #
    # plt.imshow(out_label, cmap=cm.gray, interpolation='nearest')
    # ax = plt.gca()
    # plt.axis("equal")
    #
    # class PosWrapper(dict):
    #     def __getitem__(self, key):
    #         return [key.point[1], key.point[0]]  # switch x and y
    #
    # for i, g in enumerate(graphs):
    #     nx.draw_networkx(g, PosWrapper(), ax=ax,
    #                      with_labels=False,
    #                      node_color="#%06x" % random.randint(0, 0xFFFFFF),
    #                      edge_color='b', node_size=20
    #                      )
    # plt.show()

    remain_pixel_count = len(line_graph.nodes)
    return remain_pixel_count


    # find the number of line segments, may not need to find the longest line segment.
    # If we can get the end point number, then the number of line segment is count of end point – 1.
    # An end point like the start and end point of a line.
    # In the 8-neighbour, end point should only have one connectted point.
    # since 1 is valid pixel, 0 is background.

    # update
    # it turns out that the ideas of end point not always work, maybe we can check cross points.
    # a cross point connecting two or more line segments should have >=3 pixels in its 8-neighbour


    # kernel = np.ones((3, 3), np.float32)
    # kernel[1,1] = 0     # set the middle one as 0
    # # print(kernel)
    # dst = cv2.filter2D(out_label, -1, kernel,borderType=cv2.BORDER_CONSTANT)
    # # only keep the line pixels
    # dst = dst*out_label
    # loc_end_points = np.where(dst==1)
    #
    # # by removing these end points to remove some short line segments (length=one pixel)
    # # remove line segment small than min_line_pixel
    # # min_line_pixel = 3
    # for idx in range(min_line_pixel):
    #     # if the line is too short or already only have two end points, then skip
    #     if len(loc_end_points[0]) <= 2:
    #         break
    #
    #     out_label[loc_end_points] = 0
    #     # calculate the end points again
    #     dst = cv2.filter2D(out_label, -1, kernel, borderType=cv2.BORDER_CONSTANT)
    #     dst = dst * out_label
    #     loc_end_points = np.where(dst == 1)

    # print(loc_end_points)

    # return number of line segment
    # num_line = len(loc_end_points[0]) -1

    # # debug
    # if num_line < 1:
    #     # plot for testing
    #     # note: the polygon is in the projection of EPSG:3413,
    #     # if we open the shapefile in QGIS but another prjection, it will look different.
    #     import matplotlib.pyplot as plt
    #     # imgplot = plt.imshow(out_label)
    #     fig, axs = plt.subplots(1,2)
    #     axs[0].imshow(out_label)
    #     axs[1].imshow(dst)
    #     plt.show()

    # return num_line


def calculate_remove_based_on_pixel_line_segments(medial_axis_shp,wkt, rm_line_segment_shp,raster_res=2,max_unwant_line_pixel=10):
    '''
    convert polygons to pixels again, then count line segments and find the longest line segments from pixels.
    :param medial_axis_shp:
    :param wkt:
    :param rm_line_segment_shp:
    :param raster_res:
    :param max_unwant_line_pixel:
    :return:
    '''

    # calculate the number of line segments
    polygons = vector_gpd.read_polygons_gpd(medial_axis_shp,b_fix_invalid_polygon=False)
    unwant_line_pixels_list = []
    for idx, poly in enumerate(polygons):
        remain_line_pixel_count = calculate_line_segment_polygon_pixels(poly,raster_res=raster_res)
        unwant_line_pixels_list.append(remain_line_pixel_count)

    # save to file
    add_attributes = {'unwant_pi':unwant_line_pixels_list}
    vector_gpd.add_attributes_to_shp(medial_axis_shp,add_attributes)

    # remove based on the number of line segments
    remain_polygons_idx = []
    # remove relative large but narrow ones.
    remove_count = 0
    for idx, lines in enumerate(unwant_line_pixels_list):
        # remove too long or too short ones
        if lines > max_unwant_line_pixel:
            remove_count += 1
            continue
        remain_polygons_idx.append(idx)

    basic.outputlogMessage('remove %d polygons based on the count of unwanted line pixels, remain %d ones saving to %s' %
                           (remove_count, len(remain_polygons_idx), rm_line_segment_shp))

    if len(remain_polygons_idx) < 1:
        return False

    vector_gpd.save_shapefile_subset_as(remain_polygons_idx,medial_axis_shp,rm_line_segment_shp)

    return rm_line_segment_shp



def extract_headwall_based_medial_axis_from_slope(idx, total, slope_tif, work_dir, save_dir,slope_threshold,
                                                  min_length, max_length, max_hole_count,max_unwanted_line_pixel,process_num):
    '''
    extract headwall from slope based on medial axis (skimage.morphology.medial_axis)
    :param idx: tif index
    :param total: total slope file count
    :param slope_tif: slope file
    :param work_dir:
    :param save_dir:
    :param slope_threshold:
    :param min_length: min length, the medial axis calculated from skimage.morphology has width of one pixel, the length is based pixel count of line segments
    :param max_length: max length (pixel count)
    :param max_hole_count: some complex line segment may end in holes when forming polygons
    :param process_num:
    :return:
    '''

    headwall_shp = os.path.splitext(os.path.basename(io_function.get_name_by_adding_tail(slope_tif, 'headwall')))[0] + '.shp'
    save_headwall_shp = os.path.join(save_dir, headwall_shp)
    if os.path.isfile(save_headwall_shp):
        print('%s exists, skip' % save_headwall_shp)
        return save_headwall_shp

    print('(%d/%d) extracting headwall from %s' % (idx, total, slope_tif))

    wkt = map_projection.get_raster_or_vector_srs_info_wkt(slope_tif)
    # binary slope
    slope_bin_path = os.path.join(work_dir, os.path.basename(io_function.get_name_by_adding_tail(slope_tif, 'bin')))
    if slope_tif_to_slope_bin(slope_tif, slope_bin_path, slope_threshold) is None:
        return False

    # get medial axis raster
    medial_axis_tif = io_function.get_name_by_adding_tail(slope_bin_path,'medial_axis')
    if slope_bin_to_medial_axis_raster(slope_bin_path, medial_axis_tif) is None:
        return False

    # get madial axis vector (polygons: width of one pixel)
    medial_axis_poly_shp = os.path.join(work_dir, io_function.get_name_no_ext(medial_axis_tif) + '_poly.shp')
    medial_axis_poly_shp_buff = medial_axis_raster_to_vector(medial_axis_tif, medial_axis_poly_shp)
    if medial_axis_poly_shp_buff is None:
        return False

    # only keep not too long or too short line segments
    rm_length_shp = io_function.get_name_by_adding_tail(medial_axis_poly_shp_buff, 'rmLength')
    if os.path.isfile(rm_length_shp):
        print('%s exists, skip removing based on length in pixels' % rm_length_shp)
    else:
        # medial_axis_shp,min_length, max_length,wkt, rm_length_shp
        if remove_based_on_length_pixel(medial_axis_poly_shp_buff, min_length, max_length, wkt, rm_length_shp) is False:
            return False

    # remove based on  hole count, if too many holes, not headwall
    rm_hole_shp = io_function.get_name_by_adding_tail(rm_length_shp, 'rmHole')
    if os.path.isfile(rm_hole_shp):
        print('%s exists, skip removing based holes' % rm_hole_shp)
    else:
        # medial_axis_shp,min_length, max_length,wkt, rm_length_shp
        if remove_based_on_hole(rm_length_shp, max_hole_count, wkt, rm_hole_shp) is False:
            return False

    # get line segments in polygons and remove based on line segments
    # eventually, this is the same to hole count and does not provide new information.
    # max_line_count = 10
    # rm_line_shp = io_function.get_name_by_adding_tail(rm_hole_shp, 'rmLine')
    # if os.path.isfile(rm_line_shp):
    #     print('%s exists, skip removing based the count of line segments' % rm_line_shp)
    # else:
    #     # medial_axis_shp,min_length, max_length,wkt, rm_length_shp
    #     if calculate_remove_based_on_line_segments(rm_hole_shp, max_line_count, wkt, rm_line_shp) is False:
    #         return False


    # calculate number of line segments based on pixels.
    rm_line_shp = io_function.get_name_by_adding_tail(rm_hole_shp, 'rmLine')
    if os.path.isfile(rm_line_shp):
        print('%s exists, skip removing based the number of line segments' % rm_line_shp)
    else:
        if calculate_remove_based_on_pixel_line_segments(rm_hole_shp, wkt, rm_line_shp,
                                                  raster_res=2, max_unwant_line_pixel=max_unwanted_line_pixel) is False:
            return False


    #
    # # copy the results.
    # io_function.copy_shape_file(rm_medialAxis_shp, save_headwall_shp)



def test_slope_bin_to_medial_axis():
    data_dir = os.path.expanduser('~/Data/dem_processing/headwall_shp_sub_6174/20080511_dem_slope')
    slope_bin_path = os.path.join(data_dir, '20080511_dem_slope_bin.tif' )
    medial_axis_tif = os.path.join(data_dir, '20080511_dem_slope_bin_medial_axis.tif')

    medial_axis_tif = slope_bin_to_medial_axis_raster(slope_bin_path,medial_axis_tif)

    medial_axis_poly_shp = os.path.join(data_dir, '20080511_dem_slope_bin_medial_axis_poly.shp')
    medial_axis_raster_to_vector(medial_axis_tif, medial_axis_poly_shp)

def test_extract_headwall_based_medial_axis_from_slope():

    data_dir = os.path.expanduser('~/Data/dem_processing/grid_6174_tmp_files/slope_sub_6174')
    slope_tif = os.path.join(data_dir,'20080511_dem_slope.tif')
    work_dir = os.path.expanduser('~/Data/dem_processing')
    save_dir = os.path.expanduser('~/Data/dem_processing/grid_6174_tmp_files')
    slope_threshold = 20
    min_length = 6
    max_length = 500
    max_hole_count = 0
    max_unwanted_line_pixel = 5
    process_num = 1

    extract_headwall_based_medial_axis_from_slope(0, 1, slope_tif, work_dir, save_dir, slope_threshold,
                                                  min_length, max_length, max_hole_count, max_unwanted_line_pixel,process_num)


def test_calculate_line_segment_polygon_pixels():
    one_line_poly_shp = os.path.expanduser('~/Data/dem_processing/test_cal_line_segment/one_line_segment.shp')
    polygons = vector_gpd.read_polygons_gpd(one_line_poly_shp, b_fix_invalid_polygon=False)

    # print(polygons[0])
    line_segment = calculate_line_segment_polygon_pixels(polygons[0])
    print('remain pixel count in lines', line_segment)



def main():
    # test_slope_bin_to_medial_axis()
    test_extract_headwall_based_medial_axis_from_slope()
    # test_calculate_line_segment_polygon_pixels()
    pass


if __name__ == '__main__':
    main()
    pass