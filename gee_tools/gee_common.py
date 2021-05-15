#!/usr/bin/env python
# Filename: gee_common.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 03 May, 2021
"""





# // reproject Polar stere
# function reproject(image) {
#   // print(image);
#   // print(image.get('system:index'));
#   var reprojected = image
#     .reproject('EPSG:3413', null, 500);
#   var dateString = ee.String(image.get('system:index'));//.getInfo().replace(/_/g,"-");
#   // dateString = dateString.replace('_','-');  // bug: only replace the first '_'
#   // print(dateString)
#   // rename to the band description to date of snow cover
#   return reprojected.rename(dateString);//.rename('ndsi'.concat(dateString)); //
# }