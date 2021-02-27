#!/bin/bash


py=~/codes/PycharmProjects/rs_data_proc/DEMscripts/ArcticDEM_unpack_registration.py

job=10
# delay 3 seconds
delay=3
tar_dir=./
save_dir=arcticdem_registration_tifs

parallel --progress --delay ${delay} -j ${job} ${py} {} -d ${save_dir} -r ::: $(ls ${tar_dir}/*.tar.gz)



