#!/usr/bin/env python
# Filename: python_parallel.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 01 March, 2021
"""
import os,sys

import time
import multiprocessing
from multiprocessing import Pool

def do_something(idx,idx_2):

    proc_id = multiprocessing.current_process().pid
    # print('proc_id',proc_id)
    # time.sleep(idx)
    time.sleep(1)
    # print('proc_id', proc_id, 'Done')


def main():
    prcess_num = 6
    t0 = time.time()

    threadpool = Pool(prcess_num)
    task_count = 30
    para_list = [(idx, idx) for idx, idx in enumerate(range(task_count))]

    for chunksize in range(1,7):

        t0 = time.time()
        # print(para_list)
        threadpool.starmap(do_something, para_list,chunksize=chunksize)
        print('chunksize',chunksize,'total cost time:', time.time() - t0, 'seconds')

    t0 = time.time()
    threadpool.starmap(do_something, para_list)
    print('chunksize', 'default', 'total cost time:', time.time() - t0, 'seconds')
    pass

# it seems larger chunksize end in lower speed.  May keep it default.
# chunksize 1 total cost time: 85.0051908493042 seconds
# chunksize 2 total cost time: 99.00416207313538 seconds
# chunksize 3 total cost time: 114.00342988967896 seconds
# chunksize 4 total cost time: 108.00419902801514 seconds
# chunksize 5 total cost time: 135.00311398506165 seconds
# chunksize 6 total cost time: 159.0044059753418 seconds
# chunksize default total cost time: 99.00354313850403 seconds

# from the above, we can see for task with different computing time, smaller chunksize run quick.
# for the one with task of 30, it seems the default chunksize is 2.

# if the computing time for each time is the same (time.sleep(1)), we got results as follows,
# so, still, chunksize=1 run faster.

# chunksize 1 total cost time: 5.007109880447388 seconds
# chunksize 2 total cost time: 6.0046021938323975 seconds
# chunksize 3 total cost time: 6.003887176513672 seconds
# chunksize 4 total cost time: 8.003803968429565 seconds
# chunksize 5 total cost time: 5.003170013427734 seconds
# chunksize 6 total cost time: 6.00536584854126 seconds
# chunksize default total cost time: 6.006118297576904 seconds



if __name__ == '__main__':
    main()
    pass