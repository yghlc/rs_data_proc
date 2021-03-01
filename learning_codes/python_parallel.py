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
    time.sleep(idx)
    # print('proc_id', proc_id, 'Done')


def main():
    prcess_num = 6
    t0 = time.time()

    for chunksize in range(1,7):

        threadpool = Pool(prcess_num)
        para_list = [(idx,idx) for idx,idx in enumerate(range(10))]
        # print(para_list)
        threadpool.starmap(do_something, para_list,chunksize=chunksize)

        print('chunksize',chunksize,'total cost time:', time.time() - t0, 'seconds')
    pass

# it seems larger chunksize end in lower speed.  May keep it default.
# chunksize 1 total cost time: 12.026546955108643 seconds
# chunksize 2 total cost time: 29.040847063064575 seconds
# chunksize 3 total cost time: 50.054343938827515 seconds
# chunksize 4 total cost time: 72.0702850818634 seconds
# chunksize 5 total cost time: 107.08576011657715 seconds
# chunksize 6 total cost time: 137.09962511062622 seconds

if __name__ == '__main__':
    main()
    pass