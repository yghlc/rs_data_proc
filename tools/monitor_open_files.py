#!/usr/bin/env python
# Filename: monitor_open_files.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 10 March, 2021
"""

import os,sys
import time

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.basic as basic
from datetime import datetime

import psutil

def output_prcess_info(proc_name_contain_str=None):
    # check all open files
    import getpass
    user_name = getpass.getuser()
    all_open_files = []
    for proc in psutil.process_iter():
        try:
            # _proc = proc.as_dict(attrs=['cpu_times', 'name', 'pid', 'status'])
            # print(proc)
            if proc.username() != user_name:
                continue
            if proc_name_contain_str is not None:
                if proc_name_contain_str not in proc.name():
                    continue

            # 'cpu time:', proc.cpu_times
            print(proc.pid, proc.name(), 'user:',proc.username(), 'memory: %.2f'%proc.memory_percent(),
                   'create_time:', datetime.fromtimestamp(proc.create_time()), 'status:',proc.status(), 'is running: ',proc.is_running())

        except psutil.NoSuchProcess:
            print('Error: not such psutil.NoSuchProcess')
            continue
        except psutil.ZombieProcess: #
            print('Error: psutil.ZombieProcess')
            continue
        except:
            print('unknown except')
            continue


def main():

    while True:

        open_files = basic.get_all_processes_openfiles('python')
        print('open_files count by process with python name: %d'%len(open_files))
        for item in open_files:
            print(item)

        output_prcess_info('python')

        time.sleep(1)



if __name__ == '__main__':
    main()
    pass