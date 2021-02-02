#!/usr/bin/env python
# Filename: learn_retry 
"""
introduction: some code to understand  the retrying model.

# https://github.com/rholder/retrying
# https://pypi.org/project/retrying/

# https://www.jb51.net/article/170781.htm

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 02 February, 2021
"""
import datetime
from retrying import retry

count = 0

# stop_max_attempt_number # defautl is 5 retry, after that, exit
# stop_max_delay # after try amount of this time, exit
# wait_fixed # interval between each time, defautl is 1000 ms

# wait_random_min,  # default is 0
# wait_random_max   # randomly wait time between min and max, default is 1000 ms

# wait_incrementing_increment, # each retry, it will increase wait time, default increase 100 ms.

# wait_exponential_multiplier,
# wait_exponential_max,     # "Wait 2^x *wait_exponential_multiplier milliseconds between each retry,
                            #  x is the previous try numbers, up to wait_exponential_max milliseconds,
                            # then wait_exponential_max seconds afterwards"


# retry_on_exception:

# retry_on_result

# wrap_exception:

# stop_func

# wait_func


@retry(stop_max_attempt_number=5, wait_random_min=1000, wait_random_max=5000)
def run():
    global count
    if count > 0:
        print(datetime.datetime.now(), "start again, %d th try"%count)
    else:
        print(datetime.datetime.now(), "start")
    count += 1
    # raise NameError
    # raise IOError
    raise ValueError

if __name__ == '__main__':
    run()