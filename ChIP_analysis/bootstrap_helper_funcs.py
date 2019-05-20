# Copyright (c) 2018 Michael Wolfe University of Michigan. All rights reserved.
#
#
#Developed by: Michael Wolfe
#University of Michigan
#http://freddolino-lab.med.umich.edu/
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal with
#the Software without restriction, including without limitation the rights to
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#the Software, and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:
#
#Redistributions of source code must retain the above copyright notice, this
#list of conditions and the following disclaimers.  Redistributions in binary
#form must reproduce the above copyright notice, this list of conditions and the
#following disclaimers in the documentation and/or other materials provided with
#the distribution.  Neither the names of Michael Wolfe, University of Michigan,
#nor the names of its contributors may be used to endorse or promote products
#derived from this Software without specific prior written permission.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS
#OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.

import numpy as np
import numba
from math import ceil, floor

# @numba.jit(nopython=True, parallel=True)
def credible_interval(array, alpha = 0.05):
    """ Take a bootstrapped based 95% credible interval

    Args:
        array (np.array): 1xn array (n number of bootstraps)
        alpha (float): alpha for the size of credible interval
    Returns:
        summary_stats (np.array): An array of stats 1 x 3
                                  pos 0: mean
                                  pos 1: min_ci at alpha
                                  pos 2: max_ci at alpha
    """
    out_array = np.zeros(3)
    mean = np.mean(array)
    out_array[0] = mean
    array_sort = bubblesort_jit(array)
    out_array[1] = array_sort[int(floor((alpha/2)*array.size))]
    out_array[2] = array_sort[int(floor(array.size - alpha/2*array.size))]
    return out_array

@numba.jit(nopython=True)
def bubblesort_jit(arr):
    N = arr.shape[0]
    for end in range(N, 1, -1):
        for i in range(end - 1):
            cur = arr[i]
            if cur > arr[i + 1]:
                tmp = arr[i]
                arr[i] = arr[i + 1]
                arr[i + 1] = tmp
    return(arr)

def least_extreme_value(stats):
    """ Take the value closest to zero for the min/max ci

    Args:
        stats (np.array): output from credible interval function
    Returns:
        lev (float): least extreme value
    """
    mean, minci, maxci = stats
    if minci <= 0.0 and maxci >= 0.0:
        return 0.0
    elif minci >= 0.0 and maxci > 0.0:
        return minci
    elif minci < 0.0 and maxci <= 0.0:
        return maxci
    else:
        #print "warning weird value encountered %s %s %s"%(mean, minci, maxci)
        return np.nan

if __name__ == "__main__":
    import sys
    infile = sys.argv[1] 
    outpre = sys.argv[2]

    inmat = np.load(infile)
    wt = inmat[:,:,0]
    ko = inmat[:,:,1]
    # if nan in the ko, don't subtract anything
    ko[~np.isfinite(ko)] = 0
    inmat = wt - ko
    stats = np.apply_along_axis(credible_interval, 1, inmat)
    lev = np.apply_along_axis(least_extreme_value, 1, stats) 
    np.save(outpre+"_mean", stats[:,0])
    np.save(outpre+"_lev", lev)
    np.save(outpre+"_minci",stats[:,1])
    np.save(outpre+"_maxci",stats[:,2])
