################################################################################
# Converts multiple numpy tracks into an average robust Z
#
# Written by Michael Wolfe
#
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
################################################################################
from __future__ import division
import sys
import numpy as np
import argparse


def calculate_robust_z(signal, rm_nonfinite=True):
    """
    Function to calculate the RobustZ of a 1 dimensional signal

    
    Args:
        signal (np.array): 1-dimensional numpy array of signal to convert
        rm_nonfinite(bool): Remove nonfinite values from the calculation. Set
                            these to 0.
    Returns:
        robustz (np.array): 1-dimensional numpy array of the calculation

    Raises:
        ValueError if signal is not 1-dimensional
    """

    if signal.ndim != 1:
        raise ValueError("Signal needs to be 1 dimensional for robustZ calculation")
    if rm_nonfinite:
        finite = np.isfinite(signal)
    else:
        finite = np.ones(signal.size,dtype=bool)
    robustz = np.zeros(signal.size)
    median = np.median(signal[finite])
    mad = np.median(np.abs(signal[finite] - median))
    robustz[finite] = (signal[finite]-median)/(1.4826*mad)
    return robustz

def calculate_zscore(signal, mads, rm_nonfinite=True):
    """
    Function to calculate a zero centered z-score given a 1-dimensional signal
    and a median absolute deviation for each point

    If rm_nonfinite is true, sets any non-finite values to 0 after the calculation
 
    Args:
        signal (np.array): 1-dimensional numpy array of signal to convert
        rm_nonfinite(bool): Remove nonfinite values from the calculation. Set
                            these to 0.
    Returns:
        zscore (np.array): 1-dimensional numpy array of the calculation

    Raises:
        ValueError if signal is not 1-dimensional
        ValueError if signal and mad are not the same size
    """
    if signal.ndim != 1 or mads.ndim != 1:
        raise ValueError("Signal needs to be 1 dimensional for zscore calculation")
    if signal.size != mads.size:
        raise ValueError("Signal and mads are not of same length")
    zscore = signal/(1.4826*mads)
    if rm_nonfinite:
        zscore[~np.isfinite(zscore)] = 0
    return zscore



if __name__ == "__main__":

    outpre = sys.argv[1]
    inputfiles = sys.argv[2:]

    arrays = [np.load(input_file) for input_file in inputfiles]
    arrays = [array.flatten() for array in arrays]

    # calculate robust z scores
    robust_zs = [calculate_robust_z(array) for array in arrays]
    
    outfinal = np.mean(robust_zs, axis=0)

    np.save(outpre+"_robust_z", outfinal)
