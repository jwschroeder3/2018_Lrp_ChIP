################################################################################
# bootstrap sampling for ChIP analysis with paired extraced and input files and
# a KO to subtract
#
# Written by Michael Wolfe
# Copyright (c) 2018 Michael Wolfe University of Michigan. All rights reserved.
#
#
#Developed by: Michael Wolfe
#University of Michigan
#http://freddolino-lab.med.umich.edu/
#Modified by : Jeremy Schroeder
#University of Wisconsin - Madison
#http://wanglab.bact.wisc.edu
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
#the distribution.  Neither the names of Michael Wolfe, Jeremy Schroeder 
#University of Michigan, University of Wisconsin - Madison,
#nor the names of their contributors may be used to endorse or promote products
#derived from this Software without specific prior written permission.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS
#OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
################################################################################

# base python
from __future__ import print_function
import argparse
import sys
import logging
import os
import multiprocessing as mp
from math import floor, ceil, log
import re
import time

# Widely available
import numpy as np
from scipy import stats
import pandas as pd
from dfply import *

# More custom
import numba

# custom
from bootstrap_sam_file import ReadSampler
from bootstrap_helper_funcs import least_extreme_value, credible_interval

def summary_stats_only_finite(array, alpha=0.05):
    """ This version of summary stats only considers values that aren't nan or
    inf. Creates summaries a single location based on bootstrap

    Args:
        array (np.array): np.array that is 1 x n (n is the number of bootstrap
                          replicates)

    Returns:
        summary_stats (np.array): np.array that is 1 x 9
                                  pos 0 average
                                  pos 1 min_ci at alpha
                                  pos 2 max_ci at alpha
                                  pos 3 variance
                                  pos 4 least extreme value
                                  pos 5 median
                                  pos 6 median absolute deviation
                                  pos 7 number of infinite values in bootstrap
                                  pos 8 number of nans in bootstrap
    """
    finite_vals = np.isfinite(array)
    num_inf = np.sum(np.isinf(array))
    num_nan = np.sum(np.isnan(array))
    if np.sum(finite_vals) == 0:
        return(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, num_inf, num_nan)
    else:
        these_stats = credible_interval(array[finite_vals], alpha)
        this_lev = least_extreme_value(these_stats)
        var = np.var(array[finite_vals])
        median = np.median(array[finite_vals])
        mad = np.median(np.abs(array[finite_vals] - median))
        return(these_stats[0], these_stats[1], these_stats[2], var, this_lev, these_stats[3], mad, num_inf, num_nan)

def summary_stats_nan_zeros(array, alpha=0.05):
    """ This version of summary stats turns nans to 0. NOTE THAT THIS MODIFIES
    THE ARRAY

    Creates summaries for a single location based on bootstraps

    Args:
        array (np.array): np.array that is 1 x n (n is the number of bootstrap
                          replicates)

    Returns:
        summary_stats (np.array): np.array that is 1 x 9
                                  pos 0 average
                                  pos 1 min_ci at alpha
                                  pos 2 max_ci at alpha
                                  pos 3 variance
                                  pos 4 least extreme value
                                  pos 5 median
                                  pos 6 median absolute deviation
    """
    array[np.isnan(array)] = 0
    these_stats = credible_interval(array, alpha)
    this_lev = least_extreme_value(these_stats)
    var = np.var(array)
    median = np.median(array)
    mad = np.sum(np.abs(array - median))/array.size
    return(these_stats[0], these_stats[1], these_stats[2], var, this_lev, median, mad)


def actual_reads(sampler, size, res=1):
    """ Take all reads from a sampler and map them to an array

    Args:
        sampler (class ReadSampler): object holding reads to sample
        size (int): size of numpy array to create
        res (optional, int): resolution the numpy array is in
    """
    ############# modify to support jit acceleration ############
    array = np.zeros(size)
    for read in sampler.reads:
        start, stop = read
        # UNCOMMENT TO ONLY EVERY DO SINGLE BP RESOLUTION
        #array[start:stop] +=1
        array[int(ceil(start/res)):int(ceil(stop/res))] += 1
    return array

def sp_sample_reads(sampler, size, i, j, res, jit=False):
    """ Sample reads from samplers using multiprocessing

    Args:
        sampler  (list): List of ReadSampler objects to sample from
        size     (int): Size of numpy array to create
        res      (float): resolution of numpy array
        start    (int): starting seed for random state
        p        (optional, int): number of processors

    Returns:
        arrays (np.array) final result of sampling
    """
    start = args.s+i+len(all_sim)
    array = sample_reads(sampler, size, np.random.RandomState(start+j), res, jit)
    return array

def sample_reads(sampler, size, prng, res=1, jit=False):
    """ Sample with replacement all reads from a sampler and map them to an
        array

    Args:
        sampler (class ReadSampler): object holding reads to sample
        size (int): size of numpy array to create
        prng (np.RandomState): random state to pull random numbers
                               from
        res (optional, int): resolution the numpy array is in
    """
    array = np.zeros(size)
    sampled_reads = sampler.pull_reads(sampler.total, prng)
    # # define heirarchy for using GPU 
    # threadsperblock = 32
    # blockspergrid = (array.size + (threadsperblock - 1)) // threadsperblock

    if jit:
        # threadsperblock = 32
        # blockspergrid = (sampled_reads.shape[0] + (threadsperblock - 1)) // threadsperblock
        # note that using cuda.jit here causes us to write directly to the array,
        #   without explicitly returning it from sum_coverage_from_sample function
        # sum_coverage_from_sample[blockspergrid, threadsperblock](sampled_reads, array, res)
        # device_array = cuda.to_device(array)
        # device_samples = cuda.to_device(sampled_reads)
        # device_res = cuda.to_device(res)
        array = numba_sum_coverage_from_sample(sampled_reads, array, res)
        # cuda_sum_coverage_from_sample[blockspergrid, threadsperblock](device_samples, device_array, res)
        # array = device_array.copy_to_host()
    else:
        array = sum_coverage_from_sample(sampled_reads, array, res)

    return array

def sum_coverage_from_sample(sampled_reads, array, res):
    """ Map sampled reads to an array

    Args:
        reads (): sampled read (start,end) positions
        array (): array to be populated with coverage
        res (float): resolution the array is in
    """

    for read in sampled_reads:
        start, stop = read
        res = float(res)
        x = int(ceil(start/res))
        y = int(ceil(stop/res))
        # UNCOMMENT TO ONLY EVERY DO SINGLE BP RESOLUTION
        #array[start:stop] +=1
        array[x:y] += 1

    return(array)

@numba.jit(nopython=True, parallel=True)
def numba_sum_coverage_from_sample(sampled_reads, array, res):
    """ Map sampled reads to an array

    Args:
        reads (): sampled read (start,end) positions
        array (): array to be populated with coverage
        res (float): resolution the array is in
    """
    # pos = cuda.grid(1)
    for i in range(sampled_reads.shape[0]):
        read = sampled_reads[i,:]
        # read = sampled_reads[pos,:]
        start, stop = read
        start = float(start)
        stop = float(stop)
        res = float(res)
        x = int(ceil(start/res))
        y = int(ceil(stop/res))
        # UNCOMMENT TO ONLY EVERY DO SINGLE BP RESOLUTION
        #array[start:stop] +=1
        for idx in range(x,y):
            array[idx] += 1
    return(array)

# @numba.cuda.jit
# def cuda_sum_coverage_from_sample(sampled_reads, array, res):
#     """ Map sampled reads to an array

#     Args:
#         reads (): sampled read (start,end) positions
#         array (): array to be populated with coverage
#         res (float): resolution the array is in
#     """
#     # pos = cuda.grid(1)
#     for i in range(sampled_reads.shape[0]):
#         read = sampled_reads[i,:]
#         # read = sampled_reads[pos,:]
#         start, stop = read
#         start = float(start)
#         stop = float(stop)
#         res = float(res)
#         x = int(ceil(start/res))
#         y = int(ceil(stop/res))
#         # UNCOMMENT TO ONLY EVERY DO SINGLE BP RESOLUTION
#         #array[start:stop] +=1
#         for idx in range(x,y):
#             array[idx] += 1

def log2_ratio(array1, array2):
    """ Take the median normalized log2 ratio of two arrays

    Args:
        array1 (np.array): numpy array holding extracted signal
        array2 (np.array): numpy array holding input signal
        cu (bool): whether to use cuda.jit to accelerate calculations

    Returns:
        ratio (np.array): median normalized log2 ratio
    """
    # lets normalize by the median
    arr1_median = float(np.median(array1))
    arr2_median = float(np.median(array2))
    ratio = (array1/arr1_median)/(array2/arr2_median)

    # if not cu:
    #     ratio = (array1/arr1_median)/(array2/arr2_median)
        # only nans should be np.log2(0/0) which should be 0
        #ratio[np.isnan(ratio)] = 0.0
    # else:
    #     threadsperblock = 32
    #     blockspergrid = (array1.size + (threadsperblock - 1)) // threadsperblock

    #     device_array1 = cuda.to_device(array1)
    #     device_arr1_med = cuda.to_device(arr1_median)
    #     device_array2 = cuda.to_device(array2)
    #     device_arr2_med = cuda.to_device(arr2_median)

    #     median_array1 = np.empty(array1.shape)
    #     device_median_array1 = cuda.to_device(median_array1)
    #     median_array2 = np.empty(array2.shape)
    #     device_median_array2 = cuda.to_device(median_array2)

    #     ratio = np.empty(array1.shape)
    #     device_ratio = cuda.to_device(ratio)
    #     cuda_ratio[blockspergrid, threadsperblock](device_array1, device_array2, 
    #                                                     arr1_median, arr2_median, 
    #                                                     device_median_array1, device_median_array2,
    #                                                     device_ratio)
    #     ratio = device_ratio.copy_to_host()
        
    log2_ratio = np.log2(ratio)

    return log2_ratio

# @cuda.jit
# def cuda_ratio(array1, array2, arr1_med, arr2_med,
#                     median_array1, median_array2, 
#                     ratio):

#     pos = cuda.grid(1)
#     if pos < ratio.shape[0]:
#         # calculate median for position pos of each array
#         median_array1[pos] = array1[pos]/arr1_med
#         median_array2[pos] = array2[pos]/arr2_med

#         # calculate ratio
#         ratio[pos] = median_array1[pos]/median_array2[pos]

    

def floored_sub(samp, control):
    """ Subtract control signal from sample signal. Only subtract if control is
        greater than 0.

    Args:
        samp (np.array): sample numpy array holding signal
        control (np.array): control numpy array holding signal

    Returns:
        sub (np.array): subtracted signal
    """
    # make copy as to not modify control array
    this_control = np.copy(control)
    this_control[this_control < 0] = 0
    return samp-this_control

def read_in_sampler(samplerfile):
    """ Take a saved numpy array file and load it back into a sampler

    Args:
        samplerfile (str): name of stored numpy array

    Returns:
        sampler (ReadSampler object): a sampler ready for sampling
    """
    sampler = ReadSampler()
    sampler.load_data(samplerfile)
    # sampler.sort_reads() # our reads are already sorted
    return sampler

def mp_sample_reads_helper(args):
    """ Helper function for sampling reads with multiprocessing
    """
    sampler, size, res, prng = args
    return sample_reads(sampler, size, prng, res)

def mp_sample_reads(samplers, size, res, start, p=1):
    """ Sample reads from samplers using multiprocessing

    Args:
        samplers (list): List of ReadSampler objects to sample from
        size     (int): Size of numpy array to create
        res      (int): resolution of numpy array
        start    (int): starting seed for random state
        p        (optional, int): number of processors

    Returns:
        arrays (np.array) final result of sampling
    """
    pool = mp.Pool(processes=p)
    arrays = pool.map(mp_sample_reads_helper,
            ((sampler, size, res, np.random.RandomState(start+x)) for x,sampler in enumerate(samplers)))
    pool.close()
    pool.join()
    return arrays

def mp_actual_reads_helper(args):
    """ Helper function for mapping reads with multiprocessing
    """
    sampler, size, res = args
    np.random.RandomState()
    return actual_reads(sampler, size, res)

def mp_actual_reads(samplers, size, res, p=1):
    """ Map reads from samplers using multiprocessing

    Args:
        samplers (list): List of ReadSampler objects to sample from
        size     (int): Size of numpy array to create
        res      (int): resolution of numpy array
        p        (optional, int): number of processors

    Returns:
        arrays (np.array) final result of mapping
    """
    pool = mp.Pool(processes=p)
    arrays = pool.map(mp_actual_reads_helper,
            ((sampler, size, res) for sampler in samplers))
    pool.close()
    pool.join()
    return arrays

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Given sampling objects,\
            goes through the entire chip-seq pipeline with x bootstraps.\
            outputs a numpy g by x matrix where g is the size of the genome\
            and x is the number of bootstrap replicates to do.")
    parser.add_argument("--genome_size", type=int, help="size of genome\
            controls how large to build the arrays", required=True)
    parser.add_argument('--out_prefix', help="Output prefix for final matrix.",
                        required=False, default='')
    parser.add_argument("--sample_name_luts", type=str, nargs="+",
            help="Paths to files containing the\
            lookup run information provided by Illumina. The\
            file will be used to match sample IDs\
            with sample names. Multiple file names should\
            be separated by spaces.", required=True)
    parser.add_argument('--ChIP_samps', type=str, nargs="+",
            help="Extracted read simulators for samples")
    parser.add_argument('--inp_samps', type=str, nargs="+",
            help="Input read simulators for samples")
    # parser.add_argument('--ChIP_conts', type=str, nargs="+",
    #         help="Extracted read simulators for controls")
    # parser.add_argument('--inp_conts', type=str, nargs="+",
    #         help="Input read simulators for controls")
    parser.add_argument('--num_replicates', type=int, default=5,
            help="Number of bootstrap replicates to do, default is 5")
    parser.add_argument('--identity', action="store_true",
            help="Don't sample, just use the data as is")
    parser.add_argument('-p', type=int, default=10,
            help="Number of processors to use, default is 5")
    parser.add_argument('-s', type=int, default=None,
            help="Seed for random number generation. Default is 1")
    parser.add_argument('--resolution', type=int, default=1,
                        help="resolution of data to map, default is 1bp")
    parser.add_argument('--save_summaries', type=float, default=None,
            help="Don't save all bootstrap replicates. Just save the summaries:\
                  minci, maxci, lev, mean, var. Specify the alpha level here")
    parser.add_argument('--numba', action="store_true",
            help="Adding this flag will enable jit compilation\
                  on the cpu to speed up calculations.")
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s',level=logging.INFO)
    ## TO DO:
    # allow functions to deal with more than two samples for each
    prng= np.random.RandomState(seed=args.s)
    array_size = int(floor(args.genome_size/args.resolution))

    num_ext_samp = len(args.ChIP_samps)
    num_inp_samp = len(args.inp_samps)
    # num_ext_cont = len(args.ChIP_conts)
    # num_inp_cont = len(args.inp_conts)

    if args.out_prefix == '':
        save_name_prefix = ''
    else:
        save_name_prefix = '{}_'.format(args.out_prefix)

    if args.save_summaries is not None:
        this_summary_stats_only_finite = lambda x: summary_stats_only_finite(x, alpha=args.save_summaries)
        this_summary_stats_nan_zeros = lambda x: summary_stats_nan_zeros(x, alpha=args.save_summaries)

    samp_final = np.zeros((array_size, args.num_replicates, num_ext_samp))
    # cont_final = np.zeros((array_size, args.num_replicates, num_ext_cont))

    pat = re.compile(r'Sample_\d+')
    # grab the sample names from the treatment sampler files
    samp_matches = [pat.search(s) for s in args.ChIP_samps]
    if None in samp_matches:
        print("Your sampler file names should begin with the sample id\
              Illumina gave you, i.e., Sample_NNNNNN,\
              where N is any integer. Exiting script.")
        sys.exit()
    else: samp_strings = [match.group() for match in samp_matches]
    # grab the sample names from the control sampler files
    # cont_matches = [pat.search(s) for s in args.ChIP_conts]
    # if None in cont_matches:
    #     print("Your sampler file names should begin with the sample id\
    #           Illumina gave you, i.e., Sample_NNNNNN,\
    #           where N is any integer. Exiting script.")
    #     sys.exit()
    # else: cont_strings = [match.group() for match in cont_matches]

    # read in sample info from illumina to look up descriptions from sample ids
    for i,sample_info_file_name in enumerate(args.sample_name_luts):
        if sample_info_file_name.endswith('.csv'):
            sample_info_tmp = (pd.read_csv(sample_info_file_name, header=18) >>
                               select(X.Sample_ID, X.Description))
        else:
            sample_info_tmp = (pd.read_csv(sample_info_file_name, sep='\t') >>
                               mutate(Sample_ID = 'Sample_' + X.SampleID.astype(str)) >>
                               select(X.Sample_ID, X.Description))
        if i == 0:
            sample_info = sample_info_tmp
        else:
            sample_info = sample_info.append(sample_info_tmp)

    samp_names = []
    for samp_id in samp_strings:
        sample_rows = sample_info[sample_info.Sample_ID == samp_id]
        samp_names.append(sample_rows['Description'].iloc[0]) # grab the first description for this Sample_ID

    # cont_names = []
    # for samp_id in cont_strings:
    #     cont_rows = sample_info[sample_info.Sample_ID == samp_id]
    #     cont_names.append(cont_rows['Description'].iloc[0])

    # samp_names = [os.path.basename(x).split(".")[0] for x in args.ChIP_samps]
    # cont_names = [os.path.basename(x).split(".")[0] for x in args.ChIP_conts]

    ## Read in all samplers
    all_sims = []
    all_sims.extend(args.ChIP_samps)
    all_sims.extend(args.inp_samps)
    # all_sims.extend(args.ChIP_conts)
    # all_sims.extend(args.inp_conts)
    if num_ext_samp != num_inp_samp: # or num_ext_cont != num_inp_cont:
        logging.error("Mismatch number of extracted and input samples Ext_samp: %s\
                Inp_samp: %s"%(num_ext_samp, num_inp_samp))#, num_ext_cont, num_inp_cont))

    all_sim = []
    for sampler in all_sims:
        logging.info("Reading in sampler {}".format(sampler))
        all_sim.append(read_in_sampler(sampler))

    # all_sim = [read_in_sampler(sampler) for sampler in all_sim]
    logging.info("Finished reading in samplers")

    ## sample reads
    for i in range(args.num_replicates):
        logging.info("Starting bootstrap replicate %s"%i)
        logging.info("Sampling reads %s"%i)
        begin = time.time()
        if args.numba:
            if args.identity:
                # arrays = mp_actual_reads(all_sim, array_size, args.resolution, args.p)
                ### UNCOMMENT TO GO BACK TO NOT USING MULTIPROCCESING
                arrays = [actual_reads(sampler, array_size, args.resolution) for sampler in all_sim]
            else:
                ### UNCOMMENT TO GO BACK TO NOT USING MULTIPROCCESING
                arrays = [sp_sample_reads(sampler, array_size, i, j, args.resolution, jit=args.numba) for j,sampler in enumerate(all_sim)]
                # arrays = mp_sample_reads(all_sim, array_size, args.resolution, args.s+i+len(all_sim), args.p)
        else:
            if args.identity:
                arrays = mp_actual_reads(all_sim, array_size, args.resolution, args.p)
                ### UNCOMMENT TO GO BACK TO NOT USING MULTIPROCCESING
                # arrays = [actual_reads(sampler, args.genome_size, args.resolution) for sampler in all_sim]
            else:
                ### UNCOMMENT TO GO BACK TO NOT USING MULTIPROCCESING
                # arrays = [sp_sample_reads(sampler, array_size, i, j, args.resolution, cu=args.numba) for j,sampler in enumerate(all_sim)]
                arrays = mp_sample_reads(all_sim, array_size, args.resolution, args.s+i+len(all_sim), args.p)
        finish = time.time()
        logging.info("Sampling bootstrap replicate {} took {} seconds".format(i, finish-begin))
        ## Calculate ratios
        logging.info("Calculating Ratios %s"%i)
        begin = time.time()
        for j, (ext_ind, inp_ind) in enumerate(zip(range(0, num_ext_samp),
                                                   range(num_ext_samp,
                                                         num_ext_samp+num_inp_samp))):
            samp_final[:,i,j] = log2_ratio(arrays[ext_ind], arrays[inp_ind])
        # for j, (ext_ind, inp_ind) in enumerate(zip(range(num_ext_samp + num_inp_samp,
        #                                                  num_ext_samp + num_inp_samp+num_ext_cont),
        #                                            range(num_ext_samp + num_inp_samp + num_ext_cont,
        #                                                  num_ext_samp + num_inp_samp + num_ext_cont + num_inp_cont))):
        #     cont_final[:,i,j]=log2_ratio(arrays[ext_ind], arrays[inp_ind])
        finish = time.time()
        logging.info("Calculating Ratios for replicate {} took {} seconds".format(i, finish-begin))


    # Write out final output
    if args.identity and args.save_summaries:
        # save sample summaries
        for j, name in enumerate(samp_names):
            save_name = "{}{}_actual".format(save_name_prefix,name)
            np.save(save_name, samp_final[:,:,j])
        # save control summaries
        # for j, name in enumerate(cont_names):
        #     np.save(args.out_prefix+"_%s_actual"%name, cont_final[:,:,j])
        # save each combination of sample - control summary
        # for j, samp_name in enumerate(samp_names):
        #     for k, cont_name in enumerate(cont_names):
        #         # note that floored sub MODIFIES the control array. Since we have already written the control array
        #         # that is not a big deal but be aware that this modification happens
        #         np.save(args.out_prefix+"_%s_Sub_%s_actual"%(samp_name,cont_name), floored_sub(samp_final[:,:,j], cont_final[:,:,k]))

    elif args.save_summaries:
        # ALL OF THIS IS HARDCODED QUICK CODING. Saves a lot of output files to be
        # used downstream
        logging.info("Calculating Summary Stats for finite values")
        begin = time.time()
        # save sample summaries
        for j, name in enumerate(samp_names):
            these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, samp_final[:,:,j])
            np.save("{}{}_mean".format(save_name_prefix,name), these_stats[:,0])
            np.save("{}{}_minci".format(save_name_prefix,name), these_stats[:,1])
            np.save("{}{}_maxci".format(save_name_prefix,name), these_stats[:,2])
            np.save("{}{}_var".format(save_name_prefix,name), these_stats[:,3])
            np.save("{}{}_lev".format(save_name_prefix,name), these_stats[:,4])
            np.save("{}{}_median".format(save_name_prefix,name), these_stats[:,5])
            np.save("{}{}_mad".format(save_name_prefix,name), these_stats[:,6])
            np.save("{}{}_num_inf".format(save_name_prefix,name), these_stats[:,7])
            np.save("{}{}_num_nan".format(save_name_prefix,name), these_stats[:,8])
        # save control summaries
        # for j, name in enumerate(cont_names):
        #     these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, cont_final[:,:,j])
        #     np.save("{}{}_mean".format(save_name_prefix,name), these_stats[:,0])
        #     np.save("{}{}_minci".format(save_name_prefix,name), these_stats[:,1])
        #     np.save("{}{}_maxci".format(save_name_prefix,name), these_stats[:,2])
        #     np.save("{}{}_var".format(save_name_prefix,name), these_stats[:,3])
        #     np.save("{}{}_lev".format(save_name_prefix,name), these_stats[:,4])
        #     np.save("{}{}_median".format(save_name_prefix,name), these_stats[:,5])
        #     np.save("{}{}_mad".format(save_name_prefix,name), these_stats[:,6])
        #     np.save("{}{}_num_inf".format(save_name_prefix,name), these_stats[:,7])
        #     np.save("{}{}_num_nan".format(save_name_prefix,name), these_stats[:,8])
        # # save each combination of sample - control summary
        # for j, samp_name in enumerate(samp_names):
        #     for k, cont_name in enumerate(cont_names):
        #         these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, floored_sub(samp_final[:,:,j], cont_final[:,:,k]))
        #         np.save(args.out_prefix+"_%s_Sub_%s_mean"%(samp_name,cont_name), these_stats[:,0])
        #         np.save(args.out_prefix+"_%s_Sub_%s_minci"%(samp_name,cont_name), these_stats[:,1])
        #         np.save(args.out_prefix+"_%s_Sub_%s_maxci"%(samp_name, cont_name), these_stats[:,2])
        #         np.save(args.out_prefix+"_%s_Sub_%s_var"%(samp_name,cont_name), these_stats[:,3])
        #         np.save(args.out_prefix+"_%s_Sub_%s_lev"%(samp_name,cont_name), these_stats[:,4])
        #         np.save(args.out_prefix+"_%s_Sub_%s_median"%(samp_name,cont_name), these_stats[:,5])
        #         np.save(args.out_prefix+"_%s_Sub_%s_mad"%(samp_name,cont_name), these_stats[:,6])
        #         np.save(args.out_prefix+"_%s_Sub_%s_num_inf"%(samp_name,cont_name), these_stats[:,7])
        #         np.save(args.out_prefix+"_%s_Sub_%s_num_nan"%(samp_name,cont_name), these_stats[:,8])
        # finish = time.time()
        # logging.info("Calculating Summary Stats took {} seconds".format(finish-begin))

#        logging.info("Calculating Summary Stats for nans as zeros")
#        # save sample summaries
#        for j, name in enumerate(samp_names):
#            these_stats = np.apply_along_axis(this_summary_stats_nan_zeros, 1, samp_final[:,:,j])
#            np.save("{}{}_mean".format(save_name_prefix,name), these_stats[:,0])
#            np.save("{}{}_minci".format(save_name_prefix,name), these_stats[:,1])
#            np.save("{}{}_maxci".format(save_name_prefix,name), these_stats[:,2])
#            np.save("{}{}_var".format(save_name_prefix,name), these_stats[:,3])
#            np.save("{}{}_lev".format(save_name_prefix,name), these_stats[:,4])
#            np.save("{}{}_median".format(save_name_prefix,name), these_stats[:,5])
#            np.save("{}{}_mad".format(save_name_prefix,name), these_stats[:,6])
#        # save control summaries
#        for j, name in enumerate(cont_names):
#            these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, cont_final[:,:,j])
#            np.save("{}{}_mean".format(save_name_prefix,name), these_stats[:,0])
#            np.save("{}{}_minci".format(save_name_prefix,name), these_stats[:,1])
#            np.save("{}{}_maxci".format(save_name_prefix,name), these_stats[:,2])
#            np.save("{}{}_var".format(save_name_prefix,name), these_stats[:,3])
#            np.save("{}{}_lev".format(save_name_prefix,name), these_stats[:,4])
#            np.save("{}{}_median".format(save_name_prefix,name), these_stats[:,5])
#            np.save("{}{}_mad".format(save_name_prefix,name), these_stats[:,6])
#        # save each combination of sample - control summary
#        for j, samp_name in enumerate(samp_names):
#            for k, cont_name in enumerate(cont_names):
#                these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, floored_sub(samp_final[:,:,j], cont_final[:,:,k]))
#                np.save(args.out_prefix+"_%s_Sub_%s_mean"%(samp_name,cont_name), these_stats[:,0])
#                np.save(args.out_prefix+"_%s_Sub_%s_minci"%(samp_name,cont_name), these_stats[:,1])
#                np.save(args.out_prefix+"_%s_Sub_%s_maxci"%(samp_name, cont_name), these_stats[:,2])
#                np.save(args.out_prefix+"_%s_Sub_%s_var"%(samp_name,cont_name), these_stats[:,3])
#                np.save(args.out_prefix+"_%s_Sub_%s_lev"%(samp_name,cont_name), these_stats[:,4])
#                np.save(args.out_prefix+"_%s_Sub_%s_median"%(samp_name,cont_name), these_stats[:,5])
#                np.save(args.out_prefix+"_%s_Sub_%s_mad"%(samp_name,cont_name), these_stats[:,6])

    else:
        # write final array
        logging.info("Writing final array")
        np.save("{}samp".format(save_name_prefix), samp_final)
        # np.save(args.out_prefix+"_cont", cont_final)
