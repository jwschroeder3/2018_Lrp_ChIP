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

# base python
import argparse
import sys
import logging
import os
import multiprocessing as mp
from math import floor, ceil

# Widely available
import numpy as np
from scipy import stats

# custom
from bootstrap_sam_file import ReadSampler
from bootstrap_helper_funcs import credible_interval, least_extreme_value


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
        return(these_stats[0], these_stats[1], these_stats[2], var, this_lev, median, mad, num_inf, num_nan)

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


def actual_reads(sampler, size, res=1.0):
    """ Take all reads from a sampler and map them to an array

    Args:
        sampler (class ReadSampler): object holding reads to sample
        size (int): size of numpy array to create
        res (optional, float): resolution the numpy array is in
    """
    array = np.zeros(size)
    for read in sampler.reads:
        start, stop = read
        # UNCOMMENT TO ONLY EVERY DO SINGLE BP RESOLUTION
        #array[start:stop] +=1
        array[int(ceil(start/res)):int(ceil(stop/res))] += 1
    return array

def sample_reads(sampler, size, prng, res=1.0):
    """ Sample with replacement all reads from a sampler and map them to an
        array

    Args:
        sampler (class ReadSampler): object holding reads to sample
        size (int): size of numpy array to create
        prng (np.RandomState): random state to pull random numbers 
                               from
        res (optional, float): resolution the numpy array is in
    """
    array = np.zeros(size)
    for read in sampler.pull_reads(sampler.total, prng):
        start, stop = read
        # UNCOMMENT TO ONLY EVERY DO SINGLE BP RESOLUTION
        #array[start:stop] +=1
        array[int(ceil(start/res)):int(ceil(stop/res))] += 1
    return array

def log2_ratio(array1, array2):
    """ Take the median normalized log2 ratio of two arrays

    Args:
        array1 (np.array): numpy array holding extracted signal
        array2 (np.array): numpy array holding input signal

    Returns:
        ratio (np.array): median normalized log2 ratio
    """
    # lets normalize by the median
    ratio = (array1/(np.median(array1)+0.0))/(array2/(np.median(array2)+0.0))
    ratio = np.log2(ratio)
    # only nans should be np.log2(0/0) which should be 0
    #ratio[np.isnan(ratio)] = 0.0
    return ratio
 
def floored_sub(samp, control):
    """ Subtract control signal from sample signal. Only subtract is control is
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
    sampler.sort_reads()
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
        res      (float): resolution of numpy array
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
        res      (float): resolution of numpy array
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
    parser.add_argument("genome_size", type=int, help="size of genome\
            controls how large to build the arrays")
    parser.add_argument('outpre', help="Output prefix for final matrix.")
    parser.add_argument('--ext_samps', type=str, nargs="+", 
            help="Extracted read simulators for samples")
    parser.add_argument('--inp_samps', type=str, nargs="+", 
            help="Input read simulators for samples")
    parser.add_argument('--ext_conts', type=str, nargs="+", 
            help="Extracted read simulators for controls")
    parser.add_argument('--inp_conts', type=str, nargs="+", 
            help="Input read simulators for controls")
    parser.add_argument('--num_replicates', type=int, default=5, 
            help="Number of bootstrap replicates to do, default is 5")
    parser.add_argument('--identity', action="store_true",
            help="Don't sample, just use the data as is")
    parser.add_argument('-p', type=int, default=5, 
            help="Number of processors to use, default is 5")
    parser.add_argument('-s', type=int, default=None, 
            help="Seed for random number generation. Default is 1")
    parser.add_argument('--resolution', type=float, default=1.0,
                        help="resolution of data to map, default is 1bp")
    parser.add_argument('--save_summaries', type=float, default=None,
            help="Don't save all bootstrap replicates. Just save the summaries:\
                  minci, maxci, lev, mean, var. Specify the alpha level here")
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s',level=logging.INFO)    
    ## TO DO:
    # allow functions to deal with more than two samples for each
    prng= np.random.RandomState(seed=args.s)
    array_size = int(floor(args.genome_size/args.resolution))

    num_ext_samp = len(args.ext_samps)
    num_inp_samp = len(args.inp_samps)
    num_ext_cont = len(args.ext_conts)
    num_inp_cont = len(args.inp_conts)

    if args.save_summaries is not None:
        this_summary_stats_only_finite = lambda x: summary_stats_only_finite(x, alpha=args.save_summaries)
        this_summary_stats_nan_zeros = lambda x: summary_stats_nan_zeros(x, alpha=args.save_summaries)
        
    samp_final = np.zeros((array_size, args.num_replicates, num_ext_samp))
    cont_final = np.zeros((array_size, args.num_replicates, num_ext_cont))

    samp_names = [os.path.basename(x).split(".")[0] for x in args.ext_samps]
    cont_names = [os.path.basename(x).split(".")[0] for x in args.ext_conts]

    ## Read in all samplers
    all_sim = []
    all_sim.extend(args.ext_samps)
    all_sim.extend(args.inp_samps)
    all_sim.extend(args.ext_conts)
    all_sim.extend(args.inp_conts)
    if num_ext_samp != num_inp_samp or num_ext_cont != num_inp_cont:
        logging.error("Mismatch number of extracted and input samples Ext_samp: %s\
                Inp_samp: %s, Ext_cont: %s, Inp_cont: %s"%(num_ext_samp, num_inp_samp, num_ext_cont, num_inp_cont))

    logging.info("Reading in samplers")
    all_sim = [read_in_sampler(sampler) for sampler in all_sim]
    logging.info("Finished reading in samplers")

    ## sample reads
    for i in xrange(args.num_replicates):
        logging.info("Starting bootstrap replicate %s"%i)
        logging.info("Sampling reads %s"%i)
        if args.identity:
            arrays = mp_actual_reads(all_sim, array_size, args.resolution, args.p)
            ### UNCOMMENT TO GO BACK TO NOT USING MULTIPROCCESING
#            arrays = [actual_reads(sampler, args.genome_size, args.resolution) for sampler in all_sim]
        else:
            ### UNCOMMENT TO GO BACK TO NOT USING MULTIPROCCESING
            #arrays = [sample_reads(sampler, args.genome_size, args.resolution) for sampler in all_sim]
            arrays = mp_sample_reads(all_sim, array_size, args.resolution, args.s+i+len(all_sim), args.p)
        ## Calculate ratios
        logging.info("Calculating Ratios %s"%i)
        for j, (ext_ind, inp_ind) in enumerate(zip(range(0, num_ext_samp), 
                                                   range(num_ext_samp, 
                                                         num_ext_samp+num_inp_samp))):
            samp_final[:,i,j] = log2_ratio(arrays[ext_ind], arrays[inp_ind])
        for j, (ext_ind, inp_ind) in enumerate(zip(range(num_ext_samp + num_inp_samp, 
                                                         num_ext_samp + num_inp_samp+num_ext_cont), 
                                                   range(num_ext_samp + num_inp_samp + num_ext_cont, 
                                                         num_ext_samp + num_inp_samp + num_ext_cont + num_inp_cont))):
            cont_final[:,i,j]=log2_ratio(arrays[ext_ind], arrays[inp_ind])


    # Write out final output
    if args.identity and args.save_summaries:
        for j, name in enumerate(samp_names):
            np.save(args.outpre+"_%s_actual"%name, samp_final[:,:,j])
        # save control summaries
        for j, name in enumerate(cont_names):
            np.save(args.outpre+"_%s_actual"%name, cont_final[:,:,j])
        # save each combination of sample - control summary
        for j, samp_name in enumerate(samp_names):
            for k, cont_name in enumerate(cont_names):
                # note that floored sub MODIFIES the control array. Since we have already written the control array
                # that is not a big deal but be aware that this modification happens
                np.save(args.outpre+"_%s_Sub_%s_actual"%(samp_name,cont_name), floored_sub(samp_final[:,:,j], cont_final[:,:,k]))

    elif args.save_summaries:
        # ALL OF THIS IS HARDCODED QUICK CODING. Saves a lot of output files to be
        # used downstream
        logging.info("Calculating Summary Stats for finite values")
        # save sample summaries
        for j, name in enumerate(samp_names):
            these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, samp_final[:,:,j])
            np.save(args.outpre+"_%s_mean"%name, these_stats[:,0])
            np.save(args.outpre+"_%s_minci"%name, these_stats[:,1])
            np.save(args.outpre+"_%s_maxci"%name, these_stats[:,2])
            np.save(args.outpre+"_%s_var"%name, these_stats[:,3])
            np.save(args.outpre+"_%s_lev"%name, these_stats[:,4])
            np.save(args.outpre+"_%s_median"%name, these_stats[:,5])
            np.save(args.outpre+"_%s_mad"%name, these_stats[:,6])
            np.save(args.outpre+"_%s_num_inf"%name, these_stats[:,7])
            np.save(args.outpre+"_%s_num_nan"%name, these_stats[:,8])
        # save control summaries
        for j, name in enumerate(cont_names):
            these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, cont_final[:,:,j])
            np.save(args.outpre+"_%s_mean"%name, these_stats[:,0])
            np.save(args.outpre+"_%s_minci"%name, these_stats[:,1])
            np.save(args.outpre+"_%s_maxci"%name, these_stats[:,2])
            np.save(args.outpre+"_%s_var"%name, these_stats[:,3])
            np.save(args.outpre+"_%s_lev"%name, these_stats[:,4])
            np.save(args.outpre+"_%s_median"%name, these_stats[:,5])
            np.save(args.outpre+"_%s_mad"%name, these_stats[:,6])
            np.save(args.outpre+"_%s_num_inf"%name, these_stats[:,7])
            np.save(args.outpre+"_%s_num_nan"%name, these_stats[:,8])
        # save each combination of sample - control summary
        for j, samp_name in enumerate(samp_names):
            for k, cont_name in enumerate(cont_names):
                these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, floored_sub(samp_final[:,:,j], cont_final[:,:,k]))
                np.save(args.outpre+"_%s_Sub_%s_mean"%(samp_name,cont_name), these_stats[:,0])
                np.save(args.outpre+"_%s_Sub_%s_minci"%(samp_name,cont_name), these_stats[:,1])
                np.save(args.outpre+"_%s_Sub_%s_maxci"%(samp_name, cont_name), these_stats[:,2])
                np.save(args.outpre+"_%s_Sub_%s_var"%(samp_name,cont_name), these_stats[:,3])
                np.save(args.outpre+"_%s_Sub_%s_lev"%(samp_name,cont_name), these_stats[:,4])
                np.save(args.outpre+"_%s_Sub_%s_median"%(samp_name,cont_name), these_stats[:,5])
                np.save(args.outpre+"_%s_Sub_%s_mad"%(samp_name,cont_name), these_stats[:,6])
                np.save(args.outpre+"_%s_Sub_%s_num_inf"%(samp_name,cont_name), these_stats[:,7])
                np.save(args.outpre+"_%s_Sub_%s_num_nan"%(samp_name,cont_name), these_stats[:,8])

#        logging.info("Calculating Summary Stats for nans as zeros")
#        # save sample summaries
#        for j, name in enumerate(samp_names):
#            these_stats = np.apply_along_axis(this_summary_stats_nan_zeros, 1, samp_final[:,:,j])
#            np.save(args.outpre+"_%s_mean"%name, these_stats[:,0])
#            np.save(args.outpre+"_%s_minci"%name, these_stats[:,1])
#            np.save(args.outpre+"_%s_maxci"%name, these_stats[:,2])
#            np.save(args.outpre+"_%s_var"%name, these_stats[:,3])
#            np.save(args.outpre+"_%s_lev"%name, these_stats[:,4])
#            np.save(args.outpre+"_%s_median"%name, these_stats[:,5])
#            np.save(args.outpre+"_%s_mad"%name, these_stats[:,6])
#        # save control summaries
#        for j, name in enumerate(cont_names):
#            these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, cont_final[:,:,j])
#            np.save(args.outpre+"_%s_mean"%name, these_stats[:,0])
#            np.save(args.outpre+"_%s_minci"%name, these_stats[:,1])
#            np.save(args.outpre+"_%s_maxci"%name, these_stats[:,2])
#            np.save(args.outpre+"_%s_var"%name, these_stats[:,3])
#            np.save(args.outpre+"_%s_lev"%name, these_stats[:,4])
#            np.save(args.outpre+"_%s_median"%name, these_stats[:,5])
#            np.save(args.outpre+"_%s_mad"%name, these_stats[:,6])
#        # save each combination of sample - control summary
#        for j, samp_name in enumerate(samp_names):
#            for k, cont_name in enumerate(cont_names):
#                these_stats = np.apply_along_axis(this_summary_stats_only_finite, 1, floored_sub(samp_final[:,:,j], cont_final[:,:,k]))
#                np.save(args.outpre+"_%s_Sub_%s_mean"%(samp_name,cont_name), these_stats[:,0])
#                np.save(args.outpre+"_%s_Sub_%s_minci"%(samp_name,cont_name), these_stats[:,1])
#                np.save(args.outpre+"_%s_Sub_%s_maxci"%(samp_name, cont_name), these_stats[:,2])
#                np.save(args.outpre+"_%s_Sub_%s_var"%(samp_name,cont_name), these_stats[:,3])
#                np.save(args.outpre+"_%s_Sub_%s_lev"%(samp_name,cont_name), these_stats[:,4])
#                np.save(args.outpre+"_%s_Sub_%s_median"%(samp_name,cont_name), these_stats[:,5])
#                np.save(args.outpre+"_%s_Sub_%s_mad"%(samp_name,cont_name), these_stats[:,6])

    else:
        # write final array 
        logging.info("Writing final array")
        np.save(args.outpre+"_samp", samp_final)
        np.save(args.outpre+"_cont", cont_final)
