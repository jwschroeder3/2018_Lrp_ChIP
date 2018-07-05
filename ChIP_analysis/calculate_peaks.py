################################################################################
# Calculate binding regions through combination of IDR, tech and bio filters
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

# base
import sys
import logging
import argparse
# widely available
import numpy as np
import scipy.stats as stats
# custom
import peak
import convert_multiple_to_average_robustz as rz

def neg_logp_to_neg_logq(array):
    """ Convert -log_10(pvalues) to -log_10(qvalues) using the benjamini-hochberg
    procedure. This gives answers identical to R's p.adjust (if it dealt with log
    transformed p values)

    Be careful with large arrays. Two additional array.size are 
    created internally, one is returned.

    Args:
        array (np.array) - contains -log_10(pvalues), unsorted.
    Returns:
        qvalues (np.array) - -log_10(qvalues) in the same order as p values
    """
    # get an index for the sorted array, This is sorted from largest p value
    # to smallest (since -log(p) higher means smaller p)
    idx = np.argsort(array)
    asize = array.size
    qvalues = np.zeros(asize)
    M = np.log10(asize)
    qval_prev = 0
    for i, idx_val in enumerate(idx):
        qval = array[idx_val] + (np.log10(asize-i) - M)
        qval = max(qval, qval_prev)
#        print idx_val,array[idx_val], np.log(asize-i), M, qval_this,qval,np.e**(-1*qval)
        qvalues[idx_val] = qval
        qval_prev = qval
    return qvalues

def call_peaks(array, consolidate = 0):
    """ Take a logical array and find the start and stop of ones across
    the array. Consolidate any peaks that are within consolidate.

    TODO: Deal with peaks that go over the end of the array

    Args:
        array (np.array) - logical array to call peaks from.
        consolidate (int) - number of bins to consolidate peaks over
    Returns:
        peak_indices (list of lists) - a list of [start, stop] indices in the
                                        array
    """
    # first find all the places where there is a change from 1 to 0 or vice
    # versa, here we pad the incoming array with a zero on each end, then take
    # the difference along the array and finally take the absolute value to flip
    # the places where the difference was 0 - 1
    changes = np.abs(np.diff(np.concatenate(([0], array.view(np.int8), [0]))))
    # changes is now a logical array, I then find all the indices where changes
    # happen and reshape them into an ndarray of start, stop locations
    start_stops = np.where(changes)[0].reshape(-1, 2)
    if start_stops.size == 0:
        logging.warning("No bp was considered to be within a peak.")
        return []

    # now lets consolidate any peaks that are within consolidate
    consolidate_peaks = [[start_stops[0][0], start_stops[0][1]]]
    consolidated_peaks = 0
    removed_peaks = 0
    for start, stop in start_stops[1:]:
        if start - consolidate_peaks[-1][1] < consolidate:
            consolidate_peaks[-1][1] = stop
            consolidated_peaks += 1
        else:
            if stop-start > consolidate:
                consolidate_peaks.append([start, stop])
            else:
                removed_peaks += 1
    logging.info("Consolidated %s peaks within %s bps"%(consolidated_peaks, consolidate))
    logging.info("Removed %s peaks < %s bps"%(removed_peaks, consolidate))
    return consolidate_peaks

def convert_technical_to_pvalue(signal, mad):
    """ Take signal and the mad from bootstrapping and generate a pvalue

    Args:
        signal (np.array): input signal
        mad    (np.array): median absolute deviation from bootstrapping

    Returns:
        pvalue (np.array): a pvalue for every location
    """
    zscore = rz.calculate_zscore(signal, mad)
    return 1-stats.norm.cdf(zscore)

def convert_biological_to_pvalue(signal):
    """ Take robustz signal and convert to pvalue

    Args:
        signal (np.array): input robustz

    Returns:
        pvalue (np.array): pvalues from robustZ
    """
    robustz = rz.calculate_robust_z(signal)
    return 1-stats.norm.cdf(robustz)

def convert_p_to_q(pvalues):
    """ Convert pvalues to qvalues using Benjamini Hochberg procedure

    Args:
        pvalues (np.array): input pvalues

    Returns:
        qvalues (np.array)
    """
    qvalues = np.ones(pvalues.size)
    these_pvalues = -np.log10(pvalues)
    # truncate any infinite pvalues to the last finite pvalue
    these_pvalues[these_pvalues == np.inf] = np.max(these_pvalues[np.isfinite(these_pvalues)])
    qvalues = neg_logp_to_neg_logq(these_pvalues)
    return qvalues

def bio_and_tech_filter(bio_qs, tech_qs):
    """ Apply biological and technical filters

    Args:
        bio_qs (np.array): biological qvalues
        tech_qs (np.array): technical qvalues

    Returns:
        temp_peaks (np.array) boolean array where TRUE indicates peak region
    """
    temp_peaks = np.zeros(bio_qs[0].size)
    this_rep = np.zeros(bio_qs[0].size)
    for bioq, techq in zip(bio_qs, tech_qs):
        this_rep += 1*np.logical_and(bioq > -np.log10(args.bioalpha), techq > -np.log10(args.techalpha))
        temp_peaks = np.logical_or(temp_peaks, this_rep > 0)
    print "Fractions passing bio and tech filter"
    print np.sum(this_rep == 1)/(bio_qs[0].size + 0.0)
    print np.sum(this_rep == 2)/(bio_qs[0].size + 0.0)
    print np.sum(this_rep == 3)/(bio_qs[0].size + 0.0)
    print np.sum(this_rep == 4)/(bio_qs[0].size + 0.0)
    return temp_peaks

def idr_filter(temp_peaks, idrs):
    """ Apply IDR filter

    Args:
        temp_peaks (np.array): output from bio and tech filter function
        idrs (np.array): results from IDR calculations

    Returns:
        temp_peaks (np.array) boolean array where TRUE indicates peak region
    """
    new_peaks = temp_peaks.copy()
    for i,idr in enumerate(idrs):
        # only consider locations that are already a peak
        idr[temp_peaks] = convert_p_to_q(idr[temp_peaks])
        # set all other q values to 1
        idr[~temp_peaks] = -np.log10(1)
        new_peaks = np.logical_and(new_peaks, idr > -np.log10(args.idralpha))
    return new_peaks

def find_max_value(start, stop, vallist, ret_index=False):
    max_val = 0
    for val in vallist:
        this_max = np.nanmax(val[start:stop])
        if this_max > max_val:
            max_val = this_max
            if ret_index:
                maxidx = np.where(val[start:stop] == max_val)[0][0]
    if ret_index:
        return max_val, maxidx
    else:
        return max_val

def find_min_value(start, stop, vallist, ret_index=False):
    max_val = 0
    for val in vallist:
        this_max = np.nanmax(val[start:stop])
        if this_max < max_val:
            max_val = this_max
            if ret_index:
                maxidx = np.where(val[start:stop] == max_val)[0][0]
    if ret_index:
        return max_val, maxidx
    else:
        return max_val

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate peaks from bio, tech, and IDR filters")
    parser.add_argument('--log2ratios', type=str, nargs="+", help="log2 ratio data")
    parser.add_argument('--mad', type=str, nargs="+", help="mad from bootstrapping for each log2 ratio")
    parser.add_argument('--idr', type=str, nargs="+", help="idr between samples")
    parser.add_argument('--bioalpha', type=float, help="alpha for log2 ratio", default=0.001)
    parser.add_argument('--techalpha', type=float, help="alpha for technical reproducibility", default=0.05)
    parser.add_argument('--idralpha', type=float, help="alpha for idr", default=0.05)
    parser.add_argument('--resolution', type=float, help="resolution of data", default=1)
    parser.add_argument('--bins', type=int, help="bins to consolidate over", default=30)
    parser.add_argument('--outpre', type=str, help="output prefix")
    args = parser.parse_args()

    actual_sigs = [np.load(x) for x in args.log2ratios]
    actual_sigs = [sig.flatten() for sig in actual_sigs]
    mad_sigs = [np.load(x) for x in args.mad]
    mad_sigs = [sig.flatten() for sig in mad_sigs]
    idrs = [np.load(x) for x in args.idr]
    idrs = [sig.flatten() for sig in idrs]
    bio_ps = [convert_biological_to_pvalue(sig) for sig in actual_sigs]
    tech_ps = [convert_technical_to_pvalue(sig, mad) for sig, mad in zip(actual_sigs, mad_sigs)]
    bio_qs = [convert_p_to_q(bio_p) for bio_p in bio_ps]
    tech_qs= [convert_p_to_q(tech_p) for tech_p in tech_ps]
    bio_tech_peaks = bio_and_tech_filter(bio_qs, tech_qs)
    peaks = call_peaks(idr_filter(bio_tech_peaks, idrs), consolidate=args.bins)
    averagerobustz = np.mean([rz.calculate_robust_z(sig) for sig in actual_sigs], 0)


    # Make peaks
    all_peaks = peak.PeakList()
    for i, entry in enumerate(peaks):
        start = entry[0]*args.resolution
        end = entry[1]*args.resolution
        max_value, max_idx = find_max_value(entry[0], entry[1], [averagerobustz], True)
        max_valbio = find_max_value(entry[0], entry[1], bio_qs)
        max_validr = find_max_value(entry[0], entry[1], idrs, False)
        max_valtech = find_max_value(entry[0], entry[1], tech_qs, False)
        this_peak = peak.Peak("gi|48994873|gb|U00096.2|mod|ATCC.47076|", start=int(start), end=int(end), 
                              name = "%s_%s"%(args.outpre,i),
                              score = max_value,
                              signalval = max_valtech,
                              pval = max_valbio,
                              qval = max_validr, peak = int(max_idx*args.resolution))
        if end-start >= args.bins*args.resolution:
            all_peaks.add_Peak(this_peak)
        else:
            continue
    # Write peaks
    all_peaks.write_narrowPeak_file(args.outpre + ".narrowPeak")
