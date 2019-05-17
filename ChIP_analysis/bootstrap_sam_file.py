################################################################################
# bootstrap_sam_file
# A script for sampling reads from a sam/bam file
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

## Python base imports ##
from __future__ import division, print_function
import argparse
import sys
import random
import logging
from math import floor, ceil
# Wildely available
import numpy as np
# custom module
import sam_utils

## TO DO:
# add support for multiple chromosomed organisms
# add support for strandedness

class ReadSampler(object):
    """Class to hold and sample from processed reads. Reads are stored internally
       as [0,end) coordinates.

    Attributes:
        reads (list): holds (start, end) tuples of reads. Converted to numpy
                      array when sampling
        total (int): number of reads in the list
        sampling (boolean): True if ready for sampling, False otherwise
    """
    def __init__(self):
        self.reads = []
        self.total = 0
        self.sampling = False

    def add_read(self, new_read):
        """ Method to add a new read to the sampler.

        If sampler is in sampling mode, then converts self.reads to a list
        and sets the sampler flag back to False.


        Args:
            new_read (tuple): new (start, end) read
        """
        if self.sampling:
            self.convert_to_list()
        self.reads.append(new_read)
        self.total+=1

    def add_reads(self, new_reads):
        """ Method to add a multiple new reads to the sampler.

        If sampler is in sampling mode, then converts self.reads to a list
        and sets the sampler flag back to False.

        Args:
            new_reads (list): new list of (start, end) reads
        """
        if self.sampling:
            self.convert_to_list()
        self.reads.extend(new_reads)

    def convert_to_array(self):
        """ Method to convert reads to a numpy array for efficient sampling
        """
        self.reads = np.asarray(self.reads, dtype="int64")
        self.sampling=True

    def convert_to_list(self):
        """ Method to convert reads to a list for efficient addition of reads
        """
        self.reads = list(self.reads)
        self.sampling = False

    def pull_read(self, prng):
        """ Method to sample a single read from the sampler.

        Will convert sampler to sampling mode if not already.

        Args:
            prng (np.RandomState object): random state to use for sampling
        """
        if not self.sampling:
            self.convert_to_array()
        index = prng.random.randint(0, self.total)
        return self.reads[index, :]

    def pull_reads(self, n, prng):
        """ Method to sample a multiple reads from the sampler.

        Will convert sampler to sampling mode if not already.

        Args:
            n (int): number of reads to sample
            prng (np.RandomState object): random state to use for sampling
        """
        if not self.sampling:
            self.convert_to_array()
        index = prng.randint(0, self.total, size=n)
        index = np.sort(index)
        return self.reads[index,:]

    def save_data(self, f):
        """ Method to save sampler reads for later use

        Will convert sampler to sampling mode if not already. Only saves
        the reads attribute as a binary numpy array

        Args:
            f (fhandle): a file handle to save to
        """
        if not self.sampling:
            self.convert_to_array()
        np.save(f, self.reads)

    def sort_reads(self):
        """ Method to sort reads by starting location

        Will convert sampler to sampling mode if not already.
        """
        if not self.sampling:
            self.convert_to_array()
        self.reads = self.reads[self.reads[:,0].argsort()]
    def load_data(self, f):
        """ Method to load data from a saved sampler object

        Will overwrite any attributes that were previously in object
        """
        self.sampling = True
        self.reads = np.load(f)
        self.total = self.reads.shape[0]



def merge(intervals):
    """ Merge several individual intervals into one (start, stop) interval

    Used to take two paired end reads and convert them to a single start stop
    interval. Not appropriate for RNA-seq reads, only DNA reads (i.e. ChIP-Seq).

    Args:
        intervals (list): list of [start, stop] intervals to merge

    Returns:
        merged (list): list containing final [start, stop] interval
    """
    intervals.sort(key=lambda x: x[0])
    # take the first interval
    merged = [intervals[0]]
    # loop through all the intervals
    for this_interval in intervals:
        if this_interval[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], this_interval[1]))
        else:
            merged.append(this_interval)
    return merged

def get_paired_blocks(r1, r2):
    """ Take two individual read alignments and return the bases that the
    reads overlap in a single [start, stop] interval

    Args:
        r1 (sam_utils.SamAlignment): read 1 of a pair (order doesn't matter)
        r2 (sam_utils.SamAlignment): read 2 of a pair (order doesn't matter)

    Returns:
        total_blocks: list containing final [start, stop] interval

    Raises:
        ValueError if reads aren't a proper pair
        RuntimeError if a gapped read is found
    """
    if r1.TLEN > 0 and r2.TLEN < 0:
        left = r1.get_aligned_blocks()
        right = r2.get_aligned_blocks()
    elif r1.TLEN < 0 and r2.TLEN > 0:
        left = r2.get_aligned_blocks()
        right = r1.get_aligned_blocks()
    elif r1.POS == r2.POS:
        left = r1.get_aligned_blocks()
        right = r2.get_aligned_blocks()
    else:
        raise ValueError("Pair not consistent %s %s"%(r1.QNAME, r2.QNAME))
    total_blocks = []
    if left[-1][1] < right[0][0]:
        total_blocks.append((left[-1][1], right[0][0]))
    total_blocks.extend(left)
    total_blocks.extend(right)
    total_blocks = merge(total_blocks)
    if len(total_blocks) > 1:
        raise RuntimeError("Gapped read found %s %s %s"%(r1.QNAME, r2.QNAME, str(total_blocks)))
    return total_blocks[0]

def create_read_list(samfile):
    """ Read in a samfile and convert it to a list of reads for the sampler
    object

    This function is for single end reads only. Skips any reads that are
    gapped reads

    Args:
        samfile (fhandle): an open filehandle for reading of a samfile

    Returns:
        read_sampler(obj(ReadSampler)): final read sampler for the samfile

    """
    read_sampler = ReadSampler()
    for line in samfile:
        line = sam_utils.SamAlignment(line)
        vals = line.get_aligned_blocks()
        if len(vals) > 1:
            logging.info("Skipping gapped read %s %s"%(line.QNAME, str(vals)))
        read_sampler.add_read(vals[0])
    return read_sampler

def create_read_list_paired(samfile):
    """ Read in a samfile and convert it to a list of reads for the sampler
    object.

    This function is for paired end reads only. Skips any reads that are gapped
    reads or are not properly paired. Assumes samfile is sorted by readname and
    only one alignment per pair is present in the file. If these assumptions
    are not met than this function will yield nonsense.

    Args:
        samfile (fhandle): an open filehandle for reading of a samfile

    Returns:
        read_sampler(obj(ReadSampler)): final read sampler for the samfile

    Raises:
        ValueError: If pair readnames don't match. Not considered a failsafe
        for catching violations of assumptions above but should catch most
        mistakes.
    """
    read_sampler = ReadSampler()
    while True:
        line1 = samfile.readline()
        line2 = samfile.readline()
        if not line2:
            break
        line1 = sam_utils.SamAlignment(line1)
        line2 = sam_utils.SamAlignment(line2)
        if line1.QNAME != line2.QNAME:
            raise ValueError("Unpaired read or read with more than one pair\
                              encountered. Check your input file. File must\
                              be sorted by read name, every read must have\
                              a single pair and each pair must have one\
                              mapping. %s %s"%(line1.QNAME, line2.QNAME))
        try:
            read_sampler.add_read(get_paired_blocks(line1,line2))
        except ValueError as err:
            logging.error("Skipping pair %s"%err)
        except RuntimeError as err:
            logging.error("Skipping pair %s"%err)
    return read_sampler

def map_read(array, read, res=1.0):
    """ Take in a [start, stop] read and map it to an array.

    Read must always be in single basepair coordinates [0,end). Array can
    be any resolution desired through control of res parameter. Modifies
    array in place.

    Args:
        array (1d np.array): array to store coverage in
        read (list): Single [start, stop) read to map
        res (optional, float): resolution the numpy array is in

    """
    start, stop = read
    array[int(ceil(start/res)):int(ceil(stop/res))] += 1

def sample(read_sampler, n, array, res=1.0, prng = np.random.RandomState()):
    """ Sample reads with replacement from a sampler and map them to an array

    Modifies array in place

    Args:
        read_sampler (class ReadSampler): object holding reads to sample
        n (int): number of reads to sample
        array (1d np.array): array to store coverage in
        res (optional, float): resolution the numpy array is in
        prng (optional, np.RandomState): random state to pull random numbers
                                         from
    """
    for read in read_sampler.pull_reads(n, prng):
        map_read(array, read, res)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='commands', dest='command')

    # parse
    parse_parser = subparsers.add_parser('parse', help="create a sampling\
            object from a sam file")
    parse_parser.add_argument('samfile', help="Input samfile, this tool does no\
            filtering and will consider every line in the file. Accepts input\
            from stdin if '-' is specified here.")
    parse_parser.add_argument('outpre', help="output prefix to np.save the\
            sampling object data that is created")
    parse_parser.add_argument('--paired', action="store_true", help="Consider\
            the sam file as paired. If this flag is specified then the sam file\
            MUST be pre-filtered to have only ONE alignment per pair. Further,\
            there must be NO unpaired reads in the file and the reads must be\
            sorted by read name.")

    # sample
    sample_parser = subparsers.add_parser('sample', help="sample coverage from a\
            sampling object.")
    sample_parser.add_argument('samplerfile', help="output file from using parse\
            on the samfile of interest")
    sample_parser.add_argument('outpre', help="output file to np.save the\
            numpy array that is created")
    sample_parser.add_argument('array_size',type=int, help="length of genome")
    sample_parser.add_argument('--num_samples', type=int, default=1,
    help="number of full samples to pull from the sampler, default is 1")
    sample_parser.add_argument('--num_reads', type=int, default=None,
    help="number of reads to pull for each sample. Default is the size of\
            sampling object.")
    sample_parser.add_argument('--identity', action="store_true",
            help="write an array of the actual sample without sampling, ignores\
                    num_reads, num_samples and seed options")
    sample_parser.add_argument('--resolution', type=int, default=1,
            help="only report every x basepairs, default=1")
    sample_parser.add_argument('--seed', type=int, default=1234,
            help="psuedo-random number generator seed, default=1234")
    args = parser.parse_args()


    if args.command == "parse":
        if args.samfile == "-":
            f = sys.stdin
        else:
            f = open(args.samfile, mode="r")
        if args.paired:
            sampler = create_read_list_paired(f)
        else:
            sampler = create_read_list(f)
        f.close()
        sampler.sort_reads()
        sampler.save_data(args.outpre)
    elif args.command == "sample":
        prng = np.random.RandomState(args.seed)
        array = np.zeros((int(floor(args.array_size/args.resolution)), args.num_samples))
        sampler = ReadSampler()
        sampler.load_data(args.samplerfile)
        if args.identity:
            for read in sampler.reads:
                map_read(array, read, args.resolution)
            np.save(args.outpre, array)
        else:
            for i in xrange(args.num_samples):
                if args.num_reads:
                    num_reads = args.num_reads
                else:
                    num_reads = sampler.total
                sample(sampler, num_reads, array[:,i], args.resolution, prng)
            np.save(args.outpre, array)
