################################################################################
# Script to parse an alignment for reads that overlap locations we want to
# eliminate. Writes a text file of read names to eliminate
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

import sam_utils
import sys
import gfftools
from collections import defaultdict


## FUNCTIONS INVOLVED IN DETERMINING READ OVERLAP FOR ELIMINATION ##
def interval_overlap(range1, range2):
    """Function to determine if two intervals overlap. These are python
    intervals. I.e. [start, end)
    
    Inputs: range1 - tuple (x y) where x <= y
            range2 - tuple (x y) where x <= y

    Modifies: nothing
    Returns: True if there is an overlap, False otherwise

    Tests:
    Non-overlapping but next to each other
    >>> interval_overlap((1,2), (2,3))
    False

    Non-overlapping but next to each other (flipped)
    >>> interval_overlap((2,3), (1,2))
    False

    Overlapping by one
    >>> interval_overlap((1,3), (2,3))
    True

    Same exact range
    >>> interval_overlap((1,3), (1,3))
    True

    range2 is entirely within range1
    >>> interval_overlap((1,4), (2,3))
    True

    range1 is entirely within range2
    >>> interval_overlap((1,2), (0,3))
    True
    """
    return(range1[0] < range2[1] and range2[0] < range1[1])

def read_overlap(locs, start_end):
    """Function to determine whether a start_end tuple overlaps with
    locations of interest.

    Inputs: locs      - list. of locations of the following structure 
                        [(start, stop), ... (start, stop)]
            start_end - tuple. of start and end of the fragment or read.
                        i.e. (start, end)
    Modifies: nothing
    Returns: True if the read overlaps, False if it does not and raises
             an error if sense_to_ref could not be determined
    """
    # If the locs dictionary is empty, then no locs need to be checked
    if locs == None:
        return False
    # otherwise for each start end passed to the function
    for read_interval in start_end:
        # loop through each testing start and testing end
        for test_interval in locs:
            # test if the read is outside the testing locations
            if interval_overlap(read_interval, test_interval):
                return True
    return False

def read_list_of_locations(locations_file, chrom_names):
    """Takes in a .txt file with tab delimitted locations by format:
    chrom start end strand. Ex. chr1 100 200 + would be a location
    that started at 100 ended at 200 and was on the plus strand
    
    Inputs: locations_file - string. .txt file  with tab delimited locations, or
                             a .gff file with features to remove, or a .gtf
                             file with features to remove.
            chrom_names    - list. of chromosome names
    Modifies: nothing
    Returns: A dictionary of dictionaries where the internal dictionary has
    keys of '+' or '-' holding the locations for either and the external
    dictionary has keys of each chromosome name.
    """
    # create default dictionary of type dictionary to hold everything
    locs = defaultdict(dict)
    # create an empty list for each chromosome and strand type
    for chrom in chrom_names:
        locs[chrom]['+'] = []
        locs[chrom]['-'] = []
        locs[chrom]['.'] = []
    # open up and read the input file
    with open(locations_file, 'r') as locs_file:
        # go through each line and split it
        for line in locs_file:
            # if the file is a txt file than parse it appropriately
            if locations_file.endswith('.txt'):
                line_arr = line.rstrip().split()
                # chrom should be first field
                chrom = line_arr[0]
                # start should be the min of the next two fields
                start = min(int(line_arr[1]), int(line_arr[2]))
                # end should be the max of the same fields
                end = max(int(line_arr[1]), int(line_arr[2]))
                # strand should be the last field
                strand = line_arr[3]
            # if the field is a .gff or .gtf file, treat it differently
            elif locations_file.endswith('.gff') or \
                 locations_file.endswith('.gtf'):
                line_arr = line.rstrip().split("\t")
                # first field will be the chromosome
                chrom = line_arr[0]
                # start will be the min of the 3rd and 4th fields
                start = min(int(line_arr[3]), int(line_arr[4]))
                # end will be the max of the same fields
                end = max(int(line_arr[3]), int(line_arr[4]))
                # strand will be in field 5 (all of these assume 0-based field)
                strand = line_arr[6]
            else:
                # if not one of those file types, raise a value error
                raise ValueError('File must end in .gtf, .gff, or .txt')
            # finally append it to the appropriate place in the dictionary
            locs[chrom][strand].append((start, end))
    # sort the resulting
    return(locs)
## END FUNCTIONS INVOLVED IN DETERMINING READ OVERLAP FOR ELIMINATION ##

    
def main():
    prefix = sys.argv[1]
    gff = sys.argv[2]
    locs = read_list_of_locations(gff, ["gi|48994873|gb|U00096.2|mod|ATCC.47076|"])
    read_set = set()
    with open(prefix+"_bad_reads.txt", mode="w") as f:
        for line in sys.stdin:
            read = sam_utils.SamAlignment(line)
            try:
                if read.sense_to_ref(True, "R2"):
                    strand = "+"
                else:
                    strand = "-"
            except:
                continue
            start,end,gaps = read.start_end_gaps(True)
            if read_overlap(locs[read.RNAME][strand], [(start, end)]):
               f.write(read.QNAME+"\n")
            else:
                continue

main()

        
