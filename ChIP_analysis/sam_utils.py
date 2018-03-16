################################################################################
# sam_utils module
# A module for storing a manipulating sam/bam alignments
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
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.#
################################################################################

## Python base imports ##
import re
import itertools
import logging

## Widely available ##
import numpy as np
import pandas as pd
import statsmodels.sandbox.stats.multicomp as multicomp

## LOGGING SET UP ##
FORMAT='%(levelname)s %(asctime)s: %(message)s'
logging.basicConfig(format=FORMAT, level=logging.INFO)

## HELPER FUNCTIONS ##
def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks. Taken from
    itertools documentation

    Test:
    >>> out = grouper('ABCDEFG', 3, 'x')
    >>> for group in out:
    ...    print group
    ('A', 'B', 'C')
    ('D', 'E', 'F')
    ('G', 'x', 'x')
    """
    args = [iter(iterable)] * n
    return itertools.izip_longest(fillvalue=fillvalue, *args)

def complement(sequence):
    """Complement a nucleotide sequence
    >>> complement("AGTC")
    'TCAG'
    >>> complement("AGNT")
    'TCNA'
    >>> complement("AG-T")
    'TC-A'
    """
    # create a dictionary to act as a mapper
    comp_dict = {'A': 'T', 'G':'C', 'C':'G', 'T': 'A', 'N':'N', '-':'-'}
    # turn the sequence into a list
    sequence = list(sequence)
    # remap it to the compelmentary sequence using the mapping dict
    sequence = [comp_dict[base] for base in sequence]
    # join the new complemented sequence list into a string
    sequence = ''.join(sequence)
    return sequence
## END HELPER FUNCTIONS ##

## CLASSES ##
class SamAlignment:
    """ Class for holding a single line out of a SAM file. Along with methods
    to manipulate and extract data from the alignment.
    """
        
    def __init__(self, line, sep = '\t'):
        """ Class initalizer. Takes in a single line from a Sam alignment and
        stores each and every field as a class attribute.

        Input: line - str Raw Sam input file line.
        Modifies: self.QNAME - str Name of the read
                  self.FLAG  - int Bitwise flag of read information
                  self.RNAME - str Name of the reference chromosome
                  self.POS   - int 0-based left most mapping of first base
                  self.MAPQ  - int mapping quality
                  self.CIGAR - str contains information on mapping with gaps
                  self.RNEXT - str name of reference chromosome of read pair
                  self.PNEXT - int POS of other read in pair
                  self.TLEN  - int number of bases from left to right in pair
                  self.SEQ   - str contains sequence of the read
                  self.QUAL  - str contains ASCII base quality information
                  self.OPT   - str optional fields specific to each aligner
        Returns: None

        Tests:

        Test reading in a typical SAM entry
        >>> line = "testQ.1.testR.20.30.3M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.QNAME == "testQ"
        True
        >>> read.FLAG == 1
        True
        >>> read.RNAME == "testR"
        True
        >>> read.POS == 19
        True
        >>> read.MAPQ == 30
        True
        >>> read.CIGAR == "3M"
        True
        >>> read.RNEXT == "testR"
        True
        >>> read.PNEXT == 24
        True
        >>> read.TLEN == 9
        True
        >>> read.SEQ == "AAA"
        True
        >>> read.QUAL == "((("
        True
        >>> read.OPT["NM"] == 0
        True
        """
        linearr = line.rstrip().split(sep)
        # All information comes from SAM Specificationv1 18 Nov 2015
        # This is the name of the read, two reads from the same
        # template share the same QNAME. * indicates the information
        # is unknown.
        self.QNAME = linearr[0]
        # This is a combination of bitwise flags. Each bit is as follows:
        # Bit    Hex    Description
        # 1      0x1    template having multiple segments
        # 2      0x2    each segment properly aligned
        # 4      0x4    segment unmapped
        # 8      0x8    next segment in template unmapped
        # 16     0x10   SEQ being reverse complemented
        # 32     0x20   SEQ of the next segment reverse complemented
        # 64     0x40   the first segment in the template
        # 128    0x80   the last segment in the template
        # 256    0x100  secondary alignment
        # 512    0x200  not passing filters(platform/vendor)
        # 1024   0x400  PCR or optical duplicate
        # 2048   0x800  supplementary alignment
        self.FLAG = int(linearr[1])
        # Reference sequence name. Would be name of the chromosome in an
        # organism with multiple chromosomes
        self.RNAME = linearr[2]
        # 0-based left most mapping of the first base. -1 if unmapped. If
        # -1 no assumptions can be made about RNAME and CIGAR
        self.POS = int(linearr[3]) - 1
        # Mapping quality. Equal to -10log10(P(mapping position is wrong)).
        # Rounded to nearest integer. 255 means it is unavailable
        self.MAPQ = int(linearr[4])
        # CIGAR string, set to * if unavailable.
        # Value  Description
        #   M     alignment match
        #   I     insertion to the reference
        #   D     deletion from the reference
        #   N     skipped region from the reference
        #   S     soft clipping (clipped sequences present in SEQ)
        #   H     hard clipping (clipped sequences NOT in SEQ)
        #   P     padding (silent deletion from padded reference)
        #   =     sequence match
        #   X     sequence mismatch
        # H can only be present as first or last operation
        # S may only have H operations between them and the ends
        # N operation represents an intron for mRNA-Genome
        # Sum of the lengths shall equal the length of SEQ
        self.CIGAR = linearr[5]
        # Reference sequence name for other read in template, set to
        # = if they are identical and * if unavailable.
        self.RNEXT = linearr[6]
        # POS of the other read in template. -1 when information is unavail.
        # If -1 no assumptions can be made on RNEXT and bit 0x20
        # Subtract one to keep the value in 0-based coordinates
        self.PNEXT = int(linearr[7]) - 1
        # If all mapped to same reference. Then equal number of bases
        # from leftmost mapped base to righmost mapped base. Set to 0
        # for single-segment template or when information is unavail.
        self.TLEN = int(linearr[8])
        # Sequence for the segment
        self.SEQ = linearr[9]
        # ASCII quality score plus 33 (Phred 33)
        self.QUAL = linearr[10]
        # Optional fields. See SAM spec for more details
        self.OPT = self.parse_opt_fields(linearr[11:])

        ## Attributes that act as caches
        self.cigar_tuples = None
        self.aligned_blocks = None
        self.aligned_seq = None
        self.aligned_phred = None
        self.aligned_locs = None
        self.aligned_reference = None
        self.aligned_muts = None

    def parse_opt_fields(self, fields_list):
        """ Parses each field from the optional fields in the 
        SAM file. 

        Inputs: List. Contains all the optional fields in a sam file
        Modifies: nothing
        Returns: dictionary. Has field tags as keys and correctly converted
                 values. 'H' and 'B' field types are converted to strings
                 for now. Could be implemented to arrays in the future
        Tests:
        >>> line = "testQ.4.testR.20.30.3M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> ans = read.parse_opt_fields(["NM:i:0", "JK:f:0", "HT:Z:0", "QR:A:0",
        ...                              "AB:H:0", "CD:B:0"])
        >>> type(ans["NM"]) is int
        True
        >>> type(ans["JK"]) is float
        True
        >>> type(ans["HT"]) is str
        True
        >>> type(ans["QR"]) is str
        True
        >>> type(ans["AB"]) is str
        True
        >>> type(ans["CD"]) is str
        True
        """
        # Technically, H and B are arrays, but I don't have that implemented
        # right now, could be implemented in the future
        # Have a dictionary to hold the functions to convert each type of
        # field to the type it indicates
        d_type_func = {'A': str, 'i': int, 'f':float, 'Z':str, 
                       'H': str, 'B': str}
        # initialize a dictionary to hold the parsed fields
        field_dict = {}
        # loop through the additional fields
        for field in fields_list:
            # get name, type and value
            name, d_type, value  = field.split(':', 3)
            # convert value to appropriate type and put in dictionary for
            # holding fields
            field_dict[name] = d_type_func[d_type](value)
        # return this dictionary
        return(field_dict)

    def quality_filter(self, good_flags=[], bad_flags=[], mapq=None, stats=None):
        """ Function to filter a read based on flag and mapq.  Checks the
        self.FLAG field for the 4 (unmapped) and 2 (properly aligned) bits.
        Inputs allow additional, filters for determining if a read is good or
        not.

        Input: good_flags - flags that are required to be present in self.FLAG 
                            in order to pass the filter
               bad_flags - flags that are required to be absent from self.FLAG 
                           in order to pass the filter
               [mapq]  - int. If positive, only return true if self.MAPQ is
                         >= mapq. If negative only return true if self.MAPQ
                         is < mapq. Default is None, which will not filter by
                         mapq value.
               [stats] - AlignmentStats object. Allows for gathering stats
                         on how many reads are filtered for what issues.

        Modifies: stats
        Returns: True if the read passes filters and False otherwise. 

        Tests:

        Test a bad flag
        >>> line = "testQ.4.testR.20.30.3M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.quality_filter(bad_flags = [4])
        False
        >>> read.quality_filter()
        True

        Test a good flag
        >>> read.FLAG = 2
        >>> read.quality_filter(good_flags = [2])
        True
        >>> read.FLAG = 0
        >>> read.quality_filter(good_flags = [2])
        False

        Test multiple flags
        >>> read.FLAG = 2 | 1 | 256
        >>> read.quality_filter(good_flags = [2, 1])
        True
        >>> read.quality_filter(good_flags = [2, 1], bad_flags = [256])
        False
        >>> read.FLAG = 2
        >>> read.quality_filter(good_flags = [2, 1], bad_flags = [256])
        False

        Test must be above a mapq value
        >>> read.FLAG = 0
        >>> read.quality_filter(mapq=30)
        True
        >>> read.quality_filter(mapq=31)
        False

        Test must be below a mapq value
        >>> read.quality_filter(mapq=-30)
        False
        >>> read.quality_filter(mapq=-31)
        True

        """
        good_al = True
        # 1. Check bad flags
        for flag in bad_flags:
            if flag & self.FLAG:
                good_al = False
                if stats:
                    stats.bad_flag_reads[flag] = \
                            stats.bad_flag_reads.get(flag, 0) + 1
        # 2. Check good flags
        for flag in good_flags:
            if flag & self.FLAG:
                pass
            else:
                good_al = False
                if stats:
                    stats.missed_good_flag_reads[flag] = \
                            stats.missed_good_flag_reads.get(flag, 0) + 1

        # 3. Check for mapq conditions if required
        # If no mapq then we definitely passed the condition.
        if mapq == None:
            mapq_condition = True
        elif mapq >= 0:
            # if the mapq is greater than 0 then we are looking for reads
            # ABOVE the mapq score specified
            mapq_condition = self.MAPQ >= mapq 
        elif mapq < 0:
            # if the mapq is less than 0 then we are looking for reads
            # BELOW the mapq score specified
            mapq_condition = self.MAPQ < abs(mapq) 
        else:
            raise ValueError("mapq must be a positive or negative integer "\
                             "value.") 
        good_al = good_al and mapq_condition
        if stats and not mapq_condition:
            stats.total_mapq += 1
        # finally return whether the alignment was good or not
        return bool(good_al)
    
    def is_gapped(self):
        """ Determine if the alignment has any gaps in it as determined by the
        CIGAR field.

        Inputs: nothing
        Modifies: nothing
        Returns: True if 'N' is in the CIGAR field. 
                 False if not. Raises error if the CIGAR field = "*"
        Tests:
        >>> line = "testQ.1.testR.20.30.3M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.is_gapped()
        False
        >>> line = "testQ.1.testR.20.30.1M1N1M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.is_gapped()
        True
        >>> line = "testQ.1.testR.20.30.*.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.is_gapped()#doctest: +ELLIPSIS
        Traceback (most recent call last):
        RuntimeError: CIGAR field is empty,...for read testQ
        """

        # if N is present then there is a gap, if not then there isnt
        if('N' in self.CIGAR):
            return(True)
        # however if the cigar field doesn't exist, then we can't tell
        elif self.CIGAR == '*':
            raise RuntimeError("CIGAR field is empty, can not determine if " +
                               "the read is gapped or not " +
                               "for read %s"%self.QNAME)
        else:
            return(False)

    def is_rc(self, paired):

        """ Function to determine if the read is reverse complemented to the
        reference. If paired, forces each read of the pair to have opposite
        orientations, otherwise it raises an error

        Inputs: paired - boolean. if true, checks reads mate to see if it is 
                         in a consistent orientation. If not, raises error.
                         if false, just checks if read is rc to reference.
        Modifies: nothing
        Returns: True if reverse complemented. False otherwise

        Tests:

        Test a single end read
        >>> line = "testQ.1.testR.20.30.3M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.FLAG = 16
        >>> read.is_rc(False)
        True
        >>> read.FLAG = 1
        >>> read.is_rc(False)
        False

        Test a paired end read both in same direction
        >>> read.FLAG = 16 | 32
        >>> read.is_rc(True) #doctest: +ELLIPSIS
        Traceback (most recent call last):
        RuntimeError: Read mate ...
        >>> read.FLAG = 1
        >>> read.is_rc(True) #doctest: +ELLIPSIS
        Traceback (most recent call last):
        RuntimeError: Read mate ...

        Test a paired end read that is not reverse complemented
        >>> read.FLAG = 32
        >>> read.is_rc(True)
        False
        >>> read.FLAG = 16
        >>> read.is_rc(True)
        True

        """

        # check if the read is reverse complement to the reference
        this_is_rc = eval(bin(self.FLAG)) & 0x10
        # check if the read's pair is reverse complement to the reference
        other_is_rc = eval(bin(self.FLAG)) & 0x20
        
        # if paired, force the read's pair to have the opposite orientation
        # of the read. Raise error if not
        if paired:
            if ( this_is_rc and (not(other_is_rc)) ):
                return(True)
            elif ( other_is_rc and (not(this_is_rc)) ):
                return(False)
            else:
                raise RuntimeError("Read mate is in the same orientation as "\
                                   "read. Could not determine whether fragment "\
                                   "was reverse complemented or not")
        else:
        # if the read isn't paired, then we can just return whether it was
        # reverse complemented to the reference or not
            return(bool(this_is_rc))

    def sense_to_ref(self, paired, library):
        """ Function to determine whether the ORIGINAL RNA molecule was sense
        to the reference or antisense to the reference. Only implemented for 
        pairs when each read is in a different orientation
        
        Inputs: paired  - boolean. If True than additional checking will be done
                          to make sure the paired read is consistent with this
                          read. If False, then it is considered an individual
                          read.
                library - str. ["R1" | "R2" | "unstranded"]. R1 indicates that 
                          the first
                          read sequenced (*_R1.fasta) is sense to the original
                          RNA strand. R2 indicates that the 2nd read sequenced
                          (*_R2.fasta) is sense to the original RNA strand.
                          unstranded library types follow the the R1 convention.
        Modifies: nothing
        Returns: True if the original read sequence is the same as the reference
                 False if the original read sequence is not.
        Tests:

        Test a single end read
        >>> line = "testQ.1.testR.20.30.3M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.FLAG = 16
        >>> read.sense_to_ref(False, "R1")
        False
        >>> read.sense_to_ref(False, "R2")
        True
        >>> read.FLAG = 1
        >>> read.sense_to_ref(False, "R1")
        True
        >>> read.sense_to_ref(False, "R2")
        False

        Test a paired end read both in same direction
        >>> read.FLAG = 16 | 32
        >>> read.sense_to_ref(True, "R1") #doctest: +ELLIPSIS
        Traceback (most recent call last):
        RuntimeError: Read mate ...
        >>> read.FLAG = 1
        >>> read.sense_to_ref(True, "R1") #doctest: +ELLIPSIS
        Traceback (most recent call last):
        RuntimeError: Read mate ...

        Test a paired end read 
        >>> read.FLAG = 32 | 64
        >>> read.sense_to_ref(True, "R1")
        True
        >>> read.sense_to_ref(True, "R2")
        False
        >>> read.FLAG = 32 | 128
        >>> read.sense_to_ref(True, "R1")
        False
        >>> read.sense_to_ref(True, "R2")
        True
        >>> read.FLAG = 16 | 64
        >>> read.sense_to_ref(True, "R1")
        False
        >>> read.sense_to_ref(True, "R2")
        True
        >>> read.FLAG = 16 | 128
        >>> read.sense_to_ref(True, "R1")
        True
        >>> read.sense_to_ref(True, "R2")
        False

        """
        # if the read is paired, we have to consider the orientation of
        # the reads
        if paired: 
            # first segment in the template
            is_first = eval(bin(self.FLAG)) & 0x40
            # last segment in the template
            is_second = eval(bin(self.FLAG)) & 0x80
            # first try to figure out with this read is reverse complemented
            try:
                is_rc = self.is_rc(True)
            # raise the error if this can't be determined
            except RuntimeError:
                raise 
            # determine what orientation that entire segment is in
            if ( (is_first and not(is_rc)) or (is_second and is_rc) ):
                if library == "R1" or library == "unstranded":
                    return(True)
                elif library == "R2":
                    return(False)
                else:
                    raise ValueError("library must be R1, R2, or unstranded")
            elif ( (is_second and not(is_rc)) or (is_first and is_rc) ):
                if library == "R1" or library == "unstranded":
                    return(False)
                elif library == "R2":
                    return(True)
                else:
                    raise ValueError("library must be R1, R2, or unstranded")
            else:
                raise RuntimeError("Read pairs are not in the "+
                                   " the proper orientation, could " +
                                   "not determine whether read is sense "+
                                   " or not: %s."%self.QNAME)
        else:
            # if the read is not paired, then we determine if it is sense to
            # the reference based soley on its reverse complementarity
            is_rc = self.is_rc(False)
            if is_rc:
                if library == "R1" or library == "unstranded":
                    return(False)
                elif library == "R2":
                    return(True)
                else:
                    raise ValueError("library must be R1, R2, or unstranded")
            elif not(is_rc):
                if library == "R1" or library == "unstranded":
                    return(True)
                elif library == "R2":
                    return(False)
                else:
                    raise ValueError("library must be R1, R2, or unstranded")
            else:
                raise RuntimeError("Couldn't determine if single end read "+
                                   "was reverse complemented or not.")

    def get_cigar_tuples(self):

        """Function to parse the cigar field and turn it into tuples
        using regular expressions. Will also cache the return value into
        self.cigar_tuples so that the calculation only has to be performed
        once, when handling the record.

        Input: nothing
        Modifies: self.cigar_tuples
        Returns: a tuple of the cigar field broken up by separators.

        Tests:

        All M test
        >>> line = "testQ.1.testR.20.30.7M.testR.25.9.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_cigar_tuples()
        [(7, 'M')]

        More complex test
        >>> read.CIGAR = "3M4D2M1I1M" 
        >>> read.cigar_tuples = None
        >>> read.get_cigar_tuples()
        [(3, 'M'), (4, 'D'), (2, 'M'), (1, 'I'), (1, 'M')]

        """
        # if the value returned by this function has already been cached, then
        # just return that value and skip the calculation
        if self.cigar_tuples:
            return self.cigar_tuples

        # otherwise initalize a list to store each of the tuples
        cigar_tuples = []
        # split using regular expressions. The parenthesis in the re gives us
        # the seperators as well
        cigar_list = re.split('([MIDNSHP=X])', self.CIGAR)
        
        # loop through by twos using itertools grouper recipe from the 
        # python itertools documentation.
        # The cigar string always starts with a number and ends in a char,
        # so we cut the list short by one since the last value by twos will
        # end up being a None.
        for number, char in grouper(cigar_list[:-1], 2):
            cigar_tuples.append((int(number), char))
        # set the value into the cache so the calculation doesn't have to be
        # performed again
        self.cigar_tuples = cigar_tuples
        # return the now cached cigar_tuples value
        return self.cigar_tuples


    def get_aligned_blocks(self):

        """ Function to take the cigar field and determine the locations
        where there is continuous mapping coverage. 
        
        Inputs: nothing
        Modifies: nothing
        Returns: a list of (start, end) locations where continuous mapping 
                 coverage occurs. Continuous coverage includes locations where 
                 there is continuous 'M', '=', 'D', or 'X' in the CIGAR field. 
                 Breaks occur at 'S', 'I', or 'N' in the CIGAR field.
        Tests:
        All matching test, read.POS is at 1, cigar is 7M
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_aligned_blocks()
        [(1, 8)]

        Test internal deletions
        >>> read.CIGAR = "2M3D5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 11)]

        Test internal introns
        >>> read.CIGAR = "2M3N5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 3), (6, 11)]

        Test internal insertions
        >>> read.CIGAR = "2M3I2M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 5)]

        Test soft clipping on either end
        >>> read.CIGAR = "2S5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 6)]
        >>> read.CIGAR = "5M2S"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 6)]
        """
        if self.aligned_blocks:
            return self.aligned_blocks
        # initialize variables needed in the loop
        end = None
        tuple_list = []
        start = self.POS
        # parse the cigar string into tuples
        cigar_tuples = self.get_cigar_tuples()
       
        # go through each cigar number and character
        for val, char in cigar_tuples:
            # if it is an alignment match of any sort,
            # just add the value to determine the end. If the end
            # had not been determined previously, add it. otherwise,
            # just add the value to the previous end.
            if char in ['M', '=', 'X', 'D']:
                if end is not None:
                    end += val 
                else:
                    end = start + val
            # If there is an intron, go ahead and append
            # the previous (start, end) pair to the list and reset the
            # next end to being unknown. Additionally, push the next start
            # to the end of the intron. If an end had not been 
            # determined yet, do not add the (start, end) pair to the list
            # as this is a continuing insertion of some sort.
            elif char == 'N':
                if end is not None:
                    tuple_list.append((start, end))
                    start = end + val
                    end = None
                else:
                    start += val
            elif char in ['S', 'I']:
                continue
        # Finally, once all the way through the cigar field. Append the last
        # (start, stop) pair as long as it doesn't end in an insertion. If
        # it ends in an insertion then don't append the last pair.
        if end is not None:
            tuple_list.append((start, end))

        # last add the final list to the cache in the class and return
        # the value.
        self.aligned_blocks = tuple_list
        return self.aligned_blocks

    def get_aligned_seq_and_phred(self):

        """This function uses the CIGAR field of the read to determine
        what the sequence that aligns to the reference looks like after
        removing any insertions and accounting for deletions. It likewise 
        returns the corresponding quality scores associated with the sequence.

        Inputs: nothing
        Modifies: nothing
        Returns: seq       - read sequence with the insertions and deletions
                             removed
                 phred     - phred scores corresponding to each base in seq
        Tests:
        
        Test an all matching sequence
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_aligned_seq_and_phred()
        ('AGTCGCT', '!#%()+,')

        Test internal deletions
        >>> read.CIGAR = "2M3D5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_seq = None
        >>> read.aligned_phred = None
        >>> read.get_aligned_seq_and_phred()
        ('AG---TCGCT', '!#---%()+,')

        Test internal introns
        >>> read.CIGAR = "2M3N5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_seq = None
        >>> read.aligned_phred = None
        >>> read.get_aligned_seq_and_phred()
        ('AGTCGCT', '!#%()+,')

        Test internal insertions
        >>> read.CIGAR = "2M3I2M"
        >>> read.cigar_tuples = None
        >>> read.aligned_seq = None
        >>> read.aligned_phred = None
        >>> read.get_aligned_seq_and_phred()
        ('AGCT', '!#+,')

        Test soft clipping on either end
        >>> read.CIGAR = "2S5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_seq = None
        >>> read.aligned_phred = None
        >>> read.get_aligned_seq_and_phred()
        ('TCGCT', '%()+,')
        >>> read.CIGAR = "5M2S"
        >>> read.cigar_tuples = None
        >>> read.aligned_seq = None
        >>> read.aligned_phred = None
        >>> read.get_aligned_seq_and_phred()
        ('AGTCG', '!#%()')

        """

        # check to see if the cache is populated, if so return the cache
        if self.aligned_seq and self.aligned_phred:
            return(self.aligned_seq, self.aligned_phred)
        # first parse the cigar field to get the tuples
        cigar_tuples = self.get_cigar_tuples()
        expanded_cigar = ''
        # expand out the tuples to make a big cigar string
        for val, char in cigar_tuples:
            expanded_cigar += val*char
        
        # remove the elements that would not be included in the SEQ field
        expanded_cigar = expanded_cigar.replace('N', '')
        expanded_cigar = expanded_cigar.replace('H', '')
    
        #Initialize a new seq variable
        new_seq = ''
        new_phred = ''
        i = 0
        # loop through the expanded cigar string to create the new string
        for cigar in expanded_cigar:
            # If the alignment is a match add it to the new seq and increment
            # the index
            if cigar in ['M', '=', 'X']:
                new_seq += self.SEQ[i]
                new_phred += self.QUAL[i]
                i += 1
            # if there is a deletion, DON'T
            # increment the counter
            elif cigar == 'D':
                new_seq += '-'
                new_phred += '-'
                continue
            # If there is an insertion, increment the counter WITHOUT
            # adding anything to the new seq
            elif cigar in ['I', 'S', 'N']:
                i+= 1
                continue
            else:
                # I still don't know how to handle P in the CIGAR field
                raise ValueError("Can't handle cigar value %s"%cigar)

        # add the newly created seq and phred to the cache
        self.aligned_seq = new_seq
        self.aligned_phred = new_phred

        # return the cache
        return self.aligned_seq, self.aligned_phred

    def get_gap_blocks(self):
        """ Function to take the cigar field and determine the locations
        where there is a gap in coverage
        
        Inputs: nothing
        Modifies: nothing
        Returns: a list of (start, end) locations where gaps occur in mapping 
                 coverage.
        Tests:

        All matching test, read.POS is at 1, cigar is 7M
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_gap_blocks()
        []

        Test internal deletions
        >>> read.CIGAR = "2M3D5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_gap_blocks()
        []

        Test internal introns
        >>> read.CIGAR = "2M3N5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_gap_blocks()
        [(3, 6)]

        Test internal insertions
        >>> read.CIGAR = "2M3I2M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_gap_blocks()
        []

        Test soft clipping on either end
        >>> read.CIGAR = "2S5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_gap_blocks()
        []
        >>> read.CIGAR = "5M2S"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_gap_blocks()
        []
        """
        # first get the aligned blocks
        aligned_blocks = self.get_aligned_blocks()
        # next take the list of tuples of blocks and split it into two
        # lists, one with the starts and one with the ends
        start, end = zip(*aligned_blocks)
        # take all the starts from 0 to the second to last start
        new_start = end[0:-1]
        # take all the ends from the second to the last end
        new_end = start[1:]
        # zip them together into a new list of tuples and return
        return zip(new_start, new_end)
        
    def get_aligned_locs(self):

        """ Function to get the corresponding locations with
        the aligned sequence of the read.

        Inputs: nothing
        Modifies: nothing
        Returns: seq       - read sequence with the insertions and deletions
                             removed
                 phred     - phred scores corresponding to each base in seq
                 locs_list - locations corresponding to each base in seq
        Tests:
        All matching test, read.POS is at 1, cigar is 7M
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_aligned_locs()
        [1, 2, 3, 4, 5, 6, 7]

        Test internal deletions
        >>> read.CIGAR = "2M3D5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.aligned_locs = None
        >>> read.get_aligned_locs()
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        Test internal introns
        >>> read.CIGAR = "2M3N5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.aligned_locs = None
        >>> read.get_aligned_locs()
        [1, 2, 6, 7, 8, 9, 10]

        Test internal insertions
        >>> read.CIGAR = "2M3I2M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.aligned_locs = None
        >>> read.get_aligned_locs()
        [1, 2, 3, 4]

        Test soft clipping on either end
        >>> read.CIGAR = "2S5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.aligned_locs = None
        >>> read.get_aligned_locs()
        [1, 2, 3, 4, 5]
        >>> read.CIGAR = "5M2S"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.aligned_locs = None
        >>> read.get_aligned_locs()
        [1, 2, 3, 4, 5]
        """
    
        # get both the seq with no gaps and the start and stop
        # locations
        if self.aligned_locs:
            return self.aligned_locs
        locs_list = []
        locs = self.get_aligned_blocks()
        # Loop through to make a list of locations that corresponds
        # to the location of each base in the final gapped read with the
        # gaps removed
        for start, end in locs:
            locs_list.extend(range(start, end))
        self.aligned_locs = locs_list
        return self.aligned_locs


    def start_end_gaps(self, paired):
        """ Function to get the start, end, and locations within gaps within
        a sequence. Only returns the gaps for the individual read under
        consideration. Thus, if paired is True, it will return the start
        and end for the entire pair, but only the locations within gaps for
        the individual read it was called on.
        
        Inputs: paired - boolean. Treat read as SE or PE. If paired is True,
                         returns the pos + TLEN as the end.
                         If false then returns the right most mapped location
                         of the read itself

        Modifies: nothing

        Returns: start  - integer. left most start of the read or fragment
                 end    - integer. right most end of the read or fragment
                 gaps   - list. locations that reside within gaps

        Tests:
        
        All matching test, read.POS is at 1, cigar is 7M, read.TLEN is 12
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.start_end_gaps(False)
        (1, 8, [])
        >>> read.start_end_gaps(True)
        (1, 13, [])


        Test internal deletions, read.POS is at 1, read.TLEN is 12
        >>> line = "testQ.1.testR.2.30.2M3D5M.=.4.12.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.start_end_gaps(False)
        (1, 11, [])
        >>> read.start_end_gaps(True)
        (1, 13, [])

        Test internal introns
        >>> line = "testQ.1.testR.2.30.2M3N5M.=.4.12.AGTCGCT.!#%()+,.NM:i:0" 
        >>> read = SamAlignment(line, sep = '.')
        >>> read.start_end_gaps(False)
        (1, 11, [(3, 6)])
        >>> read.start_end_gaps(True)
        (1, 13, [(3, 6)])

        Test internal insertions
        >>> line = "testQ.1.testR.2.30.2M3I2M.=.4.12.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.start_end_gaps(False)
        (1, 5, [])
        >>> read.start_end_gaps(True)
        (1, 13, [])

        Test soft clipping on either end
        >>> line = "testQ.1.testR.4.30.2S5M.=.4.12.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.start_end_gaps(False)
        (3, 8, [])
        >>> read.start_end_gaps(True)
        (3, 15, [])
        >>> line = "testQ.1.testR.2.30.5M2S.=.4.12.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.start_end_gaps(False)
        (1, 6, [])
        >>> read.start_end_gaps(True)
        (1, 13, [])
        """
        # if the data is paired, then determining the start and the end of
        # the segement is simple
        if paired:
            # make sure that the pairs mapped on the same chromosome
            if self.RNEXT == "=":
                # take the minimum left-most start of the two reads
                start = min(self.POS, self.PNEXT)
                # find the end by adding the abs value of the TLEN
                end = start + abs(self.TLEN)
            else:
                raise ValueError("Cannot determine start and end for pair," +
                                 " the pair is chimeric. RNEXT != '=' %s"%(
                                     self.QNAME))
        else:
            start = self.POS
            # grab the end value from the last tuple that is returned from
            # self.get_aligned_blocks to determine the end of the read
            end = (self.get_aligned_blocks())[-1][-1]

        return(start, end, self.get_gap_blocks())


    def reconstruct_reference(self):

        """Function to parse the MD field and use the information in it
        to recreate the reference sequence. Simultaneously determines where
        mutations in the read from the reference are and returns them as a
        list of tuples containing relevant information.

        Input: nothing
        Modifies: nothing
        Returns: string. Genome sequence reconstructed from MD field
                 list. Tuples containing mutations in the read from the 
                 reference, [(loc, ref_base, mut_base, phred, loc)]

        Tests:
        All matching test
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,.MD:Z:7"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        ('AGTCGCT', [])

        Missing MD field
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        Traceback (most recent call last):
        KeyError: 'MD field not found for read testQ'

        Every other mis-matching test
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,.MD:Z:1A1A1A1"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        ('AATAGAT', [(2, 'A', 'G', '#'), (4, 'A', 'C', '('), (6, 'A', 'C', '+')])

        Mis-matches on ends
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,.MD:Z:G5G"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        ('GGTCGCG', [(1, 'G', 'A', '!'), (7, 'G', 'T', ',')])

        Test internal deletions, 
        >>> line = "testQ.1.testR.2.30.2M3D5M.=.4.12.AGTCGCT.!#%()+,.MD:Z:2^AAA5"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        ('AGAAATCGCT', [(3, 'A', '-', '-'), (4, 'A', '-', '-'), (5, 'A', '-', '-')])

        Test internal introns
        >>> line = "testQ.1.testR.2.30.2M3N5M.=.4.12.AGTCGCT.!#%()+,.MD:Z:7" 
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        ('AGTCGCT', [])

        Test internal insertions
        >>> line = "testQ.1.testR.2.30.2M3I2M.=.4.12.AGTCGCT.!#%()+,.MD:Z:4"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        ('AGCT', [])

        Test soft clipping on either end
        >>> line = "testQ.1.testR.4.30.2S5M.=.4.12.AGTCGCT.!#%()+,.MD:Z:5"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        ('TCGCT', [])
        >>> line = "testQ.1.testR.2.30.5M2S.=.4.12.AGTCGCT.!#%()+,.MD:Z:5"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.reconstruct_reference()
        ('AGTCG', [])
        """

        # If cached already, return cached values
        if self.aligned_reference and self.aligned_muts:
            return self.aligned_reference, self.aligned_muts
        # get the aligned sequence and phred scores
        aligned_seq, aligned_phred = self.get_aligned_seq_and_phred() 
        # get the aligned locations
        aligned_locs = self.get_aligned_locs()
        # create a list to hold the split up MD values
        MD_list = []
        # compile a list of possible values to split on
        splitters = set(['A','G','T','C','^', ''])
        # split using regular expressions. The parenthesis give us the
        # seperators as well
        try:
            MD_list = re.split('([AGTC^])', self.OPT['MD']) 
        except KeyError:
            raise KeyError("MD field not found for read %s"%self.QNAME)
        # make a list to hold of the mutations that are found
        mut_list = []
        # make a list to hold the reference sequence
        reference_seq = []
        # initialize values to hold the locations that we are at
        start = 0
        end = None
        # loop through all the values that were in the MD field
        for val in MD_list:
            # test to see if the value is a splitter, if it is one of the
            # deletion or empty string splitters, then we just continue
            if val in splitters:
                if val == '^' or val == '':
                    continue
                else:
                    # however if it is not a deletion or empty string, we
                    # append this value to the reference sequence
                    reference_seq.append(val)
                    # we also append mutation information to the mutation list
                    mut_list.append((aligned_locs[start], val, 
                                     aligned_seq[start],
                                     aligned_phred[start]))
                    # finally we increment the index of where we are on the
                    # read
                    start = start + 1
            else:
                # if the value is not a splitter, then it is suggesting that the
                # self.SEQ has the correct information for the sequence
                # We then set the end to equal the current index location plus
                # the number of bases that match
                end = start + int(val)
                # we append all of this to the reference sequence
                reference_seq.append(aligned_seq[start:end])
                # and change the end to become the start
                start = end 
        # finally we join the list of bases in the reference sequence to a
        # string and cache it along with the mutation list
        self.aligned_reference = ''.join(reference_seq)
        self.aligned_muts = mut_list
        # we finally return the cached values
        return self.aligned_reference, self.aligned_muts

    def record_muts(self, strand, phred_encoding=33):
        """ Function to find mutations found in the read as compared to the
        reference. Considers all reads individually (doesn't matter if paired 
        or not).  This depends on the presence of both a CIGAR field and an
        MD field in the alignment. If these are not present then this cannot
        be used. It also requires the XM field to be present

        Inputs: strand           - char. ["+" | "-" | "."] what strand was the
                                   read from? If "-" will complement the
                                   mutation and its base.
                [phred_encoding] - int. what is the phred-encoding for the read?
                                   the default is 33.
        Modifies: nothing
        Returns: A list of mutation records in the form: 
        (chrm, loc, strand, base, mut_base, phred)

        Tests:
        All matching test
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,.MD:Z:7.XM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.record_muts('+')
        []

        Every other mis-matching test
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,.MD:Z:1A1A1A1.XM:i:3"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.record_muts("+")
        [('testR', 2, '+', 'A', 'G', 2), ('testR', 4, '+', 'A', 'C', 7), ('testR', 6, '+', 'A', 'C', 10)]
        >>> read.record_muts("-")
        [('testR', 2, '-', 'T', 'C', 2), ('testR', 4, '-', 'T', 'G', 7), ('testR', 6, '-', 'T', 'G', 10)]

        Mis-matches on ends
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,.MD:Z:G5G.XM:i:2"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.record_muts("+")
        [('testR', 1, '+', 'G', 'A', 0), ('testR', 7, '+', 'G', 'T', 11)]
        >>> read.record_muts("-")
        [('testR', 1, '-', 'C', 'T', 0), ('testR', 7, '-', 'C', 'A', 11)]

        Test internal deletions
        >>> line = "testQ.1.testR.2.30.2M3D5M.=.4.12.AGTCGCT.!#%()+,.MD:Z:2^AAA5.XM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.record_muts("+")
        []

        Missing XM field
        >>> line = "testQ.1.testR.2.30.7M.=.4.12.AGTCGCT.!#%()+,.MD:Z:7"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.record_muts("+")
        Traceback (most recent call last):
        KeyError: 'XM field not found for read testQ'
        """

        mutation_list = []
        # test if there are any mutations according to the XM field of the
        # SAM record
#        try:
#            if self.OPT['XM'] == 0:
#                return mutation_list
#        except KeyError:
#            raise KeyError("XM field not found for read %s"%self.QNAME)

        # muts is coming in the form [(loc, base, mut_base, phred)...]
        muts = self.reconstruct_reference()[1]
        # Go through each mutation and build a mutation record
        for mut in muts:
            # skip any deletions where the read had no information
            if mut[2] == "-":
                continue
            # gather general information that does not depend on strand
            chrm = self.RNAME
            loc = mut[0]
            phred = ord(mut[3])-phred_encoding
            # figure out what the mutation actually was using strand information
            if strand == "-":
                base = complement(mut[1])
                mut_base = complement(mut[2])
            else:
                base = mut[1]
                mut_base = mut[2]
            # finally we add the new mutation record a list to return
            mutation_list.append((chrm, loc, strand, base, mut_base, phred))
        return(mutation_list)

    def get_ref_base_locs(self, base, strand):

        """ Function to return the location of where a particular base
        is in the sequence of the read. Considers only the locations
        that actually have sequence for the read, not its pair or the
        interpolated region between the pair and itself.

        Inputs: base    - char. The identity of the reference base 
                strand  - char. ["+" | "-" | "."] If + treats the base as is,
                          if - searches for the bases complement. If . it 
                          searches for both the base and its complement (
                          essentially search for base pairs)
        Modifies: nothing
        Returns: a numpy array of locations where the base occurs on the read.

        Tests:
        
        All base of interest
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.TTTTTTT.!#%()+,.NM:i:0.MD:Z:7"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_ref_base_locs('T', '+')
        [1, 2, 3, 4, 5, 6, 7]
        >>> read.get_ref_base_locs('A', '-')
        [1, 2, 3, 4, 5, 6, 7]

        Mix of bases and complement
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.TATATAT.!#%()+,.NM:i:0.MD:Z:7"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_ref_base_locs('T', '+')
        [1, 3, 5, 7]
        >>> read.get_ref_base_locs('T', '-')
        [2, 4, 6]
        >>> read.get_ref_base_locs('T', '.')
        [1, 3, 5, 7, 2, 4, 6]

        Base at end
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.TTTTTTG.!#%()+,.NM:i:0.MD:Z:7"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_ref_base_locs('G', '+')
        [7]
        >>> read.get_ref_base_locs('C', '-')
        [7]

        Base at beginning
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.TGGGGGG.!#%()+,.NM:i:0.MD:Z:7"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_ref_base_locs('T', '+')
        [1]
        >>> read.get_ref_base_locs('A', '-')
        [1]
        """

        # get the reference sequence and the associated locations
        genome_sequence = self.reconstruct_reference()[0]
        aligned_locs = self.get_aligned_locs()

        if strand == "-":
            # NOTE only complementing! not reverse complementing for
            # minus strand read 
            base = complement(base) 

        # use numpy arrays to quickly find where bases are in the reference
        # sequence
        base_index = np.where(np.asarray(list(genome_sequence)) == base)[0]
        locs_list = [aligned_locs[i] for i in base_index]
        # if we don't have stranded information, then we have to consider that
        # the "base" could be from either strand. We therefore search for a 
        # "base-pair" i.e. find all the locations of both the base and its
        # complement
        if strand == ".":
            base_index2 = np.where(np.asarray(list(
                                   genome_sequence)) == complement(base))[0]
            locs_list.extend([aligned_locs[i] for i in base_index2])
        return(locs_list)


class MutDataFrame:
    """ Class for holding mutations in a pandas dataframe, includes
        methods for combining dataframes and calculating q values
        using benjamini-hochberg procedures. Initially mutation information
        is stored in a list of tuples.
    """

    def __init__(self):
        """Initialize the data frame. It starts out being empty.
        Mutations should be added to the mutation list like this:

        self.mutation_list.append((chrom, loc, strand, base, mut, phred))

        Inputs: None

        Modifies: Creates self.df holding the empty df
        and self.mutation_list holding a tuple of the mutation
        information

        Returns: None

        Tests:
        
        Test initialization
        >>> x = MutDataFrame()
        >>> x.mutation_list
        []
        >>> x.df
        Empty DataFrame
        Columns: [chrom, loc, strand, base, mut, phred, pvalue, qvalue]
        Index: []

        """
        # initialize the columns to sit in the data frame
        column_names = ['chrom', 'loc', 'strand', 'base', 'mut', 'phred',
                        'pvalue','qvalue']
        # intialize the data frame
        self.df = pd.DataFrame(columns = column_names)
        # initialize the mutations list that will be later turned into
        # a dataframe
        self.mutation_list = []


    def create_df_from_list(self):
        """Takes the list of mutation data and enters it into the dataframe all
        at once. Only to be done once all mutations are found. Mutation_list
        has to be in the form of: 

        [(chrom1, loc1, strand1, base1, mut1, phred1), ...,(chromN, locN,
        strandN, baseN, mutN, phredN)]
        
        Inputs: None

        Modifies: self.df

        Returns: True if mutations exist, False otherwise

        Tests:

        Test no mutations
        >>> x = MutDataFrame()
        >>> x.create_df_from_list()
        False

        Test mutations
        >>> x.mutation_list = [("test", 1, "+", "A", "T", 30),("test", 2, "-","G","C",2)]
        >>> x.create_df_from_list()
        True
        >>> x.df
          chrom  loc strand base mut  phred pvalue qvalue
        0  test    1      +    A   T     30    NaN    NaN
        1  test    2      -    G   C      2    NaN    NaN
        """
        # first make sure there is at least one mutation
        if len(self.mutation_list) >= 1:
            # unzip the tuples into seperate lists
            chrom, loc, strand, base, mut, phred = zip(*self.mutation_list)
            # convert each list to a pd series and add to the data frame
            self.df.chrom = pd.Series(chrom)
            self.df.loc = pd.Series(loc)
            self.df.strand = pd.Series(strand)
            self.df.base = pd.Series(base)
            self.df.mut = pd.Series(mut)
            self.df.phred = pd.Series(phred)
            return True
        else:
            # otherwise warn the user that no mutations were found
            logging.warning("No mutations found.")
            return False

    def phred_to_p(self):
        """Creates the pvalue column from the phred column in place.

        Inputs: None

        Modifies: self.df['pvalue']

        Returns: None

        Tests:

        Test mutations
        >>> x = MutDataFrame()
        >>> x.mutation_list = [("test", 1, "+", "A", "T", 30),("test", 2, "-","G","C",2)]
        >>> success = x.create_df_from_list()
        >>> x.phred_to_p()
        >>> x.df
          chrom  loc strand base mut  phred    pvalue qvalue
        0  test    1      +    A   T     30  0.001000    NaN
        1  test    2      -    G   C      2  0.630957    NaN
        """
        # create an anonymous function to convert phred scores to p values
        convert_to_p = lambda x: 10**(-x/10.0)
        # apply that across the entire phred series and put the answer in the
        # pvalue series
        self.df['pvalue'] = self.df['phred'].apply(convert_to_p)

    def do_fdr(self, alpha=0.05):
        """ Uses statsmodels module to perform benjamini-hochberg FDR
        correction. Fills in the qvalue column of the data frame. Filters
        using the qvalue > alpha

        Inputs: alpha - float the desired alpha level for FDR correction.

        Modifies: self.df["qvalue"]

        Returns: None

        Tests:

        Test a passing value
        >>> x = MutDataFrame()
        >>> x.mutation_list = [("test", 1, "+", "A", "T", 30),("test", 2, "-","G","C",2)]
        >>> success = x.create_df_from_list()
        >>> x.do_fdr(0.05)
        >>> x.df
          chrom  loc strand base mut  phred  pvalue  qvalue
        0  test    1      +    A   T     30   0.001   0.002
 
        Test a non-passing value
        >>> x = MutDataFrame()
        >>> x.mutation_list = [("test", 2, "-","G","C",2)]
        >>> success = x.create_df_from_list()
        >>> x.do_fdr(0.05)
        >>> x.df
        Empty DataFrame
        Columns: [chrom, loc, strand, base, mut, phred, pvalue, qvalue]
        Index: []
        """
        # defaults to assuming phred scores havent been converted to p values
        self.phred_to_p()
        # calculate the qvalues, alpha does nothing in this regard here
        self.df["qvalue"] = multicomp.fdrcorrection0(self.df["pvalue"], 
                                                     alpha)[1]
        # filter by alpha
        self.df = self.df[self.df["qvalue"] < alpha]

    def filter_mut_type(self, base, mut, data_frame=None):
        """Function to grab mutations of a certain type. Can be coupled
        with filter_strand or filter_region by specifying the output data_frame
        from either function

        Inputs: base         - char reference base you are looking for
                mut          - char mutation base you are looking for
                [data_frame] - dataframe from a previous filter. Default = None
                               When None, will use self.df
        Modifies: nothing

        Returns: pandas dataframe object with results of the filter

        Example: get all the T to C mutations from region chr1:1-100 that were
        from the "+" strand:
        
        filtered_df = Mut_Df.filter_strand('+', 
                            data_frame= Mut_Df.filter_region("chr1", (1,100), 
                            data_frame=Mut_Df.filter_mut_type(T, C)))
        Tests:
        >>> x = MutDataFrame()
        >>> x.mutation_list = [("test", 1, "+", "A", "T", 30),("test", 2, "-","G","C",2)]
        >>> success = x.create_df_from_list()
        >>> x.filter_mut_type("A", "T")
          chrom  loc strand base mut  phred pvalue qvalue
        0  test    1      +    A   T     30    NaN    NaN
        >>> x.filter_mut_type("G", "C")
          chrom  loc strand base mut  phred pvalue qvalue
        1  test    2      -    G   C      2    NaN    NaN

        """
        # if they are passing in a dataframe use that, otherwise use
        # self.df
        if data_frame is not None:
            df = data_frame
        else:
            df = self.df
        # if you want all bases that mutate to the mut then just return the
        # locations that match the mutation
        if base == "*":
            return df[df['mut'] == mut]

        # if you want all mutations of a base  then just return the locations
        # that match that base
        elif mut == "*":
            return df[df['base'] == base]
        else:
        # otherwise match the mut and the base
            df = df[(df['base'] == base) & (df['mut'] == mut)]
        return df

    def filter_region(self, chrom, locs=None, data_frame=None):

        """Function to grab mutations within a certain region. Can be coupled
        with filter_mut_type or filter_strand by specifying 
        data_frame = output_from <self.mut_type> on each subsequent function

        Inputs: chrom        - char chromosome you are looking for
                locs         - tuple (start, stop) tuple containing the region
                               you want to find mutations in. Region is
                               inclusive.
                [data_frame] - dataframe from a previous filter. Default = None
                               When None, will use self.df
        Modifies: nothing

        Returns: pandas dataframe object with results of the filter

        Example: get all the T to C mutations from region chr1:1-100 that were
        from the "+" strand:
        
        filtered_df = Mut_Df.filter_strand('+', 
                            data_frame= Mut_Df.filter_region("chr1", (1,100), 
                            data_frame=Mut_Df.filter_mut_type(T, C)))
        Tests:
        >>> x = MutDataFrame()
        >>> x.mutation_list = [("test", 1, "+", "A", "T", 30),("test", 2, "-","G","C",2)]
        >>> success = x.create_df_from_list()
        >>> x.filter_region("test", (0,1))
          chrom  loc strand base mut  phred pvalue qvalue
        0  test    1      +    A   T     30    NaN    NaN
        >>> x.filter_region("test", (2,3))
          chrom  loc strand base mut  phred pvalue qvalue
        1  test    2      -    G   C      2    NaN    NaN
        """

        if data_frame is not None:
            df = data_frame
        else:
            df = self.df
        if locs:
            df = df[(df['chrom'] == chrom) & 
                 ((df['loc'] >= locs[0]) & (df['loc'] <= locs[1]))]
        else:
            df = df[(df['chrom']==chrom)]

        return df

    def filter_strand(self, strand, data_frame=None):

        """Function to grab mutations from a certain strand. Can be coupled
        with filter_mut_type or filter_region by specifing the output data_frame
        from either function

        Inputs: strand       - char strand you want mutations from. Can be
                               '+', '-' or '.'.

                [data_frame] - dataframe from a previous filter. Default = None
                               When None, will use self.df
        Modifies: nothing

        Returns: pandas dataframe object with results of the filter

        Example: get all the T to C mutations from region chr1:1-100 that were
        from the "+" strand:
        
        filtered_df = Mut_Df.filter_strand('+', 
                            data_frame= Mut_Df.filter_region("chr1", (1,100), 
                            data_frame=Mut_Df.filter_mut_type(T, C)))
        Tests:
        >>> x = MutDataFrame()
        >>> x.mutation_list = [("test", 1, "+", "A", "T", 30),("test", 2, "-","G","C",2)]
        >>> success = x.create_df_from_list()
        >>> x.filter_strand("+")
          chrom  loc strand base mut  phred pvalue qvalue
        0  test    1      +    A   T     30    NaN    NaN
        >>> x.filter_strand("-")
          chrom  loc strand base mut  phred pvalue qvalue
        1  test    2      -    G   C      2    NaN    NaN
        """

        if data_frame is not None:
            df = data_frame
        else:
            df = self.df

        df = df[df['strand']==strand]

        return df


    def count_muts(self, data_frame=None):
        """Increments counters for each type of mutation.
        
        Inputs: [data_frame] - dataframe from the results of a filter. If None
                               (default), then uses self.df
        Modifies: nothing

        Returns: dictionary of mutations with mutation type as the key. For 
        example, dict['AT']: 30, if there were 30 mutations that were A->T

        Tests:
        >>> x = MutDataFrame()
        >>> x.mutation_list = [("test", 1, "+", "A", "T", 30),("test", 2, "-","G","C",2)]
        >>> success = x.create_df_from_list()
        >>> count_dict = x.count_muts()
        >>> count_dict['AT'] == 1
        True
        >>> count_dict['GC'] == 1
        True
        >>> count_dict['TA'] == 0
        True
        >>> count_dict['CG'] == 0
        True
        """
        bases = ['A','G','C','T']
        out_dict = {}
        if data_frame is not None:
            df = data_frame
        else:
            df = self.df
        for pair in itertools.permutations(bases, 2):
            base = pair[0]
            mut = pair[1]
            filtered = self.filter_mut_type(base, mut, data_frame=df)
            out_dict[base+mut] = len(filtered)

        return out_dict


class AlignmentStats:

    def __init__(self):
        """ A class to hold various stats about many reads in a samfile
        Initialization just creates a bunch of attributes that act as counters
        for various stats that need to be gathered.
        Tests:
        >>> x = AlignmentStats()
        >>> x.base_counts['A'] == 0
        True
        >>> x.base_counts['G'] == 0
        True
        >>> x.base_counts['C'] == 0
        True
        >>> x.base_counts['T'] == 0
        True
        """
        ## CHECKED FOR ALL READS
        ## Holds the total number of reads that were in the file
        self.total_reads = 0
        # total that overlapped with off-limits positions
        self.total_overlap = 0
        # dictionary holding values for reads filtered by flag:
        self.bad_flag_reads = {}
        # dictionary hold values for reads filtered for not having flag
        self.missed_good_flag_reads = {}
        # mapping for flags to their names:
        self.flag_mapping = {1:"template having multiple segments in sequencing",
                        2:"each segment properly aligned according to the aligner",
                        4:"segment unmapped",
                        8:"next segment in the template unmapped",
                       16:"SEQ being reverse complemented",
                       32:"SEQ of the next segment in the template being reverse complemented",
                       64:"the first segment in the template",
                      128:"the last segment in the template",
                      256:"secondary alignment",
                      512:"not passing filters, such as platform/vendor quality controls",
                     1024:"PCR or optical duplicate",
                     2048:"supplementary alignment"
                       }
        # eliminated by mapq filter
        self.total_mapq = 0
        # couldnt determine if read was sense or not
        self.total_bad_sense = 0
        ## CHECKED ONLY FOR READS THAT PASS FILTERS
        # total that were mapped
        self.total_mapped = 0
        # Paired end counters, considering things as SEGMENTS
        self.paired = {"+": 0,
                       "-": 0,
                       ".": 0}
        self.total_paired_bp = 0
        # considering each READ
        self.unpaired = {"+": 0,
                         "-": 0,
                         ".": 0}
        self.total_read_bp = 0
        ## FOR MUTATION INFORMATION
        bases = ['A', 'G', 'T', 'C']
        self.base_counts = {key: 0 for key in bases}
        # number of mutations before FDR
        self.muts_before_filt = 0
        # total mutation counts (after FDR) for each type of mutation
        self.mut_counts = {''.join(key): 0 for key in itertools.permutations(bases,2)}

    def __add__(self, other):
        """ Create a new AlignmentStats object with the sum of all the attributes
        of the self and other.
        """
        new = AlignmentStats()
        ## CHECKED FOR ALL READS
        new.total_reads = self.total_reads + other.total_reads
        new.total_overlap = self.total_overlap + other.total_overlap
        for key in set(self.bad_flag_reads.keys(), other.bad_flag_reads.keys()):
            new.bad_flag_reads[key] = self.bad_flag_reads.get(key, 0) +\
                                      other.bad_flag_reads.get(key, 0)
        for key in set(self.missed_good_flag_reads.keys(), 
                       other.missed_good_flag_reads.keys()):
            new.missed_good_flag_reads[key] = self.missed_good_flag_reads.get(key, 0) +\
                                              other.missed_good_flag_reads.get(key, 0)
        new.total_mapq = self.total_mapq + other.total_mapq
        new.total_bad_sense = self.total_bad_sense + other.total_bad_sense
        ## CHECKED ONLY FOR READS THAT PASS FILTERS
        new.total_mapped = self.total_mapped + other.total_mapped
        # Paired end
        for key in new.paired.keys():
            new.paired[key] = self.paired[key] + other.paired[key]
        new.total_paired_bp = self.total_paired_bp + other.total_paired_bp
        # Single end
        for key in new.unpaired.keys():
            new.unpaired[key] = self.unpaired[key] + other.unpaired[key]
        new.total_read_bp = self.total_read_bp + other.total_read_bp
        ## FOR MUTATION INFORMATION
        for key in new.base_counts.keys():
            new.base_counts[key] = self.base_counts[key] + other.base_counts[key]
        new.muts_before_filt = other.must_before_filt + self.muts_before_filt
        for key in self.mut_counts.keys():
            new.mut_counts[key] = self.mut_counts[key] + other.mut_counts[key]
        return new

    def increment_bases(self, read, strand):
        """ Takes in SamAlignment read and increments the base_counts
        dictionary with the bases present within the read.

        Tests:
        >>> line = "testQ.1.testR.20.30.5M.testR.25.9.AGTCC.(((((.MD:Z:5"
        >>> read = SamAlignment(line, sep = '.')
        >>> x = AlignmentStats()
        >>> x.increment_bases(read, "+")
        >>> x.base_counts["A"] == 1
        True
        >>> x.base_counts["G"] == 1
        True
        >>> x.base_counts["C"] == 2
        True
        >>> x.base_counts["T"] == 1
        True
        """
        # we want to rebuild the reference first
        seq = read.reconstruct_reference()[0]
        # next we figure out if it should be complemented or not
        if strand == "-":
            seq = complement(seq)
        # count the amount of each base
        self.base_counts['T'] += seq.count('T')
        self.base_counts['G'] += seq.count('G')
        self.base_counts['C'] += seq.count('C')
        self.base_counts['A'] += seq.count('A')
        
    def print_map_stats(self, section, mapq, paired, locs=False):
        """ This prints out each of the counters in an easy to read format
        """
        print "MAPPING STATS FOR %s"%section
        if paired:
            print "Paired-End Mode"
        else:
            print "Single-End Mode"
        print 20*"-"
        print "Total reads in file: %s" %self.total_reads
        print ""
        print "Quality Filtering"
        print "Total reads filtered: %s"%(self.total_reads - self.total_mapped)
        print "Filtered for:"
        for flag in self.missed_good_flag_reads.keys():
            print "\tMissing %s: %s"%(self.flag_mapping[flag], self.missed_good_flag_reads[flag])
        for flag in self.bad_flag_reads.keys():
            print "\t%s: %s"%(self.flag_mapping[flag], self.bad_flag_reads[flag])
        if mapq is None:
            pass
        elif mapq >= 0:
            print "\tMAPQ < %s : %s"%(mapq, self.total_mapq)
        elif mapq < 0:
            print "\tMAPQ > %s : %s"%(mapq, self.total_mapq)
        if locs:
            print "\tOverlap with locations: %s"%(self.total_overlap)
        print "\tStrandedness could not be determined: %s"%(self.total_bad_sense)
        if paired:
            print "Filtered pairs that mapped to plus strand: %s"%(
                                                               self.paired["+"])
            print "Filtered pairs that mapped to minus strand: %s"%(
                                                               self.paired["-"])
            print "Filtered pairs that mapped to either strand: %s"%(
                                                               self.paired["."])
            print "Total pair basepairs considered: %s"%(self.total_paired_bp)
        else:
            print "Filtered reads that mapped to plus strand: %s"%(
                                                        self.unpaired["+"])
            print "Filtered reads that mapped to minus strand: %s"%(
                                                        self.unpaired["-"])
            print "Filtered reads that mapped to either strand: %s"%(
                                                        self.unpaired["."])
            print "Total read basepairs considered: %s"%(self.total_read_bp)

    def calc_mut_rates(self, stranded):
        """ Calculates mutation rate for each type of mutation

        Tests:
        >>> x = AlignmentStats()
        >>> x.base_counts['G'] = 10
        >>> x.base_counts['C'] = 5
        >>> x.mut_counts['GC'] = 2
        >>> x.mut_counts['CG'] = 4
        >>> y = x.calc_mut_rates(True)
        >>> y['GC']
        0.2
        >>> y['CG']
        0.8
        >>> y=x.calc_mut_rates(False)
        >>> y['GC']
        0.4
        >>> y['CG']
        0.4
        """

        bases = ['A', 'G', 'C', 'T']
        out_dict = {}
        for pair in itertools.permutations(bases, 2):
            base = pair[0]
            mut = pair[1]
            try:
                if stranded:
                    mut_rate = (self.mut_counts[base+mut]+0.0)/\
                                self.base_counts[base]
                else:
                    mut_rate = (self.mut_counts[base+mut] + 0.0 +
                                self.mut_counts[complement(base)
                                                +complement(mut)])/(
                                self.base_counts[base] + 
                                self.base_counts[complement(base)])
            except ZeroDivisionError:
                mut_rate = 0.0
            out_dict[base+mut] = mut_rate
        return out_dict

    def print_mut_stats(self, stranded): 

        print "Mutation Stats"
        print 20*"-"
        print "Total mutations before FDR: %s"%self.muts_before_filt
        print "Total mutations after FDR: %s"%sum(self.mut_counts.values())
        for mut_type in self.mut_counts.keys():
            print "Total %s->%s mutations after FDR: %s"%(mut_type[0],
                                                          mut_type[1], 
                                                      self.mut_counts[mut_type])
        for base in self.base_counts.keys():
            print "Total reference %s sequenced: %s"%(base, self.base_counts[base])
        print "Mutation Rate Matrix (after filtering):"
        print "Original Base on left, Mutation base on top"
        print "\tA\tG\tT\tC"
        mut_rates = self.calc_mut_rates(stranded)
        print "A\tX\t%.3e\t%.3e\t%.3e"%(mut_rates["AG"], mut_rates["AT"],
                                  mut_rates["AC"])
        print "G\t%.3e\tX\t%.3e\t%.3e"%(mut_rates["GA"], mut_rates["GT"],
                                  mut_rates["GC"])
        print "T\t%.3e\t%.3e\tX\t%.3e"%(mut_rates["TA"], mut_rates["TG"],
                                  mut_rates["TC"])
        print "C\t%.3e\t%.3e\t%.3e\tX"%(mut_rates["CA"], mut_rates["CG"],
                                  mut_rates["CT"])

    def write_mut_stats(self, out_prefix, stranded):
        mut_rates = self.calc_mut_rates(stranded)
        file_out = out_prefix + "_mut_rates.txt"
        with open(file_out, 'w') as fn:
            for key in mut_rates.keys():
                fn.write("%s\t%f\n"%(key, mut_rates[key]))

    def write_exposure_df(self, out_prefix, paired=False):
        file_out = out_prefix + "_exposure.txt"
        with open(file_out, 'w') as fn:
            if paired:
                fn.write("%s\t%s\n" %((self.total_mapped/2), 
                                     self.total_paired_bp))
            else:
                fn.write("%s\t%s\n" %((self.total_mapped), self.total_read_bp))
