################################################################################
# Script to read in a fastq file a remove reads that match a set of read names
# to eliminate. Writes a new fastq without these reads
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
import sys
import gzip
import os
import shutil

def read_in_bad_reads(infile):
    x = set()
    with open(infile) as f:
        for line in f:
            x.add(line.rstrip())
    return x

def write_out_read(read_list, fhandle):
    for line in read_list:
        fhandle.write(line)


def parse_file(filename, n=4):
    with open(filename, mode="r") as f:
        for i in xrange(n):
            lines.append(f.readline())
        yield lines

def process_file(filename, outfilename, bad_reads):
    with open(outfilename, mode="w") as f:
        with gzip.open(filename, mode="r") as infile:
            read = [1]
            while '' not in read:
                read = []
                for i in xrange(4):
                    read.append(infile.readline())
                try:
                    readname = read[0].rstrip().split()[0][1:]
                except IndexError:
                    print("End of File reached %s"%(read))
                    break
                if readname in bad_reads:
                    continue
                else:
                    write_out_read(read, f)
    print("compressing %s"%outfilename)
    with open(outfilename, mode="rb") as fin:
        with gzip.open(outfilename+".gz", "wb") as fout:
            shutil.copyfileobj(fin, fout)
    print("done")
    os.unlink(outfilename)




def main():
    bad_reads = sys.argv[1]
    in_file_R1 = sys.argv[2]
    in_file_R2 = sys.argv[3]
    out_pre = sys.argv[4]
    bad_reads = read_in_bad_reads(bad_reads)
    process_file(in_file_R1, out_pre+"_R1.fastq", bad_reads)
    process_file(in_file_R2, out_pre+"_R2.fastq", bad_reads)

main()

    
