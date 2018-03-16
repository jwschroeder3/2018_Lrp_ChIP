################################################################################
# Script to estimate log2(WT/KO) TPM ratios from kallisto bootstrap replicates
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
import numpy as np
import os
import glob
import sys

# we have two WT and two KO samples for each replicate
samp_dir1 = sys.argv[1]
samp_dir2 = sys.argv[2]

cont_dir1 = sys.argv[3]
cont_dir2 = sys.argv[4]

outf = sys.argv[5]

# take an unpacked kallisto directory and read it in as a table
def read_dir_as_table(in_dir):
    columns = []
    genes = []
    for this_file in glob.glob(in_dir+"/bs*.tsv"):
        column = []
        with open(this_file) as f:
            f.readline()
            for line in f:
                column.append(float(line.rstrip().split("\t")[4]))
            columns.append(column)
    return np.array(columns)

# read in the actual values presented in the abundance file
def read_actual(in_dir):
    column=[]
    genes = []
    with open(in_dir+"/abundance.tsv") as f:
        f.readline()
        for line in f:
            linearr = line.rstrip().split("\t")
            column.append(float(linearr[4]))
            genes.append(linearr[0])
    return(genes, column)

# calculate the credible interval for the bootstrapping
def credible_interval(in_vector, alpha=0.05):
    size = in_vector.size
    lower = int(np.floor((alpha/2)*size))
    upper = size-lower
    in_vec_sort = np.sort(in_vector)
    return (in_vec_sort[lower], in_vec_sort[upper])

# read in all the data
s1_genes, s1_actual = read_actual(samp_dir1)
s2_genes, s2_actual = read_actual(samp_dir2)
c1_genes, c1_actual = read_actual(cont_dir1)
c2_genes, c2_actual = read_actual(cont_dir2)
# get the average actual ratio log2(average(WT)/average(KO))
actual = np.log2(np.mean(np.array([s1_actual, s2_actual]), axis=0)/np.mean(np.array([c1_actual, c2_actual]), axis=0))
# read in all the bootstrap replicates
samp1 = read_dir_as_table(samp_dir1)
samp2 = read_dir_as_table(samp_dir2)
cont1 = read_dir_as_table(cont_dir1)
cont2 = read_dir_as_table(cont_dir2)

# calculate the average ratio of the bootstrap replicates
average_ratio = np.log2(np.mean(np.array([samp1,samp2]), axis=0)/np.mean(np.array([cont1,cont2]), axis=0))
# get the credible interval of the bootstraps
lower = np.apply_along_axis(credible_interval, 0, average_ratio)
# create the final numpy array
final = np.column_stack((actual, lower.transpose()))
# write out the final array as a text file, add the gene names in the first row
with open(outf, mode="w") as f:
    for gene, row in zip(s1_genes,final):
        f.write(gene+"\t"+"\t".join([str(val) for val in row])+"\n")
