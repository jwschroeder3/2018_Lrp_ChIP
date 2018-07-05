################################################################################
# Script to estimate log2(WT/KO) count ratios from kallisto bootstrap replicates
#
# Written by Michael Wolfe. Modified by Peter Freddolino
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


# take an unpacked kallisto directory and read it in as a table
# note that the values that we read are now the estimated counts, NOT the tpm
def read_dir_as_table(in_dir):
    columns = []
    genes = []
    for this_file in glob.glob(in_dir+"/bootstraps/bs*.tsv"):
        column = []
        with open(this_file) as f:
            f.readline()
            for line in f:
                column.append(float(line.rstrip().split("\t")[3]))
            columns.append(column)
    return np.array(columns)

# read in the actual values presented in the abundance file
# note that we now read in the estimated counts, NOT the tpm
def read_actual(in_dir):
    column=[]
    genes = []
    with open(in_dir+"/abundance.tsv") as f:
        f.readline()
        for line in f:
            linearr = line.rstrip().split("\t")
            column.append(float(linearr[3]))
            genes.append(linearr[0])
    return(genes, np.array(column))

# calculate the credible interval for the bootstrapping
def credible_interval(in_vector, alpha=0.05):
    size = in_vector.size
    lower = int(np.floor((alpha/2)*size))
    upper = size-lower
    in_vec_sort = np.sort(in_vector)
    return (in_vec_sort[lower], in_vec_sort[upper])

def calc_size_factors( dat_mat ):
    # given a data matrix with samples by row and genes by column, 
    # calculate and return the size factor for each sample
    # dat_mat should contain the LOG COUNTS

    #print dat_mat.shape
    ref_vals = np.mean(dat_mat, axis=0)
    norm_consts = np.median( dat_mat - ref_vals, axis=1)
    #print norm_consts
    #print np.mean(norm_consts)
    return np.exp(norm_consts)


if __name__ == "__main__":

    # we have two WT and two KO samples for each replicate
    samp_dir1 = sys.argv[1]
    samp_dir2 = sys.argv[2]
    
    cont_dir1 = sys.argv[3]
    cont_dir2 = sys.argv[4]
    
    outf = sys.argv[5]

    # read in all the data
    s1_genes, s1_actual = read_actual(samp_dir1)
    s2_genes, s2_actual = read_actual(samp_dir2)
    c1_genes, c1_actual = read_actual(cont_dir1)
    c2_genes, c2_actual = read_actual(cont_dir2)

    # get the size factor for each sample
    size_facs = calc_size_factors(np.log(np.vstack( (s1_actual, s2_actual, c1_actual, c2_actual) ) + 0.1 ))
    for i, (sample, dirname) in enumerate(zip([s1_actual, s2_actual, c1_actual, c2_actual],
                                   [samp_dir1, samp_dir2, cont_dir1, cont_dir2])):
        np.savetxt(dirname+"_expr.txt", sample/size_facs[i], fmt="%.4e")


    #print size_facs
    # get the average actual ratio log2(average(WT)/average(KO))
    actual = np.mean(np.column_stack([np.log2(s1_actual/size_facs[0]), np.log2(s2_actual/size_facs[1])]), axis=1) - np.mean(np.column_stack([np.log2(c1_actual/size_facs[2]), np.log2(c2_actual/size_facs[3])]),axis=1)
    # read in all the bootstrap replicates, correcting for size factors
    samp1 = read_dir_as_table(samp_dir1) / size_facs[0]
    samp2 = read_dir_as_table(samp_dir2) / size_facs[1]
    cont1 = read_dir_as_table(cont_dir1) / size_facs[2]
    cont2 = read_dir_as_table(cont_dir2) / size_facs[3]

    # calculate the average ratio of the bootstrap replicates
    ratio1 = np.log2(samp1/cont1)
    ratio2 = np.log2(samp2/cont2)
    ratio3 = np.log2(samp1/cont2)
    ratio4 = np.log2(samp2/cont1)
    # get the credible interval of the bootstraps
    lower1 = np.apply_along_axis(credible_interval, 0, ratio1)
    lower2 = np.apply_along_axis(credible_interval, 0, ratio2)
    lower3 = np.apply_along_axis(credible_interval, 0, ratio3)
    lower4 = np.apply_along_axis(credible_interval, 0, ratio4)
    lower = np.dstack([lower1.transpose(), lower2.transpose(), lower3.transpose(),
                       lower4.transpose()])
    overall_lower = np.apply_along_axis(min, 2, lower)[:,0]
    overall_upper = np.apply_along_axis(max, 2, lower)[:,1]
    final = np.column_stack((actual, overall_lower, overall_upper))
    
    ## UNCOMMENT HERE TO SWITCH TO THE OLD MEAN METHOD OF CALCULATING THE CI
    #average_ratio = np.mean(np.array([np.log2(samp1),np.log2(samp2)]), axis=0)-np.mean(np.array([np.log2(cont1),np.log2(cont2)]), axis=0)
    #lower = np.apply_along_axis(credible_interval, 0, average_ratio)
    # create the final np array
    #final = np.column_stack((actual, lower.transpose()))
    ##


    # write out the final array as a text file, add the gene names in the first row
    with open(outf, mode="w") as f:
        for gene, row in zip(s1_genes,final):
            f.write(gene+"\t"+"\t".join([str(val) for val in row])+"\n")
