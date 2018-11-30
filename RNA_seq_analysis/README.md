## Introduction ##
This directory contains custom python and R scripts used in the analysis of
RNA-seq data. Additionally, this file contains instructions for a pipeline
using these scripts, and available software, to determine Lrp-dependent RNA
expression

## Pre-processing of RNA-Seq data ##
In order to remove adapters and low quality reads we used a combination of
cutadapt and trimmomatic. Fastqc was used for visualization of read quality

Here
arg0 is a sample prefix
arg1 represents the fastq for R1 of the paired end sequencing
arg2 represents the fastq for R2 of the paired end sequencing
```
fastqc arg1 -f fastq -o fastqc_before/	5
fastqc arg2 -f fastq -o fastqc_before/	5
cutadapt --quality-base=33 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -n 3 -m 20 --mask-adapter --match-read-wildcards -o arg0_R1_cutadapt.fastq.gz -p arg0_R2_cutadapt.fastq.gz arg1 arg2 > arg0_cutadapt.log 2> arg0_cutadapt.err
TrimmomaticPE -phred33 arg0_R1_cutadapt.fastq.gz arg0_R2_cutadapt.fastq.gz arg0_R1_trim_paired.fastq.gz arg0_R1_trim_unpaired.fastq.gz arg0_R2_trim_paired.fastq.gz arg0_R2_trim_unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 > arg0_trim.log 2> arg0_trim.err
fastqc arg0_R1_trim_paired.fastq.gz -f fastq -o fastqc_after/
fastqc arg0_R2_trim_paired.fastq.gz -f fastq -o fastqc_after/
```

## Alignment for removal of high abundance RNAs ##
Next, each pair of processed reads was aligned with bowtie2 in order to find
reads to remove for downstream analysis

Here
arg0 is substituted for a samples prefix
arg1 represents the trimmed paired fastq for R1 of the paired end sequencing
arg2 represents the trimmed unpaired fastq for R1 of the paired end sequencing
arg3 represents the trimmed paired fastq for R2 of the paired end sequencing
arg4 represents the trimmed unpaired fastq for R2 of the paired end sequencing

```
bowtie2 -x /home/mbwolfe/genomes/bowtie2indices/ATCC_47076 -1 arg1 -2 arg3 -U arg2,arg4 -q --end-to-end --very-sensitive -p 5 --phred33 --dovetail 2> arg0_bow.log | samtools view -bSh - > arg0.bam 
```

## Filtering of high abundance RNAs ##
In order to find reads that overlapped with our high abundance RNA set we used
the custom python script `parse_for_overlaps.py` to find reads that
overlapped with the locations of interest. NOTE: You will need the
`sam_utils.py` file in your path (or copied over to the directory
`parse_for_overlaps.py` is being run from). This file is in `../ChIP_analysis/`.

Here
arg0 represents the .bam of interest

```
samtools view arg0.bam | python parse_for_overlaps.py arg0 abundantRNASet_ATCC_47076.gff
```

Next, new fasta were written from the trim_paired.fastq files in the
preprocessing step

Here
arg0 is an outprefix
arg1 is the txt file containing read names to filter
arg2 is the trim_paired.fastq for R1
arg3 is the trim_paried.fastq for R2

```
python ../rewrite_filtered_fastq.py arg1 arg2 arg3 arg0
```

## Using filtered fastqs with kallisto ##

Each filtered fastq pair was then quantified with kallisto using the following
command.

Here
arg0 is an output prefix
arg1 is the filtered fastq for R1
arg2 is the filtered fastq for R2

```
kallisto quant -i ATCC_47076_transcriptome.idx -o arg0 -t 4 -b 100 --rf-stranded arg1 arg2 > arg0_kallisto.log 2> arg0_kallisto.err
```

## Estimating log2(WT/KO) Expression ratio ##

Next the .hd5 files were unpacked for use with the next script.
These files must be unpacked into a subdirectory called `/bootstraps/`
in the kallisto directories for each sample. In order to do this, one
can write a simple bash script that looks something like this:

```
# this assumes you ran kallisto all in one directory and you are in
# that directory
for dir in *kallisto_dirs; # <- *kallisto_dirs is a wild card that encompasses all your sample names

do
    # change into a kallisto data directory
    cd $dir
    # dump the .h5 file into a bootstraps directory
    /path/to/kallisto/binary/kallisto hd5dump abundance.h5 -o bootstraps
    # change back up to the main directory
    cd -
done;
```

Finally, the log2(WT/KO) expression ratio estimates were calculated
from each kallisto output directory using the custom script `calculated_expr.py` 
after unpacking the .hd5 files for each sample 

Here
arg0 indicates a sample name
arg1 is WT rep 1
arg2 is WT rep 2
arg3 is KO rep 1
arg4 is KO rep 2

```
python calculate_expr.py arg1 arg2 arg3 arg4 arg0_log2_wt_ko_ratio.txt
```

## Wald test results ##
Wald test results were obtained using the Sleuth companion R package as
described in the methods.

For any questions or clarifications please contact Michael Wolfe at
mbwolfe@umich.edu
