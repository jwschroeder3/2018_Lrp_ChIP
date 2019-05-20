## Introduction ##
This directory contains custom python and R scripts used in the analysis of
ChIP-seq data. Additionally, this file contains instructions for a pipeline
using these scripts, and available software,  to call final high-confidence Lrp
binding sites

## Pre-processing of ChIP data ##
For each fastq file, the following steps were taken to preprocess the data

Reads were concatenated from multiple files into a single large fastq.gz file

Here:  
arg0 is substituted for a samples prefix

```bash
cat *R1*.fastq.gz > arg0_all_R1.fastq.gz
cat *R2*.fastq.gz > arg0_all_R2.fastq.gz
```

Illumina adapters were trimmed from reads.

Here:
arg0 is substituted for a samples prefix
arg1 represents the fastq for R1 of the paired end sequencing
arg2 represents the fastq for R2 of the paired end sequencing

```bash
fastqc arg1 -f fastq -o fastqc_before/
fastqc arg2 -f fastq -o fastqc_before/
cutadapt --quality-base=33 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -n 3 -m 20 --mask-adapter --match-read-wildcards -o arg0_R1_cutadapt.fastq.gz -p arg0_R2_cutadapt.fastq.gz arg1 arg2 > arg0_cutadapt.log 2> arg0_cutadapt.err
java -jar trimmomatic-0.39.jar PE -phred33 arg0_R1_cutadapt.fastq.gz arg0_R2_cutadapt.fastq.gz arg0_R1_trim_paired.fastq.gz arg0_R1_trim_unpaired.fastq.gz arg0_R2_trim_paired.fastq.gz arg0_R2_trim_unpaired.fastq.gz TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 > arg0_trim.log 2> arg0_trim.err
fastqc arg0_R1_trim_paired.fastq.gz -f fastq -o fastqc_after/
fastqc arg0_R2_trim_paired.fastq.gz -f fastq -o fastqc_after/
```

If you have single end data:

```bash
fastqc arg1 -f fastq -o fastqc_before/
cutadapt --quality-base=33 -a AGATCGGAAGAGC -a CGGAAGAGCACAC -n 3 -m 20 --mask-adapter --match-read-wildcards -o arg0_R1_cutadapt.fastq.gz arg1 > arg0_cutadapt.log 2> arg0_cutadapt.err
java -jar trimmomatic-0.39.jar SE -phred33 arg0_R1_cutadapt.fastq.gz arg0_R1_trim_unpaired.fastq.gz TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 > arg0_trim.log 2> arg0_trim.err
fastqc arg0_R1_trim_paired.fastq.gz -f fastq -o fastqc_after/
```

## Alignment of processed reads ##
Next, each pair of processed reads was aligned with bowtie2.

Here
arg0 is substituted for a samples prefix
arg1 represents the trimmed paired fastq for R1 of the paired end sequencing
arg2 represents the trimmed paired fastq for R2 of the paired end sequencing
arg3 represents the trimmed unpaired fastq for R1 of the paired end sequencing
arg4 represents the trimmed unpaired fastq for R2 of the paired end sequencing

```bash
bowtie2 -x /home/mbwolfe/genomes/bowtie2indices/ATCC_47076 -1 arg1 -2 arg2 -U arg3,arg4 -X 2000 -q --end-to-end --very-sensitive -p 5 --phred33 --dovetail 2> arg0_bow.log | samtools view -bSh - > arg0.bam
```

## Mapping of coverage and preparation for bootstrapping ##
Next, reads were filtered with samtools and made into a sampling object using
the custom script `bootstrap_sam_file.py`'s `parse` option.
### Note: if you have single-end reads, DO NOT include the `-f 3` flag, as it will exclude all reads without a proper pair, and you ain't got 'em.

Here
arg0 is the sample prefix
arg1 is the input bam file

```bash
samtools view -f 3 -F 2828 -q 30 arg1 | python3 bootstrap_sam_file.py parse - arg0.ob --paired 2> arg0_sampler.err
```

## Obtaining summary tracks at 10 bp resolution ##
To obtain RE and RSE summary signals we used the samplers from the previous step
as input the the `bootstrapped_chip_no_consolidation.py` script.

Here
arg0 is the sample prefix
arg1 and arg2 are samplers for the WT extracted samples
arg3 and arg4 are samplers for the WT input samples
arg5 and arg6 are samplers for the KO extracted samples
arg7 and arg8 are samplers for the KO input samples

Note that arg1 and arg3 are paired, arg2 and arg4 are paired etc.

```bash
python3 bootstrapped_chip_no_consolidation.py --sample_name_luts run_info --genome_size 4215607 --out_prefix arg0 --ChIP_samps arg1 arg2 --inp_samps arg3 arg4 --ChIP_conts arg5 arg6 --inp_conts arg7 arg8 --num_replicates 1 --identity -s 1234 -p 8 --save_summaries 0.05 --resolution 10 2> arg0.log
```

## Obtaining bootstrap replicate summary statistics ##
To obtain the bootstrap MAD stat (as well as additional statistics) from 1000
bootstrap replicates for each 10 bp location in the genome we also used the
`bootstrapped_chip_no_consolidation.py` script.

Here
arg0 is the sample prefix
arg1 and arg2 are samplers for the WT extracted samples
arg3 and arg4 are samplers for the WT input samples
arg5 and arg6 are samplers for the KO extracted samples
arg7 and arg8 are samplers for the KO input samples

Note that arg1 and arg3 are paired, arg2 and arg4 are paired etc.

```bash
python2.7 bootstrapped_chip_no_consolidation.py 4215607 arg0 --ext_samps arg1 arg2 --inp_samps arg3 arg4 --ext_conts arg5 arg6 --inp_conts arg7 arg8 --num_replicates 1000 -s 1234 -p 8 --save_summaries 0.05 --resolution 10 2> arg0_1000_bootstrap.log
```

## Calculating IDR statistic ##
To calculate the IDR statistic. Each combination of WT-KO replicates were input
to the `calculate_IDR.R` script

Here
arg0 is the sample prefix
arg1 is one WT-KO replicate
arg2 is a second WT-KO replicate
arg3 is the estimated number of Lrp octamers for that condition

```bash
Rscript calculate_idr.R arg1 arg2 arg3 arg0 > arg0.log 2> arg0.err
```

## Calculating final Lrp binding regions ##
Finally, final Lrp binding regions were called using a combination of actual
binding coverage (Summary tracks at 10 bp resolution), bootstrap MAD (bootstrap
replicate summary statistics), and IDR (Calculating IDR statistic). Using the
`calculate_peaks.py` script

Here
sub1-4.npy represent the RSE replicates
mad1-4.npy represent the corresponding bootstrap mad for each RSE replicate
idr_combo1-2.npy represent each combination of subtracted replicates for the idr
                 calculation
out_peaks is the output prefix

```bash
python2.7 calculate_peaks.py --log2ratios actual_sub1.npy actual_sub2.npy actual_sub3.npy actual_sub4.npy --mad mad_sub1.npy mad_sub2.npy mad_sub3.npy mad_sub4.npy --idr idr_combo1.npy idr_combo2.npy --resolution 10 --bins 3 --outpre out_peaks --idralpha 0.01 --bioalpha 0.001 --techalpha 0.001
```

For any questions or clarifications please contact Michael Wolfe at
mbwolfe@umich.edu
