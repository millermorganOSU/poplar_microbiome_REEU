#!/bin/sh

## Install the following programs (brew install)
# bowtie2
# samtools 
# bedtools

# Start by indexing your reference sequence
# $ bowtie2-build Reference_sequence.fasta Species
# Indexing only needs to happen once

# This script will take care of one sample (FWD and REV reads)

# Arguments
host_i=$1 # host index 		Static
read1_i=$2 # read 1 input	variable
read2_i=$3 # read 2 input	variable
host_o_s=$4 # sam output	Static
host_o_b=$5 # bam output	Static
name=$6 # sample name		variable


# Expected Arguments

# ./align.sh reference/Ptrichocarpa myco.03.01.a.R1.filt.fastq.gz myco.03.01.a.R2.filt.fastq.gz sam/Ptrichocarpa bam/Ptrichocarpa nh/myco.03.01

# Aligning the reads to the host reference
bowtie2 -x ${host_i} -1 ${read1_i} -2 ${read2_i} -S ${host_o_s}_mapped_and_unmapped.sam
samtools view -S -b ${host_o_s}_mapped_and_unmapped.sam > ${host_o_b}_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 ${host_o_b}_mapped_and_unmapped.bam > ${host_o_b}_bothEndsUnmapped.bam
samtools sort -n ${host_o_b}_bothEndsUnmapped.bam -o ${host_o_b}_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i ${host_o_b}_bothEndsUnmapped_sorted.bam -fq ${name}.R1.nh.fastq -fq2 ${name}.R2.nh.fastq

date;
