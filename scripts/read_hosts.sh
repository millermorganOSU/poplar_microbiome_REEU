#!/bin/sh

## Install the following programs
# bowtie2  # it can be in brew for mac or apt install for linux
# samtools 
# bedtools

# Remove host working
# Tutorial http://www.metagenomics.wiki/tools/samtools/remove-host-sequences
host=$1 # host index
read1=$2 # read 1
read2=$3 # read 2

# Start indexing your reference_sequence
#$ bowtie2-build Reference_sequence.fasta Species

# Aligning the reads to the host reference 
bowtie2 -x $host -1 {$read1} -2 {$read2} -S {$host}_mapped_and_unmapped.sam
samtools view -bS {$hots}_mapped_and_unmapped.sam > {$host}_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 {$host}_mapped_and_unmapped.bam > {$host}_bothEndsUnmapped.bam
samtools sort -n {$host}_bothEndsUnmapped.bam {$host}_bothEndsUnmapped_sorted
bedtools bamtofastq -i {$host}_bothEndsUnmapped_sorted.bam -fq {$read1}_r1_nh.fastq -fq2 {$read2}_r2_nh.fastq

date;
