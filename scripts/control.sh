#!/bin/sh

# This script will control align.sh

# Static variables

host_i=reference/Ptrichocarpa # host index
host_o_s=sam/Ptrichocarpa # sam output
host_o_b=bam/Ptrichocarpa # bam output 

# Create names text file

rm inames.txt
rm names.txt

for FILE in ../data/filtered_outputs/* ; do
echo "${FILE:25:12}" >> inames.txt
done

sort -u inames.txt >> names.txt

rm inames.txt

# Loop through filtered reads

while read name; do
read1_i="../data/filtered_outputs/${name}.R1.filt.fastq.gz"
read2_i="../data/filtered_outputs/${name}.R2.filt.fastq.gz"
# echo "${read1_i}"
# echo "${read2_i}"
./align.sh ${host_i} ${read1_i} ${read2_i} ${host_o_s} ${host_o_b} nh/${name}
done <names.txt

# After running some files contained zero length sequences
# This causes issued with dada2

#nh/myco.04.01.h.R1.nh.fastq
#nh/myco.04.01.h.R2.nh.fastq
#nh/myco.04.06.f.R1.nh.fastq
#nh/myco.04.06.f.R2.nh.fastq
#nh/myco.05.01.h.R1.nh.fastq
#nh/myco.05.01.h.R2.nh.fastq
#nh/myco.06.06.d.R1.nh.fastq
#nh/myco.06.06.d.R2.nh.fastq
#nh/myco.06.09.h.R1.nh.fastq
#nh/myco.06.09.h.R2.nh.fastq
#nh/myco.07.05.b.R1.nh.fastq
#nh/myco.07.05.b.R2.nh.fastq

# removed using
# $find nh -size 0c -delete


