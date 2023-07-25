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
done <namesTEST.txt




