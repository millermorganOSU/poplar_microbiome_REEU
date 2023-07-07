#!/bin/bash

# Make directory in output directory only if it does not already exist

echo "Make Sure FastQC is downloaded on your machine"

postfilter_FQC_dir="outputs/postfilter_fastQC_outputs"

if [ ! -d "$postfilter_FQC_dir" ]; then
    mkdir -p "$postfilter_FQC_dir"
        echo "Directory '$postfilter_FQC_dir' created."
else
    echo "Directory '$postfilter_FQC_dir' already exists."
fi

echo "FastQC-ing filtered reads"

fastqc outputs/filtered_outputs/myco*gz -o outputs/postfilter_fastQC_outputs

