#!/bin/bash

# Make directory in output directory only if it does not already exist

echo "Make Sure FastQC is downloaded on your machine"

raw_FQC_dir="outputs/raw_fastQC_outputs"

if [ ! -d "$raw_FQC_dir" ]; then
    mkdir -p "$raw_FQC_dir"
        echo "Directory '$raw_FQC_dir' created."
else
    echo "Directory '$raw_FQC_dir' already exists."
fi

echo "FastQC-ing all raw reads"

fastqc data/reads/myco*gz -o outputs/raw_fastQC_outputs

