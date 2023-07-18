#!/bin/bash

# Unzipping files after fastQC

# Make directory if it does not already exist

unzip_FQC_dir="outputs/unzip_fastQC_outputs"

if [ ! -d "$unzip_FQC_dir" ]; then
    mkdir -p "$unzip_FQC_dir"
        echo "Directory '$unzip_FQC_dir' created."
else
    echo "Directory '$unzip_FQC_dir' already exists."
fi

# Unzip

for zipfile in outputs/raw_fastQC_outputs/*.zip
do
unzip $zipfile -d outputs/unzip_fastQC_outputs
