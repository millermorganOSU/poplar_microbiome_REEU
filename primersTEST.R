# This script will:

## Remove primer sequences from forward and reverse reads using cutadapt
## ---(primarily reverse compliment primer sequences from read-through)
## Filter and trim reads ##### THIS SCRIPT IS FOR TRUNCQ = 28 #####
## Determine and plot error rates
## Merge ASVs into sequence table
## Remove chimeras
## Create tracking violin plots

# Import packages
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(glue); packageVersion("glue")
library(reshape2); packageVersion("reshape2")
library(ggplot2); packageVersion("ggplot2")

# Set WD for or local machine
setwd('~/LeBoldus/local_git/poplar_microbiome_REEU')

# Define path for data and outputs
data_path <- "data/reads"
outputs_path <- "outputs"
filtered_path <- "data/filtered_outputs"
filtered_path_Q20 <-"data/filtered_outputs_Q20"
n_path <- "data/n_filtered_reads"
cut_path <- "data/cut_reads"
Q28_path <- "outputs/Q28"
Q20_path <- "outputs/Q20"

# Create directories if they do not all ready exist

## outputs dir
if (!dir.exists(outputs_path)) {
  dir.create(outputs_path)
}else{
  glue("Output directory {outputs_path} allready exists")
}

## n filtered reads
if (!dir.exists(n_path)) {
  dir.create(n_path)
}else{
  glue("Output directory {n_path} allready exists")
}

## cut reads
if (!dir.exists(cut_path)) {
  dir.create(cut_path)
}else{
  glue("Output directory {cut_path} allready exists")
}

## filtered outputs (these ones are the Q28)
if (!dir.exists(filtered_path)) {
  dir.create(filtered_path)
}else{
  glue("Output directory {filtered_path} allready exists")
}

## filtered outputs (these will be Q20)
if (!dir.exists(filtered_path_Q20 )) {
  dir.create(filtered_path_Q20 )
}else{
  glue("Output directory {filtered_path_Q20 } allready exists")
}

## summary outputs as we tweak filter (tight filter: trunQ = 28)
if (!dir.exists(Q28_path)) {
  dir.create(Q28_path)
}else{
  glue("Output directory {Q28_path} allready exists")
}

## summary outputs as we tweak filter (lax filter: trunQ = 20)
if (!dir.exists(Q20_path)) {
  dir.create(Q20_path)
}else{
  glue("Output directory {Q20_path} allready exists")
}

# Subset and sort forward and reverse reads

## read the relative paths of all of the raw files
## sort so they are in order
## separate into forward and reverse
fnFs <-sort(list.files(data_path, pattern=".R1.fq.gz$", full.names = TRUE))
fnRs <-sort(list.files(data_path, pattern=".R2.fq.gz$", full.names = TRUE))

# List sample names in sorted order
sample_names <- sapply(strsplit(basename(fnFs), ".R[1-2].fq.gz$"), `[`, 1)

## at this point we have 528 samples each with a F and R read for a total of 1056 files

# Assign primer variables with correct strings
ITS4 <- "TCCTCCGCTTATTGATATGC"
REV <- ITS4
REV.RC <- dada2:::rc(REV)
ITS3_KYO1 <- "AHCGATGAAGAACRYAG"
FWD <- ITS3_KYO1
FWD.RC <- dada2:::rc(FWD)

# Define function that creates orients of primers
allOrients <- function(primer) {
  require(Biostrings) ## require() just outputs a warning if package is not installed
  dna <- DNAString(primer) ## turns input sequence (string) into DNAString object
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna)) ## create orients by manipulating DNAString
  return(sapply(orients, toString))  ## convert back to string
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Pre-filter out Ns so that primer sequences can be mapped effectively
fnFs_filtN <- file.path(n_path, basename(fnFs))
fnRs_filtN <- file.path(n_path, basename(fnRs))

####### commented out so as to not accidentally run again
# nout <- filterAndTrim(fnFs, fnFs_filtN, fnRs, fnRs_filtN, maxN = 0, multithread = TRUE)
# nout
####### commented out so as to not accidentally run again

# Define function to count number of reads containing primer sequence
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Test files to see how many primer sequences are present
sample_numbers = seq(1, 528, 25)
for (x in sample_numbers) {
  capture.output(
    print(
      rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs_filtN[[x]]), 
            FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs_filtN[[x]]), 
            REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs_filtN[[x]]),
            REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs_filtN[[x]]))
    ), file = "outputs/withPrimersTEST.txt", append = TRUE)
}

## mainly reverse compliments present 
## ---indicates forward primers were already removed 
## ---read through errors still present

# Create file names for cut files
fnFs_cut <- file.path(cut_path, basename(fnFs))
fnRs_cut <- file.path(cut_path, basename(fnRs))

# Create arguments for cutadapt

## capital argument means effective for the reverse read
## -a means regular three prime adapter
## -g means regular five prime adapter
## in English: looking for forward and reverse RC on both the 3' and 5' ends of the forward reads
## in English: looking for reverse and forward RC on both the 3' and 5' ends of the reverse reads
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

# Run cutadapt from R script
# cutadapt path only will work for local machine, replace path if replicating code

####### commented out so as to not accidentally run again
# cutadapt <- "/Users/morganmiller/miniconda3/envs/cutadaptenv/bin/cutadapt"
# for(i in seq_along(fnFs)) {
#   system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
#                              "-o", fnFs_cut[i], "-p", fnRs_cut[i],
#                              fnFs_filtN[i], fnRs_filtN[i]))
# }
####### commented out so as to not accidentally run again

# Test samples to see if primer sequences were removed
sample_numbers = seq(1, 528, 25)
for (i in sample_numbers) {
  capture.output(
    print(
      rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs_cut[[i]]), 
            FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs_cut[[i]]), 
            REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs_cut[[i]]),
            REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs_cut[[i]]))
    ), file = "outputs/noPrimers.txt", append = TRUE)
}