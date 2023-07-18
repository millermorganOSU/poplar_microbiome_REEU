# This script will:

## Remove primer sequences from forward and reverse reads using cutadapt
## ---(primarily reverse compliment primer sequences from read-through)
## Filter and trim reads ##### THIS SCRIPT IS FOR TRUNCQ = 20 #####
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
FWD <- ITS4
FWD.RC <- dada2:::rc(FWD)
ITS3_KYO1 <- "AHCGATGAAGAACRYAG"
REV <- ITS3_KYO1
REV.RC <- dada2:::rc(REV)

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
      ), file = "outputs/withPrimers.txt", append = TRUE)
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

# All zeros as expected

# Filter and trim
?filterAndTrim()

# Create file names for filtered output files
filtFs <- file.path(filtered_path_Q20, paste0(sample_names, ".R1.filt.fastq.gz"))
filtRs <- file.path(filtered_path_Q20, paste0(sample_names, ".R2.filt.fastq.gz"))

# Max N = zero shouldn't matter (we allready filtered out Ns)
# Max EE = maximum allowed errors
# TruncQ = minimum quality allowed (sequences that dip to 28 get cut)
# minLength = sequences shorter than 90 are filtered out

####### commented out so as to not accidentally run again
# fout <- filterAndTrim(fnFs_cut, filtFs, fnRs_cut, filtRs, maxN=0, 
#                      maxEE=c(2,2), truncQ=20, compress=TRUE, multithread=TRUE,
#                      minLen=90)
# fout
####### commented out so as to not accidentally run again

# Check which sample got filtered out
fout_df <- as.data.frame(fout)
Q20_filtered_out <- subset(fout_df, reads.out == 0)
Q20_filtered_out

# One file completely filtered out (1 sample, 2 files), started with only one sequence
#myco.04.12.a.R1.fq.gz        1         0
#myco.04.12.a.R2.fq.gz        1         0 (?)

# Explore quality plots
plotQualityProfile(fnRs_cut[53])
plotQualityProfile(filtRs[53])

# Rename files so that missing files don't mess with the error algorithm 
post_filter_fnFs <-sort(list.files(filtered_path_Q20, pattern=".R1.filt.fastq.gz$", full.names = TRUE))
post_filter_fnRs <-sort(list.files(filtered_path_Q20, pattern=".R2.filt.fastq.gz$", full.names = TRUE))
post_filter_sample_names <- sapply(strsplit(basename(post_filter_fnFs), ".R[1-2].filt.fastq.gz$"), `[`, 1)

# Filter out files not present in new file list
"data/filtered_outputs_Q20/myco.04.12.a.R1.filt.fastq.gz" %in% post_filter_fnFs # should be F
"data/filtered_outputs_Q20/myco.03.12.d.R2.filt.fastq.gz" %in% post_filter_fnRs # should be T
"data/filtered_outputs_Q20/myco.04.12.a.R1.filt.fastq.gz" %in% post_filter_fnFs # should be F
"data/filtered_outputs_Q20/myco.03.12.d.R2.filt.fastq.gz" %in% post_filter_fnRs # should be T

# Learn error rates
errF <- learnErrors(post_filter_fnFs, multithread=TRUE)
errR <- learnErrors(post_filter_fnRs, multithread=TRUE)

# Plot error rates and create image files
err_F_plot <- plotErrors(errF, nominalQ = TRUE)
jpeg(filename="outputs/Q20/Error_F_plots_Q20.jpeg", width = 4, height = 4, units = 'in', res = 300)
err_F_plot
dev.off()
err_R_plot <- plotErrors(errR, nominalQ = TRUE)
jpeg(filename="outputs/Q20/Error_R_plots_Q20.jpeg", width = 4, height = 4, units = 'in', res = 300)
err_R_plot
dev.off()

# Deprecation: combines unique sequences and defines abundance rather than perform redundant computations

####### commented out so as to not accidentally run again
# derepFs <- derepFastq(post_filter_fnFs, verbose=TRUE)
# derepRs <- derepFastq(post_filter_fnRs, verbose=TRUE)
####### commented out so as to not accidentally run again

# Name the derep-class objects by the sample names rather than entire file name
names(derepFs) <- post_filter_sample_names
names(derepRs) <- post_filter_sample_names

#Sample inference

####### commented out so as to not accidentally run again
# dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
# dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
####### commented out so as to not accidentally run again

# Merging paired-end (forward and reverse) sequences 

##objects are a list of dfs for each sample with sequence and abundance
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Make sequence table, a single matrix with samples as rows and ASVs as columns
seqtab <- makeSequenceTable(mergers)

# 526 samples made it through the filter and 1074 ASVs were detected
dim(seqtab)

# Wide range of seq length distribution, spikes at 361bp
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Post chimera removal stats and distribution
head(seqtab.nochim)
class(seqtab.nochim)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #chimeric sequence frequency
table(nchar(getSequences(seqtab.nochim)))

# Make sequence table a data frame and export as CSV
seqtab.nochim.df.20 <- as.data.frame(seqtab.nochim)
write.csv(seqtab.nochim.df.20, "outputs/Q20/SequenceTableQ20.csv", row.names=TRUE)

# Tracking table
options(max.print=100000)

## row 92 myco.03.12.d
## row 185 myco.04.12.a

# Sample that was filtered out needs to be removed from filter output
head(fout,185)
track_fout <- fout[-c(185),]

# Define function that counts unique sequences (?)
getN <- function(x) sum(getUniques(x))

# Create tracking table
track <- cbind(track_fout, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("nFilt", "qFilt", "DenoisedF", "DenoisedR", "MergedASV", "remChim")
rownames(track) <- post_filter_sample_names

# As data frame
track_df <- as.data.frame(track)

# Put head 20 of tracking table in stats file
Q20_track_head <- capture.output(head(track_df,100), file = "outputs/Q20/Q20_stats.txt", append = TRUE)
Q20_track_head

# Melt data frame for violin plots
track_df_melt <-reshape2::melt(track_df, id.vars = NULL)
track_df_melt

# Create violin plot
track_plot <- ggplot(track_df_melt, aes(x=variable, y=value)) + 
  geom_violin(trim = FALSE) +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  stat_summary(fun=mean, geom="point", size=2, color="blue") +
  xlab("Process Step") +
  ylab("Number of Sequences") +
  ggtitle("Tracking Sequences Through Filtering Process")

# Export plot as jpeg
jpeg(filename="outputs/Q20/Tracking_Q20.jpeg", width = 8, height = 4, units = 'in', res = 300)
track_plot
dev.off()