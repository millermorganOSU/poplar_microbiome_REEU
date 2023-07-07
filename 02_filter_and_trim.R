# This script will:

## Remove primer sequences from forward and reverse reads (including reverse compliment primer sequences from read-through)
## Filter and trim reads
## Determine error rates

# Import packages
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(glue); packageVersion("glue")

# Set WD for or local machine
setwd('~/LeBoldus/local_git/poplar_microbiome_REEU')

# Define path for data and outputs
data_path <- "data/reads"
outputs_path <- "outputs"
filtered_path <- "data/filtered_outputs"
n_path <- "data/n_filtered_reads"
cut_path <- "data/cut_reads"

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

## filtered outputs
if (!dir.exists(filtered_path)) {
  dir.create(filtered_path)
}else{
  glue("Output directory {filtered_path} allready exists")
}

# Subset and sort forward and reverse reads

## basically this is reading the relative paths of all of the raw R1 files into the F list and all of the raw R2 files into the R list
fnFs <-sort(list.files(data_path, pattern=".R1.fq.gz$", full.names = TRUE))
fnRs <-sort(list.files(data_path, pattern=".R2.fq.gz$", full.names = TRUE))

# List sample names in sorted order
sample_names <- sapply(strsplit(basename(fnFs), ".R[1-2].fq.gz$"), `[`, 1)

## at this point we have 528 samples each with a F and R read for a total of 1056 read files

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
nout <- filterAndTrim(fnFs, fnFs_filtN, fnRs, fnRs_filtN, maxN = 0, multithread = TRUE)

nout

# Define function to count number of reads containing primer sequence
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Test a file to see how many primer sequences are present
print(fnFs_filtN[[33]])
print(fnRs_filtN[[33]])
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs_filtN[[33]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs_filtN[[33]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs_filtN[[33]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs_filtN[[33]]))

## only RCs found indicating that forward primers were allready removed but readthrough errors still present

# Create file names for filtered output files
filtFs <- file.path(filtered_path, paste0(sample_names, ".R1.filt.fastq.gz"))
filtRs <- file.path(filtered_path, paste0(sample_names, ".R2.filt.fastq.gz"))

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
# Cutadapt path only will work for local machine, replace path if replicating code

cutadapt <- "/Users/morganmiller/miniconda3/envs/cutadaptenv/bin/cutadapt"
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs_cut[i], "-p", fnRs_cut[i],
                             fnFs_filtN[i], fnRs_filtN[i]))
}

# Test samples to see if primer sequences were removed

sample_numbers = list(109, 44, 16, 14, 33, 199, 350)
for (i in sample_numbers) {
  print(rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs_cut[[i]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs_cut[[i]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs_cut[[i]]),
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs_cut[[i]])))
}

# All zeros as expected

# Filter and trim
?filterAndTrim()

# Max N = zero shouldn't matter (we allready filtered out Ns)
# Max EE = maximum allowed errors
# TruncQ = minimum quality allowed (sequences that dip to 28 get cut)
# minLength = sequences shorter than 90 are filtered out

fout <- filterAndTrim(fnFs_cut, filtFs, fnRs_cut, filtRs, maxN=0, 
                     maxEE=c(2,2), truncQ=28, compress=TRUE, multithread=TRUE,
                     minLen=90)

fout

#some files completely filtered out (2 samples, 4 files), started with allmost no sequences
#myco.04.12.a.R1.fq.gz        1         0
#myco.03.12.d.R1.fq.gz        5         0
#myco.04.12.a.R2.fq.gz        1         0 (?)
#myco.03.12.d.R2.fq.gz        5         0 (?)
#reverse reads were also discarded, at least their files do not appear in the output directory

# Explore quality plots
plotQualityProfile(fnRs_cut[53])
plotQualityProfile(filtRs[53])

# Rename files so that missing files don't mess with the error algorithm 
post_filter_fnFs <-sort(list.files(filtered_path, pattern=".R1.filt.fastq.gz$", full.names = TRUE))
post_filter_fnRs <-sort(list.files(filtered_path, pattern=".R2.filt.fastq.gz$", full.names = TRUE))
post_filter_sample_names <- sapply(strsplit(basename(post_filter_fnFs), ".R[1-2].filt.fastq.gz$"), `[`, 1)

# Filter out files not present in new file list
"outputs/filtered_outputs/myco.04.12.a_F_filt.fastq.gz" %in% post_filter_fnFs
"outputs/filtered_outputs/myco.03.12.d_F_filt.fastq.gz" %in% post_filter_fnFs
"outputs/filtered_outputs/myco.04.12.a_R_filt.fastq.gz" %in% post_filter_fnFs
"outputs/filtered_outputs/myco.03.12.d_R_filt.fastq.gz" %in% post_filter_fnFs

# Learn error rates
errF <- learnErrors(post_filter_fnFs, multithread=TRUE)
errR <- learnErrors(post_filter_fnRs, multithread=TRUE)

# Plot error rates
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Deprecation: in simple terms combines unique sequences and defines abundance rather than perform redundant computations
derepFs <- derepFastq(post_filter_fnFs, verbose=TRUE)
derepRs <- derepFastq(post_filter_fnRs, verbose=TRUE)

# Name the derep-class objects by the sample names rather than entire file name
names(derepFs) <- post_filter_sample_names
names(derepRs) <- post_filter_sample_names

#Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merging paired-end (FW and reverse) sequences 

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

# Stats and Distribution
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(nout, fout, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("raw", "nfilt","input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
head(track)
class(track)
