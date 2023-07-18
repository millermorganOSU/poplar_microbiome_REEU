
# This script will:

## Determine and plot error rates ##### AFTER HOST SEQUENCE REMOVAL #####
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
nh_path <- "genome/nh"
outputs_path <- "outputs"
nh_o_path <- "outputs/nh"

# Create directories if they do not all ready exist

## outputs directory
if (!dir.exists(outputs_path)) {
  dir.create(outputs_path)
}else{
  glue("Output directory {outputs_path} allready exists")
}

## nh outputs directory
if (!dir.exists(nh_o_path)) {
  dir.create(nh_o_path)
}else{
  glue("Output directory {nh_o_path} allready exists")
}

# Rename files so that missing files (filtered out) don't mess with the error algorithm 
post_filter_fnFs <-sort(list.files(nh_path, pattern=".R1.nh.fastq$", full.names = TRUE))
post_filter_fnRs <-sort(list.files(nh_path, pattern=".R2.nh.fastq$", full.names = TRUE))
post_filter_sample_names <- sapply(strsplit(basename(post_filter_fnFs), ".R[1-2].nh.fastq$"), `[`, 1)

post_filter_fnFs

# Filtered out files not present in new file list
"genome/nh/myco.04.12.a.R1.nh.fastq" %in% post_filter_fnFs # should be F
"genome/nh/myco.03.12.d.R2.nh.fastq" %in% post_filter_fnRs # should be F
"genome/nh/myco.04.12.a.R1.nh.fastq" %in% post_filter_fnFs # should be F
"genome/nh/myco.03.12.d.R2.nh.fastq" %in% post_filter_fnRs # should be F
"genome/nh/myco.04.01.h.R1.nh.fastq" %in% post_filter_fnFs # should be F
"genome/nh/myco.04.01.h.R2.nh.fastq" %in% post_filter_fnRs # should be F
"genome/nh/myco.04.06.f.R1.nh.fastq" %in% post_filter_fnFs # should be F
"genome/nh/myco.04.06.f.R2.nh.fastq" %in% post_filter_fnRs # should be F
"genome/nh/myco.05.01.h.R1.nh.fastq" %in% post_filter_fnFs # should be F
"genome/nh/myco.05.01.h.R2.nh.fastq" %in% post_filter_fnRs # should be F
"genome/nh/myco.06.06.d.R1.nh.fastq" %in% post_filter_fnFs # should be F
"genome/nh/myco.06.06.d.R2.nh.fastq" %in% post_filter_fnRs # should be F
"genome/nh/myco.06.09.h.R1.nh.fastq" %in% post_filter_fnFs # should be F
"genome/nh/myco.06.09.h.R2.nh.fastq" %in% post_filter_fnRs # should be F
"genome/nh/myco.07.05.b.R1.nh.fastq" %in% post_filter_fnFs # should be F
"genome/nh/myco.07.05.b.R2.nh.fastq" %in% post_filter_fnRs # should be F
"genome/nh/myco.03.06.a.R2.nh.fastq" %in% post_filter_fnRs # should be T
"genome/nh/myco.03.06.a.R1.nh.fastq" %in% post_filter_fnFs # should be T

# Learn error rates
errF <- learnErrors(post_filter_fnFs, multithread=TRUE)
errR <- learnErrors(post_filter_fnRs, multithread=TRUE)

# Plot error rates and create image files
err_F_plot <- plotErrors(errF, nominalQ = TRUE)
jpeg(filename="outputs/nh/Error_F_plots_nh.jpeg", width = 4, height = 4, units = 'in', res = 300)
err_F_plot
dev.off()
err_R_plot <- plotErrors(errR, nominalQ = TRUE)
jpeg(filename="outputs/nh/Error_R_plots_nh.jpeg", width = 4, height = 4, units = 'in', res = 300)
err_R_plot
dev.off()

# Derep combines unique sequences and defines abundance rather than perform redundant computations

####### commented out so as to not accidentally run again
derepFs <- derepFastq(post_filter_fnFs, verbose=TRUE)
derepRs <- derepFastq(post_filter_fnRs, verbose=TRUE)
####### commented out so as to not accidentally run again

# Name the derep-class objects by the sample names rather than entire file name
names(derepFs) <- post_filter_sample_names
names(derepRs) <- post_filter_sample_names

#Sample inference

####### commented out so as to not accidentally run again
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
####### commented out so as to not accidentally run again

# Merging paired-end (forward and reverse) sequences 

##objects are a list of dfs for each sample with sequence and abundance
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Make sequence table, a single matrix with samples as rows and ASVs as columns
seqtab <- makeSequenceTable(mergers)

# 526 samples made it through the filter and 1074 ASVs were detected
dim(seqtab)

# Wide range of seq length distribution, spikes at 303bp
seq.dist.tab <- table(nchar(getSequences(seqtab)))

# Make distribution plot

# Make data frame with Length_pb field
seq.dist.tab.df <- as.data.frame(seq.dist.tab)
colnames(seq.dist.tab.df) <- c("Length_bp", "Frequency")

# Make Length_pb index
index <- as.data.frame(seq(100,450))
colnames(index) <- "Length_bp"

# Merge so independent is continuous
seq.dist.tab.df <- merge(index, seq.dist.tab.df, by = "Length_bp", all.x = TRUE, all.y = TRUE)

# Simple plot
nh.dist.plot <- ggplot(seq.dist.tab.df, aes(x=Length_bp, y=Frequency)) +
  geom_point(colour = "red", size = 3) +
  labs(title = "Sequence Length Distribution (host removal pipeline)", xlab = "Length (bp)", ylab = "Frequency") +
  ggsave("outputs/nh/Sequence_Length_Distribution_nh.png")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Post chimera removal stats and distribution
head(seqtab.nochim)
class(seqtab.nochim)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #chimeric sequence frequency
table(nchar(getSequences(seqtab.nochim)))

# Save sequence table as Rdata
save(seqtab.nochim, file = "outputs/nh/SeqTab_nh.RData")