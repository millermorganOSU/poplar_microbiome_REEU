## Process Taxonomy

# Import packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(readxl); packageVersion("readxl")
library(reshape2); packageVersion("reshape2")
library(seqinr); packageVersion("seqinr")

# Set WD for or local machine
setwd('~/LeBoldus/local_git/poplar_microbiome_REEU')

# Load assigned taxonomy
load("R/AssignedTaxa_nh.RData")

# Load sequence table
load("R/SeqTab_nh.RData")

# Load metadata
disease.scores <- read_xlsx("data/sample_diseasescores.xlsx")
wood.meta <- read.csv("data/wood_meta2.csv")

# Some ASVs were not assigned beyond the Fungi Kingdom
# -indicating that they might not be Fungi at all

# Index taxa table
taxa.nh.i <- cbind(taxa.nh, seq(1,613))
colnames(taxa.nh.i) <- c("Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species", "i")

# Separate with Phylum
taxa.nh.phy <- taxa.nh.i[!is.na(taxa.nh.i[,2]),]

# Separate with out Phylum
taxa.nh.nophy <- taxa.nh.i[is.na(taxa.nh.i[,2]),]

# Dimensions
dim(taxa.nh.i)
dim(taxa.nh.phy)
dim(taxa.nh.nophy)

# ASV stats
tnh <- nrow(taxa.nh.i) # 613 ASVs including ones with no Phylum
pnh <- nrow(taxa.nh.phy) # 516 ASVs with Phylum
nnh <- nrow(taxa.nh.nophy) # 97 ASVs with no Phylum

# About 16 percent of out data has no Phylum assignment
percent.nophy.nh <- ((nnh/tnh)*100)

# Separate out ASVs with no Phylum and export to .fasta
nophy.ASVs <- as.list(rownames(taxa.nh.nophy))
length(nophy.ASVs)
class(nophy.ASVs)

nophy.index <- taxa.nh.nophy[,8]
length(nophy.index)
class(nophy.index)

write.fasta(nophy.ASVs, nophy.index, "outputs/nophy_ASVs.fasta", open = "w", as.string = TRUE)
