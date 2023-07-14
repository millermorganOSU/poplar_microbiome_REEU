##### Assigning Taxonomy to the Q28 sequence table
## seqtab.nochim needs to be:
# -- Large matrix (488654 elements, 2.3 MB)
# -- Dimensions 526 samples by 929 ASVs
## This is how you know you have the correct matrix loaded

dim(seqtab.nochim)
save(seqtab.nochim, file = "outputs/Q28/SeqTabQ28.RData")

# Import packages
library(dada2); packageVersion("dada2")

# Set WD for or local machine
setwd('~/LeBoldus/local_git/poplar_microbiome_REEU')

# Define UNITE reference file path
unite.ref <- "UNITE/sh_general_release_29.11.2022/sh_general_release_dynamic_29.11.2022.fasta"

#Assign Taxonomy
taxa.28 <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE, verbose = TRUE)

