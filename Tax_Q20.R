##### Assigning Taxonomy to the Q20 sequence table
## seqtab.nochim needs to be:
# -- Large matrix (787338 elements, 3.7 MB)
# -- Dimensions 527 samples by 1494 ASVs
## This is how you know you have the correct matrix loaded

# Import packages
library(dada2); packageVersion("dada2")

# Set WD for or local machine
setwd('~/LeBoldus/local_git/poplar_microbiome_REEU')

# Define UNITE reference file path
unite.ref <- "UNITE/sh_general_release_29.11.2022/sh_general_release_dynamic_29.11.2022.fasta"

#Assign Taxonomy
taxa.28 <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE, verbose = TRUE)

