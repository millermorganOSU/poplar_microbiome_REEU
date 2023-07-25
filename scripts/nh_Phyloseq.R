## Phyloseq

# Import packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(readxl); packageVersion("readxl")
library(reshape2); packageVersion("reshape2")
library(seqinr); packageVersion("seqinr")
library(tidyverse)
library(glue)
library(vegan)

# Set WD for or local machine
setwd('~/LeBoldus/local_git/poplar_microbiome_REEU')

# Load sequence table
load("R/SeqTab_nh.RData")

# Load complete Taxa table
load("R/FilledTaxa.RData")

# Load metadata
disease.scores <- read_xlsx("data/sample_diseasescores.xlsx")
wood.meta <- read.csv("data/wood_meta2.csv", row.names = 1)

##### Trimming data and working with Phyloseq #####

# Separate out fungi
fungal.taxa <- as.matrix(taxa.V2[taxa.V2[,1]=="k__Fungi",])
dim(fungal.taxa)

# 520 samples with 613 ASVs as columns
dim(seqtab.nochim)
class(seqtab.nochim)

# List ASVs
all.ASVs <- as.list(colnames(seqtab.nochim))
fungal.ASVs <- as.list(rownames(fungal.taxa))

# Subset columns (ONLY FUNGAL ASVs)
seqtab.nochim.fungal <- seqtab.nochim[,unlist(fungal.ASVs)]
dim(seqtab.nochim.fungal)

# Test
test.names.ASVs <- as.list(colnames(seqtab.nochim.fungal))
identical(test.names.ASVs, fungal.ASVs) # Worked!

# 520 samples with 470 fungal ASVs as columns
dim(seqtab.nochim.fungal)
class(seqtab.nochim.fungal)

# List samples with metadata
sample.names.wmeta <- as.list(rownames(wood.meta))

# Subset rows (ONLY SAMPLES WITH METADATA)
seqtab.nochim.fungal.wmeta <- seqtab.nochim.fungal[unlist(sample.names.wmeta),]
dim(seqtab.nochim.fungal.wmeta)

# 236 samples with 470 fungal ASVs as columns
dim(seqtab.nochim.fungal.wmeta)
class(seqtab.nochim.fungal.wmeta)

# Test
test.names.samples <- as.list(rownames(seqtab.nochim.fungal.wmeta))
identical(test.names.samples, sample.names.wmeta) # Worked!

# Remove zero value columns and rows

zero.col <- as.matrix(which(colSums(seqtab.nochim.fungal.wmeta)==0))
dim(zero.col)
non.zero.col <- as.matrix(which(colSums(seqtab.nochim.fungal.wmeta)!=0))
dim(non.zero.col)
zero.row <- as.matrix(which(rowSums(seqtab.nochim.fungal.wmeta)==0))
dim(zero.row)
non.zero.row <- as.matrix(which(rowSums(seqtab.nochim.fungal.wmeta)!=0))
dim(non.zero.row)

non.zero.col.keep <- as.list(non.zero.col)
non.zero.row.keep <- as.list(non.zero.row)

clean.seq.tab <- seqtab.nochim.fungal.wmeta[,unlist(non.zero.col.keep)]
dim(clean.seq.tab)
clean.seq.tab <- clean.seq.tab[unlist(non.zero.row.keep),]
dim(clean.seq.tab)

# Now the seq tab has
# - only fungal ASVs
# - only samples that have metdata
# - only samples and ASVs with abundance > 0
# 235 samples by 359 ASVs

# 236 samples with 10 metadata columns
dim(wood.meta)
class(wood.meta)

# Remove sample with no fungal ASVs
wood.meta.clean <- wood.meta[-38,]

# 470 ASVs with taxa (and index and accession)
dim(fungal.taxa)
class(fungal.taxa)

# Keep only ASVs that have abundance in seq table
non.zero.ASVs <- as.list(colnames(clean.seq.tab))

fungal.taxa.clean <- fungal.taxa[unlist(non.zero.ASVs),]
fungal.taxa.clean <- fungal.taxa[,1:7]
dim(fungal.taxa.clean)

# ### Add total occurrence
# # The code below was a test that wasn't very useful
# 
# identical(colnames(seqtab.nochim.fungal.wmeta), rownames(fungal.taxa)) # Seemingly ordered the same
# 
# # Create column of total occurrence
# ASV.overall.dist <- as.matrix(colSums(seqtab.nochim.fungal.wmeta))
# ASV.overall.dist[,1] <- as.character(ASV.overall.dist[,1])
# class(ASV.overall.dist[,1])
# colnames(ASV.overall.dist) <- "Total"
# dim(ASV.overall.dist)
# 
# # Create data column of ASVs
# fungal.ASVs.mtrx <- as.matrix(unlist(fungal.ASVs))
# fungal.ASVs.mtrx[,1] <- as.character(fungal.ASVs.mtrx[,1])
# class(fungal.ASVs.mtrx[,1])
# colnames(fungal.ASVs.mtrx) <- "ASV"
# dim(fungal.ASVs.mtrx)
# 
# # Create new ID matrix for ASVs
# id.mtrx <- as.matrix(as.character(seq(1,470)))
# id.mtrx[,1] <- as.character(id.mtrx[,1])
# class(id.mtrx[,1])
# colnames(id.mtrx) <- "id"
# dim(id.mtrx)
# 
# # Character matrix for phyloseq
# fungal.taxa.wtotals.cbind <- cbind(fungal.taxa, ASV.overall.dist)
# fungal.taxa.wtotals.cbind <- cbind(fungal.taxa.wtotals.cbind, fungal.ASVs.mtrx)
# fungal.taxa.wtotals.cbind <- cbind(fungal.taxa.wtotals.cbind, id.mtrx)
# dim(fungal.taxa.wtotals.cbind)
# 
# # Test Totals
# test1 <- "TTAAGTTCAGCGGGTAGTCCTACCTGATTTGAGGTCAATGATAAAAGATGGGGGTTGTGAGCAGGCATCACACCAGGGCTTGACGAAACTTATTACGTCCAACCGGTGCTTATCCCACTGACGCTTTTAAGGTGAGCCGAGCGGCAGCACCCAAGTCCAACCTCCCCAGATCAGAAACCTGAGGGGTTGAGATTACATGACACTCAAACAGGCATGCCTTTCGGAATACCAAAAGGCGCAAGGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCG"
# ASV.overall.dist[test1,]
# fungal.taxa.wtotals.cbind[test1,10]
# test2 <- "ATAAGTTCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTAGAAAAATAAAGGTTTCAGTCGGCAGAGTTCCTCTCCTTTGACAGACGTTCGAATAAATTCTACTACGCCTAAAGCCGGAGTGGCCTCGCCGAGGTCTTTAAGGCGCGCCCAACTAAGGACGACGCCCAATACCAAGCATAGCTTGAGTGGTGTAATGACGCTCGAACAGGCATGCCCCTCGGAATACCAAGGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCG"
# ASV.overall.dist[test2,]
# fungal.taxa.wtotals.cbind[test2,10]
# test3 <- "TTAAGTTCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTTGAGGTGTTGGGTTTTGGAGGCGGGCGACCACAGACTCTAGAAACGGGAAGTATTACTAACGCTTAAGAGGCCTGAGCCACCGCCGACCGGTTTGAAGCGCGCCCACCCGGTGAGGGGGAGGTGACGCTCAATGCCAAGCAGAGCTTGATGGTTGATAATGACGCTCGAACAGGCATGCCCACCGGAATACCAGTGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATAACTTATCGCATTTCG"
# ASV.overall.dist[test3,]
# fungal.taxa.wtotals.cbind[test3,10]
# 
# # Test ASVs
# test4 <- "TTAAGTTCAGCGGGTAGTCCTACCTGATTTGAGGTCAATGATAAAAGATGGGGGTTGTGAGCAGGCATCACACCAGGGCTTGACGAAACTTATTACGTCCAACCGGTGCTTATCCCACTGACGCTTTTAAGGTGAGCCGAGCGGCAGCACCCAAGTCCAACCTCCCCAGATCAGAAACCTGAGGGGTTGAGATTACATGACACTCAAACAGGCATGCCTTTCGGAATACCAAAAGGCGCAAGGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCG"
# identical(unlist(fungal.taxa.wtotals.cbind[test4,11]), test4)
# test5 <- "ATAAGTTCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTAGAAAAATAAAGGTTTCAGTCGGCAGAGTTCCTCTCCTTTGACAGACGTTCGAATAAATTCTACTACGCCTAAAGCCGGAGTGGCCTCGCCGAGGTCTTTAAGGCGCGCCCAACTAAGGACGACGCCCAATACCAAGCATAGCTTGAGTGGTGTAATGACGCTCGAACAGGCATGCCCCTCGGAATACCAAGGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCG"
# identical(unlist(fungal.taxa.wtotals.cbind[test5,11]), test5)
# test6 <- "TTAAGTTCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTTGAGGTGTTGGGTTTTGGAGGCGGGCGACCACAGACTCTAGAAACGGGAAGTATTACTAACGCTTAAGAGGCCTGAGCCACCGCCGACCGGTTTGAAGCGCGCCCACCCGGTGAGGGGGAGGTGACGCTCAATGCCAAGCAGAGCTTGATGGTTGATAATGACGCTCGAACAGGCATGCCCACCGGAATACCAGTGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATAACTTATCGCATTTCG"
# identical(unlist(fungal.taxa.wtotals.cbind[test6,11]), test6)
# 
# Examine ASVs with the most abundance (single sample)
# 
# high_rep <- as.list(tail(sort(seqtab.nochim.fungal.wmeta),10))
# rc <- matrix(, nrow = 0, ncol = 17)
# for ( i in high_rep) {
#   add <- which(seqtab.nochim.fungal.wmeta == i, arr.ind = TRUE)
#   row <- add[1,1]
#   col <- add[1,2]
#   sample <- as.vector(test.names.samples[row])
#   asv <- as.vector(test.names.ASVs[col])
#   value <- as.vector(i)
#   lin <- t(as.matrix(fungal.taxa.wtotals.cbind[unlist(asv),]))
#   add <- cbind(add, asv)
#   add <- cbind(add, sample)
#   add <- cbind(add, value)
#   add <- cbind(add, lin)
#   rc <- rbind(rc, add)
# }
# 
# # Another test
# for (i in 1:(nrow(rc))) {
#   print(identical(rc[i,3], rc[i,16]))
# }

# # Create phyloseq object
# ps <- phyloseq(otu_table(seqtab.nochim.fungal.wmeta, taxa_are_rows=FALSE), sample_data(wood.meta), tax_table(fungal.taxa.no.index))

# Create phyloseq object clean
ps.clean <- phyloseq(otu_table(clean.seq.tab, taxa_are_rows=FALSE), sample_data(wood.meta.clean), tax_table(fungal.taxa.clean))

# Plot ASV abundance per sample

ASV.per.sample.clean <- plot_bar(ps.clean, fill = "Abundance") +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  scale_fill_gradient(low="lightgoldenrodyellow", high="red", name = "Individual ASV \n Abundance") +
  labs(title = "ASV Abundance per Sample") +
  xlab("Sample") +
  ylab("ASV Abundance")
ASV.per.sample.clean

# Plot ASV abundance per sample (sorted)

ASV.per.sample.sorted.clean <- plot_bar(ps.clean, fill = "Abundance") +
  aes(reorder(Sample, -Abundance)) +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  scale_fill_gradient(low="lightgoldenrodyellow", high="red", name = "Individual ASV \n Abundance") +
  labs(title = "ASV Abundance per Sample") +
  xlab("Sample") +
  ylab("ASV Abundance")
ASV.per.sample.sorted.clean

# Plot ordination NMDS
ps.ord.nmds <- ordinate(ps.clean, "NMDS", "bray")
ASV.ordination.nmds <- plot_ordination(ps.clean, ps.ord.nmds, type="taxa", color="Phylum", title="taxa")
ASV.ordination.nmds

# Plot PCoA
ps.ord.pcoa <- ordinate(ps.clean, "PCoA", "bray")
ASV.ordination.pcoa <- plot_ordination(ps.clean, ps.ord.pcoa, type="taxa", color="Phylum", title="taxa")
ASV.ordination.pcoa

#### ordination not working with phyloseq, trying with vegan #####

ntaxa(ps.clean)
MDS <- metaMDS(clean.seq.tab, distance = "bray", weakties = FALSE)
MDS
stressplot(MDS)
dist_bc <- distance(ps.clean, method = "bray")
stressplot(dist_bc)
