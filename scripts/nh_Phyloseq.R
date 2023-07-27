## Phyloseq

# Import packages
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion("microbiome")
library(ggplot2); packageVersion("ggplot2")
library(readxl); packageVersion("readxl")
library(reshape2); packageVersion("reshape2")
library(seqinr); packageVersion("seqinr")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(glue); packageVersion("glue")
library(vegan); packageVersion("vegan")

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

# Separate out fungi from taxa table
fungal.taxa <- as.matrix(taxa.V2[taxa.V2[,1]=="k__Fungi",])
# 470 fungal ASVs
dim(fungal.taxa)

# 520 samples with 613 ASVs
dim(seqtab.nochim)

# List ASVs
all.ASVs <- as.list(colnames(seqtab.nochim))
fungal.ASVs <- as.list(rownames(fungal.taxa))

# Subset columns (ONLY FUNGAL ASVs)
seqtab.nochim.fungal <- seqtab.nochim[,unlist(fungal.ASVs)]
dim(seqtab.nochim.fungal)

# Test
test.names.ASVs <- as.list(colnames(seqtab.nochim.fungal))
identical(test.names.ASVs, fungal.ASVs) # Worked!

# 520 samples with 470 FUNGAL ASVs
dim(seqtab.nochim.fungal)

# List samples with metadata
sample.names.wmeta <- as.list(rownames(wood.meta))

# Subset rows (ONLY SAMPLES WITH METADATA)
seqtab.nochim.fungal.wmeta <- seqtab.nochim.fungal[unlist(sample.names.wmeta),]

# 236 samples with 470 fungal ASVs as columns
dim(seqtab.nochim.fungal.wmeta)

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

# Row 104 could be messing with NMDS
hist(clean.seq.tab[104,])
sum(clean.seq.tab[104,])

# Remove it
clean.seq.tab <- clean.seq.tab[-104,]
dim(clean.seq.tab)

# # test
# test <- rownames(clean.seq.tab)
# which(test == "myco.04.05.f") #83
# which(test == "myco.04.05.g") #84
# hist(clean.seq.tab[83,])
# sum(clean.seq.tab[83,])
# hist(clean.seq.tab[84,])
# sum(clean.seq.tab[84,])
# clean.seq.tab <- clean.seq.tab[-(83:84),]
# dim(clean.seq.tab)

# Any other low sum rows?
rowsum <- as.matrix(rowSums(clean.seq.tab))
hist(rowsum)
lowsumR <- as.matrix(tail(sort(rowSums(clean.seq.tab)), n = 20))
hist(lowsumR)

# Any other low sum cols?
colsum <- as.matrix(colSums(clean.seq.tab))
hist(colsum)
lowsumC <- as.matrix(tail(sort(colSums(clean.seq.tab)), n = 20))
hist(lowsumC)

# Looks ok to me - no unreasonable outliers

# Now the seq tab has
# - only fungal ASVs
# - only samples that have metdata
# - only samples and ASVs with abundance > 0
# - removed potential outlier row 104
# 234 samples by 359 ASVs

# 236 samples with 10 metadata columns
dim(wood.meta)
class(wood.meta)

# Remove sample with no fungal ASVs
wood.meta.clean <- wood.meta[-38,]
# Remove sample with virtually no fungal ASVs
wood.meta.clean <- wood.meta.clean[-104,]
dim(wood.meta.clean)

# # test
#wood.meta.clean <- wood.meta.clean[-(83:84),]

# 470 ASVs with taxa (and index and accession)
dim(fungal.taxa)

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

ASV.per.sample <- plot_bar(ps.clean, fill = "Abundance") +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  scale_fill_gradient(low="lightgoldenrodyellow", high="red", name = "Individual ASV \n Abundance") +
  labs(title = "ASV Abundance per Sample") +
  xlab("Sample") +
  ylab("ASV Abundance")
ASV.per.sample
ggsave("outputs/phyloseq_figures/ASV_abundance_per_sample_unsorted.png", width = 9, height = 6)

# Plot ASV abundance per sample (sorted)

ASV.per.sample.sorted <- plot_bar(ps.clean, fill = "Abundance") +
  aes(reorder(Sample, -Abundance)) +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  scale_fill_gradient(low="lightgoldenrodyellow", high="red", name = "Individual ASV \n Abundance") +
  labs(title = "ASV Abundance per Sample") +
  xlab("Sample") +
  ylab("ASV Abundance")
ASV.per.sample.sorted
ggsave("outputs/phyloseq_figures/ASV_abundance_per_sample_sorted.png", width = 9, height = 6)

# Plot ASV abundance per sample (sorted) (fill = phylum)

ASV.per.sample.sorted.phy <- plot_bar(ps.clean, fill = "Phylum") +
  aes(reorder(Sample, -Abundance)) +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  labs(title = "ASV Abundance per Sample") +
  scale_fill_discrete(labels = c("Ascomycota", "Basidiomycota","Chytridiomycota","Fungi Incertae sedis")) +
  xlab("Sample") +
  ylab("ASV Abundance")
ASV.per.sample.sorted.phy
ggsave("outputs/phyloseq_figures/ASV_abundance_per_sample_sorted_phylum_fill.png", width = 9, height = 6)

# Plot ASV abundance per sample (unsorted) (fill = phylum)

ASV.per.sample.unsorted.phy <- plot_bar(ps.clean, fill = "Phylum") +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  labs(title = "ASV Abundance per Sample") +
  scale_fill_discrete(labels = c("Ascomycota", "Basidiomycota","Chytridiomycota", "Fungi Incertae sedis")) +
  xlab("Sample") +
  ylab("ASV Abundance")
ASV.per.sample.unsorted.phy
ggsave("outputs/phyloseq_figures/ASV_abundance_per_sample_unsorted_phylum_fill.png", width = 9, height = 6)

# Plot ASV abundance by Phylum

ASV.per.phy <- plot_bar(ps.clean, x = "Phylum", fill = "Phylum") +
aes(reorder(Phylum, -Abundance)) +
geom_bar(aes(color= Phylum , fill= Phylum), stat="identity", position="stack") + theme(legend.position = "none") +
scale_x_discrete(labels = c("Ascomycota", "Basidiomycota", "Fungi Incertae sedis","Chytridiomycota")) +
scale_fill_discrete(labels = c("Ascomycota", "Basidiomycota","Chytridiomycota", "Fungi Incertae sedis")) +
labs(title = "Total ASV Abundance by Phylum") +
theme(axis.text.x = element_text(angle = 0, hjust= .5)) +
xlab("Phylum") +
ylab("Total ASV Abundance")
ASV.per.phy
ggsave("outputs/phyloseq_figures/ASV_abundance_total_phylum_fill.png", width = 9, height = 6)

# Plot ASV abundance by Class
x_labs <- ggplot_build(ASV.per.cls)$layout$panel_params[[1]]$x$get_labels()
x_labs.2 <- sub("c__", "", x_labs)
x_labs.2 <- sub("c_", "", x_labs.2)
x_labs.2 <- sub("_cls_Incertae_sedis", " Incertae sedis", x_labs.2)
print(x_labs)
print(x_labs.2)

ASV.per.cls <- plot_bar(ps.clean, x = "Class", fill = "Class") +
  aes(reorder(Class, -Abundance, sum)) +
  geom_bar(aes(color= Class, fill= Class), stat="identity", position="stack") + theme(legend.position = "none") +
  scale_x_discrete(labels = x_labs.2) +
  labs(title = "Total ASV Abundance by Class") +
  xlab("Class") +
  ylab("Total ASV Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ASV.per.cls
ggsave("outputs/phyloseq_figures/ASV_abundance_total_class_fill.png", width = 9, height = 6)

# Code below might be used later
# # Plot ASV abundance by Genus
# x_labs.3 <- ggplot_build(ASV.per.gen)$layout$panel_params[[1]]$x$get_labels()
# x_labs.4 <- sub("g__", "", x_labs.3)
# x_labs.4 <- sub("g_", "", x_labs.4)
# x_labs.4 <- sub("_gen_Incertae_sedis", " Incertae sedis", x_labs.4)
# print(x_labs.3)
# print(x_labs.4)
# 
# ASV.per.gen <- plot_bar(ps.clean, x = "Genus", fill = "Genus") +
#   aes(reorder(Genus, -Abundance, sum)) +
#   geom_bar(aes(color= Genus, fill= Genus), stat="identity", position="stack") + theme(legend.position = "none") +
#   scale_x_discrete(labels = x_labs.4) +
#   labs(title = "Total ASV Abundance by Genus") +
#   xlab("Genus") +
#   ylab("Total ASV Abundance") +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# ASV.per.gen
# ggsave("outputs/phyloseq_figures/ASV_abundance_total_class_fill.png", width = 9, height = 6)
# 
# # Subset out top 20 species
# top20 <- names(sort(taxa_sums(ps.clean), decreasing=TRUE))[1:20]
# ps.top20 <- transform_sample_counts(ps.clean, function(OTU) OTU/sum(OTU))
# ps.top20 <- prune_taxa(top20, ps.top20)
# 
# ASV.per.sp <- plot_bar(ps.top20, x = "Species", fill = "Disease") +
#   geom_bar(aes(color= Disease, fill= Disease), stat="identity", position="stack") +
#   aes(reorder(Species, -Abundance, sum))
# ASV.per.sp
# ggsave("outputs/phyloseq_figures/ASV_Abundance_per_Sample_sorted.png", width = 9, height = 6)


# Plot ordination NMDS
ps.ord.nmds <- ordinate(ps.clean, "NMDS", "bray")
ASV.ordination.nmds.phy <- plot_ordination(ps.clean, ps.ord.nmds, type="taxa", color="Phylum", title="ASV Ordination (NMDS)") +
  scale_color_hue(labels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Fungi Incertae sedis"))
ASV.ordination.nmds.phy
ggsave("outputs/phyloseq_figures/ASV_ord_NMDS_phylum_fill.png", width = 8, height = 6)

# Plot ordination NMDS (by Disease)
ps.ord.nmds <- ordinate(ps.clean, "NMDS", "bray")
SAM.ordination.nmds.dis <- plot_ordination(ps.clean, ps.ord.nmds, type="sample", color="Disease", title="Sample Ordination (NMDS)") # +scale_color_hue(labels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Fungi Incertae sedis"))
SAM.ordination.nmds.dis
ggsave("outputs/phyloseq_figures/SAM_ord_NMDS_disease_fill.png", width = 8, height = 6)

# Plot PCoA by phylum
ps.ord.pcoa <- ordinate(ps.clean, "PCoA", "bray")
ASV.ordination.pcoa.tax <- plot_ordination(ps.clean, ps.ord.pcoa, type="taxa", color="Phylum", title="ASV Ordination (PCoA)") +
  scale_color_hue(labels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Fungi Incertae sedis"))
ASV.ordination.pcoa.tax
ggsave("outputs/phyloseq_figures/ASV_ord_PCoA_phylum_fill.png", width = 8, height = 6)

# Plot PCoA by disease
ps.ord.pcoa <- ordinate(ps.clean, "PCoA", "bray")
ASV.ordination.pcoa.dis <- plot_ordination(ps.clean, ps.ord.pcoa, type="sample", color="Disease", title="ASV Ordination (PCoA)")
ASV.ordination.pcoa.dis
ggsave("outputs/phyloseq_figures/ASV_ord_PCoA_disease_fill.png", width = 8, height = 6)

# Plot PCA ASV + Phylum
ps.ord.pca <- ordinate(ps.clean, "RDA", "bray")
ASV.ordination.pca <- plot_ordination(ps.clean, ps.ord.pca, type="taxa", color="Phylum", title="ASV Ordination (PCA)") +
scale_color_hue(labels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Fungi Incertae sedis"))
ASV.ordination.pca
ggsave("outputs/phyloseq_figures/ASV_ord_PCA_phylum_fill.png", width = 8, height = 6)

Total_rep_per_ASV <- as.matrix(colSums(clean.seq.tab))
hist(Total_rep_per_ASV, 100)

sample_data(ps.clean)['sample_id'] <- row.names(sample_data(ps.clean))
# Plot PCA sample + Disease
ps.ord.pca <- ordinate(ps.clean, "RDA", "bray")
SAM.ordination.pca <- plot_ordination(ps.clean, ps.ord.pca, type="samples", color="Disease", title="Sample Ordination (PCA)", label = "sample_id") # + theme(legend.position = "none")
SAM.ordination.pca
ggsave("outputs/phyloseq_figures/SAMPLE_ord_PCA_disease_fill.png", width = 8, height = 6)

#### Investigating NMDS ordination to see if we actually have low variability
#### Turns out it was just an outlier, a sample with virtually no ASVs

### using phyloseq's distance() -> vegan metaMDS
dist_1 <- as.matrix(distance(ps.clean, method = "bray", type = "samples"))
NMDS1 <- metaMDS(dist_1, k = 2, trymax = 100, trace = F)
stressplot(NMDS1)
jpeg(filename="outputs/phyloseq_figures/Stress_plot_phyloseq.png", width = 7, height = 4, units = 'in', res = 300)
stressplot(NMDS1)
dev.off()

dist_1[1:5,1:5]
plot(NMDS1, type = "t")
hist(dist_1)

### using vegan
dist_2 <-as.matrix(vegdist(clean.seq.tab, method = "bray"))
NMDS2 <- metaMDS(dist_2, k = 2, trymax = 100, trace = F)
stressplot(NMDS2)
jpeg(filename="outputs/phyloseq_figures/Stress_plot_vegan.png", width = 7, height = 4, units = 'in', res = 300)
stressplot(NMDS2)
dev.off()

dist_2[1:5,1:5]
plot(NMDS2, type = "t")
hist(dist_2)

### is number of dimensions an issue? Even with one dimension, stress = zero
### after removal of outlier, stress plot behaving as expected although a little high
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}
#NMDS.scree(dist_1)
#NMDS.scree(dist_2)

# RDA still clumped?
rda1 <- rda(otu_table(ps.clean))
rda2 <- rda(clean.seq.tab)
rda
plot(rda1, display = "sites")
plot(rda2)
text(rda1, labels=rownames(otu_table(ps.clean)))
text(rda1, labels=rownames(tax_table(ps.clean)), cex = .1)

# Relative abundance by disease - not amazing
ps.t <- transform_sample_counts(ps.clean, function(x) x / sum(x) * 100)
ps.t.p <- aggregate_taxa(ps.t, "Phylum")
relative <- plot_composition(ps.t.p, group_by = "Disease", average_by = "Disease") +
scale_fill_hue(labels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Fungi Incertae sedis"), name = "Phylum") +
labs(title = "Relative Phylum Abundance by Averaged by Disease Class") +
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
ylab("Relative Phylum Abundance")
relative
ggsave("outputs/phyloseq_figures/relative_abundance_disease_facet.png", width = 8, height = 6)

Sph <- fungal.taxa.clean[fungal.taxa.clean[,6]=="g__Sphaerulina",]
print(Sph)
which(fungal.taxa[,6]=="g__Sphaerulina")
fungal.taxa[340,6]
# microViz is cool

#ps.clean%>%
    
#phyloseq::merge_samples(group = "Disease") %>%

test <-comp_barplot(
    ps.clean,
    tax_level = "Family",
    order_with_all_taxa = TRUE,
    taxon_renamer = function(x) stringr::str_replace_all(x, "f__", ""),
    tax_order = prev,
    facet_by = "Disease",
    label = "Genotype",
    n_taxa = 25,
    merge_other = FALSE,
    bar_outline_colour = NA) +
    guides(fill = guide_legend(ncol = 1)) +
    ggtitle("Relative Endophyte Family Abundance by Host Genotype and Disease Class") +
    xlab("Host Genotype") +
    ylab("Relative Abundance (Family)") +
    theme(axis.text.y=element_text(size=6)) +
    coord_flip()
test
ggsave("outputs/phyloseq_figures/relative_abundance_viz.png", width = 12, height = 9)


# Aplha Diversity Boxplot
plot_richness(ps.clean,
              x="Disease",
              measures=c("Shannon", "Simpson", "InvSimpson"),
              color="Disease") +
  geom_boxplot(alpha=0.6) +
  xlab("Disease Class") +
  theme(axis.text.x = element_text(angle = 0, hjust= .5)) +
  ggtitle("Alpha Diversity by Disease Class") +
  theme(legend.position = "none") 
ggsave("outputs/phyloseq_figures/alpha_diversity_boxplot.png", width = 8, height = 6)
