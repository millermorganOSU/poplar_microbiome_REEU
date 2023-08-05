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
library(scales); packageVersion("scales")
library(patchwork); packageVersion("patchwork")
library(gridExtra); packageVersion("gridExtra")
library(ggpubr); packageVersion("ggpubr")
library(forcats); packageVersion("forcats")
library(RColorBrewer); packageVersion("RColorBrewer")

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
unique(fungal.taxa[,2])
fungal.taxa <- as.matrix(fungal.taxa[fungal.taxa[,2]!="p__Chytridiomycota",])
# 466 fungal ASVs ## remove Chytrids that were likely incorrectly classified by UNITE
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

# 520 samples with 466 FUNGAL ASVs
dim(seqtab.nochim.fungal)

# List samples with metadata
sample.names.wmeta <- as.list(rownames(wood.meta))

# Subset rows (ONLY SAMPLES WITH METADATA)
seqtab.nochim.fungal.wmeta <- seqtab.nochim.fungal[unlist(sample.names.wmeta),]

# 236 samples with 466 fungal ASVs as columns
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
# 234 samples by 356 ASVs

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

# 466 ASVs with taxa (and index and accession)
dim(fungal.taxa)

# Keep only ASVs that have abundance in seq table
non.zero.ASVs <- as.list(colnames(clean.seq.tab))

fungal.taxa.clean <- fungal.taxa[unlist(non.zero.ASVs),]
fungal.taxa.clean <- fungal.taxa.clean[,1:7]
dim(fungal.taxa.clean)
unique(fungal.taxa.clean[,2])

# Replace Mycosphaerellaceae with musiva
view(fungal.taxa.clean)
fungal.taxa.clean[fungal.taxa.clean == "g__Mycosphaerellaceae_gen_Incertae_sedis"] <- "g__Sphaerulina"
fungal.taxa.clean[fungal.taxa.clean == "s__Mycosphaerellaceae_sp_Incertae_sedis"] <- "s__Sphaerulina_musiva"
fungal.taxa.clean[fungal.taxa.clean == "s__Sphaerulina_musiva"] <- "s__musiva"

# High rep
Total_rep_per_ASV <- as.matrix(sort(colSums(clean.seq.tab), decreasing = TRUE))
hist(Total_rep_per_ASV, 50)
high_rep <- head(Total_rep_per_ASV, 20)

rc <- matrix(, nrow = 0, ncol = 9)
for ( i in rownames(high_rep)) {
  ASV <- as.matrix(i)
  lin <- t(as.matrix(fungal.taxa.clean[i,]))
  add <- as.matrix(high_rep[i,])
  add <- cbind(add, ASV)
  add <- cbind(add, lin)
  rc <- rbind(rc, add)
}

colnames(rc) <- c("Total Abundance","ASV","Kingdom","Phylum","Class","Order","Family","Genus","Species")

# Create phyloseq object clean
ps.clean <- phyloseq(otu_table(clean.seq.tab, taxa_are_rows=FALSE), sample_data(wood.meta.clean), tax_table(fungal.taxa.clean))

# Plot ASV abundance per sample

options(scipen = 999)

ASV.per.sample <- plot_bar(ps.clean, fill = "Abundance") +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  scale_fill_gradient(low="lightgoldenrodyellow", high="red", name = "Individual ASV \n Abundance") +
  labs(title = "ASV Abundance per Sample") +
  xlab("Sample") +
  ylab("ASV Abundance") +
  scale_y_continuous(labels = scales::label_number_si())
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
  ylab("ASV Abundance") +
  scale_y_continuous(labels = scales::label_number_si())
ASV.per.sample.sorted
ggsave("outputs/phyloseq_figures/ASV_abundance_per_sample_sorted.png", width = 13.333, height = 7.5)

# Plot ASV abundance per sample (sorted) (fill = phylum)

phy.colors <- c("#E41A1C", "#377EB8", "green3")
phy.colors <- c("#0072B2", "#F0E442", "maroon1")

"#999999", "#E69F00", "#56B4E9", "#009E73", 
"#F0E442", "#0072B2", "#D55E00", "#CC79A7"
"#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
"#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

ASV.per.sample.sorted.phy <- plot_bar(ps.clean, fill = "Phylum") +
  aes(reorder(Sample, -Abundance)) +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  labs(title = "ASV Abundance per Sample") +
  scale_fill_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors) +
  xlab("Sample") +
  ylab("ASV Abundance") +
  scale_y_continuous(labels = scales::label_number_si())
ASV.per.sample.sorted.phy
ggsave("outputs/phyloseq_figures/ASV_abundance_per_sample_sorted_phylum_fill.png", width = 13.333, height = 7.5)

# Plot ASV abundance per sample (unsorted) (fill = phylum)

ASV.per.sample.unsorted.phy <- plot_bar(ps.clean, fill = "Phylum") +
  scale_x_discrete(labels=seq(1,236)) +
  theme(axis.text.x = element_blank()) +
  labs(title = "ASV Abundance per Sample") +
  scale_fill_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors) +
  xlab("Sample") +
  ylab("ASV Abundance") +
  scale_y_continuous(labels = scales::label_number_si())
ASV.per.sample.unsorted.phy
ggsave("outputs/phyloseq_figures/ASV_abundance_per_sample_unsorted_phylum_fill.png", width = 9, height = 6)

# Plot ASV abundance by Phylum

my_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_grid = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

ASV.per.phy <- my_plot_bar(ps.clean, x = "Phylum", fill = "Phylum")+
scale_fill_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors) +
aes(reorder(Phylum, -Abundance)) +
#geom_bar(aes(color = Phylum , fill= Phylum), stat="identity", position="stack") + 
theme(legend.position = "none") +
scale_x_discrete(labels = c("Ascomycota", "Basidiomycota", "Fungi Incertae sedis")) +
labs(title = "Total ASV Abundance by Phylum") +
theme(axis.text.x = element_text(angle = 0, hjust= .5)) +
xlab("\nPhylum") +
ylab("Total ASV Abundance\n") +
scale_y_continuous(labels = scales::label_number_si())
ASV.per.phy
ggsave("outputs/phyloseq_figures/ASV_abundance_total_phylum_fill.png", width = 13.333, height = 7.5)

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
  xlab("\nClass") +
  ylab("Total ASV Abundance\n") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(labels = scales::label_number_si())
ASV.per.cls
ggsave("outputs/phyloseq_figures/ASV_abundance_total_class_fill.png", width = 13.333, height = 7.5)

# Make ords 
ps.ord.nmds.b <- ordinate(ps.clean, "NMDS", "bray")
ps.ord.nmds.j <- ordinate(ps.clean, "NMDS", "jaccard")
ps.ord.pcoa.b <- ordinate(ps.clean, "PCoA", "bray")
ps.ord.pcoa.j <- ordinate(ps.clean, "PCoA", "jaccard")
ps.ord.pca.b <- ordinate(ps.clean, "RDA", "bray")
ps.ord.pca.j <- ordinate(ps.clean, "RDA", "jaccard")

###

# Plot ordination NMDS
ASV.ordination.nmds.phy.b <- plot_ordination(ps.clean, ps.ord.nmds.b, type="taxa", color="Phylum", title="ASV Ordination (NMDS, Bray-Curtis)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  scale_color_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors)
  ASV.ordination.nmds.phy.b
ggsave("outputs/phyloseq_figures/ORD_B_NMDS_ASV.png", width = 13.333, height = 7.5)

# Plot ordination NMDS (by Disease)
SAM.ordination.nmds.dis.b <- plot_ordination(ps.clean, ps.ord.nmds.b, type="sample", color="Disease", title="Sample Ordination (NMDS, Bray-Curtis)") + 
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black") +
SAM.ordination.nmds.dis.b
ggsave("outputs/phyloseq_figures/ORD_B_NMDS_SAM.png", width = 13.333, height = 7.5)

# Plot ordination NMDS
ASV.ordination.nmds.phy.j <- plot_ordination(ps.clean, ps.ord.nmds.j, type="taxa", color="Phylum", title="ASV Ordination (NMDS, Jaccard)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  scale_color_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors)
ASV.ordination.nmds.phy.j
ggsave("outputs/phyloseq_figures/ORD_J_NMDS_ASV.png", width = 13.333, height = 7.5)

# Plot ordination NMDS (by Disease)
SAM.ordination.nmds.dis.j <- plot_ordination(ps.clean, ps.ord.nmds.j, type="sample", color="Disease", title="Sample Ordination (NMDS, Jaccard)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black")
SAM.ordination.nmds.dis.j
ggsave("outputs/phyloseq_figures/ORD_J_NMDS_SAM.png", width = 13.333, height = 7.5)

patchNMDS <- patchwork::wrap_plots(ASV.ordination.nmds.phy.b, 
                                SAM.ordination.nmds.dis.b, 
                                ASV.ordination.nmds.phy.j,
                                SAM.ordination.nmds.dis.j,
                                nrow = 2, guides = "collect")
patchNMDS

ggsave("outputs/phyloseq_figures/ORD_NMDS_PATCH.png", width = 13.333, height = 7.5)

###

# Plot PCoA by phylum
ASV.ordination.pcoa.tax.b <- plot_ordination(ps.clean, ps.ord.pcoa.b, type="taxa", color="Phylum", title="ASV Ordination (PCoA, Bray-Curtis)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  scale_color_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors)
ASV.ordination.pcoa.tax.b 
ggsave("outputs/phyloseq_figures/ORD_B_PCoA_ASV.png", width = 13.333, height = 7.5)

# Plot PCoA by disease
ASV.ordination.pcoa.dis.b <- plot_ordination(ps.clean, ps.ord.pcoa.b, type="sample", color="Disease", title="Sample Ordination (PCoA, Bray-Curtis)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black")
ASV.ordination.pcoa.dis.b
ggsave("outputs/phyloseq_figures/ORD_B_PCoA_SAM.png", width = 13.333, height = 7.5)

# Plot PCoA by phylum
ASV.ordination.pcoa.tax.j <- plot_ordination(ps.clean, ps.ord.pcoa.j, type="taxa", color="Phylum", title="ASV Ordination (PCoA, Jaccard)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  scale_color_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors)
ASV.ordination.pcoa.tax.j 
ggsave("outputs/phyloseq_figures/ORD_J_PCoA_ASV.png", width = 13.333, height = 7.5)

# Plot PCoA by disease
ASV.ordination.pcoa.dis.j <- plot_ordination(ps.clean, ps.ord.pcoa.j, type="sample", color="Disease", title="Sample Ordination (PCoA, Jaccard)") +
geom_point(size=2) +
geom_point(shape = 1,size = 2,colour = "black")
ASV.ordination.pcoa.dis.j
ggsave("outputs/phyloseq_figures/ORD_J_PCoA_SAM.png", width = 13.333, height = 7.5)

###

patchPCoA <- patchwork::wrap_plots(ASV.ordination.pcoa.tax.b, 
                                   ASV.ordination.pcoa.dis.b, 
                                   ASV.ordination.pcoa.tax.j,
                                   ASV.ordination.pcoa.dis.j,
                                   nrow = 2, guides = "collect")
patchPCoA

ggsave("outputs/phyloseq_figures/ORD_PCoA_PATCH.png", width = 13.333, height = 7.5)

#sample_data(ps.clean)['sample_id'] <- row.names(sample_data(ps.clean))

# Plot PCA ASV + Phylum
ASV.ordination.pca.b <- plot_ordination(ps.clean, ps.ord.pca.b, type="taxa", color="Phylum", title="ASV Ordination (PCA, Bray-Curtis)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  scale_color_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors)
ASV.ordination.pca.b
ggsave("outputs/phyloseq_figures/ORD_B_PCA_ASV.png", width = 13.333, height = 7.5)

# Plot PCA sample + Disease
SAM.ordination.pca.b <- plot_ordination(ps.clean, ps.ord.pca.b, type="samples", color="Disease", title="Sample Ordination (PCA, Bray-Curtis)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black")
SAM.ordination.pca.b
ggsave("outputs/phyloseq_figures/ORD_B_PCA_SAM.png", width = 13.333, height = 7.5)

# Plot PCA ASV + Phylum
ASV.ordination.pca.j <- plot_ordination(ps.clean, ps.ord.pca.j, type="taxa", color="Phylum", title="ASV Ordination (PCA, Jaccard)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  scale_color_manual(labels = c("Ascomycota", "Basidiomycota","Fungi Incertae sedis"), values = phy.colors)
ASV.ordination.pca.j
ggsave("outputs/phyloseq_figures/ORD_J_PCA_ASV.png", width = 13.333, height = 7.5)

# Plot PCA sample + Disease
SAM.ordination.pca.j <- plot_ordination(ps.clean, ps.ord.pca.j, type="samples", color="Disease", title="Sample Ordination (PCA, Jaccard)") +
  geom_point(size=2) +
  geom_point(shape = 1,size = 2,colour = "black")
SAM.ordination.pca.j
ggsave("outputs/phyloseq_figures/ORD_J_PCA_SAM.png", width = 13.333, height = 7.5)

patchPCA <- patchwork::wrap_plots(ASV.ordination.pca.b, 
                                  SAM.ordination.pca.b, 
                                  ASV.ordination.pca.j,
                                  SAM.ordination.pca.j,
                                   nrow = 2, guides = "collect")
patchPCA

ggsave("outputs/phyloseq_figures/ORD_PCA_PATCH.png", width = 13.333, height = 7.5)


# Plot ordination NMDS with rings
SAM.ordination.nmds.dis.b.e <- plot_ordination(ps.clean, ps.ord.nmds.b, type="sample", color="Disease", title="Sample Ordination (NMDS, Bray-Curtis)") + 
  geom_point(size=4) +
  geom_point(shape = 1,size = 4,colour = "black") +
  stat_ellipse(type = "t", linewidth = 2) +
  theme_bw() +
  labs(caption = "Ellipses created using multivariate 't' distribution\nFinal Stress: 0.205") +
  theme(plot.caption.position = "plot",
        plot.caption = element_text(hjust = 1))
  SAM.ordination.nmds.dis.b.e
ggsave("outputs/phyloseq_figures/ORD_B_NMDS_SAM_ELLIPSES.png", width = 13.333, height = 7.5)

# Relative abundance by disease - not amazing
ps.t <- transform_sample_counts(ps.clean, function(x) x / sum(x) * 100)
ps.t.p <- aggregate_taxa(ps.t, "Phylum")
relative.p <- plot_composition(ps.t.p, group_by = "Disease", average_by = "Disease") +
scale_fill_hue(labels = c("Ascomycota", "Basidiomycota", "Fungi Incertae sedis"), name = "Phylum") +
labs(title = "Relative Phylum Abundance Averaged by Disease Class") +
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
ylab("Relative Phylum Abundance")
relative.p
ggsave("outputs/phyloseq_figures/relative_abundance_disease_facet.png", width = 8, height = 6)

# Relative abundance by disease - not amazing
ps.t <- transform_sample_counts(ps.clean, function(x) x / sum(x) * 100)
ps.t.c <- aggregate_taxa(ps.t, "Class")
relative.c <- plot_composition(ps.t.c, group_by = "Disease", average_by = "Disease") +
  scale_fill_hue(name = "Class") +
  labs(title = "Relative Class Abundance Averaged by Disease Class") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Relative Class Abundance")
relative.c
ggsave("outputs/phyloseq_figures/relative_abundance_disease_facet.png", width = 8, height = 6)

unique(fungal.taxa.clean[,5])
Sph <- fungal.taxa.clean[fungal.taxa.clean[,6]=="g__Sphaerulina",]
print(Sph)
ch <- as.matrix(fungal.taxa.clean[fungal.taxa.clean[,2]=="p__Chytridiomycota",])
fungal.taxa[340,6]

# microViz is cool

ps.clean %>%
  phyloseq::merge_samples(group = "Disease") %>%
  comp_barplot(tax_level = "Phylum", n_taxa = 3, bar_width = 0.8,
               taxon_renamer = function(x) stringr::str_replace_all(x, c("p__" = "", "_phy_Incertae_sedis" = " Incertae sedis"))) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggtitle("Average Endophyte Phylum Abundance by Disease Class") +
  xlab("\nDisease Class") +
  ylab("Relative Abundance\n") 
ggsave("outputs/phyloseq_figures/relative_abundance_phy_average.png", width = 13.333, height = 7.5)

ps.clean %>%
  phyloseq::merge_samples(group = "Disease") %>%
  comp_barplot(tax_level = "Class", n_taxa = 10, bar_width = 0.8,
               taxon_renamer = function(x) stringr::str_replace_all(x, c("c__" = "", "_cls_Incertae_sedis" = " Incertae sedis"))) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggtitle("Average Endophyte Class Abundance by Disease Class") +
  xlab("\nDisease Class") +
  ylab("Relative Abundance\n") 
ggsave("outputs/phyloseq_figures/relative_abundance_class_average.png", width = 13.333, height = 7.5)

ps.clean %>%
  phyloseq::merge_samples(group = "Disease") %>%
  comp_barplot(tax_level = "Order", n_taxa = 10, bar_width = 0.8,
               taxon_renamer = function(x) stringr::str_replace_all(x, c("o__" = "", "_ord_Incertae_sedis" = " Incertae sedis"))) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggtitle("Average Endophyte Order Abundance by Disease Class") +
  xlab("\nDisease Class") +
  ylab("Relative Abundance\n") 
ggsave("outputs/phyloseq_figures/relative_abundance_ord_average.png", width = 13.333, height = 7.5)

ps.clean %>%
  phyloseq::merge_samples(group = "Disease") %>%
  comp_barplot(tax_level = "Family", n_taxa = 10, bar_width = 0.8,
               taxon_renamer = function(x) stringr::str_replace_all(x, c("f__" = "", "_fam_Incertae_sedis" = " Incertae sedis"))) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggtitle("Average Endophyte Family Abundance by Disease Class") +
  xlab("\nDisease Class") +
  ylab("Relative Abundance\n") 
ggsave("outputs/phyloseq_figures/relative_abundance_fam_average.png", width = 13.333, height = 7.5)

ps.clean %>%
  phyloseq::merge_samples(group = "Disease") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 10, bar_width = 0.8,
               taxon_renamer = function(x) stringr::str_replace_all(x, c("g__" = "", "_gen_Incertae_sedis" = " Incertae sedis"))) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggtitle("Average Endophyte Genus Abundance by Disease Class") +
  xlab("\nDisease Class") +
  ylab("Relative Abundance\n") 
ggsave("outputs/phyloseq_figures/relative_abundance_gen_average.png", width = 13.333, height = 7.5)

ps.clean %>% 
  tax_fix(unknowns = c("s__penicillioides")) %>%
  phyloseq::merge_samples(group = "Disease") %>%
  comp_barplot(tax_level = "Species", n_taxa = 10, bar_width = 0.8,
               taxon_renamer = function(x) stringr::str_replace_all(x, c("s__" = "", "_sp_Incertae_sedis" = " Incertae sedis"))) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggtitle("Average Endophyte Species Abundance by Disease Class") +
  xlab("\nDisease Class") +
  ylab("Relative Abundance\n") 
ggsave("outputs/phyloseq_figures/relative_abundance_sp_average.png", width = 13.333, height = 7.5)

C.f <- ps.clean %>%
  ps_filter(Disease == "Canker") %>%
  comp_barplot(tax_level = "Family",
               taxon_renamer = function(x) stringr::str_replace_all(x, c("f__" = "", "_fam_Incertae_sedis" = " Incertae sedis")),
               label = "Genotype",
               n_taxa = 10,
               bar_outline_colour = NA) +
  theme(axis.text.y=element_text(size=7),
        #axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        #legend.position = "None",
        plot.title = element_text(hjust = .5)) +
  xlab("Host Genotype\n") +
  ggtitle("Canker")
C.f
H.f <- ps.clean %>%
  ps_filter(Disease == "Healthy") %>%
  comp_barplot(tax_level = "Family",
               taxon_renamer = function(x) stringr::str_replace_all(x, c("f__" = "", "_fam_Incertae_sedis" = " Incertae sedis")),
               label = "Genotype",
               n_taxa = 10,
               bar_outline_colour = NA) +
  theme(axis.text.y=element_text(size=4),
        axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        legend.position = "None",
        plot.title = element_text(hjust = .5)) +
  ylab("\nRelative Abundance") +
  ggtitle("Healthy")
H.f
nC.f <- ps.clean %>%
  ps_filter(Disease == "Non-canker") %>%
  comp_barplot(tax_level = "Family",
               taxon_renamer = function(x) stringr::str_replace_all(x, c("f__" = "", "_fam_Incertae_sedis" = " Incertae sedis")),
               label = "Genotype",
               n_taxa = 10,
               bar_outline_colour = NA) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "None",
        plot.title = element_text(hjust = .5)) +
  
  ggtitle("Non-Canker")
nC.f

patch <- patchwork::wrap_plots(C.f, H.f, nC.f, nrow = 1, guides = "collect") & coord_flip() & plot_annotation(
    title = "Relative Endophyte Family Abundance by Host Genotype and Disease Class",
    theme = theme(plot.title = element_text(size = 14, face = "bold",)))
patch

ggsave("outputs/phyloseq_figures/relative_abundance_fam_VIZ.png", width = 13.333, height = 7.5)
  
C.s <- ps.clean %>%
  tax_fix(unknowns = c("s__penicillioides")) %>%
  ps_filter(Disease == "Canker") %>%
  comp_barplot(tax_level = "Species",
               taxon_renamer = function(x) stringr::str_replace_all(x, c("s__" = "", "_sp_Incertae_sedis" = " Incertae sedis")),
               label = "Genotype",
               n_taxa = 10,
               bar_outline_colour = NA) +
  theme(axis.text.y=element_text(size=7),
        #axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        #legend.position = "None",
        plot.title = element_text(hjust = .5)) +
  xlab("Host Genotype\n") +
  ggtitle("Canker")
C.s
H.s <- ps.clean %>%
  tax_fix(unknowns = c("s__penicillioides")) %>%
  ps_filter(Disease == "Healthy") %>%
  comp_barplot(tax_level = "Species",
               taxon_renamer = function(x) stringr::str_replace_all(x, c("s__" = "", "_sp_Incertae_sedis" = " Incertae sedis")),
               label = "Genotype",
               n_taxa = 10,
               bar_outline_colour = NA) +
  theme(axis.text.y=element_text(size=4),
        axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        legend.position = "None",
        plot.title = element_text(hjust = .5)) +
  ylab("\nRelative Abundance") +
  ggtitle("Healthy")
H.s
nC.s <- ps.clean %>%
  tax_fix(unknowns = c("s__penicillioides")) %>%
  ps_filter(Disease == "Non-canker") %>%
  comp_barplot(tax_level = "Species",
               taxon_renamer = function(x) stringr::str_replace_all(x, c("s__" = "", "_sp_Incertae_sedis" = " Incertae sedis")),
               label = "Genotype",
               n_taxa = 10,
               bar_outline_colour = NA) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "None",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Non-Canker")
nC.s

patch2 <- patchwork::wrap_plots(C.s, H.s, nC.s, nrow = 1, guides = "collect") & coord_flip() & plot_annotation(
  title = "Relative Endophyte Species Abundance by Host Genotype and Disease Class",
  theme = theme(plot.title = element_text(size = 14, face = "bold",)))
patch2

ggsave("outputs/phyloseq_figures/relative_abundance_sp_VIZ.png", width = 13.333, height = 7.5)

# Alpha Diversity Boxplot

my.comparison <- list( c("Canker", "Healthy"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))


Alpha.plot <- plot_richness(ps.clean,
              x="Disease",
              measures=c("Shannon", "Simpson", "InvSimpson", "Observed"),
              color="Disease") +
  geom_boxplot(alpha=0.6) +
  xlab("\nDisease Class") +
  ylab("Alpha Diversity Measure\n") +
  theme(axis.text.x = element_text(angle = 0, hjust= .5)) +
  ggtitle("Alpha Diversity by Disease Class") +
  theme(legend.position = "none") +
  labs(caption = "Unpaired Wilcoxon Test (p-values <= 0.001)") +
  stat_compare_means(method = "wilcox.test", comparisons = my.comparison, label = "p.signif", symnum.args = symnum.args)
ggsave("outputs/phyloseq_figures/alpha_diversity_boxplot_observed.png", units = "in", width = 13.333, height = 7.5)


# Beta diversity

# Permanova test using the vegan package
metadata <- as(sample_data(ps.clean), "data.frame")

bray <- distance(ps.clean, method="bray")
jaccard <- distance(ps.clean, method="jaccard")

pre.bray <- adonis2(bray ~ Disease, data = metadata)
pre.jaccard <- adonis2(jaccard ~ Disease, data = metadata)

pre.bray
pre.jaccard

beta.bray <- betadisper(bray, metadata$Disease)
beta.jaccard <- betadisper(jaccard, metadata$Disease)
permutest(beta.bray)
permutest(beta.jaccard)

plot(beta.bray, identify.or)
plot(beta.jaccard)
boxplot(beta.bray)
boxplot(beta.jaccard)
eigenvals(beta.bray)
eigenvals(beta.jaccard)
TukeyHSD(beta.bray, which = "group", ordered = FALSE,
         conf.level = 0.95)
TukeyHSD(beta.jaccard, which = "group", ordered = FALSE,
         conf.level = 0.95)
# Misc.

unique(fungal.taxa.clean[,2])
chy <- fungal.taxa.clean[fungal.taxa.clean[,2]=="p__Chytridiomycota",]

unique(fungal.taxa.clean[,5])
myco <- fungal.taxa.clean[fungal.taxa.clean[,5]=="f__Mycosphaerellaceae",]
sph <- fungal.taxa.clean[fungal.taxa.clean[,6]=="g__Sphaerulina",]

raw.stats <- read.table("outputs/raw_stats.tsv", header = TRUE)
mean(raw.stats[,4])
sum(raw.stats[,4])

n.stats <- read.table("outputs/n_stats.tsv", header = TRUE)
mean(n.stats[,4])
sum(n.stats[,4])

cut.stats <- read.table("outputs/cut_stats.tsv", header = TRUE)
mean(cut.stats[,4])
sum(cut.stats[,4])

filt.stats <- read.table("outputs/filtered_stats.tsv", header = TRUE)
mean(filt.stats[,4])
sum(filt.stats[,4])
nrow((filt.stats))

dim(seqtab.nochim)

dis.class.dis <- as.matrix(table(wood.meta.clean$Disease))
gen.dis <- as.matrix(table(wood.meta.clean$Genotype))
which(gen.dis[,1]==2)

canker.library <- as.matrix(ps.clean %>%
                              ps_filter(Disease == "Canker") %>%
                              sample_sums())

mean(canker.library[,1])

healthy.library <- as.matrix(ps.clean %>%
                              ps_filter(Disease == "Healthy") %>%
                              sample_sums())

mean(healthy.library[,1])

non.library <- as.matrix(ps.clean %>%
                               ps_filter(Disease == "Non-canker") %>%
                               sample_sums())

mean(non.library[,1])

dist_ex <- dist_1[1:10,1:10]

write.csv(dist_ex, "outputs/distance_example.csv")

richness <- estimate_richness(ps.clean, measures="Shannon")
richness_ex <- head(richness, 10)

write.csv(richness_ex, "outputs/richness_example.csv")

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

# test <-comp_barplot(
#   ps.clean,
#   tax_level = "Family",
#   taxon_renamer = function(x) stringr::str_replace_all(x, c("f__" = "", "_fam_Incertae_sedis" = " Incertae sedis")),
#   facet_by = "Disease",
#   label = "Genotype",
#   n_taxa = 10,
#   bar_outline_colour = NA) +
#   guides(fill = guide_legend(ncol = 1)) +
#   ggtitle("Relative Endophyte Family Abundance by Host Genotype and Disease Class") +
#   xlab("Host Genotype") +
#   ylab("Relative Abundance") +
#   theme(axis.text.y=element_text(size=6)) +
#   coord_flip()
# test
# ggsave("outputs/phyloseq_figures/relative_abundance_viz_1.png", width = 13.333, height = 7.5)


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

#### Investigating NMDS ordination to see if we actually have low variability
#### Turns out it was just an outlier, a sample with virtually no ASVs

# ### using phyloseq's distance() -> vegan metaMDS
#dist_1 <- as.matrix(distance(ps.clean, method = "bray", type = "samples"))
#NMDS1 <- metaMDS(dist_1, k = 2, trymax = 100, trace = F)
#NMDS1
#stressplot(NMDS1, main = "NMDS Sample Ordination Stress Plot")
#mtext(stp, "Final Stress: 0.205")
#NMDS1$stress
# jpeg(filename="outputs/phyloseq_figures/Stress_plot_phyloseq.png", width = 7, height = 4, units = 'in', res = 300)
# stressplot(NMDS1)
# dev.off()
# 
# dist_1[1:5,1:5]
# plot(NMDS1, type = "t")
# hist(dist_1)
# 
# ### using vegan
#dist_2 <-as.matrix(vegdist(clean.seq.tab, method = "jaccard"))
#NMDS2 <- metaMDS(dist_2, k = 2, trymax = 100, trace = F)
#stressplot(NMDS2)
# jpeg(filename="outputs/phyloseq_figures/Stress_plot_vegan.png", width = 7, height = 4, units = 'in', res = 300)
# stressplot(NMDS2)
# dev.off()
# 
# dist_2[1:5,1:5]
# plot(NMDS2, type = "t")
# hist(dist_2)
# 
# ### is number of dimensions an issue? Even with one dimension, stress = zero
# ### after removal of outlier, stress plot behaving as expected although a little high
# NMDS.scree <- function(x) { #where x is the name of the data frame variable
#   plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
#   for (i in 1:10) {
#     points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
#   }
# }
# #NMDS.scree(dist_1)
# #NMDS.scree(dist_2)

# # RDA still clumped?
# rda1 <- rda(otu_table(ps.clean))
# rda2 <- rda(clean.seq.tab)
# rda
# plot(rda1, display = "sites")
# plot(rda2)
# text(rda2, labels=colSums(clean.seq.tab))
# text(rda1, labels=(tax_table(ps.clean))[,7], cex = 1)

