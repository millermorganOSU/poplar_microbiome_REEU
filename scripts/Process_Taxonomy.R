## Process Taxonomy

# Import packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(readxl); packageVersion("readxl")
library(reshape2); packageVersion("reshape2")

# Set WD for or local machine
setwd('~/LeBoldus/local_git/poplar_microbiome_REEU')

# Load assigned taxonomy
load("outputs/Q20/AssignedTaxa_Q20.RData")
load("outputs/Q28/AssignedTaxa_Q28.RData")

# Load sequence tables
load("outputs/Q20/SeqTabQ20.RData")
SeqTabQ20 <- seqtab.nochim
load("outputs/Q28/SeqTabQ28.RData")
SeqTabQ28 <- seqtab.nochim

# Load metadata
disease.scores <- read_xlsx("data/sample_diseasescores.xlsx")
wood.meta <- read.csv("data/wood_meta2.csv")

##### ----- Q20

# Read sample names
samples.out.Q20 <- rownames(SeqTabQ20)
samples.out.Q20.df <- as.data.frame(samples.out.Q20)

# Rename sample name column for metadata merge
colnames(samples.out.Q20.df) <- "id"
colnames(wood.meta) <-  c("id", "Plate", "Col", "Row", "Control", "Genotype",
                         "Disease", "Region", "MID.fwd.i5", "MID.rev.i7",
                         "MID.rev.i7.rc")

# Merge metadata
meta.Q20 <- merge(samples.out.Q20.df, wood.meta, by = "id", all.x = TRUE, all.y = TRUE)
rownames(meta.Q20) <- samples.out.Q20
meta.Q20 <- meta.Q20[,-1]

# Create phyloseq object
ps.Q20 <- phyloseq(otu_table(SeqTabQ20, taxa_are_rows=FALSE), sample_data(meta.Q20), tax_table(taxa.20))

# Create bar plot

Q20_plot <- plot_bar(ps.Q20, x="Phylum", facet_grid=~Disease) +
  labs(title = "Q20")

##### ----- Q28

# Read sample names
samples.out.Q28 <- rownames(SeqTabQ28)
samples.out.Q28.df <- as.data.frame(samples.out.Q28)

# Rename sample name column for metadata merge
colnames(samples.out.Q28.df) <- "id"
colnames(wood.meta) <-  c("id", "Plate", "Col", "Row", "Control", "Genotype",
                          "Disease", "Region", "MID.fwd.i5", "MID.rev.i7",
                          "MID.rev.i7.rc")

# Merge metadata
meta.Q28 <- merge(samples.out.Q28.df, wood.meta, by = "id", all.x = TRUE, all.y = TRUE)
rownames(meta.Q28) <- samples.out.Q28
meta.Q28 <- meta.Q28[,-1]

# Create phyloseq object
ps.Q28 <- phyloseq(otu_table(SeqTabQ28, taxa_are_rows=FALSE), sample_data(meta.Q28), tax_table(taxa.28))

# Create bar plot

Q28_plot <- plot_bar(ps.Q28, x="Phylum", facet_grid=~Disease) +
  labs(title = "Q28")

# View plots

Q20_plot
ggsave("outputs/Q20/phy_abundance_by_disease_Q20.png")
Q28_plot
ggsave("outputs/Q28/phy_abundance_by_disease_Q28.png")
general_plot <- plot_bar(ps.Q20, x="Phylum") +
  labs(title = "general")
general_plot

# Some ASVs were not assigned beyond the Fungi Kingdom
# -indicating that they might not be Fungi at all

# Separate with Phylum
taxa.20.phy <- taxa.20[!is.na(taxa.20[,2]),]
taxa.28.phy <- taxa.28[!is.na(taxa.28[,2]),]

# Separate with out Phylum
taxa.20.nophy <- taxa.20[is.na(taxa.20[,2]),]
taxa.28.nophy <- taxa.28[is.na(taxa.28[,2]),]

# Dimensions
dim(taxa.20)
dim(taxa.28)
dim(taxa.20.phy)
dim(taxa.28.phy)
dim(taxa.20.nophy)
dim(taxa.28.nophy)

# ASV stats
t20 <- nrow(taxa.20) # 1494 ASVs including ones with no Phylum
t28 <- nrow(taxa.28) # 929 ASVs including ones with no Phylum
p20 <- nrow(taxa.20.phy) # 1063 ASVs with Phylum
p28 <- nrow(taxa.28.phy) # 694 ASVs with Phylum
n20 <- nrow(taxa.20.nophy) # 431 ASVs with no Phylum
n28 <- nrow(taxa.28.nophy) # 235 ASVs with no Phylum

# About a quarter of our data has no Phylum assignment
percent.nophy.20 <- ((n20/t20)*100)
percent.nophy.28 <- ((n28/t28)*100)

# --- Work below this was just exploration of the data / r practice, will remove before script is complete
### Fungal stats Q20

# Unique clades represented
rep.phy.20 <- gsub("p__", "", (unique(taxa.20.phy[,2])), )
rep.cls.20 <- gsub("c__", "", (unique(taxa.20.phy[,3])), )
rep.ord.20 <- gsub("o__", "", (unique(taxa.20.phy[,4])), )
rep.fam.20 <- gsub("f__", "", (unique(taxa.20.phy[,5])), )
rep.gen.20 <- gsub("g__", "", (unique(taxa.20.phy[,6])), )
rep.spe.20 <- gsub("s__", "", (unique(taxa.20.phy[,7])), )

# Print
print(rep.phy.20)
print(rep.cls.20)
print(rep.ord.20)
print(rep.fam.20)
print(rep.gen.20)
print(rep.spe.20)

# Number of clades represented
Q20.clade.num <- list()
Q20.clade.num <- append(Q20.clade.num, length(rep.phy.20))
Q20.clade.num <- append(Q20.clade.num, length(rep.cls.20))
Q20.clade.num <- append(Q20.clade.num, length(rep.ord.20))
Q20.clade.num <- append(Q20.clade.num, length(rep.fam.20))
Q20.clade.num <- append(Q20.clade.num, length(rep.gen.20))
Q20.clade.num <- append(Q20.clade.num, length(rep.spe.20))
Q20.clade.num <- as.numeric(unlist(Q20.clade.num))
print(Q20.clade.num)

# Phylum Distribution
Q20.phy.dist <- list()
for (p in seq(1, Q20.clade.num[1])) {
  name <- paste0("p__", rep.phy.20[p])
  rcount <- sum(taxa.20.phy[,2] == name)
  print(rcount)
  Q20.phy.dist <- append(Q20.phy.dist, rcount)
}

Q20.phy.dist <- as.numeric(unlist(Q20.phy.dist))
phy.tab.dist <- rbind(rep.phy.20, Q20.phy.dist)
print(phy.tab.dist)

### Fungal stats Q28

# Unique clades represented
rep.phy.28 <- gsub("p__", "", (unique(taxa.28.phy[,2])), )
rep.cls.28 <- gsub("c__", "", (unique(taxa.28.phy[,3])), )
rep.ord.28 <- gsub("o__", "", (unique(taxa.28.phy[,4])), )
rep.fam.28 <- gsub("f__", "", (unique(taxa.28.phy[,5])), )
rep.gen.28 <- gsub("g__", "", (unique(taxa.28.phy[,6])), )
rep.spe.28 <- gsub("s__", "", (unique(taxa.28.phy[,7])), )

# Print
print(rep.phy.28)
print(rep.cls.28)
print(rep.ord.28)
print(rep.fam.28)
print(rep.gen.28)
print(rep.spe.28)

# Number of clades represented
Q28.clade.num <- list()
Q28.clade.num <- append(Q28.clade.num, length(rep.phy.28))
Q28.clade.num <- append(Q28.clade.num, length(rep.cls.28))
Q28.clade.num <- append(Q28.clade.num, length(rep.ord.28))
Q28.clade.num <- append(Q28.clade.num, length(rep.fam.28))
Q28.clade.num <- append(Q28.clade.num, length(rep.gen.28))
Q28.clade.num <- append(Q28.clade.num, length(rep.spe.28))
Q28.clade.num <- as.numeric(unlist(Q28.clade.num))
print(Q28.clade.num)

colnames()
index <- seq(1,6)
index <- append(index, seq(1,6))
index
index <- as.data.frame(index)
index
Q28.clade.num.df <- as.data.frame(Q28.clade.num)
Q20.clade.num.df <- as.data.frame(Q20.clade.num)
tax.compare <- data.frame(Q28.clade.num.df, Q20.clade.num.df)
tax.compare <- melt(tax.compare)
tax.compare <- data.frame(index, tax.compare)
tax.compare
ggplot(data = tax.compare, aes(x = index, y = value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs()
