## Process Taxonomy

# Some ASVs were not assigned beyond the Fungi Kingdom
# -indicating that they might not be Fungi at all

# Set WD for or local machine
setwd('~/LeBoldus/local_git/poplar_microbiome_REEU')

# Load data
load("outputs/Q20/AssignedTaxa_Q20.RData")
load("outputs/Q28/AssignedTaxa_Q28.RData")

# Separate with Phylum
taxa.20.phy <- taxa.20[!is.na(taxa.20[,2]),]
taxa.28.phy <- taxa.28[!is.na(taxa.28[,2]),]

# Separate with out Phylum
taxa.20.nophy <- taxa.20[is.na(taxa.20[,2]),]
taxa.28.nophy <- taxa.28[is.na(taxa.28[,2]),]

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



