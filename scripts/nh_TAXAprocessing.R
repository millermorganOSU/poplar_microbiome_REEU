## Process Taxonomy

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

# Load assigned taxonomy (after host genome alignment nh = NO HOST *supposedly)
load("R/AssignedTaxa_nh.RData")

# Load sequence table
load("R/SeqTab_nh.RData")

# Load metadata
disease.scores <- read_xlsx("data/sample_diseasescores.xlsx")
wood.meta <- read.csv("data/wood_meta2.csv", row.names = 1)

# Some ASVs were not assigned beyond the Fungi Kingdom
# -indicating that they might not be Fungi at all

# Index taxa matrix
taxa.nh<- cbind(taxa.nh, seq(1,613))
colnames(taxa.nh) <- c("Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species", "i")

# Separate with Phylum
taxa.nh.phy <- taxa.nh[!is.na(taxa.nh[,2]),]

# Separate with out Phylum
taxa.nh.nophy <- taxa.nh[is.na(taxa.nh[,2]),]

# Dimensions
dim(taxa.nh)
dim(taxa.nh.phy)
dim(taxa.nh.nophy)

# ASV stats
tnh <- nrow(taxa.nh) # 613 ASVs including ones with no Phylum
pnh <- nrow(taxa.nh.phy) # 516 ASVs with Phylum
nnh <- nrow(taxa.nh.nophy) # 97 ASVs with no Phylum

# About 16 percent of out data has no Phylum assignment
percent.nophy.nh <- ((nnh/tnh)*100)

# Separate out ASVs with no Phylum
nophy.ASVs <- as.list(rownames(taxa.nh.nophy))
length(nophy.ASVs)
class(nophy.ASVs)

nophy.index <- taxa.nh.nophy[,8]
length(nophy.index)
class(nophy.index)

# Export to .fasta for clustering

### write.fasta(nophy.ASVs, nophy.index, "outputs/nophy_ASVs.fasta", open = "w", as.string = TRUE)

# run clustal-o.sh from clustal directory 

### $./clustal-o.sh -i ../outputs/align/nophy_ASVs.fasta -o ../outputs/align/clust_nophy_tree.fasta -v --output-order=tree-output

# Subset data based on clustal-o output

# Import likely host sequences
likely_poplar <- read.table(file = "outputs/align/likely_poplar_seq.txt")

# Create data frame for easy subsetting
taxa.nh.df <- as.data.frame(taxa.nh)

# Change index to numeric
glimpse(taxa.nh.df)
taxa.nh.df$i <- as.numeric(taxa.nh.df$i)
class(taxa.nh.df$i)

# Subset based on likely poplar sequence index
# 74 ASVs are likely host DNA
taxa.nh.probpop <- taxa.nh.df[taxa.nh.df$i %in% as.list(as.numeric((likely_poplar$V1))),]

prob.pop.ASVs <- as.list(rownames(taxa.nh.probpop))
prob.pop.index <- taxa.nh.probpop[,8] 

### write.fasta(prob.pop.ASVs, prob.pop.index, "outputs/prob_pop_ASVs.fasta", open = "w", as.string = TRUE)

# Unable to come up with consensus sequence
# Add 364 to UNITE basal on phylo-tree of likely poplar DNA
# Assign taxonomy again

load("R/AssignedTaxa_ptri.RData")
dim(taxa.ptri)

# Index taxa table
taxa.ptri<- cbind(taxa.ptri, seq(1,613))
colnames(taxa.ptri) <- c("Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species", "i")

# Separate with Phylum
taxa.ptri.phy <- taxa.ptri[!is.na(taxa.ptri[,2]),]
dim(taxa.ptri.phy) # 592 with Phylum

# Separate with out Phylum
taxa.ptri.nophy <- taxa.ptri[is.na(taxa.ptri[,2]),]
dim(taxa.ptri.nophy) # 21 with no Phylum

# ASV stats
ptri.tnh <- nrow(taxa.ptri) # 613 ASVs including ones with no Phylum
ptri.pnh <- nrow(taxa.ptri.phy) # 592 ASVs with Phylum
ptri.nnh <- nrow(taxa.ptri.nophy) # 21 ASVs with no Phylum

# About 3 percent of out data has no Phylum assignment
percent.nophy.ptri <- ((ptri.nnh/ptri.tnh)*100)

# 584 has no Kingdom assignment but blastable with blastn
# 301 and 365 have no Phylum and are only blastable with blastx

plantea <- taxa.ptri[taxa.ptri[,1]=="k__Plantae",]
plantea <- plantea[!is.na(plantea[,1]),] # Remove 584 with no assignment
dim(plantea) # 74 ASVs were expected to become Ptri, now 135 are assigned as Ptri

no.plantea <- taxa.ptri[taxa.ptri[,1]!="k__Plantae",]
no.plantea <- no.plantea[!is.na(no.plantea[,1]),] # Remove 584 with no assignment
dim(no.plantea) # 477 ASVs non-Plantea

# Separate out ASVs with no Phylum and export to .fasta

blast.ASVs <- as.list(rownames(taxa.ptri.nophy))
length(blast.ASVs)
class(blast.ASVs)

blast.index <- taxa.ptri.nophy[,8]
length(nophy.index)
class(nophy.index)

blast.ASVs
blast.index

### write.fasta(blast.ASVs, blast.index, "outputs/blast/blast.ASVs.fasta", open = "w", as.string = TRUE)

### Blast remotely $blastn -db nt -query blast.ASVs.fasta -remote -out results.bls -outfmt 6

blasted <- read.table("outputs/blast/results.bls")
colnames(blasted) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                       "qend", "sstart", "send", "evalue", "bitscore")

blasted.index <- unique(blasted[,1]) # 301 and 365 unassigned by blastn (need blastx manually)
blasted.index
taxa.ptri.nophy.test <- taxa.ptri.nophy
taxa.ptri.nophy.test <- cbind(taxa.ptri.nophy.test, Accession=NA)
dim(taxa.ptri.nophy.test)


taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==92,] = c("k__Viridiplantae", 
                                                        "p__Chlorophyta",
                                                        "c__Trebouxiophyceae",
                                                        "o__Trebouxiales",
                                                        "f__Trebouxiaceae", 
                                                        "g__Trebouxia",
                                                        "s__Trebouxia_sp_Incertae_sedis",
                                                        92,
                                                        "KJ754204.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==140,] = c("k__Fungi", 
                                                        "p__Fungi_phy_Incertae_sedis",
                                                        "c__Fungi_cls_Incertae_sedis",
                                                        "o__Fungi_ord_Incertae_sedis",
                                                        "f__Fungi_fam_Incertae_sedis", 
                                                        "g__Fungi_gen_Incertae_sedis",
                                                        "s__Fungi_ps__Incertae_sedis",
                                                        140,
                                                        "MW472276.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==152,] = c("k__Fungi", 
                                                         "p__Basidiomycota",
                                                         "c_Basidiomycota_cls_Incertae_sedis",
                                                         "o_Basidiomycota_ord_Incertae_sedis",
                                                         "f_Basidiomycota_fam_Incertae_sedis", 
                                                         "g_Basidiomycota_gen_Incertae_sedis",
                                                         "s_Basidiomycota_sp_Incertae_sedis",
                                                         152,
                                                         "AM901784.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==182,] = c("k__Fungi", 
                                                         "p__Basidiomycota",
                                                         "c_Basidiomycota_cls_Incertae_sedis",
                                                         "o_Basidiomycota_ord_Incertae_sedis",
                                                         "f_Basidiomycota_fam_Incertae_sedis", 
                                                         "g_Basidiomycota_gen_Incertae_sedis",
                                                         "s_Basidiomycota_sp_Incertae_sedis",
                                                         182,
                                                         "AM901784.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==190,] = c("k__Fungi", 
                                                         "p__Fungi_phy_Incertae_sedis",
                                                         "c__Fungi_cls_Incertae_sedis",
                                                         "o__Fungi_ord_Incertae_sedis",
                                                         "f__Fungi_fam_Incertae_sedis", 
                                                         "g__Fungi_gen_Incertae_sedis",
                                                         "s__Fungi_ps__Incertae_sedis",
                                                         190,
                                                         "MW164317.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==228,] = c("k__Fungi", 
                                                         "p__Ascomycota",
                                                         "c__Ascomycota_cls_Incertae_sedis",
                                                         "o__Ascomycota_ord_Incertae_sedis",
                                                         "f__Ascomycota_fam_Incertae_sedis", 
                                                         "g__Ascomycota_gen_Incertae_sedis",
                                                         "s__Ascomycota_sp_Incertae_sedis",
                                                         228,
                                                         "NR_166261.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==266,] = c("k__Fungi", 
                                                         "p__Fungi_phy_Incertae_sedis",
                                                         "c__Fungi_cls_Incertae_sedis",
                                                         "o__Fungi_ord_Incertae_sedis",
                                                         "f__Fungi_fam_Incertae_sedis", 
                                                         "g__Fungi_gen_Incertae_sedis",
                                                         "s__Fungi_ps__Incertae_sedis",
                                                         266,
                                                         "MT236898.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==319,] = c("k__Fungi", 
                                                         "p__Fungi_phy_Incertae_sedis",
                                                         "c__Fungi_cls_Incertae_sedis",
                                                         "o__Fungi_ord_Incertae_sedis",
                                                         "f__Fungi_fam_Incertae_sedis", 
                                                         "g__Fungi_gen_Incertae_sedis",
                                                         "s__Fungi_ps__Incertae_sedis",
                                                         319,
                                                         "MT236898.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==380,] = c("k__Fungi", 
                                                         "p__Fungi_phy_Incertae_sedis",
                                                         "c__Fungi_cls_Incertae_sedis",
                                                         "o__Fungi_ord_Incertae_sedis",
                                                         "f__Fungi_fam_Incertae_sedis", 
                                                         "g__Fungi_gen_Incertae_sedis",
                                                         "s__Fungi_ps__Incertae_sedis",
                                                         380,
                                                         "KP897592.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==381,] = c("k__Fungi", 
                                                         "p__Fungi_phy_Incertae_sedis",
                                                         "c__Fungi_cls_Incertae_sedis",
                                                         "o__Fungi_ord_Incertae_sedis",
                                                         "f__Fungi_fam_Incertae_sedis", 
                                                         "g__Fungi_gen_Incertae_sedis",
                                                         "s__Fungi_ps__Incertae_sedis",
                                                         381,
                                                         "MT236847.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==435,] = c("k__Fungi", 
                                                         "p__Fungi_phy_Incertae_sedis",
                                                         "c__Fungi_cls_Incertae_sedis",
                                                         "o__Fungi_ord_Incertae_sedis",
                                                         "f__Fungi_fam_Incertae_sedis", 
                                                         "g__Fungi_gen_Incertae_sedis",
                                                         "s__Fungi_ps__Incertae_sedis",
                                                         435,
                                                         "MK040530.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==436,] = c("k__Fungi", 
                                                         "p__Fungi_phy_Incertae_sedis",
                                                         "c__Fungi_cls_Incertae_sedis",
                                                         "o__Fungi_ord_Incertae_sedis",
                                                         "f__Fungi_fam_Incertae_sedis", 
                                                         "g__Fungi_gen_Incertae_sedis",
                                                         "s__Fungi_ps__Incertae_sedis",
                                                         436,
                                                         "KJ182679.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==447,] = c("k__Fungi", 
                                                         "p__Basidiomycota",
                                                         "c_Basidiomycota_cls_Incertae_sedis",
                                                         "o_Basidiomycota_ord_Incertae_sedis",
                                                         "f_Basidiomycota_fam_Incertae_sedis", 
                                                         "g_Basidiomycota_gen_Incertae_sedis",
                                                         "s_Basidiomycota_sp_Incertae_sedis",
                                                         447,
                                                         "AM901784.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==496,] = c("k__Fungi", 
                                                         "p__Basidiomycota",
                                                         "c_Basidiomycota_cls_Incertae_sedis",
                                                         "o_Basidiomycota_ord_Incertae_sedis",
                                                         "f_Basidiomycota_fam_Incertae_sedis", 
                                                         "g_Basidiomycota_gen_Incertae_sedis",
                                                         "s_Basidiomycota_sp_Incertae_sedis",
                                                         496,
                                                         "MT236898.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==554,] = c("k__Plantae", 
                                                         "p__Tracheophyta",
                                                         "c__Angiospermae",
                                                         "o__Malpighiales",
                                                         "f__Salicaceae", 
                                                         "g_Populus",
                                                         "s__Trichocarpa",
                                                         554,
                                                         "AC216848.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==566,] = c("k__Plantae", 
                                                         "p__Tracheophyta",
                                                         "c__Angiospermae",
                                                         "o__Malpighiales",
                                                         "f__Salicaceae", 
                                                         "g_Populus",
                                                         "s__Trichocarpa",
                                                         566,
                                                         "OX637686.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==570,] = c("k__Bacteria", 
                                                         "p__Bacteria_phy_Incertae_sedis",
                                                         "c__Bacteria_cls_Incertae_sedis",
                                                         "o__Bacteria_ord_Incertae_sedis",
                                                         "f__Bacteria_fam_Incertae_sedis", 
                                                         "g__Bacteria_gen_Incertae_sedis",
                                                         "s__Bacteria_sp_Incertae_sedis",
                                                         570,
                                                         "CP015243.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==575,] = c("k__Plantae", 
                                                         "p__Tracheophyta",
                                                         "c__Angiospermae",
                                                         "o__Malpighiales",
                                                         "f__Salicaceae", 
                                                         "g_Populus",
                                                         "s__Trichocarpa",
                                                         575,
                                                         "OX637690.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==584,] = c("k__Viridiplantae", 
                                                         "p__Streptophyta",
                                                         "c__Magnoliopsida",
                                                         "o__Malvales",
                                                         "f__Malvaceae", 
                                                         "g_Hibiscus",
                                                         "s__Hibiscus_tridactylites",
                                                         584,
                                                         "OY288291.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==301,] = c("k__Bacteria", 
                                                         "p__Pseudomonadota",
                                                         "c__Gammaproteobacteria",
                                                         "o__Enterobacterales",
                                                         "f__Pectobacteriaceae", 
                                                         "g_Brenneria",
                                                         "s__Brenneria_salicis",
                                                         301,
                                                         "WP_113866948.1")
taxa.ptri.nophy.test[taxa.ptri.nophy.test[,8]==365,] = c("k__Bacteria", 
                                                         "p__Pseudomonadota",
                                                         "c__Gammaproteobacteria",
                                                         "o__Enterobacterales",
                                                         "f__Pectobacteriaceae", 
                                                         "g_Brenneria",
                                                         "s__Brenneria_salicis",
                                                         365,
                                                         "WP_113866948.1")

# Manual commands for workflow

blasted.index
working <- blasted[blasted[,1]=="584",]
head(working)

# Now merge back with taxa.ptri to see where we stand

taxa.ptri <- cbind(taxa.ptri, Accession=NA)
taxa.ptri.df <- as.data.frame(taxa.ptri)
dim(taxa.ptri.df)

addphy.df <- as.data.frame(taxa.ptri.nophy.test)
dim(addphy.df)

taxa.V1 <- rows_update(taxa.ptri.df, addphy.df, by = "i")

# Separate with Phylum
V1.phy <- taxa.V1[!is.na(taxa.V1[,2]),]

# Separate with out Phylum
V1.nophy <- taxa.V1[is.na(taxa.V1[,2]),]

# Subset species NA for filling

wNA <- taxa.V1[is.na(taxa.V1[,7]),]
wNA.BU <- wNA
iNA <- wNA[,8]

# Filling loop #####

for (i in 1:nrow(wNA)) { # loop repeated for every row without species assignment
  if (is.na(wNA[i,2])) { # if row i phy. is empty, don't worry about it
    print("301 and 365 have no phylum assignment") # retroactivaley blastx-ed
  } else if (is.na(wNA[i,3])) { # if row i cls. is empty
    kng <- wNA[i,1] # king = king of row i
    kng <- sub('^k__', '', kng) # remove prefix
    phy <- wNA[i,2] # phy = phy of row i
    phy <- sub('^p__', '', phy) # remove prefix
    id <- wNA[i,8] # id = id of row i
    print(phy)
    print(kng)
    print(id)
    # replace the row with the id index
    wNA[wNA[,8]==id,] =  c( glue("k__{kng}"), 
                            glue("p__{phy}"),
                            glue("c__{phy}_cls_Incertae_sedis"),
                            glue("o__{phy}_ord_Incertae_sedis"),
                            glue("f__{phy}_fam_Incertae_sedis"), 
                            glue("g__{phy}_gen_Incertae_sedis"),
                            glue("s__{phy}_sp_Incertae_sedis"),
                            id,
                            "Filled")
  } else if (is.na(wNA[i,4])) { # if row i ord. is empty
    kng <- wNA[i,1] # king = king of row i
    kng <- sub('^k__', '', kng) # remove prefix
    phy <- wNA[i,2] # phy = phy of row i
    phy <- sub('^p__', '', phy) # remove prefix
    cls <- wNA[i,3] # cls = cls of row i
    cls <- sub('^c__', '', cls) # remove prefix
    id <- wNA[i,8] # id = id of row i
    print(phy)
    print(kng)
    print(cls)
    print(id)
    # replace the row with the id index
    wNA[wNA[,8]==id,] =  c( glue("k__{kng}"), 
                            glue("p__{phy}"),
                            glue("c__{cls}"),
                            glue("o__{cls}_ord_Incertae_sedis"),
                            glue("f__{cls}_fam_Incertae_sedis"), 
                            glue("g__{cls}_gen_Incertae_sedis"),
                            glue("s__{cls}_sp_Incertae_sedis"),
                            id,
                            "Filled")
  } else if (is.na(wNA[i,5])) { # if row i fam. is empty
    kng <- wNA[i,1] # king = king of row i
    kng <- sub('^k__', '', kng) # remove prefix
    phy <- wNA[i,2] # phy = phy of row i
    phy <- sub('^p__', '', phy) # remove prefix
    cls <- wNA[i,3] # cls = cls of row i
    cls <- sub('^c__', '', cls) # remove prefix
    ord <- wNA[i,4] # ord = ord of row i
    ord <- sub('^o__', '', ord) # remove prefix
    id <- wNA[i,8] # id = id of row i
    print(phy)
    print(kng)
    print(cls)
    print(ord)
    print(id)
    # replace the row with the id index
    wNA[wNA[,8]==id,] =  c( glue("k__{kng}"), 
                            glue("p__{phy}"),
                            glue("c__{cls}"),
                            glue("o__{ord}"),
                            glue("f__{ord}_fam_Incertae_sedis"), 
                            glue("g__{ord}_gen_Incertae_sedis"),
                            glue("s__{ord}_sp_Incertae_sedis"),
                            id,
                            "Filled")
  } else if (is.na(wNA[i,6])) { # if row i gen. is empty
    kng <- wNA[i,1] # king = king of row i
    kng <- sub('^k__', '', kng) # remove prefix
    phy <- wNA[i,2] # phy = phy of row i
    phy <- sub('^p__', '', phy) # remove prefix
    cls <- wNA[i,3] # cls = cls of row i
    cls <- sub('^c__', '', cls) # remove prefix
    ord <- wNA[i,4] # ord = ord of row i
    ord <- sub('^o__', '', ord) # remove prefix
    fam <- wNA[i,5] # fam = fam of row i
    fam <- sub('^f__', '', fam) # remove prefix
    id <- wNA[i,8] # id = id of row i
    print(phy)
    print(kng)
    print(cls)
    print(ord)
    print(fam)
    print(id)
    # replace the row with the id index
    wNA[wNA[,8]==id,] =  c( glue("k__{kng}"), 
                            glue("p__{phy}"),
                            glue("c__{cls}"),
                            glue("o__{ord}"),
                            glue("f__{fam}"), 
                            glue("g__{fam}_gen_Incertae_sedis"),
                            glue("s__{fam}_sp_Incertae_sedis"),
                            id,
                            "Filled")
  } else if (is.na(wNA[i,7])) { # if row i sp. is empty
    kng <- wNA[i,1] # king = king of row i
    kng <- sub('^k__', '', kng) # remove prefix
    phy <- wNA[i,2] # phy = phy of row i
    phy <- sub('^p__', '', phy) # remove prefix
    cls <- wNA[i,3] # cls = cls of row i
    cls <- sub('^c__', '', cls) # remove prefix
    ord <- wNA[i,4] # ord = ord of row i
    ord <- sub('^o__', '', ord) # remove prefix
    fam <- wNA[i,5] # fam = fam of row i
    fam <- sub('^f__', '', fam) # remove prefix
    gen <- wNA[i,6] # gen = gen of row i
    gen <- sub('^g__', '', gen) # remove prefix
    gen <- sub('_gen_Incertae_sedis$', '', gen) # remove suffix
    id <- wNA[i,8] # id = id of row i
    print(phy)
    print(kng)
    print(cls)
    print(ord)
    print(fam)
    print(gen)
    print(id)
    # replace the row with the id index
    wNA[wNA[,8]==id,] =  c( glue("k__{kng}"), 
                            glue("p__{phy}"),
                            glue("c__{cls}"),
                            glue("o__{ord}"),
                            glue("f__{fam}"), 
                            glue("g__{gen}_gen_Incertae_sedis"),
                            glue("s__{gen}_sp_Incertae_sedis"),
                            id,
                            "Filled")
}
}

### Create matrix with all assigned #####
taxa.V2 <- rows_patch(taxa.V1, wNA, by = "i")
#save(taxa.V2, file = "R/FilledTaxa.RData")
#load("R/FilledTaxa.RData")

### Check to make sure all ASVs assigned #####

# Separate with kingdom
V2.kng <- taxa.V2[!is.na(taxa.V2[,1]),]
dim(V2.kng)

# Separate with out class
V2.nokng <- taxa.V2[is.na(taxa.V2[,1]),]
dim(V2.nokng)

# Separate with class
V2.phy <- taxa.V2[!is.na(taxa.V2[,2]),]
dim(V2.phy)

# Separate with out class
V2.nophy <- taxa.V2[is.na(taxa.V2[,2]),]
dim(V2.nophy)

# Separate with class
V2.cls <- taxa.V2[!is.na(taxa.V2[,3]),]
dim(V2.cls)

# Separate with out class
V2.nocls <- taxa.V2[is.na(taxa.V2[,3]),]
dim(V2.nophy)

# Separate with order
V2.ord <- taxa.V2[!is.na(taxa.V2[,4]),]
dim(V2.ord)

# Separate with out order
V2.noord <- taxa.V2[is.na(taxa.V2[,4]),]
dim(V2.noord)

# Separate with family
V2.fam <- taxa.V2[!is.na(taxa.V2[,5]),]
dim(V2.fam)

# Separate with out family
V2.nofam<- taxa.V2[is.na(taxa.V2[,5]),]
dim(V2.nofam)

# Separate with genus
V2.gen <- taxa.V2[!is.na(taxa.V2[,6]),]
dim(V2.gen)

# Separate with out genus
V2.nogen<- taxa.V2[is.na(taxa.V2[,6]),]
dim(V2.nogen)

# Separate with species
V2.sp <- taxa.V2[!is.na(taxa.V2[,7]),]
dim(V2.sp)

# Separate with out species
V2.nosp<- taxa.V2[is.na(taxa.V2[,7]),]
dim(V2.nosp)

### All assigned!