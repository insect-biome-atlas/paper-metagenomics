# Script for generating table matching metagenomic and
# metabarcoding results
library(data.table)
library(reshape2)

# Read in the metagenomics results
G <- read.delim("metagenomics_genus_occurrence.tsv")

# Get metabarcoding reads for the same 40 samples
# TODO: change directory to your local copy of the IBA utils repo
cur_dir <- getwd()
setwd("~/dev/ms-repos-iba/utils/")
source("get_iba_co1_data_fxn.R")
setwd(cur_dir)

# Set IBA data paths
# TODO: Change to your local paths to the IBA data
data_path <- "~/dev/figshare-repos/iba/processed_data/v2/"
metadata_path <- "~/dev/figshare-repos/iba/raw_data/v4/"

# Get IBA lysate and homogenate data
# Note that these are data tables and not standard R data frames
L <- get_iba_co1_data(data_path, metadata_path, "SE", "lysate")
H <- get_iba_co1_data(data_path, metadata_path, "SE", "homogenate")

# Extract the 40 samples we are interested in
# Convert the data to standard R data frames in the process
M <- read.delim(paste0(metadata_path,"CO1_sequencing_metadata_SE.tsv"))
field_samples <- (colnames(G))[-1]

lysate_samples <- M$sampleID_NGI[M$sampleID_FIELD %in% field_samples & M$dataset=="CO1_lysate_2019_SE"]
homogenate_samples <- M$sampleID_NGI[M$sampleID_FIELD %in% field_samples & M$dataset=="CO1_homogenate_2019_SE"]

index <- which(colnames(L) %in% lysate_samples)
index <- c(1,2,6:13,index)
L <- L[,..index]
index <- 11:ncol(L)
row_sums <- rowSums(L[,..index])
L <- L[row_sums>0,]
L <- data.frame(L)
x <- colnames(L)[index]
x <- M$sampleID_FIELD[match(x,M$sampleID_NGI)]
colnames(L) <- c(colnames(L)[1:10],x)

index <- which(colnames(H) %in% homogenate_samples)
index <- c(1,2,6:13,index)
H <- H[,..index]
index <- 11:ncol(H)
row_sums <- rowSums(H[,..index])
H <- H[row_sums>0,]
H <- data.frame(H)
x <- colnames(H)[index]
x <- M$sampleID_FIELD[match(x,M$sampleID_NGI)]
colnames(H) <- c(colnames(H)[1:10],x)

# Correct some genus names
# See Iwaszkiewicz et al. for background on these cluster
# annotations
L$Genus[L$cluster=="Noctuidae_cluster19"] <- "Charanyca"
H$Genus[H$cluster=="Noctuidae_cluster19"] <- "Charanyca"

L$Genus[L$cluster=="Crambidae_cluster3"] <- "Agriphila"
H$Genus[H$cluster=="Crambidae_cluster3"] <- "Agriphila"

# Write the resulting tables
write.table(L, "iba_lysate_cluster_counts.tsv", sep="\t", row.names=FALSE)
write.table(H, "iba_homogenate_cluster_counts.tsv", sep="\t", row.names=FALSE)

# Compute genus counts tables
lysate_genus <- unique(L$Genus)
homogenate_genus <- unique(H$Genus)

L1 <- data.frame(L[FALSE,c(1:8,11:50)])
for (genus in lysate_genus) {
    temp <- L[L$Genus==genus,c(1:8,11:50)]
    x <- colSums(temp[,9:48])
    L1 <- rbind(L1, c(temp[1,1:8],x))
}

H1 <- data.frame(H[FALSE,])
for (genus in homogenate_genus) {
    temp <- H[H$Genus==genus,c(1:8,11:50)]
    x <- colSums(temp[,9:48])
    H1 <- rbind(H1, c(temp[1,1:8],x))
}

# Write the resulting tables
write.table(L1, "iba_lysate_genus_counts.tsv", sep="\t", row.names=FALSE)
write.table(H1, "iba_homogenate_genus_counts.tsv", sep="\t", row.names=FALSE)

# Make box plots
L2 <- reshape2::melt(L1,id.vars=1:8,measure.vars=9:48,value.name="reads")
H2 <- reshape2::melt(H1,id.vars=1:8,measure.vars=9:48,value.name="reads")
G2 <- reshape2::melt(G,id.vars=1:5,measure.vars=6:45,value.name="present")

DL <- merge(G2,L2)
DL <- DL[DL$reads > 0,]     # If not detected with metabarcoding, assume it is an error
DL <- DL[DL$Genus!="Drosophila",]       # This corresponds to spikeins removed in metabarcoding
DH <- merge(G2,H2)
DH <- DH[DH$reads > 0,]     # If not detected with metabarcoding, assume it is an error
DH <- DH[DH$Genus!="Poecilobothrus",]   # This is an obvious mismatch
DH <- DH[DH$Genus!="Drosophila",]       # This corresponds to spikeins removed in metabarcoding

#DL$reads <- DL$reads + 1    # Avoid -inf log value for 0 reads (if 0-read entries not removed above)
#DH$reads <- DH$reads + 1    # Avoid -inf log value for 0 reads (if 0-read entries not removed above)

# Write the resulting tables
write.table(DL, "lysate_match_counts.tsv", sep="\t", row.names=FALSE)
write.table(DH, "homogenate_match_counts.tsv", sep="\t", row.names=FALSE)

# We focus below on the homogenate data, as the metagenomic analysis
# used the homogenate from these samples

# Run a T test on homogenate data
t.test(x=DH$reads[DH$present==1],y=DH$reads[DH$present==0])

# Generate a box plot for homogenate data
pdf("Fig_SX.pdf")
boxplot(DH$reads~DH$present,log="y")
dev.off()
