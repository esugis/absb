# This script reads in coexpression datasets in CA1, CA2, CA3, CA4, DG, SptN and subiculum brain regions
# Set working directory
setwd("~/absb/results/allenbrain/")

# Used liraries
library(reshape2)
library(biomaRt)

# Read in probeset original names from the dataset
probes <- read.csv("~/absb/data/allenbrain/178236545_ds/Probes.csv")

# Keep probeset_id and gene_symbol
p <- probes[,c(1,4)]
gene <- as.character(p$gene_symbol)
gene <- unique(gene)
length(gene) #current number of genes is 29131

# Convert gene ID using biomart
mart.g <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "ensembl.org")
g2ensg <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),filters=c("external_gene_name"), values = gene, mart = mart.g)

# Rename the columns
colnames(g2ensg) <- c(".id", "Target")

# Remove duplicates
g2ensg <- g2ensg [!duplicated(g2ensg), ]

# Non-converted IDs
gene_nc<-gene[!gene%in%g2ensg$.id]

# Merge with the original data
genes_aba_ensg <- merge(p, g2ensg, by.x = "gene_symbol", by.y = ".id", all = F)

# For debuging
#save(genes_aba_ensg, file = "genes_aba_ensg.RData")

# Add ensg IDs for DG co-expression dataset


