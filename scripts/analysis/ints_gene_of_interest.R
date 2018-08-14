# This script selects the interactions between genes of interest and their neighbours.
# The list of genes of interest was provided by the domain expert.

# Load the data for the analysis
# Integrated dataset
integrated_int <- read.table(file = "~/absb/results/integration/integrated_int.txt", sep = "\t", stringsAsFactors = F)

# Load genes of interst
genes <- read.table(file = "~/absb/data/analysis_data/Int_nodes_gname.txt", sep = "\t", stringsAsFactors = F)
genes <- genes$V1

# Convert genes of interest to ensg IDs
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="ensembl.org")
name2ensg <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "gene_biotype"), filters=c("external_gene_name"), values=genes, mart=mart)
colnames(name2ensg) <- c("ensg", "gene_name", "biotype")
write.table(name2ensg, file = "~/absb/data/analysis_data/name2ensg_int_genes.txt", sep="\t", row.names=F, quote=F)

# Select only ensg IDs
ensg_int_genes <- unique(name2ensg$ensg)

# Subset interactions realted to the genes of interest

interactions_int_genes <- integrated_int[integrated_int$ensg.A%in%ensg_int_genes|integrated_int$ensg.B%in%ensg_int_genes, ]
save(interactions_int_genes, file="~/absb/results/integration/interactions_int_genes.RData")

# Add gene names
ensg_selected <- unique(c(interactions_int_genes$ensg.A,interactions_int_genes$ensg.B))
ensg2names_selected <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "gene_biotype"), filters=c("ensembl_gene_id"), values=ensg_selected, mart=mart)
colnames(ensg2names_selected) <- c("ensg", "gene_name", "biotype")

# Merge gene names with the selected interactions
selected_interactions <- merge(interactions_int_genes, ensg2names_selected, by.x = "ensg.A", by.y = "ensg", all.x=T)
selected_interactions <- merge(selected_interactions, ensg2names_selected, by.x = "ensg.B", by.y = "ensg", all.x=T)
head(selected_interactions)

# Rename columns.
selected_interactions <- selected_interactions[,c("gene_name.x","gene_name.y","score", "interaction_type", "data_source","ensg.A", "ensg.B", "biotype.x","biotype.y")] 
colnames(selected_interactions) <- c("gene_name.A","gene_name.B","score", "interaction_type", "data_source","ensg.A", "ensg.B", "biotype.A","biotype.B")

# Save to file
save(selected_interactions, file="~/absb/results/integration/selected_interactions_genes_of_interest.RData")
write.table(selected_interactions, file="~/absb/results/integration/selected_interactions_genes_of_interest.txt", sep="\t", quote=F, row.names=F)
