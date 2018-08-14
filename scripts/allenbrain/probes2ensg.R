# Convert the probe names to ensg.

# Set created directory as working dirrectory
setwd("~/absb/results/allenbrain/")

# Used liraries
library(reshape2)
library(gProfileR)

# Ontologies and probe annotations are the same for 6 microarray data sets, so read files from one as below
onto <- read.csv(file = "~/absb/data/allenbrain/178236545_ds/Ontology.csv")
probes <- read.csv("~/absb/data/allenbrain/178236545_ds/Probes.csv")

# Keep probe_id and entrez_id
p <- probes[,c(1,6)]
gene <- as.numeric(p$entrez_id)
gene <- unique(gene)
length(gene)

# Convert probes to ENSG IDs
p2ensg1 <- gconvert(gene[1:5000])
save(p2ensg1, file = "p2ensg1.RData")
p2ensg2 <- gconvert(gene[5001:10000])
save(p2ensg2, file = "p2ensg2.RData")
p2ensg3 <- gconvert(gene[10001:15000])
save(p2ensg3, file = "p2ensg3.RData")
p2ensg4 <- gconvert(gene[15001:20000])
save(p2ensg4, file = "p2ensg4.RData")
p2ensg5 <- gconvert(gene[20001:20788])
save(p2ensg5, file = "p2ensg5.RData")

# Combine all converted names
p2ensg <- rbind(p2ensg1,p2ensg2,p2ensg3,p2ensg4,p2ensg5)
p2ensg <- p2ensg[, c(2,4)] 
colnames(p2ensg) <- c(".id", "Target")
# Size of the converted data frame
dim(p2ensg)

# Remove duplicates
p2ensg <- p2ensg[!duplicated(p2ensg),]
dim(p2ensg)

p$gene_symbol <- as.character(p$entrez_id)
p2ensg$Target <- as.character(p2ensg$Target)
p2ensg$.id <- as.character(p2ensg$.id)
p2ensg_orig <- merge(p, p2ensg, by.x = "entrez_id", by.y = ".id", all = F)
p2ensg_orig <- p2ensg_orig[!duplicated(p2ensg_orig),]

# Save probes and corresponding ENSG IDs to the file
save(p2ensg, file = "p2ensg_all.RData")


# Load z-score matrix (cols-tissues, rows-probe ids)
load(file = "~/absb/results/allenbrain/tissues_zscores_mtx.RData")

# Add ENSG ids to tissues_zscores_mtx
tissues_zscores_mtx <- cbind(probe_id = row.names(tissues_zscores_mtx),tissues_zscores_mtx)
tissues_zscores_mtx <- data.frame(tissues_zscores_mtx)
tissues_zscores_mtx <- data.frame(sapply(tissues_zscores_mtx,function(x) as.numeric(as.vector(x))))
# Size of the data
dim(tissues_zscores_mtx)
tissues_zscores_mtx_ensg <- merge(tissues_zscores_mtx, p2ensg_orig, by.x = "probe_id", by.y = "probe_id", all = F)
tissues_zscores_mtx_ensg <- tissues_zscores_mtx_ensg[!duplicated(tissues_zscores_mtx_ensg),]

# Save to the file 
save(tissues_zscores_mtx_ensg, file = "tissues_zscores_mtx_ensg.RData")
write.table(tissues_zscores_mtx_ensg, file = "tissues_zscores_mtx_ensg.txt", sep = "\t",quote = F, row.names = F)

# Size of the filnal data
dim(tissues_zscores_mtx_ensg)

