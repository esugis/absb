# This script computes zscores for each probe in each tissue based on global expression

# Read in the data 
# Set working directory
setwd("~/absb/results/allenbrain/")

# Used libaries
library(reshape)
#library(foreach)
#library(doParallel)
#registerDoParallel(cores=20)

# Path to the results
pathRdata <- "~/absb/results/allenbrain/tissues_rdata/acronym/"
dir.create(file.path(pathRdata), showWarnings = FALSE, recursive = TRUE)

# Load melted ds probe_id sample_id value  for 178236545
load(file = "~/absb/results/allenbrain/178236545_ds/maexp.probes_178236545_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178236545 <- maexp.probes

# Load processed data for 178238266
load(file = "~/absb/results/allenbrain/178238266_ds/maexp.probes_178238266_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238266 <- maexp.probes

# Load processed data for 178238316
load(file = "~/absb/results/allenbrain/178238316_ds/maexp.probes_178238316_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238316 <- maexp.probes

# Load processed data for 178238359
load(file = "~/absb/results/allenbrain/178238359_ds/maexp.probes_178238359_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238359 <- maexp.probes

# Load processed data for 178238373
load(file = "~/absb/results/allenbrain/178238373_ds/maexp.probes_178238373_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238373 <- maexp.probes

# Load processed data for 178238387
load(file = "~/absb/results/allenbrain/178238387_ds/maexp.probes_178238387_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238387 <- maexp.probes

# Select unique tissues in all 6 data sets
s1 <- read.csv(file = "~/absb/data/allenbrain/178236545_ds/SampleAnnot.csv")
s2 <- read.csv(file = "~/absb/data/allenbrain/178238266_ds/SampleAnnot.csv")
s3 <- read.csv(file = "~/absb/data/allenbrain/178238316_ds/SampleAnnot.csv")
s4 <- read.csv(file = "~/absb/data/allenbrain/178238359_ds/SampleAnnot.csv")
s5 <- read.csv(file = "~/absb/data/allenbrain/178238373_ds/SampleAnnot.csv")
s6 <- read.csv(file = "~/absb/data/allenbrain/178238387_ds/SampleAnnot.csv")

t11 <- s1[,c("structure_id", "structure_acronym")]
t22 <- s2[,c("structure_id", "structure_acronym")]
t33 <- s3[,c("structure_id", "structure_acronym")]
t44 <- s4[,c("structure_id", "structure_acronym")]
t55 <- s5[,c("structure_id", "structure_acronym")]
t66 <- s6[,c("structure_id", "structure_acronym")]

# Combine all tissues
tissues <- rbind(t11, t22, t33, t44, t55, t66)
tissues <- tissues[!duplicated(tissues),]

# List of unique tissues in all 6 data sets
tissues_all <- unique(c(as.character(t11$structure_acronym),as.character(t22$structure_acronym),as.character(t33$structure_acronym),as.character(t44$structure_acronym),as.character(t55$structure_acronym),as.character(t66$structure_acronym)))

# Number of tissues
length(tissues_all)#231

# Merge datasets and tissue acronyms by tissue id.
maexp_178236545 <- merge(maexp_178236545, t11, by.x="structure_id", by.y="structure_id", all=F)
#maexp_178236545 <- maexp_178236545[, -1]
maexp_178238266 <- merge(maexp_178238266, t22, by.x="structure_id", by.y="structure_id", all=F)
#maexp_178238266 <- maexp_178238266[, -1]
maexp_178238316 <- merge(maexp_178238316, t33, by.x="structure_id", by.y="structure_id", all=F)
#maexp_178238316 <- maexp_178238316[, -1]
maexp_178238359 <- merge(maexp_178238359, t44, by.x="structure_id", by.y="structure_id", all=F)
#maexp_178238359 <- maexp_178238359[, -1]
maexp_178238373 <- merge(maexp_178238373, t55, by.x="structure_id", by.y="structure_id", all=F)
#maexp_178238373 <- maexp_178238373[, -1]
maexp_178238387 <- merge(maexp_178238387, t66, by.x="structure_id", by.y="structure_id", all=F)
#maexp_178238387 <- maexp_178238387[, -1]

# Combine datasets
brain_data <- rbind(maexp_178236545,maexp_178238266, maexp_178238316, maexp_178238359, maexp_178238373, maexp_178238387)
brain_data$structure_acronym <- as.character(brain_data$structure_acronym)
# Remove duplicates 
brain_data <- brain_data[!duplicated(brain_data),]
# Save joint dataset
save(brain_data, file="brain_data.RData")
brain_data_tis=brain_data
library(reshape2))
# Apply quantile normalization for the joint dataset
# Add tissue id for casting the matrix for quanlite normalization
#brain_data_tis <- merge(brain_data, tissues, by.x= "structure_acronym", by.y = "structure_acronym", all=T)
# Remove duplicates
#brain_data_tis$structure_acronym <- as.character(brain_data_tis$structure_acronym)
#brain_data_tis <- brain_data_tis[!duplicated(brain_data_tis),]


# Cast a matrix
brain_mtrx <- dcast(brain_data_tis, probe_id ~ sample_num + structure_acronym + structure_id,mean,value.var="value")

library(preprocessCore)
norm_brain_mtrx <- normalize.quantiles(as.matrix(brain_mtrx[,2:length(colnames(brain_mtrx))]))
colnames(norm_brain_mtrx) <- colnames(brain_mtrx)[-1]
norm_brain_mtrx <- cbind(probe_id=brain_mtrx$probe_id, norm_brain_mtrx)

# Melt
norm_brain_mtrx <- data.frame(norm_brain_mtrx)
rownames(norm_brain_mtrx)= norm_brain_mtrx$probe_id
melt_brain <- melt(t(norm_brain_mtrx[,-1]))

# check the colnames
colnames(melt_brain) <- c("structure_acronym","probe_id","value")

melt_brain$structure_acronym <- gsub("^[^\\s_]+_", "\\1", melt_brain$structure_acronym, perl=T)
melt_brain$structure_acronym <- gsub( "_.*", "", melt_brain$structure_acronym, perl=T)

# Compute mean tissue expression in the relevant samples
# check how to aggregate over 2 columns
brain_tissue_mean <- aggregate(value ~ probe_id + structure_acronym, melt_brain, mean)

#Alltissues mean
all_mean <- mean(brain_tissue_mean[,2])
all_sd <-sd(brain_tissue_mean[,2])

# Check the formula!
# Compute zscores based on these average values. i.e compute mean over mean values per tissue and corresponding SD.
brain_tissue_mean <- cbind(brain_tissue_mean, zscore = ((brain_tissue_mean[,2] - all_mean) / all_sd))

# Save to file
save(brain_tissue_mean, file="all_brain_tissue_zscores.RData")
write.table(brain_tissue_mean,file="all_brain_tissue_zscores", sep="\t", row.names=F, quote=F)
