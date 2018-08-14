# This script calculates z-scores in Allen Brain Atlas datsets based on tissue acronym.
# Tissue IDs from left and right hemisphear corresponds to the same tissue acronym.

# Set working directory
setwd("~/absb/results/allenbrain/")

# Used libaries
library(reshape)
library(foreach)
library(doParallel)
registerDoParallel(cores=2)
	
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

# List of unique tissues in all 6 data sets
tissues_all <- unique(c(as.character(t11$structure_acronym),as.character(t22$structure_acronym),as.character(t33$structure_acronym),as.character(t44$structure_acronym),as.character(t55$structure_acronym),as.character(t66$structure_acronym)))

# Number of tissues
length(tissues_all)#231

# Merge datasets and tissue acronyms by tissue id.
maexp_178236545 <- merge(maexp_178236545, t11, by.x="structure_id", by.y="structure_id", all=F)
maexp_178236545 <- maexp_178236545[, -1]

maexp_178238266 <- merge(maexp_178238266, t22, by.x="structure_id", by.y="structure_id", all=F)
maexp_178238266 <- maexp_178238266[, -1]

maexp_178238316 <- merge(maexp_178238316, t33, by.x="structure_id", by.y="structure_id", all=F)
maexp_178238316 <- maexp_178238316[, -1]

maexp_178238359 <- merge(maexp_178238359, t44, by.x="structure_id", by.y="structure_id", all=F)
maexp_178238359 <- maexp_178238359[, -1]

maexp_178238373 <- merge(maexp_178238373, t55, by.x="structure_id", by.y="structure_id", all=F)
maexp_178238373 <- maexp_178238373[, -1]

maexp_178238387 <- merge(maexp_178238387, t66, by.x="structure_id", by.y="structure_id", all=F)
maexp_178238387 <- maexp_178238387[, -1]



# Save the list of tissues
save(tissues_all,file = "all_tissues_acronym.RData")

# Compute z-scores in each of the tissues
#foreach(i = 1:length(tissues_all))%dopar%{

for(i in 1:length(tissues_all)){
# Subset one tissue
tissue <- tissues_all[i]
print("current tissue")
print(tissue)

# Subset tissue related data from each dataset
onetissue1 <- subset(maexp_178236545, structure_acronym==tissue)
onetissue2 <- subset(maexp_178238266, structure_acronym==tissue)
onetissue3 <- subset(maexp_178238316, structure_acronym==tissue)
onetissue4 <- subset(maexp_178238359, structure_acronym==tissue)
onetissue5 <- subset(maexp_178238373, structure_acronym==tissue)
onetissue6 <- subset(maexp_178238387, structure_acronym==tissue)

# Bind tissue related data from each ds together
onetissue <- rbind(onetissue1,onetissue2,onetissue3,onetissue4,onetissue5,onetissue6)

# Aggregate based on probe name, take mean
onetissue_av <- merge(aggregate(value ~ probe_id, onetissue, mean),tissue)

# Calculate z-score over one tissue
onetissue_av <- cbind(onetissue_av, zscore = ((onetissue_av[,2] - mean(onetissue_av[,2])) / sd(onetissue_av[,2])))

# Create path to to save results for one tissue
filedata <- sprintf("%s.RData",tissue);
pathdata <- file.path(pathRdata, filedata);

# Write file for each tissue
save(onetissue_av,file = pathdata)

}



