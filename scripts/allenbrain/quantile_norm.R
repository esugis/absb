# This script is merging all 6 datasets based on the probesets names
# Then quantile normalization is applied over combined data.

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv 
maexp_178236545 <- read.csv(file = "~/absb/data/allenbrain/178236545_ds/MicroarrayExpression.csv", header = F)
dim(maexp_178236545)

# Read ontology
onto_178236545 <- read.csv(file = "~/absb/data/allenbrain/178236545_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178236545 <- read.csv("~/absb/data/allenbrain/178236545_ds/Probes.csv")

# Read in sample annotations
sampleann_178236545 <- read.csv("~/absb/data/allenbrain/178236545_ds/SampleAnnot.csv")
dim(sampleann_178236545)
rownames(maexp_178236545) <- maexp_178236545[,1]
maexp_178236545 <- maexp_178236545[,-1]
colnames(maexp_178236545) <- sampleann_178236545$structure_acronym

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238266 <- read.csv(file = "~/absb/data/allenbrain/178238266_ds/MicroarrayExpression.csv", header = F)
dim(maexp_178238266)

# Read ontology
onto_178238266 <- read.csv(file = "~/absb/data/allenbrain/178238266_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238266 <- read.csv("~/absb/data/allenbrain/178238266_ds/Probes.csv")

# Read in sample annotations
sampleann_178238266 <- read.csv("~/absb/data/allenbrain/178238266_ds/SampleAnnot.csv")
dim(sampleann_178238266)
rownames(maexp_178238266) <- maexp_178238266[,1]
maexp_178238266 <- maexp_178238266[,-1]
colnames(maexp_178238266) <- sampleann_178238266$structure_acrony 

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238316 <- read.csv(file = "~/absb/data/allenbrain/178238316_ds/MicroarrayExpression.csv", header = F)
dim(maexp_178238316)

# Read ontology
onto_178238316 <- read.csv(file = "~/absb/data/allenbrain/178238316_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238316 <- read.csv("~/absb/data/allenbrain/178238316_ds/Probes.csv")

# Read in sample annotations
sampleann_178238316 <- read.csv("~/absb/data/allenbrain/178238316_ds/SampleAnnot.csv")
dim(sampleann_178238316)
rownames(maexp_178238316) <- maexp_178238316[,1]
maexp_178238316 <- maexp_178238316[,-1]
colnames(maexp_178238316) <- sampleann_178238316$structure_acronym

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238359 <- read.csv(file = "~/absb/data/allenbrain/178238359_ds/MicroarrayExpression.csv", header = F)
dim(maexp_178238359)

# Read ontology
onto_178238359 <- read.csv(file = "~/absb/data/allenbrain/178238359_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238359 <- read.csv("~/absb/data/allenbrain/178238359_ds/Probes.csv")

# Read in sample annotations
sampleann_178238359 <- read.csv("~/absb/data/allenbrain/178238359_ds/SampleAnnot.csv")
dim(sampleann_178238359)
rownames(maexp_178238359) <- maexp_178238359[,1]
maexp_178238359 <- maexp_178238359[,-1]
colnames(maexp_178238359) <- sampleann_178238359$structure_acronym

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238373<- read.csv(file = "~/absb/data/allenbrain/178238373_ds/MicroarrayExpression.csv", header = F)
dim(maexp_178238373)

# Read ontology
onto_178238373<- read.csv(file = "~/absb/data/allenbrain/178238373_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238373<- read.csv("~/absb/data/allenbrain/178238373_ds/Probes.csv")

# Read in sample annotations
sampleann_178238373<- read.csv("~/absb/data/allenbrain/178238373_ds/SampleAnnot.csv")
dim(sampleann_178238373)
rownames(maexp_178238373) <- maexp_178238373[,1]
maexp_178238373<- maexp_178238373[,-1]
colnames(maexp_178238373) <- sampleann_178238373$structure_acronym

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv
maexp_178238387<- read.csv(file = "~/absb/data/allenbrain/178238387_ds/MicroarrayExpression.csv", header = F)
dim(maexp_178238387)

# Read ontology
onto_178238387<- read.csv(file = "~/absb/data/allenbrain/178238387_ds/MicroarrayExpression.csv", header = F)

# Read in probes names
probes_178238387<- read.csv("~/absb/data/allenbrain/178238387_ds/Probes.csv")

# Read in sample annotations
sampleann_178238387<- read.csv("~/absb/data/allenbrain/178238387_ds/SampleAnnot.csv")
dim(sampleann_178238387)
rownames(maexp_178238387) <- maexp_178238387[,1]
maexp_178238387<- maexp_178238387[,-1]
colnames(maexp_178238387) <- sampleann_178238387$structure_acronym

# Merge datasets based on the probesests
# Merge maexp_178236545 and maexp_178238266
ds1_2 <- merge(maexp_178236545[,1:501], maexp_178238266[,1:407], by="row.names", all=T)
row.names(ds1_2) <- ds1_2$Row.names
ds1_2 <- ds1_2[, -1]

# Merge result with maexp_178238316
ds1_2_3 <- merge(ds1_2[,1:908], maexp_178238316[,1:529], by="row.names", all=T)
row.names(ds1_2_3) <- ds1_2_3$Row.names
ds1_2_3 <- ds1_2_3[, -1]

# Merge result with maexp_178238359 
ds1_2_3_4 <- merge(ds1_2_3[,1:1437], maexp_178238359[,1:363], by="row.names", all=T)
row.names(ds1_2_3_4) <- ds1_2_3_4$Row.names
ds1_2_3_4 <- ds1_2_3_4[, -1]

# Merge result with maexp_178238373 
ds1_2_3_4_5 <- merge(ds1_2_3_4[,1:1800], maexp_178238373[,1:893], by="row.names", all=T)
row.names(ds1_2_3_4_5) <- ds1_2_3_4_5$Row.names
ds1_2_3_4_5 <- ds1_2_3_4_5[, -1]

# Merge with maexp_178238387
ds1_2_3_4_5_6 <- merge(ds1_2_3_4_5[,1:2693], maexp_178238387[,1:946], by="row.names", all=T)
row.names(ds1_2_3_4_5_6) <- ds1_2_3_4_5_6$Row.names
ds1_2_3_4_5_6 <- ds1_2_3_4_5_6[, -1]

dim(ds1_2_3_4_5_6) #58692  3639
save(ds1_2_3_4_5_6, file = "~/absb/results/allenbrain/combined_ds123456.Rdata")

# Convert expression to log2
ds1_2_3_4_5_6_log <- log2(ds1_2_3_4_5_6)
save(ds1_2_3_4_5_6_log, file = "~/absb/results/allenbrain/combined_ds123456_log2.Rdata")

# Apply quantile normalization over probesets
#Copy from the local computer file
# Addede part with checking for batch effect
#library(preprocessCore)
#norm_brain_mtrx <- normalize.quantiles(as.matrix(ds1_2_3_4_5_6_log))
#colnames(norm_brain_mtrx) <- colnames(ds1_2_3_4_5_6_log)
#norm_brain_mtrx <- cbind(probe_id=row.names(ds1_2_3_4_5_6_log), norm_brain_mtrx)

## Remove .x, .y, etc from the column names.
#colnames(norm_brain_mtrx) <- gsub("\\..*", "", colnames(norm_brain_mtrx))

## Melt structure
#norm_brain_mtrx <- data.frame(norm_brain_mtrx)
#colnames(norm_brain_mtrx) <- gsub("\\..*", "", colnames(norm_brain_mtrx))
#colnames(norm_brain_mtrx) <- gsub("^X", "", colnames(norm_brain_mtrx))
#rownames(norm_brain_mtrx)= norm_brain_mtrx$probe_id
#melt_brain <- melt(t(norm_brain_mtrx[,-1]))

## Assign colnames
#colnames(melt_brain) <- c("structure_acronym","probe_id","value")

#melt_brain$structure_acronym <- as.character(melt_brain$structure_acronym)
#melt_brain$value <- as.numeric(as.character(melt_brain$value))
#melt_brain$probe_id <- as.character(melt_brain$probe_id)
#melt_brain$structure_acronym <- gsub("\\..*", "", melt_brain$structure_acronym)



############################### From local file
library(preprocessCore)
norm_brain_mtrx <- normalize.quantiles(as.matrix(ds1_2_3_4_5_6_log))
colnames(norm_brain_mtrx) <- colnames(ds1_2_3_4_5_6_log)
norm_brain_mtrx <- cbind(probe_id=row.names(ds1_2_3_4_5_6_log), norm_brain_mtrx)

# Remove .x, .y, etc from the column names.
colnames(norm_brain_mtrx) <- gsub("\\..*", "", colnames(norm_brain_mtrx))
norm_brain_mtrx <- apply(norm_brain_mtrx,2, as.numeric)
tdata_nb <- t(norm_brain_mtrx)
tdata_nb <- cbind(study=row.names(tdata_nb),tdata_nb)
colnames(tdata_nb) <- tdata_nb[1,]
tdata_nb <- tdata_nb[-1,]
colnames(tdata_nb)[1] <- "study"
dim(tdata_nb)#[1]  3639 58693
tdata_nb[1:501,1] <- "1"
tdata_nb[502:908,1] <- "2"
tdata_nb[909:1437,1] <- "3"
tdata_nb[1438:1800,1] <- "4"
tdata_nb[1801:2693,1] <- "5"
tdata_nb[2694:3639,1] <- "6"
study <- tdata_nb[,1]
data.na.omit <- tdata_nb[,-1]
data.na.omit <- na.omit(data.na.omit)
data.na.omit <- apply(data.na.omit,2, as.numeric)
row.names(data.na.omit) <- row.names(tdata_nb)

# PCA
pca.all <- prcomp(data.na.omit, scale=T, center=T)
summary(pca.all)
samples <- rownames(data.na.omit)

plotdata<- as.data.frame(pca.all$x, stringsAsFactors = F)
library(ggplot2)
ggplot(plotdata,aes(x=PC1,y=PC2,label=samples, colour=study))+geom_point()+geom_text()+theme_bw()

# Conclusion - judging by the plot
# There is no obvious batch effect between the datasets

# Compute the average expression per brain region
# Melt ds
library(reshape2)
norm_brain_mtrx <- data.na.omit
melt_brain <- melt(norm_brain_mtrx)
colnames(melt_brain) <- c("structure_acronym","probe_id","value")

# Change variable types
melt_brain$structure_acronym <- as.character(melt_brain$structure_acronym)
melt_brain$value <- as.numeric(as.character(melt_brain$value))
melt_brain$probe_id <- as.character(melt_brain$probe_id)
length(unique(melt_brain$structure_acronym)) # should be 231
save(melt_brain, file="~/absb/results/melt_brain.RData")


library(tidyr)
#tst <- brain_tissue_mean[1:100,]
data_wide <- spread( brain_tissue_mean, structure_acronym,value)

# Make one-sided  t-test or wilcoxon test to find 
#genes that are expressed HIGHER than the mean over all brain regions
         
###wilcoxon test
p.val=c()
selcols <- c("CA1", "CA2", "CA3", "CA4", "SptN", "DG", "S")
data <- data_wide[, -1]
rownames(data) <- data_wide[,1]

for(i in 1: nrow(data)){
  u <- wilcox.test(unlist(data[i,c(colnames(data)%in%selcols)]),unlist(data[i,-c(colnames(data)%in%selcols)]), alternative = "greater")$p.value  
  p.val <- c(p.val, u)
}
p.val
p.val.adj <- p.adjust(p.val, method = "fdr")  
genes <- cbind(rownames(data), p.val.adj)
genes <- as.data.frame(genes)
genes[,1] <- as.character(genes[,1])
genes[,2] <- as.numeric(as.character(genes[,2]))
genes <- genes[genes$p.val.adj < 0.05,]
wlx_genes_means <- genes
save(wlx_genes_means, file="~/absb/results/allenbrain/wlx_genes_on_means.RData")
dim(wlx_genes_means)#  7095    2

# Compute z-scores per probe per tissue.
# Compute mean tissue expression in the relevant samples
brain_tissue_mean <- aggregate(value ~ probe_id + structure_acronym, melt_brain, FUN=mean, na.rm=TRUE, na.action="na.pass")

#Alltissues mean
all_mean <- mean(brain_tissue_mean$value, na.rm=TRUE)
all_sd <-sd(brain_tissue_mean$value, na.rm=TRUE)

# Compute zscores based on these average values. i.e compute mean over mean values per tissue and corresponding SD.
brain_tissue_mean <- cbind(brain_tissue_mean, zscore = ((brain_tissue_mean$value - all_mean) / all_sd))

# Save to file
save(brain_tissue_mean, file="~/absb/results/allenbrain/all_brain_tissue_zscores.RData")
write.table(brain_tissue_mean,file="~/absb/results/allenbrain/all_brain_tissue_zscores", sep="\t", row.names=F, quote=F)


