# This script computes co-expression in >= 400 brain tissues from Allen Brain Atlas.
# Co-expression is computed over 6 datasets of gene expression in healthy brains.

# Set working directory
setwd("~/absb/results/allenbrain/")

# Used libaries
library(reshape)
library(foreach)
        
# Path to the results  
pathRdata <- "~/absb/results/allenbrain/tissues_rdata"
dir.create(file.path(pathRdata), showWarnings = FALSE, recursive = TRUE)


# Load melted ds probe_id sample_id value  for 178236545
load(file = "~/absb/results/allenbrain/178236545_ds/maexp.probes_178236545_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178236545 <- maexp.probes
t1 <- unique(maexp_178236545[,4])
print("tissue 1")
length(t1)

# Load processed data for 178238266
load(file = "~/absb/results/allenbrain/178238266_ds/maexp.probes_178238266_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238266 <- maexp.probes
t2 <- unique(maexp_178238266[,4])
print("tissue 2")
length(t2)

# Load processed data for 178238316
load(file = "~/absb/results/allenbrain/178238316_ds/maexp.probes_178238316_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238316 <- maexp.probes
t3 <- unique(maexp_178238316[,4])
print("tissue 3")
length(t3)

# Load processed data for 178238359
load(file = "~/absb/results/allenbrain/178238359_ds/maexp.probes_178238359_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238359 <- maexp.probes
t4 <- unique(maexp_178238359[,4])
print("tissue 4")
length(t4)

# Load processed data for 178238373
load(file = "~/absb/results/allenbrain/178238373_ds/maexp.probes_178238373_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238373 <- maexp.probes
t5 <- unique(maexp_178238373[,4])
print("tissue 5")
length(t5)

# Load processed data for 178238387
load(file = "~/absb/results/allenbrain/178238387_ds/maexp.probes_178238387_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238387 <- maexp.probes
t6 <- unique(maexp_178238387[,4])
print("tissue 6")
length(t6)


# Select unique tissues in all 6 data sets
s1 <- read.csv(file = "~/absb/data/allenbrain/178236545_ds/SampleAnnot.csv")
s2 <- read.csv(file = "~/absb/data/allenbrain/178238266_ds/SampleAnnot.csv")
s3 <- read.csv(file = "~/absb/data/allenbrain/178238316_ds/SampleAnnot.csv")
s4 <- read.csv(file = "~/absb/data/allenbrain/178238359_ds/SampleAnnot.csv")
s5 <- read.csv(file = "~/absb/data/allenbrain/178238373_ds/SampleAnnot.csv")
s6 <- read.csv(file = "~/absb/data/allenbrain/178238387_ds/SampleAnnot.csv")

# A number of tissues relevant to the disease were selected based on the domain expert input
# Allen Brain atlas brain tissue ontology was used for defning the regions and their IDs
# DG, CA1, CA2, CA3, CA4, S (subiculum), SptN (septal nuclei)
# Tissue IDs for left and right part are separated.
# They will be grouped together in individual correlation analysis.

selected_tissues<-c(4258, 4267,4254, 4263, 4255, 4264, 4256, 4265, 4257, 4266, 4251, 4260, 4301, 4304)

tst1 <- maexp_178236545
tst1 <- subset(tst1, structure_id%in%selected_tissues)
tst2 <- maexp_178238266
tst2 <- subset(tst2, structure_id%in%selected_tissues)
tst3 <- maexp_178238316
tst3 <- subset(tst3, structure_id%in%selected_tissues)
tst4 <- maexp_178238359
tst4 <- subset(tst4, structure_id%in%selected_tissues)
tst5 <- maexp_178238373
tst5 <- subset(tst5, structure_id%in%selected_tissues)
tst6 <- maexp_178238387
tst6 <- subset(tst6, structure_id%in%selected_tissues)

# Temporarily
rm(maexp_178236545)
rm(maexp_178238266)
rm(maexp_178238316)
rm(maexp_178238359)
rm(maexp_178238373)
rm(maexp_178238387)

# Save preprocessed datasets
save(tst1, tst2, tst3, tst4, tst5, tst6, file="preprocessed_data.RData")

# Subset one tissue
#tissue <- selected_tissues[i]

# DG
tissue <- "DG"
tissue_id <-c(4258, 4267)
print("current tissue")
print(tissue)

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)


# Subset tissue related data from each dataset
tst_onetissue1 <- subset(tst1, structure_id%in%tissue_id)
tst_onetissue2 <- subset(tst2, structure_id%in%tissue_id)
tst_onetissue3 <- subset(tst3, structure_id%in%tissue_id)
tst_onetissue4 <- subset(tst4, structure_id%in%tissue_id)
tst_onetissue5 <- subset(tst5, structure_id%in%tissue_id)
tst_onetissue6 <- subset(tst6, structure_id%in%tissue_id)

cast_tst1 <- data.frame(cast(tst_onetissue1, probe_id~sample_num))
cast_tst2 <- data.frame(cast(tst_onetissue2, probe_id~sample_num))
cast_tst3 <- data.frame(cast(tst_onetissue3, probe_id~sample_num))
cast_tst4 <- data.frame(cast(tst_onetissue4, probe_id~sample_num))
cast_tst5 <- data.frame(cast(tst_onetissue5, probe_id~sample_num))
cast_tst6 <- data.frame(cast(tst_onetissue6, probe_id~sample_num))

# Merge individual datasets
cast12 <- merge(cast_tst1, cast_tst2, by.x="probe_id", by.y="probe_id", all=T)
cast123 <- merge(cast12, cast_tst3, by.x="probe_id", by.y="probe_id", all=T)
cast1234 <- merge(cast123, cast_tst4, by.x="probe_id", by.y="probe_id", all=T)
cast12345 <- merge(cast1234, cast_tst5, by.x="probe_id", by.y="probe_id", all=T)
cast123456 <- merge(cast12345, cast_tst6, by.x="probe_id", by.y="probe_id", all=T)

# Exclude probesets that are not changing
rownames_cast123456 <- cast123456$probe_id 
dim(cast123456)
cast123456 <-cast123456[, -1]

library(preprocessCore)
# Exclude probesets that are not changing
cast123456 <- normalize.quantiles(as.matrix(cast123456))
rownames(cast123456) <- rownames_cast123456
SD <- apply(cast123456, 1, sd, na.rm = T)
cast123456_filt <- cast123456[SD >= 0.29, ]
dim(cast123456_filt)

# Remove unnesessary data structures 
rm(cast_tst1,cast_tst2, cast_tst3,
     cast_tst4, cast_tst5, cast_tst6, 
     cast12, cast123, cast1234, cast12345, cast123456)

# Transpose matrix and compute coexpression and correlation p-value
m <- t(cast123456_filt)

# Write example data
file_name <- "DG_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))

# CA1
tissue <- "CA1"
tissue_id <-c(4254,4263)
print("current tissue")
print(tissue)

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)

# Subset tissue related data from each dataset
tst_onetissue1 <- subset(tst1, structure_id%in%tissue_id)
tst_onetissue2 <- subset(tst2, structure_id%in%tissue_id)
tst_onetissue3 <- subset(tst3, structure_id%in%tissue_id)
tst_onetissue4 <- subset(tst4, structure_id%in%tissue_id)
tst_onetissue5 <- subset(tst5, structure_id%in%tissue_id)
tst_onetissue6 <- subset(tst6, structure_id%in%tissue_id)

cast_tst1 <- data.frame(cast(tst_onetissue1, probe_id~sample_num))
cast_tst2 <- data.frame(cast(tst_onetissue2, probe_id~sample_num))
cast_tst3 <- data.frame(cast(tst_onetissue3, probe_id~sample_num))
cast_tst4 <- data.frame(cast(tst_onetissue4, probe_id~sample_num))
cast_tst5 <- data.frame(cast(tst_onetissue5, probe_id~sample_num))
cast_tst6 <- data.frame(cast(tst_onetissue6, probe_id~sample_num))

# Merge individual datasets
cast12 <- merge(cast_tst1, cast_tst2, by.x="probe_id", by.y="probe_id", all=T)
cast123 <- merge(cast12, cast_tst3, by.x="probe_id", by.y="probe_id", all=T)
cast1234 <- merge(cast123, cast_tst4, by.x="probe_id", by.y="probe_id", all=T)
cast12345 <- merge(cast1234, cast_tst5, by.x="probe_id", by.y="probe_id", all=T)
cast123456 <- merge(cast12345, cast_tst6, by.x="probe_id", by.y="probe_id", all=T)

# Exclude probesets that are not changing
rownames_cast123456 <- cast123456$probe_id
dim(cast123456)
cast123456 <-cast123456[, -1]

library(preprocessCore)
# Exclude probesets that are not changing
cast123456 <- normalize.quantiles(as.matrix(cast123456))
rownames(cast123456) <- rownames_cast123456
SD <- apply(cast123456, 1, sd, na.rm = T)
cast123456_filt <- cast123456[SD >= 0.29, ]
dim(cast123456_filt)

# Remove unnesessary data structures
rm(cast_tst1,cast_tst2, cast_tst3,
     cast_tst4, cast_tst5, cast_tst6,
     cast12, cast123, cast1234, cast12345, cast123456)

# Transpose matrix and compute coexpression and correlation p-value
m <- t(cast123456_filt)

# Write example data
file_name <- "CA1_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))


# CA2
tissue <- "CA2"
tissue_id <-c(4255, 4264)

# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)

# Subset tissue related data from each dataset
tst_onetissue1 <- subset(tst1, structure_id%in%tissue_id)
tst_onetissue2 <- subset(tst2, structure_id%in%tissue_id)
tst_onetissue3 <- subset(tst3, structure_id%in%tissue_id)
tst_onetissue4 <- subset(tst4, structure_id%in%tissue_id)
tst_onetissue5 <- subset(tst5, structure_id%in%tissue_id)
tst_onetissue6 <- subset(tst6, structure_id%in%tissue_id)

cast_tst1 <- data.frame(cast(tst_onetissue1, probe_id~sample_num))
cast_tst2 <- data.frame(cast(tst_onetissue2, probe_id~sample_num))
cast_tst3 <- data.frame(cast(tst_onetissue3, probe_id~sample_num))
cast_tst4 <- data.frame(cast(tst_onetissue4, probe_id~sample_num))
cast_tst5 <- data.frame(cast(tst_onetissue5, probe_id~sample_num))
cast_tst6 <- data.frame(cast(tst_onetissue6, probe_id~sample_num))

# Merge individual datasets
cast12 <- merge(cast_tst1, cast_tst2, by.x="probe_id", by.y="probe_id", all=T)
cast123 <- merge(cast12, cast_tst3, by.x="probe_id", by.y="probe_id", all=T)
cast1234 <- merge(cast123, cast_tst4, by.x="probe_id", by.y="probe_id", all=T)
cast12345 <- merge(cast1234, cast_tst5, by.x="probe_id", by.y="probe_id", all=T)
cast123456 <- merge(cast12345, cast_tst6, by.x="probe_id", by.y="probe_id", all=T)

# Exclude probesets that are not changing
rownames_cast123456 <- cast123456$probe_id
dim(cast123456)
cast123456 <-cast123456[, -1]

library(preprocessCore)
# Exclude probesets that are not changing
cast123456 <- normalize.quantiles(as.matrix(cast123456))
rownames(cast123456) <- rownames_cast123456
SD <- apply(cast123456, 1, sd, na.rm = T)
cast123456_filt <- cast123456[SD >= 0.29, ]
dim(cast123456_filt)

# Remove unnesessary data structures
rm(cast_tst1,cast_tst2, cast_tst3,
     cast_tst4, cast_tst5, cast_tst6,
     cast12, cast123, cast1234, cast12345, cast123456)

# Transpose matrix and compute coexpression and correlation p-value
m <- t(cast123456_filt)
# Write example data
file_name <- "CA2_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))

# CA3
tissue <- "CA3"
tissue_id <-c(4256, 4265)
# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)

# Subset tissue related data from each dataset
tst_onetissue1 <- subset(tst1, structure_id%in%tissue_id)
tst_onetissue2 <- subset(tst2, structure_id%in%tissue_id)
tst_onetissue3 <- subset(tst3, structure_id%in%tissue_id)
tst_onetissue4 <- subset(tst4, structure_id%in%tissue_id)
tst_onetissue5 <- subset(tst5, structure_id%in%tissue_id)
tst_onetissue6 <- subset(tst6, structure_id%in%tissue_id)

cast_tst1 <- data.frame(cast(tst_onetissue1, probe_id~sample_num))
cast_tst2 <- data.frame(cast(tst_onetissue2, probe_id~sample_num))
cast_tst3 <- data.frame(cast(tst_onetissue3, probe_id~sample_num))
cast_tst4 <- data.frame(cast(tst_onetissue4, probe_id~sample_num))
cast_tst5 <- data.frame(cast(tst_onetissue5, probe_id~sample_num))
cast_tst6 <- data.frame(cast(tst_onetissue6, probe_id~sample_num))

# Merge individual datasets
cast12 <- merge(cast_tst1, cast_tst2, by.x="probe_id", by.y="probe_id", all=T)
cast123 <- merge(cast12, cast_tst3, by.x="probe_id", by.y="probe_id", all=T)
cast1234 <- merge(cast123, cast_tst4, by.x="probe_id", by.y="probe_id", all=T)
cast12345 <- merge(cast1234, cast_tst5, by.x="probe_id", by.y="probe_id", all=T)
cast123456 <- merge(cast12345, cast_tst6, by.x="probe_id", by.y="probe_id", all=T)

# Exclude probesets that are not changing
rownames_cast123456 <- cast123456$probe_id
dim(cast123456)
cast123456 <-cast123456[, -1]

library(preprocessCore)
# Exclude probesets that are not changing
cast123456 <- normalize.quantiles(as.matrix(cast123456))
rownames(cast123456) <- rownames_cast123456
SD <- apply(cast123456, 1, sd, na.rm = T)
cast123456_filt <- cast123456[SD >= 0.29, ]
dim(cast123456_filt)

# Remove unnesessary data structures
rm(cast_tst1,cast_tst2, cast_tst3,
     cast_tst4, cast_tst5, cast_tst6,
     cast12, cast123, cast1234, cast12345, cast123456)

# Transpose matrix and compute coexpression and correlation p-value
m <- t(cast123456_filt)

# Write example data
file_name <- "CA3_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))

# CA4
tissue <- "CA4"
tissue_id <-c(4257, 4266)
# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)

# Subset tissue related data from each dataset
tst_onetissue1 <- subset(tst1, structure_id%in%tissue_id)
tst_onetissue2 <- subset(tst2, structure_id%in%tissue_id)
tst_onetissue3 <- subset(tst3, structure_id%in%tissue_id)
tst_onetissue4 <- subset(tst4, structure_id%in%tissue_id)
tst_onetissue5 <- subset(tst5, structure_id%in%tissue_id)
tst_onetissue6 <- subset(tst6, structure_id%in%tissue_id)

cast_tst1 <- data.frame(cast(tst_onetissue1, probe_id~sample_num))
cast_tst2 <- data.frame(cast(tst_onetissue2, probe_id~sample_num))
cast_tst3 <- data.frame(cast(tst_onetissue3, probe_id~sample_num))
cast_tst4 <- data.frame(cast(tst_onetissue4, probe_id~sample_num))
cast_tst5 <- data.frame(cast(tst_onetissue5, probe_id~sample_num))
cast_tst6 <- data.frame(cast(tst_onetissue6, probe_id~sample_num))

# Merge individual datasets
cast12 <- merge(cast_tst1, cast_tst2, by.x="probe_id", by.y="probe_id", all=T)
cast123 <- merge(cast12, cast_tst3, by.x="probe_id", by.y="probe_id", all=T)
cast1234 <- merge(cast123, cast_tst4, by.x="probe_id", by.y="probe_id", all=T)
cast12345 <- merge(cast1234, cast_tst5, by.x="probe_id", by.y="probe_id", all=T)
cast123456 <- merge(cast12345, cast_tst6, by.x="probe_id", by.y="probe_id", all=T)

# Exclude probesets that are not changing
rownames_cast123456 <- cast123456$probe_id
dim(cast123456)
cast123456 <-cast123456[, -1]

library(preprocessCore)
# Exclude probesets that are not changing
cast123456 <- normalize.quantiles(as.matrix(cast123456))
rownames(cast123456) <- rownames_cast123456
SD <- apply(cast123456, 1, sd, na.rm = T)
cast123456_filt <- cast123456[SD >= 0.29, ]
dim(cast123456_filt)

# Remove unnesessary data structures
rm(cast_tst1,cast_tst2, cast_tst3,
     cast_tst4, cast_tst5, cast_tst6,
     cast12, cast123, cast1234, cast12345, cast123456)

# Transpose matrix and compute coexpression and correlation p-value
m <- t(cast123456_filt)

# Write example data
file_name <- "CA4_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))


# subiculum
tissue <- "subiculum"
tissue_id <-c(4251,4260)
# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)

# Subset tissue related data from each dataset
tst_onetissue1 <- subset(tst1, structure_id%in%tissue_id)
tst_onetissue2 <- subset(tst2, structure_id%in%tissue_id)
tst_onetissue3 <- subset(tst3, structure_id%in%tissue_id)
tst_onetissue4 <- subset(tst4, structure_id%in%tissue_id)
tst_onetissue5 <- subset(tst5, structure_id%in%tissue_id)
tst_onetissue6 <- subset(tst6, structure_id%in%tissue_id)

cast_tst1 <- data.frame(cast(tst_onetissue1, probe_id~sample_num))
cast_tst2 <- data.frame(cast(tst_onetissue2, probe_id~sample_num))
cast_tst3 <- data.frame(cast(tst_onetissue3, probe_id~sample_num))
cast_tst4 <- data.frame(cast(tst_onetissue4, probe_id~sample_num))
cast_tst5 <- data.frame(cast(tst_onetissue5, probe_id~sample_num))
cast_tst6 <- data.frame(cast(tst_onetissue6, probe_id~sample_num))

# Merge individual datasets
cast12 <- merge(cast_tst1, cast_tst2, by.x="probe_id", by.y="probe_id", all=T)
cast123 <- merge(cast12, cast_tst3, by.x="probe_id", by.y="probe_id", all=T)
cast1234 <- merge(cast123, cast_tst4, by.x="probe_id", by.y="probe_id", all=T)
cast12345 <- merge(cast1234, cast_tst5, by.x="probe_id", by.y="probe_id", all=T)
cast123456 <- merge(cast12345, cast_tst6, by.x="probe_id", by.y="probe_id", all=T)

# Exclude probesets that are not changing
rownames_cast123456 <- cast123456$probe_id
dim(cast123456)
cast123456 <-cast123456[, -1]

library(preprocessCore)
# Exclude probesets that are not changing
cast123456 <- normalize.quantiles(as.matrix(cast123456))
rownames(cast123456) <- rownames_cast123456
SD <- apply(cast123456, 1, sd, na.rm = T)
cast123456_filt <- cast123456[SD >= 0.29, ]
dim(cast123456_filt)

# Remove unnesessary data structures
rm(cast_tst1,cast_tst2, cast_tst3,
     cast_tst4, cast_tst5, cast_tst6,
     cast12, cast123, cast1234, cast12345, cast123456)

# Transpose matrix and compute coexpression and correlation p-value
m <- t(cast123456_filt)

# Write example data
file_name <- "subiculum_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))


#SptN
tissue <- "SptN"
tissue_id <-c(4301, 4304)
# Create folder for storeing the results
pathRdataTissue <- sprintf("%s", tissue);
dir.create(file.path(pathRdata,pathRdataTissue),showWarnings = FALSE, recursive = TRUE)

# Subset tissue related data from each dataset
tst_onetissue1 <- subset(tst1, structure_id%in%tissue_id)
tst_onetissue2 <- subset(tst2, structure_id%in%tissue_id)
tst_onetissue3 <- subset(tst3, structure_id%in%tissue_id)
tst_onetissue4 <- subset(tst4, structure_id%in%tissue_id)
tst_onetissue5 <- subset(tst5, structure_id%in%tissue_id)
tst_onetissue6 <- subset(tst6, structure_id%in%tissue_id)

cast_tst1 <- data.frame(cast(tst_onetissue1, probe_id~sample_num))
cast_tst2 <- data.frame(cast(tst_onetissue2, probe_id~sample_num))
cast_tst3 <- data.frame(cast(tst_onetissue3, probe_id~sample_num))
cast_tst4 <- data.frame(cast(tst_onetissue4, probe_id~sample_num))
cast_tst5 <- data.frame(cast(tst_onetissue5, probe_id~sample_num))
cast_tst6 <- data.frame(cast(tst_onetissue6, probe_id~sample_num))

# Merge individual datasets
cast12 <- merge(cast_tst1, cast_tst2, by.x="probe_id", by.y="probe_id", all=T)
cast123 <- merge(cast12, cast_tst3, by.x="probe_id", by.y="probe_id", all=T)
cast1234 <- merge(cast123, cast_tst4, by.x="probe_id", by.y="probe_id", all=T)
cast12345 <- merge(cast1234, cast_tst5, by.x="probe_id", by.y="probe_id", all=T)
cast123456 <- merge(cast12345, cast_tst6, by.x="probe_id", by.y="probe_id", all=T)

# Exclude probesets that are not changing
rownames_cast123456 <- cast123456$probe_id
dim(cast123456)
cast123456 <-cast123456[, -1]

library(preprocessCore)
# Exclude probesets that are not changing
cast123456 <- normalize.quantiles(as.matrix(cast123456))
rownames(cast123456) <- rownames_cast123456
SD <- apply(cast123456, 1, sd, na.rm = T)
cast123456_filt <- cast123456[SD >= 0.29, ]
dim(cast123456_filt)

# Remove unnesessary data structures
rm(cast_tst1,cast_tst2, cast_tst3,
     cast_tst4, cast_tst5, cast_tst6,
     cast12, cast123, cast1234, cast12345, cast123456)

# Transpose matrix and compute coexpression and correlation p-value
m <- t(cast123456_filt)

# Write example data
file_name <- "SptN_preprocessed.RData"
save(m,file=file.path(pathRdata,pathRdataTissue,file_name))



  



 
