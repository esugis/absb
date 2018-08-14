# This script combines all indiviadual interaction datasets into one and removes duplicated interaction of one type

# Create the folder where current results will be written
resdir<-paste("~/absb/results","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Load individual interaction datasets

# PPIs associated with brain ageing (PBA)
load(file="~/absb/results/pba/pba_int.RData")

# Combined epistatic interactions ADNI ventrical volume, ADNI cognitive traits, TGEN, HBTRC
load(file="~/absb/results/epistasis/epistasis_all_int.RData")

# Intact hyman PPIs with MIscore >=0.45
load(file="~/absb/results/intact/intact_int.RData")

## IntAct Alzheimer's related manually curated PPI dataset with MIscore >=0.45
load(file="~/absb/results/intact/alz_intact_int.RData")

## IntAct Alzheimer's related manually curated PPI dataset with all MIscores
#load(file="~/absb/results/intact/alz_intact_int_no_score_filtering.RData")

# IntAct Synapse related automatically curated PPI dataset with MIscore >=0.45
load(file="~/absb/results/intact/syn_intact_int.RData")

# Co-expression dataset with removed self-loops RRA score <= 0.00001 
load(file="~/absb/results/adn/integration/adn_coexp_int.RData")

# Co-expression in DG
load(file="~/absb/results/allenbrain/DG_coexp_int.RData")

# Co-expression	in CA1
load(file="~/absb/results/allenbrain/CA1_coexp_int.RData")   

# Co-expression in CA2
load(file="~/absb/results/allenbrain/CA2_coexp_int.RData")

# Co-expression in CA3
load(file="~/absb/results/allenbrain/CA3_coexp_int.RData")

# Co-expression in CA4
load(file="~/absb/results/allenbrain/CA4_coexp_int.RData")

# Co-expression in subiculum (S)
load(file="~/absb/results/allenbrain/subiculum_coexp_int.RData")                 

# Co-expression in SptN
load(file="~/absb/results/allenbrain/SptN_coexp_int.RData")

# Create one DF from separate datasets
integrated_int <- rbind(pba_int, epistasis_all_int, intact_int, alz_intact_int, syn_intact_int, adn_coexp_int, DG_coexp_int, CA1_coexp_int, CA2_coexp_int, CA3_coexp_int, CA4_coexp_int, subiculum_coexp_int, SptN_coexp_int)
integrated_int <- integrated_int[!duplicated(integrated_int),]
colnames(integrated_int)[1:2]<-c("ensg.A","ensg.B")

# Size of integrated dataset
dim(integrated_int) #307036837     5

# Sort rows ensg.A ensg.B in alpahbetic order
#integrated_int[,1] <- as.character(integrated_int[,1])
#integrated_int[,2] <- as.character(integrated_int[,2])

# Sort ensg in ensg.A and ensg.B row wise
#ensgAB <- data.frame(t(apply(integrated_int[,1:2], 1, sort)))

# Substitute columns for the sorted 
#integrated_int[,1:2] <- ensgAB  

# Save sorted integrated dataset to the files.
save(integrated_int, file = "integrated_int.RData")
write.table(integrated_int, file = "integrated_int.txt", sep="\t", quote=F, row.names=F)

dim(integrated_int)
