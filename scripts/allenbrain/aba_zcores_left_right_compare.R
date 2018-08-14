# This script compares z-scores in the same brain region from right and left hemispheares i.e CA1 region with tissue ID 4254(L) and 4263(R). 

# Set working directory
setwd("~/absb/results/allenbrain/")

# Used libraries
library(reshape2)
library(foreach)

# Path to the results saved in RData format for individual tissues
pathRdata <- "~/absb/results/allenbrain/tissues_rdata/"

# Load the list of tissues
load(file = "~/absb/results/allenbrain/all_tissues.RData" )
tissues_all <-sort(tissues_all)


# Select tissue IDs in left and right hemisphears for CA1, CA2, CA3, CA4, DG, SptN, subiculum
CA1_left_right <-c("4254","4263")

# load matrix with z-scores in left and right hemisphears called tissues_zscores_mtx
load(file = "tissues_zscores_mtx.RData")
# Test on CA1 resgion if right and left expression not differ based on z-scores
