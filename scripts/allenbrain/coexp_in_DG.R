# This script computes co-expression in DG brain tissue from Allen Brain Atlas.
# Co-expression is computed over 6 datasets of gene expression in healthy brains.

# Set working directory
setwd("~/absb/results/allenbrain/tissues_rdata/DG")
tissue_path <- "~/absb/results/allenbrain/tissues_rdata/DG"

# Used libaries
library(reshape)
library(foreach)
        

# Load DG dataset
load(file="DG_preprocessed.RData")

##### Try with foreach per one probeset
tissue <- "DG"
filenametissue<- sprintf("%s.txt",tissue)

#filenametissue<- sprintf("%s.txt", "foreach_per_probeset")
sink_filenametissue<-sprintf("%s.txt", "foreach_sink_per_probeset_DG")

probesets <- colnames(m)
library("psych")

library(doParallel)
registerDoParallel(cores=4)
system.time(foreach(i = 1:length(probesets)) %dopar%{
  #system.time(foreach(i = 1:3) %dopar%{
  cor1gds <- c() # Corelation for one gene in one ds
  probeset <- probesets[i]
  print("probeset")
  print(probeset)
  print(i)
  vect <- m[,i]
  mtrx <- m[,-c(i)]
  cor1gds <- t(cor(vect, mtrx, method = "spearman", use = "pairwise.complete.obs"))
  # Compute p-values for correlation matrix
  cor_pval <- corr.p(cor1gds,length(vect),adjust="fdr",alpha=.05)$p
  # Combine probesets' names, corresponding correlations and p-values.
  cor_one_probeset <- cbind(rep(paste(probeset),length(rownames(cor1gds))),rownames(cor1gds),cor1gds,cor_pval)
  colnames(cor_one_probeset) <- c("probeset.A", "probeset.B", "rho", "pvalue")
  # Filter our rows with p-values > 0.05
  cor_one_probeset <- data.frame(cor_one_probeset, stringsAsFactors = F)
  cor_one_probeset$rho <- as.numeric(cor_one_probeset$rho)
  cor_one_probeset$pvalue <- as.numeric(cor_one_probeset$pvalue)
  cor_one_probeset <- cor_one_probeset[cor_one_probeset$pvalue<=0.05,]
  print("dim cor_one_probeset")
  print(dim(cor_one_probeset))
  print("head")
  print(head(cor_one_probeset))

# Append write to txt file
   write.table(cor_one_probeset, file=file.path(tissue_path,filenametissue), append = T, quote=F,sep="\t", row.names=F, col.names=F)
    # Sink output to the file
    sink(file.path(tissue_path,sink_filenametissue), append=T)
    print(cor_one_probeset,right=F)
    sink()
}
)

# Read the computed correlations and remove duplicated pares AB BA
tmp <- read.table(file=file.path(tissue_path,filenametissue), header=F, sep = "\t")  
  
  # Convert factors to characters
  df2string<-function(df){
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    df[,3]<-as.numeric(df[,3])
    return (df)}

# Remove duplicated values like AB BA  
  tmp <- df2string(tmp)
  str(tmp)
  tmp<- tmp[!duplicated(tmp), ]
  dim(tmp)
  tmp <- tmp[!duplicated(data.frame(t(apply(tmp[1:2], 1, sort)), tmp$V3)),]
 # New size
  dim(tmp)
 filtered_coexpression <- tmp
 colnames(filtered_coexpression) <- c("probeset.A", "probeset.B", "cor", "cor_pval")


# Save to file
filenametissue_filt<- sprintf("%s.txt", paste(tissue,"filt", sep="_"))
write.table(filtered_coexpression, file=file.path(tissue_path,filenametissue_filt), append = F, quote=F,sep="\t", row.names=F, col.names=F)



# Load converted probesets IDs
load(file = "~/absb/results/allenbrain/genes_aba_ensg.RData")
# Add ensg IDs for DG co-expression dataset
DG_filt <- filtered_coexpression
DG_filt_ensg1 <- merge(DG_filt, genes_aba_ensg, by.x = "probeset.A", by.y = "probe_id", all = F)
DG_filt_ensg12 <- merge(DG_filt_ensg1, genes_aba_ensg, by.x = "probeset.B", by.y = "probe_id", all = F)

# Select only ensg1 ensg2 score and add type of interaction and data source
library(plyr)
colnames(DG_filt_ensg12)[c(6,8)] <- c("ensg1", "ensg2")
DG <- aggregate(cor ~ ensg1 + ensg2, data = DG_filt_ensg12, max)
DG <- DG[, c("ensg1", "ensg2", "cor")]
DG <- cbind(DG, interaction_type = "coexpression")
# Add data source name Allen Brain Atlas (ABA)
DG <- cbind(DG, data_source = "ABA")

# Rename the columns
colnames(DG) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

## Remove the duplicated undirrescted edges with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5

# Convert factors to characters
df2string<-function(df){
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  df[,3]<-as.numeric(df[,3])
  return (df)}

# Co-expression in DG
DG <- df2string(DG)
str(DG)
DG <- DG[!duplicated(DG), ]
dim(DG)
DG <- DG[!duplicated(data.frame(t(apply(DG[1:2], 1, sort)), DG$score)),]

# Size after self-loops and duplicated interactions were removed
dim(DG)

head(DG)
DG_coexp_int <- DG
DG_coexp_int$data_source <- "ABA_DG"

save(DG_coexp_int, file = "~/absb/results/allenbrain/DG_coexp_int.RData")
write.table(DG_coexp_int,file = "~/absb/results/allenbrain/DG_coexp_int.txt", quote = F, sep = "\t",row.names = F)








 
