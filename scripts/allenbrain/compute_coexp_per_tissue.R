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


#t11 <- unique(s1$structure_id)
#t22 <- unique(s2$structure_id)
#t33 <- unique(s3$structure_id)
#t44 <- unique(s4$structure_id)
#t55 <- unique(s5$structure_id)
#t66 <- unique(s6$structure_id)


# List of unique tissues in all 6 data sets
#tissues_all <- unique(c(t11,t22,t33,t44,t55,t66))

# Number of tissues
#length(tissues_all) #414

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


##################################################### Below example for one tissue ########################
# Compute coexpression in each of the tissues
#foreach(i = 1:length(selected_tissues))%do%{

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
tst_onetissue1 <- subset(tst1, structure_id==tissue)
tst_onetissue2 <- subset(tst2, structure_id==tissue)
tst_onetissue3 <- subset(tst3, structure_id==tissue)
tst_onetissue4 <- subset(tst4, structure_id==tissue)
tst_onetissue5 <- subset(tst5, structure_id==tissue)
tst_onetissue6 <- subset(tst6, structure_id==tissue)

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
save(m,file="m_tissue.RData")

test<-m[,1:20]

  
#library(forward)
#  correl.stats=function(X, method = "spearman", use = "pairwise.complete.obs"){
#    require(forward)
#    combs=t(fwd.combn(colnames(X), 2))
#    temp=t(apply(combs,1, function(x){
#      Y=X[,as.character(x)]
#      res=cor.test(Y[,1],Y[,2], use = use, method = method)
#      temp2=c(res$estimate, res$p.value)
#      names(temp2)=c('rho','pvalue')
#      rm(res)
#      temp2
#    }
#    )
#    )
#    temp <- data.frame(cbind(gene.A=combs[,1],gene.B=combs[,2],temp), stringsAsFactors = F)
#    temp$rho <- as.numeric(temp$rho)
#    temp$pvalue <- as.numeric(temp$pvalue)
#    rownames(temp)=paste(combs[,1],combs[,2],sep="")
#    temp <- temp[temp$pvalue<=0.05,]
#    temp
#  }
  
# corel_test_res_one_tissue <- correl.stats(test)

#correl.stats2=function(X, method = "spearman", use = "pairwise.complete.obs"){
#    require(forward)
#    combs <- t(fwd.combn(colnames(X), 2))
#    temp <- t(apply(combs,1, function(x){
#      Y <- X[,as.character(x)]
#      res <- cor.test(Y[,1],Y[,2], use = "pairwise.complete.obs", method = "spearman")
#      temp2=c(colnames(Y)[1], colnames(Y)[2],res$estimate, res$p.value)
#      #temp2 <- c(res$estimate, res$p.value)
#      names(temp2) <- c("gene.A", "gene.B","rho","pvalue")
#      temp2 <- as.data.frame(t(as.matrix(temp2)), stringsAsFactors = F)
#    #  temp2$pvalue<- as.numeric(temp2$pvalue)
#    #  temp2$rho<- as.numeric(temp2$rho)
#      temp2 <- temp2[temp2$pvalue<=0.05,]
#      rm(res)
#      temp2_new <- as.character(temp2[1,])
#      names(temp2_new) <- colnames(temp2)
#     #   temp2_new <- na.omit(temp2_new)
#        temp2_new
#        }
#    )
#    )
#   # temp <- data.frame(cbind(gene.A=combs[,1],gene.B=combs[,2],temp), stringsAsFactors = F)
#   # temp$rho <- as.numeric(temp$rho)
#   # temp$pvalue <- as.numeric(temp$pvalue)
#   # rownames(temp)=paste(combs[,1],combs[,2],sep="")
#   # temp <- temp[temp$pvalue<=0.05,]
#
#   temp=temp[complete.cases(temp),]
#    temp
#  }
# system.time(corel_test_res_one_tissue <- correl.stats2(test))
# corel_test_res_one_tissue 

#system.time(corel_test_res_one_tissue_m <- correl.stats2(m))
# head(corel_test_res_one_tissue_m)
#dim(corel_test_res_one_tissue_m)




#### Append to the file

# tissue file where the data will be saved
#filenametissue<- sprintf("%s.txt", format(tissue, scientific=F))
#sink_filenametissue<-sprintf("%s.txt", "sink_test")
#correl.stats3=function(X, method = "spearman", use = "pairwise.complete.obs"){
#  require(forward)
#  combs <- t(fwd.combn(colnames(X), 2))
#    apply(combs,1, function(x){
#    Y <- X[,as.character(x)]
#    res <- cor.test(Y[,1],Y[,2], use = "pairwise.complete.obs", method = "spearman")
#    temp2=c(colnames(Y)[1], colnames(Y)[2],res$estimate, 
#res$p.value)
#    names(temp2) <- c("gene.A", "gene.B","rho","pvalue")
#    temp2 <- as.data.frame(t(as.matrix(temp2)), stringsAsFactors = F)
#    temp2 <- temp2[temp2$pvalue<=0.05,]
#    write.table(temp2, file=file.path(pathRdata,pathRdataTissue,filenametissue),
#                append = T, quote=F,sep="\t", row.names=F, col.names=F)
#     sink(file.path(pathRdata,pathRdataTissue,sink_filenametissue), append=T)
#    print(temp2,right=F)
#    sink()
#    rm(res)
#    temp2_new <- as.character(temp2[1,])
#    names(temp2_new) <- colnames(temp2)
#  }
#  )
#}

#system.time(correl.stats3(test))
#system.time(correl.stats3(m))






##### Try with foreach per one probeset
filenametissue<- sprintf("%s.txt", format(tissue, scientific=F))

#filenametissue<- sprintf("%s.txt", "foreach_per_probeset")
sink_filenametissue<-sprintf("%s.txt", "foreach_sink_per_probeset")

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
   write.table(cor_one_probeset, file=file.path(pathRdata,pathRdataTissue,filenametissue), append = T, quote=F,sep="\t", row.names=F, col.names=F)
    # Sink output to the file
    sink(file.path(pathRdata,pathRdataTissue,sink_filenametissue), append=T)
    print(cor_one_probeset,right=F)
    sink()

}
)



tmp <- read.table(file=file.path(pathRdata,pathRdataTissue,filenametissue), header=F, sep = "\t")  
  
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
filenametissue_filt<- sprintf("%s.txt", paste(format(tissue, scientific=F),"filt", sep="_"))
write.table(filtered_coexpression, file=file.path(pathRdata,pathRdataTissue,filenametissue_filt), append = F, quote=F,sep="\t", row.names=F, col.names=F)









# Remove unnesessary data structures 
rm(cast_tst1,cast_tst2, cast_tst3,
     cast_tst4, cast_tst5, cast_tst6, 
     castcast12, cast123, cast1234, cast12345, cast123456)

# Assign m to a dataframe
df <- m

# Since the length of all data matrix in one tissue is too 
# For the sake of computation it is split into 3 parts
# Determimne the whole number of probesets of matrix m
n <- length(colnames(m))

# Define the starting points of intervals for computaion
n1 <- round(n/3)
n2 <- round(n/3*2)


# Compute cor.test
# Define function to compute cor p-value
corpij <- function(i,j,data) {cor.test(data[,i],data[,j],  method = c("spearman"))$p.value}
corp <- Vectorize(corpij, vectorize.args=list("i","j"))

# Compute the fist part of all against all correlation p-values in m
corp_m1 <- outer(1:n,1:n1,corp,data=df)
colnames(corp_m1) <- colnames(df)[1:n1]
rownames(corp_m1) <- colnames(df)[1:n]
corp_m1 <- melt(corp_m1)
colnames(corp_m1) <- c("probe1","probe2","pval")
corp_m1$pval <- p.adjust(corp_m1$pval, method ="BH", n = length(corp_m1$pval))
corp_m1$probe1 <- as.character(corp_m1$probe1)
corp_m1$probe2 <- as.character(corp_m1$probe2)
corp_m1 <- corp_m1[!corp_m1$probe1==corp_m1$probe2,]

# Compute the second part of all against all correlation p-values in m
corp_m2 <- outer(1:n,(n1+1):n2,corp,data=df)
colnames(corp_m2) <- colnames(df)[(n1+1):n2]
rownames(corp_m2) <- colnames(df)[1:n]
corp_m2 <- melt(corp_m2)
colnames(corp_m2) <- c("probe1","probe2","pval")
corp_m2$pval <- p.adjust(corp_m2$pval, method ="BH", n = length(corp_m2$pval))
corp_m2$probe1 <- as.character(corp_m2$probe1)
corp_m2$probe2 <- as.character(corp_m2$probe2)
corp_m2 <- corp_m2[!corp_m2$probe1==corp_m2$probe2,]

# Compute the third part of all against all correlation p-values in m
corp_m3 <- outer(1:n,(n2+1):n,corp,data=df)
colnames(corp_m3) <- colnames(df)[(n2+1):n]
rownames(corp_m3) <- colnames(df)[1:n]
corp_m3 <- melt(corp_m3)
colnames(corp_m3) <- c("probe1","probe2","pval")
corp_m3$pval <- p.adjust(corp_m3$pval, method ="BH", n = length(corp_m3$pval))
corp_m3$probe1 <- as.character(corp_m3$probe1)
corp_m3$probe2 <- as.character(corp_m3$probe2)
corp_m3 <- corp_m3[!corp_m3$probe1==corp_m3$probe2,]



# Compute correlation
# Define correlation function
corij <- function(i,j,data) {cor(data[,i],data[,j],method = "spearman", use="pairwise.complete.obs")}
corm <- Vectorize(corij, vectorize.args=list("i","j"))

# Compute the fist part of all against all correlation in m
cor_m1 <- outer(1:n,1:n1, corm, data=df)
colnames(cor_m1) <- colnames(df)[1:n1]
rownames(cor_m1) <- colnames(df)[1:n]
cor_m1 <- melt(cor_m1)
colnames(cor_m1) <- c("probe1","probe2","cor")
cor_m1$probe1 <- as.character(cor_m1$probe1)
cor_m1$probe2 <- as.character(cor_m1$probe2)
cor_m1 <- cor_m1[!cor_m1$probe1==cor_m1$probe2,]

# Compute the second part of all against all correlation in m
cor_m2 <- outer(1:n,(n1+1):n2, corm, data=df)
colnames(cor_m2) <- colnames(df)[(n1+1):n2]
rownames(cor_m2) <- colnames(df)[1:n]
cor_m2 <- melt(cor_m2)
colnames(cor_m2) <- c("probe1","probe2","cor")
cor_m2$probe1 <- as.character(cor_m2$probe1)
cor_m2$probe2 <- as.character(cor_m2$probe2)
cor_m2 <- cor_m2[!cor_m2$probe1==cor_m2$probe2,]

# Compute the fird part of all against all correlation in m
cor_m3 <- outer(1:n,(n2+1):n, corm, data=df)
colnames(cor_m3) <- colnames(df)[(n2+1):n]
rownames(cor_m3) <- colnames(df)[1:n]
cor_m3 <- melt(cor_m3)
colnames(cor_m3) <- c("probe1","probe2","cor")
cor_m3$probe1 <- as.character(cor_m3$probe1)
cor_m3$probe2 <- as.character(cor_m3$probe2)
cor_m3 <- cor_m3[!cor_m3$probe1==cor_m3$probe2,]

# Display numbers in the format 10^
options("scipen"=-100)

# Combine results of cor.test and cor
cor_one_tissue1 <- cbind(cor_m1,pval=corp_m1$pval)
cor_one_tissue1$pval <- as.numeric(as.character(cor_one_tissue1$pval))
cor_one_tissue1_sig <- cor_one_tissue1[cor_one_tissue1$pval<=0.05,]
cor_one_tissue2 <- cbind(cor_m2,pval=corp_m2$pval)
cor_one_tissue2$pval <- as.numeric(as.character(cor_one_tissue2$pval))
cor_one_tissue2_sig <- cor_one_tissue2[cor_one_tissue2$pval<=0.05,]
cor_one_tissue3 <- cbind(cor_m3,pval=corp_m3$pval)
cor_one_tissue3$pval <- as.numeric(as.character(cor_one_tissue3$pval))
cor_one_tissue3_sig <- cor_one_tissue3[cor_one_tissue3$pval<=0.05,]

# Save intemediate results
# Create file name based on the tissue name
filename<- sprintf("%s.RData", paste(format(tissue, scientific=F), "1", sep="_"))
# Combine filepath and filename
pathdata <- file.path(pathRdata,pathRdataTissue,filename)
# Save I part of correlation matrix based on the 
save(cor_one_tissue1, file = pathdata)

filename<- sprintf("%s.RData", paste(format(tissue, scientific=F), "2", sep="_"))
# Combine filepath and filename
pathdata <- file.path(pathRdata,pathRdataTissue,filename)
# Save II part of correlation matrix based on the 
save(cor_one_tissue2, file = pathdata)

filename<- sprintf("%s.RData", paste(format(tissue, scientific=F), "3", sep="_"))
# Combine filepath and filename
pathdata <- file.path(pathRdata,pathRdataTissue,filename)
 # Save III part correlation matrix based on the 
save(cor_one_tissue3, file = pathdata)


# Combine all significant correlations into one datasfame 
signif_cor_tissue <- rbind(cor_one_tissue1_sig,cor_one_tissue2_sig,cor_one_tissue3_sig )
# Remove the duplicated undirrescted edges with the same score.
# For example probeset1-probeset2 0.5 and probeset2-probeset1 0.5
# Convert factors to characters
df2string<-function(df){
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  df[,3]<-as.numeric(df[,3])
  return (df)}

signif_cor_tissue=df2string(signif_cor_tissue)
signif_cor_tissue <- signif_cor_tissue[!duplicated(signif_cor_tissue), ]
signif_cor_tissue <- signif_cor_tissue[!duplicated(data.frame(t(apply(signif_cor_tissue[1:2], 1, sort)), signif_cor_tissue$cor)),]



# Create file name based on the tissue name
filename<- sprintf("%s.RData", format(tissue, scientific=F));
# Combine filepath and filename
pathdata <- file.path(pathRdata,pathRdataTissue,filename);
# Save significant correlation matrix based on the 
save(signif_cor_tissue, file = pathdata)


#}







 
