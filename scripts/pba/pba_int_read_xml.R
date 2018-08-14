# This script processes dataset of protein-protein inteacrtions related to brain ageing (PBA) 

## Create the folder where current results will be written
resdir <- paste("~/absb/results", "pba", sep = "/")
dir.create(file.path(resdir), showWarnings  =  FALSE,  recursive  =  TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file in xml format and save the interactions as RData
# Please indicate the path to the downloaded data
library(XML)
hybxml <- xmlParse("~/absb/data/pba/agedbrain_psimi25_data.xml")

# Convert xml to list
xml_data <- xmlToList(hybxml)

# Length of the data
length(xml_data$entry$interactionList)

# Run example for one interaction in the data
xml_data$entry$interactionList[[1]]$participantList$participant$interactorRef ##1st interactor ID
#[1] "1"
#xml_data$entry$interactionList[[1]]$confidenceList$confidence$score # interaction score
xml_data$entry$interactionList[[1]]$confidenceList$confidence$value # interaction value defined by Hybrigenics A,B,C,D
xml_data$entry$interactionList[[1]]$participantList[[2]]$interactorRef #2nd interactor ID
xml_data$entry$experimentList$experimentDescription$bibref$xref[[1]][3] #primary ref
xml_data$entry$experimentList$experimentDescription$bibref$xref[[2]][3] #secondary ref
xml_data$entry$experimentList$experimentDescription$bibref$xref[[3]][3] #secondary ref
xml_data$entry$interactionList[[1]]$interactionType$names$shortLabel #interaction type(physical association)
xml_data$entry$interactorList[[1]]$organism$names$shortLabel #organizm

## Get the part of the data related to the interactors (IDs, corresponding name, organism, reference ID)

# Number of all interactors
length(xml_data$entry$interactorList) 

# Toy example for one interactor
# Short name
xml_data$entry$interactorList[1]$interactor$names$shortLabel 
# Long name
xml_data$entry$interactorList[1]$interactor$names$fullName 
# ENTREZ gene Identifier
xml_data$entry$interactorList[1]$interactor$xref$secondaryRef[[3]]
# Organism ncbi taxomony ID
xml_data$entry$interactorList[1]$interactor$organism$.attrs[[1]]
# Organism
xml_data$entry$interactorList[1]$interactor$organism$names$fullName
# Interactor;s ID
xml_data$entry$interactorList[1]$interactor$.attrs[[1]]

## Create a data frame that contains info about the interactors

# Allocate ampty vectors
name.h <- c()
fullname.h <- c()
intid.h <- c()
organism.h <- c()
entrez.h <-c()

for (i in 1:length(xml_data$entry$interactorList)){
  name.h <- c(name.h, xml_data$entry$interactorList[i]$interactor$names$shortLabel)
  fullname.h <- c(fullname.h, xml_data$entry$interactorList[i]$interactor$names$fullName)
  intid.h <- c(intid.h, xml_data$entry$interactorList[i]$interactor$.attrs[[1]])
  organism.h <- c(organism.h, xml_data$entry$interactorList[i]$interactor$organism$names$fullName) 
  entrez.h <-c(entrez.h, xml_data$entry$interactorList[i]$interactor$xref$secondaryRef[[3]])
}

# Combine individual vectors about interactors into dataframe
interactors <- cbind(intid.h, entrez.h, name.h, fullname.h, organism.h)
interactors <- as.data.frame(interactors)
interactors[,1] <- as.character(interactors[,1])
interactors[,2] <- as.character(interactors[,2])
interactors[,3] <- as.character(interactors[,3])
interactors[,4] <- as.character(interactors[,4])
interactors[,5] <- as.character(interactors[,5])

# Create data frame about interactions
# Allocate empty vectors
intid1.h <- c()
intid2.h <- c()
confval.h <- c()

# Extract needed information
for (i in 1:length(xml_data$entry$interactionList)){
  intid1.h <- c(intid1.h, xml_data$entry$interactionList[[i]]$participantList$participant$interactorRef)
  intid2.h <- c(intid2.h, xml_data$entry$interactionList[[i]]$participantList[[2]]$interactorRef)
  confval.h <- c(confval.h, xml_data$entry$interactionList[[i]]$confidenceList$confidence$value)
}

# Combine firts interactors IS , second interactor's ID ans interaction score into dataframe  
interactions <- cbind(intid1.h, intid2.h, confval.h)
interactions <- as.data.frame(interactions)
interactions[,1] <- as.character(interactions[,1])
interactions[,2] <- as.character(interactions[,2])

#### Merge interaction data with interactors data
pba=merge(interactions,interactors,by.x="intid1.h",by.y="intid.h", all=F)
pbafinal=merge(pba,interactors,by.x="intid2.h",by.y="intid.h", all=F)
pba_ppi <-pbafinal[,c(4,8,6,10,3,7,2,1)]

# Rename the columns
colnames(pba_ppi)=c("entrez.p1","entrez.p2","name.p1","name.p2","confidence.value","organism","ref.id.1","ref.id.2")


# Save as R data 
save(pba_ppi, file =  "pba_ppi_original.RData")

# Assemble the file in the required format ensg1 ensg2 score interaction_type data_source
# Load the data if needed
# load(file <- "pba_ppi_original.RData")
# See how the data looks
head(pba_ppi)

# Control the dimentions of the data
dim(pba_ppi)

# Split df by organism

# Print out the organisms present in the data
unique(pba_ppi$organism)
#[1] "Human"                   "Mouse"                  
#[3] "Homo sapiens"            "Mus musculus"           
#[5] "Drosophila melanogaster" "Drosophila"

# Substitute "Human with "Homo Sapiens"
pba_ppi$organism <- gsub("Human", "Homo sapiens", pba_ppi$organism)

# Substitute "Mouse" with "Mus Musculus"
pba_ppi$organism <- gsub("Mouse", "Mus musculus", pba_ppi$organism)

# Substitute "Drosophila" with "Drosophila melanogaster"
pba_ppi$organism <- gsub("Drosophila", "Drosophila melanogaster", pba_ppi$organism)
pba_ppi$organism <- gsub("Drosophila melanogaster melanogaster", "Drosophila melanogaster", pba_ppi$organism)

# Extract part of the data related to Human
pba_ppi.hs <- pba_ppi[pba_ppi$organism%in%"Homo sapiens",]

# Save human data separately
pba_ppi.hs<-pba_ppi.hs[,1:5]
save(pba_ppi.hs, file="PBA_PPI_HS.RData")
write.table(pba_ppi.hs, file = "~/absb/data/pba/PBA_PPI_HS.txt", sep = "\t", quote = F, row.names = F)


# Extract part of the data related to Mouse
pba_ppi.mm <- pba_ppi[pba_ppi$organism%in%"Mus musculus",]

# Extract part of the data related to Drosophila
pba_ppi.dm <- pba_ppi[pba_ppi$organism%in%"Drosophila melanogaster",]

# Control the size of the data related to each of the organisms
dim(pba_ppi.hs)
dim(pba_ppi.mm)
dim(pba_ppi.dm)



# Read in the data
pba_ppi.hs <- read.table(file = "~/absb/data/pba/PBA_PPI_HS.txt", sep = "\t", header = T, stringsAsFactors = F) 



# Convert pba_ppi.hs protein names to ENSG and bing them to the dataframe.
length(unique(c(pba_ppi.hs[,1], pba_ppi.hs[,2])))#1250



