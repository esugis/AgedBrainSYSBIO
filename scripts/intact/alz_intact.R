# This script processes inteacrtions from curated datasets related to Alzheimer's disease from IntAct with mi score >0.45

# Create the folder where current results will be written
resdir<-paste("~/AgedBrainSYSBIO/results,"intact",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file of all interactions with miscore >0.45 downloaded from IntAct and save them as RData
# Please indicate the path to the downloaded data
alz_intact=read.delim("~/AgedBrainSYSBIO/data/intact/alzheimers_intact_v_4_2_6.txt", stringsAsFactors=F, sep = "\t", header=T)
save(alz_intact, file = "alz_intact.RData")

# Get the column names
colnames(alz_intact)

# Select only interactions related to human
alz_intact <- alz_intact[alz_intact$Taxid.interactor.A%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
dim(alz_intact) 
alz_intact <- alz_intact[alz_intact$Taxid.interactor.B%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
dim(alz_intact) 

# Filter out interactions with the scores < 0.45
alz_intact$Confidence.value.s. <- as.numeric(sapply(strsplit(alz_intact$Confidence.value.s., split = "miscore:"),'[',2))
alz_intact <- alz_intact[alz_intact$Confidence.value.s. >= 0.45,]
dim(alz_intact[alz_intact$Confidence.value.s.>= 0.45,])

# Convert interactors' names to ensg
# Get interactors' names
A <- unique(as.character(alz_intact[,1]))
B <- unique(as.character(alz_intact[,2]))

# Show the structure
str(A)
A_ID <- sapply(strsplit(A,split = ":"),'[',2)
str(B)
B_ID <- sapply(strsplit(B,split = ":"),'[',2)
AB <- unique(c(A_ID,B_ID))

# Length of the interactors
length(AB)

# Convert the names to ENSG ids
library(gProfileR)
AB2ensg <- gconvert(AB)
AB2ensg <- AB2ensg[!duplicated(AB2ensg), ]
dim(AB2ensg)
head(AB2ensg)

# Protein Complex Interaction(PCI) part of IntAct dataset
# Place interactions like Protein Complex Interaction(PCI) into data structure

# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from alz_intact
alz_intact_ABVAL <- alz_intact[,c(1,2,15)]
str(alz_intact_ABVAL)
alz_intact_ABVAL[,1] <- as.character(alz_intact_ABVAL[,1])
alz_intact_ABVAL[,2] <- as.character(alz_intact_ABVAL[,2])
alz_intact_ABVAL[,3] <- as.character(alz_intact_ABVAL[,3])
colnames(alz_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
alz_intact_ABVAL$A <- sapply(strsplit(alz_intact_ABVAL$A, split = ":"),'[',2)
alz_intact_ABVAL$A <- sapply(strsplit(alz_intact_ABVAL$A, split = "-P"),'[',1)
alz_intact_ABVAL$B <- sapply(strsplit(alz_intact_ABVAL$B, split = ":"),'[',2)
alz_intact_ABVAL$B <- sapply(strsplit(alz_intact_ABVAL$B, split = "-P"),'[',1)

# Get the rows containing protein complexes' IDs
alz_intact_pci_a <- alz_intact_ABVAL[grep("EBI-", alz_intact_ABVAL$A), ]
alz_intact_pci_b <- alz_intact_ABVAL[grep("EBI-",alz_intact_ABVAL$B), ]
alz_intact_pci <- unique(merge(alz_intact_pci_a, alz_intact_pci_b, by.row = T, all = T))

# Convert protein names(where present) to ensg
# Check the dimention
dim(merge(alz_intact_pci, AB2ensg, by.x = "A", by.y = ".id", all.x = T))

# Merge for the first interactor
alz_intact_pci_A <- merge(alz_intact_pci, AB2ensg, by.x = "A", by.y = ".id", all.x = T)

# Merge for the second interactor
alz_intact_pci_AB <- merge(alz_intact_pci_A, AB2ensg, by.x = "B", by.y = ".id", all.x = T)

# Substitute complex names that were not converted with original PCI IDs
alz_intact_pci_AB$Target.x[is.na(alz_intact_pci_AB$Target.x)] <- as.character(alz_intact_pci_AB$A[is.na(alz_intact_pci_AB$Target.x)])
alz_intact_pci_AB$Target.y[is.na(alz_intact_pci_AB$Target.y)] <- as.character(alz_intact_pci_AB$B[is.na(alz_intact_pci_AB$Target.y)])

# Rename the columns
colnames(alz_intact_pci_AB) <- c("B", "A", "score", "ensg1", "ensg2")
head(alz_intact_pci_AB)

# Save converted data with original names and corresponding ENSG IDs
save(alz_intact_pci_AB, file = "alz_intact_pci_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"
save(alz_intact_pci_AB, file = "alz_intact_pci_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"

# Select columns with converted names and score
alz_intact_int_pci <- alz_intact_pci_AB[, c(4,5,3)]

# Rename the columns
colnames(alz_intact_int_pci) <- c("ensg1", "ensg2", "score")
alz_intact_int_pci <- cbind(alz_intact_int_pci, interaction_type = "PCI")
alz_intact_int_pci <- cbind(alz_intact_int_pci, data_source = "ADIA") # source name of the interacions in homo sapiens with miscore >0.45 IntAct
dim(alz_intact_int_pci)  

#Save the part of the integrated dataset from Alzheimer's IntAct pci for human
save(alz_intact_int_pci, file = "alz_intact_int_pci.RData")
write.table(alz_intact_int_pci, file = "alz_intact_int_pci.txt", sep = "\t", quote = F, row.names = F)

### Protein-protein interaction(PPI) part of IntAct dataset
# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from alz_intact
dim(alz_intact)   
alz_intact_ABVAL <- alz_intact[, c(1,2,15)]
str(alz_intact_ABVAL)
alz_intact_ABVAL[,1] <- as.character(alz_intact_ABVAL[,1])
alz_intact_ABVAL[,2] <- as.character(alz_intact_ABVAL[,2])
alz_intact_ABVAL[,3] <- as.character(alz_intact_ABVAL[,3])
colnames(alz_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
alz_intact_ABVAL$A <- sapply(strsplit(alz_intact_ABVAL$A, split = ":"),'[',2)
alz_intact_ABVAL$B <- sapply(strsplit(alz_intact_ABVAL$B, split = ":"),'[',2)
head(alz_intact_ABVAL)

# Merge for the first interactor
alz_intact_ppi_A <- merge(alz_intact_ABVAL, AB2ensg, by.x = "A", by.y = ".id", all = F)

# Merge for the second interactor
alz_intact_ppi_AB <- merge(alz_intact_ppi_A, AB2ensg, by.x = "B", by.y = ".id", all = F)

# Rename the columns
colnames(alz_intact_ppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")
head(alz_intact_ppi_AB)

# Save converted data with original names and corresponding ENSG IDs
save(alz_intact_ppi_AB, file = "alz_intact_ppi_AB.RData") #file that has fields "B","A","score","ensg1","ensg2"

# Select columns "ensg1", "ensg2", "score"
alz_intact_int_ppi <- alz_intact_ppi_AB[, c(4,5,3)]

# Rename the columns
colnames(alz_intact_int_ppi) <- c("ensg1", "ensg2", "score")
alz_intact_int_ppi <- cbind(alz_intact_int_ppi, interaction_type = "PPI")
alz_intact_int_ppi <- cbind(alz_intact_int_ppi, data_source = "ADIA") # source name of the interacions in homo sapiens with miscore >0.45 alz_intact 
dim(alz_intact_int_ppi) 

#Save the part of the integrated dataset from alz_intact ppi for human 
save(alz_intact_int_ppi, file = "alz_intact_int_ppi.RData")
write.table(alz_intact_int_ppi, file = "alz_intact_int_ppi.txt", sep = "\t", quote = F, row.names = F)

### Combined PCI and PPI
# Filter out interactions with the scores < 0.45

# Merge converted PPI and PCI part of alz_intact
alz_intact_int <- rbind(alz_intact_int_ppi, alz_intact_int_pci)
dim(alz_intact_int) 

#Save the part of the integrated dataset from alz_intact PPI and PCI for human
save(alz_intact_int, file = "alz_intact_int.RData")
write.table(alz_intact_int, file = "alz_intact_int.txt", sep = "\t", quote = F, row.names = F)
