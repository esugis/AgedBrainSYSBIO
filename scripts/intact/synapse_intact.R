# This script processes manually curated synapse realted  dataset from IntAct

# Create the folder where current results will be written
resdir<-paste("~/AgedBrainSYSBIO/results","intact",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file of all interactions with miscore >0.45 downloaded from IntAct and save them as RData
# Please indicate the path to the downloaded data
syn_intact <- read.delim("~/AgedBrainSYSBIO/data/intact/synapse_intact_v_4_2_6.txt", stringsAsFactors=F, sep = "\t", header=T)
save(syn_intact, file = "syn_intact.RData")

# Get the column names
colnames(syn_intact)

# Select only interactions related to human
syn_intact <- syn_intact[syn_intact$Taxid.interactor.A%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
dim(syn_intact)     
syn_intact <- syn_intact[syn_intact$Taxid.interactor.B%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
dim(syn_intact)

# Filter out interactions with the scores <0.45
syn_intact$Confidence.value.s. <- as.numeric(sapply(strsplit(syn_intact$Confidence.value.s., split = "miscore:"),'[',2))
syn_intact <- syn_intact[syn_intact$Confidence.value.s. >= 0.45,]
dim(syn_intact[syn_intact$Confidence.value.s.>= 0.45,]) #863 15

# Convert interactors' names to ensg
# Get interactors' names
A <- unique(as.character(syn_intact[,1]))
B <- unique(as.character(syn_intact[,2]))

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

# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from syn_intact
syn_intact_ABVAL <- syn_intact[,c(1,2,15)]
str(syn_intact_ABVAL)
syn_intact_ABVAL[,1] <- as.character(syn_intact_ABVAL[,1])
syn_intact_ABVAL[,2] <- as.character(syn_intact_ABVAL[,2])
syn_intact_ABVAL[,3] <- as.character(syn_intact_ABVAL[,3])
colnames(syn_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
syn_intact_ABVAL$A <- sapply(strsplit(syn_intact_ABVAL$A, split = ":"),'[',2)
syn_intact_ABVAL$A <- sapply(strsplit(syn_intact_ABVAL$A, split = "-P"),'[',1)
syn_intact_ABVAL$B <- sapply(strsplit(syn_intact_ABVAL$B, split = ":"),'[',2)
syn_intact_ABVAL$B <- sapply(strsplit(syn_intact_ABVAL$B, split = "-P"),'[',1)

# Get the rows containing protein complexes' IDs
syn_intact_pci_a <- syn_intact_ABVAL[grep("EBI-", syn_intact_ABVAL$A), ]
syn_intact_pci_b <- syn_intact_ABVAL[grep("EBI-",syn_intact_ABVAL$B), ]
syn_intact_pci <- unique(merge(syn_intact_pci_a, syn_intact_pci_b, by.row = T, all = T))

# Convert protein names(where present) to ensg
# Check the dimention
dim(merge(syn_intact_pci, AB2ensg, by.x = "A", by.y = ".id", all.x = T))

# Merge for the first interactor
syn_intact_pci_A <- merge(syn_intact_pci, AB2ensg, by.x = "A", by.y = ".id", all.x = T)

# Merge for the second interactor
syn_intact_pci_AB <- merge(syn_intact_pci_A, AB2ensg, by.x = "B", by.y = ".id", all.x = T)

# Substitute complex names that were not converted with original PCI IDs
syn_intact_pci_AB$Target.x[is.na(syn_intact_pci_AB$Target.x)] <- as.character(syn_intact_pci_AB$A[is.na(syn_intact_pci_AB$Target.x)])
syn_intact_pci_AB$Target.y[is.na(syn_intact_pci_AB$Target.y)] <- as.character(syn_intact_pci_AB$B[is.na(syn_intact_pci_AB$Target.y)])

# Rename the columns
colnames(syn_intact_pci_AB) <- c("B", "A", "score", "ensg1", "ensg2")
head(syn_intact_pci_AB)

# Save converted data with original names and corresponding ENSG IDs
save(syn_intact_pci_AB, file = "syn_intact_pci_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"
save(syn_intact_pci_AB, file = "syn_intact_pci_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"

# Select columns with converted names and score
syn_intact_int_pci <- syn_intact_pci_AB[, c(4,5,3)]

# Rename the columns
colnames(syn_intact_int_pci) <- c("ensg1", "ensg2", "score")
syn_intact_int_pci <- cbind(syn_intact_int_pci, interaction_type = "PCI")
syn_intact_int_pci <- cbind(syn_intact_int_pci, data_source = "SIA") # source name of the interacions in homo sapiens with miscore >0.45 IntAct
dim(syn_intact_int_pci)  

#Save the part of the integrated dataset from synapse IntAct pci for human
save(syn_intact_int_pci, file = "syn_intact_int_pci.RData")
write.table(syn_intact_int_pci, file = "syn_intact_int_pci.txt", sep = "\t", quote = F, row.names = F)

### Protein-protein interaction(PPI) part of IntAct dataset
# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from syn_intact
dim(syn_intact)     
syn_intact_ABVAL <- syn_intact[, c(1,2,15)]
str(syn_intact_ABVAL)
syn_intact_ABVAL[,1] <- as.character(syn_intact_ABVAL[,1])
syn_intact_ABVAL[,2] <- as.character(syn_intact_ABVAL[,2])
syn_intact_ABVAL[,3] <- as.character(syn_intact_ABVAL[,3])
colnames(syn_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
syn_intact_ABVAL$A <- sapply(strsplit(syn_intact_ABVAL$A, split = ":"),'[',2)
syn_intact_ABVAL$B <- sapply(strsplit(syn_intact_ABVAL$B, split = ":"),'[',2)
head(syn_intact_ABVAL)

# Merge for the first interactor
syn_intact_ppi_A <- merge(syn_intact_ABVAL, AB2ensg, by.x = "A", by.y = ".id", all = F)

# Merge for the second interactor
syn_intact_ppi_AB <- merge(syn_intact_ppi_A, AB2ensg, by.x = "B", by.y = ".id", all = F)

# Rename the columns
colnames(syn_intact_ppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")
head(syn_intact_ppi_AB)

# Save converted data with original names and corresponding ENSG IDs
save(syn_intact_ppi_AB, file = "syn_intact_ppi_AB.RData") #file that has fields "B","A","score","ensg1","ensg2"

# Select columns "ensg1", "ensg2", "score"
syn_intact_int_ppi <- syn_intact_ppi_AB[, c(4,5,3)]

# Rename the columns
colnames(syn_intact_int_ppi) <- c("ensg1", "ensg2", "score")
syn_intact_int_ppi <- cbind(syn_intact_int_ppi, interaction_type = "PPI")
syn_intact_int_ppi <- cbind(syn_intact_int_ppi, data_source = "SIA") # source name of the interacions in homo sapiens with miscore >0.45 syn_intact 
dim(syn_intact_int_ppi)

#Save the part of the integrated dataset from syn_intact ppi for human 
save(syn_intact_int_ppi, file = "syn_intact_int_ppi.RData")
write.table(syn_intact_int_ppi, file = "syn_intact_int_ppi.txt", sep = "\t", quote = F, row.names = F)

### Combined PCI and PPI
# Filter out interactions with the scores < 0.45

# Merge converted PPI and PCI part of syn_intact
syn_intact_int <- rbind(syn_intact_int_ppi, syn_intact_int_pci)
dim(syn_intact_int) 

#Save the part of the integrated dataset from syn_intact PPI and PCI for human
save(syn_intact_int, file = "syn_intact_int.RData")
write.table(syn_intact_int, file = "syn_intact_int.txt", sep = "\t", quote = F, row.names = F)
