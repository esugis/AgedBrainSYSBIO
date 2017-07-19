# This script processes inteacrtion from IntAct with mi score >0.45

# Create the folder where current results will be written
resdir<-paste("~/absb/results","intact",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file of all interactions with miscore >0.45 downloaded from IntAct and save them as RData
# Please indicate the path to the downloaded data
intact=read.delim("~/AgedBrain/ppi2_data/intact/intact_hs_v_4_2_6.txt", stringsAsFactors = F, sep = "\t", header=T)
save(intact, file = "intact.RData")

# Get the column names
colnames(intact)

# Select only interactions related to human
intact <- intact[intact$Taxid.interactor.A%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
dim(intact) #167470 15   
intact <- intact[intact$Taxid.interactor.B%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
dim(intact) #162228  15 

# Convert interactors' names to ensg
# Get interactors' names
A <- unique(as.character(intact[,1]))
B <- unique(as.character(intact[,2]))

# Show the structure
str(A)
A_ID <- sapply(strsplit(A,split = ":"),'[',2)
str(B)
B_ID <- sapply(strsplit(B,split = ":"),'[',2)
AB <- unique(c(A_ID,B_ID))

# Length of the interactors
length(AB) #12695

# Convert the names to ENSG ids
library(gProfileR)
AB2ensg <- gconvert(AB)
AB2ensg <- AB2ensg[!duplicated(AB2ensg), ]
dim(AB2ensg) #9942    2
head(AB2ensg)

# Protein Complex Interaction(PCI) part of IntAct dataset
# Place interactions like Protein Complex Interaction(PCI) into data structure

# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from intact
intact_ABVAL <- intact[,c(1,2,15)]
str(intact_ABVAL)
intact_ABVAL[,1] <- as.character(intact_ABVAL[,1])
intact_ABVAL[,2] <- as.character(intact_ABVAL[,2])
intact_ABVAL[,3] <- as.character(intact_ABVAL[,3])
colnames(intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
intact_ABVAL$A <- sapply(strsplit(intact_ABVAL$A, split = ":"),'[',2)
intact_ABVAL$A <- sapply(strsplit(intact_ABVAL$A, split = "-P"),'[',1)
intact_ABVAL$B <- sapply(strsplit(intact_ABVAL$B, split = ":"),'[',2)
intact_ABVAL$B <- sapply(strsplit(intact_ABVAL$B, split = "-P"),'[',1)
intact_ABVAL$score <- sapply(strsplit(intact_ABVAL$score, split = "miscore:"),'[',2)

# Get the rows containing protein complexes' IDs
intact_pci_a <- intact_ABVAL[grep("EBI-", intact_ABVAL$A), ]
intact_pci_b <- intact_ABVAL[grep("EBI-",intact_ABVAL$B), ]
intact_pci <- unique(merge(intact_pci_a, intact_pci_b, by.row = T, all = T))

# Convert protein names(where present) to ensg
# Check the dimention
dim(merge(intact_pci, AB2ensg, by.x = "A", by.y = ".id", all.x = T))

# Merge for the first interactor
intact_pci_A <- merge(intact_pci, AB2ensg, by.x = "A", by.y = ".id", all.x = T)

# Merge for the second interactor
intact_pci_AB <- merge(intact_pci_A, AB2ensg, by.x = "B", by.y = ".id", all.x = T)

# Substitute complex names that were not converted with original PCI IDs
intact_pci_AB$Target.x[is.na(intact_pci_AB$Target.x)] <- as.character(intact_pci_AB$A[is.na(intact_pci_AB$Target.x)])
intact_pci_AB$Target.y[is.na(intact_pci_AB$Target.y)] <- as.character(intact_pci_AB$B[is.na(intact_pci_AB$Target.y)])

# Rename the columns
colnames(intact_pci_AB) <- c("B", "A", "score", "ensg1", "ensg2")
head(intact_pci_AB)

# Save converted data with original names and corresponding ENSG IDs
save(intact_pci_AB, file = "intact_pci_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"
save(intact_pci_AB, file = "intact_pci_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"

# Select columns with converted names and score
intact_int_pci <- intact_pci_AB[, c(4,5,3)]

# Rename the columns
colnames(intact_int_pci) <- c("ensg1", "ensg2", "score")
intact_int_pci <- cbind(intact_int_pci, interaction_type = "PCI")
intact_int_pci <- cbind(intact_int_pci, data_source = "IAH") # source name of the interacions in homo sapiens with miscore >0.45 IntAct
dim(intact_int_pci) # 4216  5

#Save the part of the integrated dataset from IntAct ppi for human
save(intact_int_pci, file = "intact_int_pci.RData")
write.table(intact_int_pci, file = "intact_int_pci.txt", sep = "\t", quote = F, row.names = F)

### Protein-protein interaction(PPI) part of IntAct dataset
# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from intact
dim(intact) #162228   15    
intact_ABVAL <- intact[, c(1,2,15)]
str(intact_ABVAL)
intact_ABVAL[,1] <- as.character(intact_ABVAL[,1])
intact_ABVAL[,2] <- as.character(intact_ABVAL[,2])
intact_ABVAL[,3] <- as.character(intact_ABVAL[,3])
colnames(intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
intact_ABVAL$A <- sapply(strsplit(intact_ABVAL$A, split = ":"),'[',2)
intact_ABVAL$B <- sapply(strsplit(intact_ABVAL$B, split = ":"),'[',2)
intact_ABVAL$score <- sapply(strsplit(intact_ABVAL$score, split = "miscore:"),'[',2)
head(intact_ABVAL)

# Merge for the first interactor
intact_ppi_A <- merge(intact_ABVAL, AB2ensg, by.x = "A", by.y = ".id", all = F)

# Merge for the second interactor
intact_ppi_AB <- merge(intact_ppi_A, AB2ensg, by.x = "B", by.y = ".id", all = F)

# Rename the columns
colnames(intact_ppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")
head(intact_ppi_AB)

# Save converted data with original names and corresponding ENSG IDs
save(intact_ppi_AB, file = "intact_ppi_AB.RData") #file that has fields "B","A","score","ensg1","ensg2"

# Select columns 4,5,3
intact_int_ppi <- intact_ppi_AB[, c(4,5,3)]

# Rename the columns
colnames(intact_int_ppi) <- c("ensg1", "ensg2", "score")
intact_int_ppi <- cbind(intact_int_ppi, interaction_type = "PPI")
intact_int_ppi <- cbind(intact_int_ppi, data_source = "IAH") # source name of the interacions in homo sapiens with miscore >0.45 IntAct 
dim(intact_int_ppi) # 118502    5

#Save the part of the integrated dataset from IntAct ppi for human 
save(intact_int_ppi, file = "intact_int_ppi.RData")
write.table(intact_int_ppi, file = "intact_int_ppi.txt", sep = "\t", quote = F, row.names = F)

### Combined PCI and PPI
# Merge converted PPI and PCI part of IntAct
intact_int <- rbind(intact_int_ppi, intact_int_pci)
dim(intact_int) #122718

#Save the part of the integrated dataset from IntAct PPI and PCI for human
save(intact_int, file = "intact_int.RData")
write.table(intact_int, file = "intact_int.txt", sep = "\t", quote = F, row.names = F)



