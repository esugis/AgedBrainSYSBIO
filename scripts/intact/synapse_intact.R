# This script processes manually curated synapse realted  dataset from IntAct

# Create the folder where current results will be written
resdir<-paste("~/absb/results","intact",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file of all interactions with miscore >0.45 downloaded from IntAct and save them as RData
# Please indicate the path to the downloaded data
syn_intact <- read.delim("~/absb/data/intact/synapse_intact_v_4_2_6.txt", stringsAsFactors=F, sep = "\t", header=T)
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

# Convert protein ids to ensg ids using BioMart
library(biomaRt)
#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)
mart.pr <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "ensembl.org")
AB2ensg <- getBM(attributes = c("uniprotswissprot", "ensembl_gene_id"),filters=c("uniprotswissprot"), values = AB, mart = mart.pr)

AB2ensg<-AB2ensg[,c(1,2)]
head(AB2ensg)

# Remove duplicates
AB2ensg<-AB2ensg[!duplicated(AB2ensg), ]
colnames(AB2ensg) <- c(".id", "Target")

# Non-converted IDs
AB_nc<-AB[!AB%in%AB2ensg$.id]

# Convert unconverted isoforms that have "-" in the name
AB_nc <- sapply(strsplit(AB_nc,split = "-"),'[',1)
AB_nc_original <- data.frame(cbind(AB[!AB%in%AB2ensg$.id],AB_nc))

# Comvert non-converted IDs using BioMart
AB_nc2ensg <- getBM(attributes = c("uniprotswissprot", "ensembl_gene_id"),filters=c("uniprotswissprot"), values = AB_nc, mart = mart.pr)

AB_nc2ensg<-AB_nc2ensg[,c(1,2)]
head(AB_nc2ensg)

# Remove duplicates
AB_nc2ensg<-AB_nc2ensg[!duplicated(AB_nc2ensg), ]
colnames(AB_nc2ensg) <- c(".id", "Target")

AB_nc2ensg_orig <- merge(AB_nc2ensg, AB_nc_original, by.x=".id", by.y="AB_nc", all=F)
AB_nc2ensg_orig <-AB_nc2ensg_orig[, 2:3]
colnames(AB_nc2ensg_orig)<-c("Target", ".id")
AB_nc2ensg_orig<-AB_nc2ensg_orig[,c(2,1)]

# Combine AB2ensg and AB_nc

AB2ensg<- rbind(AB2ensg, AB_nc2ensg_orig) 

# Interactions with uncharacterised proteins (UCPPI) part of IntAct dataset
# Place interactions with uncharacterised proteins into data structure
# Select all unconverted proteins
UP <-AB[!AB%in%AB2ensg$.id]

# Part of Synapse IntAct dataset describing PPIs with uncharacterised proteins

# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from syn_intact
syn_intact_ABVAL <- syn_intact[,c(1,2,15)]
str(syn_intact_ABVAL)
syn_intact_ABVAL[,1] <- as.character(syn_intact_ABVAL[,1])
syn_intact_ABVAL[,2] <- as.character(syn_intact_ABVAL[,2])
syn_intact_ABVAL[,3] <- as.character(syn_intact_ABVAL[,3])
colnames(syn_intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
syn_intact_ABVAL$A <- sapply(strsplit(syn_intact_ABVAL$A, split = ":"),'[',2)
syn_intact_ABVAL$B <- sapply(strsplit(syn_intact_ABVAL$B, split = ":"),'[',2)

# Select the rows where there are unconverted proteins in "A" or "B" columns
rn_A_UP <- rownames(syn_intact_ABVAL[syn_intact_ABVAL$A%in%UP,])
rn_B_UP <- rownames(syn_intact_ABVAL[syn_intact_ABVAL$B%in%UP,])
rn_AB_UP <- unique(c(rn_A_UP,rn_B_UP))

# Select the part of the data with uncharacterised proteins
syn_intact_ucppi <- syn_intact_ABVAL[rownames(syn_intact_ABVAL)%in%rn_AB_UP,]
syn_intact_ucppi$score<-as.numeric(syn_intact_ucppi$score)

# Convert protein names(where present) to ensg
# Check the dimention
dim(merge(syn_intact_ucppi, AB2ensg, by.x = "A", by.y = ".id", all.x = T))

# Merge for the first interactor
syn_intact_ucppi_A <- merge(syn_intact_ucppi, AB2ensg, by.x = "A", by.y = ".id", all.x = T)

# Merge for the second interactor
syn_intact_ucppi_AB <- merge(syn_intact_ucppi_A, AB2ensg, by.x = "B", by.y = ".id", all.x = T)

# Convert to character
syn_intact_ucppi_AB$Target.x[is.na(syn_intact_ucppi_AB$Target.x)] <- as.character(syn_intact_ucppi_AB$A[is.na(syn_intact_ucppi_AB$Target.x)])
syn_intact_ucppi_AB$Target.y[is.na(syn_intact_ucppi_AB$Target.y)] <- as.character(syn_intact_ucppi_AB$B[is.na(syn_intact_ucppi_AB$Target.y)])

# Rename the columns
colnames(syn_intact_ucppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")
head(syn_intact_ucppi_AB)

# Save converted data with original names and corresponding ENSG IDs
save(syn_intact_ucppi_AB, file = "syn_intact_ucppi_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"
save(syn_intact_ucppi_AB, file = "syn_intact_ucppi_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"

# Select columns with converted names and score
syn_intact_int_ucppi <- syn_intact_ucppi_AB[, c(4,5,3)]

# Rename the columns
colnames(syn_intact_int_ucppi) <- c("ensg1", "ensg2", "score")
syn_intact_int_ucppi <- cbind(syn_intact_int_ucppi, interaction_type = "UCPPI")
syn_intact_int_ucppi <- cbind(syn_intact_int_ucppi, data_source = "SIA") # source name of the interacions in homo sapiens with miscore >0.45 IntAct
dim(syn_intact_int_ucppi)

#Save the part of the integrated dataset from IntAct ppi for human
save(syn_intact_int_ucppi, file = "syn_intact_int_ucppi.RData")
write.table(syn_intact_int_ucppi, file = "syn_intact_int_ucppi.txt", sep = "\t", quote = F, row.names = F)

dim(syn_intact_int_ucppi)

#Save the part of the integrated dataset from IntAct ppi for human
save(syn_intact_int_ucppi, file = "syn_intact_int_ucppi.RData")
write.table(syn_intact_int_ucppi, file = "syn_intact_int_ucppi.txt", sep = "\t", quote = F, row.names = F)

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

# Select columns 4,5,3
syn_intact_int_ppi <- syn_intact_ppi_AB[, c(4,5,3)]

# Rename the columns
colnames(syn_intact_int_ppi) <- c("ensg1", "ensg2", "score")
syn_intact_int_ppi <- cbind(syn_intact_int_ppi, interaction_type = "PPI")
syn_intact_int_ppi <- cbind(syn_intact_int_ppi, data_source = "SIA") # source name of the interacions in homo sapiens with miscore >0.45 IntAct 
dim(syn_intact_int_ppi)

#Save the part of the integrated dataset from Synapse IntAct ppi for human 
save(syn_intact_int_ppi, file = "syn_intact_int_ppi.RData")
write.table(syn_intact_int_ppi, file = "syn_intact_int_ppi.txt", sep = "\t", quote = F, row.names = F)

### Combined UCPPI and PPI
# Merge converted PPI and PCI part of IntAct
syn_intact_int <- rbind(syn_intact_int_ppi, syn_intact_int_ucppi)
dim(syn_intact_int) 

# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG 0.5
# PPIs from IntAct
syn_intact_int=df2string(syn_intact_int)
str(syn_intact_int)
syn_intact_int <- syn_intact_int[!duplicated(syn_intact_int), ]
dim(syn_intact_int)
syn_intact_int <- syn_intact_int[!duplicated(data.frame(t(apply(syn_intact_int[1:2], 1, sort)), syn_intact_int$score)),]
# New size
dim(syn_intact_int)#  [1] 372   5

 Exclude interaction_type = "UCPPI" from the integrated dataset
syn_intact_int <- syn_intact_int[!syn_intact_int$interaction_type%in%"UCPPI",]
dim(syn_intact_int)     #362  5

#Save the part of the integrated dataset from Synapse IntAct PPI and UCPPI for human
save(syn_intact_int, file = "syn_intact_int.RData")
write.table(syn_intact_int, file = "syn_intact_int.txt", sep = "\t", quote = F, row.names = F)

