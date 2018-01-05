# This script processes inteacrtion from IntAct with mi score >0.45

# Create the folder where current results will be written
resdir<-paste("~/absb/results","intact",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Read in text file of all interactions with miscore >0.45 downloaded from IntAct and save them as RData
# Please indicate the path to the downloaded data
intact <- read.delim("~/absb/data/intact/intact_hs_v_4_2_6.txt", stringsAsFactors = F, sep = "\t", header=T)
save(intact, file = "intact.RData")

# Get the column names
colnames(intact)

# Select only interactions related to human
intact <- intact[intact$Taxid.interactor.A%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
dim(intact)    
intact <- intact[intact$Taxid.interactor.B%in%"taxid:9606(human)|taxid:9606(Homo sapiens)",]
dim(intact) 

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
length(AB)

# Convert the names to ENSG ids
#library(gProfileR)
#AB2ensg <- gconvert(AB)
#AB2ensg <- AB2ensg[!duplicated(AB2ensg), ]
#dim(AB2ensg)
#head(AB2ensg)
#AB2ensg <- AB2ensg[,c(2,4)]
#colnames(AB2ensg)<-c(".id", "Target")

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

## Convert with gProfileR
#AB_nc2ensg <- gconvert(AB_nc)
#AB_nc2ensg <- AB_nc2ensg[!duplicated(AB_nc2ensg), ]
#dim(AB_nc2ensg)
#head(AB_nc2ensg)
#AB_nc2ensg <- AB_nc2ensg[,c(2,4)]
#colnames(AB_nc2ensg)<-c(".id", "Target")

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



# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from intact
intact_ABVAL <- intact[,c(1,2,15)]
str(intact_ABVAL)
intact_ABVAL[,1] <- as.character(intact_ABVAL[,1])
intact_ABVAL[,2] <- as.character(intact_ABVAL[,2])
intact_ABVAL[,3] <- as.character(intact_ABVAL[,3])
colnames(intact_ABVAL) <- c("A", "B", "score")

# Split column content on ":"
intact_ABVAL$A <- sapply(strsplit(intact_ABVAL$A, split = ":"),'[',2)
#intact_ABVAL$A <- sapply(strsplit(intact_ABVAL$A, split = "-P"),'[',1)
intact_ABVAL$B <- sapply(strsplit(intact_ABVAL$B, split = ":"),'[',2)
#intact_ABVAL$B <- sapply(strsplit(intact_ABVAL$B, split = "-P"),'[',1)
intact_ABVAL$score <- sapply(strsplit(intact_ABVAL$score, split = "miscore:"),'[',2)

# Get the rows containing protein complexes' IDs
#intact_pci_a <- intact_ABVAL[grep("EBI-", intact_ABVAL$A), ]
#intact_pci_b <- intact_ABVAL[grep("EBI-",intact_ABVAL$B), ]
#intact_pci <- unique(merge(intact_pci_a, intact_pci_b, by.row = T, all = T))

# Select the rows where there are unconverted proteins in "A" or "B" columns
rn_A_UP <- rownames(intact_ABVAL[intact_ABVAL$A%in%UP,])
rn_B_UP <- rownames(intact_ABVAL[intact_ABVAL$B%in%UP,])
rn_AB_UP <- unique(c(rn_A_UP,rn_B_UP))

# Select the part of the data with uncharacterised proteins
intact_ucppi <- intact_ABVAL[rownames(intact_ABVAL)%in%rn_AB_UP,]
intact_ucppi$score<-as.numeric(intact_ucppi$score)

# Convert protein names(where present) to ensg
# Check the dimention
dim(merge(intact_ucppi, AB2ensg, by.x = "A", by.y = ".id", all.x = T))

# Merge for the first interactor
intact_ucppi_A <- merge(intact_ucppi, AB2ensg, by.x = "A", by.y = ".id", all.x = T)

# Merge for the second interactor
intact_ucppi_AB <- merge(intact_ucppi_A, AB2ensg, by.x = "B", by.y = ".id", all.x = T)

# Convert to character
intact_ucppi_AB$Target.x[is.na(intact_ucppi_AB$Target.x)] <- as.character(intact_ucppi_AB$A[is.na(intact_ucppi_AB$Target.x)])
intact_ucppi_AB$Target.y[is.na(intact_ucppi_AB$Target.y)] <- as.character(intact_ucppi_AB$B[is.na(intact_ucppi_AB$Target.y)])

# Rename the columns
colnames(intact_ucppi_AB) <- c("B", "A", "score", "ensg1", "ensg2")
head(intact_ucppi_AB)

# Save converted data with original names and corresponding ENSG IDs
save(intact_ucppi_AB, file = "intact_ucppi_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"
save(intact_ucppi_AB, file = "intact_ucppi_AB.RData")#file that has fields "B","A","score","ensg1","ensg2"

# Select columns with converted names and score
intact_int_ucppi <- intact_ucppi_AB[, c(4,5,3)]

# Rename the columns
colnames(intact_int_ucppi) <- c("ensg1", "ensg2", "score")
intact_int_ucppi <- cbind(intact_int_ucppi, interaction_type = "UCPPI")
intact_int_ucppi <- cbind(intact_int_ucppi, data_source = "IAH") # source name of the interacions in homo sapiens with miscore >0.45 IntAct
dim(intact_int_ucppi)

#Save the part of the integrated dataset from IntAct ppi for human
save(intact_int_ucppi, file = "intact_int_ucppi.RData")
write.table(intact_int_ucppi, file = "intact_int_ucppi.txt", sep = "\t", quote = F, row.names = F)


### Protein-protein interaction(PPI) part of IntAct dataset
# Select columns "X.ID.s..interactor.A",  "ID.s..interactor.B"  "Confidence.value.s." from intact
dim(intact)    
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
dim(intact_int_ppi)

#Save the part of the integrated dataset from IntAct ppi for human 
save(intact_int_ppi, file = "intact_int_ppi.RData")
write.table(intact_int_ppi, file = "intact_int_ppi.txt", sep = "\t", quote = F, row.names = F)

### Combined UCPPI and PPI
# Merge converted PPI and PCI part of IntAct
intact_int <- rbind(intact_int_ppi, intact_int_ucppi)
dim(intact_int) 

# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG 0.5
# PPIs from IntAct
intact_int=df2string(intact_int)
str(intact_int)
intact_int <- intact_int[!duplicated(intact_int), ]
dim(intact_int)
intact_int <- intact_int[!duplicated(data.frame(t(apply(intact_int[1:2], 1, sort)), intact_int$score)),]
# New size
dim(intact_int)# !!!!74979     5


# Exclude interaction_type = "UCPPI" from the integrated dataset
intact_int <- intact_int[!intact_int$interaction_type%in%"UCPPI",]
dim(intact_int) #  64774     5

#Save the part of the integrated dataset from IntAct PPI and PCI for human
save(intact_int, file = "intact_int.RData")
write.table(intact_int, file = "intact_int.txt", sep = "\t", quote = F, row.names = F)





