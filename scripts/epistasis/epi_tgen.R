# Script preprocesses file containing epistatic interactions in TGEN cohort.
library(gProfileR)

# Read the file
epi <- read.table(file = "~/absb/data/epistasis/TGEN_epistasis.tsv", header = T)

# Create the folder where current results will be written
resdir <- "~/absb/results/epistasis/"
dir.create(file.path(resdir),showWarnings  =  FALSE, recursive  =  TRUE)

# Data size
dim(epi)

# Create the folder where current results will be written
resdir <- "~/absb/results/epistasis/"
dir.create(file.path(resdir),showWarnings  =  FALSE, recursive  =  TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Find LRG genes
epi_a <- epi[grep("LRG_", epi$ENSG_A), ]
epi_b <- epi[grep("LRG_", epi$ENSG_B), ]

# Number of LRG genes among interactors A
dim(epi_a)

# Number of LRG genes among interactors B
dim(epi_b)#
rows_cut <- row.names(epi_b)

# epi_b contains part of epi df that contains LRG genes 

# Convert LRG to ENSG
B <- as.character(epi_b$ENSG_B)
LRG2ensg_tgen <- gconvert(B)
LRG2ensg_tgen <- LRG2ensg_tgen[, c(2,4)]
colnames(LRG2ensg_tgen) <- c(".id", "Target")

# Save to the file in RData and txt formats 
#save(LRG2ensg_tgen, file = "LRG2ensg_tgen.RData")
#write.table(LRG2ensg_tgen, file = "LRG2ensg_tgen.txt", quote = F,sep = "\t",row.names = F)

# Remove duplicates
LRG2ensg_tgen <- LRG2ensg_tgen[!duplicated(LRG2ensg_tgen), ]
dim(LRG2ensg_tgen)

# Merge with epi_b
colnames(epi_b)

dim(merge(epi_b,LRG2ensg_tgen, by.x = "ENSG_B", by.y = ".id", all = F))
epi_tgen_lrg <- merge(epi_b,LRG2ensg_tgen, by.x = "ENSG_B", by.y = ".id", all = F)

# Select only ensg1 ensg2 and score
epi_tgen_lrg_fin <- epi_tgen_lrg[,c(2,8,3)]

# Print head of the data
head(epi_tgen_lrg_fin) 
 
# Bind columns with interaction_type, data_source.
epi_tgen_lrg_fin <- cbind(epi_tgen_lrg_fin, interaction_type = "epistasis")

epi_tgen_lrg_fin <- cbind(epi_tgen_lrg_fin, data_source = "TGEN")
colnames(epi_tgen_lrg_fin) <- c("ensg1","ensg2","score","interaction_type","data_source")

# Save data
#save(epi_tgen_lrg_fin, file  =  "epi_tgen_lrg.RData")
#write.table(epi_tgen_lrg_fin,file = "epi_tgen_lrg_fin.txt",sep = "\t", quote = F, row.names = F)

# Combine with the main data frame 
epi_cut <- epi[!row.names(epi)%in%rows_cut,]
dim(epi_cut)
dim(epi)

# Take only ensg1 ensg2 and score
epi_cut_fin <- epi_cut[,c(1,2,3)]
epi_cut_fin <- cbind(epi_cut_fin,interaction_type = "epistasis")
epi_cut_fin <- cbind(epi_cut_fin, data_source = "TGEN")
colnames(epi_cut_fin) <- c("ensg1","ensg2","score","interaction_type","data_source")
epi_tgen <- epi_cut_fin

# Convert gene ids and ensg id to tha latest Ensembl version
epi_tgen_ensg12ensg <- gconvert(epi_tgen$ensg1)
epi_tgen_ensg12ensg <- epi_tgen_ensg12ensg[, c(2,4)]
colnames(epi_tgen_ensg12ensg) <- c(".id", "Target") 
dim(epi_tgen_ensg12ensg)

# Remove duplicates
epi_tgen_ensg12ensg <- epi_tgen_ensg12ensg[!duplicated(epi_tgen_ensg12ensg), ]
dim(epi_tgen_ensg12ensg)

# Convert second column of interactors
epi_tgen_ensg22ensg <- gconvert(epi_tgen$ensg2)
epi_tgen_ensg22ensg <- epi_tgen_ensg22ensg[, c(2,4)]
colnames(epi_tgen_ensg22ensg) <- c(".id", "Target")
dim(epi_tgen_ensg22ensg)

# Remove duplicates
epi_tgen_ensg22ensg <- epi_tgen_ensg22ensg[!duplicated(epi_tgen_ensg22ensg), ]
dim(epi_tgen_ensg22ensg)

# Merge by ensg1
epi_tgen_2ensg <- merge(epi_tgen,epi_tgen_ensg12ensg, by.x = "ensg1", by.y = ".id", all = F)
dim(epi_tgen_2ensg)

# Merge by ensg2
epi_tgen_2ensg <- merge(epi_tgen_2ensg,epi_tgen_ensg22ensg, by.x = "ensg2", by.y = ".id", all = F)

# Size of the dataset with old ENSG IDs(ver 74) converted to ver 90
dim(epi_tgen_2ensg)

# Size of the dataset with old ENSG IDs(ver74)
dim(epi_tgen)

# Find differences between ENSG IDs in Ensembl ver 90 and Ensembl ver 74 
# All ENSG IDs in ver 90
ensg_tgen_ver90<-unique(c(as.character(epi_tgen_2ensg$ensg1),as.character(epi_tgen_2ensg$ensg2)))
length(ensg_tgen_ver90)

# All ENSG IDs in ver 74
ensg_tgen_ver74<-unique(c(as.character(epi_tgen$ensg1),as.character(epi_tgen$ensg2)))
length(ensg_tgen_ver74)

# ENSG IDs present in ver 74 and not present in ver 90
ensg_90_vs_74_tgen<-ensg_tgen_ver74[!ensg_tgen_ver74%in%ensg_tgen_ver90]
length(ensg_90_vs_74_tgen)#

# Save to file
#write.table(ensg_90_vs_74_tgen, file = "ensg_90_vs_74_tgen.txt", quote = F, row.names = F, sep = "\t")


epi_tgen<-epi_tgen_2ensg[,c(6,7,3,4,5)]
colnames(epi_tgen) <- c("ensg1","ensg2","score","interaction_type","data_source")
epi_tgen <- epi_tgen[!duplicated(epi_tgen), ]
dim(epi_tgen)

# Write to files
#save(epi_tgen, file = "epi_tgen.RData")
#write.table(epi_tgen,file = "epi_tgen.txt",sep = "\t",quote = F, row.names = F)

# Combine dataframes of LRG genes interactions and the rest
epi_tgen_int <- rbind(epi_tgen,epi_tgen_lrg_fin)

# Remove the duplicated undirrescted edges with the same score.
# For example ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5

# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

epi_tgen_int <- df2string(epi_tgen_int)
str(epi_tgen_int)
dim(epi_tgen_int)
epi_tgen_int <- epi_tgen_int[!duplicated(data.frame(t(apply(epi_tgen_int[1:2], 1, sort)), epi_tgen_int[,c(3,5)])),]
# New size
dim(epi_tgen_int)




#Save the part of the integrated dataset related to TGEN cohort
save(epi_tgen_int, file = "epi_tgen_int.RData")
write.table(epi_tgen_int,file = "epi_tgen_int.txt",sep = "\t",quote = F, row.names = F)




