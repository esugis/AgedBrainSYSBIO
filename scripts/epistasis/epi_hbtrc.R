# Script splits complex ensg1-ensg2 names and formats the interactions to ensg1/ensg2/p-value/int
# it also converts LRG ids to ENSG ids
library(gProfileR)

# Read the file related to HTBRC cohort
epi=read.table(file="~/AgedBrainSYSBIO/data/epistatis/HBTRC_epistatic_gene_pairs_Braak_1e-8.tsv", header=T)

# Create the folder where current results will be written
resdir <- "~/AgedBrainSYSBIO/results/epistasis/"
dir.create(file.path(resdir),showWarnings  =  FALSE, recursive  =  TRUE)

# Data size
dim(epi)

# Set created directory as working dirrectory
setwd(resdir)

epi_a<- epi[grep("LRG_", epi$ENSG_A), ]
epi_b<- epi[grep("LRG_", epi$ENSG_B), ]

# Number of LRG genes among interactors A
dim(epi_a)

# Number of LRG genes among interactors B
dim(epi_b)
rows_cut=row.names(epi_b)

# Convert LRG to ENSG
B=as.character(epi_b$ENSG_B)
LRG2ensg_hbtrc = gconvert(B)

# Write results
save(LRG2ensg_hbtrc, file="LRG2ensg_hbtrc.RData")
write.table(LRG2ensg_hbtrc, file="LRG2ensg_hbtrc.txt", quote=F,sep="\t", row.names=F)

# Remove duplicates
LRG2ensg_hbtrc=LRG2ensg_hbtrc[!duplicated(LRG2ensg_hbtrc), ]
# Size
dim(LRG2ensg_hbtrc)

# Column names
colnames(epi_b)

# Megre
dim(merge(epi_b,LRG2ensg_hbtrc, by.x="ENSG_B", by.y=".id", all=F))
epi_hbtrc_lrg=merge(epi_b,LRG2ensg_hbtrc, by.x="ENSG_B", by.y=".id", all=F)

# Select only ensg1 ensg2 and score
epi_hbtrc_lrg_fin=epi_hbtrc_lrg[,c(2,8,3)]

# Data head
epi_hbtrc_lrg_fin 

# Bind columns with interaction_type, data_source
epi_hbtrc_lrg_fin=cbind(epi_hbtrc_lrg_fin, interaction_type="epistasis")
epi_hbtrc_lrg_fin=cbind(epi_hbtrc_lrg_fin, data_source="HBTRC")

# Rename the coluns
colnames(epi_hbtrc_lrg_fin)=c("ensg1","ensg2","score","interaction_type","data_source")

# Write results
save(epi_hbtrc_lrg_fin, file="epi_hbtrc_lrg.RData")
write.table(epi_hbtrc_lrg_fin,file="epi_hbtrc_lrg.txt",sep="\t", quote=F, row.names=F)

# Combine with the main data frame 
epi_cut=epi[!row.names(epi)%in%rows_cut,]
dim(epi_cut)
dim(epi)
epi_cut_fin=epi_cut[,c(1,2,3)]
epi_cut_fin=cbind(epi_cut_fin,interaction_type="epistasis")
epi_cut_fin=cbind(epi_cut_fin, data_source="HBTRC")
colnames(epi_cut_fin)=c("ensg1","ensg2","score","interaction_type","data_source")
epi_hbtrc=epi_cut_fin

# Convert gene ids and ensg id to tha latest Ensembl version
epi_hbtrc_ensg12ensg = gconvert(epi_hbtrc$ensg1)
dim(epi_hbtrc_ensg12ensg)

# Remove duplicates
epi_hbtrc_ensg12ensg=epi_hbtrc_ensg12ensg[!duplicated(epi_hbtrc_ensg12ensg), ]
dim(epi_hbtrc_ensg12ensg)

# Convert ids
epi_hbtrc_ensg22ensg = gconvert(epi_hbtrc$ensg2)
dim(epi_hbtrc_ensg22ensg)

# Remove duplicates
epi_hbtrc_ensg22ensg=epi_hbtrc_ensg22ensg[!duplicated(epi_hbtrc_ensg22ensg), ]
dim(epi_hbtrc_ensg22ensg)

# Merge by ensg1
epi_hbtrc_2ensg=merge(epi_hbtrc,epi_hbtrc_ensg12ensg, by.x="ensg1", by.y=".id", all=F)
dim(epi_hbtrc_2ensg)

# Merge by ensg2
epi_hbtrc_2ensg=merge(epi_hbtrc_2ensg,epi_hbtrc_ensg22ensg, by.x="ensg2", by.y=".id", all=F)

# Size of the dataset with old ENSG IDs(ver 74) converted to ver 87
dim(epi_hbtrc_2ensg)

# Size of the dataset with old ENSG IDs(ver74)
dim(epi_hbtrc)

# Find differences between ENSG IDs in Ensembl ver 87 and Ensembl ver 74
# All ENSG IDs in ver 87
ensg_hbtrc_ver87 <- unique(c(as.character(epi_hbtrc_2ensg$ensg1),as.character(epi_hbtrc_2ensg$ensg2)))
length(ensg_hbtrc_ver87)

# All ENSG IDs in ver 74
ensg_hbtrc_ver74<-unique(c(as.character(epi_hbtrc$ensg1),as.character(epi_hbtrc$ensg2)))
length(ensg_hbtrc_ver74)    

# ENSG IDs present in ver 74 and not present in ver 87
ensg_87_vs_74_hbtrc<-ensg_hbtrc_ver74[!ensg_hbtrc_ver74%in%ensg_hbtrc_ver87]
length(ensg_87_vs_74_hbtrc)#

# Save to file
write.table(ensg_87_vs_74_hbtrc, file="ensg_87_vs_74_hbtrc.txt", quote=F, row.names=F, sep="\t")


epi_hbtrc<-epi_hbtrc_2ensg[,c(6,7,3,4,5)]
colnames(epi_hbtrc)=c("ensg1","ensg2","score","interaction_type","data_source")
epi_hbtrc=epi_hbtrc[!duplicated(epi_hbtrc), ]
dim(epi_hbtrc)

# Write results
save(epi_hbtrc, file="epi_hbtrc.RData")
write.table(epi_hbtrc,file="epi_hbtrc.txt",sep="\t", quote=F, row.names=F)

# Combine dataframes of LRG genes interactions and the rest
epi_hbtrc_int=rbind(epi_hbtrc,epi_hbtrc_lrg_fin)


# Save the part of the integrated dataset related to hbtrc cohort
save(epi_hbtrc_int, file="epi_hbtrc_int.RData")
write.table(epi_hbtrc_int,file="epi_hbtrc_int.txt",sep="\t", quote=F, row.names=F)




