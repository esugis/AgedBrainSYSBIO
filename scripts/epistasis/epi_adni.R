# Script splits complex ensg1-ensg2 names and formats the interactions to ensg1/ensg2/p-value/int
# it also converts LRG ids to ENSG ids

library(gProfileR)

# Read the file related to HTBRC cohort
epi <- read.table(file = "~/AgedBrainSYSBIO/data/epistasis/epistatic_bin-bin_interactions_at_MAF_prod_GE_0.01_bonf_pvalue_LE_1.0E-11.tsv", header=T)

# Create the folder where current results will be written
resdir <- "~/AgedBrainSYSBIO/results/epistasis/"
dir.create(file.path(resdir),showWarnings  =  FALSE, recursive  =  TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Convert to character
epi$bin_1 <- as.character(epi$bin_1)

# Grep all the interactions that has intergenic regions
epi_a <- epi[grep("-ENSG", epi$bin_1), ]
epi_b <- epi[grep("-ENSG", epi$bin_2), ]

# Number of LRG genes among interactors 
dim(epi_a)
dim(epi_b)

rows_cut_a <- row.names(epi_a)
rows_cut_b <- row.names(epi_b)
rows_cut <- unique(c(rows_cut_a, rows_cut_b))
length(rows_cut)

# Intergenic regions
epi_adni_igri <- epi[rownames(epi)%in%rows_cut,c(1,2,3)]

# Bind columns with interaction_type, data_source.
epi_adni_igri <- cbind(epi_adni_igri, interaction_type="IGRI")
epi_adni_igri <- cbind(epi_adni_igri, data_source="ADNI_VER")

# df that has intercators as intergenic regions in one of the columns and correcponding scores
colnames(epi_adni_igri) <- c("ensg1","ensg2","score","interaction_type","data_source")

# Save data for intergenic regions interactions(IGRI)
save(epi_adni_igri, file = "epi_adni_igri.RData")
write.table(epi_adni_igri,file = "epi_adni_igri.txt",sep = "\t", quote = F, row.names = F)

# Size
dim(epi_adni_igri)

# Combine with the main data frame 
epi_adni <- epi[!row.names(epi)%in%rows_cut,c(1,2,3)]
dim(epi_adni)

# Bind columns with interaction_type, data_source
epi_adni <- cbind(epi_adni, interaction_type = "epistasis")
epi_adni <- cbind(epi_adni, data_source = "ADNI_VER")
colnames(epi_adni) <- c("ensg1","ensg2","score","interaction_type","data_source")

# Convert gene ids and ensg id to tha latest Ensembl version
# First interactor
epi_adni_ensg12ensg <- gconvert(epi_adni$ensg1)
dim(epi_adni_ensg12ensg)
epi_adni_ensg12ensg <- epi_adni_ensg12ensg[!duplicated(epi_adni_ensg12ensg), ]
dim(epi_adni_ensg12ensg)

# Second interactor
epi_adni_ensg22ensg <- gconvert(epi_adni$ensg2)
dim(epi_adni_ensg22ensg)
epi_adni_ensg22ensg <- epi_adni_ensg22ensg[!duplicated(epi_adni_ensg22ensg), ]
dim(epi_adni_ensg22ensg)

# Merge by ensg1
epi_adni_2ensg <- merge(epi_adni,epi_adni_ensg12ensg, by.x="ensg1", by.y=".id", all=F)
dim(epi_adni_2ensg)

#Merge by ensg2
epi_adni_2ensg <- merge(epi_adni_2ensg,epi_adni_ensg22ensg, by.x="ensg2", by.y=".id", all=F)

# Size of the dataset with old ENSG IDs(ver 74) converted to ver 887
dim(epi_adni_2ensg)
# Size of the dataset with old ENSG IDs(ver74)
dim(epi_adni)

# Find differences between ENSG IDs in Ensembl ver 87 and Ensembl ver 74
# All ENSG IDs in ver 87
ensg_adni_ver87 <- unique(c(as.character(epi_adni_2ensg$ensg1),as.character(epi_adni_2ensg$ensg2)))
length(ensg_adni_ver87)

# All ENSG IDs in ver 74
ensg_adni_ver74 <- unique(c(as.character(epi_adni$ensg1),as.character(epi_adni$ensg2)))
length(ensg_adni_ver74)       

# ENSG IDs present in ver 74 and not present in ver 87
ensg_87_vs_74_adni <- ensg_adni_ver74[!ensg_adni_ver74%in%ensg_adni_ver87]
length(ensg_87_vs_74_adni)

# Save to file
write.table(ensg_87_vs_74_adni, file="ensg_87_vs_74_adni.txt", quote=F, row.names=F, sep="\t")

# Find differences in IGRI
epi_adni_igri_ensg <- unique(c(as.character(epi_adni_igri$ensg1),as.character(epi_adni_igri$ensg2)))
length(epi_adni_igri_ensg)
library(stringr)
epi_adni_igri_ensg<-unique(unlist(str_split(epi_adni_igri_ensg,"-")))
length(epi_adni_igri_ensg)

# Convert gene ids and ensg id to tha latest Ensembl version
epi_adni_igri_ensg12ensg <- gconvert(epi_adni_igri_ensg)
dim(epi_adni_igri_ensg12ensg)

epi_adni_igri_ensg12ensg <- epi_adni_igri_ensg12ensg[!duplicated(epi_adni_igri_ensg12ensg), ]
dim(epi_adni_igri_ensg12ensg)

epi_adni_igri_ensg_ver87 <- unique(as.character(epi_adni_igri_ensg12ensg$Target))
epi_adni_igri_ensg_ver74 <- epi_adni_igri_ensg

ensg_87_vs_74_adni_igri <- epi_adni_igri_ensg_ver74[!epi_adni_igri_ensg_ver74%in%epi_adni_igri_ensg_ver87]
length(ensg_87_vs_74_adni_igri)

# Save 
write.table(ensg_87_vs_74_adni_igri, file="ensg_87_vs_74_adni_igri.txt", quote=F, row.names=F, sep="\t")

# Find overlap between missing IDs in IGRI and normal gene IDs
ensg_87_vs_74_adni_all <- unique(c(ensg_87_vs_74_adni,ensg_87_vs_74_adni_igri))
length(ensg_87_vs_74_adni_all)
write.table(ensg_87_vs_74_adni_all, file = "ensg_87_vs_74_adni_all.txt", quote=F, row.names=F, sep="\t")

# Exclude the rows with old ids that are not mapped to Ensembl  ver 87
#epi_adni_igri
excluderow <- c()
for(i in 1:length(rownames(epi_adni_igri))){
 pair1 <- unique(unlist(str_split(as.character(epi_adni_igri[i,1]),"-")))
 p1 <- c(pair1%in%ensg_87_vs_74_adni_all)
 c1 <- length(p1[p1==T]) 
 pair2 <- unique(unlist(str_split(as.character(epi_adni_igri[i,2]),"-")))
 p2 <- c(pair2%in%ensg_87_vs_74_adni_all)
 c2 <- length(p2[p2==T])
if(c1>0){ 
  excluderow <- c(excluderow,i)
  }else if(c2>0){ 
        excluderow <- c(excluderow,i)
        }
}
epi_adni_igri_ver87 <- epi_adni_igri[!rownames(epi_adni_igri)%in%excluderow,]
save(epi_adni_igri_ver87, file = "epi_adni_igri_ver87.RData")
write.table(epi_adni_igri_ver87,file = "epi_adni_igri_ver87.txt", sep = "\t", quote = F, row.names = F)
dim(epi_adni_igri_ver87)
dim(epi_adni_igri)

epi_adni <- epi_adni_2ensg[,c(6,7,3,4,5)]
colnames(epi_adni)=c("ensg1","ensg2","score","interaction_type", "data_source")
epi_adni <- epi_adni[!duplicated(epi_adni), ]
dim(epi_adni)

# Write result
save(epi_adni, file="epi_adni.RData")
write.table(epi_adni,file="epi_adni.txt",sep="\t", quote=F, row.names=F)

# Combine dataframes of IGRI and the rest genes interactions. Write final ds to the file
epi_adni_ver_int <- rbind(epi_adni,epi_adni_igri_ver87)

# Save the part of the integrated dataset related to hbtrc cohort
save(epi_adni_ver_int, file="epi_adni_ver_int.RData")
write.table(epi_adni_ver_int,file="epi_adni_ver_int.txt",sep="\t", quote=F, row.names=F)




