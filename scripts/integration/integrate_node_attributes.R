# This script assembles the nodes attribute file from GWAS, Allen Brain atlas expression, 
# and  adds biotype of all the genes(gwas, allen brain atlas, IntAct, pba, epistasis, co-expression)

# Create the folder where current results will be written
resdir<-paste("~/absb/results","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Load individual  datasets of gene attributes

load(file="~/absb/results/gwas/gwas_ensg.RData")
load(file="~/absb/results/allenbrain/mtx_zscores_max_sorted.RData")
load(file="~/absb/results/ps/ps.RData")

# Load integrated interactions
load(file="~/absb/results/integration/integrated_int.RData")

# Add gene names and biotype to the interactors
rownames(integrated_int) <- 1:nrow(integrated_int)

library(biomaRt)
library(tidyr)
library(dplyr)
# Split by ENSG*-ENSG*
integrated_int_split <-  integrated_int %>% separate(ensg.A, c("ENSG.A1", "ENSG.A2"), "-") %>% separate(ensg.B, c("ENSG.B1", "ENSG.B2"), "-")
integrated_int  <- cbind(integrated_int,integrated_int_split[,c(1:4)])

# Select all individual ENSG
ensg_integrated_int <- unique(c(integrated_int$ENSG.A1,integrated_int$ENSG.A2,integrated_int$ENSG.B1,integrated_int$ENSG.B2))

# Convert to gene names
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "ensembl.org")
ensg2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name","gene_biotype"), filters=c("ensembl_gene_id"),values=ensg_integrated_int, mart=mart)
colnames(ensg2name) <- c("ensg", "gene_name", "biotype")

rm_ensg <- c(integrated_int$ensg.A,integrated_int$ensg.B)
integrated_int <- integrated_int[!integrated_int$ensg.A%in%rm_ensg[grep("\\ENSG.*-$",rm_ensg)],]
integrated_int <- integrated_int[!integrated_int$ensg.B%in%rm_ensg[grep("\\ENSG.*-$",rm_ensg)],]
integrated_int <- integrated_int[!integrated_int$ensg.A%in%rm_ensg[grep("^-.*",rm_ensg)],]
integrated_int <- integrated_int[!integrated_int$ensg.B%in%rm_ensg[grep("^-.*",rm_ensg)],]


# Combined gene names with the original data
integrated_int <- merge(integrated_int, ensg2name, by.x = "ENSG.A1", by.y = "ensg", all.x=T)
colnames(integrated_int) <- gsub("gene_name","gene_name.A1", colnames(integrated_int))  
colnames(integrated_int) <- gsub("biotype","biotype.A1", colnames(integrated_int)) 
integrated_int <- merge(integrated_int, ensg2name, by.x = "ENSG.A2", by.y = "ensg", all.x=T)
colnames(integrated_int) <- gsub("\\gene_name$","gene_name.A2", colnames(integrated_int))   
colnames(integrated_int) <- gsub("biotype$","biotype.A2", colnames(integrated_int))   

integrated_int <- merge(integrated_int, ensg2name, by.x = "ENSG.B1", by.y = "ensg", all.x=T)
colnames(integrated_int) <- gsub("\\gene_name$","gene_name.B1", colnames(integrated_int))      
colnames(integrated_int) <- gsub("biotype$","biotype.B1", colnames(integrated_int)) 
integrated_int <- merge(integrated_int, ensg2name, by.x = "ENSG.B2", by.y = "ensg", all.x=T)
colnames(integrated_int) <- gsub("\\gene_name$","gene_name.B2", colnames(integrated_int))   
colnames(integrated_int) <- gsub("biotype$","biotype.B2", colnames(integrated_int)) 

# Create colunms corresponding to gene names
integrated_int <- cbind(integrated_int, gene_name.A=paste(integrated_int$gene_name.A1, integrated_int$gene_name.A2, sep="-"), 
gene_name.B=paste(integrated_int$gene_name.B1, integrated_int$gene_name.B2, sep="-"))
integrated_int <- cbind(integrated_int, biotype.A=paste(integrated_int$biotype.A1, integrated_int$biotype.A2, sep="-"),
 biotype.B=paste(integrated_int$biotype.B1, integrated_int$biotype.B2, sep="-"))

# Remove part with NA
integrated_int$gene_name.A <- gsub("NA-", "", integrated_int$gene_name.A)
integrated_int$gene_name.A <- gsub("-NA", "", integrated_int$gene_name.A)
integrated_int$gene_name.B <- gsub("NA-", "", integrated_int$gene_name.B)
integrated_int$gene_name.B <- gsub("-NA", "", integrated_int$gene_name.B)

integrated_int$biotype.A <- gsub("NA-", "", integrated_int$biotype.A)
integrated_int$biotype.A <- gsub("-NA", "", integrated_int$biotype.A)
integrated_int$biotype.B <- gsub("NA-", "", integrated_int$biotype.B)
integrated_int$biotype.B <- gsub("-NA", "", integrated_int$biotype.B)

integrated_int <- integrated_int[,colnames(integrated_int)%in%c("ensg.A", "ensg.B","score", "interaction_type", "data_source", 
                                  "gene_name.A", "gene_name.B",
                                  "biotype.A", "biotype.B")]
dim(integrated_int)

# Extract ensg ids

# GWAS
#ensg_gwas <- unique(as.character(as.vector(gwas_ensg$ensg)))
#length(ensg_gwas)

# AllenBrain atlas data
colnames(mtx_zscores_max_sorted)[2] <- "ensg"
#ensg_allen <- unique(as.character(as.vector(mtx_zscores_max_sorted$ensg)))
allen <- mtx_zscores_max_sorted
#length(ensg_allen)

# Positve selection
#ensg_ps <- unique(as.character(as.vector(ps$ensg)))
#length(ensg_ps)

# Integrated interactions
#int_ensgs <- integrated_int
#ensg_int <- unique(as.character(c(int_ensgs[,1], int_ensgs[,2])))

# Combine all ensg ids
#ensg <- unique(c(ensg_int,ensg_ps, ensg_allen, ensg_gwas))
#length(ensg) 

# Get gene names and biotype
# Biomart
#library(biomaRt)
#mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="ensembl.org")
#ensg2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "gene_biotype"), filters=c("ensembl_gene_id"),values=ensg, mart=mart)
#colnames(ensg2name) <- c("ensg", "gene_name", "biotype")

# Save converted ensg ids and corresponding name
#write.table(ensg2name, file = "ensg2name.txt", sep="\t", quote=F, row.names=F)
#save(ensg2name,file = "ensg2name.RData")





# Merge individual attribute data together based on ensg

# Combine gwas and positive selection data 
node_attr <- merge(gwas_ensg, ps, by.x="ensg",by.y="ensg", all=T)

# Combine with allen brain expression
node_attr <- merge(node_attr, allen, by.x="ensg",by.y="ensg", all=T)

# Combine with gene names and biotype
#node_attr <- merge(ensg2name, node_attr, by.x="ensg",by.y="ensg", all=T)

# Remove duplicates
node_attr<-node_attr[!duplicated(node_attr), ]
dim(node_attr)# [1] 61490   420

# Extract gene names and biotype from integrated_int

ensg_int <- merge(integrated_int, integrated_int, by.x=ensg.A, by.y=ensg.B, all=T)




# Save to file
write.table(node_attr, file="node_attributes.txt", sep="\t", quote=F, row.names=F )
save(node_attr, file="node_attributes.RData")

# Combine interactions and attributes
integrated_int_1attr <- merge(integrated_int, node_attr, by.x="ensg.A", by.y="ensg", all.x=T)
integrated_int_attr <-  merge(integrated_int_1attr, node_attr, by.x="ensg.B", by.y="ensg", all.x=T)
integrated_int_attr <-integrated_int_attr[!duplicated(integrated_int_attr), ]
integrated_int_attr <- cbind(integrated_int_attr$ensg.A, integrated_int_attr)
integrated_int_attr <- integrated_int_attr[, -3]
colnames(integrated_int_attr)[1] <- "ensg.A"

# Rename columns .x corresponds to the second interactor, .y corresponds to the first interactor
colnames(integrated_int_attr)<-gsub("\\.x", ".A", colnames(integrated_int_attr))
colnames(integrated_int_attr)<-gsub("\\.y", ".B", colnames(integrated_int_attr))
colnames(integrated_int_attr)[1:2]<- c("ensg.A", "ensg.B")

# Save to file
write.table(integrated_int_attr, file="integrated_int_attributes.txt", sep="\t", quote=F, row.names=F )
save(integrated_int_attr, file="integrated_int_attributes.RData")


# Network for the web
# Select brain regions of interest
# DG, CA1, CA2, CA

selected_columns <- c("gene_name.A","ensg.A", "gene_name.B","ensg.B", "score","interaction_type","data_source","biotype.A", "SNP_id.A", "GWAS_pvalue.A", "ps_pval.A","biotype.B",
"SNP_id.B", "GWAS_pvalue.B", "ps_pval.B", "DG.A", "CA1.A","CA2.A", "CA3.A", "DG.B", "CA1.B", "CA2.B", "CA3.B")
integrated_int_attr_web <- integrated_int_attr[, colnames(integrated_int_attr)%in%selected_columns]
integrated_int_attr_web <- integrated_int_attr_web[,c("gene_name.A","gene_name.B","ensg.A","ensg.B","score","interaction_type","data_source","biotype.A","biotype.B","SNP_id.A",
"GWAS_pvalue.A","ps_pval.A","SNP_id.B","GWAS_pvalue.B", "ps_pval.B","DG.A", "CA1.A","CA2.A", "CA3.A","DG.B", "CA1.B", "CA2.B", "CA3.B")]

write.table(integrated_int_attr_web, file="integrated_int_attributes_web.txt", sep="\t", quote=F, row.names=F )
save(integrated_int_attr_web, file="integrated_int_attributes_web.RData")

# R Shiny requires factors
net  <- integrated_int_attr_web
character_vars <- lapply(net, class) == "character"
net[, character_vars] <- lapply(net[, character_vars], as.factor)

write.table(net, file="integrated_int_attributes_web_net.txt", sep="\t", quote=F, row.names=F )
save(net, file="integrated_int_attributes_web_net.RData")



