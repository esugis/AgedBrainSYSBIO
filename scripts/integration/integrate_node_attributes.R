# This script assembles the nodes attribute file from GWAS, Allen Brain atlas expression, 
# and  adds biotype of all the genes(gwas, allen brain atlas, IntAct, pba, epistasis, co-expression)

# Set working directory
setwd("~/absb/results/integration/")

# Load individual  datasets of gene attributes

load(file="~/absb/results/gwas/gwas_ensg.RData")
load(file="~/absb/results/allenbrain/mtx_zscores_max_sorted.RData")
load(file="~/absb/results/ps/ps.RData")

# Load integrated interactions
load(file="~/absb/results/integration/integrated_int.RData")


# Extract ensg ids

# GWAS
ensg_gwas <- unique(as.character(as.vector(gwas_ensg$ensg)))
length(ensg_gwas)
#[1] 2034

# AllenBrain atlas data
colnames(mtx_zscores_max_sorted)[2] <- "ensg"
ensg_allen <- unique(as.character(as.vector(mtx_zscores_max_sorted$ensg)))
allen <- mtx_zscores_max_sorted
length(ensg_allen)
#[1] 19849

# Positve selection
ensg_ps <- unique(as.character(as.vector(ps$ensg)))
length(ensg_ps)
# 22

# Integrated interactions
int_ensgs <- integrated_int
ensg_int <- unique(as.character(c(int_ensgs[,1], int_ensgs[,2])))

# Size
length(ensg_int)
#37567

# Combine all ensg ids

ensg <- unique(c(ensg_int,ensg_ps, ensg_allen, ensg_gwas))
length(ensg) #37567

# Get gene names and biotype
# Biomart
library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl", host="ensembl.org")
ensg2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "gene_biotype"), filters=c("ensembl_gene_id"),values=ensg, mart=mart)
colnames(ensg2name) <- c("ensg", "gene_name", "biotype")

# Save converted ensg ids and corresponding name
write.table(ensg2name, file = "ensg2name.txt", sep="\t", quote=F, row.names=F)
save(ensg2name,file = "ensg2name.RData")


# Merge individual attribute data together based on ensg

# Combine gwas and positive selection data 
node_attr <- merge(gwas_ensg, ps, by.x="ensg",by.y="ensg", all=T)

# Combine with allen brain expression
node_attr <- merge(node_attr, allen, by.x="ensg",by.y="ensg", all=T)

# Remove probe names from Allen brain atlas
node_attr<-node_attr[ , -which(names(node_attr) %in% c("probe_id"))]
  
# Combine with gene names and biotype
node_attr <- merge(ensg2name, node_attr, by.x="ensg",by.y="ensg", all=T)

# Remove duplicates
node_attr<-node_attr[!duplicated(node_attr), ]


# Save to file
write.table(node_attr, file="node_attributes.txt", sep="\t", quote=F, row.names=F )
save(node_attr, file="node_attributes.RData")

# Combine interactions and attributes
integrated_int_1attr <- merge(integrated_int, node_attr, by.x="ensg1", by.y="ensg", all.x=T)
integrated_int_attr <-	merge(integrated_int_1attr, node_attr, by.x="ensg2", by.y="ensg", all.x=T)
integrated_int_attr <-integrated_int_attr[!duplicated(integrated_int_attr), ]
integrated_int_attr <- cbind(integrated_int_attr$ensg1, integrated_int_attr)
integrated_int_attr <- integrated_int_attr[, -3]
colnames(integrated_int_attr)[1] <- "ensg1"

# Rename columns .x corresponds to the second interactor, .y corresponds to the first interactor
colnames(integrated_int_attr)<-gsub(".y", ".1", colnames(integrated_int_attr))
colnames(integrated_int_attr)<-gsub(".x", ".2", colnames(integrated_int_attr))
colnames(integrated_int_attr)[colnames(integrated_int_attr)%in%"bio.1pe.1"] <- "biotype.1"
colnames(integrated_int_attr)[colnames(integrated_int_attr)%in%"bio.1pe.2"] <-	"biotype.2"
colnames(integrated_int_attr)[colnames(integrated_int_attr)%in%"interaction_.1pe"] <- "interaction_type"

# Save to file
write.table(integrated_int_attr, file="integrated_int_attributes.txt", sep="\t", quote=F, row.names=F )
save(integrated_int_attr, file="integrated_int_attributes.RData")


# Network for the web
# Select brain regions of interest
# DG, CA1, CA2, CA	

selected_columns <- c("gene_name.1","ensg1", "gene_name.2","ensg2", "score","interaction_type","data_source","biotype.1", "SNP_id.1", "GWAS_pvalue.1", "ps_pval.1","biotype.2", 
"SNP_id.2", "GWAS_pvalue.2", "ps_pval.2", "DG.1", "CA1.1","CA2.1", "CA3.1", "DG.2", "CA1.2", "CA2.2", "CA3.2")
integrated_int_attr_web <- integrated_int_attr[, colnames(integrated_int_attr)%%selected_columns]
integrated_int_attr_web <- integrated_int_attr_web[,c("gene_name.1","ensg1","SNP_id.1", "biotype.1","GWAS_pvalue.1","ps_pval.1",
"DG.1", "CA1.1","CA2.1", "CA3.1","gene_name.2","ensg2","biotype.2", "SNP_id.2", "GWAS_pvalue.2", "ps_pval.2","DG.2", "CA1.2", "CA2.2", "CA3.2")]

write.table(integrated_int_attr_web, file="integrated_int_attributes_web.txt", sep="\t", quote=F, row.names=F )
save(integrated_int_attr_web, file="integrated_int_attributes_web.RData")
