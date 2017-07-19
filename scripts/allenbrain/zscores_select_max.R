# Script selects maximum z-scores for probes corresponding to the same Ensembl gene name

# Set working dirrectory
setwd("~/absb/results/allenbrain/")

# Load z-score matrix
load(file="tissues_zscores_mtx_ensg.RData")

# Define function to select max z-score for probes corrsponding to the same ENSG id
absmax<-function(x){m=x[which.max( abs(x) )];return(m)}

# Subset only tissue columns from the dataframe
tissues_zscores_mtx_ensg_a=tissues_zscores_mtx_ensg[,c(1:222,224)]

# Aggregare based on ENSG id, take absolute max of z-scores 
mtx_zscores_max=aggregate(. ~ Target,tissues_zscores_mtx_ensg_a, absmax)

# Save aggregated data
save(mtx_zscores_max, file="mtx_zscores_max.RData")
write.table(mtx_zscores_max, file="mtx_zscores_max.txt", sep="\t",quote=F, row.names=F)

# Dimentions of the data
dim(mtx_zscores_max)

# Add tissue short names as colnames.
onto=read.csv(file="~/absb/data/allenbrain/Ontology.csv")
t_names=colnames(mtx_zscores_max)

# Take only colnames corresponding to  tissue names
t_names=t_names[3:223]
t_names_sub=gsub("X","",t_names)
onto_selected=onto[onto$id%in%t_names_sub,]

# Save selected ontologies
save(onto_selected, file="onto_selected.RData")

# Acronyms
onto_acr=onto_selected[,1:2]
onto_acr$id=paste("X",onto_acr$id, sep="")

# Map the values from the dataframe to colnames
# Order cols by name
mtx_zscores_max_sorted=mtx_zscores_max[,order(names(mtx_zscores_max))] 
onto_acr_sorted=onto_acr[order(onto_acr$id),]

# Rename the coluns
colnames(mtx_zscores_max_sorted)[3:223]=onto_acr_sorted$acronym
acr=onto_acr_sorted$acronym
acr=as.character(acr)

# Rename the column
colnames(mtx_zscores_max_sorted)[3:223]=acr  

# Save the results
save(mtx_zscores_max_sorted, file="mtx_zscores_max_sorted.RData")
write.table(mtx_zscores_max_sorted, file="mtx_zscores_max_sorted.txt", quote=F, row.names=F, sep="\t")
