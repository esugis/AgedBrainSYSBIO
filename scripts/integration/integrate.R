# This script combines all indiviadual interaction datasets into one and removes duplicated interaction of one type

# Create the folder where current results will be written
resdir<-paste("~/AgedBrainSYSBIO/results","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Libraries
library(gProfileR)

# Load individual interaction datasets

# PPIs associated with brain ageing (PBA)
load(file="~/AgedBrainSYSBIO/results/pba/pba_int.RData")

# Combined epistatic interactions ADNI ventrical volume, ADNI cognitive traits, TGEN, HBTRC
load(file="~/AgedBrainSYSBIO/results/epistasis/epistasis_all_int.RData")

# Intact hyman PPIs with MIscore >=0.45
load(file="~/AgedBrainSYSBIO/results/intact/intact_int.RData")

# IntAct Alzheimer's related manually curated PPI dataset with MIscore >=0.45
load(file="~/AgedBrainSYSBIO/results/intact/alz_intact_int.RData")

# IntAct Synapse related automatically curated PPI dataset with MIscore >=0.45
load(file="~/AgedBrainSYSBIO/results/intact/syn_intact_int.RData")

# Co-expression dataset with removed self-loops RRA score <= 0.0001 
load(file="~/AgedBrainSYSBIO/results/adn/integration/adn_coexp_int.RData")

# Convert factors to characters

df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# PBA
pba_int <- df2string(pba_int)
# Structure
str(pba_int)
#Size
dim(pba_int)  

# Epistatic interactions
epistasis_all_int <- df2string(epistasis_all_int)
str(epistasis_all_int)
dim(epistasis_all_int) 

# PPIs from IntAct
intact_int=df2string(intact_int)
str(intact_int)
intact_int <- intact_int[!duplicated(intact_int), ]
size(intact_int)

# Alzheimer's related PPIs from IntAct
alz_intact_int=df2string(alz_intact_int)
str(alz_intact_int)
alz_intact_int <- alz_intact_int[!duplicated(alz_intact_int), ]
dim(alz_intact_int)

# Synaptic interactions from IntAct
syn_intact_int=df2string(syn_intact_int)
str(syn_intact_int)
syn_intact_int <- syn_intact_int[!duplicated(syn_intact_int), ]
dim(syn_intact_int)

# Co-expression in Alzheimer's and normal brain
adn_coexp_int=df2string(adn_coexp_int)
str(adn_coexp_int)
adn_coexp_int <- adn_coexp_int[!duplicated(adn_coexp_int), ]
dim(adn_coexp_int)

# Create one DF from separate datasets
integrated_int <- rbind(pba_int, epistasis_all_int, intact_int, alz_intact_int,syn_intact_int,adn_coexp_int)
integrated_int <- integrated_int[!duplicated(integrated_int),]
colnames(integrated_int)<-c("ensg.A","ensg.B")

save(integrated_int, file = "integrated_int.RData")
write.table(integrated_int, file = "integrated_int.txt", sep="\t", quote=F, row.names=F)


# Size of integrated dataset
dim(integrated_int)

