# This script combines all indiviadual interaction datasets into one and removes duplicated interaction of one type

# Create the folder where current results will be written
resdir<-paste("~/absb/results","integration",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Load individual interaction datasets

# PPIs associated with brain ageing (PBA)
load(file="~/absb/results/pba/pba_int.RData")

# Combined epistatic interactions ADNI ventrical volume, ADNI cognitive traits, TGEN, HBTRC
load(file="~/absb/results/epistasis/epistasis_all_int.RData")

# Intact hyman PPIs with MIscore >=0.45
load(file="~/absb/results/intact/intact_int.RData")

# IntAct Alzheimer's related manually curated PPI dataset with MIscore >=0.45
load(file="~/absb/results/intact/alz_intact_int.RData")

# IntAct Synapse related automatically curated PPI dataset with MIscore >=0.45
load(file="~/absb/results/intact/syn_intact_int.RData")

# Co-expression dataset with removed self-loops RRA score <= 0.00001 
load(file="~/absb/results/adn/integration/adn_coexp_int.RData")

# Size of PBA dataset
dim(pba_int) 

# Size of all epistatic interactions dataset
dim(epistasis_all_int)

# Size of PPIs from IntAct
dim(intact_int) 

# Size of Alzheimer's related PPIs from IntAct
dim(alz_intact_int) 

# Size of Synaptic interactions from IntAct
dim(syn_intact_int)

# Size of co-expression in Alzheimer's and normal brain
dim(adn_coexp_int) 

# Create one DF from separate datasets
integrated_int <- rbind(pba_int, epistasis_all_int, intact_int, alz_intact_int, syn_intact_int, adn_coexp_int)
integrated_int <- integrated_int[!duplicated(integrated_int),]
colnames(integrated_int)[1:2]<-c("ensg.A","ensg.B")

save(integrated_int, file = "integrated_int.RData")
write.table(integrated_int, file = "integrated_int.txt", sep="\t", quote=F, row.names=F)


# Size of integrated dataset
dim(integrated_int) #144015     5

