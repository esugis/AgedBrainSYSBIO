# This script does the sanity check for the coexpression data set.
# Parameters that can vary are the p-value cut off(currently 1e-5)
# Currently we used consrvative approach to select interactions based on RRA score. 
# In case of 2 probeset  pairs AB BA, we keep the one with the larger p-value.
# We select one of 2 pars because RRA results of AB and BA are not symmetrical.
# Adds missing columns interaction_type="coexpression" and data_source="ADN"
# The filnal coexpression part of the integrated dataset is stored in 
# file="~/absb/results/adn/integration/adn_coexp_int.RData"
# and file="~/absb/results/adn/integration/adn_coexp_int.txt" 

# Load the dataset of interest
# Load ds with cut off 0.00001
load(file = "~/absb/results/adn/integration/adn_coexp_pairs_int.RData")
dim(adn_coexp_pairs_int)
alzcoexp_int<-adn_coexp_pairs_int

# Remove columns 4, 5
alzcoexp_int<-alzcoexp_int[,1:3]

# Add columns interaction_type, data_source
alzcoexp_int=cbind(alzcoexp_int, interaction_type="coexpression")
alzcoexp_int=cbind(alzcoexp_int, data_source="ADN") # Stands for alzheimer's disease and normal samples
colnames(alzcoexp_int)=c("ensg1","ensg2", "score", "interaction_type","data_source")

# Check the size
dim(alzcoexp_int)

# Sanity check for self-loops 
alzcoexp_int <- subset(alzcoexp_int, !ensg1 == ensg2)

# Size after self-loops were removed
dim(alzcoexp_int)

head(alzcoexp_int)
adn_coexp_int <- alzcoexp_int
save(adn_coexp_int, file="~/absb/results/adn/integration/adn_coexp_int.RData")
write.table(adn_coexp_int,file="~/absb/results/adn/integration/adn_coexp_int.txt", quote=F, sep="\t",row.names=F)


