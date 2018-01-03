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
alzcoexp_int <- adn_coexp_pairs_int

# Remove columns 4, 5
alzcoexp_int <- alzcoexp_int[,1:3]

# Add columns interaction_type, data_source
alzcoexp_int <- cbind(alzcoexp_int, interaction_type = "coexpression")
alzcoexp_int <- cbind(alzcoexp_int, data_source = "ADN") # Stands for alzheimer's disease and normal samples
colnames(alzcoexp_int) <- c("ensg1","ensg2", "score", "interaction_type","data_source")

# Check the size
dim(alzcoexp_int)

# Remove the duplicated undirrescted edges with the same score.
# For example, ENSG1-ENSG2 0.5 and ENSG2-ENSG1 0.5

# Convert factors to characters
df2string<-function(df){
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df[,3]<-as.numeric(df[,3])
return (df)}

# Co-expression in Alzheimer's and normal brain
alzcoexp_int <- df2string(alzcoexp_int)
str(alzcoexp_int)
alzcoexp_int <- alzcoexp_int[!duplicated(alzcoexp_int), ]
dim(alzcoexp_int)
alzcoexp_int <- alzcoexp_int[!duplicated(data.frame(t(apply(alzcoexp_int[1:2], 1, sort)), alzcoexp_int$score)),]
# New size
dim(alzcoexp_int)

# Size after self-loops and duplicated interactions were removed
dim(alzcoexp_int)

head(alzcoexp_int)
adn_coexp_int <- alzcoexp_int
save(adn_coexp_int, file = "~/absb/results/adn/integration/adn_coexp_int.RData")
write.table(adn_coexp_int,file = "~/absb/results/adn/integration/adn_coexp_int.txt", quote = F, sep = "\t",row.names = F)


