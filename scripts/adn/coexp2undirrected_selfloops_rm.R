# This script takes the all to all probes interactions
# From 2 corresponding pairs ENSG1-ENSG2 and ENSG2-ENSG1 selects the conservative RRA score
# It also removes interaction with gene itself(self-loop)
# The initial dataset is stored in ~/absb/results/adn/integration/


# Set working directory

resdir <- "~/absb/results/adn/integration/"

# Create directories
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)
setwd(resdir)

# Load the coexpression interactions dataset
load(file = "~/absb/results/adn/integration/alzcoexp_int.RData")

# Correct the colnames

# Add colnames
colnames(alzcoexp_int) <- c("ensg1","ensg2","score","interaction_type","data_source")

# Create new object to work with
alzcoexp_int_ud <-  alzcoexp_int

# Add the pair combination to the df
alzcoexp_int_ud <- alzcoexp_int_ud[,1:3]
alzcoexp_int_ud$ensg12 <- paste(alzcoexp_int_ud$ensg1, alzcoexp_int_ud$ensg2, sep = "_")
alzcoexp_int_ud$ensg21 <- paste(alzcoexp_int_ud$ensg2, alzcoexp_int_ud$ensg1, sep = "_")
#save(alzcoexp_int_ud, file = "alzcoexp_int_ud.RData")

# Create a temporary dataframe
test <- alzcoexp_int_ud
# Summary of the data
summary(test)

# Remove self loops(interactions between the same genes)
test  <-  subset(test, !ensg1  ==  ensg2)

# Compare with full data 
# Size of filtered data
dim(test)
# Size of the full data
dim(alzcoexp_int_ud)

# Search for the rows with the same interactors A_B B_A and select the most conservative p-value
outDF <- data.frame(test[1,])
for (i in 1:length(test$ensg12)){
   	AB <- test[i,4]
   	BA <- test[i,5]
   	s_AB <- test[i,3]
	s_BA <- test[test$ensg12%in%BA,3]
	
	if(length(rownames(test[test$ensg12%in%BA,]))>0){
  		if(s_AB > s_BA ){
   			outDF <- rbind(outDF, test[i,])
   		} else if(s_AB  ==  s_BA ){
        		rBA <- as.numeric(as.character(rownames(test[test$ensg12%in%BA,])))
        		rAB <- as.numeric(as.character(rownames(test[test$ensg12%in%AB,])))
        		maxr <- max(rBA,rAB)
        		outDF <- rbind(outDF,test[rownames(test)%in%maxr,])
        	} else if(s_AB < s_BA ){
        		outDF <- rbind(outDF, test[test$AB%in%BA,])
        		}
	 } else if(length(rownames(test[test$ensg12%in%BA,])) == 0) {
        	outDF <- rbind(outDF, test[i,])
        	} 
}

outDF <- outDF[-1,]
outDF <- outDF[!duplicated(outDF),]
adn_coexp_pairs_int<- outDF
save(adn_coexp_pairs_int , file = "~/absb/results/adn/integration/adn_coexp_pairs_int.RData")
write.table(adn_coexp_pairs_int, file = "~/absb/results/adn/integration/adn_coexp_pairs_int.txt", quote = F, sep = "\t",row.names = F)
dim(adn_coexp_pairs_int) 
