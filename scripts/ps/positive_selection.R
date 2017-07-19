# This script expracts ensg name, score from positiv selection dataset.

# Set working directory
setwd("~/absb/results/ps/")

# Load the data
ps <- read.csv(file = "~/absb/data/ps/positive_selection_data.csv", sep  = ",", header = T, stringsAsFactors = F)

# Extract the columns containing ENSG id and p-value
ps <- ps[,c(1,8)]

# Rename columns
colnames(ps) <- c("ensg", "ps_pval")

# Remove NAs
ps <- na.omit(ps)

# Save the data
save(ps, file = "ps.RData")
write.table(ps, file = "ps.txt",sep = "\t", row.names = F, col.names = T, quote = F)
