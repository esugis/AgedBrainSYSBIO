# This script expracts ensg name, score from positiv selection dataset.

# Create the folder where current results will be written
resdir<-paste("~/AgedBrainSYSBIO/results","ps",sep="/")
dir.create(file.path(resdir),showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Load the data
ps <- read.csv(file = "~/AgedBrainSYSBIO/data/ps/positive_selection_data.csv", sep  = ",", header = T, stringsAsFactors = F)

# Extract the columns containing ENSG id and p-value
ps <- ps[,c(1,8)]

# Rename columns
colnames(ps) <- c("ensg", "ps_pval")

# Remove NAs
ps <- na.omit(ps)

# Save the data
save(ps, file = "ps.RData")
write.table(ps, file = "ps.txt",sep = "\t", row.names = F, col.names = T, quote = F)
