# This script combines expression data with tissue annotations in 178238373 data set
# Processed data is a data frame with columns probe_id, sample_num, value, structure_id

# Create the folder where current results will be written
resdir <- "~/absb/results/allenbrain/178238373_ds"
dir.create(file.path(resdir), showWarnings = FALSE, recursive = TRUE)

# Set created directory as working dirrectory
setwd(resdir)

# Used librarries
library(reshape)

# Read the expression from MicroarrayExpression.csv  Ontology.csv  PACall.csv  Probes.csv  Readme.txt  SampleAnnot.csv 
maexp <- read.csv(file = "~/absb/data/allenbrain/178238373_ds/MicroarrayExpression.csv", header = F)
dim(maexp)

# Read ontology
onto <- read.csv(file = "~/absb/data/allenbrain/178238373_ds/Ontology.csv")

# Read in probes names
probes <- read.csv("~/absb/data/allenbrain/178238373_ds/Probes.csv")
dim(probes)

# Read in sample annotations
sampleann <- read.csv("~/absb/data/allenbrain/178238373_ds/SampleAnnot.csv")
dim(sampleann)
rownames(maexp) <- maexp[,1]
maexp <- maexp[,-1]
maexp.melt <- melt(t(maexp))
maexp.melt$X1 <- gsub("V","",maexp.melt$X1)

# Add sample numbers to merge by it
sampleann <- cbind(sampleann, rownames(sampleann))
colnames(sampleann)[14] <- "samplenum"
sampleann.short <- sampleann[,c(1,5,6,14)]
 i <- sapply(sampleann.short, is.factor)
sampleann.short[i] <- lapply(sampleann.short[i], as.character)
sampleann.short$samplenum <- as.numeric(sampleann.short$samplenum)

# Merge with sample annotations
sample_tissueid <- sampleann.short[,c(1,4)]
maexp.merged <- merge(maexp.melt,sample_tissueid, by.x = "X1", by.y = "samplenum", all = F)
maexp.sort <- maexp.merged[order(maexp.merged$X2),]
maexp.probes <- merge(maexp.sort, probes, by.x = "X2", by.y = "probe_id", all = F)
colnames(maexp.probes)[1:2] <- c("probe_id","sample_num")
maexp.probes <- maexp.probes[,c(1,2,3,4)]
head(maexp.probes)

# Save the preprocessed data
save(maexp.probes,file = "maexp.probes_178238373_ds.RData")
write.table(maexp.probes, file = "maexp.probes_178238373_ds.txt", sep = "\t", quote = F, row.names = F)


