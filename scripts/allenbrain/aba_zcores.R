y'# This script calculates z-scores in Allen Brain Atlas datsets

# Set working directory
setwd("~/AgedBrainSYSBIO/results/allenbrain/")

# Used libaries
library(reshape)
library(foreach)
	
# Path to the results  
pathRdata <- "~/AgedBrainSYSBIO/results/allenbrainatlas/tissues_rdata/"

# Load melted ds probe_id sample_id value  for 178236545
load(file = "~/AgedBrainSYSBIO/results/allenbrain/178236545_ds/maexp.probes_178236545_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178236545 <- maexp.probes
t1 <- unique(maexp_178236545[,4])
print("tissue 1")
length(t1)

# Load processed data for 178238266
load(file = "~/AgedBrainSYSBIO/results/allenbrain/178238266_ds/maexp.probes_178238266_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238266 <- maexp.probes
t2 <- unique(maexp_178238266[,4])
print("tissue 2")
length(t2)

# Load processed data for 178238316
load(file = "~/AgedBrainSYSBIO/results/allenbrain/178238316_ds/maexp.probes_178238316_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238316 <- maexp.probes
t3 <- unique(maexp_178238316[,4])
print("tissue 3")
length(t3)

# Load processed data for 178238359
load(file = "~/AgedBrainSYSBIO/results/allenbrain/178238359_ds/maexp.probes_178238359_ds.RData")
maexp.probes <- maexp.probes[order(maexp.probes[, c(4,1)]),]
maexp_178238359 <- maexp.probes
t4 <- unique(maexp_178238359[,4])
print("tissue4 ")
length(t4)

# Select unique tissues in all 4 data sets
s1 <- read.csv(file = "~/AgedBrainSYSBIO/data/allenbrain/178236545_ds/SampleAnnot.csv")
s2 <- read.csv(file = "~/AgedBrainSYSBIO/data/allenbrain/178238266_ds/SampleAnnot.csv")
s3 <- read.csv(file = "~/AgedBrainSYSBIO/data/allenbrain/178238316_ds/SampleAnnot.csv")
s4 <- read.csv(file = "~/AgedBrainSYSBIO/data/allenbrain/178238359_ds/SampleAnnot.csv")
t11 <- unique(s1$structure_id)
t22 <- unique(s2$structure_id)
t33 <- unique(s3$structure_id)
t44 <- unique(s4$structure_id)

# List of unique tissues in all 4 data sets
tissues_all <- unique(c(t11,t22,t33,t44))

# Number of tissues
length(tissues_all)


# Save the list of tissues
save(tissues_all,file = "all_tissues.RData" )

# Compute z-scores in each of the tissues
tissues_z <- foreach(i = 1:length(tissues_all), .combine = rbind)%do%{

# Subset one tissue
tissue <- tissues[i]
print("current tissue")
print(tissue)

# Subset tissue related data from each ds
onetissue1 <- subset(maexp_178236545, structure_id==tissue)
onetissue2 <- subset(maexp_178238266, structure_id==tissue)
onetissue3 <- subset(maexp_178238316, structure_id==tissue)
onetissue4 <- subset(maexp_178238359, structure_id==tissue)

# Bind tissue related data from each ds together
onetissue <- rbind(onetissue1,onetissue2,onetissue3,onetissue4)

# Aggregate based on probe name, take mean
onetissue_av <- merge(aggregate(value ~ probe_id, onetissue, mean),tissue)

# Calculate z-score over one tissue
onetissue_av <- cbind(onetissue_av, zscore = ((onetissue_av[,2] - mean(onetissue_av[,2])) / sd(onetissue_av[,2])))

# Create path to to save results for one tissue
filedata <- sprintf("%s.RData",tissue);
pathdata <- file.path(pathRdata, filedata);

# Write file for each tissue
save(onetissue_av,file = pathdata)

}


