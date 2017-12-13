# This script assembles matrix of z-scores for all probes in the datasets from individual results for each tissue  

# Set working directory
setwd("~/absb/results/allenbrain/")

# Used libraries
library(reshape2)
library(foreach)

# Path to the results saved in RData format for individual tissues
pathRdata <- "~/absb/results/allenbrain/tissues_rdata/"

# Load the list of tissues
load(file = "~/absb/results/allenbrain/all_tissues.RData" )
tissues_all <-sort(tissues_all)
# Combine z-scores in individual tissues for all probes in the form of one matrix

tissues_z <- foreach(i = 1:length(tissues_all), .combine = rbind)%do%{
tissue=tissues_all[i]
print("current tissue")
print(tissue)
print("current iteration")
print(i)
# Create path to the file
filedata <- sprintf("%s.RData",tissue);
pathdata <- file.path(pathRdata, filedata);
print("path to file")
print(pathdata)
load(file = pathdata)
# Load onetissue_av
tst<-head(onetissue_av)
print("head")
print(tst)
onetissue_av <-onetissue_av[,1:4]
return(onetissue_av)
}

save(tissues_z, file = "tissues_z.RData")

# Cast a matrix
tissues_zscores_mtx <- acast(tissues_z, probe_id~y, value.var = "zscore")

# Save to the file 
save(tissues_zscores_mtx,file = "tissues_zscores_mtx.RData")
write.table(tissues_zscores_mtx, file = "tissues_zscores_mtx.txt", quote = F, row.names = F, sep = "\t")
