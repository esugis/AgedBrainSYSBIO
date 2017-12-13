# This script installs required packages

# gProfiler
install.packages("gProfileR") 
library(gProfileR) 

install.packages("reshape2")
library(reshape2)

install.packages("foreach")
library(foreach)

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
