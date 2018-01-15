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

install.packages("stringr")
library(stringr)

install.packages("R.utils")
library("R.utils")

install.packages("ncdf4")
library(ncdf4)

install.packages("RobustRankAggreg")
library(RobustRankAggreg)

install.packages("doMC")
library(doMC)

install.packages("tidyr")
libray(tidyr)

install.packaged("dplyr")
library(dplyr)
