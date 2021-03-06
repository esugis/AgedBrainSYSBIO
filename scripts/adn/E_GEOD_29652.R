# This script:
# Extracts Alzheimer’s related and healthy samples from the data sets.
# Filters out probesets with SD < 0.29.
# Calculating the co-expression between all probesets in each of the data sets using Spearman correlation coefficient.
# For each individual probeset in each of the datasets script creates 2 separate files in .txt and .RData formats.
# Files are named after the probeset.
# Created files contain the names of the correlated probesets and the corresponding Spearman coefficient.

# Used libarries
library(foreach); library(doMC); cores=10 ; registerDoMC(cores);
library(stringr);
library("R.utils");
library(ncdf4);

# Open dataset in NetCDF format
# Please indicate the path to the saved .nc file, e.g. as shown below
E_GEOD_29652 <- nc_open("~/absb/data/adn/E-GEOD-29652.nc"); 

# Extract only Alzheimer and healthy saples from E_GEOD_29652
# List of variable in E_GEOD_29652
names(E_GEOD_29652$var)

# Select only Alzheimer's disease and healthy samples
data_E_GEOD_29652 <- ncvar_get(E_GEOD_29652, "data")
rownames(data_E_GEOD_29652) <- ncvar_get(E_GEOD_29652, "MetadataOrder")
colnames(data_E_GEOD_29652)<- ncvar_get(E_GEOD_29652,"gene")
data_E_GEOD_29652 <- t(data_E_GEOD_29652)

# Dimentions of the data
dim(data_E_GEOD_29652)

# Close NetCDF file
nc_lose(E_GEOD_29652)

# Rename the samples
colnames(data_E_GEOD_29652)[1:18]="alz"

# Filter out rows with SD values less then 0.29
SD <- apply(data_E_GEOD_29652, 1, sd, na.rm = T)
data_E_GEOD_29652_filt <- data_E_GEOD_29652[SD >= 0.29, ]
dim(data_E_GEOD_29652_filt)                          
m <- t(data_E_GEOD_29652_filt)

# Path to the foldet where results will be saved
pathRdata <- "~/absb/results/adn/all_probes/rdata/E_GEOD_29652/"
pathtxt <- "~/absb/results/adn/all_probes/txt/E_GEOD_29652/"

# Create directories
dir.create(file.path(pathRdata),showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(pathtxt),showWarnings = FALSE, recursive = TRUE)

# Substitute "-" and "/" in probesets' names with "_"
colnames(m) <- gsub("-", "_",colnames(m))
colnames(m) <- gsub("/", "_",colnames(m))

## Probesets
ds_genes <- colnames(m)

# Compute Sperman correlation between expression profiles of probesets "all against all".
foreach(i = 1:length(ds_genes)) %dopar%{
  cor1gds <- c() # Corelation for one gene in one ds
   gene <- ds_genes[i];
   vect <- m[,i];
   cor1gds <- t(cor(vect, m, method = "spearman"))

   # Save result as txt  and  RData
   filename <- sprintf("%s.txt", gene);
   filedata <- sprintf("%s.RData", gene);
   pathname <- file.path(pathtxt, filename);
   pathdata <- file.path(pathRdata, filedata);
   write.table(cor1gds, file = pathname, sep = "\t", quote=F, row.names = F);
   save(cor1gds,file = pathdata)
 }

# Exract the list of probesets
E_GEOD_29652_pr <- colnames(m)
save(E_GEOD_29652_pr, file = "~/absb/results/adn/all_probes/E_GEOD_29652_all_probes.RData")


