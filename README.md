# AgedBrainSYSBIO
Bringing together heterogeneous datasets related to Alzheimer's disease.

This set of  R scripts were used to create integrated dataset described in the article:

Bringing together heterogeneous datasets related to Alzheimer's disease
Elena Sügis1,2† , Jerome Dauvillier3† , Priit Adler1,2, Valerie Hindie4 , Thomas Moncion4 , Vincent Collura4 , Rachel Daudin5,6, Yann Loe-Mie5,6,Yann Herault7 , Jean-Charles Lambert8 , Henning Hermjakob9 ,Tal Pupko10, Jean-Christophe Rain4 , Ioannis Xenarios3,11, Jaak Vilo1,2, Michel Simonneau5,6*, Hedi Peterson1,2*
*corresponding author(s): Michel Simonneau (michel.simonneau@inserm.fr), Hedi Peterson (peterson@quretec.com) †these authors contributed equally to this work

For detailed description of the methods, please see the publication.

This file contains information about the structure of the files and some useful notes.

## CLONING REPOSITORY
To clone repository to your local machine use the following command:  
git clone https://github.com/esugis/AgedBrainSYSBIO.git  
By default the project will be cloned to the folder "AgedBrainSYSBIO".  
The following project structure will be created:  
AgedBrainSYSBIO  
AgedBrainSYSBIO  
              |── README.md  
              |── data  
              │   |── adn  
              │   |── allenbrain  
              │   |── epistasis  
              │   |── gwas  
              │   |── intact  
              │   |── pba  
              │   |── ps  
              |── data_downloader.sh  
              └── scripts  
                  |── adn  
                  │   |── E_GEOD_18309.R  
                  │   |── E_GEOD_28146.R  
                  │   |── E_GEOD_29652.R  
                  │   |── E_GEOD_4757.R  
                  │   |── E_GEOD_5281.R  
                  │   |── E_MEXP_2280.R  
                  │   |── RRA_probesets.R  
                  │   |── affy2ensg.R  
                  │   |── coexp2undirrected_selfloops_rm.R  
                  │   |── coexp_int.R  
                  │   |── final_adn.R  
                  │   └── tmp_coexp_int.R.save  
                  |── allenbrain  
                  │   |── aba_178236545_ds.R  
                  │   |── aba_178238266_ds.R  
                  │   |── aba_178238316_ds.R  
                  │   |── aba_178238359_ds.R  
                  │   |── aba_178238373_ds.R  
                  │   |── aba_178238387_ds.R  
                  │   |── aba_zcores.R  
                  │   |── aba_zcores_assemble.R  
                  │   |── expressed_regions.R  
                  │   |── probes2ensg.R  
                  │   └── zscores_select_max.R  
                  |── epistasis  
                  │   |── combine.R  
                  │   |── epi_adni.R  
                  │   |── epi_adni_cog.R  
                  │   |── epi_hbtrc.R  
                  │   └── epi_tgen.R  
                  |── gwas  
                  │   └── gwas2ensg.R  
                  |── intact  
                  │   |── alz_intact.R  
                  │   |── intact.R  
                  │   └── synapse_intact.R  
                  |── integration  
                  │   |── integrate.R  
                  │   └── integrate_node_attributes.R  
                  |── libraries.R  
                  |── pba  
                  │   └── pba_int.R  
                  └── ps  
  
You can also clone the project to the folder of your choice as following:  

git clone https://github.com/esugis/AgedBrainSYSBIO.git <DestinationDirectory>  

## STRUCTURE

The following projects is organised in the following structure that you can use in your analysis:
The main catalogue of the project is ~/AgedBrainSYSBIO (where ~ stands for your home dirrectory).
Project is divided into 3 major parts:
 ~/AgedBrainSYSBIO/data  a place to store downloaded datasets
 ~/AgedBrainSYSBIO/scripts a place to store analysis scripts
 ~/AgedBrainSYSBIO/results a place to store the results

## DATA

Data is stored in the folder ~/AgedBrainSYSBIO/data
Folder ~/absb/data contains sub-folders for individual data types used in the analysis.
If you follow the structure of the project , please use the folders to upload the corresponding data sets.
 ~/AgedBrainSYSBIO/data/adn is for microarrays in .nc format
 ~/AgedBrainSYSBIO/data/epistasis is for epistatic datasets
 ~/AgedBrainSYSBIO/data/pba is for PPI dataset related to brain ageing
 ~/AgedBrainSYSBIO/data/intact is for PPI datasets downloaded from IntAct database
 ~/AgedBrainSYSBIO/data/gwas is for GWAS dataset
 ~/AgedBrainSYSBIO/data/allenbrain is for microarray datasets downloaded from Allen brain atlas
 ~/AgedBrainSYSBIO/data/ps is for positive selection data set  

## SCRIPTS

### 1. Data location
Please specify in the beginning of the script the path to the dataset(s) that you are going to process.
For easier organisation of your project you can save the uploaded data into the folder “data” with the corresponding 

### 2. Scripts related to the calculation of interactions are located in the following folders:

#### 2.2 Calculation of gene co-expression using microarray data. 
Loction: " ~/AgedBrainSYSBIO/scripts/adn” 
2.2.1 Filtering out probesets with SD < 0.29. Extracting Alzheimer’s related and healthy samples from the data sets. Calculating the co-expression between  all probesets in each of the data sets using  Spearman correlation coefficient. For each individual probeset in each of the datasets script creates 2 separate files in .txt and .RData formats. Files are named after the probeset. Created files contain the names of the correlated probesets and the corresponding Spearman coefficient.
E_GEOD_28146.R
E_GEOD_4757.R
E_MEXP_2280.R
E_GEOD_18309.R
E_GEOD_29652.R
E_GEOD_5281.R

Results are written into files for individual probesets

pathRdata <- "~/AgedBrainSYSBIO/results/adn/all_probes/rdata/E_GEOD_18309/"

pathtxt <- "~/AgedBrainSYSBIO/results/adn/all_probes/txt/E_GEOD_18309/"


2.1.2 Ranking the co-expressed values for each probeset in each of the datasets. Aggregating the ranks in all the datasets.
RRA_probesets.R

2.1.3 Converting the values to ensg ids. The script also handles the cases when affymetrix probeset id was not recognised by the gconvert() function due to the presence of the unrecognised sexual symbols in the probeset name.
affy2ensg.R

2.1.4 Assembling together calculated RRA scores for all affy probes is performed in coexp_int.R  

2.1.5 Removing “self loops”(co-expression of gene with itself) from the co-expression dataset.
coexp2undirrected_selfloops_rm.R 

2.1.6 final_adn.R  Adding columns interaction_type, data_source.
Script adds missing columns to match the integrated dataset format interaction_type="co-expression" and data_source="ADN”. The filnal co-expression part of the integrated dataset is stored in 

"~/AgedBrainSYSBIO/results/adn/integration/adn_coexp_int.RData" 

"~/AgedBrainSYSBIO/results/adn/integration/adn_coexp_int.txt" 

#### 2.2 Scripts processing epistatic interactions.
Location “~/AgedBrainSYSBIO/scripts/epistasis”

2.2.1 epi_adni.R Script processes interactions related to ADNI cohort taking to account change in ventrical volume.

2.2.2 epi_hbtrc.R Script processes interactions related to HBTRC cohort.

2.2.3 epi_tgen.R Script processes interactions related to TGEN cohort.

2.2.3 epi_adni_cog.R Script processes interactions related to ADNI cohort taking to account cognitive traits.

2.3.4 combine.R Script combines all epistatic data into one ds.

#### 2.3 Reading in and processing datasets dowloaded from IntAct database.
Location : “scripts/intact” 

#### 2.4 Reading in and processing of the protein-protein interaction dataset realted to brainageing generated by AgedBrainSYSBIO consortium.
Location: “ ~/AgedBrainSYSBIO/scripts/pba” 
 
### 3. Scripts related to the calculation of interactions are located in the following folders:

#### 3.1 Processing of GWAS data.
Location: “ ~/AgedBrainSYSBIO/scripts/gwas”

#### 3.2 Processing of the gene expression datasets from Allen Brain Atlas.
Location: “ ~/AgedBrainSYSBIO/scripts/allenbrain”

3.2.1 Scripts for the preprocessing of the individual datasets.
Read in initial data and prepare data for the analysis.Combine expression values, sample annotations and probe annotations into matrix.
"aba_178236545_ds.R", "aba_178238266_ds.R", "aba_178238316_ds.R", "aba_178238359_ds.R".

NB! Please make sure you have placed microarray datasets under "~/AgedBrainSYSBIO/data/allenbrain/" into the individual folders with the following names: /178236545_ds, /178238266_ds, /178238316_ds, /178238359_ds.

3.2.2 Calculate z scores in each of the tissue in all data sets "aba_zcores.R".

3.2.3 Assemble matrix of z-scores for all probes in the datasets from individual results for each tissue "aba_scores_assemble.R".

3.2.4 Convert probe names to Ensembl namespace "probes2ensg.R".

3.2.5 Select maximum z-scores for probes corresponding to the same Ensembl gene name "zscores_select_max.R".


#### 3.3 Reading in and processing of the dataset of positively selected genes generated by AgedBrainSYSBIO consortium.
Location:  “ ~/AgedBrainSYSBIO/scripts/ps” 

### 4. Integrate individual results for interactions and interactors.
Location : “ ~/AgedBrainSYSBIO/scripts/integration” 

4.1 All types of interactions are combines using integrate.R script.

4.2 All the nodes attributes are combined using integrate_node_attributes.R script.

## RESULTS
#### Intermediate results are saved as .RData objects and as .txt files with the corresponding names in the folder  "~/AgedBrainSYSBIO/results/".

Microarrays in .nc format "~/AgedBrainSYSBIO/results/adn/"  

Intermediate results of affy2ensg.R, E_GEOD_28146.R, E_GEOD_4757.R, E_MEXP_2280.R, E_GEOD_18309.R, E_GEOD_29652.R, E_GEOD_5281.R  "~/AgedBrainSYSBIO/results/adn/all_probes" 

Spearman correlation results for each probe in RData format  "~/AgedBrainSYSBIO/results/adn/all_probes/rdata" 

Spearman correlation results for each probe in txt format "~/AgedBrainSYSBIO/results/adn/all_probes/txt" 

Aggregated ranks for each probeset in RData format  "~/AgedBrainSYSBIO/results/adn/all_probes/scores/rdata" 

Aggregated ranks for each probeset in txt format "~/AgedBrainSYSBIO/results/adn/all_probes/scores/txt"

Epistatic datasets "~/AgedBrainSYSBIO/results/epistasis" 

PPI dataset related to brain ageing "~/AgedBrainSYSBIO/results/pba" 

PPI datasets downloaded from IntAct database "~/AgedBrainSYSBIO/results/intact"

GWAS dataset "~/AgedBrainSYSBIO/results/gwas"

Microarray datasets downloaded from Allen brain atlas "~/AgedBrainSYSBIO/results/allenbrain"

Positive selection data set "~/AgedBrainSYSBIO/results/ps" 

Integrated dataset "~/AgedBrainSYSBIO/results/integration"

