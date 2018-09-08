# AgedBrainSYSBIO
The scripts in this repository were used to create integrated dataset described in the article:
*Bringing together heterogeneous datasets related to Alzheimer's disease, Sügis et. al., 2018*.      
For detailed description of the methods please see the publication.

**NB! This file contains information about the structure of the files and some useful notes on how to repeat the analysis.**

## CLONING REPOSITORY
To clone repository to your local machine use the following command:  
```
git clone https://github.com/esugis/AgedBrainSYSBIO.git  
```
By default the project will be cloned to the folder "AgedBrainSYSBIO".  
The following project structure will be created:   
```
AgedBrainSYSBIO
|── README.md
|── data
│   |── adn
│   |── allenbrain
│   |── epistasis
│   |── gwas
│   |── intact
│   |── pba
│   |── ps
|── data_downloader.sh
└── scripts
|   |── libraries.R
|   |── adn
|   │   |── E_GEOD_18309.R
|   │   |── E_GEOD_28146.R
|   │   |── E_GEOD_29652.R
|   │   |── E_GEOD_4757.R
|   │   |── E_GEOD_5281.R
|   │   |── E_MEXP_2280.R
|   │   |── RRA_probesets.R
|   │   |── affy2ensg.R
|   │   |── coexp2undirrected_selfloops_rm.R
|   │   |── coexp_int.R
|   │   └── final_adn.R
|   |── epistasis
|   │   |── combine.R
|   │   |── epi_adni.R
|   │   |── epi_adni_cog.R
|   │   |── epi_hbtrc.R
|   │   └── epi_tgen.R
|   |── intact
|   │   |── intact.R
|   │   |── alz_intact.R
|   │   └── synapse_intact.R
|   |── pba
|   │   └── pba_int.R
|   |── gwas
|   │   └── gwas2ensg.R
|   |── allenbrain
|   │   |── aba_178236545_ds.R
|   │   |── aba_178238266_ds.R
|   │   |── aba_178238316_ds.R
|   │   |── aba_178238359_ds.R
|   │   |── aba_178238373_ds.R
|   │   |── aba_178238387_ds.R
|   │   |── aba_zcores.R
|   │   |── aba_zcores_assemble.R
|   │   |── expressed_regions.R
|   │   |── probes2ensg.R
|   │   └── zscores_select_max.R
|   |── ps
|   │   └── positive_selection.R
|   └── integration
|       |── integrate.R
|       └── integrate_node_attributes.R
└── results
    |── adn
    |   |── all_probes
    |   └── integration
    |── allenbrain 
    |   |── all_probes
    |   └── integration
    |── epistasis 
    |── gwas 
    |── intact 
    |── pba 
    |── ps 
    └──integration
```


## STRUCTURE

The main default catalogue of the project is ~/AgedBrainSYSBIO (where ~ stands for your home dirrectory).  
If you have cloned the repository into other directory, then project route folder is ~/DestinationDirectory.    
Project is divided into 3 major parts:  
 - **data**  a place to store downloaded datasets  
 - **scripts** a place to store analysis scripts  
 - **results** a place to store the results  

## DATA

Data is stored in the folder ~/AgedBrainSYSBIO/**data**  
Folder **data** contains sub-folders for individual data types used in the analysis.  
For automatic data download run data_downloader.sh script in your terminal:  
```
|── data_downloader.sh 
```

```
cd AgedBrainSYSBIO
./data_downloader.sh
```
All necessary datasets will be downloaded into the corresponding locations:

```
|── data
│   |── adn contains microarray datasets in .nc format
│   |── allenbrain contains microarray datasets downloaded from Allen Brain Atlas
│   |── epistasis contains epistatic datasets computed based on ADNI, TGEN and HBTRC cohorts data
│   |── gwas contains GWAS dataset
│   |── intact contains PPI datasets downloaded from IntAct database
│   |── pba contains PPI dataset related to brain ageing
│   |── ps contains positive selection dataset
```

## SCRIPTS

### General information.
Project scripts are essentially serve for three purposes.      
1. The first category of scripts are used to preprocess and analyse the datasets related to the experimental (Y2H PPIs, PPIs from IntAct) and computational (epistasis, co-expression) interactions.      
2. The second category of scripts are used to preprocess and analyse the datasets related to the genes and proteins present in the interactions from clause 1, i.e. GWAS, expression in the brain regions, etc.         
3. The third group of scripts are used to assemble (integrate) the individual resulting datasets  about the interactions from 2.1.1 and the individual genes and proteins from clause 2.      
4. Attach all R libraries used in the project.

```
|── libraries.R
```


#### 1.1 Calculation of gene co-expression using microarray data. 
1.1.1 Filtering out probesets with SD < 0.29. Extracting Alzheimer’s related and healthy samples from the data sets. Calculating the co-expression between  all probesets in each of the data sets using  Spearman correlation coefficient. For each individual probeset in each of the datasets script creates 2 separate files in .txt and .RData formats. Files are named after the probeset. Created files contain the names of the correlated probesets and the corresponding Spearman coefficient.
```
|── adn
│   |── E_GEOD_18309.R
│   |── E_GEOD_28146.R
│   |── E_GEOD_29652.R
│   |── E_GEOD_4757.R
│   |── E_GEOD_5281.R
│   |── E_MEXP_2280.R
```
Results are written into files for individual probesets
```
pathRdata <- "~/AgedBrainSYSBIO/results/adn/all_probes/rdata/E_GEOD_18309/"
pathtxt <- "~/AgedBrainSYSBIO/results/adn/all_probes/txt/E_GEOD_18309/"
```
1.1.2 Ranking the co-expressed values for each probeset in each of the datasets. Aggregating the ranks in all the datasets.
```
    │   |── RRA_probesets.R  
```
1.1.3 Converting the values to ensg ids. The script also handles the cases when affymetrix probeset id was not recognised by the gconvert() function due to the presence of the unrecognised symbols in the probeset name.
```
    │   |── affy2ensg.R  
```
1.1.4 Assembling together calculated RRA scores for all the probes.
```
    │   |── coexp_int.R 
```

1.1.5 Removing “self loops”(co-expression of gene with itself) from the co-expression dataset.
```
    │   |── coexp2undirrected_selfloops_rm.R  
```

1.1.6 Adding columns interaction_type, data_source.
```
    │   └──  final_adn.R  
```

Script adds missing columns to match the integrated dataset format interaction_type="co-expression" and data_source="ADN”. The filnal co-expression part of the integrated dataset is stored in 
```
"~/AgedBrainSYSBIO/results/adn/integration/adn_coexp_int.RData" 
"~/AgedBrainSYSBIO/results/adn/integration/adn_coexp_int.txt" 
```
#### 1.2 Scripts processing epistatic interactions.

1.2.1 epi_adni.R Script processes epistatic interactions associated with change in ventrical volume in ADNI cohort.
```
    |── epistasis             
    │   |── epi_adni.R  
```
1.2.2 Preprocessing epistatic interactions in HBTRC cohort dataset.
```
│   |── epi_hbtrc.R  
```
1.2.3 Preprocessing epistatic interactions in TGEN cohort dataset.
```
│   └── epi_tgen.R    
```

1.2.4 Preprocessing epistatic interactions associated with a set of cognitive traits in ADNI cohort dataset.
```
│   |── epi_adni_cog.R  
```

1.2.5 Coombining all epistatic data into one dataset.
```
│   |── combine.R 
```
#### 1.3 Reading in and processing datasets dowloaded from IntAct database.
1.3.1 Preprocessing all human protein-protein interactions from IntAct database. 
```
|── intact  
│   |── intact.R  
```

1.3.2 Preprocessing protein-protein interactions from IntAct database where one interacting partner is associated with Alzheimre's disease.
```
│   |── alz_intact.R   
```
1.3.3 Preprocessing synaptic protein-protein interactions from IntAct database.
```
│   └── synapse_intact.R   
```
#### 1.4 Reading in and processing of the protein-protein interaction dataset realted to brainageing generated by AgedBrainSYSBIO consortium.
```
|── pba  
│   └── pba_int.R     
``` 
### 2. Scripts related to the calculation of interactions are located in the following folders:

#### 2.1 Processing of GWAS data.   
```
|── gwas  
│   └── gwas2ensg.R      
```
#### 2.2 Analysing gene expression datasets from Allen Brain Atlas. 
2.2.1 Preprocess the individual datasets.       
```
|── allenbrain
│   |── aba_178236545_ds.R
│   |── aba_178238266_ds.R
│   |── aba_178238316_ds.R
│   |── aba_178238359_ds.R
│   |── aba_178238373_ds.R
│   |── aba_178238387_ds.R
```
2.2.2 Calculate z scores in each of the tissue in all data set.     
```
|   |── aba_zcores.R    
```
2.2.3 Assemble matrix of z-scores for all probes in the datasets from individual results for each tissue.    
```   
|   |── aba_zcores_assemble.R    
```
2.2.4 Convert probe names to Ensembl namespace.   
```
|   |── probes2ensg.R     
```
2.2.5 Select maximum z-scores for probes corresponding to the same Ensembl gene name.
```
|   └── zscores_select_max.R    
```
#### 2.3 Reading in and processing of the dataset of positively selected genes generated by AgedBrainSYSBIO consortium.
```
└── ps  
    └── positive_selection.R     
```
### 3. Integrate individual results for interactions and interactors.
```
└── integration  
```
3.1 Combine interactions.
```
    |── integrate.R  
```
3.2 Combine node (gene) attributes.
```
    └── integrate_node_attributes.R  
```
## RESULTS
4.1 The results of the individual data type analyses are save as .RData and .txt formats. They are placed into the newly created folder **AgedBrainSYSBIO/results**.      
Folder has similar structure as **AgedBrainSYSBIO/data**

```
└── results
    |── adn contains resultig co-expression interactions
    |── allenbrain contains co-expression intercation in brain regions related to Alzheimer's disease and aggregated expression in 221 brain regions
    |── epistasis contains preprocessed epistatic interactions from ADNI, TGEN and HBTRC cohorts
    |── gwas contains preprocessed data from GWAS studies
    |── intact contains preprcessed PPI from IntAct database
    |── pba contains preperocessed PPI related to brain ageing
    |── ps contains preprocessed positive selection p-values
```
4.2 The resulted integrated dataset cosists of two parts: integrated interactions and node attributes. These parts are stored in the sub-folder *integration* correspondingly in the files *integrated_int.txt* and *integrated_int_attributes.txt* . 
```
    └──integration contains integrated interactions of individual data types and node attributes
```
