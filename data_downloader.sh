#!/bin/bash

mkdir data
mkdir data/epistasis
wget https://ndownloader.figshare.com/files/9918205?private_link=65859545441cbc5b70ac -O data/epistasis/ADNI_CT_epistasis.txt
wget https://ndownloader.figshare.com/files/9918220?private_link=65b137e4198c1849177f -O data/epistasis/ADNI_VER_epistasis.tsv
wget https://ndownloader.figshare.com/files/9918214?private_link=65da7ebf68a1938b4df0 -O data/epistasis/HBTRC_epistasis.tsv
wget https://ndownloader.figshare.com/files/9918217?private_link=687b22e229f32f6aa302 -O data/epistasis/TGEN_epistasis.tsv

mkdir data/ps
wget https://ndownloader.figshare.com/files/8943232?private_link=7d26e27526f9524fc5f8 -O data/ps/positive_darwinian_selection.csv

mkdir data/gwas
wget https://ndownloader.figshare.com/files/9919579?private_link=d2e8de1569c2c953f84f -O data/gwas/IGAP_stage_1_2_combined.txt

mkdir data/pba
wget https://ndownloader.figshare.com/files/9982702?private_link=2aec9f67eb45a362c2af -O data/pba/PBA_PPI_HS.txt

mkdir data/adn
wget https://ndownloader.figshare.com/files/8753989?private_link=e469a93aad71bba4a9b5 -O data/adn/E-GEOD-18309.nc
wget https://ndownloader.figshare.com/files/8753992?private_link=e469a93aad71bba4a9b5 -O data/adn/E-GEOD-28146.nc
wget https://ndownloader.figshare.com/files/8753983?private_link=e469a93aad71bba4a9b5 -O data/adn/E-GEOD-4757.nc
wget https://ndownloader.figshare.com/files/8753977?private_link=e469a93aad71bba4a9b5 -O data/adn/E-MEXP-2280.nc
wget https://ndownloader.figshare.com/files/8753986?private_link=e469a93aad71bba4a9b5 -O data/adn/E-GEOD-29652.nc
wget https://ndownloader.figshare.com/files/8753980?private_link=e469a93aad71bba4a9b5 -O data/adn/E-GEOD-5281.nc

mkdir data/intact
wget https://ndownloader.figshare.com/files/9924304?private_link=56a99af6889415a71264 -O data/intact/intact_hs_v_4_2_6.txt
wget https://ndownloader.figshare.com/files/9924376?private_link=4c9b36ad30e4919fca5d -O data/intact/alzheimers_intact_v_4_2_6.txt
wget https://ndownloader.figshare.com/files/9924400?private_link=7dd03a13f95a8991eaa6 -O data/intact/synapse_intact_v_4_2_6.txt

mkdir data/allenbrain
mkdir data/allenbrain/178236545_ds
mkdir data/allenbrain/178238266_ds
mkdir data/allenbrain/178238316_ds
mkdir data/allenbrain/178238359_ds
mkdir data/allenbrain/178238387_ds
mkdir data/allenbrain/178238373_ds

wget http://human.brain-map.org/api/v2/well_known_file_download/178236545 -O data/allenbrain/normalized_microarray_donor15697.zip
unzip data/allenbrain/normalized_microarray_donor15697.zip -d data/allenbrain/178236545_ds
wget http://human.brain-map.org/api/v2/well_known_file_download/178238359 -O data/allenbrain/normalized_microarray_donor12876.zip
unzip data/allenbrain/normalized_microarray_donor12876.zip -d data/allenbrain/178238359_ds
wget http://human.brain-map.org/api/v2/well_known_file_download/178238316 -O data/allenbrain/normalized_microarray_donor14380.zip
unzip data/allenbrain/normalized_microarray_donor14380.zip -d data/allenbrain/178238316_ds 
wget http://human.brain-map.org/api/v2/well_known_file_download/178238266 -O data/allenbrain/normalized_microarray_donor15496.zip
unzip data/allenbrain/normalized_microarray_donor15496.zip -d data/allenbrain/178238266_ds
wget http://human.brain-map.org/api/v2/well_known_file_download/178238387 -O normalized_microarray_donor9861.zip
unzip data/allenbrain/normalized_microarray_donor9861.zip -d data/allenbrain/178238387_ds
wget http://human.brain-map.org/api/v2/well_known_file_download/178238373 -O normalized_microarray_donor10021.zip
unzip data/allenbrain/normalized_microarray_donor10021.zip -d data/allenbrain/178238373_ds
