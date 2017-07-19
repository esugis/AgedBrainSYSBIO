# This script calculate the relative expression score of each probe.
# We have 4 brain expression dataset. 
# To each datasets corresponds a matrix of 0/1 values indicating if the probe is expressed or not in the specific tissue.
# We aggregate the values of all 4 datasets by calculating a relative expression score.
# A relative expression score indicates in how many out of 4 datasets the probe is expressed.
# example 4/4 gets weight 1, 3/4 gets weight 0.75


# Load libraries
library(foreach)

 
#Read data for each of the dataset and check the details

# 178236545_ds
onto_178236545_ds <- read.csv(file="/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178236545_ds/Ontology.csv",stringsAsFactors=F)
dim(onto_178236545_ds)
#[1] 1839    8
PACall_178236545_ds <- read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178236545_ds/PACall.csv", header=F, row.names = 1,stringsAsFactors=F)
dim(PACall_178236545_ds)
#[1] 58692   501
probes_178236545_ds <- read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178236545_ds/Probes.csv",stringsAsFactors=F)
dim(probes_178236545_ds)
#[1] 58692     7
sampleann_178236545_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178236545_ds/SampleAnnot.csv",stringsAsFactors=F)
dim(sampleann_178236545_ds)
#[1] 501  13


#178238266_ds
onto_178238266_ds=read.csv(file="/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238266_ds/Ontology.csv",stringsAsFactors=F)
dim(onto_178238266_ds)
#[1] 1839    8
PACall_178238266_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238266_ds/PACall.csv", header=F, row.names = 1,stringsAsFactors=F)
#first column has probes' numbers
dim(PACall_178238266_ds)
#[1] 58692   470
probes_178238266_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238266_ds/Probes.csv",stringsAsFactors=F)
dim(probes_178238266_ds)
#[1] 58692     7
sampleann_178238266_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238266_ds/SampleAnnot.csv",stringsAsFactors=F)
dim(sampleann_178238266_ds)
#[1] 470  13


onto_178238316_ds=read.csv(file="/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238316_ds/Ontology.csv",stringsAsFactors=F)
dim(onto_178238316_ds)
#[1] 1839    8
PACall_178238316_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238316_ds/PACall.csv", header=F,row.names = 1,stringsAsFactors=F)
dim(PACall_178238316_ds)
#[1] 58692   529
probes_178238316_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238316_ds/Probes.csv",stringsAsFactors=F)
dim(probes_178238316_ds)
#[1] 58692     7
sampleann_178238316_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238316_ds/SampleAnnot.csv",stringsAsFactors=F)
dim(sampleann_178238316_ds)
#[1] 529  13


onto_178238359_ds=read.csv(file="/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238359_ds/Ontology.csv",stringsAsFactors=F)
dim(onto_178238359_ds)
#[1] 1839    8
PACall_178238359_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238359_ds/PACall.csv", header=F,row.names = 1,stringsAsFactors=F)
dim(PACall_178238359_ds)
#[1] 58692   363
probes_178238359_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238359_ds/Probes.csv",stringsAsFactors=F)
dim(probes_178238359_ds)
#[1] 58692     7
sampleann_178238359_ds=read.csv("/home/nikolaeva/AgedBrain/ppi2_data/allenbrainatlas/178238359_ds/SampleAnnot.csv",stringsAsFactors=F)
dim(sampleann_178238359_ds)
#[1] 363  13


# Attach tissue IDs as column names
colnames(PACall_178238359_ds)<-sampleann_178238359_ds$structure_id
colnames(PACall_178238316_ds)<-sampleann_178238316_ds$structure_id
colnames(PACall_178238266_ds)<-sampleann_178238266_ds$structure_id
colnames(PACall_178236545_ds)<-sampleann_178236545_ds$structure_id

# Find tissues, that are present in all datasets
tissues_all<-unique(Reduce(intersect, 
list(colnames(PACall_178238359_ds),colnames(PACall_178238316_ds),colnames(PACall_178238266_ds),colnames(PACall_178236545_ds))))

# For each probe find the weighted sum of scores(1s and 0s)
# Add column probe_id for each   

PACall_178238359_ds<-cbind(PACall_178238359_ds, probe_id=rownames(PACall_178238359_ds))
PACall_178238316_ds<-cbind(PACall_178238316_ds, probe_id=rownames(PACall_178238316_ds))
PACall_178238266_ds<-cbind(PACall_178238266_ds, probe_id=rownames(PACall_178238266_ds))
PACall_178236545_ds<-cbind(PACall_178236545_ds, probe_id=rownames(PACall_178236545_ds))



binary_expression=foreach(i = 1:length(tissues_all), .combine = rbind)%do%{
#subset for one tissue
tissue=tissues_all[i]
print("current tissue")
print(tissue)
#subset tissue related data from each ds
tissue1=subset(PACall_178238359_ds, colname(PACall_178238359_ds)%in%tissue)
onetissue2=subset(PACall_178238316_ds, structure_id==tissue)
onetissue3=subset(PACall_178238266_ds, structure_id==tissue)
onetissue4=subset(PACall_178236545_ds, structure_id==tissue)
onetissue=rbind(onetissue1,onetissue2,onetissue3,onetissue4)
onetissue_av <- merge(aggregate(value ~ probe_id, onetissue, mean),tissue)
#calculate z-score over one tissue
#onetissue_av=cbind(onetissue_av, zscore=((onetissue_av[,2] - mean(onetissue_av[,2])) / sd(onetissue_av[,2])))
#write file for each tissue
filename=sprintf("%s.txt",tissue);
filedata=sprintf("%s.RData",tissue);
pathname= file.path(pathtxt, filename);
pathdata=file.path(pathRdata, filedata);
write.table(onetissue_av, file=pathname, sep="\t", quote=F);
save(onetissue_av,file=pathdata)
}
save(binary_expression,file="binary_expression.RData")
head(tissues_z[,1:6])
write.table(binary_expression, file="binary_expression.txt", sep="\t", quote=F);
