## HCA Africa: Single-Cell RNA-Seq Computational Workshop, March 7, 2023
library(Seurat)
###Project 4 Ebola
setwd("Computational_Projects/Project4")

##Now lets load the mouse data 
matrix <- read.csv(file = "counts.csv.gz",sep = ",",row.names = 1)
matrix[1:5,1:5]
metadata <- read.csv(file = "metadata.csv",sep = ",",row.names = 1,header = T)
head(metadata)
##Lets create the seurat object for Mouse
seuratObject <- CreateSeuratObject(counts = matrix[,rownames(metadata) %in% colnames(matrix)], meta.data = metadata[rownames(metadata) %in% colnames(matrix),], project = "Project4")
seuratObject
mito.genes<-rownames(seuratObject)[rownames(seuratObject) %in% c('ENSMMUG00000028704','ENSMMUG00000028703','ENSMMUG00000028702','ENSMMUG00000028701','ENSMMUG00000028700','ND1',
  'ENSMMUG00000028698','ENSMMUG00000028697','ENSMMUG00000028696','ND2','ENSMMUG00000028694','ENSMMUG00000028693','ENSMMUG00000028692','ENSMMUG00000028691',
  'ENSMMUG00000028690','COX1','ENSMMUG00000028688','ENSMMUG00000028687','COX2','ENSMMUG00000028685','ATP8','ATP6','COX3','ENSMMUG00000028681','ND3','ENSMMUG00000028679','ND4L',
  'ND4','ENSMMUG00000028676','ENSMMUG00000028675','ENSMMUG00000028674','ND5','ND6','ENSMMUG00000028671','CYTB','ENSMMUG00000028669','ENSMMUG00000028668')]
seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, features = mito.genes)

##How many cells do we have in each sample? hbu each time point (DPI)?

## Lets do QC! Use the practical from the previous day and
#check the seurat tutorial! https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

## Once We QC and filter lets do normalization, scaling , PCA, clustering and umap!

## Identify Clusters! Hint: Check for  immunes cells and intestinal cells for example B cells,
#Dendritic Cells, , Macrophages,Platelets, Monocytes , Neutrophils , NK Cells and T Cells
#Use websites like https://singlecell.broadinstitute.org/single_cell and
#https://panglaodb.se/ (if cant connect use an online proxy) as well a protein tissue atlas 

## How many cell do we have in every condition and sample? Is there any batch effect? Do we need to integrate? 
#Tip use the practicals and the seurat tutorial https://satijalab.org/seurat/articles/integration_introduction.html

##How about per cluster/Cell type

## Where can you identify infected cells, check the ebola genes
## "EBOV-GENOME" "EBOV-GP"     "EBOV-L"      "EBOV-NP"     "EBOV-VP24"   "EBOV-VP30"   "EBOV-VP35"   "EBOV-VP40" 

## How do cytokines and ISGs change across time
##Check genes like "STAT1", "ISG15, "MX1" among others! 

##Lets play with Differential exppression

## Lets make pretty plots! 
#Check visualization from Seurat https://satijalab.org/seurat/articles/visualization_vignette.html and GGPlot2
