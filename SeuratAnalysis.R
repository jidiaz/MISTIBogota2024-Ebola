library(dplyr)
library(Matrix)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)
library(scDblFinder)
library(future) #paralelization

######################### Paralelization function ######################################
# check the current active plan
plan()
# change the current plan to access parallelization
plan("multicore", workers = 4)
availableCores()
plan()
########################################################################################

setwd("~/Documents/scRNA-Seq-MIT-Training/MISTIBogota2024-Ebola/Projects/Project4-OriginalMaterial/")

#Now lets load the data 
matrix <- read.csv(file = "counts.csv.gz",sep = ",",row.names = 1)
matrix[1:5,1:5]
metadata <- read.csv(file = "metadata.csv",sep = ",",row.names = 1,header = T)
head(metadata)
##Lets create the seurat object
seuratObject <- CreateSeuratObject(counts = matrix[,rownames(metadata) %in% colnames(matrix)], meta.data = metadata[rownames(metadata) %in% colnames(matrix),], project = "Project4")
seuratObject
mito.genes<-rownames(seuratObject)[rownames(seuratObject) %in% c('ENSMMUG00000028704','ENSMMUG00000028703','ENSMMUG00000028702','ENSMMUG00000028701','ENSMMUG00000028700','ND1',
                                                                 'ENSMMUG00000028698','ENSMMUG00000028697','ENSMMUG00000028696','ND2','ENSMMUG00000028694','ENSMMUG00000028693','ENSMMUG00000028692','ENSMMUG00000028691',
                                                                 'ENSMMUG00000028690','COX1','ENSMMUG00000028688','ENSMMUG00000028687','COX2','ENSMMUG00000028685','ATP8','ATP6','COX3','ENSMMUG00000028681','ND3','ENSMMUG00000028679','ND4L',
                                                                 'ND4','ENSMMUG00000028676','ENSMMUG00000028675','ENSMMUG00000028674','ND5','ND6','ENSMMUG00000028671','CYTB','ENSMMUG00000028669','ENSMMUG00000028668')]
seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, features = mito.genes)

##How many cells do we have in each sample? hbu each time point (DPI)?
metadata1 <-seuratObject@meta.data
cell_counts <- table(metadata1$orig.ident)
cell_counts_df <- as.data.frame(cell_counts)
cell_counts_df
write.csv(cell_counts_df, file = "cell_counts_per_sample.csv", row.names = FALSE)

## Lets do QC! Use the practical from the previous day and
# Violin plots
vln_plot <- VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)
vln_plot
plot1 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position = "none") + labs(x = "n_UMIs", y = "mithocondrial_Content")
plot2 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position = "right") + labs(x = "n_UMIs", y = "n_Genes")
plot1 + plot2 + plot_layout(guides = "collect")

## Let's do doublets detection
sce <- as.SingleCellExperiment(seuratObject)
sce
set.seed(123)
results <- scDblFinder(sce, returnType = 'table') %>%
  as.data.frame() %>%
  filter(type == 'real')
head(results)

results %>%
  dplyr::count(class)

outfile = file.path('Ebola_doubletFile.txt')
write.table(results, outfile, sep='\t', quote=F,
            col.names=TRUE, row.names=TRUE)

keep = results %>%
  dplyr::filter(class == "singlet") %>%
  rownames()
seuratObject = seuratObject[, keep]
seuratObject

## Let's normalize the dataset

seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

## Let's look for HVF
seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(seuratObject), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seuratObject) + theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + theme(legend.position="none")
plot1 + plot2

# Scaling the data
all.genes <- rownames(seuratObject)
seuratObject <- ScaleData(seuratObject, features = all.genes)

# Perform linear dimensional reduction
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
print(seuratObject[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seuratObject, dims = 1:2, reduction = "pca")
DimPlot(seuratObject, reduction = "pca")



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
