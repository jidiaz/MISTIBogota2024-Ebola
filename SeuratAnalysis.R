library(dplyr)
library(Matrix)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)
##Libraries dublets
library(scDblFinder)
## Libraries integration and batch correction
library(harmony)
##Libraries cell annotation
library(celldex)
library(SingleR)
library(pheatmap)

library(future) #paralelization

######################### Paralelization function ######################################
# check the current active plan
plan()
# change the current plan to access parallelization
plan("multicore", workers = 4)
availableCores()
plan()
########################################################################################

setwd("~/MISTIBogota2024-Ebola/Projects/Project4-OriginalMaterial/")

#Now lets load the data 
matrix <- read.csv(file = "counts.csv.gz",sep = ",",row.names = 1)
matrix[1:5,1:5]
metadata <- read.csv(file = "metadata.csv",sep = ",",row.names = 1,header = T)
rownames(metadata)<-gsub("-",".",rownames(metadata))
head(metadata)
##Lets create the seurat object
SO <- CreateSeuratObject(counts = matrix[,rownames(metadata) %in% colnames(matrix)], meta.data = metadata[rownames(metadata) %in% colnames(matrix),], project = "Project4")
mito.genes<-rownames(SO)[rownames(SO) %in% c('ENSMMUG00000028704','ENSMMUG00000028703','ENSMMUG00000028702','ENSMMUG00000028701','ENSMMUG00000028700','ND1',
                                                                 'ENSMMUG00000028698','ENSMMUG00000028697','ENSMMUG00000028696','ND2','ENSMMUG00000028694','ENSMMUG00000028693','ENSMMUG00000028692','ENSMMUG00000028691',
                                                                 'ENSMMUG00000028690','COX1','ENSMMUG00000028688','ENSMMUG00000028687','COX2','ENSMMUG00000028685','ATP8','ATP6','COX3','ENSMMUG00000028681','ND3','ENSMMUG00000028679','ND4L',
                                                                 'ND4','ENSMMUG00000028676','ENSMMUG00000028675','ENSMMUG00000028674','ND5','ND6','ENSMMUG00000028671','CYTB','ENSMMUG00000028669','ENSMMUG00000028668')]
SO[["percent.mt"]] <- PercentageFeatureSet(SO, features = mito.genes)
SO

##How many cells do we have in each sample? hbu each time point (DPI)?
metadata1 <-SO@meta.data
cell_counts <- table(metadata1$orig.ident)
cell_counts_df <- as.data.frame(cell_counts)
cell_counts
write.csv(cell_counts_df, file = "cell_counts_per_sample.csv", row.names = FALSE)

cell_countsDPI <- table(metadata1$DPI)
cell_countsDPI_df <- as.data.frame(cell_countsDPI)
cell_countsDPI
write.csv(cell_countsDPI_df, file = "cell_counts_per_DPI.csv", row.names = FALSE)


## Lets do QC!
vln_plot_feat_animal <- VlnPlot(SO, features = "nFeature_RNA", pt.size = 0.0001, group.by = "animal") +
  theme(legend.position = "none")
vln_plot_count_animal <- VlnPlot(SO, features = "nCount_RNA", pt.size = 0.0001, group.by = "animal") +
  theme(legend.position = "none")
vln_plot_mt_animal <- VlnPlot(SO, features = "percent.mt", pt.size = 0.0001, group.by = "animal") +
  theme(legend.position = "none")
vln_plot_feat_DPI <- VlnPlot(SO, features = "nFeature_RNA", pt.size = 0.0001, group.by = "DPI") +
  theme(legend.position = "none")
vln_plot_count_DPI <- VlnPlot(SO, features = "nCount_RNA", pt.size = 0.0001, group.by = "DPI") +
  theme(legend.position = "none")
vln_plot_mt_DPI <- VlnPlot(SO, features = "percent.mt", pt.size = 0.0001, group.by = "DPI") +
  theme(legend.position = "none")

ggsave("01A.QC_filtered.png", plot = vln_plot_feat_animal, width = 24, height = 12, dpi = 300)
ggsave("01B.QC_filtered.png", plot = vln_plot_count_animal, width = 12, height = 4, dpi = 300)
ggsave("01C.QC_filtered.png", plot = vln_plot_mt_animal, width = 12, height = 4, dpi = 300)
ggsave("01D.QC_filtered.png", plot = vln_plot_feat_DPI, width = 24, height = 12, dpi = 300)
ggsave("01E.QC_filtered.png", plot = vln_plot_count_DPI, width = 24, height = 12, dpi = 300)
ggsave("01F.QC_filtered.png", plot = vln_plot_mt_DPI, width = 24, height = 12, dpi = 300)

plot1 <- FeatureScatter(SO, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position = "none") + labs(x = "n_UMIs", y = "mithocondrial_Content")
plot2 <- FeatureScatter(SO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position = "none") + labs(x = "n_UMIs", y = "n_Genes")
combined_plot <- plot1 + plot2 + plot_layout(guides = "collect")
ggsave("02.feature_scatter_high_res.png", plot = combined_plot, width = 12, height = 6, dpi = 300)

## Let's do doublets detection
sce <- as.SingleCellExperiment(SO)
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
SO= SO[, keep]
SO

##How many cells do we have in each sample? hbu each time point (DPI)? AFTER removing duplets
metadata2 <-SO@meta.data
cell_counts <- table(metadata2$orig.ident)
cell_counts_df <- as.data.frame(cell_counts)
write.csv(cell_counts_df, file = "cell_counts_per_sample_nodups.csv", row.names = FALSE)

cell_countsDPI <- table(metadata2$DPI)
cell_countsDPI_df <- as.data.frame(cell_countsDPI)
write.csv(cell_countsDPI_df, file = "cell_counts_per_DPI_nodups.csv", row.names = FALSE)

## Let's normalize the dataset
SO <- NormalizeData(SO, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

## Let's look for HVF
SO <- FindVariableFeatures(SO, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(SO), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SO) + theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + theme(legend.position="none")
HighVariableFeatures <- plot1 + plot2
ggsave("03.Scatter_HVF.png", plot = HighVariableFeatures, width = 12, height = 6, dpi = 300)

# Scaling the data
all.genes <- rownames(SO)
SO <- ScaleData(SO, features = all.genes)

# Perform linear dimensional reduction
SO <- RunPCA(SO, features = VariableFeatures(object = SO))
print(SO[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SO, dims = 1:2, reduction = "pca")
PCA <- DimPlot(SO, reduction = "pca")
ggsave("04.PCA.png", plot = PCA, width = 12, height = 6, dpi = 300)

pca = SO[["pca"]]

## get the eigenvalues
evs = pca@stdev^2
total.var = pca@misc$total.variance
varExplained = evs/total.var
pca.data = data.frame(PC=factor(1:length(evs)),
                      percVar=varExplained*100)
pca.data$cumulVar = cumsum(pca.data$percVar)

head(pca.data, 20)

scPlot <- pca.data[1:50,] %>%
  ggplot(aes(x=PC, y=percVar)) +
  geom_bar(stat='identity') +
  geom_hline(yintercept = 1, colour="red", linetype=3) +
  labs(title="Variance Explanation by PCA") +
  xlab("Principal Components") +
  ylab("Percentage of Explained Variance") +
  theme_bw()
scPlot
ggsave("04A.PCA_geom_bar.png",plot = scPlot, bg = 'white')

scPlot <- pca.data[1:50,] %>%
  ggplot(aes(x=PC, y=cumulVar)) +
  geom_bar(stat='identity') +
  geom_hline(yintercept = 50, colour="red", linetype=3) +
  labs(title="Cumulative Variance Explanation by PCA") +
  xlab("Principal Components") +
  ylab("Cumulative Percentage of Explained Variance") +
  theme_bw()
scPlot
ggsave("04B.PCA_geom_bar.png",plot = scPlot, bg = 'white')

# Cluster the cells
SO <- FindNeighbors(SO,  dims = 1:20)
SO <- FindClusters(SO, resolution=0.1)
head(Idents(SO), 5)
identities <- Idents(SO)
clusterorder <-SO$seurat_clusters
write.csv(identities, file = "Cluster_identities.csv")

# Run non-linear dimensional reduction (UMAP)
SO <- RunUMAP(SO, dims = 1:20)
umap_plot <- DimPlot(SO, reduction = "umap")
ggsave("05A.umap_plot_high_res.png", plot = umap_plot, width = 10, height = 8, dpi = 300)

# Run non-linear dimensional reduction (tSNE)
SO <- RunTSNE(SO, dims = 1:20)
tsne_plot <- DimPlot(SO, reduction = "tsne")
ggsave("05B.tsne_plot_high_res.png", plot = tsne_plot, width = 10, height = 8, dpi = 300)

#Plot clusters against metadata for batch errors
umap_plot_dpi <- DimPlot(SO, reduction = "umap",group.by = "DPI")
umap_plot_cell <- DimPlot(SO, reduction = "umap",group.by = "cell")
umap_plot_animal <- DimPlot(SO, reduction = "umap",group.by = "animal")
umap_plot_stype <- DimPlot(SO, reduction = "umap",group.by = "sample_type")
ggsave("06A.umap_plot_high_res_dpi.png", plot = umap_plot_dpi, width = 10, height = 8, dpi = 300)
ggsave("06B.umap_plot_high_res_cell.png", plot = umap_plot_cell, width = 10, height = 8, dpi = 300)
ggsave("06C.umap_plot_high_res_animal.png", plot = umap_plot_animal, width = 10, height = 8, dpi = 300)
ggsave("06D.umap_plot_high_res_stype.png", plot = umap_plot_stype, width = 10, height = 8, dpi = 300)

# Finding differentially expressed features (cluster markers)
levels(SO)
SO.markers <- FindAllMarkers(SO, only.pos = TRUE, min.pct = 0.25, test.use="negbinom", slot="counts")
#seuratObject.markers <- FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.25, test.use="wilcox")

SO.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
filtered_markers <- SO.markers[SO.markers$avg_log2FC > 2 & SO.markers$p_val_adj < 0.05, ]
SO.markers
write.csv(filtered_markers, file = "Significative_markers.csv", row.names = FALSE)

SO.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10

# Generate the heatmap for raw counts
heatmap_raw_counts <- DoHeatmap(SO, features = top10$gene, slot = "counts")
ggsave("07A.heatmap_raw_counts_high_res01.png", plot = heatmap_raw_counts, width = 10, height = 8, dpi = 300)

# Generate the heatmap for normalized expression
heatmap_normalized_expression <- DoHeatmap(SO, features = top10$gene, slot = "data")
ggsave("07B.heatmap_normalized_expression_high_res01.png", plot = heatmap_normalized_expression, width = 10, height = 8, dpi = 300)

#Cell annotation
ref = BlueprintEncodeData()
ref

sce = as.SingleCellExperiment(SO)
pred = SingleR(sce, ref=ref, labels=ref$label.main)
table(pred$labels)
cluster_markers <- print(top10, n = Inf)

table_celltypes_bycluster = table(Assigned=pred$pruned.labels,
                 cluster=sce$seurat_clusters)
table_celltypes_bycluster

pheatmap(log2(table_celltypes_bycluster + 1), filename = "08.CellAnnotation_pheatmap.png")

write.csv(cluster_markers, file = "Top10Markers_perCluster.csv", row.names = FALSE)

scPlot <- FeaturePlot(SO, features=c("COL5A2"))


#Batch correction
options(repr.plot.height = 4, repr.plot.width = 6)
SOharmony <- SO %>%
  RunHarmony("cell", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(SOharmony, 'harmony')
harmony_embeddings[1:5, 1:5]

#Clustering 
SOharmony <- RunUMAP(SOharmony, reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

#Plot clusters against metadata after batch correction
umap_plot_dpi <- DimPlot(SOharmony, reduction = "umap",group.by = "DPI")
umap_plot_cell <- DimPlot(SOharmony, reduction = "umap",group.by = "cell")
umap_plot_animal <- DimPlot(SOharmony, reduction = "umap",group.by = "animal")
umap_plot_stype <- DimPlot(SOharmony, reduction = "umap",group.by = "sample_type")
ggsave("07A.umap_plot_high_res_dpi.png", plot = umap_plot_dpi, width = 10, height = 8, dpi = 300)
ggsave("07B.umap_plot_high_res_cell.png", plot = umap_plot_cell, width = 10, height = 8, dpi = 300)
ggsave("07C.umap_plot_high_res_animal.png", plot = umap_plot_animal, width = 10, height = 8, dpi = 300)
ggsave("07D.umap_plot_high_res_stype.png", plot = umap_plot_stype, width = 10, height = 8, dpi = 300)

#Integration metrics
plot_integrated_clusters = function (srat, batchcolumn) {
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  
  count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data[[batchcolumn]])
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)
  
  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "dataset"
  
  
  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") +
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")
  
  p2 + p1 + plot_layout(widths = c(3,1))
}

plot_integrated_clusters(SOharmony, 'cell')

# Finding differentially expressed features (cluster markers)
levels(SOharmony)
seuratObject.markers <- FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.25, test.use="negbinom", slot="counts")
#seuratObject.markers <- FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.25, test.use="wilcox")
seuratObject.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
filtered_markers <- seuratObject.markers[seuratObject.markers$avg_log2FC > 2 & seuratObject.markers$p_val_adj < 0.05, ]
seuratObject.markers
write.csv(filtered_markers, file = "Significative_markers.csv", row.names = FALSE)

seuratObject.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10

# Generate the heatmap for raw counts
heatmap_raw_counts <- DoHeatmap(seuratObject, features = top10$gene, slot = "counts")
ggsave("heatmap_raw_counts_high_res.png", plot = heatmap_raw_counts, width = 10, height = 8, dpi = 300)

# Generate the heatmap for normalized expression
heatmap_normalized_expression <- DoHeatmap(seuratObject, features = top10$gene, slot = "data")
ggsave("heatmap_normalized_expression_high_res.png", plot = heatmap_normalized_expression, width = 10, height = 8, dpi = 300)


cluster_markers <- print(top10, n = Inf)

write.csv(cluster_markers, file = "Top10Markers_perCluster.csv", row.names = FALSE)


################### SAVE RDS OBJECT ###########################################
saveRDS(seuratObject, file = "seurat.RDS")
################### READ RDS OBJECT ###########################################
#seuratObject <- readRDS(file = "seuratObject.RDS")
###############################################################################

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
