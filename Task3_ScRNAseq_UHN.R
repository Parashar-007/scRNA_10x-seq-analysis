# Title: Single-cell RNA  analysis on RNA Counts data
# File: Task2_DeSeq2_UHN
# Project: UHN_Assignment

# Working Dir: setwd("C:/Users/surab/Desktop/UHN_Tasks")

# Install required packages
install.packages("Seurat")

# Load libraries
library(Seurat)
library(tidyverse)

# Load the ScRNA seq dataset
hbm.sc10x <- 
  Read10X(data.dir = 
            "C:/Users/surab/Desktop/UHN_Tasks/scRNA10x-20240221T021714Z-001/scRNA10x/frozen_bmmc_healthy_donor1_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/")
str(hbm.sc10x)


# Initialize the Seurat object with the raw (non-normalized data).
hbm.sc10x.seurat.obj <- CreateSeuratObject(counts = hbm.sc10x, project = "ScRNA", min.cells = 3, min.features = 200)
str(hbm.sc10x.seurat.obj)


# 1. QC -------
View(hbm.sc10x.seurat.obj@meta.data)
# % MT reads
hbm.sc10x.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(hbm.sc10x.seurat.obj, pattern = "^MT-")
View(hbm.sc10x.seurat.obj@meta.data)

VlnPlot(hbm.sc10x.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(hbm.sc10x.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
hbm.sc10x.seurat.obj <- subset(hbm.sc10x.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)

# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
hbm.sc10x.seurat.obj <- NormalizeData(hbm.sc10x.seurat.obj)
str(hbm.sc10x.seurat.obj)

# 4. Identify highly variable features --------------
hbm.sc10x.seurat.obj <- FindVariableFeatures(hbm.sc10x.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hbm.sc10x.seurat.obj), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hbm.sc10x.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling -------------
all.genes <- rownames(hbm.sc10x.seurat.obj)
hbm.sc10x.seurat.obj <- ScaleData(hbm.sc10x.seurat.obj, features = all.genes)

str(hbm.sc10x.seurat.obj)

# 6. Perform Linear dimensionality reduction --------------
hbm.sc10x.seurat.obj <- RunPCA(hbm.sc10x.seurat.obj, features = VariableFeatures(object = hbm.sc10x.seurat.obj))

# visualize PCA results
print(hbm.sc10x.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(hbm.sc10x.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

# determine dimensionality of the data
ElbowPlot(hbm.sc10x.seurat.obj)

# 7. Clustering ------------
hbm.sc10x.seurat.obj <- FindNeighbors(hbm.sc10x.seurat.obj, dims = 1:15)

# understanding resolution
hbm.sc10x.seurat.obj <- FindClusters(hbm.sc10x.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(hbm.sc10x.seurat.obj@meta.data)

DimPlot(hbm.sc10x.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(hbm.sc10x.seurat.obj)
Idents(hbm.sc10x.seurat.obj) <- "RNA_snn_res.0.1"
Idents(hbm.sc10x.seurat.obj)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
hbm.sc10x.seurat.obj <- RunUMAP(hbm.sc10x.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(hbm.sc10x.seurat.obj, reduction = "umap")

save.image(file = "ScRNA.10x.RData")
load("ScRNA.10x.RData")

# Single cell Annotation
# install packages

BiocManager::install("celldex")
BiocManager::install("SingleR")
install.packages("pheatmap")

# Load libraries
library(SingleR)
library(celldex)
library(pheatmap)
library(Seurat)
library(tidyverse)

# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

# expression values are log counts (log normalized counts)


# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

hbm.sc10x.seurat.obj
hbm_counts <- GetAssayData(hbm.sc10x.seurat.obj, layer = 'counts')

pred <- SingleR(test = hbm_counts,
                ref = ref,
                labels = ref$label.main)

pred

hbm.sc10x.seurat.obj$singleR.labels <- pred$labels[match(rownames(hbm.sc10x.seurat.obj@meta.data), rownames(pred))]
DimPlot(hbm.sc10x.seurat.obj, reduction = 'umap', group.by = 'singleR.labels', label=TRUE)


# Annotation diagnostics ----------
install.packages("viridis")
library(viridis)
library(viridisLite)
# ...Based on the scores within cells -----------
pred
pred$scores

plotScoreHeatmap(pred)


# ...Based on deltas across cells ----------

plotDeltaDistribution(pred)

save.image(file = "ScRNA.10x.RData")
load("ScRNA.10x.RData")
