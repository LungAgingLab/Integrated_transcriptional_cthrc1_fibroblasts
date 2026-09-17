############################################################################
# Project: Mayo (WL tissue and cultured fibroblasts)
# Purpose: Address reviewer comments for Mayo paper
# Description: Fine tuning differential expression between adventitial fibroblasts and inflammatory fibroblasts
# Author: Mahnoor Hamid
# Date: 2026-07-15
############################################################################


# ==============================================================================
# 1. ENVIRONMENT SETUP & LIBRARIES (ESSENTIALS ONLY)
# ==============================================================================
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(pdftools)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggpubr)
library(readxl)
library(dplyr)
library(stringr)
library(qs2)
library(UCell)
library(rstatix)

# Directories
datadir   <- "/fs/ess/PAS2148/Projects_Mora_Rojas_Lab/Mayo/"
plotdir   <- "/fs/ess/PAS2148/Projects_Mora_Rojas_Lab/Mayo/Plots/"

#The lung tissue data path: /fs/ess/PAS2148/huchen/fibroblasts_data/WL_fibroblasts.h5ad
#Cultured data path: /fs/ess/PAS2148/huchen/fibroblasts_data/Cultured_fibroblasts.h5ad

# 1. Install the matrix reading wrapper package if you don't have it
if (!requireNamespace("anndata", quietly = TRUE)) install.packages("anndata")

library(anndata)
library(Seurat)

# 2. Define the explicit cluster data paths
lung_tissue_path <- "/fs/ess/PAS2148/huchen/fibroblasts_data/WL_fibroblasts.h5ad"
cultured_path    <- "/fs/ess/PAS2148/huchen/fibroblasts_data/Cultured_fibroblasts.h5ad"

# 3. Read files directly into memory space
adata_lung     <- read_h5ad(lung_tissue_path)
adata_cultured <- read_h5ad(cultured_path)
# ==============================================================================
# MANUALLY PARSE LUNG TISSUE DATA
# ==============================================================================
# Extract the sparse expression matrix and transpose it (Seurat requires Features x Cells)
lung_counts <- t(adata_lung$X)

# Assign feature (gene) names and cell barcodes cleanly
rownames(lung_counts) <- adata_lung$var_names
colnames(lung_counts) <- adata_lung$obs_names

# Extract cell-level metadata as a clean dataframe
lung_metadata <- as.data.frame(adata_lung$obs)

# Create native Seurat Object
seurat_lung <- CreateSeuratObject(counts = lung_counts, meta.data = lung_metadata)

# ==============================================================================
# MANUALLY PARSE CULTURED DATA
# ==============================================================================
# Repeat same steps for the cultured fibroblast tracking matrix
cultured_counts <- t(adata_cultured$X)
rownames(cultured_counts) <- adata_cultured$var_names
colnames(cultured_counts) <- adata_cultured$obs_names
cultured_metadata <- as.data.frame(adata_cultured$obs)

seurat_cultured <- CreateSeuratObject(counts = cultured_counts, meta.data = cultured_metadata)

# ==============================================================================
# VERIFY SUCCESS
# ==============================================================================
print(seurat_lung)
print(seurat_cultured)
message("Both datasets successfully built as native Seurat objects!")


head(seurat_lung[[]])
head(seurat_cultured[[]])

# 1. Update identity for the Lung Tissue object
# (Replace "annotation" with the exact metadata column name if different)
Idents(seurat_lung) <- "annotation"

# 2. Update identity for the Cultured Fibroblast object
Idents(seurat_cultured) <- "annotation"

# 3. Quickly verify that the active identities updated successfully
print(table(Idents(seurat_lung)))
print(table(Idents(seurat_cultured)))

print(unique(seurat_lung$Sample_ID))
print(unique(seurat_cultured$Sample_ID))


#okay, using whole long downstream

names(seurat_lung@reductions)


# Set random seed to ensure exact mathematical reproducibility 
set.seed(42)

# ==============================================================================
# 1. SCANPY-STYLE REGULARIZED LOG2 NORMALIZATION
# ==============================================================================
message("Performing library size normalization to 10k...")
# Scanpy's sc.pp.normalize_total(target_sum=1e4)
seurat_lung <- NormalizeData(
  seurat_lung, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

# CRITICAL CONVERSION: Convert Seurat's natural log (ln) to Scanpy's base-2 log (log2)
# Formula derivation: log2(x) = ln(x) / ln(2)
seurat_lung@assays$RNA$data <- seurat_lung@assays$RNA$data / log(2)

# ==============================================================================
# 2. HIGHLY VARIABLE FEATURE SELECTION (3,000 GENES)
# ==============================================================================
message("Identifying top 3000 highly variable genes...")
seurat_lung <- FindVariableFeatures(
  seurat_lung, 
  selection.method = "vst", 
  nfeatures = 3000
)

# ==============================================================================
# 3. REGRESSION & CAP MAX MATRIX SCALING
# ==============================================================================
message("Regressing tech variation & capping features at 10...")
# Regresses out total UMI counts (nCount_RNA) and percentage of mitochondrial counts (pct_counts_mt)
# Scanpy's sc.pp.scale(max_value=10) is replicated by setting clip.max = 10
seurat_lung <- ScaleData(
  seurat_lung, 
  features = VariableFeatures(seurat_lung),
  vars.to.regress = c("nCount_RNA", "pct_counts_mt"),
  clip.max = 10
)

# ==============================================================================
# 4. PRINCIPAL COMPONENT ANALYSIS (PCA)
# ==============================================================================
message("Running PCA using arpack approximation solver...")
# Seurat uses 'irlba' by default, which is equivalent to the 'arpack' solver behavior
seurat_lung <- RunPCA(
  seurat_lung, 
  npcs = 30, 
  approx = TRUE, 
  verbose = FALSE
)

# ==============================================================================
# 5. BATCH EFFECT CORRECTION (HARMONY)
# ==============================================================================
message("Applying Harmony batch correction over Sample groups...")

# Explicitly calling the package via harmony:: bypasses the Rcpp namespace lock
seurat_lung <- harmony::RunHarmony(
  object = seurat_lung,
  group.by.vars = "Sample", 
  reduction.use = "pca",
  dims.use = 1:30,
  verbose = FALSE
)


# ==============================================================================
# 6. NEIGHBORHOOD GRAPH CONSTRUCTION & UMAP (LUNG CONFIG: 10 PCs, 20 NEIGHBORS)
# ==============================================================================
message("Constructing neighbor maps and UMAP embedding coordinates...")
# Whole lung specifications: 10 corrected PCs, 20 neighbors
seurat_lung <- FindNeighbors(
  seurat_lung, 
  reduction = "harmony", 
  dims = 1:10, 
  k.param = 20, 
  verbose = FALSE
)

seurat_lung <- RunUMAP(
  seurat_lung, 
  reduction = "harmony", 
  dims = 1:10, 
  verbose = FALSE
)

# ==============================================================================
# 7. GENERATE REPLICATED VISUALIZATION
# ==============================================================================
umap_plot <- DimPlot(
  seurat_lung, 
  reduction = "umap", 
  label = TRUE, 
  repel = TRUE, 
  label.size = 3.5
) +
  labs(
    title = "Whole Lung Fibroblast UMAP", 
    x = "UMAP 1", 
    y = "UMAP 2"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 9),
    axis.text = element_text(size = 8, color = "black")
  )

print(umap_plot)
