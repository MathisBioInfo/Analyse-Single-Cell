library(Seurat)
library(ggplot2)
library(sctransform)

setwd("/home/stagiaire/Documents/Mathis")

load("data/seurat_phase.RData")

### TEST SCTransform

pbmc <- seurat_phase

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = c("mitoRatio","G2M.Score"), verbose = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("glmGamPoi")
pbmc <- SCTransform(pbmc, method = "glmGamPoi", vars.to.regress = c("mitoRatio","G2M.Score"), verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
jpeg("DimPlot_SCT.jpg")
DimPlot(pbmc, label = TRUE) + NoLegend()
dev.off()

### VRAIE PARTIE

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

split_seurat <- split_seurat[c("nn_strv", "strv")]

options(future.globals.maxSize = 4000 * 1024^2)


for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio","G2M.Score"))
}

# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")

# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")

### TEST INTEGRATION

## DO NOT RUN (you have already run this)
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

## DO NOT RUN (you have already run this)
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

## DO NOT RUN (you have already run this)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

## DO NOT RUN (you have already run this)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
jpeg("PCAPlot.jpg")
PCAPlot(seurat_integrated,
        split.by = "sample")  
dev.off()

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP  
jpeg("UMAPIntegrated.jpg")# Plot UMAP split by sample
DimPlot(seurat_integrated,
        split.by = "sample")  
DimPlot(seurat_integrated)
dev.off()

# Plot UMAP split by sample
jpeg("UMAPIntegratedSplited.jpg")
DimPlot(seurat_integrated,
        split.by = "sample")  
dev.off()

# Save integrated seurat object
saveRDS(seurat_integrated, "data/integrated_seurat.rds")
