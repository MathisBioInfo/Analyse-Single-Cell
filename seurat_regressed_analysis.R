# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

library(AnnotationHub)
library(ensembldb)

setwd("/home/stagiaire/Documents/Mathis/script")
load("data/filtered_seurat.RData")

if (file.exists("regression/regressed_seurat.rds") == FALSE) {
  print("integrating seurat object")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  seurat_phase <- SCTransform(filtered_seurat, method = "glmGamPoi", verbose = FALSE)
  
  seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase))
  
  seurat_phase <- CellCycleScoring(seurat_phase, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  # Split seurat object by condition to perform cell cycle scoring and SCT on all samples
  split_seurat <- SplitObject(seurat_phase, split.by = "sample")
  
  split_seurat <- split_seurat[c("nn_strv", "strv")]
  
  options(future.globals.maxSize = 4000 * 1024^2)
  
  for (i in 1:length(split_seurat)) {
    split_seurat[[i]]$CC.Difference <- split_seurat[[i]]$S.Score - split_seurat[[i]]$G2M.Score
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio","S.Score", "G2M.Score"))
  }
  
  # Select the most variable features to use for integration
  integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                              nfeatures = 3000)
  
  # Prepare the SCT list object for integration
  split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                     anchor.features = integ_features)
  
  # Find best buddies - can take a while to run
  integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                          normalization.method = "SCT", 
                                          anchor.features = integ_features)
  
  # Integrate across conditions
  seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                     normalization.method = "SCT")
  
  # Run PCA
  seurat_integrated <- RunPCA(object = seurat_integrated)
  
  # Plot PCA
  jpeg("regression/PCAPlot.jpg")
  PCAPlot(seurat_integrated,
          split.by = "sample")  
  dev.off()
  
  # Run UMAP
  seurat_integrated <- RunUMAP(seurat_integrated, 
                               dims = 1:40,
                               reduction = "pca")
  
  # Plot UMAP  
  jpeg("regression/UMAPIntegrated.jpg") 
  DimPlot(seurat_integrated)
  dev.off()
  
  # Plot UMAP split by sample
  jpeg("regression/UMAPIntegratedSplited.jpg")
  DimPlot(seurat_integrated,
          split.by = "sample")  
  dev.off()
  
  # Defining the information in the seurat object of interest
  columns <- c(paste0("PC_", 1:16),
               "ident",
               "UMAP_1", "UMAP_2")
  
  # Extracting this data from the seurat object
  pc_data <- FetchData(seurat_integrated, 
                       vars = columns)
  # Adding cluster label to center of cluster on UMAP
  umap_label <- FetchData(seurat_integrated, 
                          vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
    group_by(ident) %>%
    summarise(x=mean(UMAP_1), y=mean(UMAP_2))
  
  # Plotting a UMAP plot for each of the PCs
  jpeg("regressionUMAP_PC.jpg")
  map(paste0("PC_", 1:16), function(pc){
    ggplot(pc_data, 
           aes(UMAP_1, UMAP_2)) +
      geom_point(aes_string(color=pc), 
                 alpha = 0.7) +
      scale_color_gradient(guide = "none", 
                           low = "grey90", 
                           high = "blue")  +
      geom_text(data=umap_label, 
                aes(label=ident, x, y)) +
      ggtitle(pc)
  }) %>% 
    plot_grid(plotlist = .)
  dev.off()
  
  # Save integrated seurat object
  saveRDS(seurat_integrated, "regression/regressed_seurat.rds")
} else {
  print("loading integrated seurat object")
  seurat_integrated <- readRDS("regression/regressed_seurat.rds")
}
