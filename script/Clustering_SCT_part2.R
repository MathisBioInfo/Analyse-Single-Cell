# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

setwd("/home/stagiaire/Documents/Mathis")
seurat_integrated <- readRDS("data/integrated_seurat.rds")

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Explore resolutions
seurat_integrated@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
jpeg("UMAPIntegratedClusterised.jpeg")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# UMAP of cells in each cluster by sample
jpeg("UMAPIntegratedClusterisedSplited.jpeg")
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
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
jpeg("UMAP_PC.jpg")
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
dev.off()

# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

jpeg("CellMarker.jpg")
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()
dev.off()

jpeg("FeaturePLotCell.jpg")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("SERPINE1", "TUBA1B"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()

# find all markers of cluster 3
cluster3.markers <- FindMarkers(seurat_integrated, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)

# find all markers of cluster 3
cluster1.markers <- FindMarkers(seurat_integrated, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)