# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

setwd("/home/stagiaire/Documents/Mathis")

load("data/seurat_filtered.RData")

filtered_seurat <- NormalizeData(filtered_seurat)

filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(filtered_seurat), 10)

plot1 <- VariableFeaturePlot(filtered_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2