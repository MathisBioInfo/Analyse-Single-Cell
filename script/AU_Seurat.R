library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

setwd("/home/stagiaire/Documents/Mathis")

dirAU <- "/media/HDDstorage/InVitro scRNAseq 2020/Fichiers_SC_Pour_Analyse/AU/filtered_feature_bc_matrix/"
dirAU_strv <- "/media/HDDstorage/InVitro scRNAseq 2020/Fichiers_SC_Pour_Analyse/AU_strv//filtered_feature_bc_matrix/"
adj.matrix <- Read10X(c(dirAU,dirAU_strv))
srat <- CreateSeuratObject(adj.matrix,project = "AU")
srat
adj.matrix <- NULL
str(srat)
meta <- srat@meta.data
dim(meta) #1409 rows 3 columns
head(meta) #columns index: [orig.ident,nCount_RNA,nFeature_RNA]
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

jpeg('figure/Vln.jpg')
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()
