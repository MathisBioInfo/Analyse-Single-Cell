library(Seurat)
library(patchwork)

AU_dir <- "/media/HDDstorage/InVitro scRNAseq 2020/Fichiers_SC_Pour_Analyse/AU/filtered_feature_bc_matrix/"
AU_strv_dir <- "/media/HDDstorage/InVitro scRNAseq 2020/Fichiers_SC_Pour_Analyse/AU_strv//filtered_feature_bc_matrix/"
AU_data <- Read10X(data.dir = AU_dir)
AU_strv_data <- Read10X(data.dir = AU_strv_dir)
AU_srat = CreateSeuratObject(AU_data,project = "AU")
AU_strv_srat = CreateSeuratObject(AU_strv_data,project = "AU_strv")
meta_ctrl <- AU_strat@meta.data
meta_strv <- AU_strv_strat@meta.data
dim(meta_ctrl)
dim(meta_strv)
summary(meta_ctrl$nFeature_RNA)
summary(meta_strv$nFeature_RNA)
summary(meta_ctrl$nCount_RNA)
summary(meta_strv$nCount_RNA)
AU_srat[["percent.mt"]] <- PercentageFeatureSet(AU_srat, pattern = "^MT-")
AU_srat[["percent.rb"]] <- PercentageFeatureSet(AU_srat, pattern = "^RP[SL]")
AU_strv_srat[["percent.mt"]] <- PercentageFeatureSet(AU_strv_srat, pattern = "^MT-")
AU_strv_srat[["percent.rb"]] <- PercentageFeatureSet(AU_strv_srat, pattern = "^RP[SL]")
jpeg('figure/AU_Vln.jpg')
VlnPlot(AU_srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4)
dev.off()
jpeg('figure/AU_strv_Vln.jpg')
VlnPlot(AU_strv_srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4)
dev.off()

###Test###
library(ggplot2)
meta_ctrl %>% 
  ggplot(aes(color='red', x=nFeature_RNA, fill=nFeature_RNA)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
