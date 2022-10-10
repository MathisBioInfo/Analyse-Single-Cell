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

setwd("/home/stagiaire/Documents/Mathis")

load("data/seurat_filtered.RData")

# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seurat_phase)
seurat_phase <- ScaleData(seurat_phase, features = all.genes)

# Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle. Read it in with:

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Perform cell cycle scoring
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

# Perform PCA and color by cell cycle phase
seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(object = seurat_phase))

# Visualize the PCA, grouping by cell cycle phase
jpeg("test_acp.jpg")
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")
dev.off()

seurat_phase <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
seurat_phase <- RunPCA(seurat_phase, features = c(s_genes, g2m_genes))

jpeg("test_acp_reg.jpg")
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")
dev.off()

