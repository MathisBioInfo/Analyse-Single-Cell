# Bring in Seurat object
seurat <- readRDS("data/regressed_seurat.rds")

library(SingleCellExperiment)


# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts 

metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample")]

# Explore the raw counts for the dataset

## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))

counts(sce)[1:6, 1:6]

## Explore the cellular metadata for the dataset
dim(colData(sce))

head(colData(sce))
