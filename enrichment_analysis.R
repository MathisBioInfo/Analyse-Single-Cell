if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}

if (!("enrichplot" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("enrichplot", update = FALSE)
}

if (!("pathview" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("pathview", update = FALSE)
}


#Package Used to perform RNAseq analysis
library(Seurat)
#Genome wide annotation for Human, mostly using Entrez Gene identifiers
library(org.Hs.eg.db)
#Enrichment tool for gene annotation
library(clusterProfiler)
#Visualization methods for enrichment results
library(enrichplot)
#Tool to visualize enrichment results over pathway graphs
library(pathview)

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Read
seurat_integrated <- readRDS("data/regressed_seurat.rds")

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

cluster.markers <- FindAllMarkers(seurat_integrated)
df <- subset(cluster.markers,select=c('gene','avg_log2FC'))

# Extract the top 10 labels based on adjusted p-value
top_10 <- head(cluster.markers[order(cluster.markers$p_val_adj),], 10)
top_10

cluster3.markers <- FindMarkers(seurat_integrated, ident.1 = 3)
df <- subset(cluster3.markers,select=c('avg_log2FC'))

hs <- org.Hs.eg.db
my.symbols <- row.names(df)
gene <- AnnotationDbi::select(hs, 
                      keys=my.symbols, 
                      columns=c("SYMBOL", "ENTREZID"), 
                      keytype="SYMBOL")

ego <- enrichGO(gene          = gene[,2],
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

egox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

jpeg("cnetplot_MF_ALL.jpg", height = 1200, width = 1200, quality = 100)
cnetplot(egox,
         showCategory = 10,
         colorEdge = TRUE,
         node_label="category",
         cex_lab_category = 1.2,
         cex_label_gene = 0.8)
dev.off()

egox2 <- pairwise_termsim(egox)

jpeg("treeplot_c3.jpg", height = 1200, width = 1200, quality = 100)
treeplot(egox2)
dev.off()

ego <- pairwise_termsim(ego)
jpeg("emapplot_c3.jpg", height = 1200, width = 1200, quality = 100)
emapplot(ego)
dev.off()

library(EnhancedVolcano)
top10 <- head(cluster3.markers[order(cluster3.markers$p_val_adj),], 10)
jpeg("volcano_c3.jpg")
EnhancedVolcano(cluster3.markers,
                lab = rownames(cluster3.markers),
                labels = rownames(top10),
                x = 'avg_log2FC',
                y = 'p_val_adj')
dev.off()

library(ggplot2)

# Extract the top 10 labels based on adjusted p-value
top_10 <- head(cluster3.markers[order(cluster3.markers$p_val_adj),], 10)

# Create a new column indicating whether the gene meets the cut-off criteria
cluster3.markers$color <- ifelse(abs(cluster3.markers$avg_log2FC) > 2 & cluster3.markers$p_val_adj < 10^-6, "red", "black")
cluster3.markers$color <- ifelse(abs(cluster3.markers$avg_log2FC) > 2 & cluster3.markers$p_val_adj < 10^-6, "red", 
                                 ifelse(cluster3.markers$p_val_adj < 10^-6, "blue", 
                                        ifelse(abs(cluster3.markers$avg_log2FC) > 2, "green", "black")))
# Plot all the data points
p <- ggplot(cluster3.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
  geom_point() +
  xlab("avg_log2FC") +
  ylab("-log10(p_val_adj)") +
  scale_color_identity()

# Manually annotate only the top 10 labels
jpeg("testplot.jpg")
p + geom_text(data = top_10,
              aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(top_10)),
              hjust = 0, vjust = 0,
              size = 3)
dev.off()


?jpeg("dotplot.jpg")
dotplot(ego,showCategory = 20)
dev.off

function(){
  
}

#test comparaison entre sample
Idents(seurat_integrated) <- "sample"
conserved.markers <- FindAllMarkers(seurat_integrated)
df <- subset(conserved.markers,select=c('avg_log2FC'))

jpeg("volcano_strv.jpg")
EnhancedVolcano(conserved.markers,
                lab = rownames(conserved.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj')
dev.off()

# Extract the top 10 labels based on adjusted p-value
top_10 <- head(conserved.markers[order(conserved.markers$p_val_adj),], 10)

# Create a new column indicating whether the gene meets the cut-off criteria
conserved.markers$color <- ifelse(abs(conserved.markers$avg_log2FC) > 2 & conserved.markers$p_val_adj < 10^-6, "red", "black")
conserved.markers$color <- ifelse(abs(conserved.markers$avg_log2FC) > 2 & conserved.markers$p_val_adj < 10^-6, "red", 
                                  ifelse(conserved.markers$p_val_adj < 10^-6, "blue", 
                                         ifelse(abs(conserved.markers$avg_log2FC) > 2, "green", "black")))
# Plot all the data points
p <- ggplot(conserved.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
  geom_point() +
  xlab("avg_log2FC") +
  ylab("-log10(p_val_adj)") +
  scale_color_identity()

# Manually annotate only the top 10 labels
jpeg("testvolcanoplot.jpg")
p + geom_text(data = top_10,
              aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(top_10)),
              hjust = 0, vjust = 0,
              size = 3)
dev.off()


df <- subset(conserved.markers,select=c('avg_log2FC'))

hs <- org.Hs.eg.db
my.symbols <- row.names(df)
gene <- AnnotationDbi::select(hs, 
                              keys=my.symbols, 
                              columns=c("SYMBOL", "ENTREZID"), 
                              keytype="SYMBOL")

ego <- enrichGO(gene          = gene[,2],
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

egox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

jpeg("cnetplot_BP_cond.jpg", height = 1200, width = 1200, quality = 100)
cnetplot(egox,
         showCategory = 20,
         colorEdge = TRUE,
         node_label="category",
         cex_lab_category = 1.2,
         cex_label_gene = 0.8)
dev.off()

egox2 <- pairwise_termsim(egox)

jpeg("treeplot_cond.jpg", height = 1200, width = 1200, quality = 100)
treeplot(egox2)
dev.off()
