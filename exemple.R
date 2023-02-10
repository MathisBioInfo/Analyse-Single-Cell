#Package Used to perform RNAseq analysis
library(Seurat)
#Package Used to draw plot
library(ggplot2)

# Read
seurat_integrated <- readRDS("data/regressed_seurat.rds")
head(seurat_integrated)
tail(seurat_integrated)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

#Split the data by sample condition
Idents(seurat_integrated) <- "sample"
#ŧest
conserved.markers <- FindConservedMarkers(seurat_integrated, ident.1 = 3, grouping.var = 'sample')
cluster.markers <- FindAllMarkers(seurat_integrated)
head(conserved.markers)

# Plot all the data points
p <- ggplot(object_marker, aes(x = avg_log2FC, y = -log10(p_val_adj)) +
  geom_point() +
  xlab("avg_log2FC") +
  ylab("-log10(p_val_adj)")

jpeg("VolcanoPlot.jpg")
p
dev.off()

# Extract the top 10 labels based on adjusted p-value
top_10 <- head(all.markers[order(all.markers$p_val_adj),], 10)
top_10

# Manually annotate only the top 10 labels
p <- p + geom_text(data = top_10,
              aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(top_10)),
              hjust = 0, vjust = 0,
              size = 3)

jpeg("labelVolcanoPlot.jpg")
p + geom_text(data = top_10,
              aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(top_10)),
              hjust = 0, vjust = 0,
              size = 3)
dev.off()

# Create a new column indicating whether the gene meets the cut-off criteria
all.markers$color <- ifelse(abs(all.markers$avg_log2FC) > 1 & all.markers$p_val_adj < 10^-6, "red", "black")
all.markers$color <- ifelse(abs(all.markers$avg_log2FC) > 1 & all.markers$p_val_adj < 10^-6, "red", 
                                  ifelse(all.markers$p_val_adj < 10^-6, "blue", 
                                         ifelse(abs(all.markers$avg_log2FC) > 1, "green", "black")))

# Extract the top 10 labels based on adjusted p-value
top_10 <- head(conserved.markers[order(conserved.markers$p_val_adj),], 10)

# Make the final plot
p <- ggplot(all.markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
  geom_point() +
  geom_text(data = top_10,
            aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(top_10)),
            hjust = 0, vjust = 0,
            size = 3) +
  xlab("avg_log2FC") +
  ylab("-log10(p_val_adj)") +
  ggtitle("Volcano Plot of control vs starved") +
  scale_color_identity() +
  guides(color = guide_legend(title = "Color"))

jpeg("finalVolcanoPlot.jpeg")
p
dev.off()



volcano_Seurat <- function(object_marker, n, FC_cut_off, p_value_cut_off){
  # Create a new column indicating whether the gene meets the cut-off criteria
  object_marker$color <- ifelse(abs(object_marker$avg_log2FC) > FC_cut_off & object_marker$p_val_adj < p_value_cut_off, "red", "black")
  object_marker$color <- ifelse(abs(object_marker$avg_log2FC) > FC_cut_off & object_marker$p_val_adj < p_value_cut_off, "red", 
                              ifelse(object_marker$p_val_adj < p_value_cut_off, "blue", 
                                     ifelse(abs(object_marker$avg_log2FC) > FC_cut_off, "green", "black")))
  
  # Extract the top 10 labels based on adjusted p-value
  top_n <- head(object_marker[order(object_marker$p_val_adj),], n)
  
  # Make the final plot
  p <- ggplot(object_marker, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point() +
    geom_text(data = top_n,
              aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(top_n)),
              hjust = 0, vjust = 0,
              size = 3) +
    xlab("avg_log2FC") +
    ylab("-log10(p_val_adj)") +
    ggtitle("Volcano Plot of control vs starved") +
    scale_color_identity() +
    guides(color = guide_legend(title = "Color"))
  ggsave("functionVolcano.jpg")
  
   gene_list <- row.names(object_marker[object_marker$color=="red",])
   print(gene_list)
   return(gene_list)
}



converter_Symbol_2_Entrez <- function(gene_list){
  hs <- org.Hs.eg.db
  my.symbols <- gene_list
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
  print(head(egox))
  return(egox)
}

egox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

jpeg("cnetplot_exemple_MF_ALL.jpg", height = 1200, width = 1200, quality = 100)
cnetplot(egox,
         showCategory = 10,
         colorEdge = TRUE,
         node_label="category",
         cex_lab_category = 1.2,
         cex_label_gene = 0.8)
dev.off()

osteo <- subset(object_marker, idents = 3)
Idents(osteo) <- "sample"
avg.osteo <- log1p(AverageExpression(osteo, verbose = FALSE)$RNA)
avg.osteo <- as.data.frame(avg.osteo)
avg.osteo$gene <- row.names(avg.osteo)
avg.osteo$outlier <- abs(avg.osteo$nn_strv - avg.osteo$strv)
gene.outlier <- head(avg.osteo[order(-avg.osteo$outlier),],15)
jpeg("test_DE.jpg")
ggplot(avg.osteo, aes(nn_strv,strv)) +
  geom_point() +
  geom_text(data = gene.outlier, 
            aes(nn_strv,strv, label = row.names(gene.outlier)),
            hjust = 0, vjust = 0,
            size = 3) +
  ggtitle("Osteo")
dev.off()
