# R code
# Lukas Simon
# This code will generate parts of Strunz et al Fig 2abcde

# Load R libs ####
library(Seurat)
library(Matrix)
library(readxl)

# Load Whole Lung Alveolar Epithelium Seurat object ####
load('../data/Whole_lung_seurat_alveolar_epithelium.RData')

# Identify DE genes ####
markers <- FindAllMarkers(seu)

# Generate Fig 2a ####
DimPlot(seu, reduction.use = "umap", do.label = T)

# Generate Fig 2b ####
day <- as.character(seu@meta.data$grouping)
day <- gsub("PBS", "d0", day)
day <- factor(day, levels = c("d0", "d3", "d7", "d10", "d14", "d21", "d28"))
seu@meta.data$day <- day
DimPlot(seu, reduction.use = "umap",group.by = "day")

# Generate Fig 2e ####
clusters <- c(1, 12, 15, 39)
top_genes<- unlist(lapply(clusters, function(x){
  subm <- markers[which(markers$cluster == x),]
  subm <- subm[which(as.numeric(subm$p_val_adj) < 0.25),]
  tail(subm$gene[order(subm$avg_logFC)], 50)
}))
cells <- unlist(lapply(clusters, function(x){
  subm <- seu@meta.data[which(seu@meta.data$res_2 == x),]
  tail(rownames(subm)[order(subm$nGene)], 35)
}))
DoHeatmap(seu, group.by = 'res_2', genes.use = top_genes, cells.use = cells, slim.col.label = T, use.scaled = T, disp.min = -2, disp.max = 2)

# Generate Fig d ####
FeaturePlot(seu, c("Sftpc", "Lcn2", "Pdpn", "Krt8"), reduction.use = "umap", cols.use = c("grey", "red"))
