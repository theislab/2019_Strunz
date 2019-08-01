## Generates all Figures in Supplementary Figure 3, based on whole lung data set
## taken, adapted from Bleo/IM_AM_Score/with_published_data/try_limma_DGE_V2.R
## and Bleo/IM_AM_Score/calculate_IM_AM_V2.R

## Generate Subset of Macrophages
library(dbscan)
library(Seurat)
library(limma)

library(data.table)
library(ggplot2)
library(plyr)

## whole lung data set
load("../data/Whole_lung_seurat.RData")

## Subset to Macrophage Cluster
cluster <- c(0, 6, 2, 14, 19, 21, 3, 8)
subset <- SubsetData(seu.ica, accept.value = cluster, subset.name = "res.2")

## Add cell type names
cell.types <- c("activated AM", "AM", "Arg1+ M2", "CD 163- IM", "CD 163+ IM",
                "Cebpb+ AM", "Mfge8+ resolution IM", "Monocytes")
names(cell.types) <- cluster
subset@meta.data$cell.type <- cell.types[subset@meta.data$res.2]

## Clean subset with dbscan and recalculate PCs with new set of variable genes
subset@meta.data$dbscan <- dbscan(subset@dr$umap@cell.embeddings, 0.2, 4)$cluster
subset <- SubsetData(subset, accept.value = "1", subset.name = "dbscan")

subset <- FindVariableGenes(subset, y.cutoff = 0, do.plot = F)
subset <- RunPCA(subset, pcs.compute = 50, do.print = F)
subset <- RunUMAP(subset, reduction.use = "pca", dims.use = 1:20, n_neighbors = 20)
DimPlot(object = subset, reduction.use = "umap", do.label = T, group.by = "res.2", no.legend = T)

## Scoring to bulk count data provided by Misharin et. al
mat <- read.delim(file = "../data/relevant_bulk_all.txt", sep = "\t", header = T)

## DGE between these 4 conditions
mat.cur <- mat[, -which(colnames(mat) %in% c("Symbol", "gene"))]
all.conds = c("mono", "im", "mo_am", "tr_am")
conds = colnames(mat.cur)
for(c in all.conds){
  conds[grepl(c, conds)] = c
}

## Create Model Matrix and Contrasts
design <- model.matrix(~0+conds)
colnames(design) <- gsub("conds", "", colnames(design))

contr.matrix <- makeContrasts(imVSrest = im - (mo_am + mono + tr_am)/3,
                              monoVSrest = mono - (mo_am + im + tr_am)/3,
                              tr_amVSrest = tr_am - (mo_am + im + mono)/3,
                              mo_amVSrest = mo_am - (tr_am + mono + im)/3,
                              levels = colnames(design))
contr.matrix
rownames(mat.cur) <- mat$gene

## Logarithmize Matrix and fit Model
mat.cur <- log(mat.cur)
vfit <- lmFit(mat.cur, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

im.vs.rest <- topTreat(efit, coef = "imVSrest", n = Inf, sort.by = "P")
mono.vs.rest <- topTreat(efit, coef = "monoVSrest", n = Inf, sort.by = "P")
tr_am.vs.rest <- topTreat(efit, coef = "tr_amVSrest", n = Inf, sort.by = "P")
mo_am.vs.rest <- topTreat(efit, coef = "mo_amVSrest", n = Inf, sort.by = "P")

## Correlate gene expression in Macrophage subset to the log foldchanges from bulk data
calc_mish_score <- function(seu, dge.data, scorename){
  topmarker <- data.frame(logFC = dge.data$logFC, genes = rownames(dge.data))
  
  topmarker = topmarker[topmarker$genes %in% rownames(seu@data),]
  topmarker <- topmarker[1:500,]
  genExpr <- FetchData(seu, use.scaled = F, vars.all = topmarker$genes) 
  
  res <- cor(topmarker$logFC, t(genExpr))
  clust <- FetchData(seu, vars.all = "res.2", cells.use = colnames(res))
  res <- data.frame(corr = res[1,], cluster = clust)
  colnames(res) <- c(scorename, "cluster")
  res
}

plot_score <- function(seu, col, scorename){
  colors <- colorRampPalette(c("grey", col))(length(unique(seu@meta.data[,scorename])))
  p <- DimPlot(seu, reduction.use = "umap", do.label = F, plot.title = scorename,
               pt.size = 1, group.by = scorename, no.legend = T, cols.use = colors, do.return = T)
  p
}

subset <- AddMetaData(subset, calc_mish_score(subset, im.vs.rest, "IM_score"))
subset <- AddMetaData(subset, calc_mish_score(subset, mono.vs.rest, "mono_score"))
subset <- AddMetaData(subset, calc_mish_score(subset, mo_am.vs.rest, "mo_am_score"))
subset <- AddMetaData(subset, calc_mish_score(subset, tr_am.vs.rest, "AM_score"))

## Assign cells to either IM or AM category, if difference is greater thresh
thresh = 0.05
scores <- subset@meta.data[, c("AM_score", "IM_score")]
scores$origin <- "unassigned"
scores$origin[scores$AM_score > scores$IM_score & abs(scores$AM_score - scores$IM_score) > thresh ] <- "AM"
scores$origin[scores$AM_score < scores$IM_score & abs(scores$AM_score - scores$IM_score) > thresh] <- "IM"

subset <- AddMetaData(subset, data.frame(row.names = rownames(scores), origin = scores$origin))

### Generate Figures ###

## Supplementary Figure A
DimPlot(subset, reduction.use = "umap", group.by = "cell.type", pt.size = 2,
             vector.friendly = T, do.return = T)

## Supplementary Figure B - Feature Plots
genes = c("Cd68", "Mrc1", "Cd163", "Cebpb", "Arg1", "Mfge8", "Mmp12", "Ear2")
for(gene in genes){
  FeaturePlot(subset, reduction.use = "umap", features.plot = gene, pt.size = 2, cols.use = c("gray", "red"),
                   vector.friendly = T, do.return = T)[[1]]
}

## Supplementary Figure C
DimPlot(subset, reduction.use = "umap", group.by = "grouping", pt.size = 2,
             vector.friendly = T, do.return = T)

## Supplementary Figure D
DimPlot(subset, reduction.use = "umap", group.by = "origin", pt.size = 1.5,
             cols.use = c("#0101DF", "#FAE316", "gray"), do.return = T)

## Supplementary Figure E
thresh = 0.1
values <- rownames(subset@meta.data[subset@meta.data["mo_am_score"] > thresh,])
DimPlot(object = subset, cells.highlight = values, reduction.use = "umap",
             do.label = F, sizes.highlight = 1.5, pt.size = 1.5, no.legend = T, 
             cols.use = "gray", cols.highlight = "#B40404", do.return = T,
             plot.title = paste("Monocyte derived (lineage trace) - Threshold", thresh))

## Figure F, G, H and I
## Taken from IM_AM_score/subset_kinetics_plots and avg_expression_gene_in_cluster_V2
genLinePlot <- function(subset, gene, col = "blue", lev = 0.95, save = F){
  meta <- subset@meta.data
  identifier <- meta$identifier
  logCPM_ok <- subset@data
  
  expr_tmp <- logCPM_ok[gene, ]
  means <- unlist(lapply(split(expr_tmp, identifier), mean))
  
  timepoint <- unique(meta[,c("identifier", "grouping")])
  timepoint$grouping <- gsub("d", "", as.character(timepoint$grouping))
  timepoint$grouping <- gsub("PBS", "0", timepoint$grouping)
  data <- data.frame(expression = means, day = timepoint$grouping[match(names(means), timepoint$identifier)])
  data$day <- as.numeric(as.character(data$day))
  
  ggplot(data, aes(y = expression, x = day)) + 
      geom_smooth(se = T, level = lev, col = col) +
      ylab(paste0("mean expression")) + xlab("Days") + ggtitle(gene) +
      theme(plot.title = element_text(face = "italic", hjust = 0.5))
}

gene = "Arg1"
genLinePlot(subset, gene, col = "blue")

gene = "Mfge8"
genLinePlot(subset, gene, col = "blue")

gene = "Cebpb"
genLinePlot(subset, gene, col = "blue")

gene = "Ear2"
genLinePlot(subset, gene, col = "blue")



