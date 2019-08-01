# R code
# Lukas Simon
# This code will generate parts of Strunz et al Fig 1

# Load R libs ####
library(Seurat)
library(Matrix)
library(limma)
library(preprocessCore)
library(reshape)
library(data.table)
library(ggplot2)

# Load Seurat object ####
load('../data/Whole_lung_Seurat.RData')

# Generate Fig 1a ####
DimPlot(seu.ica, reduction.use = "umap", group.by = "grouping", vector.friendly = T)

# Load bulk and in silico bulk RNAseq expression matrices ###
insilico <- read.csv('../data/MultiomicIntegration/Bulk_singleCell.csv', row.names = 1)
bulk <- read.csv('../data/MultiomicIntegration/Bulk_transcriptome.csv', row.names = 1)
ok <- intersect(rownames(insilico), rownames(bulk))
insilico <- insilico[ok, ]
bulk <- bulk[ok,]

# Load bulk protein expression matrix ####
protein <- read.csv('../data/MultiomicIntegration/Bulk_protein.csv', row.names = 1)

# Match genes for all three matrices ####
ok <- intersect(rownames(protein), rownames(insilico))
protein_ok <- protein[ok,]
bulk_ok <- bulk[ok,]
insilico_ok <- insilico[ok, ]
rownames(bulk_ok) <- rownames(protein_ok) <- rownames(insilico_ok) <- ok

# Normalize and merge all three matrices ####
bulk_ok <- voom(bulk_ok, normalize.method = 'none')$E
insilico_ok <- voom(insilico_ok)$E
merged <- data.matrix(cbind(bulk_ok, insilico_ok, protein_ok))
merged <- normalize.quantiles(merged)
colnames(merged) <- c(colnames(bulk_ok), colnames(insilico_ok), colnames(protein_ok))
rownames(merged) <- rownames(bulk_ok)

# Get the day information ####
tmp <- do.call(rbind, lapply(colnames(insilico), function(x) strsplit(x, '_', fixed = T)[[1]]))
day <- tmp[,3]
day[which(tmp[,2] == 'PBS')] <- 'd0'

# Define sample attributes ####
batch <- c(rep("bulk", ncol(bulk_ok)), rep("insilico", ncol(insilico_ok)), rep("protein", ncol(protein_ok)))
treat_bulk <- unlist(lapply(colnames(bulk_ok), function(x) strsplit(x, '_', fixed = T)[[1]][1]))
treat_prot <- c(rep('PBS', 3), rep('Bleo', 3))
treat_insilico <- unlist(lapply(colnames(insilico_ok), function(x) strsplit(x, '_', fixed = T)[[1]][2]))
treat <- c(treat_bulk, treat_insilico, treat_prot)

farben <- c("blue", "red")
names(farben) <- c('PBS', 'Bleo')

shape <- c(1, 2, 3)
names(shape) <- c("bulk", "insilico", "protein")

# Calculate PCA ####
pca <- prcomp(t(merged))

# Generate Fig 1b ####
par(mfrow = c(1, 2))
plot(pca$x[,1:2], col = farben[treat], pch = shape[batch])
legend("bottomright", c(names(farben), names(shape)), pch = c(16, 16, shape), bty = "n", col = c("blue", "red", "black", "black", "black"))

plot(pca$x[,2:3], col = farben[treat], pch = shape[batch])
legend("bottomleft", c(names(farben), names(shape)), pch = c(16, 16, shape), bty = "n", col = c("blue", "red", "black", "black", "black"))

# Generate Fig 1c ####
par(mfrow = c(1, 1))
loadings <- pca$rotation[,3]
top_loadings <- c(head(sort(loadings), 20), tail(sort(loadings), 20))
farben <- c(colorRampPalette(c("darkblue", "blue"))(20), colorRampPalette(c("red", "darkred"))(20))
barplot(top_loadings, las = 2, ylab = 'Loadings PC3', col = farben, horiz = T)

# Generate Fig 1d ####
ok <- grep('muc', colnames(merged))
boxplot(split(pca$x[ok, 3], day)[c('d0', 'd3', 'd7', 'd10', 'd14', 'd21', 'd28')], las = 2, ylab = 'PC3')

# Define functions for plotting relative freqeuncies ####
get_colours <- function(seu, groupby = "grouping"){
  p <- Seurat::DimPlot(seu, group.by = groupby, do.return = T, reduction.use = "umap")
  pdata <- ggplot2::ggplot_build(p)$data[[1]]
  cols <- unique(pdata[,c("colour", "group")])
  newCols <- cols$colour
  names(newCols) <- levels(as.factor(seu@meta.data[,groupby]))[cols$group]
  newCols
}
calcFreq <- function(seu, by1, by2) table(as.character(seu@meta.data[,by1]), as.character(seu@meta.data[,by2]))
get_relFreqs <- function(seu, grouping = "batch"){
  freqByGrouping <- as.data.frame.matrix(calcFreq(seu, by1 = 'identifier', by2 = 'res.2'))
  
  relFreqs <- data.frame(freqByGrouping / apply(freqByGrouping, 1, sum), check.names = F)
  relFreqs$id <- rownames(freqByGrouping)
  
  temp <- unique(seu@meta.data[, c("identifier", grouping)])
  map <- temp[, grouping]
  names(map) <- temp$identifier
  relFreqs$group <- map[relFreqs$id]
  
  freqs <- melt(relFreqs, id.vars = c("id", "group"), variable.name = "cluster")
  freqs
}
plot_relFreqs <- function(freqs, clust, cols, order, title = F, save = F, 
                          path = "Plots/rel_frequencies/"){
  freqs$cluster <- paste("Cluster", freqs$cluster)
  freqs$group <- factor(freqs$group, order)
  
  clust = paste("Cluster", clust)
  if(title == F){
    title = clust
  }
  cols <- cols[order]
  p <- ggplot(subset(freqs, cluster %in% clust)) + geom_boxplot(aes(x = group, y = value), fill= cols) +
    labs(y = "rel. freqency", x = "group", title = title) + guides(fill=FALSE) +
    scale_x_discrete(labels = order)
  
  if(save == T){
    ggsave(paste0(gsub(" ", "_", title), ".pdf"), plot = p, device = "pdf", width = 6, height = 5, units = "in", path = path)
  }
  else{
    plot(p)
  }
}


# Define the order ####
order <- c("PBS", "d3", "d7", "d10", "d14", "d21", "d28")
freqs <- get_relFreqs(seu.ica, grouping = "grouping")
cols <- get_colours(seu.ica, groupby = "grouping")

# Generate Fig 1e ####
clust = "9"
plot_relFreqs(freqs, clust, cols, order, title = "Mki67+ T cells (proliferating)", save = F, path = "<filepath>")

# Generate Fig 1f ####
clust = "8"
plot_relFreqs(freqs, clust, cols, order, title = "Ly6c2+ monocytes", save = F, path = "<filepath>")

# Generate Fig 1g ####
clust = "2"
plot_relFreqs(freqs, clust, cols, order, title = "Arg1+ M2-macrophages", save = F, path = "<filepath>")

# Generate Fig 1h ####
clust = "18"
plot_relFreqs(freqs, clust, cols, order, title = "Acta2+ myofibroblasts", save = F, path = "<filepath>")

# Generate Fig 1i ####
clust = "3"
plot_relFreqs(freqs, clust, cols, order, title = "Mfge8+ resolution macrophages", save = F, path = "<filepath>")

