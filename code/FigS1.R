# R code
# Lukas Simon
# This code will generate parts of Strunz et al Fig S1

# Load R libs ####
library(Seurat)
library(Matrix)
library(cluster)
library(dplyr)
library(tidyr)

# Load Seurat object ####
load('../data/Whole_lung_seurat.RData')

# Calculate Silhouette coefficient for mouse clustering ####
dist_euc <- dist(seu.ica@dr$ica@cell.embeddings, method = "euclidean")
identifier <- as.numeric(as.factor(seu.ica@meta.data$identifier))
sil <- silhouette(identifier, dist_euc)
summary(sil)[["avg.width"]]

# Generate Fig S1a ####
farben <- c(rgb(0,0,1,0.1), rgb(1,0,0,0.5))
seu.ica@meta.data$pbs <- seu.ica@meta.data$grouping == "PBS"
r <- sample(rownames(seu.ica@meta.data))
seu.ica@meta.data <- seu.ica@meta.data[r,]
seu.ica@dr$umap@cell.embeddings <- seu.ica@dr$umap@cell.embeddings[r,]
DimPlot(seu.ica, reduction.use = "umap", group.by = "pbs", vector.friendly = T, cols.use = farben, pt.size = 0.5)

# Generate Fig S1b ####
p <- DimPlot(seu.ica, reduction.use = "umap", group.by = "orig.ident", vector.friendly = T)
p
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]]
farben <- unlist(lapply(split(pdata$colour, seu.ica@meta.data$identifier), unique))
DimPlot(seu.ica, reduction.use = "umap", group.by = "identifier", vector.friendly = T, pt.size = 0.5)

# Generate Fig S1c ####
DimPlot(seu.ica, reduction.use = "umap", group.by = "res.2", vector.friendly = T, do.label = T)

# Generate Fig S1g ####
genes <- c("Acoxl", "Fabp1", "Lgi3", "Slc34a2", "Gzma", "Arhgef38", "Aox3", "Lamp3", "Ppp1r14c", "Lyz1",
           "Chi3l1", "Itih4", "Cxcl15", "Acot1", "Enpep", "Sftpa1", "Atp6v1c2", "Slc4a5", "Ank3", "Tinag",
           "Pgam1", "Lox", "Ecm1", "Tpm2", "Cxcl14", "C1qc", "Thbs1", "C1qa", "C1qb", "Mmp19", "Mmp12",
           "Emilin2", "S100a4", "Serpina3n", "Gpnmb", "Fn1", "Msr1", "Tnc", "Spp1", "Arg1")

noms <- rbind(c(11, "Ciliated cells"),
  c(37, "Ciliated cell subet"),
  c(36, "low quality cells"),
  c(27, "Goblet cells"),
  c(7, "Club cells"),
  c(1, "AT2 cells"),
  c(12, "activated AT2"),
  c(15, "Krt8+ ADI"),
  c(39, "AT1 cells"),
  c(28, "LECs"),
  c(32, "CECs"),
  c(13, "VECs"),
  c(24, "Vcam1+ VECs"),
  c(23, "Mesothelial cells"),
  c(25, "activated Mesothelial cells"),
  c(35, "SMCs"),
  c(18, "Myofibroblasts"),
  c(22, "Fibroblasts"),
  c(26, "Ccl17+ DCs"),
  c(16, "CD103+ DCs"),
  c(4, "DCs"),
  c(30, "Plasma cells"),
  c(31, "NK cells"),
  c(10, "B-lymphocytes"),
  c(5, "T-lymphocytes"),
  c(2, "M2 macrophages"),
  c(3, "resolution macrophage"),
  c(8, "recruited monocytes"),
  c(14, "CD163-/CD11c+ IMs"),
  c(19, "CD163+/CD11c- IMs"),
  c(0, "AM (Bleo)"),
  c(21, "AM (Bleo)"),
  c(6, "AM (PBS)"),
  c(20, "Neutrophils"),
  c(29, "Fn1+ macrophages"),
  c(33, "non-classical monocytes (Ly6c2-)"),
  c(34, "Themis+ T-lymphocytes"),
  c(9, "Mki67+ Proliferating cells"),
  c(38, "T cell subset"),
  c(17, "Mki67+/Top2a+ Proliferting cells"),
  c(40, "Mki67+ Proliferting cells"))
noms <- data.frame(noms)
colnames(noms) <- c("res.2", "name")
noms$res.2 <- as.numeric(as.character(noms$res.2))

seu.ica@meta.data$tmp <- NA
lapply(1:nrow(noms), function(x) seu.ica@meta.data$tmp[which(seu.ica@meta.data$res.2 == noms$res.2[x])] <<- as.character(paste(noms$res.2[x], noms$name[x], sep = " - ")))

cells <- rownames(seu.ica@meta.data)[which(!is.na(seu.ica@meta.data$tmp))]
seu.ica <- SubsetData(seu.ica, cells.use = cells)

PercentAbove <- function (x, threshold) length(x = x[x > threshold])/length(x = x)

data.to.plot <- data.frame(FetchData(object = seu.ica, vars.all = genes))
data.to.plot$cell <- rownames(x = data.to.plot)
data.to.plot$id <- factor(seu.ica@meta.data$tmp, levels = rev(paste(noms$res.2, noms$name, sep = " - ")))
data.to.plot <- data.to.plot %>% gather(key = genes.plot, value = expression, -c(cell, id))
data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression, threshold = 0))
data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, max = 2.5, min = -2.5))
data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, levels = rev(x = sub(pattern = "-", replacement = ".", x = genes)))
data.to.plot$pct.exp[data.to.plot$pct.exp < 0] <- NA

p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
  geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) + 
  scale_radius(range = c(0, 6)) + scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p