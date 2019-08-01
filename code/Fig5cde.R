# R code
# Lukas Simon
# This code will generate parts of Strunz et al Fig4

# load R libs ####
library(plotrix)
library(mgcv)
library(slingshot)
library(SingleCellExperiment)
library(dbscan)
library(mclust)
library(RColorBrewer)
library(enrichR)
library(gam)
library(viridis)
load('../data/BlueYellowColormaps_V1.RData')

# load data ####
load("../data/EpiHiRes_seurat.RData")
louvain <- read.table('../data/EpiHiRes_louvain_labels.txt', row.names = 1)
coord <- read.csv('../data/EpiHiRes_diffusion_map_coordinates_convergence.txt', row.names = 1)
dm <- coord[,2:3]

# subset to common cells ####
ok <- intersect(rownames(dm), colnames(subset@data))
dm <- dm[ok,]
cm <- subset@raw.data
cm <- cm[, ok]
md <- subset@meta.data
md <- md[ok,]
md$louvain <- louvain[rownames(md),]

# Clean with DBSCAN ####
rd1 <- data.matrix(dm)
db <- dbscan(rd1, eps = 0.001)$cluster
ok <- which(db == names(which(table(db) == max(table(db)))))

# Fit trajectory using slingshot ####
rd1 <- rd1[ok,]
sim <- SingleCellExperiment(assays = List(counts = cm[,ok], norm = subset@data[, colnames(cm)[ok]]))
reducedDims(sim) <- SimpleList(DiffMap = rd1)
cl1 <- Mclust(rd1, G = 5)$classification
colData(sim)$GMM <- cl1
sce <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'DiffMap')

# Define pseudotime ####
t <- sce$slingPseudotime_1
t <- (t - min(t))/(max(t) - min(t))

# Define count object ####
cm <- assays(sce)$counts
md <- md[colnames(cm),]

# Define plotting function ####
plot_traj_single_gene <- function(gene = 'Krt8'){
  
  scale <- function(x) (x - min(x))/(max(x) - min(x))
  
  expr <- scale(data.matrix(log(cm[gene,] + 1)))
  
  aframe <- data.frame(expr, pt = t)
  
  melted <- melt(data = aframe, id.vars = "pt", measure.vars = setdiff(colnames(aframe), 'pt'))
  
  ggplot(melted, aes(y = value, x = pt, group = variable, color = variable)) + geom_smooth(method = "loess") +
    ylab(paste("Relative expression")) + xlab("Pseudotime") + ggtitle(gene)
}

# Generate Fig 5c ####
plot_traj_single_gene("Scgb1a1")
plot_traj_single_gene("Krt8")
plot_traj_single_gene("Sftpc")

# Generate Fig 5e ####
genes <- c("Nkx2-1", "Cebpa", "Foxp2",
           "Etv5", "Epas1", "Yeats4",
           "Gata6", "Hint1", "Kdm5c",
           "Rbpms", "Nupr1", "Sox4",
           "Sra1", "Nfkb2", "Nfkb1")
lapply(genes, plot_traj_single_gene)

# Restrict analysis to genes with more than 10 non-zero counts in more than 5 mice ####
asplit <- split(colnames(cm), as.character(md$identifier))
sums <- do.call(cbind, lapply(asplit, function(x){
  if(length(x) == 1) return(rep(0, nrow(cm)))
  apply(cm[,x], 1, function(y) sum(y > 0)) 
}))
genes <- rownames(cm)[which(apply(sums, 1, function(x) sum(x > 10)) > 5)]

# Run GAM regression - this may take some time (~ 1.5hrs) ####
tmp <- data.matrix(cm[genes,])
fits <- apply(tmp, 1, function(x){
  aframe <- data.frame(outcome = (x > 0), pt = t, nUMI = log(md$nUMI))
  gam::gam(outcome ~ nUMI + lo(pt), data = aframe, family = binomial(link = "logit"))
})
names(fits) <- genes
pvals <- unlist(lapply(fits, function(x) summary(x)[4][[1]][2, 5]))
names(pvals) <- genes

# Calculate fitted values ####
fitted_values <- do.call(rbind, lapply(fits, function(x){
  xpred <- seq(min(t), max(t), length = 50)
  aframe <- data.frame(pt = xpred, nUMI = mean(md$nUMI))
  predict(x, aframe)
}))

# Define significant genes ####
sig <- which(pvals < 1e-4)
tmp <- fitted_values[sig,]

# Load terminal state likelihoods ####
endstate <- read.csv('../data/EpiHiRes_convergence_terminal _state_likelihood.csv', row.names = 1)
endstate$id <- paste(endstate$identifier, rownames(endstate), sep = '_')
rownames(endstate) <- endstate$id
endstate <- endstate[colnames(cm),]
aframe <- data.frame(endstate = endstate$end_points, pt = t)
fit <- gam(endstate ~ lo(pt), data = aframe)
preds <- predict(fit, data.frame(pt = seq(min(t), max(t), length = 50)))
anno_col = data.frame(endstate = preds)
rownames(anno_col) <- colnames(tmp)
anno_color <- list(endstate = viridis(50))

peaks <- apply(tmp, 1, function(x) which(x == max(x)))
rowOrd <- order(peaks)
tmp <- tmp[rowOrd,]
rownames(anno_row) <- rownames(tmp)

# Generate Fig 5d ####
pheatmap::pheatmap(tmp, show_colnames = F, scale = 'row', cluster_cols = F, cluster_rows = F,
                   breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue, show_rownames = F,
                   annotation_col = anno_col, annotation_colors = anno_color)

