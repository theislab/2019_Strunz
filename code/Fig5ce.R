# load R libs ####
library(plotrix)
library(mgcv)
library(slingshot)
library(SingleCellExperiment)
library(mclust)
library(dbscan)
library(RColorBrewer)
library(enrichR)
library(gam)

# load data ####
load("../data/EpiHiRes_seurat.RData")
louvain <- read.table('../data/EpiHiRes_louvain_labels.txt', row.names = 1)
load("../data/Convergence_trajectory.RData")
dm <- convergence[,1:2]

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

# Subset to genes expressed in at least 10 cells in more than 5 mice ####
cm <- assays(sce)$counts
md <- md[colnames(cm),]

plot_traj_single_gene <- function(gene = 'Krt8'){
  
  scale <- function(x) (x - min(x))/(max(x) - min(x))
  
  expr <- scale(data.matrix(log(cm[gene,] + 1)))
  
  aframe <- data.frame(expr, pt = t)
  
  melted <- melt(data = aframe, id.vars = "pt", measure.vars = setdiff(colnames(aframe), 'pt'))
  
  ggplot(melted, aes(y = value, x = pt, group = variable, color = variable)) + geom_smooth(method = "loess") +
    ylab(paste("Relative expression")) + xlab("Pseudotime") + ggtitle(gene)
}

# Generate panel c ####
plot_traj_single_gene("Krt8")
plot_traj_single_gene("Scgb1a1")
plot_traj_single_gene("Sftpc")

# Generate panel e ####
plot_traj_single_gene("Nkx2-1")
plot_traj_single_gene("Cebpa")
plot_traj_single_gene("Foxp2")

plot_traj_single_gene("Etv5")
plot_traj_single_gene("Epas1")
plot_traj_single_gene("Yeats4")

plot_traj_single_gene("Gata6")
plot_traj_single_gene("Hint1")
plot_traj_single_gene("Kdm5c")

plot_traj_single_gene("Rbpms")
plot_traj_single_gene("Nupr1")
plot_traj_single_gene("Sox4")

plot_traj_single_gene("Sra1")
plot_traj_single_gene("Nfkb2")
plot_traj_single_gene("Nfkb1")
