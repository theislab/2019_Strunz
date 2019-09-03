# R code
# Lukas Simon
# This code will generate Strun et al Fig 1j

# Load R libs ####
library(readxl)
library(igraph)
library(plotrix)

# Read receptor-ligand table ####
merged <- read_excel("../data/Table_S6_(Receptor_ligand_pairs).xlsx", sheet = 1)
merged <- merged[-which(merged$cluster.lig %in% c("plasma_cells", "unassigned")),]
merged <- merged[-which(merged$cluster.rec %in% c("plasma_cells", "unassigned")),]

# Add metacelltype groups ####
assignments <- list(
  T_cells = c("t_cells"),
  B_cells = c("b_cells"),
  NK_cells = "nk_cells",
  Dendritic_cells = c("ccl17_dc", "cd103_dc", "dc"),
  Monocytes = c("nonclassical_mo", "classical_mo"),
  Granulocytes = "granulocytes",
  Macrophages = c("am", "fn1_mcpg", "mono_im", "im", "cd163_im", 14, 19, 29),
  Endothelial_cells = c("cec", "lec", "vec", "vcam1_vec"),
  Mesothelial_cells = c("mesothelial_cells"),
  Fibroblasts = c("fibroblasts"),
  Smooth_muscle_cells = "smooth_muscle_cells",
  Alveolar_epithelium = c("alv_epithelium"),
  Club_cells = "club_cells",
  Goblet_cells = "goblet_cells",
  Ciliated_cells = c("ciliated_cells")
)
merged$metacelltype.lig <- NA
merged$metacelltype.rec <- NA
lapply(names(assignments), function(x){
  print(x)
  tmp <- intersect(assignments[[x]], merged$cluster.lig)
  merged$metacelltype.lig[which(merged$cluster.lig %in% tmp)] <<- x
  tmp <- intersect(assignments[[x]], merged$cluster.rec)
  merged$metacelltype.rec[which(merged$cluster.rec %in% tmp)] <<- x
})
   
# Generate adjacency table of all pairs ####
tmp <- table(paste(merged$metacelltype.lig, merged$metacelltype.rec, sep = "|"))
tmp <- data.frame(do.call(rbind, lapply(names(tmp),function(x) strsplit(x, '|', fixed = T)[[1]])), as.numeric(tmp))
colnames(tmp) <- c('receptor', 'ligand', 'freq')
tmp$receptor <- as.character(tmp$receptor)
tmp$ligand <- as.character(tmp$ligand)
celltypes <- setdiff(unique(c(tmp[,1], tmp[,2])), "NA")
matr <- matrix(0, length(celltypes), length(celltypes))
colnames(matr) <- rownames(matr) <- celltypes
lapply(1:nrow(tmp), function(x){
  try(matr[tmp[x,1], tmp[x,2]] <<- tmp[x, 3]  )
})
adj_all <- matr

# Generate plots ####
node_size <- 10
farben <- "lightblue"
vertex_order <- c("Ciliated_cells", "Club_cells", "Goblet_cells", "Alveolar_epithelium",
                  "B_cells", "Dendritic_cells", "Granulocytes", "Macrophages", "Monocytes", "NK_cells", "T_cells",
                  "Mesothelial_cells", "Smooth_muscle_cells", "Fibroblasts", "Endothelial_cells")

G <- as.undirected(graph.adjacency(adj_all, weighted = TRUE, mode = "undirected", diag = F))
E(G)$width <- E(G)$weight/20
E(G)$color <- color.scale(E(G)$weight, extremes = c("grey", "red"))
coords <- layout_in_circle(G, order = vertex_order)
plot(G, edge.width = E(G)$width, layout = coords,
     vertex.size = node_size, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.color = farben)
