setwd("Cell-type analysis")

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)


####################################################################################################################
##                                     #  Create initial cell-type object   #                                     ##
####################################################################################################################

# Required RAM: 12GB
source("Cell-type analysis/utils/subset.cell.type.R")
obj <- subset.cell.type("endo")
obj@misc$graph.path <- "vascular.niche/graphs"
lapply(c("vascular.niche", "vascular.niche/graphs", "vascular.niche/data"), dir.create)


# Basic analysis to visualize object doublets
# Required RAM: 12GB
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 30, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:30, verbose = F)
obj <- FindClusters(obj, resolution = .2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:30, verbose = F)
InitialDoubletsStatus(obj)
SaveH5Seurat(obj, "vascular.niche/data/vascular.niche.h5Seurat")


####################################################################################################################
##                                     #  Search For Low Quality Clusters   #                                     ##
####################################################################################################################

obj <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat")
obj <- subset(obj, is.doublet == F)
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 2.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.2.5", filename=file.path(obj@misc$graph.path, "QCs1.pdf"))
SaveH5Seurat(obj, "vascular.niche/data/vascular.niche.h5Seurat", overwrite = T, verbose = F)


# Remove low-quality/doublet clusters
# 18,24,26,31 - low-qual/multiples; 38 - Exc-Ast-Endom; 25 - multiples; 30 - Endo-Olig
obj <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[obj$percent.mt <= 10 & (! obj$SCT_snn_res.2.5 %in% c(18,24,25,26,30,31,38))])
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 2.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.2.5", filename=file.path(obj@misc$graph.path, "QCs2.pdf"))
SaveH5Seurat(obj, "vascular.niche/data/vascular.niche.h5Seurat", overwrite = T, verbose = F)



obj <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat")
obj <- FindClusters(obj, resolution = c(.3, .5, .8, 1.2, 2.5), algorithm = 4, method = "igraph", verbose = F)
PlotClusteringResolutions(obj)
SaveH5Seurat(obj, "vascular.niche/data/vascular.niche.h5Seurat", overwrite = T, verbose = F)



####################################################################################################################
##                                #  Sub-clustering Analysis of Vascular Subset    #                              ##
####################################################################################################################
main <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat")

# ------------------------------------------------ #
#         Search For Low Quality Clusters          #
# ------------------------------------------------ #
obj <- subset(main, SCT_snn_res.0.3 %in% c(4,9,10,11,12))
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1", filename="vascular.niche/graphs/vascular.QCs1.pdf")
SaveH5Seurat(obj, "vascular.niche/data/vascular.h5Seurat", overwrite = T, verbose = F)
rm(main)

obj <- LoadH5Seurat("vascular.niche/data/vascular.h5Seurat")
obj <- FindClusters(obj, resolution = c(.1, .3, .5, .7, 1), algorithm = 4, method = "igraph", verbose = F)
PlotClusteringResolutions(obj, path="vascular.niche/graphs/vascular.clustering.resolutions.pdf")
SaveH5Seurat(obj, "vascular.niche/data/vascular.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#               Determine clustering               #
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/vascular.h5Seurat")
resolutions <- c("0.3", "0.5", "0.7")

for (res in resolutions) {
  # Building cluster tree
  Idents(obj) <- obj@meta.data[paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)
  
  # Differential expression (DE) + pathway analysis (PA)
  de <- FindMarkersWrapper(obj, test.use="negbinom")
  pa <- EnrichmentAnalysis(pa %>% dplyr::filter(avg_log2FC > 0))
  
  saveRDS(de, paste0("vascular.niche/data/vascular.de.", res, ".rds"))
  saveRDS(pa, paste0("vascular.niche/data/vascular.pa.", res, ".rds"))

  # Preliminary trait associations functionality
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res), filename= paste0("vascular.niche/graphs/vascular.trait.association.", ident, ".pdf"))
}
SaveH5Seurat(obj, "vascular.niche/data/vascular.h5Seurat", overwrite = T, verbose = F)



# Plot top differential expressed genes as well as expression of known signatures
pdf("vascular.niche/graphs/vascular.signatures.pdf", width=24, height = 10)

# UMAP of clustering resolutions *after* reordering according to cluster tree
DimPlot(obj, group.by = paste0("SCT_snn_res.", resolutions), label=T, raster=T) & NoAxes() & NoLegend()

source("utils/signatures.R")
for (res in resolutions) {
  de <- readRDS(paste0("vascular.niche/data/vascular.de.", res, ".rds"))
  
  print(DotPlot(obj, 
                group.by = paste0("SCT_snn_res.", res), 
                features = de %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC) %>% pull(gene) %>% unique) + 
          scale_color_viridis_c() + theme(legend.position = "bottom") + 
          RotatedAxis() + btheme +
          labs(x=NULL, y=NULL, title=paste0("Differential Expression - ", paste0("SCT_snn_res.", res))))
  
  for(publication in c("Garcia F. J. et al., (2021)", "Yang A. C. et al., (2021)", "Junyi H. et al., (2020)")) {
    print(DotPlot(obj, 
                  features=signatures$Perivascular[[publication]],
                  group.by = paste0("SCT_snn_res.", res)) + 
            RotatedAxis() + 
            scale_color_viridis_c() + 
            theme(legend.position = "bottom") + 
            labs(X=NULL, y=NULL, title=paste0(publication, " - Cluster Resolution:", res)))
  }
}
while (!is.null(dev.list()))  dev.off()


# ------------------------------------------------ #
#        Atlas Identities of Sub-Populations       #
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/vascular.h5Seurat")
Idents(obj) <- obj$v1.base <- plyr::mapvalues(obj$SCT_snn_res.0.7, from=c(4,6,8), to=c(3,5,7))
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
saveRDS(Tool(obj, "BuildClusterTree"), file.path(dirname(obj@misc$data.path), paste0("tree.v1.idents.rds")))
obj$v1.base <- Idents(obj)

obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "Endothelial",
                            class     = "Vascular Niche",
                            cell.type = "Endothelial",
                            state     = case_when(v1.base == 1 ~ "Arteriole",
                                                  v1.base == 2 ~ "Venule",
                                                  v1.base == 3 ~ "End.1",
                                                  v1.base == 4 ~ "End.2",
                                                  v1.base == 5 ~ "End.3",
                                                  v1.base == 6 ~ "End.4",
                                                  v1.base == 7 ~ "End.5"),
                            sub.population = case_when(v1.base == 4 & SCT_snn_res.0.7 == 5 ~ "End.2.1",
                                                       v1.base == 4 & SCT_snn_res.0.7 != 5 ~ "End.2.2",
                                                       v1.base == 5 & SCT_snn_res.0.7 == 7 ~ "End.3.1",
                                                       v1.base == 5 & SCT_snn_res.0.7 != 7 ~ "End.3.2",
                                                       v1.base == 3 & SCT_snn_res.0.7 == 3 ~ "End.1.1",
                                                       v1.base == 3 & SCT_snn_res.0.7 != 3 ~ "End.1.2"),
                            annotation = case_when(v1.base == 6 ~ "Pericyte-like")) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))

Idents(obj) <- obj$state
SaveH5Seurat(obj, "vascular.niche/data/vascular.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#                  Compute DE Genes                # 
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/vascular.h5Seurat")
de <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 4)
saveRDS(de, "vascular.niche/data/vascular.de.rds")

de.pairwise <- do.call(rbind, combn(unique(Idents(obj)), 2, simplify = F, FUN = function(pair) {
  pair = as.character(pair)
  message(paste(pair, collapse = " - "))
  comps = c(paste(pair, collapse = " vs. "), paste(rev(pair), collapse = " vs. "))  
  FindMarkers(obj, ident.1 = pair[[1]], ident.2 = pair[[2]], test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = TRUE) %>%
    mutate(comparison = if_else(avg_log2FC >= 0, comps[[1]], comps[[2]]),
           cluster = if_else(avg_log2FC >=0, pair[[1]], pair[[2]]),
           avg_log2FC = abs(avg_log2FC),
           gene = rownames(.),
           id=GeneIdMapping()$ids[gene])
}))
saveRDS(de.pairwise, "vascular.niche/data/vascular.de.pairwise.rds")



# Required RAM: 4GB
obj <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat", assays=list(SCT=c("data")), misc=T, graphs=F, neighbors=F, reductions=F, verbose=F)

saveRDS(EnrichmentAnalysis(readRDS("vascular.niche/data/vascular.de.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           cluster = gsub("\\.","_", cluster)),
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "vascular.niche/data/vascular.pa.rds")

saveRDS(EnrichmentAnalysis(readRDS("vascular.niche/data/vascular.de.pairwise.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           comparison = gsub("\\.", "_", gsub(" vs\\. ", "_vs_", comparison))), 
                           formula="id~comparison+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "vascular.niche/data/vascular.pa.pairwise.rds")


####################################################################################################################
##                                  #  Sub-clustering Analysis of Mural Subset    #                               ##
####################################################################################################################
main <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat")

# ------------------------------------------------ #
#         Search For Low Quality Clusters          #
# ------------------------------------------------ #
obj <- subset(main, SCT_snn_res.0.3 %in% c(6,7))
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1", filename="vascular.niche/graphs/mural.QCs1.pdf")
SaveH5Seurat(obj, "vascular.niche/data/mural.h5Seurat", overwrite = T, verbose = F)
rm(main)


obj <- LoadH5Seurat("vascular.niche/data/mural.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[obj$nCount_RNA <= 20000 & obj$nFeature_RNA <= 5000 & obj$SCT_snn_res.1 != 7])
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = c(.3, .5, .6, .8, 1), algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1", filename="vascular.niche/graphs/mural.QCs2.pdf")
PlotClusteringResolutions(obj, path="vascular.niche/graphs/mural.clustering.resolutions.pdf")
SaveH5Seurat(obj, "vascular.niche/data/mural.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#               Determine clustering               #
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/mural.h5Seurat")
resolutions <- c("0.3", "0.6", "0.8", "1")

for (res in resolutions) {
  # Building cluster tree
  Idents(obj) <- obj@meta.data[paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)
  
  # Differential expression (DE) + pathway analysis (PA)
  de <- FindMarkersWrapper(obj, test.use="negbinom")
  pa <- EnrichmentAnalysis(pa %>% dplyr::filter(avg_log2FC > 0))
  
  saveRDS(de, paste0("vascular.niche/data/mural.de.", res, ".rds"))
  saveRDS(pa, paste0("vascular.niche/data/mural.pa.", res, ".rds"))
  
  # Preliminary trait associations functionality
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res), filename= paste0("vascular.niche/graphs/mural.trait.association.", ident, ".pdf"))
}
SaveH5Seurat(obj, "vascular.niche/data/mural.h5Seurat", overwrite = T, verbose = F)



# Plot top differential expressed genes as well as expression of known signatures
pdf("vascular.niche/graphs/mural.signatures.pdf", width=24, height = 10)

# UMAP of clustering resolutions *after* reordering according to cluster tree
DimPlot(obj, group.by = paste0("SCT_snn_res.", resolutions), label=T, raster=T) & NoAxes() & NoLegend()

source("utils/signatures.R")
for (res in resolutions) {
  de <- readRDS(paste0("vascular.niche/data/mural.de.", res, ".rds"))
  
  print(DotPlot(obj, 
                group.by = paste0("SCT_snn_res.", res), 
                features = de %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC) %>% pull(gene) %>% unique) + 
          scale_color_viridis_c() + theme(legend.position = "bottom") + 
          RotatedAxis() + btheme +
          labs(x=NULL, y=NULL, title=paste0("Differential Expression - ", paste0("SCT_snn_res.", res))))
  
  for(publication in c("Garcia F. J. et al., (2021)", "Yang A. C. et al., (2021)", "Junyi H. et al., (2020)")) {
    print(DotPlot(obj, 
                  features=signatures$Perivascular[[publication]],
                  group.by = paste0("SCT_snn_res.", res)) + 
            RotatedAxis() + 
            scale_color_viridis_c() + 
            theme(legend.position = "bottom") + 
            labs(X=NULL, y=NULL, title=paste0(publication, " - Cluster Resolution:", res)))
  }
}
while (!is.null(dev.list()))  dev.off()



# ------------------------------------------------ #
#        Atlas Identities of Sub-Populations       #
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/vascular.h5Seurat")

obj$v1.idents <- plyr::mapvalues(as.character(obj$SCT_snn_res.0.6), from=c(2,3,4), to=c(2,2,2))
Idents(obj) <- obj$v1.idents
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
obj$v1.idents <- Idents(obj)

obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "Endothelial",
                            class     = "Vascular Niche",
                            cell.type = "Mural",
                            state     = plyr::mapvalues(v1.idents, from=1:4, to=c("Peri.1", "Peri.2", "SMC.1", "SMC.2")),
                            sub.population = state,
                            annotation = plyr::mapvalues(v1.idents, from=1:4, to=c("M-Pericyte", "T-Pericyte", "aSMC", "vSMC"))) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))

Idents(obj) <- obj$state
SaveH5Seurat(obj, "vascular.niche/data/mural.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#                  Compute DE Genes                # 
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/mural.h5Seurat")
de <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 4)
saveRDS(de, "vascular.niche/data/mural.de.rds")

de.pariwise <- do.call(rbind, combn(unique(Idents(obj)), 2, simplify = F, FUN = function(pair) {
  pair = as.character(pair)
  message(paste(pair, collapse = " - "))
  comps = c(paste(pair, collapse = " vs. "), paste(rev(pair), collapse = " vs. "))  
  FindMarkers(obj, ident.1 = pair[[1]], ident.2 = pair[[2]], test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = TRUE) %>%
    mutate(comparison = if_else(avg_log2FC >= 0, comps[[1]], comps[[2]]),
           cluster = if_else(avg_log2FC >=0, pair[[1]], pair[[2]]),
           avg_log2FC = abs(avg_log2FC),
           gene = rownames(.),
           id=GeneIdMapping()$ids[gene])
}))
saveRDS(de.pairwise, "vascular.niche/data/mural.de.pairwise.rds")


# Required RAM: 3GB
obj <- LoadH5Seurat("vascular.niche/data/mural.h5Seurat", assays=list(SCT=c("data")), misc=T, graphs=F, neighbors=F, reductions=F,verbose=T)

saveRDS(EnrichmentAnalysis(readRDS("vascular.niche/data/mural.de.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           cluster = gsub("\\.","_", cluster)),
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "vascular.niche/data/mural.pa.rds")

saveRDS(EnrichmentAnalysis(readRDS("vascular.niche/data/mural.de.pairwise.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           comparison = gsub("\\.", "_", gsub(" vs\\. ", "_vs_", comparison))), 
                           formula="id~comparison+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "vascular.niche/data/mural.pa.pairwise.rds")


####################################################################################################################
##                               #  Sub-clustering Analysis of Fibroblast Subset    #                             ##
####################################################################################################################
main <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat")

# ------------------------------------------------ #
#         Search For Low Quality Clusters          #
# ------------------------------------------------ #
obj <- subset(main, SCT_snn_res.0.3 %in% c(1,3))
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1", filename="vascular.niche/graphs/fibroblast.QCs1.pdf")
SaveH5Seurat(obj, "vascular.niche/data/fibroblast.h5Seurat")
rm(main)

# Remove multiplet cluster 11
obj <- LoadH5Seurat("vascular.niche/data/fibroblast.h5Seurat")
obj <- subset(obj, SCT_snn_res.1 != 11)
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = c(.5, 1), algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1", filename="vascular.niche/graphs/fibroblast.QCs2.pdf")
SaveH5Seurat(obj, "vascular.niche/data/fibroblast.h5Seurat", overwrite = T, verbose = F)


# Plot top differential expressed genes as well as expression of known signatures
pdf("vascular.niche/graphs/fibroblast.signatures.pdf", width=24, height = 10)

# UMAP of clustering resolutions *after* reordering according to cluster tree
DimPlot(obj, group.by = paste0("SCT_snn_res.", resolutions), label=T, raster=T) & NoAxes() & NoLegend()

source("utils/signatures.R")
for (res in resolutions) {
  de <- readRDS(paste0("microglia/data/de.", res, ".rds"))
  
  print(DotPlot(obj, 
                group.by = paste0("SCT_snn_res.", res), 
                features = de %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC) %>% pull(gene) %>% unique) + 
          scale_color_viridis_c() + theme(legend.position = "bottom") + 
          RotatedAxis() + btheme +
          labs(x=NULL, y=NULL, title=paste0("Differential Expression - ", paste0("SCT_snn_res.", res))))
  
  for(publication in c("Garcia F. J. et al., (2021)", "Yang A. C. et al., (2021)","Junyi H. et al., (2020)")) {
    print(DotPlot(obj, 
                  features=signatures$Perivascular[[publication]],
                  group.by = paste0("SCT_snn_res.", res)) + 
            RotatedAxis() + 
            scale_color_viridis_c() + 
            theme(legend.position = "bottom") + 
            labs(X=NULL, y=NULL, title=paste0(publication, " - Cluster Resolution:", res)))
  }
}
while (!is.null(dev.list()))  dev.off()



# ------------------------------------------------ #
#        Atlas Identities of Sub-Populations       #
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/fibroblast.h5Seurat")
Idents(obj) <- obj$v1.base <- obj$SCT_snn_res.0.5
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
obj$v1.base <- Idents(obj)

obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "Endothelial",
                            class     = "Vascular Niche",
                            cell.type = "Fibroblast",
                            state     = paste0("Fib.", dplyr::recode(v1.base, "3"="3", "4"="3","5"="3","6"="3")),
                            sub.population = state,
                            annotation = case_when(v1.base == 1 ~ "Meningeal Fibroblast",
                                                   v1.base != 1 ~ "Perivascular Fibroblast-like")) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))
Idents(obj) <- obj$state
SaveH5Seurat(obj, "vascular.niche/data/fibroblast.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#                  Compute DE Genes                # 
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/fibroblast.h5Seurat")
de <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 4)
saveRDS(de, "vascular.niche/data/fibroblast.de.rds")

de.pariwise <- do.call(rbind, combn(unique(Idents(obj)), 2, simplify = F, FUN = function(pair) {
  pair = as.character(pair)
  message(paste(pair, collapse = " - "))
  comps = c(paste(pair, collapse = " vs. "), paste(rev(pair), collapse = " vs. "))  
  FindMarkers(obj, ident.1 = pair[[1]], ident.2 = pair[[2]], test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = TRUE) %>%
    mutate(comparison = if_else(avg_log2FC >= 0, comps[[1]], comps[[2]]),
           cluster = if_else(avg_log2FC >=0, pair[[1]], pair[[2]]),
           avg_log2FC = abs(avg_log2FC),
           gene = rownames(.),
           id=GeneIdMapping()$ids[gene])
}))
saveRDS(de.pairwise, "vascular.niche/data/fibroblast.de.pairwise.rds")


# Required RAM: 3GB
obj <- LoadH5Seurat("vascular.niche/data/fibroblast.h5Seurat", assays=list(SCT=c("data")), misc=T, graphs=F, neighbors=F, reductions=F,verbose=T)

saveRDS(EnrichmentAnalysis(readRDS("vascular.niche/data/fibroblast.de.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           cluster = gsub("\\.","_", cluster)),
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "vascular.niche/data/fibroblast.pa.rds")

saveRDS(EnrichmentAnalysis(readRDS("vascular.niche/data/fibroblast.de.pairwise.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           comparison = gsub("\\.", "_", gsub(" vs\\. ", "_vs_", comparison))), 
                           formula="id~comparison+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "vascular.niche/data/fibroblast.pa.pairwise.rds")



####################################################################################################################
##                                 #  Sub-clustering Analysis of Immune Subset    #                               ##
####################################################################################################################
main <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat")

# ------------------------------------------------ #
#         Search For Low Quality Clusters          #
# ------------------------------------------------ #
obj <- subset(main, SCT_snn_res.0.3 %in% c(2,5,8))
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = c(1,1.4), algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1", filename="vascular.niche/graphs/immune.QCs1.pdf")
SaveH5Seurat(obj, "vascular.niche/data/immune.h5Seurat")
rm(main)


# Plot expression of known signatures
pdf("vascular.niche/graphs/immune.signatures.pdf", width=24, height = 10)
DimPlot(obj, group.by = "SCT_snn_res.1.4", label=T, raster=T) & NoAxes() & NoLegend()

source("utils/signatures.R")
features = list(unlist(signatures$Perivascular$`Garcia F. J. et al., (2021)`, recursive = F),
                unlist(signatures$Perivascular$`Yang A. C. et al., (2021)`, recursive = F),
                signatures$Perivascular$`Junyi H. et al., (2020)`,
                signatures$Machrophages$`Milich L. M. et al., (2021)`,
                signatures$Microglia$`Milich L. M. et al., (2021)`,
                signatures$`Cell Types`$`Milich L. M. et al., (2021)`)
titles = c("Garcia F.J. et al., (2021)", "Yang A. C. et al., (2021)", "Junyi H. et al., (2020)", 
           "Macrophages - Milich L. M. et al., (2021)", "Microglia - Milich L. M. et al., (2021)", "General - Milich L. M. et al., (2021)")

for(i in seq_along(features)) {
  f = setNames(features[[i]], gsub("\\.","\n", names(features[[i]])))
  print(DotPlot(obj, features=f, group.by = "SCT_snn_res.1.4", assay = "RNA") + RotatedAxis() + 
          scale_color_viridis_c() + theme(legend.position = "bottom", strip.text.x = element_text(size=7)) + 
          btheme + labs(x=NULL, y=NULL, title=paste(titles[[i]], "- Cluster Resolution: SCT_snn_res.1.4")))
}
while (!is.null(dev.list()))  dev.off()



# ------------------------------------------------ #
#        Atlas Identities of Sub-Populations       #
# ------------------------------------------------ #
obj <- LoadH5Seurat("vascular.niche/data/immune.h5Seurat")
Idents(obj) <- obj$v1.base <- obj$SCT_snn_res.1.4
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
obj$v1.base <- Idents(obj)

obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     dplyr::mutate(subset    = "Endothelial",
                                   class     = "Immune",
                                   cell.type = case_when(v1.base == 1 ~ "NK Cells",
                                                         v1.base == 2 ~ "CD8+ T Cells",
                                                         v1.base == 3 ~ "Macrophages",
                                                         v1.base == 8 ~ "Neutrophils",
                                                         T ~ "Erythrocytes"),
                                   state = cell.type,
                                   sub.population = cell.type,
                                   annotation = cell.type) %>%
                     dplyr::select(subset, class, cell.type, state, sub.population, annotation, v1.base))
Idents(obj) <- obj$state
SaveH5Seurat(obj, "vascular.niche/data/immune.h5Seurat", overwrite = T, verbose = F)



####################################################################################################################
##                                       #  Merge endo branches and analysis    #                                 ##
####################################################################################################################

obj <- LoadH5Seurat("vascular.niche/data/vascular.niche.h5Seurat")

df <- lapply(paste0("vascular.niche/data/", c("vascular", "mural", "fibroblast","immune"), ".h5Seurat"), function(o)
  LoadH5Seurat(o, assays=list(SCT=c("data")), misc=F, graphs=F, reductions=F, verbose=F)@meta.data %>% 
    dplyr::select(subset, class, cell.type, state, sub.population, annotation, v1.base)) %>%
  do.call(rbind, .)

obj <- subset(obj, cells = rownames(df))
obj <- FindNeighbors(obj, dims=1:50)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)

obj <- AddMetaData(obj, df)
pdf("vascular.niche/graphs/final.idents.pdf", width=9, height = 6)
for(g in c("class","cell.type","state","sub.population", "annotation")) {
  print(LabelClusters(DimPlot(obj, label=F, pt.size=.1, raster=T, group.by = g), id = g,
                      position = "median", box = T, repel=T, size=5, box.padding=1, min.segment.length=1, fill = rgb(1, 1, 1, alpha = 0.8),colour = '#0C0C0C') + NoAxes() + labs(title=NULL))
}
while (!is.null(dev.list()))  dev.off()
SaveH5Seurat(obj, "vascular.niche/data/vascular.h5Seurat", overwrite = T, verbose = F)

