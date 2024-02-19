setwd("2. Cell-type analysis")

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
suppressPackageStartupMessages(source("utils/analysis.reports.R"))

####################################################################################################################
##                                     #  Create initial cell-type object   #                                     ##
####################################################################################################################

# Required RAM: 20GB
source("utils/subset.cell.type.R")
obj <- subset.cell.type("microglia")
obj@misc$graph.path <- file.path(name, "graphs")
lapply(c("microglia", "microglia/graphs", "microglia/data"), dir.create)


# Basic analysis to visualize object doublets 
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 30, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:30, verbose = F)
obj <- FindClusters(obj, resolution = .2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:30, verbose = F)
InitialDoubletsStatus(obj)
SaveH5Seurat(obj, "microglia/data/microglia.h5Seurat")


####################################################################################################################
##                                     #  Search For Low Quality Clusters   #                                     ##
####################################################################################################################
obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat")
obj <- subset(obj, is.doublet == F)
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F) %>%
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs1.pdf"))
SaveH5Seurat(obj, "microglia/data/microglia.h5Seurat", overwrite = T, verbose = F)


# Required RAM: 40GB
# 34,35,37,33,29,14,11 - multiples
obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[obj$percent.mt <= 10 & !obj$SCT_snn_res.1.5 %in% c(11,14,29,33,34,35,37)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs2.pdf"))
SaveH5Seurat(obj, "microglia/data/microglia.h5Seurat", overwrite = T, verbose = F)


# Required RAM: 42GB
#Clusters: #30 - High MT, #33 - Low-qual
obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[!obj$SCT_snn_res.1.5 %in% c(30, 33)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs3.pdf"))
SaveH5Seurat(obj, "microglia/data/microglia.h5Seurat", overwrite = T, verbose = F)


obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat")
obj$SCT_snn_res.1.5 <- NULL
obj <- FindClusters(obj, resolution = c(.1,.3,.5,.8,1.5), algorithm = 4, method = "igraph", verbose = F)
PlotClusteringResolutions(obj)
SaveH5Seurat(obj, "microglia/data/microglia.h5Seurat", overwrite = T, verbose = F)




####################################################################################################################
##                                            #  Determine Clustering   #                                         ##
####################################################################################################################
# Investigated resolutions
resolutions <- c("0.3", "0.5", "0.8")

# Required RAM: Clustering Trees: 18BG; DE+PA: `0.3`=25GB, `0.5`=28GB, `0.8`=28GB
obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat")
for (res in resolutions) {
  # Building cluster tree
  Idents(obj) <- obj@meta.data[, paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)
  
  # Differential expression (DE) + pathway analysis (PA)
  de <- obj@misc[[paste0("de.", res)]] <- FindMarkersWrapper(obj, test.use="negbinom")
  pa <- EnrichmentAnalysis(de %>% dplyr::filter(avg_log2FC > 0))
  
  saveRDS(de, paste0("microglia/data/de.", res, ".rds"))
  saveRDS(pa, paste0("microglia/data/pa.", res, ".rds"))
  
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res))
}
SaveH5Seurat(obj, "microglia/data/microglia.h5Seurat", overwrite = T, verbose = F)

# Plot top differential expressed genes as well as expression of known signatures
pdf(file.path(obj@misc$graph.path, "signatures.pdf"), width=24, height = 10)

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
  
  for(publication in c("Cain A. et al., (2021)", "Gerrits E. et al., (2021)")) {
    print(DotPlot(obj, 
                  features=signatures$Microglia[[publication]],
                  group.by = paste0("SCT_snn_res.", res)) + 
            RotatedAxis() + 
            scale_color_viridis_c() + 
            theme(legend.position = "bottom") + 
            labs(X=NULL, y=NULL, title=paste0(publication, " - Cluster Resolution:", res)))
  }
}
while (!is.null(dev.list()))  dev.off()



####################################################################################################################
##                                 #  Subletting and refining microglial clusters  #                              ##
####################################################################################################################

obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat")
for (res in c("0.3", "0.5", "0.8"))
  main[[paste0("main.clusters.res.", res)]] <- main[[paste0("SCT_snn_res.", res)]]
main@meta.data <- main@meta.data[, !grepl("SCT_snn_res", colnames(main@meta.data))]

# ------------------------------------------------ #
#                  Redox Microglia                 #
# ------------------------------------------------ #
obj <- subset(main, main.clusters.res.0.5 == 7)
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = c(.1, .2, .3), algorithm = 4, verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.0.3", filename=file.path(main@misc$graph.path, "redox.QCs.pdf"))
PlotClusteringResolutions(obj, path=file.path(main@misc$graph.path, "redox.clustering.resolutions.pdf"))


obj$merged.idents.0.2 <- plyr::mapvalues(obj$SCT_snn_res.0.2, from = c(4,5,6), to=c(2,2,1))
obj$merged.idents.0.3 <- plyr::mapvalues(obj$SCT_snn_res.0.3, from = c(7,8,9,10), to=c(2,3,1,1))
obj$less.merged.idents.0.3 <- plyr::mapvalues(obj$SCT_snn_res.0.3, from = c(7,9), to=c(2,1))

# Differential expression (DE) + pathway analysis (PA) for investigated resolutions
idents <- c("merged.idents.0.2", "merged.idents.0.3", "less.merged.idents.0.3")
for(ident in idents) {
  Idents(obj) <- obj@meta.data[,ident]
  obj@misc[[paste0("de.", ident)]] <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 2)
  obj@misc[[paste0("pa.", ident)]] <- EnrichmentAnalysis(obj@misc[[paste0("de.", ident)]] %>% dplyr::filter(avg_log2FC > 0))
}
saveRDS(obj, "microglia/data/microglia.redox.rds")

# Plot top differential expressed genes as well as expression of known signatures
pdf(file.path(main@misc$graph.path,"redox.signatures.pdf"), width=24, height = 6)
DimPlot(obj, group.by = idents, label=T, raster=T, ncol = 4) & NoAxes() & NoLegend()

for(ident in idents) {
  print(DotPlot(obj, 
                group.by = ident, 
                features = unique((obj@misc[[paste0("de.", ident)]] %>% group_by(cluster) %>% top_n(25, wt=avg_log2FC))$gene)) + 
          scale_color_viridis_c() + 
          theme(legend.position = "bottom") + 
          RotatedAxis() + 
          btheme +
          labs(x=NULL, y=NULL, title=paste0("Differential Expression - ", ident)))
}
while (!is.null(dev.list()))  dev.off()



# ------------------------------------------------ #
#              Disease-like Microglia              #
# ------------------------------------------------ #
obj <- subset(main, main.clusters.res.0.5 %in% c(3,4,13))
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = c(.2,.3,.4,.5), algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
PlotClusteringResolutions(obj, path=file.path(main@misc$graph.path, "redox.clustering.resolutions.pdf"))

# Differential expression (DE) + pathway analysis (PA) for investigated resolutions
idents <- paste0("SCT_snn_res.0.", c(2,3,4))
for(ident in idents) {
  Idents(obj) <- obj@meta.data[,ident]
  obj@misc[[paste0("de.", ident)]] <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 2)
  obj@misc[[paste0("pa.", ident)]] <- EnrichmentAnalysis(obj@misc[[paste0("de.", ident)]] %>% dplyr::filter(avg_log2FC > 0))
}

# Plot top differential expressed genes as well as expression of known signatures
pdf(file.path(main@misc$graph.path,"disease.like.signatures.pdf"), width=24, height = 6)
DimPlot(obj, group.by = idents, label=T, raster=T, ncol = 4) & NoAxes() & NoLegend()

for(ident in idents) {
  print(DotPlot(obj, 
                group.by = ident, 
                features = unique((obj@misc[[paste0("de.", ident)]] %>% group_by(cluster) %>% top_n(25, wt=avg_log2FC))$gene)) + 
          scale_color_viridis_c() + 
          theme(legend.position = "bottom") + 
          RotatedAxis() + 
          btheme +
          labs(x=NULL, y=NULL, title=paste0("Differential Expression - ", ident)))
}
while (!is.null(dev.list()))  dev.off()

saveRDS("microglia/data/microglia.disease.like.rds")



####################################################################################################################
##                                    #  Atlas Identities of Sub-populations    #                                 ##
####################################################################################################################
obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat")

# Create base identities
Idents(obj) <- obj$v1.base <- plyr::mapvalues(obj$SCT_snn_res.0.5, from=c(5,8), to=c(18,9))
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
obj$v1.base <- Idents(obj)

redox        <- saveRDS("microglia/data/microglia.redox.rds")
disease.like <- saveRDS("microglia/data/microglia.disease.like.rds")

obj <- AddMetaData(obj,
                   obj@meta.data %>%
                     select(v1.base) %>%
                     mutate(subset    = "Microglia",
                            class     = recode(v1.base, "1" = "Immune", "3" = "Immune", .default = "Glia"),
                            cell.type = recode(v1.base, "1" = "Monocytes", "3" = "Macrophages", .default = "Microglia"),
                            state     = recode(v1.base, 
                                               "2"  = "Mic.1",
                                               # Homeostatic clusters
                                               "10" = "Mic.2",
                                               "9"  = "Mic.3",
                                               "11" = "Mic.4",
                                               "8"  = "Mic.5",
                                               # Reactive clusters
                                               "14" = "Mic.6",
                                               "16" = "Mic.7",
                                               "15" = "Mic.8",
                                               # Other clusters
                                               "6"  = "Mic.9/10",
                                               "7"  = "Mic.11",
                                               # Disease clusters
                                               "13" = "Mic.12/13",
                                               # Other clusters
                                               "12"  = "Mic.14",
                                               "4"   = "Mic.15",
                                               "5"   = "Mic.16",
                                               # Non-microglia
                                               "1"   = "Monocytes",
                                               "3"   = "Macrophages")) %>%
                     mutate(state     = case_when(state == "Mic.9/10" & rownames(.) %in% colnames(redox)[redox$SCT_snn_res.0.2 %in% c(2,3)] ~ "Mic.9", # Homeostatic-Redox
                                                  state == "Mic.9/10" & (!rownames(.) %in% colnames(redox)[redox$SCT_snn_res.0.2 %in% c(2,3)]) ~ "Mic.10", # Reactive-Redox
                                                  state == "Mic.12/13" & rownames(.) %in% colnames(disease.like)[disease.like$SCT_snn_res.0.4 %in% c(2)] ~ "Mic.13",
                                                  state == "Mic.12/13" & (!rownames(.) %in% colnames(disease.like)[disease.like$SCT_snn_res.0.4 %in% c(2)])  ~ "Mic.12",
                                                  T ~ state),
                            sub.population = NA_character_,
                            annotation = case_when(state %in% c("Monocytes","Macropgages") ~ state,
                                                   state %in% paste0("Mic.",c(2,3,4,5)) ~ "Homeostatic",
                                                   state %in% paste0("Mic.",c(6,7,8)) ~ "Reactive",
                                                   state == "Mic.1" ~ "Proliferating",
                                                   state == "Mic.9" ~ "Homeostatic-Redox",
                                                   state == "Mic.10" ~ "Reactive-Redox",
                                                   state == "Mic.11" ~ "Stress-Response",
                                                   state == "Mic.14" ~ "Interferon-Response",
                                                   state %in% c("Mic.12","Mic.13") ~ "Disease-Elevated",
                                                   T ~ state)
                     ) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))

Idents(obj) <- obj$state
SaveH5Seurat(obj, "microglia/data/microglia.h5Seurat", overwrite = T, verbose = F)
rm(redox, disease.like)


pdf(file.path(obj@misc$graph.path,"final.idents.pdf"), width=9, height = 6)
for(g in c("v1.base","class","cell.type","state", "annotation")) {
  print(LabelClusters(DimPlot(obj, label=F, pt.size=.1, raster=T, group.by = g), id = g,
                      position = "median", box = T, repel=T, size=5, box.padding=1, min.segment.length=1, 
                      fill = rgb(1, 1, 1, alpha = 0.8),colour = '#0C0C0C') + 
          NoAxes() + 
          labs(title=g))
}
while (!is.null(dev.list()))  dev.off()


# ----------------------------------------------------- #
#                     Compute DE Genes                  #
# ----------------------------------------------------- #
obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat")

# Required RAM: 104G
options(future.globals.maxSize= 5*1024^3)
obj@misc$de.v1.idents <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 4)
saveRDS(obj@misc$de.v1.idents, "microglia/data/de.rds")

obj@misc$de.v1.idents.pairwise <- do.call(rbind, combn(unique(Idents(obj)), 2, simplify = F, FUN = function(pair) {
  pair = as.character(pair)
  plan(multicore, workers = 4)
  FindMarkers(obj, ident.1 = pair[[1]], ident.2 = pair[[2]], test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = F) %>%
    mutate(comparison = paste(if(avg_log2FC >= 0) pair else rev(pair), collapse = " vs. "),
           cluster = if_else(avg_log2FC >=0, pair[[1]], pair[[2]]),
           avg_log2FC = abs(avg_log2FC),
           gene = rownames(.),
           id=GeneIdMapping()$ids[gene])
}))
saveRDS(obj@misc$de.v1.idents.pairwise, "microglia/data/de.pairwise.rds")


# Required RAM: 9GB
obj <- LoadH5Seurat("microglia/data/microglia.h5Seurat", assays=list(SCT=c("data")), misc=T, graphs=F, neighbors=F, reductions=F, verbose=F)

saveRDS(EnrichmentAnalysis(readRDS("microglia/data/de.rds") %>% 
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"), 
                                           cluster = gsub("\\.","_", cluster)), 
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "microglia/data/pa.rds")

saveRDS(EnrichmentAnalysis(readRDS("microglia/data/de.pairwise.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           comparison = gsub("\\.", "_", gsub(" vs\\. ", "_vs_", comparison))), 
                           formula="id~comparison+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "microglia/data/pa.pairwise.rds")

