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
obj <- subset.cell.type("inhibitory")
obj@misc$graph.path <- file.path(name, "graphs")
lapply(c("inhibitory", "inhibitory/graphs", "inhibitory/data"), dir.create)


# Basic analysis to visualize object doublets 
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 30, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:30, verbose = F)
obj <- FindClusters(obj, resolution = .2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:30, verbose = F)
InitialDoubletsStatus(obj)
SaveH5Seurat(obj, "inhibitory/data/inhibitory.h5Seurat")


####################################################################################################################
##                                     #  Search For Low Quality Clusters   #                                     ##
####################################################################################################################
neuronal.markers = c("MEG3", "PVALB","NPY","SST","VIP","KIT","GAD2","SLC17A7", "RORB", "TOX","FOXP2", "CUX2")

obj <- LoadH5Seurat("inhibitory/data/inhibitory.h5Seurat")
obj <- subset(obj, is.doublet == F)
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose=F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs1.pdf"))
SaveH5Seurat(obj, "inhibitory/data/inhibitory.h5Seurat", overwrite = T, verbose = F)


pdf(file.path(obj@misc$graph.path,"neuronal.markers1.pdf"), width=20, height = 20)
FeaturePlot(obj, features=neuronal.markers, order=T, raster=T) & NoAxes() & scale_color_viridis_c()
while (!is.null(dev.list()))  dev.off()


# Required RAM: 199GB
# Remove low-quality/doublet clusters
# 41,42,44,45,46 - Olig-Inh; 43 - Micr-Endo-Astr; 40 - Astr-Inh; 38 - Micr; 33 - Low-qual
obj <- LoadH5Seurat("inhibitory/data/inhibitory.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[obj$percent.mt <= 10 & !obj$SCT_snn_res.1.5 %in% c(33,38,40,41,42,43,44,45,46)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj@meta.data <- obj@meta.data[,!grepl("SCT_snn", colnames(obj@meta.data))]
obj <- FindClusters(obj, resolution = 1.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
SaveH5Seurat(obj, "inhibitory/data/inhibitory.h5Seurat", overwrite = T, verbose = F)

LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs2.pdf"))
pdf(file.path(obj@misc$graph.path,"neuronal.markers2.pdf"), width=20, height = 20)
FeaturePlot(obj, features=neuronal.markers, order=T, raster=T) & NoAxes() & scale_color_viridis_c()
while (!is.null(dev.list()))  dev.off()


# Required RAM: 218GB
obj <- LoadH5Seurat("inhibitory/data/inhibitory.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[! obj$SCT_snn_res.1.5 %in% c(12,34)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj@meta.data <- obj@meta.data[,!grepl("SCT_snn", colnames(obj@meta.data))]
obj <- FindClusters(obj, resolution = 1.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs3.pdf"))
SaveH5Seurat(obj, "inhibitory/data/inhibitory.h5Seurat", overwrite = T, verbose = F)


# Required RAM: 192GB
obj <- LoadH5Seurat("inhibitory/data/inhibitory.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[obj$percent.mt <= 10 & !obj$SCT_snn_res.1.5 %in% c(35,36,37)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj@meta.data <- obj@meta.data[,!grepl("SCT_snn", colnames(obj@meta.data))]
obj <- FindClusters(obj, resolution = c(.3, .5, .7, .9, 1.5), algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs4.pdf"))
PlotClusteringResolutions(obj)
SaveH5Seurat(obj, "inhibitory/data/inhibitory.h5Seurat", overwrite = T, verbose = F)



####################################################################################################################
##                                            #  Determine Clustering   #                                         ##
####################################################################################################################
# Investigated resolutions
resolutions <- c("0.3", "0.5", "0.9")

# Required RAM: `0.3`=421GB, `0.5`=422GB, `0.9`=423GB
obj <- LoadH5Seurat("inhibitory/data/inhibitory.h5Seurat")
for (res in resolutions) {
  # Building cluster tree
  Idents(obj) <- obj@meta.data[, paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)

  # Differential expression (DE) + pathway analysis (PA)
  options(future.globals.maxSize= 14*1024^3)
  de <- FindMarkersWrapper(obj, test.use="negbinom", workers = 3)
  pa <- EnrichmentAnalysis(de %>% dplyr::filter(avg_log2FC > 0))

  saveRDS(de, paste0("inhibitory/data/de.", res, ".rds"))
  saveRDS(pa, paste0("inhibitory/data/pa.", res, ".rds"))
  
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res))
}
SaveH5Seurat(obj, "inhibitory/data/inhibitory.h5Seurat", overwrite = T, verbose = F)



# Plot top differential expressed genes as well as expression of known signatures
pdf(file.path(obj@misc$graph.path, "signatures.pdf"), width=24, height = 10)

# UMAP of clustering resolutions *after* reordering according to cluster tree
DimPlot(obj, group.by = paste0("SCT_snn_res.", resolutions), label=T, raster=T) & NoAxes() & NoLegend()

source("utils/signatures.R")
for (res in resolutions) {
  de <- readRDS(paste0("inhibitory/data/de.", res, ".rds"))
  
  print(DotPlot(obj, 
                group.by = paste0("SCT_snn_res.", res), 
                features = de %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC) %>% pull(gene) %>% unique) + 
          scale_color_viridis_c() + theme(legend.position = "bottom") + 
          RotatedAxis() + btheme +
          labs(x=NULL, y=NULL, title=paste0("Differential Expression - ", paste0("SCT_snn_res.", res))))
  
  for(publication in c("Cain A. et al., (2021)")) {
    print(DotPlot(obj, 
                  features=signatures$`Inhibitory Neurons`[[publication]],
                  group.by = paste0("SCT_snn_res.", res)) + 
            RotatedAxis() + 
            scale_color_viridis_c() + 
            theme(legend.position = "bottom") + 
            labs(X=NULL, y=NULL, title=paste0(publication, " - Cluster Resolution:", res)))
  }
}
while (!is.null(dev.list()))  dev.off()




####################################################################################################################
##                                    #  Atlas Identities of Sub-populations    #                                 ##
####################################################################################################################
# Required RAM: Load object: 95GB, rest of analysis 425GB (17.5h)
obj <- LoadH5Seurat("inhibitory/data/inhibitory.h5Seurat")

# Create base identities
obj$v1.base <- plyr::mapvalues(obj$SCT_snn_res.0.5, from=c(12,18,9), to=c(11,6,8))
Idents(obj) <- obj$v1.base
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
obj$v1.base <- Idents(obj)

obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "Inhibitory",
                            class     = "Neuronal",
                            cell.type = "Inhibitory Neurons",
                            state     = paste0("Inh.", v1.base),
                            sub.population = case_when(v1.base == 6 & SCT_snn_res.0.5 == 6 ~ "Inh.6.1",
                                                       v1.base == 6 & SCT_snn_res.0.5 != 6 ~ "Inh.6.2",
                                                       v1.base == 10 & SCT_snn_res.0.5 == 11 ~ "Inh.10.1",
                                                       v1.base == 10 & SCT_snn_res.0.5 != 11 ~ "Inh.10.2"),
                            annotation = case_when(v1.base %in% c(9,10,11) ~ "VIP",
                                                   v1.base == 2 ~ "KIT",
                                                   v1.base %in% c(6,7) ~ "SST",
                                                   v1.base %in% c(13,14,15,16) ~ "PVALB",
                                                   v1.base == 5 ~ "NPY")
                     ) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))
Idents(obj) <- obj$state
SaveH5Seurat(obj, "inhibitory/data/inhibitory.h5Seurat", overwrite = T, verbose = F)


pdf(file.path(obj@misc$graph.path,"final.idents.pdf"), width=9, height = 6)
print(DimPlot(obj, label=T, pt.size=.1, raster=T) + NoAxes() + NoLegend() + labs(title=NULL))

for(g in c("class","cell.type","state","sub.population", "annotation")) {
  print(LabelClusters(DimPlot(obj, label=F, pt.size=.1, raster=T, group.by = g), id = g,
                      position = "median", box = T, repel=T, size=5, box.padding=1, min.segment.length=1, fill = rgb(1, 1, 1, alpha = 0.8),colour = '#0C0C0C') + NoAxes() + labs(title=NULL))
}
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------- #
#                     Compute DE Genes                  # 
# ----------------------------------------------------- #
obj <- LoadH5Seurat("inibitory/data/inibitory.h5Seurat")

options(future.globals.maxSize= 14*1024^3)
de <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 4)
saveRDS(de, "inhibitory/data/de.rds")

de.pairwise <- do.call(rbind, combn(unique(Idents(obj)), 2, simplify = F, FUN = function(pair) {
  pair = as.character(pair)
  plan(multicore, workers = 4)
  FindMarkers(obj, ident.1 = pair[[1]], ident.2 = pair[[2]], test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = F) %>%
    mutate(comparison = paste(if(avg_log2FC >= 0) pair else rev(pair), collapse = " vs. "),
           cluster = if_else(avg_log2FC >=0, pair[[1]], pair[[2]]),
           avg_log2FC = abs(avg_log2FC),
           gene = rownames(.),
           id=GeneIdMapping()$ids[gene])
}))
saveRDS(de.pairwise, "inhibitory/datade.pairwise.rds")



obj <- LoadH5Seurat("inibitory/data/inibitory.h5Seurat", assays=list(SCT=c("data")), misc=T, graphs=F, neighbors=F, reductions=F, verbose=F)
saveRDS(EnrichmentAnalysis(readRDS("inhibitory/data/de.rds") %>% 
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"), 
                                           cluster = gsub("\\.","_", cluster)), 
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "pa.rds")

saveRDS(EnrichmentAnalysis(readRDS("de.pairwise.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           comparison = gsub("\\.", "_", gsub(" vs\\. ", "_vs_", comparison))), 
                           formula="id~comparison+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "pa.pairwise.rds")
