setwd("2. Cell-type analysis")

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
suppressPackageStartupMessages(source("utils/analysis.reports.R"))

####################################################################################################################
##                                  #  Create Layer Specific Neurons Objects   #                                  ##
####################################################################################################################
# At the time of analysis, using the available Seurat package, we were not able to create a single 
# excitatory neurons object. As such we randomly split excitatory neuron cells into 3 subsets for 
# the initial QCs. Then, we merged cells into a CUX2+ and CUX2- subsets of excitatory neurons. Sub-
# clustering analysis was performed within each of the two subsets. 

# ----------------------------------------------------------- #
# Create temporary objects of multiple batches                #
# ----------------------------------------------------------- #
# Required RAM: 216GB
source("utils/subset.cell.type.R")
lapply(c("excitatory", "excitatory/graphs", "excitatory/data"), dir.create)

libs <- subset.cell.type("excitatory", remove.doublets=T, return.library.lists=T)
all.lib.names <- names(libs$counts)

cells <- sample(unlist(unname(lapply(names(libs$counts), function(n) colnames(libs$counts[[n]]))))) %>% split(., 1:3)
pb <- progress_bar$new(format=paste0("Splitting excitatory neurons into ", length(cells), " objects :current/:total [:bar] :percent in :elapsed. ETA :eta"),
                       total = length(cells), clear=F, width=100, force = T)
for(i in seq_along(cells)) {
  obj <- lapply(names(libs$counts), function(n) {
    ids <- dplyr::intersect(colnames(libs$counts[[n]]), cells[[i]])
    CreateSeuratObject(counts = libs$counts[[n]][, ids], meta.data = libs$meta.data[[n]][ids, ])
  })
  obj <- merge(obj[[1]], obj[-1])
  
  Project(obj) <- paste0("excitatory.subset.", i)
  obj@misc$mmc <- data.frame(mmc=do.call(rbind, libs$mmc))
  obj@misc$cell.type.predictions <- do.call(rbind, unname(lapply(libs$cell.type.predictions, function(df) df[dplyr::intersect(rownames(df), colnames(obj)),] )))
  obj@misc$doub.parent.distribution <- do.call(rbind, unname(lapply(names(libs$doub.parent.distribution), function(n) {
    df <- libs$doub.parent.distribution[[n]]
    data.frame(batch=n, df[df$cell %in% colnames(obj),])
  }))) 
  
  obj@misc$graph.path <- "excitatory/graphs"

  capture.output(SaveH5Seurat(obj, filename = paste0("excitatory/data/", Project(obj), ".h5Seurat"), verbose=F, overwrite = T), file=NULL)
  obj <- NULL; pb$tick()
}



# ----------------------------------------------------------- #
# Basic analysis and QCs for each subset                      #
# ----------------------------------------------------------- #
# Required RAM: `subset.1` = 163GB, `subset.2` = 165GB, `subset.3` = 162GB
for(f in c(list.files("excitatory/data", "excitatory.subset.*.h5Seurat", full.names = T, recursive = F))) {
  message(f)
  obj <- LoadH5Seurat(f, verbose=F)
  obj <- subset(obj, cells = colnames(obj)[obj$percent.mt <= 10 & obj$nCount_RNA <= 80000])
  obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
  obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
  obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
  obj <- FindClusters(obj, resolution = 1, algorithm = 4, method = "igraph", verbose = F)
  obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
  LowQualityStatus(obj, "SCT_snn_res.1", filename=file.path(obj@misc$graph.path, paste0(Project(obj), ".QCs1.pdf")))
  SaveH5Seurat(obj, f, overwrite = T, verbose = F)
}


# ----------------------------------------------------------- #
# Removing low quality & additional doublets                  #
# ----------------------------------------------------------- #
# Required RAM: `subset.1` = 200GB, `subset.2` = 222GB, `subset.3` = 218GB
to.remove <- list(`excitatory.subset.1`  = c(13,14,25,26,27,28,30,32), 
                  `excitatory.subset.2`  = c(13,23,24,25,26,27,31),
                  `excitatory.subset.3`  = c(14,23,24,25,26,28,31))
for(f in c(list.files("excitatory/data", "excitatory.subset.*.h5Seurat", full.names = T, recursive = F))) {
  message("File: ", f)
  n   <- gsub(".h5Seurat", "", basename(f))
  
  obj <- LoadH5Seurat(f, verbose=F)
  obj <- subset(obj, cells = colnames(obj)[!obj$SCT_snn_res.1 %in% to.remove[[n]]])
  obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
  obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
  obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
  obj <- FindClusters(obj, resolution = 1, algorithm = 4, method = "igraph", verbose = F)
  obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
  LowQualityStatus(obj, "SCT_snn_res.1", filename=file.path(obj@misc$graph.path, paste0(Project(obj), ".QCs2.pdf")))
  SaveH5Seurat(obj, gsub(".h5Seurat", "", f), overwrite = T, verbose = F)
  
  
  markers <- list(
    GABAergic=c("MEG3", "PVALB","SST","VIP", "KIT","GAD2"),
    Glutamatergic=c("SLC17A7", "RORB", "TOX","FOXP2", "CUX2"),
    OPCs=c("PCDH15", "MEGF11", "VCAN", "PDGFRA", "CSPG4"),
    Microglia=c("C3","PTPRC", "TREM2"),
    Endothelial=c("FLT1", "CLDN5","ABCB1","ATP10A"),
    Pericytes=c("PDGFRB","DCN"),
    Oligodendrocytes=c("MBP", "MOG", "MAG"),
    Astrocytes=c("CD44","SLC1A2","SLC1A3","APOE","GJA1","GFAP", "ALDH1L1")
  )
  
  pdf(paste0("excitatory/graphs/", Project(obj), ".cell.type.markers.pdf"), width=20, height = 20)
  print(FeaturePlot(obj, features=unlist(markers), order=T, raster = T) & NoAxes() & scale_color_viridis_c())

  print(plot_grid(plot_grid(DimPlot(obj, label=T, raster=T) + NoAxes() + NoLegend(),
                            DotPlot(obj, features=unlist(markers$Glutamatergic)) + btheme + scale_color_viridis_c(),
                            ncol=1),
                  FeaturePlot(obj, features=unlist(markers$Glutamatergic), order=T, raster = T) & NoAxes() & scale_color_viridis_c()))
  while (!is.null(dev.list()))  dev.off()
}



####################################################################################################################
##                                #  Sub-clustering analysis for CUX2+ cells   #                                  ##
####################################################################################################################
# In the case of CUX2+ cells, we again were not able to create a single Seurat object due Seurat's use of 32-bit 
# matrices (at time of analysis). We therefore created 2 CUX2+ subsets, performed QCs and filtering for each, and
# only then merged into a single CUX2+ object. Final sub-clustering analysis was done over the joint CUX2+ object.

# ------------------------------------------------ #
#            Create CUX2+ Seurat object            #
# ------------------------------------------------ #
# Required RAM: 140GB
subsets <- paste0("excitatory/data/excitatory.subset.", 1:3, ".h5Seurat")
clusters = list(`excitatory.subset.1`  = c(1,2,3,5,6,7,10,11,21,22,23,25,26,27), 
                `excitatory.subset.2`  = c(1,2,4,5,6,10,11,12,14,19,23,24,26,27), # Excluding #12 
                `excitatory.subset.3`  = c(4,5,6,7,8,9,11,12,13,14,21,22,24,25,26)) # Excluding #9

obj = list(); pred.df = list(); doub.par = list()
for(f in subsets) {
  n <- gsub(".h5Seurat", "", basename(f))
  o <- LoadH5Seurat(f, assays=list(RNA=c("counts")), misc=T, verbose = F)
  
  cell.ids       <- colnames(o)[(o$SCT_snn_res.1 %in% clusters[[n]]) & (o$nCount_RNA < 100000)]
  obj[[n]]       <- subset(o, cells = cell.ids)
  pred.df[[n]]   <- o@misc$cell.type.predictions[cell.ids, ]
  doub.par[[n]]  <- o@misc$doub.parent.distribution[o@misc$doub.parent.distribution$cell %in% cell.ids, ]
  
  rm(o); gc()
}
pred.df <- do.call(rbind, unname(pred.df))
doub.par <- do.call(rbind, doub.par)


message("Total number of CUX2+ cells is:", sum(sapply(obj, ncol)), ". Splitting into 2 objects.")
cells <- split(sample(unlist(unname(lapply(obj, colnames)))), 1:2)
name <- "cux2+"
for(i in seq_along(cells)) {
  message("Creating Seurat object of size: ", length(cells[[i]]))
  o <- lapply(names(obj), function(n) subset(obj[[n]], cells = dplyr::intersect(colnames(obj[[n]]), cells[[i]])))
  o <- merge(o[[1]], o[-1])
  print(o)
  
  Project(o) <- paste0("cux2+.", i)
  o@misc$cell.type.predictions <- pred.df[colnames(o),]
  o@misc$doub.parent.distribution <- doub.par[doub.par$cell %in% colnames(o),]
  
  SaveH5Seurat(o, paste0("excitatory/data/excitatory.cux2+.", i, ".h5Seurat"))
  rm(o); gc()
}



# ------------------------------------------------ #
#         Search For Low Quality Clusters          #
# ------------------------------------------------ #
# Required RAM: `cux2+.1`=137GB, `cux2+.2`=137GB
for (path in paste0("excitatory/data/excitatory.cux2+.", 1:2, ".h5Seurat")) {
  message(path)
  obj <- LoadH5Seurat(path, verbose = F)
  obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
  obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
  obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
  obj@meta.data <- obj@meta.data[,!grepl("SCT_snn_res", colnames(obj@meta.data))]
  obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method = "igraph", verbose = F)
  obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model=T,  verbose = F)
  LowQualityStatus(obj, "SCT_snn_res.1.2", filename=paste0("excitatory/graphs/", gsub(".h5Seurat", "", basename(path)), ".QCs1.pdf"))
  SaveH5Seurat(obj, path, overwrite = T, verbose = T)
}


# Required RAM: `cux2+.1`=178GB, `cux2+.2`=171GB
clusters = list(`cux2+.1` = c(7), `cux2+.2`=c(12,13))
for (path in paste0("excitatory/data/excitatory.cux2+.", 1:2, ".h5Seurat")) {
  obj <- LoadH5Seurat(path, verbose = F)
  obj <- subset(obj, cells = colnames(obj)[!obj$SCT_snn_res.1.2 %in% clusters[[gsub(".h5Seurat", "", basename(path))]]])
  obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
  obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
  obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
  obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method = "igraph", verbose = F)
  obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model=T,  verbose = F)
  LowQualityStatus(obj, "SCT_snn_res.1.2", filename=paste0("excitatory/graphs/", gsub(".h5Seurat", "", basename(path)), ".QCs1.pdf"))
  SaveH5Seurat(obj, path, overwrite = T, verbose = T)
}


# ------------------------------------------------ #
#       Create Signle Seurat CUX2+ Object          #
# ------------------------------------------------ #
# Required RAM: 125GB
obj <- list(); pred.df <- list(); doub.par <- list()
for (path in paste0("excitatory/data/excitatory.cux2+.", 1:2, ".h5Seurat")) {
  o <- LoadH5Seurat(path, assays = list(RNA=c("counts","data")), misc=T, reductions=F, graphs=F, verbose = T)
  o@misc$sim.umap <- NULL
  
  pred.df[[path]]  <- o@misc$cell.type.predictions
  doub.par[[path]] <- o@misc$doub.parent.distribution
  obj[[path]]      <- o
  rm(o); gc()
}
obj <- merge(obj[[1]], obj[[2]])

Project(obj) <- "cux2+"
obj@misc$cell.type.predictions    <- do.call(rbind, unname(pred.df))
obj@misc$doub.parent.distribution <- do.call(rbind, doub.par)
obj@misc$graph.path     <- "excitatory/graphs"

obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj@meta.data <- obj@meta.data[,!grepl("SCT_snn_res", colnames(obj@meta.data))]
obj <- FindClusters(obj, resolution = c(.1, .2, .3, 1.2), algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, return.model=T,  verbose = F)

LowQualityStatus(obj, "SCT_snn_res.1.2", filename="excitatory/graphs/cux2+.QCs1.pdf")
PlotClusteringResolutions(obj, path="excitatory/graphs/cux2+.clustering.resolutions.pdf")
SaveH5Seurat(obj, "excitatory/data/cux2+.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#               Determine clustering               #
# ------------------------------------------------ #

# To avoid clusters splitting by sex (with PCs driven only by sex-chromosome genes) we rerun pipeline using only the first 10 PCs
obj <- LoadH5Seurat("excitatory/data/cux2+.h5Seurat")
obj <- FindNeighbors(obj, dims=1:10, verbose=F)
obj <- FindClusters(obj, resolution = c(.1, .2, .3), algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:10, n.components = 2, min.dist = .1, return.model=T,  verbose = F)
PlotClusteringResolutions(obj, path="excitatory/graphs/cux2+.clustering.resolutions.10pcs.pdf")
SaveH5Seurat(obj, "excitatory/data/cux2+.h5Seurat", overwrite = T, verbose = F)


# Required RAM: `0.1` = 248GB (13h), `0.2` = 277GB (13h), `0.3` = 297GB
obj <- LoadH5Seurat("excitatory/data/cux2+.h5Seurat")
for (res in c("0.1", "0.2", "0.3")) {
  # Building cluster tree
  Idents(obj) <- obj@meta.data[, paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)
  
  # Differential expression (DE) + pathway analysis (PA)
  de <- FindMarkersWrapper(obj, test.use="negbinom")
  pa <- EnrichmentAnalysis(de %>% dplyr::filter(avg_log2FC > 0))

  saveRDS(de, paste0("excitatory/data/cux2+.de.", res, ".rds"))
  saveRDS(pa, paste0("excitatory/data/cux2+.pa.", res, ".rds"))

  # Preliminary trait associations functionality
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res))
}
SaveH5Seurat(obj, "excitatory/data/cux2+.h5Seurat", overwrite = T, verbose = F)


# Plot top differential expressed genes 
pdf(file.path(obj@misc$graph.path,".signatures.pdf"), width=24, height = 6)
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



####################################################################################################################
##                                #  Sub-clustering analysis for CUX2- cells   #                                  ##
####################################################################################################################

# ------------------------------------------------ #
#            Create CUX2- Seurat object            #
# ------------------------------------------------ #
# Required RAM: 100GB
library(progress)
subsets <- list.files("excitatory/data", "excitatory.subset.*.h5Seurat", full.names = T, recursive = F)
obj = list(); pred.df = list(); doub.par = list()
clusters = list(`excitatory.subset.1`  = c(4,8,9,12,13,14,15,16,17,18,19,20,24), 
                `excitatory.subset.2`  = c(3,7,8,9,13,15,16,17,18,20,21,22,25), 
                `excitatory.subset.3`  = c(1,2,3,10,15,16,17,18,19,20,23,27))

# ----------------------------------------------------------- #
# Extracting desired clusters from each Ex. neurons subset    #
# ----------------------------------------------------------- #
pb <- progress_bar$new(format="Splitting CUX2- excitatory neuron subsets :current/:total [:bar] :percent in :elapsed. ETA :eta",
                       total = length(subsets), clear=F, width=100, force = T)
for(f in subsets) {
  n <- gsub(".h5Seurat", "", basename(f))
  o <- LoadH5Seurat(f, assays=list(RNA=c("counts")), misc=T, verbose = F)
  
  cell.ids      <- colnames(o)[o$SCT_snn_res.1 %in% clusters[[n]]]
  obj[[n]]      <- CreateSeuratObject(o@assays$RNA@counts[, cell.ids], meta.data = o@meta.data[cell.ids,])
  pred.df[[n]]  <- o@misc$cell.type.predictions[cell.ids, ]
  doub.par[[n]] <- o@misc$doub.parent.distribution[o@misc$doub.parent.distribution$cell %in% cell.ids, ]
  
  rm(o); gc()
  pb$tick()
}
rm(pb, f, n, cell.ids, subsets, clusters)


message(paste0("Mering cells to create CUX2- branch: ", do.call(sum, lapply(obj, ncol)), " cells"))
obj <- merge(obj[[1]], obj[-1]); gc()
Project(obj) <- "cux2-"
obj@misc$cell.type.predictions <- do.call(rbind, unname(pred.df))
obj@misc$doub.parent.distribution <- do.call(rbind, unname(doub.par))

obj@misc$graph.path <- file.path("excitatory/graphs")
SaveH5Seurat(obj, "excitatory/data/cux2-.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#         Search For Low Quality Clusters          #
# ------------------------------------------------ #
# Required RAM: 188GB
obj <- LoadH5Seurat("excitatory/data/cux2-.h5Seurat")
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, dims = 1:50, verbose = F)
obj@meta.data <- obj@meta.data[,!grepl("SCT_snn_res", colnames(obj@meta.data))]
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj, dims=1:50, n.components = 2, min.dist = .1, return.model=T,  verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename="excitatory/graphs/cux2-.QCs1.pdf")
SaveH5Seurat(obj, "excitatory/data/cux2-.h5Seurat", overwrite = T, verbose = F)



# Remove low-quality/doublet clusters
# Required RAM: 282GB
obj <- LoadH5Seurat("excitatory/data/cux2-.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[!obj$SCT_snn_res.1.2 %in% c(33,32,31,30,29,28,27,24)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj, dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename="excitatory/graphs/cux2-.QCs2.pdf")
SaveH5Seurat(obj, "excitatory/data/cux2-.h5Seurat", overwrite = T, verbose = F)



# Required RAM: 312GB
obj <- LoadH5Seurat("excitatory/data/cux2-.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[!obj$SCT_snn_res.1.2 %in% c(19, 24)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj, dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename="excitatory/graphs/cux2-.QCs3.pdf")
SaveH5Seurat(obj, "excitatory/data/cux2-.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#               Determine clustering               #
# ------------------------------------------------ #
library(future)
obj <- LoadH5Seurat("excitatory/data/cux2-.h5Seurat")
plan(multiprocess, workers = 2)
obj <- FindClusters(obj, resolution = c(.1, .3, .7), algorithm = 4, method = "igraph", verbose = F)
PlotClusteringResolutions(obj, path="excitatory/graphs/cux2-.clustering.resolutions.pdf")
SaveH5Seurat(obj, "excitatory/data/cux2-.h5Seurat", overwrite = T, verbose = F)

# Required RAM: 235GB
obj <- LoadH5Seurat("excitatory/data/cux2-.h5Seurat")
for (res in c("0.1", "0.3", "0.7")) {
  # Building cluster tree
  Idents(obj) <- obj@meta.data[, paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)
  
  # Differential expression (DE) + pathway analysis (PA)
  de <- FindMarkersWrapper(obj, test.use="negbinom")
  pa <- EnrichmentAnalysis(de %>% dplyr::filter(avg_log2FC > 0))
  
  saveRDS(de, paste0("excitatory/data/cux2-.de.", res, ".rds"))
  saveRDS(pa, paste0("excitatory/data/cux2-.pa.", res, ".rds"))
  
  # Preliminary trait associations functionality
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res))
}
SaveH5Seurat(obj, "excitatory/data/cux2-.h5Seurat", overwrite = T, verbose = F)


# Plot top differential expressed genes as well as expression of known signatures
pdf(file.path(obj@misc$graph.path, "signatures.pdf"), width=24, height = 10)

# UMAP of clustering resolutions *after* reordering according to cluster tree
DimPlot(obj, group.by = paste0("SCT_snn_res.", resolutions), label=T, raster=T) & NoAxes() & NoLegend()

source("utils/signatures.R")
for (res in resolutions) {
  de <- readRDS(paste0("excitatory/data/cux2-.de.", res, ".rds"))
  
  print(DotPlot(obj, 
                group.by = paste0("SCT_snn_res.", res), 
                features = de %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC) %>% pull(gene) %>% unique) + 
          scale_color_viridis_c() + theme(legend.position = "bottom") + 
          RotatedAxis() + btheme +
          labs(x=NULL, y=NULL, title=paste0("Differential Expression - ", paste0("SCT_snn_res.", res))))
  
  for(publication in c("Cain A. et al., (2021)")) {
    print(DotPlot(obj, 
                  features=signatures$`Excitatory Neurons`[[publication]],
                  group.by = paste0("SCT_snn_res.", res)) + 
            RotatedAxis() + 
            scale_color_viridis_c() + 
            theme(legend.position = "bottom") + 
            labs(X=NULL, y=NULL, title=paste0(publication, " - Cluster Resolution:", res)))
  }
}
while (!is.null(dev.list()))  dev.off()




####################################################################################################################
##                         #  Refining annotations of "bridge" cells between CUX2+ and CUX2- #                    ##
####################################################################################################################
# While sub-clustering CUX2+ and CUX2- cells we notice that for some "bridge" cells the separation into +/- is not
# accurate enough. To refine the overall clustering of excitatory neurons we pull these "bridge" cells from the two
# subsets, and perform a separate clustering analysis.

# Required RAM: 140GB
o1 <- LoadH5Seurat("excitatory/data/cux2-.h5Seurat", assays = list(RNA=c("counts")), misc=F, graphs=F, reductions=F, neighbors=F) %>%
  subset(., SCT_snn_res.0.1 %in% c(2,4))
o1$original.clusters <- o1$SCT_snn_res.0.1

o2 <- LoadH5Seurat("excitatory/data/cux2+.h5Seurat", assays = list(RNA=c("counts")), misc=F, graphs=F, reductions=F, neighbors=F) %>%
  subset(, SCT_snn_res.0.3 %in% c(1,2,3))
o2$original.clusters <- o2$SCT_snn_res.0.3

obj <- merge(o1, o2)
rm(o1,o2); gc()

# Seurat pipeline analysis
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, dims = 1:50, verbose = F)
obj@meta.data <- obj@meta.data[,!grepl("SCT_snn_res", colnames(obj@meta.data))]
obj <- FindClusters(obj, resolution = c(.05, .1, .2, .3), algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj, dims=1:50, n.components = 2, min.dist = .1, return.model=T,  verbose = F)
PlotClusteringResolutions(obj, path="excitatory/graphs/bridge.clustering.resolutions.pdf")
SaveH5Seurat(obj, "excitatory/data/bridge.h5Seurat", overwrite = T, verbose = F)



####################################################################################################################
##                                    #  Atlas Identities of Sub-populations    #                                 ##
####################################################################################################################
bridge <- LoadH5Seurat("excitatory/data/bridge.h5Seurat", assays = list(SCT=c("data")), reductions=F, neighbors=F, graphs=F, misc=F, verbose=F)

# ------------------------------------------------ #
#                     CUX2+ Subset                 #
# ------------------------------------------------ #
obj <- LoadH5Seurat("excitatory/data/cux2+.h5Seurat")
obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "CUX2+",
                            class     = "Neuronal",
                            cell.type = "Excitatory Neurons",
                            v1.base   = as.character(SCT_snn_res.0.3),
                            state     = case_when(v1.base == "5" ~ "Exc.1",
                                                  v1.base == "7" ~ "Exc.2",
                                                  v1.base == "6" ~ "Exc.3",
                                                  v1.base == "8" ~ "Exc.4",
                                                  v1.base == "4" ~ "Exc.5",
                                                  rownames(.) %in% colnames(bridge)[bridge$SCT_snn_res.0.1 == 2] ~ "Exc.6",
                                                  rownames(.) %in% colnames(bridge)[bridge$SCT_snn_res.0.1 == 3] ~ "Exc.7",
                                                  rownames(.) %in% colnames(bridge)[bridge$SCT_snn_res.0.1 == 1] ~ "Exc.8",
                                                  rownames(.) %in% colnames(bridge)[bridge$SCT_snn_res.0.1 == 4] ~ "Exc.11"),
                            sub.population = NA_character_,
                            annotation = NA_character_) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))

Idents(obj) <- obj$state
SaveH5Seurat(obj, "excitatory/data/cux2+.h5Seurat", overwrite = T, verbose = F)



# ------------------------------------------------ #
#                     CUX2- Subset                 #
# ------------------------------------------------ #
obj <- LoadH5Seurat("excitatory/data/cux2-.h5Seurat")
obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "CUX2-",
                            class     = "Neuronal",
                            cell.type = "Excitatory Neurons",
                            v1.base   = as.character(SCT_snn_res.0.1),
                            state     = case_when(rownames(.) %in% colnames(bridge)[bridge$SCT_snn_res.0.1 == 1] ~ "Exc.8",
                                                  rownames(.) %in% colnames(bridge)[bridge$SCT_snn_res.0.1 == 2] ~ "Exc.6",
                                                  rownames(.) %in% colnames(bridge)[bridge$SCT_snn_res.0.1 == 3] ~ "Exc.7",
                                                  rownames(.) %in% colnames(bridge)[bridge$SCT_snn_res.0.1 == 4] ~ "Exc.11",
                                                  v1.base == "1" ~ "Exc.9",
                                                  v1.base == "3" ~ "Exc.10",
                                                  v1.base == "5" ~ "Exc.12",
                                                  v1.base == "6" ~ "Exc.13",
                                                  v1.base == "7" ~ "Exc.14",
                                                  v1.base == "8" ~ "Exc.15",
                                                  v1.base == "9" ~ "Exc.16"),
                            sub.population = NA_character_,
                            annotation = NA_character_) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))
Idents(obj) <- obj$state
SaveH5Seurat(obj, "excitatory/data/cux2-.h5Seurat", overwrite = T, verbose = F)




# ----------------------------------------------------------- #
# DEGs and PA over all excitatory neurons                     #
# ----------------------------------------------------------- #
obj <- list(`cux2-`=LoadProject("excitatory", "cux2-", assays=list(SCT=c("counts")), graphs=F, misc=F, neighbors=F, reductions=F, verbose=F),
            `cux2+`=LoadProject("excitatory", "cux2+", assays=list(SCT=c("counts")), graphs=F, misc=F, neighbors=F, reductions=F, verbose=F))

for (ident in unique(obj[[1]]$state)) {
  pb <- progress_bar$new(format=paste("DEGs for", ident, " :current/:total [:bar] :percent in :elapsed. ETA :eta"),
                         total = length(genes), clear=F, width=100, force = T)
  
  genes <- split(unique(unlist(lapply(obj, rownames))), 1:6)
  de <- do.call(rbind, lapply(genes, function(gene.subset) {
    o <- merge(subset(obj[[1]], features = gene.subset),
               subset(obj[[2]], features = gene.subset))
    Idents(o) <- o$state
    
    df <- FindMarkers(o, ident.1 = ident, test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = F)
    pb$tick()
    return(df)
  }))
  saveRDS(de, paste0("excitatory/data/de.",ident, ".rds"))
}
rm(obj, ident, pb, genes, de)

# Merge results
de <- lapply(list.files("excitatory/data", "Exc.\\d+.rds", full.names = T), function(f) 
  readRDS(f) %>% 
    mutate(cluster = gsub(".rds", "", basename(f)),
           gene = gsub("\\d+\\.","", rownames(.)),
           id = GeneIdMapping()$ids[gene]) %>%
    `rownames<-`(NULL)) %>%
  do.call(rbind, .) %>% 
  mutate(p_val_adj = p.adjust(p_val, method = "BH"))
saveRDS(de, "excitatory/data/de.rds")


# Run pathways analysis
universe <- lapply(c("cux2+","cux2-"), function(b)
  LoadProject("excitatory", b, assays=list(SCT=c("data")), 
              graphs=F, misc=F, neighbors=F, reductions=F, verbose=F) %>% rownames) %>%
  unlist %>% unique()

saveRDS(EnrichmentAnalysis(readRDS("excitatory/data/de.rds") %>% 
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"), 
                                           cluster = gsub("\\.","_", cluster)), 
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[universe]),
        "excitatory/data/pa.rds")
rm(universe)

