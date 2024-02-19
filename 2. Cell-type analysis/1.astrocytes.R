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

# Required RAM: 74GB
source("utils/subset.cell.type.R")
obj <- subset.cell.type("astrocytes")
obj@misc$graph.path <- file.path(name, "graphs")
lapply(c("astrocytes", "astrocytes/graphs", "astrocytes/data"), dir.create)


# Basic analysis to visualize object doublets 
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 30, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:30, verbose = F)
obj <- FindClusters(obj, resolution = .2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:30, verbose = F)
InitialDoubletsStatus(obj)
SaveH5Seurat(obj, "astrocytes/data/astrocytes.h5Seurat")



####################################################################################################################
##                                     #  Search For Low Quality Clusters   #                                     ##
####################################################################################################################
obj <- LoadH5Seurat("astrocytes/data/astrocytes.h5Seurat")
obj <- subset(obj, is.doublet == F)
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs1.pdf"))
SaveH5Seurat(obj, "astrocytes/data/astrocytes.h5Seurat", overwrite = T, verbose = F)


# Required RAM: 149GB
# clusters: #45 - Endo-Astr; #39 - Micr-Astr; #35 - OPC-Inh-Ast; #29 - Olig-Astr; #26 - Exc-Astr;
obj <- LoadH5Seurat("astrocytes/data/astrocytes.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[obj$percent.mt <=10 & obj$nFeature_RNA <= 6000 & !obj$SCT_snn_res.1.5 %in% c(26,29,35,39,45)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.5, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj, dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, ident="SCT_snn_res.1.5", filename=file.path(obj@misc$graph.path, "QCs2.pdf"))
SaveH5Seurat(obj, "astrocytes/data/astrocytes.h5Seurat", overwrite = T, verbose = F)


# Required RAM: 74GB
obj <- LoadH5Seurat("astrocytes/data/astrocytes.h5Seurat")
obj$SCT_snn_res.1.5 <- NULL
obj <- FindClusters(obj, resolution = c(.1, .3, .5, .7, 1), algorithm = 4, method = "igraph", verbose = F)
PlotClusteringResolutions(obj)
SaveH5Seurat(obj, "astrocytes/data/astrocytes.h5Seurat", overwrite = T, verbose = F)


####################################################################################################################
##                                            #  Determine Clustering   #                                         ##
####################################################################################################################
# Investigated resolutions
resolutions <- c("0.1", "0.3"  ,"0.5")

# Executed in a separate job for each resolution, with 2 cores per job
# Required RAM: `0.1`=236GB, `0.3`=237GB, `0.5`=238GB
obj <- LoadH5Seurat("astrocytes/data/astrocytes.h5Seurat")
for (res in resolutions) {
  # Building cluster tree
  Idents(obj) <- obj@meta.data[, paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)

  # Differential expression (DE) + pathway analysis (PA)
  options(future.globals.maxSize= 10*1024^3)
  de <- FindMarkersWrapper(obj, test.use="negbinom", workers = 2)
  pa <- EnrichmentAnalysis(de %>% dplyr::filter(avg_log2FC > 0))

  saveRDS(de, paste0("astrocytes/data/de.", res, ".rds"))
  saveRDS(pa, paste0("astrocytes/data/pa.", res, ".rds"))

  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res))
}
SaveH5Seurat(obj, "astrocytes/data/astrocytes.h5Seurat", overwrite = T, verbose = F)


# Plot top differential expressed genes as well as expression of known signatures
pdf(file.path(obj@misc$graph.path, "signatures.pdf"), width=24, height = 10)

# UMAP of clustering resolutions *after* reordering according to cluster tree
DimPlot(obj, group.by = paste0("SCT_snn_res.", resolutions), label=T, raster=T) & NoAxes() & NoLegend()

source("utils/signatures.R")
for (res in resolutions) {
  de <- readRDS(paste0("astrocytes/data/de.", res, ".rds"))
  
  print(DotPlot(obj, 
                group.by = paste0("SCT_snn_res.", res), 
                features = de %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC) %>% pull(gene) %>% unique) + 
          scale_color_viridis_c() + theme(legend.position = "bottom") + 
          RotatedAxis() + btheme +
          labs(x=NULL, y=NULL, title=paste0("Differential Expression - ", paste0("SCT_snn_res.", res))))
  
  for(publication in c("Cain A. et al., (2021)")) {
    print(DotPlot(obj, 
                  features=signatures$Astrocytes[[publication]],
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
# Required RAM: 237GB
# Required RAM for loading object (60GB)
obj <- LoadH5Seurat("astrocytes/data/astrocytes.h5Seurat")

# Create base identities
obj$v1.base <- as.character(plyr::mapvalues(obj$SCT_snn_res.0.3, from=c(5,8,9,10,2,3,6,14), to=c(12,11,11,11,1,1,7,13)))

Idents(obj) <- obj$v1.base
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
saveRDS(Tool(obj, "BuildClusterTree"), file.path(dirname(obj@misc$data.path), paste0("tree.v1.idents.rds")))
obj$v1.base <- Idents(obj)

obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "Astrocytes",
                            class     = "Glia",
                            cell.type = "Astrocyte",
                            state     = paste0("Ast.", v1.base),
                            sub.population = NA_character_,
                            annotation = NA_character_) %>%
                     mutate(state = recode(state,
                                           "Ast.1"="Ast.4",
                                           "Ast.2"="Ast.7",
                                           "Ast.4"="Ast.8",
                                           "Ast.7"="Ast.2",
                                           "Ast.8"="Ast.1",
                                           .default = state)) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))
Idents(obj) <- obj$state
SaveH5Seurat(obj, "astrocytes/data/astrocytes.h5Seurat", overwrite = T, verbose = F)


pdf(file.path(obj@misc$graph.path,"final.idents.pdf"), width=9, height = 6)
for(g in c("v1.base","class","cell.type","state", "annotation")) {
  print(LabelClusters(DimPlot(obj, label=F, pt.size=.1, raster=T, group.by = g), id = g,
                      position = "median", box = T, repel=T, size=5, box.padding=1, min.segment.length=1, fill = rgb(1, 1, 1, alpha = 0.8),colour = '#0C0C0C') + NoAxes() + labs(title=g))
}
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------- #
#                     Compute DE Genes                  # 
# ----------------------------------------------------- #
# Required RAM: 230GB
obj <- LoadH5Seurat("astrocytes/data/astrocytes.h5Seurat")

options(future.globals.maxSize= 10*1024^3)
obj@misc$de.v1.idents <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 2)
saveRDS(obj@misc$de.v1.idents, "astrocytes/data/de.rds")

obj@misc$de.v1.idents.pairwise <- do.call(rbind, combn(unique(Idents(obj)), 2, simplify = F, FUN = function(pair) {
  pair = as.character(pair)
  plan(multicore, workers = 2)
  FindMarkers(obj, ident.1 = pair[[1]], ident.2 = pair[[2]], test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = F) %>%
    mutate(comparison = paste(if(avg_log2FC >= 0) pair else rev(pair), collapse = " vs. "),
           cluster = if_else(avg_log2FC >=0, pair[[1]], pair[[2]]),
           avg_log2FC = abs(avg_log2FC),
           gene = rownames(.),
           id=GeneIdMapping()$ids[gene])
}))
saveRDS(obj@misc$de.v1.idents.pairwise, "astrocytes/data/de.pairwise.rds")


# Required RAM: 38BG
obj <- LoadH5Seurat("astrocytes/data/astrocytes.h5Seurat", assays=list(SCT=c("data")), misc=T, graphs=F, neighbors=F, reductions=F,verbose=F)

saveRDS(EnrichmentAnalysis(readRDS("astrocytes/data/de.rds") %>% 
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"), 
                                           cluster = gsub("\\.","_", cluster)), 
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "pa.rds")


saveRDS(EnrichmentAnalysis(readRDS("astrocytes/data/de.rds") %>%
                             tidyr::separate("comparison", c("a","b"), " vs. ",remove = F) %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           comparison = if_else(cluster == a,paste0(a, " vs. ", b), paste0(b, " vs. ", a))) %>%
                             dplyr::mutate(comparison = gsub("\\.", "_", gsub(" vs\\. ", "_vs_", comparison))), 
                           formula="id~comparison+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "pa.pairwise.rds")

