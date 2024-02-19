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

# Required RAM: 100GB
source("utils/subset.cell.type.R")
obj <- subset.cell.type("oligodendroglia")
obj@misc$graph.path <- file.path(name, "graphs")
lapply(c("oligodendroglia", "oligodendroglia/graphs", "oligodendroglia/data"), dir.create)


# Basic analysis to visualize object doublets 
# Required RAM: 165GB
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50)
obj <- FindClusters(obj, resolution = .5, algorithm = 4, method="igraph")
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1)

InitialDoubletsStatus(obj)
SaveH5Seurat(obj, "oligodendroglia/data/oligodendroglia.h5Seurat")


# ----------------------------------------------------------- #
# Remove doubltes, keeping oligo-opc-astr doublets            #
# ----------------------------------------------------------- #
source("utils/doublet.calling.plots.R")

# Required RAM: 155GB
obj <- LoadH5Seurat("oligodendroglia/data/oligodendroglia.h5Seurat")
obj <- FindClusters(obj, resolution = 2, algorithm = 4, method="igraph")
obj$prior.doub.removal.SCT_snn_res.2 <- obj$SCT_snn_res.2
PlotDoubletDistributions(obj, "SCT_snn_res.2", remove.doublet.group=c("Olig-OPC"))

# Removing non (oli-opc-ast) DoubletFinder doublets (demux-doublets - removed)
# Removed clusters: 
#   c56 - oli-neuronal(-ast)
#   c49 endo/peri 
#   c26,48 mic
#   c39 - suspected as opc/ast - will be revisited later
#   c43,46 - suspected as opc/exc+opc/inh doublets - will be revisited later - DoubletFinder doublets are removed
#   c53 - suspected as opc-oli trajectory - not removed even if with high DoubletFinder scores
obj <- subset(obj, cells = colnames(obj)[!(obj$demux.droplet.type == "DBL" |
                                          (obj$SCT_snn_res.2 %in% c(26,48,49,56)) |
                                          (obj$is.doublet & ! obj$SCT_snn_res.2 %in% c(39,53) ))] )
SaveH5Seurat(obj, "oligodendroglia/data/oligodendroglia.h5Seurat", overwrite = T, verbose = F)


 
####################################################################################################################
##                                     #  Search For Low Quality Clusters   #                                     ##
####################################################################################################################

# Required RAM: 150GB
obj <- LoadH5Seurat("oligodendroglia/data/oligodendroglia.h5Seurat")
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj, dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename=file.path(obj@misc$graph.path,"QCs1.pdf"))
SaveH5Seurat(obj, "oligodendroglia/data/oligodendroglia.h5Seurat", overwrite = T, verbose = F)




# Required RAM: 185GB
# c30 - some mic/endothelial genes - to be revisited later
# c34 - suspected as OPC-Ast doubltes - to be revisited later
# c35 - Some inh neurons genes but low doublet/classification scores - to be revisited later
obj <- LoadH5Seurat("oligodendroglia/data/oligodendroglia.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[obj$percent.mt <= 10])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method = "igraph", verbose = F)
obj <- RunUMAP(obj, dims=1:50, n.components = 2, min.dist = .1, return.model = T, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename=file.path(obj@misc$graph.path,"QCs2.pdf"))
SaveH5Seurat(obj, "oligodendroglia/data/oligodendroglia.h5Seurat", overwrite = T, verbose = F)



####################################################################################################################
##                                  #  Split Object To Oligodendrocytes & OPCs  #                                 ##
####################################################################################################################
# While there are small doublet-suspected clusters they will be addressed as part of the oligodendrocytes or OPCs subsets
obj <- LoadH5Seurat("oligodendroglia/data/oligodendroglia.h5Seurat")
obj <- FindClusters(obj, resolution = .3, algorithm = 4, method = "igraph", verbose = F)
SaveH5Seurat(obj, "oligodendroglia/data/oligodendroglia.h5Seurat", overwrite = T, verbose = F)


subs <- list(opcs = c(2,18,19,20), oligodendrocytes = setdiff(unique(obj$SCT_snn_res.0.3), c(2,18,19,20)))
for (i in seq_along(subs)) {
  sub <- subset(obj, SCT_snn_res.0.3 %in% subs[[i]])
  
  Project(sub) <- names(subs)[[i]]
  sub@misc$graph.path <- file.path(obj@misc$graph.path, Project(sub))
  lapply(file.path(c("oligodendroglia/graphs", "oligodendroglia/data"), Project(sub)), dir.create)

  SaveH5Seurat(sub, paste0("oligodendroglia/data/",Project(sub),"/", Project(sub), ".h5Seurat"))
}
rm(obj, subs, sub, i)




####################################################################################################################
##                                       #  Sub-clustering analysis OPCs   #                                      ##
####################################################################################################################
suppressPackageStartupMessages(source("src/analysis.reports.R"))

# ------------------------------------------------ #
#         Search For Low Quality Clusters          #
# ------------------------------------------------ #

# Required RAM: 33GB
obj <- LoadH5Seurat("oligodendroglia/data/opcs/opcs.h5Seurat")
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method="igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename=file.path(obj@misc$graph.path, "QCs1.pdf"))
SaveH5Seurat(obj, "oligodendroglia/data/opcs/opcs.h5Seurat", overwrite = T, verbose = F)



# Required RAM: 34GB
# Removed clusters:
#   c11 - OPC-Oli high doublets scores - check differentiation markers
#   c13 - OPC-Exc suspected doublets
#   c17 - OPC-Ast suspected doublets
#   c22 - OPC-Inh suspected doublets
#   c26 - high for Oli genes - revisit in next iteration
obj <- LoadH5Seurat("oligodendroglia/data/opcs/opcs.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[!obj$SCT_snn_res.1.2 %in% c(13, 17, 22)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method="igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename=file.path(obj@misc$graph.path, "QCs2.pdf"))
SaveH5Seurat(obj, "oligodendroglia/data/opcs/opcs.h5Seurat", overwrite = T, verbose = F)



obj <- LoadH5Seurat("oligodendroglia/data/opcs/opcs.h5Seurat")
obj@meta.data <- obj@meta.data[,!grepl("SCT_snn_res", colnames(obj@meta.data))]
obj <- FindClusters(obj, resolution = c(.05,.1, .2, .3, .5), algorithm = 4, method = "igraph", verbose = F)
PlotClusteringResolutions(obj)
SaveH5Seurat(obj, "oligodendroglia/data/opcs/opcs.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#               Determine clustering               #
# ------------------------------------------------ #
# Ran each resolution independently
# Required RAM: 30GB
obj <- LoadH5Seurat("oligodendroglia/data/opcs/opcs.h5Seurat")
for (res in c("0.05", "0.1", "0.2", "0.3", "0.5")) {
  # Building cluster tree
  Idents(obj) <- obj@meta.data[, paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)
  
  # Differential expression (DE) + pathway analysis (PA)
  de <- FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"))
  pa <- EnrichmentAnalysis(de %>% dplyr::filter(avg_log2FC > 0))
  
  saveRDS(de, paste0("oligodendroglia/data/opcs/de.", res, ".rds"))
  saveRDS(pa, paste0("oligodendroglia/data/opcs/pa.", res, ".rds"))
  
  # Preliminary trait associations functionality
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res))
}
SaveH5Seurat(obj, "oligodendroglia/data/opcs/opcs.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#        Atlas Identities of Sub-Populations       #
# ------------------------------------------------ #
obj <- LoadH5Seurat("oligodendroglia/data/opcs/opcs.h5Seurat")

# Create base identities
Idents(obj) <- obj$v1.base <- obj$SCT_snn_res.0.1
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
obj$v1.base <- Idents(obj)

obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "OPCs",
                            class     = "Glia",
                            cell.type = recode(v1.base, "1" = "COP", "2" = "MFOL", .default = "OPCs"),
                            state     = case_when(v1.base == 1 ~ "COP",
                                                  v1.base == 2 ~ "MFOL",
                                                  v1.base == 3 ~ "OPC.1",
                                                  v1.base %in% c(5,6,7,8) ~ "OPC.2",
                                                  v1.base == 4 ~ "OPC.3"),
                            sub.population = state) %>%
                     mutate(annotation = state) %>%
                     select(subset, class, cell.type, state, sub.population, annotation, v1.base))
Idents(obj) <- obj$state
SaveH5Seurat(obj, "oligodendroglia/data/opcs/opcs.h5Seurat", overwrite = T, verbose = F)

pdf(file.path(obj@misc$graph.path,"final.idents.pdf"), width=9, height = 6)
for(g in c("v1.base","class","cell.type","state", "annotation")) {
  print(LabelClusters(DimPlot(obj, label=F, pt.size=.1, raster=T, group.by = g), id = g,
                      position = "median", box = T, repel=T, size=5, box.padding=1, min.segment.length=1, fill = rgb(1, 1, 1, alpha = 0.8),colour = '#0C0C0C') + NoAxes() + labs(title=g))
}
while (!is.null(dev.list()))  dev.off()


# ----------------------------------------------------- #
#                     Compute DE Genes                  # 
# ----------------------------------------------------- #
# Required RAM: 85GB
obj <- LoadH5Seurat("oligodendroglia/data/opcs/opcs.h5Seurat")

options(future.globals.maxSize= 5*1024^3)
saveRDS(FindMarkersWrapper(obj, test.use="negbinom", latent.vars = c("batch","pmi"), workers = 4), 
        "oligodendroglia/data/opcs/de.rds")

de.pairwise <- do.call(rbind, combn(unique(Idents(obj)), 2, simplify = F, FUN = function(pair) {
  pair = as.character(pair)
  plan(multicore, workers = 4)
  FindMarkers(obj, ident.1 = pair[[1]], ident.2 = pair[[2]], test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = F) %>%
    mutate(cluster = if_else(avg_log2FC >=0, pair[[1]], pair[[2]]),
           avg_log2FC = abs(avg_log2FC),
           gene = rownames(.),
           comparison = paste(pair, collapse = " vs. "),
           id=GeneIdMapping()$ids[gene])
}))
saveRDS(de.pairwise, "oligodendroglia/data/opcs/de.pairwise.rds")



obj <- LoadH5Seurat("oligodendroglia/data/opcs/opcs.h5Seurat", assays=list(SCT=c("data")), misc=T, graphs=F, neighbors=F, reductions=F,verbose=T)

saveRDS(EnrichmentAnalysis(readRDS("oligodendroglia/data/opcs/de.rds") %>% 
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"), 
                                           cluster = gsub("\\.","_", cluster)), 
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "oligodendroglia/data/opcs/pa.rds")

saveRDS(EnrichmentAnalysis(readRDS("oligodendroglia/data/opcs/de.pairwise.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           comparison = gsub("\\.", "_", gsub(" vs\\. ", "_vs_", comparison))), 
                           formula="id~comparison+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "oligodendroglia/data/opcs/pa.pairwise.rds")





####################################################################################################################
##                                 #  Sub-clustering analysis Oligodendrocytes  #                                 ##
####################################################################################################################
suppressPackageStartupMessages(source("src/analysis.reports.R"))

# ------------------------------------------------ #
#         Search For Low Quality Clusters          #
# ------------------------------------------------ #

# Required RAM: 120GB
obj <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat")
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method="igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename=file.path(obj@misc$graph.path, "QCs1.pdf"))
SaveH5Seurat(obj, "oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat", overwrite = T, verbose = F)



# Required RAM: 147GB
# Removed clusters:
# c21 - Oli-Exc suspected doublets
# c29 - Suspected doublets with Ast/Mic/Peri/End - probably a generally low-qual/noisy cluster
obj <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[!obj$SCT_snn_res.1.2 %in% c(21, 29)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = 1.2, algorithm = 4, method="igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)
LowQualityStatus(obj, "SCT_snn_res.1.2", filename=file.path(obj@misc$graph.path, "QCs2.pdf"))
SaveH5Seurat(obj, "oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat", overwrite = T, verbose = F)


# ------------------------------------------------ #
#               Determine clustering               #
# ------------------------------------------------ #

# Required RAM: 103GB
obj <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat")
obj@meta.data <- obj@meta.data[,!grepl("SCT_snn_res", colnames(obj@meta.data))]
obj <- FindClusters(obj, resolution = c(.1, .15, .2, .25, .3, .5), algorithm = 4, method = "igraph", verbose = F)
PlotClusteringResolutions(obj)
SaveH5Seurat(obj, "oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat", overwrite = T, verbose = F)


obj <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat")
for (res in c(.1, .2, .3, .5)) {
  Idents(obj) <- obj@meta.data[, paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)
  
  # Preliminary trait associations functionality
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res))
}
SaveH5Seurat(obj, "oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat", overwrite = T, verbose = F)


# -------------------------------------------------- #
# Sub-clustering of core cells for better separation #
# -------------------------------------------------- #
obj <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat")
obj <- subset(obj, cells = colnames(obj)[!obj$SCT_snn_res.0.2 %in% c(1,2,3,4,6,10,16)])
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:50, verbose = F)
obj <- FindClusters(obj, resolution = c(.1, .15, .2, .25 ,.3), algorithm = 4, method="igraph", verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)

for (res in c(.1, .15, .2, .25 ,.3)) {
  Idents(obj) <- obj@meta.data[, paste0("SCT_snn_res.", res)]
  obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
  obj@meta.data[, paste0("SCT_snn_res.", res)] <- Idents(obj)
  
  # Preliminary trait associations functionality
  NaitveTraitsAssociation(obj, paste0("SCT_snn_res.", res), filename = file.path(obj@misc$graph.path, paste0("core.trait.association.", res, ".pdf")))
}
PlotClusteringResolutions(obj, path=file.path(obj@misc$graph.path, "core.clustering.resolutions.pdf"))
SaveH5Seurat(obj, "oligodendroglia/data/oligodendrocytes/oligodendrocytes.core.h5Seurat")


# ------------------------------------------------ #
#        Atlas Identities of Sub-Populations       #
# ------------------------------------------------ #
obj  <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat")
core <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.core.h5Seurat",
                     assays=list(SCT=c("data")), reductions=F, neighbors=F, misc=F)

df <- data.frame(main = plyr::mapvalues(as.character(obj$SCT_snn_res.0.2), from=c(9), to=c(16)), 
                 core = NA,
                 row.names = colnames(obj))

# c16 was decided to be merged with 16 instead of being part of the "core"
core.ids <- setdiff(colnames(core), colnames(obj)[obj$SCT_snn_res.0.2 == 9])
# Merging merging clusters by cluster tree
df[core.ids,]$core <- 
  plyr::mapvalues(as.character(core@meta.data[core.ids,]$SCT_snn_res.0.25),
                  from = c(1, 2, 6,8, 9, 10,11,13,14,15,16,17,18),
                  to   = c(12,12,7,12,12,12,12,12,19,19,19,19,19)) 
rm(core)

obj <- AddMetaData(obj, df %>% mutate(v1.base = if_else(is.na(core), paste0("all.", main), paste0("core.", core))) %>% 
                     dplyr::select(v1.base))
obj <- subset(obj, cells = colnames(obj)[obj$v1.base != "all.1"])
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, min.dist = .1, verbose = F)

# Create base identities
Idents(obj) <- obj$v1.base
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T, dims = 1:50, verbose = F)
saveRDS(Tool(obj, "BuildClusterTree"), file.path(dirname(obj@misc$data.path), paste0("tree.v1.idents.rds")))
obj$v1.base <- Idents(obj)

obj <- AddMetaData(obj, 
                   obj@meta.data %>% 
                     mutate(subset    = "Oligodendrocytes",
                            class     = "Glia",
                            cell.type = "Oligodendrocytes",
                            state     = paste0("Oli.", plyr::mapvalues(as.character(v1.base), 
                                                                       from=c(1, 3, 5, 7,2,10,11,12,4,6,8,9),
                                                                       to=c(12,11,10,9,1,2, 3, 4, 5,6,7,8))),
                            sub.population = state,
                            annotation = state) %>%
                     select(subset, class, cell.type, state, v1.base, sub.population, annotation))
Idents(obj) <- obj$state
SaveH5Seurat(obj, "oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat", overwrite = T, verbose = F)


pdf(file.path(obj@misc$graph.path,"final.idents.pdf"), width=9, height = 6)
for(g in c("v1.base","class","cell.type","state")) {
  print(LabelClusters(DimPlot(obj, label=F, pt.size=.5, raster=T, group.by = g), id = g,
                      position = "median", box = T, repel=T, size=5, box.padding=1, min.segment.length=1, fill = rgb(1, 1, 1, alpha = 0.8),colour = '#0C0C0C') + NoAxes() + labs(title=g))
}
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------- #
#                     Compute DE Genes                  # 
# ----------------------------------------------------- #
# Due to long running times, computation was split into batches 
# of idents, each running as an independent job
obj  <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat")

# ------------------------ one-vs-all ----------------- #
options(future.globals.maxSize= 6*1024^3)
x <- split(unique(Idents(obj)), 1:4)
for(i in 1:4) {
  de <- lapply(x[[i]], function(v) {
    FindMarkers(obj, ident.1 = v, test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = F) %>%
      mutate(cluster = v,
             gene = rownames(.),
             id=GeneIdMapping()$ids[gene])
  }) %>% do.call(rbind, .)
  saveRDS(de, paste0("oligodendroglia/data/oligodendrocytes/de.final.", i,".rds"))
}
rm(x,i)

# Merge results 
lapply(paste0("oligodendroglia/data/oligodendrocytes/de.final.", 1:5,".rds"), readRDS) %>%
  do.call(rbind, .) %>%
  saveRDS(., "oligodendroglia/data/oligodendrocytes/de.rds")



# ------------------------ all-vs-all ----------------- #
options(future.globals.maxSize= 6*1024^3)
x <- split(combn(unique(Idents(obj)), 2, simplify = F, FUN = function(pair) pair), 1:10)
for(i in 1:10) {
  de.pairwise <- do.call(rbind, lapply(x[[i]], function(pair) {
    pair = as.character(pair)
    FindMarkers(obj, ident.1 = pair[[1]], ident.2 = pair[[2]], test.use = "negbinom", latent.vars = c("batch","pmi"), verbose = F) %>%
      mutate(gene = rownames(.), id=GeneIdMapping()$ids[gene])
  }))
  saveRDS(de.pairwise, paste0("oligodendroglia/data/oligodendrocytes/de.pairwise.", i, ".rds"))
}
rm(x,i)

# Merge results 
de <- lapply(paste0("oligodendroglia/data/oligodendrocytes/de.pairwise.", 1:10, ".rds"), readRDS) %>%
  do.call(rbind, .) %>%
  tidyr::separate(comparison, c("a", "b"), " vs\\. ") %>% 
  mutate(comparison = if_else(avg_log2FC >= 0, paste(a,b, sep=" vs. "), paste(b,a, sep=" vs. ")), 
         cluster = if_else(avg_log2FC >= 0, a, b), 
         avg_log2FC = abs(avg_log2FC)) %>% 
  dplyr::select(p_val, avg_log2FC, pct.1, pct.2, p_val_adj, comparison, cluster, gene, id)
saveRDS(de, "oligodendroglia/data/oligodendrocytes/de.pairwise.rds")


# ----------------------------------------------------- #
#                     Compute Pathways                  # 
# ----------------------------------------------------- #
obj  <- LoadH5Seurat("oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat", 
                     assays=list(SCT=c("data")), misc=T, graphs=F, neighbors=F, reductions=F,verbose=F)

saveRDS(EnrichmentAnalysis(readRDS("oligodendroglia/data/oligodendrocytes/de.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           cluster = gsub("\\.","_", cluster)),
                           formula="id~cluster+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "oligodendroglia/data/oligodendrocytes/pa.rds")

saveRDS(EnrichmentAnalysis(readRDS("oligodendroglia/data/oligodendrocytes/de.pairwise.rds") %>%
                             dplyr::mutate(direction = if_else(avg_log2FC > 0, "upregulated","downregulated"),
                                           comparison = gsub("\\.", "_", gsub(" vs\\. ", "_vs_", comparison))), 
                           formula="id~comparison+direction",
                           universe = GeneIdMapping()$ids[rownames(obj)]),
        "oligodendroglia/data/oligodendrocytes/pa.pairwise.rds")



####################################################################################################################
##                                      #  Merge data from Oli & OPC subsets  #                                   ##
####################################################################################################################

cells <- lapply(c("oligodendroglia/data/opcs/opcs.h5Seurat", "oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat"), function(o)
  LoadH5Seurat(o, assays=list(SCT=c("data")))@meta.data[,c("v1.base","subset","class","state","sub.population","annotation")]) %>%
  do.call(rbind, .)

obj <- LoadH5Seurat("oligodendroglia/data/oligodendroglia.h5Seurat")
obj <- subset(obj, cells = rownames(cells))
obj <- AddMetaData(obj, cells)
Idents(obj) <- obj$state

# Rerunning UMAP embedding after oli/opc cells' cleanup, and with different parameters to promote a larger spread 
# within each of the oli or opc, but not between
obj <- RunUMAP(obj,  reduction = "pca", dims=1:50, n.components = 2, spread=2, min.dist=.1, verbose = F)
SaveH5Seurat(obj, "oligodendroglia/data/oligodendroglia.h5Seurat", overwrite = T, verbose = F)
