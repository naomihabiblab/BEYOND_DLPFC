suppressPackageStartupMessages(source("2. Cell-type analysis/utils/analysis.reports.R"))

library(Seurat)
library(SeuratDisk)
library(progress)
library(dplyr)

#####################################################################################################################
#                                        Create Unified UMAP Space For All Libraries                                #
#####################################################################################################################

# -------------------------------------------------------- #
#         Create Seurat object with 517556 cells           #
# -------------------------------------------------------- #
# # Execution time: ~4:45h
# # Required RAM: 235 GB
libs <- sample(list.files("1. Library preprocessing/data/snRNA-seq libraries", "*.seurat.rds", full.names = T, recursive = T))[1:40]
cells  <- h5read(aggregated.data, "annotations")$cell
pb <- progress_bar$new(format="Loading libraries :current/:total [:bar] :percent in :elapsed. ETA :eta",
                       total = length(libs), clear=F, width=100, force = T)
obj <- lapply(libs, function(p) {
  o <- readRDS(p)
  t <- CreateSeuratObject(counts = o@assays$RNA@counts[rownames(o)[!grepl("^(AC\\d+{3}|AL\\d+{3}|AP\\d+{3}|LINC\\d+{3})", rownames(o))],
                                                       intersect(cells, colnames(o))])
  rm(o); pb$tick(); return(t)
})

message("Merge libraries into single Seurat object")
obj <- merge(obj[[1]], obj[-1])

message("Begin Seurat pipeline analysis")
obj <- SCTransform(obj, variable.features.n = 4000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 50, verbose = F)
obj <- FindNeighbors(obj, dims = 1:50, verbose = F)
obj <- RunUMAP(obj, dims=1:50, n.components = 2, n.neighbors = 200, return.model = T, verbose = F)
SaveH5Seurat(obj, "2. Cell-type analysis/data/unified.umap.h5Seurat", overwrite = T, verbose = F)



# ------------------------------------------------------- #
#      Predict UMAP Embedding on Grouped Libraries        #
# ------------------------------------------------------- #
reference <- LoadH5Seurat("2. Cell-type analysis/data/unified.umap.h5Seurat", assays=list(SCT=c("data")), verbose=F)
message("Loaded reference UMAP object and model. Containing ", ncol(reference), " nuclei")
ref.pca.features <- rownames(reference[["pca"]]@feature.loadings)

libs <- sample(list.files("1. Library preprocessing/data/snRNA-seq libraries", "*.seurat.rds", full.names = T, recursive = T))
pb   <- progress_bar$new(format="Projecting libraries onto reference UMAP space :current/:total [:bar] :percent in :elapsed. ETA :eta",
                         total = ceiling(length(libs) / 5), clear=F, width=100, force = T)
embedding <- lapply(seq(1, length(libs), by=5), function(i) {
  o <- lapply(libs[i:min(i+4,length(libs))], function(p) {
    .o <- readRDS(p)
    CreateSeuratObject(counts = .o@assays$RNA@counts[rownames(.o)[!grepl("^(AC\\d+{3}|AL\\d+{3}|AP\\d+{3}|LINC\\d+{3})", rownames(.o))],])
  })
  o <- merge(o[[1]], o[-1])
  
  # Project onto reference PCA space
  features <- intersect(ref.pca.features, rownames(o))
  o <- SCTransform(o, residual.features = features,  verbose = F)
  o.pca.embeddings <- t(GetAssayData(o[["SCT"]], slot = "scale.data")[features,]) %*% reference[["pca"]]@feature.loadings[features,]
  o[["pca"]] <- CreateDimReducObject(embeddings = o.pca.embeddings, key = "PC_", assay = "SCT"  )
  
  # Project onto reference UMAP space
  o <- ProjectUMAP(query = o, query.reduction = "pca",
                   reference = reference, reference.reduction = "pca",
                   reduction.model = "umap", verbose=F)
  e <- Embeddings(o, "ref.umap")
  rm(o); gc(); pb$tick(); return(e)
})
saveRDS(do.call(rbind, embedding), "2. Cell-type analysis/data/unified.umap.coordinates.rds")
rm(embedding)
