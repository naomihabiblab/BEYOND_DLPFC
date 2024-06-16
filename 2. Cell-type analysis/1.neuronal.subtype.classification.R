library(SingleR)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(dplyr)
library(tibble)


# The following script trains a SingleR classifier over the Allen Brain institute atlas
# to annotate the neuronal subtypes in the present atlas. 

ref.labels.full <- read.csv("2. Cell-type analysis/data/allen.annotations.csv") %>% column_to_rownames("sample_name")
ref.counts.full <- read.csv("2. Cell-type analysis/data/allen.counts.matrix.csv", sep = ",", row.names = 1)

neurons <- list(
  inhibitory = list(label.atlas = "GABAergic",     seurat.obj = "2. Cell-type analysis/inhibitory/data/inhibitory.h5Seurat"),
  `cux2+`    = list(label.atlas = "Glutamatergic", seurat.obj = "2. Cell-type analysis/excitatory/data/cux2+.h5Seurat"),
  `cux2-`    = list(label.atlas = "Glutamatergic", seurat.obj = "2. Cell-type analysis/excitatory/data/cux2-.h5Seurat"))

for(conf in neurons) {
  # Load query neuronal dataset
  query.counts <- GetAssayData(LoadH5Seurat(conf$seurat.obj, assays = c(SCT = "data"), graphs = FALSE, images = NULL)[["SCT"]]) %>% as.matrix()
  
  # Define shared gene space
  gene.space <- intersect(colnames(ref.counts), rownames(query.counts))
  
  # Subset reference:
  #   1. Samples -  to cells of current neuronal type
  #   2. Features - to shared gene space
  cells      <- intersect(ref.labels.full %>% filter(class_label == conf$label.atlas) %>% rownames, rownames(ref.counts))
  ref.counts <- t(ref.counts[cells, gene.space])
  ref.labels <- ref.labels.full %>% `[`(cells, "cluster_label")
  rm(cells)
  
  # Predict neuronal sub-type annotations
  preds <- SingleR(test = query.counts[gene.space,], ref = ref.counts, labels = ref.labels)
  saveRDS(preds, gsub(".h5Seurat", ".predicted.annotation.rds", conf$seurat.obj))
}

