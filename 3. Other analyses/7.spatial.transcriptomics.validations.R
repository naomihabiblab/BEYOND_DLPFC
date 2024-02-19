#' Natacha Comandante-Lou (nc3018 at columbia dot edu)

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(logger)
library(tidyverse)
library(scCustomize)
library(Rmagic)
library(viridis)

source("3. Other analyses/utils/spatial.utils.R")

markers <- list(Ast.10=c("SMTN","SLC38A2"), Ast.5=c("OSMR","SERPINA3"), Mic.13=c("TREM2","GPNMB"))
samples <- c("20260284_2", "35941263_2", "50107871_2","50301675_2", "69982533_1", 
             "50500136_3", "20242958_2", "20865035_2", "11327005_1", "20254902_2")

df <- lapply(samples, function(s) {
  # Load Seurat spatial object and 
  meta.data = read.csv(file.path("3. Other analyses/data/", sample, "spot.annotation.csv")) %>% tibble::column_to_rownames("SpotID")
  obj <- preprocess_visium_h5(file.path("Other analyses/data/", samples, "outs"), meta.data=meta.data)
  obj <- subset(obj, features = rownames(obj)[Matrix::rowSums(obj@assays$RNA@counts) > 10])
  
  # Add spatial coordinate as reduction
  obj[["coord"]] <- CreateDimReducObject(
    embeddings = obj@images[["slice1"]]@coordinates %>% select(Spatial_1 = imagerow, Spatial_2 = imagecol) %>% as.matrix, 
    key = "Spatial_", assay = "Spatial")
  
  # Run MAGIC imputation
  # TODO: library.size.normalize(data) %>% sqrt - still normalize and sqrt data
  magic.res <- magic(t(GetAssayData(obj, slot = "count",assay = "Spatial")), 
                     genes = colnames(data), 
                     init = NULL, 
                     knn = 5, 
                     t = 3)
  obj[["MAGIC"]] <- CreateAssayObject(data = magic.res[["result"]] %>% as.matrix %>% t)
  Misc(obj, slot = "magic") = marig.res$operator
  rm(magic.res)
  
  # remove WM or unknown layer
  obj <- subset(obj, !is.na(SpatialCluster) & SpatialCluster != "WM")
  
  # Obtain gene signature expression, scaling each gene expression between 0-1
  DefaultAssay(spt) <- "MAGIC"
  FetchData(obj, vars = unlist(markers), slot = "data") %>% 
    mutate(across(unlist(markers), scales::rescale)) %>%
    mutate(sample = gsub("_.*", "", s), 
           Ast.10 = rowSums(select(., matches(markers$Ast.10))),
           Ast.5 = rowSums(select(., matches(markers$Ast.5))),
           Mic.13 = rowSums(select(., matches(markers$Mic.13))))
}) %>% do.call(rbind, .) %>%
  mutate(trajectory = case_when(sample %in% c("20254902","50301675","11327005","20260284") ~ "Early",
                                sample %in% c("20865035","50500136") ~ "ABA",
                                sample %in% c("50107871","69982533","35941263","20242958") ~ "prAD")) %>%
  dplyr::select(sample, trajectory, Mic.13, Ast.10, Ast.5)

saveRDS(df, "3. Other analyses/data/ST.validation.state.signatures.rds")