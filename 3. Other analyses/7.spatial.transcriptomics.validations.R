

####################################################################################################################
##                                       #  Initial QCs and Quantification  #                                     ##
####################################################################################################################
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
samples <- c("R7594705", "R3631183", "R1777884", "R5955028", "R7876628",
             "R4077358", "R6665276", "R8882846", "R3111222", "R2006886")

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
  mutate(trajectory = case_when(sample %in% c("R3111222","R8882846","R4077358","R1777884") ~ "Early",
                                sample %in% c("R2006886","R6665276") ~ "ABA",
                                sample %in% c("R3631183","R5955028","R7594705","R7876628") ~ "prAD")) %>%
  dplyr::select(sample, trajectory, Mic.13, Ast.10, Ast.5)

saveRDS(df, "3. Other analyses/data/ST.validation.state.signatures.rds")


####################################################################################################################
##                                        #  ST Co-localization Analysis  #                                       ##
####################################################################################################################

results <- list()
df <- readRDS("3. Other analyses/data/ST.validation.state.signatures.rds") %>%
  arrange(trajectory, sample) %>% mutate(sample=paste(trajectory, sample))


# -------------------------------------------------------- #
# Within participant Mic.13-Ast.10 colocalization          #
# -------------------------------------------------------- #
vals <- df[,c("sample","trajectory")] %>% unique %>% `rownames<-`(NULL)
vals <- lapply(1:nrow(vals), function(i){ 
  t <- cor.test(df[df$sample == vals[i,"sample"], "Ast.10"], df[df$sample == vals[i,"sample"], "Mic.13"], alternative = "greater")
  vals[i,] %>% mutate(cor=t$estimate, pval=t$p.value)
}) %>% do.call(rbind,.) %>%
  mutate(adj.pval = p.adjust(pval, method="BH"),
         sig = cut(adj.pval, c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", "")),
         is.sig = adj.pval < .01) %>% 
  arrange(trajectory)

results$mic13.ast10.colocalization <- vals
rm(vals)


# -------------------------------------------------------- #
# Summarize Mic.13-Ast.10 colocalization                   #
# -------------------------------------------------------- #
corrs <- results$mic13.ast10.colocalization %>% split(., .$trajectory) %>% lapply(., \(v) v$cor)
results$mic13.ast10.colocalization.summarize <- 
  t.test(corrs$prAD, c(corrs$ABA, corrs$Early), alternative = "greater")


# -------------------------------------------------------- #
# Within participant Ast.10-Ast.5 colocalization           #
# -------------------------------------------------------- #

# Obtain mean Ast.10 and Ast.5 for all participants
vals <- df %>% 
  group_by(trajectory, sample) %>% 
  summarise(across(c(Ast.10, Ast.5), mean), .groups = "drop") %>%
  mutate(`10.vs.5` = log(Ast.10/Ast.5))

# For prAD or ABA participants test colocalization
vals <- vals %>% 
  merge(., 
        rbind(
          # For prAD participants - is mean Ast.10 > Ast.5
          lapply(vals[vals$trajectory == "prAD", ]$sample, \(s) {
            t <- t.test(df[df$sample == s, "Ast.10"], df[df$sample == s, "Ast.5"], "greater")
            data.frame(sample = s, test = "Is Ast.10 > Ast.5", t.stat = t$statistic, pval=t$p.value)
          }) %>% do.call(rbind, .),
          
          # For ABA participants - is mean Ast.55 > Ast.10
          lapply(vals[vals$trajectory == "ABA", ]$sample, \(s) {
            t <- t.test(df[df$sample == s, "Ast.5"], df[df$sample == s, "Ast.10"], "greater")
            data.frame(sample = s, test = "Is Ast.5 > Ast.10", t.stat = t$statistic, pval=t$p.value)
          }) %>% do.call(rbind, .)),
        all.x=TRUE) %>%
  mutate(adj.pval = p.adjust(pval, method="BH"),
         sig = cut(adj.pval, c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", "")),
         is.sig = adj.pval < .01) %>%
  arrange(trajectory, sample) %>%
  column_to_rownames("sample")

results$ast10.ast5.colocalization <- vals
rm(vals)


# -------------------------------------------------------- #
# Summarize Ast.10-Ast.5 colocalization                    #
# -------------------------------------------------------- #

`10.vs.5` <- results$ast10.ast5.colocalization %>% split(., .$trajectory) %>% lapply(., \(v) v$`10.vs.5`)
results$ast.10.ast5.colocalization.summarize <- 
  t.test(`10.vs.5`$prAD, `10.vs.5`$ABA, alternative = "greater")
rm(`10.vs.5`)

saveRDS(results, "3. Other analyses/data/ST.validation.rds")