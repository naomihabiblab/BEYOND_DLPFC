source("2. Cell-type analysis/load.code.env.R")


####################################################################################################################
##                                #  Creating cellular landscape representation #                                 ##
####################################################################################################################
sc <- reticulate::import("scanpy")
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

# Participant clustering based on cellular environment representation
sc$pp$neighbors(data, n_neighbors = as.integer(10), use_rep = "X", metric = "cosine")
sc$tl$leiden(data, resolution =.25)
sc$tl$leiden(data, resolution = .75, restrict_to=reticulate::tuple("leiden", reticulate::np_array(c("0"))))
data$obs["clusters"] = plyr::mapvalues(data$obs$leiden_R, levels(data$obs$leiden_R), 1:length(levels(data$obs$leiden_R)))
data$obs$leiden <- data$obs$leiden_R <- NULL

data$obs$core <- !data$obs$clusters %in% c(9,10)

# 2D visualization of landscape
sc$tl$tsne(data, n_pcs = 0, use_rep = "X", learning_rate = 100)
sc$tl$umap(data, maxiter = as.integer(3000), spread = 3)

# Compute PHATE embedding for all donors
sc$external$tl$phate(data, 
                     n_components = as.integer(3),  
                     k = as.integer(10), a = as.integer(40), 
                     knn_dist =  "euclidean", mds_dist = "euclidean", 
                     mds_solver = "smacof", verbose = F)
data$obsm$X_all_3d_phate <- data$obsm$X_phate
data$obsm$X_phate <- NULL

# Compute PHATE embedding for core donors
sub <- data[data$obs$core]
sc$external$tl$phate(sub, n_components = as.integer(2),  
                     k = as.integer(15), a = as.integer(100), 
                     knn_dist =  "euclidean", mds_dist = "correlation", 
                     mds_solver = "smacof", verbose = F)

data$obsm$X_core_phate <- 
  data.frame(sub$obsm$X_phate, row.names = sub$obs_names) %>%
  merge(data.frame(row.names = data$obs_names), 
        by.x="row.names", by.y="row.names", all.y=T) %>% 
  tibble::column_to_rownames("Row.names") %>%
  `[`(data$obs_names,)
rm(sub)

# Local similarities of participants based on cellular environments
data$obsp <- list()
for(e in c("X_all_3d_phate","X_core_phate","X_umap","X_tsne")) 
  data$obsp[[paste0("similarity_", e)]] <- embedding.similatity(data$obsm[[e]], knn = 5)
anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(e, sc)
