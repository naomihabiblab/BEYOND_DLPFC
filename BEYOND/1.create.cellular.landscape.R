source("BEYOND/utils.R")


####################################################################################################################
##                                #  Creating landscape representation object  #                                  ##
####################################################################################################################
df  <- h5read(file = "Cell-type analysis/DLPFC.Green.atlas.h5", name = "annotations") %>% filter(state != "NA")
qcs <- read.csv("BEYOND/data/donors.qc.csv") %>% filter(Final_QC == "Pass")
qcs <- qcs[qcs$projid %in% unique(df$projid),] %>% column_to_rownames("projid")
df  <- df[df$projid %in% rownames(qcs), ]


donor.batches <- df %>% group_by(projid) %>% summarise(batch = paste(unique(batch), collapse=", ")) %>% column_to_rownames("projid")
main.batch <- df %>% mutate(batch = gsub("-[A|B]$", "", batch)) %>% 
  count(projid, batch) %>% 
  group_by(projid) %>% 
  slice_max(order_by=n, n=1) %>% 
  column_to_rownames("projid") %>% 
  dplyr::select(batch)


df <- df %>% 
  dplyr::select(-one_of(c("cell","batch", "sub.population", "annotation"))) %>%
  group_by_all() %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(grouping.by, projid) %>%
  mutate(prevalence=n/sum(n)) %>%
  ungroup()
gc()

ids <- unique(df$projid)
excluded.states.df <- df[df$grouping.by == "Immune",]
df <- df[!df$state %in% unique(excluded.states.df$state),]

# transform to donor~state data frame with either counts or prevalence and in row order as specified by ids
tables <- sapply(list(df=df, excluded.df=excluded.states.df), function(df_) 
  sapply(c("n","prevalence"), function(f) {
    . <- dcast(df_, projid~state, value.var = f, fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid")
    
    (merge(data.frame(ids, row.names = ids), ., all.x = T, by = "row.names") %>% 
      dplyr::select(-ids) %>% 
      tibble::column_to_rownames("Row.names") %>% 
      replace(is.na(.), 0))[ids,]
  }, 
  simplify = F, USE.NAMES = T),
  simplify = F, USE.NAMES = T)


data <- AnnData(X = tables$df$prevalence,
                layers = list(counts = tables$df$n,
                              sqrt.prev = sqrt(tables$df$prevalence)),
                var = (df %>% dplyr::select(class, grouping.by, cell.type, state) %>% unique() %>% column_to_rownames("state"))[colnames(tables$df$prevalence),],
                obs = list(batches = donor.batches[ids,],
                           main.batch = main.batch[ids,]),
                obsm = list(QCs = qcs[ids,] %>% `rownames<-`(rownames(.) %>% as.character())),
                uns = list(excluded.states = list(counts = tables$excluded.df$n,
                                                  prev = tables$excluded.df$prevalence,
                                                  sqrt.prev = sqrt(tables$excluded.df$prevalence)))
                )
rm(qcs, ids, excluded.states.df, donor.batches, main.batch)


# -------------------------------------------------------- #
# Append donor metadata                                    #
# -------------------------------------------------------- #
meta.data <- load.metadata()
data$obsm$meta.data <- meta.data[data$obs_names,]
rm(meta.data)


# -------------------------------------------------------- #
# Cell type level of aggregation                           #
# -------------------------------------------------------- #
cell.type.df <- df %>% group_by(projid, class, cell.type) %>%
  dplyr::summarise(n = sum(n), .groups = "drop") %>%
  dplyr::group_by(projid) %>%
  dplyr::mutate(prevalence = n/sum(n)) %>%
  dplyr::group_by(projid, class) %>%
  dplyr::mutate(within.class.prevalence = n/sum(n))

ct.counts <- (dcast(cell.type.df, projid~cell.type, value.var = "n", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid"))[data$obs_names, ]
ct.prev   <- (dcast(cell.type.df, projid~cell.type, value.var = "prevalence", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid"))[data$obs_names, ]
ct.c.prev <- (dcast(cell.type.df, projid~cell.type, value.var = "within.class.prevalence", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid"))[data$obs_names, ]

data$uns$cell.types <- list(counts = ct.counts, 
                            prev = ct.prev, 
                            sqrt.prev = sqrt(ct.prev), 
                            wc.prev = ct.c.prev, 
                            sqrt.wc.prev = sqrt(ct.c.prev))
rm(cell.type.df, ct.counts, ct.prev, ct.c.prev)

anndata::write_h5ad(data, "BEYOND/data/BEYOND.DLPFC.h5ad")
rm(df, data, tables)


# -------------------------------------------------------- #
# Basic scanpy pipeline over donors in landscape           #
# -------------------------------------------------------- #
sc <- reticulate::import("scanpy")

data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")
sc$pp$neighbors(data, n_neighbors = as.integer(10), use_rep = "X", metric = "cosine")
sc$tl$leiden(data, resolution =.25)
sc$tl$leiden(data, resolution = .75, restrict_to=reticulate::tuple("leiden", reticulate::np_array(c("0"))))
data$obs["clusters"] = plyr::mapvalues(data$obs$leiden_R, levels(data$obs$leiden_R), 1:length(levels(data$obs$leiden_R)))
data$obs$leiden <- data$obs$leiden_R <- NULL

data$obs$core <- !data$obs$clusters %in% c(9,10)


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

data$obsp <- list()
for(e in c("X_all_3d_phate","X_core_phate","X_umap","X_tsne")) 
  data$obsp[[paste0("similarity_", e)]] <- embedding.similatity(data$obsm[[e]], knn = 5)
anndata::write_h5ad(data, "BEYOND/data/BEYOND.DLPFC.h5ad")
rm(e, sc)
