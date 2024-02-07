source("BEYOND/utils.R")


####################################################################################################################
##                        #  Fitting Celmod for prevalence predictions in bulk RNA-seq #                          ##
####################################################################################################################
if(!"Celmod" %in% installed.packages())
  devtools::install_github("MenonLab/Celmod")

# Load snRNAseq and bulk data
data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")
bulk <- read.table("BEYOND/data/ROSMAP.bulkRNAseq.adjusted", as.is=T, header=T, row.names = 1) %>% 
  `colnames<-`(as.character(as.numeric(gsub("X","",colnames(.)))))

# Keep only donors for which we have both bulk and snuc data
shared.donors <- intersect(colnames(bulk), data$obs_names)
sub.data <- data[shared.donors,]

# Perform many train-test splits to average Celmod predictions over
. <- lapply(1:100, function(i) {
  # Split trait- and test sets
  set.seed(i)
  idx <- sample(sub.data$obs_names, size = floor(.75 * nrow(sub.data)))
  
  # Fit CelMod models over each cell type
  start.celmod <- Sys.time()
  res <- lapply(as.character(unique(sub.data$var$grouping.by)), function(g) {
    message("\n",Sys.time(), "\tStart fitting CelMod for ", g)
    
    states <- sub.data$var_names[sub.data$var$grouping.by == g]
    
    # Fitting CelMod model - passing non-sqrt prevalences and setting usesqrt=T
    fit <- Celmod::train_model(bulk[,idx], t(sub.data[idx,]$X[,states]), numgenevec = 10:40, crossval_times = 5, usesqrt=T, method_type = "spearman")
    rm(.Random.seed) # Override Celmod's set of a seed
    
    # Predicting state sqrt prevalences, storing non-sqrt prevalences
    prd <- Celmod::predict_estimates(fit, bulk, usesqrt=T)$proportions[[fit$cv_bestgenes]] %>% `colnames<-`(colnames(bulk)) %>% t
    
    # Extract fitted model parameters and test performance
    fit.genes <- do.call(rbind, lapply(1:ncol(fit$modelgenerank), function(s) {
      used.genes <- which(fit$modelgenerank[,s] <= fit$cv_bestgenes)
      t(fit$model[[s]][,used.genes]) %>% `colnames<-`(c("intercept","coef")) %>% data.frame %>%
        rownames_to_column("gene") %>% mutate(state = states[[s]])
    }))
    
    return(list(genes=fit.genes, preds=prd, model=fit))
  })
  message("CelMod fitting took ", round(difftime(Sys.time(), start.celmod, units = "mins"), 2), " mins")
  
  # Extract results
  fitted.models <- lapply(res, "[[", "model")
  celmod <- list(genes = do.call(rbind, lapply(res, "[[", "genes")),
                 preds = do.call(cbind, lapply(res, "[[", "preds")) %>% 
                   data.frame(., row.names = colnames(bulk)) %>% 
                   dplyr::select_(.dots = colnames(data$X)),
                 shared.donors = shared.donors,
                 train.donors = idx)
  
  saveRDS(fitted.models, paste0("BEYOND/data/celmod/celmod.model", i, ".rds"))
  saveRDS(celmod, paste0("BEYOND/data/celmod/celmod.results.", i, ".rds"))
  return(celmod)
})

# -------------------------------------------------------- #
# State proportion correlation                             #
#   snRNA-seq vs. avg Celmod prediction for test samples   #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")
celmod <- lapply(list.files("BEYOND/data/celmod", "celmod.results.*rds", full.names = TRUE), readRDS)

predicted.proportions <- lapply(celmod, function(x) 
  x$preds %>% rownames_to_column("projid") %>% 
    mutate(set = case_when(!projid %in% x$shared.donors ~ "validation",
                           projid %in% x$train.donors ~ "train",
                           .default = "test") )) %>%
  do.call(rbind, .) %>% 
  melt(id.vars=c("set","projid"), variable.name="state", value.name="prev") %>%
  group_by(set, projid, state) %>%
  summarise(n=n(), across(c(prev), list(median=median, mean=mean, sd=sd)),.groups = "drop") %>%
  mutate(sqrt.prev_mean=sqrt(prev_mean))

avg.predicted.prop <- 
  predicted.proportions %>% split(., .$set) %>% 
  lapply(., function(df) dcast(df, projid~state, value.var = "sqrt.prev_mean") %>% 
           column_to_rownames("projid") %>%
           `[`(,colnames(data$X)))

corrs <- do.call(rbind, sapply(colnames(avg.predicted.prop$test), function(s) 
  cor.test(avg.predicted.prop$test[,s], 
           data$layers[["sqrt.prev"]][rownames(avg.predicted.prop$test), s],
           use = "pairwise.complete.obs", 
           method = "spearman")[c("estimate","p.value")] %>% as.data.frame(), simplify = FALSE, USE.NAMES = TRUE )) %>%
  dplyr::select(corr=estimate, pval=p.value) %>%
  dplyr::mutate(adj.pval = p.adjust(pval, method="BH"),
                sig = cut(adj.pval, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "")))

data$uns$celmod <- list(
  shared.donors = celmod[[1]]$shared.donors,
  predicted.proportions = predicted.proportions,
  avg.predicted.prop = avg.predicted.prop,
  test.corrs = corrs, # correlations between (sqrt) average over test set predictions and actual snRNA-seq proportions
  celmod.states = corrs %>% filter(adj.pval < .01 & corr > 0) %>% rownames() # list of FDR<0.01 states
)
anndata::write_h5ad(data, "BEYOND/data/BEYOND.DLPFC.h5ad")
rm(celmod, predicted.proportions, avg.predicted.prop, corrs, data)



####################################################################################################################
##                           #  Trait associations analysis over snRNA-seq & bulk predictions                     ##
####################################################################################################################

# -------------------------------------------------------- #
# State-Trait associations                                 #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")
bulk <- py_to_r(data$uns$celmod$avg.predicted.prop$train)

bulk.metadata <- merge(load.metadata()[rownames(bulk),], 
                       read.csv("BEYOND/data/ROSMAP.bulk.RIN.values.csv", row.names = 1), 
                       by.x="row.names", by.y="row.names", all.x=TRUE)

traits <- c("sqrt.amyloid", "sqrt.amyloid_mf","sqrt.tangles","sqrt.tangles_mf","cogng_demog_slope", "ceradsc","braaksc","cogdx_ad")
sets <- list(`snuc` = list(covariates = data$layers[["sqrt.prev"]] %>% `colnames<-`(data$var_names), 
                           traits     = data$obsm$meta.data[,c(traits)] %>% mutate(across(all_of(traits), ~as.numeric(as.character(.)))),
                           controls   = data.frame(data$obsm$meta.data[,c("age_death","msex","pmi")],
                                                   data$obsm$QCs[,c("Total_Genes_Detected","Estimated_Number_of_Cells")] )),
             `celmod`=list(covariates = bulk[, data$uns$celmod$celmod.states],
                           traits     = bulk.metadata[,traits] %>% mutate(across(all_of(traits), ~as.numeric(as.character(.)))),
                           controls   = bulk.metadata[,c("age_death","msex","pmi","RIN")]))

data$uns$trait.analysis <- sapply(names(sets), function(n) 
  associate.traits(sets[[n]]$traits, sets[[n]]$covariates, sets[[n]]$controls) %>% dplyr::rename("state"="covariate"),
  simplify = F, USE.NAMES = T)

anndata::write_h5ad(data, "BEYOND/data/BEYOND.DLPFC.h5ad")
rm(data, bulk, bulk.metadata, traits, sets)


# -------------------------------------------------------- #
# State-Trait associations - meta-analysis                 #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")
data$uns$trait.analysis$meta.analysis <-
  merge(py_to_r(data$uns$trait.analysis$snuc),
        py_to_r(data$uns$trait.analysis$celmod),
        by = c("trait","state"),
        suffixes = c(".sc",".b"),
        all.x = T) %>% 
  merge(., py_to_r(data$uns$celmod$test.corrs) %>% `colnames<-`(paste0(colnames(.),".celmod")),
        by.x = "state",
        by.y = "row.names") %>% arrange(-corr.celmod) %>%
  mutate(z.meta = (tstat.sc*sqrt(n.sc) + tstat.b*sqrt(n.b)) / sqrt(n.sc + n.b),
         pval.meta = 2*pnorm(-abs(z.meta)),
         adj.pval.meta = p.adjust(pval.meta, method="BH"),
         sig.meta = cut(adj.pval.meta, c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", "")))
anndata::write_h5ad(data, "BEYOND/data/BEYOND.DLPFC.h5ad")
