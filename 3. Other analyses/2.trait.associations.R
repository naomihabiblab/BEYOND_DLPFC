source("Cell-type analysis/load.code.env.R")


####################################################################################################################
##                           #  Trait associations analysis over snRNA-seq & bulk predictions                     ##
####################################################################################################################

# -------------------------------------------------------- #
# State-Trait associations                                 #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("Cell-type analysis/data/subpopulation.proportions.h5ad")
bulk <- py_to_r(data$uns$celmod$avg.predicted.prop$validation)

source("Cell-type analysis/ROSMAP.metadata.R")
bulk.metadata <- merge(load.metadata()[rownames(bulk),], 
                       read.csv("Other analyses/data/ROSMAP.bulk.RIN.values.csv", row.names = 1), 
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

anndata::write_h5ad(data, "Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(data, bulk, bulk.metadata, traits, sets)


# -------------------------------------------------------- #
# State-Trait associations - meta-analysis                 #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("Cell-type analysis/data/subpopulation.proportions.h5ad")
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
anndata::write_h5ad(data, "Cell-type analysis/data/subpopulation.proportions.h5ad")
