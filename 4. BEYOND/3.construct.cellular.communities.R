source("2. Cell-type analysis/load.code.env.R")
source("3. Other analyses/trait.association.util.R")
source("4. BEYOND/utils/utils.R")


####################################################################################################################
##                                           #  Communities Analysis #                                            ##
####################################################################################################################

# -------------------------------------------------------- #
# State-State Correlations                                 #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
data$uns$ss.cor <- list(
  names = colnames(data$X),
  corr = stats::cor(data$X, use = "pairwise.complete.obs", method = "spearman"),
  pval = outer(1:ncol(data$X), 1:ncol(data$X), Vectorize(function(i,j)
    cor.test(data$X[,i], data$X[,j], use="pairwise.complete.obs", method = "spearman")[["p.value"]]))
)
data$uns$ss.cor$adj.pval <- matrix(p.adjust(data$uns$ss.cor$pval, method = "BH"), nrow=nrow(data$uns$ss.cor$pval))
data$uns$ss.cor$sig <- matrix(cut(data$uns$ss.cor$adj.pval, c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", "")), nrow=nrow(data$uns$ss.cor$pval))

data$uns$ss.cor$params <- list(cor.method = "spearman",
                               cor.use = "pairwise.complete.obs",
                               p.adjust.method = "BH")
anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(data)


# -------------------------------------------------------------- #
# Construct communities based on state correlations and dynamics #
# -------------------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

# ----------------------------------- #
#    High Resolution Communities      # 
# Partitioning by multiplexed Leiden  #
# ----------------------------------- #
communities <- construct.communities(
  dynamics     = py_to_r(data$uns$trajectories$palantir$dynamics$pred.vals) %>% filter(feature %in% colnames(data)), 
  correlations =  data$uns$ss.cor$corr %>% `dimnames<-`(list(data$uns$ss.cor$names, data$uns$ss.cor$names)),
  k=10,
  weights = c(1,-1,1),
  partition.type = "RBERVertexPartition")
communities$membership <- recode(communities$membership, "0"="C2", "1"="C1", "2"="C3")

# ----------------------------------- #
#     Low Resolution Commuinties      # 
# Hierarchical clustering & splitting #
# ----------------------------------- #
library(dendextend)
dends <- split(data$var_names, communities$membership) %>%
  lapply(., function(states) {
    lapply(c("dynamics.adjacency", "correlation.mtx"), function(n)
      communities[[n]][states, states]) %>%
      base::Reduce(`+`, .) %>%
      dist %>%
      hclust %>%
      dendsort::dendsort(.) %>%
      as.dendrogram() %>%
      ladderize %>%
      set("labels_to_character")

  })

sub.dends <- lapply(1:2, function(i, comm, nsplit)
  cutree(dends[[comm[[i]]]], nsplit[[i]]) %>%
    split(names(.), .) %>%
    lapply(., function(s) prune(dends[[comm[[i]]]], leaves = setdiff(labels(dends[[comm[[i]]]]), s))) %>%
    `names<-`(paste0(comm[[i]], ".", names(.))),
  comm = c("C1","C2"), nsplit = c(2,3)) %>% unlist(., recursive=F)

membership <- 
  data.frame(community = communities$membership, row.names = data$var_names) %>%
  merge(., lapply(sub.dends, labels) %>% stack, by.x="row.names", by.y="values", all.x=TRUE) %>%
  column_to_rownames("Row.names") %>% 
  dplyr::select(community, sub.community=ind) %>%
  mutate(sub.community = recode(sub.community, "C2.1"="C2.3", "C2.3"="C2.1")) %>%
  `[`(data$var_names,)

# ----------------------------------- #
#     Low Resolution Communities      # 
# Hierarchical clustering & splitting #
# ----------------------------------- #
data$uns$communities <- list(
  similarities = communities %>% `[`(c("dynamics.mtx","dynamics.adjacency","correlation.mtx")) %>% 
    `names<-`(c("dynamics","dynamics.adjacency","correlation")),
  dynamics.colnames = colnames(communities$dynamics.mtx)
)
data$var <- cbind(data$var, membership)
anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(data)



# -------------------------------------------------------------- #
# Community dynamics and trait-associations                      #
# -------------------------------------------------------------- #
data <- anndata::read_h5ad("Cell-type analysis/data/subpopulation.proportions.h5ad")
data$obsm$communities <- donor.community.proportion(data$X, data$var %>% dplyr::select(community))
data$obsm$sub.communities <- donor.community.proportion(data$X, data$var %>% dplyr::select(sub.community))

data$uns$communities$dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$palantir$pseudotime,
                                              features = cbind(data$obsm$communities, data$obsm$sub.communities %>% dplyr::select(-`NA.`)),
                                              trajectory.probs = py_to_r(data$uns$trajectories$palantir$branch.probs),
                                              trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$palantir$terminals)[,1], data$uns$trajectories$palantir$terminals$index),
                                              bootstrap=F)


data$uns$communities$trait.association <- associate.traits(
  traits = data$obsm$meta.data[,c("sqrt.amyloid_mf","sqrt.tangles_mf","cogng_demog_slope")],
  covariates = cbind(data$obsm$communities, data$obsm$sub.communities %>% dplyr::select(-`NA.`)),
  controls = data.frame(data$obsm$meta.data[,c("age_death","msex","pmi")], data$obsm$QCs[,c("Estimated_Number_of_Cells")]))

anndata::write_h5ad(data, "Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(data)



# -------------------------------------------------------------- #
# Community clustering of pathways                               #
# -------------------------------------------------------------- #
source("Cell-type analysis/utils/pathways.analysis.R")

aggregated.data <- "Cell-type analysis/DLPFC.Green.atlas.h5"
mapping <- setNames(h5read(aggregated.data, "index")$values, h5read(aggregated.data, "index")$ind)

data <- anndata::read_h5ad("Cell-type analysis/data/subpopulation.proportions.h5ad")
states.all <- data$var %>% 
  mutate(grouping.by = recode(grouping.by, "Excitatory Neurons" = "excitatory","Inhibitory Neurons" = "inhibitory",
                              "Oligodendrocytes" = "oligodendrocytes", "Astrocyte" = "astrocytes","Microglia" = "microglia", "OPCs" = "opcs",
                              "Vascular Niche" = "endo") %>% as.character())

communities <- lapply(c("community", "sub.community"), function(x) 
  split(data$var_names, data$var[[x]])) %>% unlist(recursive = F)

clustering.params <- list(
  C1 = list(up=.925, down=.9),
  C2 = list(up=.9, down=.85),
  C3 = list(up=.89, down=.85),
  C4.1 = list(up=.875, down=.9),
  C4.2 = list(up=.9, down=.9),
  C4.3 = list(up=.9, down=.9),
  C5 = list(up=.825, down=.85)
)

adjacencies <- list()
pathways <- lapply(names(clustering.params), function(comm) {
  message("Pathways clustering of community ", comm)
  # Specify states to retrieve their pathways
  states <- states.all %>% filter(community == comm | sub.community == comm) %>%
    rownames_to_column() %>%
    unstack(rowname~grouping.by)
  
  pathways <- do.call(rbind, lapply(names(states), function(ct) 
    h5read(aggregated.data, file.path(mapping[[ct]], "pa")) %>% filter(state %in% states[[ct]])))
  
  grouped.pathways <- pathways %>% 
    group_by(ID, Description, direction) %>%
    summarise(n = n(),
              across(c(Count, GeneRatio, p.adjust), mean),
              geneID = paste(sort(geneID), collapse = "/"),
              gene = paste(sort(gene), collapse = "/"),
              state = paste(sort(state), collapse = "/"),
              .groups = "drop") %>%
    rowwise() %>%
    mutate(geneID = strsplit(geneID, "/")[[1]] %>% unique %>% paste(., collapse = "/")) %>%
    as.data.frame()
  
  lapply(list(list("upregulated", clustering.params[[comm]]$up), list("downregulated", clustering.params[[comm]]$down)), function(p) {
    .sub <- grouped.pathways %>% filter(direction == p[[1]])
    res <- ClusterPathways(.sub, 
                           adjacency.args = list(method="kappa"), 
                           clustering.args = list(control =list(cutoff = p[[2]])),
                           rank.pathways.args = list(n.select = function(n) 3, 
                                                     attributes = .sub %>% mutate(p.adjust=-p.adjust)  %>% dplyr::select(p.adjust)))    
    
    adjacencies[[paste0(comm, ".", p[[1]])]] <- res$adjacencies$joint
    
    grouped.pathways %>% 
      merge(., res$ranks %>% rownames_to_column("ID") %>% dplyr::select(-Description) %>% mutate(direction = p[[1]], community = comm), 
            by.x=c("ID", "direction"), by.y=c("ID", "direction"))
  }) %>% do.call(rbind, .)
}) %>%
  do.call(rbind, .) %>%
  dplyr::select(community, direction, membership, top.ranked=selected, ID, Description, n.states=n, state, gene, everything()) %>%
  arrange(community, desc(direction), membership, desc(n.states), p.adjust)


openxlsx::write.xlsx(
  list(pathways = pathways,
       args     = lapply(names(clustering.params), function(comm) data.frame(
         community = comm,
         upregulated.cutoff=clustering.params[[comm]]$up,
         downregulated.cutoff=clustering.params[[comm]]$down)) %>% do.call(rbind, .)),
  file = "BEYOND/data/community.pathways.xlsx")
