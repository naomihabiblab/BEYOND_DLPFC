source("BEYOND/utils.R")


####################################################################################################################
##                                           #  Communities Analysis #                                            ##
####################################################################################################################

# -------------------------------------------------------- #
# State-State Correlations                                 #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")
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
anndata::write_h5ad(data, "BEYOND/data/BEYOND.DLPFC.h5ad")
rm(data)


# -------------------------------------------------------------- #
# Construct communities based on state correlations and dynamics #
# -------------------------------------------------------------- #
data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")

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
anndata::write_h5ad(data, "BEYOND/data/BEYOND.DLPFC.h5ad")
rm(data)



# -------------------------------------------------------------- #
# Community dynamics and trait-associations                      #
# -------------------------------------------------------------- #
data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")
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

anndata::write_h5ad(data, "BEYOND/data/BEYOND.DLPFC.h5ad")
rm(data)



# -------------------------------------------------------------- #
# Community clustering of pathways                               #
# -------------------------------------------------------------- #
source("Cell-type analysis/utils/pathways.analysis.R")

data.extended <- "Cell-type analysis/DLPFC.Green.atlas.h5"
mapping <- setNames(h5read(data.extended, "index")$values, h5read(data.extended, "index")$ind)

data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")
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
    h5read(data.extended, file.path(mapping[[ct]], "pa")) %>% filter(state %in% states[[ct]])))
  
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



####################################################################################################################
##                              #  Landscape, dynamics & communities validation in bulk  #                        ##
####################################################################################################################

# -------------------------------------------------------- #
# Creating AnnData object with embedding for bulk data     #
# -------------------------------------------------------- #
snuc.data <- anndata::read_h5ad("all/data/500.h5ad")

# Create anndata object from celmod predictions non-overlapping donors
data <- AnnData(X = (py_to_r(snuc.data$uns$celmod$avg.predicted.prop$validation)[, snuc.data$uns$celmod$celmod.states])**2,
                obsm = list(meta.data = load.metadata()[rownames(py_to_r(snuc.data$uns$celmod$avg.predicted.prop$validation)),]))

# 2D embedding and donors' similarities in embedding
sc <- reticulate::import("scanpy")
sc$pp$neighbors(data, n_neighbors = as.integer(10), use_rep = "X", metric = "cosine")
sc$tl$leiden(data, resolution = 1.2)
names(data$obs) <- c("clusters")

sc$tl$umap(data, maxiter = as.integer(1000), min_dist = .1, spread=2)

data$obsp <- list()
for(e in c("X_umap"))
  data$obsp[[paste0("similarity_", e)]] <- embedding.similatity(data$obsm[[e]], knn = 5)

anndata::write_h5ad(data, "all/data/500.celmod.h5ad")
rm(e, sc, snuc.data, data)

# -------------------------------------------------------- #
# Trajectory Analysis and dynamics                         #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("all/data/500.celmod.h5ad")

# Fitting trajectories
data$uns$trajectories <- fit.trajectories.palantir(data, dm = 5, dm.k=30, palantir.k = 15, scale=F, root.clusters = c("3","11"), use_early_cell_as_start = F)
data$uns$trajectories$branch.probs$columns <- c("prAD.like", "ABA.like")
data$uns$trajectories$terminals$index <- c("prAD.like", "ABA.like")
# Fitting dynamics
features <- data.frame(sqrt(data$X), data$obsm$meta.data[,
  c("cogng_demog_slope","sqrt.tangles_mf","sqrt.amyloid_mf",
    "sqrt.tangles", "sqrt.amyloid","age_death","msex")])

data$uns$trajectories$dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$pseudotime,
                                           features = features,
                                           trajectory.probs = py_to_r(data$uns$trajectories$branch.probs),
                                           trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$terminals)[,1], data$uns$trajectories$terminals$index),
                                           evaluate.fit = F,
                                           bootstrap = F)

data$uns$trajectories$params <- NULL
anndata::write_h5ad(data, "all/data/500.celmod.h5ad")
rm(features)




####################################################################################################################
##                                     #  RNAscope validations & related analysis  #                              ##
####################################################################################################################
validations <- list(RNAscope = list())


# -------------------------------------------------------- #
# Pool RNAscope quantification from files                  #
# -------------------------------------------------------- #
files <- list(batch2 = "Iba1_DAPI_FINAL.csv") # batch1 = "Iba1_DAPI.csv" - batch had high background noise and therefore was not used
names <- c(Children_Cy3_F_Count = "TPRG1", Children_FilterObjects_cy3_F_Count="TPRG1", 
           Children_Cy5_F_Count = "MRC1",  Children_FilterObjects_cy5_F_Count="MRC1", 
           Children_Cy7_F_Count = "CPM",   Children_FilterObjects_cy7_F_Count="CPM")
path  <- "all/data/RNAscope quantification/Quantification using John Pipeline"

df <- lapply(names(files), function(batch) {
  lapply(list.files(file.path(path, batch), pattern = files[[batch]], recursive = TRUE, full.names = TRUE), function(f) {
    # Load quantification file of image
    parts <- stringr::str_match(f,".*\\/(.*)\\/(\\d+)")[,2:3]
    df <- read.csv(f) %>% mutate(ImageNumber = parts[[2]], 
                                 sample = gsub("_ctrl|-|_.*| .*", "", parts[[1]]))
    if(batch == "batch2") {
      morph <- read.csv(file.path(dirname(f), "MyExpt_Iba1_DAPI.csv")) %>% 
        dplyr::select(ObjectNumber, AreaShape_Compactness, AreaShape_Eccentricity, Parent_IBA1_final)
      df <- merge(df, morph, by.x = "Parent_Iba1_DAPI", by.y = "ObjectNumber", all.x = TRUE)
    }
    
    df %>% mutate(batch=batch,
                  AD = factor(case_when(sample %in% c("T5961","T5993","MP05187","MP04112", "MP05187","MP1141","MP08229","OC0746",
                                                      "OC1139","OC1850") ~ "No AD",
                                        sample %in% c("T5915","MP0671","OC1731") ~ "MCI",
                                        sample %in% c("T5939","T5918","T5899","T5869","T5818","T5972","T5668","MP1060","MP117",
                                                      "MP11189","MP1242","MP13101","MP14134","MP1791","OC1253") ~ "AD"),
                              levels = c("No AD", "MCI", "AD"))) %>%
      plyr::rename(names, warn_missing = FALSE) %>%
      dplyr::select(batch, sample, AD, ImageNumber, ObjectNumber, 
                    iba1=Parent_IBA1_final, TPRG1, MRC1, CPM,
                    compactness=AreaShape_Compactness, eccentricity=AreaShape_Eccentricity) %>%
      filter(!( is.na(TPRG1) | is.na(MRC1) | is.na(CPM)))
  }) %>% plyr::rbind.fill(.)
}) %>% plyr::rbind.fill(.)
validations$RNAscope$df <- df

rm(files, names, path)


# -------------------------------------------------------- #
# Obtain snRNAseq expression values of validation genes    #
# -------------------------------------------------------- #
library(Seurat)
library(SeuratDisk)

o <- LoadH5Seurat("microglia/data/microglia.h5Seurat", assays=list(SCT=c("counts")), reductions=F, graphs=F, neighbors=F, misc=F, verboth=F)
snuc.exp <- FetchData(o, vars = c("state","TPRG1","CPM","MRC1")) %>%
  filter(TPRG1+CPM+MRC1 != 0) %>%
  mutate(state = case_when(state %in% c("Mic.12","Mic.13","Macrophages") ~ state, T~"Other"))

validations$RNAscope$snuc.exp <- snuc.exp
rm(o, snuc.exp)


# -------------------------------------------------------- #
# Estimate RNAscope state prevalences in donors            #
# -------------------------------------------------------- #
validations$RNAscope$df <- df <- df %>% 
  mutate(state = case_when(TPRG1>5 &  MRC1<10 ~ "Mic.13",
                           TPRG1<2 & CPM>5 &  MRC1<5 ~ "Mic.12",
                           MRC1>=5 ~ "Macrophages",
                           MRC1<3 & TPRG1 < 2 & CPM < 2 ~ "Other",
                           .default = "none")) %>%
  filter(!is.na(state)) %>% 
  mutate(v=1) %>% 
  tidyr::spread(state, v, fill = 0)

validations$RNAscope$predicted.proportions <- df %>%
  group_by(sample, AD, batch) %>%
  summarise(across(c(Mic.12, Mic.13, Macrophages, Other), ~sum(.)/n()), ncells = n(),
            across(c(TPRG1, MRC1, CPM), mean),
            .groups = "drop") %>%
  column_to_rownames("sample")



# -------------------------------------------------------- #
# Append and summarise pTau quantification for donors      #
# -------------------------------------------------------- #
df <- list.dirs("all/data/pTAU validations/batch2") %>% # batch1 had high background noise and therefore was not used
  Filter(function(f) file.exists(file.path(f, "MyExpt_Image.csv")) ,.) %>%
  lapply(., function(f) {
    batch = basename(dirname(f))
    sample = basename(f) %>% gsub(".*/", "", .) %>% gsub("_ctrl|-|_.*| .*", "",.)
    
    read.csv(file.path(f, "MyExpt_Image.csv")) %>% 
      dplyr::select(ImageNumber, 
                    pTau_merged = AreaOccupied_AreaOccupied_pTau_merged,
                    pTau_tangles_plaques = AreaOccupied_AreaOccupied_PTau_Tangles_and_plaques,
                    pTau_filaments = AreaOccupied_AreaOccupied_pTau_other_than_tangles_and_plaques) %>%
      mutate(sample=sample, batch=batch)
  }) %>%
  do.call(rbind, .)
df$AD <- validations$RNAscope$predicted.proportions[df$sample,]$AD
validations$pTau = list(raw=df)

validations$pTau$summarised <- df %>% group_by(batch,sample, AD) %>% 
  summarise(across(c(pTau_merged, pTau_tangles_plaques, pTau_filaments), list(mean=mean, median=median))) %>% 
  data.frame



# -------------------------------------------------------- #
# Subpopulation-proportion morphology analysis             #
# -------------------------------------------------------- #
validations$RNAscope$morpholoy.association <-
  lapply(c("compactness", "eccentricity"), function(y)
    sapply(c("Mic.12","Mic.13","Macrophages","Other"), function(s)
      summary(lm(paste0(y, "~", s), validations$RNAscope$df))[["coefficients"]][2,]) %>% t %>%  # y, "~", s, "+ batch"
      `colnames<-`(c("beta", "se", "tstat", "pval")) %>% as.data.frame %>%
      mutate(adj.pval = p.adjust(pval, method = "BH"),
             sig = cut(adj.pval, c(-Inf, .0001, .001,.01,.05,1), c("****","***","**","*","")),
             trait = y) %>% 
      rownames_to_column("state")) %>%
  do.call(rbind, .) %>%
  dplyr::select(trait, state, everything())
  

# ------------------------------------------------------------- #
# Subpopulation-proportion morphology analysis - ROSMAP cohort  #
# ------------------------------------------------------------- #
data <- anndata::read_h5ad("all/data/500.h5ad")

others = paste0("Mic.", c(1:11,14:16) )
df <- data.frame(data$X[,data$var$grouping.by == "Microglia"],
                 data$obsm$meta.data[,c("mglia3_mf","mglia123_mf","pmi","sex","age_death")]) %>%
  filter(!is.na(mglia3_mf)) %>%
  mutate(Other = rowSums(dplyr::select(., all_of(others))),
         across(c(Mic.12,Mic.13,Macrophages,Other), ~sqrt(.)),
         pam = sqrt(mglia3_mf / mglia123_mf)) %>%
  dplyr::select(-all_of(c(others, "Monocytes")))

# Independent regressions of pam scores on states
res <- lapply(c("Mic.12","Mic.13","Macrophages","Other"), 
              function(s) summary(lm(paste0("pam~", s,"+pmi+age_death+sex"), df))$coefficients %>% 
                data.frame %>% `[`(2,)) %>% do.call(rbind, .)  


# Mic.12/13 associations while controlling for one another
controls <- c("pmi","age_death","sex", "Mic.12", "Mic.13")
res.ctrl <- lapply(c("Mic.12","Mic.13"), function(state) {
  .df <- df %>% mutate(covariate = get(state))
  fit <- lm(paste0("pam~", state,"+", paste(setdiff(controls, state), collapse = "+")), .df)
  summary(fit)$coefficients %>% data.frame() %>% `[`(2,) %>% `rownames<-`(paste0("ctrl.", state))
}) %>% do.call(rbind, .)

# Multiple hypothesis testing correction
res <- rbind(res, res.ctrl) %>% 
  mutate(adj.pval = p.adjust(Pr...t.., method = "BH"),
         sig = cut(adj.pval, c(-Inf, .0001, .001, .01, .05, 1), c("****","***","**","*","")))

validations$ROSMAP.PAM.Association <- res


saveRDS(validations, "all/data/RNAscope.rds")
h5write(validations, file = "all/data/500.h5","validations")
rm(validations, res, controls, res.ctrl, others,data, df)


