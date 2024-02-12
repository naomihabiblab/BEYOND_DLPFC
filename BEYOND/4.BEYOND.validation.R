source("Cell-type analysis/load.code.env.R")
source("BEYOND/utils/utils.R")


####################################################################################################################
##                              #  Landscape, dynamics & communities validation in bulk  #                        ##
####################################################################################################################

# -------------------------------------------------------- #
# Creating AnnData object with embedding for bulk data     #
# -------------------------------------------------------- #
snuc.data <- anndata::read_h5ad("Cell-type analysis/data/subpopulation.proportions.h5ad")
celmod <- readRDS("Other analyses/data/celmod.predictions.rds")

# Create anndata object from celmod predictions non-overlapping donors
data <- AnnData(X = (py_to_r(celmod$avg.predicted.prop$validation)[, celmod$celmod.states])**2,
                obsm = list(meta.data = load.metadata()[rownames(py_to_r(celmod$avg.predicted.prop$validation)),]))

# 2D embedding and donors' similarities in embedding
sc <- reticulate::import("scanpy")
sc$pp$neighbors(data, n_neighbors = as.integer(10), use_rep = "X", metric = "cosine")
sc$tl$leiden(data, resolution = 1.2)
names(data$obs) <- c("clusters")

sc$tl$umap(data, maxiter = as.integer(1000), min_dist = .1, spread=2)

data$obsp <- list()
for(e in c("X_umap"))
  data$obsp[[paste0("similarity_", e)]] <- embedding.similatity(data$obsm[[e]], knn = 5)

anndata::write_h5ad(data, "BEYOND/data/Celmod.subpopulation.proportion.h5ad")
rm(e, sc, snuc.data, data, celmod)

# -------------------------------------------------------- #
# Trajectory Analysis and dynamics                         #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("BEYOND/data/Celmod.subpopulation.proportion.h5ad")

# Fitting trajectories
data$uns$trajectories <- fit.trajectories.palantir(data, dm = 5, dm.k=30, palantir.k = 15, scale=F, root.clusters = c("3","11"), use_early_cell_as_start = F)
data$uns$trajectories$branch.probs$columns <- c("prAD.like", "ABA.like")
data$uns$trajectories$terminals$index <- c("prAD.like", "ABA.like")
# Fitting dynamics
features <- data.frame(sqrt(data$X), data$obsm$meta.data[,c("cogng_demog_slope","sqrt.tangles_mf","sqrt.amyloid_mf")])

data$uns$trajectories$dynamics <- fit.dynamics(pseudotime = data$uns$trajectories$pseudotime,
                                           features = features,
                                           trajectory.probs = py_to_r(data$uns$trajectories$branch.probs),
                                           trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories$terminals)[,1], data$uns$trajectories$terminals$index),
                                           evaluate.fit = F,
                                           bootstrap = F)

data$uns$trajectories$params <- NULL
anndata::write_h5ad(data, "BEYOND/data/Celmod.subpopulation.proportion.h5ad")
rm(features)
