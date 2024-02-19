source("2. Cell-type analysis/load.code.env.R")
source("4. BEYOND/utils/utils.R")


####################################################################################################################
##                           #  Trajectory Analysis Including Dynamics Fitting & Clustering                       ##
####################################################################################################################

# -------------------------------------------------------- #
# Trajectory Analysis using Palantir & VIA                 #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
data$uns$trajectories <- list()

# Run Palantir algorithm
data$uns$trajectories$palantir <- fit.trajectories.palantir(data, dm = 5, dm.k=30, 
                                                            palantir.k = 15, scale=F, 
                                                            root.clusters = c("3","4"), 
                                                            exclude.clusters=c("9","10"))
data$uns$trajectories$palantir$branch.probs$columns <- c("prAD", "ABA")
data$uns$trajectories$palantir$terminals$index <- c("prAD", "ABA")

# Run VIA algorithm
data$uns$trajectories$via <- fit.trajectories.via(data.path = "Cell-type analysis/data/subpopulation.proportions.h5ad", knn=20, 
                                                  too.big = c(.2, .075), too.small=5,
                                                  root.clusters = c("3","4"))
data$uns$trajectories$via$branch.probs$columns <- 
  data$uns$trajectories$via$branch.probs.scaled$columns <-
  data$uns$trajectories$via$terminals$index <- 
  c("ABA", "Trajectory.4", "prAD", "Trajectory.3")
anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")


# -------------------------------------------------------- #
# Fit trait/state dynamics to found trajectories           #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
features <- data.frame(data$layers[["sqrt.prev"]],
                       data$obsm$meta.data[,c("cogng_demog_slope","sqrt.tangles_mf","sqrt.amyloid_mf")], row.names = data$obs_names)

for(model in c("palantir", "via")) {
  data$uns$trajectories[[model]]$dynamics <- 
    fit.dynamics(pseudotime = data$uns$trajectories[[model]]$pseudotime,
                 features = features,
                 trajectory.probs = py_to_r(data$uns$trajectories[[model]]$branch.probs),
                 trajectory.terminal.pseudotime = setNames(py_to_r(data$uns$trajectories[[model]]$terminals)[,1], data$uns$trajectories[[model]]$terminals$index),
                 evaluate.fit = T,
                 bootstrap = F) 
}
anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(features, model)
