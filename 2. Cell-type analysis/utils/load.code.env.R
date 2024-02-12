

#####################################################################################################################
#                                           General settings and functions                                          #
#####################################################################################################################

if(!"BEYOND.env" %in% reticulate::conda_list()[,"name"]) {
  # Create conda envoriment for running python analysis (mainly scanpy+palantir) through R
  reticulate::conda_create("BEYOND.env",  packages = c("numpy","seaborn", "scikit-learn", "statsmodels", "numba", "pytables==3.6.1"))
  reticulate::conda_install("BEYOND.env", packages = c("python-igraph","leidenalg"), pip=TRUE)
  reticulate::conda_install("BEYOND.env", packages = c("scanpy","palantir","phate"), pip = TRUE)
}

reticulate::use_condaenv("BEYOND.env")
invisible(lapply(c("dplyr","tibble","reticulate","reshape2","anndata","progress","rhdf5"), library, character.only = TRUE))