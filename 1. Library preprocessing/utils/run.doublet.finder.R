

#' Search for doublets using DoubletFinder
#'
#' Run DoubletFinder on given object. Store DoubletFinder's `proportion of Artificial Nearest Neighbors` as a meta.data
#' column with specified name. If requested, the simulated doublet parent names as well as parent identity distribution
#' is stored under \code{Tool(object = object, "RunDoubletFinder")}
#'
#' @param object Seurat object to run DoubletFinder over
#' @param pN proportion of simulated doublets from object after adding these doublets. 
#' @param pK proportion of dataset size to be used as number of neighbors. If NULL DoubletFinder's \code{paramSweep_v3}
#' is executed. If `k` is the desired number of neighbors then pass `k/(ncol(object)/(1-.25))`
#' @param pcs A vector PCs to use for the run of DoubletFinder. If NULL uses the `dims` logged under the `FindNeighbors`
#' command of the specified `assay` and `reduction`
#'
RunDoubletFinder <- function(object, assay = DefaultAssay(object = object), pN = .25, pK = NULL, pcs = NULL, sct = F, 
                             score.name="doublet.score", reduction="PCA", batch.size=5000, sim.idents = NULL,
                             compute.parent.ident.distribution=T, ...) {
  require(DoubletFinder); require(data.table)
  
  if(is.null(pcs)) {
    pcs <- Command(object, paste("FindNeighbors", assay, reduction, sep = "."))$dims
  }
  if (is.null(pK)) {
    sweep.res.list <- DoubletFinder::paramSweep_v3(object, PCs = pcs, sct = sct)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <-DoubletFinder::find.pK(sweep.stats)
    pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), 
    ]$pK))
  }
  
  df <- doubletFinder_v3(object, PCs = pcs, pN = pN, pK = pK, nExp = 0, idents = sim.idents,
                         reuse.pANN = F, sct = sct, batch.size = batch.size, 
                         get.neighbor.doublets = compute.parent.ident.distribution)
  object@meta.data[,score.name] <- df@meta.data[, grep("pANN", colnames(df@meta.data))]
  
  if(compute.parent.ident.distribution) {
    message("Computing Doublet Neighbors Parents Ident Distribution...")
    tool <- df@tools$doubletFinder_v3
    
    tool$neighbor.doublets <- tool$neighbor.doublets[sapply(tool$neighbor.doublets, function(l) length(l) != 0)]
    n <- names(tool$neighbor.doublets)
    tool$neighbor.doublets <- as.data.frame(rbindlist(lapply(tool$neighbor.doublets, function(l) data.table(t(l))), fill=T))
    rownames(tool$neighbor.doublets) <- n
    
    tool$parent.ident.distribution <- ComputeDoubletNeighborParentsDistribution(object, tool, score.column = score.name, ...)
    Tool(object) <- tool
  }
  
  object <- Seurat::LogSeuratCommand(object)
  return(object)
}