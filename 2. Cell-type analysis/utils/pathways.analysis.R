#' Obtain gene name-id mappings
#'
#' @return list containing gene names given ids or gene ides given names
#' @examples 
#' GeneIdMapping()$ids[c("APOE", "TREM2","SLC38A2")]
#'
GeneIdMapping <- function() {
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  df <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db, "ENTREZID"), columns=c("SYMBOL","ENTREZID")))

  return(list(ids   = setNames(df$ENTREZID, df$SYMBOL), 
              names = setNames(df$SYMBOL, df$ENTREZID)))
}


#' Wrapper of `Seurat::FindAllMarkers` 
#' Sets multicore running pan if number of workers is specified, 
#' and appends gene ids to returned data frame returned from `Seurat::FindAllMarkers` 
FindMarkersWrapper <- function(object, idents="ident", workers=NULL, future.globals.maxSize = 1024^3, ...) {
  Idents(object) <- FetchData(object, idents)
  
  if(!is.null(workers)) {
    library(future)
    options(future.globals.maxSize = future.globals.maxSize)
    plan(multicore, workers = workers)
  }
  de <- FindAllMarkers(object, ...)
  if (nrow(de) == 0)
    return(de)
  return(de %>% dplyr::mutate(id=GeneIdMapping()$ids[gene]))
}


#' Run pathway enrichment analysis using `clusterProfiler`
#'
#' @param de Differential expression dataframe in the format required by `clusterProfiler::compareCluster`
#' @param formula Formula specifying differential gene grouping. Default is `"id~cluster"` where `id` 
#' refers to the gene id and `cluster` refers to the single cell cluster in which the gene is differentially expressed
#' @param fun List of pathways annotations to run 
#' @param universe Gene list to be used as universe for enrichment test performed by `clusterProfiler::compareCluster`
EnrichmentAnalysis <- function(
    de, 
    formula="id~cluster", 
    fun=c("enrichKEGG","enrichPathway","enrichGO:BP","enrichGO:MF","enrichGO:CC"), 
    universe=NULL) {
  if(nrow(de) == 0) return(list())
  
  suppressPackageStartupMessages(require(clusterProfiler))
  f <- as.formula(formula)
  
  pathways <- list(`enrichKEGG` = function() compareCluster(f, data = de, fun="enrichKEGG", organism="hsa", universe=universe),
                   `enrichPathway` = function() {
                     suppressPackageStartupMessages(require(ReactomePA))
                     compareCluster(f, data = de, fun="enrichPathway", universe=universe)
                   },
                   `enrichGO:BP` = function() clusterProfiler::simplify(compareCluster(f, data = de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP", universe=universe)),
                   `enrichGO:MF` = function() clusterProfiler::simplify(compareCluster(f, data = de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="MF", universe=universe)),
                   `enrichGO:CC` = function() clusterProfiler::simplify(compareCluster(f, data = de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="CC", universe=universe)))
  
  ea <- list()
  for(p in intersect(fun, names(pathways)))
    tryCatch(ea[[p]] <- pathways[[p]](), error=function(e){})
  
  for(i in seq_along(ea)) {
    # Patch to fix bug of duplicated descriptions but different IDs in Pathways package 
    problematic.descriptions <- ea[[i]]@compareClusterResult %>% 
      dplyr::select(ID, Description) %>% 
      dplyr::group_by_all() %>%
      dplyr::filter(length(unique(ID)) > 1) %>% 
      unique() %>%
      dplyr::mutate(desc_id = paste0(Description, " (ID: ", ID, ")"))
    
    ea[[i]]@compareClusterResult <- plyr::join(ea[[i]]@compareClusterResult, problematic.descriptions, by = c("Description","ID")) %>% 
      dplyr::mutate(Description = ifelse(is.na(desc_id), Description, desc_id))
  }
  return(ea)
}
