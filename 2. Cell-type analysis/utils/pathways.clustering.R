source("Cell-type analysis/utils/pathways.analysis.R")

PathwayGenes <- function(pathways = NULL, organism = c("hsa", "mmu"), funs = c("KEGG","Pathway","GO:BP","GO:MF","GO:CC")) {
  require(dplyr)
  
  organism = match.arg(organism)
  funs = match.arg(funs, several.ok = T)
  
  orgdbs = list(hsa = "org.Hs.eg.db", mmu = "org.Mm.eg.db")
  sets = list(
    KEGG = function(...) tryCatch({
      suppressMessages(clusterProfiler::download_KEGG(organism)$KEGGPATHID2EXTID)
    }, error=function(e) data.frame(a=1,b=1)[0,]),
    Pathway = function(...) tryCatch({
      suppressMessages(AnnotationDbi::select(reactome.db::reactome.db, keys=pathways, columns =c("ENTREZID"), keytype = "PATHID") %>% filter(!is.null(ENTREZID)))
    }, error=function(e) data.frame(a=1,b=1)[0,]),
    GO = function(domain, ...) {
      tryCatch({
        terms <- pathways
        if(is.null(terms))
          terms <- GO.db::GOTERM
        k = names(Filter(function(v) v == gsub("^.*:", "", domain), AnnotationDbi::Ontology(terms))) 
        suppressMessages(AnnotationDbi::mapIds(GOSemSim::load_OrgDb("org.Hs.eg.db"), 
                                               keys = k, keytype = "GOALL", 
                                               column="ENTREZID", multiVals = 'list') %>% stack() %>% dplyr::select(2,1))
      }, error=function(e) data.frame(a=1,b=1)[0,])
    })
  
  return(do.call(rbind, lapply(funs, function(n) 
    sets[[gsub(":.*$", "", n)]](n) %>% `colnames<-`(c("pathway","gene")) %>% 
      filter(is.null(pathways) | pathway %in% pathways)
    )))
}



count_words <- function(terms, ngram.n=1:2) {
  # Copied and adjusted from simplifyEnrichment::count_words and https://stackoverflow.com/a/53226662/6400526
  library(tm)
  .tokenizer <- function(x)
    unlist(lapply(ngrams(stemDocument(words(x)), ngram.n), paste, collapse = "_"), use.names = FALSE)
  
  docs = VCorpus(VectorSource(terms))
  docs = tm_map(docs, content_transformer(tolower))
  docs = tm_map(docs, removeNumbers)
  docs = tm_map(docs, removeWords, stopwords())
  docs = tm_map(docs, removePunctuation)
  docs = tm_map(docs, stripWhitespace)
  docs = tm_map(docs, removeWords, NULL)
  
  tdm = TermDocumentMatrix(
    docs,
    control = list(
      wordLengths = c(1, Inf),
      tokenize = .tokenizer,
      stemming = FALSE,
      dictionary = NULL,
      tolower = FALSE,
      weighting = function(x) weightSMART(x, "ntn")
    )
  )
  
  return(list(terms=sort(slam::col_means(tdm) %>% `names<-`(names(terms)), decreasing = T),
              words=sort(slam::row_sums(tdm) %>% `names<-`(rownames(tdm)), decreasing = T)))
}

RankClusteredPathways <- function(membership, adjacency, method=c("betweenness","closeness","hub.score","eigen_centrality","page.rank","attributes","semantics"), 
                                  attributes=NULL, descriptions=NULL, ties.method="first", 
                                  n.select = function(n) as.numeric(cut(n, c(0, 5, 10, 25, 50, Inf), c(1, 2, 4, 6, 8))),
                                  ...) {
  library(igraph)
  
  method = match.arg(method, several.ok = T)
  if(hasArg(attributes) & !"attributes" %in% method)
    method <- c(method, "attributes")
  
  if(hasArg(descriptions) & !"semantics" %in% method)
    method <- c(method, "semantics")
  
  
  g <- graph_from_adjacency_matrix(adjacency, mode="undirected", weighted = T)
  V(g)$membership <- membership
  
  # Rank pathways in each pathway cluster independently
  scores <- lapply(unique(membership), function(cluster) {
    sub  <- igraph::induced.subgraph(g, vids = names(V(g))[V(g)$membership == cluster])
    args <- modifyList(list(), list(directed=F, graph=sub))
    
    # Run different pathways' scoring methods
    do.call(cbind, sapply(method, function(m) 
      switch(m,
             "betweenness"      =, 
             "closeness"        =R.utils::doCall(match.fun(m), args=args),
             "hub.score"        =,
             "eigen_centrality" =,
             "page.rank"        =R.utils::doCall(match.fun(m), args=args)$vector,
             "attributes"       =attributes[names(V(sub)),],
             "semantics"        ={
               descs <- descriptions[names(V(sub)),] %>% `names<-`(names(V(sub)))
               count_words(descs)$terms / sapply(strsplit(descs, " "), length)
             })
    ,simplify = F)) %>% data.frame %>% mutate(membership=as.character(cluster))
  }) %>% do.call(rbind, .)
  
  score.cols <- setdiff(colnames(scores), c("membership"))
  scores %>% rownames_to_column() %>% 
    group_by(membership) %>% 
    mutate(across(score.cols, list(r= ~ base::rank(-., ties.method = ties.method)), .names="{col}.rank")) %>%
    ungroup() %>% 

    rowwise() %>%
    mutate(mean.rank = mean(c_across(paste0(score.cols,".rank")))) %>%    
    
    group_by(membership) %>% 
    mutate(selected = base::rank(mean.rank, ties.method = "first") <= n.select(n())) %>%
    ungroup() %>%
    
    column_to_rownames() %>%
    dplyr::select(membership, everything()) %>%
    `[`(rownames(adjacency),)
}


ClusterPathways <- function(pathways.gl,
                            adjacency.args = list(method = "kappa"),
                            clustering.args = list(method = "binary_cut"),
                            rank.pathways.args = list(method = c("page.rank","semantics")),
                            ...) {
  # If not provided set default argument values
  defaults <- list(
    adjacency.args = list(method = "kappa"),
    clustering.args = list(method = "binary_cut"),
    rank.pathways.args = list(method = c("page.rank","semantics"),
                              n.select = function(n) as.numeric(cut(n, c(0, 5, 10, 25, 50, Inf), c(1, 2, 4, 6, 8))))
  )
  for(n in names(defaults))
    assign(n, modifyList(defaults[[n]], get(n)))
  
  # Compute adjacency matrices
  adjacencies <- list(
    # Adjacency by DEGs
    degs = do.call(simplifyEnrichment::term_similarity, 
                   modifyList(adjacency.args, list(
                     gl = pathways.gl %>% tidyr::separate_rows(geneID, sep="/") %>% unique %>% unstack(geneID~ID))))
    ,
    # "Apriori" similarity - adjacency by pathways' genes
    genes = do.call(simplifyEnrichment::term_similarity,
                    modifyList(adjacency.args, list(
                      gl = PathwayGenes(pathways.gl$ID) %>% dplyr::select(2,1) %>%
                        filter_all(Negate(is.na)) %>% unique %>% unstack() )))
  ) %>%

    # Align rows and columns of adjacency matrices
    lapply(., function(a) {
      missing <- setdiff(pathways.gl$ID, rownames(a))
      a <- rbind(a, matrix(0, nrow = length(missing), ncol = ncol(a)) %>% `rownames<-`(missing))
      a <- cbind(a, matrix(0, ncol = length(missing), nrow = nrow(a)) %>% `colnames<-`(missing))
      
      a[is.na(a)] <- 0
      a[pathways.gl$ID, pathways.gl$ID]
    })

  # Element-wise mean of adjacency matrices
  adjacency <- base::Reduce("+", adjacencies) / 2
  adjacency[is.na(adjacency)] <- 0
  
  if(nrow(pathways.gl) == 1) {
    membership <- setNames(1:nrow(pathways.gl), rownames(pathways.gl))
    ranks <- NULL
  }
  else {
    # Cluster pathways using binary cut
    membership <- setNames(rep("-", nrow(adjacency)), rownames(adjacency))
    tryCatch(membership <- do.call(simplifyEnrichment::cluster_terms,
                                   modifyList(clustering.args, list(mat = adjacency))) %>% 
               `names<-`(rownames(adjacency)),
             error = function(e) warning(e))

    # Rank pathways within clusters
    ranks <- do.call(RankClusteredPathways,
                     modifyList(rank.pathways.args, list(membership = membership,
                                                         adjacency = adjacency,
                                                         descriptions = pathways.gl %>% dplyr::select(Description))))
    ranks <- ranks %>% mutate(Description = pathways.gl[rownames(.),"Description"])
  }

  return(list(
    membership  = membership,
    adjacencies = modifyList(adjacencies, list(joint=adjacency)),
    ranks       = ranks
  ))
}

