library(dplyr)
library(rhdf5)
library(progress)
library(tibble)

#####################################################################################################################
#                                     Create h5 object with groups for cell types                                   #
#####################################################################################################################
safe.remove <- function(n, f=file) {
  tryCatch(h5delete(f, n),
           error = function(e) {
             warning(paste0("Could not delete [", f, ":", n, "]. Reason: ", e))
           })
}

file <- "2. Cell-type analysis/data/DLPFC.Green.atlas.h5"

if(FALSE) {
  h5createFile(file)
  
  h5createGroup(file, "glia")
  h5createGroup(file, "glia/microglia")
  h5createGroup(file, "glia/astrocytes")
  h5createGroup(file, "glia/oligodendroglia")
  h5createGroup(file, "glia/oligodendroglia/opcs")
  h5createGroup(file, "glia/oligodendroglia/oligodendrocytes")
  h5createGroup(file, "neuronal")
  h5createGroup(file, "neuronal/excitatory")
  h5createGroup(file, "neuronal/inhibitory")
  h5createGroup(file, "vascular.niche")
  
  h5write(file=file, name="index", DataFrameAsCompound=F, obj=
            stack(list(vascular.niche="vascular.niche",
                       microglia="glia/microglia",
                       astrocytes="glia/astrocytes",
                       oligodendroglia="glia/oligodendroglia",
                       opcs="glia/oligodendroglia/opcs",
                       oligodendrocytes="glia/oligodendroglia/oligodendrocytes",
                       neuronal="neuronal",
                       excitatory="neuronal/excitatory",
                       inhibitory="neuronal/inhibitory")))
}
mapping <- setNames(h5read(file, "index")$values, h5read(file, "index")$ind)

#####################################################################################################################
#                                    Create unified dataframe of cell annotations                                   #
#####################################################################################################################
# Required RAM: 106GB

df <- list()

df <- lapply(c("2. Cell-type analysis/microglia/data/microglia.h5Seurat",
               "2. Cell-type analysis/astrocytes/data/astrocytes.h5Seurat",
               "2. Cell-type analysis/inhibitory/data/inhibitory.h5Seurat",
               "2. Cell-type analysis/oligodendroglia/data/oligodendrocytes/oligodendrocytes.h5Seurat",
               "2. Cell-type analysis/oligodendroglia/data/opcs/opcs.h5Seurat",
               "2. Cell-type analysis/excitatory/data/cux2+.h5Seurat",
               "2. Cell-type analysis/excitatory/data/cux2-.h5Seurat",
               "2. Cell-type analysis/vascular.niche/data/vascular.niche.h5Seurat"), 
             function(name)
               LoadH5Seurat(name, assays = list(SCT=c("data")), misc=F, graphs=F, reductions=F, neighbors=F)@meta.data %>% 
               mutate(cell=rownames(.), projid=as.character(projid)) %>%
               select(cell, subset, class, cell.type, state, sub.population, annotation, orig.ident, batch, projid, nCount_RNA, nFeature_RNA, nCount_SCT, nFeature_SCT) 
) %>% do.call(rbind, .)

df <- df %>% mutate(
  sub.population = case_when(sub.population == "NA" ~ state, T ~ sub.population),
  cell.type = case_when(cell.type %in% c("COP", "MFOL") ~ "OPCs", T ~ cell.type),
  class = case_when(cell.type %in% c("Neutrophils", "NK Cells","CD8+ T Cells") ~ "Immune", T~class),
  grouping.by = case_when(cell.type %in% c("Monocytes", "Macrophages") ~ "Microglia",
                          class %in% c("Vascular Niche","Immune") ~ class,
                          T ~ cell.type))


safe.remove("annotations")
h5write(df, file, "annotations")


# Append neuronal subtype predicted annotations
preds <- lapply(neurons, \(v) readRDS(gsub(".h5Seurat", ".predicted.annotation.rds", v$seurat.obj))["pruned.labels"] %>% as.matrix) %>% do.call(rbind, .)
preds <- merge(df %>% select(cell, subset, grouping.by), preds %>% select(pruned.labels), by.x="cell", by.y=0) %>%
  rename(allen.labels = pruned.labels)
h5write(preds, file, "neuronal/allen.annotations")
rm(preds, df)


#####################################################################################################################
#                                      Append DE results to cell-type groups                                        #
#####################################################################################################################

de <- list(de="data/de.rds", de.pairwise="data/de.pairwise.rds")
for(ct in c("microglia", "astrocytes","inhibitory", "excitatory")) {
  for(name in names(de)) {
    group <- file.path(mapping[[ct]], name)
    safe.remove(group)
    h5write(readRDS(file.path("2. Cell-type analysis", ct, de[[name]])) %>% mutate(cluster=as.character(cluster)), file, group)
  }
}

# Append branches of vascular niche or oligodendroglia while merge information from separate files
sets <- list(vascular.niche = c("vascular", "mural", "fibroblast"),
             oligodendroglia = c("oligodendrocytes", "opcs"))

for (n in names(sets)) {
  for(name in names(de)) {
    group <- file.path(mapping[[n]], name)
    safe.remove(group)
    h5write(do.call(rbind, lapply(file.path("2. Cell-type analysis", n, "data", sets[[n]], de[[name]]), readRDS)), file, group)
  }
}
rm(de, ct, name, group, paths, n, sets)



#####################################################################################################################
#                                      Append PA results to cell-type groups                                        #
#####################################################################################################################

# Pathway clustering cutoffs
params <- list(
  microglia = list(
    Mic.1=c(.87, 1),
    Mic.2=c(.9, .975),
    Mic.3=c(.85, .9825),
    Mic.4=c(1, .9),
    Mic.5=c(.87, 1),
    Mic.6=c(.89, .88),
    Mic.7=c(.94, .96),
    Mic.8=c(1, 1),
    Mic.9=c(.939, 1),
    Mic.10=c(.935, 1),
    Mic.11=c(.87, .95),
    Mic.12=c(.95, .95),
    Mic.13=c(.98, .98),
    Mic.14=c(.97, .9),
    Mic.15=c(.87, 1),
    Mic.16=c(1, 1),
    Macrophages=c(.975, .85),
    Monocytes=c(.85, .95)
  ),
  astrocytes = list(
    Ast.1=c( .9, .85),
    Ast.2=c(.9, .95),
    Ast.3=c(.9, 1),
    Ast.4=c(.97,.925),
    Ast.5=c(.85, .93),
    Ast.6=c(1, .93),
    Ast.7=c(.875, 1),
    Ast.8=c(.85, 1),
    Ast.9=c(.9,.9),
    Ast.10=c(.9, 1)
  ),
  oligodendroglia = list(
    Oli.1=c(2,2),
    Oli.2=c(2,2),
    Oli.3=c(2,2),
    Oli.4=c(2,2),
    Oli.5=c(2,2),
    Oli.6=c(2,2),
    Oli.7=c(2,2),
    Oli.8=c(2,2),
    Oli.9=c(2,2),
    Oli.10=c(2,2),
    Oli.11=c(2,2),
    Oli.12=c(2,2),
    OPC.1=c(.9,.89),
    OPC.2=c(1,.9),
    OPC.3=c(.89,.865),
    COP=c(.9,.9),
    MFOL=c(.9,.9)
  ),
  excitatory = list(
    Exc.1=c(.85,.9),
    Exc.2=c(.85,.9),
    Exc.3=c(.85,.9),
    Exc.4=c(.85,.9),
    Exc.5=c(.9,.9),
    Exc.6=c(.9,.9),
    Exc.7=c(.85,.9),
    Exc.8=c(.9,.9),
    Exc.9=c(.8,.85),
    Exc.10=c(.85,1),
    Exc.11=c(.85,.85),
    Exc.12=c(.85,.85),
    Exc.13=c(.85,.85),
    Exc.14=c(.85,.85),
    Exc.15=c(.85,.85),
    Exc.16=c(.85,.85)
  ),
  inhibitory = list(
    Inh.1=c(.92,.92),
    Inh.2=c(.92,.9),
    Inh.3=c(.87,.89),
    Inh.4=c(.865,.88),
    Inh.5=c(.91,.87),
    Inh.6=c(.89,.95),
    Inh.7=c(.88,.95),
    Inh.8=c(.93,.89),
    Inh.9=c(.9,.95),
    Inh.10=c(.9,.95),
    Inh.11=c(.92,.94),
    Inh.12=c(.95,.94),
    Inh.13=c(.95,.9),
    Inh.14=c(.95,.95),
    Inh.15=c(.9,.9),
    Inh.16=c(1,.925)
  ),
  vascular.niche = list(
    End.1=c(.935,.95),
    End.2=c(2,.925),
    End.3=c(2,.97),
    End.4=c(.92,.94),
    End.5=c(.87,.95),
    Arteriole=c(.893,.96),
    Venule=c(1,.9),
    Peri.1=c(1,.93),
    Peri.2=c(.9,.91),
    SMC.1=c(.87,.91),
    SMC.2=c(.92,.9),
    Fib.1=c(.9,.91),
    Fib.2=c(.87,.9),
    Fib.3=c(.95,.87)
  ))

source("2. Cell-type analysis/utils/pathways.clustering.R")
sets <- list(vascular.niche = c("vascular", "mural", "fibroblast"),
             oligodendroglia = c("oligodendrocytes", "opcs"))

for(ct in c("oligodendroglia","microglia", "astrocytes","inhibitory", "excitatory", "vascular.niche")) {
  group <- file.path(mapping[[ct]], "pa")
  
  # Merge pathway analysis results of different pathway datasets
  if(!ct %in% names(sets))
    pathways <- readRDS(file.path(ct, "data/pa.rds"))
  else
    pathways <- lapply(paste0("2. Cell-type analysis/", ct, "/data/", sets[[ct]], ".pa.rds"), readRDS) %>% unlist(., recursive=FALSE)
  
  pathways <- lapply(seq_along(pathways), function(i) {
    n <- names(pathways)[i]
    pathways[[i]]@compareClusterResult %>%
      rowwise() %>%
      # Evaluate string-formatted gene and background ratios
      dplyr::mutate(dataset      = gsub("enrich", "", n),
                    cluster      = gsub("_", "\\.", cluster),
                    GeneRatio.cp = GeneRatio,
                    BgRatio.cp   = BgRatio,
                    across(c(GeneRatio, BgRatio), ~eval(parse(text=.)))) %>%
      
      # Add gene names to gene IDs and arrange names and ids by names lexicographical order
      tidyr::separate_rows(geneID, sep="/") %>%
      mutate(gene = GeneIdMapping()$names[geneID]) %>%
      group_by(across(c(-gene, -geneID))) %>%
      arrange(gene, .by_group = T) %>%
      summarise(gene = paste(gene, collapse = "/"),
                geneID = paste(geneID, collapse = "/"),
                .groups = "drop")
    
  }) %>% do.call(rbind, .) %>%
    dplyr::select(state=cluster, direction, dataset, ID, Description, gene, Count, GeneRatio, 
                  BgRatio, pvalue, p.adjust, qvalue, GeneRatio.cp, BgRatio.cp, geneID, Cluster) 
  
  # Cluster pathways within each state
  clustered.pathways <- lapply(unique(pathways$state), function(s) {
    lapply(list(list("upregulated", 1), list("downregulated", 2)), function(p) {
      cutoff = params[[ct]][[s]][[p[[2]]]]
      message("\nClustering ", p[[1]], " pathways of ", s, " - using cutoff of ", cutoff)
      
      .sub <- pathways %>% filter(state == s & direction == p[[1]]) %>% mutate(rowname = ID) %>% column_to_rownames %>% as.data.frame
      
      tryCatch({
        res <- ClusterPathways(.sub, 
                               adjacency.args = list(method="kappa"), 
                               clustering.args = list(control =list(cutoff = cutoff)),
                               rank.pathways.args = list(n.select = function(n) 3, 
                                                         attributes = .sub %>% mutate(p.adjust=-p.adjust)  %>% dplyr::select(p.adjust)))
        return(res$ranks %>% mutate(direction = p[[1]], state = s) %>% 
                 rownames_to_column("ID") %>% 
                 dplyr::select(-Description, -attributes, p.adjust.rank=attributes.rank, top.ranked=selected))
      }, error = function(e){ return(NULL) })
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
  
  
  # Append pathway clustering results and arrange records    
  pathways <- merge(pathways, clustered.pathways, 
                    by.x=c("state","direction", "ID"), 
                    by.y=c("state","direction","ID"),
                    all.x = T) %>% 
    dplyr::select(state, direction, membership, top.ranked, Description, gene, everything()) %>% 
    arrange(state, desc(direction), membership, p.adjust) %>% 
    rowwise() %>%
    mutate(top.ranked = as.numeric(top.ranked)) %>%
    ungroup()
  
  pathways$row.groups <- pathways %>%
    mutate(change.group   = rowSums(across(c(state, direction, membership), ~. != lag(.)))) %>% 
    mutate(grouping.number = cumsum(replace_na(change.group != 0, 0))) %>%
    pull(grouping.number)
  
  
  h5write(pathways, file, group)
}



for(ct in c("microglia", "astrocytes","inhibitory", "vascular.niche")) {
  group <- file.path(mapping[[ct]], "pa.pairwise")
  
  # Merge pathway analysis results of different pathway datasets
  if(!ct %in% names(sets))
    pathways <- readRDS(file.path("2. Cell-type analysis", ct, "data/pa.pairwise.rds"))
  else
    pathways <- lapply(paste0("2. Cell-type analysis/", ct, "/data/", sets[[ct]], ".pa.pairwise.rds"), readRDS) %>% unlist(., recursive=FALSE)
  
  pathways <- lapply(seq_along(pathways), function(i) {
    n <- names(pathways)[i]
    pathways[[i]]@compareClusterResult %>%
      rowwise() %>%
      # Evaluate string-formatted gene and background ratios
      dplyr::mutate(dataset      = gsub("enrich", "", n),
                    comparison   = gsub("_", "\\.", gsub("_vs_", " vs\\. ", comparison)),
                    GeneRatio.cp = GeneRatio,
                    BgRatio.cp   = BgRatio,
                    across(c(GeneRatio, BgRatio), ~eval(parse(text=.)))) %>%
      
      # Add gene names to gene IDs and arrange names and ids by names lexicographical order
      tidyr::separate_rows(geneID, sep = "/") %>%
      mutate(gene = GeneIdMapping()$names[geneID]) %>%
      group_by(across(c(-gene, -geneID))) %>%
      summarise(across(c(gene, geneID), ~paste(sort(.), collapse = "/")), .groups = "drop")
  }) %>% do.call(rbind, .)
  
  h5write(pathways, file, group)
}


#####################################################################################################################
#                                                  Obtain UMAP Coordinates                                          #
#####################################################################################################################
# Required RAM: 61GB
for(name in c("vascular.niche", "microglia", "astrocytes", "inhibitory","oligodendroglia")) {
  o <- LoadH5Seurat(paste0("2. Cell-type analysis/", name, "/data/", name, ".h5Seurat"), assays = list(SCT=c("data")), misc=F, graphs=F, neighbors=F, verbose=F)
  v <- file.path(mapping[[name]], "umap")

  safe.remove(v)
  h5write(data.frame(cell=colnames(o), group=name, state=o$state,
                     o[["umap"]]@cell.embeddings %>% `colnames<-`(c("x","y"))),
          file, v)
  rm(o, v); gc()
}
rm(name)


# ------------------------------------------------------ #
# Append other umaps                                     #
# ------------------------------------------------------ #
# Following objects and UMAP embedding are generated by code in the `2. Cell-type analysis/2.unified.umaps.R` file

# Append cell-type unified umap: reference and predictions
ref <- LoadH5Seurat("2. Cell-type analysis/data/unified.umap.h5Seurat", assays=list(SCT=c("data")), graphs=F, neighbors=F, verbose=F)
df  <- data.frame(cell=colnames(ref), state=ref$state, ref[["umap"]]@cell.embeddings %>% `colnames<-`(c("x","y")))

h5write(df, file, "umap.ref")
h5write(readRDS("2. Cell-type analysis/data/unified.umap.coordinates.rds") %>% data.frame %>% tibble::rownames_to_column(), file, "umap")
rm(ref, df)

# Append excitatory neurons umap: reference and predictions
# Required RAM: 120GB
ref <- LoadH5Seurat("2. Cell-type analysis/excitatory/data/unified.umap.h5Seurat", assays=list(SCT=c("data")), graphs=F, neighbors=F, verbose=F)
df  <- data.frame(cell=colnames(ref), group="excitatory", state=ref$state, ref[["umap"]]@cell.embeddings %>% `colnames<-`(c("x","y")))

h5write(df, file, "neuronal/excitatory/umap.ref")
h5write(readRDS("2. Cell-type analysis/excitatory/data/unified.umap.embeddings.rds") %>% `colnames<-`(c("cell","group","state","x","y")), file, "neuronal/excitatory/umap")
rm(ref, df)



#####################################################################################################################
#                                    Obtain average gene expression across all state                                #
#####################################################################################################################
summarise.expression <- function(mtx,
                                 states,
                                 donors=setNames(rep(1, length(states)), names(states)),
                                 name="",
                                 gene.splits=2) {

  gene.groups <- split(rownames(mtx), 1:gene.splits)
  barcodes    <- data.frame(state=states, donor=donors)

  # Count number of nuclei per state per donor
  barcodes <- barcodes %>% tibble::rownames_to_column() %>%
    group_by(state, donor) %>%
    mutate(n=n()) %>%
    tibble::column_to_rownames()

  # Iterate over barcodes of each state separately
  pb  <- progress_bar$new(format=paste0(name, " gene expression statistics :current/:total [:bar] :percent in :elapsed. ETA :eta"),
                          total = length(unique(states))*gene.splits, clear=F, width=100, force = T)
  exp <- list()
  for(s in unique(barcodes$state)) {

    # Then, for each subset of genes
    for(i in seq_along(gene.groups)) {

      # Subset expression matrix to genes and state subsets
      sub.mtx <- mtx[gene.groups[[i]], rownames(barcodes)[barcodes$state == s]]

      # Obtain sparse representation (row,column,expression) and append
      # gene, barcode, state and donor information
      res <- Matrix::summary(sub.mtx) %>%
        dplyr::mutate(gene = rownames(sub.mtx)[i],
                      barcode=colnames(sub.mtx)[j]) %>%
        cbind(., barcodes[.$barcode,]) %>%
        dplyr::select(-i,-j,-barcode) %>%
        `rownames<-`(NULL)

      # Compute statistics by state-donor-gene groups
      res <- res %>% dplyr::group_by(state, gene, donor, n) %>%
        dplyr::summarise(mean.exp = mean(x),
                         median.exp = median(x),
                         n.exp = n(),
                         .groups = "drop") %>%
        dplyr::mutate(pct.exp = n.exp / n,
                      group=name)

      # Remove donor column if no donors were provided
      if(length(unique(donors)) == 1)
        res$donor <- NULL

      exp[[paste(s,i)]] <- res
      rm(sub.mtx, res)
      gc()
      pb$tick()
    }
  }
  return(do.call(rbind, exp))
}

# Required RAM: 105 GB
for(n in c("vascular.niche","microglia","inhibitory","astrocytes","oligodendroglia")) {
  o <- LoadH5Seurat(paste0("2. Cell-type analysis/", n, "/data/", n, ".h5Seurat"), assays = list(SCT=c("data")), misc=F, graphs=F, reductions=F, neighbors=F, verbose=F)

  safe.remove((v <- file.path(mapping[[n]], "gene.exp")))
  message("Appending gene expression to ", v)
  h5write(summarise.expression(o[[DefaultAssay(o)]]@data, o$state, name = n), file, v)

  safe.remove((v <- file.path(mapping[[n]], "gene.exp.donor")))
  message("Appending donor gene expression to ", v)
  h5write(summarise.expression(o[[DefaultAssay(o)]]@data, o$state, o$projid, n), file, v)
  rm(o, v); gc()
}
rm(n)

# Required RAM: 187GB
objs <- list(LoadH5Seurat("2. Cell-type analysis/excitatory/data/cux2-.h5Seurat", assays = list(SCT=c("data")), misc=F, graphs=F, reductions=F, neighbors=F, verbose=F),
             LoadH5Seurat("2. Cell-type analysis/excitatory/data/cux2+.h5Seurat", assays = list(SCT=c("data")), misc=F, graphs=F, reductions=F, neighbors=F, verbose=F))

states <- unique(c(objs[[1]]$state, objs[[2]]$state))
for(by.donor in c(T, F)) {
  v <- file.path(mapping["excitatory"], ifelse(by.donor, "gene.exp.donor", "gene.exp"))
  message("Appending gene expression to ", v)
  
  res <- do.call(rbind, lapply(split(states, seq_along(states) %% 4), function(s) {
    message(s)
    o <- lapply(objs, function(o) tryCatch({subset(o, cells = colnames(o)[o$state %in% s])}, error = function(e){})) %>%
      Filter(Negate(is.null), .)
    if(length(o) == 1) o <- o[[1]]
    else o <- merge(o[[1]], o[-1])
    
    donors <- o$projid
    if(!by.donor)
      donors <- setNames(rep(1, length(o$state)), names(o$state))

    r <- summarise.expression(o[[DefaultAssay(o)]]@data, o$state, donors, "excitatory")
    rm(o); gc()
    r
  }))

  h5write(res, file, v)
  rm(res); gc()
}
rm(summarise.expression, name, objs, states, by.donor)



#####################################################################################################################
#                                  Compute module signatures for different cell types                               #
#####################################################################################################################

source("2. Cell-type analysis/utils/signatures.R")

for(ct in c("microglia","astrocytes","oligodendroglia","vascular.niche")) {
  # Load object
  o <- LoadH5Seurat(paste0("2. Cell-type analysis/", ct, "/data/", ct, ".h5Seurat"), misc=F, graphs=F, reductions=F, neighbors=F, verbose=F)
  
  # Arrange gene lists of signatures to evaluate
  signature.df <- do.call(rbind, lapply(names(markers[[ct]]), function(r) 
    stack(markers[[ct]][[r]]) %>% mutate(reference=r) %>% 
      dplyr::select(reference, signature=ind, gene=values))) %>% 
    mutate(signature = as.character(signature))
  
  # Compute signatures using Seurat's AddModuleScore
  signatures <- signature.df %>% mutate(sig=paste(reference, signature, sep = "___")) %>% unstack(gene~sig)
  df <- AddModuleScore(o, features = signatures, name = "signature")@meta.data %>% 
    dplyr::select(state, contains("signature")) %>% 
    `colnames<-`(c("state", names(signatures))) %>%
    reshape2::melt(var.ids="state") #%>% 
  
  # Test for signatures' enrichment using a wilcox text for every state and signature
  significance <- lapply(unique(df$state), function(s) {
    lapply(unique(df$variable), function(v)
      wilcox.test(df[df$variable == v & df$state == s, "value"],
                  df[df$variable == v & df$state != s, "value"],
                  alternative = "greater")[c("statistic","p.value")] %>% 
        as.data.frame %>% 
        mutate(state=s, signature=v)) %>%
      do.call(rbind, .)
  }) %>%
    do.call(rbind, .) %>%
    mutate(p.adjust = p.adjust(p.value, "BH")) %>%
    tidyr::separate(signature, c("reference", "signature"), sep = "___")
  
  # Summarize state information and append significance
  df <- df %>% tidyr::separate(variable, c("reference", "signature"), sep = "___") %>%
    group_by(state, reference, signature) %>%
    summarise(mean=mean(value), median=median(value), `pct.abv.0`=mean(value>0)) %>%
    merge(., significance %>% dplyr::rename(W=statistic), all.x = T)
  
  # Store results in h5 file
  safe.remove((v <- file.path(mapping[[ct]], "signatures")))
  h5write(list(genes = signature.df, scores = df), file, v)
  
  rm(o, signature.df, signatures, df, v); gc()
}
rm(ct)
