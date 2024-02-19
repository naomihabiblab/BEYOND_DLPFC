
subset.cell.type <- function(name, proj.name=name, lib.path="1. Library preprocessing/data/snRNA-seq libraries", file.pattern="*.seurat.rds", 
                             low.exp.thr = 15, mt.exp.thr=15, remove.doublets=F, remove.features = "^(AC\\d+{3}|AL\\d+{3}|AP\\d+{3}|LINC\\d+{3})",
                             self.merge=F, return.library.lists=F) {
  cell.types = list(microglia=c("Micr"), astrocytes=c("Astr"), endo=c("Endo", "Peri"), oligodendrocytes=c("Olig"), opcs=c("OPC"),
                    inhibitory=c("Inh"), excitatory=c("Exc"))
  
  suppressPackageStartupMessages(require(Seurat))
  suppressPackageStartupMessages(library(progress))
  
  # ----------------------------------------------------------- #
  # Subset all libraries to requested cell types                #
  # ----------------------------------------------------------- #
  message(paste0(Sys.time(), " Begin subsetting libraries for cell type ", name))
  libs <- list.files(lib.path, file.pattern, full.names = T, recursive = T)
  pb <- progress_bar$new(format="Loading libraries :current/:total [:bar] :percent in :elapsed. ETA :eta",
                         total = length(libs), clear=F, width=100, force = T)
  
  missing.libs <- c(); counts <- list(); meta.data <- list(); pred.df <- list(); doub.parents <- list(); doub.mmc <- list()
  for(p in libs) {
    o <- readRDS(p)
    if(!any(o$cell.type %in% unlist(cell.types[name]))) {
      missing.libs <- c(missing.libs, p)
      next
    }
    
    # Keep cells by type
    cell.ids <- colnames(o)[o$cell.type %in% unlist(cell.types[name])]
    
    # Keep non-doublet cells (if specified to remove)
    if(remove.doublets) cell.ids <- cell.ids[(is.na(o@meta.data[cell.ids,]$is.doublet) | o@meta.data[cell.ids,]$is.doublet == F)]
    
    # Keep cells with low MT content
    o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
    cell.ids <- cell.ids[o@meta.data[cell.ids,]$percent.mt <= mt.exp.thr]
    
    # Remove columns from pre-processing
    md.columns <- colnames(o@meta.data)[!grepl("(snn_res|seurat_clusters|tree.ident|SCT|quality.text)", colnames(o@meta.data))]
    
    # Remove unwanted features (genes)
    features <- rownames(o)[!grepl(remove.features, rownames(o))]
    
    counts[[o@project.name]] <- o@assays$RNA@counts[features, cell.ids]
    meta.data[[o@project.name]] <- o@meta.data[cell.ids, md.columns]
    pred.df[[o@project.name]] <- o@tools$PredictCellType$prediction[cell.ids,]
    doub.parents[[o@project.name]] <- o@tools$RunDoubletFinder$parent.ident.distribution %>% dplyr::filter(cell %in% cell.ids)
    doub.mmc[[o@project.name]] <- o@misc$mmc
    rm(o, cell.ids, md.columns); pb$tick()
  }
  
  if(length(missing.libs) != 0) {
    message(paste0(Sys.time(), " Warning: The following libraries have no cells of type ", name, ":", paste(missing.libs, collapse = ", ")))
  }
  rm(p, libs, missing.libs)


  # ----------------------------------------------------------- #
  # Remove lowly expressed genes                                #
  # ----------------------------------------------------------- #
  if(!is.na(low.exp.thr) & !is.null(low.exp.thr) & low.exp.thr >0) {
    message(paste0(Sys.time(), " Discarding lowly expressed genes (in less than ", low.exp.thr, " cells)"))
    discard.genes <- (do.call(rbind, lapply(counts, function(o) 
      data.frame(n=Matrix::rowSums(o != 0)) %>% rownames_to_column(var="gene"))) %>%
        dplyr::group_by(gene) %>% dplyr::summarise(n=sum(n)) %>% dplyr::filter(n < low.exp.thr))$gene
    
    for(n in names(counts)) counts[[n]] <- counts[[n]][!rownames(counts[[n]]) %in% discard.genes,]
    message(paste0(Sys.time(), " Discarded ", length(discard.genes), " lowly expressed genes"))
    rm(discard.genes, n)
  }
  
  
  columns <- unique(unlist(lapply(meta.data, colnames)))
  meta.data <- lapply(meta.data, function(md) { md[setdiff(columns, colnames(md))] <- NA; md})
  if(return.library.lists) {
    return(list(name=name, counts=counts, meta.data=meta.data, mmc=doub.mmc,
                cell.type.predictions=pred.df, doub.parent.distribution=doub.parents))
  }
  
  # ----------------------------------------------------------- #
  # Merge all to a single object                                #
  # ----------------------------------------------------------- #
  nsamples <- do.call(sum, lapply(counts, ncol))
  message(paste0(Sys.time(), " Merging libraries to a single object: ", nsamples, " cells"))
  if(self.merge) {
    # Workaround for limitation in seurat-object when merging very large matrices
    # This code was adapted from https://stackoverflow.com/a/52236148/6400526
    pb <- progress_bar$new(format=paste0("Merging libraries :current/:total (:cells.prog cells) [:bar] :percent in :elapsed. ETA :eta"),
                           total = length(counts), clear=F, width=100, force = T)
    
    i <- c(); j <- c(); x <- c(); all.colnames <- c();
    all.rownames <- sort(unique(unlist(lapply(counts, rownames))))
    while(length(counts) > 0) {
      new.row.idx <- match(rownames(counts[[1]]), all.rownames)
      values.idx <- Matrix::which(counts[[1]] > 0, arr.ind = T)
      
      i <- c(i, new.row.idx[values.idx[,1]])
      j <- c(j, max(j, 0)+ unname(values.idx[,2]))
      x <- c(x, counts[[1]]@x)
      all.colnames <- c(all.colnames, colnames(counts[[1]]))
      
      counts[[1]] <- NULL
      pb$tick(tokens=list(cells.prog=paste0(length(all.colnames), "/", nsamples)))
    }
    
    print(paste0(Sys.time()," Count matrix size: (",length(all.rownames),",", 
          length(all.colnames),") with ", length(x), " non-zero entries - sparsity:",
          round(100*(length(x)/ length(all.rownames)* length(all.colnames)),3)))
    
    obj <- Matrix::sparseMatrix(i=i, j=j, x=x,
                                dimnames = list(all.rownames, all.colnames),
                                dims = c(length(all.rownames), length(all.colnames)))
    rm(all.rownames, i, j, x, new.row.idx, values.idx)
    obj <- CreateSeuratObject(counts = obj)
  } else {
    for(n in names(counts)) counts[[n]] <- CreateSeuratObject(counts = counts[[n]])
    obj <- merge(counts[[1]], counts[-1])
  }
  
  obj@project.name <- proj.name
  obj@meta.data <- do.call(rbind, unname(meta.data))
  obj@misc$mmc <- data.frame(mmc=do.call(rbind, doub.mmc))
  obj@misc$cell.type.predictions <- do.call(rbind, unname(pred.df))
  obj@misc$doub.parent.distribution <- do.call(rbind, lapply(seq_along(doub.parents), function(i) data.frame(batch=names(doub.parents)[i], doub.parents[[i]])))
  
  # Fix bug in rownames of cell.type.predictions dataframe
  rownames(obj@misc$cell.type.predictions) <- sub(".*\\.", "", rownames(obj@misc$cell.type.predictions))
  
  message(Sys.time(), " Finished merging libraries")
  return(obj)
}

