btheme <- theme(axis.text = element_text(size=7), axis.title = element_text(size=7))

#' Show dataset's initial doublet detection status
#'
#' Given an object with a basic clustering+visualizing analysis plot umap showing clusters, doublet
#' type and distribution of doublet type in each cluster. Plot is saved to file
#'
#' @param object Seurat object to use
#' @param path Path to save pdf to. Default is the `graphs` directory of the current object
#' @param width Width of generated pdf
#' @param height Height of generated pdf
#' 
#' @examples 
#' InitialDoubletsStatus(object)
#' 
InitialDoubletsStatus <- function(object, path=object@misc$graph.path, width=20, height=7, raster = T) {
  levs <- c("DoubletFinder & Demux", "Demux", "DoubletFinder", "Cluster Expanded")
  object@meta.data <- object@meta.data %>% 
    mutate(doublet.text = factor(case_when(
      as.logical(is.doublet.demux) & as.logical(is.doublet.df) ~ levs[1],
      as.logical(is.doublet.demux) ~ levs[2],
      as.logical(is.doublet.df) ~ levs[3],
      as.logical(is.doublet.expanded) ~ levs[4],
      T ~ NA_character_
  ), levels = levs))
  
  object@meta.data$doublet.text.full <- plyr::mapvalues(as.character(object$doublet.text), from=c(NA), to=c("Singlet"))
  prop <- round(prop.table((cm <- table(object@active.ident, object$doublet.text.full)), margin = 1)*100, 3)
  
  
  pdf(file.path(path, "initial.doublets.status.pdf"), width = width, height = height)
  print(plot_grid(DimPlot(object, label=T, raster=raster) + NoLegend() + labs(title="Clusters") + btheme,
                  DimPlot(object, group.by = "doublet.text", order=T, raster = raster) + scale_color_discrete(na.value="lightgrey") + btheme + labs(title="Doublets"),
                  pheatmap::pheatmap(prop, display_numbers = matrix(paste0(cm, "\n", prop, "%"), ncol=ncol(cm), dimnames = list(rownames(cm), colnames(cm))), silent = T,
                                     cluster_rows = F, cluster_cols = F, colorRampPalette(c("grey97","red3"))(50), angle_col = 315, main = "Doublet Distribution")[[4]],
                  nrow=1, rel_widths = c(1,1.6,1)))
  while (!is.null(dev.list()))  dev.off()
}



#' Investigate low quality cells in dataset (by clusters)
#'
#' Plot the nCounts, nFeature and expression of MT-genes percentage for a given clustering resolution.
#' In addition, plot the `Seurat::DotPlot` for cell-type markers, a heatmap over clusters of cell-type 
#' classifier predictions and doublet parents distribution
#'
#' @param object Seurat object to use
#' @param resolution Clustering resolution to show results for
#' @param width Width of generated pdf
#' @param height Height of generated pdf
#' @param filename Path and file name to save pdf to.
#' 
#' @examples 
#' LowQualityStatus(object, resolution="SCT_snn_res.XX")
#' 
LowQualityStatus <- function(object, ident, nCount.y.max=NULL, nFeature.y.max=NULL, mt.y.max=NULL, raster = T,
                             width=20, height=12, filename=file.path(object@misc$graph.path, "QC.pdf")) {
  df <- data.frame(cluster=as.character(object@meta.data[,ident]), 
                   object@misc$cell.type.predictions[colnames(object),1:8]) %>% 
    dplyr::group_by(cluster) %>% dplyr::summarise_all(mean) %>% ungroup()
  df <- data.frame(100*df[,-1], row.names = df$cluster)
  
  df2 <- object@misc$doub.parent.distribution %>% dplyr::filter(cell %in% colnames(object))
  df2 <- data.frame(cluster = object@meta.data[df2$cell, ident], df2)
  
  pdf(filename, width=width, height=height)
  print(plot_grid(
    DimPlot(object, group.by = ident, label = T, raster = raster) + NoLegend() + NoAxes(),
    plot_grid(VlnPlot(object, group.by = ident, features="nCount_RNA", pt.size = 0, ncol=1, y.max = nCount.y.max) + scale_y_log10() + NoLegend() + btheme + labs(x=NULL, y=NULL) + theme(plot.title = element_text(vjust = -2, size=10)),
              VlnPlot(object, group.by = ident, features="nFeature_RNA", pt.size = 0, ncol=1, y.max = nFeature.y.max) + NoLegend() + btheme + labs(x=NULL, y=NULL) + theme(plot.title = element_text(vjust=-2, size=10)),
              VlnPlot(object, group.by = ident, features="percent.mt", pt.size = 0, ncol=1, y.max = mt.y.max) + NoLegend() + btheme + labs(x=NULL, y=NULL) + theme(plot.title = element_text(vjust=-2, size=10)),
              VlnPlot(object, group.by = ident, features="doublet.score", pt.size = 0, ncol=1) + NoLegend() + btheme + labs(x=NULL, y=NULL) + theme(plot.title = element_text(vjust=-2, size=10)), 
              ncol = 1)))
  
  print(plot_grid(
    # data points not ordered
    FeaturePlot(object, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "doublet.score"), raster = raster, ncol=4) & DarkTheme() & NoAxes() & scale_color_viridis_c() &
      theme(plot.title = element_text(vjust = -2, size=10), legend.position = "bottom", legend.text = element_text(angle = 90, size=9)),
    # data points ordered
    FeaturePlot(object, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "doublet.score"), raster = raster, ncol=4, order=T) & DarkTheme() & NoAxes() & scale_color_viridis_c() &
      theme(plot.title = element_text(vjust = -2, size=10), legend.position = "bottom", legend.text = element_text(angle = 90, size=9)),
    ncol=1))
  
  
  cell.type.markers <- list(
    GABAergic=c("MEG3", "PVALB","SST","VIP", "KIT","GAD2"),
    Glutamatergic=c("SLC17A7", "RORB", "TOX","FOXP2", "CUX2"),
    OPCs=c("PCDH15", "MEGF11", "VCAN", "PDGFRA", "CSPG4"),
    Microglia=c("C3","PTPRC", "TREM2"),
    Endothelial=c("FLT1", "CLDN5","ABCB1","ATP10A"),
    Pericytes=c("PDGFRB","DCN"),
    Oligodendrocytes=c("MBP", "MOG", "MAG"),
    Astrocytes=c("CD44","SLC1A2","SLC1A3","APOE","GJA1","GFAP", "ALDH1L1"))
  
  print(plot_grid(DotPlot(object, group.by=ident, features=cell.type.markers) + scale_color_viridis_c() + btheme + RotatedAxis() + theme(legend.position = "bottom"),
                  plot_grid(pheatmap::pheatmap(df, display_numbers = round(df, 2), cluster_cols = F, cluster_rows = F,
                                               silent = T, color = colorRampPalette(c("white","darksalmon","red3","darkred"))(40)[1:37])[[4]],
                            NULL,
                            DoubletNeighborParentHeatmap(NA, df2, ident = "cluster", cluster.rows = F, cluster.cols=T, 
                                                         main = "Simulated Doublets Parent Distribution - Normalized", angle_col = 45,
                                                         treeheight_row = 5, treeheight_col=5, silent=T)[[4]],
                            rel_widths = c(1,.2,1.5), nrow=1),
                  ncol=1))
  while (!is.null(dev.list()))  dev.off()
}

#' Plot all clustering resolutions computed for given object
#'
#' `Seurat::DimPlot` of different clustering resolutions for specified visualizations
#'
#' @param object Seurat object to use
#' @param resolution.prefix Prefix of `meta.data` columns to plot
#' @param reductions Collection of reduction slots to plot
#' @param width Width of generated pdf
#' @param height Height of generated pdf
#' @param path Path to save pdf to. Default is the `graphs` directory of the current project
#' 
#' @examples 
#' PlotClusteringResolutions(object, resolution.prefix="SCT_snn", reductions=c("umap", "phate"))
#' 
PlotClusteringResolutions <- function(object, resolution.prefix="SCT_snn", reductions=c("umap"), raster=T,
                                      width=12, height=12, path=file.path(object@misc$graph.path, "clustering.resolutions.pdf"), ...) {
  
  pdf(path, width=width, height=height)
  for(reduction in reductions) {
    lims <- apply(Reductions(object, reduction)@cell.embeddings, 2, function(x) c(min(x), max(x)))
    
    print(DimPlot(object, 
                  group.by = colnames(object@meta.data)[grepl(resolution.prefix, colnames(object@meta.data))],
                  reduction=reduction, label=T, raster=raster, ...) & 
            NoLegend() & 
            NoAxes() & 
            lims(x=lims[1:2, 1], y=lims[1:2, 2]))
  }
  while (!is.null(dev.list()))  dev.off()
}


#' Donor-wise proportion of cells in each `ident`, split by `param`
#' 
#' @param object Seurat object to use
#' @param param parameter/feature to split by (e.g. donor sex, Braak stage, cognitive diagnosis, etc.)
#' @param ident Cell grouping for which to show cell proportion (e.g. seurat clustering resolution)
#'
#' @examples
#' # Grouped bar plot with showing for each cluster (x-axis) the donor-wise cluster proportion (y-axis), 
#' # split by donor sex (barplot color fill)
#' PlotIdentByParam(obj, "sex", "SCT_snn_res.0.1") 
PlotIdentByParam <- function(object, param, ident="ident", stat_compare=T) {
  
  df <- FetchData(object, c(ident, param, "projid")) %>%
    dplyr::filter_all(all_vars(!is.na(.))) %>%
    dplyr::group_by_at(c(ident, param, "projid")) %>%
    dplyr::summarise(f=n(), .groups = "keep") %>%
    dplyr::group_by_at(c(param, "projid")) %>%
    dplyr::mutate(f=f/sum(f))
  
  p <- ggplot(IdentByParam(object, param, ident), aes_string(x=ident, y="f")) +
    geom_boxplot(aes_string(fill=param), position=position_dodge(1), outlier.size = .3) +
    geom_point(aes_string(fill=param), position=position_dodge(1), size=.1, alpha=.5) +
    facet_wrap(paste0("~", ident) , nrow = 1, scales = "free_x") +
    labs(x=NULL, y="Mean %Cells", fill=param) +
    theme_classic() + btheme + 
    theme(strip.text = element_blank())
  
  if(stat_compare) {
    library(ggpubr)
    p <- p + stat_compare_means(aes_string(group = param), label = "p.signif", 
                                method="wilcox.test", label.y=-0.05, hide.ns = T)
  }
  return(p) 
}


ClusterTraitsCorrelation <- function(object, 
                                     ident, 
                                     traits = c("age_death", "pAD", "niareagansc", "braaksc","amyloid", "cognep_demog_slope",
                                                "dlbdx","nft", "tangles", "tdp_stage4", "cvda_4gp2", "caa_4gp"),
                                     correction.method="BH", 
                                     cor.args=list(method="spearman", exact=F)) {
  
  cluster.prop <- GetIdentFrequencies(object, ident, flat=F)
  traits <- data.frame(unique(FetchData(object, c("projid", traits))), row.names = NULL) %>% 
    dplyr::filter(!is.na(projid)) %>%
    dplyr::mutate_at(traits, as.numeric) %>% 
    column_to_rownames("projid")
  
  n <- intersect(rownames(cluster.prop), rownames(traits))
  cluster.prop <- cluster.prop[n,]
  traits <- traits[n, ]
  
  cors  <- stats::cor(traits, cluster.prop, use = "pairwise.complete.obs")
  pvals <- outer(colnames(traits), colnames(cluster.prop), Vectorize(function(i,j)
    do.call(cor.test, modifyList(cor.args, list(x=traits[,i], y=cluster.prop[,j], use="pairwise.complete.obs")))[["p.value"]]))
  asterisks <- matrix(cut(pvals, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "")), nrow=nrow(pvals))
  dimnames(asterisks) <- dimnames(pvals) <- dimnames(cors)
  
  if(!is.null(correction.method)) {
    corrected.pvals <- matrix(p.adjust(pvals, method = correction.method), nrow=nrow(pvals))
    asterisks <- matrix(cut(corrected.pvals, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "")), nrow=nrow(pvals))
    dimnames(asterisks) <- dimnames(corrected.pvals) <- dimnames(cors)
    return(list(cors = cors, pvals = pvals, corrected.pvals = corrected.pvals, asterisks = asterisks))
  } 
  
  return(list(cors = cors, pvals = pvals, asterisks = asterisks))
}


PlotClusterTraitsCorrelation <- function(object,
                                         ident,
                                         traits.correlation.args = list(), 
                                         display=c("*","corr","pval"),
                                         color=colorRampPalette(c("navy","white", "red3"))(60), 
                                         clustering_distance_cols = "correlation", 
                                         clustering_distance_rows = "correlation",
                                         color.lim = 0,
                                         transpose = F,
                                         ...) {
  display <- match.arg(display, several.ok = T)
  res     <- do.call(ClusterTraitsCorrelation, modifyList(traits.correlation.args, list(object=object, ident=ident)))
  
  display <- matrix(do.call(paste, c(lapply(display, function(disp) {
    if(disp == "*")    return(res$asterisks)
    if(disp == "corr") return(round(res$cors, 2))
    if(disp == "pval") { if("corrected.pvals" %in% names(res)) return(round(res$corrected.pvals, 2)) else return(round(res$pvals,2)) }
  }),  list(sep="\n"))), nrow=nrow(res$cors))
  
  dimnames(display) <- dimnames(res$cors)
  m=max(c(abs(res$cors), color.lim), na.rm = T)
  
  if(transpose) {res$cors <- Matrix::t(res$cors); display <- Matrix::t(display)}
  
  return(pheatmap::pheatmap(res$cors, display_numbers = display, color = color, 
                            clustering_distance_cols = clustering_distance_cols, clustering_distance_rows = clustering_distance_rows,
                            breaks=c(seq(-m, 0, length.out=floor(length(color)/2)) , seq(0.0001, m, length.out=floor(length(color)/2) )), ...))
}


ClusterDonorHeatmap <- function(object, col.param="ident", row.param="projid", cols.plot=NULL,
                    color = colorRampPalette(c("navy","white","red3"))(20),
                    row.annotations = c(),
                    clustering_distance_cols = "correlation", clustering_distance_rows = "correlation",
                    ...) {
  
  df <- FetchData(object, c(row.param, col.param)) %>%
    drop_na() %>%
    dcast(paste0(row.param, "~", col.param), fun.aggregate = length, value.var = row.param) %>% column_to_rownames(row.param)
  df <- df/Matrix::rowSums(df)
  
  if(!is.null(cols.plot)) df <- df[, cols.plot]
  
  annotation.row = NULL
  if (length(row.annotations) > 0) {
    annotation.row <- data.frame(unique(FetchData(object, c(row.param, row.annotations))), row.names = NULL) %>%
      dplyr::filter_at(c(row.param), all_vars(. %in% rownames(df))) %>%
      column_to_rownames(row.param)
  }
  
  pheatmap::pheatmap(df, show_rownames = F, color=color, annotation_row = annotation.row,
                     clustering_distance_cols = clustering_distance_cols, clustering_distance_rows = clustering_distance_rows, ...)
}



#' Basic statistics of associating cell type sub-clusters to different traits of interest
#' Of note, these associations were only used to better inform our decision on clustering resolution. 
#' The full trait associations used in the manuscript are found under the BEYOND code folder
NaitveTraitsAssociation <- function(object, ident, raster=T, width=24, height=10,
                               filename = file.path(object@misc$graph.path, paste0("trait.association.", ident, ".pdf"))) {
  object@meta.data[object@meta.data == "NA"] <- NA
  object$diff.projid <- DifferentNeighborIdentity(object)
  object$braaksc <- factor(object$braaksc, levels = 0:6, ordered = T)
  object$niareagansc <- factor(object$niareagansc, levels = 4:1, ordered=T)
  object$dlbdx <- factor(object$dlbdx)
  object$apoe_genotype <- factor(object$apoe_genotype)
  
  
  pdf(filename, width=width, height = height)
  print(plot_grid(
    plot_grid(DimPlot(object, group.by = ident, label=T, raster = raster) + NoLegend() + NoAxes(),
              PlotIdentByParam(object, "sex", ident, stat_compare = T) + theme(legend.margin=margin(0,0,0,0)),
              PlotIdentByParam(object, "pAD", ident, stat_compare = T) + theme(legend.margin=margin(0,0,0,0)),
              PlotIdentByParam(object, "apoe_genotype", ident, stat_compare = F) + theme(legend.margin=margin(0,0,0,0)),
              nrow=1, rel_widths = c(.5,1,1,1)),
    plot_grid(plotlist = lapply(c("niareagansc", "Cdx", "braaksc", "dlbdx"), function(p) 
      PlotIdentByParam(object, p, ident, stat_compare = F) + theme(legend.margin=margin(0,0,0,0))),
      nrow=2),
    ncol=1, rel_heights = c(1,2)))
  
  tryCatch(
    assign("hm", ClusterDonorHeatmap(object, ident, row.annotations = c("braaksc","niareagansc","sex"), 
                         color=colorRampPalette(c("white","salmon","red3"))(20), silent=T)[[4]]),
    error = function(e) {
      assign("hm", ClusterDonorHeatmap(object, ident, row.annotations = c("braaksc","niareagansc","sex"), 
                           color=colorRampPalette(c("white","salmon","red3"))(20), silent=T,
                           cluster_rows = T, cluster_cols=F, clustering_distance_rows = "euclidean")[[4]])
    })
  
  tryCatch({ cluster.cor <<- PlotClusterTraitsCorrelation(object, ident, silent=T)[[4]] }, 
           error = function(e) { cluster.cor <<- PlotClusterTraitsCorrelation(object, ident, cluster_rows=F, cluster_cols=F, silent=T)[[4]] })
  
  print(plot_grid(hm, cluster.cor))

  while (!is.null(dev.list()))  dev.off()
}