

DoubletNeighborParentHeatmap <- function(object, distribution = Tool(object, "RunDoubletFinder")$parent.ident.distribution, 
                                         measurement="freq", ident = "ident", normalize="row",
                                         color = viridisLite::inferno(30),
                                         average.score.color = viridisLite::viridis(30),
                                         cluster.size.color = colorRampPalette(c("white","red"))(30),
                                         cluster.rows=F, cluster.cols=F, border_color = NA, ...) {
  require(reshape2)
  
  dist <- reshape2::dcast(distribution, paste0(ident, "~doublet"), value.var = measurement, fun.aggregate = sum)
  dist <- data.frame(dist[,-1], row.names = dist[,1])
  annotation <- distribution %>% group_by_at(ident) %>% summarise(Average.Score=mean(score), Cluster.Size=n()) %>% as.data.frame()
  annotation <- data.frame(annotation[,-1], row.names = annotation[,1])
  
  cols <- order(colnames(dist))
  
  if(normalize == "row") dist <- dist/Matrix::rowSums(dist)
  else { if (normalize == "col") dist <- t(t(dist)/Matrix::rowSums(t(dist))) }
  
  return(pheatmap::pheatmap(dist[,cols], annotation_row = annotation,
                            annotation_colors = list(Average.Score=average.score.color, Cluster.Size=cluster.size.color),
                            color = color, border_color = border_color,
                            cluster_rows = cluster.rows, cluster_cols = cluster.cols, ...))
}


PlotDoubletDistributions <- function(object, 
                                     ident="ident", 
                                     remove.doublet.group = c(),
                                     file.name=file.path(object@misc$graph.path, "doublet.distribution.pdf")) {
  
  object$demux.doublet <- NA; obj$demux.singlet <- NA
  object@meta.data[object$demux.droplet.type == "DBL",]$demux.doublet <- as.character(object@meta.data[object$demux.droplet.type == "DBL", ident])
  object@meta.data[object$demux.droplet.type != "DBL",]$demux.singlet <- as.character(object@meta.data[object$demux.droplet.type != "DBL", ident])
  
  
  df <- object@misc$doub.parent.distribution %>% filter(cell %in% colnames(object)[object$demux.droplet.type != "DBL"])
  df <- data.frame(merged.clusters = as.character(Idents(object)[df$cell]), df)
  df.without <- df %>% filter(!doublet %in% c("Olig-OPC")) %>% group_by(cell) %>% mutate(tot.without.lineage=sum(n))
  
  pdf(file.name, width = 24, height = 26)
  print(plot_grid(plot_grid(DimPlot(object, group.by = "ident", label=T,raster=T) + NoLegend() + NoAxes() + labs(title="High Resolution Clustering"),
                            FeaturePlot(object, "doublet.score", order=T, raster=T) + NoAxes() + scale_color_viridis_c(option="inferno") + labs(title="DoubletFinder Score"), 
                            DimPlot(object, group.by = "demux.doublet", order=T, label=T, raster=T) + NoAxes() + NoLegend() + labs(title="Demux Doublet") + scale_color_discrete(na.value="lightgrey"),
                            DimPlot(object, group.by = "demux.singlet", order=T, label=T, raster=T) + NoAxes() + NoLegend() + labs(title="Demux Non-Doublet") + scale_color_discrete(na.value="lightgrey"), nrow=1),
                  DotPlot(object, features = markers$`Cell Types`$General) + scale_color_viridis_c() + btheme + RotatedAxis(),
                  plot_grid(DoubletNeighborParentHeatmap(NA, df, ident = "merged.clusters", cluster.rows = T, cluster.cols=T, 
                                                         main = "Simulated Doublets Parent Distribution - Normalized", angle_col = 315,
                                                         treeheight_row = 5, treeheight_col=5, silent=T)[[4]],
                            DoubletNeighborParentHeatmap(NA, df.without, ident = "merged.clusters", cluster.rows = T, cluster.cols=T, 
                                                         main = paste0("Simulated Doublets Parent Distribution - Without ", paste(remove.doublet.group, collapse = ",") , " - Normalized"), angle_col = 315,
                                                         treeheight_row = 5, treeheight_col=5, silent=T)[[4]]),
                  
                  plot_grid(DoubletNeighborParentHeatmap(NA, df, ident = "merged.clusters", cluster.rows = T, cluster.cols=T, 
                                                         main = "Simulated Doublets Parent Distribution", angle_col = 315, normalize = "NA",
                                                         treeheight_row = 5, treeheight_col=5, silent=T)[[4]],
                            DoubletNeighborParentHeatmap(NA, df.without, ident = "merged.clusters", cluster.rows = T, cluster.cols=T, normalize = "NA",
                                                         main = paste0("Simulated Doublets Parent Distribution - Without ", paste(remove.doublet.group, collapse = ",")), angle_col = 315,
                                                         treeheight_row = 5, treeheight_col=5, silent=T)[[4]]),
                  ncol=1))
  while (!is.null(dev.list()))  dev.off()
}

