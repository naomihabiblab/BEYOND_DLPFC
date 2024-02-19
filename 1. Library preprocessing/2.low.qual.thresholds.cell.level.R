# The purpose of the following analysis is to find sufficient (but not necessarily tight) 
# bounds for low quality counts/features thresholds specific for each cell type.

library(Seurat)
library(dplyr)
library(ggridges)


####################################################################################################################
##                       #  Initial Preprocessing & Manual calling of low quality clusters   #                    ##
####################################################################################################################
# The libraries specified below were randomly selected (out of the multiplexed libraries).
# The following basic analysis is identical to that performed for all libraries (including the current) 
# under `3.library.analysis.R`
libs <- c("191213-B7-B", "200225-B10-A", "200302-B13-B", "200309-B17-B", "200313-B22-B", 
          "200720-B36-A", "200806-B44-B", "201007-B57-B", "201021-B60-B", "200701-B28-B")

for (lib in libs) {
  obj <- CreateSeuratObject(Read10X_h5(file.path("1. Library preprocessing/data/snRNA-seq libraries", lib, paste0(lib, "_filtered.h5"))))
  obj <- subset(obj, features = rownames(obj)[Matrix::rowSums(obj@assays$RNA@counts) > 5],
                cells = rownames(obj@meta.data)[obj$nCount_RNA >= 1])
  
  obj <- PredictCellType(obj)
  obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T)
  obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 30, verbose = F)
  obj <- FindNeighbors(obj, reduction="pca", dims = 1:30, verbose = F)
  obj <- FindClusters(obj, resolution = .2, verbose = F)
  obj <- RunUMAP(obj,  dims=1:30, n.components = 2)
  
  
  # Append demultiplex information to created object
  # For some libraries, demultiplexing was performed twice with different methods. For those 
  # libraries, extract results from the `wgs_and_array` file, in which donor id is the `projid` and not `sample-id`
  demux.file <- file.path("1. Library preprocessing/data/snRNA-seq libraries", lib, c("demultiplexed_allsites_gt.best.wgs_and_array","demultiplexed_allsites_gt.best"))
  i <- which(file.exists(demux.file))[[1]]
  obj <- AddMetaData(obj, 
                     read.table(demux.file[[i]], header = T, row.names = "BARCODE") %>%
                       `[`(colnames(obj), ) %>% 
                       dplyr::select(demux.droplet.type=DROPLET.TYPE))
  
  saveRDS(obj, file.path("1. Library preprocessing/data/low.qual.thr.libs", paste0(lib,".seurat.rds")))
  
  pdf(past0("1. Library preprocessing/graphs/low.qual.report.", lib, ".pdf"), width=10, height=5)
  plot_grid(
    DimPlot(obj, group.by = "cell.type", label=T, raster=TRUE) + NoLegend() + labs(title="Cell Types") + btheme,
    FeaturePlot(obj, "cell.type.entropy", order = T, raster=TRUE) + labs(title="Prediction Entropy") + scale_color_viridis_c(option="inferno") + btheme,
    DimPlot(obj, label=T, raster=TRUE) + NoLegend() + labs(title="Clusters") + btheme,
    plot_grid(VlnPlot(obj, features=c("nCount_RNA"), pt.size = 0) + NoLegend()+ scale_y_log10() + labs(x="", title="# Counts") + btheme + theme(title = element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")),
              VlnPlot(obj, features=c("nFeature_RNA"), pt.size = 0) + NoLegend() + labs(x="", title="# Features") + btheme + theme(title = element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")),
              VlnPlot(obj, features=c("cell.type.entropy"), pt.size = 0) + NoLegend() + labs(x="", title="Prediction Entropy") + btheme + theme(title = element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")),
              ncol=1),
    nrow=1)  
  dev.off()
}



####################################################################################################################
##                             #  Selection of thresholds for low quality cells   #                               ##
####################################################################################################################
low.qual.clusters <- list(`191213-B7-B` =c(3,4),        `200225-B10-A`=c(0,1,7,14,15,18,20,21),  `200302-B13-B`=c(0,1,12,13,14,16,20),
                          `200309-B17-B`=c(0,1,11,14),  `200313-B22-B`=c(0,1,13,14,17,18,19),    `200701-B28-B`=c(0,1,13,14),
                          `200720-B36-A`=c(0,1,14),     `200806-B44-B`=c(0,1,12),                `201007-B57-B`=c(0,8,9),     
                          `201021-B60-B`=c(0,1))

libs <- list.files("1. Library preprocessing/data/low.qual.thr.libs", ".seurat.rds", full.names = TRUE, recursive = TRUE)
df <- do.call(rbind, lapply(libs, function(p) { 
  o <- readRDS(p)
  o$high.qual <- !(o$SCT_snn_res.0.2 %in% low.qual.clusters[[o@project.name]])
  return(o@meta.data[o$demux.droplet.type != "DBL",])
}))

n=20
thresholds <- list(
  Exc = data.frame(merge(data.frame(counts=seq(400,12000, length.out = n)), data.frame(features=seq(100, 7000, length.out = n)))) %>% rownames_to_column("i"),
  Inh = data.frame(merge(data.frame(counts=seq(400,8000, length.out = n)),  data.frame(features=seq(100, 6000, length.out = n)))) %>% rownames_to_column("i"),
  Astr= data.frame(merge(data.frame(counts=seq(400,8000, length.out = n)),  data.frame(features=seq(100, 5000, length.out = n)))) %>% rownames_to_column("i"),
  Micr= data.frame(merge(data.frame(counts=seq(400,6000, length.out = n)),  data.frame(features=seq(100, 3000, length.out = n)))) %>% rownames_to_column("i"),
  Olig= data.frame(merge(data.frame(counts=seq(400,6000, length.out = n)),  data.frame(features=seq(100, 3000, length.out = n)))) %>% rownames_to_column("i"),
  OPC = data.frame(merge(data.frame(counts=seq(400,6000, length.out = n)),  data.frame(features=seq(100, 3000, length.out = n)))) %>% rownames_to_column("i"),
  Endo= data.frame(merge(data.frame(counts=seq(400,6000, length.out = n)),  data.frame(features=seq(100, 3000, length.out = n)))) %>% rownames_to_column("i"),
  Peri= data.frame(merge(data.frame(counts=seq(400,6000, length.out = n)),  data.frame(features=seq(100, 3000, length.out = n)))) %>% rownames_to_column("i")
)

chosen.thresholds <- list()
for(c in names(thresholds)) {
  ct.thr <- thresholds[[c]]
  ct.df <- df %>% filter(cell.type == c)
  
  # -------------------------------------------------- #
  # Computing threshold statistics                     #
  # -------------------------------------------------- #
  message(paste0(c, ": Computing statistics"))
  ct.res <- do.call(rbind, lapply(1:nrow(ct.thr), function(i) {
    .df <- ct.df %>% select(cell.type, high.qual, batch, nCount_RNA, nFeature_RNA) %>%
      mutate(i = i,
             counts         = ct.thr[i, "counts"],
             features       = ct.thr[i, "features"]) %>%
      mutate(above.counts   = nCount_RNA >= counts,
             above.features = nFeature_RNA >= features) %>%
      mutate(pred.pos       = above.counts & above.features) %>%
      mutate(tp             =  high.qual &  pred.pos,
             fp             = !high.qual &  pred.pos,
             tn             = !high.qual & !pred.pos,
             fn             =  high.qual & !pred.pos)
    
    .df %>% group_by(i, cell.type, counts, features) %>% 
      summarise_at(c("tp", "fp", "tn", "fn"), sum) %>% 
      mutate(f1=ifelse(counts >= features, 2*tp/(2*tp + fp + fn), NA))
  }))
  
  # -------------------------------------------------- #
  # Plotting & Selecting thresholds                    #
  # -------------------------------------------------- #
  message(paste0(c, ": Plotting figures"))
  # Plot f1 scores heatmap
  f1.scores <- 
    ggplot(ct.res, aes(x=counts, y=features, fill=f1, z=f1)) + 
    geom_raster(interpolate = F) +
    scale_fill_viridis_c(option="inferno", na.value = "black") +
    theme_classic() + 
    theme(aspect.ratio=1,
          axis.text = element_text(size=7), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x="#Counts", y="#Features", fill="Avg. F1\nScore") + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0))
  
  # Locate argmax parameters and choose one of argmaxs
  d <- layer_data(f1.scores, 1) %>% filter(!is.na(z)) %>% filter(z == max(z)) %>% 
    mutate(mean_x = mean(x), mean_y = mean(y), sd_x = sd(x), sd_y = sd(y)) %>%
    mutate(dist_mean_x = ifelse(sd_x == 0 | is.na(sd_x), 0, abs(x-mean_x)/sd_x),
           dist_mean_y = ifelse(sd_y == 0 | is.na(sd_y), 0, abs(y-mean_y)/sd_y)) %>% 
    mutate(case = dist_mean_x + dist_mean_y == min(dist_mean_x + dist_mean_y)) %>%
    group_by(case) %>% mutate(case = case & row_number() == 1) %>% arrange(case) %>% as.data.frame()
  
  # Annotate heatmap with selection
  f1.scores <-  f1.scores + 
    geom_rect(data=d, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color=case),
              inherit.aes = F, fill = NA, size = .01, show.legend = F) + 
    scale_color_manual(breaks = c(F,T), values=c("grey50", "red"))
  
  # Batch-Cell-Type distributions
  dists <- lapply(c("nCount_RNA", "nFeature_RNA"), function(f)
    ggplot(ct.df, aes_string(y="batch", x=f, fill="high.qual")) + 
      geom_density_ridges(alpha=.7) + 
      geom_segment(data=data.frame(nCount_RNA=d[d$case,"x"], nFeature_RNA=d[d$case,"y"], batch=unique(ct.df$batch)), 
                   mapping=aes_string(x=f, xend=f,
                                      y="as.numeric(as.factor(batch))", yend="as.numeric(as.factor(batch))+1.1"),
                   color="black", size=1, inherit.aes = F) +
      scale_x_log10() + 
      scale_fill_manual(values=c(`FALSE`="navy",`TRUE`="red2")) + 
      theme_classic() +
      theme(legend.position = "bottom", legend.justification =c(0,1)) +
      labs(x="", y="", fill="High Quality", title=f))
  
  # Classification confusion matrices for high and low quality cells
  cms <- plot_grid(plotlist = lapply(list(c(1,F), c(2,T)), function(h) {
    i = h[1]; h=h[2]
    .df <- ct.df[ct.df$high.qual == h,]
    cm <- matrix(c(sum(.df$nCount_RNA < d[d$case,"x"]  & .df$nFeature_RNA < d[d$case,"y"]), 
                   sum(.df$nCount_RNA >= d[d$case,"x"] & .df$nFeature_RNA < d[d$case,"y"]),
                   sum(.df$nCount_RNA < d[d$case,"x"]  & .df$nFeature_RNA >= d[d$case,"y"]),
                   sum(.df$nCount_RNA >= d[d$case,"x"] & .df$nFeature_RNA >= d[d$case,"y"])),2,2)
    
    colnames(cm) <- c("Under Features Thr.", "Over Features Thr.")
    rownames(cm) <- c("Under Counts Thr.", "Over Counts Thr.")
    prop <- round(prop.table(cm)*100,2)
    
    pheatmap::pheatmap(prop, display_numbers = matrix(paste0(cm, "\n", prop, "%"), ncol=ncol(cm)), cluster_rows = F, cluster_cols = F, silent = T, 
                       show_colnames = i==2, angle_col = 315, breaks = seq(0,100,length.out = 50),
                       colorRampPalette(c("grey97","red3"))(50), fontsize_number = 10,
                       main=ifelse(h, "High Quality Nuclei", "Low Quality Nuclei"))[[4]]
  }), ncol=1, rel_heights = c(1,1.6))
  
  # Save figures to pdf file
  pdf(paste0("1. Library preprocessing/graphs/low.quality.threshold.", c, ".pdf"), width = 20, height = 5)
  print(plot_grid(
    ggdraw() + draw_label(paste("Cell Type", c), fontface='bold'),
    plot_grid(plotlist = append(append(list(f1.scores), dists),list(cms)), nrow=1),
    rel_heights = c(.2,1), ncol=1))
  while (!is.null(dev.list()))  dev.off()
  
  chosen.thresholds[[c]] <- data.frame(cell.type = c, nCount_RNA=d[d$case,"x"], nFeature_RNA=d[d$case,"y"])
}

chosen.thresholds <- round(data.frame(do.call(rbind, chosen.thresholds)[, -1]))
saveRDS(chosen.thresholds, "1. Library preprocessing/data/low.quality.thresholds.rds")
print(chosen.thresholds)

