####################################################################################################################
##                               #  Run CellBender Over CellRanger Libraries   #                                  ##
#                                                                                                                  #
# The purpose of this document is to create a Seurat object from the output of CellBender for                      #
# a single batch. This object is analyzed for the different QCs, doublet removal and cell type                     #
# identification.                                                                                                  #
####################################################################################################################
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggridges)
btheme <- theme(axis.text = element_text(size=7), axis.title = element_text(size=7))


# -------------------------------------------------- #
# Receive file name from command-line input          #
# -------------------------------------------------- #
suppressPackageStartupMessages(library(optparse))
args <- parse_args(OptionParser(usage="%prog [options]", 
                                option_list=list(
                                  make_option(c("-n","--name"), action = "store",  help="Name of library"),
                                  make_option(c("-l","--lib"), action = "store", help="Path to counts h5 file"),
                                  make_option(c("-d","--demux.file"), action = "store",
                                              help="Path to demultiplexed_allsites_gt.best file"),
                                  make_option(c("-o", "--output"), action="store",
                                              help="Folder to place output RDS file containing Seurat object"),
                                  make_option(c("-c", "--remove.lowqual.clusters"), action="store", type='logical',
                                              help="Logical if to remove clusters of low quality cells"),
                                  make_option(c("-m", "--meta.data"), action="store", 
                                              default="1. Library preprocessing/data/ROSMAP.participants.metadata.csv",
                                              help="Path to meta data file on donors"))))
# Example of args list:
# args <- list(name="191213-B7-B", 
#              lib="191213-B7-B_filtered.h5", 
#              demux.file=".", 
#              output=".",
#              remove.lowqual.clusters=T, 
#              meta.data="1. Library preprocessing/data/snrnaseq_sample_selection_20210810.csv")

unloadNamespace("optparse")

####################################################################################################################
##                                   #  Initial Preprocessing & Annotations   #                                   ##
####################################################################################################################
message("\n############################### Library Analysis ###############################")
for (i in seq_along(args)) message(paste0("\t - ", names(args)[[i]], ": ", args[[i]])); rm(i)
message("################################################################################\n")

# --------------------------------------------------- #
# Create Seurat object with basic filtering           #
# --------------------------------------------------- #
obj <- CreateSeuratObject(Read10X_h5(args$lib), project = args$name)
obj <- RenameCells(obj, add.cell.id = obj@project.name)
obj$batch <- obj@project.name

report.info <- list()
report.info["Raw Dims"] = paste("Library contains", nrow(obj), "features across" , ncol(obj), "samples")
message(paste(Sys.time(), "\t", report.info[["Raw Dims"]]))

message(paste(Sys.time(), "\t", "Removing low-represented features and cells with no counts"))
obj <- subset(obj, features = rownames(obj)[Matrix::rowSums(obj@assays$RNA@counts) > 5],
              cells = rownames(obj@meta.data)[obj$nCount_RNA >= 1])


# --------------------------------------------------- #
# Run cell type clustering and initial analysis       #
# --------------------------------------------------- #
obj <- PredictCellType(obj, "1. Library preprocessing/data/celltype.classifier.rds")
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 30, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:30, verbose = F)
obj <- FindClusters(obj, resolution = .2, verbose = F)
obj <- RunUMAP(obj,  dims=1:30, n.components = 2)

basic.qc.plot <- plot_grid(
  DimPlot(obj, group.by = "cell.type", label=T, raster=TRUE) + NoLegend() + labs(title="Cell Types") + btheme,
  FeaturePlot(obj, "cell.type.entropy", order = T, raster=TRUE) + labs(title="Prediction Entropy") + scale_color_viridis_c(option="inferno") + btheme,
  DimPlot(obj, label=T, raster=TRUE) + NoLegend() + labs(title="Clusters") + btheme,
  plot_grid(VlnPlot(obj, features=c("nCount_RNA"), pt.size = 0) + NoLegend()+ scale_y_log10() + labs(x="", title="# Counts") + btheme + theme(title = element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")),
            VlnPlot(obj, features=c("nFeature_RNA"), pt.size = 0) + NoLegend() + labs(x="", title="# Features") + btheme + theme(title = element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")),
            VlnPlot(obj, features=c("cell.type.entropy"), pt.size = 0) + NoLegend() + labs(x="", title="Prediction Entropy") + btheme + theme(title = element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")),
            ncol=1),
  nrow=1)


# -------------------------------------------------- #
# Append Demux and meta-data information             #
# -------------------------------------------------- #
source("1. Library preprocessing/utils/ROSMAP.metadata.R")
meta.data <-  load.metadata(args$meta.data)

if(grepl("^MAP\\d+", args$name)) {
  # library of a single individual and not the result of demultiplexing
  md <- meta.data[meta.data$projid == gsub("^MAP", "", args$name),]
  md <- md %>% tibble::column_to_rownames("SampleID")
  if(nrow(md) != 1) {
    warning(paste(Sys.time(), "\tIndividual with identity", args$name, "was not found in meta data file"))
  } else {
    obj$orig.ident <- rownames(md)
    obj <- AddMetaData(obj, data.frame(md[rep(1, nrow(obj@meta.data)),], row.names = rownames(obj@meta.data)))
  }
  rm(md)
  
} else {
  # Append demultiplex information to created object
  # For some libraries, demultiplexing was performed twice with different methods. For those 
  # libraries, extract results from the `wgs_and_array` file, in which donor id is the `projid` and not `sample-id`
  demux.file <- file.path(args$demux.file, c("demultiplexed_allsites_gt.best.wgs_and_array","demultiplexed_allsites_gt.best"))

  i <- which(file.exists(demux.file))[[1]]
  ident.column <- c("projid", "SampleID")[[i]]
  demux <- read.table(demux.file[[i]], header = T, row.names = "BARCODE")

  rownames(demux) <- paste(obj@project.name, rownames(demux), sep = "_")
  if(length(  missing.barcodes <- setdiff(colnames(obj), rownames(demux))  ) > 0) {
    write.csv(missing.barcodes, file = (f <- file.path(args$output, paste0(args$name,".missing.barcodes.csv"))))
    warning(paste(Sys.time(), "\tSome barcodes appear in object but are missing from specified demultiplexed_allsites_gt.best file. Missing barcodes are saved to", f))
    rm(f)
  }
  
  demux <- demux %>% `[`(colnames(obj), ) %>%
    dplyr::select(demux.droplet.type=DROPLET.TYPE, orig.ident=BEST.GUESS, next.orig.ident=NEXT.GUESS, doublet.ident=DBL.BEST.GUESS) %>%
    mutate(orig.ident = sub("\\,.*", "", orig.ident),
           next.orig.ident = gsub(".*,(.+),.*", "\\1", next.orig.ident),
           doublet.ident = sub(",[^,]*$", "", doublet.ident))
  
  if(ident.column == "projid") {
    demux <- demux %>% mutate(across(c(orig.ident, doublet.ident), function(x) as.character(as.numeric(x)) ))
    meta.data[,ident.column] <- as.character(as.numeric(meta.data[,ident.column]))
  }
  rownames(meta.data) <- meta.data[,ident.column]
  
  message(Sys.time(),"\tAppending demux metadata")
  obj <- AddMetaData(obj, demux)
  
  message(Sys.time(),"\tAppending donor metadata")  
  obj <- AddMetaData(obj, data.frame(meta.data[obj$orig.ident, ], row.names = rownames(obj@meta.data)))
  
  rm(demux, demux.file)
}
report.info["Individuals"] <- paste(unique(obj$orig.ident), collapse = ", ")

if(exists("missing.idents")) {report.info["Missing individuals"] <- paste(missing.idents, collapse = ", "); rm(missing.idents)}
if(exists("missing.barcodes")) {report.info["Missing barcodes"] <- T; rm(missing.barcodes)}

message(paste(Sys.time(), "\tAppended meta data annotations to object from", args$meta.data))
rm(meta.data)


# -------------------------------------------------- #
# Removing low quality cells                         #
# -------------------------------------------------- #
message(paste(Sys.time(), "\tRemove low quality cells (cell-type differential threshold). Thresholds used:"))
thresholds <- readRDS("1. Library preprocessing/data/low.quality.thresholds.rds")
print(thresholds)



obj@meta.data <- obj@meta.data %>%
  mutate(over_nCount_RNA   =  nCount_RNA      >= thresholds[cell.type, "nCount_RNA"],
         over_nFeature_RNA =  Feature_RNA     >= thresholds[cell.type, "nFeature_RNA"],
         under.both        = !over_nCount_RNA & !over_nFeature_RNA,
         under.features    =  over_nCount_RNA & !over_nFeature_RNA,
         under.counts      = !over_nCount_RNA &  over_nFeature_RNA,
         over.both         =  over_nCount_RNA &  over_nFeature_RNA) %>%
  mutate(quality.text = case_when(
    under.both ~ "Under Both",
    under.features ~ "Under Features",
    under.counts ~ "Under Counts",
    T ~ NA_character_
  ))

# Create graphs of low quality cells
cell.type.low.qual.table.plot <- obj@meta.data[,c("cell.type","under.both","under.features","under.counts","over.both")] %>%
  group_by(cell.type) %>% summarise_all(sum) %>% as.data.frame()
cell.type.low.qual.table.plot <- data.frame(cell.type.low.qual.table.plot[,-1],
                                            row.names = cell.type.low.qual.table.plot[,1])

cell.type.low.qual.table.plot <-
  pheatmap::pheatmap(cell.type.low.qual.table.plot/Matrix::rowSums(cell.type.low.qual.table.plot), display_numbers = cell.type.low.qual.table.plot,
                     cluster_rows = F, cluster_cols = F, colorRampPalette(c("navy","grey99","red3"))(10), angle_col = 315, silent = T)[[4]]

cell.type.quality.distribution.plot <-
  plot_grid(ggplot(obj@meta.data, aes(x=nCount_RNA, y=cell.type, fill=over.both)) +
              geom_density_ridges(alpha=.7) +
              scale_fill_manual(values=c(`FALSE`="navy",`TRUE`="red2")) +
              scale_x_log10(oob = scales::squish_infinite) +
              theme_classic() +
              theme(legend.position = "bottom", legend.justification =c(0,1)) +
              labs(x="", y="", title="#Counts", fill="High Quality"),
            ggplot(obj@meta.data, aes(x=nFeature_RNA, y=cell.type, fill=over.both)) +
              geom_density_ridges(alpha=.7) +
              scale_fill_manual(values=c(`FALSE`="navy",`TRUE`="red2")) +
              scale_x_log10(oob = scales::squish_infinite) +
              theme_classic() +
              theme(legend.position = "bottom", legend.justification =c(0,1)) +
              labs(x="", y="", title="#Features", fill="High Quality"))



# -------------------------------------------------- #
# Expand to removal of low quality clusters          #
# -------------------------------------------------- #
# Use fitted classifier to remove clusters of low quality cells based on percentage of low
# quality cells, average cluster entropy (of cell type classification) and percentage of 
# demux doublets. Pericyte and endothelial clusters are excluded from removal.
if(args$remove.lowqual.clusters == T) {
  suppressPackageStartupMessages(library(e1071))
  message(paste(Sys.time(), "\tRemove clusters of low quality cells"))
  high.qual.cluster.classifier <- readRDS("1. Library preprocessing/data/low.quality.clusters.classifier.rds")
  
  df <- obj@meta.data %>% 
    mutate(cluster = SCT_snn_res.0.2,
           demux.droplet.type = demux.droplet.type == "DBL") %>%
    group_by(cluster) %>% 
    summarise_at(c("cell.type.entropy", "demux.droplet.type", "over.both"), mean) %>% 
    ungroup() %>%
    mutate(
      x1 = 100 - round(100*over.both, 2), 
      x2 = round(cell.type.entropy, 3),
      x3 = round(100*demux.droplet.type, 2),
      cluster = as.character(cluster)) %>%
    select(x1, x2, x3, cluster) %>%
    as.data.frame
  
  high.qual.clusters <- df[predict(high.qual.cluster.classifier, df) == "TRUE", "cluster"]
  
  # Create report plots of low quality removal
  obj$low.qual.cluster <- plyr::mapvalues(obj$SCT_snn_res.0.2, high.qual.clusters, rep(" ", length(high.qual.clusters)))
  colors <- c(` `="lightgrey", setNames(scales::hue_pal()(length(unique(obj$SCT_snn_res.0.2))), levels(obj$SCT_snn_res.0.2)))
  
  low.quality.plot <-
    plot_grid(cell.type.quality.distribution.plot,
              cell.type.low.qual.table.plot,
              DimPlot(obj, group.by = "quality.text", order=T, raster=RASTER) + scale_color_discrete(na.value="lightgrey") + labs(x="UMAP 1", y="UMAP 2", title="Low Quality Cells") + btheme,
              DimPlot(obj, group.by = "low.qual.cluster", label=T, raster=RASTER) + scale_color_manual(values=colors) + NoLegend() + labs(x="UMAP 1", y="UMAP 2", title="Low Quality Clusters") + btheme,
              nrow=1, rel_widths = c(1, 1, 1.5, 1))
  
  keep <- rownames(obj@meta.data)[(obj$over.both & obj$SCT_snn_res.0.2 %in% high.qual.clusters) | obj$cell.type %in% c("Peri", "Endo")]
  rm(df, high.qual.clusters, high.qual.cluster.classifier, colors)
} else { 
  low.quality.plot <-
    plot_grid(cell.type.quality.distribution.plot,
              cell.type.low.qual.table.plot,
              DimPlot(obj, group.by = "quality.text", order=T, raster=RASTER) + scale_color_discrete(na.value="lightgrey") + labs(x="UMAP 1", y="UMAP 2", title="Low Quality Cells") + btheme,
              nrow=1, rel_widths = c(1, 1, 1.5))

  keep <- rownames(obj@meta.data)[obj$over.both | obj$cell.type %in% c("Peri", "Endo")]
}
rm(cell.type.low.qual.table.plot, cell.type.quality.distribution.plot, thresholds)

# Remove low quality cells and clusters
obj <- subset(obj, cells = keep)
rm(keep)
report.info["Filtered Dims"] = paste("Library contains", nrow(obj), "features across" , ncol(obj), "samples")
message(paste(Sys.time(), "\t", report.info["Filtered Dims"]))
####################################################################################################################





####################################################################################################################
##                                            #  Doublet Detection   #                                            ##
####################################################################################################################
DefaultAssay(obj) <- "RNA"; obj@assays$SCT <- NULL; gc();
obj <- SCTransform(obj, variable.features.n = 2000, conserve.memory = T, verbose = F)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs = 30, verbose = F)
obj <- FindNeighbors(obj, reduction="pca", dims = 1:30, verbose = F)
obj <- FindClusters(obj, resolution = c(.2, 1.5), verbose = F)
obj <- RunUMAP(obj,  reduction = "pca", dims=1:30, verbose = F)

Idents(obj) <- obj$SCT_snn_res.1.5;
obj <- BuildClusterTree(obj, reorder = T, reorder.numeric = T);
obj$SCT_snn_res.1.5 <- Idents(obj);
Idents(obj) <- obj$SCT_snn_res.0.2

# The following DoubletFinder version is found under: https://github.com/GreenGilad/DoubletFinder
# and differs from the official version as follows:
# 1. At time of writing DoubletFinder has a high memory consumption - it generates a single large pairwise distance
#    matrix of all real- and simulated cells, resulting in O(n^2) memory, for n=#real+#simulated. The used version
#    partitions the cells into batches, calculating the pairwise distances for cells within batch compared to all cells,
#    resulting in O(n*k) memory, with k the batch.size parameter of the function. Since the time of the analysis, 
#    DoubletFinder does support batching but only of a fixed size, which might still cause a high memory consumption
#    for large objects or large number of simulated doublets
# 2. DoubletFinder simulates artificial doublets by randomly selecting pairs of cells and averaging their gene expression.
#    For datasets of unbalanced cell-types this might result in many simulated doublets in which both parent cells are
#    of the same type. DoubletFinder is not designed to be able to find such sort of doublets. Therefore, the used version
#    accepts a `parent.ident` vector, and simulates doublets with parents of different identities. To this point (early 2024)
#    the officient DoubletFinder does not support this functionality.
source("1. Library preprocessing/utils.run.doublet.finder.R")
obj <- RunDoubletFinder(obj, pN = .5, pK = 75/(1.5*ncol(obj)), sct = TRUE,
                        sim.idents = obj$cell.type, ident="SCT_snn_res.1.5", parent.ident="cell.type")

obj@misc$markers <- list(
  GABAergic=c("MEG3", "PVALB","SST","VIP", "KIT","GAD2"),
  Glutamatergic=c("SLC17A7", "RORB", "TOX","FOXP2", "CUX2"),
  OPCs=c("PCDH15", "MEGF11", "VCAN", "PDGFRA", "CSPG4"),
  Microglia=c("C3","PTPRC", "TREM2"),
  Endothelial=c("FLT1", "CLDN5","ABCB1","ATP10A"),
  Pericytes=c("PDGFRB","DCN"),
  Oligodendrocytes=c("MBP", "MOG", "MAG"),
  Astrocytes=c("CD44","SLC1A2","SLC1A3","APOE","GJA1","GFAP", "ALDH1L1")
)


# -------------------------------------------------- #
# Determine threshold for doublet score              #
# -------------------------------------------------- #
if(! "demux.droplet.type" %in% colnames(obj@meta.data)) {
  obj$is.doublet.df <- obj$doublet.score >= .35
} else {
  # Based on demux singlet/doublet classification calculate Matthews correlation coefficients for different thresholds
  # Use threshold with maximum correlation as doublet classification threshold
  library(ROCit)
  
  doublets.df <- data.frame(class=obj$demux.droplet.type, score=obj$doublet.score, row.names = rownames(obj@meta.data))[obj$demux.droplet.type %in% c("SNG", "DBL"),]
  doublets.measure <- measureit(score=doublets.df$score, class=doublets.df$class, negref = "SNG", measure = c("ACC", "FSCR", "FPR", "TPR"))
  doublets.measure$MMC <- (doublets.measure$TP*doublets.measure$TN - doublets.measure$FP * doublets.measure$FN) / 
    sqrt((doublets.measure$TP+doublets.measure$FP)*(doublets.measure$TP+doublets.measure$FN)*(doublets.measure$TN+doublets.measure$FP)*(doublets.measure$TP+doublets.measure$FN))
  
  mcc.idx <- which.max(doublets.measure$MMC)
  obj@misc$mmc <- doublets.measure$Cutoff[mcc.idx]
  
  obj$is.doublet.df <- obj$doublet.score >= obj@misc$mmc
  obj$is.doublet.demux <- obj$demux.droplet.type == "DBL"
}

cluster.doublet.prob <- 100*prop.table(table(obj$SCT_snn_res.1.5, factor(obj$is.doublet.df, levels = c(F,T))), margin = 1)
colnames(cluster.doublet.prob) <- c("Predicted Singlets", "Predicted Doublets")

obj$is.doublet.expanded <- obj$SCT_snn_res.1.5 %in% rownames(cluster.doublet.prob)[cluster.doublet.prob[,2] > 70]

# Classify doublet cells by all criteria
obj$is.doublet <- obj$is.doublet.df | obj$is.doublet.expanded
if("demux.droplet.type" %in% colnames(obj@meta.data)) obj$is.doublet <- obj$is.doublet | obj$is.doublet.demux



####################################################################################################################
##                                            #  Create Library Report   #                                        ##
####################################################################################################################
tryCatch({
if("demux.droplet.type" %in% colnames(obj@meta.data)) {
  message("\n\nCreating library report")  

  title.plot <- ggdraw() + draw_label(paste(Project(obj), "Library Analysis Report"), fontface='bold')
  report.into.plot <- tableGrob(data.frame(unlist(report.info)), cols = c(), theme = ttheme_minimal(base_size = 8))
  clean.obj.plot <- plot_grid(
    DimPlot(obj, group.by = "cell.type", label=T, raster=RASTER) + NoLegend() + labs(title="Cell Types") + btheme,
    FeaturePlot(obj, "cell.type.entropy", order = T, raster=RASTER) + labs(title="Prediction Entropy") + scale_color_viridis_c(option="inferno") + btheme,
    DimPlot(obj, label=T, raster=RASTER) + NoLegend() + labs(title="Clusters") + btheme,
    plot_grid(VlnPlot(obj, features=c("nCount_RNA"), pt.size = 0) + NoLegend()+ scale_y_log10() + labs(x="", title="# Counts") + btheme + theme(title = element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")),
              VlnPlot(obj, features=c("nFeature_RNA"), pt.size = 0) + NoLegend() + labs(x="", title="# Features") + btheme + theme(title = element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")),
              ncol=1),
    nrow=1)
  
  doublet.scores.plot <- plot_grid(
    DimPlot(obj, group.by = "SCT_snn_res.1.5", label=T, raster=RASTER) + NoLegend() + btheme,
    FeaturePlot(obj, "doublet.score", order = T, raster=RASTER) + scale_color_viridis_c(option="inferno") + btheme,
    VlnPlot(obj, group.by = "SCT_snn_res.1.5", features = "doublet.score", pt.size = 0) + NoLegend() + labs(x="") + btheme,
    nrow=1, rel_widths = c(1,1,2))
  
  
  mcc <- data.frame(x=doublets.measure$Cutoff[mcc.idx], y=doublets.measure$MMC[mcc.idx], label=paste("MCC:\n", round(obj@misc$mmc, 3)))
  stats.curve <- ggplot(melt(data.frame(do.call(cbind, doublets.measure)), id.vars = "Cutoff", 
                             measure.vars = c("ACC","FSCR","MMC", "FPR", "TPR"), variable.name = "stat"),
                        aes(x=Cutoff, y=value, color=stat)) + 
    geom_line() + 
    geom_point(aes(x=x, y=y), mcc, inherit.aes = F) + 
    geom_text(aes(x=x, y=y, label=label), mcc, hjust=-.1, inherit.aes = F) + 
    labs(title="Classification Statistics", subtitle = " ", color="Statistic") +
    theme_classic() + theme(plot.title = element_text(hjust = .5, face="bold"), plot.subtitle = element_text(hjust = .5, face="italic")) + btheme
  
  doublet.classification.plot <- plot_grid(
    stats.curve,
    DimPlot(obj, group.by = "is.doublet.demux", cols = c("navy", "red2"), order=T, raster=RASTER) + NoLegend() + btheme + 
      labs(title="Demux Doublets", caption = paste(sum(obj$is.doublet.demux), "/", ncol(obj), "nuclei")),
    DimPlot(obj, group.by = "is.doublet.df", cols = c("navy", "red2"), order=T, raster=RASTER) + NoLegend() + btheme + 
      labs(title="DoubletFinder Doublets", caption = paste(sum(obj$is.doublet.df), "/", ncol(obj), "nuclei")),
    DimPlot(obj, group.by = "is.doublet.expanded", cols = c("navy", "red2"), order=T, raster=RASTER) + NoLegend() + btheme + 
      labs(title="Cluster Expantion Doublets", caption = paste(sum(obj$is.doublet.expanded), "/", ncol(obj), "nuclei")),
    pheatmap::pheatmap(cluster.doublet.prob, display_numbers = round(cluster.doublet.prob, 1), cluster_rows = T, cluster_cols = F, 
                       treeheight_row = 0, angle_col=315, fontsize = 7, silent=T, color=colorRampPalette(c("navy","white","red2"))(40))[[4]],
    nrow=1, rel_widths = c(1,1,1,1,.4)
  )
  rm(mcc, stats.curve)

  doublets.additional.data.plot <- plot_grid(PlotMarkers(obj, features = "markers", group.by = "SCT_snn_res.1.5") + btheme + theme(legend.position = "none") + scale_color_viridis_c(),
                  DoubletNeighborParentHeatmap(obj, silent = T, main="Doublet Parent Cell Type Distribution", fontsize=8, legend=F, annotation_legend = F)[[4]],
                  nrow=1)
  
  #Patch to avoid failure in plotting below
  flags <- list(projid=0, sex=0, pAD=0, Cdx=0, braaksc=0, age_death=0, apoe_genotype=0)
  if(all(is.na(obj$projid)))        { obj$projid <- 0; flags$projid=1 }
  if(all(is.na(obj$sex)))           { obj$sex <- 0; flags$sex=1 }
  if(all(is.na(obj$pAD)))           { obj$pAD <- 0; flags$pAD=1 }
  if(all(is.na(obj$Cdx)))           { obj$Cdx <- 0; flags$Cdx=1 }
  if(all(is.na(obj$braaksc)))       { obj$braaksc <- 0; flags$braaksc=1 }
  if(all(is.na(obj$age_death)))     { obj$age_death <- 0; flags$age_death=1 }
  if(all(is.na(obj$apoe_genotype))) { obj$apoe_genotype <- 0; flags$apoe_genotype=1 }
    
  biological.details.plot <- plot_grid(
    DimPlot(obj, group.by = "projid", raster=RASTER) + btheme + NoLegend() + labs(title="Donors"),
    DimPlot(obj, group.by = "sex", raster=RASTER) + btheme + labs("Sex"),
    DimPlot(obj, group.by = "pAD", raster=RASTER) + btheme + labs(title="AD Pathology"),
    DimPlot(obj, group.by = "Cdx", raster=RASTER) + btheme + labs(title="AD Cognitive State"),
    DimPlot(obj, group.by = "braaksc", raster=RASTER) + btheme + labs(title="Braak Stage"),
    FeaturePlot(obj, features="age_death",raster=RASTER, order=T) + scale_color_viridis_c(option="inferno") + btheme + labs(title="Age of Death"),
    DimPlot(obj, group.by = "apoe_genotype", raster=RASTER) + btheme + labs(title="APOE Genotype"),
    DimPlot(obj, group.by = "demux.droplet.type", raster=RASTER) + btheme + labs(title="Demux Doublet Annotation"),
    ncol=4
  )
  
  if(flags$projid == 1)        obj$projid <- NA
  if(flags$sex == 1)           obj$sex <- NA
  if(flags$pAD == 1)           obj$pAD <- NA
  if(flags$Cdx == 1)           obj$Cdx <- NA
  if(flags$braaksc == 1)       obj$braaksc <- NA
  if(flags$age_death == 1)     obj$age_death <- NA
  if(flags$apoe_genotype == 1) obj$apoe_genotype <- NA
  
  
  pdf(file.path(args$output, paste0(args$name, ".pdf")), width = 22, height = 40)
  print(plot_grid(plotlist = list(title.plot, report.into.plot, basic.qc.plot,
                            low.quality.plot, 
                            clean.obj.plot, doublet.scores.plot,
                            doublet.classification.plot, doublets.additional.data.plot,
                            biological.details.plot), 
            rel_heights = c(.1, .3, 1, 1, 1, 1, 1,1.1, 2), ncol=1))
  while (!is.null(dev.list()))  dev.off()
  warnings()
}
})

message(paste("\n\n", Sys.time(), "\tSaving library to", (save.to <- file.path(args$output, paste0(args$name, ".seurat.rds"))) ))
saveRDS(obj, save.to)
message(paste(Sys.time(),"\tCompleted analyzing library", args$name, "successfully."))


