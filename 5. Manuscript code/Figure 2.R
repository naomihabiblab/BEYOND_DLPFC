source("5. Manuscript code/utils.R")


#####################################################################################################################
#                                          Figure 2 - Cell-State Characterization                                   #
#####################################################################################################################

# ----------------------------------------------------------------------------------------------------------------- #
#                                                     Cell-type UMAPS                                               #
# ----------------------------------------------------------------------------------------------------------------- #
plot.umap <- function(name, cols=NULL) {
  embedding <- h5read(aggregated.data, name)
  st        <- embedding$state %>% unique()
  
  # If not specified randomly assign colors to states
  if(is.null(cols)) {
    cols <- unlist(lapply(split(1:length(st), 1:3), length))
    cols <- setNames(c(colorRange("darkgreen")(cols[[1]]), colorRange("darkorchid4")(cols[[2]]), colorRange("midnightblue")(cols[[3]])), sample(st))
  }
  
  # Annotate plot with number of nuclei 
  text <- annotations[embedding$cell,] %>% 
    group_by(grouping.by) %>% summarise(n=n()) %>% 
    mutate(label=paste0(grouping.by," n=", n)) %>% 
    pull(label) %>% paste0(collapse = "\n")
  pos  <- apply(embedding[c("x","y")], 2, min)
  
  return(ggplot(embedding, aes(x, y, color=state, label=state)) +
           scale_color_manual(values=cols) +
           ggrastr::geom_point_rast(size=.01, raster.dpi = 600) +
           geom_text(data = embedding %>% group_by(state) %>% summarise_at(vars(x, y), list(mean)), color="black") +
           annotate(geom = "text", label=text, x=pos[1], y=pos[2], hjust=0)+
           no.labs +
           theme_embedding +
           theme(legend.position = "none"))
}

pdf(file.path(panel.path, "2.UMAPs.pdf"), width=embed.width, height = embed.height)
for(ct in c("vascular.niche","microglia","astrocytes","oligodendroglia","inhibitory","excitatory","neuronal")) {
  plot.umap(file.path(mapping[[ct]], ifelse(ct %in% c("excitatory","neuronal"), "umap.ref","umap")), 
            cols = state.colors[[ct]]) %>% print(.)
}
while (!is.null(dev.list()))  dev.off()
rm(ct, plot.umap)


# ----------------------------------------------------------------------------------------------------------------- #
#                                                         State QCs                                                 #
# ----------------------------------------------------------------------------------------------------------------- #


pdf(file.path(panel.path, "2.QCs.pdf"), width=embed.width, height=embed.height)
for(ct in names(atlas)) {
  o <- LoadH5Seurat(paste0(ct,"/data/", ct, ".H5Seurat"), assays=list(SCT=c("data")), misc=F, graphs=F, reductions=F, neighbors=F, verbose=F)
  Idents(o) <- factor(as.character(Idents(o)), levels=atlas[[ct]]$state.order)
  
  if(any(is.na(Idents(o))))
    o <- subset(o, cells = colnames(o)[!is.na(Idents(o))])
  
  plot_grid(
    VlnPlot(o, features=c("nCount_RNA","nFeature_RNA"), pt.size = 0, ncol = 1, log = TRUE, cols = state.colors[[ct]]) &
      labs(x=NULL, y=NULL, title=NULL) &
      NoLegend() & 
      theme(axis.text.x = element_blank()),

    melt(data$X) %>% 
      `colnames<-`(c("donor","state","prev")) %>% 
      filter(state %in% levels(Idents(o))) %>% 
      mutate(state = factor(state, levels = levels(Idents(o)))) %>%
      ggplot(aes(state, prev, fill=state, color=state)) + 
      geom_boxplot(alpha=.3, outlier.shape = NA) +
      geom_jitter(size=.5, width=.2, alpha=.75) + 
      scale_color_manual(values=state.colors[[ct]]) + 
      scale_fill_manual(values=state.colors[[ct]]) +
      scale_y_sqrt(breaks = c(.01, .05, .1, .25, .5, .75), labels = paste0(as.integer(100*c(.01, .05, .1, .25, .5, .75)), "%"), expand = expansion(0)) +
      labs(x=NULL, y=NULL) + 
      theme_classic() + 
      theme(legend.position = "none"),
    ncol=1,
    rel_heights = c(1,1)) %>% print(.)
  
  rm(o); gc()
}
while (!is.null(dev.list()))  dev.off()
rm(ct)



# ----------------------------------------------------------------------------------------------------------------- #
#                                                       Signature scores                                            #
# ----------------------------------------------------------------------------------------------------------------- #

# Removing Sadick et al annotations in favor of their integrated annotations
excluded.signatures = list(dummy=c("ABC"), `Sadick J.S. et al (2022)`=c())

pdf(file.path(panel.path, "2.signatures.pdf"), height=embed.height*1.25, width = embed.width*1.5)
for(ct in c("vascular.niche","microglia","astrocytes","oligodendroglia")) {
  
  sigs <- h5read(aggregated.data, file.path(mapping[[ct]], "signatures"))
  sigs$scores <- sigs$scores %>% mutate(sig.full=paste(reference, signature)) %>%
    # Remove specified signatures
    merge(., stack(excluded.signatures) %>% `colnames<-`(c("signature","reference")) %>% mutate(to.exclude=T), all.x=T) %>% 
    filter(is.na(to.exclude)) %>%
    dplyr::select(-to.exclude) %>%
     # Remove all signatures from specified references
    filter(! reference %in% names(Filter(is.null, excluded.signatures)))
  
  labels <- sigs$scores %>% dplyr::select(sig.full, signature) %>% unique %>% `rownames<-`(NULL) %>% column_to_rownames("sig.full")
  references <- sigs$scores %>% dplyr::select(sig.full, reference) %>% unique %>% `rownames<-`(NULL) %>% column_to_rownames("sig.full")
  
  
  mtx <- sigs$scores %>% pivot_wider(id_cols = "state", names_from = "sig.full", values_from = "mean") %>% 
    column_to_rownames("state") %>%
    `[`(atlas[[ct]]$state.order, ) %>%
    scale
  
  significance <- sigs$scores %>% 
    rowwise() %>%
    mutate(sig = ifelse(p.adjust <= .01, "*", ""))  %>%
    pivot_wider(id_cols = "state", names_from = "sig.full", values_from = "sig", values_fill = "") %>%
    column_to_rownames("state") %>%
    `[`(rownames(mtx), colnames(mtx)) %>%
    as.matrix()
  
  v <- max(abs(mtx), na.rm = T)
  draw(Heatmap(mtx,
          col = circlize::colorRamp2(seq(-v, v, length.out=21), green2purple(21)),
          cluster_rows=F,
          row_split = atlas[[ct]]$main.group,
          show_column_dend = F,
          column_labels = labels[colnames(mtx),],
          column_names_rot = 45,
          row_title = " ",
          layer_fun = function(j, i, x, y, w, h, fill) grid.text(pindex(significance, i, j), x, y),
          bottom_annotation = HeatmapAnnotation(ref=references[colnames(mtx),]),
          height = unit(18, "cm"),
          width = unit(14, "cm"),
          border=T,
          row_gap = unit(0,"pt")
          ), merge_legend=T)
}
while (!is.null(dev.list()))  dev.off()
rm(ct, sigs, labels, references, mtx, significance, v,excluded.signatures)


# ----------------------------------------------------------------------------------------------------------------- #
#                                                  Signature gene expression                                        #
# ----------------------------------------------------------------------------------------------------------------- #

# Removing Sadick et al annotations in favor of their integrated annotations
excluded.signatures = list(dummy=c("ABC"), `Sadick J.S. et al (2022)`=c())

pdf(file.path(panel.path, "2.signatures.genes.pdf"), height=embed.height*3, width = embed.width*1.5)
for(ct in c("vascular.niche","microglia","astrocytes","oligodendroglia")) {
  
  # Load signatures' gene sets
  sigs <- h5read(aggregated.data, file.path(mapping[[ct]], "signatures"))
  sigs <- sigs$genes %>% mutate(section=paste(reference, signature), col=paste(reference, signature, gene)) %>%
    # Remove specified signatures
    merge(., stack(excluded.signatures) %>% `colnames<-`(c("signature","reference")) %>% mutate(to.exclude=T), all.x=T) %>% 
    filter(is.na(to.exclude)) %>%
    dplyr::select(-to.exclude) %>%
    # Remove all signatures from specified references
    filter(! reference %in% names(Filter(is.null, excluded.signatures)))
  
  de  <- h5read(aggregated.data, file.path(mapping[[ct]], "de")) %>% 
    filter(avg_log2FC > 0 & p_val_adj < .01 & gene %in% unique(sigs$gene)) %>%
    mutate(sig = "*") %>%
    dplyr::select(cluster, gene, sig)
  exp <- h5read(aggregated.data, file.path(mapping[[ct]], "gene.exp")) %>% filter(gene %in% unique(sigs$gene))
  
  
  # Append expression levels and if gene is DE in state
  df <- sigs %>% 
    merge(., exp, by.x = "gene", by.y="gene", all.x = T) %>%
    merge(., de,  by.x = c("state","gene"), by.y=c("cluster","gene"), all.x = T) %>%
    mutate(sig = replace_na(sig, " "))

  
  mtx <- pivot_wider(df %>% filter(!is.na(state)), id_cols = "state", names_from = "col", values_from = "mean.exp", values_fill = 0) %>% 
    column_to_rownames("state") %>%
    `[`(atlas[[ct]]$state.order, intersect(sigs$col,colnames(.))) %>%
    as.matrix %>%
    scale %>%
    t
  
  significance <- pivot_wider(df %>% filter(!is.na(state)), id_cols = "state", names_from = "col", values_from = "sig", values_fill = "") %>% 
    column_to_rownames("state") %>%
    `[`(atlas[[ct]]$state.order, intersect(sigs$col,colnames(.))) %>%
    t
  
  labeling <- sigs %>% column_to_rownames("col") %>% mutate(section = gsub(" .* ", " - ", section)) %>% `[`(intersect(sigs$col, rownames(mtx)),)
  gaps <- labeling %>% dplyr::select(reference, section) %>% unique() %>% mutate(gap = ifelse(reference == lag(reference), 0, 3)) %>% pull(gap) %>% Filter(Negate(is.na), .)
  draw(Heatmap(mtx, 
               col = green2purple(21),
               cluster_columns = F,
               column_labels = case_when(ct == "vascular.niche"~colnames(mtx), T~gsub("^.*\\.", "", colnames(mtx))),
               column_names_rot = 0,
               row_labels = labeling$gene,
               row_split = labeling$section,
               cluster_rows = T,
               cluster_row_slices = F,
               show_row_dend = F,
               layer_fun = function(j, i, x, y, w, h, fill) grid.text(pindex(significance, i, j), x, y),
               left_annotation = rowAnnotation(ref = labeling$reference),
               heatmap_legend_param = list(title = "scaled mean exp."),
               border=T,
               row_gap = unit(gaps, "pt"),
               width = unit(5,"cm"),
               height= unit(20, "cm")))

}
while (!is.null(dev.list()))  dev.off()
rm(ct,  sigs, de, exp, df, mtx, significance, labeling)



# ----------------------------------------------------------------------------------------------------------------- #
#                                          Dotplots of top and selected DEGs                                        #
# ----------------------------------------------------------------------------------------------------------------- #
n.genes = list(inhibitory=4, excitatory=4, endo=4, astrocytes=4, oligodendroglia=4, microglia=3)

pdf(file.path(panel.path, "2.dotplots.pdf"), height=embed.height, width = embed.width*2)
for(ct in names(atlas)) {
  
  additional.genes = c()
  if("genes" %in% names(atlas[[ct]]))
    additional.genes = unname(unlist(atlas[[ct]]$genes))
  
  genes <- h5read(aggregated.data, file.path(mapping[[ct]], "de")) %>%
    filter(avg_log2FC > 0) %>%
    group_by(gene) %>%
    slice_max(avg_log2FC, n=1) %>%
    group_by(cluster) %>% 
    arrange(desc(avg_log2FC)) %>%
    filter(row_number() <= n.genes[[ct]] | gene %in% additional.genes) %>% 
    unstack(gene~cluster) %>% 
    as.list
  
  hm <- h5read(aggregated.data, file.path(mapping[[ct]], "gene.exp")) %>%
    filter(gene %in% (genes[atlas[[ct]]$state.order] %>% unlist)) %>%
    gheatmap(color.by="mean.exp", size.by="pct.exp", 
             row.order = list(a=atlas[[ct]]$state.order),
             column.order = genes[atlas[[ct]]$state.order],
             cluster_columns=FALSE, column_gap = unit(0, "mm"),
             cols = green2purple,
             border=TRUE,
             cluster_rows=FALSE, scale="col")
  
  draw(hm[[1]], merge_legend=TRUE, annotation_legend_list=hm[[2]])
  }
while (!is.null(dev.list()))  dev.off()
rm(ct, genes, hm, additional.genes, n.genes)





# ----------------------------------------------------------------------------------------------------------------- #
#                                Dotplots of top and selected DEGs - Oligodendroglia specific                       #
# ----------------------------------------------------------------------------------------------------------------- #

# Retrieve oligodendroglia DEGs + pseudobulk gene expression
de  <- h5read(aggregated.data, file.path(mapping[["oligodendroglia"]], "de")) %>% filter(avg_log2FC > 0)
exp <- h5read(aggregated.data, file.path(mapping[["oligodendroglia"]], "gene.exp"))

# Create general OPC vs Oligo heatmap
genes <- list(OPCs=c("PCDH15", "VCAN", "PDGFRA"),
              Oligodendrocytes=c("MBP", "MOG", "MAG"))

ct.hm <- exp %>% filter(gene %in% unlist(genes)) %>%
  gheatmap(color.by="mean.exp", size.by="pct.exp", 
           name = "Cell-type marker",
           row.order = list(a=atlas$oligodendroglia$state.order),
           column.order = genes,
           cluster_columns=FALSE, column_gap = unit(0, "mm"),
           cols = green2purple,
           border=TRUE,
           cluster_rows=FALSE, scale="col")

# Create within OPCs and within Oligos heatmaps
states <- atlas$oligodendroglia$state.order %>% split(., if_else(grepl("Oli.", .),"Oligodendrocytes", "OPCs"))
additional.genes <- list(
  OPCs = c("PINK1", "TOMM20", "TOMM7", "SERPINA3", "OSMR", "CNTN2", "SOX10", "SOX6"),
  Oligodendrocytes = c("QDPR", "DPYD", "S100A6", "SEMA3C", "SLC38A2", "DNAJB1", "DNAJB6", "DNAJC1", "HSPA1A", "HSPA1B", "HSPA4L", "HSPH1", "PTGES3"))

hms <- lapply(names(states), \(ct) {
  genes <- de %>% 
    filter(cluster %in% states[[ct]]) %>%
    group_by(gene) %>%
    slice_max(avg_log2FC, n=1) %>%
    group_by(cluster) %>% 
    arrange(desc(avg_log2FC)) %>%
    filter(row_number() <= 3 | gene %in% additional.genes[[ct]]) %>% 
    unstack(gene~cluster) %>% as.list
  
  exp %>% filter(gene %in% unlist(genes)) %>%
    gheatmap(color.by="mean.exp", size.by="pct.exp", 
             name = paste0(ct, "\nmean exp."),
             row.order = list(a=atlas$oligodendroglia$state.order),
             column.order = genes[states[[ct]]],
             cluster_columns=FALSE, column_gap = unit(0, "mm"),
             cols = green2purple,
             border=TRUE,
             cluster_rows=FALSE, scale="col")
})


pdf(file.path(panel.path, "2.dotplots.oligodendroglia.specific.pdf"), height=embed.height, width = embed.width*2)
draw(ct.hm[[1]] + hm[[2]][[1]] + hm[[1]][[1]] , merge_legend=TRUE, 
     annotation_legend_list=hm[[1]][[2]][[1]])
dev.off()


# ----------------------------------------------------------------------------------------------------------------- #
#                                              Clustering-goodness confusion plots                                  #
# ----------------------------------------------------------------------------------------------------------------- #

library(SeuratDisk)
library(dplyr)
library(Seurat)

objs <- c(sapply(c("microglia","astrocytes","oligodendroglia","inhibitory"), 
                 \(ct) paste0(ct, "/data/", ct, ".h5Seurat")),
          endo = "vascular.niche/data/vascular.niche.h5Seurat",
          `cux2-`="excitatory/data/cux2-.h5Seurat",
          `cux2+`="excitatory/data/cux2+.h5Seurat")

for (ct in names(objs)) {
  message(ct)
  
  obj      <- LoadH5Seurat(objs[[ct]], assays=list(SCT=c("data")), verbose=F)
  clusters <- Idents(obj)
  
  if(! "SCT_snn" %in% names(obj@graphs))
    obj <- FindNeighbors(obj, reduction="pca", dims = 1:50)
  
  combs <- expand.grid(clusters %>% unique() %>% as.vector(), 
                       clusters %>% unique() %>% as.vector()) %>% t()
  dists <- apply(combs, 2, function(pair_){
    cells1 <- names(clusters)[clusters == pair_[[1]] %>% as.vector()]
    cells2 <- names(clusters)[clusters == pair_[[2]] %>% as.vector()]
    (rowSums(obj@graphs$SCT_snn[cells1, cells2])/rowSums(obj@graphs$SCT_snn[cells1, ])) %>% mean()
  })
  # saveRDS(list(dists = rbind(combs, dists) %>% t() %>% as.data.frame() %>% tidyr::spread(Var1, dists) %>% tibble::column_to_rownames("Var2"), 
             # clusters = clusters), paste0("5. Manuscript code/data/", ct, ".knn.rds"))
  
  saveRDS(list(dists = rbind(combs, dists) %>% t() %>% as.data.frame() %>% tidyr::spread(Var1, dists) %>% tibble::column_to_rownames("Var2"), 
               clusters = clusters), paste0(ct, ".knn.rds"))
  
  rm(obj, clusters, combs, dists)
}


pdf(file.path(panel.path, "2.cluster.goodness.confusions.pdf"), height=embed.height, width = embed.width+.25)
for(ct in names(objs)) {
  df <- readRDS(paste0("5. Manuscript code/data/", ct, ".knn.rds")) 
  mtx <- sapply(df$dists, as.numeric) %>% `rownames<-`(rownames(df$dists))
  
  if(ct %in% c("cux2+", "cux2-"))
    ord <- intersect(atlas$excitatory$state.order, rownames(mtx))
  else
    ord <- atlas[[ct]]$state.order
  
  Heatmap(mtx[ord, ord],
          col=colorRampPalette(c("white","darkgreen","#5F9E5F","#A074B6","darkorchid4"))(21),
          cluster_rows = F, cluster_columns = F, border=T) %>% draw()
  
}
while (!is.null(dev.list()))  dev.off()


# ----------------------------------------------------------------------------------------------------------------- #
#                                            Heatmaps of predicted neuronal subtypes                                #
# ----------------------------------------------------------------------------------------------------------------- #
neurons <- merge(annotations %>% rownames_to_column("cell"), h5read(aggregated.data, "neuronal/allen.annotations")) %>%
  count(grouping.by, state, allen.labels) %>% 
  group_by(grouping.by, state) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup()

neuronal.subtypes <- neurons %>%
  select(allen.labels) %>% unique %>%
  tidyr::separate(col = "allen.labels", sep = " ", into = c("type", "layer", "marker.1", "marker.2"), extra = "merge", remove = FALSE) %>%
  mutate(marker.1 = factor(marker.1, rev(c("LAMP5", "PAX6", "SST", "VIP", "PVALB","RORB", "LINC00507","THEMIS","FEZF2"))),
         across(c(layer, marker.2), ~factor(., sort(unique(.))))) %>%
  column_to_rownames("allen.labels")

ord <- c(
"Exc L2 LINC00507 GLRA3", "Exc L2 LAMP5 KCNG3", "Exc L2 LINC00507 ATP7B", "Exc L2-3 LINC00507 DSG3",
"Exc L3 LAMP5 CARM1P1", "Exc L3 THEMIS ENPEP", "Exc L2-3 RORB RTKN2", "Exc L2-3 RORB PTPN3",
"Exc L2-3 RORB CCDC68", "Exc L3-5 RORB TNNT2", "Exc L3-5 RORB LAMA4", "Exc L3 RORB OTOGL",
"Exc L3-5 RORB LINC01202", "Exc L3-5 RORB LNX2", "Exc L3-5 RORB RPRM", "Exc L5 THEMIS SLC22A18",
"Exc L6 THEMIS LINC00343", "Exc L6 THEMIS SNTG2", "Exc L5-6 THEMIS TNFAIP6", "Exc L6 THEMIS SLN",
"Exc L5 RORB MED8", "Exc L5 THEMIS FGF10", "Exc L5 THEMIS VILL", "Exc L5-6 FEZF2 IFNG-AS1",
"Exc L5 FEZF2 NREP-AS1", "Exc L5 FEZF2 RNF144A-AS1", "Exc L5 FEZF2 PKD2L1", "Exc L5-6 FEZF2 LPO",
"Exc L6 FEZF2 KLK7", "Exc L6 FEZF2 POGK", "Exc L6 FEZF2 FFAR4", "Exc L6 FEZF2 PROKR2",
"Exc L5-6 FEZF2 CFTR", "Exc L6 FEZF2 PDYN", "Exc L5-6 FEZF2 C9orf135-AS1", "Exc L5-6 FEZF2 SH2D1B",
"Exc L5 THEMIS RGPD6", "Exc L5-6 FEZF2 FILIP1L", "Exc L5-6 FEZF2 OR1L8", "Exc L5 THEMIS LINC01116",
"Exc L5-6 THEMIS SMYD1", "Exc L5 FEZF2 CSN1S1", "Exc L3-5 FEZF2 ASGR2", "Exc L3-5 FEZF2 LINC01107",
"Inh L1-6 SST NPY", "Inh L1-6 LAMP5 AARD", "Inh L1 LAMP5 RAB11FIP1", "Inh L1-6 LAMP5 NES", "Inh L5-6 LAMP5 CRABP1",
"Inh L1-6 LAMP5 CA1", "Inh L1 PAX6 MIR101-1", "Inh L3-6 PAX6 LINC01497", "Inh L5-6 SST BEAN1",
"Inh L5-6 SST DNAJC14", "Inh L5-6 SST KLHL1", "Inh L5-6 SST FBN2", "Inh L5-6 SST C4orf26",
"Inh L1-2 SST PRRT4", "Inh L3-5 SST GGTLC3", "Inh L1-2 SST CLIC6", "Inh L2 PVALB FRZB",
"Inh L2-3 SST NMU", "Inh L1-2 SST CCNJL", "Inh L1-3 SST FAM20A", "Inh L3-5 SST CDH3",
"Inh L5 SST RPL35AP11", "Inh L5-6 SST PAWR", "Inh L3-5 SST OR5AH1P", "Inh L5-6 SST PIK3CD",
"Inh L5-6 PVALB SST CRHR2", "Inh L5-6 SST ISX", "Inh L1 SST DEFB108B", "Inh L1 LAMP5 BMP2",
"Inh L1 PVALB SST ASIC4", "Inh L1 PAX6 CHRFAM7A", "Inh L1-6 VIP SLC7A6OS", "Inh L1 LAMP5 PVRL2",
"Inh L1 SST P4HA3", "Inh L1 LAMP5 NMBR", "Inh L2 PAX6 FREM2", "Inh L1-2 VIP HTR3A", "Inh L1-2 VIP WNT4",
"Inh L1-5 VIP PHLDB3", "Inh L1-5 VIP LINC01013", "Inh L3-6 VIP UG0898H09", "Inh L3-6 VIP ZIM2-AS1",
"Inh L1-5 VIP CD27-AS1", "Inh L2-5 VIP SOX11", "Inh L1-2 VIP PTGER3", "Inh L2-5 VIP BSPRY", "Inh L5-6 VIP COL4A3",
"Inh L1-2 VIP SCML4", "Inh L2 VIP SLC6A16", "Inh L1-3 VIP HSPB6", "Inh L1-5 VIP SMOC1", "Inh L1-3 VIP CHRNA2",
"Inh L3-5 VIP HS3ST3A1", "Inh L1-2 VIP EXPH5", "Inh L1-3 VIP FNDC1", "Inh L1 VIP KLHDC8B", "Inh L1-3 VIP CBLN1",
"Inh L3-5 VIP IGDCC3", "Inh L3-5 VIP TAC3", "Inh L1-6 PVALB COL15A1", "Inh L5-6 PVALB ZFPM2-AS1",
"Inh L6 SST TH", "Inh L5-6 PVALB GAPDHP60", "Inh L5-6 PVALB KCNIP2", "Inh L1-2 SST CLIC6",
"Inh L2 PVALB FRZB", "Inh L2-5 PVALB RPH3AL", "Inh L1-2 PVALB CDK20", "Inh L3 PVALB SAMD13",
"Inh L3-5 PVALB ISG20", "Inh L2-5 PVALB HHIPL1", "Inh L5-6 PVALB FAM150B", "Inh L5-6 PVALB MEPE", "Inh L5 PVALB LRIG3")

sets = list(c("Inhibitory Neurons", "inhibitory"), c("Excitatory Neurons", "excitatory"))

pdf(file.path(panel.path, "2.neuronal.confusions.pdf"), height=embed.height*2, width = embed.width*2)
for(set in sets) {
df <- pivot_wider(neurons %>% filter(grouping.by == set[[1]]), id_cols="allen.labels", names_from = "state", values_from = "prop", values_fill = 0) %>% 
  column_to_rownames("allen.labels") %>%
  `[`(intersect(v, rownames(.)), atlas[[set[[2]]]]$state.order) %>%
  as.matrix

Heatmap(df, 
        col=colorRampPalette(c("white","darkgreen","#5F9E5F","#A074B6","darkorchid4"))(21),#colorRampPalette(c("white","#eb8694","#A074B6","darkorchid4"))(21), 
        width = unit(ncol(df)*.35,"cm"),
        border=T,
        column_order = atlas[[set[[2]]]]$state.order, 
        cluster_rows = F, cluster_row_slices = F,
        right_annotation = rowAnnotation(
          `Cortical layer` = anno_text(paste(neuronal.subtypes[rownames(df), "layer"], "  ")),
          `Marker 1` = anno_text(paste(neuronal.subtypes[rownames(df), "marker.1"], "  ")),
          `Marker 2` = anno_text(neuronal.subtypes[rownames(df), "marker.2"])
        ),
        left_annotation = rowAnnotation(
          `Cortical layer` = neuronal.subtypes[rownames(df), "layer"],
          marker1 = neuronal.subtypes[rownames(df), "marker.1"]),
        show_row_names = F) %>% draw()
}
while (!is.null(dev.list()))  dev.off()
rm(df, sets, set, ord, neuronal.subtypes, neurons)
