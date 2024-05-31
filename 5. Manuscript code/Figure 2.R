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
n.genes = list(inhibitory=4, excitatory=4, endo=4, astrocytes=4, oligodendroglia=3, microglia=3)

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

