source("5. Manuscript code/utils.R")


#####################################################################################################################
#                                    Figure 3 - State-Trait Associations & Validations                              #
#####################################################################################################################

# ----------------------------------------------------------------------------------------------------------------- #
#                                  Panel B+C - joint - snuc+bulk Trait Associations                                 #
# ----------------------------------------------------------------------------------------------------------------- #

pdf(file.path(panel.path, "3B-D.pdf"), width=2*3 + 1, height=5)
plot.trait.associations.cross.cohort(
  names(AD.traits), AD.traits, fdr.thr = .01, use_raster=TRUE, raster_quality = 10)
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                     Panel B+C - snuc+bulk Trait Associations                                      #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "3B.pdf"), width=3.5, height=5)
plot.trait.associations(py_to_r(data$uns$trait.analysis$snuc), 
                        params = names(AD.traits),
                        column_title="snRNA-seq",
                        column_labels = AD.traits,
                        column_names_rot = 45,
                        column_names_centered = T,
                        row_names_side = "left",
                        use_raster=T,
                        raster_quality = 10) %>% print()
while (!is.null(dev.list()))  dev.off()


pdf(file.path(panel.path, "3C.pdf"), width=3.5, height=5)
plot.trait.associations(py_to_r(data$uns$trait.analysis$celmod), 
                        params = names(AD.traits),
                        column_title="bulk-predicted",
                        column_labels = AD.traits,
                        column_names_rot = 45,
                        column_names_centered = T,
                        row_names_side = "left",
                        use_raster=T,
                        raster_quality = 10) %>% print()
while (!is.null(dev.list()))  dev.off()


# ----------------------------------------------------------------------------------------------------------------- #
#                                      Panel D - Meta analysis of associations                                      #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "3D.pdf"), width=3.5, height=5)
plot.trait.associations(py_to_r(data$uns$trait.analysis$meta.analysis),
                        params = names(AD.traits),
                        value.by = "z.meta", pval.by = "adj.pval.meta",
                        column_title="meta-analysis",
                        column_labels = AD.traits,
                        column_names_rot = 45,
                        column_names_centered = T,
                        row_names_side = "left") %>% print()
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                             Panels F-J - Causal modeling                                          #
# ----------------------------------------------------------------------------------------------------------------- #
# Panels were designed based on the results in the Other analyses/3.causal modeling.R file


#####################################################################################################################
#                                 Supp Figure 3 - State-Trait Associations & Validations                            #
#####################################################################################################################


# ----------------------------------------------------------------------------------------------------------------- #
#                                         Panel A - CelMod correlations for states                                  #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "s6A.pdf"), width=14, height=3)
print(py_to_r(data$uns$celmod$test.corrs) %>% 
        rownames_to_column("state") %>% 
        mutate(group = factor(data$var[state, "grouping.by"], levels=names(cell.group.color))) %>%
        ggplot(aes(reorder(state, -corr),corr, fill=group, label=sig)) +     
        geom_bar(stat="identity") + 
        geom_text(nudge_y = .01, angle=90, hjust=0, vjust=.65) + 
        facet_wrap(.~group, scales = "free_x", nrow = 1) + 
        scale_fill_manual(values=scales::alpha(cell.group.color, .75)) + 
        scale_y_continuous(expand = expansion(add=c(.1, .2))) + 
        labs(x=NULL, y=NULL) + 
        theme_classic() + 
        theme(legend.position = "none",
              strip.background = element_blank(),
              axis.text.x = element_text(angle=90, hjust = 1, vjust = .5)))
while (!is.null(dev.list()))  dev.off()




# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel B - snRNAseq/bulk tstat comparison                                   #
# ----------------------------------------------------------------------------------------------------------------- #

pdf(file.path(panel.path, "s6B.pdf"), width=embed.width.small, height=embed.height.small*3)
for(traits in list(AD.traits, AD.traits.cat)) {
  df <- py_to_r(data$uns$trait.analysis$meta.analysis) %>% 
    filter(trait %in% names(traits)) %>%
    group_by(trait) %>% 
    mutate(trait = factor(trait, names(traits)),
           label = if_else(adj.pval.meta <= .05, state, NA_character_),
           color = case_when(is.na(label) ~ NA_character_, 
                             tstat.sc > 0 ~ "1",
                             tstat.sc < 0 ~ "-1"))
  
  
  corrs <- do.call(rbind, lapply(unique(df$trait), function(t) 
    data.frame(cor.test(df[df$trait == t,]$tstat.sc, 
                        df[df$trait == t, ]$tstat.b, 
                        method="spearman")[c("estimate","p.value")], row.names = t))) %>%
    mutate(adj.pval = p.adjust(p.value, method = "BH"))
  
  
  plot_grid(plotlist = lapply(levels(df$trait), function(t) {
    .df <- df[df$trait == t,]
    
    range.min <- min(.df[,c("tstat.sc","tstat.b")], na.rm = T)
    range.max <- max(.df[,c("tstat.sc","tstat.b")], na.rm = T)
    
    .label <- paste0("R=", round(corrs[t,]$estimate, 2), 
                     "\nFDR=", scales::scientific(corrs[t,]$adj.pval))
    
    ggplot(.df, aes(tstat.sc, tstat.b, label=label)) + 
      geom_abline(linetype="dashed", color="grey30") + 
      geom_point(size=.75) + 
      ggplot2::annotate(geom="text", x=range.min, y=.95*range.max, hjust=0, vjust=1,label=.label) +
      ggrepel::geom_text_repel(aes(color=color), max.overlaps = 30, show.legend = F, min.segment.length = unit(0, "pt")) + 
      scale_color_manual(values = list("1"=green2purple(3)[3], "-1"=green2purple(3)[1])) + 
      scale_x_continuous(limits = c(range.min, range.max), expand = expansion(mult = .025)) + 
      scale_y_continuous(limits = c(range.min, range.max), expand = expansion(mult = .025)) + 
      labs(x="t-stat snRNA-seq", y="t-stat bulk pred.", title=traits[[t]])
  }),
  ncol=1) %>% print()
  rm(df, corrs)
}
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel C - snuc+bulk Trait Associations                                      #
# ----------------------------------------------------------------------------------------------------------------- #

pdf(file.path(panel.path, "s6C.joint.pdf"), width=2*3 + 1, height=5)
plot.trait.associations.cross.cohort(names(AD.traits.cat), AD.traits.cat, 
                                     fdr.thr = .01, use_raster=TRUE, raster_quality = 10)
while (!is.null(dev.list()))  dev.off()


pdf(file.path(panel.path, "s6C.pdf"), width=3.5, height=5)
plot.trait.associations(py_to_r(data$uns$trait.analysis$snuc), 
                        params = names(AD.traits.cat),
                        column_title="snRNA-seq",
                        column_labels = AD.traits.cat,
                        column_names_rot = 45,
                        column_names_centered = T,
                        row_names_side = "left",
                        use_raster=T,
                        raster_quality = 10) %>% print()

plot.trait.associations(py_to_r(data$uns$trait.analysis$celmod), 
                        params = names(AD.traits.cat),
                        column_title="bulk-predicted",
                        column_labels = AD.traits.cat,
                        column_names_rot = 45,
                        column_names_centered = T,
                        row_names_side = "left",
                        use_raster=T,
                        raster_quality = 10) %>% print()

plot.trait.associations(py_to_r(data$uns$trait.analysis$meta.analysis),
                        params = names(AD.traits.cat),
                        value.by = "z.meta", pval.by = "adj.pval.meta",
                        column_title="meta-analysis",
                        column_labels = AD.traits.cat,
                        column_names_rot = 45,
                        column_names_centered = T,
                        row_names_side = "left") %>% print()
while (!is.null(dev.list()))  dev.off()


# ----------------------------------------------------------------------------------------------------------------- #
#                                             Panels D-I - Causal modeling                                          #
# ----------------------------------------------------------------------------------------------------------------- #
# Panels were designed based on the results in the `3. Other analyses/3.causal modeling.R` file

