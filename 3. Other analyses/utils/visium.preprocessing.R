#' Natacha Comandante-Lou (nc3018 at columbia dot edu)

library(ggplot2)
library(cowplot)
library(dplyr)

preprocess_visium_h5 = function(sample,
                                umi_min = 500,
                                umi_max = Inf,
                                nfeatures_min = 500,
                                mt_filter = TRUE,
                                rp_filter = FALSE,
                                max_mito_pct = 30,
                                plot_qc = TRUE,
                                figure.path = "figures")
{
  
  obj <- Load10X_Spatial(sample, filename = "filtered_feature_bc_matrix.h5",
                         assay = "Spatial", slice = "slice1", filter.matrix = TRUE, 
                         to.upper = FALSE, image = NULL) %>%
    AddMetaData(., meta.data)
  
  # ------------------------------------------- #
  # Filter high mitochondrial RNA spots         #
  # ------------------------------------------- #
  if(mt_filter) {

    obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    med_cutoff <- median(obj$percent.mt) + 2*mad(obj$percent.mt)
    IQR_cutoff <- as.numeric(quantile(obj$percent.mt)[4] + IQR(obj$percent.mt)*1.5)
    
    mitofilter <- min(med_cutoff, IQR_cutoff)
    mitofilter <- max(mitofilter, max_mito_pct)
    
    if (plot_qc) {
      pdf(file.path(figure.path, " qc_mito-pct.pdf"), width=4, height=3.5)
      ggplot(obj@meta.data, aes(percent.my, y = nCount_Spatial)) +
        geom_point(size = 1, color = "lightblue", alpha = 0.7)+
        labs(x="Percent Mitochrondrial RNA",  y="Number of Expressed Genes") +
        theme_classic()
      dev.off()
    }

    obj <- subset(obj, percent.my < mitofilter)
  }
  

  # ------------------------------------------- #
  # Filter spots based on UMI counts            #
  # ------------------------------------------- #
  keep_spots <- which(nCount_Spatial >= umi_min & nCount_Spatial <= umi_max) 

  if (plot_qc) {
    # Plotting UMI count distribution
    pdf(file.path(figure.path, "qc_umi_distr.pdf"), width=4, height=3.5)
    ggplot(obj@meta.data, aes(nCount_Spatial)) +
      geom_histogram(color = "grey")+
      geom_vline(xintercept = umi_min, linetype = "dashed",color = "blue")+
      geom_vline(xintercept = umi_max, linetype = "dashed",color = "blue")+
      theme_classic()
    dev.off()

    # Visualize variance in sequencing depth across data points 
    pdf(file.path(figure.path, "qc_feature_count_spatial.pdf"), width=6, height=6)
    plot_grid(VlnPlot(spt, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend(),
              SpatialFeaturePlot(spt, features = "nFeature_Spatial",pt.size.factor = 2) + theme(legend.position = "right"))
    dev.off()
  }

  obj <- obj[, keep_spots]

  
  # ------------------------------------------- #
  # Filter spots based on number of features    #
  # ------------------------------------------- #
  obj <- subst(obj, nFeature_Spatial > nfeatures_min) 

  return(obj)
}

