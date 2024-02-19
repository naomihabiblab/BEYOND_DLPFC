library(dplyr)
library(Seurat)

####################################################################################################################
##                            #  Fit SVM classifier for calling low quality clusters   #                          ##
####################################################################################################################
# The following libraries underwent the basic Seurat pipeline analysis, and manual calling for low quality clusters, 
# as seen in `low.qul.thresholds.cell.type.R`


# -------------------------------------------------- #
# Loading Datasets                                   #
# -------------------------------------------------- #

low.qual.clusters <- list(`191213-B7-B` =c(3,4),        `200225-B10-A`=c(0,1,7,14,15,18,20,21),  `200302-B13-B`=c(0,1,12,13,14,16,20),
                          `200309-B17-B`=c(0,1,11,14),  `200313-B22-B`=c(0,1,13,14,17,18,19),    `200701-B28-B`=c(0,1,13,14),
                          `200720-B36-A`=c(0,1,14),     `200806-B44-B`=c(0,1,12),                `201007-B57-B`=c(0,8,9),     
                          `201021-B60-B`=c(0,1))

libs <- list.files("1. Library preprocessing/data/low.qual.thr.libs", ".seurat.rds", full.names = TRUE, recursive = TRUE)
thresholds <- readRDS("1. Library preprocessing/data/low.quality.thresholds.rds")

df <- do.call(rbind, lapply(libs, function(p) { 
  o <- readRDS(p)
  
  return(
    o@meta.data %>%
      mutate(clusters  = SCT_snn_res.0.2,
             high.qual = !(SCT_snn_res.0.2 %in% low.qual.clusters[[o@project.name]]),
             
             over_nCount_RNA   =  nCount_RNA      >= thresholds[cell.type, "nCount_RNA"],
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
      )))
}))


# -------------------------------------------------- #
# Calculate Batch-Cluster-wise stats                 #
# -------------------------------------------------- #
features = c("% Under Both Thresholds", "Avg. Pred. Entropy", "% Demux Doublets")
batch.df <- df %>% 
  group_by(batch, clusters, high.qual) %>% 
  mutate(demux.droplet.type = demux.droplet.type == "DBL") %>% 
  summarise_at(c("cell.type.entropy", "demux.droplet.type", "over.both"), mean) %>% 
  ungroup() %>%
  mutate(
    x1 = 100 - round(100*over.both, 2), 
    x2 = round(cell.type.entropy, 3),
    x3 = round(100*demux.droplet.type, 2),
    y  = as.factor(high.qual),
    cluster = as.character(clusters)) %>%
  select(x1, x2, x3, y, batch, cluster) %>%
  as.data.frame


  
  
# -------------------------------------------------- #
# Fit Soft-SVM model to separate classes             #
# -------------------------------------------------- #  
library(e1071)
library(caTools)  
library(pracma)

# Cross-Validation for selecting cost and gamma parameters
svm.tune <- tune(svm, y~x1+x2+x3, data=batch.df, kernel="linear", scale=T, shrinking=T,
                 ranges=list(cost=seq(10^-2, 10^2, length.out = 6), 
                             gamma=seq(2^-2, 2^2, length.out = 6),
                             class.weights=lapply(seq(1, 10, by = 1), function(a) c(`TRUE`=a, `FALSE`=1))),
                 tunecontrol = tune.control(cross=10))

fit <- svm(y~x1+x2+x3, batch.df, method="C-classification", kernel="linear", scale = T, shrinking=T, class.weights=svm.tune$best.parameters$class.weights[[1]], 
           cost=svm.tune$best.parameters$cost, gamma=svm.tune$best.parameters$gamma)
p.coefs = t(fit$coefs) %*% fit$SV


saveRDS(fit, "1. Library preprocessing/data/low.quality.clusters.classifier.rds")


# -------------------------------------------------- #
# Plot data and SVM hyperplane                       #
# -------------------------------------------------- #  
fig.df <- cbind(batch.df, as.data.frame(setNames(lapply(1:3, function(i) 
  (batch.df[,i] - fit$x.scale$`scaled:center`[i]) / fit$x.scale$`scaled:scale`[i] ), c("x1.scaled", "x2.scaled", "x3.scaled"))))


x <- seq(min(fig.df$x1.scaled) - 5, max(fig.df$x1.scaled) + 5, length.out = 2)
y <- seq(min(fig.df$x2.scaled) - 5, max(fig.df$x2.scaled) + 5, length.out = 2)
z <- t(outer(x, y, FUN=function(x,y) -(p.coefs[1]*x + p.coefs[2]*y - fit$rho) / p.coefs[3] ))

library(plotly)
p <- plot_ly() %>% 
  add_markers(
    type="scatter3d", mode="markers", data=fig.df,
    x=~x1.scaled, y=~x2.scaled, z=~x3.scaled,
    color=~y, symbol = ~y, text=~cluster,
    symbols=c("circle", "square"), colors=c("navy","red3"), size=2,
    hovertemplate = ~paste0("<b>Batch: ",batch,"</b>",
                            "<br>Cluster #",cluster,
                            "<br>High Quality: ",y,
                            "<br>Under Thresholds: ", x1,"%",
                            "<br>Avg. Entropy: ", x2 ,
                            "<br>Demux Doublet: ", x3, "%",
                            "<extra></extra>")) %>%
  add_surface(z=z, x=x, y=y, opacity=.4, showscale=F) %>%
  layout(showlegend=F, hoverlabel=list(align="left"), 
         title="Separating Hyperplane of Low/High Quality Clusters",
         margin=list(t=60,l=0,r=0,b=0, pad=0),
         scene = list(xaxis=list(title=features[1], titlefont=list(size=10), showticklabels=F, autorange=F, range=c(min(fig.df$x1.scaled), max(fig.df$x1.scaled))),
                      yaxis=list(title=features[2], titlefont=list(size=10), showticklabels=F, autorange=F, range=c(min(fig.df$x2.scaled), max(fig.df$x2.scaled))),
                      zaxis=list(title=features[3], titlefont=list(size=10), showticklabels=F, autorange=F, range=c(min(fig.df$x3.scaled), max(fig.df$x3.scaled))),
                      camera = list(eye = list(x = -1.5, y = 2, z = .75)),
                      aspectmode="cube"))

htmlwidgets::saveWidget(as_widget(p), "1. Library preprocessing/graphs/low.quality.cluster.thresholds.html")


