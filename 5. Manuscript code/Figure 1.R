source("5. Manuscript code/utils.R")
source("1. Library preprocessing/utils/ROSMAP.metadata.R")

# Loading list of snRNA-seq participants, while correcting mistakes in tracker file
participants <- openxlsx::read.xlsx("5. Manuscript code/data/ROSMAP 10X Processing Tracker.xlsx", sheet = "Processed Batches") %>%
  filter(!Batch %in% c("B1", "B2", "B3") & (!StudyCode %in% c("MAP83034844", "MAP74718818"))) %>%
  mutate(StudyCode = as.character(as.numeric(gsub("ROS|MAP", "", StudyCode)))) %>% 
  select(StudyCode, Batch) %>%
  rbind(data.frame(StudyCode="44299049", Batch=55))


# The following loads ROSMAP participants' demographic- and endophenotypic characterization 
# Parts of these data are available in supplementary table 1
cohort <- load.metadata() %>% `[`(unique(participants$StudyCode), )

#####################################################################################################################
#                                          Figure 1 - Cohort & Atlas Overview                                       #
#####################################################################################################################

# ----------------------------------------------------------------------------------------------------------------- #
#                                 Panel B - Categorical Pathology & Cognitive Decline                               #
# ----------------------------------------------------------------------------------------------------------------- #

p1 <- ggplot(cohort, aes(sex, age_death)) + 
  geom_jitter(alpha=.4, width = .15) +
  geom_boxplot(alpha=.5, fill="darkorchid4", outlier.size = 0) + 
  stat_summary(geom = "text", fun.data = function(x) list(y=65, label=paste0("n=", length(x))))  +
  labs(x=NULL, y="Age of death", title=NULL) + 
  scale_y_continuous(breaks = c(70, 90, 110)) + 
  theme_classic()

pdf(file.path(panel.path, "1B.pdf"), width=embed.width, height=embed.height*2)
plot_grid(p1, plot_grid(plotlist = lapply(c("cogdx.grouped.txt","braak.grouped.txt","cerad.txt"), function(t)
  ggplot(cohort, aes(!!sym(t), fill=!!sym(t))) + 
    geom_bar(color="black") + 
    theme_classic() + 
    scale_fill_manual(values = colorRampPalette(c("white","darkorchid4"))(1+length(unique(cohort[,t])))[-1]) + 
    scale_y_continuous(expand = c(0,0), 
                       limits = c(0, 170), 
                       breaks = c(50,100, 150)) + 
    labs(x=NULL, y=t, title=NULL)), 
  ncol=1),
  rel_heights = c(1,3), ncol = 1)
while (!is.null(dev.list()))  dev.off()
rm(p1)



# ----------------------------------------------------------------------------------------------------------------- #
#                                  Panel C - quantative Pathology & Cognitive Decline                               #
# ----------------------------------------------------------------------------------------------------------------- #
df <- cohort[,c("sqrt.amyloid","sqrt.tangles","cogng_demog_slope","sex")] %>% 
  filter(!is.na(sqrt.amyloid) & !is.na(sqrt.tangles)) %>%
  rename(x="sqrt.amyloid", y="sqrt.tangles") 

grid <- expand.grid(x=seq(-.05, 1.05*max(df$x, na.rm = T), length.out=200),
                    y=seq(-.05, 1.05*max(df$y, na.rm = T), length.out=200))


dist <- rdist::cdist(grid, df[!is.na(df$cogng_demog_slope),1:2])^2
dens <- apply(dist, 1, function(x) sqrt(mean(sort(x, partial=10)[1:10])))

C <- as.matrix(df[!is.na(df$cogng_demog_slope),"cogng_demog_slope"])
grid$cogng_demog_slope <- exp(-dist/(.75*dens)) %*% C
  

pdf(file.path(panel.path, "1C.pdf"), width=embed.width, height=embed.height+1)
ggplot(df %>% arrange(!is.na(cogng_demog_slope),cogng_demog_slope), aes(x,y, color=cogng_demog_slope, shape=sex)) + 
  ggrastr::geom_tile_rast(aes(x,y,fill=cogng_demog_slope), grid, inherit.aes = F, alpha=.8, color=NA, raster.dpi = 10) + 
  geom_point(size=3.5, color="black", alpha=.4) +
  geom_point(size=2.7) +
  scale_fill_gradientn(colors = rev(green2purple.less.white(5))) +
  scale_color_gradientn(colors = rev(green2purple.less.white(5))) + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="amyloid",y="tangles", shae=NULL, color="Cognitive decline rate") + 
  theme_classic() + 
  theme(legend.position = "bottom")
while (!is.null(dev.list()))  dev.off()
rm(df, grid, dist, C, dens)



# ----------------------------------------------------------------------------------------------------------------- #
#                                                Panel D - UMAP of all nuclei                                       #
# ----------------------------------------------------------------------------------------------------------------- #
# 2D UMAP embedding of final atlas with cell-type annotations
embedding <- h5read(aggregated.data, "umap") %>% tibble::column_to_rownames("rowname") %>% `colnames<-`(c("x","y")) 
shared    <- intersect(rownames(embedding), rownames(annotations))
embedding <- data.frame(embedding[shared,], annotations[shared, ])


# Plot UMAP embedding of atlas
pdf(file.path(panel.path, "1D.pdf"), width=8, height=8)
ggplot(embedding, aes(x,y, color=state)) +
  ggrastr::geom_point_rast(size=.05, raster.dpi = 600) + 
  scale_color_manual(values = joint.state.colors, na.value = "lightgrey") +
  theme_classic() +
  no.labs + 
  no.axes + 
  theme_embedding + 
  theme(legend.position = "none")
while (!is.null(dev.list()))  dev.off()
rm(embedding, shared)


# ----------------------------------------------------------------------------------------------------------------- #
#                                      Panel E - Number of nuclei in cell groups                                    #
# ----------------------------------------------------------------------------------------------------------------- #
df <- annotations %>% filter(projid != "NA" & grouping.by != "Immune") %>% 
  count(grouping.by, projid) %>% 
  mutate(
    grouping.by = recode(grouping.by, "Excitatory Neurons"="Exc. Neurons", "Inhibitory Neurons"="Inh. Neurons", "Oligodendrocytes"="Oligodend."),
    grouping.by = factor(grouping.by, levels = levels(grouping.by)))

cols <- cell.group.color %>% `names<-`(recode(names(.), "Excitatory Neurons"="Exc. Neurons", "Inhibitory Neurons"="Inh. Neurons", "Oligodendrocytes"="Oligodend."))

pdf(file.path(panel.path, "1E.pdf"), width=embed.width*2/3, height=embed.height)
ggplot(df, aes(grouping.by, n, fill=grouping.by, color=grouping.by)) + 
  geom_boxplot(alpha=.3, outlier.shape = NA) +
  geom_point(alpha=.5, position = position_jitterdodge(jitter.width = 2.5, dodge.width = .7, jitter.height = 0)) + 
  labs(x=NULL, y=NULL, title="Nuclei per donor") + 
  scale_y_sqrt(expand=c(0,0), breaks=c(50, 100,250, 500, 800, 1500, 3000)) + 
  scale_color_manual(values = cols) + 
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
while (!is.null(dev.list()))  dev.off()
rm(df, cols)


# ----------------------------------------------------------------------------------------------------------------- #
#                                            Panel F - State Prevalence Heatmap                                     #
# ----------------------------------------------------------------------------------------------------------------- #
# Cluster donors by distance in cellular landscape
df <- data.frame(data$X)
tree <- hclust(rdist::rdist(df, "euclidean")^2)

# Order states first by grouping and then by their numeric order (i.e Exc.2 < Exc.10)
ord <- states %>% rownames_to_column("state") %>% 
  arrange(grouping.by, suppressWarnings(as.numeric(gsub("^.*\\.","", state))), gsub("^.*\\.","", state)) %>% 
  pull(state)

# Flatten data for plotting
df <- df %>% rownames_to_column("projid") %>% 
  melt(id.vars = "projid", variable.name="state") %>%
  mutate(state = factor(gsub(" ", ".", state), levels = gsub(" ", ".", rev(ord))),
         projid = factor(projid, levels = rownames(df)[tree$order], ordered = T))


pdf(file.path(panel.path, "1F.pdf"), width=embed.width*2/3, height=embed.height)
bars <- ggplot(df, aes(projid,value,fill=state)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=joint.state.colors) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) + 
  coord_flip() + 
  labs(x="Donors", y="Proportion of cell population") + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())

bars + ggtree::ggtree(tree, size=.25) + aplot::ylim2(bars) + patchwork::plot_layout(widths = c(1,.2))
while (!is.null(dev.list()))  dev.off()
rm(df, tree, ord, bars)




#####################################################################################################################
#                                        Supp Figure 1 - Library pre-processing                                     #
#####################################################################################################################
cols <- list(
  "Exc"  = "midnightblue",
  "Inh"  = "firebrick4",
  "Olig"  = "olivedrab4",
  "Astr" = "darkorchid4",
  "Micr" = "chocolate3",
  "OPC"  = "springgreen4",
  "Endo" = "darkgoldenrod4",
  "Peri" = "cyan4")

# ----------------------------------------------------------------------------------------------------------------- #
#                                         Panel A - batch trait distributions                                       #
# ----------------------------------------------------------------------------------------------------------------- #
df <- cohort %>% select(sex, cogdx, cerad.txt, braak.grouped.txt) %>%
  merge(participants, by.x="row.names", by.y="StudyCode") %>%
  mutate(Batch = as.numeric(gsub("B", "", Batch))) %>%
  mutate(Batch = factor(Batch, levels = unique(sort(Batch))))
  
pdf(file.path(panel.path, "s1A.pdf"), width=embed.width*1.5, height=embed.height*1.5)
lapply(c("sex","cogdx","braak.grouped.txt","cerad.txt"), function(p) {
  ggplot(df, aes_string("Batch", fill=p)) + 
    geom_bar(stat = "count") + 
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), breaks = c(2,4,6,8)) + 
    scale_fill_manual(values = colorRampPalette(c("white","darkorchid4"))(1+length(unique(df[,p])))[-1]) + 
    
    labs(x=NULL, y=NULL) +
    theme_classic() + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.border = element_rect(color="black", fill=NA, size=1))
}) %>% 
  base::Reduce(`+`, .) + patchwork::plot_layout(ncol=1)
while (!is.null(dev.list()))  dev.off()
rm(df)


# ----------------------------------------------------------------------------------------------------------------- #
#                                  Panel C - classifier predictions over example libraרies                          #
# ----------------------------------------------------------------------------------------------------------------- #

pdf(file.path(path, "s1C.pdf"), width=embed.width*2, height=embed.height*2)
lapply(c("191213-B7-B", "200720-B36-A"), function(n) {
  o <- readRDS(paste0("1. Library preprocessing/data/low.qual.thr.libs/", n, ".seurat.rds") )
  plot_grid(
    DimPlot(o, group.by = "cell.type", pt.size = 1.75, ncol = 1, cols=cols, raster = T, raster.dpi = c(1024,1024)) + 
      theme_embedding + labs(x=NULL, y=NULL, title=NULL) + theme(legend.position = "none"),
    FeaturePlot(o, features = "cell.type.entropy", pt.size = 1.75, ncol = 1, cols=viridis::inferno(11), 
                order = T, raster = T, raster.dpi = c(1024,1024)) + 
      theme_embedding + labs(x=NULL, y=NULL, title=NULL) + theme(legend.position = c(3,3)))
}) %>% plot_grid(plotlist = ., nrow=2) %>% print()
while (!is.null(dev.list()))  dev.off()


# ----------------------------------------------------------------------------------------------------------------- #
#                                    Panel D - manual curation of low quality thresholds                            #
# ----------------------------------------------------------------------------------------------------------------- #
o <- readRDS("1. Library preprocessing/data/low.qual.thr.libs/191213-B7-B.seurat.rds")
cols2 <- colorRampPalett
cols2 <- cols(length(unique(o$SCT_snn_res.0.2)))

pdf(file.path(path, "s1D.pdf"), width=embed.width, height=embed.height*2)
plot_grid(
  DimPlot(o, group.by = "SCT_snn_res.0.2", pt.size = 1.75, ncol = 1, cols = cols2,
          raster = T, raster.dpi = c(1024,1024), label=T) + 
    theme_embedding + labs(x=NULL, y=NULL, title=NULL) + theme(legend.position = "none"),
  
  plot_grid(VlnPlot(o, features=c("nCount_RNA"), pt.size = 0, ncol = 1, log = T) + 
              labs(x=NULL, y=NULL, title=NULL) + 
              scale_fill_manual(values=cols2) + 
              theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = .5)),
            VlnPlot(o, features=c("nFeature_RNA"), pt.size = 0, ncol = 1, log = T) + 
              scale_fill_manual(values=cols2) + 
              labs(x=NULL, y=NULL, title=NULL) + 
              theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = .5)),
            VlnPlot(o, features=c("cell.type.entropy"), pt.size = 0, ncol = 1) + 
              scale_fill_manual(values=cols2) + 
              scale_y_sqrt() + labs(x=NULL, y=NULL, title=NULL) + 
              theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = .5)),
            ncol=1, align = "hv"),
  
  ncol=1, rel_heights = c(1,1.3))
while (!is.null(dev.list()))  dev.off()
rm(o, cols2)


# ----------------------------------------------------------------------------------------------------------------- #
#                              Panels E-F - nUMAI & nFeature distributions in curated libraries                     #
# ----------------------------------------------------------------------------------------------------------------- #
# The code for these panels can be found under Library preprocessing/low.qual.thresholds.cell.level.R


# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel G - Example library Demux doublets UMAP                               #
# ----------------------------------------------------------------------------------------------------------------- #
o <- readRDS("5. Manuscript code/data/library.for.fig.191213-B7-B.seurat.rds")

pdf(file.path(path, "s1G.pdf"), width=embed.width, height=embed.height)
FetchData(o, c("UMAP_1","UMAP_2","is.doublet.demux")) %>% 
  `colnames<-`(c("x", "y", "demux")) %>% 
  #arrange(demux) %>%
  ggplot(aes(x,y,color=demux)) + 
  geom_point(size=.75) + 
  scale_color_manual(values = c("#468E46","#8D59A8")) + 
  labs(x=NULL, y=NULL) + 
  theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = c(3,3))
while (!is.null(dev.list()))  dev.off()
rm(o)


# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel H - Example library DoubletFinder scores                              #
# ----------------------------------------------------------------------------------------------------------------- #
o <- readRDS("5. Manuscript code/data/library.for.fig.191213-B7-B.seurat.rds")

pdf(file.path(path, "s1H.pdf"), width=embed.width, height=embed.height)
FetchData(o, c("UMAP_1","UMAP_2","doublet.score")) %>% 
  `colnames<-`(c("x", "y", "score")) %>% 
  arrange(score) %>%
  ggplot(aes(x,y,color=score)) + 
  ggrastr::geom_point_rast(size=.75, raster.dpi = 600) + 
  scale_color_viridis_c(option = "inferno") + 
  labs(x=NULL, y=NULL) + 
  theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = c(3,3))
while (!is.null(dev.list()))  dev.off()
rm(o)


# ----------------------------------------------------------------------------------------------------------------- #
#                                        Panel I - MCC distributions across libraries                               #
# ----------------------------------------------------------------------------------------------------------------- #
libs <- list.files("1. Library preprocessing/data/snRNA-seq libraries/", "*.seurat.rds", full.names = TRUE, recursive = TRUE)
pb <- progress_bar$new(format="Loading libraries :current/:total [:bar] :percent in :elapsed. ETA :eta",
                       total = length(libs), clear=FALSE, width=100, force = TRUE)

library(ROCit)
mcc <- lapply(libs, function(p) {
  obj <- readRDS(p)
  
  doublets.df          <- data.frame(class=obj$demux.droplet.type, score=obj$doublet.score, row.names = rownames(obj@meta.data))[obj$demux.droplet.type %in% c("SNG", "DBL"),]
  doublets.measure     <- measureit(score=doublets.df$score, class=doublets.df$class, negref = "SNG", measure = c("ACC", "FSCR", "FPR", "TPR"))
  doublets.measure$MCC <- (doublets.measure$TP*doublets.measure$TN - doublets.measure$FP * doublets.measure$FN) / 
    sqrt((doublets.measure$TP+doublets.measure$FP)*(doublets.measure$TP+doublets.measure$FN)*(doublets.measure$TN+doublets.measure$FP)*(doublets.measure$TP+doublets.measure$FN))
  
  pb$tick()
  return(doublets.measure)
})

pdf(file.path(path, "s1I.pdf"), width=embed.width, height=embed.height)
plot_grid(
  lapply(seq_along(mcc), function(i) data.frame(i=i, thr=mcc[[i]]$Cutoff, MCC=mcc[[i]]$MCC)) %>% 
    do.call(rbind, .) %>%
    ggplot(aes(x=thr, y=MCC, group=i)) + 
    geom_line() + 
    scale_x_continuous(expand = c(0,0)) + 
    labs(x=NULL, y=NULL) + 
    theme_classic()+ 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()),
  data.frame(v=unlist(lapply(mcc, function(x) x$Cutoff[which.max(x$MCC)]))) %>%
    ggplot(aes(v)) + 
    geom_boxplot(outlier.size = 0, fill="grey80") +
    geom_jitter(aes(y=0), height = .1, width=0) +
    scale_x_continuous(expand = c(0,0)) + 
    labs(x=NULL, y=NULL) + 
    theme_classic() + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()),
  ncol=1, 
  rel_heights = c(3,1), align = "v")
while (!is.null(dev.list()))  dev.off()
rm(libs, pb, mcc)


# ----------------------------------------------------------------------------------------------------------------- #
#                              Panel J - Example library after removing low-quality + doublets                      #
# ----------------------------------------------------------------------------------------------------------------- #
# Note that in the full analysis of all libraries doulbets and cells with high levels of mitochondrial genes 
# expression were only removed once cell-type specific objects were created. The following is to illustrate
# this process in the scope of the specific exaple library

o <- readRDS("5. Manuscript code/data/library.for.fig.191213-B7-B.seurat.rds")
o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
!(o$is.doublet | o$percent.mt >= 10)
o <- subset(o, cells = colnames(o)[!(o$is.doublet | o$percent.mt >= 10)])

# Rerunning Seurat pipeline after removing cells - parameters are as used when initially analysing the library
o <- SCTransform(o, variable.features.n = 2000, conserve.memory = T)
o <- RunPCA(o, features=VariableFeatures(o), npcs = 30, verbose = F)
o <- FindNeighbors(o, reduction="pca", dims = 1:30, verbose = F)
o <- FindClusters(o, resolution = .2, verbose = F)
o <- RunUMAP(o,  dims=1:30, n.components = 2)

pdf(file.path(path, "s1J.pdf"), width=embed.width, height=embed.height)
FetchData(o, c("UMAP_1","UMAP_2","cell.type")) %>% 
  `colnames<-`(c("x", "y", "cell.type")) %>% 
  ggplot(aes(x,y,color=cell.type)) + 
  geom_point(size=.5) + 
  scale_color_manual(values = cols, name="Cell Type") +
  labs(x=NULL, y=NULL) + 
  theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = c(3,3))
while (!is.null(dev.list()))  dev.off()
rm(o)


# ----------------------------------------------------------------------------------------------------------------- #
#                                        Panel K - Distribution of #nuclei for participants                         #
# ----------------------------------------------------------------------------------------------------------------- #
df <- annotations %>% filter(projid != "NA") %>% count(projid)

pdf(file.path(path, "s1K.pdf"), width=embed.width, height=embed.height)
ggplot(df %>% rowwise() %>% mutate(n = min(n, 7000)), aes(x = n)) + 
  geom_histogram(alpha=.75, color="black",fill="#8D59A8", bins = 45) +  
  geom_vline(xintercept = 868) + 
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x="# Nuclei", y="# Participants")
while (!is.null(dev.list()))  dev.off()
rm(df)


# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel I - Number of #UMIs in cell groups                                    #
# ----------------------------------------------------------------------------------------------------------------- #
df <- annotations %>% filter(projid != "NA" & grouping.by != "Immune") %>% 
  group_by(grouping.by, projid) %>% 
  summarise(umi=mean(nCount_RNA), .groups = "drop") %>% 
  mutate(
    grouping.by = recode(grouping.by, "Excitatory Neurons"="Exc. Neurons", "Inhibitory Neurons"="Inh. Neurons", "Oligodendrocytes"="Oligodend."),
    grouping.by = factor(grouping.by, levels = levels(grouping.by)))

cols <- cell.group.color %>% `names<-`(recode(names(.), "Excitatory Neurons"="Exc. Neurons", "Inhibitory Neurons"="Inh. Neurons", "Oligodendrocytes"="Oligodend."))

pdf(file.path(path, "s1I.pdf"), width=embed.width*2/3, height=embed.height)
ggplot(df, aes(grouping.by, umi, fill=grouping.by, color=grouping.by)) + 
  geom_boxplot(alpha=.3, outlier.shape = NA) +
  geom_point(alpha=.5, position = position_jitterdodge(jitter.width = 2.5, dodge.width = .7, jitter.height = 0)) + 
  labs(x=NULL, y=NULL, title="Average UMIs per donor") + 
  scale_y_sqrt(expand=c(0,0), breaks=c(1000, 5000, 10000, 20000)) + 
  scale_color_manual(values = cols) + 
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
while (!is.null(dev.list()))  dev.off()
rm(df, cols)



# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel M - Donor-nuclei stacked bar plot                                 #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(path, "s1M.pdf"), width=embed.width*2/3, height=embed.height)
print(ggplot(annotations %>% filter(projid != "NA" & grouping.by != "Immune") %>% count(projid, grouping.by), 
             aes(reorder(projid, -n, sum), n, fill=reorder(grouping.by, -n, sum))) + 
        geom_bar(stat = "identity") + 
        scale_y_continuous(expand=c(0,0)) + 
        scale_fill_manual(values = cell.group.color) + 
        coord_flip() + 
        labs(x=NULL, y=NULL, fill="Cell Type Group", title="# Nuclei per donor") + 
        theme_classic() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = c(.75,.35)))
while (!is.null(dev.list()))  dev.off()

