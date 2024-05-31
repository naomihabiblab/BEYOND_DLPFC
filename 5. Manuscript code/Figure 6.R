source("5. Manuscript code/utils.R")


#####################################################################################################################
#                                           Figure 6 - Communities Structure                                        #
#####################################################################################################################

# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel A - Dynamics - states of interest                                     #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "6A.pdf"), width=embed.width, height=embed.height.small)
for(states in state.dynamics.sets) {
  print(plot.dynamics(names(states), cols = states, overlap.pseudotime=.1, label=TRUE) + 
    theme(strip.text = element_blank(), 
          legend.position="none", 
          axis.line = element_line(), 
          axis.text.y = element_text()) + 
    labs(x=NULL, y=NULL, title=NULL))
}
while (!is.null(dev.list()))  dev.off()
rm(states)



# ----------------------------------------------------------------------------------------------------------------- #
#                                      Panel B - Dynamics & adjacency matrices                                      #
# ----------------------------------------------------------------------------------------------------------------- #
library(dendextend)
communities   <- data$uns$communities
membership    <- split(data$var_names, data$var$community)
dynamics      <- communities$similarities$dynamics %>% `dimnames<-`(list(data$var_names, communities$dynamics.colnames))
corr          <- communities$similarities$correlation %>% `dimnames<-`(list(data$var_names, data$var_names))
dyn.adjacency <- communities$similarities$dynamics.adjacency %>% `dimnames<-`(list(data$var_names, data$var_names))


# States to annotate
mark.states <- c("Oli.2","Oli.3","Ast.1","Mic.2","OPC.1","Ast.3","Oli.6","Mic.12","Mic.13","Ast.10","Oli.7","Ast.5","OPC.3","Mic.6","Mic.7")

# Specify dynamics column-order: zero pseudotime in the middle
dynamics <- dynamics[,c(rev(colnames(dynamics)[grepl("ABA", colnames(dynamics))]), 
                        colnames(dynamics)[grepl("prAD", colnames(dynamics))])]

# Hierarchically split community's membership, creating dynamics heatmap within each community/sub-community
hm.comms <- data$var %>% split(., as.character(.$community)) %>% 
  lapply(., function(inner) {
    
    if(any(!is.na(inner$sub.community))) 
      lst <- split(rownames(inner), as.character(inner$sub.community))
    else
      lst <- list(rownames(inner))
    
    lapply(lst, function(states) {
      mtx  <- dyn.adjacency[states,states] + corr[states,states]
      dend <- dendsort::dendsort(hclust(dist(mtx))) %>% as.dendrogram() %>% set("labels_to_character")
      
      prepare(Heatmap(dynamics[states,],
                      column_split =  factor(gsub("_.*", "", colnames(dynamics)), levels=c("ABA","prAD")),
                      cluster_rows = dend, 
                      cluster_row_slices = F,
                      cluster_columns = F,
                      show_column_names = F,
                      show_row_names = T,
                      show_row_dend = T,
                      row_dend_width = unit(10,"pt"),
                      col = circlize::colorRamp2(seq(-4,4,length.out=21), 
                                                 colorRampPalette(c("darkorchid4","white","#E65100"))(21)),
                      column_title = NULL,
                      row_names_side = "left",
                      column_gap = unit(1, "pt"),
                      right_annotation = rowAnnotation(states = anno_mark(which(rownames(mtx) %in% mark.states), intersect(rownames(mtx), mark.states)))))
      
    }) %>% unlist(., recursive=FALSE)
  })  %>% unlist(., recursive=FALSE)

hm.comms <- hm.comms %>% `[`(c("C3", "C1.C1.1", "C1.C1.2","C2.C2.2","C2.C2.3","C2.C2.1"))


# Create annotations of traits and pseudotime
df <- data$uns$trajectories$palantir$dynamics$pred.vals %>% py_to_r(.) %>% 
  filter(feature %in% names(AD.traits)) %>% 
  mutate(col=paste(trajectory, x, sep = "_"))

traits <- df %>% reshape2::dcast(col~feature, value.var="fit") %>% column_to_rownames("col")
hm.traits <- lapply(names(AD.traits), function(trait) {
  mtx <- t(traits[colnames(dynamics),trait])
  Heatmap(mtx, 
          col = circlize::colorRamp2(seq(min(mtx), max(mtx), length.out=21),
                                     colorRampPalette(c("white", AD.traits.colors[[trait]]))(21)),
          name = AD.traits[[trait]],
          cluster_columns = F,
          height = unit(.25,"cm"),
          right_annotation = rowAnnotation(n=anno_mark(1, trait)))
}) %>% `names<-`(names(AD.traits))


pseudotime <- df %>% dplyr::select(col, x) %>% unique() %>% dcast(1~col, value.var = "x") %>% column_to_rownames("1") %>% `[`(,colnames(dynamics)) %>% as.matrix() 
pseudotime <- Heatmap(pseudotime,
                      col = viridis::turbo(21),
                      name = "pseudotime",
                      show_row_names = F,
                      show_column_names = F,
                      cluster_columns = F,
                      height = unit(.25,"cm"),
                      right_annotation = rowAnnotation(n=anno_mark(1, "pseudotime")))

# Save figure of all heatmaps
pdf(file.path(panel.path, "6B.1.pdf"), width=embed.width*1.5, height=embed.height*2)
draw(modifyList(modifyList(hm.comms, hm.traits), list(ps=pseudotime)) 
     %>% base::Reduce("%v%", .), ht_gap = unit(1, "pt"))
while (!is.null(dev.list()))  dev.off()

# Retrieve state order to be used in composite adjacency matrix - Fig 5c
state.order <- do.call(rbind, lapply(hm.comms, function(hm) data.frame(comm = hm@name, state=rownames(hm@matrix)[row_order(hm)])))

rm(dynamics, hm.comms, df, traits, hm.traits, pseudotime)


# Plot dynamics adjacency and state-state correlations in lower/upper triangles
dyn.adjacency <- dyn.adjacency[state.order$state, state.order$state]
corr <- corr[state.order$state, state.order$state]

mtx                 <- dyn.adjacency
mtx[upper.tri(mtx)] <- 10 + corr[upper.tri(corr)]

breaks   <- list(seq(0,1,length.out=21), 10 + seq(-1,1,length.out=21))
cols     <- list(colorRampPalette(c("white","salmon","red","firebrick4"))(21), green2purple.less.white(21))
cols.fun <- circlize::colorRamp2(unlist(breaks), unlist(cols))


marks <- sapply(mark.states, function(s) which(rownames(mtx) == s))

hm <- Heatmap(
  mtx,
  row_split = state.order$comm,
  column_split = state.order$comm,
  cluster_rows = F, 
  cluster_columns = F,
  col = cols.fun,
  cell_fun = function(j, i, x, y, w, h, col) {
    if(i==j) grid.rect(x,y,w,h, gp = gpar(fill="black", col=NA))
  },
  column_title = NULL,
  row_title = NULL,
  show_heatmap_legend = F,
  show_row_names = F, 
  show_column_names = F,
  row_gap = unit(2, "pt"),
  column_gap = unit(2, "pt"),
  right_annotation = rowAnnotation(states = anno_mark(marks, names(marks)))
)

pdf(file.path(panel.path, "6B.2.pdf"), width=embed.width*2, height=embed.height*2)
draw(hm)
while (!is.null(dev.list()))  dev.off()

rm(dyn.adjacency, corr, mtx, breaks, cols, cols.fun, marks, mark.states, hm, state.order, communities, membership)



# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel C - Communities dynamics                                          #
# ----------------------------------------------------------------------------------------------------------------- #
comm.cols <- list(C1="springgreen4", C2.2="darkorchid4", C2.3="firebrick4", C3="firebrick4")

pdf(file.path(panel.path, "6C.pdf"), width=embed.width, height=embed.height)
plot_grid(plot.dynamics(c("C1", "C2.2","C2.3"), dynamics = data$uns$communities$dynamics, cols = comm.cols, legend.position = "none", overlap.pseudotime = .1) + labs(x=NULL, y=NULL),
          plot.dynamics(paste0("C",c(1,3)), dynamics = data$uns$communities$dynamics, cols = comm.cols, legend.position = "none", overlap.pseudotime = .1) + labs(x=NULL, y=NULL),
          ncol=1)
while (!is.null(dev.list()))  dev.off()
rm(comm.cols)



# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel D - Communities trait associations                                    #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "6D.pdf"), width=embed.width*.5, height=embed.height)
plot.trait.associations(py_to_r(data$uns$communities$trait.association) %>% 
                          rename(state=covariate), 
                        params = names(AD.traits),
                        show.only.significant = F,
                        column_labels = AD.traits,
                        column_names_rot = 45,
                        column_names_centered = T,
                        cluster_rows=F,
                        row_names_side = "left",
                        use_raster=F,
                        raster_quality = 10,
                        border=T) 
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                               Panel E - Community pathways                                        #
# ----------------------------------------------------------------------------------------------------------------- #

states <- data$var %>% 
  mutate(grouping.by = recode(grouping.by, "Excitatory Neurons" = "excitatory","Inhibitory Neurons" = "inhibitory",
                              "Oligodendrocytes" = "oligodendrocytes", "Astrocyte" = "astrocytes","Microglia" = "microglia", "OPCs" = "opcs",
                              "Vascular Niche" = "endo") %>% as.character()) %>%
  rownames_to_column() %>%
  unstack(rowname~grouping.by)

pdf(file.path(panel.path, "6F.pdf"), width=embed.width, height=embed.height)
lapply(list(c("Mic.12","Mic.13"), c("OPC.1","Ast.3"), c("Ast.9","Ast.10","Oli.13")), function(vals)
  do.call(rbind, lapply(names(states), function(ct) 
    h5read(data.extended, file.path(mapping[[ct]], "pa")) %>% 
      filter(state %in% vals & direction == "upregulated") %>%
      dplyr::select(state, Description))) %>%
    split(., .$state) %>%
    lapply(., function(x) x$Description) %>%
    `[`(rev(vals)) %>%
    make_comb_mat %>% 
    t %>%
    UpSet(set_order = vals,
          comb_order = order(comb_degree(.)),
          border=TRUE, 
          width=unit(2,"cm"),
          top_annotation = NULL,
          right_annotation = upset_right_annotation(., gp = gpar(fill="grey80"))) %>%
    draw %>% 
    grid.grabExpr) %>%
  plot_grid(plotlist = ., ncol = 1)
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                               Panel E - Community pathways                                        #
# ----------------------------------------------------------------------------------------------------------------- #
clustered.pathways  <- openxlsx::read.xlsx("5. Manuscript code/data/community.pathways.xlsx", sheet = 1)
clustering.args     <- openxlsx::read.xlsx("5. Manuscript code/data/community.pathways.xlsx", sheet = 2, rowNames = T)
pathway.annotations <- openxlsx::read.xlsx("5. Manuscript code/data/figures.data.xlsx", sheet = "community.pathways")

settings <- expand.grid(rownames(clustering.args), c("upregulated")) %>% 
  `colnames<-`(c("comm","dirc")) %>%
  filter(comm %in% c("C4.2","C4.3"))

pdf(file.path(path, "6F.pdf"), width=embed.width*2, height=embed.height*2)
for(i in 1:nrow(settings)) {
  message(i)
  comm = as.character(settings[i, "comm"])
  dirc = settings[i, "dirc"]
  
  # Load community clustered pathways as well as "raw" state-specific pathways
  states <- data$var %>% 
    mutate(grouping.by = recode(grouping.by, 
                                "Excitatory Neurons" = "excitatory","Inhibitory Neurons" = "inhibitory",
                                "Oligodendrocytes" = "oligodendrocytes", "Astrocyte" = "astrocytes",
                                "Microglia" = "microglia", "OPCs" = "opcs","Vascular Niche" = "endo") %>% as.character()) %>% 
    filter(community == comm | sub.community == comm) %>%
    `[`(c("Mic.12","Mic.13","Ast.3","OPC.1","Ast.10","Oli.7","Ast.9","Mic.11"),) %>%
    rownames_to_column()
  states <- split(states, states$grouping.by) %>% lapply(., function(x) x$rowname)
  
  pathways <- do.call(rbind, lapply(names(states), function(ct) 
    h5read(data.extended, file.path(mapping[[ct]], "pa")) %>% filter(state %in% states[[ct]])))
  
  pathways <- clustered.pathways %>% 
    dplyr::select(-gene, -geneID) %>% 
    tidyr::separate_rows(state, sep="/") %>%
    merge(., pathways, 
          by.x = c("state","direction","ID", "Description"),
          by.y = c("state","direction","ID", "Description"),
          suffixes = c(".grouped",""))
  
  
  # Prepare state~pathways data for plotting
  df <- pathways %>% filter(direction == dirc) %>% 
    dplyr::select("state","direction","Description","p.adjust","membership") %>%
    reshape2::dcast(membership+direction+Description~state, value.var = "p.adjust", fun.aggregate = mean)
  
  # Determine clustered pathways ordering in heatmap
  ord <- df[,-(1:3)] %>% mutate_all(~replace(., is.na(.), 2)) %>% as.matrix %>% 
    Heatmap(show_row_dend = F, show_column_dend = F, row_split = df$membership) %>% prepare()
  row.ord <- row_order(ord)
  col.ord <- column_order(ord)
  
  # Plot clustered heatmap with pathway annotations
  row.splitting <- stack(row.ord) %>% arrange(values) %>% pull(ind)
  
  text <- pathway.annotations %>% filter(community == comm & direction == dirc) %>% tidyr::separate_rows(pathways, sep="/") %>% unstack(pathways~cluster)
  anno <- anno_textbox(align_to = row.splitting,
                       text = text,
                       gp = gpar(fontsize=12, col="black"),
                       background_gp = gpar(fill="grey95", col="grey80"),
                       max_width=unit(6, "cm"),
                       text_space = unit(8, "pt"))
  
  cols <- green2purple(3)
  if(dirc == "downregulated") 
    cols <- cols[1:2]
  else
    cols <- rev(cols[2:3])
  draw(Heatmap(df[stack(row.ord)[[1]],-(1:3)][,col.ord] %>% as.matrix,
          col = cols,
          row_split = row.splitting,
          cluster_rows = F, cluster_columns = F, show_row_names = F,
          border=T, row_gap = unit(0,"pt"),
          na_col = "grey90",
          right_annotation = rowAnnotation(pathways=anno),
          height = unit(10,"cm"),
          width = unit(.4*length(col.ord), "cm")))
  
}
while (!is.null(dev.list()))  dev.off()




# ----------------------------------------------------------------------------------------------------------------- #
#                                               Panel G - Mic.13-Ast.10 colocalization                              #
# ----------------------------------------------------------------------------------------------------------------- #
cols <- list(`Early/ABA`=green2purple(5)[2], prAD=green2purple(5)[4])
df <- readRDS("3. Other analyses/data/ST.validation.rds")$mic13.ast10.colocalization %>% 
  mutate(trajectory = case_when(trajectory %in% c("Early","ABA") ~ "Early/ABA", .default = trajectory))

pdf(file.path(panel.path, "6F.pdf"), width=embed.width.small*1.2, height=embed.height.small)
ggplot(df, aes(trajectory, cor, fill=trajectory)) + 
  geom_boxplot(outlier.size = 0) + 
  geom_point(aes(shape=is.sig)) + 
  scale_fill_manual(values=cols) + 
  scale_shape_manual(values = list("TRUE"=16, "FALSE"=1))+
  theme_classic()
while (!is.null(dev.list()))  dev.off()

# ----------------------------------------------------------------------------------------------------------------- #
#                                                 Panel H - Ast.10-Ast.5 exclusivity                                #
# ----------------------------------------------------------------------------------------------------------------- #
cols <- list(ABA=green2purple(5)[2], prAD=green2purple(5)[4], Early = "grey70")
df <- readRDS("3. Other analyses/data/ST.validation.rds")$ast10.ast5.colocalization


pdf(file.path(panel.path, "6H.1.pdf"), width=embed.width.small, height=embed.height)
df %>% pivot_longer(cols = c("Ast.10", "Ast.5"), names_to="state") %>% 
  select(trajectory, state, value) %>%
  ggplot(aes(trajectory, value, fill=trajectory)) + 
  geom_violin() +
  geom_point() +
  facet_wrap(~state, ncol=1, scales = "free") + 
  theme_classic() +
  theme(strip.background = element_blank(), 
        legend.position = "none") + 
  scale_fill_manual(values = cols)
while (!is.null(dev.list()))  dev.off()


pdf(file.path(panel.path, "6H.2.pdf"), width=embed.width, height=embed.height)
ggplot(df %>% filter(trajectory != "Early"), aes(trajectory, `10.vs.5`, fill=trajectory)) + 
  geom_boxplot(outlier.size = 0) + 
  geom_point(aes(shape=is.sig)) + 
  scale_fill_manual(values=cols) + 
  scale_shape_manual(values = list("TRUE"=16, "FALSE"=1))+
  theme_classic()
while (!is.null(dev.list()))  dev.off()



#####################################################################################################################
#                                       Supp Figure 6 - Communities Structure                                       #
#####################################################################################################################

# ----------------------------------------------------------------------------------------------------------------- #
#                             Panel A - Dynamics - states of interest - including points                            #
# ----------------------------------------------------------------------------------------------------------------- #

pdf(file.path(panel.path, "s9A.pdf"), width=embed.width, height=embed.height.small)
for(states in state.dynamics.sets) {
  print(plot.dynamics(names(states), cols = states, overlap.pseudotime=.1, label=TRUE, include.points = TRUE) + 
          theme(strip.text = element_blank(), 
                legend.position="none", 
                axis.line = element_line(), 
                axis.text.y = element_text()) + 
          labs(x=NULL, y=NULL, title=NULL))
}
while (!is.null(dev.list()))  dev.off()
rm(states)



# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel B -  Communities node graphs                                      #
# ----------------------------------------------------------------------------------------------------------------- #
library(ggraph)
library(tidygraph)

communities <- data$uns$communities
labelled <- c("Ast.1","Mic.2","Oli.3", "Mic.6","Mic.7","OPC.3","Ast.5","Mic.12","Mic.13","OPC.1","Oli.6","Ast.3","Oli.7","Ast.10")

graph <-
  # Convert similarity matrices to undirected flat edge dataframes
  lapply(c("correlation", "dynamics.adjacency"), function(n)
    communities$similarities[[n]] %>%
      `dimnames<-`(list(data$var_names, data$var_names)) %>%
      melt %>%
      filter(value != 0 & Var1 != Var2) %>%
      mutate(type=n)) %>%
  do.call(rbind, .) %>%
  
  # Convert to tidygraph - with workaround for duplicated edges
  as_tbl_graph(directed=FALSE) %>%
  activate(edges) %>%
  group_by(from, to, value, type) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  
  # Append community information
  activate(nodes) %>%
  mutate(community = data$var[name, "community"],
         sub.community = data$var[name, "sub.community"],
         grouping.by = data$var[name,"grouping.by"],
         label = ifelse(name %in% labelled, name, NA_character_)) %>%
  activate(edges) %>% 
  mutate(cross = .N()$community[from] != .N()$community[to]) %>%
  activate(nodes)


pos <- list(C1=c(0,0), C2 = c(-2, -.5), C3=c(2, -.5))
la <- lapply(names(pos), function(comm)
  graph %>% filter(community == comm) %>%
    activate(edges) %>% 
    filter(abs(value) >= .1) %>%
    create_layout(layout = "graphopt") %>%
    mutate(across(c(x,y), ~scales::rescale(., c(-.5,.5))),
           x = x+pos[[comm]][1], y = y+pos[[comm]][2])
) %>% do.call(rbind, .) %>%
  column_to_rownames("name") %>%
  `[`(graph %>% pull(name),)


pdf(file.path(panel.path, "s9B.pdf"),width=3*1.5, height=3)
for(configuration in list(list("correlation", .2, TRUE),
                          list("correlation", .2, FALSE),
                          list("dynamics.adjacency", .2, TRUE),
                          list("dynamics.adjacency", .2, FALSE))) {
  
  # Subset edges based on configuration
  .sub <- graph %>% activate(edges) %>% 
    filter(type == configuration[[1]] & 
             abs(value) >= configuration[[2]] & 
             cross == configuration[[3]]) %>% 
    activate(nodes)
  
  p <- ggraph(.sub, layout = "manual", x=la$x, y=la$y)
  if (configuration[[3]]) # If plotting between community edges
    p <- p + geom_edge_fan(aes(colour = value, edge_width=abs(value)/50), n=10, alpha=.6)
  
  p <- p + ggforce::geom_mark_hull(
    aes(x, y, group = community),
    fill="grey90", concavity = 5, expand = unit(2, "mm"), alpha = .6)
  
  if(!configuration[[3]]) # If plotting within community edges
    p <- p + geom_edge_fan(aes(colour = value, edge_width=abs(value)/50), n=10, alpha=.6)
  
  if(configuration[[1]] == "correlation") 
    p <- p + scale_edge_colour_gradient2(low=green2purple(2)[[1]], high = green2purple(2)[[2]], limit=c(-1,1))
  else
    p <- p + scale_edge_colour_gradient2(low="grey", high = "firebrick4", limit=c(0,1))
  
  p <- p + 
    scale_edge_width_binned(range=c(.1, 1), breaks=c(.2,.4,.6,.8,1)) + 
    geom_node_point(aes(fill=grouping.by), size=2, shape=21) + 
    scale_fill_manual(values = cell.group.color) + 
    geom_node_text(aes(filter = !is.na(label), label = label), colour = "black", repel = TRUE) + 
    theme(panel.background = element_blank())
  
  print(p)
}
dev.off()
rm(p, configuration, la, pos, graph, labelled, communities)



# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel C - Communities dynamics                                          #
# ----------------------------------------------------------------------------------------------------------------- #
comm.cols <- list(C1="springgreen4", C1.1="springgreen4", C1.2="olivedrab4", C2="darkorchid4", C2.1="springgreen4", C2.2="darkorchid4", C2.3="firebrick4", C3="firebrick4")

pdf(file.path(panel.path, "s9C.pdf"), width=embed.width, height=embed.height*1.5)
plot_grid(plot.dynamics(c("C1", "C2","C3"), dynamics = data$uns$communities$dynamics, cols = comm.cols, include.points = TRUE,  legend.position = c(3,3), overlap.pseudotime = .1) + labs(x=NULL, y=NULL),
          plot.dynamics(c("C1.1","C1.2"), dynamics = data$uns$communities$dynamics, cols = comm.cols, include.points = TRUE,legend.position = c(3,3), overlap.pseudotime = .1) + labs(x=NULL, y=NULL),
          plot.dynamics(c("C2.1","C2.2","C2.3"), dynamics = data$uns$communities$dynamics, cols = comm.cols, include.points = TRUE,legend.position = c(3,3), overlap.pseudotime = .1) + labs(x=NULL, y=NULL),
          ncol=1)
while (!is.null(dev.list()))  dev.off()
rm(comm.cols)



# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel D - Replication state dynamics                                        #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "s9D.pdf"), width=embed.width, height=embed.height.small)
for(states in state.dynamics.sets) {
  print(plot.dynamics(names(states), cols = states, data. = bulk, dynamics = bulk$uns$trajectories$dynamics, overlap.pseudotime=.2, label=TRUE, include.points = TRUE) + 
          theme(strip.text = element_blank(), 
                legend.position="none", 
                axis.line = element_line(), 
                axis.text.y = element_text()) + 
          labs(x=NULL, y=NULL, title=NULL))}
while (!is.null(dev.list()))  dev.off()
rm(states)



# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel E - ST participants in landscape                                      #
# ----------------------------------------------------------------------------------------------------------------- #
cols <- list(ABA=green2purple(2)[1], prAD=green2purple(2)[2], Early = "grey30")

df <- data$obsm$X_core_phate %>% `colnames<-`(c("PHATE_1","PHATE_2")) %>% filter(!is.na(PHATE_1)) %>% 
  merge(data$X[,c("Mic.13","Ast.10","Ast.5")], by.x="row.names", by.y="row.names") %>%
  merge(data.frame(data$uns$trajectories$palantir$branch.probs %>% py_to_r, pseudotime = data$uns$trajectories$palantir$pseudotime), by.x="Row.names", by.y="row.names") %>% 
  merge(readRDS("3. Other analyses/data/ST.validation.state.signatures.rds") %>% select(trajectory, sample) %>% unique, by.x="Row.names", by.y="sample", all.x=TRUE) %>% 
  arrange(!is.na(trajectory))

pdf(file.path(panel.path, "s9E.pdf"), width=embed.width, height=embed.height)
ggplot(df, aes(PHATE_1, PHATE_2, color=trajectory)) + 
  geom_point() + 
  ggrepel::geom_label_repel(aes(label=Row.names),data = df %>% filter(!is.na(trajectory))) + 
  scale_color_manual(values = cols, na.value = "lightgrey") + 
  theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(), 
        legend.position = "none") + 
  labs(x=NULL, y=NULL)
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel F - Mic.13-Ast.10 spatial correlation                                 #
# ----------------------------------------------------------------------------------------------------------------- #
cols <- list(ABA=green2purple(5)[2], prAD=green2purple(5)[4], Early = "grey70")
df <- readRDS("3. Other analyses/data/ST.validation.state.signatures.rds") %>%
  arrange(trajectory, sample) %>% mutate(sample=paste(trajectory, sample))


pdf(file.path(panel.path, "s9F.pdf"), width=embed.height.small*length(unique(df$sample)), height=embed.height.small)
ggplot(df , aes(Mic.13, Ast.10, group="1")) + 
  geom_point(size=.25, color="grey") + 
  geom_smooth(method="lm", color="black") + 
  facet_wrap(~sample, scales = "free", nrow=1) + 
  theme_classic() + 
  theme(strip.background = element_blank())
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel G - Ast.10-Ast.5 spatial incompatibility                              #
# ----------------------------------------------------------------------------------------------------------------- #
cols <- list(Ast.5=green2purple(5)[2], Ast.10=green2purple(5)[4])
df <- readRDS("3. Other analyses/data/ST.validation.state.signatures.rds") %>%
  arrange(trajectory, sample) %>% mutate(sample=paste(trajectory, sample))

pdf(file.path(panel.path, "s9G.pdf"), width=embed.height.small*.8*length(unique(df$sample)), height=embed.height.small)
df %>% pivot_longer(cols = c("Ast.10","Ast.5"), names_to = "state") %>% 
  arrange(trajectory) %>%
  ggplot(., aes(state, value, fill=state)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  facet_wrap(~sample, scales = "free", nrow=1) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        legend.position = "none") + 
  scale_fill_manual(values = cols)
while (!is.null(dev.list()))  dev.off()
