source("5. Manuscript code/utils.R")


#####################################################################################################################
#                                  Figure 5 - Cellular Landscape & State/Trait Dynamics                             #
#####################################################################################################################


# ----------------------------------------------------------------------------------------------------------------- #
#                                            Panel B - Branch Probabilities                                         #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "5B.pdf"), width=embed.width.small*.65, height=embed.height.small)
plot.landscape(data$uns$trajectories$palantir$branch.probs %>% py_to_r %>% mutate(diff=prAD-ABA) %>% dplyr::select(diff),
               smoothened = F,
               size = 1,
               cols = green2purple.less.white,
               legend.position = c(3,3))
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel C - State Prevalence PHATE                                        #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "5C.pdf"), width=embed.width.small*.75*4, height=embed.height.small*2)
plot.landscape(c("Ast.1", "Mic.12","Ast.10","Ast.5", "Mic.2","Mic.13","Oli.7","OPC.3"), enforce.same.color.scale = FALSE, smoothened = TRUE, ncol=4) %>% print(.)
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel D - Trait Prevalence PHATE                                        #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "5D.pdf"), width=embed.width.small*.75, height=embed.height.small*3)
plot_grid(
  plot.landscape(c("sqrt.amyloid_mf","sqrt.tangles_mf"), enforce.same.color.scale = FALSE, ncol = 1),
  plot.landscape(c("cogng_demog_slope"), cols = function(n) colorRampPalette(rev(c("#30123BFF", "#4777EFFF", "#1BD0D5FF", "#62FC6BFF", "#D2E935FF", "#FE9B2DFF", "#DB3A07FF")))(n), sort.direction = -1, enforce.same.color.scale = FALSE, ncol = 1),
  ncol=1, rel_heights = c(2,1)
)
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                              Panel E - Trait Dynamics                                             #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "5E.pdf"), width=embed.width, height=embed.height*1.75)
plot_grid(plotlist = lapply(names(AD.traits), function(trait)
  plot.dynamics(c(trait), dynamics = data$uns$trajectories$palantir$dynamics, cols=AD.traits.colors,
                overlap.pseudotime=.112, ncol=2, strip.position="top", scales="free",
                label = F, legend.position = c(2.5,0), include.points=F) + 
    labs(x=NULL, y=trait, title=NULL)),
  ncol=1) %>% print()
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                              Panel F - Landscape Bulk                                             #
# ----------------------------------------------------------------------------------------------------------------- #
bulk <- anndata::read_h5ad("4. BEYOND/data/Celmod.subpopulation.proportion.h5ad")

# prAD-ABA landscape
pdf(file.path(panel.path, "5F.traj.prob.pdf"), width=embed.width.small, height=embed.height.small)
plot.landscape(bulk$uns$trajectories$branch.probs %>% py_to_r %>% mutate(diff=prAD.like-ABA.like) %>% dplyr::select(diff),
               "X_umap",
               smoothened = F,
               size = 1,
               cols = green2purple.less.white,
               legend.position = c(3,3), data. = bulk)
while (!is.null(dev.list()))  dev.off()

# Trajectories root point
pdf(file.path(panel.path, "5F.root.pdf"), width=embed.width.small, height=embed.height.small)
plot.landscape(data$obs %>% mutate(root = if_else(data$uns$trajectories$pseudotime ==0, "root", NA_character_)) %>% 
                 dplyr::select(root), "X_umap",data. = bulk)
while (!is.null(dev.list()))  dev.off()


# Landscape of key states
pdf(file.path(panel.path, "5F.states.pdf"), width=embed.width.small, height=embed.height.small)
plot.landscape(c("Ast.1","Ast.5","Ast.10","Mic.13"), "X_umap", enforce.same.color.scale = FALSE, size = 1, data. = bulk)
while (!is.null(dev.list()))  dev.off()


# AD trait dynamics
pdf(file.path(panel.path, "5F.traits.pdf"),width=embed.width, height=embed.height*1.75)
plot_grid(plotlist = lapply(names(AD.traits), function(trait)
  plot.dynamics(c(trait), dynamics = bulk$uns$trajectories$dynamics, cols=AD.traits.colors,
                overlap.pseudotime=.2, ncol=2, strip.position="top", scales="free",
                label = F, legend.position = c(2.5,0), include.points=F) + 
    labs(x=NULL, y=trait, title=NULL)),
  ncol=1) %>% print()
while (!is.null(dev.list()))  dev.off()

rm(bulk)



#####################################################################################################################
#                              Supp Figure 5 - Cellular Landscape & State/Trait Dynamics                            #
#####################################################################################################################


# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel A - Plain 3D PHATE                                                #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "s8A.pdf"), width=embed.width, height=embed.height)
plot.landscape.3D(c("clusters"), theta = 135, phi = 30, cols = scales::hue_pal(), legend.position = "left")
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                          Panel B - 3D PHATE from additional states                                #
# ----------------------------------------------------------------------------------------------------------------- #

pdf(file.path(panel.path, "s8B.pdf"), width=embed.width*2, height=embed.height)
plot.landscape.3D(c("Mic.6","Ast.9","Ast.2","Mic.7","Oli.8","Oli.4"), theta = 135, phi = 30, smoothened = TRUE)
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel C - State Prevalence UMAP+tSNE                                    #
# ----------------------------------------------------------------------------------------------------------------- #
source("4. BEYOND/utils/utils.R")
sc <- reticulate::import("scanpy")

sub <- data[,data$uns$celmod$celmod.states]
sc$external$tl$phate(sub, k = as.integer(10), n_components = as.integer(2),  a = as.integer(40), knn_dist =  "cosine", 
                     mds_dist = "cosine", mds_solver = "smacof", verbose = F)
sub$obsp[["similarity_X_phate"]] <- embedding.similarity(sub$obsm[["X_phate"]], knn = 5)

pdf(file.path(panel.path, "s8C.pdf"), width=embed.width*2.5, height=embed.height)
args <- list(features = c("Ast.1","Mic.2","Mic.12","Mic.13","Oli.7","Ast.10","OPC.3","Ast.5","Ast.2","Ast.9"),
             enforce.same.color.scale = FALSE, 
             nrow=1, size=1,
             data. = data)
p1 <- plot_grid(do.call(plot.landscape, modifyList(args, list(embedding = "X_umap"))),
                do.call(plot.landscape, modifyList(args, list(embedding = "X_tsne"))),
                do.call(plot.landscape, modifyList(args, list(embedding = "X_phate", data.=sub))),
                ncol=1)

args <- list(features = c("clusters"), cols = scales::hue_pal(), size=1, nrow=1, data. = data)
p2 <- plot_grid(do.call(plot.landscape, modifyList(args, list(embedding = "X_umap"))),
          do.call(plot.landscape, modifyList(args, list(embedding = "X_tsne"))),
          do.call(plot.landscape, modifyList(args, list(embedding = "X_phate", data.=sub))),
          ncol=1)

plot_grid(p2, p1, ncol=2, rel_widths = c(1,10)) %>% print(.)
while (!is.null(dev.list()))  dev.off()
rm(sub, sc, args, p1, p2)



# ----------------------------------------------------------------------------------------------------------------- #
#                                            Panel D - VIA model over 3D PHATE                                      #
# ----------------------------------------------------------------------------------------------------------------- #
df <- data.frame(pseudotime=data$uns$trajectories$via$pseudotime,
                 entropy=apply(py_to_r(data$uns$trajectories$via$branch.probs), 1, function(i) -sum(i*log2(i)) %>% ifelse(is.nan(.), 0, .)),
                 py_to_r(data$uns$trajectories$via$branch.probs))

pdf(file.path(panel.path, "s8D.pdf"), width=embed.width, height=embed.height*1.5)
plot.landscape.3D(df, ncol=2, theta = 135, phi = 30, smoothened = F, legend.position = c(3,3), enforce.same.color.scale = F, size=1.5)
while (!is.null(dev.list()))  dev.off()
rm(df)



# ----------------------------------------------------------------------------------------------------------------- #
#                                           Panel E - Palantir model over PHATE                                     #
# ----------------------------------------------------------------------------------------------------------------- #
df <- data.frame(pseudotime=data$uns$trajectories$palantir$pseudotime,
                 entropy=apply(py_to_r(data$uns$trajectories$palantir$branch.probs), 1, function(i) -sum(i*log2(i)) %>% ifelse(is.nan(.), 0, .)),
                 py_to_r(data$uns$trajectories$palantir$branch.probs))

pdf(file.path(panel.path, "s8E.pdf"), width=embed.width*.85, height=embed.height)
plot.landscape(df, ncol=2, smoothened = F, legend.position = c(3,3), enforce.same.color.scale = F, size=1.5)
while (!is.null(dev.list()))  dev.off()
rm(df)



# ----------------------------------------------------------------------------------------------------------------- #
#                                            Panel F - Palantir VIA comparison                                      #
# ----------------------------------------------------------------------------------------------------------------- #

m <- lapply(data$uns$trajectories, function(t) list(prob = as.matrix(py_to_r(t$branch.probs)), ps = t$pseudotime))
sig <- outer(colnames(m$palantir$prob), colnames(m$via$prob), 
      Vectorize(function(p,v) cor.test(m$palantir$prob[,p], m$via$prob[,v], use="pairwise.complete.obs")[["p.value"]])) %>% 
  p.adjust(method = "BH") %>% 
  cut(., c(-.1, .0001, .001, .01, .05, Inf), c("****","***", "**", "*", "")) %>% matrix(ncol=ncol(m$via$prob))


pdf(file.path(panel.path, "s8F1.pdf"), width=embed.width, height=embed.height)
Heatmap(cor(m$palantir$prob, m$via$prob, use = "pairwise.complete.obs"),
        col = circlize::colorRamp2(seq(-1,1, length.out=11), green2purple(11)),
        cell_fun = function(j, i, x, y, w, h, fill) grid.text(sig[i,j], x,y),
        cluster_rows = F, cluster_columns = F)
while (!is.null(dev.list()))  dev.off()

pdf(file.path(panel.path, "s8F2.pdf"), width=embed.width, height=embed.height)
  ggplot(data.frame(p=m$palantir$ps, v=m$via$ps) %>% arrange(p),aes(p,v)) + 
    geom_point(size=1) + 
    geom_abline() + 
    scale_x_continuous(expand = c(.01,.01)) + 
    scale_y_continuous(expand = c(.01,.01)) + 
    ggpubr::stat_cor(aes(label = paste("Pearson", after_stat(r.label), after_stat(p.label), sep = "~`,`~")), label.y = .95) +
    ggpubr::stat_cor(aes(label = paste("Spearman", after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method="spearman", label.y = .9) + 
    theme_classic() + 
    labs(x="Pseudotime: Palantir", y="Pseudotime: VIA")
while (!is.null(dev.list()))  dev.off()



# ----------------------------------------------------------------------------------------------------------------- #
#                                            Panel G - Entropy over pseudotime                                      #
# ----------------------------------------------------------------------------------------------------------------- #
df <- data.frame(pseudotime=data$uns$trajectories$palantir$pseudotime,
                 entropy=apply(py_to_r(data$uns$trajectories$palantir$branch.probs), 1, function(i) -sum(i*log2(i)) %>% ifelse(is.nan(.), 0, .)),
                 py_to_r(data$uns$trajectories$palantir$branch.probs))

pdf(file.path(panel.path, "s8G.pdf"), width=embed.width.small, height=embed.height)
ggplot(df, aes(pseudotime, entropy, color=prAD-ABA)) + 
  geom_ribbon(aes(x,ymax=ymax,ymin=ymin),
              data.frame(x=c(0, .1), ymax=rep(1-.01, 2), ymin=rep(0,2)), 
              fill="black", alpha=.05, inherit.aes = F) + 
  geom_vline(xintercept = .1, linetype="dashed") +
  geom_point(size=1.5) + 
  scale_color_gradientn(colors = green2purple.less.white(21),
                        seq(-1,1, length.out=21), na.value = "lightgrey", aesthetics = "color") + 
  scale_x_continuous(expand = c(0,0),
                     labels=function(x) recode(x, `0`="0", `1`="1", .default = as.character(x))) + 
  scale_y_continuous(expand = expansion(add = .01),
                     labels=function(x) recode(x, `0`="0", `1`="1", .default = as.character(x))) + 
  theme_classic() + 
  theme(legend.position=c(3,3))
while (!is.null(dev.list()))  dev.off()
rm(df)

  

# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel H - Trait Dynamics - With points                                      #
# ----------------------------------------------------------------------------------------------------------------- #
pdf(file.path(panel.path, "s8H.pdf"), width=embed.width.small, height=embed.height)
plot_grid(plotlist = lapply(names(AD.traits), function(trait)
  plot.dynamics.wrapper(data$uns$trajectories$palantir$dynamics, c(trait), cols=AD.traits.colors[trait],
                        overlap.pseudotime=.1, ncol=2, strip.position="top", scales="free", size=1,
                        label = F, legend.position = c(2.5,0), include.points=TRUE) + 
    labs(x=NULL, y=trait, title=NULL)),
  ncol=1) %>% print()
while (!is.null(dev.list()))  dev.off()

  


# ----------------------------------------------------------------------------------------------------------------- #
#                                       Panel I - Replication Palantir model                                        #
# ----------------------------------------------------------------------------------------------------------------- #

bulk <- anndata::read_h5ad("4. BEYOND/data/Celmod.subpopulation.proportion.h5ad")
df <- data.frame(pseudotime=bulk$uns$trajectories$pseudotime,
                 entropy=apply(py_to_r(bulk$uns$trajectories$branch.probs), 1, function(i) -sum(i*log2(i)) %>% ifelse(is.nan(.), 0, .)),
                 py_to_r(bulk$uns$trajectories$branch.probs))

pdf(file.path(panel.path, "s8I.clusters.pdf"), width=embed.width, height=embed.height)
plot.landscape(c("clusters"),"X_umap", legend.position = c(3,3), cols = scales::hue_pal(), data. = bulk)
while (!is.null(dev.list()))  dev.off()


pdf(file.path(panel.path, "s8I.palantir.pdf"), width=embed.width*3, height=embed.height)
plot.landscape(df %>% dplyr::select(-entropy), "X_umap", ncol=3, smoothened = F, legend.position = c(3,3), enforce.same.color.scale = F, data. = bulk)
while (!is.null(dev.list()))  dev.off()


pdf(file.path(panel.path, "s8J.entropy.drop.pdf"), width=embed.width.small, height=embed.height)
ggplot(df, aes(pseudotime, entropy, color=prAD.like-ABA.like)) + 
  geom_ribbon(aes(x,ymax=ymax,ymin=ymin),
              data.frame(x=c(0, .2), ymax=rep(1-.01, 2), ymin=rep(0,2)), 
              fill="black", alpha=.05, inherit.aes = F) + 
  geom_vline(xintercept = .2, linetype="dashed") +
  geom_point(size=1.5) + 
  scale_color_gradientn(colors = green2purple.less.white(21),
                        seq(-1,1, length.out=21), na.value = "lightgrey", aesthetics = "color") + 
  scale_x_continuous(expand = c(0,0),
                     labels=function(x) recode(x, `0`="0", `1`="1", .default = as.character(x))) + 
  scale_y_continuous(expand = expansion(add = .01),
                     labels=function(x) recode(x, `0`="0", `1`="1", .default = as.character(x))) + 
  theme_classic() + 
  theme(legend.position=c(3,3))
while (!is.null(dev.list()))  dev.off()
rm(df)
