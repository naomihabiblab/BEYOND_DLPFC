library(openxlsx)
source("5. Manuscript code/utils.R")

creator <- "Gilad Sahar Green, Habib Lab, HUJI"
state.order <- lapply(atlas, "[[", "state.order") %>% unlist %>% unname %>% recode(., "NFOL.MOL"="NFOL/MOL")

header.style <- function() {
  createStyle(fontSize = 10, halign = "left", fgFill = "#808080", textDecoration = "bold")
}

add.header.style <- function(wb, sheet, ncols) {
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncols+1, gridExpand = TRUE)
}


# Swap donor identities with a sequential index 1:n
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
core.donors <- data$obs_names
ST.donors <- readRDS("Manuscript code/data/ST.validation.state.signatures.corrected.rds") %>% 
  select(donor.id = sample, trajectory=trajectory_group) %>% unique %>% 
  mutate(donor.id = gsub("_.*","",donor.id)) %>% 
  filter(trajectory != "Early") %>% 
  `rownames<-`(.$donor.id) %>% select(trajectory)

snuc.donors <- openxlsx::read.xlsx("5. Manuscript code/data/ROSMAP 10X Processing  Tracker.xlsx", sheet = "Processed Batches") %>%
  filter(!Batch %in% c("B1", "B2", "B3") & StudyCode != "MAP83034844") %>%
  pull(StudyCode) %>% unique %>% gsub("ROS|MAP", "", .) %>% as.numeric %>% as.character

bulk.donors <- data$uns$celmod$avg.predicted.prop$validation$index

donor.mapping <- c(snuc.donors, data$obs_names, bulk.donors) %>% unique
donor.mapping <- setNames(seq_along(donor.mapping), donor.mapping)
data$obs_names <- as.character(donor.mapping[data$obs_names])

stack(donor.mapping) %>% `colnames<-`(c("supp.table.id","ROSMAP.proj.id")) %>% 
  saveRDS("5. Manuscript code/data/Supp.to.ROSMAP.id.mapping.rds")

####################################################################################################################
##                            #  Supplementary Table 1 - Clinicopathological Characteristics  #                   ##
####################################################################################################################

wb <- createWorkbook(creator = creator, title = "Participants clinicopathological characteristics")
addWorksheet(wb, "Description", gridLines = TRUE)

# -------------------------------------------------------- #
# Create sheet for characteristics of snRNA-seq cohort     #
# -------------------------------------------------------- #
df <- load.metadata() %>% `[`(snuc.donors,) %>% 
  rownames_to_column("projid") %>%
  mutate(core.donor = projid %in% core.donors,
         projid = donor.mapping[projid]) %>%
  dplyr::select(donor.id=projid, core.donor, study, msex, age_death, pmi, cogng_demog_slope, cogdx, ceradsc, braaksc, niareagansc,
                sqrt.amyloid, sqrt.amyloid_mf, sqrt.tangles, sqrt.tangles_mf, mglia3_mf, mglia123_mf) %>%
  unique

addWorksheet(wb, (sheet <- "discovery cohort"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)



# -------------------------------------------------------- #
# Create sheet for characteristics of bulk cohort          #
# -------------------------------------------------------- #
df <- merge(load.metadata()[bulk.donors,], 
            read.csv("3. Other analyses/data/ROSMAP.bulk.RIN.values.csv", row.names = 1), 
            by.x="row.names", by.y="row.names", all.x=TRUE) %>%
  mutate(donor.id = donor.mapping[as.character(Row.names)]) %>%
  dplyr::select(donor.id, study, msex, age_death, pmi, RIN, cogng_demog_slope, cogdx,
                sqrt.amyloid, sqrt.amyloid_mf, sqrt.tangles, sqrt.tangles_mf)

addWorksheet(wb, (sheet <- "replication cohort"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


saveWorkbook(wb, "5. Manuscript code/data/Supplementary Table 1 - Participants Clinicopathological Characteristics.xlsx", overwrite = TRUE)
rm(df, sheet, wb)



####################################################################################################################
##                              #  Supplementary Table 2 - Atlas Characterization  #                              ##
####################################################################################################################

format.sheet <- function(wb, sheet, df, sets, horizontal.border.after=NULL) {
  columns <- setNames(LETTERS[0:ncol(df)+1], colnames(df))
  
  # identify rows of same state in order to apply different colors in sheet
  sets  <- cumsum(tidyr::replace_na(sets, FALSE))
  
  add.header.style(wb, sheet, ncol(df))
  
  # Even groups style
  addStyle(wb, sheet = sheet, gridExpand = TRUE,
           createStyle(fontSize = 11, fgFill = "#e0dede", borderStyle = "thin"),
           rows = which(as.logical(sets %% 2)) + 1, 
           cols = 0:ncol(df)+1)
  
  if(!is.null(horizontal.border.after)) {
    addStyle(wb, sheet, createStyle(border="Left", borderStyle ="thick"), rows = 0:nrow(df)+1, 
             cols = columns[c(horizontal.border.after)], 
             stack = TRUE, gridExpand = TRUE)
  }
  return(wb)
}


wb <- createWorkbook(creator = creator, title = "Atlas characterization")
addWorksheet(wb, "Description", gridLines = TRUE)

# -------------------------------------------------------- #
# Create sheet for DEGs of all subpopulations              #
# -------------------------------------------------------- #
df <- lapply(c("microglia","astrocytes","oligodendroglia","endo","excitatory","inhibitory"), function(ct)
  h5read(data.extended, file.path(mapping[[ct]], "de")) %>% 
    mutate(cell.type = ifelse(ct == "endo", "vascular niche", ct),
           cluster = factor(cluster, state.order, ordered = TRUE)) %>%
    dplyr::select(cell.type, state=cluster, geneID = id, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj) %>%
    arrange(state, desc(avg_log2FC))) %>%
  do.call(rbind, .)

addWorksheet(wb, "DEGs", gridLines = TRUE)
writeData(wb, sheet = "DEGs", df, rowNames = FALSE)
wb <- format.sheet(wb, "DEGs", df, df %>% mutate(v = state != lag(state)) %>% pull(v))


# -------------------------------------------------------- #
# Create sheet for pathways of all subpopulations          #
# -------------------------------------------------------- #
df <- lapply(c("microglia","astrocytes","oligodendroglia","endo","excitatory","inhibitory"), function(ct)
  h5read(data.extended, file.path(mapping[[ct]], "pa")) %>% 
    mutate(cell.type = ifelse(ct == "endo", "vascular niche", ct),
           state = factor(state, state.order, ordered = TRUE)) %>%
    dplyr::select(cell.type, state, direction, group=membership, pathway=Description, gene, dataset, pathwayID=ID, 
                  Count, GeneRatio.nums=GeneRatio.cp, GeneRatio, BgRatio.nums=BgRatio.cp, BgRatio, 
                  pvalue, p.adjust, qvalue,  page.rank, semantics, page.rank.rank, semantics.rank, p.adjust.rank, mean.rank, top.ranked) %>%
    arrange(state, desc(direction), group, qvalue)) %>%
  do.call(rbind, .)

addWorksheet(wb, "pathways", gridLines = TRUE)
writeData(wb, sheet = "pathways", df, rowNames = FALSE)
wb <- format.sheet(wb, "pathways", df, 
                   sets = df %>% mutate(v = state != lag(state) | direction != lag(direction) | group != lag(group)) %>% pull(v), 
                   horizontal.border.after = c("pathwayID","page.rank"))


# -------------------------------------------------------- #
# Create sheet for pairwise DEGs of all subpopulations     #
# -------------------------------------------------------- #
df <- lapply(c("microglia","astrocytes"), function(ct)
  h5read(data.extended, file.path(mapping[[ct]], "de.pairwise")) %>% 
    tidyr::separate(comparison, into=c("a","b"), sep="( vs. )", remove = FALSE) %>% 
    mutate(cell.type = ifelse(ct == "endo", "vascular niche", ct), 
           across(c(a,b), ~factor(., state.order, ordered = TRUE))) %>% 
    arrange(a,b, desc(avg_log2FC)) %>% 
    dplyr::select(cell.type, comparison, upregulated.in=cluster, geneID=id, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)) %>%
  do.call(rbind, .)

addWorksheet(wb, "DEGs.pairwise", gridLines = TRUE)
writeData(wb, sheet = "DEGs.pairwise", df, rowNames = FALSE)
wb <- format.sheet(wb, "DEGs.pairwise", df, 
                   sets = df %>% mutate(v = comparison != lag(comparison) | upregulated.in != lag(upregulated.in)) %>% pull(v))



# -------------------------------------------------------- #
# Create sheet for pairwise pathywas of all subpopulations #
# -------------------------------------------------------- #
df <- lapply(c("microglia","astrocytes"), function(ct)
  h5read(data.extended, file.path(mapping[[ct]], "pa.pairwise")) %>% 
    tidyr::separate(comparison, into=c("a","b"), sep="( vs. )", remove = FALSE) %>% 
    mutate(cell.type = ifelse(ct == "endo", "vascular niche", ct), 
           across(c(a,b), ~factor(., state.order, ordered = TRUE))) %>% 
    arrange(a,b, qvalue) %>% 
    dplyr::select(cell.type, comparison, upregulated.in=a, pathway=Description, gene, dataset, pathwayID=ID, 
                  Count, GeneRatio.nums=GeneRatio.cp, GeneRatio, BgRatio.nums=BgRatio.cp, BgRatio, 
                  pvalue, p.adjust, qvalue)) %>%
  do.call(rbind, .)

addWorksheet(wb, "pathways.pairwise", gridLines = TRUE)
writeData(wb, sheet = "pathways.pairwise", df, rowNames = FALSE)
wb <- format.sheet(wb, "pathways.pairwise", df, 
                   sets = df %>% mutate(v = comparison != lag(comparison) | upregulated.in != lag(upregulated.in)) %>% pull(v),
                   horizontal.border.after = c("pathway"))
  

saveWorkbook(wb, "5. Manuscript code/data/Supplementary Table 2 - Atlas Characterization.xlsx", overwrite = TRUE)
rm(df, wb, format.sheet)




####################################################################################################################
##                            #  Supplementary Table 3 - Endophenotype associations  #                            ##
####################################################################################################################
wb <- createWorkbook(creator = creator, title = "Endophenotypes Associations")
addWorksheet(wb, "Description", gridLines = TRUE)


# -------------------------------------------------------- #
# Participant QCs n465 -> n437                             #
# -------------------------------------------------------- #
df <- read.csv("2. Cell-type analysis/data/donors.qc.csv") %>%
  rowwise() %>%
  mutate(projid = donor.mapping[as.character(projid)]) %>% 
  rename(individualID = projid)

addWorksheet(wb, (sheet <- "Participant inclusion QCs"), gridLines = TRUE)
writeData(wb, sheet, df, startRow = 1, rowNames = FALSE)

addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# Subpopulation proportions matrix                         #
# -------------------------------------------------------- #
addWorksheet(wb, (sheet <- "snRNA-seq proportions"), gridLines = TRUE)
writeData(wb, sheet = sheet, 
          data$X %>% as.data.frame %>% `[`(, intersect(state.order, colnames(.))) %>% rownames_to_column("donor.id"), 
          startRow = 1, rowNames = FALSE)

addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(data$X)+1, gridExpand = TRUE)
addStyle(wb, sheet = sheet, header.style(), rows = 1:nrow(data$X)+2, cols = 1, gridExpand = TRUE)


# -------------------------------------------------------- #
# CelMod predicted proportions                             #
# -------------------------------------------------------- #
df <- data$uns$celmod$predicted.proportions %>% py_to_r %>% 
  rowwise() %>%
  mutate(projid = donor.mapping[as.character(projid)]) %>% 
  rename(donor.id = projid)

addWorksheet(wb, (sheet <- "CelMod predicted proportions"), gridLines = TRUE)
writeData(wb, sheet, df, startRow = 1, rowNames = FALSE)

addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)



# -------------------------------------------------------- #
# CelMod predictions correlations                          #
# -------------------------------------------------------- #
df <- data$uns$celmod$test.corrs %>% py_to_r %>% rownames_to_column("state")

addWorksheet(wb, (sheet <- "CelMod correlations"), gridLines = TRUE)
writeData(wb, sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# Subpopulaiton-endophenotype associations - both cohorts  #
# -------------------------------------------------------- #
df <- rbind(data$uns$trait.analysis$snuc %>% py_to_r %>% mutate(cohort = "discovery"),
            data$uns$trait.analysis$celmod %>% py_to_r %>% mutate(cohort = "replication")) %>%
  dplyr::select(cohort, everything())

addWorksheet(wb, (sheet <- "Endophenoptye associations"), gridLines = TRUE)
writeData(wb, sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# Meta-analysis of endophenotype associations              #
# -------------------------------------------------------- #
df <- data$uns$trait.analysis$meta.analysis %>% py_to_r %>% 
  dplyr::select(trait, state, tstat.sc, n.sc, tstat.b, n.b, z.meta, 
                pval.meta, adj.pval.meta, sig.meta)

addWorksheet(wb, (sheet <- "Endophenoptye meta-analysis"), gridLines = TRUE)
writeData(wb, sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)

saveWorkbook(wb, "5. Manuscript code/data/Supplementary Table 3 - Endophenotypes Associations.xlsx", overwrite = TRUE)
rm(df, wb, sheet)



####################################################################################################################
##                                   #  Supplementary Table 4 - smFISH analysis  #                                ##
####################################################################################################################
validations <- readRDS("3. Other analyses/data/RNAscope.rds")

wb <- createWorkbook(creator = creator, title = "smFISH quantifications")
addWorksheet(wb, "Description", gridLines = TRUE)

validation.donor.mapping <- validations$RNAscope$df$sample %>% unique
validation.donor.mapping <- setNames(seq_along(validation.donor.mapping), validation.donor.mapping)

# -------------------------------------------------------- #
# Subpopulation proportions matrix                         #
# -------------------------------------------------------- #
df <- validations$RNAscope$df %>%
  rename(donor.id = sample) %>%
  mutate(donor.id = validation.donor.mapping[as.character(donor.id)])

addWorksheet(wb, (sheet <- "Cell-wise values"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, startRow = 2, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 2, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# Subpopulation proportions matrix                         #
# -------------------------------------------------------- #
df <- validations$RNAscope$predicted.proportions %>% 
  rownames_to_column("donor.id") %>%
  mutate(donor.id = validation.donor.mapping[as.character(donor.id)])

addWorksheet(wb, (sheet <- "Subpopulation proportion"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, startRow = 1, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# pTau quantifications                                     #
# -------------------------------------------------------- #
df <- validations$pTau$summarised %>% 
  select(donor.id = sample, everything()) %>% 
  mutate(donor.id = validation.donor.mapping[as.character(donor.id)])

addWorksheet(wb, (sheet <- "Donor pTau quantification"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, startRow = 1, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)



# -------------------------------------------------------- #
# smFISH morphology association                            #
# -------------------------------------------------------- #
df <- validations$RNAscope$morpholoy.association

addWorksheet(wb, (sheet <- "Morphology association"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, startRow = 1, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# ROSMAP microglial density analysis                       #
# -------------------------------------------------------- #
df <- validations$ROSMAP.PAM.Association %>% 
  `colnames<-`(c("beta","se", "tstat", "pval", "adj.pval", "sig")) %>%
  rownames_to_column("covariate")

addWorksheet(wb, (sheet <- "PAM score associations"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, startRow = 1, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


saveWorkbook(wb, "5. Manuscript code/data/Supplementary Table 4 - smFISH quantifications.xlsx", overwrite = TRUE)
rm(df, validations, wb, sheet)



####################################################################################################################
##                                   #  Supplementary Table 5 - BEYOND analysis  #                                ##
####################################################################################################################
wb <- createWorkbook(creator = creator, title = "BEYOND analysis results")
addWorksheet(wb, "Description", gridLines = TRUE)


# -------------------------------------------------------- #
# Cellular landscape 3d embedding                          #
# -------------------------------------------------------- #
df <- data.frame(donor.id = data$obs_names, 
                 cluster = data$obs$clusters, 
                 data$obsm$X_all_3d_phate %>% `colnames<-`(paste0("PHATE_", 1:3)))

addWorksheet(wb, (sheet <- "3D Landscape embedding"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# Cellular landscape 2d embedding                          #
# -------------------------------------------------------- #
df <- data.frame(donor.id = data$obs_names, 
                 cluster = data$obs$clusters, 
                 data$obsm$X_core_phate %>% `colnames<-`(paste0("PHATE_", 1:2)))

addWorksheet(wb, (sheet <- "2D Landscape embedding"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# Donors' pseudotime and trajectory probability            #
# -------------------------------------------------------- #
algs <- list(via="VIA", palantir="Palantir")

for(a in names(algs)) {
  traj      <- data$uns$trajectories[[a]]
  terminals <- traj$terminals %>% py_to_r %>% select(x=terminal) %>% rownames_to_column("terminal") %>% column_to_rownames("x")
  root      <- setNames(c(T), data$uns$trajectories[[a]]$user.root)
  
  df <- data.frame(pseudotime = traj$pseudotime,
                   traj$branch.probs %>% py_to_r) %>%
    rownames_to_column("individualID") %>%
    mutate(terminal = terminals[individualID,"terminal"],
           ST.include = ST.donors[individualID,"trajectory"],
           root = root[individualID],
           individualID = donor.mapping[individualID])
  
  
  addWorksheet(wb, (sheet <- paste(algs[[a]], "Trajectories")), gridLines = TRUE)
  writeData(wb, sheet = sheet, df, rowNames = FALSE)
  addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)
}
rm(algs, a, traj, terminals, df)


# -------------------------------------------------------- #
# Subpopulation community assignment                       #
# -------------------------------------------------------- #
df <- data$var %>% dplyr::select(community, sub.community) %>% 
  rownames_to_column("state") %>% 
  mutate(state = factor(state, state.order, ordered = TRUE)) %>%
  arrange(state)

addWorksheet(wb, (sheet <- "Subpopulation communities"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------- #
# Donor community proportions                              #
# -------------------------------------------------------- #
df <- cbind(data$obsm$communities, data$obsm$sub.communities) %>% rownames_to_column("donor.id")
  
addWorksheet(wb, (sheet <- "Participant comm. prop."), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)



# -------------------------------------------------------------- #
# Community-endophenotype associations                           #
# -------------------------------------------------------------- #
df <- data$uns$communities$trait.association %>% py_to_r

addWorksheet(wb, (sheet <- "Community endophe. assoc."), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)



# ------------------------------------------------------------ #
# State, community and endophenotype dynamics - fitted values  #
# ------------------------------------------------------------ #
df <- rbind(data$uns$trajectories$palantir$dynamics$fitted.vals %>% py_to_r,
            data$uns$communities$dynamics$fitted.vals %>% py_to_r %>% filter(feature %in% paste0("C", 1:5))) %>%
  dplyr::select(feature, trajectory, pseudotime=x, value=y, weight=X.weights., fit, se.fit) %>% 
  mutate(trajectory = recode(trajectory, "ABA"="Alternative Aging", "prAD" = "prAD"))

addWorksheet(wb, (sheet <- "Dynamics fitted values"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


# -------------------------------------------------------------- #
# State, community and endophenotype dynamics - predicted values #
# -------------------------------------------------------------- #
df <- rbind(data$uns$trajectories$palantir$dynamics$pred.vals %>% py_to_r,
            data$uns$communities$dynamics$pred.vals %>% py_to_r %>% filter(feature %in% paste0("C", 1:5))) %>%
  dplyr::select(feature, trajectory, pseudotime=x, fit, se.fit) %>% 
  mutate(trajectory = recode(trajectory, "ABA"="Alternative Aging", "prAD" = "prAD"))

addWorksheet(wb, (sheet <- "Dynamics predicted values"), gridLines = TRUE)
writeData(wb, sheet = sheet, df, rowNames = FALSE)
addStyle(wb, sheet = sheet, header.style(), rows = 1, cols = 0:ncol(df)+1, gridExpand = TRUE)


saveWorkbook(wb, "5. Manuscript code/data/Supplementary Table 5 - BEYOND analysis results.xlsx", overwrite = TRUE)
