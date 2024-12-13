# BEYOND - Cellular communities reveal trajectories of brain ageing and Alzheimer’s disease

<p align="justify">
To reconstruct the dynamics of the brain's cellular environment along the progression of Alzheimer's Disease, and to distinguish between AD and aging effects, we built a comprehensive cell atlas of the aged prefrontal cortex from 1.65 million single-nucleus RNA-seq profiles. We associated glial, vascular and neuronal subpopulations with AD-related traits for 437 aged individuals and aligned them along the disease cascade using causal modeling. To model the coordinated dynamics of the entire cellular environment, we devised the BEYOND methodology, which uncovered two distinct trajectories of brain aging defined by differing sequences of changes in cellular communities
</p>

<p align="center"><img src="https://github.com/user-attachments/assets/554fdf75-0160-4c8a-93e1-dd45a601cd85" width="75%"></p>

<details>
<summary>Available infomration in repository</summary>
<p align="justify"> 
1. **Library preprocessing**: Library level snRNA-seq analysis. Removal of background RNA, quality control steps, cell-type classification and doublets annotation.
2. **Cell-type analysis**: snRNA-seq analysis for each of the cell types. Creating our cell atlas. Removal of low quality cells and doublets, sub-clustering analysis, differential expression and pathway enrichment analyses.
3. **Other analyses**: Celmod bulk RNA-seq deconvolution, trait associations and meta-analysis, causal modeling, RNAscope validations and pTau validations.
4. **BEYOND**: Applying BEYOND methodology steps over our cell atlas, as well as validation over replication cohort.
5. **Manuscript code**: The complete code generating all figures, extended data figures and supplementary tables shown in the manuscript.

The code within each folder runs under the assumption that the code of the previous folder was executed and generated the necessary data files. Please refer to the readme file of each folder for more details
</details>

<details>
<summary>Data accessibility</summary>
All snRNA-seq data used in our study is accessible via [Synapse (syn53366818)](https://www.synapse.org/#!Synapse:syn53366818) and contains the raw reads, library count matrices and processed cell atlas in the format of cell-type Seurat objects. 
</details>

--- 
## Reproducibility and further analysis
To facilitate further investigation into our cellular landscape and BEYOND methodology we provide the following tutorial in which we load the manuscript's supplementary tables and reproduce the different types of plots shown in our study. Clone this GitHub repository and follow the code below.
```R
# Clone the current GitHub repository and specify the path to the following supplementary tables

library(openxlsx)
SUPP = list(
  table.1 = "Supplementary Tables/Supplementary Table 1 - Participants Clinicopathological Characteristics.xlsx",
  table.3 = "Supplementary Tables/Supplementary Table 3 - Endophenotypes Associations.xlsx",
  table.5 = "Supplementary Tables/Supplementary Table 5 - BEYOND analysis results.xlsx")
```

<details>
<summary>Subpopulation proportion matrix </summary>
The starting point for all post-atlas analyses is the participant-wise subpopulation proportions matrix. Computation of the subpopulation proportion matrix starts with the extraction of the participant-cell annotations from the Seurat objects and the calculation of the participant-wise proportions. This code can be found in the following files:

-   *2.Cell-type analysis/2.aggregate.cell.type.data.R* lines 45-77.
-   *2.Cell-type analysis/3.create.proportion.matrix.R*

To recreate the subpopulation proportion matrix from the supplementary tables use the following:

```R
# Create and load a conda environment with the necessary R and python packages
source("2. Cell-type analysis/utils/load.code.env.R")
library(anndata)
library(dplyr)
library(tibble)

proportions <- read.xlsx(SUPP$table.3, "snRNA-seq proportions") %>% column_to_rownames("individualID")
clinical.information <- read.xlsx(SUPP$table.1, "Discovery cohort") %>% column_to_rownames("individualID")
qcs <- read.xlsx(SUPP$table.3, "Participant inclusion QCs") %>% filter(Final_QC == "Pass") %>% column_to_rownames("individualID")

# We represent the cellular landscape using an AnnData object. In its core `data$X` is the subpopulation proportion matrix.
# This representation provides us a data structure into which we can load participant- (rows) and subpopulation (column) information
data <- AnnData(
  X = proportions,
  layers = list(sqrt.prev = sqrt(proportions)),
  obsm = list(
    QCs = qcs[rownames(proportions), ],
    meta.data = clinical.information[rownames(proportions), ])
)

# The following `cogdx_ad` is derived from the `cogdx` characterization and is used for subpopulation-endophenotype associations (while excluding MCI or AD with another cause of CI)
data$obsm$meta.data$cogdx_ad = as.numeric(recode(data$obsm$meta.data$cogdx, "1"="1","2"="2", "3"=NA_character_, "4"="3", "5"=NA_character_,"6"=NA_character_))

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(proportions, clinical.information, qcs)
```
</details>

<details>
<summary>Endophenotype associations & meta-analysis</summary>
The populated data-structure above is sufficient for re-running the Celmod fitting, snRNA-seq and bulk RNA-seq endophenotype associations, meta-analysis of associations, and the causal modeling (see *"3. Other analyses"* folder). To load Celmod results from the supplementary tables use the following:

```R
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

data$uns$celmod <- list()
data$uns$celmod$predicted.proportions <- read.xlsx(SUPP$table.3, "CelMod predicted proportions")
data$uns$celmod$avg.predicted.prop <- 
  data$uns$celmod$predicted.proportions %>% py_to_r %>% split(., .$set) %>% 
  lapply(., function(df) dcast(df, individualID~subpopulation, value.var = "sqrt.prev_mean") %>% 
           column_to_rownames("individualID") %>%
           `[`(,colnames(data$X)))
data$uns$celmod$test.corrs    <- read.xlsx(SUPP$table.3, "CelMod correlations") %>% column_to_rownames("subpopulation")
data$uns$celmod$celmod.states <- data$uns$celmod$test.corrs %>% py_to_r %>% filter(adj.pval < .01 & corr > 0) %>% rownames()
data$uns$celmod$shared.donors <- data$uns$celmod$avg.predicted.prop$train$index

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
```

To load endophenotype associations and their meta-analysis (discovery and replication) from the supplementary tables use the following:

```R
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
associations <- read.xlsx(SUPP$table.3, "Endophenotype associations") %>% split(., .$cohort)

data$uns$trait.analysis <- list(
  snuc = associations$discovery, 
  celmod = associations$replication)

data$uns$trait.analysis$meta.analysis <- 
  merge(py_to_r(data$uns$trait.analysis$snuc),
              py_to_r(data$uns$trait.analysis$celmod),
              by = c("trait","state"),
              suffixes = c(".sc",".b"),
              all.x = T) %>% 
    merge(., py_to_r(data$uns$celmod$test.corrs) %>% `colnames<-`(paste0(colnames(.),".celmod")),
          by.x = "state",
          by.y = "row.names") %>% arrange(-corr.celmod) %>% 
    merge(., read.xlsx(SUPP$table.3, "Endophenotype meta-analysis") %>% rename("state"="subpopulation") %>% select(-ends_with(".sc"), -ends_with(".b")), by=c("trait","state"))

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(associations)
```
</details>

<details>
<summary>BEYOND - Cellular landscape embedding</summary>
Using the participant-wise subpopulation proportions as our representation we can now visualize the cellular landscape. The code for computing the used embeddings can be found in file *4. BEYOND/1.create.cellular.landscape.R*. To load embeddings from the supplementary tables use the following:

```R
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

data$obs$clusters <- read.xlsx(SUPP$table.5, "3D Landscape embedding") %>% 
    column_to_rownames("individualID") %>%
    `[`(rownames(data), "cluster")
data$obs$core <- !data$obs$clusters %in% c(9,10)

data$obsm$X_all_3d_phate <- read.xlsx(SUPP$table.5, "3D Landscape embedding") %>% 
  column_to_rownames("individualID") %>% 
  select(-cluster) %>% `[`(rownames(data), )
data$obsm$X_core_phate <- read.xlsx(SUPP$table.5, "2D Landscape embedding") %>% 
  column_to_rownames("individualID") %>% 
  select(-cluster) %>% `[`(rownames(data), )

# Compute embedding local density for smoothened landscape plots
source("4. BEYOND/utils/utils.R")
data$obsp <- list()
for(e in c("X_all_3d_phate","X_core_phate")) 
  data$obsp[[paste0("similarity_", e)]] <- embedding.similarity(data$obsm[[e]], knn = 5)

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
```
</details>

<details>
<summary>BEYOND - Cellular trajectories & dynamics</summary>

To uncover trajectories in the cellular landscape we used two different trajectory inference tools: Palantir and VIA. The code for inferring the trajectories, as well as the parameters used in the manuscript can be found in file *4. BEYOND/2.trajectories.and.dynamics.in.landscape.R* lines 10-31. To load these trajectories from the supplementary tables use the following:

```R
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

palantir <- read.xlsx(SUPP$table.5, "Palantir Trajectories") %>% column_to_rownames("individualID") %>% `[`(rownames(data), )
via <- read.xlsx(SUPP$table.5, "VIA Trajectories") %>% column_to_rownames("individualID") %>% `[`(rownames(data), )

data$uns$trajectories <- list(
  palantir  = list(
    pseudotime = setNames(palantir$pseudotime, rownames(palantir)),
    branch.probs = palantir[, c("prAD", "ABA")],
    terminals = palantir %>% filter(!is.na(terminal)) %>% 
      select(pseudotime, traj=terminal) %>% rownames_to_column("terminal") %>% column_to_rownames("traj"),
    user.root = palantir[!is.na(palantir$root), "root"]
  ),
  
  via = list(
    pseudotime = setNames(via$pseudotime, rownames(via)),
    branch.probs = via[, c("prAD.like", "ABA.like", "Trajectory.3")] %>% `colnames<-`(gsub(".like", "", colnames(.))),
    terminals = via %>% filter(!is.na(terminal)) %>% 
      select(pseudotime, traj=terminal) %>% rownames_to_column("terminal") %>% column_to_rownames("traj"),
    user.root = via[!is.na(via$root), "root"]
  )
)

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(palantir, via)
```

To load fitted endophenotype- and subpopulation dynamics from the supplementary tables use the following:

```R
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

data$uns$trajectories$palantir$dynamics = list(
  fitted.vals = read.xlsx(SUPP$table.5, "Dynamics fitted values") %>% 
    filter(!feature %in% c("C1", "C2", "C3", "C4", "C5", "C1.1", "C1.2", "C2.3", "C2.2", "C2.1")) %>% 
    mutate(fit_sd=NA, se.fit_sd=NA) %>%
    select(y=value, x=pseudotime, X.weights.=weight, fit, se.fit, fit_sd, se.fit_sd, feature, trajectory),
  
  pred.vals = read.xlsx(SUPP$table.5, "Dynamics predicted values") %>% 
    filter(!feature %in% c("C1", "C2", "C3", "C4", "C5", "C1.1", "C1.2", "C2.3", "C2.2", "C2.1")) %>% 
    mutate(fit_sd=NA, se.fit_sd=NA) %>%
    select(x=pseudotime, fit, se.fit, fit_sd, se.fit_sd, feature, trajectory)
)

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
```
</details>

<details>
<summary>BEYOND - Communities</summary>
Subpopulations were partitioned into cellular communities using: (1) sub-population co-occurrences across individuals, and (2) similarities in their dynamics across both trajectories. The code for defining the communities, and computing their dynamics and endophenotype associations can be found in file *4. BEYOND/3.construct.cellular.communities.R*. To load these data from the supplementary tables use the following:

```R
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

# Append sub-population to community partition
data$var <- read.xlsx(SUPP$table.5, "Subpopulation communities") %>% 
  column_to_rownames("subpopulation") %>% `[`(colnames(data), )

# Append participants' community proportions
data$obsm$communities <- read.xlsx(SUPP$table.5, "Participant comm. prop.") %>% 
  column_to_rownames("individualID") %>% `[`(rownames(data), c("C1", "C2", "C3"))
data$obsm$sub.communities <- read.xlsx(SUPP$table.5, "Participant comm. prop.") %>% 
    column_to_rownames("individualID") %>% `[`(rownames(data),) %>% select(-C1, -C2, -C3)

# Append community dynamics
data$uns$communities$dynamics = list(
  fitted.vals = read.xlsx(SUPP$table.5, "Dynamics fitted values") %>% 
    filter(!feature %in% c("C1", "C2", "C3", "C4", "C5", "C1.1", "C1.2", "C2.3", "C2.2", "C2.1")) %>% 
    mutate(fit_sd=NA, se.fit_sd=NA) %>%
    select(y=value, x=pseudotime, X.weights.=weight, fit, se.fit, fit_sd, se.fit_sd, feature, trajectory),
  
  pred.vals = read.xlsx(SUPP$table.5, "Dynamics predicted values") %>% 
    filter(!feature %in% c("C1", "C2", "C3", "C4", "C5", "C1.1", "C1.2", "C2.3", "C2.2", "C2.1")) %>% 
    mutate(fit_sd=NA, se.fit_sd=NA) %>%
    select(x=pseudotime, fit, se.fit, fit_sd, se.fit_sd, feature, trajectory)
)

# Append endophenotype associations
data$uns$communities$trait.association <- read.xlsx(SUPP$table.5, "Community endophe. assoc.")

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(data)
```
</details>

<details>
<summary>Plotting manuscript figures</summary>
To plot the figures seen in the manuscript load the `anndata` object and use the different plotting util functions.

-   The sourced utils file expects the *2. Cell-type analysis/data/DLPFC.Green.atlas.h5* file to be available. This file contains atlas-level information such as subpopulation differential expressed genes, pathways, UMAP coordinates and pseudo-bulk gene expression calculated for each of the different cell-type Seurat objects. To replicate this data structure from the Synapse available Seurat objects follow *2. Cell-type analysis/2.aggregate.cell.type.data.R*

- Alternately, comment out atlas related plot settings in lines 125-196

```R
source("5. Manuscript code/utils.R")
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

landscape <- plot_grid(
  ggdraw() + draw_label("Trajectories in Cellular Landscape"),
  plot.landscape(data$uns$trajectories$palantir$branch.probs %>% 
                              py_to_r %>% mutate(diff=prAD-ABA) %>% dplyr::select(diff),
               smoothened = F,
               size = 2,
               cols = green2purple.less.white,
               legend.position = "right", show.feature.name = FALSE),
    ncol=1,
    rel_heights = c(.1, 1))

dynamics <- plot_grid(
  ggdraw() + draw_label("Subpopulation Dynamics Along Trajectories"),
  plot_grid(
      plot.dynamics(features = c("Oli.7","Ast.10"), 
                    cols = list(Oli.7="olivedrab4", Ast.10="darkorchid4"), 
                    label = T, legend.position = "none", include.points = TRUE, overlap.pseudotime = .1) +
        labs(x="Pseudotime", y="Proportion") + 
        theme(strip.background = element_rect(fill="lightgrey")),
      plot.dynamics(features = c("Ast.5","OPC.3"), 
                    cols = list(Ast.5="darkorchid4",OPC.3="springgreen4"), 
                    label = T, legend.position = "none", include.points = TRUE, overlap.pseudotime = .1) +
        labs(x="Pseudotime", y="Proportion") + 
        theme(strip.background = element_rect(fill="lightgrey")),
      ncol=1),
    ncol=1,
    rel_heights = c(.1, 1))

associations <- plot_grid(
  ggdraw() + draw_label("Endophenotype Associations"),
  plot.trait.associations(py_to_r(data$uns$trait.analysis$snuc), 
                        params = names(AD.traits),
                        column_title="Discovery (snRNA-seq)",
                        column_labels = c("A", "T", "C"),
                        column_names_rot = 0,
                        column_names_centered = T,
                        row_names_side = "left",
                        use_raster=T,
                        border=T,
                        raster_quality = 10) %>% draw %>% grid.grabExpr(),
    ncol=1,
    rel_heights = c(.1, 1))

landscape + dynamics + associations
```

<kdb><img src="https://github.com/naomihabiblab/BEYOND_DLPFC/assets/43610945/028a2ea5-71b2-4705-949c-e0fd784b1937"/><kdb/>
</details>

---

> [!IMPORTANT]
> The code available in this repository is the code used during the analysis. This includes the use of non-publicly available information such as participant identifiers `projid`. In the published supplementary tables and data participants are identified using `individualID`. It is therefore possible that some lines of code will fail when cloning and running with the publicly available information. In such cases check if a `projid` data column is used and replace with `individualID`. In addition, when using participants' demographic/clinical information (by using the `load.metadata` function) be sure to specify loading the supplementary table and providing the correct path to the file.

---

## Citation
- Green, G.S., Fujita, M., Yang, HS. et al. Cellular communities reveal trajectories of brain ageing and Alzheimer’s disease. Nature (2024). [https://doi.org/10.1038/s41586-024-07871-6](https://doi.org/10.1038/s41586-024-07871-6)

