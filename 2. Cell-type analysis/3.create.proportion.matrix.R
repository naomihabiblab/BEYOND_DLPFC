source("2. Cell-type analysis/load.code.env.R")


####################################################################################################################
##                                #  Creating subpopulation proportion matrix  #                                  ##
####################################################################################################################

# Load cell-type annotations previously pulled from the different cell-type Seurat objects
df  <- h5read(file = "2. Cell-type analysis/DLPFC.Green.atlas.h5", name = "annotations") %>% filter(state != "NA")

# Retain cell-type annotations only of participants passing QCs (Methods). This step reduces the number
# of participants from the general 465 participants whose cells form the cell atlas to 437.
qcs <- read.csv("2. Cell-type analysis/data/donors.qc.csv") %>% filter(Final_QC == "Pass")
qcs <- qcs[qcs$projid %in% unique(df$projid),] %>% column_to_rownames("projid")
df  <- df[df$projid %in% rownames(qcs), ]

# Participants' nuclei were sequenced in two technical replicates (A and B) 
donor.batches <- df %>% group_by(projid) %>% summarise(batch = paste(unique(batch), collapse=", ")) %>% column_to_rownames("projid")
main.batch <- df %>% mutate(batch = gsub("-[A|B]$", "", batch)) %>% 
  count(projid, batch) %>% 
  group_by(projid) %>% 
  slice_max(order_by=n, n=1) %>% 
  column_to_rownames("projid") %>% 
  dplyr::select(batch)


# -------------------------------------------------------- #
# Compute subpopulation proportions                        #
# -------------------------------------------------------- #

# Calculate participant-wise subpopulation-proportion (within cell-type)
df <- df %>% 
  filter(grouping.by != "Immune") %>% # remove NK cells and T-cells for which we do not have sufficient cells
  count(class, grouping.by, cell.type, projid, state) %>%
  group_by(grouping.by, projid) %>%
  mutate(prevalence=n/sum(n)) %>%
  ungroup()
gc()


# Creating proportion- and count matrices of participants' subpopulations
proportions <- dcast(df, projid~state, value.var = "prevalence", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid")
counts <- dcast(df, projid~state, value.var = "n", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid")

saveRDS(proportions, "2. Cell-type analysis/data/subpopulation.proportion.matrix.rds")
saveRDS(counts, "2. Cell-type analysis/data/subpopulation.counts.matrix.rds")

# Create base AnnData representation of subpopulation proportion matrix user for downstream trait-associations and BEYOND analyses
ids <- rownames(proportions)
data <- AnnData(
  # Subpopulation proportions: participants (rows) over subpopulations (columns)
  X = proportions,
  
  # Counts and sqrt(prop) of subpopulations
  layers = list(counts = counts[ids, colnames(proportions)],
                sqrt.prev = sqrt(proportions)),
  
  # General information about the different subpopulations
  var = df %>% dplyr::select(class, grouping.by, cell.type, state) %>% unique %>% column_to_rownames("state") %>% `[`(colnames(proportions),),
  
  # General information about the participants
  obs = list(batches = donor.batches[ids,],
             main.batch = main.batch[ids,]),
  
  # Participant-wise QCs
  obsm = list(QCs = qcs[ids,] %>% `rownames<-`(rownames(.) %>% as.character()))
)

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(qcs, ids, donor.batches, main.batch)



# -------------------------------------------------------- #
# Append participant metadata                              #
# -------------------------------------------------------- #
source("1. Library preprocessing/utils/ROSMAP.metadata.R")

data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
meta.data <- load.metadata()
data$obsm$meta.data <- meta.data[data$obs_names,]
anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(meta.data)


# -------------------------------------------------------- #
# Cell type level of aggregation                           #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
cell.type.df <- df %>% group_by(projid, class, cell.type) %>%
  dplyr::summarise(n = sum(n), .groups = "drop") %>%
  dplyr::group_by(projid) %>%
  dplyr::mutate(prevalence = n/sum(n)) %>%
  dplyr::group_by(projid, class) %>%
  dplyr::mutate(within.class.prevalence = n/sum(n))

ct.counts <- (dcast(cell.type.df, projid~cell.type, value.var = "n", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid"))[data$obs_names, ]
ct.prev   <- (dcast(cell.type.df, projid~cell.type, value.var = "prevalence", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid"))[data$obs_names, ]
ct.c.prev <- (dcast(cell.type.df, projid~cell.type, value.var = "within.class.prevalence", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid"))[data$obs_names, ]

data$uns$cell.types <- list(counts = ct.counts, 
                            prev = ct.prev, 
                            sqrt.prev = sqrt(ct.prev), 
                            wc.prev = ct.c.prev, 
                            sqrt.wc.prev = sqrt(ct.c.prev))
rm(cell.type.df, ct.counts, ct.prev, ct.c.prev)

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(df, data, tables)

