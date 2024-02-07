source("BEYOND/utils.R")


####################################################################################################################
##                                     #  RNAscope validations & related analysis  #                              ##
####################################################################################################################
validations <- list(RNAscope = list())


# -------------------------------------------------------- #
# Pool RNAscope quantification from files                  #
# -------------------------------------------------------- #
files <- list(batch2 = "Iba1_DAPI_FINAL.csv") # batch1 = "Iba1_DAPI.csv" - batch had high background noise and therefore was not used
names <- c(Children_Cy3_F_Count = "TPRG1", Children_FilterObjects_cy3_F_Count="TPRG1", 
           Children_Cy5_F_Count = "MRC1",  Children_FilterObjects_cy5_F_Count="MRC1", 
           Children_Cy7_F_Count = "CPM",   Children_FilterObjects_cy7_F_Count="CPM")
path  <- "Other analyses/data/RNAscope quantification/"

df <- lapply(names(files), function(batch) {
  lapply(list.files(file.path(path, batch), pattern = files[[batch]], recursive = TRUE, full.names = TRUE), function(f) {
    # Load quantification file of image
    parts <- stringr::str_match(f,".*\\/(.*)\\/(\\d+)")[,2:3]
    df <- read.csv(f) %>% mutate(ImageNumber = parts[[2]], 
                                 sample = gsub("_ctrl|-|_.*| .*", "", parts[[1]]))
    if(batch == "batch2") {
      morph <- read.csv(file.path(dirname(f), "MyExpt_Iba1_DAPI.csv")) %>% 
        dplyr::select(ObjectNumber, AreaShape_Compactness, AreaShape_Eccentricity, Parent_IBA1_final)
      df <- merge(df, morph, by.x = "Parent_Iba1_DAPI", by.y = "ObjectNumber", all.x = TRUE)
    }
    
    df %>% mutate(batch=batch,
                  AD = factor(case_when(sample %in% c("T5961","T5993","MP05187","MP04112", "MP05187","MP1141","MP08229","OC0746",
                                                      "OC1139","OC1850") ~ "No AD",
                                        sample %in% c("T5915","MP0671","OC1731") ~ "MCI",
                                        sample %in% c("T5939","T5918","T5899","T5869","T5818","T5972","T5668","MP1060","MP117",
                                                      "MP11189","MP1242","MP13101","MP14134","MP1791","OC1253") ~ "AD"),
                              levels = c("No AD", "MCI", "AD"))) %>%
      plyr::rename(names, warn_missing = FALSE) %>%
      dplyr::select(batch, sample, AD, ImageNumber, ObjectNumber, 
                    iba1=Parent_IBA1_final, TPRG1, MRC1, CPM,
                    compactness=AreaShape_Compactness, eccentricity=AreaShape_Eccentricity) %>%
      filter(!( is.na(TPRG1) | is.na(MRC1) | is.na(CPM)))
  }) %>% plyr::rbind.fill(.)
}) %>% plyr::rbind.fill(.)
validations$RNAscope$df <- df

rm(files, names, path)


# -------------------------------------------------------- #
# Obtain snRNAseq expression values of validation genes    #
# -------------------------------------------------------- #
library(Seurat)
library(SeuratDisk)

o <- LoadH5Seurat("Cell-type analysis/microglia/data/microglia.h5Seurat", assays=list(SCT=c("counts")), reductions=F, graphs=F, neighbors=F, misc=F, verboth=F)
snuc.exp <- FetchData(o, vars = c("state","TPRG1","CPM","MRC1")) %>%
  filter(TPRG1+CPM+MRC1 != 0) %>%
  mutate(state = case_when(state %in% c("Mic.12","Mic.13","Macrophages") ~ state, T~"Other"))

validations$RNAscope$snuc.exp <- snuc.exp
rm(o, snuc.exp)


# -------------------------------------------------------- #
# Estimate RNAscope state prevalences in donors            #
# -------------------------------------------------------- #
validations$RNAscope$df <- df <- df %>% 
  mutate(state = case_when(TPRG1>5 &  MRC1<10 ~ "Mic.13",
                           TPRG1<2 & CPM>5 &  MRC1<5 ~ "Mic.12",
                           MRC1>=5 ~ "Macrophages",
                           MRC1<3 & TPRG1 < 2 & CPM < 2 ~ "Other",
                           .default = "none")) %>%
  filter(!is.na(state)) %>% 
  mutate(v=1) %>% 
  tidyr::spread(state, v, fill = 0)

validations$RNAscope$predicted.proportions <- df %>%
  group_by(sample, AD, batch) %>%
  summarise(across(c(Mic.12, Mic.13, Macrophages, Other), ~sum(.)/n()), ncells = n(),
            across(c(TPRG1, MRC1, CPM), mean),
            .groups = "drop") %>%
  column_to_rownames("sample")



# -------------------------------------------------------- #
# Append and summarise pTau quantification for donors      #
# -------------------------------------------------------- #
df <- list.dirs("Other analyses/data/pTAU validations/batch2") %>% # batch1 had high background noise and therefore was not used
  Filter(function(f) file.exists(file.path(f, "MyExpt_Image.csv")) ,.) %>%
  lapply(., function(f) {
    batch = basename(dirname(f))
    sample = basename(f) %>% gsub(".*/", "", .) %>% gsub("_ctrl|-|_.*| .*", "",.)
    
    read.csv(file.path(f, "MyExpt_Image.csv")) %>% 
      dplyr::select(ImageNumber, 
                    pTau_merged = AreaOccupied_AreaOccupied_pTau_merged,
                    pTau_tangles_plaques = AreaOccupied_AreaOccupied_PTau_Tangles_and_plaques,
                    pTau_filaments = AreaOccupied_AreaOccupied_pTau_other_than_tangles_and_plaques) %>%
      mutate(sample=sample, batch=batch)
  }) %>%
  do.call(rbind, .)
df$AD <- validations$RNAscope$predicted.proportions[df$sample,]$AD
validations$pTau = list(raw=df)

validations$pTau$summarised <- df %>% group_by(batch,sample, AD) %>% 
  summarise(across(c(pTau_merged, pTau_tangles_plaques, pTau_filaments), list(mean=mean, median=median))) %>% 
  data.frame



# -------------------------------------------------------- #
# Subpopulation-proportion morphology analysis             #
# -------------------------------------------------------- #
validations$RNAscope$morpholoy.association <-
  lapply(c("compactness", "eccentricity"), function(y)
    sapply(c("Mic.12","Mic.13","Macrophages","Other"), function(s)
      summary(lm(paste0(y, "~", s), validations$RNAscope$df))[["coefficients"]][2,]) %>% t %>%
      `colnames<-`(c("beta", "se", "tstat", "pval")) %>% as.data.frame %>%
      mutate(adj.pval = p.adjust(pval, method = "BH"),
             sig = cut(adj.pval, c(-Inf, .0001, .001,.01,.05,1), c("****","***","**","*","")),
             trait = y) %>% 
      rownames_to_column("state")) %>%
  do.call(rbind, .) %>%
  dplyr::select(trait, state, everything())
  

# ------------------------------------------------------------- #
# Subpopulation-proportion morphology analysis - ROSMAP cohort  #
# ------------------------------------------------------------- #
data <- anndata::read_h5ad("BEYOND/data/BEYOND.DLPFC.h5ad")

others = paste0("Mic.", c(1:11,14:16) )
df <- data.frame(data$X[,data$var$grouping.by == "Microglia"],
                 data$obsm$meta.data[,c("mglia3_mf","mglia123_mf","pmi","sex","age_death")]) %>%
  filter(!is.na(mglia3_mf)) %>%
  mutate(Other = rowSums(dplyr::select(., all_of(others))),
         across(c(Mic.12,Mic.13,Macrophages,Other), ~sqrt(.)),
         pam = sqrt(mglia3_mf / mglia123_mf)) %>%
  dplyr::select(-all_of(c(others, "Monocytes")))

# Independent regressions of pam scores on states
res <- lapply(c("Mic.12","Mic.13","Macrophages","Other"), 
              function(s) summary(lm(paste0("pam~", s,"+pmi+age_death+sex"), df))$coefficients %>% 
                data.frame %>% `[`(2,)) %>% do.call(rbind, .)  


# Mic.12/13 associations while controlling for one another
controls <- c("pmi","age_death","sex", "Mic.12", "Mic.13")
res.ctrl <- lapply(c("Mic.12","Mic.13"), function(state) {
  .df <- df %>% mutate(covariate = get(state))
  fit <- lm(paste0("pam~", state,"+", paste(setdiff(controls, state), collapse = "+")), .df)
  summary(fit)$coefficients %>% data.frame() %>% `[`(2,) %>% `rownames<-`(paste0("ctrl.", state))
}) %>% do.call(rbind, .)

# Multiple hypothesis testing correction
res <- rbind(res, res.ctrl) %>% 
  mutate(adj.pval = p.adjust(Pr...t.., method = "BH"),
         sig = cut(adj.pval, c(-Inf, .0001, .001, .01, .05, 1), c("****","***","**","*","")))

validations$ROSMAP.PAM.Association <- res


saveRDS(validations, "Other analyses/data/RNAscope.rds")
rm(validations, res, controls, res.ctrl, others,data, df)


