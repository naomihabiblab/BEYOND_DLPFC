source("2. Cell-type analysis/load.code.env.R")


####################################################################################################################
##                                           #  RNAscope validations  #                                           ##
####################################################################################################################
validations <- list(RNAscope = list())

# -------------------------------------------------------- #
# Pool RNAscope quantification from files                  #
# -------------------------------------------------------- #
files <- list(batch2 = "Iba1_DAPI_FINAL.csv") # batch1 = "Iba1_DAPI.csv" - batch had high background noise and therefore was not used
names <- c(Children_Cy3_F_Count = "TPRG1", Children_FilterObjects_cy3_F_Count="TPRG1", 
           Children_Cy5_F_Count = "MRC1",  Children_FilterObjects_cy5_F_Count="MRC1", 
           Children_Cy7_F_Count = "CPM",   Children_FilterObjects_cy7_F_Count="CPM")
path  <- "3. Other analyses/data/RNAscope quantification/"

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

o <- LoadH5Seurat("2. Cell-type analysis/microglia/data/microglia.h5Seurat", assays=list(SCT=c("counts")), reductions=F, graphs=F, neighbors=F, misc=F, verboth=F)
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


saveRDS(validations, "3. Other analyses/data/RNAscope.rds")
