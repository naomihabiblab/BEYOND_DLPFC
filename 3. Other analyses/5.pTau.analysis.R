source("Cell-type analysis/load.code.env.R")


####################################################################################################################
##                                        #  RNAscope pTAU validations  #                                         ##
####################################################################################################################
validations <- readRDS("Other analyses/data/RNAscope.rds")


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

saveRDS(validations, "Other analyses/data/RNAscope.rds")
