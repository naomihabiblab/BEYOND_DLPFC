library(dplyr)


load.metadata <- function(path = "1. Library preprocessing/data/ROSMAP.participants.metadata_04-25-2022.csv") {
  read.csv(path, header = T) %>%
    dplyr::mutate(
      projid = as.character(projid),
      sex = factor(plyr::mapvalues(msex, c(F,T), c("Female", "Male"))),
      
      # Cognition state
      cogdx = factor(cogdx, levels = 1:6),
      
      cogdx_ad = as.numeric(recode(cogdx, "1"="1","2"="2", "3"=NA_character_, "4"="3", "5"=NA_character_,"6"=NA_character_)),
      cogdx_grouped = as.numeric(recode(cogdx, "3"="2", "4" = "3","5"="3")),
      
      # Semi-quantitative pathology measures
      ceradsc = factor(ceradsc, levels = 4:1, ordered=T),
      braaksc = factor(case_when(braaksc == 0 ~ NA_integer_, T~braaksc), levels = 1:6, ordered = T),
      niareagansc = factor(niareagansc, levels = 4:1, ordered=T),
      pAD = factor(plyr::mapvalues(as.numeric(as.character(niareagansc)), from = c(1,2,3,4), to=c(1,1,0,0)), labels = c(0,1)),
      dlbdx = factor(dlbdx, levels=0:3, ordered = T),
      tdp_stage4 = factor(tdp_stage4, levels=0:3, ordered = T),
      
      # Quantitative pathology measures
      across(matches("amyloid|tangles|nft|plaq"), sqrt, .names="sqrt.{.col}"),
      
      apoe_genotype = factor(apoe_genotype),
      apoe_2 = as.numeric(apoe_genotype %in% c(22, 23)),
      apoe_3 = as.numeric(apoe_genotype %in% c(33)),
      apoe_4 = as.numeric(apoe_genotype %in% c(34, 44))) %>%
    dplyr::mutate(cerad.txt = factor(ceradsc, levels = 4:1, labels = c("No AD", "Possible", "Probable", "Definite")),
                  braak.txt = factor(braaksc, levels = 1:6, labels = c("I", "II", "III", "IV", "V", "VI")),
                  braak.grouped.txt = factor(braaksc, levels = 0:6, labels = c("0","I", "II", "III", "IV", "V+VI", "V+VI")),
                  cogdx.grouped.txt = factor(cogdx, levels = c(1,2,3,4,5,6), labels = c("NCI", "MCI", "MCI+Other","AD","AD+Other","Other\nDementia")),
                  cdx.txt = factor(cogdx, levels=c(1,2,3,4,5,6), labels=c("No Cognitive\nImpairment", "Mild\nCognitive Impairment", 
                                                                          "Mild Cognitive\nImpairment\nand Another\nCause of CI", "AD", 
                                                                          "AD\nand Another\nCause of CI", "Other Dementia")),
                  niareagen.txt = factor(niareagansc, levels=4:1, labels = c("No AD", "Low", "Intermediate", "High"))) %>%
    tibble::column_to_rownames("projid")
}