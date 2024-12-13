library(dplyr)

#' 
#' Load ROSMAP participant metadata
#'
#' The following function loads participants' demographic, clinical and pathological information, of either the snRNA-seq discovery 
#' or bulk replication cohorts. During the original analysis we used RADC obtained participant information which contains participants
# of both cohorts and uses the `projid` identifier. For use with the manuscript's publicly available Supplemetary Table 1 please specify
# which of the two cohorts to load. In this case the participants' identifier is the `individualID`. 
# To obtain a mapping between the `projid` and `participantID` see the ROSMAP clinical information file available on Synapse.
load.metadata <- function(
    path = "1. Library preprocessing/data/ROSMAP.participants.metadata_04-25-2022.csv",
    from.supplementary.table = FALSE,
    sheet = c("Discovery cohort", "Replication cohort")) {
  
  if(from.supplementary.table)
    df <- openxlsx::read.xlsx(path, sheet) %>% tibble::column_to_rownames("individualID")
  else {
    # Loading from ROSMAP clinical information contains more information than what is available in the supplementary table
    df <- read.csv(path, header = TRUE) %>%
      mutate(projid = as.character(projid), 
             
             # Semi-quantitative pathology measures
             dlbdx = factor(dlbdx, levels=0:3, ordered = T),
             tdp_stage4 = factor(tdp_stage4, levels=0:3, ordered = T),
             
             # Quantitative pathology measures
             across(matches("amyloid|tangles|nft|plaq"), sqrt, .names="sqrt.{.col}"),
             
             # APOE genotype
             apoe_genotype = factor(apoe_genotype),
             apoe_2 = as.numeric(apoe_genotype %in% c(22, 23)),
             apoe_3 = as.numeric(apoe_genotype %in% c(33)),
             apoe_4 = as.numeric(apoe_genotype %in% c(34, 44))) %>%
      tibble::column_to_rownames("projid")
  }
  
  df <- df %>% dplyr::mutate(
    sex = factor(plyr::mapvalues(msex, c(F,T), c("Female", "Male"))),
    
    # Cognition state
    cogdx = factor(cogdx, levels = 1:6),
    cogdx_ad = as.numeric(recode(cogdx, "1"="1","2"="2", "3"=NA_character_, "4"="3", "5"=NA_character_,"6"=NA_character_)),
    cogdx_grouped = as.numeric(recode(cogdx, "3"="2", "4" = "3","5"="3")),
    
    # Semi-quantitative pathology measures
    ceradsc = factor(ceradsc, levels = 4:1, ordered=T),
    braaksc = factor(case_when(as.numeric(braaksc) == 0 ~ NA_integer_, TRUE~as.numeric(braaksc)), levels = 1:6, ordered = T),
    niareagansc = factor(niareagansc, levels = 4:1, ordered=T),
    pAD = factor(plyr::mapvalues(as.numeric(as.character(niareagansc)), from = c(1,2,3,4), to=c(1,1,0,0)), labels = c(0,1)),
    
    # Generating textual pathology annotations
    cerad.txt = factor(ceradsc, levels = 4:1, labels = c("No AD", "Possible", "Probable", "Definite")),
    braak.txt = factor(braaksc, levels = 1:6, labels = c("I", "II", "III", "IV", "V", "VI")),
    braak.grouped.txt = factor(braaksc, levels = 0:6, labels = c("0","I", "II", "III", "IV", "V+VI", "V+VI")),
    cogdx.grouped.txt = factor(cogdx, levels = c(1,2,3,4,5,6), labels = c("NCI", "MCI", "MCI+Other","AD","AD+Other","Other\nDementia")),
    cdx.txt = factor(cogdx, levels=c(1,2,3,4,5,6), labels=c("No Cognitive\nImpairment", "Mild\nCognitive Impairment", 
                                                            "Mild Cognitive\nImpairment\nand Another\nCause of CI", "AD", 
                                                            "AD\nand Another\nCause of CI", "Other Dementia")),
    niareagen.txt = factor(niareagansc, levels=4:1, labels = c("No AD", "Low", "Intermediate", "High")))
  
  return(df)
}
