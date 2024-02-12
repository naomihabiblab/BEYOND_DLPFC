source("Cell-type analysis/load.code.env.R")


####################################################################################################################
##                                    #  RNAscope & ROSMAP morphology validations  #                              ##
####################################################################################################################
validations <- readRDS("Other analyses/data/RNAscope.rds")


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
data <- anndata::read_h5ad("Cell-type analysis/data/subpopulation.proportions.h5ad")

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
