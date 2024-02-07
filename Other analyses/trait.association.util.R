library(dplyr)

#' Perform trait associations 
#'
#' Associate co-variates to traits will controlling for technical variables
#' 
#' @param traits dataframe of shape (n_samples, n_traits)
#' @param covariates dataframe of shape (n_samples, n_covariates)
#' @param controls dataframe of shape (n_samples, n_controls)
#'
#' @return dataframe with trait association results for every n_traits*n_covariates 
#' while correcting for multiple hypothesis testing within each of the given traits
#'
associate.traits <- function(traits, covariates, controls, p.adjust.method="BH") {
  df <- data.frame(traits, covariates, controls)
  
  df <- do.call(rbind, lapply(colnames(traits), function(trait)
    do.call(rbind, lapply(colnames(covariates), function(covariate) {
      # Set controlled covariates
      control <- colnames(controls)
      control <- case_when(trait == "cogng_demog_slope" ~ paste(setdiff(control, c("msex", "age_death")), collapse = " + "), 
                           T ~ paste(control, collapse = " + "))
      
      formula <- stringr::str_interp("${trait} ~ ${covariate} + ${control}")
      m <- summary(lm(formula, df[!is.na(df[,trait]), ]))
      data.frame(trait = trait, 
                 covariate = covariate,
                 beta  = m$coefficients[covariate, 1],
                 se = m$coefficients[covariate, 2],
                 tstat = m$coefficients[covariate, 3],
                 pval  = m$coefficients[covariate, 4],
                 r.sq  = m$adj.r.squared,
                 n     =  sum(!is.na(df[,trait])),
                 formula = formula)
    }))))
  return(df %>% dplyr::group_by(trait) %>% 
           dplyr::mutate(adj.pval = p.adjust(pval, method=p.adjust.method),
                         sig = cut(adj.pval, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))) %>% 
           dplyr::ungroup())
}