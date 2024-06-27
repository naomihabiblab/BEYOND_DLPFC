source("2. Cell-type analysis/utils/load.code.env.R")
library(mediation)
library(lavaan)
library(semPlot)



####################################################################################################################
##                                                 #  Mediation Modeling  #                                       ##
####################################################################################################################
# (ROSMAP data download 04/25/2022)
source("1. Library preprocessing/utils/ROSMAP.metadata.R")

data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
df.full <- rbind(data.frame(cohort = "discovery", sqrt(data$X[,c("Mic.12","Mic.13","Ast.10","Oli.7")])),
                 data.frame(cohort = "replication", data$uns$celmod$avg.predicted.prop$validation %>% 
                              py_to_r %>% `[`(, c("Mic.12","Mic.13","Ast.10","Oli.7")))) %>% 
  # Add snRNA-seq library quality measures (for discovery participants only)
  merge(., data$obsm$QCs[,c("Total_Genes_Detected","Estimated_Number_of_Cells")], all.x=TRUE, by.x="row.names", by.y="row.names") %>%
  # Add bulk library quality measure (RIN, used only for replication participants)
  merge(., read.csv("3. Other analyses/data/ROSMAP.bulk.RIN.values_08222023.csv", row.names = 1), all.x=TRUE, by.x="Row.names", by.y="row.names") %>%
  # Add phenotypic data of participants from both cohorts
  merge(., load.metadata(id.by.individualID = FALSE) %>% 
          dplyr::select(cogng_demog_slope, sqrt.amyloid_mf, sqrt.tangles_mf, age_death, sex, pmi, apoe_genotype) %>%
          mutate(msex = sex == "Male",
                 apoe4n = stringr::str_count(apoe_genotype, "4")) %>%
          dplyr::select(-sex, -apoe_genotype), 
        all.x=TRUE, by.x="Row.names", by.y="row.names") %>%
  column_to_rownames("Row.names")



# -------------------------------------------------------- #
# APOE e4 associations                                     #
# -------------------------------------------------------- #
df <- subset(df.full, cohort == "discovery")

## APOE e4 association:
m <- lm(Mic.12 ~ apoe4n+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m)
# apoe4n                     5.399e-03  1.025e-02   0.527 0.598622    
# n=435
m <-lm(Mic.13 ~ apoe4n+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m); confint(m)
# apoe4n                     4.199e-02  1.136e-02   3.697 0.000247 ***
# apoe4n                     1.966515e-02  6.432048e-02
# n=435
m <- lm(Ast.10 ~ apoe4n+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m); confint(m)
# apoe4n                     2.390e-02  1.825e-02   1.310    0.191
# n=435
m <- lm(Oli.7 ~ apoe4n+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m); confint(m)
# apoe4n                     2.041e-02  1.501e-02   1.360   0.1746  
# n=435


# -------------------------------------------------------- #
# Mic.13 placement in cascade                              #
# -------------------------------------------------------- #
df <- subset(df.full, cohort == "discovery")
m1 <- lm(sqrt.amyloid_mf ~ apoe4n+msex+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m1) # n=432
m2 <- lm(Mic.13 ~ sqrt.amyloid_mf+apoe4n+msex+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m2) # n=432

set.seed(7725)
m_med01_e4AbMic <- mediate(m1, m2, sims=10000, treat="apoe4n", mediator="sqrt.amyloid_mf", boot=TRUE); summary(m_med01_e4AbMic)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME             0.0319       0.0219         0.04  <2e-16 ***
# ADE              0.0104      -0.0130         0.03  0.3826    
# Total Effect     0.0422       0.0199         0.07  0.0002 ***
# Prop. Mediated   0.7545       0.4454         1.60  0.0002 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 432 
# Simulations: 10000 


df <- subset(df.full, cohort == "discovery" & !is.na(sqrt.tangles_mf))
m1 <- lm(Mic.13 ~ sqrt.amyloid_mf+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m1)
# n=432
m2 <- lm(sqrt.tangles_mf ~ Mic.13+sqrt.amyloid_mf+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m2)
# n=432
set.seed(1776)
m_med03_AbMicTau <- mediate(m1, m2, sims=10000, boot=TRUE, treat="sqrt.amyloid_mf", mediator="Mic.13")
summary(m_med03_AbMicTau)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME             0.0996       0.0493         0.16  <2e-16 ***
# ADE              0.2139       0.1381         0.29  <2e-16 ***
# Total Effect     0.3136       0.2456         0.39  <2e-16 ***
# Prop. Mediated   0.3178       0.1645         0.51  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 432
# Simulations: 10000 


df <- subset(df.full, cohort == "discovery" & !is.na(cogng_demog_slope))
m1 <- lm(sqrt.tangles_mf ~ Mic.13+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m1) # n=413
m2 <- lm(cogng_demog_slope ~ sqrt.tangles_mf+Mic.13+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m2) # n=413
set.seed(2022)
m_med02_MicTauCog <- mediate(m1, m2, sims=10000, boot=TRUE, treat="Mic.13", mediator="sqrt.tangles_mf")
summary(m_med02_MicTauCog)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            -0.1347      -0.1879        -0.09  <2e-16 ***
# ADE             -0.0835      -0.1671        -0.01   0.036 *  
# Total Effect    -0.2182      -0.3074        -0.13  <2e-16 ***
# Prop. Mediated   0.6173       0.4061         0.97  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 413
# Simulations: 10000 
plot(m_med02_MicTauCog)
m_med02_MicTauCog$z1.p
# [1] 0.0362



# -------------------------------------------------------- #
# Mic.12 placement in cascade                              #
# -------------------------------------------------------- #
df <- subset(df.full, cohort == "discovery" & !is.na(tangles_mf_sqrt))
m1 <- lm(sqrt.amyloid_mf~Mic.12+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m1)
m2 <- lm(sqrt.tangles_mf~sqrt.amyloid_mf+Mic.12+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m2)

set.seed(77364)
m_med_m12AbTau<-mediate(m1, m2, sims=10000, boot=TRUE, treat="Mic.12", mediator="sqrt.amyloid_mf")
summary(m_med_m12AbTau)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME              0.951        0.595         1.38  <2e-16 ***
# ADE               2.097        0.900         3.39   2e-04 ***
# Total Effect      3.049        1.841         4.37  <2e-16 ***
# Prop. Mediated    0.312        0.185         0.55  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 432
# Simulations: 10000 


df <- subset(df.full, cohort == "discovery" & !is.na(tangles_mf_sqrt))
m1<-lm(Mic.13 ~ Mic.12+sqrt.amyloid_mf+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m1)
m2<-lm(sqrt.tangles_mf ~ Mic.13+Mic.12+sqrt.amyloid_mf+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m2)

set.seed(23421)
m_med_m12m13Tau<-mediate(m1, m2, sims=10000, boot=TRUE, treat="Mic.12", mediator="Mic.13")
summary(m_med_m12m13Tau)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME             0.9933       0.4413         1.66  0.0004 ***
# ADE              1.1040       0.0118         2.24  0.0474 *  
# Total Effect     2.0973       0.9026         3.36  0.0004 ***
# Prop. Mediated   0.4736       0.2230         0.99  0.0008 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 432
# Simulations: 10000 




# -------------------------------------------------------- #
# Ast.10 placement in cascade                              #
# -------------------------------------------------------- #
df <- subset(df.full, cohort == "discovery" & !is.na(sqrt.tangles_mf))
m <- lm(Ast.10 ~ sqrt.amyloid_mf+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m); confint(m)
# sqrt.amyloid_mf            1.611e-02  6.334e-03   2.543   0.0113 *
# sqrt.amyloid_mf            3.656569e-03 2.855798e-02
# n=432
m1<-lm(sqrt.tangles_mf~sqrt.amyloid_mf+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m1)
# sqrt.amyloid_mf            3.136e-01  3.770e-02   8.318 1.23e-15 ***
# n=432
m2<-lm(Ast.10~sqrt.tangles_mf+sqrt.amyloid_mf+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m2)
# sqrt.tangles_mf            3.578e-02  7.973e-03   4.488 9.27e-06 ***
# sqrt.amyloid_mf            4.887e-03  6.682e-03   0.731    0.465 
# n=432

set.seed(1919)
m_med04_AbTauAst<-mediate(m1, m2, sims=10000, boot=TRUE, treat="sqrt.amyloid_mf", mediator="sqrt.tangles_mf")
summary(m_med04_AbTauAst)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.01122      0.00354         0.02  0.0006 ***
# ADE             0.00489     -0.00764         0.02  0.4428    
# Total Effect    0.01611      0.00425         0.03  0.0056 ** 
# Prop. Mediated  0.69661      0.22702         2.26  0.0062 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 432
# Simulations: 10000 


df <- subset(df.full, cohort == "discovery" & !is.na(cogng_demog_slope))
m1 <- lm(Ast.10 ~ sqrt.tangles_mf+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m1) # n=413
m2 <- lm(cogng_demog_slope ~ Ast.10+sqrt.tangles_mf+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m2) # n=413

set.seed(2018)
m_med05_TauAstCog<-mediate(m1, m2, sims=10000, boot=TRUE, treat="sqrt.tangles_mf", mediator="Ast.10")
summary(m_med05_TauAstCog)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME           -0.00310     -0.00639         0.00  0.0012 ** 
# ADE            -0.03391     -0.04126        -0.03  <2e-16 ***
# Total Effect   -0.03701     -0.04429        -0.03  <2e-16 ***
# Prop. Mediated  0.08368      0.02351         0.17  0.0012 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 413 
# Simulations: 10000 

## Additional variance explained by Ast.10
df <- subset(df.full, cohort == "discovery" & !is.na(cogng_demog_slope))
m1 <- lm(cogng_demog_slope~sqrt.amyloid_mf+sqrt.tangles_mf+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m1)
# Multiple R-squared:  0.2398,	Adjusted R-squared:  0.2304 
# n=412
m2 <- lm(cogng_demog_slope~Ast.10+sqrt.amyloid_mf+sqrt.tangles_mf+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m2)
# Multiple R-squared:  0.2663,	Adjusted R-squared:  0.2554 
# n=412
# 0.2554 - 0.2304
# [1] 0.025


df <- subset(df.full, cohort == "discovery" & !is.na(sqrt.tangles_mf))
m1<-lm(sqrt.tangles_mf~Mic.13+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected,data=df);summary(m1)
m2<-lm(Ast.10~sqrt.tangles_mf+Mic.13+age_death+msex+pmi+Estimated_Number_of_Cells+Total_Genes_Detected,data=df);summary(m2)
set.seed(1347)

m_med_m13TauA10<-mediate(m1, m2, sims=10000, boot=TRUE, treat="Mic.13", mediator="sqrt.tangles_mf")
summary(m_med_m13TauA10)
# Causal Mediation Analysis
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# Estimate 95% CI Lower 95% CI Upper p-value   
# ACME             0.1124       0.0190         0.23   0.015 * 
# ADE              0.2383       0.0646         0.41   0.006 **
# Total Effect     0.3507       0.1829         0.53  <2e-16 ***
# Prop. Mediated   0.3205       0.0605         0.71   0.015 * 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 433
# Simulations: 10000



# -------------------------------------------------------- #
# Oli.7 placement in cascade                              #
# -------------------------------------------------------- #
df <- subset(df.full, cohort == "discovery")
m1 <- lm(Ast.10 ~ sqrt.tangles_mf+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df);summary(m1) # n=433
m2 <- lm(Oli.7 ~ Ast.10+sqrt.tangles_mf+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df);summary(m2) # n=433

set.seed(8876)
m_med06_TauAstOli<-mediate(m1, m2, sims=10000, boot=TRUE, treat="sqrt.tangles_mf", mediator="Ast.10")
summary(m_med06_TauAstOli)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.02784      0.01075         0.05  <2e-16 ***
# ADE             0.00234     -0.00598         0.01  0.6116    
# Total Effect    0.03018      0.00987         0.06  0.0002 ***
# Prop. Mediated  0.92247      0.71668         1.38  0.0002 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 433 
# Simulations: 10000 

df <- subset(df.full, cohort == "discovery")
m <- lm(cogng_demog_slope ~ Oli.7+Ast.10+pmi+Estimated_Number_of_Cells+Total_Genes_Detected, df); summary(m)
# Oli.7                -4.449e-03  5.890e-02  -0.076  0.93982   
# Ast.10               -1.319e-01  4.871e-02  -2.708  0.00705 **
# n=416
# Confounded by Ast.10




####################################################################################################################
##                                            #  Structural Equation Modeling #                                   ##
####################################################################################################################

df <- df.full %>% filter(cohort == "discovery") %>% filter(if_all(c(sqrt.amyloid_mf, sqrt.tangles_mf, cogng_demog_slope, Mic.13, pmi), ~ !is.na(.x))) # n=412

df$Mic.13_resid <- resid(lm(Mic.13 ~ age_death+msex+pmi+Total_Genes_Detected+Estimated_Number_of_Cells, df))
df$Ast.10_resid <- resid(lm(Ast.10 ~ age_death+msex+pmi+Total_Genes_Detected+Estimated_Number_of_Cells, df))
df$Ab_resid     <- resid(lm(sqrt.amyloid_mf ~ age_death+msex, df))
df$Tau_resid    <- resid(lm(sqrt.tangles_mf ~ age_death+msex, df))
df$Mic.12_resid <- resid(lm(Mic.12 ~ age_death+msex+pmi+Total_Genes_Detected+Estimated_Number_of_Cells, df))
df$Oli.7_resid  <- resid(lm(Oli.7 ~ age_death+msex+pmi+Total_Genes_Detected+Estimated_Number_of_Cells, df))

model1<-'
Ab_resid ~ Mic.12_resid
Mic.13_resid ~ Ab_resid + Mic.12_resid
Tau_resid ~ Mic.13_resid + Ab_resid
Ast.10_resid ~ Tau_resid + Mic.13_resid
Oli.7_resid ~ Ast.10_resid
cogng_demog_slope ~ Ast.10_resid + Tau_resid
'
model1.fitted<-sem(model1, data=df, fixed.x=FALSE)
summary(model1.fitted, standardized=TRUE, rsquare=TRUE, fit.measures=TRUE)
# Number of observations                           412
# P-value (Chi-square)                           0.284
# Comparative Fit Index (CFI)                    0.998
# Tucker-Lewis Index (TLI)                       0.996
# RMSEA                                          0.022
# 90 Percent confidence interval - lower         0.000
# 90 Percent confidence interval - upper         0.060
# P-value H_0: RMSEA <= 0.050                    0.865
# P-value H_0: RMSEA >= 0.080                    0.003

semPaths(model1.fitted,"std",edge.label.cex=1.5, label.font=1, curvePivot = TRUE, 
         nCharNodes=0, sideMan=50, Weighted=TRUE, esize=6,asize=6,
         node.width=3, sizeMan2=6, layout="tree2",residuals=FALSE,thresholds = TRUE,
         exoCov=FALSE, edge.label.position=0.65)



df <- df.full %>% filter(cohort == "replication") %>% filter(if_all(c(sqrt.amyloid_mf, sqrt.tangles_mf, cogng_demog_slope, Mic.13, pmi, RIN), ~ !is.na(.x))) # n=605

df$amyloid_resid <- resid(lm(sqrt.amyloid_mf~age_death+msex, df))
df$tangles_resid <- resid(lm(sqrt.tangles_mf~age_death+msex, df))
df$Mic.13_resid  <- resid(lm(Mic.13~age_death+msex+pmi+RIN, df))
df$Ast.10_resid  <- resid(lm(Ast.10~age_death+msex+pmi+RIN, df))
d4$Mic.12_resid  <- resid(lm(Mic.12_celmod_sqrt~age_death+msex+pmi+RIN, df))
d4$Oli.7_resid   <- resid(lm(Oli.7_celmod_sqrt~age_death+msex+pmi+RIN, df))

model1_celmod<-'
amyloid_resid ~ Mic.12_resid
Mic.13_resid ~ amyloid_resid + Mic.12_resid
tangles_resid ~ Mic.13_resid + amyloid_resid
Ast.10_resid ~ tangles_resid + Mic.13_resid
Oli.7_resid ~ Ast.10_resid
cogng_demog_slope ~ Ast.10_resid + tangles_resid
'

model1_celmod.fitted<-sem(model1_celmod, data=df, fixed.x=FALSE)
summary(model1_celmod.fitted, standardized=TRUE, rsquare=TRUE, fit.measures=TRUE)
semPaths(model1_celmod.fitted,"std",edge.label.cex=1.5, label.font=1, curvePivot = TRUE, 
         nCharNodes=0, sideMan=50, Weighted=TRUE, esize=6,asize=6,
         node.width=3, sizeMan2=6, layout="tree2",residuals=FALSE,thresholds = TRUE,
         exoCov=FALSE, edge.label.position=0.65)
 

df <- df.full %>% filter(cohort == "replication")
m<-lm(Ast.10~sqrt.tangles_mf+msex+age_death+msex+pmi+RIN, df); summary(m)
# sqrt.tangles_mf  2.507e-02  4.337e-03   5.780 1.16e-08 ***
m<-lm(Ast.10~sqrt.tangles_mf+Mic.13+age_death+msex+pmi+RIN, df); summary(m)
# sqrt.tangles_mf     3.383e-03  2.619e-03   1.292 0.196981    
# Mic.13  1.097e+00  3.130e-02  35.028  < 2e-16 ***
