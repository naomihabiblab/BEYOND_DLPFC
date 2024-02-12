# Analyses Over Cell Atlas
<p align="justify">
Various types of analyses and validations were performed in this study

- `1.CelMod.fitting.R` - Predicting snRNA-seq subpopulation proportions within bulk RNA-seq data.
- `2.trait.associations.R` - Associating subpopulation proportions to disease traits. Associations are performed over snRNA-seq subpopulation proportions (our discovery cohort of 437 participants)
  and over Celmod predicted proportions (our replication cohort of 673 independent participants), and followed by a meta-analysis of both cohorts.
- `3.causal.modeling.R` - Mediation modeling to position subpopulation changes w.r.t AD key endophenotypes.
- `4.RNAscope.analysis.R` - Quantification and analysis of microglial subpopulations using smFISH.
- `5.pTau.analysis.R` - Quantification of pTau and association with microglical subpopulations (from smFISH data), validating pTau associations found in the snRNA-seq data
- `6.morphology.analysis.R` - Assessing microglia morphology in smFISH images and associating with subpopulation proportions, as well as associating microglial activation stage with subpopulation proportion.
- `7.spatial.transcriptomics.validations.R` - Quntification gene signatures of key subpopulations in spatial transcriptomics data.
