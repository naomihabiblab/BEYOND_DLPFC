
<p align="justify">
This folder contains additional analyses performed in this study

- `CelMod.fitting.R` - predicting in bulk RNA-seq the proportions of snRNA-seq subpopulations defined in our cell atlas.
- `trait.associations.R` - associating subpopulation proportions to disease traits. Associations are performed over snRNA-seq subpopulation proportions (our discovery cohort of 437 participants)
  and over Celmod predicted proportions (our replication cohort of 673 independent participants), and followed by a meta-analysis of both cohorts.
- `RNAscope.analysis.R` - quantification and analysis of RNAscope results used in this study.

<br>

> [!IMPORTANT]
> While not part of BEYOND, the code used above loads the subpopulation proportions calculated in `../BEYOND/1.create.cellular.landscape.R`.
</p>
