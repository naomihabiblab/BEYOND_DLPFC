# BEYOND - Uncovering trajectories of cellular change in the aged human cellular landscape

<p align="justify">
To reconstruct the dynamics of the brain's cellular environment along the progression of Alzheimer's Disease, and to distinguish between AD and aging effects, we built a comprehensive cell atlas of the aged prefrontal cortex from 1.64 million single-nucleus RNA-seq profiles. We associated glial, vascular and neuronal subpopulations with AD-related traits for 437 aged individuals and aligned them along the disease cascade using causal modeling. To model the coordinated dynamics of the entire cellular environment, we devised the BEYOND methodology, which uncovered two distinct trajectories of brain aging defined by differing sequences of changes in cellular communities
</p>
<p align="center"><img src="https://github.com/GreenGilad/BEYOND_DLPFC/assets/43610945/d5e52110-1bee-452d-bf8a-1c33e0ede755" width="75%"></p>

<p align="justify">
The following repository contains:
  
- **Library level snRNA-seq analysis**: Removal of background RNA, quality control steps, cell-type classification and doublets annotation.
- **Cell-type level snRNA-seq analysis**: Creating our cell atlas. Removal of low quality cells and doublets, sub-clustering analysis, differential expression and pathway enrichment analyses.
- **BEYOND analysis and validations**: Applying BEYOND methodology steps over our cell atlas.
- **Other analyses**: Such as bulk RNA-seq deconvolution using CelMod, trait associations and meta-analysis, and RNAscope validations.
- **Manuscript figure and supplementary tables code**: The complete code generating all figures, extended data figures and supplementary tables shown in the manuscript

> [!NOTE]
> The code within each folder runs under the assumption that the code of the previous folder was executed and generated the necessary data files. Please refer to the readme file of each folder for more details

## Citation
- Green G. S., F. M., Yang H-S., Taga M., McCabe C., Cain A., White C. C., Schmidtner A. K., Zeng L., Wang Y., Regev A., Menon V., Bennett D. A, Habib N., De Jager P. L. Cellular dynamics across aged human brains uncover a multicellular cascade leading to Alzheimer’s disease. Preprint from bioRvix (2023). <a href="https://www.biorxiv.org/content/10.1101/2023.03.07.531493v1">doi:10.1101/2023.03.07.531493</a>

## Contact Us
