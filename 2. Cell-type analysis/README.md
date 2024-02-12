# Cell-type Level snRNA-seq Analysis

<p align="justify">
  The following files contain the code used to analyze data of the specified cell type, resulting in the final set of annotated cells and functionality analysis (DE and pathway enrichment) used in out cell atlas.
The starting point for the cell-type analysis is the successful execution of the library analysis for each one of the 127 snRNA-seq libraries generated in this study (see code under Library preprocessing folder)
</p>
<p align="justify">The final cell-type objects (Seurat H5 format) are available on <a href="https://www.synapse.org/#!Synapse:syn53366818">Synapse</a></p>
<p align="center"><img src="https://github.com/GreenGilad/BEYOND_DLPFC/assets/43610945/832ade56-b4e8-46ba-a358-3535e83487b8" width="40%"></p>

Once cell-type analysis was performed for all cell-types:
- `2.aggregate.cell.type.data.R` - Cell-type annotations, UMAP coordinates, differential expression, pathway analysis and pseudobulk gene expressions were pulled into a single hierarchical file.
- `3.create.proportion.matrix.R` - Creating the subpopulation proportion matrix (whose rows represent the different 437 participants and columns the 91 subpopulations of the atlas) and storing with participant- and subpopulation metadata in an <a href="https://anndata.readthedocs.io/en/latest/">AnnData</a> format.

All downstream analysis (such as trait-associations, mediation modeling and BEYOND), as well as manuscript figures and supplementary tables, used the data present in these files.
