# BEYOND Methodology & Application To The DLPFC Atlas

<p align="justify">
Trait association analyses and causal modeling over cell proportions highlight key subpopulations and predict their causal role along the disease cascade. 
These however do not provide a complete view of the cellular changes leading to AD and how it might differ from brain aging. We are limited by the lack of 
longitudinal molecular brain measurements, as specimens are obtained at autopsy and do not provide information about temporal dynamics within each individual. 
To this end we devised the **BEYOND** methodology which enables the de-novo identification of trajectories of cellular change, and placement of individuals along them.
</p>
<p align="center"><img src="https://github.com/GreenGilad/BEYOND_DLPFC/assets/43610945/60b7ce9d-549c-44b3-8354-a50adf268496" width="75%"></p>

This framework consists of 4 steps:

### (1) Modeling the cellular landscape manifold

<p align="justify">BEYOND represents each participant by its cellular composition profile in our cell atlas to model their cellular environment. It then builds a cellular landscape manifold that captures the diversity of 
cellular environments, with each participant as a single point in the high dimensional manifold (and can be visualized in low dimension). 

File `1.create.cellular.landscape.R` contains the code to calculate participants' cellular environments and to model the DLPFC cellular landscape manifold. For visualization and exploration purposes we embedded 
participants' in a low-dimensional space using PHATE (2D and 3D), UMAP and tSNE (2D). For convenience we have arranged the participants' cellular environments, the landscape, embeddings and further information 
in an <a href="https://github.com/scverse/anndata">AnnData object</a> format.
</p>

> [!NOTE] 
> BEYOND is a conceptual framework for uncovering changes in the cellular landscape. As such it is not restricted to any of the above low-dimensioal embedding. Its application to different types of data may
> result in the use of different algorithms.

<br>

### (2) Reconstruction of trajectories of cellular change

<p align="justify">The next step of BEYOND is to model trajectories of change along the manifold, using similarities in the cellular environments across participants. In the current implementation we used two different algorithms 
of pseudotime analysis, <a href="https://github.com/dpeerlab/Palantir"> Palantir</a> and <a href="https://github.com/ShobiStassen/VIA">VIA</a>. File `2.trajectories.and.dynamics.in.landscape.R` contains the code for 
running both algorithms and storing the results within the AnnData object.</p>

<br>

### (3) Assesses dynamics of sub-populations and disease-related traits

<p align="justify">Next, BEYOND assesses whether specific subpopulations or traits of interest are associated with the trajectories. It does so by calculating the dynamics of a feature of interest along the different 
  trajectories. This is done using the participants' pseudotime assignment and trajectory probabilities, provided in the step above. File `2.trajectories.and.dynamics.in.landscape.R` contains the code for 
fitting the dynamics and storing the results within the AnnData object. </p>

<br>

### (4) Constructs multi-cellular communities

Finally, BEYOND assigns cell subpopulations to multi-cellular communities. We expand our previous work (<a href="https://www.nature.com/articles/s41593-023-01356-x">Cain A. _et al._, 2023</a>) and cluster is performed using two properties reflecting coordinated changes
- Co-occurrences of subpopulations across participants. Calculated by the Spearman correlation.
- Shared patterns of dynamics along all found trajectories. Calculated by an adaptive RBF-kernel over the predicted dynamics.

File `3.construct.cellular.communities.R` contains the code for constructing the communities, as well as evaluating their dynamics, association with disease traits and pathways shared for subpopulations within communities.

<br>

### Validating BEYOND

We validate BEYOND in an independent set of participants whose cellular environments are obtained by bulk-to-single-cell deconvolution using 
the <a href="https://github.com/MenonLab/Celmod">CelMod</a> algorithm (<a href="https://www.nature.com/articles/s41593-023-01356-x">Cain A. _et al._, 2023</a>). File `BEYOND.validations.R` contains the relevant code.

- The code for fitting the Celmod models over the used bulk RNA-seq and predicting the subpopulation proportions is found in `../Other analyses/CelMod.fitting.R`.
- While this code is not directly part of the BEYOND methodology, for convenience of implementation results were stored within the same DLPFC landscape AnnData object.
