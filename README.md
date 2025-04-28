# Background

This repository contains files for the manuscript:  

> Siepe, B. S., Kloft, M., Zhang, Y., Petersen, F., Bringmann, L. F., & Heck, D. W. (2025). Using features of dynamic networks to guide treatment selection and outcome prediction: The central role of uncertainty. PsyArXiv. https://doi.org/10.31234/osf.io/2c8xf_v1 


Please cite the associated preprint as: 

```BibTeX
@article{Siepe2025,
  title={Using features of dynamic networks to guide treatment selection and outcome prediction: The central role of uncertainty},
  author={Siepe, Bj{\"o}rn S. and Kloft, Matthias and Zhang, Yong and Petersen, Fridtjof and Bringmann, Laura F. and Heck, Daniel W.},
  url={https://doi.org/10.31234/osf.io/2c8xf_v1},
  year={2025},
  note={PsyArXiv preprint}
}
```

# Overview
### `\data\`

Contains data for our simulation study and empirical example. 

### `\scripts\`

Contains the scripts used in our project:

| Script | Description |
| --- | --- |
| `00_functions.R` | contains auxiliary functions |
| `01_centrality_simulation.qmd` | contains code to replicate the main simulation study |
| `02_empirical_example_bringmann_2016.qmd` | empirical example: Bringmann et al. (2016) |
| `04_dgps.qmd` | code for data-generating processes |
| `05_simulation_viz.qmd` | code to analyze \& visualize the simulation results |
| `06_additional_mlvar_simulations.qmd` | additional simulations for the mlVAR model |
| `07_exploratory_simulations.qmd` | exploratory simulations |

Some of these scripts have rendered `.html` versions. 
`05_simulation_viz.html` also contains information about the data-generating matrices and deviations from the preregistration.
`02_empirical_example_bringmann_2016.html` contains instructions on how to fit BmlVAR to data.

#### `\scripts\models\`
Contains the Stan models used in our project.

| Model | Description |
| --- | --- |
| `MLVAR_lkj_only_empirical_example.stan` | Stan model for the empirical example with two network features |
| `MLVAR_lkj_only.stan` | Stan model for the simulation study |


### `\figures\`
Contains the figures generated in our project.

### `\output\`
Contains the simulation and empirical example output. 
As the full simulation files are too large to be stored on GitHub, we provide the summary of simulation results in `sim_results.RDS`. 


# Reproducibility

To rerun the computation of simulation figures and summaries via a Docker container:

1. Install Docker and make on your machine. 
2. Run `make docker` from the root directory. This will install all necessary dependencies.
3. RStudio can then be opened from a browser at `localhost:8787`, and `05_simulation_viz.qmd` can be run to reproduce the simulation results.
