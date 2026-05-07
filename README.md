# Background

This repository contains files for the manuscript:  

> Siepe, B. S., Kloft, M., Zhang, Y., Petersen, F., Bringmann, L. F., & Heck, D. W. (2025). Using features of dynamic networks to guide treatment selection and outcome prediction: The central role of uncertainty. PsyArXiv. https://doi.org/10.31234/osf.io/2c8xf_v1 

A timestamped version of this reproducibility archive can be found on Zenodo (DOI: [10.5281/zenodo.18155739](https://doi.org/10.5281/zenodo.18155739)).

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

# Quick navigation

| I want to … | Go to |
| --- | --- |
| Read the online supplement (simulation results, DGP details, deviations from preregistration) | `scripts/05_simulation_viz.html` |
| See results for a heterogeneous DGP | `scripts/05_simulation_viz_heterogeneous.html` |
| See additional simulations (varying sparsity & node count) | `scripts/09_additional_simulations_viz.html` |
| Read the empirical example walkthrough (incl. how to fit BmlVAR) | `scripts/02_empirical_example_bringmann_2016.html` |
| Reproduce figures using pre-computed output | `make render-viz` |
| Explore interactively in RStudio | `make rstudio` |

# Overview
### `\data\`

Contains data for our simulation study and empirical example.

| Dataset | Description |
| --- | --- |
| `data_Bringmann2016.RDS` | Data for the empirical example. Used in `02_empirical_example_bringmann_2016.qmd` |
| `graph_semisparse_synth.RDS` | True 6-node semi-sparse graph for main simulation. Generated in `04_dgps.qmd` |
| `graph_sparse_6n_highsparse.RDS` | True 6-node high-sparsity graph for additional simulations. Generated in `04_dgps.qmd` |
| `graph_semisparse_10n.RDS` | True 10-node semi-sparse graph for additional simulations. Generated in `04_dgps.qmd` |
| `graph_sparse_10n.RDS` | True 10-node sparse graph for additional simulations. Generated in `04_dgps.qmd` |
| `true_sd_semisparse_rev1.RDS` | True SDs for main simulation (revision 1). Generated in `01_centrality_simulation.qmd` |
| `true_sd_rev2_additional.RDS` | True SDs for revision 2 additional simulations. Generated in `09_additional_simulations_no_bayes.qmd` / `09_additional_simulations_only_bayes.qmd` |
| `true_sd.RDS` | True SDs for old DGM. No longer used, only kept for legacy. |
| `true_sd_semisparse.RDS` | True SDs for old DGM. No longer used, only kept for legacy. |

### `\scripts\`

Contains the scripts used in our project:

| Script | Description |
| --- | --- |
| `00_functions.R` | Contains auxiliary functions |
| `01_centrality_simulation.qmd` | Main simulation study (computationally intensive; outputs to `output/`) |
| `02_empirical_example_bringmann_2016.qmd` | Empirical example: Bringmann et al. (2016) |
| `04_dgps.qmd` | Data-generating processes; generates graph and SD files in `data/` |
| `05_simulation_viz.qmd` | Visualizations of the main simulation results |
| `05_simulation_viz_heterogeneous.qmd` | Visualizations for the heterogeneous DGP condition |
| `06_additional_mlvar_simulation.qmd` | Additional simulations investigating mlVAR performance |
| `07_exploratory_simulations.qmd` | Exploratory simulations |
| `08_centrality_values_exploration.qmd` | Exploration of centrality values |
| `09_additional_simulations.qmd` | Additional simulations for revision 2 (node count & sparsity; full pipeline) |
| `09_additional_simulations_no_bayes.qmd` | Revision 2 additional simulations without Bayesian estimator |
| `09_additional_simulations_only_bayes.qmd` | Revision 2 additional simulations with Bayesian estimator only |
| `09_additional_simulations_viz.qmd` | Visualizations for revision 2 additional simulations |

Some of these scripts have rendered `.html` versions. 
`05_simulation_viz.html` also contains information about the data-generating matrices and deviations from the preregistration.
`02_empirical_example_bringmann_2016.html` contains instructions on how to fit BmlVAR to data.

#### `\scripts\models\`
Contains the Stan models used in our project.

| Model | Description |
| --- | --- |
| `MLVAR_lkj_only_empirical_example.stan` | Stan model for the empirical example with two network features |
| `MLVAR_lkj_only.stan` | Stan model for the main simulation study |


### `\figures\`
Contains the figures generated in our project.

### `\output\`
Contains the simulation and empirical example output. 
As the full simulation files are too large to be stored on GitHub, we provide the summary of simulation results in `sim_results_rev1.RDS`. 

Our specific output structure is as follows: 

| File | Description |
| --- | --- |
| `sim_full.rds` | Simulation results for our initial manuscript version. Only kept for reproducibility purposes. |
| `sim_results_rev1_fixed_sigma.rds` | Results of the revised simulation setup (main) |
| `sim_results_rev1_heterogeneous.rds` | Results of the revised simulation setup for a heterogeneous DGP |
| `node1_centrality_rev1_fixed_sigma.rds` | Node 1 centrality estimates (all conditions) from the fixed-sigma simulation |
| `node1_centrality_rev1_cond14_fixed_sigma.rds` | Node 1 centrality estimates for condition 14 |
| `node1_centrality_rev1_cond23_fixed_sigma.rds` | Node 1 centrality estimates for condition 23 |
| `rev2/sim_full_rev2_additional_no_bayes.rds` | Revision 2 additional simulations without Bayesian estimator |
| `rev2/sim_full_rev2_additional_only_bayes.rds` | Revision 2 additional simulations with Bayesian estimator only |
| `rev2/sim_full_rev2_additional_row1_raw.rds` | Raw single-row results for revision 2 additional simulations |

The main simulation results (`sim_results_rev1_fixed_sigma.rds`, `sim_results_rev1_heterogeneous.rds`) are sourced automatically in `05_simulation_viz.qmd`. The revision 2 results are sourced in `09_additional_simulations_viz.qmd`.

# Reproducibility

All targets below require [Docker](https://docs.docker.com/get-docker/) and `make`.

## Pipeline overview

The scripts follow this dependency order:

```
04_dgps.qmd                        → data/*.RDS (graph structures, true SDs)
    ↓
01_centrality_simulation.qmd       → output/sim_results_rev1_fixed_sigma.rds
                                      output/sim_results_rev1_heterogeneous.rds
    ↓
05_simulation_viz.qmd              → figures/rev2/*.pdf
05_simulation_viz_heterogeneous.qmd

09_additional_simulations_no_bayes.qmd  }
09_additional_simulations_only_bayes.qmd} → output/rev2/*.rds
    ↓
09_additional_simulations_viz.qmd  → figures/rev2/*.pdf (additional simulations)

02_empirical_example_bringmann_2016.qmd → output/empirical_example/*.rds, figures/
```

> **Note:** `01_centrality_simulation.qmd` and the `09_additional_simulations_*.qmd` scripts are computationally intensive (intended for HPC). Pre-computed output files are provided in `output/` so that visualization scripts can be run without re-running the simulations.

## Make targets

| Target | Description |
| --- | --- |
| `make rstudio` | Launch interactive RStudio Server at `localhost:8787` |
| `make render` | Run full visualization pipeline non-interactively (uses pre-computed output) |
| `make render-dgps` | Render `04_dgps.qmd` to regenerate `data/` files |
| `make render-sim` | Render `01_centrality_simulation.qmd` (computationally intensive) |
| `make render-viz` | Render all visualization scripts (`05_*.qmd`, `09_*_viz.qmd`) |
| `make render-additional` | Render additional simulation scripts (`09_*_no_bayes.qmd`, `09_*_only_bayes.qmd`) |

### Reproducing figures (recommended)

```bash
make render-viz   # uses pre-computed output files; fast
```

### Interactive exploration

```bash
make rstudio      # opens RStudio Server at http://localhost:8787
```
