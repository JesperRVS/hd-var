# hd-var: Lasso and Related Tools for High-Dimensional Vector Autoregression
Maintained by: Jesper Riis-Vestergaard Sørensen, University of Copenhagen

## Overview
This repository contains an ``R`` implementation of Lasso and variations thereof [such as least squares refitting following Lasso selection (post-Lasso) and
square-root Lasso (sqrt-Lasso)] with a data-driven and theoretically justifiable
tuning parameter selection method designed for high-dimensional (HD) vector
autoregression (VAR).

The code in this repository was developed for and used in the paper
"***Data-Driven Tuning Parameter Selection for High-Dimensional Vector
Autoregressions***," authored by [Anders Bredahl
Kock](https://sites.google.com/site/andersbkock/), [Rasmus Søndergaard
Pedersen](https://sites.google.com/site/rspecon/) and [Jesper Riis-Vestergaard
Sørensen](https://sites.google.com/site/jesperrvs) (henceforth: KPS), and
available at https://arxiv.org/abs/2403.06657. See the paper for the theoretical
justification of the method.

## Architecture
The [root directory](.) contains the *main* scripts and functions: 
* ``lassoVAR.R``: Contains ``lasso_var`` which fits a VAR by (weighted) Lasso
  (possibly with refitting post variable selection, i.e. post-Lasso) using the
  penalization algorithm(s) in KPS.
* ``sqrtLassoVAR.R``: Contains ``sqrt_lasso_var`` which fits a VAR by (weighted) sqrt-Lasso using the penalization algorithm in KPS.

It also contains the *auxilliary* scripts and functions:
* ``icLassoVAR.R``: Contains ``ic_lasso_var`` which fits a VAR by Lasso using an
  information criterion (Akaike, Bayes or Hannan-Quinn) to determined the
  penalty level.
* ``helper_functions.R``: Contains various helper functions drawn upon by the
  above functions. These helpers are used to unpack the data to align with a
  VAR, fit a linear model using least squares (or refit post variable
  selection), and apply Lasso (or variations thereof) equation by equation.

The subfolder [simulations](./simulations/) contains the following files:
* ``simData.R``: Script with functions which allow for simulation using the
  various designs described in Section 6.
* ``runSim_v03.R``: Script running the simulations reported on in Section 6.
* ``simulations_workspace_1000_MC_100_to_1000_n_16_to_128_p_with_num_upd.RData``:
The workspace produced by ``runSim_v03.R``.
* ``createFigs_v07.R``: Script producing Figure 6.1 (relative average estimation
errors) in Section 6 and Figure H.1 in Supplementary Appendix H (raw average
estimation errors).
* ``markupDependence_v01.R``: Script running the simulations used to investigate
the mark-up dependence reported on in Figure H.2 of Supplementary Appendix H.
* ``markup_dependence_workspace_1000_MC_200_to_1000_n_16_to_128_p_diagonal_only.RData``: The workspace produced by ``markupDependence_v01.R``.
* ``markupDependenceFigs_v02.R``: Script producing Figure H.2 in Supplementary
Appendix H (raw average estimation errors).
* Figures (3) from the various simulations can be found in/are output to
[simulations/img](./simulations/img/).

The subfolder [application/FRED](./application/FRED/) contains:
* [application/FRED/pre_processing](./application/FRED/pre_processing/):
  Contains the raw FRED-MD dataset (from May, 2022) and three Matlab scripts
  (``FREDMD_preprocess.m``, ``prepare_missing.m`` and ``remove_outliers.m``)
  used to pre-process this raw dataset to end up with the dataset used in the
  empirical illustration in Section 7. (See also the *Reproducibility Workflow*
  below.)
* ``FRED-MD_2022-05_preprocessed.csv``: The FRED-MD dataset used in the
 empirical illustration in Section 7
([application/FRED/data](./application/FRED/data/))
* ``FRED-MD_forecasting_v02.R``: Script conducting the forecasting exercise
using this FRED-MD dataset and reported on in Section 7.
* ``application_workspace_N_120_qmax_12_with_num_upd.RData``: The workspace
  produced by ``FRED-MD_forecasting_v02.R``.
* ``FRED-MD_figures_v02.R``: Script producing both parts of Figure 7.1 (average
  and 95th percentile inverse variance weighted square forecast error).
* Figures (2) from the empirical illustration can be found in/are output to
  [application/FRED/img](./application/FRED/img/)
* 

## Reproducibility Workflow

**To reproduce the simulations:** Execute ```runSim_v03.R``` to run the main
text simulations and ```markupDependence_v01.R``` for the mark-up dependence
simulations, specifically. To have access to ```lassoVAR.R``` and
```sqrtLassoVAR.R```, the former two scripts must be executed from the parent
folder. Use ```setwd("..")``` to back up, if necessary. When executed in a Linux
environment, the scripts will automatically save the workspaces in the
simulations folder. If a different environment is used, the workspace must be
saved manually (via ```save.image(...)```). The figures included in the paper
arise from calling the scripts ```createFigs_v07.R``` and
```markupDependenceFigs_v02.R``` from their folder, which load the relevant
workspaces from the simulations folder and produce figures in the
```simulations/img``` folder.

**To reproduce the empirical illustration:**
Execute ```FRED-MD_forecasting_v02.R ``` 
(from the parent folder via ```setwd("../..")```) to conduct the forecasting exercise.
The script loads the pre-processed data from the ```application/FRED/data``` folder
and automatically saves the workspace (if in Linux).
Run ```FRED-MD_figures_v02.R``` to produce the figures in the empirical illustration
in the ```application/FRED/img``` folder.

***Pre-processing:*** To pre-process the raw FRED-MD dataset (from May, 2022),
run the Matlab script ``FREDMD_preprocess.m`` (changing the necessary paths).
The data pre-processing amounts to carrying out (deterministic) stationarity
inducing transformations using ``prepare_missing.m``, then removing outliers
with ``remove_outliers.m``, and finally replacing missing values with the
unconditional average of the corresponding series as in the initialization of
``factors_em.m`` (not included here). The latter three scripts 
are from Michael W. McCracken's
[website](https://www.stlouisfed.org/research/economists/mccracken/fred-databases) (accessed in May, 2022).