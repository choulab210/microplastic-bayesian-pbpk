# Bayesian PBPK Modeling of Micro- and Nanoplastic Biodistribution

This repository contains the R code and data used to develop an uncertainty-informed physiologically based pharmacokinetic (PBPK) model for oral micro- and nanoplastic biodistribution in mice.

The workflow integrates deterministic PBPK calibration, AUC-based sensitivity analysis, Bayesian MCMC calibration, random-effects meta-regression, and external evaluation without parameter refitting.

## Analysis Workflow

```text
Experimental biodistribution data
        ↓
Deterministic PBPK calibration
        ↓
AUC-based sensitivity analysis
        ↓
Bayesian MCMC calibration
        ↓
Random-effects meta-regression
        ↓
External PBPK evaluation
```

## Data

The analysis dataset is provided in:

```text
data/data.csv
```

### Model-Development Datasets

| ID | Particle size |
|---|---:|
| 1 | 20 nm |
| 2 | 220 nm |
| 3 | 1000 nm |
| 4 | 6000 nm |
| 5 | 220 nm |

### External-Evaluation Datasets

| ID | Particle size |
|---|---:|
| 6 | 50 nm |
| 7 | 500 nm |
| 8 | 5000 nm |

External datasets are used only for evaluation and are not used for parameter refitting.

## Repository Structure

```text
data/
    data.csv

model/
    MP_PBPK.R

scripts/
    01_fit_model.R
    02_run_sensitivity.R
    03_run_mcmc.R
    04_run_meta_regression.R
    05_external_evaluation.R
    06_make_figures.R
    07_make_tables.R

reference/
    manuscript_baseline/
    refactored_10k_analysis/
```

`reference/manuscript_baseline/` contains the original manuscript-era analysis scripts used to verify reproducibility.

`reference/refactored_10k_analysis/` contains an alternative refactored analysis retained for transparency and is not used for the final manuscript results.

## Reproducing the Analysis

Run the scripts sequentially from the repository root:

```bash
Rscript scripts/01_fit_model.R
Rscript scripts/02_run_sensitivity.R
Rscript scripts/03_run_mcmc.R
Rscript scripts/04_run_meta_regression.R
Rscript scripts/05_external_evaluation.R
Rscript scripts/06_make_figures.R
Rscript scripts/07_make_tables.R
```

The MCMC step is the most computationally intensive.

## Bayesian Analysis Settings

- Three MCMC chains per development dataset
- 1,000 pilot iterations
- 3,000 main iterations per chain
- Mixed prior-centering strategy
- Random seed: 1234
- Random-effects meta-regression fitted using REML
- External uncertainty propagation using 1,000 Monte Carlo simulations
- External evaluation restricted to IDs 6–8
- External datasets are evaluated without parameter refitting

## Main R Packages

Major R packages used in the workflow include:

```text
mrgsolve
FME
coda
metafor
dplyr
tidyr
ggplot2
readr
pheatmap
ggExtra
gridExtra
```

The analysis was finalized using R 4.5.x.

## Generated Outputs

Running the workflow generates:

```text
results/
figures/
tables/
```

Generated results, figures, tables, and large MCMC objects are not tracked by Git and can be recreated by running the analysis scripts.

## Citation

Please cite the associated manuscript describing the uncertainty-informed PBPK framework for micro- and nanoplastic biodistribution.

Citation information will be added after publication.

## Contact

For questions regarding the PBPK model or analysis workflow, please contact the corresponding author.