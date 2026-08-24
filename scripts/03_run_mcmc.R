## ============================================================
## 03_run_mcmc.R
##
## FINAL MANUSCRIPT-COMPATIBLE BAYESIAN PBPK CALIBRATION
##
## Purpose
## -------
## Reproduce the exact manuscript-era MCMC analysis that was
## independently verified against the frozen manuscript RDS.
##
## Strategy
## --------
## 1. Run the frozen manuscript deterministic + MCMC scripts in
##    an isolated temporary working directory.
## 2. Preserve their original numerical behavior exactly.
## 3. Copy canonical manuscript-era outputs to results/mcmc/.
## 4. Export standardized per-ID CSV/RDS files for Scripts 04–07.
##
## Scientific settings
## -------------------
## - set.seed(1234)
## - mixed prior centers
## - PRIOR_WEIGHT = 0.75
## - pilot = 1,000 iterations
## - main = 3,000 iterations
## - 3 chains per development dataset
## - 20% burn-in handled by original FME::modMCMC workflow
## - delayed rejection (ntrydr = 2; drscale = 2)
## - no adaptive covariance (updatecov = 0)
##
## Development IDs
## ---------------
## IDs 1–5
##
## IMPORTANT
## ---------
## This script intentionally uses the frozen manuscript baseline
## implementations of:
##
##   reference/manuscript_baseline/01_run_model_fit.R
##   reference/manuscript_baseline/03_run_mcmc_per_id.R
##
## because that exact combination has already been demonstrated
## to reproduce the frozen manuscript MCMC posterior summaries
## with zero numerical difference.
##
## It does NOT use the archived 10,000-iteration refactored MCMC.
## ============================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(coda)
})

options(stringsAsFactors = FALSE)


## ============================================================
## 0) USER OPTIONS
## ============================================================

MASTER_SEED <- 1234

DEV_IDS <- 1:5

KEY_META_PARAMS <- c(
  "Ktrans",
  "Kgw2b",
  "Kshed",
  "KLRESmax",
  "Kfeces",
  "KurineC"
)

CLEAN_MAIN_MCMC_RESULTS <- TRUE

KEEP_ISOLATED_WORKDIR <- TRUE


## ============================================================
## 1) PROJECT PATHS
## ============================================================

PROJECT_ROOT <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

BASELINE_DIR <- file.path(
  PROJECT_ROOT,
  "reference",
  "manuscript_baseline"
)

BASELINE_01 <- file.path(
  BASELINE_DIR,
  "01_run_model_fit.R"
)

BASELINE_03 <- file.path(
  BASELINE_DIR,
  "03_run_mcmc_per_id.R"
)

MODEL_FILE <- file.path(
  PROJECT_ROOT,
  "model",
  "MP_PBPK.R"
)

DATA_FILE <- file.path(
  PROJECT_ROOT,
  "data",
  "data.csv"
)

RESULTS_DIR <- file.path(
  PROJECT_ROOT,
  "results"
)

MAIN_MCMC_DIR <- file.path(
  RESULTS_DIR,
  "mcmc"
)

FIGURES_DIR <- file.path(
  PROJECT_ROOT,
  "figures"
)

MAIN_CONV_DIR <- file.path(
  FIGURES_DIR,
  "mcmc_convergence"
)

MAIN_PPC_DIR <- file.path(
  FIGURES_DIR,
  "mcmc_posterior"
)

WORK_DIR <- file.path(
  RESULTS_DIR,
  "manuscript_mcmc_work"
)


## ============================================================
## 2) REQUIRED FILE CHECK
## ============================================================

required_files <- c(
  BASELINE_01,
  BASELINE_03,
  MODEL_FILE,
  DATA_FILE
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (
  length(
    missing_files
  ) > 0
) {

  stop(
    "\nMissing required file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    "\n"
  )
}


cat(
  "\n============================================================\n"
)

cat(
  "MANUSCRIPT-COMPATIBLE MCMC\n"
)

cat(
  "============================================================\n"
)

cat(
  "Project root: ",
  PROJECT_ROOT,
  "\n",
  sep = ""
)

cat(
  "Frozen Script 01: ",
  BASELINE_01,
  "\n",
  sep = ""
)

cat(
  "Frozen Script 03: ",
  BASELINE_03,
  "\n",
  sep = ""
)

cat(
  "Seed: ",
  MASTER_SEED,
  "\n",
  sep = ""
)

cat(
  "Development IDs: ",
  paste(
    DEV_IDS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


## ============================================================
## 3) PREPARE ISOLATED MANUSCRIPT WORKING DIRECTORY
## ============================================================

if (
  dir.exists(
    WORK_DIR
  )
) {

  unlink(
    WORK_DIR,
    recursive = TRUE,
    force = TRUE
  )
}

dir.create(
  WORK_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ------------------------------------------------------------
## Copy exact frozen manuscript scripts
## ------------------------------------------------------------

ok_copy <- file.copy(
  from = BASELINE_01,
  to = file.path(
    WORK_DIR,
    "01_run_model_fit.R"
  ),
  overwrite = TRUE
)

if (!ok_copy) {
  stop(
    "Failed to copy frozen Script 01."
  )
}


ok_copy <- file.copy(
  from = BASELINE_03,
  to = file.path(
    WORK_DIR,
    "03_run_mcmc_per_id.R"
  ),
  overwrite = TRUE
)

if (!ok_copy) {
  stop(
    "Failed to copy frozen Script 03."
  )
}


## ------------------------------------------------------------
## Copy model/data using filenames expected by old scripts
## ------------------------------------------------------------

ok_copy <- file.copy(
  from = MODEL_FILE,
  to = file.path(
    WORK_DIR,
    "MP_PBPK.R"
  ),
  overwrite = TRUE
)

if (!ok_copy) {
  stop(
    "Failed to copy model/MP_PBPK.R."
  )
}


ok_copy <- file.copy(
  from = DATA_FILE,
  to = file.path(
    WORK_DIR,
    "data.csv"
  ),
  overwrite = TRUE
)

if (!ok_copy) {
  stop(
    "Failed to copy data/data.csv."
  )
}


dir.create(
  file.path(
    WORK_DIR,
    "results"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(
    WORK_DIR,
    "figures"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 4) RUN EXACT MANUSCRIPT WORKFLOW
## ============================================================

set.seed(
  MASTER_SEED
)

old_wd <- getwd()

on.exit(
  setwd(
    old_wd
  ),
  add = TRUE
)

setwd(
  WORK_DIR
)


cat(
  "\n============================================================\n"
)

cat(
  "STEP 1: Recreate manuscript deterministic environment\n"
)

cat(
  "============================================================\n"
)


source(
  "01_run_model_fit.R",
  local = .GlobalEnv
)


required_objects <- c(
  "MODEL_PARAM_NAMES",
  "ids",
  "Fit_list",
  "obs_df"
)


missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]


if (
  length(
    missing_objects
  ) > 0
) {

  stop(
    "\nFrozen Script 01 did not generate:\n",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}


if (
  !exists(
    "pred",
    mode = "function",
    envir = .GlobalEnv
  )
) {

  stop(
    "Frozen deterministic pred() function was not generated."
  )
}


if (
  !exists(
    "MCresid_id",
    mode = "function",
    envir = .GlobalEnv
  )
) {

  stop(
    "Frozen deterministic MCresid_id() was not generated."
  )
}


cat(
  "\nDeterministic manuscript environment recreated successfully.\n"
)


## ============================================================
## 5) VERIFY MANUSCRIPT MCMC SETTINGS BEFORE RUNNING
## ============================================================

old03_text <- readLines(
  "03_run_mcmc_per_id.R",
  warn = FALSE
)


setting_lines <- grep(
  paste(
    c(
      "^set.seed",
      "^PRIOR_CENTER_STRATEGY",
      "^PRIOR_WEIGHT",
      "^PILOT_ITER",
      "^MAIN_ITER",
      "^NSAMP_PPC"
    ),
    collapse = "|"
  ),
  old03_text,
  value = TRUE
)


cat(
  "\nFrozen MCMC settings:\n"
)

cat(
  paste0(
    "  ",
    setting_lines,
    collapse = "\n"
  ),
  "\n"
)


## ------------------------------------------------------------
## Hard checks for settings underlying manuscript reproduction
## ------------------------------------------------------------

if (
  !any(
    grepl(
      "PILOT_ITER\\s*<-\\s*1000",
      old03_text
    )
  )
) {

  stop(
    "Frozen MCMC script no longer has PILOT_ITER <- 1000."
  )
}


if (
  !any(
    grepl(
      "MAIN_ITER\\s*<-\\s*3000",
      old03_text
    )
  )
) {

  stop(
    "Frozen MCMC script no longer has MAIN_ITER <- 3000."
  )
}


if (
  !any(
    grepl(
      'PRIOR_CENTER_STRATEGY\\s*<-\\s*"mixed"',
      old03_text
    )
  )
) {

  stop(
    "Frozen MCMC script no longer uses mixed prior centers."
  )
}


## ============================================================
## 6) RUN EXACT FROZEN MCMC
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "STEP 2: Run exact manuscript MCMC\n"
)

cat(
  "============================================================\n"
)


source(
  "03_run_mcmc_per_id.R",
  local = .GlobalEnv
)


## ============================================================
## 7) LOCATE CANONICAL MANUSCRIPT OUTPUTS
## ============================================================

LEGACY_MCMC_DIR <- file.path(
  WORK_DIR,
  "results",
  "mcmc"
)


CANONICAL_MCMC_FILE <- file.path(
  LEGACY_MCMC_DIR,
  "MCMC3_results_fullset_perID_priors_mixed.rds"
)

CANONICAL_PARAM_RDS <- file.path(
  LEGACY_MCMC_DIR,
  "ParamSummary_fullset_perID_priors_mixed.rds"
)

CANONICAL_PPC_RDS <- file.path(
  LEGACY_MCMC_DIR,
  "PPCSummary_fullset_perID_priors_mixed.rds"
)

CANONICAL_GELMAN_RDS <- file.path(
  LEGACY_MCMC_DIR,
  "GelmanDiag_fullset_perID_priors_mixed.rds"
)

CANONICAL_ESS_RDS <- file.path(
  LEGACY_MCMC_DIR,
  "ESS_fullset_perID_priors_mixed.rds"
)


required_outputs <- c(
  CANONICAL_MCMC_FILE,
  CANONICAL_PARAM_RDS,
  CANONICAL_PPC_RDS,
  CANONICAL_GELMAN_RDS,
  CANONICAL_ESS_RDS
)


missing_outputs <- required_outputs[
  !file.exists(
    required_outputs
  )
]


if (
  length(
    missing_outputs
  ) > 0
) {

  stop(
    "\nManuscript MCMC finished, but expected output(s) are missing:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}


## ============================================================
## 8) RESTORE PROJECT ROOT
## ============================================================

setwd(
  PROJECT_ROOT
)


## ============================================================
## 9) CLEAN / CREATE FINAL MAIN OUTPUT DIRECTORIES
## ============================================================

if (
  CLEAN_MAIN_MCMC_RESULTS &&
    dir.exists(
      MAIN_MCMC_DIR
    )
) {

  cat(
    "\nRemoving previous results/mcmc to prevent mixing workflows...\n"
  )

  unlink(
    MAIN_MCMC_DIR,
    recursive = TRUE,
    force = TRUE
  )
}


dir.create(
  MAIN_MCMC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  MAIN_CONV_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  MAIN_PPC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 10) COPY CANONICAL MANUSCRIPT FILES
## ============================================================

legacy_files <- list.files(
  LEGACY_MCMC_DIR,
  full.names = TRUE
)


for (
  f in legacy_files
) {

  file.copy(
    from = f,
    to = file.path(
      MAIN_MCMC_DIR,
      basename(f)
    ),
    overwrite = TRUE
  )
}


## ============================================================
## 11) LOAD CANONICAL OBJECTS
## ============================================================

MCMC3_FINAL <- readRDS(
  CANONICAL_MCMC_FILE
)

PARAM_FINAL <- readRDS(
  CANONICAL_PARAM_RDS
)

PPC_FINAL <- readRDS(
  CANONICAL_PPC_RDS
)

GELMAN_FINAL <- readRDS(
  CANONICAL_GELMAN_RDS
)

ESS_FINAL <- readRDS(
  CANONICAL_ESS_RDS
)


## ============================================================
## 12) EXPORT STANDARDIZED PER-ID RESULTS
##
## These files are convenience outputs for Scripts 04–07.
## They DO NOT alter the underlying manuscript MCMC.
## ============================================================

for (
  id in DEV_IDS
) {

  id_chr <- as.character(
    id
  )

  cat(
    "\nExporting standardized results for ID ",
    id,
    "...\n",
    sep = ""
  )


  if (
    is.null(
      MCMC3_FINAL[[id_chr]]
    )
  ) {

    stop(
      "Canonical MCMC object is missing ID ",
      id
    )
  }


  ID_DIR <- file.path(
    MAIN_MCMC_DIR,
    paste0(
      "ID",
      id
    )
  )


  dir.create(
    ID_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )


  ## ----------------------------------------------------------
  ## 12a) Save original per-ID MCMC object
  ## ----------------------------------------------------------

  saveRDS(
    MCMC3_FINAL[[id_chr]],
    file.path(
      ID_DIR,
      "mcmc_results.rds"
    )
  )


  ## ----------------------------------------------------------
  ## 12b) Exact posterior summaries from actual pooled draws
  ## ----------------------------------------------------------

  draws <- do.call(
    rbind,
    lapply(
      MCMC3_FINAL[[id_chr]]$chains,
      function(m) {

        as.data.frame(
          m$pars
        )
      }
    )
  )


  if (
    nrow(draws) == 0
  ) {

    stop(
      "No posterior draws found for ID ",
      id
    )
  }


  POST_SUM <- bind_rows(
    lapply(
      names(draws),
      function(parameter) {

        x_log <- as.numeric(
          draws[[parameter]]
        )

        x_nat <- exp(
          x_log
        )


        tibble(
          ID = id,

          Parameter = parameter,

          mean_log = mean(
            x_log,
            na.rm = TRUE
          ),

          sd_log = sd(
            x_log,
            na.rm = TRUE
          ),

          mean = mean(
            x_nat,
            na.rm = TRUE
          ),

          median = median(
            x_nat,
            na.rm = TRUE
          ),

          q025 = as.numeric(
            quantile(
              x_nat,
              0.025,
              na.rm = TRUE
            )
          ),

          q975 = as.numeric(
            quantile(
              x_nat,
              0.975,
              na.rm = TRUE
            )
          )
        )
      }
    )
  )


  write_csv(
    POST_SUM,
    file.path(
      ID_DIR,
      "posterior_parameter_summary.csv"
    )
  )


  ## ----------------------------------------------------------
  ## 12c) Posterior time-course uncertainty summary
  ## ----------------------------------------------------------

  ppc_i <- PPC_FINAL[[id_chr]]


  if (
    is.null(
      ppc_i
    )
  ) {

    stop(
      "Posterior predictive summary missing for ID ",
      id
    )
  }


  ppc_i <- as_tibble(
    ppc_i
  )


  write_csv(
    ppc_i,
    file.path(
      ID_DIR,
      "posterior_prediction_summary.csv"
    )
  )


  ## ----------------------------------------------------------
  ## 12d) Acceptance-rate table
  ## ----------------------------------------------------------

  acc_i <- as.numeric(
    MCMC3_FINAL[[id_chr]]$accept
  )


  ACC_DF <- tibble(
    ID = id,

    Chain = seq_along(
      acc_i
    ),

    acceptance_rate = acc_i,

    jump = MCMC3_FINAL[[id_chr]]$jump,

    var0 = MCMC3_FINAL[[id_chr]]$var0
  )


  write_csv(
    ACC_DF,
    file.path(
      ID_DIR,
      "acceptance.csv"
    )
  )


  ## ----------------------------------------------------------
  ## 12e) ESS
  ## ----------------------------------------------------------

  ess_i <- ESS_FINAL[[id_chr]]


  if (
    !is.null(
      ess_i
    ) &&
      !all(
        is.na(
          ess_i
        )
      )
  ) {

    ESS_DF <- tibble(
      ID = id,

      Parameter = names(
        ess_i
      ),

      ESS = as.numeric(
        ess_i
      )
    )

  } else {

    ESS_DF <- tibble(
      ID = numeric(),
      Parameter = character(),
      ESS = numeric()
    )
  }


  write_csv(
    ESS_DF,
    file.path(
      ID_DIR,
      "ess.csv"
    )
  )


  ## ----------------------------------------------------------
  ## 12f) Gelman-Rubin / Rhat
  ## ----------------------------------------------------------

  gd_i <- GELMAN_FINAL[[id_chr]]


  if (
    inherits(
      gd_i,
      "gelman.diag"
    ) ||
      (
        is.list(
          gd_i
        ) &&
          !is.null(
            gd_i$psrf
          )
      )
  ) {

    psrf <- as.matrix(
      gd_i$psrf
    )


    GEL_DF <- tibble(
      ID = id,

      Parameter = rownames(
        psrf
      ),

      Rhat = as.numeric(
        psrf[
          ,
          1
        ]
      ),

      Rhat_upper = if (
        ncol(
          psrf
        ) >= 2
      ) {

        as.numeric(
          psrf[
            ,
            2
          ]
        )

      } else {

        NA_real_
      }
    )

  } else {

    GEL_DF <- tibble(
      ID = numeric(),
      Parameter = character(),
      Rhat = numeric(),
      Rhat_upper = numeric()
    )
  }


  write_csv(
    GEL_DF,
    file.path(
      ID_DIR,
      "gelman.csv"
    )
  )
}


## ============================================================
## 13) COPY ORIGINAL CONVERGENCE / PPC FIGURES
## ============================================================

LEGACY_CONV_DIR <- file.path(
  WORK_DIR,
  "figures",
  "mcmc_convergence"
)


if (
  dir.exists(
    LEGACY_CONV_DIR
  )
) {

  conv_files <- list.files(
    LEGACY_CONV_DIR,
    full.names = TRUE
  )


  for (
    f in conv_files
  ) {

    file.copy(
      from = f,
      to = file.path(
        MAIN_CONV_DIR,
        basename(f)
      ),
      overwrite = TRUE
    )
  }
}


legacy_ppc_files <- list.files(
  file.path(
    WORK_DIR,
    "figures"
  ),
  pattern = "^PosteriorPredictive_ID.*\\.tiff$",
  full.names = TRUE
)


for (
  f in legacy_ppc_files
) {

  file.copy(
    from = f,
    to = file.path(
      MAIN_PPC_DIR,
      basename(f)
    ),
    overwrite = TRUE
  )
}


## ============================================================
## 14) FINAL DIAGNOSTIC SUMMARY
## ============================================================

DIAGNOSTIC_SUMMARY <- bind_rows(
  lapply(
    DEV_IDS,
    function(id) {

      id_chr <- as.character(
        id
      )

      acc <- as.numeric(
        MCMC3_FINAL[[id_chr]]$accept
      )


      ess <- ESS_FINAL[[id_chr]]


      gd <- GELMAN_FINAL[[id_chr]]


      rhat <- if (
        is.list(gd) &&
          !is.null(
            gd$psrf
          )
      ) {

        as.numeric(
          gd$psrf[
            ,
            1
          ]
        )

      } else {

        NA_real_
      }


      tibble(
        ID = id,

        acceptance_min = min(
          acc,
          na.rm = TRUE
        ),

        acceptance_max = max(
          acc,
          na.rm = TRUE
        ),

        acceptance_mean = mean(
          acc,
          na.rm = TRUE
        ),

        minimum_ESS = suppressWarnings(
          min(
            ess,
            na.rm = TRUE
          )
        ),

        median_ESS = suppressWarnings(
          median(
            ess,
            na.rm = TRUE
          )
        ),

        maximum_Rhat = suppressWarnings(
          max(
            rhat,
            na.rm = TRUE
          )
        )
      )
    }
  )
)


write_csv(
  DIAGNOSTIC_SUMMARY,
  file.path(
    MAIN_MCMC_DIR,
    "mcmc_diagnostic_summary.csv"
  )
)


## ============================================================
## 15) OPTIONAL EXACT REPRODUCTION CHECK
##
## If the frozen manuscript MCMC RDS is available in the
## reference folder, compare the six meta-regression parameters.
## ============================================================

FROZEN_REFERENCE_RDS <- file.path(
  BASELINE_DIR,
  "MCMC3_results_fullset_perID_priors_mixed.rds"
)


if (
  file.exists(
    FROZEN_REFERENCE_RDS
  )
) {

  FROZEN_MCMC <- readRDS(
    FROZEN_REFERENCE_RDS
  )


  compare_rows <- list()


  for (
    id in DEV_IDS
  ) {

    id_chr <- as.character(
      id
    )


    old_draws <- do.call(
      rbind,
      lapply(
        FROZEN_MCMC[[id_chr]]$chains,
        function(m) {
          as.data.frame(
            m$pars
          )
        }
      )
    )


    new_draws <- do.call(
      rbind,
      lapply(
        MCMC3_FINAL[[id_chr]]$chains,
        function(m) {
          as.data.frame(
            m$pars
          )
        }
      )
    )


    for (
      parameter in KEY_META_PARAMS
    ) {

      baseline_median <- median(
        exp(
          old_draws[[parameter]]
        ),
        na.rm = TRUE
      )


      reproduced_median <- median(
        exp(
          new_draws[[parameter]]
        ),
        na.rm = TRUE
      )


      compare_rows[[length(compare_rows) + 1]] <- tibble(
        ID = id,

        Parameter = parameter,

        baseline_median =
          baseline_median,

        reproduced_median =
          reproduced_median,

        ratio =
          reproduced_median /
          baseline_median,

        percent_difference =
          100 *
          (
            reproduced_median /
              baseline_median -
              1
          )
      )
    }
  }


  REPRO_CHECK <- bind_rows(
    compare_rows
  )


  write_csv(
    REPRO_CHECK,
    file.path(
      MAIN_MCMC_DIR,
      "manuscript_reproduction_check.csv"
    )
  )


  max_abs_diff <- max(
    abs(
      REPRO_CHECK$percent_difference
    ),
    na.rm = TRUE
  )


  cat(
    "\nMaximum absolute posterior-median difference ",
    "versus frozen manuscript RDS: ",
    format(
      max_abs_diff,
      scientific = TRUE
    ),
    "%\n",
    sep = ""
  )

} else {

  cat(
    "\nFrozen reference MCMC RDS not present; ",
    "exact posterior comparison skipped.\n",
    sep = ""
  )
}


## ============================================================
## 16) ANALYSIS MANIFEST
## ============================================================

MANIFEST <- tibble(
  Item = c(
    "analysis_version",
    "master_seed",
    "prior_center_strategy",
    "prior_weight",
    "pilot_iterations",
    "main_iterations_per_chain",
    "number_of_chains",
    "posterior_predictive_draws",
    "development_IDs",
    "source_deterministic_script",
    "source_MCMC_script"
  ),

  Value = c(
    "manuscript-compatible",
    as.character(
      MASTER_SEED
    ),
    "mixed",
    "0.75",
    "1000",
    "3000",
    "3",
    "400",
    paste(
      DEV_IDS,
      collapse = ","
    ),
    "reference/manuscript_baseline/01_run_model_fit.R",
    "reference/manuscript_baseline/03_run_mcmc_per_id.R"
  )
)


write_csv(
  MANIFEST,
  file.path(
    MAIN_MCMC_DIR,
    "mcmc_analysis_manifest.csv"
  )
)


## ============================================================
## 17) CLEAN ISOLATED WORKSPACE IF REQUESTED
## ============================================================

if (
  !KEEP_ISOLATED_WORKDIR
) {

  unlink(
    WORK_DIR,
    recursive = TRUE,
    force = TRUE
  )
}


## ============================================================
## 18) COMPLETION
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "MANUSCRIPT-COMPATIBLE MCMC COMPLETED\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nFinal diagnostic summary:\n"
)

print(
  DIAGNOSTIC_SUMMARY,
  n = Inf,
  width = Inf
)


cat(
  "\nCanonical MCMC:\n  ",
  file.path(
    MAIN_MCMC_DIR,
    "MCMC3_results_fullset_perID_priors_mixed.rds"
  ),
  "\n",
  sep = ""
)


cat(
  "\nStandardized per-ID outputs:\n  ",
  file.path(
    MAIN_MCMC_DIR,
    "ID1"
  ),
  " ... ID5\n",
  sep = ""
)


cat(
  "\nScientific settings:\n"
)

cat(
  "  Pilot iterations: 1000\n"
)

cat(
  "  Main iterations per chain: 3000\n"
)

cat(
  "  Chains per dataset: 3\n"
)

cat(
  "  Prior strategy: mixed\n"
)

cat(
  "  Prior weight: 0.75\n"
)

cat(
  "  Master seed: 1234\n"
)


cat(
  "\nNo 10,000-iteration refactored results were used.\n"
)

cat(
  "============================================================\n"
)