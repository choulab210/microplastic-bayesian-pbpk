## ============================================================
## 02_run_sensitivity.R
##
## Post-calibration AUC-based normalized sensitivity analysis
##
## Workflow:
##   1. Load deterministic outputs from Script 01
##   2. Recompile PBPK model
##   3. Simulate baseline fitted profiles for each development ID
##   4. Perturb each particle-specific parameter by +1%
##   5. Calculate tissue-specific AUCs
##   6. Calculate normalized sensitivity coefficients (NSCs)
##   7. Rank parameters by mean absolute NSC across:
##        - all modeled AUC endpoints
##        - all five development datasets
##   8. Select six highest-ranked parameters
##
## Development datasets: IDs 1–5
##
## Important:
##   - IDs 2 and 5 are treated independently for ranking
##   - 220-nm averaging is NOT performed here
##   - averaging IDs 2 and 5 is only for Figure 3 visualization
## ============================================================


## ============================================================
## 0) Packages
## ============================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

options(stringsAsFactors = FALSE)


## ============================================================
## 1) Analysis settings
## ============================================================

SEED <- 1234
set.seed(SEED)

NSC_FRAC_CHANGE <- 0.01
NSC_DT_H        <- 0.05
NSC_MIN_DENOM   <- 1e-10

N_PRIORITIZED_PARAMS <- 6

DETERMINISTIC_FILE <- file.path(
  "results",
  "deterministic",
  "deterministic_outputs.rds"
)

MODEL_FILE <- file.path(
  "model",
  "MP_PBPK.R"
)

RESULTS_DIR <- file.path(
  "results",
  "sensitivity"
)

dir.create(
  RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 2) Check required files
## ============================================================

if (!file.exists(DETERMINISTIC_FILE)) {
  stop(
    "Missing deterministic output file: ",
    DETERMINISTIC_FILE,
    "\nRun scripts/01_fit_model.R first."
  )
}

if (!file.exists(MODEL_FILE)) {
  stop(
    "Missing PBPK model file: ",
    MODEL_FILE
  )
}


## ============================================================
## 3) Load deterministic outputs
## ============================================================

det <- readRDS(
  DETERMINISTIC_FILE
)

required_objects <- c(
  "ids",
  "obs_df",
  "theta_fit_by_id",
  "theta_init"
)

missing_objects <- setdiff(
  required_objects,
  names(det)
)

if (length(missing_objects) > 0) {
  stop(
    "deterministic_outputs.rds is missing required objects: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}

ids <- as.numeric(det$ids)
obs_df <- det$obs_df
theta_fit_by_id <- det$theta_fit_by_id
theta_init <- det$theta_init


cat("\n============================================================\n")
cat("Loaded deterministic calibration outputs\n")
cat("============================================================\n")

cat(
  "Development IDs: ",
  paste(ids, collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Development observations: ",
  nrow(obs_df),
  "\n",
  sep = ""
)


## ============================================================
## 4) Load and compile PBPK model
## ============================================================

source(MODEL_FILE)

if (!exists("mPBPK_RES")) {
  stop(
    "Object 'mPBPK_RES' not found after sourcing ",
    MODEL_FILE
  )
}

mod <- mcode(
  "mPBPK_RES",
  mPBPK_RES
)

MODEL_PARAM_NAMES <- names(
  param(mod)
)

cat(
  "PBPK model compiled successfully. Parameters: ",
  length(MODEL_PARAM_NAMES),
  "\n",
  sep = ""
)


## ============================================================
## 5) PBPK prediction helper
## ============================================================

pred_sensitivity <- function(
    pars_log,
    id,
    t_end_h,
    dt_h = NSC_DT_H,
    atol = 1e-8,
    rtol = 1e-6
) {

  sub_df <- obs_df %>%
    filter(ID == id)

  if (nrow(sub_df) == 0) {
    stop(
      "No observations found for ID = ",
      id
    )
  }

  dose_values <- unique(
    as.numeric(
      sub_df$dose_mg
    )
  )

  dose_values <- dose_values[
    is.finite(
      dose_values
    )
  ]

  if (length(dose_values) != 1) {
    stop(
      "Expected one unique dose for ID ",
      id,
      "; found ",
      paste(
        dose_values,
        collapse = ", "
      )
    )
  }

  dose_mg <- dose_values[1]

  ## Convert log-scale parameters to natural scale
  pars_nat <- exp(
    pars_log
  )

  pars_nat <- pars_nat[
    names(pars_nat) %in% MODEL_PARAM_NAMES
  ]

  mod_tmp <- mod

  if (length(pars_nat) > 0) {
    mod_tmp <- param(
      mod_tmp,
      pars_nat
    )
  }

  ex <- ev(
    amt = dose_mg,
    cmt = "A_lumen",
    time = 0
  )

  ts <- tgrid(
    start = 0,
    end = t_end_h,
    delta = dt_h
  )

  out <- mrgsim_d(
    mod_tmp,
    data = ex,
    tgrid = ts,
    atol = atol,
    rtol = rtol,
    maxsteps = 1e6
  ) %>%
    as.data.frame()

  lungs_total <- if (
    "total_lung" %in% names(out)
  ) {

    out$total_lung

  } else if (
    all(
      c(
        "A_lung_b",
        "A_lung_t"
      ) %in% names(out)
    )
  ) {

    out$A_lung_b +
      out$A_lung_t

  } else {

    stop(
      "Unable to determine total lung amount."
    )
  }

  tibble(
    time = out$time,

    Blood =
      out$A_blood /
      dose_mg * 100,

    Liver =
      out$total_liver /
      dose_mg * 100,

    Kidney =
      out$total_kidney /
      dose_mg * 100,

    Spleen =
      out$total_spleen /
      dose_mg * 100,

    Lungs =
      lungs_total /
      dose_mg * 100,

    GI =
      out$total_GI /
      dose_mg * 100,

    Urine =
      out$A_urine /
      dose_mg * 100,

    Rest =
      out$total_rest /
      dose_mg * 100
  )
}


## ============================================================
## 6) Trapezoidal AUC helper
## ============================================================

trapz_auc <- function(
    x,
    y
) {

  ok <- is.finite(x) &
    is.finite(y)

  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2) {
    return(
      NA_real_
    )
  }

  ord <- order(x)

  x <- x[ord]
  y <- y[ord]

  sum(
    diff(x) *
      (
        head(
          y,
          -1
        ) +
        tail(
          y,
          -1
        )
      ) /
      2
  )
}


calc_auc_by_tissue <- function(
    pred_df
) {

  pred_df %>%
    pivot_longer(
      cols = -time,
      names_to = "Tissue",
      values_to = "Value"
    ) %>%
    group_by(
      Tissue
    ) %>%
    summarise(
      AUC = trapz_auc(
        time,
        Value
      ),
      .groups = "drop"
    )
}


## ============================================================
## 7) NSC calculation for one development dataset
## ============================================================

calc_nsc_for_id <- function(
    id,
    pars_log_base,
    param_names = names(pars_log_base),
    frac_change = NSC_FRAC_CHANGE,
    dt_h = NSC_DT_H
) {

  t_end <- max(
    obs_df$Time_h[
      obs_df$ID == id
    ],
    na.rm = TRUE
  )

  if (!is.finite(t_end)) {
    stop(
      "Unable to determine final observation time for ID ",
      id
    )
  }

  ## ----------------------------------------------------------
  ## Baseline simulation
  ## ----------------------------------------------------------

  base_pred <- pred_sensitivity(
    pars_log = pars_log_base,
    id = id,
    t_end_h = t_end,
    dt_h = dt_h
  )

  base_auc <- calc_auc_by_tissue(
    base_pred
  ) %>%
    mutate(
      ID = id,
      .before = 1
    )

  ## ----------------------------------------------------------
  ## Perturb each parameter independently
  ## ----------------------------------------------------------

  nsc_list <- lapply(
    param_names,
    function(param_name) {

      pars_pert <- pars_log_base

      ## A 1% increase on the natural scale corresponds to:
      ## log(theta_pert) = log(theta_base) + log(1.01)

      pars_pert[param_name] <-
        pars_log_base[param_name] +
        log(
          1 +
            frac_change
        )

      pert_pred <- tryCatch(
        pred_sensitivity(
          pars_log = pars_pert,
          id = id,
          t_end_h = t_end,
          dt_h = dt_h
        ),
        error = function(e) {

          warning(
            "Sensitivity simulation failed for ID ",
            id,
            ", parameter ",
            param_name,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
      )

      if (is.null(pert_pred)) {
        return(
          NULL
        )
      }

      pert_auc <- calc_auc_by_tissue(
        pert_pred
      )

      base_auc %>%
        select(
          ID,
          Tissue,
          AUC_base = AUC
        ) %>%
        left_join(
          pert_auc %>%
            rename(
              AUC_pert = AUC
            ),
          by = "Tissue"
        ) %>%
        mutate(
          Parameter = param_name,

          frac_param =
            frac_change,

          frac_auc =
            (
              AUC_pert -
                AUC_base
            ) /
            pmax(
              abs(
                AUC_base
              ),
              NSC_MIN_DENOM
            ),

          NSC =
            frac_auc /
            frac_change,

          absNSC =
            abs(
              NSC
            )
        ) %>%
        select(
          ID,
          Parameter,
          Tissue,
          AUC_base,
          AUC_pert,
          frac_param,
          frac_auc,
          NSC,
          absNSC
        )
    }
  )

  nsc_df <- bind_rows(
    nsc_list
  )

  list(
    ID = id,
    t_end_h = t_end,
    base_auc = base_auc,
    nsc = nsc_df
  )
}


## ============================================================
## 8) Run AUC-based NSC analysis for IDs 1–5
## ============================================================

cat("\n============================================================\n")
cat("Running AUC-based normalized sensitivity analysis\n")
cat("============================================================\n")

NSC_results_by_id <- setNames(
  vector(
    "list",
    length(ids)
  ),
  as.character(ids)
)

for (id in ids) {

  cat(
    "Running NSC for ID ",
    id,
    " ...\n",
    sep = ""
  )

  pars_i <- theta_fit_by_id[[as.character(id)]]

  NSC_results_by_id[[as.character(id)]] <-
    calc_nsc_for_id(
      id = id,
      pars_log_base = pars_i,
      param_names = names(pars_i),
      frac_change = NSC_FRAC_CHANGE,
      dt_h = NSC_DT_H
    )
}


## ============================================================
## 9) Combine AUC and NSC outputs
## ============================================================

base_auc_all <- bind_rows(
  lapply(
    NSC_results_by_id,
    function(x) {
      x$base_auc
    }
  )
)

nsc_all <- bind_rows(
  lapply(
    NSC_results_by_id,
    function(x) {
      x$nsc
    }
  )
)

cat(
  "\nBase AUC rows: ",
  nrow(base_auc_all),
  "\n",
  sep = ""
)

cat(
  "NSC rows: ",
  nrow(nsc_all),
  "\n",
  sep = ""
)


## ============================================================
## 10) Overall parameter ranking
##
## IMPORTANT:
## IDs 1–5 are treated independently.
##
## Mean |NSC| is calculated across:
##   parameter × tissue × dataset
##
## This is the ranking used for sensitivity prioritization.
## ============================================================

nsc_parameter_ranking <- nsc_all %>%
  filter(
    is.finite(
      absNSC
    )
  ) %>%
  group_by(
    Parameter
  ) %>%
  summarise(
    mean_abs_NSC = mean(
      absNSC,
      na.rm = TRUE
    ),

    median_abs_NSC = median(
      absNSC,
      na.rm = TRUE
    ),

    max_abs_NSC = max(
      absNSC,
      na.rm = TRUE
    ),

    min_abs_NSC = min(
      absNSC,
      na.rm = TRUE
    ),

    n_NSC = sum(
      is.finite(
        absNSC
      )
    ),

    n_datasets = n_distinct(
      ID
    ),

    n_tissues = n_distinct(
      Tissue
    ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      mean_abs_NSC
    )
  ) %>%
  mutate(
    Rank = row_number(),
    .before = 1
  )


## ============================================================
## 11) Select six sensitivity-prioritized parameters
## ============================================================

sensitivity_prioritized_parameters <- nsc_parameter_ranking %>%
  slice_head(
    n = N_PRIORITIZED_PARAMS
  ) %>%
  mutate(
    Prioritized = TRUE
  )

PRIORITIZED_PARAMS <-
  sensitivity_prioritized_parameters$Parameter


cat("\n============================================================\n")
cat("Overall NSC parameter ranking\n")
cat("============================================================\n")

print(
  nsc_parameter_ranking
)


cat("\n============================================================\n")
cat("Sensitivity-prioritized parameters\n")
cat("============================================================\n")

print(
  sensitivity_prioritized_parameters
)

cat(
  "\nSelected parameters: ",
  paste(
    PRIORITIZED_PARAMS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


## ============================================================
## 12) Expected manuscript parameters check
## ============================================================

EXPECTED_PRIORITIZED_PARAMS <- c(
  "Kfeces",
  "Ktrans",
  "Kgw2b",
  "Kshed",
  "KLRESmax",
  "KurineC"
)

if (!setequal(
  PRIORITIZED_PARAMS,
  EXPECTED_PRIORITIZED_PARAMS
)) {

  warning(
    "\nSensitivity-prioritized parameters differ from the ",
    "current manuscript baseline.\n",
    "Expected: ",
    paste(
      EXPECTED_PRIORITIZED_PARAMS,
      collapse = ", "
    ),
    "\nObserved: ",
    paste(
      PRIORITIZED_PARAMS,
      collapse = ", "
    )
  )

} else {

  cat(
    "\nPASS: prioritized parameter set matches manuscript baseline.\n"
  )
}


## ============================================================
## 13) Tissue-specific ranking summary
## ============================================================

nsc_by_tissue_summary <- nsc_all %>%
  filter(
    is.finite(
      absNSC
    )
  ) %>%
  group_by(
    Parameter,
    Tissue
  ) %>%
  summarise(
    mean_abs_NSC = mean(
      absNSC,
      na.rm = TRUE
    ),

    max_abs_NSC = max(
      absNSC,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  arrange(
    Parameter,
    desc(
      mean_abs_NSC
    )
  )


## ============================================================
## 14) Dataset-specific ranking summary
## ============================================================

nsc_by_id_summary <- nsc_all %>%
  filter(
    is.finite(
      absNSC
    )
  ) %>%
  group_by(
    ID,
    Parameter
  ) %>%
  summarise(
    mean_abs_NSC = mean(
      absNSC,
      na.rm = TRUE
    ),

    max_abs_NSC = max(
      absNSC,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  arrange(
    ID,
    desc(
      mean_abs_NSC
    )
  )


## ============================================================
## 15) Save CSV outputs
## ============================================================

write.csv(
  base_auc_all,
  file.path(
    RESULTS_DIR,
    "base_auc_by_id.csv"
  ),
  row.names = FALSE
)

write.csv(
  nsc_all,
  file.path(
    RESULTS_DIR,
    "nsc_by_id_tissue_parameter.csv"
  ),
  row.names = FALSE
)

write.csv(
  nsc_parameter_ranking,
  file.path(
    RESULTS_DIR,
    "nsc_parameter_ranking.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity_prioritized_parameters,
  file.path(
    RESULTS_DIR,
    "sensitivity_prioritized_parameters.csv"
  ),
  row.names = FALSE
)

write.csv(
  nsc_by_tissue_summary,
  file.path(
    RESULTS_DIR,
    "nsc_by_tissue_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  nsc_by_id_summary,
  file.path(
    RESULTS_DIR,
    "nsc_by_id_summary.csv"
  ),
  row.names = FALSE
)


## ============================================================
## 16) Backward-compatible outputs
##
## These names match the older manuscript workflow.
## Useful for existing artwork scripts or manual checks.
## ============================================================

write.csv(
  nsc_all,
  file.path(
    RESULTS_DIR,
    "NSC_amount_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  base_auc_all,
  file.path(
    RESULTS_DIR,
    "base_auc_all.csv"
  ),
  row.names = FALSE
)


## ============================================================
## 17) Save full sensitivity object
## ============================================================

sensitivity_outputs <- list(
  seed = SEED,

  ids = ids,

  frac_change = NSC_FRAC_CHANGE,

  dt_h = NSC_DT_H,

  min_denom = NSC_MIN_DENOM,

  n_prioritized = N_PRIORITIZED_PARAMS,

  prioritized_params = PRIORITIZED_PARAMS,

  expected_prioritized_params =
    EXPECTED_PRIORITIZED_PARAMS,

  base_auc_all = base_auc_all,

  nsc_all = nsc_all,

  nsc_parameter_ranking =
    nsc_parameter_ranking,

  sensitivity_prioritized_parameters =
    sensitivity_prioritized_parameters,

  nsc_by_tissue_summary =
    nsc_by_tissue_summary,

  nsc_by_id_summary =
    nsc_by_id_summary,

  NSC_results_by_id =
    NSC_results_by_id
)

saveRDS(
  sensitivity_outputs,
  file = file.path(
    RESULTS_DIR,
    "sensitivity_outputs.rds"
  )
)


## ============================================================
## 18) Final summary
## ============================================================

cat("\n============================================================\n")
cat("AUC-based sensitivity analysis completed successfully\n")
cat("============================================================\n")

cat(
  "Development IDs: ",
  paste(
    ids,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Parameters evaluated: ",
  dplyr::n_distinct(
    nsc_all$Parameter
  ),
  "\n",
  sep = ""
)

cat(
  "AUC endpoints evaluated: ",
  paste(
    sort(
      unique(
        nsc_all$Tissue
      )
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Prioritized parameters: ",
  paste(
    PRIORITIZED_PARAMS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Outputs saved to: ",
  RESULTS_DIR,
  "\n",
  sep = ""
)

cat("============================================================\n")