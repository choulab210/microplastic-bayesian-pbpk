## ============================================================
## 05_external_evaluation.R
##
## External evaluation of particle-size-informed PBPK
## predictions without parameter refitting.
##
## External datasets:
##   ID 6 = 50 nm
##   ID 7 = 500 nm
##   ID 8 = 5000 nm
##
## POINT PREDICTION:
##   Deterministic PBPK simulation using the central predicted
##   parameter vector.
##
## UNCERTAINTY:
##   1000 Monte Carlo parameter draws are propagated through
##   the PBPK model to obtain 95% parameter-propagated
##   prediction intervals.
##
## Inputs:
##   data/data.csv
##   model/MP_PBPK.R
##   results/deterministic/deterministic_outputs.rds
##   results/meta/meta_analysis_outputs.rds
##
## Outputs:
##   results/external/
##
## No external observations are used for parameter fitting.
## ============================================================


## ============================================================
## 0) Packages
## ============================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(tools)
})

options(stringsAsFactors = FALSE)


## ============================================================
## 1) Analysis settings
## ============================================================

SEED <- 12345
set.seed(SEED)

TEST_IDS <- 6:8

N_MC <- 1000

SIM_DT_H <- 0.05

META_PRED_TYPE <- "new_id"


DATA_FILE <- file.path(
  "data",
  "data.csv"
)

MODEL_FILE <- file.path(
  "model",
  "MP_PBPK.R"
)

DETERMINISTIC_FILE <- file.path(
  "results",
  "deterministic",
  "deterministic_outputs.rds"
)

META_FILE <- file.path(
  "results",
  "meta",
  "meta_analysis_outputs.rds"
)

RESULTS_DIR <- file.path(
  "results",
  "external"
)

dir.create(
  RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 2) Check files
## ============================================================

required_files <- c(
  DATA_FILE,
  MODEL_FILE,
  DETERMINISTIC_FILE,
  META_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


## ============================================================
## 3) Load analysis objects
## ============================================================

det <- readRDS(
  DETERMINISTIC_FILE
)

meta <- readRDS(
  META_FILE
)


if (is.null(det$theta_init)) {
  stop("theta_init missing from deterministic outputs.")
}

if (is.null(meta$meta_params)) {
  stop("meta_params missing from meta-analysis outputs.")
}

if (is.null(meta$meta_summary)) {
  stop("meta_summary missing from meta-analysis outputs.")
}

if (is.null(meta$global_summary)) {
  stop("global_summary missing from meta-analysis outputs.")
}


theta_init <- det$theta_init

META_PARAMS <- meta$meta_params

META_SUMMARY <- meta$meta_summary

GLOBAL_SUMMARY <- meta$global_summary


cat("\n============================================================\n")
cat("External PBPK evaluation\n")
cat("============================================================\n")

cat(
  "External IDs: ",
  paste(TEST_IDS, collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Monte Carlo simulations per dataset: ",
  N_MC,
  "\n",
  sep = ""
)

cat(
  "Meta-regression parameters: ",
  paste(META_PARAMS, collapse = ", "),
  "\n",
  sep = ""
)


## ============================================================
## 4) Load and harmonize data
## ============================================================

obs_raw <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_columns <- c(
  "ID",
  "Size_nm",
  "dose_mg",
  "Time_h",
  "Tissue",
  "Conc"
)

missing_columns <- setdiff(
  required_columns,
  names(obs_raw)
)

if (length(missing_columns) > 0) {

  stop(
    "data.csv is missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


if (!"SD" %in% names(obs_raw)) {
  obs_raw$SD <- NA_real_
}


harmonize_tissue <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x <- tolower(x)

  x <- dplyr::recode(
    x,

    "lung" = "lungs",

    "gi tract" = "gi",
    "g.i." = "gi",
    "gitract" = "gi",

    "stomach" = "gi",
    "intestine" = "gi",
    "colon" = "gi",

    .default = x
  )

  x <- tools::toTitleCase(x)

  x <- dplyr::recode(
    x,

    "Lung" = "Lungs",
    "Gi" = "GI",

    .default = x
  )

  x
}


EVALUATION_TISSUES <- c(
  "Blood",
  "Liver",
  "Kidney",
  "Spleen",
  "Lungs",
  "GI",
  "Urine"
)


obs_all <- obs_raw %>%
  mutate(

    ID = as.numeric(ID),

    Size_nm = as.numeric(Size_nm),

    dose_mg = as.numeric(dose_mg),

    Time_h = as.numeric(Time_h),

    Conc = as.numeric(Conc),

    SD = as.numeric(SD),

    Tissue = harmonize_tissue(
      Tissue
    ),

    SD = ifelse(
      is.na(SD),
      pmax(
        Conc * 0.1,
        1e-6
      ),
      SD
    )
  )


test_all <- obs_all %>%
  filter(
    ID %in% TEST_IDS
  )


if (nrow(test_all) == 0) {
  stop("No external observations found for IDs 6–8.")
}


## ------------------------------------------------------------
## Explicitly report unsupported endpoints
## ------------------------------------------------------------

excluded_external_tissues <- test_all %>%
  filter(
    !Tissue %in% EVALUATION_TISSUES
  ) %>%
  distinct(
    Tissue
  ) %>%
  pull(
    Tissue
  )


if (length(excluded_external_tissues) > 0) {

  cat(
    "\nExternal endpoints excluded because they are not ",
    "represented as explicit evaluation compartments:\n"
  )

  print(
    excluded_external_tissues
  )
}


test_df <- test_all %>%
  filter(
    Tissue %in% EVALUATION_TISSUES
  )


cat(
  "\nExternal observations retained: ",
  nrow(test_df),
  "\n",
  sep = ""
)


## ============================================================
## 5) External dataset information
## ============================================================

test_info <- test_df %>%
  group_by(ID) %>%
  summarise(

    Size_nm = first(
      Size_nm[
        is.finite(Size_nm)
      ]
    ),

    dose_mg = first(
      dose_mg[
        is.finite(dose_mg)
      ]
    ),

    t_end_h = max(
      Time_h,
      na.rm = TRUE
    ),

    n_observations = n(),

    n_times = n_distinct(
      Time_h
    ),

    .groups = "drop"
  )


cat("\nExternal dataset characteristics:\n")

print(
  test_info
)


## ============================================================
## 6) Compile PBPK model
## ============================================================

source(
  MODEL_FILE
)


if (!exists("mPBPK_RES")) {

  stop(
    "mPBPK_RES not found in ",
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
  "\nPBPK model compiled successfully. Parameters: ",
  length(MODEL_PARAM_NAMES),
  "\n",
  sep = ""
)


## ============================================================
## 7) Meta-regression distribution helper
##
## Returns log-scale mean and variance for a specified
## particle size.
## ============================================================

get_meta_distribution <- function(
    parameter,
    size_nm,
    type = c(
      "population",
      "new_id"
    )
) {

  type <- match.arg(
    type
  )


  row <- META_SUMMARY %>%
    filter(
      Parameter ==
        parameter
    )


  if (nrow(row) != 1) {

    stop(
      "Meta-regression row not found for: ",
      parameter
    )
  }


  beta0 <- row$beta0
  beta1 <- row$beta1


  Vbeta <- matrix(
    c(
      row$vc11,
      row$vc12,
      row$vc12,
      row$vc22
    ),
    nrow = 2,
    byrow = TRUE
  )


  x <- c(
    1,
    log(
      size_nm
    )
  )


  mu_log <- beta0 +
    beta1 *
    log(
      size_nm
    )


  var_log <- as.numeric(
    t(x) %*%
      Vbeta %*%
      x
  )


  if (type == "new_id") {

    var_log <- var_log +
      row$tau2
  }


  list(
    mu_log = as.numeric(
      mu_log
    ),

    var_log = as.numeric(
      pmax(
        var_log,
        0
      )
    )
  )
}


## ============================================================
## 8) Central meta-regression parameter value
##
## Central point prediction uses the expected value of the
## predicted parameter distribution:
##
##   E(theta) = exp(mu + 0.5 * variance)
##
## This reproduces the central parameter-vector logic of
## the previous manuscript workflow.
## ============================================================

central_meta_parameter <- function(
    parameter,
    size_nm,
    type = META_PRED_TYPE
) {

  dist <- get_meta_distribution(
    parameter = parameter,
    size_nm = size_nm,
    type = type
  )


  exp(
    dist$mu_log +
      0.5 *
      dist$var_log
  )
}


## ============================================================
## 9) Random meta-regression parameter draw
## ============================================================

sample_meta_parameter <- function(
    parameter,
    size_nm,
    type = META_PRED_TYPE
) {

  dist <- get_meta_distribution(
    parameter = parameter,
    size_nm = size_nm,
    type = type
  )


  exp(
    rnorm(
      1,
      mean = dist$mu_log,
      sd = sqrt(
        pmax(
          dist$var_log,
          1e-12
        )
      )
    )
  )
}


## ============================================================
## 10) Central pooled global parameter
##
## Use posterior mean on natural scale for deterministic
## central PBPK prediction.
## ============================================================

central_global_parameter <- function(
    parameter
) {

  row <- GLOBAL_SUMMARY %>%
    filter(
      Parameter ==
        parameter
    )


  if (nrow(row) != 1) {

    stop(
      "Global parameter summary not found for: ",
      parameter
    )
  }


  as.numeric(
    row$mean
  )
}


## ============================================================
## 11) Random pooled global parameter draw
##
## Global parameters are approximated using pooled posterior
## mean and SD on the log scale.
## ============================================================

sample_global_parameter <- function(
    parameter
) {

  row <- GLOBAL_SUMMARY %>%
    filter(
      Parameter ==
        parameter
    )


  if (nrow(row) != 1) {

    stop(
      "Global parameter summary not found for: ",
      parameter
    )
  }


  if (
    !is.finite(row$mean_log) ||
    !is.finite(row$sd_log)
  ) {

    return(
      as.numeric(
        row$mean
      )
    )
  }


  exp(
    rnorm(
      1,
      mean = row$mean_log,
      sd = pmax(
        row$sd_log,
        1e-8
      )
    )
  )
}


## ============================================================
## 12) Build deterministic central parameter vector
## ============================================================

build_central_parameter_vector <- function(
    size_nm
) {

  ## Start from Script 01 baseline values
  pars_nat <- exp(
    theta_init
  )


  ## ----------------------------------------------------------
  ## Six particle-size-informed parameters
  ## ----------------------------------------------------------

  for (parameter in META_PARAMS) {

    pars_nat[parameter] <-
      central_meta_parameter(
        parameter = parameter,
        size_nm = size_nm,
        type = META_PRED_TYPE
      )
  }


  ## ----------------------------------------------------------
  ## Remaining Bayesian parameters
  ## ----------------------------------------------------------

  if (
    !is.null(GLOBAL_SUMMARY) &&
    nrow(GLOBAL_SUMMARY) > 0
  ) {

    for (parameter in GLOBAL_SUMMARY$Parameter) {

      pars_nat[parameter] <-
        central_global_parameter(
          parameter
        )
    }
  }


  pars_nat
}


## ============================================================
## 13) Build one random parameter draw
## ============================================================

build_random_parameter_vector <- function(
    size_nm
) {

  pars_nat <- exp(
    theta_init
  )


  for (parameter in META_PARAMS) {

    pars_nat[parameter] <-
      sample_meta_parameter(
        parameter = parameter,
        size_nm = size_nm,
        type = META_PRED_TYPE
      )
  }


  if (
    !is.null(GLOBAL_SUMMARY) &&
    nrow(GLOBAL_SUMMARY) > 0
  ) {

    for (parameter in GLOBAL_SUMMARY$Parameter) {

      pars_nat[parameter] <-
        sample_global_parameter(
          parameter
        )
    }
  }


  pars_nat
}


## ============================================================
## 14) PBPK simulation helper
## ============================================================

simulate_pbpk <- function(
    pars_nat,
    dose_mg,
    t_end_h,
    observation_times,
    dt_h = SIM_DT_H,
    atol = 1e-8,
    rtol = 1e-6
) {

  pars_nat <- pars_nat[
    names(pars_nat) %in%
      MODEL_PARAM_NAMES
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
    delta = dt_h,
    add = observation_times
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
    "total_lung" %in%
      names(out)
  ) {

    out$total_lung

  } else if (
    all(
      c(
        "A_lung_b",
        "A_lung_t"
      ) %in%
        names(out)
    )
  ) {

    out$A_lung_b +
      out$A_lung_t

  } else {

    stop(
      "Unable to identify total lung amount."
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
      dose_mg * 100
  )
}


## ============================================================
## 15) Deterministic point prediction for one external ID
## ============================================================

run_central_prediction <- function(
    id,
    size_nm,
    dose_mg,
    t_end_h,
    observation_times
) {

  pars_central <- build_central_parameter_vector(
    size_nm = size_nm
  )


  simulation <- simulate_pbpk(
    pars_nat = pars_central,
    dose_mg = dose_mg,
    t_end_h = t_end_h,
    observation_times = observation_times,
    dt_h = SIM_DT_H
  ) %>%
    pivot_longer(
      cols = -time,
      names_to = "Tissue",
      values_to = "pred"
    ) %>%
    mutate(
      ID = id,
      Size_nm = size_nm,
      .before = 1
    )


  parameter_table <- tibble(
    ID = id,
    Size_nm = size_nm,
    Parameter = names(
      pars_central
    ),
    Value = as.numeric(
      pars_central
    )
  )


  list(
    simulation =
      simulation,

    parameters =
      parameter_table
  )
}


## ============================================================
## 16) Monte Carlo uncertainty propagation for one ID
## ============================================================

run_uncertainty_mc <- function(
    id,
    size_nm,
    dose_mg,
    t_end_h,
    observation_times,
    n_mc = N_MC
) {

  cat(
    "Running uncertainty propagation for ID ",
    id,
    " (",
    size_nm,
    " nm; N=",
    n_mc,
    ") ...\n",
    sep = ""
  )


  set.seed(
    SEED +
      id *
      1000
  )


  simulations <- vector(
    "list",
    n_mc
  )


  for (j in seq_len(n_mc)) {

    pars_j <- build_random_parameter_vector(
      size_nm =
        size_nm
    )


    simulations[[j]] <- simulate_pbpk(
      pars_nat =
        pars_j,

      dose_mg =
        dose_mg,

      t_end_h =
        t_end_h,

      observation_times =
        observation_times,

      dt_h =
        SIM_DT_H
    ) %>%
      pivot_longer(
        cols = -time,
        names_to = "Tissue",
        values_to = "pred"
      ) %>%
      mutate(
        ID = id,
        Size_nm = size_nm,
        Draw = j,
        .before = 1
      )
  }


  simulation_draws <- bind_rows(
    simulations
  )


  interval_summary <- simulation_draws %>%
    group_by(
      ID,
      Size_nm,
      time,
      Tissue
    ) %>%
    summarise(

      mc_mean = mean(
        pred,
        na.rm = TRUE
      ),

      mc_median = median(
        pred,
        na.rm = TRUE
      ),

      pred_lo = as.numeric(
        quantile(
          pred,
          0.025,
          na.rm = TRUE
        )
      ),

      pred_hi = as.numeric(
        quantile(
          pred,
          0.975,
          na.rm = TRUE
        )
      ),

      .groups = "drop"
    )


  list(
    draws =
      simulation_draws,

    interval =
      interval_summary
  )
}


## ============================================================
## 17) Run IDs 6–8
## ============================================================

CENTRAL_RESULTS <- setNames(
  vector(
    "list",
    length(TEST_IDS)
  ),
  as.character(TEST_IDS)
)


MC_RESULTS <- setNames(
  vector(
    "list",
    length(TEST_IDS)
  ),
  as.character(TEST_IDS)
)


for (id in TEST_IDS) {

  info_i <- test_info %>%
    filter(
      ID == id
    )


  if (nrow(info_i) != 1) {

    stop(
      "External dataset information not uniquely defined for ID ",
      id
    )
  }


  observation_times <- sort(
    unique(
      test_df$Time_h[
        test_df$ID == id
      ]
    )
  )


  cat(
    "\nID ",
    id,
    ": generating deterministic central prediction\n",
    sep = ""
  )


  CENTRAL_RESULTS[[as.character(id)]] <-
    run_central_prediction(
      id =
        id,

      size_nm =
        info_i$Size_nm,

      dose_mg =
        info_i$dose_mg,

      t_end_h =
        info_i$t_end_h,

      observation_times =
        observation_times
    )


  MC_RESULTS[[as.character(id)]] <-
    run_uncertainty_mc(
      id =
        id,

      size_nm =
        info_i$Size_nm,

      dose_mg =
        info_i$dose_mg,

      t_end_h =
        info_i$t_end_h,

      observation_times =
        observation_times,

      n_mc =
        N_MC
    )
}


## ============================================================
## 18) Combine deterministic predictions
## ============================================================

central_predictions <- bind_rows(
  lapply(
    CENTRAL_RESULTS,
    function(x) {
      x$simulation
    }
  )
)


central_parameters <- bind_rows(
  lapply(
    CENTRAL_RESULTS,
    function(x) {
      x$parameters
    }
  )
)


interval_predictions <- bind_rows(
  lapply(
    MC_RESULTS,
    function(x) {
      x$interval
    }
  )
)


write_csv(
  central_predictions,
  file.path(
    RESULTS_DIR,
    "external_central_predictions.csv"
  )
)


write_csv(
  central_parameters,
  file.path(
    RESULTS_DIR,
    "external_central_parameter_vectors.csv"
  )
)


write_csv(
  interval_predictions,
  file.path(
    RESULTS_DIR,
    "external_prediction_intervals.csv"
  )
)


## ============================================================
## 19) Match observations with central predictions + intervals
## ============================================================

obs_eval <- test_df %>%
  transmute(
    ID = as.numeric(ID),
    Size_nm = as.numeric(Size_nm),
    time = as.numeric(Time_h),
    Tissue = as.character(Tissue),
    obs = as.numeric(Conc),
    SD = as.numeric(SD)
  )


evaluation_df <- obs_eval %>%

  inner_join(
    central_predictions,
    by = c(
      "ID",
      "Size_nm",
      "time",
      "Tissue"
    )
  ) %>%

  left_join(
    interval_predictions,
    by = c(
      "ID",
      "Size_nm",
      "time",
      "Tissue"
    )
  )


cat(
  "\nMatched external observations: ",
  nrow(evaluation_df),
  "\n",
  sep = ""
)


if (nrow(evaluation_df) != nrow(obs_eval)) {

  warning(
    "Not all external observations matched predictions."
  )
}


## ============================================================
## 20) Performance metrics
##
## R2 and RMSE:
##   all finite observations
##
## MAPE and fold error:
##   positive observed + predicted values
##
## Coverage:
##   positive observations
##
## IMPORTANT:
##   Metrics use deterministic central PBPK prediction "pred",
##   not Monte Carlo mean.
## ============================================================

calculate_metrics <- function(
    df
) {

  finite_df <- df %>%
    filter(
      is.finite(obs),
      is.finite(pred)
    )


  positive_df <- finite_df %>%
    filter(
      obs > 0,
      pred > 0
    )


  ## ----------------------------------------------------------
  ## R2
  ## ----------------------------------------------------------

  r2 <- if (
    nrow(finite_df) > 1 &&
      length(
        unique(
          finite_df$obs
        )
      ) > 1 &&
      length(
        unique(
          finite_df$pred
        )
      ) > 1
  ) {

    cor(
      finite_df$obs,
      finite_df$pred,
      use = "complete.obs"
    )^2

  } else {

    NA_real_
  }


  ## ----------------------------------------------------------
  ## RMSE
  ## ----------------------------------------------------------

  rmse <- if (
    nrow(finite_df) > 0
  ) {

    sqrt(
      mean(
        (
          finite_df$pred -
            finite_df$obs
        )^2,
        na.rm = TRUE
      )
    )

  } else {

    NA_real_
  }


  ## ----------------------------------------------------------
  ## Positive-value metrics
  ## ----------------------------------------------------------

  if (nrow(positive_df) > 0) {

    ratio <- positive_df$pred /
      positive_df$obs


    mape <- mean(
      abs(
        positive_df$pred -
          positive_df$obs
      ) /
        positive_df$obs,
      na.rm = TRUE
    ) *
      100


    within_2fold_n <- sum(
      ratio >=
        0.5 &
        ratio <=
        2
    )


    within_2fold_pct <- mean(
      ratio >=
        0.5 &
        ratio <=
        2
    ) *
      100


    within_3fold_n <- sum(
      ratio >=
        1 / 3 &
        ratio <=
        3
    )


    within_3fold_pct <- mean(
      ratio >=
        1 / 3 &
        ratio <=
        3
    ) *
      100


  } else {

    mape <- NA_real_

    within_2fold_n <- 0
    within_2fold_pct <- NA_real_

    within_3fold_n <- 0
    within_3fold_pct <- NA_real_
  }


  ## ----------------------------------------------------------
  ## Interval coverage
  ## ----------------------------------------------------------

  coverage_df <- finite_df %>%
    filter(
      obs > 0,
      is.finite(pred_lo),
      is.finite(pred_hi)
    )


  if (nrow(coverage_df) > 0) {

    covered <- coverage_df$obs >=
      coverage_df$pred_lo &
      coverage_df$obs <=
      coverage_df$pred_hi


    coverage_n <- sum(
      covered
    )


    coverage_pct <- mean(
      covered
    ) *
      100


  } else {

    coverage_n <- 0
    coverage_pct <- NA_real_
  }


  tibble(

    n_all =
      nrow(
        finite_df
      ),

    n_positive =
      nrow(
        positive_df
      ),

    R2 =
      r2,

    RMSE =
      rmse,

    MAPE =
      mape,

    within_2fold_n =
      within_2fold_n,

    within_2fold_pct =
      within_2fold_pct,

    within_3fold_n =
      within_3fold_n,

    within_3fold_pct =
      within_3fold_pct,

    coverage_n =
      coverage_n,

    coverage_pct =
      coverage_pct
  )
}


## ============================================================
## 21) Overall metrics
## ============================================================

metrics_overall <- calculate_metrics(
  evaluation_df
)


cat("\n============================================================\n")
cat("External evaluation: overall\n")
cat("============================================================\n")

print(
  metrics_overall
)


## ============================================================
## 22) Metrics by ID
## ============================================================

metrics_by_id <- bind_rows(
  lapply(
    TEST_IDS,
    function(id) {

      bind_cols(
        tibble(
          ID = id
        ),

        calculate_metrics(
          evaluation_df %>%
            filter(
              ID == id
            )
        )
      )
    }
  )
)


cat("\nExternal evaluation by ID:\n")

print(
  metrics_by_id
)


## ============================================================
## 23) Metrics by tissue
## ============================================================

metrics_by_tissue <- bind_rows(
  lapply(
    sort(
      unique(
        evaluation_df$Tissue
      )
    ),
    function(tissue) {

      bind_cols(
        tibble(
          Tissue =
            tissue
        ),

        calculate_metrics(
          evaluation_df %>%
            filter(
              Tissue ==
                tissue
            )
        )
      )
    }
  )
)


cat("\nExternal evaluation by tissue:\n")

print(
  metrics_by_tissue
)


## ============================================================
## 24) Observation-level metrics
## ============================================================

observation_metrics <- evaluation_df %>%
  mutate(

    ratio = ifelse(
      obs > 0 &
        pred > 0,

      pred /
        obs,

      NA_real_
    ),

    absolute_error =
      abs(
        pred -
          obs
      ),

    relative_error_pct = ifelse(
      obs > 0,

      absolute_error /
        obs *
        100,

      NA_real_
    ),

    within_2fold = ifelse(
      is.finite(
        ratio
      ),

      ratio >=
        0.5 &
        ratio <=
        2,

      NA
    ),

    within_3fold = ifelse(
      is.finite(
        ratio
      ),

      ratio >=
        1 / 3 &
        ratio <=
        3,

      NA
    ),

    covered_95 = ifelse(
      obs > 0,

      obs >=
        pred_lo &
        obs <=
        pred_hi,

      NA
    )
  )


## ============================================================
## 25) Save metrics
## ============================================================

write_csv(
  evaluation_df,
  file.path(
    RESULTS_DIR,
    "external_observed_vs_predicted.csv"
  )
)


write_csv(
  observation_metrics,
  file.path(
    RESULTS_DIR,
    "external_observation_level_metrics.csv"
  )
)


write_csv(
  metrics_overall,
  file.path(
    RESULTS_DIR,
    "external_metrics_overall.csv"
  )
)


write_csv(
  metrics_by_id,
  file.path(
    RESULTS_DIR,
    "external_metrics_by_id.csv"
  )
)


write_csv(
  metrics_by_tissue,
  file.path(
    RESULTS_DIR,
    "external_metrics_by_tissue.csv"
  )
)


## ============================================================
## 26) Table S10 source
## ============================================================

TABLE_S10_OVERALL <- metrics_overall %>%
  mutate(
    Group =
      "Overall",
    .before =
      1
  )


TABLE_S10_BY_ID <- metrics_by_id %>%
  mutate(
    Group = paste0(
      "ID ",
      ID
    ),
    .before =
      1
  )


TABLE_S10 <- bind_rows(
  TABLE_S10_OVERALL,
  TABLE_S10_BY_ID
)


write_csv(
  TABLE_S10,
  file.path(
    RESULTS_DIR,
    "Table_S10_external_evaluation.csv"
  )
)


## ============================================================
## 27) Save consolidated R object
## ============================================================

external_outputs <- list(

  seed =
    SEED,

  test_ids =
    TEST_IDS,

  n_mc =
    N_MC,

  meta_prediction_type =
    META_PRED_TYPE,

  test_info =
    test_info,

  excluded_external_tissues =
    excluded_external_tissues,

  central_parameters =
    central_parameters,

  central_predictions =
    central_predictions,

  prediction_intervals =
    interval_predictions,

  evaluation_df =
    evaluation_df,

  observation_metrics =
    observation_metrics,

  metrics_overall =
    metrics_overall,

  metrics_by_id =
    metrics_by_id,

  metrics_by_tissue =
    metrics_by_tissue,

  table_s10 =
    TABLE_S10
)


saveRDS(
  external_outputs,
  file.path(
    RESULTS_DIR,
    "external_evaluation_outputs.rds"
  )
)


## ============================================================
## 28) Final summary
## ============================================================

cat("\n============================================================\n")
cat("External evaluation completed successfully\n")
cat("============================================================\n")

cat(
  "External IDs: ",
  paste(
    TEST_IDS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "External observations evaluated: ",
  nrow(
    evaluation_df
  ),
  "\n",
  sep = ""
)

cat(
  "Monte Carlo simulations per dataset: ",
  N_MC,
  "\n",
  sep = ""
)

cat(
  "Point prediction: deterministic PBPK using ",
  "central predicted parameter vector\n"
)

cat(
  "95% interval: 2.5th–97.5th percentiles of ",
  N_MC,
  " Monte Carlo PBPK simulations\n",
  sep = ""
)


cat("\nOverall performance:\n")

print(
  metrics_overall
)


cat(
  "\nTable S10 source: ",
  file.path(
    RESULTS_DIR,
    "Table_S10_external_evaluation.csv"
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