## ============================================================
## 05_external_evaluation.R
##
## FINAL MANUSCRIPT-COMPATIBLE EXTERNAL PBPK EVALUATION
##
## External evaluation datasets:
##   ID 6 =   50 nm
##   ID 7 =  500 nm
##   ID 8 = 5000 nm
##
## IMPORTANT
## ---------
## - No refitting to external data
## - Uses manuscript-compatible MCMC/meta-regression outputs
## - Uses new-ID random-effects prediction
## - Fixed Monte Carlo seed for reproducibility
## - N = 1000 uncertainty-propagation simulations per dataset
##
## Reported point-performance metrics:
##   Observed vs Monte Carlo mean across ALL 18 observations
##
## Parameter-propagated 95% interval coverage:
##   Positive observations only (n = 14)
##
## Expected reproducible external result:
##   R2       ~ 0.7856
##   RMSE     ~ 3.022 %ID
##   Coverage = 10 / 14 = 71.4%
##
## NOTE
## ----
## The 95% intervals represent propagated parameter uncertainty.
## They do not include an additional residual-error model.
## ============================================================


suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(tools)
})

options(stringsAsFactors = FALSE)


## ============================================================
## 0) REPRODUCIBILITY / USER SETTINGS
## ============================================================

MASTER_SEED <- 1234

set.seed(
  MASTER_SEED
)


TEST_IDS <- c(
  6,
  7,
  8
)


stopifnot(
  identical(
    sort(
      as.integer(
        TEST_IDS
      )
    ),
    c(
      6L,
      7L,
      8L
    )
  )
)


META_PRED_TYPE <- "new_id"

META_STATISTIC <- "mean"

N_SIM_UNC <- 1000

SIM_DT_H <- 0.05

USE_LOG_Y <- FALSE

Y_FLOOR <- 1e-8


## ============================================================
## 1) DIRECTORIES
## ============================================================

DATA_FILE <- file.path(
  "data",
  "data.csv"
)

MODEL_FILE <- file.path(
  "model",
  "MP_PBPK.R"
)

META_FILE <- file.path(
  "results",
  "meta",
  "meta_analysis_outputs.rds"
)


META_DIR <- file.path(
  "results",
  "meta"
)

EXTERNAL_DIR <- file.path(
  "results",
  "external"
)

FIGURE_DIR <- file.path(
  "figures",
  "external"
)


dir.create(
  EXTERNAL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


required_files <- c(
  DATA_FILE,
  MODEL_FILE,
  META_FILE
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
    )
  )
}


## ============================================================
## 2) LOAD META-REGRESSION OUTPUTS
## ============================================================

META_OBJ <- readRDS(
  META_FILE
)


required_meta_objects <- c(
  "SIZE_PARAMS",
  "GLOBAL_SUM",
  "META_FITS"
)


missing_meta_objects <- setdiff(
  required_meta_objects,
  names(
    META_OBJ
  )
)


if (
  length(
    missing_meta_objects
  ) > 0
) {

  stop(
    "meta_analysis_outputs.rds missing: ",
    paste(
      missing_meta_objects,
      collapse = ", "
    )
  )
}


SIZE_PARAMS <- META_OBJ$SIZE_PARAMS

GLOBAL_SUM <- META_OBJ$GLOBAL_SUM

META_FITS <- META_OBJ$META_FITS


cat(
  "\n============================================================\n"
)

cat(
  "EXTERNAL PBPK EVALUATION\n"
)

cat(
  "============================================================\n"
)

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
  "Monte Carlo simulations per dataset: ",
  N_SIM_UNC,
  "\n",
  sep = ""
)

cat(
  "Meta-regression parameters: ",
  paste(
    SIZE_PARAMS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Master seed: ",
  MASTER_SEED,
  "\n",
  sep = ""
)


## ============================================================
## 3) LOAD EXTERNAL DATA
## ============================================================

harmonize_tissue <- function(x) {

  y <- trimws(
    as.character(x)
  )


  yl <- tolower(
    y
  )


  yl <- dplyr::recode(
    yl,

    "lung" =
      "lungs",

    "gi tract" =
      "gi",

    "g.i." =
      "gi",

    "gitract" =
      "gi",

    "intestine" =
      "gi",

    "colon" =
      "gi",

    .default =
      yl
  )


  out <- tools::toTitleCase(
    yl
  )


  dplyr::recode(
    out,

    "Lung" =
      "Lungs",

    "Gi" =
      "GI",

    .default =
      out
  )
}


ALLOWED_EXTERNAL_TISSUES <- c(
  "Blood",
  "Liver",
  "Kidney",
  "Spleen",
  "Lungs",
  "GI",
  "Urine"
)


OBS_ALL <- read.csv(
  DATA_FILE,
  encoding = "latin1",
  check.names = FALSE
) %>%

  mutate(
    ID =
      as.numeric(
        ID
      ),

    Size_nm =
      as.numeric(
        Size_nm
      ),

    Time_h =
      as.numeric(
        Time_h
      ),

    dose_mg =
      as.numeric(
        dose_mg
      ),

    Conc =
      as.numeric(
        Conc
      ),

    SD =
      as.numeric(
        SD
      ),

    Tissue =
      harmonize_tissue(
        Tissue
      ),

    SD =
      ifelse(
        is.na(
          SD
        ),

        pmax(
          Conc *
            0.1,
          1e-6
        ),

        SD
      )
  )


## ------------------------------------------------------------
## Explicitly identify endpoints excluded from evaluation
## ------------------------------------------------------------

EXCLUDED_ENDPOINTS <- OBS_ALL %>%

  filter(
    ID %in%
      TEST_IDS,

    !Tissue %in%
      ALLOWED_EXTERNAL_TISSUES
  ) %>%

  distinct(
    Tissue
  ) %>%

  pull(
    Tissue
  )


if (
  length(
    EXCLUDED_ENDPOINTS
  ) > 0
) {

  cat(
    "\nExternal endpoints excluded because they are not ",
    "explicit evaluation compartments:\n",
    sep = ""
  )

  print(
    EXCLUDED_ENDPOINTS
  )
}


TEST_DF <- OBS_ALL %>%

  filter(
    ID %in%
      TEST_IDS,

    Tissue %in%
      ALLOWED_EXTERNAL_TISSUES
  )


if (
  nrow(
    TEST_DF
  ) == 0
) {

  stop(
    "No external observations found for IDs 6–8."
  )
}


## ============================================================
## 4) CHARACTERIZE EXTERNAL DATASETS
## ============================================================

TEST_SIZES_TBL <- TEST_DF %>%

  group_by(
    ID
  ) %>%

  summarise(
    Size_nm =
      first(
        na.omit(
          Size_nm
        )
      ),

    dose_mg =
      first(
        na.omit(
          dose_mg
        )
      ),

    t_end_h =
      max(
        Time_h,
        na.rm = TRUE
      ),

    n_observations =
      n(),

    n_times =
      n_distinct(
        Time_h
      ),

    .groups =
      "drop"
  ) %>%

  arrange(
    ID
  )


if (
  !identical(
    TEST_SIZES_TBL$ID,
    as.numeric(
      TEST_IDS
    )
  )
) {

  stop(
    "External dataset IDs do not exactly match 6, 7, and 8."
  )
}


if (
  any(
    TEST_SIZES_TBL$n_times !=
      1
  )
) {

  stop(
    "IDs 6–8 are expected to be single-time-point external datasets."
  )
}


cat(
  "\nExternal dataset characteristics:\n"
)

print(
  TEST_SIZES_TBL,
  n = Inf
)


cat(
  "\nExternal observations retained: ",
  nrow(
    TEST_DF
  ),
  "\n",
  sep = ""
)


## ============================================================
## 5) LOAD / COMPILE PBPK MODEL
## ============================================================

source(
  MODEL_FILE
)


MOD <- mcode(
  "mPBPK_RES",
  mPBPK_RES
)


MODEL_PARAM_NAMES <- names(
  param(
    MOD
  )
)


cat(
  "\nPBPK model compiled successfully. Parameters: ",
  length(
    MODEL_PARAM_NAMES
  ),
  "\n",
  sep = ""
)


## ============================================================
## 6) MANUSCRIPT-ERA BASELINE PARAMETER VECTOR
##
## This is the theta.int used by the original deterministic /
## external-analysis workflow.
##
## Size-informed and pooled MCMC parameters are overwritten
## below. Parameters not represented in MCMC/meta-regression
## retain these baseline values.
## ============================================================

THETA_BASE_LOG <- log(
  c(
    KLRESmax =
      0.4,

    KLRESrelease =
      7.5e-3,

    KSRESmax =
      1,

    KSRESrelease =
      3e-3,

    KKRESmax =
      0.1,

    KKRESrelease =
      1e-2,

    KLURESmax =
      0.1,

    KLURESrelease =
      5e-3,

    KRRESmax =
      0.2,

    KRRESrelease =
      2e-3,

    PAL =
      1e-3,

    PAK =
      1e-3,

    PAS =
      1e-3,

    PALU =
      1e-3,

    PAGI =
      1e-3,

    PAR =
      1e-6,

    KpL =
      0.01,

    KpK =
      0.15,

    KpS =
      0.15,

    KpR =
      0.15,

    KpLU =
      0.15,

    KbileC =
      2.2e-3,

    KurineC =
      3.2e-4,

    Kfeces =
      0.53,

    Ktrans =
      5e-1,

    Kgw2b =
      3e-4,

    Kshed =
      5e-1
  )
)


THETA_BASE_LOG <- THETA_BASE_LOG[
  names(
    THETA_BASE_LOG
  ) %in%
    MODEL_PARAM_NAMES
]


## ============================================================
## 7) META-REGRESSION PARAMETER HELPERS
## ============================================================

get_meta_row <- function(parameter) {

  row <- META_FITS %>%

    filter(
      .data$param ==
        parameter
    )


  if (
    nrow(
      row
    ) != 1
  ) {

    stop(
      "Meta fit not found for parameter: ",
      parameter
    )
  }


  row
}


predict_param_by_size <- function(
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


  row <- get_meta_row(
    parameter
  )


  b0 <- as.numeric(
    row$coef0
  )

  b1 <- as.numeric(
    row$coef1
  )


  Vb <- matrix(
    c(
      row$vc11,
      row$vc12,
      row$vc12,
      row$vc22
    ),
    nrow = 2,
    ncol = 2
  )


  tau2 <- as.numeric(
    row$tau2
  )


  x <- c(
    1,
    log(
      size_nm
    )
  )


  mu <- as.numeric(
    b0 +
      b1 *
      log(
        size_nm
      )
  )


  var_mu <- as.numeric(
    t(x) %*%
      Vb %*%
      x
  )


  var_mu <- pmax(
    var_mu,
    0
  )


  var_total <- var_mu


  if (
    type ==
      "new_id"
  ) {

    var_total <-
      var_total +
      pmax(
        tau2,
        0
      )
  }


  tibble(
    Parameter =
      parameter,

    Size_nm =
      size_nm,

    type =
      type,

    mu_log =
      mu,

    var_mu_log =
      var_mu,

    tau2 =
      tau2,

    var_total_log =
      var_total,

    mean =
      exp(
        mu +
          0.5 *
          var_total
      ),

    median =
      exp(
        mu
      ),

    q025 =
      exp(
        mu -
          1.96 *
          sqrt(
            var_total
          )
      ),

    q975 =
      exp(
        mu +
          1.96 *
          sqrt(
            var_total
          )
      )
  )
}


predict_params_by_size <- function(
    parameters,
    size_nm,
    type = c(
      "population",
      "new_id"
    )
) {

  type <- match.arg(
    type
  )


  bind_rows(
    lapply(
      parameters,
      function(parameter) {

        predict_param_by_size(
          parameter =
            parameter,

          size_nm =
            size_nm,

          type =
            type
        )
      }
    )
  )
}


## ============================================================
## 8) CENTRAL PARAMETER VECTOR
##
## This preserves the original manuscript external-analysis
## definition:
##
## six meta parameters:
##   natural-scale new-ID mean
##
## pooled global parameters:
##   pooled natural-scale posterior mean
## ============================================================

build_central_param_vector <- function(
    size_nm
) {

  pars_nat <- exp(
    THETA_BASE_LOG
  )


  SZ <- predict_params_by_size(
    parameters =
      SIZE_PARAMS,

    size_nm =
      size_nm,

    type =
      META_PRED_TYPE
  )


  pars_nat[
    SZ$Parameter
  ] <- SZ[[META_STATISTIC]]


  if (
    nrow(
      GLOBAL_SUM
    ) > 0
  ) {

    pars_nat[
      GLOBAL_SUM$Parameter
    ] <- GLOBAL_SUM[[META_STATISTIC]]
  }


  pars_nat
}


## ============================================================
## 9) UNCERTAINTY SAMPLING
##
## Important:
## These functions intentionally preserve the original
## manuscript Script-05 Monte Carlo parameter-sampling logic.
## ============================================================

sample_size_param_one <- function(
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


  row <- get_meta_row(
    parameter
  )


  b0 <- as.numeric(
    row$coef0
  )

  b1 <- as.numeric(
    row$coef1
  )


  Vb <- matrix(
    c(
      row$vc11,
      row$vc12,
      row$vc12,
      row$vc22
    ),
    nrow = 2,
    ncol = 2
  )


  tau2 <- as.numeric(
    row$tau2
  )


  x <- c(
    1,
    log(
      size_nm
    )
  )


  mu <- as.numeric(
    b0 +
      b1 *
      log(
        size_nm
      )
  )


  var_mu <- as.numeric(
    t(x) %*%
      Vb %*%
      x
  )


  if (
    type ==
      "new_id"
  ) {

    var_total <-
      var_mu +
      tau2

  } else {

    var_total <-
      var_mu
  }


  exp(
    rnorm(
      1,

      mean =
        mu,

      sd =
        sqrt(
          pmax(
            var_total,
            1e-12
          )
        )
    )
  )
}


sample_global_param_one <- function(
    parameter
) {

  row <- GLOBAL_SUM %>%

    filter(
      .data$Parameter ==
        parameter
    )


  if (
    nrow(
      row
    ) != 1
  ) {

    stop(
      "Global parameter summary not found for: ",
      parameter
    )
  }


  ## ----------------------------------------------------------
  ## Preserve manuscript-era reconstruction of a lognormal
  ## distribution from natural-scale mean and 95% interval.
  ## ----------------------------------------------------------

  m <- as.numeric(
    row$mean
  )


  ql <- as.numeric(
    row$q025
  )


  qu <- as.numeric(
    row$q975
  )


  if (
    !is.finite(
      ql
    ) ||
      !is.finite(
        qu
      ) ||
      ql <= 0 ||
      qu <= 0 ||
      qu <= ql
  ) {

    return(
      m
    )
  }


  sdlog <- (
    log(
      qu
    ) -
      log(
        ql
      )
  ) /
    (
      2 *
        1.96
    )


  sdlog <- pmax(
    sdlog,
    1e-8
  )


  mulog <- log(
    pmax(
      m,
      1e-12
    )
  ) -
    0.5 *
    sdlog^2


  exp(
    rnorm(
      1,
      mean =
        mulog,
      sd =
        sdlog
    )
  )
}


build_parameter_draw <- function(
    size_nm
) {

  pars_nat <- exp(
    THETA_BASE_LOG
  )


  for (
    parameter in SIZE_PARAMS
  ) {

    pars_nat[[parameter]] <-
      sample_size_param_one(
        parameter =
          parameter,

        size_nm =
          size_nm,

        type =
          META_PRED_TYPE
      )
  }


  if (
    nrow(
      GLOBAL_SUM
    ) > 0
  ) {

    for (
      parameter in GLOBAL_SUM$Parameter
    ) {

      pars_nat[[parameter]] <-
        sample_global_param_one(
          parameter
        )
    }
  }


  pars_nat
}


## ============================================================
## 10) PBPK SIMULATION HELPER
## ============================================================

simulate_pbpk <- function(
    pars_nat,
    dose_mg,
    t_end_h,
    add_times = NULL
) {

  pars_nat <- pars_nat[
    names(
      pars_nat
    ) %in%
      MODEL_PARAM_NAMES
  ]


  mod_tmp <- MOD


  if (
    length(
      pars_nat
    ) > 0
  ) {

    mod_tmp <- param(
      mod_tmp,
      pars_nat
    )
  }


  event <- ev(
    amt =
      dose_mg,

    cmt =
      "A_lumen",

    time =
      0
  )


  grid <- tgrid(
    start =
      0,

    end =
      t_end_h,

    delta =
      SIM_DT_H,

    add =
      add_times
  )


  out <- mrgsim_d(
    mod_tmp,

    data =
      event,

    tgrid =
      grid,

    atol =
      1e-8,

    rtol =
      1e-6,

    maxsteps =
      1e6
  ) %>%

    as.data.frame()


  lungs_total <- if (
    "total_lung" %in%
      names(
        out
      )
  ) {

    out$total_lung

  } else if (
    all(
      c(
        "A_lung_b",
        "A_lung_t"
      ) %in%
        names(
          out
        )
    )
  ) {

    out$A_lung_b +
      out$A_lung_t

  } else {

    0
  }


  tibble(
    time =
      out$time,

    Blood =
      (
        out$A_blood /
          dose_mg
      ) *
      100,

    Liver =
      (
        out$total_liver /
          dose_mg
      ) *
      100,

    Kidney =
      (
        out$total_kidney /
          dose_mg
      ) *
      100,

    Spleen =
      (
        out$total_spleen /
          dose_mg
      ) *
      100,

    Lungs =
      (
        lungs_total /
          dose_mg
      ) *
      100,

    GI =
      (
        out$total_GI /
          dose_mg
      ) *
      100,

    Urine =
      (
        out$A_urine /
          dose_mg
      ) *
      100,

    Rest =
      (
        out$total_rest /
          dose_mg
      ) *
      100
  )
}


## ============================================================
## 11) DETERMINISTIC CENTRAL PREDICTIONS
## ============================================================

CENTRAL_SIM_ROWS <- list()

CENTRAL_VECTOR_ROWS <- list()


for (
  i in seq_len(
    nrow(
      TEST_SIZES_TBL
    )
  )
) {

  row_i <- TEST_SIZES_TBL[
    i,
    ,
    drop = FALSE
  ]


  id_i <- row_i$ID

  size_i <- row_i$Size_nm

  dose_i <- row_i$dose_mg

  tend_i <- row_i$t_end_h


  add_times_i <- sort(
    unique(
      TEST_DF$Time_h[
        TEST_DF$ID ==
          id_i
      ]
    )
  )


  cat(
    "\nID ",
    id_i,
    ": generating deterministic central prediction\n",
    sep = ""
  )


  pars_i <- build_central_param_vector(
    size_nm =
      size_i
  )


  CENTRAL_VECTOR_ROWS[[length(CENTRAL_VECTOR_ROWS) + 1]] <- tibble(
    ID =
      id_i,

    Size_nm =
      size_i,

    Parameter =
      names(
        pars_i
      ),

    Value =
      as.numeric(
        pars_i
      )
  )


  sim_i <- simulate_pbpk(
    pars_nat =
      pars_i,

    dose_mg =
      dose_i,

    t_end_h =
      tend_i,

    add_times =
      add_times_i
  ) %>%

    mutate(
      ID =
        id_i,

      Size_nm =
        size_i,

      dose_mg =
        dose_i,

      .before =
        1
    )


  CENTRAL_SIM_ROWS[[length(CENTRAL_SIM_ROWS) + 1]] <- sim_i
}


CENTRAL_SIM <- bind_rows(
  CENTRAL_SIM_ROWS
)


CENTRAL_PARAMETER_VECTORS <- bind_rows(
  CENTRAL_VECTOR_ROWS
)


write_csv(
  CENTRAL_SIM,
  file.path(
    EXTERNAL_DIR,
    "external_central_predictions.csv"
  )
)


write_csv(
  CENTRAL_PARAMETER_VECTORS,
  file.path(
    EXTERNAL_DIR,
    "external_central_parameter_vectors.csv"
  )
)


## ============================================================
## 12) MONTE CARLO PARAMETER-UNCERTAINTY PROPAGATION
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Running uncertainty propagation\n"
)

cat(
  "============================================================\n"
)


INTERVAL_ROWS <- list()


for (
  i in seq_len(
    nrow(
      TEST_SIZES_TBL
    )
  )
) {

  row_i <- TEST_SIZES_TBL[
    i,
    ,
    drop = FALSE
  ]


  id_i <- row_i$ID

  size_i <- row_i$Size_nm

  dose_i <- row_i$dose_mg

  tend_i <- row_i$t_end_h


  add_times_i <- sort(
    unique(
      TEST_DF$Time_h[
        TEST_DF$ID ==
          id_i
      ]
    )
  )


  cat(
    "Running uncertainty propagation for ID ",
    id_i,
    " (",
    size_i,
    " nm; N=",
    N_SIM_UNC,
    ") ...\n",
    sep = ""
  )


  DRAW_SIM_ROWS <- vector(
    "list",
    N_SIM_UNC
  )


  for (
    draw_i in seq_len(
      N_SIM_UNC
    )
  ) {

    pars_draw <- build_parameter_draw(
      size_nm =
        size_i
    )


    sim_draw <- simulate_pbpk(
      pars_nat =
        pars_draw,

      dose_mg =
        dose_i,

      t_end_h =
        tend_i,

      add_times =
        add_times_i
    ) %>%

      pivot_longer(
        cols =
          all_of(
            ALLOWED_EXTERNAL_TISSUES
          ),

        names_to =
          "Tissue",

        values_to =
          "prediction"
      ) %>%

      mutate(
        ID =
          id_i,

        Size_nm =
          size_i,

        draw =
          draw_i,

        .before =
          1
      )


    DRAW_SIM_ROWS[[draw_i]] <- sim_draw
  }


  DRAW_SIM <- bind_rows(
    DRAW_SIM_ROWS
  )


  interval_i <- DRAW_SIM %>%

    group_by(
      ID,
      Size_nm,
      time,
      Tissue
    ) %>%

    summarise(
      mc_mean =
        mean(
          prediction,
          na.rm = TRUE
        ),

      mc_median =
        median(
          prediction,
          na.rm = TRUE
        ),

      pred_lo =
        as.numeric(
          quantile(
            prediction,
            0.025,
            na.rm = TRUE
          )
        ),

      pred_hi =
        as.numeric(
          quantile(
            prediction,
            0.975,
            na.rm = TRUE
          )
        ),

      .groups =
        "drop"
    )


  INTERVAL_ROWS[[length(INTERVAL_ROWS) + 1]] <- interval_i
}


PRED_INTERVALS <- bind_rows(
  INTERVAL_ROWS
)


write_csv(
  PRED_INTERVALS,
  file.path(
    EXTERNAL_DIR,
    "external_prediction_intervals.csv"
  )
)


## ============================================================
## 13) MATCH CENTRAL / MONTE CARLO PREDICTIONS TO OBSERVATIONS
## ============================================================

OBS_EXTERNAL <- TEST_DF %>%

  transmute(
    ID =
      ID,

    Size_nm =
      Size_nm,

    time =
      Time_h,

    Tissue =
      Tissue,

    obs =
      Conc,

    SD =
      SD
  )


CENTRAL_LONG <- CENTRAL_SIM %>%

  pivot_longer(
    cols =
      all_of(
        ALLOWED_EXTERNAL_TISSUES
      ),

    names_to =
      "Tissue",

    values_to =
      "pred"
  ) %>%

  transmute(
    ID =
      ID,

    Size_nm =
      Size_nm,

    time =
      time,

    Tissue =
      Tissue,

    pred =
      pred
  )


EXTERNAL_COMPARE <- OBS_EXTERNAL %>%

  inner_join(
    CENTRAL_LONG,

    by =
      c(
        "ID",
        "Size_nm",
        "time",
        "Tissue"
      )
  ) %>%

  left_join(
    PRED_INTERVALS,

    by =
      c(
        "ID",
        "Size_nm",
        "time",
        "Tissue"
      )
  ) %>%

  arrange(
    ID,
    Tissue
  )


cat(
  "\nMatched external observations: ",
  nrow(
    EXTERNAL_COMPARE
  ),
  "\n",
  sep = ""
)


if (
  nrow(
    EXTERNAL_COMPARE
  ) !=
    18
) {

  warning(
    "Expected 18 external observations but matched ",
    nrow(
      EXTERNAL_COMPARE
    ),
    "."
  )
}


## ============================================================
## 14) FINAL PERFORMANCE METRICS
##
## Point-performance:
##   Observed vs Monte Carlo mean
##   ALL finite observations
##
## Coverage:
##   Positive observations only
## ============================================================

calc_external_metrics <- function(df) {

  all_df <- df %>%

    filter(
      is.finite(
        obs
      ),

      is.finite(
        mc_mean
      )
    )


  positive_df <- all_df %>%

    filter(
      obs > 0
    )


  R2 <- if (
    nrow(
      all_df
    ) > 1 &&
      length(
        unique(
          all_df$obs
        )
      ) > 1 &&
      length(
        unique(
          all_df$mc_mean
        )
      ) > 1
  ) {

    cor(
      all_df$obs,
      all_df$mc_mean
    )^2

  } else {

    NA_real_
  }


  RMSE <- if (
    nrow(
      all_df
    ) > 0
  ) {

    sqrt(
      mean(
        (
          all_df$mc_mean -
            all_df$obs
        )^2
      )
    )

  } else {

    NA_real_
  }


  MAPE <- if (
    nrow(
      positive_df
    ) > 0
  ) {

    mean(
      abs(
        positive_df$mc_mean -
          positive_df$obs
      ) /
        positive_df$obs
    ) *
      100

  } else {

    NA_real_
  }


  covered <- positive_df$obs >=
    positive_df$pred_lo &
    positive_df$obs <=
    positive_df$pred_hi


  tibble(
    n_all =
      nrow(
        all_df
      ),

    n_positive =
      nrow(
        positive_df
      ),

    R2 =
      R2,

    RMSE =
      RMSE,

    MAPE =
      MAPE,

    coverage_n =
      sum(
        covered,
        na.rm = TRUE
      ),

    coverage_pct =
      if (
        length(
          covered
        ) > 0
      ) {

        100 *
          mean(
            covered,
            na.rm = TRUE
          )

      } else {

        NA_real_
      }
  )
}


EXTERNAL_METRICS_OVERALL <- calc_external_metrics(
  EXTERNAL_COMPARE
)


EXTERNAL_METRICS_BY_ID <- bind_rows(
  lapply(
    TEST_IDS,
    function(id_i) {

      bind_cols(
        tibble(
          ID =
            id_i
        ),

        calc_external_metrics(
          EXTERNAL_COMPARE %>%
            filter(
              ID ==
                id_i
            )
        )
      )
    }
  )
)


EXTERNAL_METRICS_BY_TISSUE <- bind_rows(
  lapply(
    sort(
      unique(
        EXTERNAL_COMPARE$Tissue
      )
    ),
    function(tissue_i) {

      bind_cols(
        tibble(
          Tissue =
            tissue_i
        ),

        calc_external_metrics(
          EXTERNAL_COMPARE %>%
            filter(
              Tissue ==
                tissue_i
            )
        )
      )
    }
  )
)


## ============================================================
## 15) OBSERVATION-LEVEL METRICS / COVERAGE
## ============================================================

EXTERNAL_OBSERVATION_LEVEL <- EXTERNAL_COMPARE %>%

  mutate(
    positive_observation =
      obs > 0,

    covered_95 =
      ifelse(
        positive_observation &
          is.finite(
            pred_lo
          ) &
          is.finite(
            pred_hi
          ),

        obs >=
          pred_lo &
          obs <=
          pred_hi,

        NA
      ),

    AE_mc_mean =
      abs(
        mc_mean -
          obs
      ),

    RE_mc_mean_pct =
      ifelse(
        obs > 0,

        100 *
          abs(
            mc_mean -
              obs
          ) /
          obs,

        NA_real_
      )
  )


## ============================================================
## 16) SAVE FINAL STANDARDIZED OUTPUTS
## ============================================================

write_csv(
  EXTERNAL_COMPARE,
  file.path(
    EXTERNAL_DIR,
    "external_observed_vs_predicted.csv"
  )
)


write_csv(
  EXTERNAL_OBSERVATION_LEVEL,
  file.path(
    EXTERNAL_DIR,
    "external_observation_level_metrics.csv"
  )
)


write_csv(
  EXTERNAL_METRICS_OVERALL,
  file.path(
    EXTERNAL_DIR,
    "external_metrics_overall.csv"
  )
)


write_csv(
  EXTERNAL_METRICS_BY_ID,
  file.path(
    EXTERNAL_DIR,
    "external_metrics_by_id.csv"
  )
)


write_csv(
  EXTERNAL_METRICS_BY_TISSUE,
  file.path(
    EXTERNAL_DIR,
    "external_metrics_by_tissue.csv"
  )
)


## ------------------------------------------------------------
## Table S10 source:
## Overall + external datasets only.
##
## Tissue-level metrics remain available in a separate file but
## are not included in the primary S10 source because n=3 per
## tissue gives unstable tissue-specific R2 values.
## ------------------------------------------------------------

TABLE_S10 <- bind_rows(

  EXTERNAL_METRICS_OVERALL %>%

    mutate(
      Group =
        "Overall",

      `Size (nm)` =
        "50–5000",

      .before =
        1
    ),


  EXTERNAL_METRICS_BY_ID %>%

    mutate(
      Group =
        paste0(
          "ID ",
          ID
        ),

      `Size (nm)` =
        case_when(
          ID == 6 ~
            "50",

          ID == 7 ~
            "500",

          ID == 8 ~
            "5000",

          TRUE ~
            NA_character_
        ),

      .before =
        1
    ) %>%

    select(
      -ID
    )
)


write_csv(
  TABLE_S10,
  file.path(
    EXTERNAL_DIR,
    "Table_S10_external_evaluation.csv"
  )
)


## ============================================================
## 17) SAVE OLD-COMPATIBLE SINGLE-POINT TABLE
##
## Useful for comparison with the manuscript-era Script 05.
## ============================================================

OLD_COMPAT_SINGLEPOINT <- EXTERNAL_COMPARE %>%

  transmute(
    ID =
      ID,

    Size_nm =
      Size_nm,

    time =
      time,

    Tissue =
      Tissue,

    Obs =
      obs,

    Pred =
      pred,

    PredMean =
      mc_mean,

    PredLo =
      pred_lo,

    PredHi =
      pred_hi,

    AE =
      abs(
        pred -
          obs
      ),

    RE_percent =
      ifelse(
        obs > 0,

        100 *
          abs(
            pred -
              obs
          ) /
          obs,

        NA_real_
      ),

    Ratio =
      ifelse(
        obs > 0,

        pred /
          obs,

        NA_real_
      )
  )


write_csv(
  OLD_COMPAT_SINGLEPOINT,
  file.path(
    EXTERNAL_DIR,
    "pbpk_external_test_singlepoint_metrics.csv"
  )
)


## ============================================================
## 18) SAVE FULL EXTERNAL ANALYSIS OBJECT
## ============================================================

saveRDS(
  list(
    analysis_version =
      "manuscript-compatible",

    settings =
      list(
        MASTER_SEED =
          MASTER_SEED,

        TEST_IDS =
          TEST_IDS,

        META_PRED_TYPE =
          META_PRED_TYPE,

        META_STATISTIC =
          META_STATISTIC,

        SIM_DT_H =
          SIM_DT_H,

        N_SIM_UNC =
          N_SIM_UNC
      ),

    test_sizes_tbl =
      TEST_SIZES_TBL,

    central_parameter_vectors =
      CENTRAL_PARAMETER_VECTORS,

    central_predictions =
      CENTRAL_SIM,

    prediction_intervals =
      PRED_INTERVALS,

    external_compare =
      EXTERNAL_COMPARE,

    observation_level =
      EXTERNAL_OBSERVATION_LEVEL,

    external_metrics_overall =
      EXTERNAL_METRICS_OVERALL,

    external_metrics_by_id =
      EXTERNAL_METRICS_BY_ID,

    external_metrics_by_tissue =
      EXTERNAL_METRICS_BY_TISSUE
  ),

  file.path(
    EXTERNAL_DIR,
    "external_evaluation_outputs.rds"
  )
)


## ============================================================
## 19) ANALYSIS MANIFEST
## ============================================================

EXTERNAL_MANIFEST <- tibble(
  Item = c(
    "analysis_version",
    "external_IDs",
    "external_sizes_nm",
    "master_seed",
    "MC_draws_per_dataset",
    "meta_prediction_type",
    "central_parameter_statistic",
    "reported_point_prediction",
    "point_metric_population",
    "coverage_population",
    "coverage_interval"
  ),

  Value = c(
    "manuscript-compatible",
    "6,7,8",
    "50,500,5000",
    as.character(
      MASTER_SEED
    ),
    as.character(
      N_SIM_UNC
    ),
    META_PRED_TYPE,
    META_STATISTIC,
    "Monte Carlo mean",
    "all 18 matched observations",
    "14 positive observations",
    "2.5th–97.5th percentiles of parameter-propagated simulations"
  )
)


write_csv(
  EXTERNAL_MANIFEST,
  file.path(
    EXTERNAL_DIR,
    "external_evaluation_manifest.csv"
  )
)


## ============================================================
## 20) REPRODUCIBILITY CHECK AGAINST VERIFIED RUN
## ============================================================

EXPECTED_R2 <- 0.7855971

EXPECTED_RMSE <- 3.022362

EXPECTED_COVERAGE_N <- 10

EXPECTED_POSITIVE_N <- 14


R2_DIFF <- abs(
  EXTERNAL_METRICS_OVERALL$R2 -
    EXPECTED_R2
)


RMSE_DIFF <- abs(
  EXTERNAL_METRICS_OVERALL$RMSE -
    EXPECTED_RMSE
)


if (
  R2_DIFF >
    0.01
) {

  warning(
    "Final R2 differs from the previously verified run by ",
    R2_DIFF,
    "."
  )
}


if (
  RMSE_DIFF >
    0.05
) {

  warning(
    "Final RMSE differs from the previously verified run by ",
    RMSE_DIFF,
    "."
  )
}


if (
  EXTERNAL_METRICS_OVERALL$n_positive !=
    EXPECTED_POSITIVE_N
) {

  warning(
    "Expected 14 positive observations but found ",
    EXTERNAL_METRICS_OVERALL$n_positive,
    "."
  )
}


if (
  EXTERNAL_METRICS_OVERALL$coverage_n !=
    EXPECTED_COVERAGE_N
) {

  warning(
    "Expected 10 covered positive observations but found ",
    EXTERNAL_METRICS_OVERALL$coverage_n,
    "."
  )
}


## ============================================================
## 21) CONSOLE SUMMARY
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "External evaluation: overall\n"
)

cat(
  "============================================================\n"
)


print(
  EXTERNAL_METRICS_OVERALL,
  n = Inf,
  width = Inf
)


cat(
  "\nExternal evaluation by ID:\n"
)


print(
  EXTERNAL_METRICS_BY_ID,
  n = Inf,
  width = Inf
)


cat(
  "\n95% coverage by ID:\n"
)


print(
  EXTERNAL_OBSERVATION_LEVEL %>%

    filter(
      positive_observation
    ) %>%

    group_by(
      ID
    ) %>%

    summarise(
      n_positive =
        n(),

      coverage_n =
        sum(
          covered_95,
          na.rm = TRUE
        ),

      coverage_pct =
        100 *
        mean(
          covered_95,
          na.rm = TRUE
        ),

      .groups =
        "drop"
    ),
  n = Inf,
  width = Inf
)


cat(
  "\n============================================================\n"
)

cat(
  "FINAL MANUSCRIPT-COMPATIBLE EXTERNAL EVALUATION\n"
)

cat(
  "============================================================\n"
)


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
  "N all: ",
  EXTERNAL_METRICS_OVERALL$n_all,
  "\n",
  sep = ""
)


cat(
  "N positive: ",
  EXTERNAL_METRICS_OVERALL$n_positive,
  "\n",
  sep = ""
)


cat(
  "R2 (Observed vs MC mean; all observations): ",
  sprintf(
    "%.7f",
    EXTERNAL_METRICS_OVERALL$R2
  ),
  "\n",
  sep = ""
)


cat(
  "RMSE (Observed vs MC mean; all observations): ",
  sprintf(
    "%.6f",
    EXTERNAL_METRICS_OVERALL$RMSE
  ),
  " %ID\n",
  sep = ""
)


cat(
  "95% parameter-propagated interval coverage: ",
  EXTERNAL_METRICS_OVERALL$coverage_n,
  "/",
  EXTERNAL_METRICS_OVERALL$n_positive,
  " = ",
  sprintf(
    "%.1f",
    EXTERNAL_METRICS_OVERALL$coverage_pct
  ),
  "%\n",
  sep = ""
)


cat(
  "Monte Carlo simulations per dataset: ",
  N_SIM_UNC,
  "\n",
  sep = ""
)


cat(
  "Master seed: ",
  MASTER_SEED,
  "\n",
  sep = ""
)


cat(
  "\nOutputs saved to:\n  ",
  normalizePath(
    EXTERNAL_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat(
  "============================================================\n"
)