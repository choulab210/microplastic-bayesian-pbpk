## ============================================================
## 04_run_meta_regression.R
##
## FINAL MANUSCRIPT-COMPATIBLE RANDOM-EFFECTS META-REGRESSION
##
## Reads:
##   results/mcmc/MCMC3_results_fullset_perID_priors_mixed.rds
##
## Does NOT rerun MCMC.
##
## Model for parameter k:
##
##   log(theta_i,k) =
##       beta0_k + beta1_k * log(Size_i) + u_i,k
##
##   u_i,k ~ N(0, tau_k^2)
##
## Within-dataset variance:
##   vi = SD[log(theta_i,k)]^2
##
## Estimation:
##   metafor::rma.uni(..., method = "REML")
##
## Development datasets:
##   ID 1 =   20 nm
##   ID 2 =  220 nm
##   ID 3 = 1000 nm
##   ID 4 = 6000 nm
##   ID 5 =  220 nm
##
## IDs 2 and 5 remain independent meta-regression observations.
## ============================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(metafor)
})

options(stringsAsFactors = FALSE)


## ============================================================
## 0) DIRECTORIES
## ============================================================

MCMC_DIR <- file.path(
  "results",
  "mcmc"
)

META_DIR <- file.path(
  "results",
  "meta"
)

dir.create(
  META_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 1) ANALYSIS SETTINGS
## ============================================================

DEV_IDS <- 1:5

SIZE_PARAMS <- c(
  "Ktrans",
  "Kgw2b",
  "Kshed",
  "KLRESmax",
  "Kfeces",
  "KurineC"
)


## Remaining sampled parameters pooled across development IDs
## exactly as in the manuscript-era workflow.

GLOBAL_PARAMS_INTENDED <- c(
  "KLRESrelease",
  "KSRESmax",
  "KSRESrelease",
  "KKRESmax",
  "KKRESrelease",
  "KLURESmax",
  "KLURESrelease",
  "KRRESmax",
  "KRRESrelease",
  "PAL",
  "PAK",
  "PAS",
  "PALU",
  "PAGI",
  "PAR",
  "KbileC"
)


ID_SIZE_MAP <- c(
  "1" = 20,
  "2" = 220,
  "3" = 1000,
  "4" = 6000,
  "5" = 220
)


SIZE_TABLE <- tibble(
  ID = as.character(
    DEV_IDS
  ),

  Size_nm = as.numeric(
    ID_SIZE_MAP[
      as.character(DEV_IDS)
    ]
  )
) %>%

  mutate(
    logSize = log(
      Size_nm
    )
  )


## ============================================================
## 2) LOAD CANONICAL MANUSCRIPT MCMC
## ============================================================

MCMC_FILE <- file.path(
  MCMC_DIR,
  "MCMC3_results_fullset_perID_priors_mixed.rds"
)


if (!file.exists(MCMC_FILE)) {

  stop(
    "\nMissing manuscript-compatible MCMC file:\n",
    MCMC_FILE,
    "\n"
  )
}


cat(
  "\n============================================================\n"
)

cat(
  "MANUSCRIPT-COMPATIBLE RANDOM-EFFECTS META-REGRESSION\n"
)

cat(
  "============================================================\n"
)

cat(
  "MCMC source: ",
  MCMC_FILE,
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

cat(
  "Meta-regression parameters: ",
  paste(
    SIZE_PARAMS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


MCMC3 <- readRDS(
  MCMC_FILE
)


missing_ids <- DEV_IDS[
  !as.character(DEV_IDS) %in%
    names(MCMC3)
]


if (
  length(missing_ids) > 0
) {

  stop(
    "Canonical MCMC object missing IDs: ",
    paste(
      missing_ids,
      collapse = ", "
    )
  )
}


cat(
  "\nDevelopment dataset sizes:\n"
)

print(
  SIZE_TABLE,
  n = Inf
)


write_csv(
  SIZE_TABLE,
  file.path(
    META_DIR,
    "development_id_size_table.csv"
  )
)


## ============================================================
## 3) HELPER: COMBINE THREE CHAINS FOR ONE ID
## ============================================================

combine_id_draws <- function(id) {

  id_chr <- as.character(
    id
  )


  obj <- MCMC3[[id_chr]]


  if (
    is.null(obj)
  ) {

    stop(
      "MCMC object missing ID ",
      id
    )
  }


  if (
    is.null(obj$chains) ||
      length(obj$chains) == 0
  ) {

    stop(
      "No MCMC chains found for ID ",
      id
    )
  }


  draws <- do.call(
    rbind,
    lapply(
      obj$chains,
      function(m) {

        as.data.frame(
          m$pars
        )
      }
    )
  )


  if (
    is.null(draws) ||
      nrow(draws) == 0
  ) {

    stop(
      "No posterior draws recovered for ID ",
      id
    )
  }


  draws
}


FIRST_DRAWS <- combine_id_draws(
  DEV_IDS[1]
)


MCMC_PARAMETER_NAMES <- colnames(
  FIRST_DRAWS
)


missing_params <- setdiff(
  SIZE_PARAMS,
  MCMC_PARAMETER_NAMES
)


if (
  length(missing_params) > 0
) {

  stop(
    "Meta-regression parameters missing from MCMC: ",
    paste(
      missing_params,
      collapse = ", "
    )
  )
}


GLOBAL_PARAMS <- intersect(
  GLOBAL_PARAMS_INTENDED,
  MCMC_PARAMETER_NAMES
)


## ============================================================
## 4) EXTRACT PER-ID POSTERIOR LOG MOMENTS
## ============================================================

LOG_MOM_ROWS <- list()


for (
  parameter in SIZE_PARAMS
) {

  for (
    id in DEV_IDS
  ) {

    D <- combine_id_draws(
      id
    )


    x <- as.numeric(
      D[[parameter]]
    )


    LOG_MOM_ROWS[[length(LOG_MOM_ROWS) + 1]] <- tibble(
      ID = as.character(
        id
      ),

      Parameter = parameter,

      mean_log = mean(
        x,
        na.rm = TRUE
      ),

      sd_log = sd(
        x,
        na.rm = TRUE
      )
    )
  }
}


LOG_MOM <- bind_rows(
  LOG_MOM_ROWS
) %>%

  left_join(
    SIZE_TABLE,
    by = "ID"
  ) %>%

  mutate(
    vi = pmax(
      sd_log^2,
      1e-10
    )
  ) %>%

  arrange(
    match(
      Parameter,
      SIZE_PARAMS
    ),
    as.numeric(ID)
  )


expected_rows <- length(
  SIZE_PARAMS
) *
  length(
    DEV_IDS
  )


if (
  nrow(LOG_MOM) !=
    expected_rows
) {

  stop(
    "Expected ",
    expected_rows,
    " meta-regression rows but found ",
    nrow(LOG_MOM)
  )
}


write_csv(
  LOG_MOM,
  file.path(
    META_DIR,
    "meta_log_moments_by_id.csv"
  )
)


cat(
  "\nPosterior log moments extracted successfully: ",
  nrow(LOG_MOM),
  " rows\n",
  sep = ""
)


## ============================================================
## 5) FIT ONE RANDOM-EFFECTS META-REGRESSION PER PARAMETER
## ============================================================

META_FIT_OBJECTS <- list()

META_SUMMARY_ROWS <- list()


for (
  parameter in SIZE_PARAMS
) {

  df <- LOG_MOM %>%

    filter(
      Parameter ==
        parameter
    )


  if (
    nrow(df) !=
      length(DEV_IDS)
  ) {

    stop(
      parameter,
      " has ",
      nrow(df),
      " datasets; expected ",
      length(DEV_IDS)
    )
  }


  fit <- metafor::rma.uni(
    yi =
      mean_log,

    vi =
      vi,

    mods =
      ~ logSize,

    method =
      "REML",

    data =
      df,

    control =
      list(
        stepadj = 0.5,
        maxiter = 200
      )
  )


  META_FIT_OBJECTS[[parameter]] <- fit


  beta <- as.numeric(
    coef(fit)
  )


  V <- as.matrix(
    vcov(fit)
  )


  se_beta1 <- sqrt(
    V[2, 2]
  )


  p_beta1 <- 2 *
    pnorm(
      -abs(
        beta[2] /
          se_beta1
      )
    )


  META_SUMMARY_ROWS[[length(META_SUMMARY_ROWS) + 1]] <- tibble(
    Parameter =
      parameter,

    N =
      nrow(df),

    beta0 =
      beta[1],

    beta1 =
      beta[2],

    se_beta1 =
      se_beta1,

    beta1_ci_low =
      beta[2] -
      1.96 *
      se_beta1,

    beta1_ci_high =
      beta[2] +
      1.96 *
      se_beta1,

    p_beta1 =
      p_beta1,

    tau2 =
      as.numeric(
        fit$tau2
      ),

    vc11 =
      V[1, 1],

    vc22 =
      V[2, 2],

    vc12 =
      V[1, 2]
  )
}


META_SUMMARY <- bind_rows(
  META_SUMMARY_ROWS
) %>%

  mutate(
    Parameter = factor(
      Parameter,
      levels =
        SIZE_PARAMS
    )
  ) %>%

  arrange(
    Parameter
  ) %>%

  mutate(
    Parameter =
      as.character(
        Parameter
      )
  )


cat(
  "\n============================================================\n"
)

cat(
  "Random-effects meta-regression summary\n"
)

cat(
  "============================================================\n"
)


print(
  META_SUMMARY %>%

    select(
      Parameter,
      N,
      beta1,
      se_beta1,
      beta1_ci_low,
      beta1_ci_high,
      p_beta1,
      tau2
    ),
  n = Inf,
  width = Inf
)


## ============================================================
## 6) OLD-SCRIPT-COMPATIBLE META_FITS OBJECT
##
## Required fields:
##   param
##   engine
##   n
##   tau2
##   coef0
##   coef1
##   vc11
##   vc22
##   vc12
##   obj
## ============================================================

META_OLD_ROWS <- list()


for (
  parameter in SIZE_PARAMS
) {

  row <- META_SUMMARY %>%

    filter(
      Parameter ==
        parameter
    )


  META_OLD_ROWS[[length(META_OLD_ROWS) + 1]] <- tibble(
    param =
      parameter,

    engine =
      "metafor",

    n =
      row$N,

    tau2 =
      row$tau2,

    coef0 =
      row$beta0,

    coef1 =
      row$beta1,

    vc11 =
      row$vc11,

    vc22 =
      row$vc22,

    vc12 =
      row$vc12,

    obj =
      list(
        META_FIT_OBJECTS[[parameter]]
      )
  )
}


META_FITS <- bind_rows(
  META_OLD_ROWS
)


saveRDS(
  META_FITS,
  file.path(
    META_DIR,
    "META_FITS_metafor_REML.rds"
  )
)


write_csv(
  META_FITS %>%
    select(
      -obj
    ),
  file.path(
    META_DIR,
    "META_FITS_summary.csv"
  )
)


## ============================================================
## 7) TABLE S9 SOURCE
## ============================================================

TABLE_S9 <- META_SUMMARY %>%

  transmute(
    Parameter,

    N,

    `Size slope, beta1` =
      beta1,

    SE =
      se_beta1,

    `95% CI lower` =
      beta1_ci_low,

    `95% CI upper` =
      beta1_ci_high,

    `p-value` =
      p_beta1,

    tau2 =
      tau2
  )


write_csv(
  TABLE_S9,
  file.path(
    META_DIR,
    "Table_S9_meta_regression.csv"
  )
)


## ============================================================
## 8) OLD-COMPATIBLE SLOPE TEST TABLE
## ============================================================

SLOPE_TESTS <- META_SUMMARY %>%

  transmute(
    param =
      Parameter,

    slope =
      beta1,

    se =
      se_beta1,

    lo95 =
      beta1_ci_low,

    hi95 =
      beta1_ci_high,

    pval =
      p_beta1
  ) %>%

  arrange(
    pval
  )


write_csv(
  SLOPE_TESTS,
  file.path(
    META_DIR,
    "META_slope_tests.csv"
  )
)


## ============================================================
## 9) META-PREDICTION HELPERS
## ============================================================

get_meta_summary_row <- function(
    parameter
) {

  row <- META_SUMMARY %>%

    filter(
      Parameter ==
        parameter
    )


  if (
    nrow(row) != 1
  ) {

    stop(
      "Meta-regression summary missing for ",
      parameter
    )
  }


  row
}


predict_param_by_size <- function(
    parameter,
    sizes_nm,
    type = c(
      "population",
      "new_id"
    )
) {

  type <- match.arg(
    type
  )


  row <- get_meta_summary_row(
    parameter
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


  output_rows <- list()


  for (
    S in sizes_nm
  ) {

    x <- c(
      1,
      log(S)
    )


    mu <- as.numeric(
      row$beta0 +
        row$beta1 *
        log(S)
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

      var_total <- var_total +
        pmax(
          row$tau2,
          0
        )
    }


    output_rows[[length(output_rows) + 1]] <- tibble(
      Parameter =
        parameter,

      Size_nm =
        S,

      type =
        type,

      mu_log =
        mu,

      var_mu_log =
        var_mu,

      tau2 =
        row$tau2,

      var_total_log =
        var_total,

      ## Manuscript-era natural-scale mean
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


  bind_rows(
    output_rows
  )
}


predict_params_by_size <- function(
    parameters,
    sizes_nm,
    type = c(
      "population",
      "new_id"
    )
) {

  type <- match.arg(
    type
  )


  out <- list()


  for (
    parameter in parameters
  ) {

    out[[length(out) + 1]] <- predict_param_by_size(
      parameter =
        parameter,

      sizes_nm =
        sizes_nm,

      type =
        type
    )
  }


  bind_rows(
    out
  )
}


## ============================================================
## 10) META CURVES FOR FIGURE 5
## ============================================================

SIZE_GRID <- exp(
  seq(
    log(20),
    log(6000),
    length.out = 200
  )
)


CURVE_ROWS <- list()


for (
  parameter in SIZE_PARAMS
) {

  row <- get_meta_summary_row(
    parameter
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


  for (
    S in SIZE_GRID
  ) {

    x <- c(
      1,
      log(S)
    )


    mu <- as.numeric(
      row$beta0 +
        row$beta1 *
        log(S)
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


    var_new <- var_mu +
      pmax(
        row$tau2,
        0
      )


    CURVE_ROWS[[length(CURVE_ROWS) + 1]] <- tibble(
      Parameter =
        parameter,

      Size_nm =
        S,

      mu_log =
        mu,

      fitted_mean =
        exp(
          mu +
            0.5 *
            var_mu
        ),

      fitted_median =
        exp(
          mu
        ),

      CI_lower =
        exp(
          mu -
            1.96 *
            sqrt(
              var_mu
            )
        ),

      CI_upper =
        exp(
          mu +
            1.96 *
            sqrt(
              var_mu
            )
        ),

      PI_lower =
        exp(
          mu -
            1.96 *
            sqrt(
              var_new
            )
        ),

      PI_upper =
        exp(
          mu +
            1.96 *
            sqrt(
              var_new
            )
        )
    )
  }
}


META_CURVES <- bind_rows(
  CURVE_ROWS
)


write_csv(
  META_CURVES,
  file.path(
    META_DIR,
    "meta_regression_curves.csv"
  )
)


## ============================================================
## 11) POOL NON-META PARAMETERS ACROSS ALL DEVELOPMENT IDs
## ============================================================

GLOBAL_DRAW_ROWS <- list()


for (
  id in DEV_IDS
) {

  D <- combine_id_draws(
    id
  )


  keep <- intersect(
    GLOBAL_PARAMS,
    colnames(D)
  )


  if (
    length(keep) == 0
  ) {
    next
  }


  DD <- D[
    ,
    keep,
    drop = FALSE
  ]


  DD$ID <- id


  GLOBAL_DRAW_ROWS[[length(GLOBAL_DRAW_ROWS) + 1]] <- DD
}


GLOBAL_DRAWS <- bind_rows(
  GLOBAL_DRAW_ROWS
)


GLOBAL_SUM_ROWS <- list()


for (
  parameter in GLOBAL_PARAMS
) {

  if (
    !parameter %in%
      colnames(
        GLOBAL_DRAWS
      )
  ) {
    next
  }


  x_log <- as.numeric(
    GLOBAL_DRAWS[[parameter]]
  )


  x_nat <- exp(
    x_log
  )


  GLOBAL_SUM_ROWS[[length(GLOBAL_SUM_ROWS) + 1]] <- tibble(
    Parameter =
      parameter,

    mean_log =
      mean(
        x_log,
        na.rm = TRUE
      ),

    sd_log =
      sd(
        x_log,
        na.rm = TRUE
      ),

    mean =
      mean(
        x_nat,
        na.rm = TRUE
      ),

    median =
      median(
        x_nat,
        na.rm = TRUE
      ),

    q025 =
      as.numeric(
        quantile(
          x_nat,
          0.025,
          na.rm = TRUE
        )
      ),

    q975 =
      as.numeric(
        quantile(
          x_nat,
          0.975,
          na.rm = TRUE
        )
      )
  )
}


GLOBAL_SUM <- bind_rows(
  GLOBAL_SUM_ROWS
)


write_csv(
  GLOBAL_SUM,
  file.path(
    META_DIR,
    "meta_global_params_metafor.csv"
  )
)


cat(
  "\nGlobal pooled parameter summaries generated: ",
  nrow(GLOBAL_SUM),
  "\n",
  sep = ""
)


## ============================================================
## 12) SAVE EXAMPLE POPULATION / NEW-ID PREDICTIONS
## ============================================================

EXAMPLE_SIZES <- c(
  20,
  50,
  100,
  200,
  500,
  1000,
  5000,
  6000
)


PRED_POP <- predict_params_by_size(
  parameters =
    SIZE_PARAMS,

  sizes_nm =
    EXAMPLE_SIZES,

  type =
    "population"
)


PRED_NEW <- predict_params_by_size(
  parameters =
    SIZE_PARAMS,

  sizes_nm =
    EXAMPLE_SIZES,

  type =
    "new_id"
)


write_csv(
  PRED_POP,
  file.path(
    META_DIR,
    "meta_pred_population_metafor.csv"
  )
)


write_csv(
  PRED_NEW,
  file.path(
    META_DIR,
    "meta_pred_newid_metafor.csv"
  )
)


## ============================================================
## 13) CENTRAL EXTERNAL PARAMETER VECTORS
## ============================================================

build_param_vector_for_size <- function(
    size_nm,
    type = c(
      "population",
      "new_id"
    ),
    statistic = c(
      "mean",
      "median",
      "q025",
      "q975"
    )
) {

  type <- match.arg(
    type
  )

  statistic <- match.arg(
    statistic
  )


  SZ <- predict_params_by_size(
    parameters =
      SIZE_PARAMS,

    sizes_nm =
      size_nm,

    type =
      type
  )


  sp <- setNames(
    SZ[[statistic]],
    SZ$Parameter
  )


  if (
    nrow(GLOBAL_SUM) > 0
  ) {

    gp <- setNames(
      GLOBAL_SUM[[statistic]],
      GLOBAL_SUM$Parameter
    )


    return(
      c(
        sp,
        gp
      )
    )
  }


  sp
}


EXTERNAL_SIZES <- c(
  50,
  500,
  5000
)


EXT_VECTOR_ROWS <- list()


for (
  S in EXTERNAL_SIZES
) {

  v <- build_param_vector_for_size(
    size_nm =
      S,

    type =
      "new_id",

    statistic =
      "mean"
  )


  EXT_VECTOR_ROWS[[length(EXT_VECTOR_ROWS) + 1]] <- tibble(
    Size_nm =
      S,

    Parameter =
      names(v),

    Value =
      as.numeric(v)
  )
}


EXTERNAL_CENTRAL_VECTORS <- bind_rows(
  EXT_VECTOR_ROWS
)


write_csv(
  EXTERNAL_CENTRAL_VECTORS,
  file.path(
    META_DIR,
    "pbpk_external_central_parameter_vectors.csv"
  )
)


## ============================================================
## 14) SAVE AUTHORITATIVE META OBJECT
## ============================================================

saveRDS(
  list(
    analysis_version =
      "manuscript-compatible",

    SIZE_PARAMS =
      SIZE_PARAMS,

    GLOBAL_PARAMS =
      GLOBAL_PARAMS,

    sizes_tbl =
      SIZE_TABLE,

    LOG_MOM =
      LOG_MOM,

    ## Old Script 05-compatible object
    META_FITS =
      META_FITS,

    GLOBAL_SUM =
      GLOBAL_SUM,

    ## Standardized downstream objects
    META_SUMMARY =
      META_SUMMARY,

    META_CURVES =
      META_CURVES,

    pred_pop =
      PRED_POP,

    pred_new =
      PRED_NEW,

    external_central_vectors =
      EXTERNAL_CENTRAL_VECTORS
  ),

  file.path(
    META_DIR,
    "meta_analysis_outputs.rds"
  )
)


## ============================================================
## 15) META ANALYSIS MANIFEST
## ============================================================

META_MANIFEST <- tibble(
  Item = c(
    "analysis_version",
    "development_IDs",
    "development_sizes_nm",
    "meta_parameters",
    "model",
    "estimation",
    "within_dataset_variance",
    "duplicate_size_handling"
  ),

  Value = c(
    "manuscript-compatible",
    "1,2,3,4,5",
    "20,220,1000,6000,220",
    paste(
      SIZE_PARAMS,
      collapse = ","
    ),
    "log(theta) ~ log(Size_nm)",
    "metafor REML",
    "posterior SD(log(theta))^2",
    "IDs 2 and 5 retained as independent datasets"
  )
)


write_csv(
  META_MANIFEST,
  file.path(
    META_DIR,
    "meta_analysis_manifest.csv"
  )
)


## ============================================================
## 16) COMPLETION
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "META-REGRESSION COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nFinal Table S9 source:\n"
)


print(
  TABLE_S9,
  n = Inf,
  width = Inf
)


cat(
  "\nBetween-dataset heterogeneity:\n"
)


print(
  META_SUMMARY %>%
    select(
      Parameter,
      tau2
    ),
  n = Inf
)


cat(
  "\nKey outputs:\n"
)

cat(
  "  results/meta/meta_analysis_outputs.rds\n"
)

cat(
  "  results/meta/Table_S9_meta_regression.csv\n"
)

cat(
  "  results/meta/META_FITS_metafor_REML.rds\n"
)

cat(
  "  results/meta/meta_global_params_metafor.csv\n"
)

cat(
  "  results/meta/meta_regression_curves.csv\n"
)

cat(
  "  results/meta/pbpk_external_central_parameter_vectors.csv\n"
)


cat(
  "\nNo MCMC was rerun.\n"
)

cat(
  "============================================================\n"
)