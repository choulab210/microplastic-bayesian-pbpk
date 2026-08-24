## ============================================================
## 04_run_meta_regression.R
##
## Random-effects meta-regression for particle-size-informed
## PBPK parameter prediction.
##
## Inputs:
##   results/deterministic/deterministic_outputs.rds
##   results/sensitivity/sensitivity_outputs.rds
##   results/mcmc/ID1/mcmc_results.rds
##   ...
##   results/mcmc/ID5/mcmc_results.rds
##
## Outputs:
##   results/meta/
##
## Model for each prioritized parameter:
##
##   meanlog_i =
##     beta0 + beta1 * log(Size_i) + u_i + epsilon_i
##
##   u_i ~ N(0, tau^2)
##
## Within-dataset variance:
##
##   vi = posterior SD(log parameter)^2
##
## Fitting:
##   metafor::rma.uni(..., method = "REML")
##
## Important:
##   - IDs 1–5 remain separate dataset-level inputs
##   - IDs 2 and 5 are both 220 nm but are NOT averaged
##   - particle size is natural-log transformed
## ============================================================


## ============================================================
## 0) Packages
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(metafor)
})

options(stringsAsFactors = FALSE)


## ============================================================
## 1) Paths and settings
## ============================================================

DETERMINISTIC_FILE <- file.path(
  "results",
  "deterministic",
  "deterministic_outputs.rds"
)

SENSITIVITY_FILE <- file.path(
  "results",
  "sensitivity",
  "sensitivity_outputs.rds"
)

MCMC_DIR <- file.path(
  "results",
  "mcmc"
)

META_DIR <- file.path(
  "results",
  "meta"
)

FIGURE_DIR <- file.path(
  "figures",
  "meta"
)

dir.create(
  META_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ------------------------------------------------------------
## Canonical parameter order for manuscript Figure 5 / Table S9
## ------------------------------------------------------------

META_PARAM_ORDER <- c(
  "Ktrans",
  "Kgw2b",
  "Kshed",
  "KLRESmax",
  "Kfeces",
  "KurineC"
)


## ------------------------------------------------------------
## Size grid for meta-regression prediction curves
## ------------------------------------------------------------

SIZE_GRID <- exp(
  seq(
    log(20),
    log(6000),
    length.out = 200
  )
)


## ------------------------------------------------------------
## Fallback ID -> particle size mapping
## ------------------------------------------------------------

ID_SIZE_MAP <- c(
  "1" = 20,
  "2" = 220,
  "3" = 1000,
  "4" = 6000,
  "5" = 220
)


## ============================================================
## 2) Check required files
## ============================================================

if (!file.exists(DETERMINISTIC_FILE)) {
  stop(
    "Missing deterministic output: ",
    DETERMINISTIC_FILE
  )
}

if (!file.exists(SENSITIVITY_FILE)) {
  stop(
    "Missing sensitivity output: ",
    SENSITIVITY_FILE
  )
}


## ============================================================
## 3) Load deterministic and sensitivity outputs
## ============================================================

det <- readRDS(
  DETERMINISTIC_FILE
)

sens <- readRDS(
  SENSITIVITY_FILE
)


if (is.null(det$ids)) {
  stop("ids not found in deterministic_outputs.rds.")
}

if (is.null(det$obs_df)) {
  stop("obs_df not found in deterministic_outputs.rds.")
}

if (is.null(sens$prioritized_params)) {
  stop("prioritized_params not found in sensitivity_outputs.rds.")
}


ids <- as.numeric(
  det$ids
)

obs_df <- det$obs_df

PRIORITIZED_PARAMS <- sens$prioritized_params


## ------------------------------------------------------------
## Verify prioritized parameters
## ------------------------------------------------------------

if (!setequal(
  PRIORITIZED_PARAMS,
  META_PARAM_ORDER
)) {

  stop(
    "Sensitivity-prioritized parameters do not match ",
    "the intended six meta-regression parameters.\n",
    "Sensitivity output: ",
    paste(
      PRIORITIZED_PARAMS,
      collapse = ", "
    ),
    "\nExpected: ",
    paste(
      META_PARAM_ORDER,
      collapse = ", "
    )
  )
}


META_PARAMS <- META_PARAM_ORDER


cat("\n============================================================\n")
cat("Random-effects PBPK meta-regression\n")
cat("============================================================\n")

cat(
  "Development IDs: ",
  paste(ids, collapse = ", "),
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
## 4) Build development ID-size table
## ============================================================

build_size_table <- function(
    obs_df,
    ids,
    fallback_map = ID_SIZE_MAP
) {

  if ("Size_nm" %in% names(obs_df)) {

    size_tbl <- obs_df %>%
      mutate(
        ID = as.character(ID),
        Size_nm = as.numeric(Size_nm)
      ) %>%
      filter(
        ID %in% as.character(ids)
      ) %>%
      group_by(ID) %>%
      summarise(
        Size_nm = first(
          Size_nm[
            is.finite(Size_nm)
          ]
        ),
        .groups = "drop"
      )

  } else {

    size_tbl <- tibble(
      ID = as.character(ids),
      Size_nm = as.numeric(
        fallback_map[
          as.character(ids)
        ]
      )
    )
  }


  size_tbl <- size_tbl %>%
    mutate(
      Size_nm = ifelse(
        is.finite(Size_nm),
        Size_nm,
        as.numeric(
          fallback_map[ID]
        )
      ),
      logSize = log(Size_nm)
    )


  if (
    any(
      !is.finite(
        size_tbl$Size_nm
      )
    )
  ) {

    bad_ids <- size_tbl$ID[
      !is.finite(
        size_tbl$Size_nm
      )
    ]

    stop(
      "Unable to determine particle size for IDs: ",
      paste(
        bad_ids,
        collapse = ", "
      )
    )
  }


  size_tbl
}


size_tbl <- build_size_table(
  obs_df = obs_df,
  ids = ids
)


cat("\nDevelopment dataset sizes:\n")
print(
  size_tbl
)


## ============================================================
## 5) Load final MCMC outputs for IDs 1–5
## ============================================================

MCMC_RESULTS <- setNames(
  vector(
    "list",
    length(ids)
  ),
  as.character(ids)
)


for (id in ids) {

  mcmc_file <- file.path(
    MCMC_DIR,
    paste0(
      "ID",
      id
    ),
    "mcmc_results.rds"
  )


  if (!file.exists(mcmc_file)) {

    stop(
      "Missing final MCMC output for ID ",
      id,
      ": ",
      mcmc_file
    )
  }


  MCMC_RESULTS[[as.character(id)]] <- readRDS(
    mcmc_file
  )
}


cat(
  "\nLoaded final MCMC outputs for IDs: ",
  paste(
    ids,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


## ============================================================
## 6) Extract posterior log-scale moments
##
## yi = posterior mean(log parameter)
## vi = posterior variance(log parameter)
## ============================================================

extract_log_moments <- function(
    id,
    mcmc_result,
    parameters
) {

  posterior <- mcmc_result$posterior_summary


  if (is.null(posterior)) {
    stop(
      "posterior_summary missing for ID ",
      id
    )
  }


  required_cols <- c(
    "Parameter",
    "mean_log",
    "sd_log"
  )


  missing_cols <- setdiff(
    required_cols,
    names(posterior)
  )


  if (length(missing_cols) > 0) {

    stop(
      "Posterior summary for ID ",
      id,
      " is missing columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }


  out <- posterior %>%
    filter(
      Parameter %in% parameters
    ) %>%
    transmute(
      ID = as.character(id),
      Parameter = as.character(Parameter),
      meanlog = as.numeric(mean_log),
      sdlog = as.numeric(sd_log),
      vi = pmax(
        as.numeric(sd_log)^2,
        1e-10
      )
    )


  missing_parameters <- setdiff(
    parameters,
    out$Parameter
  )


  if (length(missing_parameters) > 0) {

    stop(
      "ID ",
      id,
      " is missing meta-regression parameters: ",
      paste(
        missing_parameters,
        collapse = ", "
      )
    )
  }


  out
}


LOG_MOMENTS <- bind_rows(
  lapply(
    ids,
    function(id) {

      extract_log_moments(
        id = id,
        mcmc_result = MCMC_RESULTS[[as.character(id)]],
        parameters = META_PARAMS
      )
    }
  )
) %>%
  left_join(
    size_tbl,
    by = "ID"
  ) %>%
  mutate(
    Parameter = factor(
      Parameter,
      levels = META_PARAM_ORDER
    )
  ) %>%
  arrange(
    Parameter,
    as.numeric(ID)
  )


write_csv(
  LOG_MOMENTS,
  file.path(
    META_DIR,
    "meta_log_moments_by_id.csv"
  )
)


cat("\nPosterior moments used in meta-regression:\n")
print(
  LOG_MOMENTS
)


## ============================================================
## 7) Fit random-effects meta-regression
## ============================================================

fit_meta_parameter <- function(
    parameter,
    log_moments
) {

  df <- log_moments %>%
    filter(
      as.character(Parameter) ==
        parameter
    ) %>%
    mutate(
      vi = pmax(
        vi,
        1e-10
      )
    )


  if (nrow(df) != length(ids)) {

    stop(
      "Expected ",
      length(ids),
      " datasets for ",
      parameter,
      " but found ",
      nrow(df),
      "."
    )
  }


  fit <- metafor::rma.uni(
    yi = meanlog,
    vi = vi,
    mods = ~ logSize,
    method = "REML",
    data = df,
    control = list(
      stepadj = 0.5,
      maxiter = 500
    )
  )


  beta <- as.numeric(
    coef(
      fit
    )
  )


  Vbeta <- as.matrix(
    vcov(
      fit
    )
  )


  list(
    parameter = parameter,
    data = df,
    fit = fit,
    beta0 = beta[1],
    beta1 = beta[2],
    vcov = Vbeta,
    tau2 = as.numeric(
      fit$tau2
    ),
    n = nrow(df)
  )
}


META_FIT_LIST <- setNames(
  lapply(
    META_PARAMS,
    function(parameter) {

      fit_meta_parameter(
        parameter = parameter,
        log_moments = LOG_MOMENTS
      )
    }
  ),
  META_PARAMS
)


## ============================================================
## 8) Create canonical meta-regression summary
##
## Everything is generated by parameter name.
## ============================================================

META_SUMMARY <- bind_rows(
  lapply(
    META_PARAMS,
    function(parameter) {

      x <- META_FIT_LIST[[parameter]]


      se_beta1 <- sqrt(
        x$vcov[
          2,
          2
        ]
      )


      z_beta1 <- x$beta1 /
        se_beta1


      tibble(
        Parameter = parameter,

        N = x$n,

        beta0 = x$beta0,

        beta1 = x$beta1,

        SE_beta1 = se_beta1,

        CI_lower = x$beta1 -
          1.96 *
          se_beta1,

        CI_upper = x$beta1 +
          1.96 *
          se_beta1,

        p_value = 2 *
          pnorm(
            -abs(
              z_beta1
            )
          ),

        tau2 = x$tau2,

        vc11 = x$vcov[
          1,
          1
        ],

        vc22 = x$vcov[
          2,
          2
        ],

        vc12 = x$vcov[
          1,
          2
        ]
      )
    }
  )
) %>%
  mutate(
    Parameter = factor(
      Parameter,
      levels = META_PARAM_ORDER
    )
  ) %>%
  arrange(
    Parameter
  ) %>%
  mutate(
    Parameter = as.character(
      Parameter
    )
  )


cat("\n============================================================\n")
cat("Random-effects meta-regression results\n")
cat("============================================================\n")

print(
  META_SUMMARY
)


write_csv(
  META_SUMMARY,
  file.path(
    META_DIR,
    "meta_regression_summary.csv"
  )
)


## ============================================================
## 9) Generate Table S9 source data
## ============================================================

TABLE_S9 <- META_SUMMARY %>%
  transmute(
    Parameter,
    N,
    `Size slope, beta1` = beta1,
    SE = SE_beta1,
    `95% CI lower` = CI_lower,
    `95% CI upper` = CI_upper,
    `p-value` = p_value,
    tau2 = tau2
  )


write_csv(
  TABLE_S9,
  file.path(
    META_DIR,
    "Table_S9_meta_regression.csv"
  )
)


cat("\nTable S9 source data:\n")
print(
  TABLE_S9
)


## ============================================================
## 10) Meta-regression prediction helper
##
## type = population
##   coefficient uncertainty only
##
## type = new_id
##   coefficient uncertainty + tau2
## ============================================================

predict_meta_parameter <- function(
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


  meta_row <- META_SUMMARY %>%
    filter(
      Parameter ==
        parameter
    )


  if (nrow(meta_row) != 1) {

    stop(
      "Meta-regression result not found for parameter: ",
      parameter
    )
  }


  beta0 <- meta_row$beta0
  beta1 <- meta_row$beta1


  Vbeta <- matrix(
    c(
      meta_row$vc11,
      meta_row$vc12,
      meta_row$vc12,
      meta_row$vc22
    ),
    nrow = 2,
    byrow = TRUE
  )


  tau2 <- meta_row$tau2


  bind_rows(
    lapply(
      sizes_nm,
      function(size_nm) {

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


        var_mean <- as.numeric(
          t(x) %*%
            Vbeta %*%
            x
        )


        var_prediction <- var_mean


        if (type == "new_id") {

          var_prediction <-
            var_prediction +
            tau2
        }


        sd_prediction <- sqrt(
          pmax(
            var_prediction,
            0
          )
        )


        median_nat <- exp(
          mu_log
        )


        mean_nat <- exp(
          mu_log +
            0.5 *
            var_prediction
        )


        tibble(
          Parameter = parameter,

          Size_nm = size_nm,

          prediction_type = type,

          mean_log = mu_log,

          variance_log =
            var_prediction,

          median =
            median_nat,

          mean =
            mean_nat,

          q025 =
            exp(
              mu_log -
                1.96 *
                sd_prediction
            ),

          q975 =
            exp(
              mu_log +
                1.96 *
                sd_prediction
            )
        )
      }
    )
  )
}


predict_all_meta_parameters <- function(
    sizes_nm,
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
      META_PARAMS,
      function(parameter) {

        predict_meta_parameter(
          parameter = parameter,
          sizes_nm = sizes_nm,
          type = type
        )
      }
    )
  )
}


## ============================================================
## 11) Generate size-prediction curves
## ============================================================

META_PRED_POPULATION <- predict_all_meta_parameters(
  sizes_nm = SIZE_GRID,
  type = "population"
)


META_PRED_NEW_ID <- predict_all_meta_parameters(
  sizes_nm = SIZE_GRID,
  type = "new_id"
)


write_csv(
  META_PRED_POPULATION,
  file.path(
    META_DIR,
    "meta_prediction_curves_population.csv"
  )
)


write_csv(
  META_PRED_NEW_ID,
  file.path(
    META_DIR,
    "meta_prediction_curves_new_id.csv"
  )
)


## ============================================================
## 12) Predictions at selected particle sizes
## ============================================================

PREDICTION_SIZES <- c(
  20,
  50,
  220,
  500,
  1000,
  5000,
  6000
)


META_PRED_EVAL_POPULATION <-
  predict_all_meta_parameters(
    sizes_nm = PREDICTION_SIZES,
    type = "population"
  )


META_PRED_EVAL_NEW_ID <-
  predict_all_meta_parameters(
    sizes_nm = PREDICTION_SIZES,
    type = "new_id"
  )


write_csv(
  META_PRED_EVAL_POPULATION,
  file.path(
    META_DIR,
    "meta_predictions_selected_sizes_population.csv"
  )
)


write_csv(
  META_PRED_EVAL_NEW_ID,
  file.path(
    META_DIR,
    "meta_predictions_selected_sizes_new_id.csv"
  )
)


## ============================================================
## 13) Pool non-meta MCMC parameters across development IDs
##
## Posterior draws are pooled across IDs 1–5.
## ============================================================

combine_draws_for_id <- function(
    mcmc_result,
    id
) {

  chains <- mcmc_result$chains


  if (is.null(chains)) {

    stop(
      "Chains missing from MCMC result for ID ",
      id
    )
  }


  draws <- do.call(
    rbind,
    lapply(
      chains,
      function(chain) {

        as.data.frame(
          chain$pars
        )
      }
    )
  )


  draws$ID <- id

  draws
}


ALL_MCMC_DRAWS <- bind_rows(
  lapply(
    ids,
    function(id) {

      combine_draws_for_id(
        mcmc_result =
          MCMC_RESULTS[[as.character(id)]],
        id = id
      )
    }
  )
)


ALL_SAMPLED_PARAMS <- setdiff(
  names(
    ALL_MCMC_DRAWS
  ),
  "ID"
)


GLOBAL_PARAMS <- setdiff(
  ALL_SAMPLED_PARAMS,
  META_PARAMS
)


GLOBAL_SUMMARY <- bind_rows(
  lapply(
    GLOBAL_PARAMS,
    function(parameter) {

      log_draws <- ALL_MCMC_DRAWS[[parameter]]

      nat_draws <- exp(
        log_draws
      )


      tibble(
        Parameter = parameter,

        mean_log = mean(
          log_draws,
          na.rm = TRUE
        ),

        sd_log = sd(
          log_draws,
          na.rm = TRUE
        ),

        mean = mean(
          nat_draws,
          na.rm = TRUE
        ),

        median = median(
          nat_draws,
          na.rm = TRUE
        ),

        q025 = as.numeric(
          quantile(
            nat_draws,
            0.025,
            na.rm = TRUE
          )
        ),

        q975 = as.numeric(
          quantile(
            nat_draws,
            0.975,
            na.rm = TRUE
          )
        )
      )
    }
  )
)


write_csv(
  GLOBAL_SUMMARY,
  file.path(
    META_DIR,
    "pooled_global_parameter_summary.csv"
  )
)


cat("\nPooled non-meta parameters:\n")
print(
  GLOBAL_SUMMARY
)


## ============================================================
## 14) Save consolidated meta-analysis object
## ============================================================

META_OUTPUTS <- list(

  ids =
    ids,

  size_tbl =
    size_tbl,

  meta_params =
    META_PARAMS,

  global_params =
    GLOBAL_PARAMS,

  log_moments =
    LOG_MOMENTS,

  meta_summary =
    META_SUMMARY,

  meta_fit_list =
    META_FIT_LIST,

  global_summary =
    GLOBAL_SUMMARY,

  prediction_population =
    META_PRED_POPULATION,

  prediction_new_id =
    META_PRED_NEW_ID,

  selected_size_population =
    META_PRED_EVAL_POPULATION,

  selected_size_new_id =
    META_PRED_EVAL_NEW_ID
)


saveRDS(
  META_OUTPUTS,
  file.path(
    META_DIR,
    "meta_analysis_outputs.rds"
  )
)


## ============================================================
## 15) Diagnostic meta-regression plots
##
## These are analysis diagnostics.
## Final manuscript figures can be generated later by
## a dedicated figure script.
## ============================================================

for (parameter in META_PARAMS) {

  observed_parameter <- LOG_MOMENTS %>%
    filter(
      as.character(
        Parameter
      ) ==
        parameter
    ) %>%
    mutate(
      posterior_median =
        exp(
          meanlog
        ),

      posterior_lo =
        exp(
          meanlog -
            1.96 *
            sdlog
        ),

      posterior_hi =
        exp(
          meanlog +
            1.96 *
            sdlog
        )
    )


  population_curve <- META_PRED_POPULATION %>%
    filter(
      Parameter ==
        parameter
    )


  prediction_curve <- META_PRED_NEW_ID %>%
    filter(
      Parameter ==
        parameter
    )


  p <- ggplot() +

    geom_ribbon(
      data =
        prediction_curve,
      aes(
        x = Size_nm,
        ymin = q025,
        ymax = q975
      ),
      alpha = 0.15
    ) +

    geom_ribbon(
      data =
        population_curve,
      aes(
        x = Size_nm,
        ymin = q025,
        ymax = q975
      ),
      alpha = 0.30
    ) +

    geom_line(
      data =
        population_curve,
      aes(
        x = Size_nm,
        y = median
      ),
      linewidth = 1
    ) +

    geom_errorbar(
      data =
        observed_parameter,
      aes(
        x = Size_nm,
        ymin =
          posterior_lo,
        ymax =
          posterior_hi
      ),
      width = 0
    ) +

    geom_point(
      data =
        observed_parameter,
      aes(
        x = Size_nm,
        y =
          posterior_median
      ),
      size = 2.5
    ) +

    scale_x_log10() +

    scale_y_log10() +

    labs(
      title =
        parameter,
      x =
        "Particle size (nm)",
      y =
        parameter
    ) +

    theme_bw(
      base_size = 12
    )


  ggsave(
    filename = file.path(
      FIGURE_DIR,
      paste0(
        "meta_regression_",
        parameter,
        ".tiff"
      )
    ),

    plot =
      p,

    device =
      "tiff",

    width =
      6,

    height =
      5,

    units =
      "in",

    dpi =
      600,

    compression =
      "lzw"
  )
}


## ============================================================
## 16) Final summary
## ============================================================

cat("\n============================================================\n")
cat("Random-effects meta-regression completed successfully\n")
cat("============================================================\n")

cat(
  "Development datasets: ",
  length(ids),
  "\n",
  sep = ""
)

cat(
  "Meta-regression parameters: ",
  paste(
    META_PARAMS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Table S9 source: ",
  file.path(
    META_DIR,
    "Table_S9_meta_regression.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "Meta-analysis object: ",
  file.path(
    META_DIR,
    "meta_analysis_outputs.rds"
  ),
  "\n",
  sep = ""
)

cat(
  "Outputs saved to: ",
  META_DIR,
  "\n",
  sep = ""
)

cat("============================================================\n")