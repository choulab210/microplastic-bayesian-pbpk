## ============================================================
## 05_simulate_pbpk_from_meta.R
## Use meta-regression outputs to build PBPK parameter vectors
## and simulate concentration-time profiles for new particle sizes.
##
## External test set:
## - IDs 6, 7, 8 from data.csv are used for validation only
##
## Updated features:
## - deterministic simulation for test IDs
## - uncertainty propagation for predicted values using meta parameter uncertainty
## - bar chart comparison with:
##     * Observed bars ± SD
##     * Predicted bars ± 95% CI
## - GI separated from other tissues
##
## Inputs:
## - results/fitted_objects.rds
## - results/meta/meta_analysis_outputs.rds
## - MP_PBPK.R
## - data.csv
##
## Outputs:
## - results/meta/pbpk_meta_predictions_testIDs_*.csv
## - results/meta/pbpk_meta_pred_intervals_testIDs_*.csv
## - results/meta/pbpk_external_test_metrics_*.csv
## - results/meta/pbpk_external_test_singlepoint_metrics.csv
## - figures/PBPK_meta_profiles_testIDs_*.tiff
## - figures/PBPK_meta_bar_nonGI_testIDs_*.tiff
## - figures/PBPK_meta_bar_GI_testIDs_*.tiff
## ============================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tools)
  library(readr)
  library(tibble)
})

options(stringsAsFactors = FALSE)

## ============================================================
## 0) USER SETTINGS
## ============================================================
results_dir <- "results"
meta_dir    <- file.path(results_dir, "meta")
figures_dir <- "figures"

dir.create(results_dir, showWarnings = FALSE)
dir.create(meta_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# external test IDs
TEST_IDS <- 6:8

# if NULL, inferred from test IDs' Size_nm
SIM_SIZES_NM <- NULL

# Which meta prediction to use for size-dependent parameters
# "population" = fixed-effect uncertainty only
# "new_id"     = includes between-ID heterogeneity tau2
META_PRED_TYPE <- "new_id"

# Which summary statistic to use for deterministic parameter vector
# "mean", "q025", "q975"
META_STATISTIC <- "mean"

# simulation dose/time
DEFAULT_DOSE_MG <- NULL   # NULL = use dose from each test ID if available
DEFAULT_TEND_H  <- NULL   # NULL = use max time from each test ID if available

# simulation step
SIM_DT_H <- 0.05

# number of simulations for uncertainty propagation
N_SIM_UNC <- 300

# plotting
USE_LOG_Y <- FALSE
Y_FLOOR   <- 1e-8

## ============================================================
## 1) Load required objects
## ============================================================
fit_file  <- file.path(results_dir, "fitted_objects.rds")
meta_file <- file.path(meta_dir, "meta_analysis_outputs.rds")

if (!file.exists(fit_file)) {
  stop("Missing file: ", fit_file, "\nPlease run 01_run_model_fit.R first.")
}
if (!file.exists(meta_file)) {
  stop("Missing file: ", meta_file, "\nPlease run 04_run_meta_analysis.R first.")
}

fit_obj  <- readRDS(fit_file)
meta_obj <- readRDS(meta_file)

required_fit <- c("theta.int")
missing_fit <- setdiff(required_fit, names(fit_obj))
if (length(missing_fit) > 0) {
  stop("fitted_objects.rds is missing: ", paste(missing_fit, collapse = ", "))
}

required_meta <- c("SIZE_PARAMS", "GLOBAL_SUM", "META_FITS")
missing_meta <- setdiff(required_meta, names(meta_obj))
if (length(missing_meta) > 0) {
  stop("meta_analysis_outputs.rds is missing: ", paste(missing_meta, collapse = ", "))
}

theta.int    <- fit_obj$theta.int
SIZE_PARAMS  <- meta_obj$SIZE_PARAMS
GLOBAL_SUM   <- meta_obj$GLOBAL_SUM
META_FITS    <- meta_obj$META_FITS

## ============================================================
## 2) Load full dataset for external test IDs
## ============================================================
## ============================================================
## 2) Load full dataset for external test IDs
## ============================================================
obs_all <- read.csv("data.csv", encoding = "latin1") %>%
  mutate(
    Tissue = trimws(as.character(Tissue)),
    Tissue = tolower(Tissue),
    Tissue = dplyr::recode(
      Tissue,
      "lung" = "lungs",
      "gi tract" = "gi",
      "g.i." = "gi",
      "gitract" = "gi",
      #"stomach" = "gi",
      "intestine" = "gi",
      "colon" = "gi",
      .default = Tissue
    ),
    Tissue = tools::toTitleCase(Tissue),
    Tissue = dplyr::recode(
      Tissue,
      "Lung" = "Lungs",
      "Gi" = "GI"
    ),
    Tissue = factor(
      Tissue,
      levels = c("Blood","Liver","Kidney","Spleen","Lungs","GI","Urine")
    ),
    SD = ifelse(is.na(SD), pmax(Conc * 0.1, 1e-6), SD)
  )

TEST_IDS <- 6:17

test_df <- obs_all %>%
  filter(ID %in% TEST_IDS)

if (nrow(test_df) == 0) {
  stop("No rows found for TEST_IDS = ", paste(TEST_IDS, collapse = ", "))
}

if (!"Size_nm" %in% names(test_df)) {
  stop("data.csv must contain a Size_nm column for test IDs.")
}

test_sizes_tbl <- test_df %>%
  group_by(ID) %>%
  summarise(
    Size_nm = first(na.omit(as.numeric(Size_nm))),
    dose_mg = first(na.omit(as.numeric(dose_mg))),
    t_end_h = max(as.numeric(Time_h), na.rm = TRUE),
    n_time  = n_distinct(as.numeric(Time_h)),
    .groups = "drop"
  )

# keep only single-time-point held-out IDs
test_sizes_tbl <- test_sizes_tbl %>%
  filter(n_time == 1)

TEST_IDS <- sort(test_sizes_tbl$ID)

test_df <- test_df %>%
  filter(ID %in% TEST_IDS)

if (nrow(test_df) == 0) {
  stop("No single-time-point held-out IDs found in TEST_IDS.")
}

if (any(!is.finite(test_sizes_tbl$Size_nm))) {
  stop("Could not determine Size_nm for one or more test IDs.")
}

if (is.null(SIM_SIZES_NM)) {
  SIM_SIZES_NM <- sort(unique(test_sizes_tbl$Size_nm))
}

cat("\n=== External test IDs and sizes ===\n")
print(test_sizes_tbl)
cat("Held-out evaluation IDs:", paste(TEST_IDS, collapse = ", "), "\n")

## ============================================================
## 3) Meta helper functions
## ============================================================
predict_param_by_size <- function(param, sizes_nm, meta_fits = META_FITS,
                                  type = c("population", "new_id")) {
  type <- match.arg(type)
  
  row <- meta_fits %>% filter(param == !!param)
  if (nrow(row) != 1) stop("Meta fit not found for param: ", param)
  
  b0   <- row$coef0
  b1   <- row$coef1
  Vb   <- matrix(c(row$vc11, row$vc12, row$vc12, row$vc22), 2, 2)
  tau2 <- row$tau2
  
  bind_rows(lapply(sizes_nm, function(S) {
    x   <- c(1, log(S))
    mu  <- as.numeric(b0 + b1 * log(S))
    var <- as.numeric(t(x) %*% Vb %*% x)
    
    if (type == "new_id") var <- var + tau2
    
    tibble(
      Parameter = param,
      Size_nm   = S,
      mean      = exp(mu + 0.5 * var),
      q025      = exp(mu - 1.96 * sqrt(var)),
      q975      = exp(mu + 1.96 * sqrt(var))
    )
  }))
}

predict_params_by_size <- function(params, sizes_nm, type = c("population", "new_id")) {
  type <- match.arg(type)
  bind_rows(lapply(params, function(p) {
    predict_param_by_size(p, sizes_nm = sizes_nm, meta_fits = META_FITS, type = type)
  }))
}

build_param_vector_for_size <- function(size_nm,
                                        size_params = SIZE_PARAMS,
                                        global_sum  = GLOBAL_SUM,
                                        theta_base_log = theta.int,
                                        type = c("population", "new_id"),
                                        statistic = c("mean", "q025", "q975")) {
  type      <- match.arg(type)
  statistic <- match.arg(statistic)
  
  stat_col <- switch(match.arg(statistic),
                     mean = "mean",
                     q025 = "q025",
                     q975 = "q975")
  
  pars_nat <- exp(theta_base_log)
  
  SZ <- predict_params_by_size(size_params, sizes_nm = size_nm, type = type) %>%
    filter(Size_nm == size_nm)
  
  if (nrow(SZ) > 0) {
    pars_nat[SZ$Parameter] <- SZ[[stat_col]]
  }
  
  if (nrow(global_sum) > 0) {
    pars_nat[global_sum$Parameter] <- global_sum[[stat_col]]
  }
  
  pars_nat
}

## ============================================================
## 3b) Parameter uncertainty sampling helpers
## ============================================================
sample_size_param_one <- function(param, size_nm, meta_fits = META_FITS,
                                  type = c("population", "new_id")) {
  type <- match.arg(type)
  
  row <- meta_fits %>% filter(param == !!param)
  if (nrow(row) != 1) stop("Meta fit not found for param: ", param)
  
  b0   <- row$coef0
  b1   <- row$coef1
  Vb   <- matrix(c(row$vc11, row$vc12, row$vc12, row$vc22), 2, 2)
  tau2 <- row$tau2
  
  x  <- c(1, log(size_nm))
  mu <- as.numeric(b0 + b1 * log(size_nm))
  var_mu <- as.numeric(t(x) %*% Vb %*% x)
  
  if (type == "new_id") {
    var_tot <- var_mu + tau2
  } else {
    var_tot <- var_mu
  }
  
  exp(rnorm(1, mean = mu, sd = sqrt(pmax(var_tot, 1e-12))))
}

sample_global_param_one <- function(param, global_sum = GLOBAL_SUM) {
  row <- global_sum %>% filter(Parameter == !!param)
  if (nrow(row) != 1) stop("Global summary not found for param: ", param)
  
  m  <- row$mean
  ql <- row$q025
  qu <- row$q975
  
  if (!is.finite(ql) || !is.finite(qu) || ql <= 0 || qu <= 0 || qu <= ql) {
    return(as.numeric(m))
  }
  
  sdlog <- (log(qu) - log(ql)) / (2 * 1.96)
  sdlog <- pmax(sdlog, 1e-8)
  mulog <- log(pmax(m, 1e-12)) - 0.5 * sdlog^2
  
  exp(rnorm(1, mean = mulog, sd = sdlog))
}

build_param_draw_for_size <- function(size_nm,
                                      size_params = SIZE_PARAMS,
                                      global_sum = GLOBAL_SUM,
                                      theta_base_log = theta.int,
                                      type = c("population", "new_id")) {
  type <- match.arg(type)
  
  pars_nat <- exp(theta_base_log)
  
  for (p in size_params) {
    pars_nat[p] <- sample_size_param_one(
      param = p,
      size_nm = size_nm,
      meta_fits = META_FITS,
      type = type
    )
  }
  
  if (nrow(global_sum) > 0) {
    for (p in global_sum$Parameter) {
      pars_nat[p] <- sample_global_param_one(
        param = p,
        global_sum = global_sum
      )
    }
  }
  
  pars_nat
}

## ============================================================
## 4) Load and compile PBPK model
## ============================================================
source("MP_PBPK.R")
mod <- mcode("mPBPK_RES", mPBPK_RES)
MODEL_PARAM_NAMES <- names(param(mod))

## ============================================================
## 5) PBPK simulation helpers
## ============================================================
pred_from_natpars <- function(pars_nat,
                              dose_mg,
                              t_end_h,
                              dt_h = SIM_DT_H,
                              atol = 1e-8,
                              rtol = 1e-6,
                              add_times = NULL) {
  pars_nat <- pars_nat[names(pars_nat) %in% MODEL_PARAM_NAMES]
  
  mod_tmp <- mod
  if (length(pars_nat) > 0) {
    mod_tmp <- param(mod_tmp, pars_nat)
  }
  
  ex <- ev(amt = dose_mg, cmt = "A_lumen", time = 0)
  
  ts <- tgrid(
    start = 0,
    end   = t_end_h,
    delta = dt_h,
    add   = add_times
  )
  
  out <- mrgsim_d(
    mod_tmp, data = ex, tgrid = ts,
    atol = atol, rtol = rtol, maxsteps = 1e6
  ) %>% as.data.frame()
  
  dose_total <- dose_mg
  
  lungs_total <- if ("total_lung" %in% names(out)) {
    out$total_lung
  } else if (all(c("A_lung_b", "A_lung_t") %in% names(out))) {
    out$A_lung_b + out$A_lung_t
  } else {
    0
  }
  
  tibble(
    time   = out$time,
    Blood  = (out$A_blood      / dose_total) * 100,
    Liver  = (out$total_liver  / dose_total) * 100,
    Kidney = (out$total_kidney / dose_total) * 100,
    Spleen = (out$total_spleen / dose_total) * 100,
    Lungs  = (lungs_total      / dose_total) * 100,
    GI     = (out$total_GI     / dose_total) * 100,
    Urine  = (out$A_urine      / dose_total) * 100,
    Rest   = (out$total_rest   / dose_total) * 100
  )
}

calc_gof_metrics <- function(df) {
  df <- df %>%
    filter(is.finite(obs), is.finite(pred), obs > 0, pred > 0)
  
  if (nrow(df) == 0) {
    return(data.frame(R2 = NA_real_, RMSE = NA_real_, MAPE = NA_real_, n = 0))
  }
  
  r2 <- if (length(unique(df$obs)) > 1 && length(unique(df$pred)) > 1) {
    cor(df$obs, df$pred, use = "complete.obs")^2
  } else {
    NA_real_
  }
  
  rmse <- sqrt(mean((df$pred - df$obs)^2, na.rm = TRUE))
  mape <- mean(abs(df$pred - df$obs) / pmax(df$obs, 1e-8), na.rm = TRUE) * 100
  
  data.frame(R2 = r2, RMSE = rmse, MAPE = mape, n = nrow(df))
}

predictive_interval_for_test_id <- function(id_i, size_nm, dose_i, tend_i, add_times_i,
                                            n_sim = N_SIM_UNC,
                                            type = META_PRED_TYPE) {
  sim_list <- vector("list", n_sim)
  
  for (j in seq_len(n_sim)) {
    pars_j <- build_param_draw_for_size(
      size_nm = size_nm,
      size_params = SIZE_PARAMS,
      global_sum = GLOBAL_SUM,
      theta_base_log = theta.int,
      type = type
    )
    
    sim_j <- pred_from_natpars(
      pars_nat = pars_j,
      dose_mg = dose_i,
      t_end_h = tend_i,
      dt_h = SIM_DT_H,
      add_times = add_times_i
    ) %>%
      pivot_longer(
        cols = c("Blood", "Liver", "Kidney", "Spleen", "Lungs", "GI", "Urine"),
        names_to = "Tissue",
        values_to = "pred"
      ) %>%
      mutate(
        ID = id_i,
        Size_nm = size_nm,
        draw = j,
        .before = 1
      )
    
    sim_list[[j]] <- sim_j
  }
  
  bind_rows(sim_list) %>%
    group_by(ID, Size_nm, time, Tissue) %>%
    summarise(
      pred_mean = mean(pred, na.rm = TRUE),
      pred_lo   = quantile(pred, 0.025, na.rm = TRUE),
      pred_hi   = quantile(pred, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

## ============================================================
## 6) Deterministic simulation by external test ID
## ============================================================
sim_test_list <- lapply(seq_len(nrow(test_sizes_tbl)), function(i) {
  row_i <- test_sizes_tbl[i, ]
  id_i  <- row_i$ID
  sz_i  <- row_i$Size_nm
  
  dose_i <- if (!is.null(DEFAULT_DOSE_MG)) DEFAULT_DOSE_MG else row_i$dose_mg
  tend_i <- if (!is.null(DEFAULT_TEND_H))  DEFAULT_TEND_H  else row_i$t_end_h
  
  add_times_i <- sort(unique(as.numeric(test_df$Time_h[test_df$ID == id_i])))
  
  pars_i <- build_param_vector_for_size(
    size_nm = sz_i,
    size_params = SIZE_PARAMS,
    global_sum = GLOBAL_SUM,
    theta_base_log = theta.int,
    type = META_PRED_TYPE,
    statistic = META_STATISTIC
  )
  
  sim_i <- pred_from_natpars(
    pars_nat = pars_i,
    dose_mg  = dose_i,
    t_end_h  = tend_i,
    dt_h     = SIM_DT_H,
    add_times = add_times_i
  ) %>%
    mutate(ID = id_i, Size_nm = sz_i, dose_mg = dose_i, .before = 1)
  
  sim_i
})

sim_test_all <- bind_rows(sim_test_list)

write_csv(
  sim_test_all,
  file.path(meta_dir, paste0("pbpk_meta_predictions_testIDs_", META_PRED_TYPE, "_", META_STATISTIC, ".csv"))
)

param_vectors_test <- bind_rows(lapply(seq_len(nrow(test_sizes_tbl)), function(i) {
  row_i <- test_sizes_tbl[i, ]
  pv <- build_param_vector_for_size(
    size_nm = row_i$Size_nm,
    size_params = SIZE_PARAMS,
    global_sum = GLOBAL_SUM,
    theta_base_log = theta.int,
    type = META_PRED_TYPE,
    statistic = META_STATISTIC
  )
  
  tibble(
    ID = row_i$ID,
    Size_nm = row_i$Size_nm,
    Parameter = names(pv),
    Value = as.numeric(pv)
  )
}))

write_csv(
  param_vectors_test,
  file.path(meta_dir, paste0("pbpk_param_vectors_testIDs_", META_PRED_TYPE, "_", META_STATISTIC, ".csv"))
)

## ============================================================
## 6b) Uncertainty propagation for predicted values
## ============================================================
cat("\n=== Running uncertainty propagation for predicted bars ===\n")

pred_interval_list <- lapply(seq_len(nrow(test_sizes_tbl)), function(i) {
  row_i <- test_sizes_tbl[i, ]
  id_i  <- row_i$ID
  sz_i  <- row_i$Size_nm
  
  dose_i <- if (!is.null(DEFAULT_DOSE_MG)) DEFAULT_DOSE_MG else row_i$dose_mg
  tend_i <- if (!is.null(DEFAULT_TEND_H))  DEFAULT_TEND_H  else row_i$t_end_h
  add_times_i <- sort(unique(as.numeric(test_df$Time_h[test_df$ID == id_i])))
  
  cat("  Uncertainty propagation for ID ", id_i, " (size=", sz_i, " nm)\n", sep = "")
  
  predictive_interval_for_test_id(
    id_i = id_i,
    size_nm = sz_i,
    dose_i = dose_i,
    tend_i = tend_i,
    add_times_i = add_times_i,
    n_sim = N_SIM_UNC,
    type = META_PRED_TYPE
  )
})

pred_interval_df <- bind_rows(pred_interval_list)

write_csv(
  pred_interval_df,
  file.path(meta_dir, paste0("pbpk_meta_pred_intervals_testIDs_", META_PRED_TYPE, ".csv"))
)

## ============================================================
## 7) Compare meta-PBPK predictions vs observed test IDs
## ============================================================
obs_test_long <- test_df %>%
  transmute(
    ID = as.numeric(ID),
    Size_nm = as.numeric(Size_nm),
    time = as.numeric(Time_h),
    Tissue = as.character(Tissue),
    obs = as.numeric(Conc),
    SD = as.numeric(SD)
  )

sim_test_long <- sim_test_all %>%
  pivot_longer(
    cols = c("Blood", "Liver", "Kidney", "Spleen", "Lungs", "GI", "Urine"),
    names_to = "Tissue",
    values_to = "pred"
  ) %>%
  transmute(
    ID = as.numeric(ID),
    Size_nm = as.numeric(Size_nm),
    time = as.numeric(time),
    Tissue = as.character(Tissue),
    pred = as.numeric(pred)
  )

test_compare <- obs_test_long %>%
  inner_join(sim_test_long, by = c("ID", "Size_nm", "time", "Tissue")) %>%
  left_join(
    pred_interval_df,
    by = c("ID", "Size_nm", "time", "Tissue")
  )

external_metrics_overall <- calc_gof_metrics(test_compare)
external_metrics_by_id <- bind_rows(lapply(sort(unique(test_compare$ID)), function(id_i) {
  cbind(data.frame(ID = id_i), calc_gof_metrics(filter(test_compare, ID == id_i)))
}))
external_metrics_by_tissue <- bind_rows(lapply(sort(unique(test_compare$Tissue)), function(tt) {
  cbind(data.frame(Tissue = tt), calc_gof_metrics(filter(test_compare, Tissue == tt)))
}))

write_csv(external_metrics_overall,  file.path(meta_dir, "pbpk_external_test_metrics_overall.csv"))
write_csv(external_metrics_by_id,    file.path(meta_dir, "pbpk_external_test_metrics_by_id.csv"))
write_csv(external_metrics_by_tissue,file.path(meta_dir, "pbpk_external_test_metrics_by_tissue.csv"))

cat("\n=== External test metrics: overall ===\n")
print(external_metrics_overall)

cat("\n=== External test metrics: by ID ===\n")
print(external_metrics_by_id)

cat("\n=== External test metrics: by Tissue ===\n")
print(external_metrics_by_tissue)

external_singlepoint_metrics <- test_compare %>%
  mutate(
    AE = abs(pred - obs),
    RE = abs(pred - obs) / pmax(obs, 1e-8),
    Ratio = pred / pmax(obs, 1e-8)
  ) %>%
  group_by(ID, Size_nm, time, Tissue) %>%
  summarise(
    Obs = mean(obs, na.rm = TRUE),
    Pred = mean(pred, na.rm = TRUE),
    PredMean = mean(pred_mean, na.rm = TRUE),
    PredLo = mean(pred_lo, na.rm = TRUE),
    PredHi = mean(pred_hi, na.rm = TRUE),
    AE = mean(AE, na.rm = TRUE),
    RE_percent = mean(RE, na.rm = TRUE) * 100,
    Ratio = mean(Ratio, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  external_singlepoint_metrics,
  file.path(meta_dir, "pbpk_external_test_singlepoint_metrics.csv")
)

cat("\n=== External single-point metrics ===\n")
print(external_singlepoint_metrics)

## ============================================================
## 8) Plot deterministic simulated profiles for test IDs
## ============================================================
p_profiles <- ggplot(
  sim_test_long,
  aes(x = time, y = pmax(pred, Y_FLOOR), color = factor(ID))
) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Tissue, scales = "free_y") +
  labs(
    title = paste0(
      "PBPK predictions from meta-regression parameters\n",
      "External test IDs: ", paste(TEST_IDS, collapse = ", "),
      " | Type = ", META_PRED_TYPE,
      " | statistic = ", META_STATISTIC
    ),
    x = "Time (h)",
    y = "%ID",
    color = "Test ID"
  ) +
  theme_bw(base_size = 13)

if (USE_LOG_Y) {
  p_profiles <- p_profiles + scale_y_log10()
}

ggsave(
  filename = file.path(figures_dir, paste0("PBPK_meta_profiles_testIDs_", META_PRED_TYPE, "_", META_STATISTIC, ".tiff")),
  plot = p_profiles,
  device = "tiff",
  width = 12,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

## ============================================================
## 9) Bar chart comparison (GI separated; predicted bars get 95% CI)
##    Upper error bars only
## ============================================================

## ----- Observed values -----
obs_bar_df <- test_compare %>%
  transmute(
    ID, Size_nm, time, Tissue,
    Source = "Observed",
    Value  = obs,
    ErrHi  = obs + SD
  )

## ----- Predicted values -----
pred_bar_df <- test_compare %>%
  transmute(
    ID, Size_nm, time, Tissue,
    Source = "Predicted",
    Value  = pred_mean,
    ErrHi  = pred_hi
  )

## ----- Combine -----
compare_bar_df <- bind_rows(obs_bar_df, pred_bar_df) %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = c("Blood","Liver","Kidney","Spleen","Lungs","GI","Urine")
    ),
    Source = factor(Source, levels = c("Observed","Predicted"))
  ) %>%
  group_by(ID, Size_nm, time, Tissue, Source) %>%
  summarise(
    Value = mean(Value, na.rm = TRUE),
    ErrHi = mean(ErrHi, na.rm = TRUE),
    .groups = "drop"
  )

## ============================================================
## Non-GI tissues
## ============================================================

compare_bar_nonGI <- compare_bar_df %>%
  filter(Tissue != "GI") %>%
  group_by(Tissue) %>%
  filter(max(Value, na.rm = TRUE) > 0) %>%
  ungroup()

pd <- position_dodge(width = 0.8)

p_bar_nonGI <- ggplot(compare_bar_nonGI,
                      aes(x = Tissue, y = Value, fill = Source)) +
  geom_col(position = pd, width = 0.7) +
  geom_errorbar(
    aes(ymin = Value, ymax = ErrHi, group = Source),
    position = pd,
    width = 0.25,
    linewidth = 0.7
  ) +
  facet_wrap(~ID, nrow = 1, scales = "free_y") +
  labs(
    title = paste0(
      "Observed vs meta-PBPK predicted values for external test IDs (",
      paste(TEST_IDS, collapse = ", "),
      ") — Non-GI tissues"
    ),
    x = "Tissue",
    y = "%ID",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  )

if (USE_LOG_Y) {
  p_bar_nonGI <- p_bar_nonGI + scale_y_log10()
}

ggsave(
  filename = file.path(
    figures_dir,
    paste0(
      "PBPK_meta_bar_nonGI_testIDs_",
      paste(TEST_IDS, collapse = "_"),
      "_",
      META_PRED_TYPE,
      "_",
      META_STATISTIC,
      ".tiff"
    )
  ),
  plot = p_bar_nonGI,
  device = "tiff",
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

## ============================================================
## GI tissue only
## ============================================================

compare_bar_GI <- compare_bar_df %>%
  filter(Tissue == "GI")

p_bar_GI <- ggplot(compare_bar_GI,
                   aes(x = Source, y = Value, fill = Source)) +
  geom_col(width = 0.65) +
  geom_errorbar(
    aes(ymin = Value, ymax = ErrHi),
    width = 0.20,
    linewidth = 0.7
  ) +
  facet_wrap(~ID, nrow = 1, scales = "free_y") +
  labs(
    title = paste0(
      "Observed vs meta-PBPK predicted values for external test IDs (",
      paste(TEST_IDS, collapse = ", "),
      ") — GI"
    ),
    x = NULL,
    y = "%ID",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

if (USE_LOG_Y) {
  p_bar_GI <- p_bar_GI + scale_y_log10()
}

ggsave(
  filename = file.path(
    figures_dir,
    paste0(
      "PBPK_meta_bar_GI_testIDs_",
      paste(TEST_IDS, collapse = "_"),
      "_",
      META_PRED_TYPE,
      "_",
      META_STATISTIC,
      ".tiff"
    )
  ),
  plot = p_bar_GI,
  device = "tiff",
  width = 8,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

## ============================================================
## 10) Save summary object
## ============================================================
saveRDS(
  list(
    settings = list(
      TEST_IDS = TEST_IDS,
      META_PRED_TYPE = META_PRED_TYPE,
      META_STATISTIC = META_STATISTIC,
      SIM_DT_H = SIM_DT_H,
      N_SIM_UNC = N_SIM_UNC
    ),
    test_sizes_tbl = test_sizes_tbl,
    param_vectors_test = param_vectors_test,
    sim_test_all = sim_test_all,
    pred_interval_df = pred_interval_df,
    test_compare = test_compare,
    compare_bar_df = compare_bar_df,
    external_metrics_overall = external_metrics_overall,
    external_metrics_by_id = external_metrics_by_id,
    external_metrics_by_tissue = external_metrics_by_tissue,
    external_singlepoint_metrics = external_singlepoint_metrics
  ),
  file = file.path(meta_dir, paste0("pbpk_meta_external_test_outputs_", META_PRED_TYPE, "_", META_STATISTIC, ".rds"))
)

cat("\n============================================================\n")
cat("PBPK-from-meta external test completed.\n")
cat("Test IDs used: ", paste(TEST_IDS, collapse = ", "), "\n", sep = "")
cat("Saved outputs to:\n")
cat(" - ", normalizePath(meta_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat(" - ", normalizePath(figures_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")