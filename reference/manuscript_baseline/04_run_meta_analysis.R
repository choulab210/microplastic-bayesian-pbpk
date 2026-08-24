## ============================================================
## 04_run_meta_analysis.R
## Random-effects meta-regression (metafor, REML) from per-ID MCMC3
##
## Size model per parameter k:
##   log(theta_{i,k}) = eta0_k + eta1_k * log(Size_i) + u_{i,k},
##   u_{i,k} ~ N(0, tau_k^2)
##
## Outputs:
## - meta fit summaries
## - slope CI / p-value table
## - predicted parameter-size curves
## - PBPK-ready parameter vectors at requested sizes
## - trend plots
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(metafor)
})

options(stringsAsFactors = FALSE)

## ============================================================
## 0) Directories
## ============================================================
results_dir <- "results"
figures_dir <- "figures"
meta_dir    <- file.path(results_dir, "meta")

dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)
dir.create(meta_dir, showWarnings = FALSE)

## ============================================================
## 1) User options
## ============================================================
# intended size-dependent set
INTENDED_SIZE_PARAMS <- c("Ktrans","Kshed","Kgw2b","KLRESmax","Kfeces","KurineC")

# intended global params (pooled across IDs)
INTENDED_GLOBALS <- c(
  "KLRESrelease","KSRESmax","KSRESrelease",
  "KKRESmax","KKRESrelease","KLURESmax","KLURESrelease",
  "KRRESmax","KRRESrelease","PAL","PAK","PAS","PALU","PAGI","PAR","KbileC"
)

# default plotting grid
SIZE_GRID_DEFAULT <- exp(seq(log(20), log(6000), length.out = 120))

# fallback ID -> size mapping if obs_df lacks Size_nm
ID_SIZE_MAP_DEFAULT <- c(
  "1" = 20,
  "2" = 220,
  "3" = 1000,
  "4" = 6000,
  "5" = 220
)

## ============================================================
## 2) Load fitted objects and per-ID MCMC3
## ============================================================
get_fitted_objects <- function() {
  f <- file.path(results_dir, "fitted_objects.rds")
  if (!file.exists(f)) {
    stop("Missing file: ", f,
         "\nPlease run 01_run_model_fit.R first.")
  }
  readRDS(f)
}

get_MCMC3 <- function() {
  # prefer in-memory object
  if (exists("MCMC3", inherits = TRUE)) {
    message("Using MCMC3 from current environment.")
    return(get("MCMC3", inherits = TRUE))
  }
  
  # then look under results/mcmc
  cand1 <- list.files(
    path = file.path(results_dir, "mcmc"),
    pattern = "^MCMC3_results_fullset_perID_priors_.*\\.rds$",
    full.names = TRUE
  )
  
  # then current working directory
  cand2 <- list.files(
    pattern = "^MCMC3_results_fullset_perID_priors_.*\\.rds$",
    full.names = TRUE
  )
  
  candidates <- unique(c(cand1, cand2))
  
  if (length(candidates) == 0) {
    stop("No 'MCMC3_results_fullset_perID_priors_*.rds' file found.",
         "\nPlease run/save your per-ID MCMC first.")
  }
  
  message("Loading: ", candidates[1])
  readRDS(candidates[1])
}

fit_obj <- get_fitted_objects()
MCMC3   <- get_MCMC3()

if (!"obs_df" %in% names(fit_obj)) {
  stop("obs_df not found in fitted_objects.rds.")
}
obs_df <- fit_obj$obs_df

## ============================================================
## 3) IDs and sizes per ID
## ============================================================
ids <- sort(unique(obs_df$ID))
ids <- as.character(ids)

build_sizes_tbl <- function(obs_df, ids_chr, fallback_map = ID_SIZE_MAP_DEFAULT) {
  has_size_col <- "Size_nm" %in% names(obs_df)
  
  if (has_size_col) {
    sizes_tbl <- obs_df %>%
      mutate(ID = as.character(ID)) %>%
      group_by(ID) %>%
      summarise(Size_nm = suppressWarnings(first(na.omit(as.numeric(Size_nm)))),
                .groups = "drop")
  } else {
    sizes_tbl <- tibble(
      ID = ids_chr,
      Size_nm = as.numeric(fallback_map[ids_chr])
    )
  }
  
  # fill missing size values with fallback map if available
  sizes_tbl <- sizes_tbl %>%
    mutate(
      Size_nm = ifelse(
        is.na(Size_nm),
        as.numeric(fallback_map[ID]),
        Size_nm
      )
    ) %>%
    mutate(logSize = log(Size_nm))
  
  if (any(!is.finite(sizes_tbl$Size_nm))) {
    bad_ids <- sizes_tbl$ID[!is.finite(sizes_tbl$Size_nm)]
    stop("Could not determine Size_nm for IDs: ", paste(bad_ids, collapse = ", "),
         "\nPlease add Size_nm to obs_df or update ID_SIZE_MAP_DEFAULT.")
  }
  
  sizes_tbl
}

sizes_tbl <- build_sizes_tbl(obs_df, ids)

stopifnot(all(ids %in% sizes_tbl$ID))

cat("\n=== ID-size table ===\n")
print(sizes_tbl)

## ============================================================
## 4) Helper to combine chains for one ID
## ============================================================
combine_id_draws <- function(M3, id_chr) {
  id_chr <- as.character(id_chr)
  if (is.null(M3[[id_chr]])) stop("MCMC3 entry not found for ID: ", id_chr)
  
  out <- do.call(rbind, lapply(M3[[id_chr]]$chains, function(m) as.data.frame(m$pars)))
  if (is.null(out) || nrow(out) == 0) {
    stop("No MCMC draws found for ID: ", id_chr)
  }
  out
}

first_mat <- combine_id_draws(MCMC3, ids[1])
all_cols  <- colnames(first_mat)

SIZE_PARAMS   <- intersect(INTENDED_SIZE_PARAMS, all_cols)
GLOBAL_PARAMS <- intersect(INTENDED_GLOBALS, all_cols)

if (!length(SIZE_PARAMS)) {
  stop("No intended size-dependent params found in MCMC3 columns.")
}

cat("\n=== Size-dependent parameters used ===\n")
print(SIZE_PARAMS)

cat("\n=== Global pooled parameters used ===\n")
print(GLOBAL_PARAMS)

## ============================================================
## 5) Per-ID moments (log-scale) from chains
## ============================================================
LOG_MOM <- bind_rows(lapply(SIZE_PARAMS, function(p) {
  bind_rows(lapply(ids, function(id) {
    D <- combine_id_draws(MCMC3, id)
    if (!p %in% colnames(D)) return(NULL)
    
    tibble(
      ID      = as.character(id),
      param   = p,
      meanlog = mean(D[[p]], na.rm = TRUE),
      sdlog   = sd(D[[p]],   na.rm = TRUE)
    )
  }))
})) %>%
  mutate(ID = as.character(ID)) %>%
  left_join(sizes_tbl, by = "ID") %>%
  filter(is.finite(meanlog), is.finite(sdlog), is.finite(Size_nm), is.finite(logSize))

if (nrow(LOG_MOM) == 0) {
  stop("LOG_MOM is empty after filtering. Check MCMC3 content and size mapping.")
}

write_csv(LOG_MOM, file.path(meta_dir, "meta_log_moments_by_id.csv"))

## ============================================================
## 6) Meta-regression per parameter (metafor REML)
## ============================================================
fit_meta_param <- function(df_param) {
  df_param <- df_param %>%
    mutate(vi = pmax(sdlog^2, 1e-10))  # sampling variance proxy
  
  fit <- metafor::rma.uni(
    yi = meanlog,
    vi = vi,
    mods = ~ logSize,
    method = "REML",
    data = df_param,
    control = list(stepadj = 0.5, maxiter = 200)
  )
  
  list(
    engine = "metafor",
    coef   = as.numeric(coef(fit)),   # intercept, slope
    vcov   = as.matrix(vcov(fit)),
    tau2   = as.numeric(fit$tau2),
    n      = nrow(df_param),
    obj    = fit
  )
}

SPLIT <- split(LOG_MOM, LOG_MOM$param)

fit_list <- lapply(names(SPLIT), function(p) {
  res <- fit_meta_param(SPLIT[[p]])
  
  tibble(
    param  = p,
    engine = res$engine,
    n      = res$n,
    tau2   = res$tau2,
    coef0  = res$coef[1],
    coef1  = res$coef[2],
    vc11   = res$vcov[1, 1],
    vc22   = res$vcov[2, 2],
    vc12   = res$vcov[1, 2],
    obj    = list(res$obj)
  )
})

META_FITS <- bind_rows(fit_list)

cat("\n=== Meta-regression (REML) summary ===\n")
print(META_FITS %>% select(param, engine, n, tau2, coef0, coef1))

write_csv(
  META_FITS %>% select(param, engine, n, tau2, coef0, coef1, vc11, vc22, vc12),
  file.path(meta_dir, "META_FITS_summary.csv")
)

saveRDS(META_FITS, file.path(meta_dir, "META_FITS_metafor_REML.rds"))

## ============================================================
## 7) Slope CIs and p-values
## ============================================================
SLP <- META_FITS %>%
  transmute(
    param,
    slope = coef1,
    se    = sqrt(vc22),
    lo95  = coef1 - 1.96 * se,
    hi95  = coef1 + 1.96 * se,
    pval  = 2 * pnorm(-abs(coef1 / se))
  ) %>%
  arrange(pval)

cat("\n=== Slope tests (eta1) ===\n")
print(SLP)

write_csv(SLP, file.path(meta_dir, "META_slope_tests.csv"))

## ============================================================
## 8) Predict by size (population mean vs new-ID)
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

## ============================================================
## 9) Summaries for global params pooled across IDs
## ============================================================
summarize_globals_from_MCMC3 <- function(M3, ids_chr, params = GLOBAL_PARAMS) {
  if (!length(params)) return(tibble())
  
  ALL <- bind_rows(lapply(ids_chr, function(id) {
    D <- combine_id_draws(M3, id)
    keep <- intersect(params, colnames(D))
    if (!length(keep)) return(NULL)
    
    D[, keep, drop = FALSE] %>%
      mutate(.draw = row_number(), ID = id)
  }))
  
  bind_rows(lapply(params, function(g) {
    if (!g %in% colnames(ALL)) return(NULL)
    
    th <- exp(ALL[[g]])  # natural scale
    tibble(
      Parameter = g,
      mean = mean(th, na.rm = TRUE),
      q025 = quantile(th, 0.025, na.rm = TRUE),
      q975 = quantile(th, 0.975, na.rm = TRUE)
    )
  }))
}

GLOBAL_SUM <- summarize_globals_from_MCMC3(MCMC3, ids)

cat("\n=== Global pooled parameters (natural scale) ===\n")
print(GLOBAL_SUM)

write_csv(GLOBAL_SUM, file.path(meta_dir, "meta_global_params_metafor.csv"))

## ============================================================
## 10) PBPK-ready vector at a given size (natural scale)
## ============================================================
build_param_vector_for_size <- function(size_nm,
                                        size_params = SIZE_PARAMS,
                                        global_sum  = GLOBAL_SUM,
                                        type = c("population", "new_id"),
                                        statistic = c("mean", "q025", "q975")) {
  type      <- match.arg(type)
  statistic <- match.arg(statistic)
  
  SZ <- predict_params_by_size(size_params, sizes_nm = size_nm, type = type) %>%
    filter(Size_nm == size_nm) %>%
    select(Parameter, mean, q025, q975)
  
  stat_col <- switch(statistic,
                     mean = "mean",
                     q025 = "q025",
                     q975 = "q975")
  
  sp <- setNames(SZ[[stat_col]], SZ$Parameter)
  
  if (nrow(global_sum) > 0) {
    gg <- setNames(global_sum$mean, global_sum$Parameter)
    c(sp, gg)
  } else {
    sp
  }
}

## ============================================================
## 11) Plot size-trend curves
## ============================================================
plot_size_trends <- function(params = SIZE_PARAMS,
                             sizes_grid = SIZE_GRID_DEFAULT) {
  lapply(params, function(p) {
    PR <- predict_param_by_size(p, sizes_grid, META_FITS, type = "population")
    
    ggplot(PR, aes(x = Size_nm)) +
      geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2) +
      geom_line(aes(y = mean), linewidth = 1) +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = paste0("Size trend — ", p, " (population mean)"),
        x = "Size (nm, log scale)",
        y = paste0(p, " (natural scale)")
      ) +
      theme_bw(base_size = 13)
  })
}

## ============================================================
## 12) Example outputs and saves
## ============================================================
cat("\n--- Meta coefficients & heterogeneity ---\n")
print(META_FITS %>% select(param, coef0, coef1, tau2))

cat("\n--- Slope CIs & p-values ---\n")
print(SLP)

sizes_example <- c(20, 50, 100, 200, 500, 1000, 6000)

pred_pop <- predict_params_by_size(SIZE_PARAMS, sizes_nm = sizes_example, type = "population")
pred_new <- predict_params_by_size(SIZE_PARAMS, sizes_nm = sizes_example, type = "new_id")

write_csv(pred_pop, file.path(meta_dir, "meta_pred_population_metafor.csv"))
write_csv(pred_new, file.path(meta_dir, "meta_pred_newid_metafor.csv"))

cat("\n--- Example predictions: population ---\n")
print(head(pred_pop))

cat("\n--- Example predictions: new ID ---\n")
print(head(pred_new))

pv_100 <- build_param_vector_for_size(100, type = "population", statistic = "mean")
cat("\n--- PBPK parameter vector (100 nm, population mean) ---\n")
print(pv_100)

write_csv(
  tibble(Parameter = names(pv_100), Value = as.numeric(pv_100)),
  file.path(meta_dir, "PBPK_vector_100nm_population_mean.csv")
)

plots <- plot_size_trends(SIZE_PARAMS)

for (i in seq_along(plots)) {
  p <- plots[[i]]
  param_name <- SIZE_PARAMS[i]
  
  ggsave(
    filename = file.path(figures_dir, paste0("meta_size_trend_", param_name, ".tiff")),
    plot = p,
    device = "tiff",
    width = 6,
    height = 5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

## ============================================================
## 13) Save workspace summary objects
## ============================================================
saveRDS(
  list(
    sizes_tbl   = sizes_tbl,
    LOG_MOM     = LOG_MOM,
    META_FITS   = META_FITS,
    SLP         = SLP,
    GLOBAL_SUM  = GLOBAL_SUM,
    pred_pop    = pred_pop,
    pred_new    = pred_new,
    pv_100      = pv_100,
    SIZE_PARAMS = SIZE_PARAMS,
    GLOBAL_PARAMS = GLOBAL_PARAMS
  ),
  file = file.path(meta_dir, "meta_analysis_outputs.rds")
)

cat("\n============================================================\n")
cat("Meta-analysis completed.\n")
cat("Saved outputs to:\n")
cat(" - ", normalizePath(meta_dir, winslash = "/", mustWork = FALSE), "\n")
cat(" - ", normalizePath(figures_dir, winslash = "/", mustWork = FALSE), "\n")
cat("============================================================\n")