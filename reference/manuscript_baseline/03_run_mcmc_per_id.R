## ============================================================
## 03_run_mcmc_per_id.R
## Per-ID DRAM MCMC with Priors, var0 Calibration, Pilot Tuning
## Posterior summaries + 95% CI PPC vs observations
## Includes convergence plots:
## - trace plots
## - density plots
## - running mean plots
## - autocorrelation plots
## ============================================================

suppressPackageStartupMessages({
  library(FME)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(coda)
  library(purrr)
  library(tibble)
})

set.seed(1234)
options(stringsAsFactors = FALSE)

## ------------------------------------------------------------
## 0) USER OPTIONS
## ------------------------------------------------------------
# Prior center strategy:
#  - "per_id": center every parameter prior at that ID's deterministic fit
#  - "global": center every prior at global baseline (same for all IDs)
#  - "mixed" : per-ID centers only for key levers; global centers for others
PRIOR_CENTER_STRATEGY <- "mixed"

# Strength multiplier for all priors (1.0 = as specified, <1 weaker, >1 stronger)
PRIOR_WEIGHT <- 0.75

# MCMC lengths
PILOT_ITER <- 1000
MAIN_ITER  <- 3000

# Posterior predictive sampling
NSAMP_PPC  <- 400
DT_PPC     <- 0.05    # hours
Y_FLOOR    <- 1e-6
LOG_SCALE  <- TRUE

# Convergence plot settings
TRACE_NCOL <- 4
DENSITY_NCOL <- 4
RUNMEAN_NCOL <- 4
ACF_MAX_LAG <- 50

# Optional: plot only a subset of parameters in convergence plots
# Set to NULL to plot all sampled parameters
KEY_PLOT_PARAMS <- NULL
# Example:
# KEY_PLOT_PARAMS <- c("Ktrans","Kshed","Kgw2b","KLRESmax","Kfeces","KurineC")

# folders
results_dir <- "results"
figures_dir <- "figures"
mcmc_dir    <- file.path(results_dir, "mcmc")
conv_dir    <- file.path(figures_dir, "mcmc_convergence")

dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)
dir.create(mcmc_dir, showWarnings = FALSE)
dir.create(conv_dir, showWarnings = FALSE)

## ------------------------------------------------------------
## 1) Load prerequisites
## ------------------------------------------------------------
# This script assumes the deterministic fit script has already been run.
# It expects objects/functions such as:
#   MODEL_PARAM_NAMES, ids, Fit_list, obs_df, pred(), MCresid_id()
#
# Recommended usage:
#   source("01_run_model_fit.R")
#   source("03_run_mcmc_per_id.R")

required_objects <- c("MODEL_PARAM_NAMES", "ids", "Fit_list", "obs_df")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ", paste(missing_objects, collapse = ", "),
    "\nPlease run/source 01_run_model_fit.R first so these objects/functions are available."
  )
}

if (!exists("pred", mode = "function")) {
  stop("Function 'pred()' not found. Please run/source 01_run_model_fit.R first.")
}

if (!exists("MCresid_id", mode = "function")) {
  stop("Function 'MCresid_id()' not found. Please run/source 01_run_model_fit.R first.")
}

ids <- as.numeric(ids)

## ------------------------------------------------------------
## 2) Baseline (log) and per-ID inits from deterministic fits
## ------------------------------------------------------------
theta_base <- log(c(
  # RES capture/release
  KLRESmax = 0.4,  KLRESrelease = 7.5e-3,
  KSRESmax = 1,    KSRESrelease = 3e-3,
  KKRESmax = 0.1,  KKRESrelease = 1e-2,
  KLURESmax = 0.1, KLURESrelease = 5e-3,
  KRRESmax = 0.2,  KRRESrelease = 2e-3,
  
  # Permeability
  PAL = 1e-3, PAK = 1e-3, PAS = 1e-3, PALU = 1e-3, PAGI = 1e-3, PAR = 1e-6,
  
  # Excretion
  KbileC = 2.2e-3, KurineC = 3.2e-4, Kfeces = 0.53,
  
  # Gut-wall
  Ktrans = 5e-1, Kgw2b = 3e-4, Kshed = 5e-1
))
theta_base <- theta_base[names(theta_base) %in% MODEL_PARAM_NAMES]

if (length(theta_base) == 0) {
  stop("theta_base has length 0 after matching MODEL_PARAM_NAMES.")
}

# Full parameter set sampled for every ID
mcmc_params_all <- names(theta_base)

# Per-ID initial vectors: baseline overwritten by that ID's deterministic (log) fit
theta_by_id <- setNames(vector("list", length(ids)), as.character(ids))
for (i in seq_along(ids)) {
  id <- ids[i]
  th <- theta_base
  if (!is.null(Fit_list[[i]]$par)) {
    fit_pars <- Fit_list[[i]]$par   # assumed log-scale from modFit
    ov <- intersect(names(th), names(fit_pars))
    th[ov] <- fit_pars[ov]
  }
  theta_by_id[[as.character(id)]] <- th
}

## ------------------------------------------------------------
## 3) Priors (log-scale): SD map and center builders
## ------------------------------------------------------------
PRIOR_SD_MAP <- c(
  # Gut-wall (moderately informative)
  Ktrans = 0.6, Kshed = 0.6, Kgw2b = 0.8,
  
  # RES capture/release (weak-to-moderate)
  KLRESmax = 0.7, KSRESmax = 0.7, KKRESmax = 0.8, KLURESmax = 0.8, KRRESmax = 0.8,
  KLRESrelease = 0.5, KSRESrelease = 0.5, KKRESrelease = 0.6, KLURESrelease = 0.6, KRRESrelease = 0.6,
  
  # Permeabilities (often weakly informed)
  PAL = 1.0, PAK = 1.0, PAS = 1.0, PALU = 1.0, PAGI = 1.0, PAR = 1.0,
  
  # Excretion
  KbileC = 0.6, Kfeces = 0.6, KurineC = 0.9
)
PRIOR_SD_DEFAULT <- 0.9

# Parameters that get per-ID centers under "mixed"
KEY_PER_ID <- c("Ktrans","Kshed","Kgw2b","KLRESmax","Kfeces","KurineC")

build_prior_for_id <- function(id, theta_by_id, theta_base, param_names,
                               strategy = c("mixed","per_id","global")) {
  strategy <- match.arg(strategy)
  
  if (strategy == "per_id") {
    mu <- theta_by_id[[as.character(id)]][param_names]
  } else if (strategy == "global") {
    mu <- theta_base[param_names]
  } else { # "mixed"
    mu <- theta_base[param_names]
    mu[intersect(param_names, KEY_PER_ID)] <-
      theta_by_id[[as.character(id)]][intersect(param_names, KEY_PER_ID)]
  }
  
  mu <- mu[param_names]
  sd <- PRIOR_SD_MAP[param_names]
  sd[is.na(sd)] <- PRIOR_SD_DEFAULT
  
  list(mu = mu, sd = sd)
}

## ------------------------------------------------------------
## 4) Utilities: acceptance, summaries, PPC & plots
## ------------------------------------------------------------
acc_rate_from_pars <- function(m, tol = 1e-10) {
  P <- as.matrix(m$pars)
  if (is.null(P) || nrow(P) < 2) return(NA_real_)
  moves <- c(FALSE, rowSums(abs(diff(P)) > tol) > 0)
  mean(moves, na.rm = TRUE)
}

param_summary_for_id <- function(id, mcmc3) {
  P <- do.call(rbind, lapply(mcmc3$chains, function(m) as.data.frame(m$pars)))
  tibble::as_tibble(P) %>%
    summarise(across(
      everything(),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE),
        q025 = ~quantile(.x, 0.025, na.rm = TRUE),
        q050 = ~quantile(.x, 0.5,   na.rm = TRUE),
        q975 = ~quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )) %>%
    pivot_longer(
      everything(),
      names_to = c("Parameter", ".value"),
      names_sep = "_"
    ) %>%
    arrange(Parameter) %>%
    mutate(ID = as.numeric(id), .before = 1)
}

posterior_summary_for_id <- function(id, mcmc3, nsamp = NSAMP_PPC, ppc_dt = DT_PPC) {
  P <- do.call(rbind, lapply(mcmc3$chains, function(m) as.data.frame(m$pars)))
  stopifnot(nrow(P) > 0)
  
  nsamp <- min(nsamp, nrow(P))
  set.seed(1000 + as.integer(id))
  idx <- sample(seq_len(nrow(P)), nsamp, replace = FALSE)
  
  t_end <- max(obs_df$Time_h[obs_df$ID == as.numeric(id)], na.rm = TRUE)
  
  preds <- vector("list", length(idx))
  for (k in seq_along(idx)) {
    pars_log <- unlist(P[idx[k], , drop = TRUE])
    
    simk <- pred(
      pars = pars_log,
      id = as.numeric(id),
      logscale = TRUE,
      t_end_h = t_end,
      dt_h = ppc_dt,
      atol = 1e-8,
      rtol = 1e-6
    )
    
    preds[[k]] <- simk %>%
      tidyr::pivot_longer(-time, names_to = "Tissue", values_to = "pred") %>%
      mutate(draw = k)
  }
  
  pred_all <- bind_rows(preds)
  
  pred_all %>%
    group_by(Tissue, time) %>%
    summarise(
      pred_med = median(pred, na.rm = TRUE),
      pred_lo  = quantile(pred, 0.025, na.rm = TRUE),
      pred_hi  = quantile(pred, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(ID = as.numeric(id), .before = 1)
}

plot_id_posterior_ci <- function(id, mcmc3, y_floor = Y_FLOOR, log_scale = LOG_SCALE) {
  post <- posterior_summary_for_id(id, mcmc3, nsamp = NSAMP_PPC, ppc_dt = DT_PPC)
  
  obs_i <- obs_df %>%
    filter(ID == as.numeric(id)) %>%
    transmute(
      time = as.numeric(Time_h),
      Tissue = as.character(Tissue),
      obs = as.numeric(Conc),
      SD  = as.numeric(SD)
    )
  
  post <- post %>%
    mutate(
      pred_lo  = pmax(pred_lo,  y_floor),
      pred_med = pmax(pred_med, y_floor),
      pred_hi  = pmax(pred_hi,  y_floor)
    )
  
  p <- ggplot(post, aes(x = time)) +
    geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.20) +
    geom_line(aes(y = pred_med), linewidth = 1) +
    geom_point(
      data = obs_i,
      aes(x = time, y = obs),
      inherit.aes = FALSE,
      size = 2
    ) +
    geom_errorbar(
      data = obs_i,
      aes(
        x = time,
        ymin = pmax(obs - SD, y_floor),
        ymax = obs + SD
      ),
      inherit.aes = FALSE,
      width = 0.4,
      alpha = 0.6
    ) +
    facet_wrap(~ Tissue, scales = "free_y") +
    labs(
      title = paste("Posterior predictive (95% CI) vs Observations — ID", id),
      x = "Time (h)",
      y = "%ID"
    ) +
    theme_bw(base_size = 13)
  
  if (log_scale) p <- p + scale_y_log10()
  
  p
}

## ------------------------------------------------------------
## 4b) Convergence plotting helpers
## ------------------------------------------------------------
make_trace_df <- function(mcmc3, id, keep_params = KEY_PLOT_PARAMS) {
  df <- bind_rows(lapply(seq_along(mcmc3$chains), function(ch) {
    P <- as.data.frame(mcmc3$chains[[ch]]$pars)
    P$Iter <- seq_len(nrow(P))
    P$Chain <- paste0("Chain", ch)
    P$ID <- as.numeric(id)
    P
  })) %>%
    pivot_longer(
      cols = -c(Iter, Chain, ID),
      names_to = "Parameter",
      values_to = "Value"
    )
  
  if (!is.null(keep_params)) {
    df <- df %>% filter(Parameter %in% keep_params)
  }
  
  df
}

plot_mcmc_trace <- function(id, mcmc3, ncol = TRACE_NCOL, keep_params = KEY_PLOT_PARAMS) {
  df <- make_trace_df(mcmc3, id, keep_params = keep_params)
  
  ggplot(df, aes(x = Iter, y = Value, color = Chain)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_wrap(~ Parameter, scales = "free_y", ncol = ncol) +
    labs(
      title = paste("MCMC Trace Plot — ID", id),
      x = "Iteration",
      y = "Log-parameter value",
      color = "Chain"
    ) +
    theme_bw(base_size = 12)
}

plot_mcmc_density <- function(id, mcmc3, ncol = DENSITY_NCOL, keep_params = KEY_PLOT_PARAMS) {
  df <- make_trace_df(mcmc3, id, keep_params = keep_params)
  
  ggplot(df, aes(x = Value, color = Chain, fill = Chain)) +
    geom_density(alpha = 0.15, linewidth = 0.6) +
    facet_wrap(~ Parameter, scales = "free", ncol = ncol) +
    labs(
      title = paste("Posterior Density by Chain — ID", id),
      x = "Log-parameter value",
      y = "Density",
      color = "Chain",
      fill = "Chain"
    ) +
    theme_bw(base_size = 12)
}

plot_mcmc_running_mean <- function(id, mcmc3, ncol = RUNMEAN_NCOL, keep_params = KEY_PLOT_PARAMS) {
  df <- make_trace_df(mcmc3, id, keep_params = keep_params) %>%
    group_by(Chain, Parameter) %>%
    arrange(Iter, .by_group = TRUE) %>%
    mutate(RunningMean = cumsum(Value) / seq_along(Value)) %>%
    ungroup()
  
  ggplot(df, aes(x = Iter, y = RunningMean, color = Chain)) +
    geom_line(alpha = 0.8, linewidth = 0.4) +
    facet_wrap(~ Parameter, scales = "free_y", ncol = ncol) +
    labs(
      title = paste("Running Mean Plot — ID", id),
      x = "Iteration",
      y = "Running mean",
      color = "Chain"
    ) +
    theme_bw(base_size = 12)
}

plot_mcmc_acf <- function(id, mcmc3, max_lag = ACF_MAX_LAG, keep_params = KEY_PLOT_PARAMS) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  coda_list <- mcmc3$coda
  
  if (!is.null(keep_params)) {
    coda_list <- coda::mcmc.list(lapply(coda_list, function(mc) {
      mat <- as.matrix(mc)
      keep <- intersect(colnames(mat), keep_params)
      coda::mcmc(mat[, keep, drop = FALSE])
    }))
  }
  
  coda::acfplot(coda_list, lag.max = max_lag,
                main = paste("Autocorrelation Plot — ID", id))
}

## ------------------------------------------------------------
## 5) Residual builder + var0 calibration
## ------------------------------------------------------------
build_resid_fun <- function(id,
                            prior_mu, prior_sd,
                            prior_weight = 1.0) {
  force(id); force(prior_mu); force(prior_sd); force(prior_weight)
  
  function(p) {
    data_res <- tryCatch(MCresid_id(p, id), error = function(e) rep(NA_real_, 1))
    if (any(!is.finite(data_res))) data_res <- rep(0, length.out = 1)
    
    pp <- p[names(prior_mu)]
    pr <- ((pp - prior_mu) / prior_sd) * prior_weight
    
    c(data_res, pr)
  }
}

calibrate_var0 <- function(resid_fun,
                           p_start,
                           n_probe = 5,
                           jitter_sd = 0.05) {
  probes <- vector("list", n_probe + 1L)
  probes[[1]] <- p_start
  
  if (n_probe > 0) {
    for (k in 1:n_probe) {
      probes[[k + 1]] <- p_start + rnorm(length(p_start), 0, jitter_sd)
    }
  }
  
  ss <- sapply(probes, function(p) {
    r <- resid_fun(p)
    sum(r * r, na.rm = TRUE) / length(r)
  })
  
  max(1.0, median(ss, na.rm = TRUE))
}

## ------------------------------------------------------------
## 6) Single-chain DRAM (no adaptive covariance) WITH PRIORS
## ------------------------------------------------------------
run_mcmc_one_id_chain <- function(id,
                                  start_pars_log,
                                  prior_mu, prior_sd,
                                  n_iter = MAIN_ITER,
                                  jump_scale = 0.24,
                                  seed = 1,
                                  var0_scalar = NULL,
                                  ntrydr = 2, drscale = 2,
                                  updatecov = 0, wvar0 = 1,
                                  prior_weight = PRIOR_WEIGHT) {
  set.seed(seed)
  
  pars_start <- start_pars_log[mcmc_params_all]
  pars_start <- pars_start[!is.na(pars_start)]
  
  resid_fun <- build_resid_fun(id, prior_mu, prior_sd, prior_weight)
  
  if (is.null(var0_scalar)) {
    var0_scalar <- calibrate_var0(resid_fun, pars_start, n_probe = 5, jitter_sd = 0.05)
    cat(sprintf("  Calibrated var0 (ID %s): %.3f\n", id, var0_scalar))
  }
  
  out <- FME::modMCMC(
    f = resid_fun,
    p = pars_start,
    niter = n_iter,
    jump  = jump_scale,
    var0  = var0_scalar,
    lower = log(exp(pars_start) / 100),
    upper = log(exp(pars_start) * 100),
    updatecov = updatecov,
    wvar0 = wvar0,
    burninlength = floor(n_iter / 5),
    ntrydr = ntrydr,
    drscale = drscale
  )
  
  if (is.null(colnames(out$pars))) colnames(out$pars) <- names(pars_start)
  out
}

## ------------------------------------------------------------
## 7) Pilot tuning of jump_scale per ID (with priors)
## ------------------------------------------------------------
pilot_tune_jump <- function(id, start_pars_log,
                            prior_mu, prior_sd,
                            pilot_iter = PILOT_ITER,
                            init_jump = 0.34,
                            target = c(0.20, 0.35)) {
  base <- start_pars_log[mcmc_params_all]
  base <- base[!is.na(base)]
  base <- base + rnorm(length(base), 0, 0.05)
  
  resid_fun <- build_resid_fun(id, prior_mu, prior_sd, PRIOR_WEIGHT)
  var0_cal  <- calibrate_var0(resid_fun, base, n_probe = 5, jitter_sd = 0.05)
  
  cat(sprintf("  [Pilot] initial var0 (ID %s): %.3f\n", id, var0_cal))
  
  jump <- init_jump
  
  for (k in 1:4) {
    m <- FME::modMCMC(
      f = resid_fun,
      p = base,
      niter = pilot_iter,
      jump  = jump,
      var0  = var0_cal,
      lower = log(exp(base) / 100),
      upper = log(exp(base) * 100),
      updatecov = 0,
      wvar0 = 1,
      burninlength = floor(pilot_iter / 5),
      ntrydr = 2,
      drscale = 2
    )
    
    P <- as.matrix(m$pars)
    acc <- if (is.null(P) || nrow(P) < 2) {
      NA_real_
    } else {
      moves <- c(FALSE, rowSums(abs(diff(P)) > 1e-10) > 0)
      mean(moves, na.rm = TRUE)
    }
    
    cat(sprintf(
      "  Pilot %d (ID %s): jump=%.3f var0=%.3f -> acc=%.3f\n",
      k, id, jump, var0_cal, acc
    ))
    
    if (is.na(acc)) acc <- 0
    
    if (acc < 0.05) {
      var0_cal <- var0_cal * 5
      next
    } else if (acc > target[2]) {
      jump <- jump * 1.6
    } else if (acc < target[1]) {
      jump <- jump * 0.7
    } else {
      break
    }
  }
  
  list(jump = jump, var0 = var0_cal, start = base)
}

## ------------------------------------------------------------
## 8) Three chains per ID
## ------------------------------------------------------------
run_mcmc_3chains_id <- function(id,
                                start_pars_log,
                                n_iter = MAIN_ITER,
                                strategy = PRIOR_CENTER_STRATEGY) {
  pri <- build_prior_for_id(
    id,
    theta_by_id,
    theta_base,
    mcmc_params_all,
    strategy = strategy
  )
  
  prior_mu <- pri$mu
  prior_sd <- pri$sd
  
  tuned <- pilot_tune_jump(
    id,
    start_pars_log,
    prior_mu,
    prior_sd,
    pilot_iter = PILOT_ITER,
    init_jump = 0.28,
    target = c(0.20, 0.35)
  )
  
  cat(sprintf("ID %s tuned jump=%.3f, var0=%.3f\n", id, tuned$jump, tuned$var0))
  
  base <- tuned$start
  seeds <- c(11, 22, 33)
  
  chains <- lapply(seq_along(seeds), function(k) {
    p0 <- base + rnorm(length(base), 0, 0.10)
    run_mcmc_one_id_chain(
      id, p0,
      prior_mu = prior_mu,
      prior_sd = prior_sd,
      n_iter = n_iter,
      jump_scale = tuned$jump,
      seed = seeds[k],
      var0_scalar = tuned$var0,
      updatecov = 0,
      wvar0 = 1,
      prior_weight = PRIOR_WEIGHT
    )
  })
  
  acc <- sapply(chains, function(m) {
    P <- as.matrix(m$pars)
    if (is.null(P) || nrow(P) < 2) {
      NA_real_
    } else {
      moves <- c(FALSE, rowSums(abs(diff(P)) > 1e-10) > 0)
      mean(moves, na.rm = TRUE)
    }
  })
  
  cat(sprintf(
    "ID %s acceptance per chain: %s\n",
    id, paste(sprintf("%.3f", acc), collapse = ", ")
  ))
  
  coda_list <- coda::mcmc.list(
    lapply(chains, function(m) coda::mcmc(as.matrix(m$pars)))
  )
  
  list(
    chains = chains,
    coda = coda_list,
    accept = acc,
    jump = tuned$jump,
    var0 = tuned$var0
  )
}

## ------------------------------------------------------------
## 9) Run all IDs
## ------------------------------------------------------------
MCMC3 <- setNames(vector("list", length(ids)), as.character(ids))

for (i in seq_along(ids)) {
  id <- ids[i]
  cat("\n=== Running DRAM MCMC x3 (priors; no adapt) for ID =", id, "===\n")
  
  MCMC3[[as.character(id)]] <- run_mcmc_3chains_id(
    id,
    theta_by_id[[as.character(id)]],
    n_iter = MAIN_ITER,
    strategy = PRIOR_CENTER_STRATEGY
  )
}

## ------------------------------------------------------------
## 10) Convergence diagnostics
## ------------------------------------------------------------
gelman_results <- list()
ess_results <- list()

for (id in names(MCMC3)) {
  cat("\n=== Convergence diagnostics: ID", id, "===\n")
  
  gd <- tryCatch(
    coda::gelman.diag(MCMC3[[id]]$coda, autoburnin = FALSE),
    error = function(e) e
  )
  
  ess <- tryCatch(
    coda::effectiveSize(MCMC3[[id]]$coda),
    error = function(e) NA
  )
  
  gelman_results[[id]] <- gd
  ess_results[[id]] <- ess
  
  print(gd)
  cat("Min ESS across parameters:", suppressWarnings(min(ess, na.rm = TRUE)), "\n")
}

## ------------------------------------------------------------
## 10b) Convergence plots
## ------------------------------------------------------------
for (id in names(MCMC3)) {
  cat("\n=== Saving convergence plots: ID", id, "===\n")
  
  p_trace <- plot_mcmc_trace(id, MCMC3[[id]], ncol = TRACE_NCOL, keep_params = KEY_PLOT_PARAMS)
  ggsave(
    filename = file.path(conv_dir, paste0("MCMC_trace_ID", id, ".tiff")),
    plot = p_trace,
    device = "tiff",
    width = 14,
    height = 10,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  p_dens <- plot_mcmc_density(id, MCMC3[[id]], ncol = DENSITY_NCOL, keep_params = KEY_PLOT_PARAMS)
  ggsave(
    filename = file.path(conv_dir, paste0("MCMC_density_ID", id, ".tiff")),
    plot = p_dens,
    device = "tiff",
    width = 14,
    height = 10,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  p_runmean <- plot_mcmc_running_mean(id, MCMC3[[id]], ncol = RUNMEAN_NCOL, keep_params = KEY_PLOT_PARAMS)
  ggsave(
    filename = file.path(conv_dir, paste0("MCMC_runningMean_ID", id, ".tiff")),
    plot = p_runmean,
    device = "tiff",
    width = 14,
    height = 10,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  tiff(
    filename = file.path(conv_dir, paste0("MCMC_acf_ID", id, ".tiff")),
    width = 14,
    height = 10,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  plot_mcmc_acf(id, MCMC3[[id]], max_lag = ACF_MAX_LAG, keep_params = KEY_PLOT_PARAMS)
  dev.off()
}

## ------------------------------------------------------------
## 11) Posterior parameter summaries
## ------------------------------------------------------------
POSTERIOR_PARAM_SUMMARY <- list()
POSTERIOR_PARAM_SUMMARY_DF <- list()
POSTERIOR_PPC_SUMMARY <- list()

for (id in names(MCMC3)) {
  cat("\n--- Posterior parameter summary (log-scale) for ID", id, "---\n")
  
  S <- param_summary_for_id(id, MCMC3[[id]])
  print(S)
  
  POSTERIOR_PARAM_SUMMARY[[id]] <- S
  POSTERIOR_PARAM_SUMMARY_DF[[id]] <- S
  
  ppc_sum <- posterior_summary_for_id(id, MCMC3[[id]], nsamp = NSAMP_PPC, ppc_dt = DT_PPC)
  POSTERIOR_PPC_SUMMARY[[id]] <- ppc_sum
  
  p_ppc <- plot_id_posterior_ci(id, MCMC3[[id]], y_floor = Y_FLOOR, log_scale = LOG_SCALE)
  
  ggsave(
    filename = file.path(figures_dir, paste0("PosteriorPredictive_ID", id, ".tiff")),
    plot = p_ppc,
    device = "tiff",
    width = 10,
    height = 7,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

posterior_param_summary_all <- bind_rows(POSTERIOR_PARAM_SUMMARY_DF)
posterior_ppc_summary_all   <- bind_rows(POSTERIOR_PPC_SUMMARY)

## ------------------------------------------------------------
## 12) Save results
## ------------------------------------------------------------
saveRDS(
  MCMC3,
  file = file.path(
    mcmc_dir,
    sprintf("MCMC3_results_fullset_perID_priors_%s.rds", PRIOR_CENTER_STRATEGY)
  )
)

saveRDS(
  POSTERIOR_PARAM_SUMMARY,
  file = file.path(
    mcmc_dir,
    sprintf("ParamSummary_fullset_perID_priors_%s.rds", PRIOR_CENTER_STRATEGY)
  )
)

saveRDS(
  POSTERIOR_PPC_SUMMARY,
  file = file.path(
    mcmc_dir,
    sprintf("PPCSummary_fullset_perID_priors_%s.rds", PRIOR_CENTER_STRATEGY)
  )
)

saveRDS(
  gelman_results,
  file = file.path(
    mcmc_dir,
    sprintf("GelmanDiag_fullset_perID_priors_%s.rds", PRIOR_CENTER_STRATEGY)
  )
)

saveRDS(
  ess_results,
  file = file.path(
    mcmc_dir,
    sprintf("ESS_fullset_perID_priors_%s.rds", PRIOR_CENTER_STRATEGY)
  )
)

write.csv(
  posterior_param_summary_all,
  file = file.path(
    mcmc_dir,
    sprintf("ParamSummary_fullset_perID_priors_%s.csv", PRIOR_CENTER_STRATEGY)
  ),
  row.names = FALSE
)

write.csv(
  posterior_ppc_summary_all,
  file = file.path(
    mcmc_dir,
    sprintf("PPCSummary_fullset_perID_priors_%s.csv", PRIOR_CENTER_STRATEGY)
  ),
  row.names = FALSE
)

# acceptance summary
acceptance_df <- bind_rows(lapply(names(MCMC3), function(id) {
  data.frame(
    ID = as.numeric(id),
    chain = seq_along(MCMC3[[id]]$accept),
    acceptance = MCMC3[[id]]$accept,
    jump = MCMC3[[id]]$jump,
    var0 = MCMC3[[id]]$var0
  )
}))

write.csv(
  acceptance_df,
  file = file.path(
    mcmc_dir,
    sprintf("MCMC_acceptance_fullset_perID_priors_%s.csv", PRIOR_CENTER_STRATEGY)
  ),
  row.names = FALSE
)

# ESS summary table
ess_df <- bind_rows(lapply(names(ess_results), function(id) {
  ess_i <- ess_results[[id]]
  if (all(is.na(ess_i))) return(NULL)
  
  data.frame(
    ID = as.numeric(id),
    Parameter = names(ess_i),
    ESS = as.numeric(ess_i)
  )
}))

if (nrow(ess_df) > 0) {
  write.csv(
    ess_df,
    file = file.path(
      mcmc_dir,
      sprintf("ESS_fullset_perID_priors_%s.csv", PRIOR_CENTER_STRATEGY)
    ),
    row.names = FALSE
  )
}

cat("\n============================================================\n")
cat("MCMC run completed.\n")
cat("Saved outputs to:\n")
cat("  ", normalizePath(mcmc_dir, winslash = "/", mustWork = FALSE), "\n")
cat("  ", normalizePath(figures_dir, winslash = "/", mustWork = FALSE), "\n")
cat("  ", normalizePath(conv_dir, winslash = "/", mustWork = FALSE), "\n")
cat("============================================================\n")