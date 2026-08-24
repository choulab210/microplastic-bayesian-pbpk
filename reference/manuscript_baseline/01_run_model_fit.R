## ============================================================
## 01_run_model_fit.R
## Main analysis script: model fitting + GOF + NSC outputs
## ============================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(FME)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(numDeriv)
  library(tools)
  library(tibble)
})

options(stringsAsFactors = FALSE)
set.seed(123)

VERBOSE <- FALSE
dbg <- function(...) if (isTRUE(VERBOSE)) cat(...)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

## ============================================================
## 0) USER SETTINGS
## ============================================================
NSC_FRAC_CHANGE <- 0.01     # 1% perturbation
NSC_DT_H        <- 0.05     # simulation time step for AUC calculation
NSC_MIN_DENOM   <- 1e-10    # numerical floor to avoid divide-by-zero

## ============================================================
## 1) Load & compile model
## ============================================================
source("MP_PBPK.R")
mod <- mcode("mPBPK_RES", mPBPK_RES)
MODEL_PARAM_NAMES <- names(param(mod))

cat("Model parameters (N=", length(MODEL_PARAM_NAMES), ")\n", sep = "")
print(MODEL_PARAM_NAMES)

## ============================================================
## 2) Mass balance check
## ============================================================
ev1 <- ev(amt = 1, cmt = "A_lumen", time = 0)
ev2 <- ev(amt = 1, cmt = "A_dose", time = 0)
ex <- ev1 + ev2

sim <- mrgsim_d(mod, data = ex, end = 120, delta = 0.1, atol = 1e-16, rtol = 1e-12)
sim_df <- as.data.frame(sim)

png("figures/mass_balance_check.png", width = 1800, height = 1200, res = 200)
sim_out <- sim_df %>% filter(time > 0)
plot(sim_out$time, sim_out$mass_balance, type = "l", lwd = 2,
     ylab = "Mass balance (mg)",
     xlab = "Time (h)",
     main = "Mass balance check: A_dose - (A_body + A_excreted)")
abline(h = 0, col = "gray50", lty = 2)
dev.off()

## ============================================================
## 3) Load & normalize merged data
## ============================================================
obs_df <- read.csv("data.csv") %>%
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
  ) %>%
  filter(ID %in% 1:5)

stopifnot(all(c("ID","dose_mg","Time_h","Tissue","Conc") %in% names(obs_df)))

ids <- sort(unique(obs_df$ID))
cat("IDs found:", paste(ids, collapse = ", "), "\n")

## ============================================================
## 4) Prediction function
## ============================================================
pred <- function(pars = NULL,
                 id = 1,
                 logscale = TRUE,
                 t_end_h = NULL,
                 dt_h = 0.01,
                 atol = 1e-8,
                 rtol = 1e-6) {
  
  pars_nat <- if (is.null(pars)) numeric(0) else if (logscale) exp(pars) else pars
  pars_nat <- pars_nat[names(pars_nat) %in% MODEL_PARAM_NAMES]
  
  mod_tmp <- mod
  if (length(pars_nat) > 0) mod_tmp <- param(mod_tmp, pars_nat)
  
  sub_df <- obs_df %>% filter(ID == id)
  if (nrow(sub_df) == 0) stop("No observations for ID=", id)
  
  dose_mg <- unique(sub_df$dose_mg)
  if (length(dose_mg) != 1) dose_mg <- dose_mg[1]
  if (is.null(t_end_h)) t_end_h <- max(sub_df$Time_h)
  
  ex <- ev(amt = dose_mg, cmt = "A_lumen", time = 0)
  add_times <- sort(unique(as.numeric(sub_df$Time_h)))
  ts <- tgrid(start = 0, end = t_end_h, delta = dt_h, add = add_times)
  
  out <- mrgsim_d(mod_tmp, data = ex, tgrid = ts,
                  atol = atol, rtol = rtol, maxsteps = 1e6) %>%
    as.data.frame()
  
  dose_total <- dose_mg
  
  lungs_total <- if ("total_lung" %in% names(out)) {
    out$total_lung
  } else if (all(c("A_lung_b","A_lung_t") %in% names(out))) {
    out$A_lung_b + out$A_lung_t
  } else 0
  
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

## ============================================================
## 5) Cost/residual functions
## ============================================================
.allowed_tissues <- c("Blood","Liver","Kidney","Spleen","Lungs","GI","Urine")

.to_num_df <- function(df) {
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(z) {
        z <- unlist(z)
        if (length(z) == 0) return(NA_real_)
        as.numeric(z[1])
      }, numeric(1))
    } else {
      df[[nm]] <- as.numeric(df[[nm]])
    }
  }
  as.data.frame(df)
}

MCcost_id <- function(pars_log, id) {
  sub_df <- obs_df %>% filter(ID == id)
  if (nrow(sub_df) == 0) return(list(ssr = NA_real_, residuals = NULL))
  
  sim <- tryCatch(
    pred(pars = pars_log, id = id, logscale = TRUE,
         t_end_h = max(sub_df$Time_h), dt_h = 0.01,
         atol = 1e-12, rtol = 1e-10),
    error = function(e) NULL
  )
  if (is.null(sim)) return(list(ssr = NA_real_, residuals = NULL))
  
  obs_wide <- sub_df %>%
    transmute(time = as.numeric(Time_h),
              Tissue = as.character(Tissue),
              obs    = as.numeric(Conc)) %>%
    group_by(time, Tissue) %>%
    summarise(obs = mean(obs, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Tissue, values_from = obs, values_fill = NA_real_) %>%
    arrange(time)
  
  sim      <- as.data.frame(sim)
  obs_wide <- as.data.frame(obs_wide)
  
  vars <- intersect(.allowed_tissues, intersect(names(sim), names(obs_wide)))
  if (length(vars) == 0) return(list(ssr = NA_real_, residuals = NULL))
  
  sim_small <- sim[, c("time", vars), drop = FALSE] %>%
    arrange(time) %>%
    distinct(time, .keep_all = TRUE) %>%
    .to_num_df()
  
  obs_small <- obs_wide[, c("time", vars), drop = FALSE] %>%
    arrange(time) %>%
    distinct(time, .keep_all = TRUE) %>%
    .to_num_df()
  
  vars <- vars[sapply(vars, function(v) any(is.finite(obs_small[[v]])))]
  if (length(vars) == 0) return(list(ssr = NA_real_, residuals = NULL))
  
  sim_small <- sim_small[, c("time", vars), drop = FALSE]
  obs_small <- obs_small[, c("time", vars), drop = FALSE]
  
  sim_small[, vars] <- lapply(sim_small[, vars, drop = FALSE], function(x) log1p(pmax(x, 0)))
  obs_small[, vars] <- lapply(obs_small[, vars, drop = FALSE], function(x) log1p(pmax(x, 0)))
  
  cost_i <- tryCatch(
    suppressWarnings(FME::modCost(model = sim_small, obs = obs_small, x = "time")),
    error = function(e) NULL
  )
  if (is.null(cost_i)) return(list(ssr = NA_real_, residuals = NULL))
  
  ssr_i <- if (!is.null(cost_i$ssr)) {
    cost_i$ssr
  } else if (!is.null(cost_i$residuals)) {
    sum(cost_i$residuals$res^2, na.rm = TRUE)
  } else 0
  
  list(ssr = ssr_i, residuals = as.data.frame(cost_i$residuals))
}

MCcost_scalar_id <- function(pars_log, id) {
  out <- tryCatch(MCcost_id(pars_log, id), error = function(e) NULL)
  if (is.null(out) || !is.finite(out$ssr) || out$ssr <= 0) return(c(logSSR = 1))
  c(logSSR = log10(out$ssr))
}

MCresid_id <- function(pars_log, id) {
  out <- MCcost_id(pars_log, id)
  if (is.null(out$residuals) || nrow(out$residuals) == 0) return(numeric(0))
  as.numeric(out$residuals$res)
}

## ============================================================
## 6) Initial parameter vector
## ============================================================
theta.int <- log(c(
  KLRESmax = 0.4,  KLRESrelease = 7.5e-3,
  KSRESmax = 1,    KSRESrelease = 3e-3,
  KKRESmax = 0.1,  KKRESrelease = 1e-2,
  KLURESmax = 0.1, KLURESrelease = 5e-3,
  KRRESmax = 0.2,  KRRESrelease = 2e-3,
  PAL = 1e-3, PAK = 1e-3, PAS = 1e-3, PALU = 1e-3, PAGI = 1e-3, PAR = 1e-6,
  KpL = 0.01, KpK = 0.15, KpS = 0.15, KpR = 0.15, KpLU = 0.15,
  KbileC = 2.2e-3, KurineC = 3.2e-4, Kfeces = 0.53,
  Ktrans = 5e-1, Kgw2b = 3e-4, Kshed = 5e-1
))

theta.int <- theta.int[names(theta.int) %in% MODEL_PARAM_NAMES]
stopifnot(length(theta.int) > 0)

## ============================================================
## 7) Sensitivity / consensus parameter selection
## ============================================================
sens_for_id <- function(theta_log, id, eps = 1e-4) {
  g <- tryCatch(
    numDeriv::grad(
      func = function(p) MCcost_scalar_id(p, id = id)[["logSSR"]],
      x = theta_log,
      method.args = list(eps = eps)
    ),
    error = function(e) rep(NA_real_, length(theta_log))
  )
  
  tibble(
    ID = id,
    Parameter = names(theta_log),
    Sensitivity = as.numeric(g),
    AbsSensitivity = abs(as.numeric(g))
  ) %>% arrange(desc(AbsSensitivity))
}

compute_consensus_params <- function(theta_log, ids, topK_each = 12, norm_thresh = 0.5) {
  sens_list <- lapply(ids, function(id) sens_for_id(theta_log, id))
  S_all <- bind_rows(sens_list)
  
  S_all <- S_all %>%
    group_by(ID) %>%
    mutate(
      Rank = rank(-AbsSensitivity, ties.method = "min"),
      Norm = ifelse(max(AbsSensitivity, na.rm = TRUE) > 0,
                    AbsSensitivity / max(AbsSensitivity, na.rm = TRUE), 0)
    ) %>%
    ungroup()
  
  split_list <- split(S_all, S_all$ID)
  topK_sets <- lapply(split_list, function(df) {
    df <- df[order(-df$AbsSensitivity), , drop = FALSE]
    head(df$Parameter, min(topK_each, nrow(df)))
  })
  common_topK <- Reduce(intersect, topK_sets)
  
  S_aggr <- S_all %>%
    group_by(Parameter) %>%
    summarise(
      MeanAbs = mean(AbsSensitivity, na.rm = TRUE),
      GeoNorm = exp(mean(log(pmax(Norm, 1e-12)), na.rm = TRUE)),
      MinNorm = min(Norm, na.rm = TRUE),
      n_ids = n_distinct(ID),
      .groups = "drop"
    ) %>%
    arrange(desc(GeoNorm), desc(MeanAbs))
  
  if (length(common_topK) == 0) {
    common_topK <- unique(unlist(topK_sets, use.names = FALSE))
  }
  
  list(
    SensitivityTable = S_all,
    AggregateTable   = S_aggr,
    Common_TopK      = common_topK
  )
}

cons <- compute_consensus_params(theta.int, ids, topK_each = 15, norm_thresh = 0.5)
common_params <- cons$Common_TopK

cat("\nUsing COMMON parameter set:\n")
print(common_params)

## ============================================================
## 8) Fit each ID
## ============================================================
fit_one_id <- function(id, theta_log_init, parm_names, widen = 1e4, trace = 0) {
  flag_params <- c("Ktrans", "Kshed", "Kgw2b")
  
  tight_bounds <- data.frame(
    lower = c(1e-3, 0.05, 1e-4),
    upper = c(0.9,  3,    1e-2),
    row.names = flag_params
  )
  
  theta_sel <- theta_log_init[parm_names]
  theta_nat <- exp(theta_sel)
  
  lower_nat <- theta_nat / widen
  upper_nat <- theta_nat * widen
  
  for (nm in parm_names) {
    if (nm %in% flag_params) {
      lower_nat[nm] <- tight_bounds[nm, "lower"]
      upper_nat[nm] <- tight_bounds[nm, "upper"]
    }
  }
  
  lower_log <- log(lower_nat)
  upper_log <- log(upper_nat)
  
  theta_sel <- pmin(pmax(theta_sel, lower_log + 1e-9), upper_log - 1e-9)
  
  FME::modFit(
    f = function(p) MCresid_id(p, id = id),
    p = theta_sel,
    lower = lower_log,
    upper = upper_log,
    method = "Port",
    control = list(eval.max = 250, iter.max = 250, trace = trace)
  )
}

Fit_list <- lapply(ids, function(id) fit_one_id(id, theta.int, common_params, trace = 0))

## Build full per-ID parameter vectors on log scale
theta_fit_by_id <- setNames(vector("list", length(ids)), as.character(ids))
for (i in seq_along(ids)) {
  id <- ids[i]
  th <- theta.int
  fit_pars <- Fit_list[[i]]$par
  ov <- intersect(names(th), names(fit_pars))
  th[ov] <- fit_pars[ov]
  theta_fit_by_id[[as.character(id)]] <- th
}

## ============================================================
## 9) GOF data for plotting
## ============================================================
collect_gof_data <- function(ids, Fit_list) {
  bind_rows(lapply(seq_along(ids), function(i) {
    id <- ids[i]
    
    sim_fit <- pred(pars = Fit_list[[i]]$par, id = id, logscale = TRUE) %>%
      pivot_longer(-time, names_to = "Tissue", values_to = "pred")
    
    obs_plot <- obs_df %>%
      filter(ID == id) %>%
      transmute(
        ID = ID,
        time = as.numeric(Time_h),
        Tissue = as.character(Tissue),
        obs = as.numeric(Conc),
        SD = as.numeric(SD)
      )
    
    sim_fit %>%
      inner_join(obs_plot, by = c("time", "Tissue")) %>%
      mutate(ID = id)
  }))
}

calc_gof_metrics <- function(df) {
  df <- df %>% filter(is.finite(obs), is.finite(pred), obs > 0, pred > 0)
  
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

gof_all <- collect_gof_data(ids, Fit_list)

overall_metrics <- calc_gof_metrics(gof_all)

per_id_metrics <- bind_rows(lapply(ids, function(id) {
  cbind(data.frame(ID = id), calc_gof_metrics(gof_all %>% filter(ID == id)))
}))

per_tissue_metrics <- bind_rows(lapply(unique(gof_all$Tissue), function(tt) {
  cbind(data.frame(Tissue = tt), calc_gof_metrics(gof_all %>% filter(Tissue == tt)))
}))

## ============================================================
## 10) NSC workflow (AUC-based)
## ============================================================
trapz_auc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

calc_auc_by_tissue <- function(pred_df) {
  pred_df %>%
    pivot_longer(-time, names_to = "Tissue", values_to = "Value") %>%
    group_by(Tissue) %>%
    summarise(
      AUC = trapz_auc(time, Value),
      .groups = "drop"
    )
}

calc_nsc_for_id <- function(id,
                            pars_log_base,
                            param_names = names(pars_log_base),
                            frac_change = NSC_FRAC_CHANGE,
                            dt_h = NSC_DT_H) {
  
  t_end <- max(obs_df$Time_h[obs_df$ID == id], na.rm = TRUE)
  
  base_pred <- pred(
    pars = pars_log_base,
    id = id,
    logscale = TRUE,
    t_end_h = t_end,
    dt_h = dt_h
  )
  
  base_auc <- calc_auc_by_tissue(base_pred) %>%
    mutate(ID = id)
  
  nsc_list <- lapply(param_names, function(pn) {
    pars_pert <- pars_log_base
    
    # log-scale perturbation corresponding to multiplicative (1 + frac_change)
    pars_pert[pn] <- pars_log_base[pn] + log(1 + frac_change)
    
    pert_pred <- tryCatch(
      pred(
        pars = pars_pert,
        id = id,
        logscale = TRUE,
        t_end_h = t_end,
        dt_h = dt_h
      ),
      error = function(e) NULL
    )
    
    if (is.null(pert_pred)) return(NULL)
    
    pert_auc <- calc_auc_by_tissue(pert_pred)
    
    out <- base_auc %>%
      select(Tissue, AUC_base = AUC) %>%
      left_join(pert_auc %>% rename(AUC_pert = AUC), by = "Tissue") %>%
      mutate(
        ID = id,
        Parameter = pn,
        frac_param = frac_change,
        frac_auc = (AUC_pert - AUC_base) / pmax(abs(AUC_base), NSC_MIN_DENOM),
        NSC = frac_auc / frac_change,
        absNSC = abs(NSC)
      ) %>%
      select(ID, Parameter, Tissue, AUC_base, AUC_pert, frac_param, frac_auc, NSC, absNSC)
    
    out
  })
  
  NSC_id <- bind_rows(nsc_list)
  
  list(
    base_auc = base_auc,
    nsc = NSC_id
  )
}

cat("\n============================================================\n")
cat("Running AUC-based NSC analysis\n")
cat("============================================================\n")

NSC_results_by_id <- lapply(ids, function(id) {
  cat("Running NSC for ID =", id, "\n")
  calc_nsc_for_id(
    id = id,
    pars_log_base = theta_fit_by_id[[as.character(id)]],
    param_names = names(theta_fit_by_id[[as.character(id)]])
  )
})
names(NSC_results_by_id) <- as.character(ids)

base_auc_all <- bind_rows(lapply(NSC_results_by_id, `[[`, "base_auc"))
NSC_amount_all <- bind_rows(lapply(NSC_results_by_id, `[[`, "nsc"))

cat("\nNSC rows generated:", nrow(NSC_amount_all), "\n")
cat("Base AUC rows generated:", nrow(base_auc_all), "\n")

## ============================================================
## 11) Save analysis outputs
## ============================================================
write.csv(gof_all, "results/gof_all.csv", row.names = FALSE)
write.csv(per_id_metrics, "results/per_id_metrics.csv", row.names = FALSE)
write.csv(per_tissue_metrics, "results/per_tissue_metrics.csv", row.names = FALSE)

write.csv(base_auc_all, "results/base_auc_all.csv", row.names = FALSE)
write.csv(NSC_amount_all, "results/NSC_amount_all.csv", row.names = FALSE)

save_list <- list(
  ids = ids,
  Fit_list = Fit_list,
  gof_all = gof_all,
  per_id_metrics = per_id_metrics,
  per_tissue_metrics = per_tissue_metrics,
  obs_df = obs_df,
  common_params = common_params,
  theta.int = theta.int,
  theta_fit_by_id = theta_fit_by_id,
  base_auc_all = base_auc_all,
  NSC_amount_all = NSC_amount_all
)

saveRDS(save_list, file = "results/fitted_objects.rds")

cat("\nSaved:\n")
cat(" - results/gof_all.csv\n")
cat(" - results/per_id_metrics.csv\n")
cat(" - results/per_tissue_metrics.csv\n")
cat(" - results/base_auc_all.csv\n")
cat(" - results/NSC_amount_all.csv\n")
cat(" - results/fitted_objects.rds\n")

cat("\nAnalysis completed and saved to results/ folder.\n")