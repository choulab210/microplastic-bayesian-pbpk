## ============================================================
## 03_run_mcmc.R
##
## Bayesian PBPK calibration by development dataset
##
## Usage:
##
## Run all development IDs:
##   Rscript scripts/03_run_mcmc.R --all
##
## Run one ID:
##   Rscript scripts/03_run_mcmc.R --id 1
##
## If no argument is supplied:
##   all development IDs are run.
##
## Inputs:
##   results/deterministic/deterministic_outputs.rds
##   results/sensitivity/sensitivity_outputs.rds
##   model/MP_PBPK.R
##
## Outputs:
##   results/mcmc/ID1/
##   results/mcmc/ID2/
##   ...
##
## Algorithm:
##   Delayed-rejection Metropolis MCMC using FME::modMCMC
##
## NOTE:
##   updatecov = 0; therefore adaptive covariance is not used.
## ============================================================


## ============================================================
## 0) Packages
## ============================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(FME)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(coda)
})

options(stringsAsFactors = FALSE)


## ============================================================
## 1) Analysis settings
## ============================================================

PRIOR_CENTER_STRATEGY <- "mixed"
PRIOR_WEIGHT <- 0.75

PILOT_ITER <- 1000
MAIN_ITER  <- 10000

N_CHAINS <- 3

CHAIN_SEEDS <- c(
  11,
  22,
  33
)

BURN_FRACTION <- 0.20


## ------------------------------------------------------------
## Posterior parameter-uncertainty propagation
## ------------------------------------------------------------

NSAMP_PPC <- 400
PPC_DT_H  <- 0.05


## ------------------------------------------------------------
## Proposal settings
## ------------------------------------------------------------

INITIAL_JUMP <- 0.28

TARGET_ACCEPTANCE <- c(
  0.20,
  0.35
)

NTRY_DR  <- 2
DR_SCALE <- 2


## ------------------------------------------------------------
## Input / output paths
## ------------------------------------------------------------

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

MODEL_FILE <- file.path(
  "model",
  "MP_PBPK.R"
)

RESULTS_DIR <- file.path(
  "results",
  "mcmc"
)

dir.create(
  RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 2) Parse command-line arguments
## ============================================================

args <- commandArgs(
  trailingOnly = TRUE
)

RUN_ALL <- FALSE
RUN_ID  <- NA_integer_


if (length(args) == 0) {

  RUN_ALL <- TRUE

} else if ("--all" %in% args) {

  RUN_ALL <- TRUE

} else if ("--id" %in% args) {

  id_pos <- match(
    "--id",
    args
  )

  if (id_pos == length(args)) {
    stop("--id requires a dataset number.")
  }

  RUN_ID <- suppressWarnings(
    as.integer(
      args[id_pos + 1]
    )
  )

  if (!is.finite(RUN_ID)) {
    stop("Invalid dataset ID supplied after --id.")
  }

} else {

  stop(
    "Unknown command-line arguments.\n\n",
    "Use:\n",
    "  Rscript scripts/03_run_mcmc.R --all\n",
    "or:\n",
    "  Rscript scripts/03_run_mcmc.R --id 1"
  )
}


## ============================================================
## 3) Check required files
## ============================================================

required_files <- c(
  DETERMINISTIC_FILE,
  SENSITIVITY_FILE,
  MODEL_FILE
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
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
## 4) Load deterministic and sensitivity outputs
## ============================================================

det <- readRDS(
  DETERMINISTIC_FILE
)

sens <- readRDS(
  SENSITIVITY_FILE
)


required_det <- c(
  "ids",
  "obs_df",
  "theta_init",
  "theta_fit_by_id"
)

missing_det <- setdiff(
  required_det,
  names(det)
)

if (length(missing_det) > 0) {

  stop(
    "Deterministic output is missing: ",
    paste(
      missing_det,
      collapse = ", "
    )
  )
}


ids <- as.numeric(
  det$ids
)

obs_df <- det$obs_df

theta_init <- det$theta_init

theta_fit_by_id <- det$theta_fit_by_id


PRIORITIZED_PARAMS <- sens$prioritized_params

if (is.null(PRIORITIZED_PARAMS)) {

  stop(
    "Sensitivity-prioritized parameters were not found in ",
    SENSITIVITY_FILE
  )
}


cat("\n============================================================\n")
cat("Bayesian PBPK calibration\n")
cat("============================================================\n")

cat(
  "Development IDs: ",
  paste(ids, collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Sensitivity-prioritized parameters: ",
  paste(PRIORITIZED_PARAMS, collapse = ", "),
  "\n",
  sep = ""
)


## ============================================================
## 5) Determine IDs to run
## ============================================================

if (RUN_ALL) {

  RUN_IDS <- ids

} else {

  if (!RUN_ID %in% ids) {

    stop(
      "Requested ID ",
      RUN_ID,
      " is not a development dataset."
    )
  }

  RUN_IDS <- RUN_ID
}


cat(
  "IDs to run: ",
  paste(RUN_IDS, collapse = ", "),
  "\n",
  sep = ""
)


## ============================================================
## 6) Compile PBPK model
## ============================================================

source(
  MODEL_FILE
)

if (!exists("mPBPK_RES")) {

  stop(
    "Object 'mPBPK_RES' was not found in ",
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
  "PBPK model compiled. Model parameters: ",
  length(MODEL_PARAM_NAMES),
  "\n",
  sep = ""
)


## ============================================================
## 7) Parameters sampled by Bayesian MCMC
##
## This closely follows the manuscript-baseline Bayesian
## parameterization.
##
## Inactive Rest-RES parameters are excluded.
##
## Partition coefficients are not sampled here.
## Parameters not sampled by MCMC remain at the baseline
## values specified in model/MP_PBPK.R.
## ============================================================

MCMC_PARAMS <- c(

  ## RES capture / release
  "KLRESmax",
  "KLRESrelease",

  "KSRESmax",
  "KSRESrelease",

  "KKRESmax",
  "KKRESrelease",

  "KLURESmax",
  "KLURESrelease",

  ## Permeability
  "PAL",
  "PAK",
  "PAS",
  "PALU",
  "PAGI",
  "PAR",

  ## Elimination
  "KbileC",
  "KurineC",
  "Kfeces",

  ## GI transport
  "Ktrans",
  "Kgw2b",
  "Kshed"
)


MCMC_PARAMS <- intersect(
  MCMC_PARAMS,
  names(theta_init)
)


if (length(MCMC_PARAMS) == 0) {
  stop("No MCMC parameters matched theta_init.")
}


missing_prioritized <- setdiff(
  PRIORITIZED_PARAMS,
  MCMC_PARAMS
)

if (length(missing_prioritized) > 0) {

  stop(
    "Sensitivity-prioritized parameters absent from ",
    "the MCMC parameter set: ",
    paste(
      missing_prioritized,
      collapse = ", "
    )
  )
}


cat(
  "MCMC parameters (N=",
  length(MCMC_PARAMS),
  "): ",
  paste(MCMC_PARAMS, collapse = ", "),
  "\n",
  sep = ""
)


## ============================================================
## 8) Global baseline parameter vector
##
## Global prior centers come from theta_init.
## ============================================================

theta_base <- theta_init[
  MCMC_PARAMS
]


## ============================================================
## 9) Dataset-specific deterministic starting vectors
##
## Deterministic estimates are used to:
##   - initialize MCMC
##   - define dataset-specific centers for prioritized
##     parameters under the mixed prior strategy
## ============================================================

theta_by_id <- setNames(
  vector(
    "list",
    length(ids)
  ),
  as.character(ids)
)


for (id in ids) {

  th <- theta_base

  deterministic_i <- theta_fit_by_id[[as.character(id)]]

  overlap <- intersect(
    names(th),
    names(deterministic_i)
  )

  th[overlap] <- deterministic_i[
    overlap
  ]

  theta_by_id[[as.character(id)]] <- th
}


## ============================================================
## 10) Bayesian priors
##
## All priors are defined on the log parameter scale.
##
## mixed:
##   sensitivity-prioritized parameters
##     -> dataset-specific deterministic centers
##
##   all other sampled parameters
##     -> global baseline centers
## ============================================================

PRIOR_SD_MAP <- c(

  ## GI transport
  Ktrans = 0.6,
  Kshed  = 0.6,
  Kgw2b  = 0.8,

  ## RES uptake
  KLRESmax  = 0.7,
  KSRESmax  = 0.7,
  KKRESmax  = 0.8,
  KLURESmax = 0.8,

  ## RES release
  KLRESrelease  = 0.5,
  KSRESrelease  = 0.5,
  KKRESrelease  = 0.6,
  KLURESrelease = 0.6,

  ## Permeability
  PAL  = 1.0,
  PAK  = 1.0,
  PAS  = 1.0,
  PALU = 1.0,
  PAGI = 1.0,
  PAR  = 1.0,

  ## Elimination
  KbileC  = 0.6,
  Kfeces  = 0.6,
  KurineC = 0.9
)

PRIOR_SD_DEFAULT <- 0.9


build_prior_for_id <- function(
    id,
    strategy = PRIOR_CENTER_STRATEGY
) {

  param_names <- MCMC_PARAMS


  if (strategy == "global") {

    mu <- theta_base[
      param_names
    ]

  } else if (strategy == "per_id") {

    mu <- theta_by_id[[as.character(id)]][
      param_names
    ]

  } else if (strategy == "mixed") {

    mu <- theta_base[
      param_names
    ]

    id_center <- theta_by_id[[as.character(id)]]

    use_id_center <- intersect(
      PRIORITIZED_PARAMS,
      param_names
    )

    mu[
      use_id_center
    ] <- id_center[
      use_id_center
    ]

  } else {

    stop(
      "Unknown prior strategy: ",
      strategy
    )
  }


  sd <- PRIOR_SD_MAP[
    param_names
  ]

  sd[
    is.na(sd)
  ] <- PRIOR_SD_DEFAULT


  list(
    mu = mu,
    sd = sd
  )
}


## ============================================================
## 11) PBPK prediction function for MCMC
##
## IMPORTANT:
##
## Only parameters sampled by MCMC are supplied to mrgsolve.
##
## Parameters not sampled by MCMC remain at the baseline
## values defined in model/MP_PBPK.R.
##
## This matches the logic of the manuscript-baseline MCMC
## workflow more closely than fixing non-sampled parameters
## at ID-specific deterministic estimates.
## ============================================================

pred_mcmc <- function(
    pars_log,
    id,
    t_end_h = NULL,
    dt_h = 0.01,
    atol = 1e-8,
    rtol = 1e-6
) {

  sub_df <- obs_df %>%
    filter(
      ID == id
    )


  if (nrow(sub_df) == 0) {

    stop(
      "No observations for ID = ",
      id
    )
  }


  ## ----------------------------------------------------------
  ## Convert MCMC parameters to natural scale
  ## ----------------------------------------------------------

  pars_nat <- exp(
    pars_log
  )


  pars_nat <- pars_nat[
    names(pars_nat) %in%
      MODEL_PARAM_NAMES
  ]


  ## ----------------------------------------------------------
  ## Begin with compiled baseline model
  ## ----------------------------------------------------------

  mod_tmp <- mod


  if (length(pars_nat) > 0) {

    mod_tmp <- param(
      mod_tmp,
      pars_nat
    )
  }


  ## ----------------------------------------------------------
  ## Dose
  ## ----------------------------------------------------------

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
      "; found: ",
      paste(
        dose_values,
        collapse = ", "
      )
    )
  }


  dose_mg <- dose_values[1]


  ## ----------------------------------------------------------
  ## Simulation time
  ## ----------------------------------------------------------

  if (is.null(t_end_h)) {

    t_end_h <- max(
      sub_df$Time_h,
      na.rm = TRUE
    )
  }


  add_times <- sort(
    unique(
      as.numeric(
        sub_df$Time_h
      )
    )
  )


  ## ----------------------------------------------------------
  ## Oral dose
  ## ----------------------------------------------------------

  ex <- ev(
    amt = dose_mg,
    cmt = "A_lumen",
    time = 0
  )


  ts <- tgrid(
    start = 0,
    end   = t_end_h,
    delta = dt_h,
    add   = add_times
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


  ## ----------------------------------------------------------
  ## Lung total
  ## ----------------------------------------------------------

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
      "Unable to determine total lung amount."
    )
  }


  ## ----------------------------------------------------------
  ## Outputs as % administered dose
  ## ----------------------------------------------------------

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
## 12) Data residual function
##
## Residuals are calculated on log(1 + %ID), consistent
## with deterministic calibration.
## ============================================================

MCresid_id <- function(
    pars_log,
    id
) {

  sub_df <- obs_df %>%
    filter(
      ID == id
    )


  n_expected <- nrow(
    sub_df
  )


  sim <- tryCatch(

    pred_mcmc(
      pars_log = pars_log,
      id = id,
      t_end_h = max(
        sub_df$Time_h,
        na.rm = TRUE
      ),
      dt_h = 0.01,
      atol = 1e-12,
      rtol = 1e-10
    ),

    error = function(e) NULL
  )


  ## ----------------------------------------------------------
  ## Invalid simulation receives a large penalty
  ## ----------------------------------------------------------

  if (is.null(sim)) {

    return(
      rep(
        1e6,
        n_expected
      )
    )
  }


  sim_long <- sim %>%
    pivot_longer(
      cols = -time,
      names_to = "Tissue",
      values_to = "pred"
    )


  compare <- sub_df %>%
    transmute(
      time = as.numeric(
        Time_h
      ),
      Tissue = as.character(
        Tissue
      ),
      obs = as.numeric(
        Conc
      )
    ) %>%
    inner_join(
      sim_long,
      by = c(
        "time",
        "Tissue"
      )
    )


  if (nrow(compare) != n_expected) {

    return(
      rep(
        1e6,
        n_expected
      )
    )
  }


  residuals <- log1p(
    pmax(
      compare$pred,
      0
    )
  ) -
    log1p(
      pmax(
        compare$obs,
        0
      )
    )


  if (
    any(
      !is.finite(
        residuals
      )
    )
  ) {

    return(
      rep(
        1e6,
        n_expected
      )
    )
  }


  residuals
}


## ============================================================
## 13) Data + prior residual function
## ============================================================

build_resid_fun <- function(
    id,
    prior_mu,
    prior_sd,
    prior_weight = PRIOR_WEIGHT
) {

  force(id)
  force(prior_mu)
  force(prior_sd)
  force(prior_weight)


  function(p) {

    data_res <- MCresid_id(
      p,
      id
    )


    pp <- p[
      names(prior_mu)
    ]


    prior_res <- (
      pp -
        prior_mu
    ) /
      prior_sd *
      prior_weight


    c(
      data_res,
      prior_res
    )
  }
}


## ============================================================
## 14) Initial residual variance calibration
## ============================================================

calibrate_var0 <- function(
    resid_fun,
    p_start,
    n_probe = 5,
    jitter_sd = 0.05
) {

  probes <- vector(
    "list",
    n_probe + 1
  )


  probes[[1]] <- p_start


  if (n_probe > 0) {

    for (k in seq_len(n_probe)) {

      probes[[k + 1]] <- p_start +
        rnorm(
          length(p_start),
          mean = 0,
          sd = jitter_sd
        )
    }
  }


  ss <- vapply(
    probes,
    function(p) {

      r <- resid_fun(
        p
      )

      sum(
        r^2,
        na.rm = TRUE
      ) /
        length(
          r
        )
    },
    numeric(1)
  )


  max(
    1,
    median(
      ss,
      na.rm = TRUE
    )
  )
}


## ============================================================
## 15) Acceptance-rate helper
## ============================================================

calc_acceptance <- function(
    m,
    tol = 1e-10
) {

  P <- as.matrix(
    m$pars
  )


  if (
    is.null(P) ||
    nrow(P) < 2
  ) {

    return(
      NA_real_
    )
  }


  moves <- c(
    FALSE,
    rowSums(
      abs(
        diff(P)
      ) >
        tol
    ) >
      0
  )


  mean(
    moves,
    na.rm = TRUE
  )
}


## ============================================================
## 16) Pilot proposal tuning
## ============================================================

pilot_tune <- function(
    id,
    start_pars,
    prior_mu,
    prior_sd,
    pilot_iter = PILOT_ITER
) {

  ## Reproducible ID-specific pilot seed
  set.seed(
    10000 +
      id
  )


  base <- start_pars +
    rnorm(
      length(start_pars),
      mean = 0,
      sd = 0.05
    )


  resid_fun <- build_resid_fun(
    id = id,
    prior_mu = prior_mu,
    prior_sd = prior_sd
  )


  var0 <- calibrate_var0(
    resid_fun = resid_fun,
    p_start = base,
    n_probe = 5,
    jitter_sd = 0.05
  )


  jump <- INITIAL_JUMP


  for (iteration in 1:4) {

    pilot <- FME::modMCMC(
      f = resid_fun,
      p = base,
      niter = pilot_iter,
      jump = jump,
      var0 = var0,

      lower = log(
        exp(base) /
          100
      ),

      upper = log(
        exp(base) *
          100
      ),

      updatecov = 0,
      wvar0 = 1,

      burninlength = floor(
        pilot_iter *
          BURN_FRACTION
      ),

      ntrydr = NTRY_DR,
      drscale = DR_SCALE
    )


    acc <- calc_acceptance(
      pilot
    )


    cat(
      "  Pilot ",
      iteration,
      ": jump=",
      sprintf("%.3f", jump),
      ", var0=",
      sprintf("%.3f", var0),
      ", acceptance=",
      sprintf("%.3f", acc),
      "\n",
      sep = ""
    )


    if (!is.finite(acc)) {
      acc <- 0
    }


    if (acc < 0.05) {

      var0 <- var0 *
        5

    } else if (
      acc >
        TARGET_ACCEPTANCE[2]
    ) {

      jump <- jump *
        1.6

    } else if (
      acc <
        TARGET_ACCEPTANCE[1]
    ) {

      jump <- jump *
        0.7

    } else {

      break
    }
  }


  list(
    start = base,
    jump = jump,
    var0 = var0
  )
}


## ============================================================
## 17) Run one MCMC chain
## ============================================================

run_chain <- function(
    id,
    start,
    prior_mu,
    prior_sd,
    jump,
    var0,
    seed
) {

  set.seed(
    seed
  )


  resid_fun <- build_resid_fun(
    id = id,
    prior_mu = prior_mu,
    prior_sd = prior_sd
  )


  out <- FME::modMCMC(
    f = resid_fun,
    p = start,
    niter = MAIN_ITER,
    jump = jump,
    var0 = var0,

    lower = log(
      exp(start) /
        100
    ),

    upper = log(
      exp(start) *
        100
    ),

    updatecov = 0,
    wvar0 = 1,

    burninlength = floor(
      MAIN_ITER *
        BURN_FRACTION
    ),

    ntrydr = NTRY_DR,
    drscale = DR_SCALE
  )


  if (is.null(colnames(out$pars))) {

    colnames(
      out$pars
    ) <- names(
      start
    )
  }


  out
}


## ============================================================
## 18) Posterior parameter summary
## ============================================================

summarize_parameters <- function(
    id,
    chains
) {

  draws <- do.call(
    rbind,
    lapply(
      chains,
      function(m) {

        as.data.frame(
          m$pars
        )
      }
    )
  )


  bind_rows(
    lapply(
      names(draws),
      function(parameter) {

        x <- draws[[parameter]]


        tibble(
          ID = id,
          Parameter = parameter,

          mean_log = mean(
            x,
            na.rm = TRUE
          ),

          sd_log = sd(
            x,
            na.rm = TRUE
          ),

          q025_log = as.numeric(
            quantile(
              x,
              0.025,
              na.rm = TRUE
            )
          ),

          median_log = median(
            x,
            na.rm = TRUE
          ),

          q975_log = as.numeric(
            quantile(
              x,
              0.975,
              na.rm = TRUE
            )
          ),

          mean = mean(
            exp(x),
            na.rm = TRUE
          ),

          median = median(
            exp(x),
            na.rm = TRUE
          ),

          q025 = as.numeric(
            quantile(
              exp(x),
              0.025,
              na.rm = TRUE
            )
          ),

          q975 = as.numeric(
            quantile(
              exp(x),
              0.975,
              na.rm = TRUE
            )
          )
        )
      }
    )
  )
}


## ============================================================
## 19) Posterior PBPK parameter-uncertainty summary
## ============================================================

summarize_posterior_predictions <- function(
    id,
    chains
) {

  draws <- do.call(
    rbind,
    lapply(
      chains,
      function(m) {

        as.data.frame(
          m$pars
        )
      }
    )
  )


  n_use <- min(
    NSAMP_PPC,
    nrow(
      draws
    )
  )


  set.seed(
    20000 +
      id
  )


  selected <- sample(
    seq_len(
      nrow(
        draws
      )
    ),
    size = n_use,
    replace = FALSE
  )


  t_end <- max(
    obs_df$Time_h[
      obs_df$ID == id
    ],
    na.rm = TRUE
  )


  simulations <- vector(
    "list",
    n_use
  )


  for (j in seq_along(selected)) {

    pars_j <- unlist(
      draws[
        selected[j],
        ,
        drop = TRUE
      ]
    )


    sim_j <- pred_mcmc(
      pars_log = pars_j,
      id = id,
      t_end_h = t_end,
      dt_h = PPC_DT_H
    ) %>%
      pivot_longer(
        cols = -time,
        names_to = "Tissue",
        values_to = "pred"
      ) %>%
      mutate(
        draw = j
      )


    simulations[[j]] <- sim_j
  }


  bind_rows(
    simulations
  ) %>%
    group_by(
      Tissue,
      time
    ) %>%
    summarise(
      pred_median = median(
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
    ) %>%
    mutate(
      ID = id,
      .before = 1
    )
}


## ============================================================
## 20) Run one development ID
## ============================================================

run_mcmc_id <- function(
    id
) {

  cat("\n============================================================\n")
  cat(
    "Running Bayesian calibration for ID ",
    id,
    "\n",
    sep = ""
  )
  cat("============================================================\n")


  output_dir <- file.path(
    RESULTS_DIR,
    paste0(
      "ID",
      id
    )
  )


  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  prior <- build_prior_for_id(
    id = id
  )


  start <- theta_by_id[[as.character(id)]]


  tuned <- pilot_tune(
    id = id,
    start_pars = start,
    prior_mu = prior$mu,
    prior_sd = prior$sd,
    pilot_iter = PILOT_ITER
  )


  cat(
    "Final pilot settings: jump=",
    sprintf(
      "%.3f",
      tuned$jump
    ),
    ", var0=",
    sprintf(
      "%.3f",
      tuned$var0
    ),
    "\n",
    sep = ""
  )


  chains <- vector(
    "list",
    N_CHAINS
  )


  for (chain_id in seq_len(N_CHAINS)) {

    initialization_seed <- 30000 +
      id *
      100 +
      chain_id


    set.seed(
      initialization_seed
    )


    chain_start <- tuned$start +
      rnorm(
        length(
          tuned$start
        ),
        mean = 0,
        sd = 0.10
      )


    chain_seed <- 40000 +
      id *
      100 +
      CHAIN_SEEDS[
        chain_id
      ]


    cat(
      "Running chain ",
      chain_id,
      " ...\n",
      sep = ""
    )


    chains[[chain_id]] <- run_chain(
      id = id,
      start = chain_start,
      prior_mu = prior$mu,
      prior_sd = prior$sd,
      jump = tuned$jump,
      var0 = tuned$var0,
      seed = chain_seed
    )
  }


  ## ----------------------------------------------------------
  ## Acceptance
  ## ----------------------------------------------------------

  acceptance <- vapply(
    chains,
    calc_acceptance,
    numeric(1)
  )


  acceptance_df <- tibble(
    ID = id,
    Chain = seq_len(
      N_CHAINS
    ),
    Acceptance = acceptance,
    Jump = tuned$jump,
    Var0 = tuned$var0
  )


  ## ----------------------------------------------------------
  ## CODA object
  ## ----------------------------------------------------------

  coda_list <- coda::mcmc.list(
    lapply(
      chains,
      function(m) {

        coda::mcmc(
          as.matrix(
            m$pars
          )
        )
      }
    )
  )


  ## ----------------------------------------------------------
  ## Gelman-Rubin
  ## ----------------------------------------------------------

  gelman <- tryCatch(
    coda::gelman.diag(
      coda_list,
      autoburnin = FALSE
    ),
    error = function(e) NULL
  )


  if (!is.null(gelman)) {

    gelman_df <- tibble(
      ID = id,

      Parameter = rownames(
        gelman$psrf
      ),

      Rhat = gelman$psrf[
        ,
        1
      ],

      Rhat_upper = gelman$psrf[
        ,
        2
      ]
    )

  } else {

    gelman_df <- tibble()
  }


  ## ----------------------------------------------------------
  ## Effective sample size
  ## ----------------------------------------------------------

  ess <- tryCatch(
    coda::effectiveSize(
      coda_list
    ),
    error = function(e) NULL
  )


  if (!is.null(ess)) {

    ess_df <- tibble(
      ID = id,

      Parameter = names(
        ess
      ),

      ESS = as.numeric(
        ess
      )
    )

  } else {

    ess_df <- tibble()
  }


  ## ----------------------------------------------------------
  ## Posterior summaries
  ## ----------------------------------------------------------

  posterior_summary <- summarize_parameters(
    id = id,
    chains = chains
  )


  posterior_predictions <- summarize_posterior_predictions(
    id = id,
    chains = chains
  )


  ## ----------------------------------------------------------
  ## Save outputs
  ## ----------------------------------------------------------

  write.csv(
    posterior_summary,
    file.path(
      output_dir,
      "posterior_parameter_summary.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    posterior_predictions,
    file.path(
      output_dir,
      "posterior_prediction_summary.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    acceptance_df,
    file.path(
      output_dir,
      "acceptance.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    ess_df,
    file.path(
      output_dir,
      "ess.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    gelman_df,
    file.path(
      output_dir,
      "gelman.csv"
    ),
    row.names = FALSE
  )


  result <- list(

    ID = id,

    settings = list(
      prior_strategy = PRIOR_CENTER_STRATEGY,
      prior_weight = PRIOR_WEIGHT,
      pilot_iter = PILOT_ITER,
      main_iter = MAIN_ITER,
      burn_fraction = BURN_FRACTION,
      n_chains = N_CHAINS,
      jump = tuned$jump,
      var0 = tuned$var0
    ),

    prior = prior,

    chains = chains,

    coda = coda_list,

    acceptance = acceptance_df,

    gelman = gelman_df,

    ess = ess_df,

    posterior_summary = posterior_summary,

    posterior_predictions = posterior_predictions
  )


  saveRDS(
    result,
    file.path(
      output_dir,
      "mcmc_results.rds"
    )
  )


  ## ----------------------------------------------------------
  ## Console diagnostics
  ## ----------------------------------------------------------

  cat("\nAcceptance rates:\n")

  print(
    acceptance_df
  )


  if (nrow(ess_df) > 0) {

    cat(
      "Minimum ESS: ",
      round(
        min(
          ess_df$ESS,
          na.rm = TRUE
        ),
        1
      ),
      "\n",
      sep = ""
    )


    cat(
      "Median ESS: ",
      round(
        median(
          ess_df$ESS,
          na.rm = TRUE
        ),
        1
      ),
      "\n",
      sep = ""
    )
  }


  if (nrow(gelman_df) > 0) {

    cat(
      "Maximum Rhat: ",
      round(
        max(
          gelman_df$Rhat,
          na.rm = TRUE
        ),
        3
      ),
      "\n",
      sep = ""
    )
  }


  cat(
    "Saved ID ",
    id,
    " results to ",
    output_dir,
    "\n",
    sep = ""
  )


  result
}


## ============================================================
## 21) Run requested IDs
## ============================================================

RUN_RESULTS <- setNames(
  vector(
    "list",
    length(
      RUN_IDS
    )
  ),
  as.character(
    RUN_IDS
  )
)


for (id in RUN_IDS) {

  RUN_RESULTS[[as.character(id)]] <- run_mcmc_id(
    id = id
  )
}


## ============================================================
## 22) Combined outputs when all IDs are run
## ============================================================

if (RUN_ALL) {

  combined_dir <- file.path(
    RESULTS_DIR,
    "combined"
  )


  dir.create(
    combined_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  posterior_all <- bind_rows(
    lapply(
      RUN_RESULTS,
      function(x) {
        x$posterior_summary
      }
    )
  )


  predictions_all <- bind_rows(
    lapply(
      RUN_RESULTS,
      function(x) {
        x$posterior_predictions
      }
    )
  )


  acceptance_all <- bind_rows(
    lapply(
      RUN_RESULTS,
      function(x) {
        x$acceptance
      }
    )
  )


  ess_all <- bind_rows(
    lapply(
      RUN_RESULTS,
      function(x) {
        x$ess
      }
    )
  )


  gelman_all <- bind_rows(
    lapply(
      RUN_RESULTS,
      function(x) {
        x$gelman
      }
    )
  )


  write.csv(
    posterior_all,
    file.path(
      combined_dir,
      "posterior_parameter_summary_all.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    predictions_all,
    file.path(
      combined_dir,
      "posterior_prediction_summary_all.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    acceptance_all,
    file.path(
      combined_dir,
      "acceptance_all.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    ess_all,
    file.path(
      combined_dir,
      "ess_all.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    gelman_all,
    file.path(
      combined_dir,
      "gelman_all.csv"
    ),
    row.names = FALSE
  )


  saveRDS(
    RUN_RESULTS,
    file.path(
      combined_dir,
      "mcmc_all_ids.rds"
    )
  )
}


## ============================================================
## 23) Final summary
## ============================================================

cat("\n============================================================\n")
cat("Bayesian PBPK calibration completed\n")
cat("============================================================\n")

cat(
  "IDs completed: ",
  paste(
    RUN_IDS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Prior strategy: ",
  PRIOR_CENTER_STRATEGY,
  "\n",
  sep = ""
)

cat(
  "Pilot iterations: ",
  PILOT_ITER,
  "\n",
  sep = ""
)

cat(
  "Main iterations per chain: ",
  MAIN_ITER,
  "\n",
  sep = ""
)

cat(
  "Chains per ID: ",
  N_CHAINS,
  "\n",
  sep = ""
)

cat(
  "Results directory: ",
  RESULTS_DIR,
  "\n",
  sep = ""
)

cat("============================================================\n")