## ============================================================
## 01_fit_model.R
##
## Deterministic PBPK calibration for development datasets
##
## Workflow:
##   1. Load and harmonize development data
##   2. Compile PBPK model
##   3. Perform pre-calibration objective-function sensitivity
##   4. Select common influential parameter set
##   5. Fit each development dataset independently
##   6. Calculate goodness-of-fit metrics
##   7. Save deterministic outputs for downstream analyses
##
## Development datasets: IDs 1–5
##
## This script does NOT perform:
##   - post-calibration AUC-based normalized sensitivity analysis
##   - Bayesian MCMC
##   - random-effects meta-regression
##   - external evaluation
##   - manuscript figure generation
## ============================================================


## ============================================================
## 0) Packages
## ============================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(FME)
  library(dplyr)
  library(tidyr)
  library(numDeriv)
  library(tools)
  library(tibble)
})

options(stringsAsFactors = FALSE)

SEED <- 123
set.seed(SEED)


## ============================================================
## 1) Project settings
## ============================================================

CALIBRATION_IDS <- 1:5

DATA_FILE  <- file.path("data", "data.csv")
MODEL_FILE <- file.path("model", "MP_PBPK.R")

RESULTS_DIR <- file.path("results", "deterministic")

SCREENING_TOP_K <- 15

dir.create(
  RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 2) Check required files
## ============================================================

if (!file.exists(DATA_FILE)) {
  stop(
    "Data file not found: ", DATA_FILE,
    "\nRun this script from the project root directory."
  )
}

if (!file.exists(MODEL_FILE)) {
  stop(
    "PBPK model file not found: ", MODEL_FILE,
    "\nRun this script from the project root directory."
  )
}


## ============================================================
## 3) Load and compile PBPK model
## ============================================================

source(MODEL_FILE)

if (!exists("mPBPK_RES")) {
  stop(
    "Object 'mPBPK_RES' was not found after sourcing ",
    MODEL_FILE
  )
}

mod <- mcode(
  "mPBPK_RES",
  mPBPK_RES
)

MODEL_PARAM_NAMES <- names(param(mod))

cat("\n============================================================\n")
cat("PBPK model compiled successfully\n")
cat("Number of model parameters: ", length(MODEL_PARAM_NAMES), "\n", sep = "")
cat("============================================================\n")


## ============================================================
## 4) Load and harmonize development data
## ============================================================

obs_raw <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "ID",
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
    "Missing required data columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


## ------------------------------------------------------------
## Tissue harmonization helper
## ------------------------------------------------------------

harmonize_tissue <- function(x) {

  x <- trimws(as.character(x))
  x <- tolower(x)

  x <- dplyr::recode(
    x,
    "lung"      = "lungs",
    "gi tract"  = "gi",
    "g.i."      = "gi",
    "gitract"   = "gi",
    "stomach"   = "gi",
    "intestine" = "gi",
    "colon"     = "gi",
    .default    = x
  )

  x <- tools::toTitleCase(x)

  x <- dplyr::recode(
    x,
    "Lung" = "Lungs",
    "Gi"   = "GI",
    .default = x
  )

  x
}


ALLOWED_TISSUES <- c(
  "Blood",
  "Liver",
  "Kidney",
  "Spleen",
  "Lungs",
  "GI",
  "Urine"
)


## ------------------------------------------------------------
## Create SD column if absent
## ------------------------------------------------------------

if (!"SD" %in% names(obs_raw)) {
  obs_raw$SD <- NA_real_
}


obs_df <- obs_raw %>%
  mutate(
    ID      = as.numeric(ID),
    dose_mg = as.numeric(dose_mg),
    Time_h  = as.numeric(Time_h),
    Conc    = as.numeric(Conc),
    Tissue  = harmonize_tissue(Tissue),
    SD      = as.numeric(SD),
    SD = ifelse(
      is.na(SD),
      pmax(Conc * 0.1, 1e-6),
      SD
    )
  ) %>%
  filter(ID %in% CALIBRATION_IDS)


## ------------------------------------------------------------
## Report unsupported tissues
## ------------------------------------------------------------

excluded_tissues <- obs_df %>%
  filter(!Tissue %in% ALLOWED_TISSUES) %>%
  distinct(Tissue) %>%
  pull(Tissue)

if (length(excluded_tissues) > 0) {

  cat(
    "\nExcluded tissues not represented as explicit calibration targets:\n"
  )

  print(excluded_tissues)
}


## ------------------------------------------------------------
## Keep explicit calibration tissues only
## ------------------------------------------------------------

obs_df <- obs_df %>%
  filter(Tissue %in% ALLOWED_TISSUES) %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = ALLOWED_TISSUES
    )
  )


ids <- sort(unique(obs_df$ID))


cat("\nDevelopment IDs found:\n")
print(ids)

cat(
  "\nNumber of development observations: ",
  nrow(obs_df),
  "\n",
  sep = ""
)


if (!all(CALIBRATION_IDS %in% ids)) {
  warning(
    "Not all expected development IDs were found.\n",
    "Expected: ", paste(CALIBRATION_IDS, collapse = ", "),
    "\nFound: ", paste(ids, collapse = ", ")
  )
}


## ============================================================
## 5) PBPK prediction function
## ============================================================

pred <- function(
    pars = NULL,
    id,
    logscale = TRUE,
    t_end_h = NULL,
    dt_h = 0.01,
    atol = 1e-8,
    rtol = 1e-6
) {

  sub_df <- obs_df %>%
    filter(ID == id)

  if (nrow(sub_df) == 0) {
    stop("No observations found for ID = ", id)
  }


  ## Convert log-scale parameters back to natural scale
  pars_nat <- if (is.null(pars)) {

    numeric(0)

  } else if (logscale) {

    exp(pars)

  } else {

    pars
  }


  ## Keep parameters that exist in the compiled model
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


  ## ----------------------------------------------------------
  ## Dose
  ## ----------------------------------------------------------

  dose_values <- unique(
    sub_df$dose_mg[
      is.finite(sub_df$dose_mg)
    ]
  )

  if (length(dose_values) != 1) {
    stop(
      "Expected one unique dose for ID ",
      id,
      "; found: ",
      paste(dose_values, collapse = ", ")
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
      sub_df$Time_h[
        is.finite(sub_df$Time_h)
      ]
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
  ## Lung amount compatibility
  ## ----------------------------------------------------------

  lungs_total <- if ("total_lung" %in% names(out)) {

    out$total_lung

  } else if (
    all(
      c("A_lung_b", "A_lung_t") %in% names(out)
    )
  ) {

    out$A_lung_b + out$A_lung_t

  } else {

    stop(
      "Unable to identify total lung amount from PBPK output."
    )
  }


  ## ----------------------------------------------------------
  ## Convert amounts to % administered dose
  ## ----------------------------------------------------------

  tibble(
    time   = out$time,
    Blood  = (out$A_blood      / dose_mg) * 100,
    Liver  = (out$total_liver  / dose_mg) * 100,
    Kidney = (out$total_kidney / dose_mg) * 100,
    Spleen = (out$total_spleen / dose_mg) * 100,
    Lungs  = (lungs_total      / dose_mg) * 100,
    GI     = (out$total_GI     / dose_mg) * 100,
    Urine  = (out$A_urine      / dose_mg) * 100,
    Rest   = (out$total_rest   / dose_mg) * 100
  )
}


## ============================================================
## 6) Cost and residual functions
## ============================================================

to_numeric_df <- function(df) {

  for (nm in names(df)) {

    if (is.list(df[[nm]])) {

      df[[nm]] <- vapply(
        df[[nm]],
        function(z) {

          z <- unlist(z)

          if (length(z) == 0) {
            return(NA_real_)
          }

          as.numeric(z[1])
        },
        numeric(1)
      )

    } else {

      df[[nm]] <- as.numeric(df[[nm]])
    }
  }

  as.data.frame(df)
}


MCcost_id <- function(
    pars_log,
    id
) {

  sub_df <- obs_df %>%
    filter(ID == id)

  if (nrow(sub_df) == 0) {
    return(
      list(
        ssr = NA_real_,
        residuals = NULL
      )
    )
  }


  ## ----------------------------------------------------------
  ## Simulate model
  ## ----------------------------------------------------------

  sim <- tryCatch(

    pred(
      pars = pars_log,
      id = id,
      logscale = TRUE,
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


  if (is.null(sim)) {
    return(
      list(
        ssr = NA_real_,
        residuals = NULL
      )
    )
  }


  ## ----------------------------------------------------------
  ## Reshape observations
  ## ----------------------------------------------------------

  obs_wide <- sub_df %>%
    transmute(
      time   = as.numeric(Time_h),
      Tissue = as.character(Tissue),
      obs    = as.numeric(Conc)
    ) %>%
    group_by(
      time,
      Tissue
    ) %>%
    summarise(
      obs = mean(
        obs,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Tissue,
      values_from = obs,
      values_fill = NA_real_
    ) %>%
    arrange(time)


  sim <- as.data.frame(sim)
  obs_wide <- as.data.frame(obs_wide)


  vars <- intersect(
    ALLOWED_TISSUES,
    intersect(
      names(sim),
      names(obs_wide)
    )
  )


  if (length(vars) == 0) {
    return(
      list(
        ssr = NA_real_,
        residuals = NULL
      )
    )
  }


  sim_small <- sim[
    ,
    c("time", vars),
    drop = FALSE
  ] %>%
    arrange(time) %>%
    distinct(
      time,
      .keep_all = TRUE
    ) %>%
    to_numeric_df()


  obs_small <- obs_wide[
    ,
    c("time", vars),
    drop = FALSE
  ] %>%
    arrange(time) %>%
    distinct(
      time,
      .keep_all = TRUE
    ) %>%
    to_numeric_df()


  ## ----------------------------------------------------------
  ## Keep tissues with finite observations
  ## ----------------------------------------------------------

  vars <- vars[
    vapply(
      vars,
      function(v) {
        any(
          is.finite(
            obs_small[[v]]
          )
        )
      },
      logical(1)
    )
  ]


  if (length(vars) == 0) {
    return(
      list(
        ssr = NA_real_,
        residuals = NULL
      )
    )
  }


  sim_small <- sim_small[
    ,
    c("time", vars),
    drop = FALSE
  ]

  obs_small <- obs_small[
    ,
    c("time", vars),
    drop = FALSE
  ]


  ## ----------------------------------------------------------
  ## Deterministic fitting objective:
  ## log(1 + %ID)
  ## ----------------------------------------------------------

  sim_small[, vars] <- lapply(
    sim_small[
      ,
      vars,
      drop = FALSE
    ],
    function(x) {
      log1p(
        pmax(
          x,
          0
        )
      )
    }
  )


  obs_small[, vars] <- lapply(
    obs_small[
      ,
      vars,
      drop = FALSE
    ],
    function(x) {
      log1p(
        pmax(
          x,
          0
        )
      )
    }
  )


  ## ----------------------------------------------------------
  ## FME model cost
  ## ----------------------------------------------------------

  cost_i <- tryCatch(

    suppressWarnings(
      FME::modCost(
        model = sim_small,
        obs = obs_small,
        x = "time"
      )
    ),

    error = function(e) NULL
  )


  if (is.null(cost_i)) {
    return(
      list(
        ssr = NA_real_,
        residuals = NULL
      )
    )
  }


  ssr_i <- if (!is.null(cost_i$ssr)) {

    cost_i$ssr

  } else if (!is.null(cost_i$residuals)) {

    sum(
      cost_i$residuals$res^2,
      na.rm = TRUE
    )

  } else {

    NA_real_
  }


  list(
    ssr = ssr_i,
    residuals = as.data.frame(
      cost_i$residuals
    )
  )
}


MCcost_scalar_id <- function(
    pars_log,
    id
) {

  out <- tryCatch(
    MCcost_id(
      pars_log,
      id
    ),
    error = function(e) NULL
  )


  if (
    is.null(out) ||
    !is.finite(out$ssr) ||
    out$ssr <= 0
  ) {

    return(
      c(logSSR = NA_real_)
    )
  }


  c(
    logSSR = log10(
      out$ssr
    )
  )
}


MCresid_id <- function(
    pars_log,
    id
) {

  out <- MCcost_id(
    pars_log,
    id
  )


  if (
    is.null(out$residuals) ||
    nrow(out$residuals) == 0
  ) {

    return(
      numeric(0)
    )
  }


  as.numeric(
    out$residuals$res
  )
}


## ============================================================
## 7) Initial particle-specific parameter vector
## ============================================================

theta_init <- log(
  c(

    ## --------------------------------------------------------
    ## RES capture / release
    ## --------------------------------------------------------
    KLRESmax     = 0.4,
    KLRESrelease = 7.5e-3,

    KSRESmax     = 1,
    KSRESrelease = 3e-3,

    KKRESmax     = 0.1,
    KKRESrelease = 1e-2,

    KLURESmax     = 0.1,
    KLURESrelease = 5e-3,


    ## --------------------------------------------------------
    ## Permeability
    ## --------------------------------------------------------
    PAL  = 1e-3,
    PAK  = 1e-3,
    PAS  = 1e-3,
    PALU = 1e-3,
    PAGI = 1e-3,
    PAR  = 1e-6,


    ## --------------------------------------------------------
    ## Partition coefficients
    ## --------------------------------------------------------
    KpL  = 0.01,
    KpK  = 0.15,
    KpS  = 0.15,
    KpR  = 0.15,
    KpLU = 0.15,


    ## --------------------------------------------------------
    ## Elimination
    ## --------------------------------------------------------
    KbileC  = 2.2e-3,
    KurineC = 3.2e-4,
    Kfeces  = 0.53,


    ## --------------------------------------------------------
    ## GI transport
    ## --------------------------------------------------------
    Ktrans = 5e-1,
    Kgw2b  = 3e-4,
    Kshed  = 5e-1

  )
)


theta_init <- theta_init[
  names(theta_init) %in% MODEL_PARAM_NAMES
]


if (length(theta_init) == 0) {
  stop(
    "No initial parameters matched the PBPK model."
  )
}


cat("\nInitial particle-specific parameters:\n")
print(names(theta_init))


## ============================================================
## 8) Pre-calibration objective-function sensitivity screening
## ============================================================

sens_for_id <- function(
    theta_log,
    id,
    eps = 1e-4
) {

  g <- tryCatch(

    numDeriv::grad(
      func = function(p) {

        MCcost_scalar_id(
          p,
          id
        )[["logSSR"]]

      },
      x = theta_log,
      method.args = list(
        eps = eps
      )
    ),

    error = function(e) {
      rep(
        NA_real_,
        length(theta_log)
      )
    }
  )


  tibble(
    ID = id,
    Parameter = names(theta_log),
    Sensitivity = as.numeric(g),
    AbsSensitivity = abs(
      as.numeric(g)
    )
  ) %>%
    arrange(
      desc(
        AbsSensitivity
      )
    )
}


compute_consensus_params <- function(
    theta_log,
    ids,
    topK_each = 15
) {

  sensitivity_table <- bind_rows(
    lapply(
      ids,
      function(id) {
        sens_for_id(
          theta_log,
          id
        )
      }
    )
  )


  sensitivity_table <- sensitivity_table %>%
    group_by(ID) %>%
    mutate(
      Rank = rank(
        -AbsSensitivity,
        ties.method = "min"
      )
    ) %>%
    ungroup()


  top_sets <- split(
    sensitivity_table,
    sensitivity_table$ID
  )


  top_sets <- lapply(
    top_sets,
    function(df) {

      df %>%
        arrange(
          desc(
            AbsSensitivity
          )
        ) %>%
        slice_head(
          n = topK_each
        ) %>%
        pull(Parameter)
    }
  )


  common_params <- Reduce(
    intersect,
    top_sets
  )


  if (length(common_params) == 0) {
    stop(
      "No parameters were shared among the top ",
      topK_each,
      " parameters across all development datasets."
    )
  }


  aggregate_table <- sensitivity_table %>%
    group_by(Parameter) %>%
    summarise(
      MeanAbsSensitivity = mean(
        AbsSensitivity,
        na.rm = TRUE
      ),
      MaxAbsSensitivity = max(
        AbsSensitivity,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    arrange(
      desc(
        MeanAbsSensitivity
      )
    )


  list(
    sensitivity_table = sensitivity_table,
    aggregate_table = aggregate_table,
    common_params = common_params
  )
}


screening <- compute_consensus_params(
  theta_log = theta_init,
  ids = ids,
  topK_each = SCREENING_TOP_K
)


common_params <- screening$common_params


cat("\n============================================================\n")
cat("Common deterministic calibration parameters\n")
cat("============================================================\n")

print(common_params)

cat(
  "\nNumber of fitted parameters: ",
  length(common_params),
  "\n",
  sep = ""
)


if (length(common_params) != 10) {
  warning(
    "Expected approximately 10 common fitted parameters based on ",
    "the manuscript baseline, but found ",
    length(common_params),
    "."
  )
}


## ============================================================
## 9) Deterministic fitting
## ============================================================

fit_one_id <- function(
    id,
    theta_log_init,
    parm_names,
    widen = 1e4,
    trace = 0
) {

  GI_BOUND_PARAMS <- c(
    "Ktrans",
    "Kshed",
    "Kgw2b"
  )


  GI_BOUNDS <- data.frame(
    lower = c(
      1e-3,
      0.05,
      1e-4
    ),
    upper = c(
      0.9,
      3,
      1e-2
    ),
    row.names = GI_BOUND_PARAMS
  )


  theta_sel <- theta_log_init[
    parm_names
  ]

  theta_nat <- exp(
    theta_sel
  )


  lower_nat <- theta_nat /
    widen

  upper_nat <- theta_nat *
    widen


  for (nm in parm_names) {

    if (nm %in% GI_BOUND_PARAMS) {

      lower_nat[nm] <- GI_BOUNDS[
        nm,
        "lower"
      ]

      upper_nat[nm] <- GI_BOUNDS[
        nm,
        "upper"
      ]
    }
  }


  lower_log <- log(
    lower_nat
  )

  upper_log <- log(
    upper_nat
  )


  theta_sel <- pmin(
    pmax(
      theta_sel,
      lower_log + 1e-9
    ),
    upper_log - 1e-9
  )


  FME::modFit(
    f = function(p) {
      MCresid_id(
        p,
        id
      )
    },
    p = theta_sel,
    lower = lower_log,
    upper = upper_log,
    method = "Port",
    control = list(
      eval.max = 250,
      iter.max = 250,
      trace = trace
    )
  )
}


cat("\n============================================================\n")
cat("Running deterministic calibration\n")
cat("============================================================\n")


Fit_list <- lapply(
  ids,
  function(id) {

    cat(
      "Fitting ID ",
      id,
      " ...\n",
      sep = ""
    )

    fit_one_id(
      id = id,
      theta_log_init = theta_init,
      parm_names = common_params,
      trace = 0
    )
  }
)


## Name by actual dataset ID rather than list position
names(Fit_list) <- as.character(ids)


## ============================================================
## 10) Build full fitted parameter vectors
## ============================================================

theta_fit_by_id <- setNames(
  vector(
    "list",
    length(ids)
  ),
  as.character(ids)
)


for (id in ids) {

  th <- theta_init

  fit_pars <- Fit_list[[as.character(id)]]$par

  overlap <- intersect(
    names(th),
    names(fit_pars)
  )

  th[overlap] <- fit_pars[
    overlap
  ]

  theta_fit_by_id[[as.character(id)]] <- th
}


## ------------------------------------------------------------
## Fitted parameter table
## ------------------------------------------------------------

fit_parameter_table <- bind_rows(
  lapply(
    names(theta_fit_by_id),
    function(id_chr) {

      pars_i <- theta_fit_by_id[[id_chr]]

      tibble(
        ID = as.numeric(id_chr),
        Parameter = names(pars_i),
        log_value = as.numeric(pars_i),
        value = exp(
          as.numeric(pars_i)
        ),
        Estimated = names(pars_i) %in%
          common_params
      )
    }
  )
)


cat("\n============================================================\n")
cat("Fitted parameter estimates\n")
cat("============================================================\n")

print(
  fit_parameter_table %>%
    filter(Estimated)
)


## ============================================================
## 11) Collect observed versus fitted predictions
## ============================================================

collect_gof_data <- function(
    ids,
    Fit_list
) {

  bind_rows(
    lapply(
      ids,
      function(id) {

        fit_obj <- Fit_list[[as.character(id)]]


        sim_fit <- pred(
          pars = fit_obj$par,
          id = id,
          logscale = TRUE
        ) %>%
          pivot_longer(
            -time,
            names_to = "Tissue",
            values_to = "pred"
          )


        obs_plot <- obs_df %>%
          filter(
            ID == id
          ) %>%
          transmute(
            ID = ID,
            time = as.numeric(
              Time_h
            ),
            Tissue = as.character(
              Tissue
            ),
            obs = as.numeric(
              Conc
            ),
            SD = as.numeric(
              SD
            )
          )


        sim_fit %>%
          inner_join(
            obs_plot,
            by = c(
              "time",
              "Tissue"
            )
          ) %>%
          mutate(
            ID = id
          )
      }
    )
  )
}


gof_all <- collect_gof_data(
  ids = ids,
  Fit_list = Fit_list
)


cat(
  "\nNumber of observed-predicted pairs: ",
  nrow(gof_all),
  "\n",
  sep = ""
)


## ============================================================
## 12) Goodness-of-fit metrics
## ============================================================

## ------------------------------------------------------------
## Primary deterministic GOF
##
## To preserve consistency with the manuscript baseline,
## R2, RMSE, MAPE, and fold error are calculated using
## positive observed and predicted values.
##
## We retain n_all separately for transparency.
## ------------------------------------------------------------

calc_gof_metrics <- function(df) {

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


  if (nrow(positive_df) == 0) {

    return(
      tibble(
        R2 = NA_real_,
        RMSE = NA_real_,
        MAPE = NA_real_,
        n_all = nrow(finite_df),
        n_positive = 0
      )
    )
  }


  r2 <- if (
    length(
      unique(
        positive_df$obs
      )
    ) > 1 &&
    length(
      unique(
        positive_df$pred
      )
    ) > 1
  ) {

    cor(
      positive_df$obs,
      positive_df$pred,
      use = "complete.obs"
    )^2

  } else {

    NA_real_
  }


  rmse <- sqrt(
    mean(
      (
        positive_df$pred -
        positive_df$obs
      )^2,
      na.rm = TRUE
    )
  )


  mape <- mean(
    abs(
      positive_df$pred -
      positive_df$obs
    ) /
      positive_df$obs,
    na.rm = TRUE
  ) * 100


  tibble(
    R2 = r2,
    RMSE = rmse,
    MAPE = mape,
    n_all = nrow(finite_df),
    n_positive = nrow(positive_df)
  )
}


calc_fold_metrics <- function(df) {

  valid <- df %>%
    filter(
      is.finite(obs),
      is.finite(pred),
      obs > 0,
      pred > 0
    ) %>%
    mutate(
      ratio = pred / obs
    )


  if (nrow(valid) == 0) {

    return(
      tibble(
        n_positive = 0,
        within_2fold_n = 0,
        within_2fold_pct = NA_real_,
        within_3fold_n = 0,
        within_3fold_pct = NA_real_
      )
    )
  }


  tibble(
    n_positive = nrow(valid),

    within_2fold_n = sum(
      valid$ratio >= 0.5 &
      valid$ratio <= 2
    ),

    within_2fold_pct = mean(
      valid$ratio >= 0.5 &
      valid$ratio <= 2
    ) * 100,

    within_3fold_n = sum(
      valid$ratio >= 1 / 3 &
      valid$ratio <= 3
    ),

    within_3fold_pct = mean(
      valid$ratio >= 1 / 3 &
      valid$ratio <= 3
    ) * 100
  )
}


## ------------------------------------------------------------
## Overall GOF
## ------------------------------------------------------------

overall_gof <- calc_gof_metrics(
  gof_all
)

overall_fold <- calc_fold_metrics(
  gof_all
)


overall_summary <- bind_cols(
  overall_gof,
  overall_fold %>%
    select(
      -n_positive
    )
)


## ------------------------------------------------------------
## Per-ID GOF
## ------------------------------------------------------------

per_id_metrics <- bind_rows(
  lapply(
    ids,
    function(id) {

      bind_cols(
        tibble(
          ID = id
        ),
        calc_gof_metrics(
          gof_all %>%
            filter(
              ID == id
            )
        )
      )
    }
  )
)


## ------------------------------------------------------------
## Per-tissue GOF
## ------------------------------------------------------------

per_tissue_metrics <- bind_rows(
  lapply(
    sort(
      unique(
        gof_all$Tissue
      )
    ),
    function(tt) {

      bind_cols(
        tibble(
          Tissue = tt
        ),
        calc_gof_metrics(
          gof_all %>%
            filter(
              Tissue == tt
            )
        )
      )
    }
  )
)


cat("\n============================================================\n")
cat("Deterministic goodness-of-fit\n")
cat("============================================================\n")

cat("\nOverall:\n")
print(overall_summary)

cat("\nBy development dataset:\n")
print(per_id_metrics)

cat("\nBy tissue:\n")
print(per_tissue_metrics)


## ============================================================
## 13) Save outputs
## ============================================================

## ------------------------------------------------------------
## Pre-calibration sensitivity
## ------------------------------------------------------------

write.csv(
  screening$sensitivity_table,
  file.path(
    RESULTS_DIR,
    "precalibration_sensitivity_by_id.csv"
  ),
  row.names = FALSE
)


write.csv(
  screening$aggregate_table,
  file.path(
    RESULTS_DIR,
    "precalibration_sensitivity_summary.csv"
  ),
  row.names = FALSE
)


## ------------------------------------------------------------
## Fitted parameters
## ------------------------------------------------------------

write.csv(
  fit_parameter_table,
  file.path(
    RESULTS_DIR,
    "fitted_parameter_table.csv"
  ),
  row.names = FALSE
)


## ------------------------------------------------------------
## GOF source data
## ------------------------------------------------------------

write.csv(
  gof_all,
  file.path(
    RESULTS_DIR,
    "gof_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  overall_summary,
  file.path(
    RESULTS_DIR,
    "overall_gof_metrics.csv"
  ),
  row.names = FALSE
)


write.csv(
  per_id_metrics,
  file.path(
    RESULTS_DIR,
    "per_id_gof_metrics.csv"
  ),
  row.names = FALSE
)


write.csv(
  per_tissue_metrics,
  file.path(
    RESULTS_DIR,
    "per_tissue_gof_metrics.csv"
  ),
  row.names = FALSE
)


## ------------------------------------------------------------
## Full R object for downstream scripts
## ------------------------------------------------------------

deterministic_outputs <- list(

  seed = SEED,

  calibration_ids = CALIBRATION_IDS,

  ids = ids,

  obs_df = obs_df,

  model_param_names = MODEL_PARAM_NAMES,

  theta_init = theta_init,

  screening_top_k = SCREENING_TOP_K,

  common_params = common_params,

  screening = screening,

  Fit_list = Fit_list,

  theta_fit_by_id = theta_fit_by_id,

  fit_parameter_table = fit_parameter_table,

  gof_all = gof_all,

  overall_summary = overall_summary,

  per_id_metrics = per_id_metrics,

  per_tissue_metrics = per_tissue_metrics
)


saveRDS(
  deterministic_outputs,
  file = file.path(
    RESULTS_DIR,
    "deterministic_outputs.rds"
  )
)


## ============================================================
## 14) Final run summary
## ============================================================

cat("\n============================================================\n")
cat("Deterministic calibration completed successfully\n")
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
  "Development observations: ",
  nrow(obs_df),
  "\n",
  sep = ""
)

cat(
  "Observed-predicted pairs: ",
  nrow(gof_all),
  "\n",
  sep = ""
)

cat(
  "Common fitted parameters (N=",
  length(common_params),
  "): ",
  paste(
    common_params,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat("\nOverall deterministic GOF:\n")
print(overall_summary)

cat(
  "\nOutputs saved to: ",
  RESULTS_DIR,
  "\n",
  sep = ""
)

cat("============================================================\n")