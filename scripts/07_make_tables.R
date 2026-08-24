## ============================================================
## 07_make_tables.R
##
## Generate manuscript/SI analysis tables from FINAL saved
## outputs produced by Scripts 01–05.
##
## IMPORTANT
## - Does NOT rerun deterministic fitting
## - Does NOT rerun sensitivity analysis
## - Does NOT rerun MCMC
## - Does NOT rerun meta-regression
## - Does NOT rerun external PBPK simulations
##
## Generated:
##
## Main manuscript
##   Table 2  Deterministic fitted parameters, if recoverable
##
## Supporting Information
##   Table S3  Deterministic GOF
##   Table S4  Detailed AUC-based NSCs
##   Table S5  Overall NSC ranking
##   Table S6  MCMC diagnostics
##   Table S7  Posterior estimates of prioritized parameters
##   Table S8  Development-set 95% interval coverage
##   Table S9  Random-effects meta-regression
##   Table S10 External evaluation
##
## Tables S1 and S2 are model-definition/reference tables
## and are intentionally not regenerated here.
## ============================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
})

options(stringsAsFactors = FALSE)


## ============================================================
## 0) DIRECTORIES
## ============================================================

DET_DIR  <- file.path("results", "deterministic")
SENS_DIR <- file.path("results", "sensitivity")
MCMC_DIR <- file.path("results", "mcmc")
META_DIR <- file.path("results", "meta")
EXT_DIR  <- file.path("results", "external")

TABLE_DIR <- "tables"

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 1) CONSTANTS
## ============================================================

DEV_IDS <- 1:5

PRIORITIZED_PARAMS <- c(
  "Kfeces",
  "Ktrans",
  "Kgw2b",
  "Kshed",
  "KLRESmax",
  "KurineC"
)

META_PARAM_ORDER <- c(
  "Ktrans",
  "Kgw2b",
  "Kshed",
  "KLRESmax",
  "Kfeces",
  "KurineC"
)

DETERMINISTIC_PARAMS <- c(
  "Kfeces",
  "KLRESmax",
  "Kshed",
  "Kgw2b",
  "Ktrans",
  "KLRESrelease",
  "KSRESmax",
  "KbileC",
  "KurineC",
  "KpK"
)

ID_SIZE_MAP <- c(
  "1" = 20,
  "2" = 220,
  "3" = 1000,
  "4" = 6000,
  "5" = 220,
  "6" = 50,
  "7" = 500,
  "8" = 5000
)

TISSUE_ORDER <- c(
  "Blood",
  "GI",
  "Spleen",
  "Kidney",
  "Liver",
  "Lungs",
  "Urine",
  "Rest"
)

AUC_LABEL_MAP <- c(
  "Blood"  = "AUCB",
  "GI"     = "AUCGI",
  "Spleen" = "AUCS",
  "Kidney" = "AUCK",
  "Liver"  = "AUCL",
  "Lungs"  = "AUCLu",
  "Urine"  = "AUCU",
  "Rest"   = "AUCRt"
)


## ============================================================
## 2) HELPERS
## ============================================================

first_existing <- function(paths) {

  hit <- paths[
    file.exists(paths)
  ]

  if (length(hit) == 0) {
    return(NULL)
  }

  hit[1]
}


find_col <- function(
    df,
    candidates,
    required = TRUE
) {

  nm <- trimws(
    names(df)
  )

  names(df) <- nm

  hit <- candidates[
    candidates %in% nm
  ]

  if (length(hit) > 0) {
    return(hit[1])
  }

  if (required) {

    stop(
      "Could not find required column. Tried: ",
      paste(
        candidates,
        collapse = ", "
      ),
      "\nAvailable columns: ",
      paste(
        nm,
        collapse = ", "
      )
    )
  }

  NA_character_
}


harmonize_tissue <- function(x) {

  y <- trimws(
    as.character(x)
  )

  yl <- tolower(y)

  yl <- dplyr::recode(
    yl,
    "lung" = "lungs",
    "gi tract" = "gi",
    "g.i." = "gi",
    "gitract" = "gi",
    "stomach" = "gi",
    "intestine" = "gi",
    "small intestine" = "gi",
    "large intestine" = "gi",
    "colon" = "gi",
    .default = yl
  )

  out <- tools::toTitleCase(
    yl
  )

  dplyr::recode(
    out,
    "Gi" = "GI",
    "Lung" = "Lungs",
    .default = out
  )
}


size_label <- function(id) {

  size_nm <- as.numeric(
    ID_SIZE_MAP[
      as.character(id)
    ]
  )

  case_when(
    size_nm == 20   ~ "20 nm",
    size_nm == 220  ~ "220 nm",
    size_nm == 1000 ~ "1 µm",
    size_nm == 6000 ~ "6 µm",
    size_nm == 50   ~ "50 nm",
    size_nm == 500  ~ "500 nm",
    size_nm == 5000 ~ "5 µm",
    TRUE ~ paste0(
      size_nm,
      " nm"
    )
  )
}


format_num <- function(
    x,
    digits = 3
) {

  ifelse(
    is.na(x),
    NA_character_,
    ifelse(
      x == 0,
      "0",
      ifelse(
        abs(x) < 0.001 |
          abs(x) >= 10000,
        formatC(
          x,
          format = "e",
          digits = 2
        ),
        formatC(
          x,
          format = "fg",
          digits = digits,
          flag = "#"
        )
      )
    )
  )
}


format_ci <- function(
    median,
    lower,
    upper
) {

  paste0(
    format_num(median),
    " (",
    format_num(lower),
    "–",
    format_num(upper),
    ")"
  )
}


safe_r2 <- function(
    obs,
    pred
) {

  ok <- is.finite(obs) &
    is.finite(pred)

  obs <- obs[ok]
  pred <- pred[ok]

  if (
    length(obs) < 2 ||
      length(unique(obs)) < 2 ||
      length(unique(pred)) < 2
  ) {

    return(NA_real_)
  }

  cor(
    obs,
    pred
  )^2
}


safe_rmse <- function(
    obs,
    pred
) {

  ok <- is.finite(obs) &
    is.finite(pred)

  if (!any(ok)) {
    return(NA_real_)
  }

  sqrt(
    mean(
      (
        pred[ok] -
          obs[ok]
      )^2
    )
  )
}


safe_mape <- function(
    obs,
    pred
) {

  ok <- is.finite(obs) &
    is.finite(pred) &
    obs > 0

  if (!any(ok)) {
    return(NA_real_)
  }

  mean(
    abs(
      pred[ok] -
        obs[ok]
    ) /
      obs[ok]
  ) * 100
}


first_numeric_column <- function(
    df,
    exclude = character()
) {

  candidates <- names(df)[
    vapply(
      df,
      is.numeric,
      logical(1)
    )
  ]

  candidates <- setdiff(
    candidates,
    exclude
  )

  if (length(candidates) == 0) {
    return(NA_character_)
  }

  candidates[1]
}


## ============================================================
## 3) LOAD DETERMINISTIC GOF
## ============================================================

find_gof_df_recursive <- function(x) {

  if (is.data.frame(x)) {

    nm <- names(x)

    has_obs <- any(
      c(
        "obs",
        "Obs",
        "Observed",
        "observed"
      ) %in% nm
    )

    has_pred <- any(
      c(
        "pred",
        "Pred",
        "Predicted",
        "predicted"
      ) %in% nm
    )

    has_id <- any(
      c(
        "ID",
        "id"
      ) %in% nm
    )

    has_tissue <- any(
      c(
        "Tissue",
        "tissue"
      ) %in% nm
    )

    if (
      has_obs &&
        has_pred &&
        has_id &&
        has_tissue
    ) {
      return(x)
    }
  }

  if (is.list(x)) {

    for (i in seq_along(x)) {

      ans <- find_gof_df_recursive(
        x[[i]]
      )

      if (!is.null(ans)) {
        return(ans)
      }
    }
  }

  NULL
}


load_deterministic_gof <- function() {

  rds_file <- first_existing(
    c(
      file.path(
        DET_DIR,
        "deterministic_outputs.rds"
      ),
      file.path(
        DET_DIR,
        "fitted_objects.rds"
      )
    )
  )

  if (!is.null(rds_file)) {

    obj <- readRDS(
      rds_file
    )

    ans <- find_gof_df_recursive(
      obj
    )

    if (!is.null(ans)) {

      message(
        "Deterministic GOF source: ",
        rds_file
      )

      return(ans)
    }
  }


  csv_files <- list.files(
    DET_DIR,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )


  for (f in csv_files) {

    x <- tryCatch(
      read_csv(
        f,
        show_col_types = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(x)) {
      next
    }

    nm <- names(x)

    if (
      any(
        c(
          "obs",
          "Observed",
          "observed"
        ) %in% nm
      ) &&
        any(
          c(
            "pred",
            "Predicted",
            "predicted"
          ) %in% nm
        ) &&
        any(
          c(
            "ID",
            "id"
          ) %in% nm
        )
    ) {
      return(x)
    }
  }

  stop(
    "Could not recover deterministic GOF data."
  )
}


DET_GOF_RAW <- load_deterministic_gof()


det_obs_col <- find_col(
  DET_GOF_RAW,
  c(
    "obs",
    "Obs",
    "Observed",
    "observed"
  )
)

det_pred_col <- find_col(
  DET_GOF_RAW,
  c(
    "pred",
    "Pred",
    "Predicted",
    "predicted"
  )
)

det_id_col <- find_col(
  DET_GOF_RAW,
  c(
    "ID",
    "id"
  )
)

det_tissue_col <- find_col(
  DET_GOF_RAW,
  c(
    "Tissue",
    "tissue"
  )
)


DET_GOF <- DET_GOF_RAW %>%

  transmute(
    ID = as.numeric(
      .data[[det_id_col]]
    ),

    Tissue = harmonize_tissue(
      .data[[det_tissue_col]]
    ),

    obs = as.numeric(
      .data[[det_obs_col]]
    ),

    pred = as.numeric(
      .data[[det_pred_col]]
    )
  ) %>%

  filter(
    ID %in% DEV_IDS
  )


## ============================================================
## 4) TABLE S3
## DETERMINISTIC GOODNESS-OF-FIT
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table S3: deterministic GOF\n"
)

cat(
  "============================================================\n"
)


S3_OVERALL <- tibble(
  Section =
    "Overall and dataset-specific",

  Group =
    "Overall",

  Particle_size =
    "20 nm–6 µm",

  N =
    nrow(
      DET_GOF
    ),

  R2 =
    safe_r2(
      DET_GOF$obs,
      DET_GOF$pred
    ),

  RMSE =
    safe_rmse(
      DET_GOF$obs,
      DET_GOF$pred
    ),

  MAPE =
    safe_mape(
      DET_GOF$obs,
      DET_GOF$pred
    )
)


S3_BY_ID <- bind_rows(
  lapply(
    DEV_IDS,
    function(id) {

      x <- DET_GOF %>%
        filter(
          ID == id
        )

      tibble(
        Section =
          "Overall and dataset-specific",

        Group =
          paste0(
            "ID ",
            id
          ),

        Particle_size =
          size_label(id),

        N =
          nrow(x),

        R2 =
          safe_r2(
            x$obs,
            x$pred
          ),

        RMSE =
          safe_rmse(
            x$obs,
            x$pred
          ),

        MAPE =
          safe_mape(
            x$obs,
            x$pred
          )
      )
    }
  )
)


S3_BY_TISSUE <- DET_GOF %>%

  group_by(
    Tissue
  ) %>%

  summarise(
    Section =
      "Tissue- and excretion-specific",

    Group =
      first(Tissue),

    Particle_size =
      NA_character_,

    N =
      n(),

    R2 =
      safe_r2(
        obs,
        pred
      ),

    RMSE =
      safe_rmse(
        obs,
        pred
      ),

    MAPE =
      safe_mape(
        obs,
        pred
      ),

    .groups =
      "drop"
  ) %>%

  select(
    Section,
    Group,
    Particle_size,
    N,
    R2,
    RMSE,
    MAPE
  )


TABLE_S3 <- bind_rows(
  S3_OVERALL,
  S3_BY_ID,
  S3_BY_TISSUE
)


write_csv(
  TABLE_S3,
  file.path(
    TABLE_DIR,
    "Table_S3_deterministic_GOF.csv"
  )
)


print(
  TABLE_S3
)


## ============================================================
## 5) TABLE S4
## DETAILED NORMALIZED SENSITIVITY COEFFICIENTS
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table S4: detailed NSCs\n"
)

cat(
  "============================================================\n"
)


nsc_file <- first_existing(
  c(
    file.path(
      SENS_DIR,
      "nsc_by_id_tissue_parameter.csv"
    ),
    file.path(
      SENS_DIR,
      "NSC_amount_all.csv"
    )
  )
)


if (is.null(nsc_file)) {
  stop(
    "Missing detailed sensitivity output."
  )
}


NSC_RAW <- read_csv(
  nsc_file,
  show_col_types = FALSE
)


nsc_id_col <- find_col(
  NSC_RAW,
  c(
    "ID",
    "id"
  )
)

nsc_param_col <- find_col(
  NSC_RAW,
  c(
    "Parameter",
    "parameter",
    "param"
  )
)

nsc_tissue_col <- find_col(
  NSC_RAW,
  c(
    "Tissue",
    "tissue",
    "Endpoint",
    "endpoint"
  )
)

nsc_value_col <- find_col(
  NSC_RAW,
  c(
    "NSC",
    "nsc"
  )
)


TABLE_S4_LONG <- NSC_RAW %>%

  transmute(
    ID =
      as.numeric(
        .data[[nsc_id_col]]
      ),

    Particle_size_nm =
      as.numeric(
        ID_SIZE_MAP[
          as.character(ID)
        ]
      ),

    Parameter =
      as.character(
        .data[[nsc_param_col]]
      ),

    Tissue =
      harmonize_tissue(
        .data[[nsc_tissue_col]]
      ),

    AUC_endpoint =
      unname(
        AUC_LABEL_MAP[
          Tissue
        ]
      ),

    NSC =
      as.numeric(
        .data[[nsc_value_col]]
      )
  ) %>%

  arrange(
    ID,
    Parameter,
    factor(
      Tissue,
      levels = TISSUE_ORDER
    )
  )


TABLE_S4_WIDE <- TABLE_S4_LONG %>%

  select(
    ID,
    Particle_size_nm,
    Parameter,
    AUC_endpoint,
    NSC
  ) %>%

  pivot_wider(
    names_from =
      AUC_endpoint,

    values_from =
      NSC
  ) %>%

  arrange(
    ID,
    Parameter
  )


write_csv(
  TABLE_S4_LONG,
  file.path(
    TABLE_DIR,
    "Table_S4_NSC_long.csv"
  )
)


write_csv(
  TABLE_S4_WIDE,
  file.path(
    TABLE_DIR,
    "Table_S4_NSC_wide.csv"
  )
)


cat(
  "Table S4 rows: ",
  nrow(
    TABLE_S4_LONG
  ),
  "\n",
  sep = ""
)


## ============================================================
## 6) TABLE S5
## OVERALL NSC PARAMETER RANKING
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table S5: sensitivity ranking\n"
)

cat(
  "============================================================\n"
)


rank_file <- file.path(
  SENS_DIR,
  "nsc_parameter_ranking.csv"
)


if (!file.exists(rank_file)) {
  stop(
    "Missing: ",
    rank_file
  )
}


RANK_RAW <- read_csv(
  rank_file,
  show_col_types = FALSE
)


rank_param_col <- find_col(
  RANK_RAW,
  c(
    "Parameter",
    "parameter",
    "param"
  )
)

rank_mean_col <- find_col(
  RANK_RAW,
  c(
    "MeanAbsNSC",
    "mean_abs_nsc",
    "mean_abs_NSC",
    "mean_abs",
    "Mean_abs_NSC"
  )
)

rank_max_col <- find_col(
  RANK_RAW,
  c(
    "MaxAbsNSC",
    "max_abs_nsc",
    "max_abs_NSC",
    "max_abs",
    "MaximumAbsNSC"
  )
)


TABLE_S5 <- RANK_RAW %>%

  transmute(
    Parameter =
      as.character(
        .data[[rank_param_col]]
      ),

    `Mean |NSC|` =
      as.numeric(
        .data[[rank_mean_col]]
      ),

    `Maximum |NSC|` =
      as.numeric(
        .data[[rank_max_col]]
      )
  ) %>%

  arrange(
    desc(
      `Mean |NSC|`
    )
  ) %>%

  mutate(
    Rank =
      row_number(),

    Prioritized =
      Parameter %in%
        PRIORITIZED_PARAMS
  ) %>%

  select(
    Rank,
    Parameter,
    `Mean |NSC|`,
    `Maximum |NSC|`,
    Prioritized
  )


write_csv(
  TABLE_S5,
  file.path(
    TABLE_DIR,
    "Table_S5_NSC_parameter_ranking.csv"
  )
)


print(
  TABLE_S5 %>%
    slice_head(
      n = 10
    )
)


## ============================================================
## 7) TABLE S6
## FINAL MCMC DIAGNOSTICS
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table S6: MCMC diagnostics\n"
)

cat(
  "============================================================\n"
)


load_mcmc_diag_id <- function(id) {

  id_dir <- file.path(
    MCMC_DIR,
    paste0(
      "ID",
      id
    )
  )


  acc_file <- file.path(
    id_dir,
    "acceptance.csv"
  )

  ess_file <- file.path(
    id_dir,
    "ess.csv"
  )

  gel_file <- file.path(
    id_dir,
    "gelman.csv"
  )


  if (
    !file.exists(acc_file) ||
      !file.exists(ess_file)
  ) {

    stop(
      "Missing MCMC diagnostic output for ID ",
      id
    )
  }


  ACC <- read_csv(
    acc_file,
    show_col_types = FALSE
  )

  ESS <- read_csv(
    ess_file,
    show_col_types = FALSE
  )

  GEL <- if (
    file.exists(gel_file)
  ) {

    read_csv(
      gel_file,
      show_col_types = FALSE
    )

  } else {

    NULL
  }


  ## ----------------------------------------------------------
  ## Acceptance
  ## ----------------------------------------------------------

  acc_col <- find_col(
    ACC,
    c(
      "acceptance_rate",
      "AcceptanceRate",
      "Acceptance",
      "acceptance",
      "rate"
    ),
    required = FALSE
  )


  if (is.na(acc_col)) {

    acc_col <- first_numeric_column(
      ACC,
      exclude = c(
        "Chain",
        "chain"
      )
    )
  }


  acc <- as.numeric(
    ACC[[acc_col]]
  )


  if (
    max(
      acc,
      na.rm = TRUE
    ) > 1
  ) {

    acc <- acc / 100
  }


  ## ----------------------------------------------------------
  ## ESS
  ## ----------------------------------------------------------

  ess_col <- find_col(
    ESS,
    c(
      "ESS",
      "ess",
      "EffectiveSampleSize",
      "effective_sample_size"
    ),
    required = FALSE
  )


  if (is.na(ess_col)) {
    ess_col <- first_numeric_column(
      ESS
    )
  }


  ess <- as.numeric(
    ESS[[ess_col]]
  )


  ## ----------------------------------------------------------
  ## Gelman-Rubin
  ## ----------------------------------------------------------

  max_rhat <- NA_real_
  max_rhat_prioritized <- NA_real_


  if (!is.null(GEL)) {

    rhat_col <- find_col(
      GEL,
      c(
        "Rhat",
        "R_hat",
        "rhat",
        "Point.est.",
        "PointEstimate",
        "PSRF"
      ),
      required = FALSE
    )


    if (is.na(rhat_col)) {
      rhat_col <- first_numeric_column(
        GEL
      )
    }


    rhat <- as.numeric(
      GEL[[rhat_col]]
    )


    if (
      any(
        is.finite(rhat)
      )
    ) {

      max_rhat <- max(
        rhat,
        na.rm = TRUE
      )
    }


    param_col <- find_col(
      GEL,
      c(
        "Parameter",
        "parameter",
        "param"
      ),
      required = FALSE
    )


    if (!is.na(param_col)) {

      keep <- as.character(
        GEL[[param_col]]
      ) %in%
        PRIORITIZED_PARAMS


      if (
        any(keep) &&
          any(
            is.finite(
              rhat[keep]
            )
          )
      ) {

        max_rhat_prioritized <- max(
          rhat[keep],
          na.rm = TRUE
        )
      }
    }
  }


  ## ----------------------------------------------------------
  ## Prioritized ESS
  ## ----------------------------------------------------------

  min_ess_prioritized <- NA_real_


  ess_param_col <- find_col(
    ESS,
    c(
      "Parameter",
      "parameter",
      "param"
    ),
    required = FALSE
  )


  if (!is.na(ess_param_col)) {

    keep <- as.character(
      ESS[[ess_param_col]]
    ) %in%
      PRIORITIZED_PARAMS


    if (any(keep)) {

      min_ess_prioritized <- min(
        ess[keep],
        na.rm = TRUE
      )
    }
  }


  tibble(
    Dataset =
      paste0(
        "ID ",
        id
      ),

    ID =
      id,

    Particle_size =
      size_label(id),

    Acceptance_min =
      min(
        acc,
        na.rm = TRUE
      ),

    Acceptance_max =
      max(
        acc,
        na.rm = TRUE
      ),

    Mean_acceptance =
      mean(
        acc,
        na.rm = TRUE
      ),

    Minimum_ESS =
      min(
        ess,
        na.rm = TRUE
      ),

    Median_ESS =
      median(
        ess,
        na.rm = TRUE
      ),

    Maximum_Rhat =
      max_rhat,

    Minimum_ESS_prioritized =
      min_ess_prioritized,

    Maximum_Rhat_prioritized =
      max_rhat_prioritized
  )
}


TABLE_S6_MACHINE <- bind_rows(
  lapply(
    DEV_IDS,
    load_mcmc_diag_id
  )
)


TABLE_S6 <- TABLE_S6_MACHINE %>%

  transmute(
    Dataset,

    `Particle size` =
      Particle_size,

    `Acceptance rate across chains` =
      paste0(
        sprintf(
          "%.3f",
          Acceptance_min
        ),
        "–",
        sprintf(
          "%.3f",
          Acceptance_max
        )
      ),

    `Mean acceptance rate` =
      round(
        Mean_acceptance,
        3
      ),

    `Minimum ESS` =
      round(
        Minimum_ESS,
        1
      ),

    `Median ESS` =
      round(
        Median_ESS,
        1
      ),

    `Maximum Rhat` =
      round(
        Maximum_Rhat,
        3
      )
  )


write_csv(
  TABLE_S6,
  file.path(
    TABLE_DIR,
    "Table_S6_MCMC_diagnostics.csv"
  )
)


write_csv(
  TABLE_S6_MACHINE %>%
    select(
      Dataset,
      Particle_size,
      Minimum_ESS_prioritized,
      Maximum_Rhat_prioritized
    ),
  file.path(
    TABLE_DIR,
    "Table_S6_internal_prioritized_diagnostics.csv"
  )
)


print(
  TABLE_S6
)


## ============================================================
## 8) TABLE S7
## POSTERIOR ESTIMATES
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table S7: posterior parameter summaries\n"
)

cat(
  "============================================================\n"
)


load_posterior_summary <- function(id) {

  f <- file.path(
    MCMC_DIR,
    paste0(
      "ID",
      id
    ),
    "posterior_parameter_summary.csv"
  )


  if (!file.exists(f)) {
    stop(
      "Missing: ",
      f
    )
  }


  x <- read_csv(
    f,
    show_col_types = FALSE
  )


  pcol <- find_col(
    x,
    c(
      "Parameter",
      "parameter",
      "param"
    )
  )

  medcol <- find_col(
    x,
    c(
      "median",
      "Median",
      "q50",
      "q050"
    )
  )

  locol <- find_col(
    x,
    c(
      "q025",
      "Q025",
      "lower",
      "Lower"
    )
  )

  hicol <- find_col(
    x,
    c(
      "q975",
      "Q975",
      "upper",
      "Upper"
    )
  )


  tibble(
    ID =
      id,

    Parameter =
      as.character(
        x[[pcol]]
      ),

    Median =
      as.numeric(
        x[[medcol]]
      ),

    Q025 =
      as.numeric(
        x[[locol]]
      ),

    Q975 =
      as.numeric(
        x[[hicol]]
      )
  )
}


POST_ALL <- bind_rows(
  lapply(
    DEV_IDS,
    load_posterior_summary
  )
)


TABLE_S7_LONG <- POST_ALL %>%

  filter(
    Parameter %in%
      META_PARAM_ORDER
  ) %>%

  mutate(
    Parameter = factor(
      Parameter,
      levels =
        META_PARAM_ORDER
    ),

    Unit = ifelse(
      Parameter ==
        "KurineC",
      "L h^-1 kg^-0.75",
      "h^-1"
    ),

    Particle_size =
      size_label(ID),

    Estimate =
      format_ci(
        Median,
        Q025,
        Q975
      )
  ) %>%

  arrange(
    Parameter,
    ID
  )


TABLE_S7 <- TABLE_S7_LONG %>%

  transmute(
    Parameter =
      as.character(
        Parameter
      ),

    Unit,

    ID =
      paste0(
        "ID ",
        ID,
        " (",
        Particle_size,
        ")"
      ),

    Estimate
  ) %>%

  pivot_wider(
    names_from =
      ID,

    values_from =
      Estimate
  )


write_csv(
  TABLE_S7_LONG,
  file.path(
    TABLE_DIR,
    "Table_S7_posterior_parameters_long.csv"
  )
)


write_csv(
  TABLE_S7,
  file.path(
    TABLE_DIR,
    "Table_S7_posterior_parameters.csv"
  )
)


print(
  TABLE_S7
)


## ============================================================
## 9) TABLE S8
## DEVELOPMENT-DATASET 95% INTERVAL COVERAGE
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table S8: development interval coverage\n"
)

cat(
  "============================================================\n"
)


DATA_FILE <- first_existing(
  c(
    file.path(
      "data",
      "data.csv"
    ),
    "data.csv"
  )
)


if (is.null(DATA_FILE)) {
  stop(
    "Could not find data/data.csv."
  )
}


OBS_DATA <- read.csv(
  DATA_FILE,
  encoding = "latin1",
  check.names = FALSE
) %>%

  mutate(
    ID =
      as.numeric(ID),

    Tissue =
      harmonize_tissue(
        Tissue
      ),

    Time_h =
      as.numeric(
        Time_h
      ),

    Conc =
      as.numeric(
        Conc
      )
  ) %>%

  filter(
    ID %in%
      DEV_IDS,

    Tissue %in%
      c(
        "Blood",
        "Liver",
        "Kidney",
        "Spleen",
        "Lungs",
        "GI",
        "Urine"
      )
  )


load_ppc_for_coverage <- function(id) {

  f <- file.path(
    MCMC_DIR,
    paste0(
      "ID",
      id
    ),
    "posterior_prediction_summary.csv"
  )


  if (!file.exists(f)) {
    stop(
      "Missing: ",
      f
    )
  }


  x <- read_csv(
    f,
    show_col_types = FALSE
  )


  time_col <- find_col(
    x,
    c(
      "time",
      "Time",
      "time_h",
      "Time_h"
    )
  )

  tissue_col <- find_col(
    x,
    c(
      "Tissue",
      "tissue"
    )
  )

  low_col <- find_col(
    x,
    c(
      "pred_lo",
      "q025",
      "lower",
      "Lower"
    )
  )

  high_col <- find_col(
    x,
    c(
      "pred_hi",
      "q975",
      "upper",
      "Upper"
    )
  )


  tibble(
    ID =
      id,

    time =
      as.numeric(
        x[[time_col]]
      ),

    Tissue =
      harmonize_tissue(
        x[[tissue_col]]
      ),

    Lower =
      as.numeric(
        x[[low_col]]
      ),

    Upper =
      as.numeric(
        x[[high_col]]
      )
  )
}


PPC_COVERAGE <- bind_rows(
  lapply(
    DEV_IDS,
    load_ppc_for_coverage
  )
)


interpolate_interval <- function(
    id,
    tissue,
    time
) {

  x <- PPC_COVERAGE %>%

    filter(
      ID == id,
      Tissue == tissue,
      is.finite(.data$time),
      is.finite(Lower),
      is.finite(Upper)
    ) %>%

    arrange(
      .data$time
    )


  if (nrow(x) == 0) {

    return(
      c(
        Lower = NA_real_,
        Upper = NA_real_
      )
    )
  }


  if (nrow(x) == 1) {

    return(
      c(
        Lower = x$Lower[1],
        Upper = x$Upper[1]
      )
    )
  }


  lower <- approx(
    x =
      x$time,

    y =
      x$Lower,

    xout =
      time,

    rule =
      2,

    ties =
      mean
  )$y


  upper <- approx(
    x =
      x$time,

    y =
      x$Upper,

    xout =
      time,

    rule =
      2,

    ties =
      mean
  )$y


  c(
    Lower =
      lower,
    Upper =
      upper
  )
}


coverage_intervals <- t(
  mapply(
    FUN =
      interpolate_interval,

    id =
      OBS_DATA$ID,

    tissue =
      OBS_DATA$Tissue,

    time =
      OBS_DATA$Time_h
  )
)


OBS_COVERAGE <- OBS_DATA %>%

  mutate(
    Lower =
      coverage_intervals[
        ,
        "Lower"
      ],

    Upper =
      coverage_intervals[
        ,
        "Upper"
      ],

    Covered =
      is.finite(Conc) &
      is.finite(Lower) &
      is.finite(Upper) &
      Conc >= Lower &
      Conc <= Upper
  )


TABLE_S8_BY_ID <- OBS_COVERAGE %>%

  group_by(
    ID
  ) %>%

  summarise(
    `Particle size (nm)` =
      as.numeric(
        ID_SIZE_MAP[
          as.character(
            first(ID)
          )
        ]
      ),

    Observations =
      sum(
        is.finite(Conc) &
          is.finite(Lower) &
          is.finite(Upper)
      ),

    `Within 95% interval` =
      sum(
        Covered,
        na.rm = TRUE
      ),

    Coverage =
      100 *
      `Within 95% interval` /
      Observations,

    .groups =
      "drop"
  ) %>%

  mutate(
    ID =
      as.character(ID)
  )


TABLE_S8_OVERALL <- tibble(
  ID =
    "Overall",

  `Particle size (nm)` =
    NA_real_,

  Observations =
    sum(
      is.finite(
        OBS_COVERAGE$Conc
      ) &
        is.finite(
          OBS_COVERAGE$Lower
        ) &
        is.finite(
          OBS_COVERAGE$Upper
        )
    ),

  `Within 95% interval` =
    sum(
      OBS_COVERAGE$Covered,
      na.rm = TRUE
    ),

  Coverage =
    100 *
    `Within 95% interval` /
    Observations
)


TABLE_S8 <- bind_rows(
  TABLE_S8_BY_ID,
  TABLE_S8_OVERALL
)


write_csv(
  TABLE_S8,
  file.path(
    TABLE_DIR,
    "Table_S8_development_95PI_coverage.csv"
  )
)


write_csv(
  OBS_COVERAGE,
  file.path(
    TABLE_DIR,
    "Table_S8_observation_level_coverage_source.csv"
  )
)


print(
  TABLE_S8
)


## ============================================================
## 10) TABLE S9
## RANDOM-EFFECTS META-REGRESSION
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table S9: random-effects meta-regression\n"
)

cat(
  "============================================================\n"
)


meta_csv <- file.path(
  META_DIR,
  "Table_S9_meta_regression.csv"
)


if (!file.exists(meta_csv)) {
  stop(
    "Missing: ",
    meta_csv
  )
}


META_RAW <- read_csv(
  meta_csv,
  show_col_types = FALSE
)


names(META_RAW) <- trimws(
  names(META_RAW)
)


cat(
  "Meta-regression source columns:\n"
)

print(
  names(META_RAW)
)


meta_param_col <- find_col(
  META_RAW,
  c(
    "Parameter",
    "parameter",
    "param"
  )
)


meta_n_col <- find_col(
  META_RAW,
  c(
    "N",
    "n"
  )
)


## ------------------------------------------------------------
## IMPORTANT:
## Supports the actual current Script 04 header:
## "Size slope, beta1"
## ------------------------------------------------------------

meta_slope_col <- find_col(
  META_RAW,
  c(
    "Size slope, beta1",
    "beta1",
    "Size slope",
    "coef1",
    "slope",
    "Size_slope"
  )
)


meta_se_col <- find_col(
  META_RAW,
  c(
    "SE",
    "se",
    "se_beta1",
    "beta1_se",
    "slope_se"
  )
)


meta_lo_col <- find_col(
  META_RAW,
  c(
    "95% CI lower",
    "beta1_ci_low",
    "ci_beta1_low",
    "ci_low",
    "CI_lower",
    "lower"
  )
)


meta_hi_col <- find_col(
  META_RAW,
  c(
    "95% CI upper",
    "beta1_ci_high",
    "ci_beta1_high",
    "ci_high",
    "CI_upper",
    "upper"
  )
)


meta_p_col <- find_col(
  META_RAW,
  c(
    "p-value",
    "p_beta1",
    "beta1_p",
    "p_value",
    "pvalue",
    "p"
  )
)


meta_tau_col <- find_col(
  META_RAW,
  c(
    "tau2",
    "tau_sq",
    "tau_squared"
  )
)


TABLE_S9 <- META_RAW %>%

  transmute(
    Parameter =
      as.character(
        .data[[meta_param_col]]
      ),

    N =
      as.numeric(
        .data[[meta_n_col]]
      ),

    `Size slope, beta1` =
      as.numeric(
        .data[[meta_slope_col]]
      ),

    SE =
      as.numeric(
        .data[[meta_se_col]]
      ),

    `95% CI lower` =
      as.numeric(
        .data[[meta_lo_col]]
      ),

    `95% CI upper` =
      as.numeric(
        .data[[meta_hi_col]]
      ),

    `p-value` =
      as.numeric(
        .data[[meta_p_col]]
      ),

    tau2 =
      as.numeric(
        .data[[meta_tau_col]]
      )
  ) %>%

  mutate(
    Parameter =
      factor(
        Parameter,
        levels =
          META_PARAM_ORDER
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


write_csv(
  TABLE_S9,
  file.path(
    TABLE_DIR,
    "Table_S9_meta_regression.csv"
  )
)


print(
  TABLE_S9
)


## ============================================================
## 11) TABLE S10
## EXTERNAL EVALUATION
##
## Use the authoritative four-row table produced directly
## by Script 05.
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table S10: external evaluation\n"
)

cat(
  "============================================================\n"
)

s10_source <- file.path(
  EXT_DIR,
  "Table_S10_external_evaluation.csv"
)

if (!file.exists(s10_source)) {

  stop(
    "Missing authoritative Table S10 source: ",
    s10_source
  )
}

TABLE_S10_RAW <- read_csv(
  s10_source,
  show_col_types = FALSE
)

TABLE_S10 <- TABLE_S10_RAW %>%

  transmute(
    Group = Group,

    `Size (nm)` = `Size (nm)`,

    `N (all)` = n_all,

    `N (positive)` = n_positive,

    R2 = R2,

    `RMSE (%ID)` = RMSE,

    `95% uncertainty covered (n)` =
      coverage_n,

    `95% uncertainty coverage (%)` =
      coverage_pct
  )

write_csv(
  TABLE_S10,
  file.path(
    TABLE_DIR,
    "Table_S10_external_evaluation.csv"
  )
)

print(
  TABLE_S10,
  n = Inf,
  width = Inf
)


## ============================================================
## 12) MAIN MANUSCRIPT TABLE 2
## DETERMINISTIC FITTED PARAMETERS
##
## Search saved outputs only.
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Table 2: deterministic fitted parameter estimates\n"
)

cat(
  "============================================================\n"
)


find_parameter_df_recursive <- function(x) {

  if (is.data.frame(x)) {

    nm <- names(x)

    has_id <- any(
      c(
        "ID",
        "id"
      ) %in% nm
    )

    has_param <- any(
      c(
        "Parameter",
        "parameter",
        "param"
      ) %in% nm
    )

    has_value <- any(
      c(
        "Estimate",
        "estimate",
        "Value",
        "value",
        "Fitted",
        "fitted",
        "fitted_value"
      ) %in% nm
    )


    if (
      has_id &&
        has_param &&
        has_value
    ) {
      return(x)
    }
  }


  if (is.list(x)) {

    for (i in seq_along(x)) {

      ans <- find_parameter_df_recursive(
        x[[i]]
      )

      if (!is.null(ans)) {
        return(ans)
      }
    }
  }


  NULL
}


DET_PARAM_DF <- NULL


## ------------------------------------------------------------
## Search CSVs
## ------------------------------------------------------------

det_csvs <- list.files(
  DET_DIR,
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)


for (f in det_csvs) {

  x <- tryCatch(
    read_csv(
      f,
      show_col_types = FALSE
    ),
    error = function(e) NULL
  )


  if (is.null(x)) {
    next
  }


  ans <- find_parameter_df_recursive(
    x
  )


  if (!is.null(ans)) {

    DET_PARAM_DF <- ans

    message(
      "Table 2 parameter source: ",
      f
    )

    break
  }
}


## ------------------------------------------------------------
## Search RDS if CSV search failed
## ------------------------------------------------------------

if (is.null(DET_PARAM_DF)) {

  det_rds <- first_existing(
    c(
      file.path(
        DET_DIR,
        "deterministic_outputs.rds"
      ),
      file.path(
        DET_DIR,
        "fitted_objects.rds"
      )
    )
  )


  if (!is.null(det_rds)) {

    obj <- readRDS(
      det_rds
    )

    DET_PARAM_DF <- find_parameter_df_recursive(
      obj
    )


    if (!is.null(DET_PARAM_DF)) {

      message(
        "Table 2 parameter source: ",
        det_rds
      )
    }
  }
}


if (!is.null(DET_PARAM_DF)) {

  p_id_col <- find_col(
    DET_PARAM_DF,
    c(
      "ID",
      "id"
    )
  )

  p_param_col <- find_col(
    DET_PARAM_DF,
    c(
      "Parameter",
      "parameter",
      "param"
    )
  )

  p_value_col <- find_col(
    DET_PARAM_DF,
    c(
      "Estimate",
      "estimate",
      "Value",
      "value",
      "Fitted",
      "fitted",
      "fitted_value"
    )
  )


  TABLE2_LONG <- DET_PARAM_DF %>%

    transmute(
      ID =
        as.numeric(
          .data[[p_id_col]]
        ),

      Parameter =
        as.character(
          .data[[p_param_col]]
        ),

      Estimate =
        as.numeric(
          .data[[p_value_col]]
        )
    ) %>%

    filter(
      ID %in%
        DEV_IDS,

      Parameter %in%
        DETERMINISTIC_PARAMS
    ) %>%

    arrange(
      match(
        Parameter,
        DETERMINISTIC_PARAMS
      ),
      ID
    )


  TABLE2 <- TABLE2_LONG %>%

    mutate(
      Dataset =
        paste0(
          "ID ",
          ID,
          " (",
          size_label(ID),
          ")"
        )
    ) %>%

    select(
      Parameter,
      Dataset,
      Estimate
    ) %>%

    pivot_wider(
      names_from =
        Dataset,

      values_from =
        Estimate
    )


  write_csv(
    TABLE2_LONG,
    file.path(
      TABLE_DIR,
      "Table_2_deterministic_parameters_long.csv"
    )
  )


  write_csv(
    TABLE2,
    file.path(
      TABLE_DIR,
      "Table_2_deterministic_parameters.csv"
    )
  )


  print(
    TABLE2
  )

} else {

  warning(
    paste0(
      "Could not automatically identify a saved deterministic ",
      "parameter table for manuscript Table 2. ",
      "Tables S3–S10 are unaffected."
    )
  )
}


## ============================================================
## 13) FINAL VALIDATION SUMMARY
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "TABLE GENERATION COMPLETED\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nGenerated analysis-derived tables:\n"
)


generated_files <- list.files(
  TABLE_DIR,
  pattern = "\\.csv$",
  full.names = FALSE
)


for (f in generated_files) {

  cat(
    "  - ",
    f,
    "\n",
    sep = ""
  )
}


cat(
  "\nValidation checks:\n"
)


cat(
  "  Deterministic observations: ",
  nrow(
    DET_GOF
  ),
  "\n",
  sep = ""
)


cat(
  "  Detailed NSC rows: ",
  nrow(
    TABLE_S4_LONG
  ),
  "\n",
  sep = ""
)


cat(
  "  Ranked NSC parameters: ",
  nrow(
    TABLE_S5
  ),
  "\n",
  sep = ""
)


cat(
  "  MCMC development datasets: ",
  nrow(
    TABLE_S6
  ),
  "\n",
  sep = ""
)


cat(
  "  Prioritized posterior parameters: ",
  length(
    unique(
      TABLE_S7_LONG$Parameter
    )
  ),
  "\n",
  sep = ""
)


cat(
  "  Development coverage observations: ",
  TABLE_S8$Observations[
    TABLE_S8$ID ==
      "Overall"
  ],
  "\n",
  sep = ""
)


cat(
  "  Development observations within 95% interval: ",
  TABLE_S8$`Within 95% interval`[
    TABLE_S8$ID ==
      "Overall"
  ],
  "\n",
  sep = ""
)


cat(
  "  Overall development coverage (%): ",
  round(
    TABLE_S8$Coverage[
      TABLE_S8$ID ==
        "Overall"
    ],
    1
  ),
  "\n",
  sep = ""
)


cat(
  "  Meta-regression parameters: ",
  nrow(
    TABLE_S9
  ),
  "\n",
  sep = ""
)


cat(
  "  External evaluation table rows: ",
  nrow(
    TABLE_S10
  ),
  "\n",
  sep = ""
)


cat(
  "\nNo model fitting or simulation was rerun.\n"
)


cat(
  "Tables saved to: ",
  normalizePath(
    TABLE_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat(
  "============================================================\n"
)