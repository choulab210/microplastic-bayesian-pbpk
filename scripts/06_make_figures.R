## ============================================================
## 06_make_figures.R
##
## Reproduce manuscript figures while reading FINAL saved outputs from Scripts 01–05.
##
## IMPORTANT
## - Does NOT rerun deterministic fitting
## - Does NOT rerun sensitivity analysis
## - Does NOT rerun MCMC
## - Does NOT rerun meta-regression
## - Does NOT rerun external PBPK simulations
##
## Main outputs:
##   Figure 2: Deterministic PBPK GOF
##   Figure 3: Combined NSC heatmap by particle size
##   Figure 4: Bayesian posterior time-course predictions
##   Figure 5: Random-effects meta-regression
##   Figure 6: External evaluation
##
## Supplement:
##   Figure S6: Posterior predictive time course for ID 5
## ============================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggExtra)
  library(gridExtra)
  library(grid)
  library(scales)
  library(pheatmap)
  library(readr)
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

FIG_MAIN <- file.path("figures", "main")
FIG_SUPP <- file.path("figures", "supplement")

dir.create(
  FIG_MAIN,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_SUPP,
  recursive = TRUE,
  showWarnings = FALSE
)

FIG_DPI <- 600


## ============================================================
## 1) SHARED HELPERS
## ============================================================

find_col <- function(
    df,
    candidates,
    required = TRUE
) {

  hit <- candidates[
    candidates %in% names(df)
  ]

  if (length(hit) > 0) {
    return(hit[1])
  }

  if (required) {

    stop(
      "Could not identify required column. Tried: ",
      paste(
        candidates,
        collapse = ", "
      ),
      "\nAvailable columns: ",
      paste(
        names(df),
        collapse = ", "
      )
    )
  }

  NA_character_
}


first_existing <- function(paths) {

  hit <- paths[
    file.exists(paths)
  ]

  if (length(hit) == 0) {
    return(NULL)
  }

  hit[1]
}


save_tiff <- function(
    plot,
    path,
    width,
    height
) {

  ggsave(
    filename = path,
    plot = plot,
    device = "tiff",
    width = width,
    height = height,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw"
  )

  message(
    "Saved: ",
    path
  )
}


save_gtable_tiff <- function(
    gtable,
    path,
    width,
    height
) {

  ggsave(
    filename = path,
    plot = gtable,
    device = "tiff",
    width = width,
    height = height,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw"
  )

  message(
    "Saved: ",
    path
  )
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
    "Lung" = "Lungs",
    "Gi" = "GI",
    .default = out
  )
}


ID_SIZE_MAP <- c(
  "1" = 20,
  "2" = 220,
  "3" = 1000,
  "4" = 6000,
  "5" = 220
)


size_label_fun <- function(size_nm) {

  dplyr::case_when(
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


## ============================================================
## 2) FIGURE 2
## ORIGINAL DETERMINISTIC GOF ARTWORK
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Figure 2: original GOF artwork\n"
)

cat(
  "============================================================\n"
)


theme_gof_poster <- function(
    show_legend = TRUE
) {

  theme_bw(
    base_size = 13
  ) +

    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 2
      ),

      axis.ticks = element_line(
        linewidth = 1.5
      ),

      axis.title = element_text(
        size = 18,
        colour = "black",
        face = "bold"
      ),

      axis.text = element_text(
        size = 16,
        colour = "black",
        face = "bold"
      ),

      legend.title = element_text(
        size = 16,
        face = "bold"
      ),

      legend.text = element_text(
        size = 14,
        face = "bold"
      ),

      legend.position = if (
        show_legend
      ) {
        "right"
      } else {
        "none"
      }
    )
}


log10_formatter <- trans_format(
  "log10",
  math_format(10^.x)
)


calc_gof_metrics_old <- function(df) {

  df <- df %>%
    dplyr::filter(
      is.finite(obs),
      is.finite(pred),
      obs > 0,
      pred > 0
    )

  if (nrow(df) == 0) {

    return(
      data.frame(
        R2 = NA_real_,
        RMSE = NA_real_,
        MAPE = NA_real_,
        n = 0
      )
    )
  }


  r2 <- if (
    length(unique(df$obs)) > 1 &&
      length(unique(df$pred)) > 1
  ) {

    cor(
      df$obs,
      df$pred,
      use = "complete.obs"
    )^2

  } else {

    NA_real_
  }


  rmse <- sqrt(
    mean(
      (
        df$pred -
          df$obs
      )^2,
      na.rm = TRUE
    )
  )


  mape <- mean(
    abs(
      df$pred -
        df$obs
    ) /
      pmax(
        df$obs,
        1e-8
      ),
    na.rm = TRUE
  ) * 100


  data.frame(
    R2 = r2,
    RMSE = rmse,
    MAPE = mape,
    n = nrow(df)
  )
}


## ------------------------------------------------------------
## Find deterministic GOF dataframe from saved RDS
## ------------------------------------------------------------

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

    for (
      i in seq_along(x)
    ) {

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

  rds_candidates <- c(
    file.path(
      DET_DIR,
      "deterministic_outputs.rds"
    ),
    file.path(
      DET_DIR,
      "fitted_objects.rds"
    )
  )


  rds_file <- first_existing(
    rds_candidates
  )


  if (!is.null(rds_file)) {

    obj <- readRDS(
      rds_file
    )

    df <- find_gof_df_recursive(
      obj
    )

    if (!is.null(df)) {

      message(
        "Figure 2 GOF source: ",
        rds_file
      )

      return(df)
    }
  }


  csvs <- list.files(
    DET_DIR,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )


  for (
    f in csvs
  ) {

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
          "Obs",
          "Observed",
          "observed"
        ) %in% nm
      ) &&
        any(
          c(
            "pred",
            "Pred",
            "Predicted",
            "predicted"
          ) %in% nm
        ) &&
        any(
          c(
            "ID",
            "id"
          ) %in% nm
        ) &&
        any(
          c(
            "Tissue",
            "tissue"
          ) %in% nm
        )
    ) {

      message(
        "Figure 2 GOF source: ",
        f
      )

      return(x)
    }
  }


  stop(
    "Could not locate deterministic GOF data in ",
    DET_DIR,
    ".\nExpected saved data with ID, Tissue, obs, and pred."
  )
}


gof_raw <- load_deterministic_gof()


obs_col <- find_col(
  gof_raw,
  c(
    "obs",
    "Obs",
    "Observed",
    "observed"
  )
)


pred_col <- find_col(
  gof_raw,
  c(
    "pred",
    "Pred",
    "Predicted",
    "predicted"
  )
)


id_col <- find_col(
  gof_raw,
  c(
    "ID",
    "id"
  )
)


tissue_col <- find_col(
  gof_raw,
  c(
    "Tissue",
    "tissue"
  )
)


gof_all <- gof_raw %>%

  transmute(
    ID = as.numeric(
      .data[[id_col]]
    ),

    Tissue = harmonize_tissue(
      .data[[tissue_col]]
    ),

    obs = as.numeric(
      .data[[obs_col]]
    ),

    pred = as.numeric(
      .data[[pred_col]]
    )
  ) %>%

  filter(
    ID %in% 1:5
  )


plot_overall_gof <- function(
    gof_all,
    log_scale = TRUE,
    show_legend = FALSE
) {

  gof_all <- gof_all %>%

    dplyr::filter(
      is.finite(obs),
      is.finite(pred),
      obs > 0,
      pred > 0
    )


  met <- calc_gof_metrics_old(
    gof_all
  )


  xmin <- min(
    c(
      gof_all$obs,
      gof_all$pred
    ),
    na.rm = TRUE
  ) / 1.5


  xmax <- max(
    c(
      gof_all$obs,
      gof_all$pred
    ),
    na.rm = TRUE
  ) * 1.5


  if (
    !is.finite(xmin) ||
      xmin <= 0
  ) {

    xmin <- 1e-6
  }


  if (
    !is.finite(xmax) ||
      xmax <= xmin
  ) {

    xmax <- xmin * 10
  }


  p <- ggplot(
    gof_all,
    aes(
      x = obs,
      y = pred,
      color = Tissue,
      shape = factor(ID)
    )
  ) +

    geom_point(
      size = 2.8,
      alpha = 0.80
    ) +

    ## Identity line
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +

    ## 2-fold upper line
    geom_abline(
      slope = 2,
      intercept = 0,
      linetype = "dotted",
      color = "gray40",
      linewidth = 0.7
    ) +

    ## 2-fold lower line
    geom_abline(
      slope = 0.5,
      intercept = 0,
      linetype = "dotted",
      color = "gray40",
      linewidth = 0.7
    ) +

    labs(
      x = "Observed %ID",
      y = "Predicted %ID",
      color = "Tissue",
      shape = "ID"
    ) +

    theme_gof_poster(
      show_legend = show_legend
    )


  if (log_scale) {

    p <- p +

      scale_x_log10(
        limits = c(
          xmin,
          xmax
        ),
        labels = log10_formatter
      ) +

      scale_y_log10(
        limits = c(
          xmin,
          xmax
        ),
        labels = log10_formatter
      )

  } else {

    p <- p +

      coord_cartesian(
        xlim = c(
          xmin,
          xmax
        ),
        ylim = c(
          xmin,
          xmax
        )
      )
  }


  invisible(
    list(
      plot = p,
      met = met
    )
  )
}


plot_ratio_gof <- function(
    gof_all,
    log_scale = TRUE,
    show_legend = FALSE
) {

  df <- gof_all %>%

    dplyr::filter(
      is.finite(obs),
      is.finite(pred),
      obs > 0,
      pred > 0
    ) %>%

    dplyr::mutate(
      ratio = pred / obs
    )


  p <- ggplot(
    df,
    aes(
      x = pred,
      y = ratio,
      color = Tissue,
      shape = factor(ID)
    )
  ) +

    geom_point(
      size = 2.8,
      alpha = 0.80
    ) +

    geom_hline(
      yintercept = 2,
      linetype = "dotted",
      color = "red3",
      linewidth = 0.9
    ) +

    geom_hline(
      yintercept = 0.5,
      linetype = "dotted",
      color = "red3",
      linewidth = 0.9
    ) +

    geom_hline(
      yintercept = 3,
      linetype = "dashed",
      color = "black",
      linewidth = 0.9
    ) +

    geom_hline(
      yintercept = 1 / 3,
      linetype = "dashed",
      color = "black",
      linewidth = 0.9
    ) +

    labs(
      x = "Predicted %ID",
      y = "Predicted/Observed",
      color = "Tissue",
      shape = "ID"
    ) +

    theme_gof_poster(
      show_legend = show_legend
    )


  if (log_scale) {

    p <- p +

      scale_x_log10(
        labels = log10_formatter
      ) +

      scale_y_log10(
        labels = log10_formatter
      )
  }


  p
}


p2A <- plot_overall_gof(
  gof_all,
  log_scale = TRUE,
  show_legend = FALSE
)$plot


p2B0 <- plot_ratio_gof(
  gof_all,
  log_scale = TRUE,
  show_legend = FALSE
)


p2B <- ggMarginal(
  p2B0,
  type = "histogram",
  margins = "y",
  yparams = list(
    binwidth = 0.1,
    fill = "#FFA488",
    color = "black"
  )
)


FIG2 <- gridExtra::arrangeGrob(
  p2A,
  p2B,
  nrow = 1
)


save_gtable_tiff(
  FIG2,
  file.path(
    FIG_MAIN,
    "Figure_2_PBPK_GOF_plot.tiff"
  ),
  width = 12,
  height = 5
)


## ============================================================
## 3) FIGURE 3
## ORIGINAL COMBINED NSC PHEATMAP
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Figure 3: original combined NSC heatmap\n"
)

cat(
  "============================================================\n"
)


nsc_file <- first_existing(
  c(
    file.path(
      SENS_DIR,
      "NSC_amount_all.csv"
    ),
    file.path(
      SENS_DIR,
      "nsc_by_id_tissue_parameter.csv"
    )
  )
)


if (is.null(nsc_file)) {

  stop(
    "Missing NSC output in ",
    SENS_DIR
  )
}


NSC_amount_all <- read_csv(
  nsc_file,
  show_col_types = FALSE
)


required_nsc_cols <- c(
  "ID",
  "Parameter",
  "Tissue",
  "NSC"
)


missing_nsc_cols <- setdiff(
  required_nsc_cols,
  names(NSC_amount_all)
)


if (
  length(
    missing_nsc_cols
  ) > 0
) {

  stop(
    "NSC output is missing: ",
    paste(
      missing_nsc_cols,
      collapse = ", "
    )
  )
}


nsc_cutoff_present_combined <- 0.1
mask_weak_cells_combined <- FALSE
n_row_clusters_combined <- 4
heatmap_levels_combined <- 400


tissue_label_map <- c(
  "Blood" = "AUCB",
  "GI" = "AUCGI",
  "Spleen" = "AUCS",
  "Kidney" = "AUCK",
  "Liver" = "AUCL",
  "Lungs" = "AUCLu",
  "Urine" = "AUCU",
  "Rest" = "AUCRt"
)


NSC_amount_all2 <- NSC_amount_all %>%

  mutate(
    ID = as.numeric(ID),

    Tissue = harmonize_tissue(
      Tissue
    ),

    Size_nm = as.numeric(
      ID_SIZE_MAP[
        as.character(ID)
      ]
    ),

    absNSC = abs(NSC)
  ) %>%

  filter(
    !is.na(Size_nm)
  )


## ------------------------------------------------------------
## IDs 2 and 5 averaged only for 220-nm visualization,
## matching the original artwork.
## ------------------------------------------------------------

NSC_plot_df <- bind_rows(

  NSC_amount_all2 %>%

    filter(
      ID %in% c(
        1,
        3,
        4
      )
    ) %>%

    mutate(
      SizeGroup = Size_nm
    ),


  NSC_amount_all2 %>%

    filter(
      ID %in% c(
        2,
        5
      )
    ) %>%

    group_by(
      Size_nm,
      Parameter,
      Tissue
    ) %>%

    summarise(
      NSC = mean(
        NSC,
        na.rm = TRUE
      ),

      .groups = "drop"
    ) %>%

    mutate(
      ID = 25,
      SizeGroup = Size_nm
    )
) %>%

  mutate(
    TissueLabel = recode(
      Tissue,
      !!!tissue_label_map
    ),

    SizeLabel = size_label_fun(
      SizeGroup
    ),

    ColLabel = paste(
      SizeLabel,
      TissueLabel,
      sep = " | "
    )
  ) %>%

  filter(
    !is.na(TissueLabel)
  )


selected_pars_combined <- NSC_plot_df %>%

  group_by(
    Parameter
  ) %>%

  summarise(
    max_abs = max(
      abs(NSC),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  filter(
    is.finite(max_abs),
    max_abs >
      nsc_cutoff_present_combined
  ) %>%

  arrange(
    desc(max_abs)
  ) %>%

  pull(
    Parameter
  )


combined_df2 <- NSC_plot_df %>%

  filter(
    Parameter %in%
      selected_pars_combined
  )


mat_combined <- combined_df2 %>%

  select(
    Parameter,
    ColLabel,
    NSC
  ) %>%

  group_by(
    Parameter,
    ColLabel
  ) %>%

  summarise(
    NSC = mean(
      NSC,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  pivot_wider(
    names_from = ColLabel,
    values_from = NSC
  ) %>%

  as.data.frame()


rownames(
  mat_combined
) <- mat_combined$Parameter


mat_combined$Parameter <- NULL


mat_combined <- as.matrix(
  mat_combined
)


desired_cols_combined <- c(
  "20 nm | AUCB",
  "20 nm | AUCGI",
  "20 nm | AUCS",
  "20 nm | AUCK",
  "20 nm | AUCL",
  "20 nm | AUCLu",
  "20 nm | AUCU",
  "20 nm | AUCRt",

  "220 nm | AUCB",
  "220 nm | AUCGI",
  "220 nm | AUCS",
  "220 nm | AUCK",
  "220 nm | AUCL",
  "220 nm | AUCLu",
  "220 nm | AUCU",
  "220 nm | AUCRt",

  "1 µm | AUCB",
  "1 µm | AUCGI",
  "1 µm | AUCS",
  "1 µm | AUCK",
  "1 µm | AUCL",
  "1 µm | AUCLu",
  "1 µm | AUCU",
  "1 µm | AUCRt",

  "6 µm | AUCB",
  "6 µm | AUCGI",
  "6 µm | AUCS",
  "6 µm | AUCK",
  "6 µm | AUCL",
  "6 µm | AUCLu",
  "6 µm | AUCU",
  "6 µm | AUCRt"
)


desired_cols_combined <- desired_cols_combined[
  desired_cols_combined %in%
    colnames(mat_combined)
]


mat_combined <- mat_combined[
  ,
  desired_cols_combined,
  drop = FALSE
]


mat_plot_combined <- mat_combined


if (
  mask_weak_cells_combined
) {

  mat_plot_combined[
    abs(
      mat_plot_combined
    ) <=
      nsc_cutoff_present_combined
  ] <- NA_real_
}


col_ann <- data.frame(
  SizeGroup = sub(
    " \\|.*$",
    "",
    colnames(mat_combined)
  )
)


rownames(
  col_ann
) <- colnames(
  mat_combined
)


col_ann$SizeGroup <- factor(
  col_ann$SizeGroup,
  levels = c(
    "20 nm",
    "220 nm",
    "1 µm",
    "6 µm"
  )
)


col_ann_colors <- list(
  SizeGroup = c(
    "20 nm" = "#66C2A5",
    "220 nm" = "#FC8D62",
    "1 µm" = "#8DA0CB",
    "6 µm" = "#E78AC3"
  )
)


mat_for_clust_combined <- mat_combined


mat_for_clust_combined[
  !is.finite(
    mat_for_clust_combined
  )
] <- 0


hc_rows_combined <- hclust(
  dist(
    mat_for_clust_combined
  )
)


row_clusters_combined <- cutree(
  hc_rows_combined,
  k = min(
    n_row_clusters_combined,
    nrow(
      mat_for_clust_combined
    )
  )
)


row_ann <- data.frame(
  Cluster = factor(
    row_clusters_combined
  )
)


rownames(
  row_ann
) <- rownames(
  mat_for_clust_combined
)


row_ann_colors <- list(
  Cluster = c(
    "1" = "#9C6ADE",
    "2" = "#D9D9E8",
    "3" = "#F0E442",
    "4" = "#E6C300"
  )[
    levels(
      row_ann$Cluster
    )
  ]
)


global_lim_combined <- max(
  abs(mat_combined),
  na.rm = TRUE
)


if (
  !is.finite(
    global_lim_combined
  ) ||
    global_lim_combined == 0
) {

  global_lim_combined <- 1
}


p3 <- pheatmap(
  mat_plot_combined,

  cluster_rows = hc_rows_combined,

  cluster_cols = FALSE,

  annotation_row = row_ann,

  annotation_col = col_ann,

  annotation_colors = c(
    row_ann_colors,
    col_ann_colors
  ),

  color = colorRampPalette(
    c(
      "#D7191C",
      "#F7F7F7",
      "#2B83BA"
    )
  )(
    heatmap_levels_combined
  ),

  breaks = seq(
    -global_lim_combined,
    global_lim_combined,
    length.out =
      heatmap_levels_combined + 1
  ),

  border_color = "grey40",

  cellwidth = 18,

  cellheight = 14,

  fontsize = 9,

  fontsize_row = 10,

  fontsize_col = 8,

  angle_col = 90,

  main = paste0(
    "Combined NSC heatmap across particle sizes",

    ifelse(
      mask_weak_cells_combined,

      paste0(
        " (|NSC| > ",
        nsc_cutoff_present_combined,
        " shown)"
      ),

      paste0(
        " (parameters with any |NSC| > ",
        nsc_cutoff_present_combined,
        ")"
      )
    )
  ),

  silent = TRUE,

  na_col = "white"
)


save_gtable_tiff(
  p3$gtable,
  file.path(
    FIG_MAIN,
    "Figure_3_NSC_combined_heatmap_by_size.tiff"
  ),
  width = 14,
  height = 8
)


## ============================================================
## 4) FIGURE 4
## ORIGINAL BAYESIAN POSTERIOR PREDICTIVE STYLE
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Figure 4: Bayesian posterior time-course profiles\n"
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
    "Could not find data/data.csv or data.csv"
  )
}


obs_df <- read.csv(
  DATA_FILE,
  encoding = "latin1",
  check.names = FALSE
) %>%

  mutate(
    ID = as.numeric(ID),

    Tissue = harmonize_tissue(
      Tissue
    ),

    SD = ifelse(
      is.na(SD),
      pmax(
        as.numeric(Conc) * 0.1,
        1e-6
      ),
      as.numeric(SD)
    ),

    Time_h = as.numeric(
      Time_h
    ),

    Conc = as.numeric(
      Conc
    )
  ) %>%

  filter(
    ID %in% 1:5
  )


load_ppc_summary <- function(id) {

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
      "Missing posterior prediction summary: ",
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


  med_col <- find_col(
    x,
    c(
      "pred_med",
      "pred_median",
      "median",
      "Median",
      "q50",
      "mc_median"
    )
  )


  lo_col <- find_col(
    x,
    c(
      "pred_lo",
      "q025",
      "lower",
      "Lower",
      "lower95"
    )
  )


  hi_col <- find_col(
    x,
    c(
      "pred_hi",
      "q975",
      "upper",
      "Upper",
      "upper95"
    )
  )


  tibble(
    ID = id,

    time = as.numeric(
      x[[time_col]]
    ),

    Tissue = harmonize_tissue(
      x[[tissue_col]]
    ),

    pred_med = as.numeric(
      x[[med_col]]
    ),

    pred_lo = as.numeric(
      x[[lo_col]]
    ),

    pred_hi = as.numeric(
      x[[hi_col]]
    )
  )
}


PPC_ALL <- bind_rows(
  lapply(
    1:5,
    load_ppc_summary
  )
)


PPC_TISSUES <- c(
  "Blood",
  "Liver",
  "Kidney",
  "Spleen",
  "Lungs",
  "GI",
  "Urine"
)


Y_FLOOR <- 1e-6


make_ppc_row <- function(
    id,
    show_title = TRUE
) {

  post <- PPC_ALL %>%

    filter(
      ID == id,
      Tissue %in%
        PPC_TISSUES
    ) %>%

    mutate(
      Tissue = factor(
        Tissue,
        levels = PPC_TISSUES
      ),

      pred_lo = pmax(
        pred_lo,
        Y_FLOOR
      ),

      pred_med = pmax(
        pred_med,
        Y_FLOOR
      ),

      pred_hi = pmax(
        pred_hi,
        Y_FLOOR
      )
    )


  obs_i <- obs_df %>%

    filter(
      ID == id,
      Tissue %in%
        PPC_TISSUES
    ) %>%

    transmute(
      time = Time_h,

      Tissue = factor(
        Tissue,
        levels = PPC_TISSUES
      ),

      obs = Conc,

      SD = SD
    )


  size_nm <- as.numeric(
    ID_SIZE_MAP[
      as.character(id)
    ]
  )


  row_title <- paste0(
    size_label_fun(
      size_nm
    ),
    " (ID ",
    id,
    ")"
  )


  ggplot(
    post,
    aes(
      x = time
    )
  ) +

    geom_ribbon(
      aes(
        ymin = pred_lo,
        ymax = pred_hi
      ),
      alpha = 0.20
    ) +

    geom_line(
      aes(
        y = pred_med
      ),
      linewidth = 1
    ) +

    geom_point(
      data = obs_i,
      aes(
        x = time,
        y = obs
      ),
      inherit.aes = FALSE,
      size = 2
    ) +

    geom_errorbar(
      data = obs_i,

      aes(
        x = time,

        ymin = pmax(
          obs - SD,
          Y_FLOOR
        ),

        ymax = obs + SD
      ),

      inherit.aes = FALSE,

      width = 0.4,

      alpha = 0.6
    ) +

    facet_wrap(
      ~ Tissue,
      scales = "free_y",
      nrow = 1,
      drop = FALSE
    ) +

    scale_y_log10() +

    labs(
      title = if (
        show_title
      ) {
        row_title
      } else {
        NULL
      },

      x = "Time (h)",

      y = "%ID"
    ) +

    theme_bw(
      base_size = 13
    )
}


PPC_ROWS <- lapply(
  1:4,
  make_ppc_row
)


FIG4 <- gridExtra::arrangeGrob(
  grobs = PPC_ROWS,
  ncol = 1
)


save_gtable_tiff(
  FIG4,

  file.path(
    FIG_MAIN,
    "Figure_4_Bayesian_timecourse_ID1_ID4.tiff"
  ),

  width = 16,

  height = 14
)


## ------------------------------------------------------------
## Figure S6 — ID 5
## ------------------------------------------------------------

post5 <- PPC_ALL %>%

  filter(
    ID == 5,
    Tissue %in%
      PPC_TISSUES
  ) %>%

  mutate(
    Tissue = factor(
      Tissue,
      levels = PPC_TISSUES
    ),

    pred_lo = pmax(
      pred_lo,
      Y_FLOOR
    ),

    pred_med = pmax(
      pred_med,
      Y_FLOOR
    ),

    pred_hi = pmax(
      pred_hi,
      Y_FLOOR
    )
  )


obs5 <- obs_df %>%

  filter(
    ID == 5,
    Tissue %in%
      PPC_TISSUES
  ) %>%

  transmute(
    time = Time_h,

    Tissue = factor(
      Tissue,
      levels = PPC_TISSUES
    ),

    obs = Conc,

    SD = SD
  )


pS6 <- ggplot(
  post5,
  aes(
    x = time
  )
) +

  geom_ribbon(
    aes(
      ymin = pred_lo,
      ymax = pred_hi
    ),
    alpha = 0.20
  ) +

  geom_line(
    aes(
      y = pred_med
    ),
    linewidth = 1
  ) +

  geom_point(
    data = obs5,
    aes(
      x = time,
      y = obs
    ),
    inherit.aes = FALSE,
    size = 2
  ) +

  geom_errorbar(
    data = obs5,

    aes(
      x = time,

      ymin = pmax(
        obs - SD,
        Y_FLOOR
      ),

      ymax = obs + SD
    ),

    inherit.aes = FALSE,

    width = 0.4,

    alpha = 0.6
  ) +

  facet_wrap(
    ~ Tissue,
    scales = "free_y"
  ) +

  labs(
    title =
      "Posterior predictive (95% CI) vs Observations — ID 5",

    x = "Time (h)",

    y = "%ID"
  ) +

  theme_bw(
    base_size = 13
  ) +

  scale_y_log10()


save_tiff(
  pS6,

  file.path(
    FIG_SUPP,
    "Figure_S6_PosteriorPredictive_ID5.tiff"
  ),

  width = 10,

  height = 7
)


## ============================================================
## 5) FIGURE 5
## RANDOM-EFFECTS META-REGRESSION
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Figure 5: random-effects meta-regression\n"
)

cat(
  "============================================================\n"
)


meta_file <- file.path(
  META_DIR,
  "meta_analysis_outputs.rds"
)


if (!file.exists(meta_file)) {

  stop(
    "Missing: ",
    meta_file
  )
}


META_OBJ <- readRDS(
  meta_file
)


META_PARAMS <- c(
  "Ktrans",
  "Kgw2b",
  "Kshed",
  "KLRESmax",
  "Kfeces",
  "KurineC"
)


META_TAGS <- setNames(
  LETTERS[1:6],
  META_PARAMS
)


## ------------------------------------------------------------
## Extract meta-regression summary table
## ------------------------------------------------------------

get_meta_table <- function(obj) {

  candidates <- c(
    "META_SUMMARY",
    "meta_summary",
    "META_FITS",
    "meta_fits"
  )


  for (
    nm in candidates
  ) {

    if (
      !is.null(
        obj[[nm]]
      ) &&
        is.data.frame(
          obj[[nm]]
        )
    ) {

      return(
        obj[[nm]]
      )
    }
  }


  stop(
    "Could not find META_SUMMARY/META_FITS in ",
    "meta_analysis_outputs.rds.\n",
    "Top-level names: ",
    paste(
      names(obj),
      collapse = ", "
    )
  )
}


META_TAB <- get_meta_table(
  META_OBJ
)


meta_param_col <- find_col(
  META_TAB,
  c(
    "Parameter",
    "parameter",
    "param"
  )
)


get_scalar <- function(
    row,
    candidates,
    required = TRUE
) {

  hit <- candidates[
    candidates %in% names(row)
  ]


  if (
    length(hit) > 0
  ) {

    return(
      as.numeric(
        row[[hit[1]]][1]
      )
    )
  }


  if (required) {

    stop(
      "Missing meta column; tried: ",
      paste(
        candidates,
        collapse = ", "
      )
    )
  }


  NA_real_
}


extract_vcov <- function(
    row,
    obj = META_OBJ,
    parameter = NULL
) {

  ## ----------------------------------------------------------
  ## 1) Try list-column covariance matrix
  ## ----------------------------------------------------------

  for (
    nm in c(
      "vcov",
      "V",
      "coef_vcov",
      "VCOV"
    )
  ) {

    if (
      nm %in%
        names(row)
    ) {

      z <- row[[nm]]

      if (
        is.list(z)
      ) {

        z <- z[[1]]
      }


      if (
        is.matrix(z) &&
          all(
            dim(z) >=
              c(
                2,
                2
              )
          )
      ) {

        return(
          z[
            1:2,
            1:2,
            drop = FALSE
          ]
        )
      }
    }
  }


  ## ----------------------------------------------------------
  ## 2) Try scalar covariance columns
  ## ----------------------------------------------------------

  vc11 <- get_scalar(
    row,
    c(
      "vc11",
      "vcov11",
      "var_beta0",
      "beta0_var",
      "V11"
    ),
    FALSE
  )


  vc22 <- get_scalar(
    row,
    c(
      "vc22",
      "vcov22",
      "var_beta1",
      "beta1_var",
      "V22"
    ),
    FALSE
  )


  vc12 <- get_scalar(
    row,
    c(
      "vc12",
      "vcov12",
      "cov_beta01",
      "beta01_cov",
      "V12"
    ),
    FALSE
  )


  if (
    all(
      is.finite(
        c(
          vc11,
          vc22,
          vc12
        )
      )
    )
  ) {

    return(
      matrix(
        c(
          vc11,
          vc12,
          vc12,
          vc22
        ),
        nrow = 2,
        ncol = 2
      )
    )
  }


  ## ----------------------------------------------------------
  ## 3) Fallback: recover covariance from META_FITS object
  ## ----------------------------------------------------------

  for (
    nm in c(
      "META_FITS",
      "meta_fits"
    )
  ) {

    tab2 <- obj[[nm]]


    if (
      is.null(tab2) ||
        !is.data.frame(tab2)
    ) {

      next
    }


    pcol2 <- intersect(
      c(
        "Parameter",
        "parameter",
        "param"
      ),
      names(tab2)
    )[1]


    if (
      is.na(pcol2)
    ) {

      next
    }


    rr <- tab2[
      as.character(
        tab2[[pcol2]]
      ) == parameter,
      ,
      drop = FALSE
    ]


    if (
      nrow(rr) != 1
    ) {

      next
    }


    if (
      "obj" %in% names(rr) &&
        is.list(
          rr$obj
        )
    ) {

      fit <- rr$obj[[1]]


      V <- tryCatch(
        as.matrix(
          vcov(fit)
        ),
        error = function(e) NULL
      )


      if (
        !is.null(V) &&
          all(
            dim(V) >=
              c(
                2,
                2
              )
          )
      ) {

        return(
          V[
            1:2,
            1:2,
            drop = FALSE
          ]
        )
      }
    }


    v11 <- get_scalar(
      rr,
      c(
        "vc11",
        "vcov11",
        "var_beta0"
      ),
      FALSE
    )


    v22 <- get_scalar(
      rr,
      c(
        "vc22",
        "vcov22",
        "var_beta1"
      ),
      FALSE
    )


    v12 <- get_scalar(
      rr,
      c(
        "vc12",
        "vcov12",
        "cov_beta01"
      ),
      FALSE
    )


    if (
      all(
        is.finite(
          c(
            v11,
            v22,
            v12
          )
        )
      )
    ) {

      return(
        matrix(
          c(
            v11,
            v12,
            v12,
            v22
          ),
          nrow = 2,
          ncol = 2
        )
      )
    }
  }


  stop(
    "Could not recover coefficient covariance matrix for ",
    parameter,
    ".\nMETA table columns: ",
    paste(
      names(row),
      collapse = ", "
    )
  )
}


build_meta_curve <- function(
    parameter
) {

  row <- META_TAB[
    as.character(
      META_TAB[[meta_param_col]]
    ) == parameter,
    ,
    drop = FALSE
  ]


  if (
    nrow(row) != 1
  ) {

    stop(
      "Expected one meta row for ",
      parameter
    )
  }


  b0 <- get_scalar(
    row,
    c(
      "beta0",
      "coef0",
      "intercept"
    )
  )


  b1 <- get_scalar(
    row,
    c(
      "beta1",
      "coef1",
      "slope"
    )
  )


  tau2 <- get_scalar(
    row,
    c(
      "tau2",
      "tau_sq",
      "tau_squared"
    )
  )


  Vb <- extract_vcov(
    row,
    META_OBJ,
    parameter
  )


  sizes <- exp(
    seq(
      log(20),
      log(6000),
      length.out = 200
    )
  )


  bind_rows(
    lapply(
      sizes,
      function(S) {

        x <- c(
          1,
          log(S)
        )


        mu <- as.numeric(
          b0 +
            b1 *
            log(S)
        )


        var_mu <- as.numeric(
          t(x) %*%
            Vb %*%
            x
        )


        var_mu <- pmax(
          var_mu,
          1e-12
        )


        var_new <- var_mu +
          pmax(
            tau2,
            0
          )


        tibble(
          Parameter = parameter,

          Size_nm = S,

          FitMean = exp(
            mu +
              0.5 *
              var_mu
          ),

          CILo = exp(
            mu -
              1.96 *
              sqrt(var_mu)
          ),

          CIHi = exp(
            mu +
              1.96 *
              sqrt(var_mu)
          ),

          PILo = exp(
            mu -
              1.96 *
              sqrt(var_new)
          ),

          PIHi = exp(
            mu +
              1.96 *
              sqrt(var_new)
          )
        )
      }
    )
  )
}


META_CURVES <- bind_rows(
  lapply(
    META_PARAMS,
    build_meta_curve
  )
)


## ------------------------------------------------------------
## Per-ID posterior medians and 95% credible intervals
## ------------------------------------------------------------

load_param_summary <- function(
    id
) {

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
      "q050",
      "q50",
      "Median"
    )
  )


  locol <- find_col(
    x,
    c(
      "q025",
      "lower",
      "Lower"
    )
  )


  hicol <- find_col(
    x,
    c(
      "q975",
      "upper",
      "Upper"
    )
  )


  tibble(
    ID = id,

    Parameter = as.character(
      x[[pcol]]
    ),

    Median = as.numeric(
      x[[medcol]]
    ),

    Q025 = as.numeric(
      x[[locol]]
    ),

    Q975 = as.numeric(
      x[[hicol]]
    )
  )
}


META_POINTS <- bind_rows(
  lapply(
    1:5,
    load_param_summary
  )
) %>%

  filter(
    Parameter %in%
      META_PARAMS
  ) %>%

  mutate(
    Size_nm = as.numeric(
      ID_SIZE_MAP[
        as.character(ID)
      ]
    ),

    ID = factor(ID)
  )


PARAM_YLAB <- list(

  Ktrans =
    expression(
      K[trans] ~
        "(h"^{-1} * ")"
    ),

  Kgw2b =
    expression(
      K[gw2b] ~
        "(h"^{-1} * ")"
    ),

  Kshed =
    expression(
      K[shed] ~
        "(h"^{-1} * ")"
    ),

  KLRESmax =
    expression(
      K[LRESmax] ~
        "(h"^{-1} * ")"
    ),

  Kfeces =
    expression(
      K[feces] ~
        "(h"^{-1} * ")"
    ),

  KurineC =
    expression(
      K[urineC] ~
        "(L h"^{-1} ~
        "kg"^{-0.75} * ")"
    )
)


make_meta_panel <- function(
    parameter
) {

  cc <- META_CURVES %>%

    filter(
      Parameter ==
        parameter
    )


  pp <- META_POINTS %>%

    filter(
      Parameter ==
        parameter
    )


  ggplot() +

    ## 95% prediction interval
    geom_ribbon(
      data = cc,

      aes(
        x = Size_nm,
        ymin = PILo,
        ymax = PIHi
      ),

      fill = "grey85",

      alpha = 1
    ) +

    ## 95% confidence interval
    geom_ribbon(
      data = cc,

      aes(
        x = Size_nm,
        ymin = CILo,
        ymax = CIHi
      ),

      fill = "grey60",

      alpha = 0.8
    ) +

    ## Meta-regression fitted curve
    geom_line(
      data = cc,

      aes(
        x = Size_nm,
        y = FitMean
      ),

      linewidth = 1,

      color = "black"
    ) +

    ## Posterior 95% CI for each development dataset
    geom_errorbar(
      data = pp,

      aes(
        x = Size_nm,
        ymin = Q025,
        ymax = Q975,
        group = ID
      ),

      width = 0.08,

      linewidth = 0.6,

      position =
        position_dodge(
          width = 0.04
        )
    ) +

    ## Posterior median
    geom_point(
      data = pp,

      aes(
        x = Size_nm,
        y = Median,
        shape = ID
      ),

      size = 2.5
    ) +

    scale_x_log10(
      breaks = c(
        20,
        220,
        1000,
        6000
      ),

      labels = c(
        "20",
        "220",
        "1000",
        "6000"
      )
    ) +

    scale_y_log10() +

    labs(
      tag = META_TAGS[[parameter]],

      x = "Particle size (nm)",

      y = PARAM_YLAB[[parameter]],

      shape = "Dataset ID"
    ) +

    theme_bw(
      base_size = 13
    ) +

    theme(
      plot.tag =
        element_text(
          face = "bold",
          size = 14
        ),

      legend.position =
        "none"
    )
}


META_PLOTS <- lapply(
  META_PARAMS,
  make_meta_panel
)


FIG5 <- gridExtra::arrangeGrob(
  grobs = META_PLOTS,
  ncol = 3
)


save_gtable_tiff(
  FIG5,

  file.path(
    FIG_MAIN,
    "Figure_5_random_effects_meta_regression.tiff"
  ),

  width = 12,

  height = 8
)


## ============================================================
## 6) FIGURE 6
## ORIGINAL EXTERNAL-EVALUATION BAR CHARTS
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Figure 6: original external evaluation bar charts\n"
)

cat(
  "============================================================\n"
)


ext_file <- file.path(
  EXT_DIR,
  "external_observed_vs_predicted.csv"
)


if (!file.exists(ext_file)) {

  stop(
    "Missing: ",
    ext_file
  )
}


EXT <- read_csv(
  ext_file,
  show_col_types = FALSE
)


required_ext <- c(
  "ID",
  "Size_nm",
  "time",
  "Tissue",
  "obs",
  "SD",
  "mc_mean",
  "pred_hi"
)


missing_ext <- setdiff(
  required_ext,
  names(EXT)
)


if (
  length(
    missing_ext
  ) > 0
) {

  stop(
    "External output is missing required columns: ",
    paste(
      missing_ext,
      collapse = ", "
    ),
    "\nAvailable: ",
    paste(
      names(EXT),
      collapse = ", "
    )
  )
}


TEST_IDS <- c(
  6,
  7,
  8
)


test_compare <- EXT %>%

  mutate(
    ID = as.numeric(ID),

    Tissue =
      harmonize_tissue(
        Tissue
      )
  ) %>%

  filter(
    ID %in%
      TEST_IDS
  )


## ------------------------------------------------------------
## Observed:
## bar = observed mean
## upper whisker = obs + SD
##
## Predicted:
## bar = MC mean
## upper whisker = 97.5% Monte Carlo percentile
## ------------------------------------------------------------

obs_bar_df <- test_compare %>%

  transmute(
    ID,
    Size_nm,
    time,
    Tissue,

    Source =
      "Observed",

    Value =
      obs,

    ErrHi =
      obs + SD
  )


pred_bar_df <- test_compare %>%

  transmute(
    ID,
    Size_nm,
    time,
    Tissue,

    Source =
      "Predicted",

    Value =
      mc_mean,

    ErrHi =
      pred_hi
  )


compare_bar_df <- bind_rows(
  obs_bar_df,
  pred_bar_df
) %>%

  mutate(
    Tissue = factor(
      Tissue,

      levels = c(
        "Blood",
        "Liver",
        "Kidney",
        "Spleen",
        "Lungs",
        "GI",
        "Urine"
      )
    ),

    Source = factor(
      Source,

      levels = c(
        "Observed",
        "Predicted"
      )
    )
  ) %>%

  group_by(
    ID,
    Size_nm,
    time,
    Tissue,
    Source
  ) %>%

  summarise(
    Value = mean(
      Value,
      na.rm = TRUE
    ),

    ErrHi = mean(
      ErrHi,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


## ============================================================
## Figure 6 A–C: non-GI tissues
## ============================================================

compare_bar_nonGI <- compare_bar_df %>%

  filter(
    Tissue != "GI"
  ) %>%

  group_by(
    Tissue
  ) %>%

  filter(
    max(
      Value,
      na.rm = TRUE
    ) > 0
  ) %>%

  ungroup()


pd <- position_dodge(
  width = 0.8
)


p6_top <- ggplot(
  compare_bar_nonGI,

  aes(
    x = Tissue,
    y = Value,
    fill = Source
  )
) +

  geom_col(
    position = pd,
    width = 0.7
  ) +

  ## Upper error bars only
  geom_errorbar(
    aes(
      ymin = Value,
      ymax = ErrHi,
      group = Source
    ),

    position = pd,

    width = 0.25,

    linewidth = 0.7
  ) +

  facet_wrap(
    ~ ID,
    nrow = 1,
    scales = "free_y"
  ) +

  labs(
    x = "Tissue",
    y = "%ID",
    fill = NULL
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    panel.grid =
      element_blank(),

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),

    strip.background =
      element_blank(),

    legend.position =
      "top"
  )


## ============================================================
## Figure 6 D–F: GI
## ============================================================

compare_bar_GI <- compare_bar_df %>%

  filter(
    Tissue == "GI"
  )


p6_bottom <- ggplot(
  compare_bar_GI,

  aes(
    x = Source,
    y = Value,
    fill = Source
  )
) +

  geom_col(
    width = 0.65
  ) +

  ## Upper error bars only
  geom_errorbar(
    aes(
      ymin = Value,
      ymax = ErrHi
    ),

    width = 0.20,

    linewidth = 0.7
  ) +

  facet_wrap(
    ~ ID,
    nrow = 1,
    scales = "free_y"
  ) +

  labs(
    x = NULL,
    y = "%ID",
    fill = NULL
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    panel.grid =
      element_blank(),

    strip.background =
      element_blank(),

    legend.position =
      "none"
  )


## ------------------------------------------------------------
## Save Figure 6 components
## ------------------------------------------------------------

save_tiff(
  p6_top,

  file.path(
    FIG_MAIN,
    "Figure_6_ABC_external_nonGI.tiff"
  ),

  width = 12,

  height = 6
)


save_tiff(
  p6_bottom,

  file.path(
    FIG_MAIN,
    "Figure_6_DEF_external_GI.tiff"
  ),

  width = 8,

  height = 5
)


## ------------------------------------------------------------
## Complete Figure 6
## ------------------------------------------------------------

FIG6 <- gridExtra::arrangeGrob(
  p6_top,
  p6_bottom,

  ncol = 1,

  heights = c(
    1.1,
    0.9
  )
)


save_gtable_tiff(
  FIG6,

  file.path(
    FIG_MAIN,
    "Figure_6_external_evaluation.tiff"
  ),

  width = 12,

  height = 10
)


## ============================================================
## 7) COMPLETION
## ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "Original-style figure reproduction completed.\n"
)

cat(
  "No scientific analysis was rerun.\n"
)


cat(
  "Main figures: ",
  normalizePath(
    FIG_MAIN,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat(
  "Supplement: ",
  normalizePath(
    FIG_SUPP,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat(
  "============================================================\n"
)