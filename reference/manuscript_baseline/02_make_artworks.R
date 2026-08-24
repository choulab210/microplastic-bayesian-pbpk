## ============================================================
## 02_make_artworks.R
## Artwork-only script
## - reads saved outputs from 01_run_model_fit.R
## - generates GOF figures
## - generates individual NSC heatmaps by particle size (if available)
## - generates combined NSC heatmap across particle sizes (if available)
## - gracefully skips NSC artwork if NSC_amount_all.csv is missing
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
})

options(stringsAsFactors = FALSE)

## ============================================================
## 0) Directories and input files
## ============================================================
results_dir <- "results"
figures_dir <- "figures"

dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

fitted_rds_file <- file.path(results_dir, "fitted_objects.rds")
gof_csv_file    <- file.path(results_dir, "gof_all.csv")
nsc_csv_file    <- file.path(results_dir, "NSC_amount_all.csv")

if (!file.exists(fitted_rds_file)) {
  stop("Missing file: ", fitted_rds_file,
       "\nPlease run 01_run_model_fit.R first.")
}

res <- readRDS(fitted_rds_file)

if (!"gof_all" %in% names(res)) {
  if (file.exists(gof_csv_file)) {
    res$gof_all <- read.csv(gof_csv_file, check.names = FALSE)
  } else {
    stop("No GOF data found in fitted_objects.rds or ", gof_csv_file)
  }
}

gof_all <- as.data.frame(res$gof_all)
ids     <- if ("ids" %in% names(res)) res$ids else sort(unique(gof_all$ID))

has_nsc <- file.exists(nsc_csv_file)

if (!has_nsc) {
  warning(
    "Missing file: ", nsc_csv_file,
    "\nNSC heatmaps will be skipped. GOF artwork will still be generated."
  )
  NSC_amount_all <- NULL
} else {
  NSC_amount_all <- read.csv(nsc_csv_file, check.names = FALSE)
}

## ============================================================
## 1) Shared plotting utilities
## ============================================================
theme_gof_poster <- function(show_legend = TRUE) {
  theme_bw(base_size = 13) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      axis.ticks = element_line(linewidth = 1.5),
      axis.title = element_text(size = 18, colour = "black", face = "bold"),
      axis.text  = element_text(size = 16, colour = "black", face = "bold"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14, face = "bold"),
      legend.position = if (show_legend) "right" else "none"
    )
}

log10_formatter <- trans_format("log10", math_format(10^.x))

calc_gof_metrics <- function(df) {
  df <- df %>%
    dplyr::filter(is.finite(obs), is.finite(pred), obs > 0, pred > 0)
  
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

## ============================================================
## 2) GOF artwork
## ============================================================
plot_overall_gof <- function(gof_all, log_scale = TRUE, show_legend = TRUE) {
  gof_all <- gof_all %>%
    dplyr::filter(is.finite(obs), is.finite(pred), obs > 0, pred > 0)
  
  met <- calc_gof_metrics(gof_all)
  
  xmin <- min(c(gof_all$obs, gof_all$pred), na.rm = TRUE) / 1.5
  xmax <- max(c(gof_all$obs, gof_all$pred), na.rm = TRUE) * 1.5
  
  if (!is.finite(xmin) || xmin <= 0) xmin <- 1e-6
  if (!is.finite(xmax) || xmax <= xmin) xmax <- xmin * 10
  
  line_df <- data.frame(
    x = c(xmin, xmax),
    y1 = c(xmin, xmax),
    y_upper = 2 * c(xmin, xmax),
    y_lower = 0.5 * c(xmin, xmax)
  )
  
  p <- ggplot(gof_all, aes(x = obs, y = pred, color = Tissue, shape = factor(ID))) +
    geom_point(size = 2.8, alpha = 0.80) +
    geom_line(data = line_df, aes(x = x, y = y1),
              inherit.aes = FALSE, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_line(data = line_df, aes(x = x, y = y_upper),
              inherit.aes = FALSE, linetype = "dotted", color = "gray40", linewidth = 0.7) +
    geom_line(data = line_df, aes(x = x, y = y_lower),
              inherit.aes = FALSE, linetype = "dotted", color = "gray40", linewidth = 0.7) +
    labs(
      x = "Observed %ID",
      y = "Predicted %ID",
      color = "Tissue",
      shape = "ID"
    ) +
    theme_gof_poster(show_legend = show_legend)
  
  if (log_scale) {
    p <- p +
      scale_x_log10(limits = c(xmin, xmax), labels = log10_formatter) +
      scale_y_log10(limits = c(xmin, xmax), labels = log10_formatter)
  } else {
    p <- p + coord_cartesian(xlim = c(xmin, xmax), ylim = c(xmin, xmax))
  }
  
  invisible(list(plot = p, met = met))
}

plot_ratio_gof <- function(gof_all, log_scale = TRUE, show_legend = FALSE) {
  df <- gof_all %>%
    dplyr::filter(is.finite(obs), is.finite(pred), obs > 0, pred > 0) %>%
    dplyr::mutate(ratio = pred / obs)
  
  if (nrow(df) == 0) {
    warning("No valid data for ratio GOF plot.")
    return(invisible(NULL))
  }
  
  pct_2fold <- mean(df$ratio >= 0.5 & df$ratio <= 2, na.rm = TRUE) * 100
  pct_3fold <- mean(df$ratio >= 1/3 & df$ratio <= 3, na.rm = TRUE) * 100
  
  p <- ggplot(df, aes(x = pred, y = ratio, color = Tissue, shape = factor(ID))) +
    geom_point(size = 2.8, alpha = 0.80) +
    geom_hline(yintercept = 2,   linetype = "dotted", color = "red3", linewidth = 0.9) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "red3", linewidth = 0.9) +
    geom_hline(yintercept = 3,   linetype = "dashed", color = "black", linewidth = 0.9) +
    geom_hline(yintercept = 1/3, linetype = "dashed", color = "black", linewidth = 0.9) +
    labs(
      x = "Predicted %ID",
      y = "Predicted/Observed",
      color = "Tissue",
      shape = "ID"
    ) +
    theme_gof_poster(show_legend = show_legend)
  
  if (log_scale) {
    p <- p +
      scale_x_log10(labels = log10_formatter) +
      scale_y_log10(labels = log10_formatter)
  }
  
  invisible(list(
    plot = p,
    metrics = data.frame(
      pct_2fold = pct_2fold,
      pct_3fold = pct_3fold,
      n = nrow(df)
    )
  ))
}

cat("\n============================================================\n")
cat("GOF metrics\n")
cat("============================================================\n")

overall_metrics <- calc_gof_metrics(gof_all)
cat("\n--- Overall goodness-of-fit metrics ---\n")
print(overall_metrics)

per_id_metrics <- bind_rows(lapply(sort(unique(gof_all$ID)), function(id) {
  df_i <- gof_all %>% filter(ID == id)
  cbind(data.frame(ID = id), calc_gof_metrics(df_i))
}))
cat("\n--- Per-ID goodness-of-fit metrics ---\n")
print(per_id_metrics)

per_tissue_metrics <- bind_rows(lapply(sort(unique(gof_all$Tissue)), function(tt) {
  df_t <- gof_all %>% filter(Tissue == tt)
  cbind(data.frame(Tissue = tt), calc_gof_metrics(df_t))
}))
cat("\n--- Per-tissue goodness-of-fit metrics ---\n")
print(per_tissue_metrics)

write.csv(overall_metrics,    file.path(results_dir, "overall_gof_metrics.csv"), row.names = FALSE)
write.csv(per_id_metrics,     file.path(results_dir, "per_id_gof_metrics_from_artworks.csv"), row.names = FALSE)
write.csv(per_tissue_metrics, file.path(results_dir, "per_tissue_gof_metrics_from_artworks.csv"), row.names = FALSE)

p1 <- plot_overall_gof(gof_all, log_scale = TRUE, show_legend = FALSE)
ratio_out <- plot_ratio_gof(gof_all, log_scale = TRUE, show_legend = FALSE)

if (!is.null(ratio_out)) {
  p2 <- ratio_out$plot
  
  p3 <- ggMarginal(
    p2,
    type = "histogram",
    margins = "y",
    yparams = list(binwidth = 0.1, fill = "#FFA488", color = "black")
  )
  
  final_gof_plot <- grid.arrange(p1$plot, p3, nrow = 1)
  
  ggsave(
    filename = file.path(figures_dir, "PBPK_GOF_plot.tiff"),
    plot = final_gof_plot,
    device = "tiff",
    width = 12,
    height = 5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

## ============================================================
## 3) NSC artwork
## ============================================================
if (!is.null(NSC_amount_all)) {
  
  ## USER SETTINGS: individual heatmaps
  nsc_top_n_heatmap_individual <- 15
  nsc_n_clusters_individual    <- 4
  heatmap_levels_individual    <- 400
  
  ## USER SETTINGS: combined heatmap
  nsc_cutoff_present_combined  <- 0.1
  mask_weak_cells_combined     <- FALSE
  n_row_clusters_combined      <- 4
  heatmap_levels_combined      <- 400
  
  ## Explicit ID -> particle size mapping
  id_size_map <- c(
    "1" = 20,
    "2" = 220,
    "3" = 1000,
    "4" = 6000,
    "5" = 220
  )
  
  ## Tissue labels
  tissue_label_map <- c(
    "Blood"  = "AUCB",
    "GI"     = "AUCGI",
    "Spleen" = "AUCS",
    "Kidney" = "AUCK",
    "Liver"  = "AUCL",
    "Lungs"  = "AUCLu",
    "Urine"  = "AUCU",
    "Rest"   = "AUCRt"
  )
  
  size_label_fun <- function(size_nm) {
    dplyr::case_when(
      size_nm == 20   ~ "20 nm",
      size_nm == 220  ~ "220 nm",
      size_nm == 1000 ~ "1 µm",
      size_nm == 6000 ~ "6 µm",
      TRUE ~ paste0(size_nm, " nm")
    )
  }
  
  required_nsc_cols <- c("ID", "Parameter", "Tissue", "NSC")
  missing_nsc_cols <- setdiff(required_nsc_cols, names(NSC_amount_all))
  if (length(missing_nsc_cols) > 0) {
    stop("NSC_amount_all is missing required columns: ",
         paste(missing_nsc_cols, collapse = ", "))
  }
  
  NSC_amount_all2 <- NSC_amount_all %>%
    mutate(
      ID = as.numeric(ID),
      Size_nm = as.numeric(id_size_map[as.character(ID)])
    ) %>%
    filter(!is.na(Size_nm))
  
  if (!"absNSC" %in% names(NSC_amount_all2)) {
    NSC_amount_all2 <- NSC_amount_all2 %>%
      mutate(absNSC = abs(NSC))
  }
  
  ## Build shared plotting dataframe
  NSC_plot_df <- bind_rows(
    NSC_amount_all2 %>%
      filter(ID %in% c(1, 3, 4)) %>%
      mutate(SizeGroup = Size_nm),
    
    NSC_amount_all2 %>%
      filter(ID %in% c(2, 5)) %>%
      group_by(Size_nm, Parameter, Tissue) %>%
      summarise(
        NSC    = mean(NSC, na.rm = TRUE),
        absNSC = mean(abs(NSC), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        ID = 25,
        SizeGroup = Size_nm
      )
  ) %>%
    mutate(
      Tissue      = as.character(Tissue),
      TissueLabel = recode(Tissue, !!!tissue_label_map),
      SizeLabel   = size_label_fun(SizeGroup),
      ColLabel    = paste(SizeLabel, TissueLabel, sep = " | ")
    ) %>%
    filter(!is.na(TissueLabel))
  
  cat("\n============================================================\n")
  cat("NSC artwork\n")
  cat("============================================================\n")
  
  ## ============================================================
  ## 4) Individual NSC heatmaps by particle size
  ## ============================================================
  global_top_pars <- NSC_plot_df %>%
    filter(is.finite(NSC)) %>%
    group_by(Parameter) %>%
    summarise(
      global_mean_abs = mean(abs(NSC), na.rm = TRUE),
      global_max_abs  = max(abs(NSC), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(global_mean_abs), desc(global_max_abs)) %>%
    slice_head(n = nsc_top_n_heatmap_individual) %>%
    pull(Parameter)
  
  cat("\n--- SAME parameters used in individual heatmaps ---\n")
  print(global_top_pars)
  
  global_color_lim_individual <- max(abs(NSC_plot_df$NSC), na.rm = TRUE)
  if (!is.finite(global_color_lim_individual) || global_color_lim_individual == 0) {
    global_color_lim_individual <- 1
  }
  
  cat("\n--- SAME color limit used in individual heatmaps ---\n")
  print(global_color_lim_individual)
  
  make_nsc_pheatmap_size <- function(nsc_df,
                                     size_nm,
                                     selected_pars,
                                     n_clusters,
                                     fixed_lim,
                                     heatmap_levels) {
    
    df_size <- nsc_df %>%
      filter(SizeGroup == size_nm, is.finite(NSC), Parameter %in% selected_pars)
    
    if (nrow(df_size) == 0) return(NULL)
    
    plot_df <- df_size %>%
      group_by(Parameter, TissueLabel) %>%
      summarise(NSC = mean(NSC, na.rm = TRUE), .groups = "drop")
    
    mat <- plot_df %>%
      pivot_wider(names_from = TissueLabel, values_from = NSC) %>%
      as.data.frame()
    
    rownames(mat) <- mat$Parameter
    mat$Parameter <- NULL
    mat <- as.matrix(mat)
    
    desired_cols <- c("AUCB","AUCGI","AUCS","AUCK","AUCL","AUCLu","AUCU","AUCRt")
    desired_cols <- desired_cols[desired_cols %in% colnames(mat)]
    mat <- mat[, desired_cols, drop = FALSE]
    
    final_row_order <- selected_pars[selected_pars %in% rownames(mat)]
    mat <- mat[final_row_order, , drop = FALSE]
    
    if (nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
    
    mat_for_clust <- mat
    mat_for_clust[!is.finite(mat_for_clust)] <- 0
    
    hc_rows <- hclust(dist(mat_for_clust))
    row_clusters <- cutree(hc_rows, k = min(n_clusters, nrow(mat_for_clust)))
    
    ann_row <- data.frame(Cluster = factor(row_clusters))
    rownames(ann_row) <- rownames(mat_for_clust)
    
    ann_colors_local <- list(
      Cluster = c(
        "1" = "#9C6ADE",
        "2" = "#D9D9E8",
        "3" = "#F0E442",
        "4" = "#E6C300"
      )[levels(ann_row$Cluster)]
    )
    
    title_lab <- paste0("NSC-", size_label_fun(size_nm))
    
    pheatmap(
      mat,
      cluster_rows = hc_rows,
      cluster_cols = FALSE,
      annotation_row = ann_row,
      annotation_colors = ann_colors_local,
      color = colorRampPalette(c("#D7191C", "#F7F7F7", "#2B83BA"))(heatmap_levels),
      breaks = seq(-fixed_lim, fixed_lim, length.out = heatmap_levels + 1),
      border_color = "grey40",
      cellwidth = 28,
      cellheight = 16,
      fontsize = 10,
      fontsize_row = 10,
      fontsize_col = 10,
      angle_col = 270,
      main = title_lab,
      silent = TRUE,
      na_col = "white"
    )
  }
  
  final_sizes <- c(20, 220, 1000, 6000)
  
  for (sz in final_sizes) {
    p_sz <- make_nsc_pheatmap_size(
      nsc_df         = NSC_plot_df,
      size_nm        = sz,
      selected_pars  = global_top_pars,
      n_clusters     = nsc_n_clusters_individual,
      fixed_lim      = global_color_lim_individual,
      heatmap_levels = heatmap_levels_individual
    )
    
    if (is.null(p_sz)) next
    
    size_file <- ifelse(sz >= 1000, paste0(sz/1000, "um"), paste0(sz, "nm"))
    
    ggsave(
      filename = file.path(figures_dir, paste0("NSC_heatmap_", size_file, ".tiff")),
      plot = p_sz$gtable,
      device = "tiff",
      width = 6,
      height = 5,
      units = "in",
      dpi = 600,
      compression = "lzw"
    )
  }
  
  ## ============================================================
  ## 5) Combined NSC heatmap across particle sizes
  ## ============================================================
  selected_pars_combined <- NSC_plot_df %>%
    group_by(Parameter) %>%
    summarise(
      max_abs = max(abs(NSC), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(max_abs), max_abs > nsc_cutoff_present_combined) %>%
    arrange(desc(max_abs)) %>%
    pull(Parameter)
  
  cat("\n--- Parameters retained in combined heatmap (any |NSC| >",
      nsc_cutoff_present_combined, ") ---\n")
  print(selected_pars_combined)
  
  combined_df2 <- NSC_plot_df %>%
    filter(Parameter %in% selected_pars_combined)
  
  mat_combined <- combined_df2 %>%
    select(Parameter, ColLabel, NSC) %>%
    group_by(Parameter, ColLabel) %>%
    summarise(NSC = mean(NSC, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = ColLabel, values_from = NSC) %>%
    as.data.frame()
  
  if (nrow(mat_combined) > 0) {
    
    rownames(mat_combined) <- mat_combined$Parameter
    mat_combined$Parameter <- NULL
    mat_combined <- as.matrix(mat_combined)
    
    desired_cols_combined <- c(
      "20 nm | AUCB",  "20 nm | AUCGI",  "20 nm | AUCS",  "20 nm | AUCK",
      "20 nm | AUCL",  "20 nm | AUCLu",  "20 nm | AUCU",  "20 nm | AUCRt",
      "220 nm | AUCB", "220 nm | AUCGI", "220 nm | AUCS", "220 nm | AUCK",
      "220 nm | AUCL", "220 nm | AUCLu", "220 nm | AUCU", "220 nm | AUCRt",
      "1 µm | AUCB",   "1 µm | AUCGI",   "1 µm | AUCS",   "1 µm | AUCK",
      "1 µm | AUCL",   "1 µm | AUCLu",   "1 µm | AUCU",   "1 µm | AUCRt",
      "6 µm | AUCB",   "6 µm | AUCGI",   "6 µm | AUCS",   "6 µm | AUCK",
      "6 µm | AUCL",   "6 µm | AUCLu",   "6 µm | AUCU",   "6 µm | AUCRt"
    )
    
    desired_cols_combined <- desired_cols_combined[desired_cols_combined %in% colnames(mat_combined)]
    mat_combined <- mat_combined[, desired_cols_combined, drop = FALSE]
    
    mat_plot_combined <- mat_combined
    if (mask_weak_cells_combined) {
      mat_plot_combined[abs(mat_plot_combined) <= nsc_cutoff_present_combined] <- NA_real_
    }
    
    col_ann <- data.frame(
      SizeGroup = sub(" \\|.*$", "", colnames(mat_combined))
    )
    rownames(col_ann) <- colnames(mat_combined)
    col_ann$SizeGroup <- factor(
      col_ann$SizeGroup,
      levels = c("20 nm", "220 nm", "1 µm", "6 µm")
    )
    
    col_ann_colors <- list(
      SizeGroup = c(
        "20 nm"  = "#66C2A5",
        "220 nm" = "#FC8D62",
        "1 µm"   = "#8DA0CB",
        "6 µm"   = "#E78AC3"
      )
    )
    
    mat_for_clust_combined <- mat_combined
    mat_for_clust_combined[!is.finite(mat_for_clust_combined)] <- 0
    
    hc_rows_combined <- hclust(dist(mat_for_clust_combined))
    row_clusters_combined <- cutree(
      hc_rows_combined,
      k = min(n_row_clusters_combined, nrow(mat_for_clust_combined))
    )
    
    row_ann <- data.frame(
      Cluster = factor(row_clusters_combined)
    )
    rownames(row_ann) <- rownames(mat_for_clust_combined)
    
    row_ann_colors <- list(
      Cluster = c(
        "1" = "#9C6ADE",
        "2" = "#D9D9E8",
        "3" = "#F0E442",
        "4" = "#E6C300"
      )[levels(row_ann$Cluster)]
    )
    
    global_lim_combined <- max(abs(mat_combined), na.rm = TRUE)
    if (!is.finite(global_lim_combined) || global_lim_combined == 0) global_lim_combined <- 1
    
    cat("\n--- Combined heatmap color limit ---\n")
    print(global_lim_combined)
    
    p_combined <- pheatmap(
      mat_plot_combined,
      cluster_rows = hc_rows_combined,
      cluster_cols = FALSE,
      annotation_row = row_ann,
      annotation_col = col_ann,
      annotation_colors = c(row_ann_colors, col_ann_colors),
      color = colorRampPalette(c("#D7191C", "#F7F7F7", "#2B83BA"))(heatmap_levels_combined),
      breaks = seq(-global_lim_combined, global_lim_combined, length.out = heatmap_levels_combined + 1),
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
          paste0(" (|NSC| > ", nsc_cutoff_present_combined, " shown)"),
          paste0(" (parameters with any |NSC| > ", nsc_cutoff_present_combined, ")")
        )
      ),
      silent = TRUE,
      na_col = "white"
    )
    
    ggsave(
      filename = file.path(figures_dir, "NSC_combined_heatmap_by_size.tiff"),
      plot = p_combined$gtable,
      device = "tiff",
      width = 14,
      height = 8,
      units = "in",
      dpi = 600,
      compression = "lzw"
    )
  }
  
} else {
  cat("\n============================================================\n")
  cat("NSC artwork skipped\n")
  cat("============================================================\n")
  cat("results/NSC_amount_all.csv was not found.\n")
  cat("GOF artwork has still been generated.\n")
}

cat("\nArtwork generation completed.\n")