ground_truth_cache <- new.env(parent = emptyenv())

normalize_resolution <- function(x) {
  x <- gsub("\u00b5m", "um", as.character(x), fixed = TRUE)
  factor(x, levels = resolution_levels)
}

format_p <- function(x) {
  ifelse(is.na(x), "P=NA", ifelse(x < 0.001, "P<0.001", paste0("P=", signif(x, 3))))
}

theme_figure <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.title = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.2, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size, color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, color = "grey25", hjust = 0),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(4, "pt"),
      panel.spacing = unit(1.1, "lines"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

check_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
}

read_tsv <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

relative_path <- function(path, root) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  path <- normalizePath(path, mustWork = TRUE)
  root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(path, prefix)) stop("File is outside its input root: ", path)
  substring(path, nchar(prefix) + 1L)
}

read_ground_truth <- function(root, bin_size, sample) {
  key <- paste(bin_size, sample, sep = "::")
  if (exists(key, ground_truth_cache, inherits = FALSE)) {
    return(get(key, ground_truth_cache, inherits = FALSE))
  }
  path <- file.path(root, paste0("bin", bin_size), paste0(sample, ".tsv"))
  if (!file.exists(path)) stop("Ground-truth file not found: ", path)
  value <- read_tsv(path)
  check_columns(value, c("Hclone", "pair_name"), path)
  result <- list(data = value, path = path)
  assign(key, result, ground_truth_cache)
  result
}

summarise_sipeac_file <- function(path, ground_truth_root) {
  bin_name <- basename(dirname(path))
  bin_size <- suppressWarnings(as.integer(sub("^bin", "", bin_name)))
  sample <- sub("_final_e100\\.tsv$", "", basename(path))
  if (!is.finite(bin_size) || !bin_size %in% bin_sizes || !nzchar(sample)) {
    stop("Unexpected SIPEAC path: ", path)
  }

  prediction <- read_tsv(path)
  check_columns(prediction, c("Hclone", "pair_name", "bin_num"), path)
  prediction <- prediction[!duplicated(prediction$pair_name), , drop = FALSE]
  prediction$bin_num <- as.numeric(prediction$bin_num)
  ground_truth <- read_ground_truth(ground_truth_root, bin_size, sample)

  rows <- lapply(clone_size_thresholds, function(threshold) {
    threshold_data <- prediction[
      is.finite(prediction$bin_num) & prediction$bin_num >= threshold,
      ,
      drop = FALSE
    ]
    evaluated <- threshold_data[
      threshold_data$Hclone %in% ground_truth$data$Hclone,
      ,
      drop = FALSE
    ]
    n_correct <- sum(evaluated$pair_name %in% ground_truth$data$pair_name)
    n_eval <- nrow(evaluated)
    data.frame(
      cur_sample = sample,
      bin_size = bin_size,
      threshold = threshold,
      SIPEAC_n_pred_total = nrow(prediction),
      SIPEAC_n_after_threshold = nrow(threshold_data),
      SIPEAC_n_eval = n_eval,
      SIPEAC_n_correct = n_correct,
      SIPEAC = if (n_eval) n_correct / n_eval else NA_real_,
      source_sipeac = relative_path(path, dirname(dirname(path))),
      source_ground_truth = relative_path(ground_truth$path, ground_truth_root),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

empty_repair_result <- function() {
  list(
    n_pred_total = NA_integer_,
    n_eval = NA_integer_,
    n_correct = NA_integer_,
    accuracy = NA_real_,
    source = NA_character_
  )
}

repair_fallback_root <- function(root) {
  getOption("r4q7.repair_fallbacks", list())[[normalizePath(root, mustWork = TRUE)]]
}

register_repair_recovery <- function(roots, fallback, manifest_path) {
  if (!file.exists(manifest_path)) return(invisible(FALSE))
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  check_columns(manifest, c("status", "base_wiring", "n_steps", "relative_result", "result_md5"), manifest_path)
  if (!all(manifest$status == "completed" & manifest$base_wiring == "ignored" & manifest$n_steps == 2000L)) {
    stop("Recovery provenance is not a completed legacy run")
  }
  paths <- file.path(fallback, manifest$relative_result)
  stopifnot(all(file.exists(paths)))
  stopifnot(all(startsWith(normalizePath(paths), paste0(normalizePath(fallback), .Platform$file.sep))))
  stopifnot(identical(unname(tools::md5sum(paths)), manifest$result_md5))
  fallbacks <- getOption("r4q7.repair_fallbacks", list())
  fallbacks[[normalizePath(roots$repair_b100, mustWork = TRUE)]] <- normalizePath(fallback)
  options(r4q7.repair_fallbacks = fallbacks)
  invisible(TRUE)
}

summarise_repair <- function(root, sample, bin_size, threshold, ground_truth) {
  result_dir <- file.path(root, paste0("bin", bin_size), paste0(sample, "_", threshold))
  files <- sort(list.files(result_dir, pattern = "analysis_result.*\\.tsv$", full.names = TRUE))
  source_root <- root
  source_prefix <- ""
  fallback <- repair_fallback_root(root)
  if (!length(files) && !is.null(fallback)) {
    result_dir <- file.path(fallback, paste0("bin", bin_size), paste0(sample, "_", threshold))
    files <- sort(list.files(result_dir, pattern = "analysis_result.*\\.tsv$", full.names = TRUE))
    source_root <- fallback
    source_prefix <- "recovery:"
  }
  if (!length(files)) return(empty_repair_result())
  if (length(files) > 1L) stop("Multiple repair result files found in: ", result_dir)

  prediction <- read_tsv(files[[1]])
  check_columns(prediction, c("chainA", "chainB"), files[[1]])
  prediction$pair_prediction <- paste0(prediction$chainA, "_", prediction$chainB)
  prediction <- prediction[!duplicated(prediction$pair_prediction), , drop = FALSE]
  evaluated <- prediction[prediction$chainA %in% ground_truth$Hclone, , drop = FALSE]
  n_correct <- sum(evaluated$pair_prediction %in% ground_truth$pair_name)
  n_eval <- nrow(evaluated)

  list(
    n_pred_total = nrow(prediction),
    n_eval = n_eval,
    n_correct = n_correct,
    accuracy = if (n_eval) n_correct / n_eval else NA_real_,
    source = paste0(source_prefix, relative_path(files[[1]], source_root))
  )
}

paired_wilcox <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L) {
    return(c(n = length(x), mean_A = if (length(x)) mean(x) else NA,
             mean_B = if (length(y)) mean(y) else NA,
             delta = if (length(x)) mean(x - y) else NA, statistic = NA, p = NA))
  }
  delta <- x - y
  delta[abs(delta) < sqrt(.Machine$double.eps)] <- 0
  if (all(delta == 0)) {
    return(c(n = length(x), mean_A = mean(x), mean_B = mean(y), delta = 0, statistic = 0, p = 1))
  }
  test <- stats::wilcox.test(
    delta,
    mu = 0,
    exact = FALSE,
    correct = TRUE,
    digits.rank = 12
  )
  c(
    n = length(x),
    mean_A = mean(x),
    mean_B = mean(y),
    delta = mean(x - y),
    statistic = unname(test$statistic),
    p = test$p.value
  )
}

unpaired_wilcox <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  output <- c(n_A = length(x), n_B = length(y),
              mean_A = if (length(x)) mean(x) else NA_real_,
              mean_B = if (length(y)) mean(y) else NA_real_,
              delta = if (length(x) && length(y)) mean(x) - mean(y) else NA_real_,
              statistic = NA_real_, rank_biserial = NA_real_, p = NA_real_)
  if (min(length(x), length(y)) < 3L || length(unique(signif(c(x, y), 12))) < 2L) {
    return(output)
  }
  test <- stats::wilcox.test(x, y, paired = FALSE, alternative = "two.sided",
                            exact = FALSE, correct = TRUE, digits.rank = 12)
  output["statistic"] <- unname(test$statistic)
  output["rank_biserial"] <- 2 * unname(test$statistic) / (length(x) * length(y)) - 1
  output["p"] <- test$p.value
  output
}

calculate_pairwise_tests <- function(sample_data, threshold = 10L, paired = FALSE) {
  comparisons <- list(
    c("SIPEAC", "repair_b100"),
    c("SIPEAC", "repair_default"),
    c("repair_b100", "repair_default")
  )
  rows <- list()
  index <- 1L
  for (bin_size in bin_sizes) {
    subset_data <- sample_data[
      sample_data$bin_size == bin_size & sample_data$threshold == threshold,
      ,
      drop = FALSE
    ]
    for (comparison in comparisons) {
      if (paired) {
        fit <- paired_wilcox(subset_data[[comparison[[1]]]], subset_data[[comparison[[2]]]])
        result <- c(n_A = fit[["n"]], n_B = fit[["n"]],
                    fit[c("mean_A", "mean_B", "delta", "statistic", "p")],
                    rank_biserial = NA_real_)
      } else {
        result <- unpaired_wilcox(subset_data[[comparison[[1]]]], subset_data[[comparison[[2]]]])
      }
      rows[[index]] <- data.frame(
        bin_size = bin_size,
        resolution = paste0(bin_size / 2, "um"),
        comparison = paste(comparison, collapse = " vs "),
        n_A = as.integer(result[["n_A"]]),
        n_B = as.integer(result[["n_B"]]),
        mean_A = result[["mean_A"]],
        mean_B = result[["mean_B"]],
        delta = result[["delta"]],
        wilcoxon_statistic = result[["statistic"]],
        rank_biserial = result[["rank_biserial"]],
        p_value = result[["p"]],
        test_method = if (paired) "Wilcoxon signed-rank" else "Wilcoxon rank-sum",
        paired = paired,
        alternative = "two.sided",
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  bind_rows(rows)
}

calculate_trend_tests <- function(sample_data) {
  output <- sample_data %>%
    group_by(bin_size, resolution, threshold) %>%
    group_modify(~ as.data.frame(as.list(paired_wilcox(.x$SIPEAC, .x$repair_default)))) %>%
    ungroup() %>%
    rename(n_paired = n, wilcoxon_statistic = statistic, paired_p = p) %>%
    mutate(comparison = "SIPEAC vs repair_default")
  output
}

collapse_patient_pairs <- function(sample_data, patient_level = FALSE) {
  if (!patient_level) return(sample_data)
  sample_data %>%
    group_by(bin_size, resolution, threshold, patient_group) %>%
    summarise(
      across(c(SIPEAC, repair_b100, repair_default),
             ~ if (any(is.finite(.x))) mean(.x[is.finite(.x)]) else NA_real_),
      .groups = "drop"
    )
}

cluster_mean_se <- function(x, cluster) {
  n_groups <- length(unique(cluster))
  if (n_groups < 2L) return(NA_real_)
  residual_sums <- tapply(x - mean(x), cluster, sum)
  sqrt(n_groups / (n_groups - 1) * sum(residual_sums^2)) / length(x)
}

build_input_manifest <- function(sample_data, roots) {
  sources <- list(
    source_sipeac = c("SIPEAC prediction", roots$sipeac),
    source_ground_truth = c("Ground truth", roots$ground_truth),
    source_repair_b100 = c("repair CLI -b 100", roots$repair_b100),
    source_repair_default = c("repair default -b e", roots$repair_default)
  )
  rows <- lapply(names(sources), function(column) {
    relative <- sort(unique(stats::na.omit(sample_data[[column]])))
    full_path <- file.path(sources[[column]][[2]], relative)
    recovered <- startsWith(relative, "recovery:")
    if (any(recovered)) {
      fallback <- repair_fallback_root(sources[[column]][[2]])
      if (is.null(fallback)) stop("Recovery input root is not registered")
      full_path[recovered] <- file.path(fallback, substring(relative[recovered], nchar("recovery:") + 1L))
    }
    info <- file.info(full_path)
    data.frame(
      role = sources[[column]][[1]],
      relative_path = relative,
      source_batch = ifelse(recovered, "legacy_completion", "original"),
      bytes = info$size,
      md5 = unname(tools::md5sum(full_path)),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

save_r4q7_ab <- function(result, figure_dir, table_dir = figure_dir, prefix = "R4_Q7_AB") {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list(
    figure = list(plot = result$figure, width = 8.2, height = 7.0),
    panel_A = list(plot = result$panel_a, width = 8.2, height = 3.2),
    panel_B = list(plot = result$panel_b, width = 8.2, height = 3.4)
  )

  for (name in names(plots)) {
    item <- plots[[name]]
    stem <- file.path(figure_dir, paste0(prefix, "_", name))
    ggsave(
      paste0(stem, ".pdf"),
      item$plot,
      width = item$width,
      height = item$height,
      device = cairo_pdf,
      bg = "white"
    )
    ggsave(
      paste0(stem, ".svg"),
      item$plot,
      width = item$width,
      height = item$height,
      device = svglite::svglite,
      bg = "white"
    )
    ggsave(
      paste0(stem, ".png"),
      item$plot,
      width = item$width,
      height = item$height,
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
  }

  write.csv(
    result$sample_data,
    file.path(table_dir, paste0(prefix, "_three_way_from_raw.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    result$trend_summary,
    file.path(table_dir, paste0(prefix, "_clone_size_summary.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    result$matched_summary,
    file.path(table_dir, paste0(prefix, "_mxo10_matched_summary.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    result$pairwise_tests,
    file.path(table_dir, paste0(prefix, "_mxo10_unpaired_tests.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    result$pairwise_annotations,
    file.path(table_dir, paste0(prefix, "_mxo10_plot_annotations.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    result$paired_sensitivity_tests,
    file.path(table_dir, paste0(prefix, "_mxo10_paired_sensitivity_tests.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    result$trend_tests,
    file.path(table_dir, paste0(prefix, "_threshold_paired_tests.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    result$specimen_trend_tests,
    file.path(table_dir, paste0(prefix, "_sample_level_sensitivity_tests.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    result$matched_data,
    file.path(table_dir, paste0(prefix, "_mxo10_matched_samples.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  if (!is.null(result$input_manifest)) {
    write.csv(
      result$input_manifest,
      file.path(table_dir, paste0(prefix, "_input_manifest.csv")),
      row.names = FALSE,
      quote = TRUE
    )
  }
  writeLines(capture.output(sessionInfo()), file.path(table_dir, paste0(prefix, "_sessionInfo.txt")))
}
