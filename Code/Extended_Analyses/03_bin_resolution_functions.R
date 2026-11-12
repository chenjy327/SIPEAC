normalize_molecule_id_r1q3 <- function(x, sample_id) {
  x <- as.character(x)
  prefix <- paste0(sample_id, "_")
  has_prefix <- startsWith(x, prefix)
  x[has_prefix] <- substring(x[has_prefix], nchar(prefix) + 1L)
  x
}

assert_unique_map_r1q3 <- function(x, key, values, label) {
  stopifnot(is.data.frame(x), key %in% names(x), all(values %in% names(x)))
  check <- as.data.table(x)[, unique(.SD), .SDcols = c(key, values)]
  ambiguous <- check[, .N, by = key][N > 1L]
  if (nrow(ambiguous) > 0L) {
    stop(
      label, " has ambiguous ", key, " mappings (n = ",
      nrow(ambiguous), ").", call. = FALSE
    )
  }
  invisible(TRUE)
}

required_input_paths_r1q3 <- function(sample_id) {
  c(
    dedup = file.path(
      dirs$bin1, paste0(sample_id, "_dedup.info.BinID.1.csv")
    ),
    bin50_airr = file.path(
      BIN50_INPUT_DIR, paste0(sample_id, "_info.txt")
    ),
    score = file.path(SCORE_DIR, paste0(sample_id, "_cdr3aa.csv")),
    ground_truth = file.path(GT_DIR, paste0(sample_id, ".tsv"))
  )
}

prepare_canonical_molecules_r1q3 <- function(sample_id) {
  paths <- required_input_paths_r1q3(sample_id)
  missing <- paths[!file.exists(paths)]
  empty <- paths[file.exists(paths) & file.info(paths)$size <= 0L]
  if (length(missing) > 0L || length(empty) > 0L) {
    stop(
      sample_id, ": invalid input(s): ",
      paste(c(names(missing), names(empty)), collapse = ", "),
      call. = FALSE
    )
  }

  airr <- fread(
    paths[["bin50_airr"]],
    select = c(
      "sequence_id", "sample", "type", "locus", "clone_id",
      "clone_name", "binid"
    ),
    colClasses = list(character = c(
      "sequence_id", "sample", "type", "locus", "clone_id",
      "clone_name", "binid"
    )),
    showProgress = FALSE
  )
  airr <- airr[
    sample == sample_id & type == "Spatial" &
      locus %chin% c("IGH", "IGK", "IGL") &
      !is.na(binid) & nzchar(binid)
  ]
  if (nrow(airr) == 0L) {
    stop(sample_id, ": no spatial Ig records in Bin50 AIRR input.",
         call. = FALSE)
  }
  airr[, expected_clone_name := paste(locus, clone_id, sep = "_")]
  if (any(is.na(airr$clone_name)) ||
      any(airr$clone_name != airr$expected_clone_name)) {
    stop(sample_id, ": clone_name is inconsistent with locus + clone_id.",
         call. = FALSE)
  }
  airr[, molecule_id := normalize_molecule_id_r1q3(sequence_id, sample_id)]
  airr <- unique(airr[, .(
    molecule_id, sequence_id, locus, clone_id, clone_name,
    bin50_id = binid
  )])
  assert_unique_map_r1q3(
    airr, "molecule_id",
    c("sequence_id", "locus", "clone_id", "clone_name", "bin50_id"),
    paste0(sample_id, " Bin50 AIRR")
  )
  if (anyDuplicated(airr$molecule_id)) {
    stop(sample_id, ": duplicated spatial AIRR molecule_id.", call. = FALSE)
  }

  dedup <- fread(
    paths[["dedup"]],
    select = c("id", "binid"),
    colClasses = list(character = c("id", "binid")),
    showProgress = FALSE
  )
  dedup[, molecule_id := normalize_molecule_id_r1q3(id, sample_id)]
  dedup <- unique(dedup[
    !is.na(molecule_id) & nzchar(molecule_id) &
      !is.na(binid) & nzchar(binid),
    .(molecule_id, bin1_id = binid)
  ])
  assert_unique_map_r1q3(
    dedup, "molecule_id", "bin1_id", paste0(sample_id, " bin1 dedup")
  )
  if (anyDuplicated(dedup$molecule_id)) {
    stop(sample_id, ": duplicated dedup molecule_id.", call. = FALSE)
  }

  canonical <- merge(airr, dedup, by = "molecule_id", all = FALSE, sort = FALSE)
  if (nrow(canonical) != nrow(airr)) {
    missing_ids <- setdiff(airr$molecule_id, canonical$molecule_id)
    stop(
      sample_id, ": not all Bin50 spatial Ig molecules map to bin1 (missing n = ",
      length(missing_ids), ").", call. = FALSE
    )
  }
  if (anyDuplicated(canonical$molecule_id)) {
    stop(sample_id, ": canonical molecule table is not one-to-one.",
         call. = FALSE)
  }
  setcolorder(canonical, c(
    "molecule_id", "sequence_id", "locus", "clone_id", "clone_name",
    "bin1_id", "bin50_id"
  ))
  canonical[, sample := sample_id]
  setcolorder(canonical, c("sample", setdiff(names(canonical), "sample")))

  audit <- data.frame(
    sample = sample_id,
    bin50_spatial_ig_molecules = nrow(airr),
    matched_molecules = nrow(canonical),
    match_rate = nrow(canonical) / nrow(airr),
    n_H_clones = uniqueN(canonical[locus == "IGH", clone_name]),
    n_L_clones = uniqueN(canonical[locus %chin% c("IGK", "IGL"), clone_name]),
    stringsAsFactors = FALSE
  )
  list(canonical = canonical, audit = audit, paths = paths)
}

build_resolution_input_r1q3 <- function(canonical, resolution) {
  stopifnot(resolution %in% RESOLUTION_LEVELS)
  bin_col <- if (resolution == "Bin1") "bin1_id" else "bin50_id"
  spatial <- canonical[, .(
    binid = get(bin_col), clone_name, locus
  )]
  if (anyNA(spatial$binid) || any(!nzchar(spatial$binid))) {
    stop("Missing bin IDs in ", resolution, " canonical input.",
         call. = FALSE)
  }

  clone_num <- spatial[, .(clone_num = .N), by = clone_name]
  bin_num <- unique(spatial[, .(binid, clone_name)])[
    , .(bin_num = .N), by = clone_name
  ]
  df_clone_num <- merge(clone_num, bin_num, by = "clone_name", all = TRUE)
  if (anyNA(df_clone_num)) {
    stop("Incomplete clone counts for ", resolution, ".", call. = FALSE)
  }

  unit_flags <- unique(spatial[, .(binid, locus)])[
    , .(
      has_H = any(locus == "IGH"),
      has_L = any(locus %chin% c("IGK", "IGL"))
    ),
    by = binid
  ]
  hl_bins <- unit_flags[has_H & has_L, binid]

  active <- spatial[binid %chin% hl_bins]
  input_mat <- if (nrow(active) > 0L) {
    table(active$binid, active$clone_name)
  } else {
    matrix(integer(), nrow = 0L, ncol = 0L)
  }

  list(
    input_mat = input_mat,
    df_clone_num = as.data.frame(df_clone_num),
    n_units = uniqueN(spatial$binid),
    n_HL_units = length(hl_bins),
    n_H_clones = uniqueN(spatial[locus == "IGH", clone_name]),
    n_L_clones = uniqueN(spatial[locus %chin% c("IGK", "IGL"), clone_name]),
    n_active_clones = uniqueN(active$clone_name)
  )
}

load_fixed_ground_truth_r1q3 <- function(gt_path, canonical) {
  gt_raw <- load_ground_truth(gt_path)
  if (nrow(gt_raw) == 0L) {
    stop("Ground truth is empty: ", gt_path, call. = FALSE)
  }
  common_clones <- unique(canonical$clone_name)
  gt <- gt_raw[
    gt_raw$Hclone %in% common_clones & gt_raw$Lclone %in% common_clones,
    , drop = FALSE
  ]
  gt <- distinct(gt, pair_name, .keep_all = TRUE)
  if (nrow(gt) == 0L) {
    stop("No ground-truth pair remains in the canonical clone universe: ",
         gt_path, call. = FALSE)
  }
  list(
    gt = gt,
    n_raw_pairs = nrow(distinct(gt_raw, pair_name, .keep_all = TRUE)),
    n_fixed_pairs = nrow(gt),
    n_fixed_H = n_distinct(gt$Hclone)
  )
}

validate_score_coverage_r1q3 <- function(score_mat, canonical, sample_id) {
  h_ids <- paste0(
    "clone_", sub("^IGH_", "", unique(canonical[locus == "IGH", clone_name]))
  )
  l_ids <- paste0(
    "clone_",
    sub("^(IGK|IGL)_", "", unique(canonical[locus %chin% c("IGK", "IGL"), clone_name]))
  )
  missing_h <- setdiff(h_ids, rownames(score_mat))
  missing_l <- setdiff(l_ids, colnames(score_mat))
  if (length(missing_h) > 0L || length(missing_l) > 0L) {
    stop(
      sample_id, ": CDR3 score coverage is incomplete (missing H = ",
      length(missing_h), ", missing L = ", length(missing_l), ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

empty_prediction_table_r1q3 <- function() {
  data.frame(
    sample = character(), resolution = character(),
    Hclone = character(), Lclone = character(), pair_name = character(),
    evaluable_prediction = logical(), correct_prediction = logical(),
    stringsAsFactors = FALSE
  )
}

run_resolution_r1q3 <- function(sample_id, resolution, resolution_input,
                                 score_mat, fixed_gt) {
  gt <- fixed_gt$gt
  if (resolution_input$n_HL_units == 0L) {
    output <- data.frame()
  } else {
    output <- sipeac_pair_greedy(
      resolution_input$input_mat,
      score_mat,
      gt,
      resolution_input$df_clone_num
    )
  }

  if (nrow(output) > 0L) {
    output <- distinct(output, pair_name, .keep_all = TRUE)
    output$sample <- sample_id
    output$resolution <- resolution
    output$evaluable_prediction <- output$Hclone %in% gt$Hclone
    output$correct_prediction <- output$pair_name %in% gt$pair_name
    output <- output[, c(
      "sample", "resolution",
      setdiff(names(output), c("sample", "resolution"))
    ), drop = FALSE]
  }

  evaluable <- if (nrow(output) > 0L) {
    output[output$evaluable_prediction, , drop = FALSE]
  } else {
    data.frame()
  }
  n_predictions <- nrow(output)
  n_evaluable <- nrow(evaluable)
  n_correct <- if (n_evaluable > 0L) {
    sum(evaluable$correct_prediction)
  } else {
    0L
  }
  recovered_pairs <- if (nrow(output) > 0L) {
    intersect(unique(output$pair_name), gt$pair_name)
  } else {
    character()
  }
  recovered_H <- unique(gt$Hclone[gt$pair_name %in% recovered_pairs])
  predicted_gt_H <- if (nrow(output) > 0L) {
    intersect(unique(output$Hclone), unique(gt$Hclone))
  } else {
    character()
  }

  status <- if (resolution_input$n_HL_units == 0L) {
    "no_HL_cooccurrence"
  } else if (n_predictions == 0L) {
    "no_predictions"
  } else if (n_evaluable == 0L) {
    "no_evaluable_predictions"
  } else {
    "evaluable"
  }

  summary <- data.frame(
    sample = sample_id,
    resolution = resolution,
    status = status,
    threshold = 1L,
    n_spatial_units = resolution_input$n_units,
    n_HL_units = resolution_input$n_HL_units,
    n_H_clones = resolution_input$n_H_clones,
    n_L_clones = resolution_input$n_L_clones,
    n_active_clones = resolution_input$n_active_clones,
    n_predictions = n_predictions,
    n_evaluable_predictions = n_evaluable,
    n_correct_predictions = n_correct,
    precision = if (n_evaluable > 0L) n_correct / n_evaluable else NA_real_,
    n_fixed_GT_pairs = fixed_gt$n_fixed_pairs,
    n_GT_pairs_recovered = length(recovered_pairs),
    pair_recall = length(recovered_pairs) / fixed_gt$n_fixed_pairs,
    n_fixed_GT_H = fixed_gt$n_fixed_H,
    n_GT_H_recovered = length(recovered_H),
    H_recovery = length(recovered_H) / fixed_gt$n_fixed_H,
    n_GT_H_with_prediction = length(predicted_gt_H),
    H_prediction_coverage = length(predicted_gt_H) / fixed_gt$n_fixed_H,
    stringsAsFactors = FALSE
  )

  if (nrow(output) == 0L) {
    output <- empty_prediction_table_r1q3()
  }
  list(summary = summary, predictions = output)
}

compare_legacy_bin50_r1q3 <- function(sample_id, predictions) {
  legacy_path <- file.path(
    LEGACY_BIN50_DIR, paste0(sample_id, "_final_e100.tsv")
  )
  if (!file.exists(legacy_path)) {
    return(data.frame(
      sample = sample_id,
      legacy_bin50_n_predictions = NA_integer_,
      rerun_bin50_n_predictions = nrow(predictions),
      bin50_pair_set_matches_legacy = NA,
      stringsAsFactors = FALSE
    ))
  }
  legacy <- fread(legacy_path, showProgress = FALSE)
  legacy_pairs <- if ("pair_name" %in% names(legacy)) {
    unique(as.character(legacy$pair_name))
  } else {
    character()
  }
  rerun_pairs <- unique(as.character(predictions$pair_name))
  data.frame(
    sample = sample_id,
    legacy_bin50_n_predictions = length(legacy_pairs),
    rerun_bin50_n_predictions = length(rerun_pairs),
    bin50_pair_set_matches_legacy = setequal(legacy_pairs, rerun_pairs),
    stringsAsFactors = FALSE
  )
}

paired_metric_test_r1q3 <- function(summary_df, metric) {
  wide <- summary_df %>%
    select(sample, resolution, value = all_of(metric)) %>%
    pivot_wider(names_from = resolution, values_from = value) %>%
    filter(is.finite(Bin1), is.finite(Bin50))
  if (nrow(wide) < 2L) {
    return(data.frame(
      metric = metric, n_paired = nrow(wide), statistic = NA_real_,
      p_value = NA_real_, median_Bin1 = NA_real_, median_Bin50 = NA_real_,
      median_paired_difference_Bin50_minus_Bin1 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  test <- tryCatch(
    wilcox.test(
      wide$Bin50, wide$Bin1,
      paired = TRUE, alternative = "two.sided", exact = FALSE,
      correct = TRUE
    ),
    error = function(e) NULL
  )
  data.frame(
    metric = metric,
    n_paired = nrow(wide),
    statistic = if (is.null(test)) NA_real_ else unname(test$statistic),
    p_value = if (is.null(test)) NA_real_ else test$p.value,
    median_Bin1 = median(wide$Bin1),
    median_Bin50 = median(wide$Bin50),
    median_paired_difference_Bin50_minus_Bin1 = median(wide$Bin50 - wide$Bin1),
    stringsAsFactors = FALSE
  )
}

patient_metric_test_r1q3 <- function(summary_df, metric) {
  sample_wide <- summary_df %>%
    select(sample, resolution, value = all_of(metric)) %>%
    pivot_wider(names_from = resolution, values_from = value) %>%
    filter(is.finite(Bin1), is.finite(Bin50))
  complete_samples <- sample_wide$sample
  patient_wide <- summary_df %>%
    filter(sample %in% complete_samples) %>%
    mutate(patient = sub("[BPT]$", "", sample)) %>%
    select(patient, sample, resolution, value = all_of(metric)) %>%
    group_by(patient, resolution) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    pivot_wider(names_from = resolution, values_from = value) %>%
    filter(is.finite(Bin1), is.finite(Bin50))

  if (nrow(patient_wide) < 2L) {
    return(data.frame(
      metric = metric,
      aggregation = "mean_of_complete_paired_sections",
      n_complete_sections = length(complete_samples),
      n_patients = nrow(patient_wide), n_nonzero_differences = 0L,
      statistic = NA_real_, p_value = NA_real_,
      median_patient_mean_Bin1 = NA_real_,
      median_patient_mean_Bin50 = NA_real_,
      median_patient_mean_difference_Bin50_minus_Bin1 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  patient_test <- coin::wilcoxsign_test(
    Bin50 ~ Bin1,
    data = patient_wide,
    distribution = "exact",
    zero.method = "Wilcoxon"
  )
  data.frame(
    metric = metric,
    aggregation = "mean_of_complete_paired_sections",
    n_complete_sections = length(complete_samples),
    n_patients = nrow(patient_wide),
    n_nonzero_differences = sum(patient_wide$Bin50 != patient_wide$Bin1),
    statistic = unname(coin::statistic(patient_test, type = "linear")),
    p_value = unname(coin::pvalue(patient_test)),
    median_patient_mean_Bin1 = median(patient_wide$Bin1),
    median_patient_mean_Bin50 = median(patient_wide$Bin50),
    median_patient_mean_difference_Bin50_minus_Bin1 =
      median(patient_wide$Bin50 - patient_wide$Bin1),
    stringsAsFactors = FALSE
  )
}

summarise_resolution_r1q3 <- function(summary_df) {
  summary_df %>%
    group_by(resolution) %>%
    summarise(
      n_samples = n(),
      n_samples_with_HL_units = sum(n_HL_units > 0L),
      n_samples_with_predictions = sum(n_predictions > 0L),
      n_samples_with_evaluable_precision = sum(is.finite(precision)),
      pooled_correct_predictions = sum(n_correct_predictions),
      pooled_evaluable_predictions = sum(n_evaluable_predictions),
      pooled_precision = ifelse(
        pooled_evaluable_predictions > 0L,
        pooled_correct_predictions / pooled_evaluable_predictions,
        NA_real_
      ),
      pooled_GT_pairs_recovered = sum(n_GT_pairs_recovered),
      pooled_fixed_GT_pairs = sum(n_fixed_GT_pairs),
      pooled_pair_recall = pooled_GT_pairs_recovered / pooled_fixed_GT_pairs,
      pooled_GT_H_recovered = sum(n_GT_H_recovered),
      pooled_fixed_GT_H = sum(n_fixed_GT_H),
      pooled_H_recovery = pooled_GT_H_recovered / pooled_fixed_GT_H,
      median_pair_recall = median(pair_recall),
      q1_pair_recall = unname(quantile(pair_recall, 0.25)),
      q3_pair_recall = unname(quantile(pair_recall, 0.75)),
      median_H_recovery = median(H_recovery),
      median_precision_evaluable = ifelse(
        any(is.finite(precision)), median(precision, na.rm = TRUE), NA_real_
      ),
      q1_precision_evaluable = ifelse(
        any(is.finite(precision)),
        unname(quantile(precision, 0.25, na.rm = TRUE)), NA_real_
      ),
      q3_precision_evaluable = ifelse(
        any(is.finite(precision)),
        unname(quantile(precision, 0.75, na.rm = TRUE)), NA_real_
      ),
      .groups = "drop"
    )
}

format_p_r1q3 <- function(p) {
  if (!is.finite(p)) return("P not estimable")
  if (p < 0.001) {
    return(paste0("P = ", formatC(p, format = "e", digits = 2)))
  }
  paste0("P = ", formatC(p, format = "f", digits = 3))
}

build_pairing_figure_r1q3 <- function(summary_df, tests_df, group_df,
                                      patient_tests_df = NULL) {
  summary_df <- summary_df %>%
    mutate(
      resolution = factor(resolution, levels = RESOLUTION_LEVELS),
      x = as.numeric(resolution)
    )
  recall_test <- tests_df %>% filter(metric == "pair_recall")
  precision_test <- tests_df %>% filter(metric == "precision")
  recall_display_test <- if (is.null(patient_tests_df)) {
    recall_test %>% transmute(n_display = n_paired, p_value)
  } else {
    patient_tests_df %>%
      filter(metric == "pair_recall") %>%
      transmute(n_display = n_patients, p_value)
  }
  precision_display_test <- if (is.null(patient_tests_df)) {
    precision_test %>% transmute(n_display = n_paired, p_value)
  } else {
    patient_tests_df %>%
      filter(metric == "precision") %>%
      transmute(n_display = n_patients, p_value)
  }
  display_unit <- if (is.null(patient_tests_df)) {
    "Paired Wilcoxon"
  } else {
    "Patient-level exact Wilcoxon"
  }

  recall_summary <- summary_df %>%
    group_by(resolution, x) %>%
    summarise(
      median = median(pair_recall * 100),
      q1 = unname(quantile(pair_recall * 100, 0.25)),
      q3 = unname(quantile(pair_recall * 100, 0.75)),
      .groups = "drop"
    )
  recall_pooled <- group_df %>%
    mutate(
      resolution = factor(resolution, levels = RESOLUTION_LEVELS),
      x = as.numeric(resolution),
      pooled_label = sprintf(
        "%s/%s (%.1f%%)",
        format(pooled_GT_pairs_recovered, big.mark = ",", trim = TRUE),
        format(pooled_fixed_GT_pairs, big.mark = ",", trim = TRUE),
        pooled_pair_recall * 100
      )
    )

  p_recall <- ggplot(
    summary_df,
    aes(x = resolution, y = pair_recall * 100, group = sample)
  ) +
    geom_line(color = "#B9C1C9", alpha = 0.45, linewidth = 0.35) +
    geom_point(
      aes(color = resolution), size = 1.85, alpha = 0.88,
      show.legend = FALSE
    ) +
    geom_errorbar(
      data = recall_summary,
      aes(x = resolution, ymin = q1, ymax = q3),
      inherit.aes = FALSE, width = 0.07, linewidth = 0.55,
      color = "#252525"
    ) +
    geom_segment(
      data = recall_summary,
      aes(x = x - 0.13, xend = x + 0.13, y = median, yend = median),
      inherit.aes = FALSE, linewidth = 1.0, color = "#252525",
      lineend = "round"
    ) +
    geom_text(
      data = recall_pooled,
      aes(x = x, y = 104.5, label = pooled_label),
      inherit.aes = FALSE, size = 2.25, color = "#252525"
    ) +
    annotate(
      "text", x = 1.5, y = 113,
      label = paste0(
        display_unit, ", n = ", recall_display_test$n_display, ", ",
        format_p_r1q3(recall_display_test$p_value)
      ),
      size = 2.35, fontface = "plain"
    ) +
    scale_color_manual(values = RESOLUTION_COLORS) +
    scale_x_discrete(labels = RESOLUTION_LABELS) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 116), expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = "B   Ground-truth pair recall",
      subtitle = "Fixed sample-specific denominator; 34 matched samples",
      x = NULL, y = "GT-pair recall (%)"
    ) +
    theme_csl(base_size = 7.5) +
    theme(
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 6.8, color = "#4A4A4A"),
      axis.text.x = element_text(size = 7.0),
      plot.margin = margin(7, 5, 4, 5)
    )

  precision_df <- summary_df %>%
    mutate(precision_y = ifelse(is.finite(precision), precision * 100, -9))
  complete_precision_samples <- precision_df %>%
    filter(is.finite(precision)) %>%
    count(sample) %>%
    filter(n == 2L) %>%
    pull(sample)
  precision_matched <- precision_df %>%
    filter(sample %in% complete_precision_samples, is.finite(precision))
  precision_summary <- precision_matched %>%
    group_by(resolution, x) %>%
    summarise(
      median = median(precision * 100),
      q1 = unname(quantile(precision * 100, 0.25)),
      q3 = unname(quantile(precision * 100, 0.75)),
      .groups = "drop"
    )
  precision_pooled <- precision_matched %>%
    group_by(resolution, x) %>%
    summarise(
      n_correct = sum(n_correct_predictions),
      n_evaluable = sum(n_evaluable_predictions),
      .groups = "drop"
    ) %>%
    mutate(
      pooled_label = ifelse(
        n_evaluable > 0L,
        sprintf(
          "%s/%s (%.1f%%)",
          format(n_correct, big.mark = ",", trim = TRUE),
          format(n_evaluable, big.mark = ",", trim = TRUE),
          n_correct / n_evaluable * 100
        ),
        "NE"
      )
    )
  ne_counts <- precision_df %>%
    group_by(resolution, x) %>%
    summarise(n_ne = sum(!is.finite(precision)), .groups = "drop")

  p_precision <- ggplot(precision_df, aes(x = resolution, group = sample)) +
    geom_line(
      data = precision_df %>%
        filter(sample %in% complete_precision_samples),
      aes(y = precision * 100),
      color = "#B9C1C9", alpha = 0.45, linewidth = 0.35
    ) +
    geom_point(
      data = precision_matched,
      aes(y = precision * 100, color = resolution),
      size = 1.85, alpha = 0.88, show.legend = FALSE
    ) +
    geom_point(
      data = precision_df %>%
        filter(is.finite(precision), !sample %in% complete_precision_samples),
      aes(y = precision * 100, color = resolution),
      shape = 21, fill = "white", size = 1.65, stroke = 0.45,
      alpha = 0.65, show.legend = FALSE
    ) +
    geom_point(
      data = precision_df %>% filter(!is.finite(precision)),
      aes(y = precision_y),
      shape = 4, size = 1.7, stroke = 0.45, color = "#7A7A7A"
    ) +
    geom_errorbar(
      data = precision_summary,
      aes(x = resolution, ymin = q1, ymax = q3),
      inherit.aes = FALSE, width = 0.07, linewidth = 0.55,
      color = "#252525"
    ) +
    geom_segment(
      data = precision_summary,
      aes(x = x - 0.13, xend = x + 0.13, y = median, yend = median),
      inherit.aes = FALSE, linewidth = 1.0, color = "#252525",
      lineend = "round"
    ) +
    geom_text(
      data = precision_pooled,
      aes(x = x, y = 104.5, label = pooled_label),
      inherit.aes = FALSE, size = 2.15, color = "#252525"
    ) +
    geom_text(
      data = ne_counts,
      aes(x = x, y = -9, label = paste0("n=", n_ne)),
      inherit.aes = FALSE, nudge_x = 0.13, size = 2.0,
      color = "#6A6A6A"
    ) +
    annotate(
      "text", x = 1.5, y = 113,
      label = paste0(
        display_unit, ", n = ", precision_display_test$n_display, ", ",
        format_p_r1q3(precision_display_test$p_value)
      ),
      size = 2.25
    ) +
    scale_color_manual(values = RESOLUTION_COLORS) +
    scale_x_discrete(labels = RESOLUTION_LABELS) +
    scale_y_continuous(
      breaks = c(-9, 0, 25, 50, 75, 100),
      labels = c("NE", "0", "25", "50", "75", "100"),
      limits = c(-13, 116), expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = "C   Prediction precision",
      subtitle = "Matched 11 sections; 10 patients",
      x = NULL, y = "Precision (%)"
    ) +
    theme_csl(base_size = 7.5) +
    theme(
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 6.8, color = "#4A4A4A"),
      axis.text.x = element_text(size = 7.0),
      plot.margin = margin(7, 5, 4, 5)
    )

  p_recall + p_precision + plot_layout(widths = c(1.65, 1))
}

write_stats_log_r1q3 <- function(path, summary_df, group_df, tests_df,
                                  patient_tests_df, audit_df) {
  lines <- c(
    "R1_Q3 matched Bin1/Bin50 SIPEAC comparison",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "Script: Code/Extended_Analyses/03_Bin1_Bin50_pairing.R",
    paste0("R: ", R.version.string),
    paste0(
      "Packages: data.table ", as.character(packageVersion("data.table")),
      "; dplyr ", as.character(packageVersion("dplyr")),
      "; tidyr ", as.character(packageVersion("tidyr")),
      "; ggplot2 ", as.character(packageVersion("ggplot2")),
      "; coin ", as.character(packageVersion("coin"))
    ),
    paste0("Samples: ", n_distinct(summary_df$sample)),
    "Clone-occurrence threshold: 1 (all observed clones)",
    "Statistical unit: sample",
    paste0(
      "Fixed GT pairs across samples: ",
      sum(summary_df$n_fixed_GT_pairs[summary_df$resolution == "Bin1"])
    ),
    paste0(
      "Fixed GT H clones across samples: ",
      sum(summary_df$n_fixed_GT_H[summary_df$resolution == "Bin1"])
    ),
    "Primary endpoint: fixed-GT pair recall.",
    "Precision is conditional on evaluable predictions and is NA when its denominator is zero.",
    "A valid run with no prediction has recall = 0, not missing.",
    "Bin50 is an approximately 25-µm spatial aggregate, not 25 µm^2.",
    "",
    "Resolution summary:",
    capture.output(print(group_df, row.names = FALSE)),
    "",
    "Paired sample-level tests:",
    capture.output(print(tests_df, row.names = FALSE)),
    "",
    "Patient-aggregated exact paired Wilcoxon sensitivity tests:",
    "B/P/T specimens from the same patient are averaged after restricting to complete section pairs.",
    capture.output(print(patient_tests_df, row.names = FALSE)),
    "",
    "Status counts:",
    capture.output(print(
      summary_df %>% count(resolution, status), row.names = FALSE
    )),
    "",
    "Input and legacy-reproduction audit:",
    capture.output(print(audit_df, row.names = FALSE))
  )
  writeLines(lines, path, useBytes = TRUE)
}

write_caption_r1q3 <- function(path, summary_df, tests_df, patient_tests_df,
                                group_df) {
  recall_test <- tests_df[tests_df$metric == "pair_recall", , drop = FALSE]
  precision_test <- tests_df[tests_df$metric == "precision", , drop = FALSE]
  recall_patient_test <- patient_tests_df[
    patient_tests_df$metric == "pair_recall", , drop = FALSE
  ]
  precision_patient_test <- patient_tests_df[
    patient_tests_df$metric == "precision", , drop = FALSE
  ]
  b1 <- group_df[group_df$resolution == "Bin1", , drop = FALSE]
  b50 <- group_df[group_df$resolution == "Bin50", , drop = FALSE]
  complete_precision_samples <- summary_df %>%
    filter(is.finite(precision)) %>%
    count(sample) %>%
    filter(n == 2L) %>%
    pull(sample)
  matched_precision <- summary_df %>%
    filter(sample %in% complete_precision_samples) %>%
    group_by(resolution) %>%
    summarise(
      n_correct = sum(n_correct_predictions),
      n_evaluable = sum(n_evaluable_predictions),
      precision = n_correct / n_evaluable,
      .groups = "drop"
    )
  mp1 <- matched_precision[matched_precision$resolution == "Bin1", ]
  mp50 <- matched_precision[matched_precision$resolution == "Bin50", ]
  caption <- sprintf(
    paste0(
      "**SIPEAC pairing performance across spatial resolutions.** ",
      "**(A)** Schematic comparison of direct same-bin1 enumeration and ",
      "SIPEAC evidence integration. Direct enumeration retains all H–L ",
      "combinations co-detected in one DNB-level unit, whereas SIPEAC scores ",
      "spatial co-occurrence and CDR3 amino-acid similarity before applying a ",
      "one-to-one filter. This schematic explains candidate generation and is ",
      "not itself a quantitative Bin1–Bin50 comparison. ",
      "**(B)** Ground-truth pair recall at Bin1 and Bin50 in 34 matched ",
      "samples. The eligible single-cell BCR-supported H–L pair set was ",
      "defined once per sample and held fixed across resolutions. Grey lines ",
      "connect the same sample; coloured points show sample-level estimates; ",
      "black bars show the median and interquartile range. Descriptively, ",
      "SIPEAC recovered %d/%d GT pairs (%.1f%%) at Bin1 and %d/%d (%.1f%%) ",
      "at Bin50. The displayed P value is from a two-sided exact paired ",
      "Wilcoxon signed-rank test after equally averaging B/P/T specimens ",
      "within each of %d patients (%s). A valid analysis with no ",
      "predicted pair contributes recall zero. **(C)** Prediction precision ",
      "is the proportion of evaluable predictions supported by ground truth. ",
      "Within the 11 matched sections, it was %d/%d (%.1f%%) at Bin1 and ",
      "%d/%d (%.1f%%) at Bin50. Precision is not estimable when there is no ",
      "evaluable prediction; these observations are shown in the NE track and ",
      "excluded from the paired analysis. Open points denote finite estimates ",
      "without a Bin1 match. The displayed precision test averages the 11 ",
      "complete sections within %d patients (%s). Bin50 denotes an ",
      "approximately 25-µm spatial aggregate."
    ),
    b1$pooled_GT_pairs_recovered, b1$pooled_fixed_GT_pairs,
    100 * b1$pooled_pair_recall,
    b50$pooled_GT_pairs_recovered, b50$pooled_fixed_GT_pairs,
    100 * b50$pooled_pair_recall,
    recall_patient_test$n_patients,
    format_p_r1q3(recall_patient_test$p_value),
    mp1$n_correct, mp1$n_evaluable, 100 * mp1$precision,
    mp50$n_correct, mp50$n_evaluable, 100 * mp50$precision,
    precision_patient_test$n_patients,
    format_p_r1q3(precision_patient_test$p_value)
  )
  writeLines(caption, path, useBytes = TRUE)
}
