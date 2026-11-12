PROJECT_ROOT <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, 'data')
source('Code/Extended_Analyses/utils_io.R')
source('Code/Extended_Analyses/utils_plot.R')
source('Code/Extended_Analyses/01_spatial_validation_functions.R')

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(tibble)

SIPEAC_ROOT <- dirs$sipeac
METADATA_DIR <- dirs$st_meta
processed_dir <- file.path(dirs$processed, '01_Spatial_pairing')
table_dir <- file.path(dirs$tables, '01_Spatial_pairing')
log_dir <- file.path(dirs$stats_logs, '01_Spatial_pairing')
figure_dir <- file.path(dirs$figures, '01_Spatial_pairing')

#### 输入与样本核对 ====

BIN_SIZE <- 50L

N_PERM <- 10000L

SEED <- 20260725L

EXCLUDED_SAMPLES <- c("4266B", "4309B")

if (!is.finite(N_PERM) || N_PERM < 99L) {
  stop("N_PERM must be at least 99", call. = FALSE)
}

input_dir <- file.path(SIPEAC_ROOT, "input_allelic", paste0("bin", BIN_SIZE))

prediction_dir <- file.path(SIPEAC_ROOT, "final_version", paste0("bin", BIN_SIZE))

truth_dir <- file.path(SIPEAC_ROOT, "grouth_truth", paste0("bin", BIN_SIZE))

compare_path <- file.path(SIPEAC_ROOT, "final_version", "compare.csv")

required_roots <- c(input_dir, prediction_dir, truth_dir, METADATA_DIR, compare_path)

missing_roots <- required_roots[!file.exists(required_roots)]

if (length(missing_roots)) {
  stop("Missing input path(s):\n", paste(missing_roots, collapse = "\n"), call. = FALSE)
}

compare <- fread(compare_path)

if (!all(c("bin_size", "cur_sample") %in% names(compare))) {
  stop("compare.csv must contain bin_size and cur_sample", call. = FALSE)
}

samples <- sort(setdiff(
  unique(as.character(compare[bin_size == BIN_SIZE, cur_sample])),
  EXCLUDED_SAMPLES
))

sample_paths <- lapply(samples, function(sample_name) {
  c(
    info = file.path(input_dir, paste0(sample_name, "_info.txt")),
    final = file.path(prediction_dir, paste0(sample_name, "_final_e100.tsv")),
    test = file.path(prediction_dir, paste0(sample_name, "_test_e100.tsv")),
    truth = file.path(truth_dir, paste0(sample_name, ".tsv")),
    metadata = file.path(METADATA_DIR, paste0(sample_name, ".txt"))
  )
})

names(sample_paths) <- samples

missing_inputs <- unlist(lapply(sample_paths, function(x) x[!file.exists(x)]), use.names = FALSE)

if (length(missing_inputs)) {
  stop("Missing sample input(s):\n", paste(missing_inputs, collapse = "\n"), call. = FALSE)
}

pair_counts <- vapply(samples, function(sample_name) {
  nrow(fread(sample_paths[[sample_name]][["test"]], select = c("Hclone", "Lclone")))
}, integer(1))

if (length(samples) != 30L || sum(pair_counts) != 1102L) {
  stop(
    sprintf(
      "The locked analysis expects 30 samples and 1,102 pairs; found %d and %d.",
      length(samples), sum(pair_counts)
    ), call. = FALSE
  )
}

input_correct <- 0L

for (sample_name in samples) {
  paths <- sample_paths[[sample_name]]
  info_names <- names(fread(paths[["info"]], nrows = 0L))
  metadata_names <- names(fread(paths[["metadata"]], nrows = 0L))
  if (!all(c("type", "locus", "clone_id", "binid") %in% info_names)) {
    stop("Unexpected info columns for ", sample_name, call. = FALSE)
  }
  if (!all(c("BinID", "row", "col", "tissue") %in% metadata_names)) {
    stop("Unexpected metadata columns for ", sample_name, call. = FALSE)
  }
  final <- fread(paths[["final"]])
  final <- final[!is.na(Lclone_t), .(H = Hclone, L = Lclone)]
  test <- fread(paths[["test"]], select = c("Hclone", "Lclone"))
  truth <- read_truth(paths[["truth"]])
  if (anyDuplicated(final) || anyDuplicated(test) ||
      nrow(final) != nrow(test) ||
      !setequal(paste(final$H, final$L), paste(test$Hclone, test$Lclone))) {
    stop("final_e100 and test_e100 disagree for ", sample_name, call. = FALSE)
  }
  truth_key <- paste(truth$Hclone, truth$Lclone, sep = "\r")
  input_correct <- input_correct + sum(
    paste(test$Hclone, test$Lclone, sep = "\r") %in% truth_key
  )
}

if (input_correct != 581L) stop("Correctness count check failed", call. = FALSE)

message("Inputs: 30 samples, 1,102 predictions; 581 correct, 521 incorrect")

invisible(lapply(c(processed_dir, table_dir, log_dir, figure_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

#### 重复空间事件 ====

RNGversion("4.2.2")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

set.seed(SEED)

pair_results <- vector("list", length(samples))

names(pair_results) <- samples

for (sample_name in samples) {
  message("Tier 1: ", sample_name)
  dat <- read_sample(sample_name)
  clone_bins <- split(dat$spatial$binid, dat$spatial$clone_name)
  clone_bins <- lapply(clone_bins, unique)
  all_bins <- unique(dat$spatial$binid)
  result <- vector("list", nrow(dat$pairs))

  for (i in seq_len(nrow(dat$pairs))) {
    h <- dat$pairs$H[[i]]
    l <- dat$pairs$L[[i]]
    h_bins <- unique(as.character(clone_bins[[h]]))
    l_bins <- unique(as.character(clone_bins[[l]]))
    if (!length(h_bins) || !length(l_bins)) {
      stop("Missing spatial clone for ", sample_name, ": ", h, " x ", l, call. = FALSE)
    }
    h_coords <- dat$metadata[match(h_bins, binid), .(row, col)]
    l_coords <- dat$metadata[match(l_bins, binid), .(row, col)]
    overlap_bins <- intersect(h_bins, l_bins)
    overlap_index <- match(overlap_bins, h_bins)
    adjacency <- setNames(lapply(1:3, function(radius) {
      make_adjacency(h_coords, radius)
    }), paste0("r", 1:3))
    observed <- vapply(adjacency, function(x) {
      count_components(overlap_index, x)
    }, integer(1))
    null <- repeat_null(
      length(h_bins), length(l_bins), length(all_bins), adjacency,
      observed, N_PERM
    )
    k_h <- length(h_bins)
    k_l <- length(l_bins)
    o <- length(overlap_bins)
    expected_overlap <- k_h * k_l / length(all_bins)
    centroid_distance <- sqrt(
      (mean(h_coords$row) - mean(l_coords$row))^2 +
        (mean(h_coords$col) - mean(l_coords$col))^2
    )
    result[[i]] <- data.frame(
      sample_name = sample_name, H = h, L = l,
      pair_id = paste(h, l, sep = "_"),
      K_H = k_h, K_L = k_l, N = length(all_bins), O_ij = o,
      E_ij = expected_overlap,
      p_hyper = phyper(o - 1L, m = k_h, n = length(all_bins) - k_h,
                       k = k_l, lower.tail = FALSE),
      R_r1 = observed[["r1"]], R_r2 = observed[["r2"]], R_r3 = observed[["r3"]],
      E_R_r1 = null$expected[["r1"]], E_R_r2 = null$expected[["r2"]],
      E_R_r3 = null$expected[["r3"]],
      p_repeat_r1 = null$p[["r1"]], p_repeat_r2 = null$p[["r2"]],
      p_repeat_r3 = null$p[["r3"]],
      centroid_distance = centroid_distance,
      stringsAsFactors = FALSE
    )
  }
  pair_results[[sample_name]] <- bind_rows(result)
}

repeated <- bind_rows(pair_results) %>%
  group_by(sample_name) %>%
  mutate(
    q_repeat_r1 = p.adjust(p_repeat_r1, method = "BH"),
    q_repeat_r2 = p.adjust(p_repeat_r2, method = "BH"),
    q_repeat_r3 = p.adjust(p_repeat_r3, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    support_r1 = R_r1 >= 2L & q_repeat_r1 < 0.05,
    support_r2 = R_r2 >= 2L & q_repeat_r2 < 0.05,
    support_r3 = R_r3 >= 2L & q_repeat_r3 < 0.05,
    repeat_excess_r2 = R_r2 - E_R_r2,
    repeat_ratio_r2 = (R_r2 + 0.5) / (E_R_r2 + 0.5)
  )

prediction_truth <- bind_rows(lapply(samples, function(sample_name) {
  test <- fread(
    sample_paths[[sample_name]][["test"]],
    select = c("Hclone", "Lclone")
  )
  truth <- read_truth(sample_paths[[sample_name]][["truth"]])
  truth_key <- paste(truth$Hclone, truth$Lclone, sep = "\r")
  data.frame(
    sample_name = sample_name,
    H = test$Hclone,
    L = test$Lclone,
    correct = paste(test$Hclone, test$Lclone, sep = "\r") %in% truth_key,
    stringsAsFactors = FALSE
  )
}))

if (anyDuplicated(prediction_truth[c("sample_name", "H", "L")])) {
  stop("Duplicated prediction key", call. = FALSE)
}

repeated <- repeated %>%
  inner_join(prediction_truth, by = c("sample_name", "H", "L"))

if (nrow(repeated) != 1102L || sum(repeated$correct) != 581L) {
  stop("Prediction/ground-truth check failed", call. = FALSE)
}

#### Ridge交叉验证 ====

analysis_data <- repeated %>%
  mutate(
    correct = as.logical(correct),
    tier1_support = as.logical(support_r2),
    log_K_H = log1p(K_H),
    log_K_L = log1p(K_L),
    log_expected_overlap = log1p(E_ij),
    log_observed_overlap = log1p(O_ij),
    log_repeat_ratio_r2 = log1p(repeat_ratio_r2),
    neglog10_p_repeat_r2 = -log10(pmax(p_repeat_r2, 1e-8)),
    log_centroid_distance = log1p(centroid_distance),
    calibration_probability = NA_real_
  )

feature_names <- c(
  "log_K_H", "log_K_L", "log_expected_overlap", "log_observed_overlap",
  "R_r1", "R_r2", "R_r3", "repeat_excess_r2",
  "log_repeat_ratio_r2", "neglog10_p_repeat_r2", "log_centroid_distance"
)

lambda_grid <- 10^seq(-4, 2, length.out = 13L)

outer_samples <- sort(unique(analysis_data$sample_name))

outer_folds <- vector("list", length(outer_samples))

for (outer_index in seq_along(outer_samples)) {
  held_out <- outer_samples[[outer_index]]
  message(sprintf("Tier 2 outer-LOSO %02d/%02d: %s",
                  outer_index, length(outer_samples), held_out))
  train <- analysis_data$sample_name != held_out & !analysis_data$tier1_support
  test <- analysis_data$sample_name == held_out & !analysis_data$tier1_support
  if (!any(test)) {
    outer_folds[[outer_index]] <- data.frame(
      held_out_sample = held_out, n_test_pairs = 0L, selected_lambda = NA_real_
    )
    next
  }
  x_train <- as.matrix(analysis_data[train, feature_names])
  x_test <- as.matrix(analysis_data[test, feature_names])
  y_train <- as.integer(analysis_data$correct[train])
  sample_train <- analysis_data$sample_name[train]
  selected <- select_lambda(
    x_train, y_train, sample_train, lambda_grid, 5L, SEED + outer_index
  )
  scaled <- standardize_train_test(x_train, x_test)
  fit <- fit_ridge(
    scaled$train, y_train, equal_sample_weights(sample_train), selected$lambda
  )
  if (!selected$converged || fit$convergence != 0L) {
    stop("Ridge fit did not converge in outer fold ", held_out, call. = FALSE)
  }
  analysis_data$calibration_probability[test] <- predict_ridge(fit, scaled$test)
  outer_folds[[outer_index]] <- data.frame(
    held_out_sample = held_out, n_test_pairs = sum(test),
    selected_lambda = selected$lambda, inner_mean_sample_log_loss = selected$loss
  )
}

if (anyNA(analysis_data$calibration_probability[!analysis_data$tier1_support])) {
  stop("Missing outer-LOSO probabilities", call. = FALSE)
}

analysis_data <- analysis_data %>%
  mutate(
    tier2_rescue = !tier1_support & calibration_probability >= 0.5,
    final_support = tier1_support | tier2_rescue,
    decision_source = case_when(
      tier1_support ~ "Tier 1 repeated spatial evidence",
      tier2_rescue ~ "Tier 2 calibrated rescue",
      TRUE ~ "No support"
    )
  )

locked_run <- N_PERM == 10000L && SEED == 20260725L

tier1_counts <- with(analysis_data, c(
  sum(tier1_support & correct), sum(tier1_support & !correct),
  sum(!tier1_support & correct), sum(!tier1_support & !correct)
))

final_counts <- with(analysis_data, c(
  sum(final_support & correct), sum(final_support & !correct),
  sum(!final_support & correct), sum(!final_support & !correct)
))

if (locked_run && (!identical(as.integer(tier1_counts), c(339L, 37L, 242L, 484L)) ||
                   !identical(as.integer(final_counts), c(459L, 146L, 122L, 375L)))) {
  stop("Locked count check failed", call. = FALSE)
}

#### 代表性配对 ====

analysis_data <- analysis_data %>%
  mutate(example_group = case_when(
    tier1_support & correct ~ "Tier 1 support / Correct",
    tier1_support & !correct ~ "Tier 1 support / Incorrect",
    tier2_rescue & correct ~ "Tier 2 rescue / Correct",
    tier2_rescue & !correct ~ "Tier 2 rescue / Incorrect",
    !final_support & correct ~ "Final no support / Correct",
    TRUE ~ "Final no support / Incorrect"
  ))

medoid_features <- c(
  "log_K_H", "log_K_L", "log_observed_overlap", "R_r2",
  "repeat_excess_r2", "log_repeat_ratio_r2",
  "neglog10_p_repeat_r2", "log_centroid_distance"
)

feature_z <- scale(analysis_data[, medoid_features])

analysis_data$distance_to_group_median <- NA_real_

for (group_name in unique(analysis_data$example_group)) {
  index <- which(analysis_data$example_group == group_name)
  group_median <- apply(feature_z[index, , drop = FALSE], 2, median)
  analysis_data$distance_to_group_median[index] <- rowSums(
    sweep(feature_z[index, , drop = FALSE], 2, group_median, "-")^2
  )
}

representatives <- analysis_data %>%
  group_by(example_group) %>%
  slice_min(distance_to_group_median, n = 1L, with_ties = FALSE) %>%
  ungroup()

example_order <- c(
  "Tier 1 support / Correct", "Tier 2 rescue / Correct",
  "Final no support / Correct"
)

examples <- representatives %>%
  filter(example_group %in% example_order) %>%
  mutate(example_group = factor(example_group, levels = example_order)) %>%
  arrange(example_group)

if (locked_run) {
  expected_examples <- c(
    "2931B\rIGH_139353\rIGL_144611",
    "2931P\rIGH_59861\rIGK_170597",
    "RH16B\rIGH_65508\rIGK_159607"
  )
  observed_examples <- paste(examples$sample_name, examples$H, examples$L, sep = "\r")
  if (!identical(observed_examples, expected_examples)) {
    stop("Representative-pair check failed", call. = FALSE)
  }
}

#### 作图 ====

font_scale <- 1.18

single_h <- data.frame(x = c(3, 4, 5, 3, 4, 5), y = c(3, 3, 3, 4, 4, 4))

single_l <- data.frame(x = c(4, 5, 4, 5, 7, 8), y = c(3, 3, 4, 4, 6, 6))

repeat_h <- data.frame(
  x = c(2, 3, 2, 6, 7, 6, 9, 8, 9),
  y = c(2, 2, 3, 3, 3, 4, 7, 8, 8)
)

repeat_l <- data.frame(
  x = c(2, 3, 6, 7, 9, 8, 4, 5),
  y = c(2, 3, 3, 4, 8, 7, 6, 6)
)

schematic <- bind_rows(
  make_schematic("One local event (R2 = 1)\nTier 1: no support", single_h, single_l),
  make_schematic("Three separated events (R2 = 3)\nTier 1: support if q < 0.05",
                 repeat_h, repeat_l)
)

schematic$scenario <- factor(
  schematic$scenario,
  levels = c("One local event (R2 = 1)\nTier 1: no support",
             "Three separated events (R2 = 3)\nTier 1: support if q < 0.05")
)

plot_schematic <- ggplot(schematic, aes(x, y, fill = state)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~scenario, nrow = 1) +
  scale_fill_manual(
    values = c("No signal" = "#EEEEEE", "H only" = "#D73027",
               "L only" = "#4575B4", "H–L overlap" = "#1A9850"),
    breaks = c("H only", "L only", "H–L overlap", "No signal")
  ) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Repeated versus single overlap",
    subtitle = "Radius-2-connected bins form one event", x = NULL, y = NULL, fill = NULL
  ) +
  theme_void(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold", size = 8.3 * font_scale),
    legend.position = "bottom", legend.text = element_text(size = 7.5 * font_scale),
    legend.key.size = grid::unit(8, "pt"),
    plot.title = element_text(face = "bold", size = 10 * font_scale),
    plot.subtitle = element_text(size = 8 * font_scale, color = "#4D4D4D"),
    plot.margin = margin(5, 5, 2, 5)
  )

flow_counts <- c(
  all = nrow(analysis_data), tier1 = sum(analysis_data$tier1_support),
  tier1_correct = sum(analysis_data$tier1_support & analysis_data$correct),
  tier1_incorrect = sum(analysis_data$tier1_support & !analysis_data$correct),
  unsupported = sum(!analysis_data$tier1_support),
  rescue = sum(analysis_data$tier2_rescue),
  rescue_correct = sum(analysis_data$tier2_rescue & analysis_data$correct),
  rescue_incorrect = sum(analysis_data$tier2_rescue & !analysis_data$correct),
  no_support = sum(!analysis_data$final_support),
  no_support_correct = sum(!analysis_data$final_support & analysis_data$correct),
  no_support_incorrect = sum(!analysis_data$final_support & !analysis_data$correct),
  final_support = sum(analysis_data$final_support),
  final_correct = sum(analysis_data$final_support & analysis_data$correct),
  final_incorrect = sum(analysis_data$final_support & !analysis_data$correct)
)

flow_nodes <- tibble::tribble(
  ~id, ~x, ~y, ~width, ~height, ~fill, ~label,
  "input", 1.0, 2.0, 1.7, 1.05, "#F2F2F2",
  sprintf("Final SIPEAC\npredictions\nn = %d", flow_counts[["all"]]),
  "tier1", 3.4, 2.9, 2.35, 1.35, "#C6DBEF",
  sprintf("Tier 1 support\nrepeated-event null\n%d: %d correct / %d incorrect",
          flow_counts[["tier1"]], flow_counts[["tier1_correct"]],
          flow_counts[["tier1_incorrect"]]),
  "unsupported", 3.4, 1.05, 2.35, 1.05, "#E5E5E5",
  sprintf("Tier 1 unsupported\nn = %d", flow_counts[["unsupported"]]),
  "tier2", 6.2, 1.55, 2.45, 1.35, "#DADAEB",
  sprintf("Tier 2 ridge rescue\nouter-LOSO, P >= 0.5\n%d: %d correct / %d incorrect",
          flow_counts[["rescue"]], flow_counts[["rescue_correct"]],
          flow_counts[["rescue_incorrect"]]),
  "final_no", 8.9, 0.45, 2.25, 1.05, "#D9D9D9",
  sprintf("Final no support\n%d: %d correct / %d incorrect",
          flow_counts[["no_support"]], flow_counts[["no_support_correct"]],
          flow_counts[["no_support_incorrect"]]),
  "final_yes", 8.9, 2.75, 2.25, 1.05, "#9ECDC1",
  sprintf("Final support\nTier 1 OR Tier 2\n%d: %d correct / %d incorrect",
          flow_counts[["final_support"]], flow_counts[["final_correct"]],
          flow_counts[["final_incorrect"]])
) %>%
  mutate(xmin = x - width / 2, xmax = x + width / 2,
         ymin = y - height / 2, ymax = y + height / 2)

flow_edges <- tibble::tribble(
  ~x, ~y, ~xend, ~yend, ~label, ~label_x, ~label_y,
  1.85, 2.15, 2.22, 2.72, "pass", 2.07, 2.62,
  1.85, 1.85, 2.22, 1.25, "fail", 2.07, 1.40,
  4.58, 1.05, 4.98, 1.43, "", 4.80, 1.42,
  7.43, 1.82, 7.78, 2.55, "rescue", 7.65, 2.32,
  7.43, 1.27, 7.78, 0.65, "not rescued", 7.70, 0.87,
  4.58, 3.00, 7.78, 2.86, "", 6.20, 3.10
)

plot_flow <- ggplot() +
  geom_segment(
    data = flow_edges, aes(x, y, xend = xend, yend = yend),
    color = "#555555", linewidth = 0.55,
    arrow = grid::arrow(length = grid::unit(5, "pt"), type = "closed")
  ) +
  geom_rect(
    data = flow_nodes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    color = "#3F3F3F", linewidth = 0.4
  ) +
  geom_text(data = flow_nodes, aes(x, y, label = label),
            size = 2.65 * font_scale, lineheight = 0.95) +
  geom_text(data = filter(flow_edges, nzchar(label)),
            aes(label_x, label_y, label = label),
            size = 2.45 * font_scale, color = "#4D4D4D") +
  annotate("text", x = 3.4, y = 3.83,
           label = "Primary evidence: no ground truth",
           size = 2.65 * font_scale, color = "#08519C", fontface = "bold") +
  annotate("text", x = 6.2, y = 0.50,
           label = "Secondary: current sample held out",
           size = 2.55 * font_scale, color = "#6A51A3", fontface = "bold") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 10.1), ylim = c(-0.15, 4.05), clip = "off") +
  labs(title = "Two-stage rule", subtitle = "Support is an evidence flag") +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10 * font_scale),
    plot.subtitle = element_text(size = 8 * font_scale, color = "#4D4D4D"),
    plot.margin = margin(3, 5, 6, 5)
  )

panel_a <- (plot_schematic / plot_flow + plot_layout(heights = c(0.95, 1.05))) +
  plot_annotation(
    title = "A  Model overview",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12 * font_scale),
      plot.background = element_rect(fill = "#FCFCFC", color = "#D9D9D9", linewidth = 0.45),
      plot.margin = margin(3, 3, 3, 3)
    )
  )

tier1_fisher <- fisher.test(matrix(tier1_counts, nrow = 2, byrow = TRUE))

final_fisher <- fisher.test(matrix(final_counts, nrow = 2, byrow = TRUE))

panel_b <- make_confusion_plot(
  tier1_counts, "B  Tier 1",
  sprintf("Primary | P %s | OR %.3f | PPV %.1f%%",
          format_p(tier1_fisher$p.value), unname(tier1_fisher$estimate),
          100 * tier1_counts[[1]] / sum(tier1_counts[1:2])),
  "Repeated-event evidence"
)

panel_c <- make_confusion_plot(
  final_counts, "C  Tier 1 + ridge",
  sprintf("Secondary | outer-LOSO | P %s | OR %.3f | PPV %.1f%%",
          format_p(final_fisher$p.value), unname(final_fisher$estimate),
          100 * final_counts[[1]] / sum(final_counts[1:2])),
  "Two-stage spatial evidence"
)

signal_colors <- c("H only" = "#D73027", "L only" = "#4575B4",
                   "H–L overlap" = "#1A9850")

panel_d <- make_spatial_panel(examples[1, ], "D")

panel_e <- make_spatial_panel(examples[2, ], "E")

panel_f <- make_spatial_panel(examples[3, ], "F")

combined <- (
  wrap_elements(full = panel_a) + wrap_elements(full = panel_b) +
    wrap_elements(full = panel_c) + wrap_elements(full = panel_d) +
    wrap_elements(full = panel_e) + wrap_elements(full = panel_f) +
    plot_layout(
      design = "\nAAAABB\nAAAACC\nDDEEFF\n",
      widths = rep(1, 6), heights = c(1, 1, 1.15)
    )
) +
  plot_annotation(
    title = "Spatial validation of SIPEAC H–L pairs",
    subtitle = "Tier 1: primary test | Tier 2: outer-LOSO calibration",
    caption = paste0(
      "D–F are locked representative correct SIPEAC predictions reconstructed from raw bin50 coordinates. ",
      "Arrows connect each dashed region to a side-by-side view of the same region; E1–E3 mark ",
      "Chebyshev-radius-2 overlap components. The zoom does not represent a second occurrence. ",
      "These panels show model–data consistency, not independent molecular proof of H–L identity."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15 * font_scale),
      plot.subtitle = element_text(size = 10 * font_scale),
      plot.caption = element_text(size = 8 * font_scale, color = "#4D4D4D", hjust = 0),
      plot.margin = margin(3, 3, 3, 3)
    )
  )

panel_paths <- c(
  A = "panel_A_model_overview", B = "panel_B_tier1",
  C = "panel_C_two_stage", D = "panel_D_tier1_in_situ",
  E = "panel_E_tier2_rescue_in_situ", F = "panel_F_no_support_in_situ"
)

panel_plots <- list(A = panel_a, B = panel_b, C = panel_c,
                    D = panel_d, E = panel_e, F = panel_f)

panel_sizes <- list(A = c(8.2, 7.6), B = c(5.8, 4.3), C = c(5.8, 4.3),
                    D = c(7.6, 4.6), E = c(7.6, 4.6), F = c(7.6, 4.6))

for (panel_name in names(panel_plots)) {
  save_figure(
    panel_plots[[panel_name]], file.path(figure_dir, panel_paths[[panel_name]]),
    panel_sizes[[panel_name]][[1]], panel_sizes[[panel_name]][[2]]
  )
}

combined_path <- file.path(figure_dir, "R1_Q2a_spatial_validation_from_raw")

save_figure(combined, combined_path, 18, 11.2)

#### 输出 ====

counts_table <- bind_rows(
  data.frame(stage = "Tier 1", status = c("Support", "Support", "No support", "No support"),
             correctness = c("Correct", "Incorrect", "Correct", "Incorrect"),
             n = tier1_counts),
  data.frame(stage = "Final", status = c("Support", "Support", "No support", "No support"),
             correctness = c("Correct", "Incorrect", "Correct", "Incorrect"),
             n = final_counts)
)

example_table <- examples %>%
  select(example_group, sample_name, H, L, correct, K_H, K_L, O_ij, R_r2,
         q_repeat_r2, calibration_probability, tier1_support, tier2_rescue,
         final_support, centroid_distance)

write.csv(repeated, file.path(table_dir, "R1_Q2a_repeated_spatial_from_raw.csv"),
          row.names = FALSE)

write.csv(analysis_data,
          file.path(table_dir, "R1_Q2a_two_stage_predictions_from_raw.csv"),
          row.names = FALSE)

write.csv(bind_rows(outer_folds),
          file.path(table_dir, "R1_Q2a_outer_loso_folds_from_raw.csv"),
          row.names = FALSE)

write.csv(counts_table,
          file.path(table_dir, "R1_Q2a_spatial_summary_counts_from_raw.csv"),
          row.names = FALSE)

write.csv(example_table,
          file.path(table_dir, "R1_Q2a_spatial_summary_examples_from_raw.csv"),
          row.names = FALSE)

saveRDS(
  list(repeated = repeated, predictions = analysis_data,
       representatives = representatives,
       parameters = list(
         project_root = PROJECT_ROOT,
         sipeac_root = SIPEAC_ROOT,
         metadata_dir = METADATA_DIR,
         bin_size = BIN_SIZE,
         n_perm = N_PERM,
         seed = SEED,
         excluded_samples = EXCLUDED_SAMPLES
       )),
  file.path(processed_dir, "R1_Q2a_spatial_validation_from_raw.rds")
)

log_path <- file.path(log_dir, "R1_Q2a_spatial_validation_from_raw_sessionInfo.txt")

sink(log_path)

cat("Samples: ", length(samples), "\n", sep = "")

cat("Pairs: ", nrow(analysis_data), "\n", sep = "")

cat("N_PERM: ", N_PERM, "\n", sep = "")

cat("Seed: ", SEED, "\n", sep = "")

cat("Tier 1 counts: ", paste(tier1_counts, collapse = ", "), "\n", sep = "")

cat("Final counts: ", paste(final_counts, collapse = ", "), "\n\n", sep = "")

print(sessionInfo())

sink()

message("Saved: ", paste0(combined_path, ".pdf"))

message("Saved: ", paste0(combined_path, ".png"))

message("Saved independent panels to: ", figure_dir)
