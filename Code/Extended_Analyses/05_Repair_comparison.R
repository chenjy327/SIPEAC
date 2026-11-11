#### Packages and paths ####
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/R/io.R")
source("Code/Extended_Analyses/R/plot.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
source("Code/Extended_Analyses/R/05_Repair_functions.R")

#### Input ####
sipeac_root <- file.path(dirs$sipeac, "final_version")
ground_truth_root <- file.path(dirs$sipeac, "grouth_truth")
repair_default_root <- file.path(dirs$sipeac, "repair_result")
repair_b100_root <- file.path(dirs$sipeac, "repair_result_b100")
repair_recovery_dir <- file.path(data_dir, "processed", "R4_Q7_missing_result_recovery", "repair_b100")
repair_recovery_manifest <- file.path(PROJECT_ROOT, "results", "stats_logs",
  "R4_Q7_missing_result_recovery", "recovery_provenance.csv")
figure_dir <- file.path(dirs$figures, "R4_Q7_AB_reproducible")
table_dir <- file.path(dirs$tables, "R4_Q7_AB_reproducible")
historical_repair <- TRUE

roots <- list(sipeac = sipeac_root, ground_truth = ground_truth_root,
  repair_default = repair_default_root, repair_b100 = repair_b100_root)
roots <- lapply(roots, normalizePath, mustWork = TRUE)
register_repair_recovery(roots, repair_recovery_dir, repair_recovery_manifest)
comparison_file <- file.path(sipeac_root, "compare.csv")
comparison <- if (file.exists(comparison_file)) read.csv(comparison_file) else NULL
cohort_samples <- if (!is.null(comparison)) unique(comparison$cur_sample) else NULL
patient_map <- if (!is.null(comparison) && "newID" %in% names(comparison)) {
  comparison %>% distinct(cur_sample, newID) %>%
    transmute(cur_sample, patient_group = sub("[BTP]$", "", newID))
} else NULL

#### Parameters ####
bin_sizes <- c(20L, 50L, 110L)
clone_size_thresholds <- c(1L, 2L, 5L, 10L, 15L, 20L)
resolution_levels <- c("10um", "25um", "55um")
resolution_labels <- c("10um" = "10 \u00b5m", "25um" = "25 \u00b5m", "55um" = "55 \u00b5m")
method_engblom <- "repair (rerun)"
method_default <- "repair (initial run)"
method_levels <- c("SIPEAC", method_engblom, method_default)
method_colors <- setNames(c("#BA2029", "#1F3557", "#F39A3D"), method_levels)
method_shapes <- setNames(c(16, 17, 15), method_levels)
method_linetypes <- setNames(c("solid", "dotted", "longdash"), method_levels)
method_axis_labels <- setNames(
  c("SIPEAC", "repair\n(rerun)", "repair\n(initial run)"),
  method_levels
)

#### Accuracy ####
input_roots <- c(sipeac_root, ground_truth_root, repair_default_root, repair_b100_root)
missing <- input_roots[!dir.exists(input_roots)]
if (length(missing)) stop("Input directories not found: ", paste(missing, collapse = ", "))
rm(list = ls(ground_truth_cache, all.names = TRUE), envir = ground_truth_cache)

sipeac_files <- sort(list.files(
  sipeac_root,
  pattern = "_final_e100\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))
if (!length(sipeac_files)) stop("No *_final_e100.tsv files found in: ", sipeac_root)

sipeac_data <- bind_rows(lapply(
  sipeac_files,
  summarise_sipeac_file,
  ground_truth_root = ground_truth_root
)) %>%
  filter(is.finite(SIPEAC))

if (any(duplicated(sipeac_data[c("cur_sample", "bin_size", "threshold")]))) {
  stop("Duplicated sample-bin-threshold rows reconstructed from SIPEAC files")
}

repair_rows <- lapply(seq_len(nrow(sipeac_data)), function(i) {
  row <- sipeac_data[i, ]
  ground_truth <- read_ground_truth(
    ground_truth_root,
    row$bin_size,
    row$cur_sample
  )$data
  engblom <- summarise_repair(
    repair_b100_root,
    row$cur_sample,
    row$bin_size,
    row$threshold,
    ground_truth
  )
  default <- summarise_repair(
    repair_default_root,
    row$cur_sample,
    row$bin_size,
    row$threshold,
    ground_truth
  )
  data.frame(
    repair_b100_n_pred_total = engblom$n_pred_total,
    repair_b100_n_eval = engblom$n_eval,
    repair_b100_n_correct = engblom$n_correct,
    repair_b100 = engblom$accuracy,
    repair_default_n_pred_total = default$n_pred_total,
    repair_default_n_eval = default$n_eval,
    repair_default_n_correct = default$n_correct,
    repair_default = default$accuracy,
    source_repair_b100 = engblom$source,
    source_repair_default = default$source,
    stringsAsFactors = FALSE
  )
})

sample_data <- bind_cols(sipeac_data, bind_rows(repair_rows)) %>%
  mutate(resolution = paste0(bin_size / 2, "um")) %>%
  arrange(bin_size, cur_sample, threshold)

#### Matched samples ####
check_columns(
  sample_data,
  c(
    "cur_sample", "bin_size", "threshold", "resolution",
    "SIPEAC", "repair_b100", "repair_default"
  ),
  "sample_data"
)

sample_data <- sample_data %>%
  mutate(
    bin_size = as.integer(bin_size),
    threshold = as.integer(threshold),
    resolution = normalize_resolution(resolution)
  )

if (!is.null(cohort_samples)) {
  sample_data <- filter(sample_data, cur_sample %in% cohort_samples)
}
if (!nrow(sample_data)) stop("No sample data remain after cohort selection")
patient_level <- !is.null(patient_map)
if (patient_level) {
  check_columns(patient_map, c("cur_sample", "patient_group"), "patient_map")
  if (anyDuplicated(patient_map$cur_sample)) stop("Duplicated samples in patient_map")
  sample_data <- sample_data %>% select(-any_of("patient_group")) %>%
    left_join(patient_map[c("cur_sample", "patient_group")], by = "cur_sample")
  if (anyNA(sample_data$patient_group) || any(!nzchar(sample_data$patient_group))) {
    stop("Missing patient groups for selected samples")
  }
} else {
  sample_data$patient_group <- sample_data$cur_sample
}
test_unit <- if (patient_level) "patient" else "sample"
if (anyNA(sample_data$resolution)) stop("Unexpected resolution label in sample_data")
if (any(duplicated(sample_data[c("cur_sample", "resolution", "threshold")]))) {
  stop("sample_data contains duplicated sample-resolution-threshold rows")
}

accuracy_values <- unlist(sample_data[c("SIPEAC", "repair_b100", "repair_default")])
if (any(is.finite(accuracy_values) & (accuracy_values < 0 | accuracy_values > 1))) {
  stop("Accuracy values must be between 0 and 1")
}

#### Statistics ####
matched_data <- sample_data %>%
  filter(
    threshold == 10,
    is.finite(SIPEAC),
    is.finite(repair_b100),
    is.finite(repair_default)
  )
test_data <- collapse_patient_pairs(matched_data, patient_level)
pairwise_tests <- calculate_pairwise_tests(test_data, paired = FALSE) %>%
  mutate(statistical_unit = test_unit)
paired_sensitivity_tests <- calculate_pairwise_tests(test_data, paired = TRUE) %>%
  mutate(statistical_unit = test_unit)
trend_input <- sample_data %>% filter(is.finite(SIPEAC), is.finite(repair_default))
specimen_trend_tests <- calculate_trend_tests(trend_input) %>% mutate(statistical_unit = "sample")
trend_tests <- calculate_trend_tests(collapse_patient_pairs(trend_input, patient_level)) %>%
  mutate(statistical_unit = test_unit)
trend_data <- trend_input %>%
  pivot_longer(
    c(SIPEAC, repair_b100, repair_default),
    names_to = "method",
    values_to = "accuracy"
  ) %>%
  filter(is.finite(accuracy)) %>%
  mutate(
    method = recode(
      method,
      "SIPEAC" = "SIPEAC",
      "repair_b100" = method_engblom,
      "repair_default" = method_default
    ),
    method = factor(method, levels = method_levels),
    clone_size = factor(threshold, levels = clone_size_thresholds)
  ) %>%
  group_by(resolution, clone_size, method) %>%
  summarise(
    n = n(),
    n_clusters = n_distinct(patient_group),
    mean_accuracy = mean(accuracy),
    se = cluster_mean_se(accuracy, patient_group),
    .groups = "drop"
  )

overlap_labels <- trend_data %>%
  filter(clone_size == "10", method %in% c(method_engblom, method_default)) %>%
  group_by(resolution, clone_size) %>%
  summarise(
    overlap_y = mean(mean_accuracy),
    spread = diff(range(mean_accuracy)),
    n_methods = n(),
    .groups = "drop"
  ) %>%
  filter(n_methods == 2, spread <= 0.002) %>%
  mutate(label_x = 2.75, label_y = 0.945, point_x = 4)

#### Figures ####
panel_a <- ggplot(
  trend_data,
  aes(
    clone_size,
    mean_accuracy,
    color = method,
    shape = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.72, position = position_dodge(width = 0.24), na.rm = TRUE) +
  geom_errorbar(
    aes(ymin = mean_accuracy - se, ymax = mean_accuracy + se),
    width = 0.1,
    linewidth = 0.3,
    position = position_dodge(width = 0.24),
    na.rm = TRUE
  ) +
  geom_point(size = 1.85, position = position_dodge(width = 0.24), na.rm = TRUE) +
  geom_segment(
    data = overlap_labels,
    aes(
      x = label_x + 0.45,
      xend = point_x - 0.08,
      y = label_y - 0.015,
      yend = overlap_y + 0.012
    ),
    inherit.aes = FALSE,
    linewidth = 0.3,
    color = "grey25",
    arrow = arrow(length = unit(0.045, "inches"), type = "closed")
  ) +
  geom_label(
    data = overlap_labels,
    aes(x = label_x, y = label_y, label = "overlapping means\nrerun / initial run"),
    inherit.aes = FALSE,
    size = 2,
    lineheight = 0.88,
    label.size = 0.25,
    label.r = unit(0.08, "lines"),
    label.padding = unit(0.1, "lines"),
    color = "grey15",
    fill = "white"
  ) +
  facet_wrap(~resolution, nrow = 1, labeller = as_labeller(resolution_labels)) +
  scale_color_manual(values = method_colors, breaks = method_levels) +
  scale_shape_manual(values = method_shapes, breaks = method_levels) +
  scale_linetype_manual(values = method_linetypes, breaks = method_levels) +
  scale_y_continuous(breaks = seq(0.2, 0.8, 0.2)) +
  coord_cartesian(ylim = c(0.19, 0.99), clip = "off") +
  labs(
    title = "A  Pairing accuracy",
    x = NULL,
    y = "Pairing accuracy \u00b1 SE"
  ) +
  theme_figure(7.2) +
  theme(
    legend.key.width = unit(0.7, "cm"),
    legend.spacing.x = unit(0.12, "cm")
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  )

p_lookup <- pairwise_tests %>%
  mutate(resolution = normalize_resolution(resolution)) %>%
  filter(comparison == "SIPEAC vs repair_b100") %>%
  transmute(resolution, n_test_units = n_A, p_value)

matched_summary <- matched_data %>%
  group_by(resolution) %>%
  summarise(
    n_matched = n(),
    n_clusters = n_distinct(patient_group),
    SIPEAC_mean = mean(SIPEAC),
    repair_b100_mean = mean(repair_b100),
    repair_default_mean = mean(repair_default),
    .groups = "drop"
  ) %>%
  left_join(p_lookup, by = "resolution")

if (any(is.na(matched_summary$p_value) & matched_summary$n_test_units >= 3L)) {
  stop("Missing SIPEAC vs repair_b100 P value with at least three test units")
}
if (any(matched_summary$n_clusters != matched_summary$n_test_units) ||
    any(pairwise_tests$n_A != pairwise_tests$n_B)) {
  stop("Complete three-way test-unit counts differ from the comparison group counts")
}

pairwise_annotations <- pairwise_tests %>%
  mutate(
    resolution = normalize_resolution(resolution),
    x_start = ifelse(comparison == "repair_b100 vs repair_default", 2, 1),
    x_end = ifelse(comparison == "SIPEAC vs repair_b100", 2, 3),
    y = case_when(
      comparison == "repair_b100 vs repair_default" ~ 1.07,
      comparison == "SIPEAC vs repair_b100" ~ 1.17,
      comparison == "SIPEAC vs repair_default" ~ 1.27
    ),
    label = format_p(p_value)
  )
stopifnot(!anyNA(pairwise_annotations$y))

box_data <- matched_data %>%
  transmute(
    cur_sample,
    resolution,
    SIPEAC,
    repair_engblom = repair_b100,
    repair_default
  ) %>%
  pivot_longer(
    c(SIPEAC, repair_engblom, repair_default),
    names_to = "method",
    values_to = "accuracy"
  ) %>%
  mutate(
    method = recode(
      method,
      "SIPEAC" = "SIPEAC",
      "repair_engblom" = method_engblom,
      "repair_default" = method_default
    ),
    method = factor(method, levels = method_levels)
  )

panel_b <- ggplot(box_data, aes(method, accuracy, color = method, fill = method, shape = method)) +
  geom_boxplot(width = 0.48, linewidth = 0.32, alpha = 0.18, outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0, seed = 2026),
    size = 1.2,
    alpha = 0.74,
    stroke = 0
  ) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 6, color = "black", show.legend = FALSE) +
  geom_segment(
    data = pairwise_annotations,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 0.25, color = "grey25"
  ) +
  geom_segment(
    data = pairwise_annotations,
    aes(x = x_start, xend = x_start, y = y, yend = y - 0.02),
    inherit.aes = FALSE, linewidth = 0.25, color = "grey25"
  ) +
  geom_segment(
    data = pairwise_annotations,
    aes(x = x_end, xend = x_end, y = y, yend = y - 0.02),
    inherit.aes = FALSE, linewidth = 0.25, color = "grey25"
  ) +
  geom_text(
    data = pairwise_annotations,
    aes(x = (x_start + x_end) / 2, y = y + 0.028, label = label),
    inherit.aes = FALSE, size = 2.05, color = "grey15"
  ) +
  facet_wrap(~resolution, nrow = 1, labeller = as_labeller(resolution_labels)) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = alpha(method_colors, 0.18)) +
  scale_shape_manual(values = method_shapes) +
  scale_x_discrete(labels = method_axis_labels) +
  scale_y_continuous(
    limits = c(0, 1.34),
    breaks = seq(0, 1, 0.25),
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    title = "B  Clones \u226510 spots/bins",
    x = NULL,
    y = "Pairing accuracy"
  ) +
  theme_figure(7.5) +
  theme(
    axis.text.x = element_text(size = 5.8, lineheight = 0.88, color = "grey12"),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

figure <- wrap_plots(panel_a, panel_b, ncol = 1, heights = c(1, 1.15))

result <- list(
  panel_a = panel_a,
  panel_b = panel_b,
  figure = figure,
  sample_data = sample_data,
  trend_summary = trend_data,
  trend_tests = trend_tests,
  specimen_trend_tests = specimen_trend_tests,
  matched_data = matched_data,
  matched_summary = matched_summary,
  pairwise_annotations = pairwise_annotations,
  pairwise_tests = pairwise_tests,
  paired_sensitivity_tests = paired_sensitivity_tests
)

#### Output ####
result$input_manifest <- build_input_manifest(sample_data, roots)
save_r4q7_ab(result, figure_dir, table_dir)
print(trend_tests)
print(pairwise_tests)
