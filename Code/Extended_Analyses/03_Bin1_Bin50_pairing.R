#### Packages and paths ####
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/utils_io.R")
source("Code/Extended_Analyses/utils_plot.R")
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
source("Code/Extended_Analyses/utils_pairing.R")
source("Code/Extended_Analyses/03_bin_resolution_functions.R")

#### Input ####
ALL_SAMPLE_PATH <- file.path(dirs$metadata, "sc_34.txt")
BIN50_INPUT_DIR <- file.path(dirs$sipeac, "input_allelic", "bin50")
GT_DIR <- file.path(dirs$sipeac, "grouth_truth", "bin50")
SCORE_DIR <- dirs$sipeac_cdr
LEGACY_BIN50_DIR <- file.path(dirs$sipeac, "final_version", "bin50")

FULL_TABLE_PREFIX <- file.path(
  dirs$tables, "R1_Q3_bin1_bin50_pairing_34samples"
)
FULL_STATS_DIR <- file.path(
  dirs$stats_logs, "R1_Q3_bin1_bin50_pairing_34samples"
)
FULL_FIG_DIR <- file.path(
  dirs$figures, "R1_Q3_220nm配对准确率", "paired_34samples"
)
FULL_PROCESSED_PATH <- file.path(
  dirs$processed, "R1_Q3_bin1_bin50_matched_ig_molecules.rds"
)

RESOLUTION_LEVELS <- c("Bin1", "Bin50")
RESOLUTION_LABELS <- c(
  Bin1 = "Bin1\n(220-nm DNB)",
  Bin50 = "Bin50\n(~25-µm aggregate)"
)
RESOLUTION_COLORS <- c(Bin1 = "#2F6FB0", Bin50 = "#D97904")

#### Samples ####
SAMPLE_IDS <- character()
OUTPUT_DIR <- ""
all_samples <- as.character(fread(ALL_SAMPLE_PATH)$sample_name)
if (length(all_samples) != 34L || anyDuplicated(all_samples)) {
  stop("sc_34.txt must contain 34 unique sample names.")
}
samples <- if (length(SAMPLE_IDS)) unique(as.character(SAMPLE_IDS)) else all_samples
unknown <- setdiff(samples, all_samples)
if (length(unknown)) stop("Unknown sample(s): ", paste(unknown, collapse = ", "))
is_full_run <- identical(samples, all_samples)
if (is_full_run && !nzchar(OUTPUT_DIR)) {
  output <- list(stats_dir = FULL_STATS_DIR, figure_dir = FULL_FIG_DIR,
    table_prefix = FULL_TABLE_PREFIX, processed_path = FULL_PROCESSED_PATH,
    persistent = TRUE)
} else {
  output_base <- if (nzchar(OUTPUT_DIR)) OUTPUT_DIR else
    file.path(tempdir(), "R1_Q3_bin1_bin50_pairing_test")
  output <- list(stats_dir = file.path(output_base, "stats"),
    figure_dir = file.path(output_base, "figures"),
    table_prefix = file.path(output_base, "tables", "R1_Q3_pairing_test"),
    processed_path = file.path(output_base, "R1_Q3_matched_molecules_test.rds"),
    persistent = FALSE)
}

dir.create(output$stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output$figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output$table_prefix), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output$processed_path), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output$stats_dir, "predictions", "Bin1"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output$stats_dir, "predictions", "Bin50"),
           recursive = TRUE, showWarnings = FALSE)

message(
  "R1_Q3 matched pairing: ", length(samples), " sample(s); output = ",
  output$stats_dir
)
#### Pairing ####
canonical_list <- list()
summary_list <- list()
audit_list <- list()
legacy_list <- list()

for (sample_id in samples) {
  message("[", sample_id, "] canonical molecule mapping")
  prepared <- prepare_canonical_molecules_r1q3(sample_id)
  canonical <- prepared$canonical
  fixed_gt <- load_fixed_ground_truth_r1q3(
    prepared$paths[["ground_truth"]], canonical
  )
  prepared$audit$n_raw_GT_pairs <- fixed_gt$n_raw_pairs
  prepared$audit$n_fixed_GT_pairs <- fixed_gt$n_fixed_pairs
  prepared$audit$n_fixed_GT_H <- fixed_gt$n_fixed_H

  message("[", sample_id, "] reading CDR3 score matrix once")
  score_result <- read_cdr3aa_scores(prepared$paths[["score"]])
  score_mat <- score_result$score_mat
  rm(score_result)
  validate_score_coverage_r1q3(score_mat, canonical, sample_id)

  sample_results <- list()
  for (resolution in RESOLUTION_LEVELS) {
    message("[", sample_id, "] SIPEAC ", resolution)
    resolution_input <- build_resolution_input_r1q3(canonical, resolution)
    result <- run_resolution_r1q3(
      sample_id, resolution, resolution_input, score_mat, fixed_gt
    )
    sample_results[[resolution]] <- result
    prediction_path <- file.path(
      output$stats_dir, "predictions", resolution,
      paste0(sample_id, "_predictions.tsv")
    )
    fwrite(result$predictions, prediction_path, sep = "\t", quote = FALSE,
           na = "NA")
    summary_list[[paste(sample_id, resolution, sep = "_")]] <- result$summary
  }

  legacy_list[[sample_id]] <- compare_legacy_bin50_r1q3(
    sample_id, sample_results$Bin50$predictions
  )
  canonical_list[[sample_id]] <- canonical
  audit_list[[sample_id]] <- prepared$audit
  rm(canonical, score_mat, fixed_gt, prepared, sample_results)
  invisible(gc(verbose = FALSE))
}

#### Statistics ####
summary_df <- bind_rows(summary_list) %>%
  mutate(resolution = factor(resolution, levels = RESOLUTION_LEVELS)) %>%
  arrange(sample, resolution) %>%
  mutate(resolution = as.character(resolution))
audit_df <- bind_rows(audit_list) %>%
  left_join(bind_rows(legacy_list), by = "sample") %>%
  arrange(sample)
canonical_df <- rbindlist(canonical_list, use.names = TRUE)

if (nrow(summary_df) != length(samples) * 2L) {
  stop("Expected exactly two summary rows per sample.", call. = FALSE)
}
paired_gt_check <- summary_df %>%
  select(sample, resolution, n_fixed_GT_pairs, n_fixed_GT_H) %>%
  pivot_wider(
    names_from = resolution,
    values_from = c(n_fixed_GT_pairs, n_fixed_GT_H)
  )
if (any(
  paired_gt_check$n_fixed_GT_pairs_Bin1 !=
    paired_gt_check$n_fixed_GT_pairs_Bin50
) || any(
  paired_gt_check$n_fixed_GT_H_Bin1 != paired_gt_check$n_fixed_GT_H_Bin50
)) {
  stop("Fixed ground-truth denominators differ across resolutions.",
       call. = FALSE)
}

group_df <- summarise_resolution_r1q3(summary_df)
tests_df <- bind_rows(
  paired_metric_test_r1q3(summary_df, "pair_recall"),
  paired_metric_test_r1q3(summary_df, "H_recovery"),
  paired_metric_test_r1q3(summary_df, "precision")
)
patient_tests_df <- bind_rows(
  patient_metric_test_r1q3(summary_df, "pair_recall"),
  patient_metric_test_r1q3(summary_df, "H_recovery"),
  patient_metric_test_r1q3(summary_df, "precision")
)

#### Output ####
summary_path <- paste0(output$table_prefix, "_per_sample.csv")
group_path <- paste0(output$table_prefix, "_group_summary.csv")
audit_path <- paste0(output$table_prefix, "_input_audit.csv")
tests_path <- paste0(output$table_prefix, "_paired_tests.csv")
patient_tests_path <- paste0(
  output$table_prefix, "_patient_sensitivity_tests.csv"
)
stats_path <- file.path(output$stats_dir, "R1_Q3_pairing_34samples_stats.txt")
caption_path <- paste0(output$table_prefix, "_figure_caption.md")
fwrite(summary_df, summary_path, na = "NA")
fwrite(group_df, group_path, na = "NA")
fwrite(audit_df, audit_path, na = "NA")
fwrite(tests_df, tests_path, na = "NA")
fwrite(patient_tests_df, patient_tests_path, na = "NA")
write_stats_log_r1q3(
  stats_path, summary_df, group_df, tests_df, patient_tests_df, audit_df
)
write_caption_r1q3(
  caption_path, summary_df, tests_df, patient_tests_df, group_df
)

if (output$persistent) {
  saveRDS(canonical_df, output$processed_path, compress = "xz")
}

#### Figure ####
figure <- build_pairing_figure_r1q3(
  summary_df, tests_df, group_df, patient_tests_df
)
figure_prefix <- file.path(
  output$figure_dir, "R1_Q3_bin1_bin50_pairing_34samples"
)
ggsave(
  paste0(figure_prefix, ".pdf"), figure,
  device = cairo_pdf, width = 7.8, height = 3.1, units = "in",
  bg = "white"
)
ggsave(
  paste0(figure_prefix, ".png"), figure,
  width = 7.8, height = 3.1, units = "in", dpi = 600,
  bg = "white", limitsize = FALSE
)

message("Completed R1_Q3 matched pairing analysis.")
message("Per-sample table: ", summary_path)
message("Group table: ", group_path)
message("Paired tests: ", tests_path)
message("Patient sensitivity tests: ", patient_tests_path)
message("Figure: ", paste0(figure_prefix, ".pdf"))
print(group_df)
print(tests_df)
