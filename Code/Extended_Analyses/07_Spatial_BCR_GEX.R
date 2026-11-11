#### Packages and paths ####
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/R/io.R")
source("Code/Extended_Analyses/R/plot.R")
library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
source("Code/Extended_Analyses/R/spatial_bcr_concordance.R")
source("Code/Extended_Analyses/R/07_Spatial_functions.R")

#### Parameters ####
PRIMARY_SAMPLE <- "2907T"
mode <- "primary"
n_perm <- 1000L
overwrite <- FALSE
BCR_DIR <- file.path(dirs$sipeac, "input_allelic", "bin50")
CACHE_DIR <- file.path(dirs$processed, "R4_Q23_spatial_concordance")
FIGURE_DIR <- file.path(dirs$figures, "R4_Q23_空间BCR信号归一化")
TABLE_DIR <- dirs$tables
stopifnot(mode %in% c("primary", "cohort", "all"), is.finite(n_perm), n_perm >= 1L)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

DAY23_PALETTE <- list(
  light = c(
    blue = "#C8DFF9", pink = "#FEE0DF",
    orange = "#FFE0BF", green = "#ADDEAD"
  ),
  saturated = c(
    blue = "#57A1FC", pink = "#FE8182",
    orange = "#FFA23E", green = "#059962"
  ),
  continuous_blue = c("#E8E8E8", "#C8DFF9", "#8FC3FA", "#57A1FC", "#285E9D"),
  continuous_orange = c("#E8E8E8", "#FFE0BF", "#FFC477", "#FFA23E", "#A94F00"),
  neutral = "#E8E8E8"
)

#### Samples ####
samples <- eligible_samples()
if (!PRIMARY_SAMPLE %in% samples) stop("Primary sample is not eligible: ", PRIMARY_SAMPLE)
sample_groups <- list()
if (mode %in% c("primary", "all")) sample_groups[["2907T"]] <- PRIMARY_SAMPLE
if (mode %in% c("cohort", "all")) sample_groups[["cohort"]] <- samples

#### Concordance ####
for (label in names(sample_groups)) {
  sample_ids <- sample_groups[[label]]
  metrics_list <- vector("list", length(sample_ids))
  qc_list <- vector("list", length(sample_ids))
  caches <- vector("list", length(sample_ids))
  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[i]
    message(sprintf("[%d/%d] Processing %s", i, length(sample_ids), sample_id))
    caches[[i]] <- build_sample_cache(sample_id, overwrite_cache = overwrite)
    qc_list[[i]] <- sample_qc(caches[[i]])
    metrics_list[[i]] <- sample_metrics(caches[[i]], n_permutations = n_perm)
  }
  metrics <- rbindlist(metrics_list, fill = TRUE)
  qc <- rbindlist(qc_list, fill = TRUE)
  metrics[, p_adj_rho := p.adjust(p_value, method = "BH")]
  metrics[, p_adj_top_dice := p.adjust(p_top_dice, method = "BH")]

  write_csv(metrics, paste0("R4_Q23_", label, "_spatial_concordance_metrics.csv"))
  write_csv(qc, paste0("R4_Q23_", label, "_input_QC.csv"))
  log_stats(
    as.data.frame(metrics),
    paste0("R4_Q23_spatial_permutation_", label)
  )
  if (label == "2907T") {
    make_primary_plots(caches[[match(PRIMARY_SAMPLE, sample_ids)]], metrics[sample == PRIMARY_SAMPLE])
  }
  result <- list(metrics = metrics, qc = qc, caches = caches)
  if (label == "2907T") primary_result <- result
  if (label == "cohort") {
    cohort_result <- result
    write_csv(summarise_cohort_metrics(cohort_result$metrics), "R4_Q23_cohort_summary.csv")
    make_cohort_plots(cohort_result$metrics)
  }
}

#### Output ####
writeLines(capture.output(sessionInfo()), file.path(dirs$stats_logs, "R4_Q23_sessionInfo.txt"))
