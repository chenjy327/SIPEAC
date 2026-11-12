#### Packages and paths ####
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/utils_io.R")
source("Code/Extended_Analyses/utils_plot.R")
library(dplyr)
library(purrr)
library(ggplot2)
library(scales)
library(patchwork)
library(readxl)
library(Matrix)
source("Code/Extended_Analyses/06_sequencing_qc_functions.R")

options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026)

#### Input ####
stage <- "all"
sample_subset <- character()
threads <- 8L
overwrite <- FALSE
project_root <- PROJECT_ROOT
cfg <- list(
  qc_xlsx = file.path(dirs$raw, "饱和曲线.xlsx"),
  bcr_rds = file.path(dirs$raw, "merged_bcr.rds"),
  gem_dir = file.path(dirs$raw, "stereo_gem"),
  smrt_fastq_dir = file.path(dirs$raw, "smrt_fastq"),
  gex_root = file.path(dirs$raw, "scRNA_gex"),
  reference_fasta = file.path(dirs$raw, "reference", "Homo_sapiens.GRCh38.cdna.all.fa"),
  transcript_map = file.path(dirs$raw, "reference", "transcript_to_gene_symbol.tsv"),
  housekeeping = file.path(dirs$metadata, "Human_Mouse_Common.csv"),
  sample_list = file.path(dirs$metadata, "sample_45.txt"),
  sc_map = file.path(dirs$metadata, "sc_gex_sample_map.tsv"),
  minimap2 = unname(Sys.which("minimap2")),
  threads = threads, force = overwrite, stage = stage,
  sample_subset = sample_subset,
  work_dir = file.path(dirs$processed, "R4_Q8_six_panel_reproducible"),
  gem_suffix = ".gem.gz", smrt_suffix = "_dedup.inST.fastq",
  figure_dir = file.path(dirs$figures, "R4_Q8_序列饱和度QC"),
  table_dir = dirs$tables
)
stopifnot(cfg$stage %in% c("all", "preprocess", "plot"), cfg$threads >= 1L)
if (length(cfg$sample_subset) && cfg$stage != "preprocess") {
  stop("sample_subset is only used for the preprocess stage")
}
cfg$stereo_dir <- file.path(cfg$work_dir, "stereo_per_gene")
cfg$smrt_dir <- file.path(cfg$work_dir, "smrt_per_gene")
cfg$sc_dir <- file.path(cfg$work_dir, "sc_per_gene")
cfg$coverage_csv <- file.path(cfg$work_dir, "bcr_coverage.csv")
cfg$log_dir <- file.path(cfg$work_dir, "logs")

#### Paths and tools ####
dir.create(cfg$figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$stereo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$smrt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$sc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$log_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  qc_xlsx = cfg$qc_xlsx,
  housekeeping = cfg$housekeeping,
  sample_list = cfg$sample_list,
  sc_map = cfg$sc_map
)
if (cfg$stage %in% c("all", "preprocess")) {
  required_files <- c(
    required_files,
    bcr_rds = cfg$bcr_rds,
    reference_fasta = cfg$reference_fasta,
    transcript_map = cfg$transcript_map
  )
}
missing_files <- names(required_files)[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(sprintf("  %s: %s", missing_files, required_files[missing_files]), collapse = "\n")
  )
}
if (cfg$stage %in% c("all", "preprocess")) {
  raw_dirs <- c(gem_dir = cfg$gem_dir, smrt_fastq_dir = cfg$smrt_fastq_dir, gex_root = cfg$gex_root)
  missing_dirs <- names(raw_dirs)[!dir.exists(raw_dirs)]
  if (length(missing_dirs) > 0) {
    stop(
      "Missing raw directories:\n",
      paste(sprintf("  %s: %s", missing_dirs, raw_dirs[missing_dirs]), collapse = "\n")
    )
  }
  if (.Platform$OS.type != "unix") stop("Raw preprocessing requires a Unix-like environment")
  tools <- c(bash = "bash", gzip = "gzip", awk = "awk", sort = "sort")
  tool_paths <- Sys.which(tools)
  if (any(!nzchar(tool_paths))) stop("Missing command-line tools: ", paste(names(tool_paths)[!nzchar(tool_paths)], collapse = ", "))
  cfg$bash <- unname(tool_paths[["bash"]])
  cfg$gzip <- unname(tool_paths[["gzip"]])
  cfg$awk <- unname(tool_paths[["awk"]])
  cfg$sort <- unname(tool_paths[["sort"]])
  if (!nzchar(cfg$minimap2) || !file.exists(cfg$minimap2)) stop("Set cfg$minimap2 to the minimap2 executable")
} else {
  required_processed <- c(
    coverage = cfg$coverage_csv,
    stereo = cfg$stereo_dir,
    smrt = cfg$smrt_dir,
    sc = cfg$sc_dir
  )
  if (!all(file.exists(required_processed))) stop("Processed work directory is incomplete: ", cfg$work_dir)
}

platform_colors <- c(
  "scRNA-seq" = "#1565C0",
  "Spatial SMRT-seq" = "#B03A2E",
  "Stereo-seq" = "#2E7D32"
)
tissue_colors <- c(B = "#4C78A8", T = "#C44E52", P = "#E69F00")
r_seq <- c(seq(0.01, 0.10, 0.01), seq(0.15, 0.90, 0.05), 1)

#### Samples and preprocessing ####
sample_table <- read.csv(cfg$sample_list, stringsAsFactors = FALSE, check.names = FALSE)
sample_45 <- trimws(as.character(sample_table[[1]]))
sample_45 <- sample_45[nzchar(sample_45)]
if (anyDuplicated(sample_45)) stop("Duplicated sample IDs in sample list")

housekeeping_table <- read.csv(
  cfg$housekeeping,
  sep = ";", stringsAsFactors = FALSE, check.names = FALSE
)
if (!"Human" %in% names(housekeeping_table)) stop("Housekeeping file requires a Human column")
housekeeping <- unique(trimws(as.character(housekeeping_table$Human)))
housekeeping <- housekeeping[!is.na(housekeeping) & nzchar(housekeeping)]

sc_map <- read.delim(cfg$sc_map, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample", "gex_dir") %in% names(sc_map))) {
  stop("scRNA-seq map requires sample and gex_dir columns")
}
duplicated_gex <- sc_map %>%
  filter(!is.na(gex_dir), gex_dir != "") %>%
  count(gex_dir, name = "n") %>%
  filter(n > 1)
if (nrow(duplicated_gex) > 0) stop("Duplicated non-empty gex_dir values in scRNA-seq map")

if (cfg$stage %in% c("all", "preprocess")) {
  selected <- if (length(cfg$sample_subset) > 0) cfg$sample_subset else sample_45
  unknown <- setdiff(selected, sample_45)
  if (length(unknown) > 0) stop("Unknown samples: ", paste(unknown, collapse = ", "))
  selected_sc <- sc_map %>%
    filter(sample %in% selected, !is.na(gex_dir), gex_dir != "")
  preprocess_bcr()
  preprocess_stereo(selected)
  preprocess_smrt(selected)
  preprocess_scrna(selected_sc)
  if (cfg$stage == "preprocess") {
    message("[done] preprocessing complete: ", cfg$work_dir)
  }
}

if (cfg$stage %in% c("all", "plot")) {
  required_intermediates <- c(
    file.path(cfg$stereo_dir, paste0(sample_45, "_per_gene.tsv")),
    file.path(cfg$smrt_dir, paste0(sample_45, "_per_gene.tsv"))
  )
  sc_samples <- unique(as.character(sc_map$sample[!is.na(sc_map$gex_dir) & sc_map$gex_dir != ""]))
  required_intermediates <- c(
    required_intermediates,
    file.path(cfg$sc_dir, paste0(sc_samples, "_per_gene.tsv")),
    cfg$coverage_csv
  )
  missing_intermediates <- required_intermediates[!file.exists(required_intermediates)]
  if (length(missing_intermediates) > 0) {
    stop("Missing generated intermediates:\n", paste(missing_intermediates, collapse = "\n"))
  }

  qc_input <- read_qc(cfg$qc_xlsx)
  qc <- qc_input$table
  coverage <- read_coverage(cfg$coverage_csv)
  if (!setequal(qc_input$raw$rawID, sample_45)) stop("QC workbook and sample list do not match")
  if (!setequal(coverage$Assemblage, sample_45)) stop("BCR coverage and sample list do not match")

  curves <- bind_rows(
    build_curves("Stereo-seq", cfg$stereo_dir, sample_45, housekeeping),
    build_curves("Spatial SMRT-seq", cfg$smrt_dir, sample_45, housekeeping),
    build_curves("scRNA-seq", cfg$sc_dir, sc_samples, housekeeping)
  ) %>%
    mutate(
      platform = factor(platform, levels = names(platform_colors)),
      layer = factor(layer, levels = c("All genes", "Housekeeping"))
    )

  saturation <- summarise_saturation(curves)
  curve_plot_data <- curves %>%
    group_by(platform, sample, layer) %>%
    mutate(recovered = detected / max(detected)) %>%
    ungroup()
  median_curve <- curve_plot_data %>%
    group_by(platform, layer, depth_frac) %>%
    summarise(
      median = median(recovered, na.rm = TRUE),
      q25 = quantile(recovered, 0.25, na.rm = TRUE),
      q75 = quantile(recovered, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  qc_plot <- qc %>%
    arrange(total_reads, newID) %>%
    mutate(tissue = factor(tissue, levels = c("B", "T", "P")), rank = row_number())

  p_a <- ggplot(qc_plot, aes(rank, total_reads / 1e6, fill = tissue)) +
    geom_hline(
      yintercept = median(qc$total_reads) / 1e6, linetype = "dashed",
      linewidth = 0.28, color = "grey45"
    ) +
    geom_line(aes(group = 1), color = "grey55", linewidth = 0.30) +
    geom_point(shape = 21, size = 1.85, stroke = 0.20, color = "white") +
    scale_fill_manual(values = tissue_colors, breaks = c("B", "T", "P")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.03))) +
    labs(x = "Libraries ranked by depth", y = "Total reads (million)", title = "Read depth") +
    theme_panel()

  coverage_plot <- coverage %>%
    arrange(SC) %>%
    mutate(
      rank = row_number(),
      coverage_group = cut(
        SC,
        breaks = c(-Inf, 0.60, 0.90, Inf),
        labels = c("<0.60", "0.60-<0.90", ">=0.90"), right = FALSE
      )
    )

  p_b <- ggplot(coverage_plot, aes(rank, SC)) +
    geom_hline(yintercept = c(0.60, 0.90), linetype = "dashed", linewidth = 0.25, color = "grey55") +
    geom_line(color = "grey35", linewidth = 0.35) +
    geom_point(aes(fill = coverage_group), shape = 21, size = 1.55, color = "white", stroke = 0.20) +
    scale_fill_manual(
      values = c("<0.60" = "#D55E00", "0.60-<0.90" = "#E69F00", ">=0.90" = "#0072B2"),
      labels = c("< 0.60", "0.60 <= SC < 0.90", ">= 0.90")
    ) +
    scale_y_continuous(
      labels = label_number(accuracy = 0.01), limits = c(0.20, 1.01),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.03))) +
    labs(x = "Libraries ranked by SC", y = "Sample coverage (SC)", title = "BCR sample coverage") +
    theme_panel()

  p_c <- plot_gene_curve("All genes", "All-gene recovery")
  p_d <- plot_gene_curve("Housekeeping", "Housekeeping recovery") +
    theme(legend.position = "none")

  saturation_plot <- saturation %>%
    mutate(
      platform = factor(platform, levels = names(platform_colors)),
      layer = factor(layer, levels = c("All genes", "Housekeeping"))
    )

  p_e <- ggplot(saturation_plot, aes(platform, saturation_pct, color = platform)) +
    geom_hline(yintercept = 95, linetype = "dashed", linewidth = 0.30, color = "grey50") +
    geom_boxplot(width = 0.46, outlier.shape = NA, linewidth = 0.36, fill = "white") +
    geom_point(
      position = position_jitter(width = 0.11, height = 0, seed = 2026),
      size = 0.82, alpha = 0.68
    ) +
    facet_wrap(~layer, nrow = 1) +
    scale_color_manual(values = platform_colors) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"), limits = c(70, 101),
      breaks = c(70, 80, 90, 95, 100), expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(x = NULL, y = "Observed / fitted Vmax (%)", title = "Recovery relative to Vmax") +
    theme_panel(base_size = 8.0) +
    theme(axis.text.x = element_text(angle = 26, hjust = 1), legend.position = "none")

  spearman <- cor.test(qc$total_reads, qc$BCR_saturation, method = "spearman", exact = FALSE)

  p_f <- ggplot(qc_plot, aes(total_reads / 1e6, BCR_saturation, fill = tissue)) +
    geom_smooth(
      aes(group = 1),
      method = "lm", formula = y ~ x,
      color = "grey35", fill = "grey80", linewidth = 0.48, alpha = 0.28
    ) +
    geom_point(shape = 21, size = 2.05, stroke = 0.25, color = "white", alpha = 0.92) +
    scale_fill_manual(values = tissue_colors, breaks = c("B", "T", "P")) +
    scale_x_continuous(breaks = c(5, 10, 20, 30), expand = expansion(mult = c(0.03, 0.05))) +
    scale_y_continuous(
      labels = label_number(accuracy = 0.1), limits = c(0.1, 1.0),
      breaks = seq(0.2, 1.0, 0.2), expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(x = "Total reads (million)", y = "BCR sequence redundancy", title = "Depth vs BCR redundancy") +
    theme_panel()

  save_figure(p_a, "R4_Q8f_panel_A_reads_depth", 3.1, 3.0)
  save_figure(p_b, "R4_Q8f_panel_B_BCR_iNEXT_coverage", 3.1, 3.0)
  save_figure(p_c, "R4_Q8f_panel_C_all_gene_recovery", 3.1, 3.0)
  save_figure(p_d, "R4_Q8f_panel_D_housekeeping_recovery", 3.1, 3.0)
  save_figure(p_e, "R4_Q8f_panel_E_gene_saturation_summary", 3.1, 3.0)
  save_figure(p_f, "R4_Q8f_panel_F_reads_vs_BCR_saturation", 3.1, 3.0)

  legend_right <- theme(
    legend.position = "right",
    legend.text = element_text(size = 6.3),
    legend.title = element_text(size = 7.0, face = "bold"),
    legend.key.size = grid::unit(7, "pt"),
    legend.spacing.y = grid::unit(0.5, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.background = element_blank(),
    legend.background = element_blank()
  )

  p_a_combined <- p_a +
    guides(fill = guide_legend(title = "Tissue", ncol = 1, override.aes = list(shape = 22, size = 3.2))) +
    legend_right
  p_b_combined <- p_b +
    guides(fill = guide_legend(title = "SC", ncol = 1, override.aes = list(shape = 22, size = 3.2))) +
    legend_right
  p_c_combined <- p_c +
    guides(color = guide_legend(title = "Platform", ncol = 1, override.aes = list(linewidth = 1.2))) +
    legend_right
  p_f_combined <- p_f +
    guides(fill = guide_legend(
      title = "Tissue", ncol = 1,
      override.aes = list(shape = 22, size = 3.2, alpha = 1)
    )) +
    legend_right

  combined <- (
    (p_a_combined | p_b_combined | p_c_combined) /
      (p_d | p_e | p_f_combined)
  ) +
    plot_layout(guides = "keep") +
    plot_annotation(
      tag_levels = "a",
      theme = theme(
        plot.tag = element_text(face = "bold", size = 10.5),
        plot.margin = margin(2, 3, 2, 3)
      )
    )

  save_figure(combined, "R4_Q8_sequencing_CNS_summary", 9.3, 5.85)

  write.csv(qc, file.path(cfg$table_dir, "QC_metrics_supplementary.csv"), row.names = FALSE)
  write.csv(coverage, file.path(cfg$table_dir, "iNEXT_coverage_saturation.csv"), row.names = FALSE)
  write.csv(
    curve_plot_data %>%
      mutate(platform = as.character(platform), layer = as.character(layer)) %>%
      arrange(platform, layer, sample, depth_frac),
    file.path(cfg$table_dir, "R4_Q8_gene_recovery_curves.csv"),
    row.names = FALSE
  )

  platform_file <- c(
    "Stereo-seq" = "stereo_saturation_summary.csv",
    "Spatial SMRT-seq" = "smrt_saturation_summary.csv",
    "scRNA-seq" = "sc_saturation_summary.csv"
  )
  for (platform in names(platform_file)) {
    out <- saturation %>%
      filter(as.character(.data$platform) == .env$platform) %>%
      mutate(
        layer = recode(as.character(layer), "All genes" = "all_genes", "Housekeeping" = "housekeeping")
      ) %>%
      arrange(factor(layer, levels = c("all_genes", "housekeeping")), sample) %>%
      select(layer, sample, tissue, Vmax, K, detected_full, saturation_pct, fit_note)
    write.csv(out, file.path(cfg$table_dir, platform_file[[platform]]), row.names = FALSE)
  }

  stats <- data.frame(
    test = "Spearman",
    x = "total_reads",
    y = "BCR_sequence_redundancy",
    rho = unname(spearman$estimate),
    p_value = spearman$p.value,
    n = nrow(qc)
  )
  write.csv(stats, file.path(cfg$table_dir, "R4_Q8_reads_vs_BCR_redundancy_spearman.csv"), row.names = FALSE)

  small_inputs <- c(cfg$qc_xlsx, cfg$housekeeping, cfg$sample_list, cfg$sc_map)
  if (all(c(dir.exists(cfg$gem_dir), dir.exists(cfg$smrt_fastq_dir), dir.exists(cfg$gex_root)))) {
    sc_raw_dirs <- file.path(
      cfg$gex_root,
      sc_map$gex_dir[!is.na(sc_map$gex_dir) & sc_map$gex_dir != ""],
      "filtered_feature_bc_matrix"
    )
    input_files <- unique(c(
      small_inputs,
      cfg$bcr_rds, cfg$reference_fasta, cfg$transcript_map,
      file.path(cfg$gem_dir, paste0(sample_45, cfg$gem_suffix)),
      file.path(cfg$smrt_fastq_dir, paste0(sample_45, cfg$smrt_suffix)),
      vapply(sc_raw_dirs, find_10x_file, character(1), stem = "matrix.mtx"),
      vapply(sc_raw_dirs, find_10x_file, character(1), stem = "features.tsv")
    ))
  } else {
    input_files <- unique(c(small_inputs, required_intermediates))
  }
  input_info <- file.info(input_files)
  small_enough_to_hash <- is.finite(input_info$size) & input_info$size <= 50 * 1024^2
  input_md5 <- rep(NA_character_, length(input_files))
  input_md5[small_enough_to_hash] <- unname(tools::md5sum(input_files[small_enough_to_hash]))
  manifest <- data.frame(
    file = vapply(input_files, relative_name, character(1)),
    bytes = input_info$size,
    modified = format(input_info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    md5 = input_md5
  )
  write.csv(manifest, file.path(cfg$table_dir, "R4_Q8_reproducibility_input_manifest.csv"), row.names = FALSE)
  minimap2_version <- if (nzchar(cfg$minimap2) && file.exists(cfg$minimap2)) {
    tryCatch(system2(cfg$minimap2, "--version", stdout = TRUE, stderr = TRUE)[1], error = function(e) NA_character_)
  } else {
    NA_character_
  }
  software <- c(
    paste("stage", cfg$stage, sep = "\t"),
    paste("R", R.version.string, sep = "\t"),
    paste("minimap2", cfg$minimap2, minimap2_version, sep = "\t"),
    paste("gzip", unname(Sys.which("gzip")), sep = "\t"),
    paste("awk", unname(Sys.which("awk")), sep = "\t"),
    paste("sort", unname(Sys.which("sort")), sep = "\t"),
    "",
    capture.output(sessionInfo())
  )
  writeLines(software, file.path(cfg$table_dir, "R4_Q8_sessionInfo.txt"))

  message(sprintf(
    "[done] samples: QC=%d, BCR=%d, Stereo=%d, SMRT=%d, scRNA=%d",
    nrow(qc), nrow(coverage), length(sample_45), length(sample_45), length(sc_samples)
  ))
  message("[done] figure: ", file.path(cfg$figure_dir, "R4_Q8_sequencing_CNS_summary.pdf"))
}
