sample_object_path <- function(sample_id) {
  file.path(
    dirs$st_root,
    paste0(sample_id, "_scRefST2934T_random100_marker20_spotlight_Spatial.rds")
  )
}

sample_meta_path <- function(sample_id) {
  file.path(dirs$st_meta, paste0(sample_id, ".txt"))
}

sample_bcr_path <- function(sample_id) {
  file.path(BCR_DIR, paste0(sample_id, "_info.txt"))
}

sample_cache_path <- function(sample_id) {
  file.path(CACHE_DIR, paste0(sample_id, "_matched_bins.rds"))
}

eligible_samples <- function() {
  bcr_samples <- sub("_info\\.txt$", "", list.files(BCR_DIR, pattern = "_info\\.txt$"))
  object_samples <- sub(
    "_scRefST2934T_random100_marker20_spotlight_Spatial\\.rds$", "",
    list.files(
      dirs$st_root,
      pattern = "_scRefST2934T_random100_marker20_spotlight_Spatial\\.rds$"
    )
  )
  meta_samples <- sub("\\.txt$", "", list.files(dirs$st_meta, pattern = "\\.txt$"))
  sort(Reduce(intersect, list(bcr_samples, object_samples, meta_samples)))
}

read_spatial_components <- function(sample_id) {
  object_path <- sample_object_path(sample_id)
  meta_path <- sample_meta_path(sample_id)
  if (!file.exists(object_path)) stop("Missing spatial object: ", object_path)
  if (!file.exists(meta_path)) stop("Missing spatial metadata: ", meta_path)

  message("[", sample_id, "] Reading spatial object")
  object <- readRDS(object_path)
  if (!"Spatial" %in% names(object@assays)) {
    stop("Spatial assay missing for sample ", sample_id, call. = FALSE)
  }
  counts <- object@assays$Spatial@counts
  if (is.null(colnames(counts)) || anyDuplicated(colnames(counts))) {
    stop("Spatial bins must have unique column names for sample ", sample_id, call. = FALSE)
  }

  meta <- fread(meta_path, sep = "\t", header = TRUE, showProgress = FALSE)
  assert_columns(meta, c("BinID", "row", "col"), paste(sample_id, "metadata"))
  meta[, binid := sub("^[^_]+_", "", as.character(BinID))]
  if (anyDuplicated(meta$binid)) {
    stop("Duplicated bin IDs in metadata for sample ", sample_id, call. = FALSE)
  }
  idx <- match(colnames(counts), meta$binid)
  if (anyNA(idx)) {
    stop(sum(is.na(idx)), " spatial bins lack metadata for sample ", sample_id, call. = FALSE)
  }
  meta <- meta[idx]
  if (!identical(as.character(meta$binid), colnames(counts))) {
    stop("Metadata could not be aligned to GEX bins for sample ", sample_id, call. = FALSE)
  }

  list(counts = counts, meta = meta, object_path = object_path, meta_path = meta_path)
}

read_unique_productive_bcr <- function(sample_id) {
  bcr_path <- sample_bcr_path(sample_id)
  if (!file.exists(bcr_path)) stop("Missing BCR table: ", bcr_path)
  message("[", sample_id, "] Reading productive BCR molecules")
  bcr <- fread(
    bcr_path,
    select = c("sequence_id", "productive", "locus", "binid"),
    showProgress = interactive()
  )
  assert_columns(
    bcr, c("sequence_id", "productive", "locus", "binid"),
    paste(sample_id, "BCR input")
  )
  bcr[, `:=`(
    sequence_id = as.character(sequence_id),
    locus = toupper(trimws(as.character(locus))),
    binid = as.character(binid),
    is_productive = as_productive_flag(productive)
  )]
  bcr <- bcr[
    is_productive & locus %in% c("IGH", "IGK", "IGL") &
      !is.na(sequence_id) & nzchar(sequence_id) &
      !is.na(binid) & nzchar(binid),
    .(sequence_id, productive = TRUE, locus, binid)
  ]
  unique(bcr)
}

build_sample_cache <- function(sample_id, overwrite_cache = FALSE) {
  cache_path <- sample_cache_path(sample_id)
  if (file.exists(cache_path) && !overwrite_cache) {
    message("[", sample_id, "] Using cache: ", cache_path)
    return(readRDS(cache_path))
  }

  spatial <- read_spatial_components(sample_id)
  counts <- spatial$counts
  meta <- spatial$meta
  gex <- extract_gex_light_chains(counts)
  total_gex <- as.numeric(Matrix::colSums(counts))

  bcr_unique <- read_unique_productive_bcr(sample_id)
  sample_total_bcr <- nrow(bcr_unique)
  if (sample_total_bcr == 0L) {
    stop("No unique productive BCR molecules for sample ", sample_id, call. = FALSE)
  }
  light <- aggregate_unique_bcr(bcr_unique)
  light <- as.data.table(light)
  light_idx <- match(gex$binid, light$binid)
  bcr_igk <- fifelse(is.na(light_idx), 0L, light$IGK[light_idx])
  bcr_igl <- fifelse(is.na(light_idx), 0L, light$IGL[light_idx])

  bins <- data.table(
    sample = sample_id,
    binid = gex$binid,
    row = as.numeric(meta$row),
    col = as.numeric(meta$col),
    total_gex = total_gex,
    gex_IGKC_raw = as.numeric(gex$IGKC),
    gex_IGLC_raw = as.numeric(gex$IGLC),
    bcr_IGK_raw = as.integer(bcr_igk),
    bcr_IGL_raw = as.integer(bcr_igl)
  )
  bins[, `:=`(
    gex_IGKC_logCP10K = normalise_gex(gex_IGKC_raw, total_gex),
    gex_IGLC_logCP10K = normalise_gex(gex_IGLC_raw, total_gex),
    bcr_IGK_logCPM = normalise_bcr_sample_cpm(bcr_IGK_raw, sample_total_bcr),
    bcr_IGL_logCPM = normalise_bcr_sample_cpm(bcr_IGL_raw, sample_total_bcr)
  )]
  bins[, `:=`(
    gex_IGKC_percentile = positive_percentile(gex_IGKC_logCP10K),
    gex_IGLC_percentile = positive_percentile(gex_IGLC_logCP10K),
    bcr_IGK_percentile = positive_percentile(bcr_IGK_logCPM),
    bcr_IGL_percentile = positive_percentile(bcr_IGL_logCPM)
  )]

  matched_light_molecules <- sum(bins$bcr_IGK_raw + bins$bcr_IGL_raw)
  cache <- list(
    sample = sample_id,
    bins = as.data.frame(bins),
    iglc_features = gex$iglc_features,
    sample_total_productive_bcr = sample_total_bcr,
    total_unique_light_molecules = sum(bcr_unique$locus %in% c("IGK", "IGL")),
    matched_light_molecules = matched_light_molecules,
    unmatched_light_molecules =
      sum(bcr_unique$locus %in% c("IGK", "IGL")) - matched_light_molecules,
    input_paths = c(
      spatial_object = spatial$object_path,
      spatial_metadata = spatial$meta_path,
      bcr = sample_bcr_path(sample_id)
    )
  )
  saveRDS(cache, cache_path, compress = "xz")
  message("[", sample_id, "] Wrote cache: ", cache_path)
  cache
}

sample_qc <- function(cache) {
  bins <- cache$bins
  data.frame(
    sample = cache$sample,
    n_spatial_bins = nrow(bins),
    n_zero_total_gex_bins = sum(bins$total_gex == 0),
    sample_total_productive_bcr = cache$sample_total_productive_bcr,
    total_unique_light_molecules = cache$total_unique_light_molecules,
    matched_light_molecules = cache$matched_light_molecules,
    unmatched_light_molecules = cache$unmatched_light_molecules,
    fraction_light_molecules_matched =
      cache$matched_light_molecules / cache$total_unique_light_molecules,
    iglc_features = paste(cache$iglc_features, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

sample_metrics <- function(cache, n_permutations = 1000L) {
  bins <- as.data.table(cache$bins)
  chain_spec <- data.frame(
    chain = c("IGKC/IGK", "IGLC/IGL"),
    gex = c("gex_IGKC_logCP10K", "gex_IGLC_logCP10K"),
    bcr = c("bcr_IGK_logCPM", "bcr_IGL_logCPM"),
    stringsAsFactors = FALSE
  )

  results <- lapply(seq_len(nrow(chain_spec)), function(i) {
    gex <- bins[[chain_spec$gex[i]]]
    bcr <- bins[[chain_spec$bcr[i]]]
    descriptive <- concordance_metrics(gex, bcr, top_probability = 0.9)
    permutation <- block_permutation_test(
      data.frame(row = bins$row, col = bins$col, gex = gex, bcr = bcr),
      n_perm = n_permutations,
      seed = 20260811L + match(cache$sample, eligible_samples()) * 10L + i,
      block_size = 10L,
      top_probability = 0.9
    )
    cbind(
      data.frame(
        sample = cache$sample,
        chain = chain_spec$chain[i],
        iglc_features = if (chain_spec$chain[i] == "IGLC/IGL") {
          paste(cache$iglc_features, collapse = ";")
        } else {
          "IGKC"
        },
        stringsAsFactors = FALSE
      ),
      descriptive,
      permutation
    )
  })
  rbindlist(results, fill = TRUE)
}

write_csv <- function(x, filename) {
  path <- file.path(TABLE_DIR, filename)
  fwrite(x, path)
  message("[table] ", path)
  invisible(path)
}

spatial_theme <- function() {
  theme_void(base_size = 8) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "bottom",
      legend.key.width = grid::unit(15, "pt"),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7),
      plot.caption = element_text(size = 6, hjust = 0)
    )
}

make_primary_plots <- function(cache, metrics) {
  bins <- as.data.table(cache$bins)
  bins[, `:=`(x = max(col) - col, y = max(row) - row)]

  make_chain_map_row <- function(chain_label, gex_col, bcr_col, colours) {
    chain_maps <- rbindlist(list(
      bins[, .(x, y, modality = "Stereo-seq GEX", value = get(gex_col))],
      bins[, .(x, y, modality = "Spatial BCR", value = get(bcr_col))],
      bins[, .(
        x, y, modality = "Shared GEX-BCR",
        value = pmin(get(gex_col), get(bcr_col))
      )]
    ))
    chain_maps[, modality := factor(
      modality,
      levels = c("Stereo-seq GEX", "Spatial BCR", "Shared GEX-BCR")
    )]
    ggplot(chain_maps, aes(x, y, colour = value)) +
      geom_point(shape = 15, size = 0.18, stroke = 0) +
      facet_wrap(~modality, nrow = 1) +
      coord_fixed() +
      scale_colour_gradientn(
        colours = colours,
        limits = c(0, 1), oob = squish,
        name = "Within-section\npositive percentile"
      ) +
      labs(title = chain_label) +
      spatial_theme()
  }

  p_map_igkc <- make_chain_map_row(
    "IGKC/IGK", "gex_IGKC_percentile", "bcr_IGK_percentile",
    DAY23_PALETTE$continuous_blue
  )
  p_map_iglc <- make_chain_map_row(
    "IGLC/IGL", "gex_IGLC_percentile", "bcr_IGL_percentile",
    DAY23_PALETTE$continuous_orange
  )
  p_map <- p_map_igkc / p_map_iglc +
    plot_annotation(
      title = paste0(cache$sample, ": matched-scale and shared light-chain spatial maps"),
      subtitle = paste(
        "Shared GEX-BCR is the per-bin minimum of the two positive-percentile signals;",
        "it is zero unless both modalities are detected."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 10, margin = margin(l = 18)),
        plot.subtitle = element_text(size = 7, margin = margin(l = 18))
      )
    )

  overlap_data <- rbindlist(list(
    bins[, .(
      x, y, chain = "IGKC/IGK",
      overlap = fifelse(
        gex_IGKC_raw > 0 & bcr_IGK_raw > 0, "Detected in both",
        fifelse(
          gex_IGKC_raw > 0, "GEX detected only",
          fifelse(bcr_IGK_raw > 0, "BCR detected only", "Neither detected")
        )
      )
    )],
    bins[, .(
      x, y, chain = "IGLC/IGL",
      overlap = fifelse(
        gex_IGLC_raw > 0 & bcr_IGL_raw > 0, "Detected in both",
        fifelse(
          gex_IGLC_raw > 0, "GEX detected only",
          fifelse(bcr_IGL_raw > 0, "BCR detected only", "Neither detected")
        )
      )
    )]
  ))
  overlap_data[, overlap := factor(
    overlap,
    levels = c(
      "Detected in both", "GEX detected only",
      "BCR detected only", "Neither detected"
    )
  )]
  overlap_colours <- c(
    "Detected in both" = unname(DAY23_PALETTE$saturated["green"]),
    "GEX detected only" = unname(DAY23_PALETTE$saturated["blue"]),
    "BCR detected only" = unname(DAY23_PALETTE$saturated["pink"]),
    "Neither detected" = DAY23_PALETTE$neutral
  )
  overlap_summary <- overlap_data[, .(n_bins = .N), by = .(chain, overlap)]
  overlap_summary[, percentage := 100 * n_bins / sum(n_bins), by = chain]
  overlap_summary[, percentage_label := fifelse(
    percentage > 0 & percentage < 0.1,
    "<0.1%", sprintf("%.1f%%", percentage)
  )]
  write_csv(
    overlap_summary,
    paste0("R4_Q23_", cache$sample, "_detection_category_percentages.csv")
  )

  p_overlap_maps <- ggplot(overlap_data, aes(x, y, colour = overlap)) +
    geom_point(shape = 15, size = 0.18, stroke = 0) +
    facet_wrap(~chain, nrow = 1) +
    coord_fixed() +
    scale_colour_manual(values = overlap_colours, drop = FALSE) +
    labs(
      title = "Spatial distribution",
      colour = "Bin category"
    ) +
    guides(colour = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(shape = 15, size = 3)
    )) +
    spatial_theme()

  p_overlap_bars <- ggplot(
    overlap_summary,
    aes(x = overlap, y = percentage, fill = overlap)
  ) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
    geom_text(
      aes(label = percentage_label),
      vjust = -0.35, size = 2.5, colour = "black"
    ) +
    facet_wrap(~chain, nrow = 1) +
    scale_fill_manual(values = overlap_colours, drop = FALSE) +
    scale_x_discrete(labels = c(
      "Detected in both" = "GEX + BCR",
      "GEX detected only" = "GEX only",
      "BCR detected only" = "BCR only",
      "Neither detected" = "Neither"
    )) +
    scale_y_continuous(
      limits = c(0, 80), breaks = seq(0, 80, by = 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = "Percentage of all tissue bins",
      x = NULL, y = "Tissue bins",
      caption = paste0("Denominator for each chain: ", format(nrow(bins), big.mark = ","), " tissue bins")
    ) +
    theme_csl(base_size = 8) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
      plot.caption = element_text(size = 6, hjust = 0)
    )

  p_overlap <- (p_overlap_maps | p_overlap_bars) +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = "Which modality detects the light-chain signal in each bin?",
      subtitle = "Each tissue bin is assigned to one mutually exclusive detection category.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 10, margin = margin(l = 18)),
        plot.subtitle = element_text(size = 7, margin = margin(l = 18))
      )
    )

  block_data <- rbindlist(lapply(seq_len(nrow(metrics)), function(i) {
    chain <- metrics$chain[i]
    gex_col <- if (chain == "IGKC/IGK") "gex_IGKC_logCP10K" else "gex_IGLC_logCP10K"
    bcr_col <- if (chain == "IGKC/IGK") "bcr_IGK_logCPM" else "bcr_IGL_logCPM"
    d <- bins[, .(row, col, gex = get(gex_col), bcr = get(bcr_col))]
    d[, block := make_spatial_blocks(row, col, block_size = 10L)]
    d[, .(gex = sum(gex), bcr = sum(bcr)), by = block][, chain := chain]
  }))
  labels <- metrics[, .(
    chain,
    label = sprintf(
      "block Spearman rho = %.3f\npermutation P = %.4g\n%d blocks",
      observed_rho, p_value, n_blocks
    )
  )]
  p_blocks <- ggplot(block_data, aes(gex, bcr, colour = chain)) +
    geom_point(size = 0.6, alpha = 0.55, show.legend = FALSE) +
    facet_wrap(~chain, scales = "free", nrow = 1) +
    scale_colour_manual(values = c(
      "IGKC/IGK" = unname(DAY23_PALETTE$saturated["blue"]),
      "IGLC/IGL" = unname(DAY23_PALETTE$saturated["orange"])
    )) +
    geom_text(
      data = labels,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.05, vjust = 1.1, size = 2.2, inherit.aes = FALSE
    ) +
    labs(
      title = "Spatial-block concordance",
      x = "Block-summed GEX log(CP10K + 1)",
      y = "Block-summed BCR log(CPM + 1)"
    ) +
    theme_csl(base_size = 8) +
    theme(strip.text = element_text(face = "bold"))

  prefix <- file.path(FIGURE_DIR, paste0("R4_Q23_", cache$sample, "_normalised_spatial_maps"))
  save_plot(p_map, prefix, width = 7.8, height = 5.8, if_png = TRUE, if_pdf = TRUE)
  save_plot(
    p_overlap,
    file.path(FIGURE_DIR, paste0("R4_Q23_", cache$sample, "_detection_overlap")),
    width = 9.4, height = 3.9, if_png = TRUE, if_pdf = TRUE
  )
  save_plot(
    p_blocks,
    file.path(FIGURE_DIR, paste0("R4_Q23_", cache$sample, "_block_concordance")),
    width = 7.8, height = 3.1, if_png = TRUE, if_pdf = TRUE
  )
  combined <- wrap_elements(full = p_map) / wrap_elements(full = p_overlap) / p_blocks +
    plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(2.0, 1.15, 1))
  save_plot(
    combined,
    file.path(FIGURE_DIR, paste0("R4_Q23_", cache$sample, "_response_figure")),
    width = 9.4, height = 11.8, if_png = TRUE, if_pdf = TRUE
  )
}

make_cohort_plots <- function(metrics) {
  metrics <- as.data.table(metrics)
  metrics[, chain := factor(chain, levels = c("IGKC/IGK", "IGLC/IGL"))]
  chain_colours <- c(
    "IGKC/IGK" = unname(DAY23_PALETTE$saturated["blue"]),
    "IGLC/IGL" = unname(DAY23_PALETTE$saturated["orange"])
  )
  chain_fills <- c(
    "IGKC/IGK" = unname(DAY23_PALETTE$light["blue"]),
    "IGLC/IGL" = unname(DAY23_PALETTE$light["orange"])
  )
  p_rho <- ggplot(metrics, aes(chain, observed_rho, colour = chain, fill = chain)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_boxplot(width = 0.45, outlier.shape = NA, linewidth = 0.35, alpha = 0.45) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.1, alpha = 0.75) +
    scale_colour_manual(values = chain_colours) +
    scale_fill_manual(values = chain_fills) +
    labs(
      title = "Cohort spatial-block concordance",
      subtitle = paste0(length(unique(metrics$sample)), " matched tissue sections; each point is one section"),
      x = NULL, y = "Spatial-block Spearman rho"
    ) +
    theme_csl(base_size = 8) + theme(legend.position = "none")

  p_dice <- ggplot(metrics, aes(chain, dice_detected, colour = chain, fill = chain)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, linewidth = 0.35, alpha = 0.45) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.1, alpha = 0.75) +
    scale_colour_manual(values = chain_colours) +
    scale_fill_manual(values = chain_fills) +
    labs(x = NULL, y = "Dice overlap of detected bins") +
    theme_csl(base_size = 8) + theme(legend.position = "none")

  combined <- p_rho | p_dice
  save_plot(
    combined,
    file.path(FIGURE_DIR, "R4_Q23_cohort_concordance"),
    width = 7.0, height = 3.2, if_png = TRUE, if_pdf = TRUE
  )
}

summarise_cohort_metrics <- function(metrics) {
  metrics <- as.data.table(metrics)
  metrics[, .(
    n_sections = uniqueN(sample),
    median_bin_rho = median(rho_all),
    min_bin_rho = min(rho_all),
    max_bin_rho = max(rho_all),
    median_block_rho = median(observed_rho),
    min_block_rho = min(observed_rho),
    max_block_rho = max(observed_rho),
    n_positive_block_rho = sum(observed_rho > 0),
    n_BH_significant_block_rho = sum(p_adj_rho < 0.05),
    median_detected_bin_dice = median(dice_detected),
    median_top_block_dice = median(observed_top_dice)
  ), by = chain]
}
