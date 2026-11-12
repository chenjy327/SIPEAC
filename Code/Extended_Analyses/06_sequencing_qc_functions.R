theme_cns <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = grid::unit(1.2, "pt"),
      plot.title = element_text(face = "bold", size = base_size),
      legend.key.size = grid::unit(4, "pt"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size)
    )
}

theme_panel <- function(base_size = 8.5) {
  theme_cns(base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.box.margin = margin(-5, 0, 0, 0),
      legend.key.width = grid::unit(10, "pt"),
      panel.grid.major.y = element_line(linewidth = 0.20, color = "grey90"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

save_pdf <- function(plot, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = width, height = height, bg = "white", family = "sans")
  } else {
    grDevices::pdf(path, width = width, height = height, bg = "white", useDingbats = FALSE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}

save_figure <- function(plot, stem, width, height, png = TRUE) {
  save_pdf(plot, file.path(cfg$figure_dir, paste0(stem, ".pdf")), width, height)
  if (png) {
    ggsave(
      file.path(cfg$figure_dir, paste0(stem, ".png")), plot,
      width = width, height = height, units = "in", dpi = 400, bg = "white"
    )
  }
  invisible(plot)
}

valid_per_gene <- function(path) {
  file.exists(path) && file.info(path)$size > 0 &&
    startsWith(readLines(path, n = 1, warn = FALSE), "#total_reads=")
}

write_per_gene <- function(path, genes, total_reads) {
  genes <- as.data.frame(genes)
  if (!all(c("gene_id", "reads") %in% names(genes)) || nrow(genes) == 0) {
    stop("No per-gene counts for: ", path)
  }
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  con <- file(tmp, open = "wt")
  writeLines(sprintf("#total_reads=%.0f", total_reads), con)
  close(con)
  write.table(
    genes[, c("gene_id", "reads")], tmp,
    sep = "\t", row.names = FALSE, col.names = FALSE,
    quote = FALSE, append = TRUE
  )
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write: ", path)
  invisible(path)
}

run_bash <- function(command) {
  status <- system2(cfg$bash, c("-o", "pipefail", "-c", shQuote(command)))
  if (!identical(status, 0L)) stop("Command failed with exit status ", status, "\n", command)
  invisible(status)
}

preprocess_stereo <- function(samples) {
  awk_file <- tempfile("stereo-", tmpdir = cfg$work_dir, fileext = ".awk")
  writeLines(c(
    "BEGIN { FS=OFS=\"\\t\" }",
    "/^#/ { next }",
    "$1 == \"geneID\" { next }",
    "{ gc[$1] += $4 }",
    "END { for (g in gc) print g, gc[g] }"
  ), awk_file)
  on.exit(unlink(awk_file), add = TRUE)

  for (sample in samples) {
    input <- file.path(cfg$gem_dir, paste0(sample, cfg$gem_suffix))
    output <- file.path(cfg$stereo_dir, paste0(sample, "_per_gene.tsv"))
    if (!file.exists(input)) stop("Stereo-seq GEM not found: ", input)
    if (!cfg$force && valid_per_gene(output)) {
      message("[Stereo-seq] cached: ", sample)
      next
    }
    tmp <- tempfile(paste0(sample, "-"), tmpdir = cfg$stereo_dir, fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    command <- sprintf(
      "LC_ALL=C %s -dc %s | %s -f %s > %s",
      shQuote(cfg$gzip), shQuote(input), shQuote(cfg$awk), shQuote(awk_file), shQuote(tmp)
    )
    message("[Stereo-seq] scan: ", sample)
    run_bash(command)
    genes <- read.delim(
      tmp,
      header = FALSE, col.names = c("gene_id", "reads"),
      colClasses = c("character", "numeric")
    ) %>%
      filter(!is.na(gene_id), gene_id != "", is.finite(reads), reads > 0) %>%
      arrange(gene_id)
    write_per_gene(output, genes, sum(genes$reads))
    unlink(tmp)
  }
}

read_transcript_map <- function(path) {
  read.delim(
    path,
    header = FALSE, col.names = c("transcript_id", "gene_id"),
    colClasses = "character", quote = "", comment.char = ""
  ) %>%
    transmute(
      transcript_id = sub("\\..*$", "", transcript_id),
      gene_id = trimws(gene_id)
    ) %>%
    filter(transcript_id != "", gene_id != "") %>%
    distinct(transcript_id, .keep_all = TRUE)
}

preprocess_smrt <- function(samples) {
  transcript_map <- read_transcript_map(cfg$transcript_map)
  select_awk <- tempfile("paf-select-", tmpdir = cfg$work_dir, fileext = ".awk")
  best_awk <- tempfile("paf-best-", tmpdir = cfg$work_dir, fileext = ".awk")
  writeLines(c(
    "BEGIN { OFS=\"\\t\" }",
    "{ print $1, $6, $10, $12, $11 }"
  ), select_awk)
  writeLines(c(
    "BEGIN { FS=OFS=\"\\t\" }",
    "!seen[$1]++ { print $1, $2 }"
  ), best_awk)
  on.exit(unlink(c(select_awk, best_awk)), add = TRUE)

  for (sample in samples) {
    input <- file.path(cfg$smrt_fastq_dir, paste0(sample, cfg$smrt_suffix))
    output <- file.path(cfg$smrt_dir, paste0(sample, "_per_gene.tsv"))
    log_file <- file.path(cfg$log_dir, paste0(sample, "_minimap2.log"))
    if (!file.exists(input)) stop("SMRT-seq FASTQ not found: ", input)
    if (!cfg$force && valid_per_gene(output)) {
      message("[SMRT-seq] cached: ", sample)
      next
    }
    best_file <- tempfile(paste0(sample, "-best-"), tmpdir = cfg$smrt_dir, fileext = ".tsv")
    on.exit(unlink(best_file), add = TRUE)
    command <- paste(
      sprintf(
        "%s -x map-ont --secondary=no -t %d %s %s 2> %s",
        shQuote(cfg$minimap2), cfg$threads, shQuote(cfg$reference_fasta),
        shQuote(input), shQuote(log_file)
      ),
      sprintf("%s -f %s", shQuote(cfg$awk), shQuote(select_awk)),
      sprintf("LC_ALL=C %s -k1,1 -k3,3nr -k4,4nr -k5,5nr", shQuote(cfg$sort)),
      sprintf("%s -f %s > %s", shQuote(cfg$awk), shQuote(best_awk), shQuote(best_file)),
      sep = " | "
    )
    message("[SMRT-seq] map: ", sample)
    run_bash(command)
    best <- read.delim(
      best_file,
      header = FALSE, col.names = c("query_id", "transcript_id"),
      colClasses = "character"
    )
    total_reads <- nrow(best)
    best$transcript_id <- sub("\\..*$", "", best$transcript_id)
    gene_id <- transcript_map$gene_id[match(best$transcript_id, transcript_map$transcript_id)]
    genes <- data.frame(gene_id = gene_id) %>%
      filter(!is.na(gene_id), gene_id != "") %>%
      count(gene_id, name = "reads") %>%
      arrange(gene_id)
    write_per_gene(output, genes, total_reads)
    unlink(best_file)
  }
}

find_10x_file <- function(directory, stem) {
  candidates <- file.path(directory, c(paste0(stem, ".gz"), stem))
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) stop("10x file not found: ", file.path(directory, stem))
  found[[1]]
}

preprocess_scrna <- function(map_table) {
  map_table <- map_table %>% filter(!is.na(gex_dir), gex_dir != "")
  for (i in seq_len(nrow(map_table))) {
    sample <- as.character(map_table$sample[[i]])
    gex_dir <- as.character(map_table$gex_dir[[i]])
    output <- file.path(cfg$sc_dir, paste0(sample, "_per_gene.tsv"))
    if (!cfg$force && valid_per_gene(output)) {
      message("[scRNA-seq] cached: ", sample)
      next
    }
    matrix_dir <- file.path(cfg$gex_root, gex_dir, "filtered_feature_bc_matrix")
    matrix_file <- find_10x_file(matrix_dir, "matrix.mtx")
    feature_file <- find_10x_file(matrix_dir, "features.tsv")
    message("[scRNA-seq] extract: ", sample)
    matrix <- Matrix::readMM(matrix_file)
    features <- read.delim(
      feature_file,
      header = FALSE, stringsAsFactors = FALSE,
      quote = "", comment.char = ""
    )
    if (nrow(matrix) != nrow(features) || ncol(features) < 2) {
      stop("10x matrix/features mismatch: ", matrix_dir)
    }
    umi <- Matrix::rowSums(matrix)
    keep <- is.finite(umi) & umi > 0
    genes <- data.frame(gene_id = as.character(features[[2]][keep]), reads = as.numeric(umi[keep]))
    write_per_gene(output, genes, sum(genes$reads))
    rm(matrix)
    gc(verbose = FALSE)
  }
}

read_qc <- function(path) {
  x <- readxl::read_excel(path)
  if (ncol(x) < 7) stop("QC workbook must contain seven columns")
  x <- as.data.frame(x[, seq_len(7)])
  names(x) <- c(
    "rawID", "newID", "Total_reads", "BCR_count",
    "BCR_clone", "TCR_count", "TCR_clone"
  )
  reads_list <- strsplit(as.character(x$Total_reads), ";", fixed = TRUE)
  parsed_reads <- lapply(reads_list, function(z) suppressWarnings(as.numeric(trimws(z))))
  bad <- vapply(parsed_reads, function(z) any(!is.finite(z)), logical(1))
  if (any(bad)) stop("Non-numeric read counts in QC workbook: ", paste(x$rawID[bad], collapse = ", "))
  x$n_sections <- lengths(parsed_reads)
  x$total_reads <- vapply(parsed_reads, sum, numeric(1))
  x$patient <- sub("^.*(Pt[0-9]+).*$", "\\1", x$newID)
  x$patient[!grepl("Pt[0-9]+", x$newID)] <- NA_character_
  x$tissue <- substr(x$newID, nchar(x$newID), nchar(x$newID))
  x$igtr_frac <- (x$BCR_count + x$TCR_count) / x$total_reads
  x$BCR_saturation <- 1 - x$BCR_clone / x$BCR_count
  x$TCR_saturation <- 1 - x$TCR_clone / x$TCR_count
  x$reads_per_BCRclone <- x$BCR_count / x$BCR_clone
  x$reads_per_TCRclone <- x$TCR_count / x$TCR_clone
  table <- x %>%
    transmute(
      patient, tissue, newID, n_sections, total_reads,
      BCR_count, BCR_clone, TCR_count, TCR_clone,
      igtr_frac = round(igtr_frac * 100, 4),
      BCR_saturation = round(BCR_saturation, 3),
      TCR_saturation = round(TCR_saturation, 3),
      reads_per_BCRclone = round(reads_per_BCRclone, 1),
      reads_per_TCRclone = round(reads_per_TCRclone, 1)
    )
  list(raw = x, table = table)
}

estimate_coverage <- function(n, f1, f2) {
  a <- ifelse(
    f2 > 0,
    (n - 1) * f1 / ((n - 1) * f1 + 2 * f2),
    (n - 1) * pmax(f1 - 1, 0) / ((n - 1) * pmax(f1 - 1, 0) + 2)
  )
  1 - f1 / n * a
}

preprocess_bcr <- function() {
  if (!cfg$force && file.exists(cfg$coverage_csv) && file.info(cfg$coverage_csv)$size > 0) {
    message("[BCR] cached coverage table")
    return(invisible(cfg$coverage_csv))
  }
  message("[BCR] calculate sample coverage")
  x <- readRDS(cfg$bcr_rds)
  needed <- c("sample", "clone_id")
  if (!all(needed %in% names(x))) {
    stop("BCR RDS must contain columns: ", paste(needed, collapse = ", "))
  }
  out <- x %>%
    filter(!is.na(sample), sample != "", !is.na(clone_id), clone_id != "") %>%
    count(sample, clone_id, name = "abundance") %>%
    group_by(sample) %>%
    summarise(
      n = sum(abundance),
      S.obs = n(),
      f1 = sum(abundance == 1),
      f2 = sum(abundance == 2),
      .groups = "drop"
    ) %>%
    mutate(SC = round(estimate_coverage(n, f1, f2), 4), .after = S.obs) %>%
    rename(Assemblage = sample) %>%
    arrange(SC)
  tmp <- paste0(cfg$coverage_csv, ".tmp-", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  write.csv(out, tmp, row.names = FALSE)
  if (file.exists(cfg$coverage_csv)) unlink(cfg$coverage_csv)
  if (!file.rename(tmp, cfg$coverage_csv)) stop("Could not write: ", cfg$coverage_csv)
  invisible(cfg$coverage_csv)
}

read_coverage <- function(path) {
  out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("Assemblage", "n", "S.obs", "SC", "f1", "f2")
  if (!all(needed %in% names(out))) {
    stop("Coverage table must contain: ", paste(needed, collapse = ", "))
  }
  out[, needed] %>% arrange(SC)
}

read_per_gene <- function(input_dir, sample) {
  path <- file.path(input_dir, paste0(sample, "_per_gene.tsv"))
  if (!file.exists(path)) stop("Missing per-gene file: ", path)
  meta <- readLines(path, n = 10, warn = FALSE)
  total_line <- meta[startsWith(meta, "#total_reads=")]
  if (length(total_line) != 1) stop("Invalid total_reads header: ", path)
  total_reads <- suppressWarnings(as.numeric(sub("#total_reads=", "", total_line)))
  genes <- read.delim(
    path,
    header = FALSE, sep = "\t", comment.char = "#",
    col.names = c("gene_id", "reads"), colClasses = c("character", "numeric")
  )
  genes <- genes %>%
    filter(!is.na(gene_id), gene_id != "", is.finite(reads), reads > 0)
  if (!is.finite(total_reads) || total_reads <= 0 || nrow(genes) == 0) {
    stop("Invalid per-gene file: ", path)
  }
  list(sample = sample, total_reads = total_reads, genes = genes, path = path)
}

rarefy_genes <- function(genes, fractions, gene_set = NULL) {
  if (!is.null(gene_set)) genes <- genes[genes$gene_id %in% gene_set, , drop = FALSE]
  counts <- genes$reads
  if (length(counts) == 0) {
    return(rep(0, length(fractions)))
  }
  vapply(fractions, function(r) sum(1 - (1 - r)^counts), numeric(1))
}

build_curves <- function(platform, input_dir, samples, housekeeping) {
  message(sprintf("[%s] reading %d samples", platform, length(samples)))
  map_dfr(samples, function(sample) {
    x <- read_per_gene(input_dir, sample)
    all_genes <- rarefy_genes(x$genes, r_seq)
    hk_genes <- rarefy_genes(x$genes, r_seq, housekeeping)
    tissue <- substr(sample, nchar(sample), nchar(sample))
    bind_rows(
      data.frame(
        platform, sample, tissue,
        layer = "All genes",
        depth_frac = r_seq, depth = r_seq * x$total_reads, detected = all_genes
      ),
      data.frame(
        platform, sample, tissue,
        layer = "Housekeeping",
        depth_frac = r_seq, depth = r_seq * x$total_reads, detected = hk_genes
      )
    )
  })
}

fit_mm <- function(depth, detected) {
  d <- data.frame(x = depth, y = detected)
  max_y <- max(d$y, na.rm = TRUE)
  fit <- suppressWarnings(tryCatch(
    nls(
      y ~ Vmax * x / (K + x),
      data = d,
      start = list(Vmax = max_y * 1.05, K = median(d$x)),
      algorithm = "port", lower = c(Vmax = max_y, K = 0),
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    ),
    error = function(e) NULL
  ))
  if (is.null(fit) || any(!is.finite(coef(fit)))) {
    return(list(Vmax = max_y, K = NA_real_, note = "fit_failed; Vmax set to observed maximum"))
  }
  list(Vmax = unname(coef(fit)["Vmax"]), K = unname(coef(fit)["K"]), note = "ok")
}

summarise_saturation <- function(curves) {
  curves %>%
    group_by(platform, sample, tissue, layer) %>%
    group_modify(function(.x, .y) {
      mm <- fit_mm(.x$depth, .x$detected)
      full <- .x$detected[which.max(.x$depth)]
      data.frame(
        Vmax = round(mm$Vmax, 1),
        K = round(mm$K, 0),
        detected_full = round(full, 1),
        saturation_pct = round(min(full / mm$Vmax * 100, 100), 1),
        fit_note = mm$note
      )
    }) %>%
    ungroup()
}

plot_gene_curve <- function(layer_name, title) {
  d <- median_curve %>% filter(layer == layer_name)
  ggplot(d, aes(depth_frac, median, color = platform, fill = platform)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.72, lineend = "round") +
    scale_color_manual(values = platform_colors) +
    scale_fill_manual(values = platform_colors, guide = "none") +
    scale_x_continuous(
      labels = label_percent(accuracy = 1), breaks = c(0, 0.25, 0.50, 0.75, 1),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = label_percent(accuracy = 1), limits = c(0, 1.02),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = "Library fraction", y = "Gene recovery (%)", title = title) +
    theme_panel()
}

relative_name <- function(path) {
  prefix <- paste0(project_root, .Platform$file.sep)
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1) else path
}
