library(glue)

proj_path <- function(...) file.path(PROJECT_ROOT, ...)

dirs <- list(
  raw = file.path(data_dir, "raw"),
  processed = file.path(PROJECT_ROOT, "data", "processed", "Extended_Analyses"),
  metadata = file.path(data_dir, "metadata"),
  figures = file.path(PROJECT_ROOT, "results", "figures", "Extended_Analyses"),
  tables = file.path(PROJECT_ROOT, "results", "tables", "Extended_Analyses"),
  stats_logs = file.path(PROJECT_ROOT, "results", "stats_logs", "Extended_Analyses"),
  spatial = file.path(data_dir, "raw", "spatial_stereoseq"),
  st_root = file.path(data_dir, "raw", "spatial_stereoseq", "samples_60TLS"),
  st_meta = file.path(data_dir, "raw", "spatial_stereoseq", "Bin50经典版注释metadata"),
  bin1 = file.path(data_dir, "raw", "spatial_stereoseq", "bin1"),
  sipeac = file.path(data_dir, "raw", "spatial_stereoseq", "SIPEAC"),
  sipeac_cdr = file.path(data_dir, "raw", "spatial_stereoseq", "SIPEAC", "2025_fold0_epoch100"),
  code = file.path(PROJECT_ROOT, "Code"),
  sipeac_model = file.path(PROJECT_ROOT, "Code", "antibody_HL_SIPEAC")
)

read_rds <- function(name, subdir = "processed") {
  base <- if (subdir == "processed") dirs$processed else file.path(data_dir, subdir)
  path <- file.path(base, name)
  if (!file.exists(path)) stop("文件不存在: ", path)
  message(sprintf("[read_rds] %s", path))
  readRDS(path)
}

write_rds <- function(obj, name, subdir = "processed") {
  base <- if (subdir == "processed") dirs$processed else file.path(data_dir, subdir)
  path <- file.path(base, name)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path, compress = "xz")
  message(sprintf(
    "[write_rds] %s (%0.1f MB)",
    path, file.info(path)$size / 1e6
  ))
  invisible(path)
}

read_meta <- function(name) {
  path <- file.path(data_dir, "metadata", name)
  if (!file.exists(path)) stop("metadata 不存在: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

log_stats <- function(df, label) {
  stopifnot(is.data.frame(df))
  dir.create(dirs$stats_logs, showWarnings = FALSE, recursive = TRUE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  path <- file.path(dirs$stats_logs, sprintf("%s_%s.csv", label, stamp))
  write.csv(df, path, row.names = FALSE)
  message(sprintf("[log_stats] %s", path))
  invisible(df)
}

save_table <- function(df, prefix, row.names = FALSE, if_csv = FALSE, col.names = TRUE) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  if (if_csv) {
    write.csv(df,
      glue("{prefix}.csv"),
      row.names = row.names,
      quote = FALSE
    )
  } else {
    write.table(df,
      glue("{prefix}.txt"),
      col.names = col.names,
      row.names = row.names,
      quote = FALSE,
      sep = "\t",
      fileEncoding = "UTF-8"
    )
  }
}

calculate_countSum <- function(data_obj, l) {
  count_mat <- data_obj@assays$Spatial@counts
  for (item in names(l)) {
    used_genes <- intersect(row.names(count_mat), l[[item]])
    pro_count_mat <- count_mat[used_genes, ]
    sum_gene_counts <- colSums(pro_count_mat)
    count_str <- glue("{item}_CountSum")
    data_obj@meta.data[[count_str]] <- sum_gene_counts
  }
  return(data_obj)
}

n_fun <- function(x) {
  return(data.frame(y = min(x), label = length(x)))
}

mean_fun <- function(x) {
  return(data.frame(y = max(x), label = round(mean(x), 3)))
}
