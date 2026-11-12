assert_columns <- function(x, required, label = deparse(substitute(x))) {
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0L) {
    stop(sprintf("%s is missing required columns: %s",
                 label, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

as_productive_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  y <- toupper(trimws(as.character(x)))
  !is.na(y) & y %in% c("TRUE", "T", "1", "YES", "Y")
}

aggregate_unique_bcr <- function(x) {
  assert_columns(x, c("sequence_id", "productive", "locus", "binid"), "BCR input")
  d <- data.table::as.data.table(x)[, .(
    sequence_id = as.character(sequence_id),
    productive = as_productive_flag(productive),
    locus = toupper(trimws(as.character(locus))),
    binid = as.character(binid)
  )]
  d <- d[productive & locus %in% c("IGK", "IGL") &
           !is.na(sequence_id) & nzchar(sequence_id) &
           !is.na(binid) & nzchar(binid)]
  d <- unique(d[, .(sequence_id, locus, binid)])

  if (nrow(d) == 0L) {
    return(data.frame(binid = character(), IGK = integer(), IGL = integer()))
  }

  counts <- d[, .(molecules = .N), by = .(binid, locus)]
  wide <- data.table::dcast(counts, binid ~ locus, value.var = "molecules", fill = 0L)
  for (chain in c("IGK", "IGL")) {
    if (!chain %in% colnames(wide)) wide[, (chain) := 0L]
  }
  data.table::setorder(wide, binid)
  as.data.frame(wide[, .(binid, IGK = as.integer(IGK), IGL = as.integer(IGL))])
}

extract_gex_light_chains <- function(count_matrix) {
  features <- rownames(count_matrix)
  bins <- colnames(count_matrix)
  if (is.null(features) || is.null(bins)) {
    stop("GEX count matrix must have feature and bin names", call. = FALSE)
  }
  if (!"IGKC" %in% features) {
    stop("GEX count matrix does not contain IGKC", call. = FALSE)
  }
  iglc_features <- grep("^IGLC[0-9]+$", features, value = TRUE)
  if (length(iglc_features) == 0L) {
    stop("GEX count matrix does not contain numbered IGLC genes", call. = FALSE)
  }

  list(
    binid = bins,
    IGKC = as.numeric(count_matrix["IGKC", , drop = TRUE]),
    IGLC = as.numeric(Matrix::colSums(count_matrix[iglc_features, , drop = FALSE])),
    iglc_features = iglc_features
  )
}

normalise_gex <- function(gene_count, total_gex, scale_factor = 1e4) {
  if (length(gene_count) != length(total_gex)) {
    stop("gene_count and total_gex must have equal lengths", call. = FALSE)
  }
  if (any(gene_count < 0, na.rm = TRUE) || any(total_gex < 0, na.rm = TRUE)) {
    stop("GEX counts and totals must be non-negative", call. = FALSE)
  }
  out <- rep(0, length(gene_count))
  valid <- is.finite(gene_count) & is.finite(total_gex) & total_gex > 0
  out[valid] <- log1p(scale_factor * gene_count[valid] / total_gex[valid])
  out[!is.finite(out)] <- 0
  out
}

normalise_bcr_sample_cpm <- function(chain_count, sample_total, scale_factor = 1e6) {
  if (length(sample_total) != 1L || !is.finite(sample_total) || sample_total <= 0) {
    stop("sample_total must be one finite positive value", call. = FALSE)
  }
  if (any(chain_count < 0, na.rm = TRUE)) {
    stop("BCR counts must be non-negative", call. = FALSE)
  }
  out <- log1p(scale_factor * chain_count / sample_total)
  out[!is.finite(out)] <- 0
  out
}

positive_percentile <- function(x) {
  out <- rep(0, length(x))
  out[is.na(x)] <- NA_real_
  idx <- which(is.finite(x) & x > 0)
  if (length(idx) > 0L) {
    out[idx] <- rank(x[idx], ties.method = "average") / length(idx)
  }
  out
}

top_positive_mask <- function(x, probability = 0.9) {
  if (!is.numeric(probability) || length(probability) != 1L ||
      probability <= 0 || probability >= 1) {
    stop("probability must be strictly between 0 and 1", call. = FALSE)
  }
  idx <- which(is.finite(x) & x > 0)
  out <- rep(FALSE, length(x))
  if (length(idx) == 0L) return(out)
  threshold <- stats::quantile(x[idx], probs = probability, names = FALSE, type = 7)
  out[idx] <- x[idx] >= threshold
  out
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

binary_overlap <- function(a, b) {
  a <- as.logical(a)
  b <- as.logical(b)
  both <- sum(a & b)
  a_n <- sum(a)
  b_n <- sum(b)
  union_n <- sum(a | b)
  data.frame(
    n_a = a_n,
    n_b = b_n,
    n_both = both,
    dice = if ((a_n + b_n) > 0L) 2 * both / (a_n + b_n) else NA_real_,
    jaccard = if (union_n > 0L) both / union_n else NA_real_
  )
}

concordance_metrics <- function(gex, bcr, top_probability = 0.9) {
  if (length(gex) != length(bcr)) stop("gex and bcr must have equal lengths", call. = FALSE)
  ok <- is.finite(gex) & is.finite(bcr)
  gex <- gex[ok]
  bcr <- bcr[ok]
  g_pos <- gex > 0
  b_pos <- bcr > 0
  union_pos <- g_pos | b_pos
  overlap <- binary_overlap(g_pos, b_pos)
  g_top <- top_positive_mask(gex, top_probability)
  b_top <- top_positive_mask(bcr, top_probability)
  top_overlap <- binary_overlap(g_top, b_top)

  data.frame(
    n_bins = length(gex),
    n_gex_positive = sum(g_pos),
    n_bcr_positive = sum(b_pos),
    n_both = sum(g_pos & b_pos),
    n_gex_only = sum(g_pos & !b_pos),
    n_bcr_only = sum(!g_pos & b_pos),
    n_neither = sum(!g_pos & !b_pos),
    rho_all = safe_spearman(gex, bcr),
    rho_positive_union = if (sum(union_pos) >= 3L) safe_spearman(gex[union_pos], bcr[union_pos]) else NA_real_,
    dice_detected = overlap$dice,
    jaccard_detected = overlap$jaccard,
    n_gex_top = top_overlap$n_a,
    n_bcr_top = top_overlap$n_b,
    n_both_top = top_overlap$n_both,
    dice_top = top_overlap$dice,
    jaccard_top = top_overlap$jaccard
  )
}

make_spatial_blocks <- function(row, col, block_size = 10L) {
  if (length(row) != length(col)) stop("row and col must have equal lengths", call. = FALSE)
  if (length(block_size) != 1L || !is.finite(block_size) || block_size <= 0) {
    stop("block_size must be one finite positive value", call. = FALSE)
  }
  if (any(!is.finite(row)) || any(!is.finite(col))) {
    stop("row and col must be finite", call. = FALSE)
  }
  block_row <- floor((row - min(row)) / block_size)
  block_col <- floor((col - min(col)) / block_size)
  paste(block_row, block_col, sep = "_")
}

block_top_dice <- function(gex, bcr, probability = 0.9) {
  binary_overlap(top_positive_mask(gex, probability),
                 top_positive_mask(bcr, probability))$dice
}

block_permutation_test <- function(x, n_perm = 1000L, seed = 20260811L,
                                   block_size = 10L, top_probability = 0.9) {
  assert_columns(x, c("row", "col", "gex", "bcr"), "block permutation input")
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 1L) {
    stop("n_perm must be one positive integer", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)
  d <- data.table::as.data.table(x)[is.finite(row) & is.finite(col) &
                                      is.finite(gex) & is.finite(bcr)]
  if (nrow(d) == 0L) stop("block permutation input has no finite rows", call. = FALSE)
  d[, block := make_spatial_blocks(row, col, block_size)]
  blocks <- d[, .(gex = sum(gex), bcr = sum(bcr)), by = block]
  n_blocks <- nrow(blocks)
  if (n_blocks < 4L) stop("block permutation requires at least four spatial blocks", call. = FALSE)

  observed_rho <- safe_spearman(blocks$gex, blocks$bcr)
  observed_top_dice <- block_top_dice(blocks$gex, blocks$bcr, top_probability)
  set.seed(seed)
  null_rho <- numeric(n_perm)
  null_top_dice <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    permuted_bcr <- sample(blocks$bcr, replace = FALSE)
    null_rho[i] <- safe_spearman(blocks$gex, permuted_bcr)
    null_top_dice[i] <- block_top_dice(blocks$gex, permuted_bcr, top_probability)
  }

  valid_rho <- is.finite(null_rho)
  valid_dice <- is.finite(null_top_dice)
  data.frame(
    block_size = block_size,
    n_blocks = n_blocks,
    n_perm = n_perm,
    observed_rho = observed_rho,
    null_rho_mean = if (any(valid_rho)) mean(null_rho[valid_rho]) else NA_real_,
    null_rho_sd = if (sum(valid_rho) > 1L) stats::sd(null_rho[valid_rho]) else NA_real_,
    p_value = if (is.finite(observed_rho) && any(valid_rho))
      (1 + sum(null_rho[valid_rho] >= observed_rho)) / (1 + sum(valid_rho)) else NA_real_,
    observed_top_dice = observed_top_dice,
    null_top_dice_mean = if (any(valid_dice)) mean(null_top_dice[valid_dice]) else NA_real_,
    null_top_dice_sd = if (sum(valid_dice) > 1L) stats::sd(null_top_dice[valid_dice]) else NA_real_,
    p_top_dice = if (is.finite(observed_top_dice) && any(valid_dice))
      (1 + sum(null_top_dice[valid_dice] >= observed_top_dice)) /
        (1 + sum(valid_dice)) else NA_real_
  )
}
