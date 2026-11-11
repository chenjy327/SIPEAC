clone_category_counts <- function(df, clone_col, unit_col) {
    pairs <- unique(df[, c(clone_col, unit_col), drop = FALSE])
    pairs <- pairs[!is.na(pairs[[clone_col]]) & !is.na(pairs[[unit_col]]), , drop = FALSE]
    clone_size <- as.integer(table(pairs[[clone_col]]))
    category <- ifelse(clone_size == 1, "Singleton", ifelse(clone_size < 10, "Expanded (2-9)", "Expanded (>=10)"))
    n_by_category <- table(factor(category, levels = EXPAND_LEVELS))
    out <- data.frame(category = EXPAND_LEVELS, n_clones = as.integer(n_by_category), stringsAsFactors = FALSE)
    out$fraction <- out$n_clones/sum(out$n_clones)
    out
}

sample_units <- function(df, unit_col, target_depth) {
    units <- unique(df[[unit_col]])
    if (target_depth > length(units)) {
        stop("target_depth exceeds available unique units")
    }
    keep <- sample(units, size = target_depth, replace = FALSE)
    df[df[[unit_col]] %in% keep, , drop = FALSE]
}

paired_stat <- function(sample_df, outcome, n_boot = 5000L, seed = MASTER_SEED + 99L) {
    wide <- sample_df %>% filter(metric == outcome) %>% select(sample, platform, estimate) %>% tidyr::pivot_wider(names_from = platform, 
        values_from = estimate) %>% filter(!is.na(Spatial), !is.na(`Single-cell`))
    delta <- wide$Spatial - wide$`Single-cell`
    set.seed(seed + match(outcome, unique(sample_df$metric)))
    boot_delta <- replicate(n_boot, median(sample(delta, length(delta), replace = TRUE)))
    wt <- suppressWarnings(wilcox.test(wide$Spatial, wide$`Single-cell`, paired = TRUE, exact = FALSE, 
        conf.int = FALSE))
    data.frame(metric = outcome, n_pairs = nrow(wide), spatial_median = median(wide$Spatial), single_cell_median = median(wide$`Single-cell`), 
        median_paired_difference = median(delta), difference_lower_95 = unname(quantile(boot_delta, 0.025)), 
        difference_upper_95 = unname(quantile(boot_delta, 0.975)), wilcoxon_V = unname(wt$statistic), 
        p_value = wt$p.value, stringsAsFactors = FALSE)
}
