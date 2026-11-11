PROJECT_ROOT <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, 'data')
source('Code/Extended_Analyses/R/io.R')
source('Code/Extended_Analyses/R/plot.R')
source('Code/Extended_Analyses/R/04_Clonotype_functions.R')

#### 参数 ====

library(dplyr)

library(ggplot2)

library(patchwork)

DATA_BCR_SMRT <- file.path(dirs$raw, "scRNAseq", "merged_bcr.rds")

DATA_BCR_SC <- file.path(dirs$raw, "scRNAseq", "merge4.rds")

DATA_TCR_SMRT <- file.path(dirs$raw, "scRNAseq", "merged_tcr.rds")

FIG_DIR <- file.path(dirs$figures, "R2_Q3_克隆型归一化")

TBL_DIR <- dirs$tables

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

N_ITER <- 500L

MIN_MATCHED_DEPTH <- 50L

MASTER_SEED <- 20260725L

EXPAND_LEVELS <- c("Singleton", "Expanded (2-9)", "Expanded (>=10)")

EXPAND_COLS <- c(Singleton = "#65A3D1", `Expanded (2-9)` = "#EDB55A", `Expanded (>=10)` = "#946625")

PLATFORM_COLS <- c(Spatial = "#C0392B", `Single-cell` = "#4C78A8")

#### BCR样本内匹配深度与重复抽样 ====

bcr_spatial <- readRDS(DATA_BCR_SMRT) %>% filter(locus == "IGH", !is.na(sample), !is.na(clone_id), !is.na(bin_id)) %>% 
    mutate(clone_key = paste(sample, clone_id, sep = "__"))

bcr_sc <- readRDS(DATA_BCR_SC) %>% filter(type == "SingleCell", locus == "IGH", !is.na(sample), !is.na(clone_id), 
    !is.na(sequence_id)) %>% distinct(sequence_id, .keep_all = TRUE) %>% mutate(clone_key = paste(sample, 
    clone_id, sep = "__"))

depth_spatial <- bcr_spatial %>% summarise(spatial_units = n_distinct(bin_id), .by = sample)

depth_sc <- bcr_sc %>% summarise(single_cell_units = n_distinct(sequence_id), .by = sample)

matched_depths <- inner_join(depth_spatial, depth_sc, by = "sample") %>% mutate(matched_depth = pmin(spatial_units, 
    single_cell_units), eligible = matched_depth >= MIN_MATCHED_DEPTH, exclusion_reason = ifelse(eligible, 
    NA_character_, paste0("matched depth < ", MIN_MATCHED_DEPTH))) %>% arrange(sample)

save_table(matched_depths, file.path(TBL_DIR, "R2_Q3_BCR_matched_sample_depths"), if_csv = TRUE)

eligible_samples <- matched_depths %>% filter(eligible)

if (nrow(eligible_samples) < 2L) {
    stop("Fewer than two eligible matched samples; paired analysis is not possible.")
}

bcr_spatial_split <- split(bcr_spatial[bcr_spatial$sample %in% eligible_samples$sample, , drop = FALSE], 
    bcr_spatial$sample[bcr_spatial$sample %in% eligible_samples$sample])

bcr_sc_split <- split(bcr_sc[bcr_sc$sample %in% eligible_samples$sample, , drop = FALSE], bcr_sc$sample[bcr_sc$sample %in% 
    eligible_samples$sample])

set.seed(MASTER_SEED)

bcr_rep_list <- vector("list", N_ITER * nrow(eligible_samples) * 2L)

k <- 1L

for (iter in seq_len(N_ITER)) {
    for (i in seq_len(nrow(eligible_samples))) {
        smp <- eligible_samples$sample[i]
        target <- eligible_samples$matched_depth[i]
        sp_i <- bcr_spatial_split[[smp]]
        sc_i <- bcr_sc_split[[smp]]
        sp_sub <- sample_units(sp_i, "bin_id", target)
        sc_sub <- sample_units(sc_i, "sequence_id", target)
        sp_cat <- clone_category_counts(sp_sub, "clone_key", "bin_id") %>% mutate(iteration = iter, sample = smp, 
            platform = "Spatial", target_depth = target)
        sc_cat <- clone_category_counts(sc_sub, "clone_key", "sequence_id") %>% mutate(iteration = iter, 
            sample = smp, platform = "Single-cell", target_depth = target)
        bcr_rep_list[[k]] <- sp_cat
        bcr_rep_list[[k + 1L]] <- sc_cat
        k <- k + 2L
    }
}

bcr_reps <- bind_rows(bcr_rep_list) %>% mutate(category = factor(category, levels = EXPAND_LEVELS), platform = factor(platform, 
    levels = c("Spatial", "Single-cell")))

save_table(bcr_reps, file.path(TBL_DIR, "R2_Q3_BCR_matched_rarefaction_replicates"), if_csv = TRUE)

bcr_summary <- bcr_reps %>% group_by(sample, platform, target_depth, category) %>% summarise(estimate = median(fraction), 
    lower_95 = quantile(fraction, 0.025), upper_95 = quantile(fraction, 0.975), .groups = "drop")

save_table(bcr_summary, file.path(TBL_DIR, "R2_Q3_BCR_matched_rarefaction_summary"), if_csv = TRUE)

bcr_sample_metrics <- bcr_summary %>% transmute(sample, platform, metric = as.character(category), estimate, 
    lower_95, upper_95)

bcr_expanded <- bcr_reps %>% mutate(is_expanded = category != "Singleton") %>% group_by(iteration, sample, 
    platform) %>% summarise(fraction = sum(fraction[is_expanded]), .groups = "drop") %>% group_by(sample, 
    platform) %>% summarise(estimate = median(fraction), lower_95 = quantile(fraction, 0.025), upper_95 = quantile(fraction, 
    0.975), .groups = "drop") %>% mutate(metric = "Expanded total")

bcr_sample_metrics <- bind_rows(bcr_sample_metrics, bcr_expanded)

metrics <- c(EXPAND_LEVELS, "Expanded total")

bcr_stats <- bind_rows(lapply(metrics, function(m) {
    paired_stat(bcr_sample_metrics, m)
})) %>% mutate(p_adjust_BH = p.adjust(p_value, method = "BH"))

save_table(bcr_stats, file.path(TBL_DIR, "R2_Q3_BCR_matched_paired_stats"), if_csv = TRUE)

bcr_plot_summary <- bcr_summary %>% group_by(platform, category) %>% summarise(median = median(estimate), 
    q1 = quantile(estimate, 0.25), q3 = quantile(estimate, 0.75), .groups = "drop")

p_bcr_cat <- ggplot(bcr_plot_summary, aes(x = category, y = median * 100, fill = platform)) + geom_col(position = position_dodge(width = 0.76), 
    width = 0.68, color = "black", linewidth = 0.25) + geom_errorbar(aes(ymin = q1 * 100, ymax = q3 * 
    100), position = position_dodge(width = 0.76), width = 0.18, linewidth = 0.4) + scale_fill_manual(values = PLATFORM_COLS) + 
    labs(x = NULL, y = "Clonotype fraction (%)", title = "Matched-depth BCR clonotype distribution", 
        subtitle = sprintf("%d matched samples; 500 rarefactions/sample; bars: sample median, error bars: IQR", 
            nrow(eligible_samples)), fill = NULL) + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 20, 
    hjust = 1), legend.position = "top", plot.title = element_text(face = "bold"))

p_bcr_pair <- bcr_sample_metrics %>% filter(metric == "Expanded total") %>% ggplot(aes(x = platform, 
    y = estimate * 100, group = sample)) + geom_line(color = "grey65", linewidth = 0.45) + geom_point(aes(color = platform), 
    size = 2.1) + scale_color_manual(values = PLATFORM_COLS) + labs(x = NULL, y = "Expanded clonotypes (%)", 
    title = "Paired sample-level effect", subtitle = "Expanded = clone size >=2; points are medians across rarefactions", 
    color = NULL) + theme_classic(base_size = 10) + theme(legend.position = "none", plot.title = element_text(face = "bold"))

p_bcr <- p_bcr_cat + p_bcr_pair + plot_annotation(title = "BCR clone expansion after within-sample depth matching", 
    caption = "Statistical unit: sample. Spatial bins and single-cell receptor observations were matched within each sample.")

save_plot(p_bcr, file.path(FIG_DIR, "R2_Q3_BCR_matched_repeated_rarefaction"), width = 12, height = 5.5, 
    if_png = TRUE, if_pdf = TRUE)

#### TCR不同深度的重复抽样 ====

tcr <- readRDS(DATA_TCR_SMRT) %>% filter(topChains == "TRB", !is.na(sample), !is.na(cloneId), !is.na(Bin50)) %>% 
    mutate(clone_key = paste(sample, cloneId, sep = "__"))

tcr_total_bins <- n_distinct(tcr$Bin50)

tcr_depths <- c(1000L, 2000L, 3000L, 4000L)

tcr_depths <- tcr_depths[tcr_depths <= tcr_total_bins]

set.seed(MASTER_SEED + 1L)

tcr_rep_list <- vector("list", length(tcr_depths) * N_ITER)

k <- 1L

for (depth in tcr_depths) {
    for (iter in seq_len(N_ITER)) {
        sub <- sample_units(tcr, "Bin50", depth)
        tcr_rep_list[[k]] <- clone_category_counts(sub, "clone_key", "Bin50") %>% mutate(target_depth = depth, 
            iteration = iter)
        k <- k + 1L
    }
}

tcr_reps <- bind_rows(tcr_rep_list) %>% mutate(category = factor(category, levels = EXPAND_LEVELS))

save_table(tcr_reps, file.path(TBL_DIR, "R2_Q3_TCR_rarefaction_replicates"), if_csv = TRUE)

tcr_summary <- tcr_reps %>% group_by(target_depth, category) %>% summarise(estimate = median(fraction), 
    lower_95 = quantile(fraction, 0.025), upper_95 = quantile(fraction, 0.975), .groups = "drop")

save_table(tcr_summary, file.path(TBL_DIR, "R2_Q3_TCR_rarefaction_summary"), if_csv = TRUE)

p_tcr <- ggplot(tcr_summary, aes(x = target_depth, y = estimate * 100, color = category, fill = category)) + 
    geom_ribbon(aes(ymin = lower_95 * 100, ymax = upper_95 * 100), alpha = 0.15, color = NA) + geom_line(linewidth = 0.85) + 
    geom_point(size = 1.8) + scale_color_manual(values = EXPAND_COLS) + scale_fill_manual(values = EXPAND_COLS) + 
    scale_x_continuous(breaks = tcr_depths) + labs(x = "Sampled spatial TRB bins", y = "Clonotype fraction (%)", 
    title = "TCR(TRB) depth-sensitivity analysis", subtitle = sprintf("Total available bins = %s; 500 rarefactions/depth", 
        format(tcr_total_bins, big.mark = ",")), caption = "Lines show medians and ribbons show 95% empirical intervals; descriptive analysis of pooled data.", 
    color = NULL, fill = NULL) + theme_classic(base_size = 10) + theme(legend.position = "top", plot.title = element_text(face = "bold"))

save_plot(p_tcr, file.path(FIG_DIR, "R2_Q3_TCR_depth_sensitivity"), width = 7.2, height = 5.5, if_png = TRUE, 
    if_pdf = TRUE)
