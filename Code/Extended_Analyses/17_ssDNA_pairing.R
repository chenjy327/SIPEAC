#### SIPEAC在ssDNA与scRNA参考下的比较 ====
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/utils_io.R")
source("Code/Extended_Analyses/utils_plot.R")
input_dir <- file.path(data_dir, "hzy")
processed_dir <- file.path(dirs$processed, "hzy")
table_dir <- file.path(dirs$tables, "hzy", "ssDNA")
figure_dir <- file.path(dirs$figures, "hzy")
library(tidyr)
library(data.table)
library(dplyr)
library(glue)
library(ggplot2)
library(ggpubr)

changeo_out <- readRDS(file.path(input_dir, "ssDNA/merge4.rds"))

df <- read.csv(file.path(input_dir, "ssDNA/a_compare.csv"), header = T, row.names = 1)

sample = unique(df$cur_sample)

base_dir <- file.path(input_dir, "ssDNA/")

bin_size <- 110

allelic_summary <- fread(file.path(input_dir, "ssDNA/summary_e100.tsv")) %>% as.data.frame()

sipeac_results <- list()

for (i in seq_len(nrow(allelic_summary))) {
    cs <- allelic_summary$cur_sample[i]
    bs <- allelic_summary$bin_size[i]
    f <- glue::glue(file.path(input_dir, "ssDNA/bin{bs}/{cs}_final_e100.tsv"))
    if (!file.exists(f)) 
        next
    df <- fread(f)
    df$cur_sample <- cs
    df$bin_size <- bs
    sipeac_results[[glue::glue("{cs}_bin{bs}")]] <- df
}

sipeac_results <- bind_rows(sipeac_results)

samples <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

results_new <- list()

results_old <- list()

for (cur_sample in sample) {
    cat("\n============================\n")
    cat("Processing:", cur_sample, "\n")
    cat("============================\n")
    bcr_path <- file.path(base_dir, cur_sample, "bcr_cell_mapping_full.csv")
    if (file.exists(bcr_path)) {
        tryCatch({
            bcr_cell_mapping_full <- read.csv(bcr_path)
            changeo_cur_sample <- changeo_out %>% filter(sample == cur_sample, locus %in% c("IGH", 
                "IGK", "IGL")) %>% select(locus, clone_id, junction, junction_aa) %>% distinct() %>% 
                mutate(clone_name = paste(locus, clone_id, sep = "_")) %>% select(junction, clone_name)
            bcr_mapped <- bcr_cell_mapping_full %>% left_join(changeo_cur_sample, by = "junction", 
                relationship = "many-to-many")
            clone_abund <- bcr_mapped %>% filter(!is.na(matched_cell_id), !is.na(clone_name)) %>% 
                mutate(chain = case_when(grepl("^IGH", clone_name) ~ "H", grepl("^(IGK|IGL)", clone_name) ~ 
                  "L", TRUE ~ NA_character_)) %>% filter(!is.na(chain)) %>% mutate(sample = cur_sample) %>% 
                filter(sample == !!cur_sample) %>% distinct(sample, matched_cell_id, chain, clone_name, 
                junction) %>% count(sample, matched_cell_id, chain, clone_name, name = "n_reads")
            dominant <- clone_abund %>% group_by(sample, matched_cell_id, chain) %>% filter(n_reads == 
                max(n_reads)) %>% ungroup()
            H_by_bin <- dominant %>% filter(chain == "H") %>% select(sample, matched_cell_id, Hclone = clone_name, 
                n_reads_H = n_reads)
            L_by_bin <- dominant %>% filter(chain == "L") %>% select(sample, matched_cell_id, Lclone = clone_name, 
                n_reads_L = n_reads)
            gold_standard <- inner_join(H_by_bin, L_by_bin, by = c("sample", "matched_cell_id"), 
                relationship = "many-to-many") %>% mutate(pair_name = paste(Hclone, Lclone, sep = "_")) %>% 
                select(sample, matched_cell_id, Hclone, Lclone, pair_name, n_reads_H, n_reads_L)
            dir.create(dirname(file.path(table_dir, cur_sample, "cellbin_gold.csv")), recursive = TRUE, 
                showWarnings = FALSE)
            write.csv(gold_standard, file.path(table_dir, cur_sample, "cellbin_gold.csv"), row.names = FALSE)
            sipeac_pairs <- sipeac_results %>% filter(cur_sample == !!cur_sample, as.character(bin_size) == 
                as.character(!!bin_size)) %>% mutate(pair_name = paste(Hclone, Lclone, sep = "_")) %>% 
                distinct(Hclone, Lclone, pair_name)
            gold_H_all <- gold_standard %>% distinct(Hclone)
            gold_L_all <- gold_standard %>% distinct(Lclone)
            gold_pair_set <- unique(gold_standard$pair_name)
            output_pair_test <- sipeac_pairs %>% filter(Hclone %in% gold_H_all$Hclone)
            n_H_clone <- nrow(gold_H_all)
            n_L_clone <- nrow(gold_L_all)
            n_H_clone_pred <- n_distinct(sipeac_pairs$Hclone)
            n_test <- gold_standard %>% filter(Hclone %in% sipeac_pairs$Hclone) %>% nrow()
            n_pred <- nrow(output_pair_test)
            n_correct <- sum(output_pair_test$pair_name %in% gold_pair_set)
            acc <- if (n_pred > 0) 
                n_correct/n_pred
            else NA_real_
            eff <- n_H_clone_pred/n_H_clone
            results_new[[cur_sample]] <- data.frame(cur_sample = cur_sample, bin_size = bin_size, 
                n_H_clone = n_H_clone, n_L_clone = n_L_clone, n_H_clone_pred = n_H_clone_pred, n_test = n_test, 
                n_pred = n_pred, n_correct = n_correct, acc = acc, eff = eff)
        }, error = function(e) {
            cat("新金标准出错:", conditionMessage(e), "\n")
        })
    }
    else {
        cat("缺少 bcr_cell_mapping_full.csv，跳新金标准\n")
    }
    old_path <- glue(file.path(input_dir, "ssDNA/bin{bin_size}/{cur_sample}.tsv"))
    if (file.exists(old_path)) {
        tryCatch({
            gold_standard <- read.table(old_path, sep = "\t", header = TRUE, row.names = 1)
            sipeac_pairs <- sipeac_results %>% filter(cur_sample == !!cur_sample, as.character(bin_size) == 
                as.character(!!bin_size)) %>% mutate(pair_name = paste(Hclone, Lclone, sep = "_")) %>% 
                distinct(Hclone, Lclone, pair_name)
            gold_H_all <- gold_standard %>% distinct(Hclone)
            gold_L_all <- gold_standard %>% distinct(Lclone)
            gold_pair_set <- unique(gold_standard$pair_name)
            output_pair_test <- sipeac_pairs %>% filter(Hclone %in% gold_H_all$Hclone)
            n_H_clone <- nrow(gold_H_all)
            n_L_clone <- nrow(gold_L_all)
            n_H_clone_pred <- n_distinct(sipeac_pairs$Hclone)
            n_test <- gold_standard %>% filter(Hclone %in% sipeac_pairs$Hclone) %>% nrow()
            n_pred <- nrow(output_pair_test)
            n_correct <- sum(output_pair_test$pair_name %in% gold_pair_set)
            acc <- if (n_pred > 0) 
                n_correct/n_pred
            else NA_real_
            eff <- n_H_clone_pred/n_H_clone
            results_old[[cur_sample]] <- data.frame(cur_sample = cur_sample, bin_size = bin_size, 
                n_H_clone = n_H_clone, n_L_clone = n_L_clone, n_H_clone_pred = n_H_clone_pred, n_test = n_test, 
                n_pred = n_pred, n_correct = n_correct, acc = acc, eff = eff)
        }, error = function(e) {
            cat("老金标准出错:", conditionMessage(e), "\n")
        })
    }
}

results_new <- bind_rows(results_new)

results_old <- bind_rows(results_old)

dir.create(dirname(file.path(table_dir, glue("acc_eff_new_gold_bin{bin_size}.csv"))), recursive = TRUE, 
    showWarnings = FALSE)

write.csv(results_new, file.path(table_dir, glue("acc_eff_new_gold_bin{bin_size}.csv")), row.names = FALSE)

dir.create(dirname(file.path(table_dir, glue("acc_eff_old_gold_bin{bin_size}.csv"))), recursive = TRUE, 
    showWarnings = FALSE)

write.csv(results_old, file.path(table_dir, glue("acc_eff_old_gold_bin{bin_size}.csv")), row.names = FALSE)

print(results_new)

print(results_old)

df_plot <- bind_rows(results_new %>% mutate(gold = "New gold standard"), results_old %>% mutate(gold = "Old gold standard")) %>% 
    filter(!is.na(acc)) %>% mutate(gold = factor(gold, levels = c("Old gold standard", "New gold standard")))

p <- ggplot(df_plot, aes(x = gold, y = acc, fill = gold)) + geom_boxplot(width = 0.55, outlier.shape = NA, 
    alpha = 0.55, color = "grey25", size = 0.5) + geom_jitter(width = 0.12, size = 1.8, alpha = 0.65, 
    color = "grey15", shape = 16) + stat_summary(fun = median, geom = "point", shape = 23, size = 3.2, 
    fill = "white", color = "black", stroke = 0.8) + stat_compare_means(method = "wilcox.test", 
    label = "p.format", label.y = max(df_plot$acc, na.rm = TRUE) * 1.05, size = 4.2) + scale_fill_manual(values = c(`Old gold standard` = "#F4A582", 
    `New gold standard` = "#92C5DE")) + labs(x = NULL, y = "Accuracy (acc)", title = "Accuracy: New vs Old gold standard", 
    subtitle = glue("mean(New) = {round(mean(results_new$acc, na.rm=TRUE),3)} | ", "mean(Old) = {round(mean(results_old$acc, na.rm=TRUE),3)}")) + 
    theme_bw(base_size = 13) + theme(legend.position = "none", plot.title = element_text(face = "bold"), 
    plot.subtitle = element_text(color = "grey35"), axis.text.x = element_text(size = 12, color = "black"), 
    panel.grid.minor = element_blank())

print(p)

dir.create(dirname(file.path(figure_dir, glue("acc_boxplot_new_vs_old_bin{bin_size}.png"))), recursive = TRUE, 
    showWarnings = FALSE)

ggsave(file.path(figure_dir, glue("acc_boxplot_new_vs_old_bin{bin_size}.png")), p, width = 6.5, 
    height = 5.5, dpi = 300)

dir.create(dirname(file.path(figure_dir, glue("acc_boxplot_new_vs_old_bin{bin_size}.pdf"))), recursive = TRUE, 
    showWarnings = FALSE)

ggsave(file.path(figure_dir, glue("acc_boxplot_new_vs_old_bin{bin_size}.pdf")), p, width = 6.5, 
    height = 5.5)

