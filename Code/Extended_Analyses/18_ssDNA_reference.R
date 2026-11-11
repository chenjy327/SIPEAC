#### ssDNA直接配对对照scRNA参考 ====
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/R/io.R")
source("Code/Extended_Analyses/R/plot.R")
input_dir <- file.path(data_dir, "hzy")
processed_dir <- file.path(dirs$processed, "hzy")
table_dir <- file.path(dirs$tables, "hzy", "ssDNA")
figure_dir <- file.path(dirs$figures, "hzy")
library(tidyr)
library(data.table)
library(dplyr)
library(glue)

changeo_out <- readRDS(file.path(input_dir, "ssDNA/merge4.rds"))

df <- read.csv(file.path(input_dir, "ssDNA/a_compare.csv"), header = T, row.names = 1)

sample <- unique(df$cur_sample)

base_dir <- file.path(input_dir, "ssDNA/")

bin_size <- 20

results_ssdna <- list()

for (cur_sample in sample) {
    cat("Processing:", cur_sample, "\n")
    old_path <- glue(file.path(input_dir, "ssDNA/bin{bin_size}/{cur_sample}.tsv"))
    if (!file.exists(old_path)) {
        cat("缺少老金标准文件，跳过\n")
        next
    }
    gold_standard <- tryCatch(read.table(old_path, sep = "\t", header = TRUE, row.names = 1), error = function(e) {
        cat("读老金标准出错:", conditionMessage(e), "\n")
        NULL
    })
    if (is.null(gold_standard) || nrow(gold_standard) == 0) 
        next
    bcr_path <- file.path(base_dir, cur_sample, "bcr_cell_mapping_full.csv")
    if (!file.exists(bcr_path)) {
        cat("缺少 bcr_cell_mapping_full.csv，跳过\n")
        next
    }
    pred_pairs <- tryCatch({
        bcr_cell_mapping_full <- read.csv(bcr_path)
        changeo_cur_sample <- changeo_out %>% filter(sample == cur_sample, locus %in% c("IGH", "IGK", 
            "IGL")) %>% select(locus, clone_id, junction, junction_aa) %>% distinct() %>% mutate(clone_name = paste(locus, 
            clone_id, sep = "_")) %>% select(junction, clone_name)
        bcr_mapped <- bcr_cell_mapping_full %>% left_join(changeo_cur_sample, by = "junction", relationship = "many-to-many")
        clone_abund <- bcr_mapped %>% filter(!is.na(matched_cell_id), !is.na(clone_name)) %>% mutate(chain = case_when(grepl("^IGH", 
            clone_name) ~ "H", grepl("^(IGK|IGL)", clone_name) ~ "L", TRUE ~ NA_character_)) %>% 
            filter(!is.na(chain)) %>% mutate(sample = cur_sample) %>% filter(sample == !!cur_sample) %>% 
            distinct(sample, matched_cell_id, chain, clone_name, junction) %>% count(sample, matched_cell_id, 
            chain, clone_name, name = "n_reads")
        dominant <- clone_abund %>% group_by(sample, matched_cell_id, chain) %>% filter(n_reads == 
            max(n_reads)) %>% ungroup()
        H_by_bin <- dominant %>% filter(chain == "H") %>% select(sample, matched_cell_id, Hclone = clone_name, 
            n_reads_H = n_reads)
        L_by_bin <- dominant %>% filter(chain == "L") %>% select(sample, matched_cell_id, Lclone = clone_name, 
            n_reads_L = n_reads)
        inner_join(H_by_bin, L_by_bin, by = c("sample", "matched_cell_id"), relationship = "many-to-many") %>% 
            mutate(pair_name = paste(Hclone, Lclone, sep = "_")) %>% select(sample, matched_cell_id, 
            Hclone, Lclone, pair_name, n_reads_H, n_reads_L)
    }, error = function(e) {
        cat("构建 ssDNA 配对出错:", conditionMessage(e), "\n")
        NULL
    })
    if (is.null(pred_pairs) || nrow(pred_pairs) == 0) 
        next
    gold_H_all <- gold_standard %>% distinct(Hclone)
    gold_L_all <- gold_standard %>% distinct(Lclone)
    gold_pair_set <- unique(gold_standard$pair_name)
    output_pair_test <- pred_pairs %>% filter(Hclone %in% gold_H_all$Hclone)
    n_H_clone <- nrow(gold_H_all)
    n_L_clone <- nrow(gold_L_all)
    n_H_clone_pred <- n_distinct(pred_pairs$Hclone)
    n_test <- gold_standard %>% filter(Hclone %in% pred_pairs$Hclone) %>% nrow()
    n_pred <- nrow(output_pair_test)
    n_correct <- sum(output_pair_test$pair_name %in% gold_pair_set)
    acc <- if (n_pred > 0) 
        n_correct/n_pred
    else NA_real_
    eff <- n_H_clone_pred/n_H_clone
    results_ssdna[[cur_sample]] <- data.frame(cur_sample = cur_sample, bin_size = bin_size, n_H_clone = n_H_clone, 
        n_L_clone = n_L_clone, n_H_clone_pred = n_H_clone_pred, n_test = n_test, n_pred = n_pred, 
        n_correct = n_correct, acc = acc, eff = eff)
}

results_ssdna <- bind_rows(results_ssdna)

dir.create(dirname(file.path(table_dir, glue("acc_eff_ssDNA_vs_scRNA_gold_bin.csv"))), recursive = TRUE, 
    showWarnings = FALSE)

write.csv(results_ssdna, file.path(table_dir, glue("acc_eff_ssDNA_vs_scRNA_gold_bin.csv")), row.names = FALSE)

print(results_ssdna)

