#### 抗原抗体邻域置换和空间图 ====
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/R/io.R")
source("Code/Extended_Analyses/R/plot.R")
input_dir <- file.path(data_dir, "hzy")
processed_dir <- file.path(dirs$processed, "hzy")
table_dir <- file.path(dirs$tables, "hzy")
figure_dir <- file.path(dirs$figures, "hzy")
library(data.table)
library(dplyr)
library(tidyr)
library(glue)
library(stringr)
library(ggplot2)
library(ggpubr)
library(scales)
library(Seurat)
library(spdep)
library(magrittr)
library(purrr)
library(cli)
source("Code/Extended_Analyses/R/expand_n_layer.R")

calc_contact_matrix <- function(data_df, neighbor_list) {
    all_groups <- sort(unique(data_df$group))
    n_group <- length(all_groups)
    count_mat <- matrix(0, nrow = n_group, ncol = n_group, dimnames = list(all_groups, all_groups))
    for (i in 1:length(neighbor_list)) {
        current_cell_group <- data_df$group[i]
        neighbors_indices <- neighbor_list[[i]]
        if (length(neighbors_indices) == 0 || neighbors_indices[1] == 0) 
            next
        neighbor_groups <- data_df$group[neighbors_indices]
        for (g in neighbor_groups) {
            count_mat[current_cell_group, g] <- count_mat[current_cell_group, g] + 1
        }
    }
    return(count_mat)
}

bcr_smrt <- readRDS(file.path(input_dir, "result/03.Ag_Ab/merged_bcr.rds"))

changeo_merge <- readRDS(file.path(input_dir, "result/03.Ag_Ab/merge4.rds"))

ab_15_clone <- read.csv(file.path(input_dir, "result/03.Ag_Ab/IFI30_Ab_merge.csv"), row.names = 1)

ab_15_clone_spatial <- tidyr::separate_rows(ab_15_clone, sample, sep = ",")

ab_15_clone_spatial <- dplyr::arrange(ab_15_clone_spatial, desc(name))

pos_matrix <- matrix(nrow = nrow(ab_15_clone_spatial), ncol = 5, dimnames = list(ab_15_clone_spatial$name, 
    c("name", "clone_name", "sample", "z_score", "p_value")))

#### 邻域统计 ====
for (cur_i in 1:nrow(ab_15_clone_spatial)) {
    print(cur_i)
    cur_clone <- ab_15_clone_spatial$clone_name[cur_i]
    cur_df <- changeo_merge[changeo_merge$clone_name == cur_clone & changeo_merge$type == "Spatial", 
        ]
    cur_seq <- cur_df$sequence_id
    cur_sample <- ab_15_clone_spatial$sample[cur_i]
    cur_smrt_anno <- bcr_smrt[bcr_smrt$sequence_id %in% cur_df$sequence_id, ]
    ig_name <- ab_15_clone_spatial$name[cur_i]
    cur_smrt_anno$mab_name <- ab_15_clone_spatial$name[cur_i]
    obj <- readRDS(glue(file.path(input_dir, "spatial_objects/{cur_sample}.rds")))
    obj@meta.data %<>% mutate(Location = case_when(Bin_Region %in% c("Paratumor_tissue", "Paratumor_side_of_Margin_area") ~ 
        "Paratumor", Bin_Region %in% c("Tumor_capsule") ~ "Invasive_zone", Bin_Region %in% c("Tumor_side_of_Margin_area", 
        "Tumor_tissue") ~ "Tumor"))
    obj@meta.data$TLS <- obj@meta.data$TLS_raw
    tmp_bcr_cloneid <- cur_smrt_anno
    tmp_bcr_cloneid <- table(tmp_bcr_cloneid$Bin50, tmp_bcr_cloneid$mab_name) %>% as.data.frame()
    colnames(tmp_bcr_cloneid) <- c("BinID", "mab_name", "clone")
    tmp_bcr_cloneid <- tmp_bcr_cloneid[tmp_bcr_cloneid$clone != 0, ]
    tmp_bcr_cloneid <- dplyr::distinct(tmp_bcr_cloneid, BinID, .keep_all = T)
    df_Ab <- filter(obj@meta.data, bin %in% tmp_bcr_cloneid$BinID) %>% dplyr::select(c(row, col)) %>% 
        mutate(group = "Ab")
    df_IFI30 <- FetchData(obj, vars = c("IFI30", "row", "col")) %>% filter(IFI30 != 0) %>% dplyr::select(c(row, 
        col)) %>% mutate(group = "Ag")
    df_bg <- filter(obj@meta.data, !bin %in% c(rownames(df_IFI30), rownames(df_Ab))) %>% dplyr::select(c(row, 
        col)) %>% mutate(group = "Background")
    df_merge <- bind_rows(list(df_Ab, df_IFI30, df_bg))
    set.seed(123)
    coords <- as.matrix(df_merge[, c("row", "col")])
    knn <- knearneigh(coords, k = 6)
    nb <- knn2nb(knn)
    observed_counts <- calc_contact_matrix(df_merge, nb)
    n_permutations <- 999
    all_groups <- rownames(observed_counts)
    n_g <- length(all_groups)
    perm_results <- array(0, dim = c(n_g, n_g, n_permutations), dimnames = list(all_groups, all_groups, 
        NULL))
    message("正在进行置换检验以构建背景分布...")
    for (p in 1:n_permutations) {
        print(p)
        perm_data <- df_merge
        perm_data$group <- sample(perm_data$group)
        perm_results[, , p] <- calc_contact_matrix(perm_data, nb)
    }
    expected_mean <- apply(perm_results, c(1, 2), mean)
    expected_sd <- apply(perm_results, c(1, 2), sd)
    z_score_matrix <- (observed_counts - expected_mean)/(expected_sd + 1e-05)
    p_value_matrix <- matrix(0, nrow = n_g, ncol = n_g, dimnames = list(all_groups, all_groups))
    for (i in all_groups) {
        for (j in all_groups) {
            higher_than_obs <- sum(perm_results[i, j, ] >= observed_counts[i, j])
            p_value_matrix[i, j] <- (higher_than_obs + 1)/(n_permutations + 1)
        }
    }
    print("--- 空间邻域富集 Z-score 矩阵 (正值代表显著接近) ---")
    print(round(z_score_matrix, 2))
    print("--- 空间邻域富集 P-value 矩阵 (越接近 0 代表越显著富集) ---")
    print(round(p_value_matrix, 4))
    pos_matrix[cur_i, ] <- c(ig_name, cur_clone, cur_sample, z_score_matrix[1, 2], p_value_matrix[1, 
        2])
}

dir.create(dirname(file.path(table_dir, "antigen_neighborhood/IFI30_Neighborhood_batch1.csv")), 
    recursive = TRUE, showWarnings = FALSE)

write.csv(pos_matrix, file.path(table_dir, "antigen_neighborhood/IFI30_Neighborhood_batch1.csv"))

#### 原位图 ====
for (cur_i in 1:nrow(ab_15_clone_spatial)) {
    print(cur_i)
    cur_clone <- ab_15_clone_spatial$clone_name[cur_i]
    cur_df <- changeo_merge[changeo_merge$clone_name == cur_clone & changeo_merge$type == "Spatial", 
        ]
    cur_seq <- cur_df$sequence_id
    cur_sample <- ab_15_clone_spatial$sample[cur_i]
    cur_smrt_anno <- bcr_smrt[bcr_smrt$sequence_id %in% cur_df$sequence_id, ]
    ig_name <- ab_15_clone_spatial$name[cur_i]
    cur_smrt_anno$mab_name <- ab_15_clone_spatial$name[cur_i]
    obj <- readRDS(glue(file.path(input_dir, "spatial_objects/{cur_sample}.rds")))
    obj@meta.data %<>% mutate(Location = case_when(Bin_Region %in% c("Paratumor_tissue", "Paratumor_side_of_Margin_area") ~ 
        "Paratumor", Bin_Region %in% c("Tumor_capsule") ~ "Invasive_zone", Bin_Region %in% c("Tumor_side_of_Margin_area", 
        "Tumor_tissue") ~ "Tumor"))
    obj@meta.data$TLS <- obj@meta.data$TLS_raw
    tmp_bcr_cloneid <- cur_smrt_anno
    tmp_bcr_cloneid <- table(tmp_bcr_cloneid$Bin50, tmp_bcr_cloneid$mab_name) %>% as.data.frame()
    colnames(tmp_bcr_cloneid) <- c("BinID", "mab_name", "clone")
    tmp_bcr_cloneid <- tmp_bcr_cloneid[tmp_bcr_cloneid$clone != 0, ]
    tmp_bcr_cloneid <- dplyr::distinct(tmp_bcr_cloneid, BinID, .keep_all = T)
    df_Ab <- filter(obj@meta.data, bin %in% tmp_bcr_cloneid$BinID) %>% dplyr::select(c(row, col))
    df_IFI30 <- FetchData(obj, vars = c("IFI30", "row", "col")) %>% filter(IFI30 != 0) %>% dplyr::select(c(row, 
        col))
    df_location <- df_Ab %>% dplyr::mutate(dplyr::across(c(col, row), as.numeric)) %>% dplyr::distinct(row, 
        col)
    Ab_expand <- df_location %>% st_expand_n(n = 1, x = "col", y = "row", expand_n = "expand_n", 
        mode = 2)
    Ab_bin <- left_join(obj@meta.data, Ab_expand, by = c("row", "col")) %>% filter(expand_n %in% 
        c(0, 1)) %>% pull(bin)
    df_location <- df_IFI30 %>% dplyr::mutate(dplyr::across(c(col, row), as.numeric)) %>% dplyr::distinct(row, 
        col)
    Ag_expand <- df_location %>% st_expand_n(n = 1, x = "col", y = "row", expand_n = "expand_n", 
        mode = 2)
    Ag_bin <- left_join(obj@meta.data, Ag_expand, by = c("row", "col")) %>% filter(expand_n %in% 
        c(0, 1)) %>% pull(bin)
    obj@meta.data %<>% mutate(label = case_when(bin %in% Ab_bin ~ ig_name, bin %in% Ag_bin ~ "IFI30", 
        TRUE ~ Location))
    dir.create(dirname(glue(file.path(figure_dir, "antigen_neighborhood/Neighborhood"))), recursive = TRUE, 
        showWarnings = FALSE)
    SpatialDimPython(obj = obj, prefix = glue(file.path(figure_dir, "antigen_neighborhood/Neighborhood")), 
        plot_item = "label", plot_order = c(ig_name, "IFI30", "Tumor", "Invasive_zone", "Paratumor"), 
        tmp_color = c("#EE3743", "#9CCCE8", "#E9ECEB", "#FFE8BF", "#F9F7EE"), name = glue("{cur_sample}_{ig_name}_bin50"))
}

