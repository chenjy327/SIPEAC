#### 抗体周围IFI30表达 ====
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/utils_io.R")
source("Code/Extended_Analyses/utils_plot.R")
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
library(patchwork)
library(readxl)
library(magrittr)
library(purrr)
library(cli)
source("Code/Extended_Analyses/utils_spatial_expand.R")

merge_metadata <- read.table(file.path(input_dir, "data/Bin50_classic_celltype_metadata_merge.txt"), 
    header = T)

bcr_smrt_merge <- readRDS(file.path(input_dir, "data/merged_bcr.rds"))

changeo_merge <- readRDS(file.path(input_dir, "data/merge4.rds"))

changeo_merge$clone_name2 <- paste(changeo_merge$sample, changeo_merge$clone_name, sep = "_")

Ab_97 <- readxl::read_xlsx(file.path(input_dir, "result/03.Ag_Ab/97Ab.xlsx"))

table(Ab_97$`IFI30 反应性抗体（ELISA>0.5且BLI>0.5）`)

IFI30_neg <- Ab_97 %>% filter(`IFI30 反应性抗体（ELISA>0.5且BLI>0.5）` == "×") %>% tidyr::separate_rows(sample_st, 
    sep = ";") %>% mutate(clone_name2 = paste(sample_st, clone_name, sep = "_")) %>% filter(clone_name2 %in% 
    changeo_merge$clone_name2)

IFI30_pos <- Ab_97 %>% filter(`IFI30 反应性抗体（ELISA>0.5且BLI>0.5）` == "√") %>% tidyr::separate_rows(sample_st, 
    sep = ";") %>% mutate(clone_name2 = paste(sample_st, clone_name, sep = "_")) %>% filter(clone_name2 %in% 
    changeo_merge$clone_name2)

dir.create(dirname(file.path(table_dir, "IFI30_distance/IFI30_neg.csv")), recursive = TRUE, showWarnings = FALSE)

write.csv(IFI30_neg, file.path(table_dir, "IFI30_distance/IFI30_neg.csv"))

dir.create(dirname(file.path(table_dir, "IFI30_distance/IFI30_pos.csv")), recursive = TRUE, showWarnings = FALSE)

write.csv(IFI30_pos, file.path(table_dir, "IFI30_distance/IFI30_pos.csv"))

#### 阴性抗体周围 ====
mean_list <- list()

pos_list <- list()

for (cur_i in 1:nrow(IFI30_neg)) {
    print(cur_i)
    cur_clone <- IFI30_neg$clone_name2[cur_i]
    cur_df <- changeo_merge[changeo_merge$clone_name2 == cur_clone & changeo_merge$type == "Spatial", 
        ]
    cur_seq <- cur_df$sequence_id
    cur_sample <- IFI30_neg$sample_st[cur_i]
    cur_smrt_anno <- bcr_smrt_merge[bcr_smrt_merge$sequence_id %in% cur_df$sequence_id, ]
    ig_name <- IFI30_neg$IG_name[cur_i]
    cur_smrt_anno$mab_name <- ig_name
    obj <- readRDS(glue(file.path(input_dir, "spatial_objects/{cur_sample}.rds")))
    df <- merge_metadata %>% filter(SampleID == cur_sample)
    obj@meta.data <- obj@meta.data %>% dplyr::select(c(BinID, bin)) %>% left_join(df, by = "BinID") %>% 
        mutate(Location = case_when(Bin_Region %in% c("Paratumor_tissue", "Paratumor_side_of_Margin_area") ~ 
            "Paratumor", Bin_Region %in% c("Tumor_capsule") ~ "Invasive_zone", Bin_Region %in% c("Tumor_side_of_Margin_area", 
            "Tumor_tissue") ~ "Tumor"))
    obj@meta.data$TLS_Final <- ifelse(obj@meta.data$TLS_maturity %in% c("Rare_Naive_GCB", "Filtered"), 
        glue("ST{cur_sample}_NA"), obj@meta.data$TLS_Final)
    rownames(obj@meta.data) <- obj@meta.data$bin
    obj@meta.data %<>% mutate(Location = case_when(Bin_Region %in% c("Paratumor_tissue", "Paratumor_side_of_Margin_area") ~ 
        "Paratumor", Bin_Region %in% c("Tumor_capsule") ~ "Invasive_zone", Bin_Region %in% c("Tumor_side_of_Margin_area", 
        "Tumor_tissue") ~ "Tumor"))
    obj@meta.data$TLS <- obj@meta.data$TLS_raw
    tmp_bcr_cloneid <- cur_smrt_anno
    tmp_bcr_cloneid <- table(tmp_bcr_cloneid$Bin50, tmp_bcr_cloneid$mab_name) %>% as.data.frame()
    colnames(tmp_bcr_cloneid) <- c("BinID", "mab_name", "clone")
    tmp_bcr_cloneid <- tmp_bcr_cloneid[tmp_bcr_cloneid$clone != 0, ]
    tmp_bcr_cloneid <- dplyr::distinct(tmp_bcr_cloneid, BinID, .keep_all = T)
    obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = c(bin = "BinID"))
    rownames(obj@meta.data) <- obj@meta.data$bin
    obj@meta.data$mab_name <- as.character(obj@meta.data$mab_name)
    obj@meta.data$mab_name <- ifelse(is.na(obj@meta.data$mab_name), "other", obj@meta.data$mab_name)
    obj@meta.data$mab_name <- ifelse(obj@meta.data$TLS != glue("ST{cur_sample}_NA") & obj@meta.data$mab_name == 
        "other", "TLS", obj@meta.data$mab_name)
    obj@meta.data$mab_name <- ifelse(obj@meta.data$mab_name == "other", obj@meta.data$Location, 
        obj@meta.data$mab_name)
    df2 <- obj@meta.data[obj@meta.data$mab_name %in% c(ig_name), ]
    df_location <- df2 %>% dplyr::mutate(dplyr::across(c(col, row), as.numeric)) %>% dplyr::distinct(row, 
        col)
    df_expand <- df_location %>% st_expand_n(n = 9, x = "col", y = "row", expand_n = "expand_n", 
        mode = 2)
    df_expand <- df_expand %>% dplyr::select(col, row, Ig_center = group.by, distance = expand_n)
    table(df_expand$Ig_center, df_expand$distance)
    obj@meta.data <- obj@meta.data %>% dplyr::mutate(dplyr::across(c(col, row), as.numeric)) %>% 
        dplyr::select(-dplyr::any_of(c("Ig_center", "distance"))) %>% dplyr::left_join(df_expand, 
        by = c("row", "col"))
    rownames(obj@meta.data) <- obj$bin
    obj@meta.data$distance <- ifelse(is.na(obj@meta.data$distance), "≥10", obj@meta.data$distance)
    df_exp <- FetchData(obj, vars = c("IFI30", "distance")) %>% setnames(c(1, 2), c("exp", "group")) %>% 
        dplyr::select(exp, group) %>% filter(group != 0)
    df_mean <- df_exp %>% group_by(group) %>% summarise(mean_exp = mean(exp))
    df_mean$group <- factor(df_mean$group, levels = c(as.character(1:9), "≥10"))
    df_mean$Path <- "IFI30 exp"
    p_mean <- ggplot(df_mean, aes_string("group", "mean_exp", color = "Path", group = "Path")) + 
        geom_point(size = 2) + geom_line(size = 1) + labs(title = "IFI30 expression", x = "", y = " ") + 
        scale_color_manual(values = mycolor_new) + theme_classic(base_line_size = 1) + theme(plot.title = element_text(size = 15, 
        colour = "black", hjust = 0.5), axis.title.y = element_text(size = 15, color = "black", 
        face = "bold", vjust = 1.9, hjust = 0.5, angle = 90), legend.position = "none", axis.text.x = element_text(size = 13, 
        color = "black", face = "bold", vjust = 1, hjust = 1, angle = 0), axis.text.y = element_text(size = 13, 
        color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))
    df_pos <- df_exp %>% group_by(group) %>% summarise(total_cells = n(), positive_cells = sum(exp > 
        0, na.rm = TRUE), proportion = positive_cells/total_cells * 100)
    df_pos$group <- factor(df_pos$group, levels = c(as.character(1:9), "≥10"))
    df_pos$Path <- "IFI30 density"
    p_pos <- ggplot(df_pos, aes_string("group", "proportion", color = "Path", group = "Path")) + 
        geom_point(size = 2) + geom_line(size = 1) + labs(title = "IFI30 ratio", x = "", y = " ") + 
        scale_color_manual(values = mycolor_new) + theme_classic(base_line_size = 1) + theme(plot.title = element_text(size = 15, 
        colour = "black", hjust = 0.5), axis.title.y = element_text(size = 15, color = "black", 
        face = "bold", vjust = 1.9, hjust = 0.5, angle = 90), legend.position = "none", axis.text.x = element_text(size = 13, 
        color = "black", face = "bold", vjust = 1, hjust = 1, angle = 0), axis.text.y = element_text(size = 13, 
        color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))
    p_pos
    mean_list[[glue("{ig_name}_{cur_sample}")]] <- df_mean %>% mutate(sample = cur_sample, igname = ig_name)
    pos_list[[glue("{ig_name}_{cur_sample}")]] <- df_pos %>% mutate(sample = cur_sample, igname = ig_name)
}

mean_merge <- bind_rows(mean_list)

pos_merge <- bind_rows(pos_list)

dir.create(dirname(file.path(table_dir, "IFI30_distance/IFI30Ab_neg/mean_merge9_97Ab.csv")), recursive = TRUE, 
    showWarnings = FALSE)

write.csv(mean_merge, file.path(table_dir, "IFI30_distance/IFI30Ab_neg/mean_merge9_97Ab.csv"))

dir.create(dirname(file.path(table_dir, "IFI30_distance/IFI30Ab_neg/neg_merge9_97Ab.csv")), recursive = TRUE, 
    showWarnings = FALSE)

write.csv(pos_merge, file.path(table_dir, "IFI30_distance/IFI30Ab_neg/neg_merge9_97Ab.csv"))

#### 汇总图 ====
mean_neg <- read.csv(file.path(table_dir, "IFI30_distance/IFI30Ab_neg/mean_merge9_97Ab.csv"))

df_mean_neg <- mean_neg %>% group_by(group) %>% summarise(mean_layer = mean(mean_exp)) %>% mutate(Ab = "Negative")

mean_pos <- read.csv(file.path(input_dir, "result/03.Ag_Ab/Ab周围IFI30表达量/IFI30Ab_pos/mean_merge9_97Ab.csv"))

df_mean_pos <- mean_pos %>% group_by(group) %>% summarise(mean_layer = mean(mean_exp)) %>% mutate(Ab = "Positive")

dat_p <- rbind(df_mean_pos, df_mean_neg) %>% filter(!group %in% c("9", "≥10"))

dat_p$group <- glue("{as.numeric(dat_p$group)*25}µm")

dat_p$group <- factor(dat_p$group, levels = c("25µm", "50µm", "75µm", "100µm", "125µm", "150µm", 
    "175µm", "200µm"))

p_mean <- ggplot(dat_p, aes_string("group", "mean_layer", color = "Ab", group = "Ab")) + geom_point(size = 2) + 
    geom_line(size = 1) + labs(title = "IFI30 expression", x = "", y = " ") + scale_color_manual(values = mycolor_new) + 
    theme_classic(base_line_size = 1) + theme(plot.title = element_text(size = 15, colour = "black", 
    hjust = 0.5), axis.title.y = element_text(size = 15, color = "black", face = "bold", vjust = 1.9, 
    hjust = 0.5, angle = 90), legend.title = element_text(color = "black", size = 15, face = "bold"), 
    legend.text = element_text(color = "black", size = 10, face = "bold"), axis.text.x = element_text(size = 13, 
        color = "black", face = "bold", vjust = 0, hjust = 1, angle = 90), axis.text.y = element_text(size = 13, 
        color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))

p_mean

pos_neg <- read.csv(file.path(table_dir, "IFI30_distance/IFI30Ab_neg/neg_merge9_97Ab.csv"))

df_pos_neg <- pos_neg %>% group_by(group) %>% summarise(mean_layer = mean(proportion)) %>% mutate(Ab = "Negative")

pos_pos <- read.csv(file.path(input_dir, "result/03.Ag_Ab/Ab周围IFI30表达量/IFI30Ab_pos/pos_merge9_97Ab.csv"))

df_pos_pos <- pos_pos %>% group_by(group) %>% summarise(mean_layer = mean(proportion)) %>% mutate(Ab = "Positive")

dat_p <- rbind(df_pos_pos, df_pos_neg) %>% filter(!group %in% c("9", "≥10"))

dat_p$group <- glue("{as.numeric(dat_p$group)*25}µm")

dat_p$group <- factor(dat_p$group, levels = c("25µm", "50µm", "75µm", "100µm", "125µm", "150µm", 
    "175µm", "200µm"))

p_pos <- ggplot(dat_p, aes_string("group", "mean_layer", color = "Ab", group = "Ab")) + geom_point(size = 2) + 
    geom_line(size = 1) + labs(title = "IFI30 density", x = "", y = " ") + scale_color_manual(values = mycolor_new) + 
    theme_classic(base_line_size = 1) + theme(plot.title = element_text(size = 15, colour = "black", 
    hjust = 0.5), axis.title.y = element_text(size = 15, color = "black", face = "bold", vjust = 1.9, 
    hjust = 0.5, angle = 90), legend.title = element_text(color = "black", size = 15, face = "bold"), 
    legend.text = element_text(color = "black", size = 10, face = "bold"), axis.text.x = element_text(size = 13, 
        color = "black", face = "bold", vjust = 0, hjust = 1, angle = 90), axis.text.y = element_text(size = 13, 
        color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))

p_pos

dir.create(dirname(file.path(figure_dir, "IFI30_distance/lineplot_distance_97")), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(wrap_plots(list(p_mean, p_pos)), file.path(figure_dir, "IFI30_distance/lineplot_distance_97"), 
    if_pdf = T, width = 12, height = 4)

