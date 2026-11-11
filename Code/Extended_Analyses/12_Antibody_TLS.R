#### 40+57抗体TLS内外bin比例 ====
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/R/io.R")
source("Code/Extended_Analyses/R/plot.R")
input_dir <- file.path(data_dir, "hzy")
processed_dir <- file.path(dirs$processed, "hzy")
table_dir <- file.path(dirs$tables, "hzy")
figure_dir <- file.path(dirs$figures, "hzy")
library(dplyr)
library(tidyr)
library(stringr)
library(glue)
library(ggplot2)
library(ggpubr)
library(scales)
library(readxl)
library(writexl)

bcr_smrt_merge <- readRDS(file.path(input_dir, "result/03.Ag_Ab/merged_bcr.rds"))

changeo_merge <- readRDS(file.path(input_dir, "result/03.Ag_Ab/merge4.rds"))

merge_metadata <- read.table(file.path(input_dir, "data/Bin50_classic_celltype_metadata_merge.txt"), 
    header = T)

#### 第一批抗体 ====
ad_info_all <- readxl::read_xlsx(file.path(input_dir, "result/03.Ag_Ab/97Ab.xlsx"))

ad_info <- readxl::read_xlsx(file.path(input_dir, "result/03.Ag_Ab/final_genescript2.xlsx"))

ad_info$sequence_id <- sub("-", "_", ad_info$raw_sc_id)

ad_info$name <- str_replace_all(ad_info$IG_name, "IG_", "IgG")

changeo_ab_125 <- changeo_merge[changeo_merge$sequence_id %in% ad_info$sequence_id, ]

changeo_ab_125 <- changeo_ab_125 %>% dplyr::left_join(ad_info[, c("sequence_id", "name")])

changeo_ab_125 <- changeo_ab_125[, c("name", "clone_name")]

changeo_ab_125$sample <- ""

for (cur_i in 1:nrow(changeo_ab_125)) {
    print(cur_i)
    cur_clone <- changeo_ab_125$clone_name[cur_i]
    cur_df <- changeo_merge[changeo_merge$clone_name == cur_clone & changeo_merge$type == "Spatial", 
        ]
    cur_seq <- cur_df$sequence_id
    cur_sample <- table(cur_df$sample) %>% sort() %>% rev()
    cur_sample <- paste(names(cur_sample), collapse = ",")
    changeo_ab_125$sample[cur_i] <- cur_sample
}

dir.create(dirname(file.path(table_dir, "antibody_TLS/40Ab/Ab_clonename.csv")), recursive = TRUE, 
    showWarnings = FALSE)

write.csv(changeo_ab_125, file.path(table_dir, "antibody_TLS/40Ab/Ab_clonename.csv"))

changeo_ab_125 <- read.csv(file.path(table_dir, "antibody_TLS/40Ab/Ab_clonename.csv"), header = T, 
    row.names = 1)

df <- changeo_ab_125 %>% filter(name %in% ad_info_all$IG_name)

dir.create(dirname(file.path(table_dir, "antibody_TLS/40Ab/Ab40_clonename.csv")), recursive = TRUE, 
    showWarnings = FALSE)

write.csv(df, file.path(table_dir, "antibody_TLS/40Ab/Ab40_clonename.csv"))

df_Ab_batch1 <- filter(df, sample != "") %>% tidyr::separate_rows(sample, sep = ",")

df_Ab_batch1$clone_name2 <- paste(df_Ab_batch1$sample, df_Ab_batch1$clone_name, sep = "_")

changeo_merge$clone_name2 <- paste(changeo_merge$sample, changeo_merge$clone_name, sep = "_")

intersect(df_Ab_batch1$clone_name2, changeo_merge$clone_name2) %>% length()

df_Ab49 <- filter(changeo_merge, clone_name2 %in% df_Ab_batch1$clone_name2 & type == "Spatial")

pos_matrix <- matrix(nrow = length(unique(df_Ab49$clone_name2)), ncol = 5, dimnames = list(unique(df_Ab49$clone_name2), 
    c("clone_name", "sample", "Ig_num", "TLSin", "TLS_ratio")))

cur_Ab <- unique(df_Ab49$clone_name2)[2]

for (cur_Ab in unique(df_Ab49$clone_name2)) {
    print(cur_Ab)
    cur_sample <- df_Ab49$sample[df_Ab49$clone_name2 == cur_Ab] %>% unique()
    cur_Abname <- df_Ab_batch1$name[df_Ab_batch1$clone_name2 == cur_Ab]
    df_Ab <- filter(bcr_smrt_merge, sequence_id %in% df_Ab49$sequence_id[df_Ab49$clone_name2 == 
        cur_Ab])
    print(table(df_Ab$locus))
    all_bin <- df_Ab %>% pull(Bin50)
    all_bin <- paste(cur_sample, all_bin, sep = "_")
    df <- merge_metadata %>% filter(SampleID == cur_sample)
    df_bin <- filter(df, BinID %in% all_bin)
    df_bin <- df_bin %>% dplyr::mutate(TLS_group = ifelse(grepl("NA", TLS_raw), "OutTLS", "InTLS"))
    fill_colors <- c("#dd3497", "#a1d99b")
    names(fill_colors) <- c("InTLS", "OutTLS")
    bcr_pie_df <- df_bin %>% group_by(TLS_group) %>% summarise(value = n(), .groups = "drop") %>% 
        mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(TLS_group, 
            " (", percent, ")"))
    TLSin_bin = nrow(df_bin[df_bin$TLS_group == "InTLS", ])
    pos_matrix[cur_Ab, ] <- c(cur_Abname, cur_sample, nrow(df_bin), TLSin_bin, TLSin_bin/nrow(df_bin))
}

dir.create(dirname(file.path(table_dir, "antibody_TLS/40Ab/TLS_ratio.csv")), recursive = TRUE, showWarnings = FALSE)

write.csv(pos_matrix, file.path(table_dir, "antibody_TLS/40Ab/TLS_ratio.csv"))

bcr_smrt_merge <- readRDS(file.path(input_dir, "data/merged_bcr.rds"))

changeo_merge <- readRDS(file.path(input_dir, "data/merge4.rds"))

#### 新增抗体 ====
new_Ab <- readxl::read_xlsx(file.path(input_dir, "result/03.Ag_Ab/抗体信息_new.xlsx"))

changeo_merge$clone_name2 <- paste(changeo_merge$sample, changeo_merge$clone_name, sep = "_")

intersect(new_Ab$clone_name, changeo_merge$clone_name2) %>% length()

df_Ab57 <- filter(changeo_merge, clone_name2 %in% new_Ab$clone_name & type == "Spatial")

pos_matrix <- matrix(nrow = length(unique(df_Ab57$clone_name2)), ncol = 5, dimnames = list(unique(df_Ab57$clone_name2), 
    c("clone_name", "sample", "Ig_num", "TLSin", "TLS_ratio")))

cur_Ab <- unique(df_Ab57$clone_name2)[1]

for (cur_Ab in unique(df_Ab57$clone_name2)) {
    print(cur_Ab)
    cur_sample <- df_Ab57$sample[df_Ab57$clone_name2 == cur_Ab] %>% unique()
    cur_Abname <- new_Ab$name[new_Ab$clone_name == cur_Ab]
    df_Ab <- filter(bcr_smrt_merge, sequence_id %in% df_Ab57$sequence_id[df_Ab57$clone_name2 == 
        cur_Ab])
    print(table(df_Ab$locus))
    all_bin <- df_Ab %>% pull(Bin50)
    all_bin <- paste(cur_sample, all_bin, sep = "_")
    df <- merge_metadata %>% filter(SampleID == cur_sample)
    df_bin <- filter(df, BinID %in% all_bin)
    df_bin <- df_bin %>% dplyr::mutate(TLS_group = ifelse(grepl("NA", TLS_raw), "OutTLS", "InTLS"))
    fill_colors <- c("#dd3497", "#a1d99b")
    names(fill_colors) <- c("InTLS", "OutTLS")
    bcr_pie_df <- df_bin %>% group_by(TLS_group) %>% summarise(value = n(), .groups = "drop") %>% 
        mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(TLS_group, 
            " (", percent, ")"))
    TLSin_bin = nrow(df_bin[df_bin$TLS_group == "InTLS", ])
    pos_matrix[cur_Ab, ] <- c(cur_Abname, cur_sample, nrow(df_bin), TLSin_bin, TLSin_bin/nrow(df_bin))
}

dir.create(dirname(file.path(table_dir, "antibody_TLS/57Ab/TLS_ratio.csv")), recursive = TRUE, showWarnings = FALSE)

write.csv(pos_matrix, file.path(table_dir, "antibody_TLS/57Ab/TLS_ratio.csv"))

#### 合并比例 ====
Ab_97 <- readxl::read_xlsx(file.path(input_dir, "result/03.Ag_Ab/97Ab.xlsx"))

TLSratio_batch1 <- read.csv(file.path(table_dir, "antibody_TLS/40Ab/TLS_ratio.csv"), header = T, 
    row.names = 1)

TLSratio_batch1 <- TLSratio_batch1 %>% group_by(clone_name) %>% summarise(mean_ratio = mean(TLS_ratio))

TLSratio_batch2 <- read.csv(file.path(table_dir, "antibody_TLS/57Ab/TLS_ratio.csv"), header = T)

TLSratio_batch2 <- TLSratio_batch2 %>% dplyr::select(-c(clone_name)) %>% left_join(Ab_97, by = c(X = "clone_name")) %>% 
    group_by(IG_name) %>% summarise(mean_ratio = mean(TLS_ratio))

colnames(TLSratio_batch2)[1] <- "clone_name"

merge_ratio <- rbind(TLSratio_batch1, TLSratio_batch2)

Ab_97_merge <- Ab_97 %>% left_join(merge_ratio, by = c(IG_name = "clone_name"))

Ab_97_merge$mean_ratio <- ifelse(is.na(Ab_97_merge$mean_ratio), 0, Ab_97_merge$mean_ratio)

colnames(Ab_97_merge)[c(12, 14)] <- c("IFI30_react", "Batch")

Ab_97_merge$IFI30_react <- ifelse(Ab_97_merge$IFI30_react == "×", "No", "Yes")

Ab_97_merge$group <- Ab_97_merge$IFI30_react

dir.create(dirname(file.path(table_dir, "antibody_TLS/97Ab_TLSratio.xlsx")), recursive = TRUE, showWarnings = FALSE)

writexl::write_xlsx(Ab_97_merge, file.path(table_dir, "antibody_TLS/97Ab_TLSratio.xlsx"))

df <- readxl::read_xlsx(file.path(input_dir, "data/97Ab_TLSratio.xlsx"))

fill_colors <- c("#dd3497", "#a1d99b")

names(fill_colors) <- c("InTLS", "OutTLS")

df_plot <- df

df_plot <- filter(df, IFI30_reaction == "Yes")

df_plot <- filter(df, `TLS+Ab` == "+")

TLS_in_ratio <- mean(df_plot$mean_ratio) %>% round(digits = 3)

TLS_out_ratio <- 1 - TLS_in_ratio

TLS_pie_df <- data.frame(ratio = c(TLS_in_ratio * 100, TLS_out_ratio * 100), TLS_group = c("InTLS", 
    "OutTLS"))

bcr_pie <- ggpie(TLS_pie_df, "ratio", label = TLS_pie_df$ratio, lab.pos = "in", lab.font = "white", 
    fill = "TLS_group", color = "white", palette = fill_colors) + theme(plot.caption = element_text(hjust = 0.5, 
    size = 12, face = "bold"))

bcr_pie

dir.create(dirname(glue(file.path(figure_dir, "antibody_TLS/35Ab_IFI30反应性抗体"))), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(bcr_pie, glue(file.path(figure_dir, "antibody_TLS/35Ab_IFI30反应性抗体")), height = 5, 
    width = 4, if_pdf = T)

