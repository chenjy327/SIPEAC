#### 跨组织亚型及扩增（NL、CLD、HCC、BC、CRC） ====
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
library(ggforce)
library(rlang)

my_countClones <- function(data, clone_col = "clone_name") {
    result <- data %>% dplyr::group_by(!!sym(clone_col)) %>% dplyr::summarise(seq_count = n(), .groups = "drop")
    return(result)
}

#### 正常肝 ====
NL_sample <- c("CRC48_LP", "CRC65_LP", "CRC66_LP", "CRC67_LP", "CRC_69_LP", "CRC70_LP")

NL_merged <- read.csv(file.path(input_dir, "result/05.BCR_TCR/shm.res_NL.csv"), header = T, check.names = F, 
    row.names = 1)

pie_df <- NL_merged %>% filter(c_call %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", 
    "IGHG3", "IGHG4", "IGHM") & orig.ident %in% NL_sample) %>% dplyr::group_by(c_call) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", pie_df$value, "\n", " (", pie_df$percent, ")")

total_num <- sum(pie_df$value)

pie_df$c_call <- factor(pie_df$c_call, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "c_call", color = NA, 
    palette = fill_colors) + labs(title = "NL")

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_NL_IGH"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_NL_IGH")), if_pdf = T, width = 4, height = 4)

df_iso <- NL_merged %>% filter(c_call %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", 
    "IGHG3", "IGHG4", "IGHM") & orig.ident %in% NL_sample)

result <- df_iso %>% group_by(c_call) %>% summarise(total_umis = sum(umi_count, na.rm = TRUE)) %>% 
    mutate(proportion = total_umis/sum(total_umis))

result$percent <- scales::percent(result$proportion, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", result$total_umis, "\n", " (", result$percent, ")")

total_num <- sum(result$total_umis)

result$c_call <- factor(result$c_call, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

p <- ggpie(result, "total_umis", label = labs, lab.pos = "in", lab.font = "black", fill = "c_call", 
    color = NA, palette = fill_colors) + labs(title = "NL")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_NL_IGH_umi"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_NL_IGH_umi")), if_pdf = T, width = 4, height = 4)

df_IGH <- filter(NL_merged, locus == "IGH" & orig.ident %in% NL_sample)

df_IGH$clone_name <- paste(df_IGH$orig.ident, df_IGH$clone_id, sep = "_")

bcr_smrt_clones <- my_countClones(df_IGH, clone = "clone_name")

df_pie <- bcr_smrt_clones %>% dplyr::mutate(expand_type = case_when(seq_count == 1 ~ "singleton", 
    seq_count >= 10 ~ "clone size ≥ 10", T ~ "clone size < 10"))

IGH_pie_df <- df_pie %>% group_by(expand_type) %>% summarise(value = n(), .groups = "drop") %>% 
    mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(expand_type, 
        "\n", " (", percent, ")"))

IGH_pie_df$expand_type <- factor(IGH_pie_df$expand_type, levels = c("singleton", "clone size < 10", 
    "clone size ≥ 10"))

fill_colors <- c("#65a3d1", "#edb55a", "#946625")

names(fill_colors) <- c("singleton", "clone size < 10", "clone size ≥ 10")

total_clones <- sum(IGH_pie_df$value)

p <- ggplot(data = IGH_pie_df) + ggforce::geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, 
    r = 2, amount = ratio, fill = expand_type)) + geom_text(aes(x = 1.5 * sin(2 * pi * (cumsum(ratio) - 
    ratio/2)/sum(ratio)), y = 1.5 * cos(2 * pi * (cumsum(ratio) - ratio/2)/sum(ratio)), label = label), 
    size = 4, color = "black") + annotate("text", x = 0, y = 0, label = glue("total = {total_clones}"), 
    size = 8, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors) + theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) + labs(title = "NL", 
    x = "", y = "")

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_NL_expand_IGH"))), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_NL_expand_IGH")), if_pdf = T, width = 6, 
    height = 6)

BCR_heavy <- readRDS(file.path(input_dir, "result/06.SMRTseq/merge_BCR_heavy.rds"))

df_list <- list()

for (cur_sample in c("NL_1", "NL_2", "NL_3")) {
    df_barcode <- read.csv(glue(file.path(input_dir, "spatial_bin/Convertbin/{cur_sample}_dedup.info.BinID.50.csv")), 
        sep = "\t")
    df_list[[cur_sample]] <- BCR_heavy %>% filter(sampleid == cur_sample & barcode %in% df_barcode$BC)
}

df_merge <- bind_rows(df_list)

df_merge <- df_merge %>% mutate(isotype = str_split_i(c_call, "\\*", 1))

df_merge$isotype[df_merge$isotype == "IGHG4A"] <- "IGHG4"

df_merge <- filter(df_merge, isotype %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

table(df_merge$isotype)

pie_data <- df_merge %>% count(isotype) %>% mutate(Percentage = n/sum(n) * 100, label = paste0(isotype, 
    "\n", round(Percentage, 1), "%"))

pie_df <- df_merge %>% dplyr::group_by(isotype) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)

pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "isotype", 
    color = NULL, palette = fill_colors) + labs(title = "NL")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/IgH_NL"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/IgH_NL")), if_pdf = T, width = 4, height = 4)

df_list <- list()

for (cur_sample in c("NL_1", "NL_2", "NL_3")) {
    print(cur_sample)
    df_barcode <- read.csv(glue(file.path(input_dir, "spatial_bin/Convertbin/{cur_sample}_dedup.info.BinID.50.csv")), 
        sep = "\t")
    df_list[[cur_sample]] <- BCR_heavy %>% filter(sampleid == cur_sample & barcode %in% df_barcode$BC)
}

df_merge <- bind_rows(df_list)

df_merge <- df_merge %>% mutate(isotype = str_split_i(c_call, "\\*", 1))

df_merge$isotype[df_merge$isotype == "IGHG4A"] <- "IGHG4"

df_merge <- filter(df_merge, isotype %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

table(df_merge$isotype)

df_merge$clone_name <- paste(df_merge$sampleid, df_merge$clone_id, sep = "_")

bcr_smrt_clones <- my_countClones(df_merge, clone = "clone_name")

df_pie <- bcr_smrt_clones %>% dplyr::mutate(expand_type = case_when(seq_count == 1 ~ "singleton", 
    seq_count >= 10 ~ "clone size ≥ 10", T ~ "clone size < 10"))

IGH_pie_df <- df_pie %>% group_by(expand_type) %>% summarise(value = n(), .groups = "drop") %>% 
    mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(expand_type, 
        "\n", " (", percent, ")"))

IGH_pie_df$expand_type <- factor(IGH_pie_df$expand_type, levels = c("singleton", "clone size < 10", 
    "clone size ≥ 10"))

fill_colors <- c("#65a3d1", "#edb55a", "#946625")

names(fill_colors) <- c("singleton", "clone size < 10", "clone size ≥ 10")

total_clones <- sum(IGH_pie_df$value)

p <- ggplot(data = IGH_pie_df) + ggforce::geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, 
    r = 2, amount = ratio, fill = expand_type)) + geom_text(aes(x = 1.5 * sin(2 * pi * (cumsum(ratio) - 
    ratio/2)/sum(ratio)), y = 1.5 * cos(2 * pi * (cumsum(ratio) - ratio/2)/sum(ratio)), label = label), 
    size = 4, color = "black") + annotate("text", x = 0, y = 0, label = glue("total = {total_clones}"), 
    size = 8, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors) + theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) + labs(title = "NL", 
    x = "", y = "")

dir.create(dirname(glue(file.path(figure_dir, "isotype/Expand_NL"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/Expand_NL")), if_pdf = T, width = 6, height = 6)

#### 肝硬化 ====
CLD_sample <- c("CLD10", "CLD11", "CLD12", "CLD13", "CLD14", "CLD2", "CLD7")

CLD_sample <- c("CLD10_Liver", "CLD11_liver", "CLD12_liver", "CLD13_liver", "CLD14_liver", "CLD2", 
    "CLD7")

CLD_merged <- read.csv(file.path(input_dir, "result/05.BCR_TCR/shm.res_CLD.csv"), header = T, check.names = F, 
    row.names = 1)

pie_df <- CLD_merged %>% filter(c_call %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", 
    "IGHG3", "IGHG4", "IGHM") & orig.ident %in% CLD_sample) %>% dplyr::group_by(c_call) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", pie_df$value, "\n", " (", pie_df$percent, ")")

total_num <- sum(pie_df$value)

pie_df$c_call <- factor(pie_df$c_call, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "c_call", color = NA, 
    palette = fill_colors) + labs(title = "CLD")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_CLD_IGH"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_CLD_IGH")), if_pdf = T, width = 4, height = 4)

df_iso <- CLD_merged %>% filter(c_call %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", 
    "IGHG3", "IGHG4", "IGHM") & orig.ident %in% CLD_sample)

result <- df_iso %>% group_by(c_call) %>% summarise(total_umis = sum(umi_count, na.rm = TRUE)) %>% 
    mutate(proportion = total_umis/sum(total_umis))

result$percent <- scales::percent(result$proportion, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", result$total_umis, "\n", " (", result$percent, ")")

total_num <- sum(result$total_umis)

result$c_call <- factor(result$c_call, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

p <- ggpie(result, "total_umis", label = labs, lab.pos = "in", lab.font = "black", fill = "c_call", 
    color = NA, palette = fill_colors) + labs(title = "CLD")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_CLD_IGH_umi"))), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_CLD_IGH_umi")), if_pdf = T, width = 4, height = 4)

df_IGH <- filter(CLD_merged, locus == "IGH" & orig.ident %in% CLD_sample)

df_IGH$clone_name <- paste(df_IGH$orig.ident, df_IGH$clone_id, sep = "_")

bcr_smrt_clones <- my_countClones(df_IGH, clone = "clone_name")

df_pie <- bcr_smrt_clones %>% dplyr::mutate(expand_type = case_when(seq_count == 1 ~ "singleton", 
    seq_count >= 10 ~ "clone size ≥ 10", T ~ "clone size < 10"))

IGH_pie_df <- df_pie %>% group_by(expand_type) %>% summarise(value = n(), .groups = "drop") %>% 
    mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(expand_type, 
        "\n", " (", percent, ")"))

IGH_pie_df$expand_type <- factor(IGH_pie_df$expand_type, levels = c("singleton", "clone size < 10", 
    "clone size ≥ 10"))

fill_colors <- c("#65a3d1", "#edb55a", "#946625")

names(fill_colors) <- c("singleton", "clone size < 10", "clone size ≥ 10")

total_clones <- sum(IGH_pie_df$value)

p <- ggplot(data = IGH_pie_df) + ggforce::geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, 
    r = 2, amount = ratio, fill = expand_type)) + geom_text(aes(x = 1.5 * sin(2 * pi * (cumsum(ratio) - 
    ratio/2)/sum(ratio)), y = 1.5 * cos(2 * pi * (cumsum(ratio) - ratio/2)/sum(ratio)), label = label), 
    size = 4, color = "black") + annotate("text", x = 0, y = 0, label = glue("total = {total_clones}"), 
    size = 8, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors) + theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) + labs(title = "CLD", 
    x = "", y = "")

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_CLD_expand_IGH"))), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_CLD_expand_IGH")), if_pdf = T, width = 6, 
    height = 6)

BCR_heavy <- readRDS(file.path(input_dir, "result/06.SMRTseq/merge_BCR_heavy.rds"))

df_list <- list()

for (cur_sample in c("CLD11", "CLD13", "CLD14")) {
    df_barcode <- read.csv(glue(file.path(input_dir, "spatial_bin/Convertbin/{cur_sample}_dedup.info.BinID.50.csv")), 
        sep = "\t")
    df_list[[cur_sample]] <- BCR_heavy %>% filter(sampleid == cur_sample & barcode %in% df_barcode$BC)
}

df_merge <- bind_rows(df_list)

df_merge <- df_merge %>% mutate(isotype = str_split_i(c_call, "\\*", 1))

df_merge$isotype[df_merge$isotype == "IGHG4A"] <- "IGHG4"

df_merge <- filter(df_merge, isotype %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

table(df_merge$isotype)

pie_data <- df_merge %>% count(isotype) %>% mutate(Percentage = n/sum(n) * 100, label = paste0(isotype, 
    "\n", round(Percentage, 1), "%"))

pie_df <- df_merge %>% dplyr::group_by(isotype) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)

pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "isotype", 
    color = NULL, palette = fill_colors) + labs(title = "CLD")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/IgH_CLD"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/IgH_CLD")), if_pdf = T, width = 4, height = 4)

df_list <- list()

for (cur_sample in c("CLD11", "CLD13", "CLD14")) {
    print(cur_sample)
    df_barcode <- read.csv(glue(file.path(input_dir, "spatial_bin/Convertbin/{cur_sample}_dedup.info.BinID.50.csv")), 
        sep = "\t")
    df_list[[cur_sample]] <- BCR_heavy %>% filter(sampleid == cur_sample & barcode %in% df_barcode$BC)
}

df_merge <- bind_rows(df_list)

df_merge <- df_merge %>% mutate(isotype = str_split_i(c_call, "\\*", 1))

df_merge$isotype[df_merge$isotype == "IGHG4A"] <- "IGHG4"

df_merge <- filter(df_merge, isotype %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

table(df_merge$isotype)

df_merge$clone_name <- paste(df_merge$sampleid, df_merge$clone_id, sep = "_")

bcr_smrt_clones <- my_countClones(df_merge, clone = "clone_name")

df_pie <- bcr_smrt_clones %>% dplyr::mutate(expand_type = case_when(seq_count == 1 ~ "singleton", 
    seq_count >= 10 ~ "clone size ≥ 10", T ~ "clone size < 10"))

IGH_pie_df <- df_pie %>% group_by(expand_type) %>% summarise(value = n(), .groups = "drop") %>% 
    mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(expand_type, 
        "\n", " (", percent, ")"))

IGH_pie_df$expand_type <- factor(IGH_pie_df$expand_type, levels = c("singleton", "clone size < 10", 
    "clone size ≥ 10"))

fill_colors <- c("#65a3d1", "#edb55a", "#946625")

names(fill_colors) <- c("singleton", "clone size < 10", "clone size ≥ 10")

total_clones <- sum(IGH_pie_df$value)

p <- ggplot(data = IGH_pie_df) + ggforce::geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, 
    r = 2, amount = ratio, fill = expand_type)) + geom_text(aes(x = 1.5 * sin(2 * pi * (cumsum(ratio) - 
    ratio/2)/sum(ratio)), y = 1.5 * cos(2 * pi * (cumsum(ratio) - ratio/2)/sum(ratio)), label = label), 
    size = 4, color = "black") + annotate("text", x = 0, y = 0, label = glue("total = {total_clones}"), 
    size = 8, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors) + theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) + labs(title = "CLD", 
    x = "", y = "")

dir.create(dirname(glue(file.path(figure_dir, "isotype/Expand_CLD"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/Expand_CLD")), if_pdf = T, width = 6, height = 6)

#### 肝细胞癌 ====
samplelist <- list.files(file.path(input_dir, "data/BCR_data/filtered_contig_annotations"))

df_list <- list()

for (cur_sample in samplelist) {
    print(cur_sample)
    df_BCR <- read.csv(glue(file.path(input_dir, "data/BCR_data/filtered_contig_annotations/{cur_sample}/filtered_contig_annotations.csv")), 
        header = T)
    df_BCR <- filter(df_BCR, chain == "IGH") %>% mutate(sample = cur_sample)
    df_list[[cur_sample]] <- df_BCR
}

df_BCR_merge1 <- bind_rows(df_list)

dir.create(dirname(file.path(table_dir, "isotype/HCC_BCR_heavy_scRNA_未去重.csv")), recursive = TRUE, 
    showWarnings = FALSE)

write.csv(df_BCR_merge1, file.path(table_dir, "isotype/HCC_BCR_heavy_scRNA_未去重.csv"))

mrna_list <- readRDS(file.path(input_dir, "data/mrna_tcr_66sample.rds"))

df_BCR_filter1 <- filter(df_BCR_merge1, sample %in% mrna_list) %>% filter(!sample %in% c("2931T_B_FS", 
    "2976T_B_FS_2", "2976T_B_FS_3", "2976T_B_FS_1", "3125_P", "3163_T")) %>% filter(c_gene != "") %>% 
    mutate(sequence_id = paste(sample, contig_id, sep = "_"))

sc_bcr <- readRDS(file.path(input_dir, "data/merged_bcr_sc.rds"))

cur_sc_bcr <- sc_bcr %>% filter(sample %in% mrna_list) %>% filter(!sample %in% c("2931T_B_FS", "2976T_B_FS_2", 
    "2976T_B_FS_3", "2976T_B_FS_1", "3125_P", "3163_T")) %>% filter(locus == "IGH")

cur_sc_bcr$isotype <- substr(cur_sc_bcr$c_call, 1, 5)

cur_sc_bcr$isotype %>% unique()

cur_sc_bcr <- subset(cur_sc_bcr, isotype != "")

cur_sc_bcr[cur_sc_bcr$isotype == "IGHM*", ]$isotype <- "IGHM"

cur_sc_bcr[cur_sc_bcr$isotype == "IGHD*", ]$isotype <- "IGHD"

cur_sc_bcr[cur_sc_bcr$isotype == "IGHE*", ]$isotype <- "IGHE"

nrow(cur_sc_bcr)

cur_meta <- left_join(cur_sc_bcr, df_BCR_filter1[, c("sequence_id", "umis")], by = "sequence_id") %>% 
    filter(!is.na(umis))

pie_df <- cur_meta %>% dplyr::group_by(isotype) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", pie_df$value, "\n", " (", pie_df$percent, ")")

total_num <- sum(pie_df$value)

pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "isotype", 
    color = NA, palette = fill_colors) + labs(title = "HCC")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_HCC_IGH"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_HCC_IGH")), if_pdf = T, width = 4, height = 4)

cur_meta <- left_join(cur_sc_bcr, df_BCR_filter1[, c("sequence_id", "umis")], by = "sequence_id") %>% 
    filter(!is.na(umis))

df_iso <- cur_meta %>% filter(isotype %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", 
    "IGHG3", "IGHG4", "IGHM"))

result <- df_iso %>% group_by(isotype) %>% summarise(total_umis = sum(umis, na.rm = TRUE)) %>% mutate(proportion = total_umis/sum(total_umis))

result$percent <- scales::percent(result$proportion, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", result$total_umis, "\n", " (", result$percent, ")")

total_num <- sum(result$total_umis)

result$isotype <- factor(result$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

p <- ggpie(result, "total_umis", label = labs, lab.pos = "in", lab.font = "black", fill = "isotype", 
    color = NA, palette = fill_colors) + labs(title = "HCC")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_HCC_IGH_umi"))), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_HCC_IGH_umi")), if_pdf = T, width = 4, height = 4)

mrna_list <- readRDS(file.path(input_dir, "data/mrna_tcr_66sample.rds"))

sc_bcr <- readRDS(file.path(input_dir, "data/merged_bcr_sc.rds"))

cur_sc_bcr <- sc_bcr %>% filter(sample %in% mrna_list) %>% filter(!sample %in% c("2931T_B_FS", "2976T_B_FS_2", 
    "2976T_B_FS_3", "2976T_B_FS_1", "3125_P", "3163_T")) %>% filter(locus == "IGH")

cur_sc_bcr$isotype <- substr(cur_sc_bcr$c_call, 1, 5)

cur_sc_bcr$isotype %>% unique()

cur_sc_bcr <- subset(cur_sc_bcr, isotype != "")

cur_sc_bcr[cur_sc_bcr$isotype == "IGHM*", ]$isotype <- "IGHM"

cur_sc_bcr[cur_sc_bcr$isotype == "IGHD*", ]$isotype <- "IGHD"

cur_sc_bcr[cur_sc_bcr$isotype == "IGHE*", ]$isotype <- "IGHE"

cur_sc_bcr <- distinct(cur_sc_bcr, sequence_id, .keep_all = TRUE)

cur_sc_bcr$clone_id <- paste(cur_sc_bcr$sample, cur_sc_bcr$clone_id, sep = "_")

bcr_clones <- my_countClones(cur_sc_bcr, clone = "clone_id")

bcr_igh <- left_join(cur_sc_bcr, bcr_clones)

bcr_bin_count <- bcr_igh %>% group_by(clone_id) %>% summarise(bin_count = n(), .groups = "drop")

bcr_clones <- left_join(bcr_clones, bcr_bin_count)

bcr_clones <- bcr_clones %>% dplyr::mutate(expand_type = dplyr::case_when(bin_count >= 10 ~ "Expanded >= 10", 
    bin_count == 1 ~ "Singleton", bin_count > 1 & bin_count < 10 ~ "Expanded < 10"))

bcr_pie_df <- bcr_clones %>% group_by(expand_type) %>% summarise(value = n(), .groups = "drop") %>% 
    mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(expand_type, 
        " (", percent, ")"))

bcr_pie_df$expand_type <- factor(bcr_pie_df$expand_type, levels = c("Singleton", "Expanded < 10", 
    "Expanded >= 10"))

fill_colors <- c("#65a3d1", "#edb55a", "#946625")

names(fill_colors) <- c("Singleton", "Expanded < 10", "Expanded >= 10")

total_clones <- sum(bcr_pie_df$value)

p <- ggplot(data = bcr_pie_df) + ggforce::geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, 
    r = 2, amount = ratio, fill = expand_type)) + geom_text(aes(x = 1.5 * sin(2 * pi * (cumsum(ratio) - 
    ratio/2)/sum(ratio)), y = 1.5 * cos(2 * pi * (cumsum(ratio) - ratio/2)/sum(ratio)), label = label), 
    size = 4, color = "black") + annotate("text", x = 0, y = 0, label = glue("total = {total_clones}"), 
    size = 8, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors) + theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) + labs(title = "HCC", 
    x = "", y = "")

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_HCC_expand_IGH"))), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_HCC_expand_IGH")), if_pdf = T, width = 6, 
    height = 6)

bcr_merge <- readRDS(file.path(input_dir, "data/merged_bcr.rds"))

bcr_merge <- bcr_merge[bcr_merge$locus == "IGH", ]

pie_df <- bcr_merge %>% dplyr::group_by(isotype) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "isotype", 
    color = "white", palette = fill_colors)

p

dir.create(dirname(file.path(figure_dir, "isotype/C-IgH")), recursive = TRUE, showWarnings = FALSE)

save_plot(p, file.path(figure_dir, "isotype/C-IgH"), if_pdf = T, width = 4, height = 4)

bcr_data <- readRDS(file.path(input_dir, "data/merged_bcr.rds"))

bcr_data$clone_id <- paste(bcr_data$sample, bcr_data$clone_id, sep = "_")

bcr_igh <- bcr_data[bcr_data$locus == "IGH", ]

bcr_clones <- my_countClones(bcr_igh, clone = "clone_id")

bcr_igh <- left_join(bcr_igh, bcr_clones)

bcr_bin_count <- bcr_igh %>% distinct(clone_id, bin_id) %>% group_by(clone_id) %>% summarise(bin_count = n(), 
    .groups = "drop")

bcr_clones <- left_join(bcr_clones, bcr_bin_count)

bcr_clones <- bcr_clones %>% dplyr::mutate(expand_type = dplyr::case_when(bin_count >= 10 ~ "Expanded >= 10", 
    bin_count == 1 ~ "Singleton", bin_count > 1 & bin_count < 10 ~ "Expanded < 10"))

bcr_pie_df <- bcr_clones %>% group_by(expand_type) %>% summarise(value = n(), .groups = "drop") %>% 
    mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(expand_type, 
        " (", percent, ")"))

fill_colors <- c("#946625", "#edb55a", "#65a3d1")

names(fill_colors) <- c("Expanded >= 10", "Expanded < 10", "Singleton")

total_clones <- sum(bcr_pie_df$value)

p <- ggplot(data = bcr_pie_df) + ggforce::geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, 
    r = 2, amount = ratio, fill = expand_type)) + geom_text(aes(x = 1.5 * sin(2 * pi * (cumsum(ratio) - 
    ratio/2)/sum(ratio)), y = 1.5 * cos(2 * pi * (cumsum(ratio) - ratio/2)/sum(ratio)), label = label), 
    size = 4, color = "black") + annotate("text", x = 0, y = 0, label = glue("total = {total_clones}"), 
    size = 8, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors) + theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) + labs(title = "HCC", 
    x = "", y = "")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/Expand_HCC"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/Expand_HCC")), if_pdf = T, width = 6, height = 6)

#### 乳腺癌 ====
BC_IgH <- read.csv(file.path(input_dir, "data/BC/20220309_BCSA23_Bcell_filtered_contig_annotations.csv"))

cur_meta <- subset(BC_IgH, chain == "IGH")

cur_meta <- subset(cur_meta, c_gene != "")

cur_meta$isotype <- cur_meta$c_gene

pie_df <- cur_meta %>% filter(isotype %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", 
    "IGHG3", "IGHG4", "IGHM")) %>% dplyr::group_by(isotype) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", pie_df$value, "\n", " (", pie_df$percent, ")")

total_num <- sum(pie_df$value)

pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "isotype", 
    color = NA, palette = fill_colors) + labs(title = "BC")

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_BC_IGH"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_BC_IGH")), if_pdf = T, width = 4, height = 4)

df_iso <- cur_meta %>% filter(isotype %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", 
    "IGHG3", "IGHG4", "IGHM"))

result <- df_iso %>% group_by(isotype) %>% summarise(total_umis = sum(umis, na.rm = TRUE)) %>% mutate(proportion = total_umis/sum(total_umis))

result$percent <- scales::percent(result$proportion, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", result$total_umis, "\n", " (", result$percent, ")")

total_num <- sum(result$total_umis)

result$isotype <- factor(result$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

p <- ggpie(result, "total_umis", label = labs, lab.pos = "in", lab.font = "black", fill = "isotype", 
    color = NA, palette = fill_colors) + labs(title = "BC")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_BC_IGH_umi"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_BC_IGH_umi")), if_pdf = T, width = 4, height = 4)

P1 <- fread(file.path(input_dir, "data/BC/P1_all_ontissue_clone_list.csv"))

P1$cloneId <- paste0("P1_", P1$cloneId)

P1 <- P1[, c("cloneId", "allCHitsWithScore", "cloneCount")]

P2 <- fread(file.path(input_dir, "data/BC/P2_all_ontissue_clone_list.csv"))

P2$cloneId <- paste0("P2_", P2$cloneId)

P2 <- P2[, c("cloneId", "allCHitsWithScore", "cloneCount")]

merge_vdj <- rbind(P1, P2)

merge_vdj$isotype <- substr(merge_vdj$allCHitsWithScore, 1, 5)

table(merge_vdj$isotype)

cur_meta <- subset(merge_vdj, isotype %in% c("IGHA1", "IGHA2", "IGHD*", "IGHG1", "IGHG2", "IGHG3", 
    "IGHG4", "IGHM*"))

pie_df <- cur_meta %>% dplyr::group_by(isotype) %>% dplyr::summarise(value = sum(cloneCount))

pie_df[pie_df$isotype == "IGHD*", ]$isotype <- "IGHD"

pie_df[pie_df$isotype == "IGHM*", ]$isotype <- "IGHM"

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "white", fill = "isotype", 
    color = NULL, palette = fill_colors)

p

#### 结直肠癌 ====
CRC_sample <- c("4272_PT", "3992_PT", "4273_PT", "3969_PT", "4265_PT", "4288_PT", "4043_PT", "4194_PT")

CRC_merged <- read.csv(file.path(input_dir, "result/05.BCR_TCR/CRC_BCR.csv"), header = T, check.names = F, 
    row.names = 1)

CRC_IgH <- subset(CRC_merged, chain %in% c("IGH"))

pie_df <- CRC_IgH %>% filter(c_gene %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3", 
    "IGHG4", "IGHM")) %>% dplyr::group_by(c_gene) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- scales::percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", pie_df$value, "\n", " (", pie_df$percent, ")")

total_num <- sum(pie_df$value)

pie_df$c_gene <- factor(pie_df$c_gene, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "c_gene", color = NA, 
    palette = fill_colors) + labs(title = "CRC")

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_CRC_IGH"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_CRC_IGH")), if_pdf = T, width = 4, height = 4)

df_iso <- CRC_IgH %>% filter(c_gene %in% c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3", 
    "IGHG4", "IGHM"))

result <- df_iso %>% group_by(c_gene) %>% summarise(total_umis = sum(umis, na.rm = TRUE)) %>% mutate(proportion = total_umis/sum(total_umis))

result$percent <- scales::percent(result$proportion, accuracy = 0.1)

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0("n=", result$total_umis, "\n", " (", result$percent, ")")

total_num <- sum(result$total_umis)

result$c_gene <- factor(result$c_gene, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

p <- ggpie(result, "total_umis", label = labs, lab.pos = "in", lab.font = "black", fill = "c_gene", 
    color = NA, palette = fill_colors) + labs(title = "CRC")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_CRC_IGH_umi"))), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_CRC_IGH_umi")), if_pdf = T, width = 4, height = 4)

CRC_IgH <- subset(CRC_merged, chain %in% c("IGH"))

CRC_IgH$clone_name <- paste(CRC_IgH$orig.ident, CRC_IgH$raw_clonotype_id, sep = "_")

bcr_smrt_clones <- my_countClones(CRC_IgH, clone = "clone_name")

df_pie <- bcr_smrt_clones %>% dplyr::mutate(expand_type = case_when(seq_count == 1 ~ "singleton", 
    seq_count >= 10 ~ "clone size ≥ 10", T ~ "clone size < 10"))

IGH_pie_df <- df_pie %>% group_by(expand_type) %>% summarise(value = n(), .groups = "drop") %>% 
    mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(expand_type, 
        "\n", " (", percent, ")"))

IGH_pie_df$expand_type <- factor(IGH_pie_df$expand_type, levels = c("singleton", "clone size < 10", 
    "clone size ≥ 10"))

fill_colors <- c("#65a3d1", "#edb55a", "#946625")

names(fill_colors) <- c("singleton", "clone size < 10", "clone size ≥ 10")

total_clones <- sum(IGH_pie_df$value)

p <- ggplot(data = IGH_pie_df) + ggforce::geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, 
    r = 2, amount = ratio, fill = expand_type)) + geom_text(aes(x = 1.5 * sin(2 * pi * (cumsum(ratio) - 
    ratio/2)/sum(ratio)), y = 1.5 * cos(2 * pi * (cumsum(ratio) - ratio/2)/sum(ratio)), label = label), 
    size = 4, color = "black") + annotate("text", x = 0, y = 0, label = glue("total = {total_clones}"), 
    size = 8, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors) + theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) + labs(title = "CRC", 
    x = "", y = "")

dir.create(dirname(glue(file.path(figure_dir, "isotype/scRNA_CRC_expand_IGH"))), recursive = TRUE, 
    showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/scRNA_CRC_expand_IGH")), if_pdf = T, width = 6, 
    height = 6)

all_sample <- c("3992B", "4272T", "4288T", "CLD11", "CLD13", "CLD14", "G617", "G620", "G637", "NL_1", 
    "NL_2", "NL_3")

BCR_heavy <- readRDS(file.path(input_dir, "result/06.SMRTseq/merge_BCR_heavy.rds"))

df_list <- list()

for (cur_sample in c("3992B", "4272T", "4288T")) {
    df_barcode <- read.csv(glue(file.path(input_dir, "spatial_bin/Convertbin/{cur_sample}_dedup.info.BinID.50.csv")), 
        sep = "\t")
    df_list[[cur_sample]] <- BCR_heavy %>% filter(sampleid == cur_sample & barcode %in% df_barcode$BC)
}

df_merge <- bind_rows(df_list)

df_merge <- df_merge %>% mutate(isotype = str_split_i(c_call, "\\*", 1))

df_merge$isotype[df_merge$isotype == "IGHG4A"] <- "IGHG4"

df_merge <- filter(df_merge, isotype %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

table(df_merge$isotype)

pie_data <- df_merge %>% count(isotype) %>% mutate(Percentage = n/sum(n) * 100, label = paste0(isotype, 
    "\n", round(Percentage, 1), "%"))

pie_df <- df_merge %>% dplyr::group_by(isotype) %>% dplyr::summarise(value = n())

pie_df$ratio <- pie_df$value/sum(pie_df$value)

pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)

pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", 
    "IGHA2", "IGHM", "IGHD", "IGHE"))

fill_colors <- c("#fdc7d4", "#ef8e89", "#f6ba59", "#9fc0d6", "#aacb89", "#5c9b4a", "#7775b6", "#caabd0", 
    "#225ea8")

names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs, lab.pos = "in", lab.font = "black", fill = "isotype", 
    color = NULL, palette = fill_colors) + labs(title = "CRC")

p

dir.create(dirname(glue(file.path(figure_dir, "isotype/IgH_CRC"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/IgH_CRC")), if_pdf = T, width = 4, height = 4)

df_list <- list()

for (cur_sample in c("3992B", "4272T", "4288T")) {
    print(cur_sample)
    df_barcode <- read.csv(glue(file.path(input_dir, "spatial_bin/Convertbin/{cur_sample}_dedup.info.BinID.50.csv")), 
        sep = "\t")
    df_list[[cur_sample]] <- BCR_heavy %>% filter(sampleid == cur_sample & barcode %in% df_barcode$BC)
}

df_merge <- bind_rows(df_list)

df_merge <- df_merge %>% mutate(isotype = str_split_i(c_call, "\\*", 1))

df_merge$isotype[df_merge$isotype == "IGHG4A"] <- "IGHG4"

df_merge <- filter(df_merge, isotype %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", 
    "IGHM", "IGHD", "IGHE"))

table(df_merge$isotype)

df_merge$clone_name <- paste(df_merge$sampleid, df_merge$clone_id, sep = "_")

bcr_smrt_clones <- my_countClones(df_merge, clone = "clone_name")

df_pie <- bcr_smrt_clones %>% dplyr::mutate(expand_type = case_when(seq_count == 1 ~ "singleton", 
    seq_count >= 10 ~ "clone size ≥ 10", T ~ "clone size < 10"))

IGH_pie_df <- df_pie %>% group_by(expand_type) %>% summarise(value = n(), .groups = "drop") %>% 
    mutate(ratio = value/sum(value), percent = scales::percent(ratio, accuracy = 0.1), label = paste0(expand_type, 
        "\n", " (", percent, ")"))

IGH_pie_df$expand_type <- factor(IGH_pie_df$expand_type, levels = c("singleton", "clone size < 10", 
    "clone size ≥ 10"))

fill_colors <- c("#65a3d1", "#edb55a", "#946625")

names(fill_colors) <- c("singleton", "clone size < 10", "clone size ≥ 10")

total_clones <- sum(IGH_pie_df$value)

p <- ggplot(data = IGH_pie_df) + ggforce::geom_arc_bar(stat = "pie", aes(x0 = 0, y0 = 0, r0 = 1, 
    r = 2, amount = ratio, fill = expand_type)) + geom_text(aes(x = 1.5 * sin(2 * pi * (cumsum(ratio) - 
    ratio/2)/sum(ratio)), y = 1.5 * cos(2 * pi * (cumsum(ratio) - ratio/2)/sum(ratio)), label = label), 
    size = 4, color = "black") + annotate("text", x = 0, y = 0, label = glue("total = {total_clones}"), 
    size = 8, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors) + theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) + labs(title = "CRC", 
    x = "", y = "")

dir.create(dirname(glue(file.path(figure_dir, "isotype/Expand_CRC"))), recursive = TRUE, showWarnings = FALSE)

save_plot(p, glue(file.path(figure_dir, "isotype/Expand_CRC")), if_pdf = T, width = 6, height = 6)

