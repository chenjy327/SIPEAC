accuracy_results <- list()
allelic_summary <- fread(glue('/Users/mac/Desktop/Rproject/findlight/version5/final_version/summary_e100.tsv')) %>% as.data.frame()
allelic_summary$repair_accuracy <- NA

for (i in 1:nrow(allelic_summary)) {
  sample_id <- allelic_summary$cur_sample[i]
  threshold_value <- allelic_summary$threshold[i]
  bin_size <- allelic_summary$bin_size[i]
  result_dir <- glue("/Users/mac/Desktop/Rproject/findlight/version5/repair_result/bin{bin_size}/{sample_id}_{threshold_value}")
  
  file_names <- list.files(result_dir)
  if (length(file_names) == 0) next
  
  result_file <- file_names[grepl('analysis_result', file_names)]
  repair_data <- fread(glue('{result_dir}/{result_file}')) %>% as.data.frame()
  repair_data$pair_prediction <- paste0(repair_data$chainA, '_', repair_data$chainB)
  repair_data <- dplyr::distinct(repair_data, pair_prediction, .keep_all = TRUE)
  
  ground_truth <- read.table(glue('/Users/mac/Desktop/Rproject/findlight/version5/grouth_truth/bin{bin_size}/{sample_id}.tsv'),
                             header = TRUE, row.names = NULL)
  
  repair_data <- repair_data[repair_data$chainA %in% ground_truth$Hclone, ]
  allelic_summary$repair_accuracy[i] <- sum(repair_data$pair_prediction %in% ground_truth$pair_name) / nrow(repair_data)
}

allelic_summary <- allelic_summary[!is.na(allelic_summary$repair_accuracy) & !is.na(allelic_summary$acc), ]
allelic_summary <- dplyr::arrange(allelic_summary, bin_size, threshold, cur_sample)
# allelic_summary <- allelic_summary[!allelic_summary$cur_sample %in% c("4266B", "4309B", "B4266", "B4309"),]
new_id <- fread("/Users/mac/Desktop/Rproject/Pacbio/newid-加上标签.csv")
allelic_summary <- dplyr::left_join(allelic_summary, new_id, by = c("cur_sample" = "rawID"))
write.csv(allelic_summary, glue('/Users/mac/Desktop/Rproject/findlight/version5/final_version/compare.csv'), quote = FALSE)

comparison_data <- read.csv(glue('/Users/mac/Desktop/Rproject/findlight/version5/final_version/compare.csv'))

df_allelic <- comparison_data
df_allelic$value <- comparison_data$acc
df_allelic$method <- "allelic"

df_repair <- comparison_data
df_repair$value <- comparison_data$repair_accuracy
df_repair$method <- "repair"

combined_df <- rbind(df_allelic, df_repair)

for (bin in unique(combined_df$bin_size)) {
  for (thresh in unique(combined_df$threshold)) {
    subset_df <- combined_df[combined_df$bin_size == bin & combined_df$threshold == thresh, ]
    acc_allelic <- subset_df$value[subset_df$method == "allelic"]
    acc_repair <- subset_df$value[subset_df$method == "repair"]
    p_val <- round(wilcox.test(acc_allelic, acc_repair, alternative = "two.sided")$p.value, 3)
    
    
    p <- ggbarplot(subset_df, x = "method", y = "value", color = "method",
                   palette = c(allelic="#B03060", repair="#0B7FAB"),
                   add = c("mean_se", "jitter")) +
      theme() +
      labs(x = "", y = "Accuracy", title = glue("p = {p_val}"))
    
    save_plot(p, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/accbar_{bin}_{thresh}"),
              width = 3, height = 4, if_pdf = TRUE)
  }
}
bin <- 110; thresh <- 1;
subset_df <- combined_df[combined_df$bin_size == bin & combined_df$threshold == thresh, ]
dplyr::group_by(subset_df, method) %>% dplyr::summarise(mean = mean(value))


# Line Plot for Accuracies across Thresholds
plot_df <- combined_df[combined_df$bin_size %in% c(20, 50, 110), ]

p_line <- ggline(plot_df, x = "threshold", y = "value", color = "method",
                 add = c("mean_se"), palette = c(allelic="#B03060", repair="#0B7FAB"),
                 facet.by = "bin_size") +
  labs(x = "", y = "") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/accline"),
          width = 10, height = 5, if_pdf = TRUE)
save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/accline"),
          width = 10, height = 5, if_png = TRUE)

# ANOVA Tests for Each Bin Size
for (bin in c(20, 50, 110)) {
  df_bin <- plot_df[plot_df$bin_size == bin, ]
  print(summary(aov(value ~ method, data = df_bin)))
}

# Individual Line Plots by Bin Size
for (bin in unique(comparison_data$bin_size)) {
  p_bin <- ggline(combined_df[combined_df$bin_size == bin, ], "threshold", "value", color = "method",
                  add = c("mean_se"), palette = c(allelic="#B03060", repair="#0B7FAB"), facet.by = "bin_size") +
    labs(x = "", y = "", title = glue("Bin {bin}")) +
    theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))
  
  save_plot(p_bin, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/accline_{bin}"),
            width = 4, height = 4, if_png = TRUE)
  save_plot(p_bin, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/accline_{bin}"),
            width = 4, height = 4, if_pdf = TRUE)
}
bin <- 20;threshold = 1;
subset_df <- combined_df[combined_df$bin_size == bin & combined_df$threshold == thresh, ]
dplyr::group_by(subset_df, method) %>% dplyr::summarise(mean = mean(value))



library(dplyr)
library(ggpubr)
library(glue)

# 1. 读取原始数据
data_raw <- read.csv('/Users/mac/Desktop/Rproject/findlight/version5/final_version/compare.csv')

# 2. 创建唯一实验组标识符（样本_阈值）
data_raw$sample_id <- paste(data_raw$cur_sample, data_raw$threshold, sep = "_")

# 3. 按阈值和eff排序；去除缺失acc值的数据
data_clean <- data_raw %>%
  arrange(threshold, eff) %>%
  filter(!is.na(acc))

# 4. 构建实验组数据（不同bin_size下的eff值）
exp_data <- data_clean %>%
  mutate(
    value = eff,
    group = paste0("Bin", bin_size)
  )

# 5. 构建“修复”组（对照组），value 设为常数 1
repair_data <- data_clean %>%
  mutate(
    value = 1,
    group = "repair"
  )

# 6. 合并数据
combined_data <- bind_rows(exp_data, repair_data)

# 7. 仅保留指定bin_size的数据
bin_list <- c(20, 50, 110)
combined_data <- combined_data %>%
  filter(bin_size %in% bin_list | group == "repair") %>%
  mutate(
    group = factor(group, levels = c("Bin20", "Bin50", "Bin110", "repair")),
    threshold = as.character(threshold)
  )

# 8. 绘图：按group显示eff值随threshold的变化
p <- ggline(
  combined_data, 
  x = "threshold", 
  y = "value", 
  color = "group",
  add = "mean_se",
  palette = c(
    "Bin110" = "#c51b7d",
    "Bin50" = "#6a3d9a",
    "Bin20" = "#EE3743",
    "repair" = "#0B7FAB"
  )
)

# 9. 保存图像为 PNG 和 PDF
output_path <- "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/eff_line_2"
# save_plot(p, glue("{output_path}"), width = 4, height = 4, if_png = TRUE)
save_plot(p, glue("{output_path}"), width = 4, height = 4, if_pdf = TRUE)


library(extrafont)
# 确保你已经导入字体
font_import() # 只需第一次运行
loadfonts()

# 设置字体
par(family = "Arial Unicode MS")  # 或其他支持μ的字体

comparison_data <- read.csv(glue('/Users/mac/Desktop/Rproject/findlight/version5/final_version/compare.csv'))

df_allelic <- comparison_data
df_allelic$value <- comparison_data$acc
df_allelic$method <- "allelic"

df_repair <- comparison_data
df_repair$value <- comparison_data$repair_accuracy
df_repair$method <- "repair"

combined_df <- rbind(df_allelic, df_repair)

plot_df <- combined_df[combined_df$bin_size %in% c(20, 50, 110), ]
plot_df$resolution <- plot_df$bin_size
plot_df$resolution[plot_df$bin_size == 20] <- "10\u03bcm"
plot_df$resolution[plot_df$bin_size == 50] <- "25\u03bcm"
plot_df$resolution[plot_df$bin_size == 110] <- "55\u03bcm"
plot_df$resolution <- factor(plot_df$resolution, levels = c("10\u03bcm", "25\u03bcm", "55\u03bcm"))


plot_df$threshold <- as.character(plot_df$threshold)
plot_df$threshold <- factor(plot_df$threshold, levels = c("1", "2", "5", "10", "15", "20"))
for(current_sample in unique(plot_df$newID)){
  p_line <- ggline(plot_df[plot_df$newID == current_sample, ], x = "threshold", y = "value", color = "method",
                   palette = c(allelic="#B03060", repair="#0B7FAB"),
                   facet.by = c("resolution")) +
    labs(x = "", y = "") +
    theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))
  
  save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/31_sample_acc/{current_sample}"),
            width = 10, height = 5, if_pdf = TRUE)

}

p_line <- ggline(plot_df, x = "threshold", y = "value", color = "method",
                 palette = c(allelic="#B03060", repair="#0B7FAB"),
                 facet.by = c( "newID", "resolution")) +
  labs(x = "", y = "") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/31_sample_acc/merge"),
          width = 10, height = 40, if_pdf = TRUE)

### 31 merge ====


# fils <- list.files("/Users/mac/Desktop/Rproject/Pacbio/Result/DedupBin")
# result_list <- list()
# 
# for (i in fils) {
#   tb <- read.table(glue("/Users/mac/Desktop/Rproject/Pacbio/Result/DedupBin/{i}"), header = TRUE)
#   
#   tb2 <- table(tb$binid) %>% as.data.frame()
#   colnames(tb2) <- c("binid", "Freq")
#   
#   
#   result_list[[i]] <- tb2
# }
# 
# combined_result <- bind_rows(result_list, .id = "filename")
# combined_result$filename <- substr(combined_result$filename,1,5)
# 
tb_info <- fread("/Users/mac/Desktop/Rproject/Pacbio/newid-加上标签.csv")
# combined_result <- merge(combined_result, tb_info, by.x = "filename", by.y = "rawID", all.x = T)
# 
# combined_result$newID <- factor(combined_result$newID, levels = tb_info$newID)
# saveRDS(combined_result, "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/combined_result.rds")
combined_result <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/combined_result.rds")

# 
# combined_result_mean <- dplyr::group_by(combined_result, newID) %>% 
#   dplyr::summarise(freq_mean  = mean(Freq))
# save_table(combined_result_mean, "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/UMI_count_per_Bin50")
  
combined_result_mean <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/UMI_count_per_Bin50.txt")

# merge_bcr <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr.rds")
# merge_bcr <- merge(merge_bcr, tb_info, by.x = "sample", by.y = "rawID", all.x = T)
# bcr_count <- table(merge_bcr$newID) %>% as.data.frame()
# colnames(bcr_count) <- c("newID", "bcr_num")
# bcr_count <- merge_bcr %>%
#   dplyr::group_by(newID) %>% 
#   summarise(
#     n_IgH = sum(locus == "IGH", na.rm = TRUE),
#     n_IgL = sum(locus != "IGH", na.rm = TRUE),
#     bcr_num = n()
#   )
# merge_bcr_clone <- merge_bcr %>% 
#   dplyr::distinct(newID, clone_id, .keep_all = T) %>% 
#   dplyr::group_by(newID) %>% 
#   summarise(
#     n_IgH_clone = sum(locus == "IGH", na.rm = TRUE),
#     n_IgL_clone = sum(locus != "IGH", na.rm = TRUE),
#     clone_num = n()
#   )
# bcr_count <- dplyr::left_join(bcr_count, merge_bcr_clone)
# save_table(bcr_count,  "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/45_sample_BCR_info")
bcr_count <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/45_sample_BCR_info.txt")
# merge_anno <- fread("/Users/mac/Desktop/01finalbackup/Bin50/Bin50_classic_celltype_metadata_merge.txt")
# merge_anno <- merge_anno[merge_anno$SampleID %in% merge_bcr$sample,]
# merge_anno <- merge(merge_anno, tb_info, by.x = "SampleID", by.y = "rawID", all.x = T)
# cell_prop <- merge_anno %>% 
#   dplyr::group_by(newID) %>% 
#   dplyr::summarize(
#     total = n(),
#     count_B_Plasma = sum(SpotLight_Anno %in% c("B", "Plasma")),
#     proportion = count_B_Plasma / total
#   )
# save_table(cell_prop, "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/cell_prop")
cell_prop <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/cell_prop.txt")

# merge_anno_count <- table(merge_anno$newID, merge_anno$SpotLight_Anno) %>% as.data.frame()
# colnames(merge_anno_count) <- c("newID", "cell_type", "Freq")
comparison_data <- read.csv(glue('/Users/mac/Desktop/Rproject/findlight/version5/final_version/compare.csv'))
comparison_data <- dplyr::left_join(comparison_data, cell_prop, by = "newID")
comparison_data <- dplyr::left_join(comparison_data, bcr_count, by = "newID")
comparison_data <- dplyr::left_join(comparison_data, combined_result_mean, by = "newID")

acc_mean <- comparison_data %>% 
  dplyr::filter(threshold == 1) %>% 
  dplyr::group_by(newID) %>% 
  dplyr::summarise(acc_mean = mean(acc)) %>% 
  dplyr::arrange(desc(acc_mean))

comparison_data$newID <- factor(comparison_data$newID, levels = acc_mean$newID)
comparison_data$n_Shared <- sapply(comparison_data$cur_sample, function(sample_id) {
  file_path <- file.path(
    "/Users/mac/Desktop/Rproject/findlight/version5/grouth_truth/bin50",
    paste0(sample_id, ".tsv")
  )
  
  if (file.exists(file_path)) {
    nrow(read.table(file_path, header = TRUE, sep = "\t"))
  } else {
    NA  # 文件不存在时返回 NA
  }
})
summary_data <- dplyr::arrange(comparison_data, newID) %>% 
  dplyr::distinct(newID, .keep_all = T)


summary_data <- summary_data[, c("newID", "total", "count_B_Plasma", "proportion", "bcr_num", "freq_mean")]
summary_data <- dplyr::left_join(summary_data, acc_mean)
save_table(summary_data, "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/summary_31_data")

df_allelic <- comparison_data
df_allelic$value <- comparison_data$acc
df_allelic$method <- "allelic"

df_repair <- comparison_data
df_repair$value <- comparison_data$repair_accuracy
df_repair$method <- "repair"

combined_df <- rbind(df_allelic, df_repair)

plot_df <- combined_df[combined_df$bin_size %in% c(20, 50, 110), ]
plot_df$resolution <- plot_df$bin_size
plot_df$resolution[plot_df$bin_size == 20] <- "10\u03bcm"
plot_df$resolution[plot_df$bin_size == 50] <- "25\u03bcm"
plot_df$resolution[plot_df$bin_size == 110] <- "55\u03bcm"
plot_df$resolution <- factor(plot_df$resolution, levels = c("10\u03bcm", "25\u03bcm", "55\u03bcm"))
plot_df$newID <- factor(plot_df$newID, levels = acc_mean$newID)

plot_df$threshold <- as.character(plot_df$threshold)
plot_df$threshold <- factor(plot_df$threshold, levels = c("1", "2", "5", "10", "15", "20"))


p_line <- ggline(plot_df, x = "threshold", y = "value", color = "method",
                 palette = c(allelic="#B03060", repair="#0B7FAB"),
                 facet.by = c( "newID", "resolution")) +
  labs(x = "", y = "") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/31_sample_acc/merge"),
          width = 10, height = 40, if_pdf = TRUE)

### D v2 =====
comparison_data <- read.csv(glue('/Users/mac/Desktop/Rproject/findlight/version5/final_version/compare.csv'))
cell_prop <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/cell_prop.txt")
combined_result <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/combined_result.rds")
combined_result_mean <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/UMI_count_per_Bin50.txt")
bcr_count <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/45_sample_BCR_info.txt")
comparison_data$n_Shared <- sapply(comparison_data$cur_sample, function(sample_id) {
  file_path <- file.path(
    "/Users/mac/Desktop/Rproject/findlight/version5/grouth_truth/bin50",
    paste0(sample_id, ".tsv")
  )
  
  if (file.exists(file_path)) {
    nrow(read.table(file_path, header = TRUE, sep = "\t"))
  } else {
    NA  # 文件不存在时返回 NA
  }
})
comparison_data <- dplyr::left_join(comparison_data, cell_prop, by = "newID")
comparison_data <- dplyr::left_join(comparison_data, bcr_count, by = "newID")
comparison_data <- dplyr::left_join(comparison_data, combined_result_mean, by = "newID")

acc_mean <- comparison_data %>% 
  dplyr::filter(threshold == 1) %>% 
  dplyr::group_by(newID) %>% 
  dplyr::summarise(acc_mean = mean(acc)) %>% 
  dplyr::arrange(desc(acc_mean))
acc_mean <- as.data.frame(acc_mean)
rows_to_remove <- c(7,10,12,13,14,16,19,21,22,23,24,25,30,31)
sample_14 <- acc_mean$newID[c(7,10,12,13,14,16,19,21,22,23,24,25,30,31)]


comparison_out <- comparison_data %>% 
  dplyr::select(newID, threshold, bin_size,
                acc, repair_accuracy,
                proportion, bcr_num, freq_mean,
                n_IgH ,n_IgL ,bcr_num,        
                n_IgH_clone, n_IgL_clone, clone_num,
                n_Shared) %>% 
  dplyr::rename(clone_size = threshold, 
                B_Plasma_proportion = proportion,
                UMI_count_per_Bin50 = freq_mean,
                allelic_accuracy = acc,
                n_bcr = bcr_num,
                n_clone = clone_num)
library(openxlsx)
comparison_out$Sample_14 <- ""
comparison_out$Sample_14[comparison_out$newID %in% sample_14] <- "True"
write.xlsx(comparison_out %>% dplyr::distinct(newID, .keep_all = T), 
           file = "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/31sample_Info.xlsx")
# save_table(comparison_out, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/Figure2D_31"))
# write.xlsx(comparison_out[comparison_out$newID %in% sample_14,] %>% dplyr::distinct(newID, .keep_all = T),
#            file = "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/14sample_Info.xlsx")

acc_mean <- acc_mean[-rows_to_remove,]
comparison_data <- comparison_data[comparison_data$newID %in% acc_mean$newID,]
comparison_data$newID <- factor(comparison_data$newID, levels = acc_mean$newID)

df_allelic <- comparison_data
df_allelic$value <- comparison_data$acc
df_allelic$method <- "allelic"

df_repair <- comparison_data
df_repair$value <- comparison_data$repair_accuracy
df_repair$method <- "repair"

combined_df <- rbind(df_allelic, df_repair)

plot_df <- combined_df[combined_df$bin_size %in% c(20, 50, 110), ]


p_line <- ggline(plot_df, x = "threshold", y = "value", color = "method",
                 add = c("mean_se"), palette = c(allelic="#B03060", repair="#0B7FAB"),
                 facet.by = "bin_size") +
  labs(x = "", y = "") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/FigureD_17sample"),
          width = 10, height = 5, if_pdf = TRUE)
save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/FigureD_17sample"),
          width = 10, height = 5, if_png = TRUE)

p_box <- ggboxplot(plot_df, x = "threshold", y = "value", color = "method",
                 palette = c(allelic="#B03060", repair="#0B7FAB"),
                 facet.by = "bin_size", width = 0.5) +
  labs(x = "", y = "") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/D_box"),
          width = 15, height = 5, if_pdf = TRUE)
save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/D_box"),
          width = 10, height = 5, if_png = TRUE)


# ANOVA Tests for Each Bin Size
for (bin in c(20, 50, 110)) {
  df_bin <- plot_df[plot_df$bin_size == bin, ]
  print(summary(aov(value ~ method, data = df_bin)))
}


#### 31个样本Umi per bin50和准确率的关系 ====
comparison_data <- read.csv(glue('/Users/mac/Desktop/Rproject/findlight/version5/final_version/compare.csv'))
cell_prop <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/cell_prop.txt")
combined_result <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/combined_result.rds")
combined_result_mean <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/UMI_count_per_Bin50.txt")
bcr_count <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/45_sample_BCR_info.txt")
comparison_data$n_Shared <- sapply(comparison_data$cur_sample, function(sample_id) {
  file_path <- file.path(
    "/Users/mac/Desktop/Rproject/findlight/version5/grouth_truth/bin50",
    paste0(sample_id, ".tsv")
  )
  
  if (file.exists(file_path)) {
    nrow(read.table(file_path, header = TRUE, sep = "\t"))
  } else {
    NA  # 文件不存在时返回 NA
  }
})
comparison_data <- dplyr::left_join(comparison_data, cell_prop, by = "newID")
comparison_data <- dplyr::left_join(comparison_data, bcr_count, by = "newID")
comparison_data <- dplyr::left_join(comparison_data, combined_result_mean, by = "newID")

acc_mean <- comparison_data %>% 
  dplyr::filter(threshold == 1) %>% 
  dplyr::group_by(newID) %>% 
  dplyr::summarise(acc_mean = mean(acc)) %>% 
  dplyr::arrange(desc(acc_mean))
acc_mean <- as.data.frame(acc_mean)
rows_to_remove <- c(7,10,12,13,14,16,19,21,22,23,24,25,30,31)
sample_14 <- acc_mean$newID[c(7,10,12,13,14,16,19,21,22,23,24,25,30,31)]


plot_df <- comparison_data %>% 
  dplyr::rename(clone_size = threshold, 
                B_Plasma_proportion = proportion,
                UMI_count_per_Bin50 = freq_mean,
                allelic_accuracy = acc,
                n_bcr = bcr_num,
                n_clone = clone_num)
plot_df <- plot_df[plot_df$clone_size == 1 & plot_df$bin_size == 50,]
correlation <- round(cor(plot_df$UMI_count_per_Bin50, plot_df$allelic_accuracy),3)
lm_model <- lm(UMI_count_per_Bin50 ~ allelic_accuracy, data = plot_df)
summary(lm_model)
p <- ggplot(plot_df, aes(x=UMI_count_per_Bin50, y=allelic_accuracy)) + 
  geom_point(color="#B03060") +
  geom_smooth(method= "lm",color="black")+
  labs(title="", x=glue("UMI count per Bin50"), y=glue("Allelic accuracy"),
       subtitle = paste("r=", correlation))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    axis.line = element_line(size = 0.6),
    axis.title = element_text(size = 15, color = "black", face = "bold"),
    legend.text = element_text(size = 10, color = "black", face = "bold"),
    legend.position="none",
  )
save_plot(p, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/UMI_accuracy"),
          width = 5, height = 5, if_png = TRUE)


#### 克隆数等于10的情况下，25um分辨率情况下这31个样品的正确率跟共享BCR数目的相关性呢 ===
comparison_data <- read.csv(glue('/Users/mac/Desktop/Rproject/findlight/version5/final_version/compare.csv'))
cell_prop <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/cell_prop.txt")
combined_result <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/combined_result.rds")
combined_result_mean <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/UMI_count_per_Bin50.txt")
bcr_count <- fread("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/45_sample_BCR_info.txt")
comparison_data$n_Shared <- sapply(comparison_data$cur_sample, function(sample_id) {
  file_path <- file.path(
    "/Users/mac/Desktop/Rproject/findlight/version5/grouth_truth/bin50",
    paste0(sample_id, ".tsv")
  )
  
  if (file.exists(file_path)) {
    nrow(read.table(file_path, header = TRUE, sep = "\t"))
  } else {
    NA  # 文件不存在时返回 NA
  }
})
comparison_data <- dplyr::left_join(comparison_data, cell_prop, by = "newID")
comparison_data <- dplyr::left_join(comparison_data, bcr_count, by = "newID")
comparison_data <- dplyr::left_join(comparison_data, combined_result_mean, by = "newID")
comparison_data <- comparison_data %>% 
  dplyr::rename(clone_size = threshold, 
                B_Plasma_proportion = proportion,
                UMI_count_per_Bin50 = freq_mean,
                allelic_accuracy = acc,
                n_bcr = bcr_num,
                n_clone = clone_num)

merge_bcr <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr.rds")
merge_bcr <- merge(merge_bcr, tb_info, by.x = "sample", by.y = "rawID", all.x = T)

summary_bind <- list()

for(current_sample in comparison_data$cur_sample %>% unique()){
  grount_truth <- fread(glue("/Users/mac/Desktop/Rproject/findlight/version5/grouth_truth/bin50/{current_sample}.tsv"))
  s_df <- comparison_data[comparison_data$cur_sample ==current_sample,]
  new_id <-  s_df$newID[1]
  for(current_bin_size in s_df$bin_size %>% unique()){
    s_df_1 <- s_df[s_df$bin_size ==current_bin_size,]
    for(current_clone_size in s_df_1$clone_size %>% unique()){
      test_out <- fread(glue("/Users/mac/Desktop/Rproject/findlight/version5/final_version/bin{current_bin_size}/{current_sample}_test_e100.tsv"))
      sub_out <- test_out[test_out$bin_num >= current_clone_size, ]
      # n_Shared <- comparison_data$n_Shared[comparison_data$cur_sample ==current_sample][1]
      n_Shared_Allelic_Predict <- nrow(sub_out)
      n_True <- sum(sub_out$pair_name %in% grount_truth$pair_name)
      n_Error <- sum(!sub_out$pair_name %in% grount_truth$pair_name)
      accuracy <- n_True/n_Shared_Allelic_Predict
      sub_df <- data.frame(
        newID = new_id,
        cur_sample = current_sample,
        clone_size = current_clone_size,
        bin_size = current_bin_size,
        n_Shared_Allelic_Predict = n_Shared_Allelic_Predict,
        n_True = n_True,
        n_Error = n_Error,
        allelic_accuracy = accuracy
      )
      summary_bind[[glue("{current_sample}_{current_clone_size}_{current_bin_size}")]] <- sub_df
    }
  }
}
summary <- dplyr::bind_rows(summary_bind)
info <- comparison_data[,c("newID", "B_Plasma_proportion", "n_IgH",
                           "n_IgL",              
                           "n_bcr",              
                           "n_IgH_clone",        
                           "n_IgL_clone",        
                           "n_clone",            
                           "UMI_count_per_Bin50")] %>% 
  dplyr::distinct(newID, .keep_all = T)
summary <- summary %>% 
  dplyr::left_join(info)
summary$Sample_14 <- ""
summary$Sample_14[summary$newID %in% sample_14] <- "True"

write.xlsx(summary, 
           file = "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/summary_shared.xlsx")

plot_df <- summary[summary$clone_size == "10" & summary$bin_size == "50",]
correlation <- round(cor(plot_df$n_Shared_Allelic_Predict, plot_df$allelic_accuracy),3)
lm_model <- lm(n_Shared_Allelic_Predict ~ allelic_accuracy, data = plot_df)
summary(lm_model)
p <- ggplot(plot_df, aes(x=n_Shared_Allelic_Predict, y=allelic_accuracy)) + 
  geom_point(color="#B03060") +
  geom_smooth(method= "lm",color="black")+
  labs(title="", x=glue("n_Shared_Allelic_Predict"), y=glue("Allelic accuracy"),
       subtitle = paste("r=", correlation))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    axis.line = element_line(size = 0.6),
    axis.title = element_text(size = 15, color = "black", face = "bold"),
    legend.text = element_text(size = 10, color = "black", face = "bold"),
    legend.position="none",
  )
save_plot(p, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_accuracy"),
          width = 5, height = 5, if_png = TRUE)

#### 共享和 acc的关系 ====

summary_shared <- read_xlsx("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/summary_shared.xlsx")
plot_df <- summary_shared[summary_shared$bin_size %in% c(20, 50, 110) & summary_shared$clone_size == 1, ]
quantile(plot_df$n_Shared_Allelic_Predict[plot_df$bin_size == 20])

plot_df <- plot_df %>%
  mutate(group = case_when(
    n_Shared_Allelic_Predict >= 16 ~ ">=16",
    n_Shared_Allelic_Predict >= 10 ~ ">=10",
    n_Shared_Allelic_Predict >= 6  ~ ">=6",
    n_Shared_Allelic_Predict >= 3  ~ ">=3",
    TRUE ~ ">=1"
  )) %>%
  arrange(factor(group, levels = c(">=1",  ">=3", ">=6", ">=10", ">=16")))

df_levels <- bind_rows(
  plot_df %>% filter(n_Shared_Allelic_Predict >= 1)  %>% mutate(group = ">=1"),
  plot_df %>% filter(n_Shared_Allelic_Predict >= 3)  %>% mutate(group = ">=3"),
  plot_df %>% filter(n_Shared_Allelic_Predict >= 6)  %>% mutate(group = ">=6"),
  plot_df %>% filter(n_Shared_Allelic_Predict >= 10)  %>% mutate(group = ">=10"),
  plot_df %>% filter(n_Shared_Allelic_Predict >= 16) %>% mutate(group = ">=16")
)


df_levels <- df_levels[df_levels$clone_size == 1,]
group_stats <- df_levels %>%
  group_by(bin_size, group) %>%
  summarise(n_sample = n(), .groups = "drop") %>%
  arrange(bin_size, factor(group, levels = c(">=1",  ">=3", ">=6", ">=10", ">=16")))
group_stats

p_line <- ggline(df_levels, x = "group", y = "allelic_accuracy", color = "#B03060",
                 add = c("mean_se"), 
                 facet.by = "bin_size") +
  labs(x = "n_shared", y = "allelic_accuracy") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_acc_line"),
          width = 10, height = 5, if_pdf = TRUE)
save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_acc_line"),
          width = 10, height = 5, if_png = TRUE)

p_box <- ggboxplot(df_levels, x = "group", y = "allelic_accuracy", color = "#B03060",
                 facet.by = "bin_size") +
  labs(x = "n_shared", y = "allelic_accuracy") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_acc_box"),
          width = 10, height = 5, if_pdf = TRUE)
save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_acc_box"),
          width = 10, height = 5, if_png = TRUE)

#### 共享和 acc的关系 v2====

summary_shared <- read_xlsx("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/summary_shared.xlsx")
plot_df <- summary_shared[summary_shared$bin_size %in% c(20, 50, 110) & summary_shared$clone_size == 1, ]
quantile(plot_df$n_Shared_Allelic_Predict[plot_df$bin_size == 20])

# plot_df <- plot_df %>%
#   mutate(group = case_when(
#     n_Shared_Allelic_Predict >= 16 ~ ">=16",
#     n_Shared_Allelic_Predict >= 10 ~ ">=10",
#     n_Shared_Allelic_Predict >= 6  ~ ">=6",
#     n_Shared_Allelic_Predict >= 3  ~ ">=3",
#     n_Shared_Allelic_Predict >= 1  ~ ">=1"
#   )) %>%
#   arrange(factor(group, levels = c(">=1",  ">=3", ">=6", ">=10", ">=16")))

df_levels <- bind_rows(
  plot_df %>% filter(n_Shared_Allelic_Predict >= 1 & n_Shared_Allelic_Predict < 3)  %>%
    mutate(group = "1-2"),
  plot_df %>% filter(n_Shared_Allelic_Predict >= 3 & n_Shared_Allelic_Predict < 6)  %>%
    mutate(group = "3-5"),
  plot_df %>% filter(n_Shared_Allelic_Predict >= 6 & n_Shared_Allelic_Predict < 10)  %>%
    mutate(group = "6-9"),
  plot_df %>% filter(n_Shared_Allelic_Predict >= 10 & n_Shared_Allelic_Predict < 16)  %>% 
    mutate(group = "10-15"),
  plot_df %>% filter(n_Shared_Allelic_Predict >= 16) %>% 
    mutate(group = "16-")
)
df_levels$group <- factor(df_levels$group, levels = c("1-2",  "3-5", "6-9", "10-15", "16-"))

df_levels <- df_levels[df_levels$clone_size == 1,]
group_stats <- df_levels %>%
  group_by(bin_size, group) %>%
  summarise(n_sample = n(), .groups = "drop") %>%
  arrange(bin_size, factor(group, levels = c("1-2",  "3-5", "6-9", "10-15", "16-")))
group_stats

p_line <- ggline(df_levels, x = "group", y = "allelic_accuracy", color = "#B03060",
                 add = c("mean_se"), 
                 facet.by = "bin_size") +
  labs(x = "n_shared", y = "allelic_accuracy") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_acc_line_v2"),
          width = 10, height = 5, if_pdf = TRUE)
save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_acc_line_v2"),
          width = 10, height = 5, if_png = TRUE)

p_box <- ggboxplot(df_levels, x = "group", y = "allelic_accuracy", color = "#B03060",
                   facet.by = "bin_size") +
  labs(x = "n_shared", y = "allelic_accuracy") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_acc_box_v2"),
          width = 10, height = 5, if_pdf = TRUE)
save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/n_shared_acc_box_v2"),
          width = 10, height = 5, if_png = TRUE)


####  UMI per Bin50 spot 和 acc的关系 ====

summary_shared <- read_xlsx("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/summary_shared.xlsx")
plot_df <- summary_shared[summary_shared$bin_size %in% c(50) & summary_shared$clone_size == 1, ]
quantile(plot_df$UMI_count_per_Bin50)
umi_group <- quantile(plot_df$UMI_count_per_Bin50)
plot_df <- plot_df %>%
  mutate(group = case_when(
    UMI_count_per_Bin50 >= 20 ~ ">=20",
    UMI_count_per_Bin50 >= 13  ~ ">=13",
    UMI_count_per_Bin50 >= 8  ~ ">=8",
    TRUE ~ ">=2"
  )) %>%
  arrange(factor(group, levels = c(">=2",  ">=8", ">=13", ">=20")))

df_levels <- bind_rows(
  plot_df %>% filter(UMI_count_per_Bin50 >= 2)  %>% mutate(group = ">=2"),
  plot_df %>% filter(UMI_count_per_Bin50 >= 8)  %>% mutate(group = ">=8"),
  plot_df %>% filter(UMI_count_per_Bin50 >= 13)  %>% mutate(group = ">=13"),
  plot_df %>% filter(UMI_count_per_Bin50 >= 20)  %>% mutate(group = ">=20")
)

df_levels <- df_levels[df_levels$clone_size == 1,]
group_stats <- df_levels %>%
  group_by(bin_size, group) %>%
  summarise(n_sample = n(), .groups = "drop") %>%
  arrange(bin_size, factor(group, levels = c(">=2",  ">=8", ">=13", ">=20")))
group_stats

p_line <- ggline(df_levels, x = "group", y = "allelic_accuracy", color = "#B03060",
                 add = c("mean_se")) +
  labs(x = "UMI_count_per_Bin50", y = "allelic_accuracy") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/umi_acc_line"),
          width = 4, height = 4, if_pdf = TRUE)
save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/umi_acc_line"),
          width = 4, height = 4, if_png = TRUE)

p_box <- ggboxplot(df_levels, x = "group", y = "allelic_accuracy", color = "#B03060") +
  labs(x = "UMI_count_per_Bin50", y = "allelic_accuracy") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/umi_acc_box"),
          width = 4, height = 4, if_pdf = TRUE)
save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/umi_acc_box"),
          width = 4, height = 4, if_png = TRUE)



####  UMI per Bin50 spot 和 acc的关系 v2====

summary_shared <- read_xlsx("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/summary_shared.xlsx")
plot_df <- summary_shared[summary_shared$bin_size %in% c(50) & summary_shared$clone_size == 1, ]
quantile(plot_df$UMI_count_per_Bin50)
umi_group <- quantile(plot_df$UMI_count_per_Bin50)
# plot_df <- plot_df %>%
#   mutate(group = case_when(
#     UMI_count_per_Bin50 >= 20 ~ ">=20",
#     UMI_count_per_Bin50 >= 13  ~ ">=13",
#     UMI_count_per_Bin50 >= 8  ~ ">=8",
#     TRUE ~ ">=2"
#   )) %>%
#   arrange(factor(group, levels = c(">=2",  ">=8", ">=13", ">=20")))

df_levels <- bind_rows(
  plot_df %>% filter(UMI_count_per_Bin50 >= 2 & UMI_count_per_Bin50 < 8)  %>% mutate(group = "2-8"),
  plot_df %>% filter(UMI_count_per_Bin50 >= 8 & UMI_count_per_Bin50 < 13)  %>% mutate(group = "8-13"),
  plot_df %>% filter(UMI_count_per_Bin50 >= 13 & UMI_count_per_Bin50 < 20)  %>% mutate(group = "13-20"),
  plot_df %>% filter(UMI_count_per_Bin50 >= 20)  %>% mutate(group = ">=20")
)

df_levels <- df_levels[df_levels$clone_size == 1,]
group_stats <- df_levels %>%
  group_by(bin_size, group) %>%
  summarise(n_sample = n(), .groups = "drop") %>%
  arrange(bin_size, factor(group, levels = c("2-8", "8-13", "13-20", ">=20")))
group_stats


wilcox.test(df_levels$allelic_accuracy[df_levels$group %in% c("2-8","8-13" )], 
              df_levels$allelic_accuracy[df_levels$group %in% c("13-20",">=20" )])
t.test(df_levels$allelic_accuracy[df_levels$group %in% c("2-8","8-13" )], 
            df_levels$allelic_accuracy[df_levels$group %in% c("13-20",">=20" )])
wilcox.test(df_levels$allelic_accuracy[df_levels$group == "8-13"],
            df_levels$allelic_accuracy[df_levels$group == "13-20"])

wilcox.test(df_levels$allelic_accuracy[df_levels$group == "13-20"],
            df_levels$allelic_accuracy[df_levels$group == ">=20"])

pairwise_result <- pairwise.wilcox.test(
  df_levels$allelic_accuracy, 
  df_levels$group,
  p.adjust.method = "BH"   # Benjamini-Hochberg 校正
)
pairwise_result


p_line <- ggline(df_levels, x = "group", y = "allelic_accuracy", color = "#B03060",
                 add = c("mean_se")) +
  labs(x = "UMI_count_per_Bin50", y = "allelic_accuracy") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/umi_acc_line_v2"),
          width = 5, height = 5, if_pdf = TRUE)
save_plot(p_line, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/umi_acc_line_v2"),
          width = 5, height = 5, if_png = TRUE)

p_box <- ggboxplot(df_levels, x = "group", y = "allelic_accuracy", color = "#B03060") +
  labs(x = "UMI_count_per_Bin50", y = "allelic_accuracy") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/umi_acc_box_v2"),
          width = 4, height = 4, if_pdf = TRUE)
save_plot(p_box, glue("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/umi_acc_box_v2"),
          width = 4, height = 4, if_png = TRUE)


library(ggpubr)

ggboxplot(df_levels, x = "group", y = "allelic_accuracy", fill = "group") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("2-8", "8-13"),
                                        c("2-8", "13-20"),
                                        c("2-8", ">=20"),
                                        c("8-13", "13-20"),
                                        c("8-13", ">=20"),
                                        c("13-20", ">=20")),
                     label = "p.signif") +
  theme_minimal()


anova_model <- aov(allelic_accuracy ~ group, data = df_levels)

# ANOVA 结果
summary(anova_model)
TukeyHSD(anova_model)



#### 给李斌师兄的数据 =====
summary_shared <- read_xlsx("/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/summary_shared.xlsx")

df_summary <- bind_rows(
  summary_shared %>% filter(UMI_count_per_Bin50 >= 2 & UMI_count_per_Bin50 < 8)  %>% mutate(umi_group = "2-8"),
  summary_shared %>% filter(UMI_count_per_Bin50 >= 8 & UMI_count_per_Bin50 < 13)  %>% mutate(umi_group = "8-13"),
  summary_shared %>% filter(UMI_count_per_Bin50 >= 13 & UMI_count_per_Bin50 < 20)  %>% mutate(umi_group = "13-20"),
  summary_shared %>% filter(UMI_count_per_Bin50 >= 20)  %>% mutate(umi_group = ">=20")
)
df_summary <- df_summary[, c("newID", "clone_size", "bin_size", "umi_group", "allelic_accuracy")] %>% 
  dplyr::rename(accuracy = allelic_accuracy)
library(openxlsx)
write.xlsx(df_summary, 
           file = "/Users/mac/Desktop/Rproject/Pacbio/文章/结果图/Fig2/Fig2-CDE/31个样本准确性.xlsx")


