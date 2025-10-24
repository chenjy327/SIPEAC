#### Figure 2E ====

accuracy_results <- list()
allelic_summary <- fread(glue('~/final_version/summary_e100.tsv')) %>% as.data.frame()
allelic_summary$repair_accuracy <- NA

for (i in 1:nrow(allelic_summary)) {
  sample_id <- allelic_summary$cur_sample[i]
  threshold_value <- allelic_summary$threshold[i]
  bin_size <- allelic_summary$bin_size[i]
  result_dir <- glue("~/repair_result/bin{bin_size}/{sample_id}_{threshold_value}")
  
  file_names <- list.files(result_dir)
  if (length(file_names) == 0) next
  
  result_file <- file_names[grepl('analysis_result', file_names)]
  repair_data <- fread(glue('{result_dir}/{result_file}')) %>% as.data.frame()
  repair_data$pair_prediction <- paste0(repair_data$chainA, '_', repair_data$chainB)
  repair_data <- dplyr::distinct(repair_data, pair_prediction, .keep_all = TRUE)
  
  ground_truth <- read.table(glue('~/grouth_truth/bin{bin_size}/{sample_id}.tsv'),
                             header = TRUE, row.names = NULL)
  
  repair_data <- repair_data[repair_data$chainA %in% ground_truth$Hclone, ]
  allelic_summary$repair_accuracy[i] <- sum(repair_data$pair_prediction %in% ground_truth$pair_name) / nrow(repair_data)
}

allelic_summary <- allelic_summary[!is.na(allelic_summary$repair_accuracy) & !is.na(allelic_summary$acc), ]
allelic_summary <- dplyr::arrange(allelic_summary, bin_size, threshold, cur_sample)
# allelic_summary <- allelic_summary[!allelic_summary$cur_sample %in% c("4266B", "4309B", "B4266", "B4309"),]
new_id <- fread("~/newid-加上标签.csv")
allelic_summary <- dplyr::left_join(allelic_summary, new_id, by = c("cur_sample" = "rawID"))
write.csv(allelic_summary, glue('~/final_version/compare.csv'), quote = FALSE)

comparison_data <- read.csv(glue('~/final_version/compare.csv'))

df_allelic <- comparison_data
df_allelic$value <- comparison_data$acc
df_allelic$method <- "allelic"

df_repair <- comparison_data
df_repair$value <- comparison_data$repair_accuracy
df_repair$method <- "repair"

combined_df <- rbind(df_allelic, df_repair)

for (bin in c(20, 50, 110)) {
  for (thresh in 1) {
    subset_df <- combined_df[combined_df$bin_size == bin & combined_df$threshold == thresh, ]
    acc_allelic <- subset_df$value[subset_df$method == "allelic"]
    acc_repair <- subset_df$value[subset_df$method == "repair"]
    p_val <- round(wilcox.test(acc_allelic, acc_repair, alternative = "two.sided")$p.value, 3)
    
    
    p <- ggbarplot(subset_df, x = "method", y = "value", color = "method",
                   palette = c(allelic="#B03060", repair="#0B7FAB"),
                   add = c("mean_se", "jitter")) +
      theme() +
      labs(x = "", y = "Accuracy", title = glue("p = {p_val}"))
    
    save_plot(p, glue("~/accbar_{bin}_{thresh}"),
              width = 3, height = 4, if_pdf = TRUE)
  }
}
bin <- 110; thresh <- 1;
subset_df <- combined_df[combined_df$bin_size == bin & combined_df$threshold == thresh, ]
dplyr::group_by(subset_df, method) %>% dplyr::summarise(mean = mean(value))

#### Figure 2F ====

# Line Plot for Accuracies across Thresholds
plot_df <- combined_df[combined_df$bin_size %in% c(20, 50, 110), ]

p_line <- ggline(plot_df, x = "threshold", y = "value", color = "method",
                 add = c("mean_se"), palette = c(allelic="#B03060", repair="#0B7FAB"),
                 facet.by = "bin_size") +
  labs(x = "", y = "") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("~/accline"),
          width = 10, height = 5, if_pdf = TRUE)
save_plot(p_line, glue("~/accline"),
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
  
  save_plot(p_bin, glue("~/accline_{bin}"),
            width = 4, height = 4, if_png = TRUE)
  save_plot(p_bin, glue("~/accline_{bin}"),
            width = 4, height = 4, if_pdf = TRUE)
}
bin <- 20;threshold = 1;
subset_df <- combined_df[combined_df$bin_size == bin & combined_df$threshold == thresh, ]
dplyr::group_by(subset_df, method) %>% dplyr::summarise(mean = mean(value))
#### Figure E3 B====
sci_p1  <- readRDS("~/SpatialVDJ_forZenodo/data/breast_cancer/0_integrated/P1/BC_P1_PartIII.rds")
sci_p1_clone <- sci_p1@assays$CLON@counts
saveRDS(sci_p1_clone, "~/sci_p1_clone_count.rds")

sci_p1_clone <- as.data.frame(sci_p1_clone)
sci_p1_clone$cloneid <- rownames(sci_p1_clone)
sci_p1_clone <- sci_p1_clone  %>% 
  tidyr::pivot_longer(!cloneid, names_to = "BinID", values_to = "count")
sci_p1_clone <- sci_p1_clone[sci_p1_clone$count != 0,]

saveRDS(sci_p1_clone, "~/sci_p1_clone.rds")


sci_p2  <- readRDS("~/SpatialVDJ_forZenodo/data/breast_cancer/0_integrated/P2/BC_P2_PartIII.rds")
sci_p2_clone <- sci_p2@assays$CLON@counts
saveRDS(sci_p2_clone, "~/sci_p2_clone_count.rds")

sci_p2_clone <- as.data.frame(sci_p2_clone)
sci_p2_clone$cloneid <- rownames(sci_p2_clone)
sci_p2_clone <- sci_p2_clone  %>% 
  tidyr::pivot_longer(!cloneid, names_to = "BinID", values_to = "count")
sci_p2_clone <- sci_p2_clone[sci_p2_clone$count != 0,]

saveRDS(sci_p2_clone, "~/sci_p2_clone.rds")

sci_p1_clone$SampleID <- "sci_p1"
sci_p2_clone$SampleID <- "sci_p2"

sci_clone <- rbind(sci_p1_clone, sci_p2_clone)
sci_clone$locus <- substr(sci_clone$cloneid, 1, 3)
sci_clone <- sci_clone[!sci_clone$locus %in%  c("TRA", "TRB", "clo"),]
sci_clone$locus2 <- dplyr::if_else(sci_clone$locus == "IGH", "Heavy", "Light")
sci_clone$Bin_size <- "10X"
sci_clone <- dplyr::rename(sci_clone, clone_id = cloneid)


sci_clone_1 <- dplyr::distinct(sci_clone, SampleID, Bin_size, BinID, locus2, clone_id) %>% 
  dplyr::group_by(SampleID, Bin_size, BinID, locus2) %>% 
  dplyr::summarise(clone_count = n()) %>% as.data.frame()

sci_clone_2 <- sci_clone_1 %>% 
  pivot_wider(
    names_from="locus2",
    values_from="clone_count", values_fill = 0)
sci_clone_2 <- subset(sci_clone_2, (Heavy != 0 ) & (Light != 0)) %>% 
  dplyr::group_by(SampleID, Bin_size) %>% 
  dplyr::summarise(mean_H_bin = mean(Heavy),
                   mean_L_bin = mean(Light),
                   n_bin_with_HL = n())

saveRDS(sci_clone_1, "~/sci_summary_1.rds")
saveRDS(sci_clone_2, "~/sci_summary_2.rds")




pacbio_merge <- readRDS("~/45sample_bcr0619.rds")
bin_list <- c(20,30,40,50,60,70,80,90,100,110)
select_col <- paste("Bin", bin_list, sep = "")
summary_1 <- list()
summary_2 <- list()

st_sc_shared_sample <- read.table("~/sc_st_share_sample1.txt", header = F)
st_sc_shared_sample <- st_sc_shared_sample$V1
for(cur_sample in st_sc_shared_sample){
  cur_result <- pacbio_merge[[cur_sample]]  %>% as.data.frame()
  cur_result <- cur_result[, c(select_col, "clone_id", "locus")]
  cur_result$clone_id <- paste(cur_sample, cur_result$clone_id, sep = "_")
  nrow(cur_result)
  group_bin <- cur_result %>% 
    tidyr::pivot_longer(cols = starts_with("Bin"),
                        names_to = "Bin_size",
                        values_to = "BinID")
  group_bin$locus2 <- dplyr::if_else(group_bin$locus == "IGH", "Heavy", "Light")
  
  group_bin_1 <- dplyr::distinct(group_bin, Bin_size, BinID, locus2, clone_id) %>% 
    dplyr::group_by(Bin_size, BinID, locus2) %>% 
    dplyr::summarise(clone_count = n()) %>% as.data.frame()
  group_bin_1$SampleID <- cur_sample
  
  group_bin_2 <- group_bin_1 %>% 
    pivot_wider(
      names_from="locus2",
      values_from="clone_count", values_fill = 0)
  group_bin_2 <- subset(group_bin_2, (Heavy != 0 ) & (Light != 0)) %>% 
    dplyr::group_by(Bin_size) %>% 
    dplyr::summarise(mean_H_bin = mean(Heavy),
                     mean_L_bin = mean(Light),
                     n_bin_with_HL = n())
  group_bin_2$SampleID <- cur_sample

  summary_1[[cur_sample]] <- group_bin_1
  
  summary_2[[cur_sample]] <- group_bin_2 %>% as.data.frame()
  
}

clone_bin_summary_1 <- dplyr::bind_rows(summary_1)
clone_bin_summary_2 <- dplyr::bind_rows(summary_2)

clone_bin_summary_1$Bin_size <- factor(clone_bin_summary_1$Bin_size, select_col)
clone_bin_summary_2$Bin_size <- factor(clone_bin_summary_2$Bin_size, select_col)

saveRDS(clone_bin_summary_1, "~/clone_bin_summary_1.rds")
saveRDS(clone_bin_summary_2, "~/clone_bin_summary_2.rds")

clone_bin_summary <- readRDS("~/clone_bin_summary_1.rds")
clone_bin_summary$Bin_size <- factor(clone_bin_summary$Bin_size, select_col)
clone_bin_summary <- clone_bin_summary[clone_bin_summary$SampleID %in% st_sc_shared_sample,]


clone_bin_summary_mean <- clone_bin_summary %>% 
  dplyr::group_by(SampleID, locus2, Bin_size) %>% 
  dplyr::summarise(mean_count = sum(clone_count)/n())

p <- ggviolin(clone_bin_summary[clone_bin_summary$locus2 == "Heavy",], "Bin_size", "clone_count", color = "Bin_size",
              palette = mycolor_merge) + 
  labs(x = "", y = "Couy")  + 
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
save_plot(p, glue("~/bin_summary_H"),
          width = 5, height = 5, if_pdf = T)

p <- ggviolin(clone_bin_summary_mean[clone_bin_summary_mean$locus2 == "Heavy",], 
              "Bin_size", "mean_count", color = "Bin_size",
              palette = mycolor_merge) + 
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
save_plot(p, glue("~/bin_summary_H_2"),
          width = 5, height = 5, if_pdf = T)

p <- ggboxplot(clone_bin_summary_mean[clone_bin_summary_mean$locus2 == "Heavy",], "Bin_size", "mean_count", color = "Bin_size",
               palette = mycolor_merge) + 
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
save_plot(p, glue("~/bin_summary_H_3"),
          width = 5, height = 5, if_png = T)

clone_bin_summary_2 <- readRDS("~/clone_bin_summary_2.rds")
sci_summary_2 <- readRDS("~/sci_summary_2.rds")
sci_summary_2 <- sci_summary_2[, colnames(clone_bin_summary_2)]
clone_bin_summary_2 <- rbind(clone_bin_summary_2, sci_summary_2)

clone_bin_summary_2_H <- clone_bin_summary_2
clone_bin_summary_2_H$type <- "Heavy"
clone_bin_summary_2_H$value <- clone_bin_summary_2_H$mean_H_bin
clone_bin_summary_2_L <- clone_bin_summary_2
clone_bin_summary_2_L$type <- "Light"
clone_bin_summary_2_L$value <- clone_bin_summary_2_L$mean_L_bin

clone_bin_summary_2_L %>% dplyr::group_by(Bin_size) %>% dplyr::summarise(mean = mean(mean_L_bin))

p <- ggbarplot(clone_bin_summary_2_H,
               x = "Bin_size", y = "value", color = "black",
               # palette = c(Our="#B03060", Repair="#0B7FAB"), 
               add = c("mean_se", "jitter"),
               size = 0.5)+ 
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  labs(x = "", y = "n clone")
save_plot(p, glue("~/bin_summary_H"),
          width = 6, height = 6, if_png = T)
save_plot(p, glue("~/bin_summary_H"),
          width = 6, height = 6, if_pdf = T)

clone_bin_summary_2_H %>% dplyr::group_by(Bin_size) %>% dplyr::summarise(mean = mean(value))

#### Figure E3 H 2879T====
library(readxl)


all_bcr <- readRDS("merge4.rds")
tb <- readRDS("merged_bcr.rds")
tb <- tb[, c("sequence_id", "molecule", "isotype", "clone_family", "bin_id", "Bin50", "TLS_Final", "TLS_raw", "Bin_Region",    
             "cell_annotation", "cell_subtype", "TLS_maturity", "Type")]
all_bcr <- merge(all_bcr, tb, by = "sequence_id")

all_bcr$mab_name <- paste(all_bcr$locus, all_bcr$clone_id, sep = "_")
all_bcr$BinID <- all_bcr$Bin50


if_tls_contour = "T"

cur_sample <- "2879T"
outdir <- glue("~")
dir.create(outdir, recursive = T)

tb_denovo <- fread(glue("~/{cur_sample}_final_e100.tsv"))
tb_denovo <- tb_denovo[is.na(tb_denovo$Lclone_t),]

all_bcr_tmp <- all_bcr[all_bcr$sample == cur_sample,]

obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <- read.table(glue("~/{cur_sample}.txt"), header = T) %>% as.data.frame()
obj$BinID <- substr(obj$BinID, 7, nchar(obj$BinID))
rownames(obj@meta.data) <- obj$BinID

obj@meta.data$bin_id <- obj@meta.data[,1] 
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region %in% c("Tumor_side_of_Margin_area", "Tumor_tissue")] <- "Tumor"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region %in% c("Paratumor_side_of_Margin_area", "Paratumor_tissue")] <- "Paratumor"
obj@meta.data$TLS <- obj@meta.data$TLS_raw
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Conforming"] <- "conforming"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Deviating"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "NotMature"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Rare_Naive_GCB"] <- "STRH05B_NA"

obj@meta.data[!obj@meta.data$TLS_maturity %in%  c("Mature", "deviating","conforming"), ]$TLS <- "STRH05B_NA"

pair1 <- list(IGH = tb_denovo$Hclone[1],
              IGL = tb_denovo$Lclone[1])

IGH_bin <- all_bcr_tmp$Bin50[all_bcr_tmp$mab_name == pair1$IGH]
IGL_bin <- all_bcr_tmp$Bin50[all_bcr_tmp$mab_name == pair1$IGL]
IGHL_bin <- intersect(IGH_bin, IGL_bin)


## 重链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGH_bin] <- pair1$IGH

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)

print(unique(obj_new@meta.data$mab_name))
print(unique(obj_new@meta.data$Bin_Region))

# obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name %in% c("Tumor","Invasive_zone"), "Tumor", obj_new@meta.data$mab_name)
# obj_new$TLS2 <- obj_new$TLS
obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 plot_item = "mab_name",
                 if_tls_contour = "T",
                 plot_order = c(pair1$IGH,"TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#C03724", "#9CCCE8", "#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGH"))

## 轻链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGL_bin] <- pair1$IGL

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)
# obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name %in% c("Tumor","Invasive_zone"), "Tumor", obj_new@meta.data$mab_name)

obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 plot_item = "mab_name",
                 plot_order = c(pair1$IGL,"TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#669BBC", "#9CCCE8","#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGL"))

## 重轻链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGL_bin] <- pair1$IGL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGH_bin] <- pair1$IGH
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGHL_bin] <- paste(pair1$IGH, pair1$IGL, sep = "_")

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)

obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 plot_item = "mab_name",
                 plot_order = c(pair1$IGH, pair1$IGL, paste(pair1$IGH, pair1$IGL, sep = "_"), "TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#C03724","#669BBC", "#3A8339","#9CCCE8","#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGHL"))




#### Figure E3 H 2907T ====

cur_sample <- "2907T"
outdir <- glue("~")
dir.create(outdir, recursive = T)

tb_denovo <- fread(glue("~/{cur_sample}_final_e100.tsv"))
tb_denovo <- tb_denovo[is.na(tb_denovo$Lclone_t),]

all_bcr_tmp <- all_bcr[all_bcr$sample == cur_sample,]

obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <- read.table(glue("~/{cur_sample}.txt"), header = T) %>% as.data.frame()
obj$BinID <- substr(obj$BinID, 7, nchar(obj$BinID))
rownames(obj@meta.data) <- obj$BinID

obj@meta.data$bin_id <- obj@meta.data[,1] 
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region %in% c("Tumor_side_of_Margin_area", "Tumor_tissue")] <- "Tumor"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region %in% c("Paratumor_side_of_Margin_area", "Paratumor_tissue")] <- "Paratumor"
obj@meta.data$TLS <- obj@meta.data$TLS_raw
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Conforming"] <- "conforming"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Deviating"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "NotMature"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Rare_Naive_GCB"] <- "STRH05B_NA"

obj@meta.data[!obj@meta.data$TLS_maturity %in%  c("Mature", "deviating","conforming"), ]$TLS <- "STRH05B_NA"

pair1 <- list(IGH = tb_denovo$Hclone[1],
              IGL = tb_denovo$Lclone[1])

IGH_bin <- all_bcr_tmp$Bin50[all_bcr_tmp$mab_name == pair1$IGH]
IGL_bin <- all_bcr_tmp$Bin50[all_bcr_tmp$mab_name == pair1$IGL]
IGHL_bin <- intersect(IGH_bin, IGL_bin)



## 重链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGH_bin] <- pair1$IGH

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)

print(unique(obj_new@meta.data$mab_name))
print(unique(obj_new@meta.data$Bin_Region))


obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 plot_item = "mab_name",
                 if_tls_contour = "T",
                 plot_order = c(pair1$IGH,"TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#C03724", "#9CCCE8", "#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGH"))

## 轻链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGL_bin] <- pair1$IGL

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)

obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 # prefix = glue("/Volumes/CLS/mess/三代投稿/Fig/Fig5/spatial_cloneid/{cur_sample}"),
                 plot_item = "mab_name",
                 plot_order = c(pair1$IGL,"TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#669BBC", "#9CCCE8","#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGL"))

## 重轻链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGL_bin] <- pair1$IGL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGH_bin] <- pair1$IGH
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGHL_bin] <- paste(pair1$IGH, pair1$IGL, sep = "_")

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)

obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 # prefix = glue("/Volumes/CLS/mess/三代投稿/Fig/Fig5/spatial_cloneid/{cur_sample}"),
                 plot_item = "mab_name",
                 plot_order = c(pair1$IGH, pair1$IGL, paste(pair1$IGH, pair1$IGL, sep = "_"), "TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#C03724","#669BBC", "#3A8339","#9CCCE8","#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGHL"))





####  Figure E3 H RH05B  ====
cur_sample <- "RH05B"
outdir <- glue("~")
dir.create(outdir, recursive = T)

tb_denovo <- fread(glue("~/{cur_sample}_final_e100.tsv"))
tb_denovo <- tb_denovo[is.na(tb_denovo$Lclone_t),]

all_bcr_tmp <- all_bcr[all_bcr$sample == cur_sample,]

obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <- read.table(glue("~/{cur_sample}.txt"), header = T) %>% as.data.frame()
obj$BinID <- substr(obj$BinID, 7, nchar(obj$BinID))
rownames(obj@meta.data) <- obj$BinID

obj@meta.data$bin_id <- obj@meta.data[,1] 
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region %in% c("Tumor_side_of_Margin_area", "Tumor_tissue")] <- "Tumor"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region %in% c("Paratumor_side_of_Margin_area", "Paratumor_tissue")] <- "Paratumor"
obj@meta.data$TLS <- obj@meta.data$TLS_raw
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Conforming"] <- "conforming"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Deviating"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "NotMature"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Rare_Naive_GCB"] <- "STRH05B_NA"

obj@meta.data[!obj@meta.data$TLS_maturity %in%  c("Mature", "deviating","conforming"), ]$TLS <- "STRH05B_NA"

pair1 <- list(IGH = tb_denovo$Hclone[1],
              IGL = tb_denovo$Lclone[1])

IGH_bin <- all_bcr_tmp$Bin50[all_bcr_tmp$mab_name == pair1$IGH]
IGL_bin <- all_bcr_tmp$Bin50[all_bcr_tmp$mab_name == pair1$IGL]
IGHL_bin <- intersect(IGH_bin, IGL_bin)




## 重链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGH_bin] <- pair1$IGH

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)

print(unique(obj_new@meta.data$mab_name))
print(unique(obj_new@meta.data$Bin_Region))


obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 plot_item = "mab_name",
                 if_tls_contour = "T",
                 plot_order = c(pair1$IGH,"TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#C03724", "#9CCCE8", "#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGH"))

## 轻链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGL_bin] <- pair1$IGL

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)
# obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name %in% c("Tumor","Invasive_zone"), "Tumor", obj_new@meta.data$mab_name)

obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 # prefix = glue("/Volumes/CLS/mess/三代投稿/Fig/Fig5/spatial_cloneid/{cur_sample}"),
                 plot_item = "mab_name",
                 plot_order = c(pair1$IGL,"TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#669BBC", "#9CCCE8","#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGL"))

## 重轻链
obj_new <- obj
obj_new@meta.data$mab_name <- NULL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGL_bin] <- pair1$IGL
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGH_bin] <- pair1$IGH
obj_new@meta.data$mab_name[obj_new@meta.data$BinID %in% IGHL_bin] <- paste(pair1$IGH, pair1$IGL, sep = "_")

obj_new@meta.data$bin_id <- rownames(obj_new@meta.data)
obj_new@meta.data$mab_name <- as.character(obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(is.na(obj_new@meta.data$mab_name), "other", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$TLS != "STRH05B_NA" & obj_new@meta.data$mab_name == "other", "TLS", obj_new@meta.data$mab_name)
obj_new@meta.data$mab_name <- ifelse(obj_new@meta.data$mab_name == "other", obj_new@meta.data$Bin_Region, obj_new@meta.data$mab_name)

obj_new$TLS <- obj_new$Bin_Region
obj_new$TLS[obj_new$TLS != "Invasive_zone"] <- "_NA"
SpatialDimPython(obj = obj_new,
                 prefix = glue("{outdir}"),
                 # prefix = glue("/Volumes/CLS/mess/三代投稿/Fig/Fig5/spatial_cloneid/{cur_sample}"),
                 plot_item = "mab_name",
                 plot_order = c(pair1$IGH, pair1$IGL, paste(pair1$IGH, pair1$IGL, sep = "_"), "TLS", "Tumor", "Invasive_zone","Paratumor"),
                 tmp_color = c("#C03724","#669BBC", "#3A8339","#9CCCE8","#E5E7E6","#F3E1B9","#F5F4EA"),
                 name = glue("{cur_sample}_pair1_IGHL"))


#### Figure E3 I ====

library(dplyr)
library(ggpubr)
library(glue)

# 1. 读取原始数据
data_raw <- read.csv('~/final_version/compare.csv')

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
output_path <- "~/eff_line_2"
# save_plot(p, glue("{output_path}"), width = 4, height = 4, if_png = TRUE)
save_plot(p, glue("{output_path}"), width = 4, height = 4, if_pdf = TRUE)


library(extrafont)
# 确保你已经导入字体
font_import() # 只需第一次运行
loadfonts()

# 设置字体
par(family = "Arial Unicode MS")  # 或其他支持μ的字体

comparison_data <- read.csv(glue('~/final_version/compare.csv'))

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
  
  save_plot(p_line, glue("~/31_sample_acc/{current_sample}"),
            width = 10, height = 5, if_pdf = TRUE)

}

p_line <- ggline(plot_df, x = "threshold", y = "value", color = "method",
                 palette = c(allelic="#B03060", repair="#0B7FAB"),
                 facet.by = c( "newID", "resolution")) +
  labs(x = "", y = "") +
  theme(plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5))

save_plot(p_line, glue("~/31_sample_acc/merge"),
          width = 10, height = 40, if_pdf = TRUE)

#### Figure E3 K ====


ours <- read.csv(glue('~/compare_4.csv'))
ours <- ours[!is.na(ours$repair_acc),]
ours <- ours[!is.na(ours$acc),]
df_1 <- ours
df_1$value <- ours$acc
df_1$group <- "allelic"
df_2 <- ours
df_2$value <- ours$repair_acc
df_2$group <- "repair"
df_p <- rbind(df_1, df_2)
bin_list <- c(20, 50, 110)
df <- df_p[df_p$bin_size %in% bin_list, ]
df <- df[df$threshold == 1,]
df$x <- paste(df$bin_size, df$type, sep = "_")
df$bin_size <- as.character(df$bin_size)
summary(df)
df_p <- dplyr::group_by(df,cur_sample, group, bin_size) %>% 
  dplyr::summarise(acc_mean = mean(value))
df_p$bin_size[df_p$cur_sample == "B4266"] <- "110-10x"
df_p$bin_size[df_p$cur_sample == "B4309"] <- "110-10x"
df_p$cur_sample[df_p$cur_sample == "B4266"] <- "4266B"
df_p$cur_sample[df_p$cur_sample == "B4309"] <- "4309B"

df_p$bin_size <- factor(df_p$bin_size, levels = c("20", "50", "110","110-10x"))
# df_p$group <- paste(df_p$cur_sample, df_p$group, sep = "_")
p <- ggline(df_p, "bin_size", "acc_mean", color = "group",
            palette = c(allelic="#B03060", repair="#0B7FAB")) +
  labs(x = "", y = "accuracy")  +
  facet_wrap(~cur_sample, scales = "free_y") +
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold",hjust = 0.5),
            legend.position = "right")

save_plot(p, glue("~/accline-v2"), 
          width = 7, height = 5, if_pdf = T, if_png = T)


aov.manu <- aov(value ~ group, data=df[df$bin_size == 20,])
summary(aov.manu)
aov.manu <- aov(value ~ group, data=df[df$bin_size == 50,])
summary(aov.manu)
aov.manu <- aov(value ~ group, data=df[df$bin_size == 110,])
summary(aov.manu)



#### Figure E3 L p1====


pred_all <- fread("~/p1_final_e100.tsv")
df_clone_num <- readRDS("~/df_clone_num_p1.rds")
df_clone_num_2 <- df_clone_num %>% dplyr::rename(clone_num_L = clone_num,
                                                 bin_num_L = bin_num) %>% 
  dplyr::filter(locus!="IGH") %>% 
  dplyr::select(-locus)
test_summary <- fread("~/p1_test_summary.tsv")

all_test <- read.table(glue("~/p1_IG_groundtruth_repair.tsv"),
                       sep = "\t", header = TRUE, row.names = 1)
all_test$pair_name <- paste(all_test$chainA, all_test$chainB, sep = "_")
pred_test <- fread("~/p1_test_e100.tsv")

compare_df <- list()
compare_df_2 <- list()
thre_list <- c(1, 5, 10, 15, 20)
for(cur_thre in thre_list){
  repair_out <- fread(glue("~/p1_{cur_thre}.tsv"))
  repair_out <- repair_out %>% 
    dplyr::left_join(df_clone_num[df_clone_num$locus == "IGH",], by = c("chainA" = "clone_name")) %>% 
    dplyr::select(-locus) %>% 
    dplyr::left_join(df_clone_num_2, by = c("chainB" = "clone_name"))
  repair_out$pair_name <- paste(repair_out$chainA, repair_out$chainB, sep = "_")
  repair_test <- repair_out[(repair_out$chainA %in% all_test$chainA) | (repair_out$chainB %in% all_test$chainB),]
  
  r <- repair_test[repair_test$bin_num >= cur_thre,]
  rt <- r[r$pair_name %in% all_test$pair_name,]
  p <- pred_test[pred_test$bin_num >= cur_thre,]
  pt <- p[p$pair_name %in% all_test$pair_name,]
  compare_df[[glue("c_{cur_thre}")]] <- data.frame(
    repair_n = nrow(r),
    repair_nt = nrow(rt),
    repair_e = nrow(r)-nrow(rt),
    repair_acc = nrow(rt)/nrow(r),
    our_n = nrow(p),
    our_nt = nrow(pt),
    our_e = nrow(p)-nrow(pt),
    our_acc = nrow(pt)/nrow(p),
    threshold = as.character(cur_thre)
  )
  repair_n = nrow(r)
  repair_nt = nrow(rt)
  repair_e = nrow(r)-nrow(rt)
  repair_acc = nrow(rt)/nrow(r)
  our_n = nrow(p)
  our_nt = nrow(pt)
  our_e = nrow(p)-nrow(pt)
  our_acc = nrow(pt)/nrow(p)
  threshold = as.character(cur_thre)
  df_2_repair <- data.frame(
    group = c("n", "true", "error", "acc"),
    value = c(repair_n, repair_nt, repair_e, repair_acc)
  )
  df_2_repair$method <- "repair"
  df_2_our <- data.frame(
    group = c("n", "true", "error", "acc"),
    value = c(our_n, our_nt, our_e, our_acc)
  )
  df_2_our$method <- "allelic"
  df_2 <- rbind(df_2_our, df_2_repair)
  df_2$threshold <- cur_thre
  compare_df_2[[glue("c_{cur_thre}")]] <- df_2
}
compare_df <- dplyr::bind_rows(compare_df)
compare_df_2 <- dplyr::bind_rows(compare_df_2)

p <- ggline(compare_df_2[compare_df_2$group %in% c("acc"),], "threshold", "value", 
            color = "method", palette = c(allelic="#B03060", repair="#0B7FAB")) +
  labs(x = "", y = "")  +
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold",hjust = 0.5))
save_plot(p, glue("~/p1_accline"), 
          width = 4, height = 4,if_pdf = T)




### Figure E3 L p2====

pred_all <- fread("~/p2_final_e100.tsv")
df_clone_num <- readRDS("~/df_clone_num_p2.rds")
df_clone_num_2 <- df_clone_num %>% dplyr::rename(clone_num_L = clone_num,
                                                 bin_num_L = bin_num) %>% 
  dplyr::filter(locus!="IGH") %>% 
  dplyr::select(-locus)
test_summary <- fread("~/p2_test_summary.tsv")

all_test <- read.table(glue("~/p2_IG_groundtruth_repair.tsv"),
                       sep = "\t", header = TRUE, row.names = 1)
all_test$pair_name <- paste(all_test$chainA, all_test$chainB, sep = "_")
pred_test <- fread("~/p2_test_e100.tsv")


compare_df <- list()
compare_df_2 <- list()
thre_list <- c(1,  5, 10, 15, 20)
for(cur_thre in thre_list){
  repair_out <- fread(glue("~/p2_{cur_thre}.tsv"))
  repair_out <- repair_out %>% 
    dplyr::left_join(df_clone_num[df_clone_num$locus == "IGH",], by = c("chainA" = "clone_name")) %>% 
    dplyr::select(-locus) %>% 
    dplyr::left_join(df_clone_num_2, by = c("chainB" = "clone_name"))
  repair_out$pair_name <- paste(repair_out$chainA, repair_out$chainB, sep = "_")
  repair_test <- repair_out[(repair_out$chainA %in% all_test$chainA) | (repair_out$chainB %in% all_test$chainB),]
  r <- repair_test[repair_test$bin_num >= cur_thre,]
  rt <- r[r$pair_name %in% all_test$pair_name,]
  p <- pred_test[pred_test$bin_num >= cur_thre,]
  pt <- p[p$pair_name %in% all_test$pair_name,]
  compare_df[[glue("c_{cur_thre}")]] <- data.frame(
    repair_n = nrow(r),
    repair_nt = nrow(rt),
    repair_e = nrow(r)-nrow(rt),
    repair_acc = nrow(rt)/nrow(r),
    our_n = nrow(p),
    our_nt = nrow(pt),
    our_e = nrow(p)-nrow(pt),
    our_acc = nrow(pt)/nrow(p),
    threshold = as.character(cur_thre)
  )
  repair_n = nrow(r)
  repair_nt = nrow(rt)
  repair_e = nrow(r)-nrow(rt)
  repair_acc = nrow(rt)/nrow(r)
  our_n = nrow(p)
  our_nt = nrow(pt)
  our_e = nrow(p)-nrow(pt)
  our_acc = nrow(pt)/nrow(p)
  threshold = as.character(cur_thre)
  df_2_repair <- data.frame(
    group = c("n", "true", "error", "acc"),
    value = c(repair_n, repair_nt, repair_e, repair_acc)
  )
  df_2_repair$method <- "repair"
  df_2_our <- data.frame(
    group = c("n", "true", "error", "acc"),
    value = c(our_n, our_nt, our_e, our_acc)
  )
  df_2_our$method <- "allelic"
  df_2 <- rbind(df_2_our, df_2_repair)
  df_2$threshold <- cur_thre
  compare_df_2[[glue("c_{cur_thre}")]] <- df_2
}
compare_df <- dplyr::bind_rows(compare_df)
compare_df_2 <- dplyr::bind_rows(compare_df_2)
save_table(compare_df_2, glue("~/compare_df_2"))

p <- ggline(compare_df_2[compare_df_2$group %in% c("acc"),], "threshold", "value", 
            color = "method", palette = c(allelic="#B03060", repair="#0B7FAB")) +
  labs(x = "", y = "")  +
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold",hjust = 0.5))
save_plot(p, glue("~/p2_accline"), 
          width = 4, height = 4,if_pdf = T)


#### data summary =====

fils <- list.files("~/Result/DedupBin")
result_list <- list()

for (i in fils) {
  tb <- read.table(glue("~/Result/DedupBin/{i}"), header = TRUE)

  tb2 <- table(tb$binid) %>% as.data.frame()
  colnames(tb2) <- c("binid", "Freq")


  result_list[[i]] <- tb2
}

combined_result <- bind_rows(result_list, .id = "filename")
combined_result$filename <- substr(combined_result$filename,1,5)

tb_info <- fread("~/newid.csv")
combined_result <- merge(combined_result, tb_info, by.x = "filename", by.y = "rawID", all.x = T)

combined_result$newID <- factor(combined_result$newID, levels = tb_info$newID)
saveRDS(combined_result, "~/combined_result.rds")
combined_result <- readRDS("~/combined_result.rds")


combined_result_mean <- dplyr::group_by(combined_result, newID) %>%
  dplyr::summarise(freq_mean  = mean(Freq))
save_table(combined_result_mean, "~/UMI_count_per_Bin50")

combined_result_mean <- fread("~/UMI_count_per_Bin50.txt")

merge_bcr <- readRDS("~/merged_bcr.rds")
merge_bcr <- merge(merge_bcr, tb_info, by.x = "sample", by.y = "rawID", all.x = T)
bcr_count <- table(merge_bcr$newID) %>% as.data.frame()
colnames(bcr_count) <- c("newID", "bcr_num")
bcr_count <- merge_bcr %>%
  dplyr::group_by(newID) %>%
  summarise(
    n_IgH = sum(locus == "IGH", na.rm = TRUE),
    n_IgL = sum(locus != "IGH", na.rm = TRUE),
    bcr_num = n()
  )
merge_bcr_clone <- merge_bcr %>%
  dplyr::distinct(newID, clone_id, .keep_all = T) %>%
  dplyr::group_by(newID) %>%
  summarise(
    n_IgH_clone = sum(locus == "IGH", na.rm = TRUE),
    n_IgL_clone = sum(locus != "IGH", na.rm = TRUE),
    clone_num = n()
  )
bcr_count <- dplyr::left_join(bcr_count, merge_bcr_clone)
save_table(bcr_count,  "~/45_sample_BCR_info")
bcr_count <- fread("~/45_sample_BCR_info.txt")
merge_anno <- fread("~/Bin50/Bin50_classic_celltype_metadata_merge.txt")
merge_anno <- merge_anno[merge_anno$SampleID %in% merge_bcr$sample,]
merge_anno <- merge(merge_anno, tb_info, by.x = "SampleID", by.y = "rawID", all.x = T)
cell_prop <- merge_anno %>%
  dplyr::group_by(newID) %>%
  dplyr::summarize(
    total = n(),
    count_B_Plasma = sum(SpotLight_Anno %in% c("B", "Plasma")),
    proportion = count_B_Plasma / total
  )
save_table(cell_prop, "~/cell_prop")
cell_prop <- fread("~/cell_prop.txt")

comparison_data <- read.csv(glue('~/final_version/compare.csv'))
cell_prop <- fread("~/cell_prop.txt")
combined_result <- readRDS("~/combined_result.rds")
combined_result_mean <- fread("~/UMI_count_per_Bin50.txt")
bcr_count <- fread("~/45_sample_BCR_info.txt")
comparison_data$n_Shared <- sapply(comparison_data$cur_sample, function(sample_id) {
  file_path <- file.path(
    "~/grouth_truth/bin50",
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

merge_bcr <- readRDS("~/merged_bcr.rds")
merge_bcr <- merge(merge_bcr, tb_info, by.x = "sample", by.y = "rawID", all.x = T)

summary_bind <- list()

for(current_sample in comparison_data$cur_sample %>% unique()){
  grount_truth <- fread(glue("~/grouth_truth/bin50/{current_sample}.tsv"))
  s_df <- comparison_data[comparison_data$cur_sample ==current_sample,]
  new_id <-  s_df$newID[1]
  for(current_bin_size in s_df$bin_size %>% unique()){
    s_df_1 <- s_df[s_df$bin_size ==current_bin_size,]
    for(current_clone_size in s_df_1$clone_size %>% unique()){
      test_out <- fread(glue("~/final_version/bin{current_bin_size}/{current_sample}_test_e100.tsv"))
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
           file = "~/summary_shared.xlsx")
summary_shared <- read_xlsx("~/summary_shared.xlsx")

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
           file = "~/31_sample.xlsx")


