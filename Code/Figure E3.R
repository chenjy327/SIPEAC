#### Figure 1B ====
bcr_merge <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr.rds")

totcounts.fig.45 <- list()
receptor.sum.45 <- list()
for(cur_sample in unique(bcr_merge$sample)){
  hcc_st <- readRDS(glue('/Users/mac/Desktop/01finalbackup/60TLS_sample_for_qianwen_backup/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
  hcc_st@meta.data <-read.table(glue("/Users/mac/Desktop/01finalbackup/Bin50/Bin50经典版注释metadata/{cur_sample}.txt"),
                                sep = "\t",
                                header = T,
                                row.names = 1) %>% as.data.frame()
  
  all_constant <- c("TRAC", "TRBC1", "IGHG1", "IGHG2", "IGHG3",
                    "IGHG4", "IGHM", "IGHA1","IGHA2", "IGHD", "IGKC", "IGLC")
  RNA1 <- hcc_st@assays$Spatial@counts
  totcounts.allgenes <- as.data.frame(rowSums(RNA1))
  totcounts.constant <- as.data.frame(rowSums(RNA1[rownames(RNA1) %in% all_constant,]))
  # Add MS4A1 to chain df
  totcounts.fig <- as.data.frame(rbind(totcounts.constant, "MS4A1" = totcounts.allgenes["MS4A1",1]))
  # And CD3E
  totcounts.fig <- as.data.frame(rbind(totcounts.fig, "CD3E" = totcounts.allgenes["CD3E",1]))
  # Add gene as column and set colnames
  totcounts.fig <- data.frame(gene=rownames(totcounts.fig),
                              counts=totcounts.fig[,1])
  totcounts.fig$sample <- cur_sample
  totcounts.fig.45[[cur_sample]] <- totcounts.fig
  
  
  IG_genes <- grep(pattern = "^IG", x=totcounts.fig$gene)
  IG_genes_sum <- sum(totcounts.fig$counts[IG_genes])
  receptor_sum_IG <- data.frame(type = "IG", counts = IG_genes_sum, ratio = IG_genes_sum/sum(totcounts.allgenes))
  TR_genes <- grep(pattern = "^TR", x=totcounts.fig$gene)
  TR_genes_sum <- sum(totcounts.fig$counts[TR_genes])
  receptor_sum_TR <- data.frame(type = "TR", counts = TR_genes_sum, ratio = TR_genes_sum/sum(totcounts.allgenes))
  others_gene_sum <- sum(totcounts.allgenes)-IG_genes_sum-TR_genes_sum
  receptor_sum_others <- data.frame(type = "others", counts = others_gene_sum, ratio = others_gene_sum/sum(totcounts.allgenes))
  
  receptor_sum <- do.call(rbind, list(receptor_sum_IG, receptor_sum_TR, receptor_sum_others))
  receptor_sum$sample <- cur_sample
  receptor.sum.45[[cur_sample]] <- receptor_sum
  
}

saveRDS(totcounts.fig.45, "/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/totcounts.fig.45.rds")
saveRDS(receptor.sum.45, "/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/receptor.sum.45.rds")



receptor.sum <- do.call(rbind, receptor.sum.45)
receptor.sum$logcount <- log10(receptor.sum$counts)
receptor.sum$type <- factor(receptor.sum$type, levels = c("IG", "TR", "others"))

mean(receptor.sum[receptor.sum$type == "IG",]$ratio)
mean(receptor.sum[receptor.sum$type == "TR",]$ratio)
mean(receptor.sum[receptor.sum$type == "IG",]$ratio)/mean(receptor.sum[receptor.sum$type == "TR",]$ratio)
receptor.sum$log_ratio <- log10(receptor.sum$ratio)
p <- ggbarplot(receptor.sum,
               x = "type", y = "log_ratio", color = "type",
               palette = c(IG="#FF4040", TR="yellow", others = "black"), 
               add = c("mean_se", "jitter"))+
  labs(x = "", y = "log10(fraction)", title = "") +
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
save_plot(p, glue("/Users/mac/Desktop/Rproject/SIPEAC/Figure1/B/B"),
          width = 3, height = 3, if_png = T, if_pdf = T)



totals <- receptor.sum %>% 
  dplyr::group_by(type) %>% 
  dplyr::summarise(mean = mean(log_ratio))
totals$mean[2]/totals$mean[1]





#### Figure 1C smrt IgH====
bcr_merge <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr.rds")
bcr_merge <- bcr_merge[bcr_merge$locus == "IGH",]

pie_df <- bcr_merge %>% dplyr::group_by(isotype) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE"))
fill_colors <- c("#fdc7d4",
                 "#ef8e89",
                 "#f6ba59",
                 "#9fc0d6",
                 "#aacb89",
                 "#5c9b4a",
                 "#7775b6",
                 "#caabd0",
                 "#225ea8")
names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "isotype", color = "white",
           palette = fill_colors)
save_plot(p, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/C-IgH", if_pdf = T,
          width = 4, height = 4)

print(sum(pie_df$value))
save_table(pie_df, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/C-IgH")

#### Figure 1C smrt IgL====

bcr_merge <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr.rds")
bcr_merge <- bcr_merge[bcr_merge$locus != "IGH",]

pie_df <- bcr_merge %>% dplyr::group_by(locus) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4",
                 "#9fc0d6")
names(fill_colors) <- c("IGK", "IGL")

labs <- paste0(pie_df$locus, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "locus", color = "white",
           palette = fill_colors)
save_plot(p, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/C-IgL", if_pdf = T,
          width = 4, height = 4)

print(sum(pie_df$value))
save_table(pie_df, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/C-IgL")


#### Figure 1C smrt TCR====
tcr_merge <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_tcr.rds")
pie_df <- tcr_merge %>% dplyr::group_by(topChains) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
pie_df$topChains <- factor(pie_df$topChains, levels = rev(c("TRB", "TRA", "TRG", "TRD")))

fill_colors <- c("#fdc7d4",
                 "#ef8e89",
                 "#5c9b4a",
                 "#9fc0d6")
names(fill_colors) <- c("TRB", "TRA", "TRG", "TRD")

labs <- paste0(pie_df$topChains, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "topChains", color = NULL,
           palette = fill_colors)
save_plot(p, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/C-tcr", if_pdf = T,
          width = 4, height = 4)

print(sum(pie_df$value))
save_table(pie_df, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/C-tcr")

#### Figure 1C singlecell IgH ====
sc_bcr <- readRDS('/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr_sc.rds')
mrna_list <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/mrna_tcr_66sample.rds")

cur_sc_bcr <- subset(sc_bcr,sample %in% mrna_list)
cur_sc_bcr <- subset(cur_sc_bcr, !sample %in% c("2931T_B_FS", "2976T_B_FS_2", "2976T_B_FS_3", "2976T_B_FS_1",
                                                "3125_P","3163_T"))
cur_meta <- subset(cur_sc_bcr,locus=='IGH')
cur_meta$isotype <- substr(cur_meta$c_call,1,5)
cur_meta$isotype %>% unique()
cur_meta <- subset(cur_meta,isotype!='')
cur_meta[cur_meta$isotype=='IGHM*',]$isotype <- 'IGHM'
cur_meta[cur_meta$isotype=='IGHD*',]$isotype <- 'IGHD'
table(cur_meta$isotype)

pie_df <- cur_meta %>% dplyr::group_by(isotype) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE"))
fill_colors <- c("#fdc7d4",
                 "#ef8e89",
                 "#f6ba59",
                 "#9fc0d6",
                 "#aacb89",
                 "#5c9b4a",
                 "#7775b6",
                 "#caabd0",
                 "#225ea8")
names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "isotype", color = NULL,
           palette = fill_colors)
save_plot(p, '/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/singcell-IgH', if_pdf = T,
          width = 4, height = 4)
save_table(pie_df,'/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/singcell-IgH')

#### Figure 1C singlecell IgL ====
sc_bcr <- readRDS('/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr_sc.rds')
mrna_list <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/mrna_tcr_66sample.rds")

cur_sc_bcr <- subset(sc_bcr,sample %in% mrna_list)
cur_sc_bcr <- subset(cur_sc_bcr, !sample %in% c("2931T_B_FS", "2976T_B_FS_2", "2976T_B_FS_3", "2976T_B_FS_1",
                                                "3125_P","3163_T"))

cur_meta <- subset(cur_sc_bcr,locus!= "IGH")
pie_df <- cur_meta %>% dplyr::group_by(locus) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4",
                 "#9fc0d6")
names(fill_colors) <- c("IGK", "IGL")

labs <- paste0(pie_df$locus, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "locus", color = NULL,
           palette = fill_colors)
save_plot(p, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/singcell-IgL", if_pdf = T,
          width = 4, height = 4)

print(sum(pie_df$value))
save_table(pie_df, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/singcell-IgL")

#### Figure 1C singlecell Tcr ====
TCR_merge <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/TCR_merge.rds")
TCR_merge <- subset(TCR_merge, !sample %in% c("2931T_B_FS", "2976T_B_FS_2", "2976T_B_FS_3", "2976T_B_FS_1",
                                              "3125_P","3163_T"))
TCR_df <- subset(TCR_merge,isotype %in% c('TRB','TRA'))
pie_df <- TCR_df %>% dplyr::group_by(isotype) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
pie_df$isotype <- factor(pie_df$isotype, levels = c("TRB", "TRA"))
fill_colors <- c("#fdc7d4",
                 "#ef8e89")
names(fill_colors) <- c("TRB", "TRA")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "isotype", color = NULL,
           palette = fill_colors)
save_plot(p, '/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/singcell-tcr', if_pdf = T,
          width = 4, height = 4)
save_table(pie_df, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/C/singcell-tcr")

#### Figure 1D smrt IgH ====
bcr_data <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr.rds")

# 构建 clone_id：sample + clone_id
bcr_data$clone_id <- paste(bcr_data$sample, bcr_data$clone_id, sep = "_")

# 只保留重链（IGH）
bcr_igh <- bcr_data[bcr_data$locus == "IGH", ]

# 计算克隆数
bcr_clones <- countClones(bcr_igh, clone = "clone_id")
bcr_igh <- left_join(bcr_igh, bcr_clones)

# 计算每个 clone_id 出现的 bin_id 数量（扩增度）
bcr_bin_count <- bcr_igh %>%
  distinct(clone_id, bin_id) %>%
  group_by(clone_id) %>%
  summarise(bin_count = n(), .groups = "drop")

bcr_clones <- left_join(bcr_clones, bcr_bin_count)

# 标记扩增类型
bcr_clones <- bcr_clones %>%
  dplyr::mutate(expand_type =  dplyr::case_when(
    bin_count >= 10 ~ "Expanded >= 10",
    bin_count == 1 ~ "Singleton",
    bin_count > 1 & bin_count < 10 ~"Expanded < 10"
  ))

# 画饼图数据
bcr_pie_df <- bcr_clones %>%
  group_by(expand_type) %>%
  summarise(value = n(), .groups = "drop") %>%
  mutate(
    ratio = value / sum(value),
    percent = percent(ratio, accuracy = 0.1),
    label = paste0(expand_type, " (", percent, ")")
  )

# 饼图颜色
fill_colors <- c("#946625", "#edb55a", "#65a3d1")
names(fill_colors) <- c("Expanded >= 10", "Expanded < 10", "Singleton")


p <- ggplot() + geom_arc_bar(data=bcr_pie_df,
                             stat = "pie",
                             aes(x0=0,y0=0,r0=1,r=2,
                                 amount=ratio, fill=expand_type
                                 # explode=c(0.05,0.1,0.05,0.05,
                                 #           0.05,0.05,0.05,0.05,0.05,0.1,0.1)
                             )
) + scale_fill_manual(values = fill_colors)+ theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')
save_plot(p, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/D/smrt-IgH", if_pdf = T,
          width = 6, height = 6)
print(sum(bcr_pie_df$value))
save_table(bcr_pie_df, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/D/smrt-IgH")




#### Figure 1D smrt tcr ====

tcr_data <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_tcr.rds")

# 构建 clone_id：sample + clone_id
tcr_data$clone_id <- paste(tcr_data$sample, tcr_data$cloneId, sep = "_")

# 只保留重链（IGH）
tcr_igh <- tcr_data[tcr_data$topChains == "TRB", ]

# 计算克隆数
tcr_clones <- countClones(tcr_igh, clone = "clone_id")
tcr_igh <- left_join(tcr_igh, tcr_clones)

# 计算每个 clone_id 出现的 bin_id 数量（扩增度）
tcr_bin_count <- tcr_igh %>%
  distinct(clone_id, Bin50) %>%
  group_by(clone_id) %>%
  summarise(bin_count = n(), .groups = "drop")

tcr_clones <- left_join(tcr_clones, tcr_bin_count)

# 标记扩增类型
tcr_clones <- tcr_clones %>%
  dplyr::mutate(expand_type =  dplyr::case_when(
    bin_count >= 10 ~ "Expanded >= 10",
    bin_count == 1 ~ "Singleton",
    bin_count > 1 & bin_count < 10 ~"Expanded < 10"
  ))

# 画饼图数据
tcr_pie_df <- tcr_clones %>%
  group_by(expand_type) %>%
  summarise(value = n(), .groups = "drop") %>%
  mutate(
    ratio = value / sum(value),
    percent = percent(ratio, accuracy = 0.1),
    label = paste0(expand_type, " (", percent, ")")
  )

# 饼图颜色
fill_colors <- c("#946625", "#edb55a", "#65a3d1")
names(fill_colors) <- c("Expanded >= 10", "Expanded < 10", "Singleton")

p <- ggplot() + geom_arc_bar(data=tcr_pie_df,
                             stat = "pie",
                             aes(x0=0,y0=0,r0=1,r=2,
                                 amount=ratio, fill=expand_type
                                 # explode=c(0.05,0.1,0.05,0.05,
                                 #           0.05,0.05,0.05,0.05,0.05,0.1,0.1)
                             )
) + scale_fill_manual(values = fill_colors)+ theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')
save_plot(p, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/D/smrt-tcr", if_pdf = T,
          width = 6, height = 6)
print(sum(tcr_pie_df$value))
save_table(tcr_pie_df, "/Users/mac/Desktop/Rproject/SIPEAC/Figure1/D/smrt-tcr")


#### Figure 1D singcell IgH ====
