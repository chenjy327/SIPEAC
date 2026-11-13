source("Code/figure_helpers.R")

# Extended Data Fig. 4j：IFI30 在肿瘤和癌旁的表达（院内 scRNA、院内 bulk、TCGA-LIHC）
# 输入文件放在 data_dir 下，文件名保持原样
data_dir <- "data/IFI30_expression"
out_dir <- "results/Figure_E4j"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##scRNA IFI30####
scRNA <- readRDS(file.path(data_dir, "HCC_annotated_with_Peri.rds"))
hep <- subset(scRNA, annotation %in% c("Tumor_Hepatocyte"))
genelist <- c("IFI30")
gene_data <- FetchData(hep, vars = c("Loci","orig.ident",genelist))  %>% data.table::setnames(c("Tissue","samplename","exp"))


df_exp <- gene_data %>% group_by(samplename,Tissue) %>%
  summarise(mean_exp = mean(exp))

p1 <- ggplot(df_exp, aes_string(x = "Tissue",y = "mean_exp" ))+
  # geom_boxplot(outliers = FALSE)+
  geom_boxplot(aes(color=Tissue),position = position_dodge(0.8),width=0.5,outlier.colour = "white")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15),shape=21, size = 0.3, stroke=0.8,aes(color=Tissue),alpha = 0.5)+
  # geom_violin(width = 0.8)+
  # geom_boxplot(width=0.1, fill="white",outlier.shape = NA)+
  # scale_fill_manual(values=c(mycolor_new))+ #颜色
  scale_color_manual(values=c(mycolor_new))+ #颜色
  guides(fill=FALSE) +
  stat_summary(fun.data = mean_fun, geom = "text", size=3 , color="black") +
  labs(title = "",
       x = '',
       y = "IFI30 expression") +
  theme_classic(base_line_size = 0.5)+# 主题设置，这个是无线条主题 # 坐标轴的粗细
  theme(plot.title = element_text(size = 10,colour = "black",hjust = 0.5),
        # axis.title.y = element_text(size = 15, color = "black",face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=10, face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, face = "bold"),
        axis.text.x = element_text(size = 10,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 0.95,
                                   angle = 45), #角度
        axis.text.y = element_text(size = 10,  color = "black",face = "bold", vjust = 0.5,
                                   hjust = 0.5, angle = 0),
        panel.border = element_rect(linetype = 'solid', size = 0.5,fill = NA)
  )+
  stat_compare_means(comparisons =  list(c("T","P")),
                     method = 'wilcox.test',hide.ns = TRUE,label = "p.format",bracket.size = 0.1,size = 2)
p1
save_plot(p1, file.path(out_dir, "boxplot_scRNA_hep"),if_pdf = T,width = 4,height = 4)


##本中心bulk数据####
df_sample <- read.csv(file.path(data_dir, "Group_1104.csv"),header = T,row.names = 1)
df_exp <- read.csv(file.path(data_dir, "rpkm_1104.csv"),header = T,row.names = 1)

IFI30_exp <- df_exp["IFI30",] %>% t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(df_sample)
IFI30_exp$IFI30 <- log2(IFI30_exp$IFI30+1)
IFI30_exp$group <- ifelse(IFI30_exp$group == "L", "P",IFI30_exp$group)
IFI30_exp <- filter(IFI30_exp, group %in% c("P","T"))

p3 <- ggplot(IFI30_exp, aes_string(x = "group",y = "IFI30" ))+
  # geom_boxplot(outliers = FALSE)+
  geom_boxplot(aes(color=group),position = position_dodge(0.8),width=0.5,outlier.colour = "white")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15),shape=21, size = 0.3, stroke=0.8,aes(color=group),alpha = 0.5)+
  # geom_violin(width = 0.8)+
  # geom_boxplot(width=0.1, fill="white",outlier.shape = NA)+
  # scale_fill_manual(values=c(mycolor_new))+ #颜色
  scale_color_manual(values=c(mycolor_new))+ #颜色
  guides(fill=FALSE) +
  stat_summary(fun.data = mean_fun, geom = "text", size=3 , color="black") +
  labs(title = "",
       x = '',
       y = "IFI30 expression") +
  theme_classic(base_line_size = 0.5)+# 主题设置，这个是无线条主题 # 坐标轴的粗细
  theme(plot.title = element_text(size = 10,colour = "black",hjust = 0.5),
        # axis.title.y = element_text(size = 15, color = "black",face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=10, face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, face = "bold"),
        axis.text.x = element_text(size = 10,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 0.95,
                                   angle = 45), #角度
        axis.text.y = element_text(size = 10,  color = "black",face = "bold", vjust = 0.5,
                                   hjust = 0.5, angle = 0),
        panel.border = element_rect(linetype = 'solid', size = 0.5,fill = NA)
  )+
  stat_compare_means(comparisons =  list(c("T","P")),
                     method = 'wilcox.test',hide.ns = TRUE,label = "p.format",bracket.size = 0.1,size = 2)
p3
save_plot(p3, file.path(out_dir, "boxplot_bulk"),if_pdf = T,width = 4,height = 4)
write.csv(IFI30_exp, file.path(out_dir, "bulk_IFI30.csv"))

##TCGA####
df_sample <- read.csv(file.path(data_dir, "TCGA_group.csv"),header = T)
df_exp <- read.csv(file.path(data_dir, "TCGA_LIHC_count_symbol.csv"),header = T,row.names = 1)

IFI30_exp <- df_exp["IFI30",] %>% t() %>% as.data.frame() %>%
  rownames_to_column(var = "samplename") %>%
  left_join(df_sample)
IFI30_exp$IFI30 <- log2(IFI30_exp$IFI30+1)
IFI30_exp$TissueType <- ifelse(IFI30_exp$TissueType == "Paratumor", "P","T")
IFI30_exp <- filter(IFI30_exp, TissueType %in% c("P","T"))

p4 <- ggplot(IFI30_exp, aes_string(x = "TissueType",y = "IFI30" ))+
  # geom_boxplot(outliers = FALSE)+
  geom_boxplot(aes(color=TissueType),position = position_dodge(0.8),width=0.5,outlier.colour = "white")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15),shape=21, size = 0.3, stroke=0.8,aes(color=TissueType),alpha = 0.5)+
  # geom_violin(width = 0.8)+
  # geom_boxplot(width=0.1, fill="white",outlier.shape = NA)+
  # scale_fill_manual(values=c(mycolor_new))+ #颜色
  scale_color_manual(values=c(mycolor_new))+ #颜色
  guides(fill=FALSE) +
  stat_summary(fun.data = mean_fun, geom = "text", size=3 , color="black") +
  labs(title = "",
       x = '',
       y = "IFI30 expression") +
  theme_classic(base_line_size = 0.5)+# 主题设置，这个是无线条主题 # 坐标轴的粗细
  theme(plot.title = element_text(size = 10,colour = "black",hjust = 0.5),
        # axis.title.y = element_text(size = 15, color = "black",face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=10, face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, face = "bold"),
        axis.text.x = element_text(size = 10,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 0.95,
                                   angle = 45), #角度
        axis.text.y = element_text(size = 10,  color = "black",face = "bold", vjust = 0.5,
                                   hjust = 0.5, angle = 0),
        panel.border = element_rect(linetype = 'solid', size = 0.5,fill = NA)
  )+
  stat_compare_means(comparisons =  list(c("T","P")),
                     method = 'wilcox.test',hide.ns = TRUE,label = "p.format",bracket.size = 0.1,size = 2)
p4
save_plot(p4, file.path(out_dir, "boxplot_TCGA"),if_pdf = T,width = 4,height = 4)
