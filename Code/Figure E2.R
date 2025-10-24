#### Figure E2A ====
bcr_merge <- readRDS("~/merged_bcr.rds")
cur_sample <- "2907T"
tmp_bcr <- bcr_merge[bcr_merge$sample == "2907T",]
tmp_bcr$bin_id <- gsub(glue("{cur_sample}_"), "", tmp_bcr$bin_id)
tmp_heavy <- subset(tmp_bcr,tmp_bcr$locus %in% c("IGH"))
tmp_heavy <- tmp_heavy$bin_id %>% table() %>% as.data.frame()
colnames(tmp_heavy) <- c("bin_id","heavy")
tmp_heavy$heavy <- "IgH"
obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                           sep = "\t",
                           header = T,
                           row.names = 1) %>% as.data.frame()
obj@meta.data$bin_id <- rownames(obj@meta.data)
obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_heavy, by = "bin_id") %>%
  column_to_rownames(var="bin_id")
obj@meta.data$heavy <- ifelse(is.na(obj@meta.data$heavy), "other", "IgH")
obj@meta.data$heavy <- ifelse(obj@meta.data$heavy == "IgH", "IgH", obj@meta.data$Bin_Region)
obj@meta.data$heavy <- ifelse(obj@meta.data$heavy == "Tumor_capsule", "Invasive_zone", obj@meta.data$heavy)
obj@meta.data$heavy[obj@meta.data$heavy == "Tumor_tissue"] <- "Tumor"
SpatialDimPython(obj = obj,
                 prefix = "~",
                 plot_item = "heavy",
                 plot_order = c("IgH","Tumor"),
                 tmp_color = c("red","#E9ECEB"),
                 name = glue("{cur_sample}_BCR_heavy_dimplot"))


bcr_merge <- readRDS("~/merged_bcr.rds")
cur_sample <- "2907T"
tmp_bcr <- bcr_merge[bcr_merge$sample == "2907T",]
tmp_bcr$bin_id <- gsub(glue("{cur_sample}_"), "", tmp_bcr$bin_id)
tmp_heavy <- subset(tmp_bcr, !tmp_bcr$locus %in% c("IGH"))
tmp_heavy <- tmp_heavy$bin_id %>% table() %>% as.data.frame()
colnames(tmp_heavy) <- c("bin_id","heavy")
tmp_heavy$heavy <- "IgL"
obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                           sep = "\t",
                           header = T,
                           row.names = 1) %>% as.data.frame()
obj@meta.data$bin_id <- rownames(obj@meta.data)
obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_heavy, by = "bin_id") %>%
  column_to_rownames(var="bin_id")
obj@meta.data$heavy <- ifelse(is.na(obj@meta.data$heavy), "other", "IgL")
obj@meta.data$heavy <- ifelse(obj@meta.data$heavy == "IgL", "IgL", obj@meta.data$Bin_Region)
obj@meta.data$heavy <- ifelse(obj@meta.data$heavy == "Tumor_capsule", "Invasive_zone", obj@meta.data$heavy)
obj@meta.data$heavy[obj@meta.data$heavy == "Tumor_tissue"] <- "Tumor"
SpatialDimPython(obj = obj,
                 prefix = "~",
                 plot_item = "heavy",
                 plot_order = c("IgL","Tumor"),
                 tmp_color = c('#4cb1d2',"#E9ECEB"),
                 name = glue("{cur_sample}_BCR_light_dimplot"))


tcr_merge <- readRDS("~/merged_tcr.rds")
tmp_tcr <- tcr_merge[tcr_merge$sample == cur_sample, ]
tmp_tcr$BinID <- gsub(glue("{cur_sample}_"), "", tmp_tcr$BinID)
tmp_tcr <- tmp_tcr[tmp_tcr$topChains == "TRB",]
tmp_tcr <- tmp_tcr$BinID %>% table() %>% as.data.frame()
colnames(tmp_tcr) <- c("bin_id","tcr")
tmp_tcr$tcr <- "TCR"

obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                           sep = "\t",
                           header = T,
                           row.names = 1) %>% as.data.frame()
obj@meta.data$bin_id <- rownames(obj@meta.data)

obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_tcr, by = "bin_id") %>%
  column_to_rownames(var="bin_id")
obj@meta.data$tcr <- ifelse(is.na(obj@meta.data$tcr), "other", "TCR")
obj@meta.data$tcr <- ifelse(obj@meta.data$tcr == "TCR", "TCR", obj@meta.data$Bin_Region)
obj@meta.data$tcr[obj@meta.data$tcr == "Tumor_tissue"] <- "Tumor"
SpatialDimPython(obj = obj, 
                 prefix = "~",
                 plot_item = "tcr",
                 plot_order = c("TCR","Tumor"), 
                 tmp_color = c("#00BF05","#E9ECEB"),
                 name = glue("{cur_sample}_TCRB_dimplot"))



tcr_merge <- readRDS("~/merged_tcr.rds")
tmp_tcr <- tcr_merge[tcr_merge$sample == cur_sample, ]
tmp_tcr$BinID <- gsub(glue("{cur_sample}_"), "", tmp_tcr$BinID)
tmp_tcr <- tmp_tcr[tmp_tcr$topChains == "TRA",]
tmp_tcr <- tmp_tcr$BinID %>% table() %>% as.data.frame()
colnames(tmp_tcr) <- c("bin_id","tcr")
tmp_tcr$tcr <- "TCR"
obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                           sep = "\t",
                           header = T,
                           row.names = 1) %>% as.data.frame()
obj@meta.data$bin_id <- rownames(obj@meta.data)
obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_tcr, by = "bin_id") %>%
  column_to_rownames(var="bin_id")
obj@meta.data$tcr <- ifelse(is.na(obj@meta.data$tcr), "other", "TCR")
obj@meta.data$tcr <- ifelse(obj@meta.data$tcr == "TCR", "TCR", obj@meta.data$Bin_Region)
obj@meta.data$tcr[obj@meta.data$tcr == "Tumor_tissue"] <- "Tumor"


SpatialDimPython(obj = obj, 
                 prefix = "~",
                 plot_item = "tcr",
                 plot_order = c("TCR","Tumor"), 
                 tmp_color = c("#C198C4","#E9ECEB"),
                 name = glue("{cur_sample}_TRA_dimplot"))




bcr_merge <- readRDS("~/merged_bcr.rds")
tmp_bcr <- bcr_merge[bcr_merge$sample == "2907B",]
tmp_heavy <- subset(tmp_bcr, tmp_bcr$locus %in% c("IGH"))
tmp_heavy$bin_id <- str_remove_all(tmp_heavy$bin_id, "2907B_")
tmp_heavy <- tmp_heavy$bin_id %>% table() %>% as.data.frame()
colnames(tmp_heavy) <- c("bin_id","heavy")
tmp_heavy$heavy <- "IgH"
cur_sample <- "2907B"
obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                           sep = "\t",
                           header = T,
                           row.names = 1) %>% as.data.frame()
obj@meta.data$bin_id <- rownames(obj@meta.data)
obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_heavy, by = "bin_id") %>%
  column_to_rownames(var="bin_id")
obj@meta.data$heavy <- ifelse(is.na(obj@meta.data$heavy), "other", "IgH")
obj@meta.data$heavy <- ifelse(obj@meta.data$heavy == "IgH", "IgH", obj@meta.data$Bin_Region)
obj@meta.data$heavy <- ifelse(obj@meta.data$heavy == "Tumor_capsule", "Invasive_zone", obj@meta.data$heavy)
obj@meta.data$heavy[obj@meta.data$heavy == "Tumor_side_of_Margin_area"] <- "Tumor"
obj@meta.data$heavy[obj@meta.data$heavy == "Paratumor_side_of_Margin_area"] <- "Paratumor"

SpatialDimPython(obj = obj,
                 prefix = "~",
                 plot_item = "heavy",
                 plot_order = c("IgH","Tumor","Invasive_zone","Paratumor"),
                 tmp_color = c("red","#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2907B_BCR_heavy_dimplot"))
obj$TLS <- obj$Bin_Region
obj$TLS[obj$TLS != "Tumor_capsule"] <- "_NA"
SpatialDimPython(obj = obj,
                 prefix = "~",
                 plot_item = "heavy",
                 if_tls_contour = "T", 
                 plot_order = c("IgH","Tumor","Invasive_zone","Paratumor"),
                 tmp_color = c("red","#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2907B_IgH_v2"))


bcr_merge <- readRDS("~/merged_bcr.rds")
tmp_bcr <- bcr_merge[bcr_merge$sample == "2907B",]
tmp_heavy <- subset(tmp_bcr, !tmp_bcr$locus %in% c("IGH"))
tmp_heavy$bin_id <- str_remove_all(tmp_heavy$bin_id, "2907B_")
tmp_heavy <- tmp_heavy$bin_id %>% table() %>% as.data.frame()
colnames(tmp_heavy) <- c("bin_id","heavy")
tmp_heavy$heavy <- "IgH"
cur_sample <- "2907B"
obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                           sep = "\t",
                           header = T,
                           row.names = 1) %>% as.data.frame()
obj@meta.data$bin_id <- rownames(obj@meta.data)
obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_heavy, by = "bin_id") %>%
  column_to_rownames(var="bin_id")
obj@meta.data$heavy <- ifelse(is.na(obj@meta.data$heavy), "other", "IgH")
obj@meta.data$heavy <- ifelse(obj@meta.data$heavy == "IgH", "IgH", obj@meta.data$Bin_Region)
obj@meta.data$heavy <- ifelse(obj@meta.data$heavy == "Tumor_capsule", "Invasive_zone", obj@meta.data$heavy)
obj@meta.data$heavy[obj@meta.data$heavy == "Tumor_side_of_Margin_area"] <- "Tumor"
obj@meta.data$heavy[obj@meta.data$heavy == "Paratumor_side_of_Margin_area"] <- "Paratumor"

SpatialDimPython(obj = obj,
                 prefix = "~",
                 plot_item = "heavy",
                 plot_order = c("IgH","Tumor","Invasive_zone","Paratumor"),
                 tmp_color = c('#4cb1d2',"#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2907B_BCR_heavy_dimplot"))
obj$TLS <- obj$Bin_Region
obj$TLS[obj$TLS != "Tumor_capsule"] <- "_NA"
SpatialDimPython(obj = obj,
                 prefix = "~",
                 plot_item = "heavy",
                 if_tls_contour = "T",
                 plot_order = c("IgH","Tumor","Invasive_zone","Paratumor"),
                 tmp_color = c('#4cb1d2',"#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2907B_IgL_v2"))

tcr_merge <- readRDS("~/merged_tcr.rds")
tmp_tcr <-  tcr_merge[tcr_merge$sample == "2907B",]
tmp_tcr <- tmp_tcr[tmp_tcr$topChains == "TRB",]
tmp_tcr$Bin50 <- str_remove_all(tmp_tcr$Bin50, "2907B_")
tmp_tcr <- tmp_tcr$Bin50 %>% table() %>% as.data.frame()
cur_sample <- "2907B"
obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                           sep = "\t",
                           header = T,
                           row.names = 1) %>% as.data.frame()
obj@meta.data$bin_id <- rownames(obj@meta.data)
colnames(tmp_tcr) <- c("bin_id","tcr")
tmp_tcr$tcr <- "TCR"
obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_tcr, by = "bin_id") %>%
  column_to_rownames(var="bin_id")
obj@meta.data$tcr <- ifelse(is.na(obj@meta.data$tcr), "other", "TCR")
obj@meta.data$tcr <- ifelse(obj@meta.data$tcr == "TCR", "TCR", obj@meta.data$Bin_Region)
obj@meta.data$tcr <- ifelse(obj@meta.data$tcr == "Tumor_capsule", "Invasive_zone", obj@meta.data$tcr)
obj@meta.data$tcr[obj@meta.data$tcr == "Tumor_side_of_Margin_area"] <- "Tumor"
obj@meta.data$tcr[obj@meta.data$tcr == "Paratumor_side_of_Margin_area"] <- "Paratumor"


SpatialDimPython(obj = obj, 
                 prefix = "~",
                 plot_item = "tcr",
                 plot_order = c("TCR","Tumor","Invasive_zone","Paratumor"), 
                 tmp_color = c("#00BF05","#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2907B_TCRB_dimplot"))
obj$TLS <- obj$Bin_Region
obj$TLS[obj$TLS != "Tumor_capsule"] <- "_NA"
SpatialDimPython(obj = obj,
                 prefix = "~",
                 plot_item = "tcr",
                 if_tls_contour = "T",
                 plot_order = c("TCR","Tumor","Invasive_zone","Paratumor"),
                 tmp_color = c("#00BF05","#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2907B_TRB_v2"))




tcr_merge <- readRDS("~/merged_tcr.rds")
tmp_tcr <-  tcr_merge[tcr_merge$sample == "2907B",]
tmp_tcr <- tmp_tcr[tmp_tcr$topChains == "TRA",]
tmp_tcr$Bin50 <- str_remove_all(tmp_tcr$Bin50, "2907B_")
tmp_tcr <- tmp_tcr$Bin50 %>% table() %>% as.data.frame()
cur_sample <- "2907B"
obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                           sep = "\t",
                           header = T,
                           row.names = 1) %>% as.data.frame()
obj@meta.data$bin_id <- rownames(obj@meta.data)
colnames(tmp_tcr) <- c("bin_id","tcr")
tmp_tcr$tcr <- "TCR"
obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_tcr, by = "bin_id") %>%
  column_to_rownames(var="bin_id")
obj@meta.data$tcr <- ifelse(is.na(obj@meta.data$tcr), "other", "TCR")
obj@meta.data$tcr <- ifelse(obj@meta.data$tcr == "TCR", "TCR", obj@meta.data$Bin_Region)
obj@meta.data$tcr <- ifelse(obj@meta.data$tcr == "Tumor_capsule", "Invasive_zone", obj@meta.data$tcr)
obj@meta.data$tcr[obj@meta.data$tcr == "Tumor_side_of_Margin_area"] <- "Tumor"
obj@meta.data$tcr[obj@meta.data$tcr == "Paratumor_side_of_Margin_area"] <- "Paratumor"


SpatialDimPython(obj = obj, 
                 prefix = "~",
                 plot_item = "tcr",
                 plot_order = c("TCR","Tumor","Invasive_zone","Paratumor"), 
                 tmp_color = c("#C198C4","#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2907B_TRA_dimplot"))
obj$TLS <- obj$Bin_Region
obj$TLS[obj$TLS != "Tumor_capsule"] <- "_NA"
SpatialDimPython(obj = obj,
                 prefix = "~",
                 plot_item = "tcr",
                 if_tls_contour = "T",
                 plot_order = c("TCR","Tumor","Invasive_zone","Paratumor"),
                 tmp_color = c("#C198C4","#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2907B_TRA_v2"))






#### Figure E2B ====
pacbio_merge <- readRDS("~/merged_bcr.rds")
pacbio_merge$Bin_Region[pacbio_merge$Bin_Region == "Tumor_side_of_Margin_area"] <- "Tumor"
pacbio_merge$Bin_Region[pacbio_merge$Bin_Region == "Tumor_tissue"] <- "Tumor"
pacbio_merge$Bin_Region[pacbio_merge$Bin_Region == "Paratumor_side_of_Margin_area"] <- "Paratumor"
pacbio_merge$Bin_Region[pacbio_merge$Bin_Region == "Paratumor_tissue"] <- "Paratumor"
pacbio_merge$Bin_Region[pacbio_merge$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"

bcr_region_density <- pacbio_merge %>% 
  
  dplyr::group_by(sample, Bin_Region) %>% 
  dplyr::summarise(bcr_count = n()) %>% 
  dplyr::ungroup()
meta_anno <- fread("~/Bin50_classic_celltype_metadata_merge.txt")
meta_anno$Bin_Region[meta_anno$Bin_Region == "Tumor_side_of_Margin_area"] <- "Tumor"
meta_anno$Bin_Region[meta_anno$Bin_Region == "Tumor_tissue"] <- "Tumor"
meta_anno$Bin_Region[meta_anno$Bin_Region == "Paratumor_side_of_Margin_area"] <- "Paratumor"
meta_anno$Bin_Region[meta_anno$Bin_Region == "Paratumor_tissue"] <- "Paratumor"
meta_anno$Bin_Region[meta_anno$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"

region_count <- meta_anno %>% 
  dplyr::group_by(SampleID, Bin_Region) %>% 
  dplyr::summarise(region_count = n()) %>% 
  dplyr::rename(sample = SampleID) %>% 
  dplyr::ungroup()
bcr_region_density <- bcr_region_density %>%
  dplyr::left_join(region_count, by = c("sample", "Bin_Region"))
bcr_region_density$bcr_density <- bcr_region_density$bcr_count/bcr_region_density$region_count
bcr_region_density$Bin_Region <- factor(bcr_region_density$Bin_Region,
                                        levels = c("Tumor" , 
                                                   "Invasive_zone", 
                                                   "Paratumor"))
wilcox.test(bcr_region_density$bcr_density[bcr_region_density$Bin_Region == "Tumor"],
            bcr_region_density$bcr_density[bcr_region_density$Bin_Region == "Invasive_zone"])
wilcox.test(bcr_region_density$bcr_density[bcr_region_density$Bin_Region == "Tumor"],
            bcr_region_density$bcr_density[bcr_region_density$Bin_Region == "Paratumor"])
wilcox.test(bcr_region_density$bcr_density[bcr_region_density$Bin_Region == "Invasive_zone"],
            bcr_region_density$bcr_density[bcr_region_density$Bin_Region == "Paratumor"])
bcr_region_density %>% dplyr::group_by(Bin_Region) %>% 
  dplyr::summarise(mean = mean(bcr_density))




p <- ggbarplot(bcr_region_density,
               x = "Bin_Region", y = "bcr_density", color = "Bin_Region",
               palette = c("Tumor" = "#6F5D96", "Invasive_zone" = "#537AAB", "Paratumor" = "#9FB9D0"),
               add = c("mean_se"),
               size = 0.5)+ 
  # stat_compare_means(comparisons = combn(c("Tumor","Invasive_zone","Paratumor"),
  #                                        2, simplify = FALSE),
  #                    method = "wilcox.test", label = "p.format",
  #                    label.y = c(0.7, 0.8, 0.9)) +
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
            legend.position = "none") +
  labs(x = "", y = "bcr density")


save_plot(p, glue("~/bcr_density_bar"),
          width = 3, height = 4, if_png = T,if_pdf = T)

#### Figure E2C ====

shared_sample_list <- c("2907", "2931", "2903", "RH17", "RH16", "RH14", "RH10", "RH08", "RH05", "2976",
                        "2893", "2846")
share_merge <- list()
share_count <- list()
for(s in shared_sample_list){
  cur_rds <- readRDS(glue('~/共享RDS/{s}_h_anno_v2.rds'))
  cur_rds$sampleid <- s
  shared_summarise <- cur_rds
  shared_summarise$shared_info <- ""
  shared_summarise$shared_info[grepl("Shared", shared_summarise$region_shared_sub)] <- "Shared" 
  shared_summarise$shared_info[grepl("Private", shared_summarise$region_shared_sub)] <- "Private"
  nrow(shared_summarise)
  # shared_summarise <- shared_summarise[shared_summarise$shared_info != shared_summarise$region,]
  nrow(shared_summarise)
  unique(shared_summarise$shared_info)
  
  shared_summarise$shared_info <- factor(shared_summarise$shared_info, levels = c(
    "Private",
    "Shared"
  ))
  shared_summarise <- shared_summarise[shared_summarise$region != "PBMC",]
  
  shared_summarise$region <- factor(shared_summarise$region, levels = c(
    "Tumor",
    "Invasive_zone",
    "Paratumor"
  ))
  input <- prop.table(table(shared_summarise$region , 
                            shared_summarise$shared_info), margin = 1)
  input <- as.data.frame(input)
  input$Var1 <- factor(input$Var1, levels = c(
    "Tumor",
    "Invasive_zone",
    "Paratumor"
  ) %>% rev())
  
  input$sampleid <- s
  # input <- input[input$Var2 == "Shared", ]
  
  spatial_count <- shared_summarise[shared_summarise$type != "SingleCell",]
  spatial_count <- table(spatial_count$sampleid, spatial_count$region) %>% as.data.frame()
  colnames(spatial_count) <- c("sampleid", "region", "spatial_bcr_count")
  input <- dplyr::left_join(input, spatial_count, by = c("Var1"= "region", "sampleid" = "sampleid"))
  
  share_merge[[s]] <- input
  
}
share_merge <- dplyr::bind_rows(share_merge)
share_merge <- share_merge[!is.nan(share_merge$Freq),]

bcr_count <- share_merge %>% 
  dplyr::distinct(sampleid, Var1, .keep_all = T) %>% 
  dplyr::group_by(sampleid) %>% 
  dplyr::summarise(bcr_count = sum(spatial_bcr_count))

high_share <- c("2907",
                "2931",
                "2903",
                "RH16",
                "RH05")


low_share <- c("RH17",
               "RH14",
               "RH10",
               "RH08",
               "2976",
               "2893",
               "2846")


bcr_count$group[bcr_count$sampleid %in% high_share] <- "high_share"
bcr_count$group[bcr_count$sampleid %in% low_share] <- "low_share"
bcr_count$group <- factor(bcr_count$group, levels = c(
  "high_share",
  "low_share"
))


share_merge <- share_merge[share_merge$Var2 == "Shared", ]

share_merge$Var1 <- factor(share_merge$Var1, levels = c(
  "Tumor",
  "Invasive_zone",
  "Paratumor"
))
share_merge <- dplyr::arrange(share_merge, Var1)
p <- ggbarplot(share_merge,
               x = "Var1", y = "Freq", color = "Var1",
               palette = c("Tumor" = "#6F5D96", "Invasive_zone" = "#537AAB", "Paratumor" = "#9FB9D0"),
               add = c("mean_se", "jitter"),
               size = 0.5,
               legend = "none" )+ 
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  labs(x = "", y = "share ratio")


save_plot(p, glue("~/share_ratio"),
          width = 4, height = 4, if_png = T)
save_plot(p, glue("~/share_ratio"),
          width = 4, height = 4, if_pdf = T)


#### Figure E2D ====
pbmc_sample <- c("RH16_PBMC",  "LC_2931", "LC_2903", "RH06_PBMC", "LC_2931after", "LC_2846" ,
                 "LC_2907_3", "LC_2976after", "LC_2907", "RH06APBMC", "LC_2893after", "LC_2893",
                 "LC_2846after", "LC_2976" , "RH08APBMC", "RH05_PBMC", "RH05APBMC",
                 "LC_2740","RH03_PBMC",    "LC_2903after", "RH17APBMC", "RH14_PBMC", "LC_2740after",
                 "RH15APBMC" , "LC_2896after", "LC_2896", "RH14APBMC",  "RH15PBMC","RH17_1_PBMC",
                 "RH16APBMC", "RH10APBMC" , "RH03APBMC",  "RH08_PBMC", "RH10_RPBMC" )
p_sample <- c("RH16_A_P", "2837P" , "2903P" , "RH05P",  "2931P", "2823P",
              "RH17_P", "2934P", "2846P", "2893P", "RH03_RP", "3099P", "RH08P",
              "2896P","RH14_A_P", "2907P", "3154_P", "3001P", "2865P", "2976P", "RH10_RP", "RH10_A_P", "2740P" ,
              "3163_P", "RH15_A_P", "3125_P")
tumor_sample <- c("RH05", "2931T_B_FS", "2931T", "RH16", "RH16_A_T", "RH17_T" ,
                  "RH05_A", "2837T",
                  "2976T_B_FS_3", "2976T_B_FS_2", "2976T_B_FS_1", 
                  "2976T",  "RH15", "RH08_A",  "2893T","RH06", "RH03","2879T",
                  "RH03_A","3001T", "3154_T", "RH10_A_T",
                  "RH17_AT1", "RH14_A_T", "RH10","RH15_A_T", "2846T", "RH06_A", "RH08", "RH03_RT",
                  "3125_T", "RH10_RT", "3099T" , "2934T" , "2865T" , "RH17_AT2", 
                  "2907T", "3163_T", "2823T")
border_sample <- c( "2931B")



pacbio_sc_info <- fread("/Users/mac/Desktop/Rproject/Pacbio/Pacbio_SingleCell_43.csv")
changeo_out <- readRDS("~/merge4.rds")

#1
sample_sc <- c("2931T", "2931B","2931P", "2931T_B_FS","LC_2931","LC_2931after")
sample_sc <- paste(sample_sc, "SingleCell", sep = "_")
#2
sample_st <- c("2931T", "2931B","2931P")
sample_st <- paste(sample_st, "Spatial", sep = "_")
sub_seq <- changeo_out[changeo_out$sample_new %in% c(sample_sc, sample_st),]
nrow(sub_seq)
sub_seq <- sub_seq %>% 
  dplyr::select(
    -"sample_shared_sub", -"sample_shared_main", -"type_shared_sub", -"type_shared_main",       
    -"sample_new_shared_sub", -"sample_new_shared_main" , - "seq_count", -"clone_num"
  )



bcr_merge <- readRDS("~/merged_bcr.rds")
bcr_merge$BinID <- bcr_merge$bin_id
meta_50 <- fread("~/Bin50_classic_celltype_metadata_merge.txt")
meta_50$Bin_Region[meta_50$Bin_Region == "Tumor_side_of_Margin_area"] <- "Tumor"
meta_50$Bin_Region[meta_50$Bin_Region == "Tumor_tissue"] <- "Tumor"
meta_50$Bin_Region[meta_50$Bin_Region == "Paratumor_side_of_Margin_area"] <- "Paratumor"
meta_50$Bin_Region[meta_50$Bin_Region == "Paratumor_tissue"] <- "Paratumor"
meta_50$Bin_Region[meta_50$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"


sub_seq <- dplyr::left_join(sub_seq, bcr_merge[, c("sequence_id", "BinID")])
sub_seq <- dplyr::left_join(sub_seq, meta_50)

sub_seq <- sub_seq %>% dplyr::mutate(region = dplyr::case_when(
  type == "Spatial" ~ Bin_Region,
  sample %in% pbmc_sample ~ "PBMC", 
  sample %in% p_sample ~ "Paratumor",
  sample %in% tumor_sample ~"Tumor",
  sample %in% border_sample ~ "Boder_Singlecell"
))
shared_info <- get_shared_info(sub_seq, source  = "clone_id",target_list =  c("region", "sample_new",  "type")
)



shared_info_h <- shared_info[shared_info$locus == "IGH",]
saveRDS(shared_info_h, "~/2931_h_anno_v2.rds")

shared_info_h <- readRDS("~/2931_h_anno_v2.rds")
shared_summarise <- shared_info_h
shared_summarise$shared_info <- ""
shared_summarise$shared_info[grepl("Shared", shared_summarise$region_shared_sub)] <- "Shared" 
shared_summarise$shared_info[grepl("Private", shared_summarise$region_shared_sub)] <- "Private"
nrow(shared_summarise)
# shared_summarise <- shared_summarise[shared_summarise$shared_info != shared_summarise$region,]
nrow(shared_summarise)
unique(shared_summarise$shared_info)

shared_summarise$shared_info <- factor(shared_summarise$shared_info, levels = c(
  "Private",
  "Shared"
))
shared_summarise$region <- factor(shared_summarise$region, levels = c(
  "Tumor",
  "Invasive_zone",
  "Paratumor",
  "PBMC"
))
input <- prop.table(table(shared_summarise$region , 
                          shared_summarise$shared_info), margin = 1)
input
input <- as.data.frame(input)
input$Var1 <- factor(input$Var1, levels = c(
  "Tumor",
  "Invasive_zone",
  "Paratumor",
  "PBMC"
) %>% rev())

names(col_clone_expanded) <- c("Shared", "Private")
p <- ggplot(input, aes_string("Var1", "Freq", fill = "Var2"))+
  geom_bar(stat = "identity", position="stack")+
  labs(x = "", y = "", title = "")+
  scale_fill_manual(values = col_clone_expanded)+
  coord_flip() +
  # scale_y_discrete(position = "right", sec.axis = sec_axis(~.)) + 
  # scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "Y-axis Label on right side")) +
  scale_y_reverse()  +
  guides(
    x = "none",
    x.sec = "axis",
    y = "none",
    y.sec = "axis"
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x =  element_text(size = 10, color = "black", face = "bold",angle = 90,vjust = 1,hjust = 0.95 ),
    axis.line = element_line(size = 0.6),
    axis.title = element_text(size = 15, color = "black", face = "bold"),
    legend.text = element_text(size = 10, color = "black", face = "bold"),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"),
    plot.title = element_text(size = 20, color = "black", face = "bold",hjust = 0.5),
    legend.title = element_blank()
  )


save_plot(p, glue("~/shared_info_bar_number_2931"),
          width = 5, height = 3, if_pdf = T)
save_table(input %>% dplyr::arrange(desc(Var1)),
           glue("~/input"))


shared_info_h <- readRDS("~/2931_h_anno_v2.rds")

c_share_info <- shared_info_h
c_share_info$region <- factor(c_share_info$region,
                              levels = c("Tumor", "Invasive_zone" , 
                                         "Paratumor"  , "PBMC" ))
lt <- split(c_share_info$clone_name, c_share_info$region)
m = make_comb_mat(lt)
plot.new()
pdf(glue("~/upset_v3_number_2931.pdf"),
    width = 6,
    height = 3)
ht_raw  = UpSet(m, set_order = c("Tumor", "Invasive_zone" , 
                                 "Paratumor", "PBMC" ))
draw(ht_raw)
dev.off()

lt_2 <- split(c_share_info$clone_name, c_share_info$region)

m = make_comb_mat(lt_2)
cs = comb_size(m)
cs <- as.data.frame(cs)
cs$name <- rownames(cs)
cs$type <- "1"
cs$cs[cs$name == "1110"] <- nrow(shared_info_h[shared_info_h$region_shared_sub == "Shared_Invasive_zone, Paratumor, Tumor",]) 
cs$cs[cs$name == "1100"] <- nrow(shared_info_h[shared_info_h$region_shared_sub == "Shared_Invasive_zone, Tumor" ,])
cs$cs[cs$name == "1010"] <- nrow(shared_info_h[shared_info_h$region_shared_sub == "Shared_Paratumor, Tumor"  ,])
cs$cs[cs$name == "1001"] <- nrow(shared_info_h[shared_info_h$region_shared_sub == "Shared_PBMC, Tumor",])
cs$cs[cs$name == "0110"] <- nrow(shared_info_h[shared_info_h$region_shared_sub == "Shared_Invasive_zone, Paratumor",])
cs$cs[cs$name == "0101"] <- nrow(shared_info_h[shared_info_h$region_shared_sub == "Shared_Invasive_zone, PBMC" ,])
cs$cs[cs$name == "1000"] <- nrow(shared_info_h[shared_info_h$region_shared_sub == "Private_Tumor" ,])
cs$cs[cs$name == "0100"] <- nrow(shared_info_h[shared_info_h$region_shared_sub ==  "Private_Invasive_zone"  ,])
cs$cs[cs$name == "0010"] <- nrow(shared_info_h[shared_info_h$region_shared_sub ==  "Private_Paratumor" ,])
cs$cs[cs$name == "0001"] <- nrow(shared_info_h[shared_info_h$region_shared_sub ==  "Private_PBMC"  ,])



p <- ggbarplot(cs, "name", "cs",
               fill = "grey", color = "grey", 
               label = TRUE, lab.col = "black", lab.pos = "out") + 
  labs(title="",x ="", y = "number")+
  scale_y_continuous(expand = c(0,0),limits= c(0, 5200))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    axis.line = element_line(size = 0.6),
    axis.title = element_text(size = 15, color = "black", face = "bold"),
    legend.text = element_text(size = 10, color = "black", face = "bold"),
    legend.position="none",
  )
p
save_plot(p, glue("~/com_size_v3_number_2931"),
          width = 6, height = 3, if_png = T)
save_table(cs, glue("~/cs"))

cur_sample <- "2931T"
count_mat_13
cur_clone <- "112117"
hcc_st <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
hcc_st@meta.data$BinID <-  rownames(hcc_st@meta.data)
hcc_st@meta.data  <-  dplyr::select(hcc_st@meta.data,"BinID")

metadata <-meta_50 %>% filter(SampleID == cur_sample) %>% 
  mutate(BinID = str_replace(BinID, "\\w{5}_", ""))

hcc_st@meta.data <- hcc_st@meta.data %>% left_join(metadata, by = "BinID") %>% 
  column_to_rownames(var = "BinID")

hcc_st@meta.data$SpotLight_Anno[is.na(hcc_st@meta.data$SpotLight_Anno)] <- "None"
hcc_st@meta.data$CellSubType[is.na(hcc_st@meta.data$CellSubType)] <- "None"
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$SpotLight_Anno == "pDC"] <- "Myeloid"
hcc_st@meta.data$BinID <- rownames(hcc_st@meta.data)


tmp_bcr <- shared_info_h[shared_info_h$sample == cur_sample & shared_info_h$clone_id == cur_clone, ]
tmp_bcr$BinID <- str_replace(tmp_bcr$BinID, "\\w{5}_", "") %>% as.character()
tmp_bcr <- tmp_bcr$BinID %>% table() %>% as.data.frame()

colnames(tmp_bcr) <- c("BinID", "bcr")
tmp_bcr$bcr <- "bcr"
hcc_st@meta.data <- dplyr::left_join(hcc_st@meta.data, tmp_bcr) %>% column_to_rownames(var="BinID")
hcc_st@meta.data$bcr <- ifelse(is.na(hcc_st@meta.data$bcr), "other", "bcr")
hcc_st@meta.data$bcr <- ifelse(hcc_st@meta.data$bcr == "bcr", "bcr", hcc_st@meta.data$Bin_Region)
hcc_st@meta.data$bcr <- ifelse(hcc_st@meta.data$bcr == "Tumor_capsule", "Invasive_zone", hcc_st@meta.data$bcr)
hcc_st_e10 <- subset(hcc_st, bcr == "bcr")
# hcc_st@meta.data$bcr_1 <- hcc_st@meta.data$bcr
# hcc_st@meta.data$bcr_1[hcc_st@meta.data$bcr_1 == "Tumor_tissue"] <- 
hcc_st_e10 <-   hcc_st_e10@meta.data %>%
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "bcr", 
    mode = 2
  )
hcc_st@meta.data$BinID <- rownames(hcc_st@meta.data)
hcc_st_1 <- hcc_st
hcc_st_1@meta.data <- hcc_st_1@meta.data %>% 
  dplyr::left_join(hcc_st_e10, by = c("row", "col"))
hcc_st_1$bcr_1 <- hcc_st_1$expand_n
hcc_st_1$bcr_1[is.na(hcc_st_1$bcr_1)] <- "Tumor_tissue"
hcc_st_1$bcr_1[hcc_st_1$bcr_1!="Tumor_tissue"] <- "bcr"

rownames(hcc_st_1@meta.data) <- hcc_st_1@meta.data$BinID

SpatialDimPython(obj = hcc_st_1,
                 prefix = "~",
                 plot_item = "bcr_1",
                 plot_order = c("bcr","Tumor_tissue"),
                 tmp_color = c("red","#E9ECEB"),
                 name = glue("2931T_clone_112117_2"))

cur_sample <- "2931B"
count_mat_13
cur_clone <- "112117"
hcc_st <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
hcc_st@meta.data$BinID <-  rownames(hcc_st@meta.data)
hcc_st@meta.data  <-  dplyr::select(hcc_st@meta.data,"BinID")

metadata <-meta_50 %>% filter(SampleID == cur_sample) %>% 
  mutate(BinID = str_replace(BinID, "\\w{5}_", ""))

hcc_st@meta.data <- hcc_st@meta.data %>% left_join(metadata, by = "BinID") %>% 
  column_to_rownames(var = "BinID")

hcc_st@meta.data$SpotLight_Anno[is.na(hcc_st@meta.data$SpotLight_Anno)] <- "None"
hcc_st@meta.data$CellSubType[is.na(hcc_st@meta.data$CellSubType)] <- "None"
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$SpotLight_Anno == "pDC"] <- "Myeloid"
hcc_st@meta.data$BinID <- rownames(hcc_st@meta.data)


tmp_bcr <- shared_info_h[shared_info_h$sample == cur_sample & shared_info_h$clone_id == cur_clone, ]
tmp_bcr$BinID <- str_replace(tmp_bcr$BinID, "\\w{5}_", "") %>% as.character()
tmp_bcr <- tmp_bcr$BinID %>% table() %>% as.data.frame()
hcc_st@meta.data$BinID <- rownames(hcc_st@meta.data)
colnames(tmp_bcr) <- c("BinID", "IgH")
tmp_bcr$bcr <- "IgH"
hcc_st@meta.data <- dplyr::left_join(hcc_st@meta.data, tmp_bcr) %>% column_to_rownames(var="BinID")
hcc_st@meta.data$bcr <- ifelse(is.na(hcc_st@meta.data$bcr), "other", "IgH")
hcc_st@meta.data$bcr <- ifelse(hcc_st@meta.data$bcr == "IgH", "IgH", hcc_st@meta.data$Bin_Region)
hcc_st@meta.data$bcr <- ifelse(hcc_st@meta.data$bcr == "Tumor_capsule", "Invasive_zone", 
                               hcc_st@meta.data$bcr)
hcc_st@meta.data$bcr[hcc_st@meta.data$bcr == "Tumor_side_of_Margin_area"] <- "Tumor"
hcc_st@meta.data$bcr[hcc_st@meta.data$bcr == "Paratumor_side_of_Margin_area"] <- "Paratumor"


hcc_st_e10 <- subset(hcc_st, bcr == "IgH")
hcc_st_e10 <-   hcc_st_e10@meta.data %>%
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "bcr", 
    mode = 2
  )
hcc_st@meta.data$BinID <- rownames(hcc_st@meta.data)
hcc_st_1 <- hcc_st
hcc_st_1@meta.data <- hcc_st_1@meta.data %>% 
  dplyr::left_join(hcc_st_e10, by = c("row", "col"))
hcc_st_1$bcr_1 <- hcc_st_1$expand_n
hcc_st_1$bcr_1[!is.na(hcc_st_1$bcr_1)] <- "bcr"
hcc_st_1$bcr_1[is.na(hcc_st_1$bcr_1)] <- hcc_st_1$bcr[is.na(hcc_st_1$bcr_1)]
rownames(hcc_st_1@meta.data) <- hcc_st_1@meta.data$BinID


hcc_st_1$TLS <- hcc_st_1$Bin_Region
hcc_st_1$TLS[hcc_st_1$TLS != "Tumor_capsule"] <- "_NA"
SpatialDimPython(obj = hcc_st_1,
                 prefix = "~",
                 plot_item = "bcr_1",
                 if_tls_contour = "T",
                 plot_order = c("bcr","Tumor","Invasive_zone","Paratumor"),
                 tmp_color = c("red","#E9ECEB","#FFE8BF","#F9F7EE"),
                 name = glue("2931B_v2"))

#### Figure E2E ====

meta_anno <- fread("~/Bin50_classic_celltype_metadata_merge.txt")
meta_anno$Bin_Region[meta_anno$Bin_Region == "Tumor_side_of_Margin_area"] <- "Tumor"
meta_anno$Bin_Region[meta_anno$Bin_Region == "Tumor_tissue"] <- "Tumor"
meta_anno$Bin_Region[meta_anno$Bin_Region == "Paratumor_side_of_Margin_area"] <- "Paratumor"
meta_anno$Bin_Region[meta_anno$Bin_Region == "Paratumor_tissue"] <- "Paratumor"
meta_anno$Bin_Region[meta_anno$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"

meta_anno_sub <- meta_anno
TLS_data <- meta_anno_sub %>%
  mutate(
    TLS_group = if_else(grepl("_NA", TLS_raw), "NoTLS", "TLS")
  )
TLS_summary <- TLS_data %>%
  group_by(SampleID, Bin_Region) %>%
  summarise(
    TLS_ratio = mean(TLS_group == "TLS")   # “TLS” 占比
  )
TLS_summary$group <- TLS_summary$Bin_Region
TLS_summary$group <- factor(TLS_summary$group, levels = c(
  "Tumor",
  "Invasive_zone",
  "Paratumor"
))
p <- ggboxplot(TLS_summary,
               x = "group", y = "TLS_ratio", color = "group",
               palette = c("Tumor" = "#6F5D96", "Invasive_zone" = "#537AAB", "Paratumor" = "#9FB9D0"),
               add = c("jitter"),
               size = 0.5,
               legend = "none" )+ 
  stat_compare_means(comparisons = combn(c("Tumor","Invasive_zone","Paratumor"),
                                         2, simplify = FALSE),
                     method = "wilcox.test", label = "p.format",
                     label.y = c(0.7, 0.8, 0.9)) +
  theme(    plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  labs(x = "", y = "TLS density")

save_plot(p, glue("~/region_TLS_ratio"),
          width = 3, height = 4, if_png = T)
save_plot(p, glue("~/region_TLS_ratio"),
          width = 3, height = 4, if_pdf = T)


#### Figure E2 F-G ====
chk_12 <- read.csv("~/12趋化因子.csv", header = F)
chk_12 <- chk_12[,1]
cur_sample <- "2907T"
out_dir <- glue("~/B")
hcc_st <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
hcc_st@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                              sep = "\t",
                              header = T,
                              row.names = 1) %>% as.data.frame()
# hcc_st <- ST_nor_sca(hcc_st)
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$CellSubType == "Tumor"] <- "Tumor"
hcc_st@meta.data <- as.data.frame(hcc_st@meta.data)
hcc_st@meta.data$BinID <- rownames(hcc_st@meta.data)
rownames(hcc_st@meta.data) <- hcc_st@meta.data$BinID
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$SpotLight_Anno == "pDC"] <- "Myeloid"
hcc_st <- ST_nor_sca(hcc_st)

genelist <- fread("~/HCC_TLS_score.csv") %>% as.data.frame()
for (i in 1:ncol(genelist)) {
  gene<-genelist[,i]
  gene<-gene[gene!='']
  gene<-list(gene)
  name<-colnames(genelist)[i]
  hcc_st <- AddModuleScore(
    object = hcc_st,
    features = gene,
    ctrl = 100,
    name = name
  )
}
colnames(hcc_st@meta.data)[(ncol(hcc_st@meta.data)-1):ncol(hcc_st@meta.data)]<-colnames(genelist)
SpatialFeaturePython(obj = hcc_st, if_gene = F,
                     if_count = T,
                     plot_item = "chemokines12", # plot的基因
                     prefix = glue("{out_dir}"), # 当前路径
                     if_tls_contour = "F",
                     name = glue("{cur_sample}_chemokines12"))
SpatialFeaturePython(obj = hcc_st, if_gene = F,
                     if_count = T,
                     plot_item = "TLSsignatures", # plot的基因
                     prefix = glue("{out_dir}"), # 当前路径
                     if_tls_contour = "F",
                     name = glue("{cur_sample}_TLSsignatures"))


hcc_st$SpotLight_Anno <- factor(
  hcc_st$SpotLight_Anno,
  levels = c("B", "Plasma", "T", "Myeloid", "Endothelial", "Fibroblast", "Hep_tumor")
)
fill_colors <- c( "#238b45", "#fdbf6f", "#ce1256", "#6a3d9a", "#1f78b4", "#df65b0", "#cfc5ea")
SpatialDimPython(obj = hcc_st, 
                 prefix = glue("{out_dir}"),
                 plot_item = "SpotLight_Anno",
                 tmp_color = c("#238b45", "#ce1256"),
                 plot_order = c("B",  "T"),
                 name = glue("T_B"))

hcc_st@meta.data$TLS <- hcc_st@meta.data$TLS_raw
hcc_st@meta.data$TLS_maturity[hcc_st@meta.data$TLS_maturity == "Conforming"] <- "conforming"
hcc_st@meta.data$TLS_maturity[hcc_st@meta.data$TLS_maturity == "Deviating"] <- "deviating"
hcc_st@meta.data$TLS_maturity[hcc_st@meta.data$TLS_maturity == "NotMature"] <- "deviating"
hcc_st@meta.data$TLS_maturity[hcc_st@meta.data$TLS_maturity == "Rare_Naive_GCB"] <- "ST2931T_NA"

hcc_st@meta.data[!hcc_st@meta.data$TLS_maturity %in%  c("Mature", "deviating","conforming"), ]$TLS <- "ST2931T_NA"
hcc_st@meta.data[hcc_st@meta.data$TLS_maturity %in%  c("Mature", "deviating","conforming"), ]$TLS <- "ST2931T_TLS"


SpatialDimPython(obj = hcc_st,
                 prefix = glue("{outdir}"),
                 plot_item = "TLS",
                 plot_order = c("ST2931T_TLS"),
                 tmp_color = c("#9CCCE8"),
                 name = glue("TLS"))

#### Figure E2 H ====
cur_sampleid <- "2907T"
out_dir <- glue("~")
hcc_st <- readRDS(glue('~/{cur_sampleid}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
hcc_st@meta.data <-fread(glue("~/{cur_sampleid}.txt", sep = "\t",header = T)) %>% as.data.frame()
rownames(hcc_st@meta.data) <- sapply(hcc_st@meta.data$BinID, function(x)strsplit(x, "_")[[1]][2]) %>% as.character()
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$CellSubType == "Tumor"] <- "Tumor"
hcc_st@meta.data <- as.data.frame(hcc_st@meta.data)
hcc_st@meta.data$BinID <- rownames(hcc_st@meta.data)
rownames(hcc_st@meta.data) <- hcc_st@meta.data$BinID
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$SpotLight_Anno == "pDC"] <- "Myeloid"
hcc_st@meta.data$TLS_maturity[!hcc_st@meta.data$TLS_maturity %in% c("Mature", "Conforming", "Deviating" )] <- glue("ST{cur_sampleid}_NA")
hcc_st@meta.data$TLS_maturity[grepl("NA", hcc_st@meta.data$TLS_raw)] <- glue("ST{cur_sampleid}_NA")

# sub_obj <- subset(hcc_st, SpotLight_Anno %in% c("B", "T",  "Endothelial"))
# sub_obj$SpotLight_Anno <- factor(sub_obj$SpotLight_Anno, levels = c("B", "T",  "Endothelial"))
sub_obj <- hcc_st
sub_obj$SpotLight_Anno <- factor(sub_obj$SpotLight_Anno, levels = c("B", "T","Myeloid", "Endothelial", "Plasma", "Fibroblast","Tumor&Hepatocyte") %>% rev())
chemokines <- c("CXCL13", "CXCR5","CCL19", "CXCR4", "CCR7",  "CCL21", "CXCL12")
pdf(glue("{out_dir}/chk.all.dot.pdf"), w=5, h=4)
p_heat <- DotPlot(sub_obj, features = chemokines,
                  group.by='SpotLight_Anno', cols='Spectral') +
  # coord_flip() + 
  geom_point(aes(size=pct.exp),
             shape = 21,
             colour="black", stroke=0.5) +
  # guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  # scale_size(breaks = c(0, 20, 40, 60, 80)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.text.y = element_text(size=8))
print(p_heat)
dev.off()

#### Figure E2 I ====
cur_sampleid <- "2907T"
out_dir <- glue("~")
hcc_st <- readRDS(glue('~/{cur_sampleid}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
hcc_st@meta.data <-fread(glue("~/{cur_sampleid}.txt", sep = "\t",header = T)) %>% as.data.frame()
rownames(hcc_st@meta.data) <- sapply(hcc_st@meta.data$BinID, function(x)strsplit(x, "_")[[1]][2]) %>% as.character()
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$CellSubType == "Tumor"] <- "Tumor"
hcc_st@meta.data <- as.data.frame(hcc_st@meta.data)
hcc_st@meta.data$BinID <- rownames(hcc_st@meta.data)
rownames(hcc_st@meta.data) <- hcc_st@meta.data$BinID
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$SpotLight_Anno == "pDC"] <- "Myeloid"
hcc_st@meta.data$TLS_maturity[!hcc_st@meta.data$TLS_maturity %in% c("Mature", "Conforming", "Deviating" )] <- glue("ST{cur_sampleid}_NA")
hcc_st@meta.data$TLS_maturity[grepl("NA", hcc_st@meta.data$TLS_raw)] <- glue("ST{cur_sampleid}_NA")

# sub_obj <- subset(hcc_st, SpotLight_Anno %in% c("B", "T",  "Endothelial"))
sub_obj <- hcc_st
sub_obj$SpotLight_Anno <- factor(sub_obj$SpotLight_Anno, levels = c("B", "T","Myeloid", "Endothelial", "Plasma", "Fibroblast","Tumor&Hepatocyte") %>% rev())
chk_count <- sub_obj@assays$Spatial@counts[chemokines, ] %>% t() %>% as.data.frame()
chk_count$binid <- rownames(chk_count)

pacbio_merge <- readRDS("~/merged_bcr.rds")
sub_smrt_bcr <- pacbio_merge[pacbio_merge$sample == cur_sampleid,]
sub_smrt_bcr_counts <- dplyr::group_by(sub_smrt_bcr, Bin50, locus) %>% 
  dplyr::summarise(count = n()) %>% 
  tidyr::pivot_wider(names_from = locus, values_from = count, values_fill = 0)
colnames(sub_smrt_bcr_counts)[1] <- "binid"
chk_count <- dplyr::left_join(chk_count, sub_smrt_bcr_counts)
chk_count$IGH[is.na(chk_count$IGH)] <- 0
chk_count$IGK[is.na(chk_count$IGK)] <- 0
chk_count$IGL[is.na(chk_count$IGL)] <- 0

# sub_smrt_bcr_counts <- table(sub_smrt_bcr$Bin50, sub_smrt_bcr$locus) %>% as.data.frame()

Mixcr_TCR_ifo <- readRDS("~/merged_tcr.rds")
sub_smrt_tcr <- Mixcr_TCR_ifo[Mixcr_TCR_ifo$sample == cur_sampleid,]
sub_smrt_tcr <- sub_smrt_tcr[sub_smrt_tcr$topChains %in% c("TRA", "TRB"), ]
sub_smrt_tcr_counts <- dplyr::group_by(sub_smrt_tcr, Bin50_BinID, topChains) %>% 
  dplyr::summarise(count = n()) %>% 
  tidyr::pivot_wider(names_from = topChains, values_from = count, values_fill = 0)
colnames(sub_smrt_tcr_counts)[1] <- "binid"
chk_count <- dplyr::left_join(chk_count, sub_smrt_tcr_counts)
chk_count$TRA[is.na(chk_count$TRA)] <- 0
chk_count$TRB[is.na(chk_count$TRB)] <- 0

chains <- c("IGH", "IGK", "IGL", "TRA", "TRB")
chain.chemokine.corr <- data.frame("chain" = chains)
chemokines <- c("CXCR5","CCL19", "CXCR4",  "CCL21", "CCR7", "CXCL13", "CXCL12")

# Run through combinations antigen receptor chains and chemokines
for(chemokine in chemokines) {
  corr.vector <- vector()
  for(chain in chains) {
    sub_chk_count <- chk_count[chk_count[,chemokine ]!=0 & chk_count[,chain ]!=0, ]
    # v1 <- chk_count[,chemokine ]
    
    cor1 <-  cor(sub_chk_count[,chemokine ], sub_chk_count[,chain ], method = "spearman")
    corr.vector <- c(corr.vector, cor1)
    
  }
  chain.chemokine.corr <- cbind(chain.chemokine.corr, corr.vector)
}

rownames(chain.chemokine.corr) <- chain.chemokine.corr$chain
chain.chemokine.corr <- chain.chemokine.corr[,-1]
colnames(chain.chemokine.corr) <- chemokines


d <- dist(chain.chemokine.corr %>% t(), method = "euclidean")
fit_c <- hclust(d , method="ward.D") #进行Ward层次聚类
d <- dist(chain.chemokine.corr, method = "euclidean")
fit_r <- hclust(d , method="ward.D") #进行Ward层次聚类
chain.chemokine.corr.c <- chain.chemokine.corr[, fit_c$order]



library(corrplot)
corr.color <- rev(COL2("RdBu"))
pdf(glue("{out_dir}/chemo_chain_heatmap_LR.pdf"), width = 6.5, height = 3.5)
corrplot(as.matrix(chain.chemokine.corr), col = corr.color,  tl.col="black", method = 'color', tl.pos = "lt")
dev.off()

#### Figure E2 J ====
bcr_merge <- readRDS("~/merged_bcr.rds")
outdir <- "~/tumor/"
bcr_heavy <- subset(bcr_merge, bcr_merge$locus == "IGH")
bcr_light <- subset(bcr_merge, bcr_merge$locus != "IGH")
sample_count <- dplyr::count(bcr_heavy, sample, sort = T)
top5 <- c("2907T", "2829T", "2879T", "2931T", "3154T")

summary_tumor <- list()
for (x in top5) {
  obj <- readRDS(glue('~/{x}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
  obj@meta.data <-read.table(glue("~/{x}.txt"),
                             sep = "\t",
                             header = T,
                             row.names = 1) %>% as.data.frame()
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  
  tmp_heavy <- subset(bcr_heavy, sample == x)
  tmp_heavy <- tmp_heavy$clone_id %>% table() %>% as.data.frame()
  cloneid <- tmp_heavy[order(-tmp_heavy$Freq),][1,1] %>% as.character()
  tmp_bcr_cloneid <- subset(bcr_heavy, bcr_heavy$sample == x & bcr_heavy$clone_id == cloneid)
  tmp_bcr_cloneid <- tmp_bcr_cloneid$Bin50 %>% table() %>% as.data.frame()
  colnames(tmp_bcr_cloneid) <- c("bin_id","clone")
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = "bin_id") %>%
    column_to_rownames(var="bin_id")
  obj@meta.data$clone <- ifelse(is.na(obj@meta.data$clone), 0, obj@meta.data$clone)
  obj <- subset(obj, tissue == 1)
  obj@meta.data$clone <- ifelse(obj@meta.data$clone != 0, "most_abundant_clone", obj@meta.data$Bin_Region)
  obj@meta.data$clone <- ifelse(obj@meta.data$clone == "Tumor_capsule", "Invasive_zone", obj@meta.data$clone)
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_side_of_Margin_area"] <- "Tumor"
  obj@meta.data$clone[obj@meta.data$clone == "Paratumor_side_of_Margin_area"] <- "Paratumor"
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_tissue"] <- "Tumor"
  
  SpatialDimPython(obj = obj,
                   prefix = glue("{outdir}"),
                   plot_item = "clone",
                   plot_order = c("most_abundant_clone","Tumor"),
                   tmp_color = c("red","#E9ECEB"),
                   # "#E9ECEB","#FFE8BF","#F9F7EE"
                   name = glue("{x}_{cloneid}_dimplot"))
  tmp_bcr_cloneid <- subset(bcr_heavy, bcr_heavy$sample == x & bcr_heavy$clone_id == cloneid)
  
  summary_tumor[[paste0(x, 1)]] <- data.frame(
    sample = x,
    cloneid = cloneid,
    cdr3aa = tmp_bcr_cloneid$junction_aa[1],
    isotype = tmp_bcr_cloneid$c_call[1],
    locus = "heavy",
    n = nrow(tmp_bcr_cloneid))
  
  obj <- readRDS(glue('~/{x}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
  obj@meta.data <-read.table(glue("~/{x}.txt"),
                             sep = "\t",
                             header = T,
                             row.names = 1) %>% as.data.frame()
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  tmp_light <- subset(bcr_light, sample == x)
  tmp_light <- tmp_light$clone_id %>% table() %>% as.data.frame()
  cloneid <- tmp_light[order(-tmp_light$Freq),][1,1] %>% as.character()
  tmp_bcr_cloneid <- subset(bcr_light, bcr_light$sample == x & bcr_light$clone_id == cloneid)
  tmp_bcr_cloneid <- tmp_bcr_cloneid$Bin50 %>% table() %>% as.data.frame()
  colnames(tmp_bcr_cloneid) <- c("bin_id","clone")
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = "bin_id") %>%
    column_to_rownames(var="bin_id")
  obj@meta.data$clone <- ifelse(is.na(obj@meta.data$clone), 0, obj@meta.data$clone)
  obj <- subset(obj, tissue == 1)
  obj@meta.data$clone <- ifelse(obj@meta.data$clone != 0, "most_abundant_clone", obj@meta.data$Bin_Region)
  obj@meta.data$clone <- ifelse(obj@meta.data$clone == "Tumor_capsule", "Invasive_zone", obj@meta.data$clone)
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_side_of_Margin_area"] <- "Tumor"
  obj@meta.data$clone[obj@meta.data$clone == "Paratumor_side_of_Margin_area"] <- "Paratumor"
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_tissue"] <- "Tumor"
  
  SpatialDimPython(obj = obj,
                   prefix = glue("{outdir}"),
                   plot_item = "clone",
                   plot_order = c("most_abundant_clone","Tumor"),
                   tmp_color = c("red","#E9ECEB"),
                   name = glue("{x}_{cloneid}_dimplot"))
  tmp_bcr_cloneid <- subset(bcr_light, bcr_light$sample == x & bcr_light$clone_id == cloneid)
  
  summary_tumor[[paste0(x, 2)]] <- data.frame(
    sample = x,
    cloneid = cloneid,
    cdr3aa = tmp_bcr_cloneid$junction_aa[1],
    isotype = tmp_bcr_cloneid$c_call[1],
    locus = "light",
    n = nrow(tmp_bcr_cloneid))
}
summary_tumor <- dplyr::bind_rows(summary_tumor)
save_table(summary_tumor, glue("{outdir}/summary_tumor"), if_csv = T)

bcr_merge <- readRDS("~/merged_bcr.rds")
outdir <- "~/border_v2/"

bcr_heavy <- subset(bcr_merge, bcr_merge$locus == "IGH")
bcr_light <- subset(bcr_merge, bcr_merge$locus != "IGH")
sample_count <- dplyr::count(bcr_heavy, sample, sort = T)
top5 <- c("2931B", "RH05B", "RH16B", "2903B", "2837B")

summary_border <- list()
for (x in top5) {
  obj <- readRDS(glue('~/{x}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
  obj@meta.data <-read.table(glue("~/{x}.txt"),
                             sep = "\t",
                             header = T,
                             row.names = 1) %>% as.data.frame()
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  
  tmp_heavy <- subset(bcr_heavy, sample == x)
  tmp_heavy <- tmp_heavy$clone_id %>% table() %>% as.data.frame()
  cloneid <- tmp_heavy[order(-tmp_heavy$Freq),][1,1] %>% as.character()
  tmp_bcr_cloneid <- subset(bcr_heavy, bcr_heavy$sample == x & bcr_heavy$clone_id == cloneid)
  tmp_bcr_cloneid <- tmp_bcr_cloneid$Bin50 %>% table() %>% as.data.frame()
  colnames(tmp_bcr_cloneid) <- c("bin_id","clone")
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = "bin_id") %>%
    column_to_rownames(var="bin_id")
  obj@meta.data$clone <- ifelse(is.na(obj@meta.data$clone), 0, obj@meta.data$clone)
  obj <- subset(obj, tissue == 1)
  # SpatialFeaturePython(obj = obj, 
  #                      prefix = glue("{outdir}"),
  #                      if_gene = F,
  #                      plot_item = "clone",
  #                      name = glue("{x}_{cloneid}_featureplot")) 
  obj@meta.data$clone <- ifelse(obj@meta.data$clone != 0, "most_abundant_clone", obj@meta.data$Bin_Region)
  obj@meta.data$clone <- ifelse(obj@meta.data$clone == "Tumor_capsule", "Invasive_zone", obj@meta.data$clone)
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_side_of_Margin_area"] <- "Tumor"
  obj@meta.data$clone[obj@meta.data$clone == "Paratumor_side_of_Margin_area"] <- "Paratumor"
  obj$TLS <- obj$Bin_Region
  obj$TLS[obj$TLS != "Tumor_capsule"] <- "_NA"
  SpatialDimPython(obj = obj,
                   prefix = glue("{outdir}"),
                   plot_item = "clone",
                   if_tls_contour = "T",
                   plot_order = c("most_abundant_clone","Tumor","Invasive_zone","Paratumor"),
                   tmp_color = c("red","#E9ECEB","#FFE8BF","#F9F7EE"),
                   name = glue("{x}_{cloneid}_v2"))
  tmp_bcr_cloneid <- subset(bcr_heavy, bcr_heavy$sample == x & bcr_heavy$clone_id == cloneid)
  
  summary_border[[paste0(x, 1)]] <- data.frame(
    sample = x,
    cloneid = cloneid,
    cdr3aa = tmp_bcr_cloneid$junction_aa[1],
    isotype = tmp_bcr_cloneid$c_call[1],
    locus = "heavy",
    n = nrow(tmp_bcr_cloneid))
  
  obj <- readRDS(glue('~/{x}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
  obj@meta.data <-read.table(glue("~/{x}.txt"),
                             sep = "\t",
                             header = T,
                             row.names = 1) %>% as.data.frame()
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  tmp_light <- subset(bcr_light, sample == x)
  tmp_light <- tmp_light$clone_id %>% table() %>% as.data.frame()
  cloneid <- tmp_light[order(-tmp_light$Freq),][1,1] %>% as.character()
  tmp_bcr_cloneid <- subset(bcr_light, bcr_light$sample == x & bcr_light$clone_id == cloneid)
  tmp_bcr_cloneid <- tmp_bcr_cloneid$Bin50 %>% table() %>% as.data.frame()
  colnames(tmp_bcr_cloneid) <- c("bin_id","clone")
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = "bin_id") %>%
    column_to_rownames(var="bin_id")
  obj@meta.data$clone <- ifelse(is.na(obj@meta.data$clone), 0, obj@meta.data$clone)
  obj <- subset(obj, tissue == 1)
  # SpatialFeaturePython(obj = obj, 
  #                      prefix = glue("{outdir}"),
  #                      if_gene = F,
  #                      plot_item = "clone",
  #                      name = glue("{x}_{cloneid}_featureplot")) 
  obj@meta.data$clone <- ifelse(obj@meta.data$clone != 0, "most_abundant_clone", obj@meta.data$Bin_Region)
  obj@meta.data$clone <- ifelse(obj@meta.data$clone == "Tumor_capsule", "Invasive_zone", obj@meta.data$clone)
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_side_of_Margin_area"] <- "Tumor"
  obj@meta.data$clone[obj@meta.data$clone == "Paratumor_side_of_Margin_area"] <- "Paratumor"
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_tissue"] <- "Tumor"
  obj$TLS <- obj$Bin_Region
  obj$TLS[obj$TLS != "Tumor_capsule"] <- "_NA"
  SpatialDimPython(obj = obj,
                   prefix = glue("{outdir}"),
                   plot_item = "clone",
                   if_tls_contour = "T",
                   plot_order = c("most_abundant_clone","Tumor","Invasive_zone","Paratumor"),
                   tmp_color = c("red","#E9ECEB","#FFE8BF","#F9F7EE"),
                   name = glue("{x}_{cloneid}_v2"))
  tmp_bcr_cloneid <- subset(bcr_light, bcr_light$sample == x & bcr_light$clone_id == cloneid)
  
  summary_border[[paste0(x, 2)]] <- data.frame(
    sample = x,
    cloneid = cloneid,
    cdr3aa = tmp_bcr_cloneid$junction_aa[1],
    isotype = tmp_bcr_cloneid$c_call[1],
    locus = "light",
    n = nrow(tmp_bcr_cloneid))
}
summary_border <- dplyr::bind_rows(summary_border)
save_table(summary_border, 
           glue("{outdir}/summary_border"), if_csv = T)

bcr_merge <- readRDS("~/merged_bcr.rds")
outdir <- "~/paratumor/"

bcr_heavy <- subset(bcr_merge, bcr_merge$locus == "IGH")
bcr_light <- subset(bcr_merge, bcr_merge$locus != "IGH")
sample_count <- dplyr::count(bcr_heavy, sample, sort = T)
top5 <- c("2903P", "2740P", "2931P", "3154P", "2896P")

summary_paratumor <- list()
for (x in top5) {
  obj <- readRDS(glue('~/{x}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
  obj@meta.data <-read.table(glue("~/{x}.txt"),
                             sep = "\t",
                             header = T,
                             row.names = 1) %>% as.data.frame()
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  
  tmp_heavy <- subset(bcr_heavy, sample == x)
  tmp_heavy <- tmp_heavy$clone_id %>% table() %>% as.data.frame()
  cloneid <- tmp_heavy[order(-tmp_heavy$Freq),][1,1] %>% as.character()
  tmp_bcr_cloneid <- subset(bcr_heavy, bcr_heavy$sample == x & bcr_heavy$clone_id == cloneid)
  tmp_bcr_cloneid <- tmp_bcr_cloneid$Bin50 %>% table() %>% as.data.frame()
  colnames(tmp_bcr_cloneid) <- c("bin_id","clone")
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = "bin_id") %>%
    column_to_rownames(var="bin_id")
  obj@meta.data$clone <- ifelse(is.na(obj@meta.data$clone), 0, obj@meta.data$clone)
  obj <- subset(obj, tissue == 1)
  # SpatialFeaturePython(obj = obj, 
  #                      prefix = glue("{outdir}"),
  #                      if_gene = F,
  #                      plot_item = "clone",
  #                      name = glue("{x}_{cloneid}_featureplot")) 
  obj@meta.data$clone <- ifelse(obj@meta.data$clone != 0, "most_abundant_clone", obj@meta.data$Bin_Region)
  obj@meta.data$clone <- ifelse(obj@meta.data$clone == "Tumor_capsule", "Invasive_zone", obj@meta.data$clone)
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_side_of_Margin_area"] <- "Tumor"
  obj@meta.data$clone[obj@meta.data$clone == "Paratumor_side_of_Margin_area"] <- "Paratumor"
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_tissue"] <- "Tumor"
  obj@meta.data$clone[obj@meta.data$clone == "Paratumor_tissue"] <- "Paratumor"
  
  SpatialDimPython(obj = obj,
                   prefix = glue("{outdir}"),
                   plot_item = "clone",
                   plot_order = c("most_abundant_clone","Paratumor"),
                   tmp_color = c("red","#F9F7EE"),
                   name = glue("{x}_{cloneid}_dimplot"))
  tmp_bcr_cloneid <- subset(bcr_heavy, bcr_heavy$sample == x & bcr_heavy$clone_id == cloneid)
  
  summary_paratumor[[paste0(x, 1)]] <- data.frame(
    sample = x,
    cloneid = cloneid,
    cdr3aa = tmp_bcr_cloneid$junction_aa[1],
    isotype = tmp_bcr_cloneid$c_call[1],
    locus = "heavy",
    n = nrow(tmp_bcr_cloneid))
  
  obj <- readRDS(glue('~/{x}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
  obj@meta.data <-read.table(glue("~/{x}.txt"),
                             sep = "\t",
                             header = T,
                             row.names = 1) %>% as.data.frame()
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  tmp_light <- subset(bcr_light, sample == x)
  tmp_light <- tmp_light$clone_id %>% table() %>% as.data.frame()
  cloneid <- tmp_light[order(-tmp_light$Freq),][1,1] %>% as.character()
  tmp_bcr_cloneid <- subset(bcr_light, bcr_light$sample == x & bcr_light$clone_id == cloneid)
  tmp_bcr_cloneid <- tmp_bcr_cloneid$Bin50 %>% table() %>% as.data.frame()
  colnames(tmp_bcr_cloneid) <- c("bin_id","clone")
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = "bin_id") %>%
    column_to_rownames(var="bin_id")
  obj@meta.data$clone <- ifelse(is.na(obj@meta.data$clone), 0, obj@meta.data$clone)
  obj <- subset(obj, tissue == 1)
  # SpatialFeaturePython(obj = obj, 
  #                      prefix = glue("{outdir}"),
  #                      if_gene = F,
  #                      plot_item = "clone",
  #                      name = glue("{x}_{cloneid}_featureplot")) 
  obj@meta.data$clone <- ifelse(obj@meta.data$clone != 0, "most_abundant_clone", obj@meta.data$Bin_Region)
  obj@meta.data$clone <- ifelse(obj@meta.data$clone == "Tumor_capsule", "Invasive_zone", obj@meta.data$clone)
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_side_of_Margin_area"] <- "Tumor"
  obj@meta.data$clone[obj@meta.data$clone == "Paratumor_side_of_Margin_area"] <- "Paratumor"
  obj@meta.data$clone[obj@meta.data$clone == "Tumor_tissue"] <- "Tumor"
  obj@meta.data$clone[obj@meta.data$clone == "Paratumor_tissue"] <- "Paratumor"  
  SpatialDimPython(obj = obj,
                   prefix = glue("{outdir}"),
                   plot_item = "clone",
                   plot_order = c("most_abundant_clone","Paratumor"),
                   tmp_color = c("red","#F9F7EE"),
                   name = glue("{x}_{cloneid}_dimplot"))
  tmp_bcr_cloneid <- subset(bcr_light, bcr_light$sample == x & bcr_light$clone_id == cloneid)
  
  summary_paratumor[[paste0(x, 2)]] <- data.frame(
    sample = x,
    cloneid = cloneid,
    cdr3aa = tmp_bcr_cloneid$junction_aa[1],
    isotype = tmp_bcr_cloneid$c_call[1],
    locus = "light",
    n = nrow(tmp_bcr_cloneid))
}
summary_paratumor <- dplyr::bind_rows(summary_paratumor)
save_table(summary_paratumor,
           glue("{outdir}/summary_paratumor"), if_csv = T)


