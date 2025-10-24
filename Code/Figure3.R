library(glue)
library(readxl)
library(readr)
library(dplyr)
library(stringr)



#### get tree ====

## 抗体mapping 45个样本  ====
pacbio_merge <-  readRDS("~/45sample_bcr0626.rds")
sample_45 <-  pacbio_merge$sample %>% unique()
seq_info_merge <- list()
pacbio_sc_info <- fread("~/Pacbio_SingleCell_43.csv")
sample_45[!sample_45 %in% pacbio_sc_info$SampleID]
pacbio_sc_info$SampleID[!pacbio_sc_info$SampleID %in% sample_45]

# pacbio_merge <- dplyr::bind_rows(pacbio_merge)
out_dir <- "~/raw_fa"
unlink(glue('{out_dir}/merge.fa'))
all_seq_id <- c()
for(cur_smrt_sample in sample_45){
  cur_patient <- substr(cur_smrt_sample, 1, 4)
  cur_smrt <- pacbio_merge[pacbio_merge$sample  %in%  cur_smrt_sample ,]
  seq_info <- list()
  
  
  dir.create(glue("{out_dir}"), recursive = T)
  
  seq_smrt <- cur_smrt$sequence
  names(seq_smrt) <- cur_smrt$sequence_id
  lapply(names(seq_smrt), function(x){
    write.table(
      glue(">{x}"),
      glue('{out_dir}/merge.fa'),
      row.names = F,
      quote = F,
      col.names = F,
      append = TRUE
    )
    write.table(seq_smrt[x],
                glue('{out_dir}/merge.fa'),
                row.names = F,
                quote = F,
                col.names = F,
                append = TRUE)
  })
  
  info <- data.frame(
    id = cur_smrt$sequence_id
  )
  info$sample <- cur_smrt_sample
  info$type <- "smrt"
  seq_info_merge[[glue("{cur_smrt_sample}_smrt")]] <- info[, c("id", "sample", "type")]
  
  
  idx <- which(pacbio_sc_info$SampleID %in%  cur_smrt_sample)
  cur_sc_list <- pacbio_sc_info$SC_Path[idx[1]]
  cur_sc_list <- strsplit(cur_sc_list, ";")[[1]]
  cur_sc_bind <- list()
  
  for(cur_path in cur_sc_list){
    subname <- basename(cur_path)
    
    if(!"filtered_contig_annotations.csv" %in% list.files(cur_path)){
      print(glue("{subname} has not filtered_contig_annotations.csv !!!!!!!!!!!!!!!!!"))
      next 
    }
    
    info <- fread(glue("{cur_path}/filtered_contig_annotations.csv"))
    fa_file <- read.table(glue("{cur_path}/filtered_contig.fasta"), col.names = F)
    fa_file[,1] <- gsub(">","",fa_file[,1] )
    
    contig_name_idx <- seq(1, nrow(fa_file), by = 2)
    seq_idx <- seq(2, nrow(fa_file), by = 2)
    seq <- fa_file[seq_idx,1]
    seq_id <- paste(subname, fa_file[contig_name_idx,1], sep = "_")
    names(seq) <- seq_id
    
    if(any(seq_id %in% all_seq_id)){
      print(glue("{cur_patient}, {cur_path}"))
      next
    }
    all_seq_id <- c(all_seq_id, seq_id)
    
    
    
    
    if(identical(info$contig_id, fa_file[contig_name_idx,1])){
      lapply(names(seq), function(x){
        write.table(
          glue(">{x}"),
          glue('{out_dir}/merge.fa'),
          row.names = F,
          quote = F,
          col.names = F,
          append = TRUE
        )
        write.table(seq[x],
                    glue('{out_dir}/merge.fa'),
                    row.names = F,
                    quote = F,
                    col.names = F,
                    append = TRUE)
      })
    }else{
      warning(glue("{cur_path} not identical contig"))
    }
    
    info$id <- seq_id
    info$sample <- subname
    info$type <- "singlecell"
    if(paste("SC",subname, sep = "_") %in% names(seq_info_merge)){
      print(glue("{cur_smrt_sample} {subname}"))
    }
    seq_info_merge[[paste("SC",subname, sep = "_")]] <- info[, c("id", "sample", "type")]
  }
  
  # seq_info_bind <- dplyr::bind_rows(seq_info)
  # save_table(seq_info_bind, glue("{out_dir}/{cur_smrt_sample}_info"))
  
  # seq_info_merge[[cur_smrt_sample]] <- seq_info_bind
  
  
}
seq_info_merge <- dplyr::bind_rows(seq_info_merge)
save_table(seq_info_merge, glue("~/tls_spatial_tree/merge_info"))

# sh ~/run_changeo_new.sh RH05 ~/raw_fa/RH05.fa ~/changeo  ~/blastn

pacbio_merge <-  readRDS("~/45sample_bcr0626.rds")
sample_45 <-  pacbio_merge$sample %>% unique()
write.table(data.frame(x = sample_45), "~/sample_45.txt",col.names = F,
            row.names = F, quote = F, sep = "\t")
# changeo_merge <- dplyr::bind_rows(changeo_merge)
changeo_merge <- list()
bcr_smrt_merge <- readRDS("~/45sample_bcr0626.rds")

for(t in sample_45){
  cur_changeo <- fread(glue("~/{t}_family_07_clone-pass.tsv"))
  cur_info <- fread(glue("~/{t}_info.txt"))
  cur_changeo$id <- cur_changeo$sequence_id
  cur_changeo$patient <- t
  cur_changeo$clone_id <- paste(t, cur_changeo$clone_id, sep = "_")
  cur_changeo$clone_family <- paste(t, cur_changeo$clone_family, sep = "_")
  cur_changeo <- dplyr::left_join(cur_changeo, cur_info)
  changeo_merge[[t]] <- cur_changeo
}
changeo_merge <- dplyr::bind_rows(changeo_merge)
changeo_merge <- get_shared_info(changeo_merge, target_list =  c("type"), source = "clone_id")
changeo_merge <- dplyr::rename(changeo_merge, share_label_clone = type_shared_main,
                               share_sublabel_clone = type_shared_sub)
table(changeo_merge$share_label_clone)
table(changeo_merge$type)
saveRDS(changeo_merge, "~/changeo_merge_07.rds")


changeo_merge <- readRDS("~/changeo_merge_07.rds")
test <- changeo_merge[changeo_merge$sample == "2907T",]
test$sequence_id[test$type == "singlecell"][1]
changeo_merge$clone_family_07 <- paste(changeo_merge$patient, changeo_merge$clone_family_07, sep = "_")
ab_info <- fread("~/changeo/final_genescript2.csv")
ab_tp53 <- c(10, 14, 18, 22, 28, 38, 39, 41, 43, 49, 55, 58, 59, 67, 72, 77, 79, 81, 86, 88, 116, 117, 118, 119, 121, 125)
ab_ifi30 <- c(10, 14, 18, 22, 23, 28, 36, 38, 39, 41, 43, 49, 55, 58, 59, 61, 77, 79, 81, 116, 117, 118, 119, 120, 121, 124, 125, 5, 30, 67, 72, 86, 88)
ab_2 <- intersect(ab_tp53, ab_ifi30)
ab_tp53_u <- setdiff(ab_tp53, ab_2)

ab_info$tp53 <- "n"
ab_info <- as.data.frame(ab_info)
ab_info$tp53[(1:nrow(ab_info)) %in% ab_tp53] <- "y"

ab_info$ifi30 <- "n"
ab_info <- as.data.frame(ab_info)
ab_info$ifi30[(1:nrow(ab_info)) %in% ab_ifi30] <- "y"

ab_info$raw_sc_id[1]
ab_info$sequence_id <- str_replace(ab_info$raw_sc_id,"-","_")
ab_info$id <- "n"
ab_info$no <- 1:nrow(ab_info)


ab_info$id[ab_info$ifi30 == "y"] <-  paste("Ig_", ab_info$no[ab_info$ifi30 == "y"], sep = "")
ab_info$id[ab_info$tp53 == "y"] <-  paste("Ig_", ab_info$no[ab_info$tp53 == "y"], sep = "")
# ab_info$id <- paste("Ig_", ab_info$id, sep = "")
# test <- changeo_merge[changeo_merge$sequence_id %in% ab_info$sequence_id,]
# table(test$sequence_id)
# nrow(test)
ab_ifi30_seq <- ab_info$sequence_id[ab_info$ifi30 == "y"]
ab_ifi30_clone <- changeo_merge[changeo_merge$sequence_id %in% ab_ifi30_seq,]
ab_ifi30_clone <- ab_ifi30_clone[ab_ifi30_clone$share_label_clone == "Shared",]
ab_ifi30_cloneid <- ab_ifi30_clone$clone_id %>% unique()
saveRDS(ab_ifi30_clone, "~/ab_ifi30_clone.rds")

ab_tp53_seq <- ab_info$sequence_id[ab_info$tp53 == "y"]
ab_tp53_clone <- changeo_merge[changeo_merge$sequence_id %in% ab_tp53_seq,]
ab_tp53_clone <- ab_tp53_clone[ab_tp53_clone$share_label_clone == "Shared",]
ab_tp53_cloneid <- ab_tp53_clone$clone_id %>% unique()

saveRDS(ab_tp53_clone, "~/ab_tp53_clone.rds")

ab_merge_seq <- ab_info$sequence_id[ab_info$id != "n"]
ab_merge_clone <- changeo_merge[changeo_merge$sequence_id %in% ab_merge_seq,]
ab_merge_clone <- ab_merge_clone[ab_merge_clone$share_label_clone == "Shared",]
ab_merge_clone <- ab_merge_clone %>% dplyr::select(-id)
ab_merge_clone <- dplyr::left_join(ab_merge_clone, ab_info[, c("sequence_id",
                                                               "ifi30", "tp53", "no",
                                                               "id")])
ab_merge_cloneid <- ab_merge_clone$clone_id %>% unique()
# saveRDS(ab_merge_clone, "~/ab_merge_clone.rds")
# save_table(ab_merge_clone, "~/ab_merge_clone", if_csv = F)


ab_merge_clonefam_id <- changeo_merge[changeo_merge$clone_id %in% ab_merge_cloneid,]$clone_family_07 %>% unique()
ab_merge_clonefam <- changeo_merge[changeo_merge$clone_family_07 %in% ab_merge_clonefam_id,]
ab_merge_clonefam$isotype <- substr(ab_merge_clonefam$c_call, 1, 5)
ab_merge_clonefam$isotype[is.na(ab_merge_clonefam$isotype)] <- "IGHG1"
ab_merge_clone$clone_subid <- ab_merge_clone$id
ab_merge_clone %>% 
  group_by(clone_id) %>% 
  summarise(clone_subid = paste(sort(unique(id)), collapse = ",")) -> ab_merge_clone_subid

ab_merge_clonefam <- dplyr::left_join(ab_merge_clonefam, ab_merge_clone_subid)
# ab_merge_clonefam$p_cloneid <- paste(ab_merge_clonefam$patient, ab_merge_clonefam$clone_id, sep = "_")
ab_merge_clonefam$clone_subid[is.na(ab_merge_clonefam$clone_subid)] <- ab_merge_clonefam$clone_id[is.na(ab_merge_clonefam$clone_subid)]
ab_merge_clonefam$isotype_group <- substr(ab_merge_clonefam$isotype, 1, 4)

# saveRDS(hbcab_map_clonefam, "~/talent/hbcab_map_clonefam.rds")
# hbcab_map_clonefam <- readRDS("~/talent/hbcab_map_clonefam.rds")
ab_merge_clonefam <- dplyr::left_join(ab_merge_clonefam, bcr_smrt_merge[, c("sequence_id",
                                                                            "SpotLight_Anno", "TLS_maturity",
                                                                            "TLS_raw", "CellSubType")])
ab_merge_clonefam$TLS_maturity[is.na(ab_merge_clonefam$TLS_maturity)] <- ab_merge_clonefam$type[is.na(ab_merge_clonefam$TLS_maturity)]
ab_merge_clonefam$SpotLight_Anno[is.na(ab_merge_clonefam$SpotLight_Anno)] <- ab_merge_clonefam$type[is.na(ab_merge_clonefam$SpotLight_Anno)]
ab_merge_clonefam$TLS_raw[is.na(ab_merge_clonefam$TLS_raw)] <- ab_merge_clonefam$type[is.na(ab_merge_clonefam$TLS_raw)]

# saveRDS(ab_merge_clonefam, "~/talent/ab_merge_clonefam.rds")
# ab_merge_clonefam <- readRDS("~/talent/ab_merge_clone_v2.rds")
ab_merge_clonefam <- dplyr::left_join(ab_merge_clonefam, bcr_smrt_merge[, c("sequence_id",
                                                                            "Bin50", "Bin20")])
# ab_merge_clonefam$Bin50[is.na(ab_merge_clonefam$Bin50)] <- ab_merge_clonefam$type[is.na(ab_merge_clonefam$Bin50)]
saveRDS(ab_merge_clonefam, "~/ab_merge_clonefam_v2.rds")



#### Figure 3A 2931T ####
library(readxl)

talent_result <- readRDS("~/ab_merge_clonefam_45_v2.rds")
talent_result$mab_name <- ""
outdir <- glue("~/A")
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}
if_tls_contour = "T"
cur_sample <- "2931T"
all_clone <- talent_result$clone_subid[grepl("Ig", talent_result$id)] %>% unique()
summary <- list()
cell_summary <- list()
# for (x in all_clone) {
#   tmp_bcr_cloneid <- subset(talent_result, talent_result$clone_subid == x)
#   mab_name <- tmp_bcr_cloneid$mab_name[grepl("Ig", unique(tmp_bcr_cloneid$mab_name))]
#   mab_name <- paste(mab_name, collapse = ",")
#   talent_result$mab_name[talent_result$clone_subid == x] <- mab_name
# }

t1 = talent_result[!is.na(Bin50)]
tmp_bcr_cloneid <- t1[grepl("Ig", t1$clone_subid),]

obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <- read.table(glue("~/{cur_sample}.txt"), header = T) %>% as.data.frame()
obj$BinID <- substr(obj$BinID, 7, nchar(obj$BinID))
rownames(obj@meta.data) <- obj$BinID

obj@meta.data$bin_id <- obj@meta.data[,1] 
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_side_of_Margin_area"] <- "Tumor"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Paratumor_side_of_Margin_area"] <- "Paratumor"
obj@meta.data$TLS <- obj@meta.data$TLS_raw
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Conforming"] <- "conforming"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Deviating"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "NotMature"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Rare_Naive_GCB"] <- "STRH05B_NA"

obj@meta.data[!obj@meta.data$TLS_maturity %in%  c("Mature", "deviating","conforming"), ]$TLS <- "STRH05B_NA"

all_mab <- unique(tmp_bcr_cloneid$clone_subid)

tmp_bcr_cloneid <- tmp_bcr_cloneid[tmp_bcr_cloneid$type == "smrt" & tmp_bcr_cloneid$sample == cur_sample ,]
tmp_bcr_cloneid <- table(tmp_bcr_cloneid$Bin50, tmp_bcr_cloneid$clone_subid) %>% as.data.frame()
colnames(tmp_bcr_cloneid) <- c("BinID", "mab_name", "clone")
tmp_bcr_cloneid <- tmp_bcr_cloneid[tmp_bcr_cloneid$clone != 0,]
tmp_bcr_cloneid <- dplyr::distinct(tmp_bcr_cloneid, BinID, .keep_all = T)
# tmp_bcr_cloneid$BinID <- paste0(cur_sample, "_", tmp_bcr_cloneid$BinID)

obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = "BinID") %>%
  column_to_rownames(var="BinID")
obj@meta.data$bin_id <- rownames(obj@meta.data)
obj@meta.data$mab_name <- as.character(obj@meta.data$mab_name)
obj@meta.data$mab_name <- ifelse(is.na(obj@meta.data$mab_name), "other", obj@meta.data$mab_name)
obj@meta.data$mab_name <- ifelse(obj@meta.data$TLS != "STRH05B_NA" & obj@meta.data$mab_name == "other", "TLS", obj@meta.data$mab_name)
obj@meta.data$mab_name <- ifelse(obj@meta.data$mab_name == "other", obj@meta.data$Bin_Region, obj@meta.data$mab_name)
obj@meta.data$mab_name <- ifelse(obj@meta.data$mab_name %in% c("Tumor","Invasive_zone"), "Tumor", obj@meta.data$mab_name)

# obj@meta.data$mab_name[obj@meta.data$mab_name %in% c("Ig_10", "Ig_14", "Ig_23", "Ig_30", "Ig_43", "Ig_5", "Ig_55", "Ig_67")] <- "anti-IFI30_antibody"
# obj@meta.data$mab_name[obj@meta.data$mab_name %in% c("Ig_38","Ig_39","Ig_49","Ig_58","Ig_88")] <- "anti-TP53_antibody"

# mab_color <- c("#E13B43", "#369945", "#1979A3", "#663D8F",  "#E0507D", "#EA7B25", "#40BCE4", "#EFA7B9", "#DB896C")
# mab_color <-  c("#33a02c")
# mab_color <-  c("#EE3743")
mab_color <-  c("#EE3743", "#33a02c", "#0B7FAB", "#6a3d9a", "#ff7f00",  "#F6A9BD", "#fdbf6f", "#c51b7d",  
                "#6F6C9E", "#01665e", "#a6cee3", "#b2df8a", "#ffff99", "#bf812d",  "#999999", "#BB7DB2", "#824615", "#ffff33", "#86DBD4",
                "#BFE2E3","#A1CFFA","#78BDAD","#D45651","#397A7F","#F0918E","#EEE8DA","#1F5392","#A0BFAF",
                "#AE98D6","#ECCBDC","#54BAD3","#8b4a4b","#DB896C","#AABAC2","#ffae3b",
                '#CCCCCC','#B5B5B5','#A092C3','#DDA0DD','#03A4C6','#7AC5CD','#A3D9ED','#00E5EE','#F9EBDA','#F5DEB3',
                '#98689E','#E84115','#FFB5C5','#00FF7F','#F9D01C','#B03060','#00ABDC','#D2691E','#03A464','#FF7F00',
                '#8968CD','#1C5B75')
# SpatialDimPython(obj = obj,
#                  prefix = glue("{outdir}"),
#                  plot_item = "mab_name",
#                  plot_order = c("mAb41", "mAb42", "mAb32,mAb40", "mAb43", "mAb34",  "mAb35",  
#                                 "mAb36",  "mAb39", "mAb33",      
#                                  "TLS","Tumor","Paratumor"),
# tmp_color = c(mab_color,"#9CCCE8", "#eaeaea","#EAD2DF"),
#                  name = glue("merge_bin20"))
obj$BinID = obj$bin_id
table(obj$mab_name)
unique(obj$mab_name)

## 外扩一圈
# df2 <- obj@meta.data[obj@meta.data$mab_name == "anti-TP53_antibody",]
# for (i in c("Ig_88", "Ig_58", "Ig_39", "Ig_49", "Ig_38")){
for (i in c("Ig_5", "Ig_10", "Ig_14", "Ig_23", "Ig_30", "Ig_43", "Ig_55", "Ig_67")){
  
  df2 <- obj@meta.data[obj@meta.data$mab_name %in% i,]
  
  df_location <- df2 %>%
    dplyr::mutate(dplyr::across(c(col, row), as.numeric)) %>%
    dplyr::distinct(row, col)
  
  df_expand <- df_location %>%
    st_expand_n(
      n = 1, 
      x = "col", 
      y = "row", 
      expand_n = "expand_n", 
      mode = 2
    )
  df_expand <- df_expand %>%
    dplyr::select(col, row, TLS_center = group.by, distance = expand_n)
  
  
  obj@meta.data <- obj@meta.data %>%
    dplyr::mutate(dplyr::across(c(col, row), as.numeric)) %>%
    dplyr::select(-dplyr::any_of(c("TLS_center", "distance"))) %>%
    dplyr::left_join(df_expand, by = c("row", "col"))
  
  obj@meta.data$mab_name[obj@meta.data$distance %in% c("0", "1")] <- i
}


# cur_sample <- "RH05B"
rownames(obj@meta.data) <- obj$bin_id


## 2931T
rownames(obj@meta.data) <- obj$bin_id
SpatialDimPython(obj = obj,
                 prefix = glue("{outdir}"),
                 plot_item = "mab_name",
                 plot_order = c("Ig_10", "Ig_14", "Ig_23",  "Ig_30", "Ig_43", "Ig_55", "Ig_67",
                                "TLS", "Tumor_tissue"),
                 tmp_color = c(mab_color[2:8],"#9CCCE8", "#eaeaea","#EAD2DF"),
                 name = glue("{cur_sample}_cloneid_bin50"))



#### Figure 3A 2837B ####
library(readxl)

talent_result <- readRDS("~/ab_merge_clonefam_45_v2.rds")
talent_result$mab_name <- ""
outdir <- glue("~/A")
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}
if_tls_contour = "T"
cur_sample <- "2837B"
all_clone <- talent_result$clone_subid[grepl("Ig", talent_result$id)] %>% unique()
summary <- list()
cell_summary <- list()
# for (x in all_clone) {
#   tmp_bcr_cloneid <- subset(talent_result, talent_result$clone_subid == x)
#   mab_name <- tmp_bcr_cloneid$mab_name[grepl("Ig", unique(tmp_bcr_cloneid$mab_name))]
#   mab_name <- paste(mab_name, collapse = ",")
#   talent_result$mab_name[talent_result$clone_subid == x] <- mab_name
# }

t1 = talent_result[!is.na(Bin50)]
tmp_bcr_cloneid <- t1[grepl("Ig", t1$clone_subid),]

obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
obj@meta.data <- read.table(glue("~/{cur_sample}.txt"), header = T) %>% as.data.frame()
obj$BinID <- substr(obj$BinID, 7, nchar(obj$BinID))
rownames(obj@meta.data) <- obj$BinID

obj@meta.data$bin_id <- obj@meta.data[,1] 
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_side_of_Margin_area"] <- "Tumor"
obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Paratumor_side_of_Margin_area"] <- "Paratumor"
obj@meta.data$TLS <- obj@meta.data$TLS_raw
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Conforming"] <- "conforming"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Deviating"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "NotMature"] <- "deviating"
obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Rare_Naive_GCB"] <- "STRH05B_NA"

obj@meta.data[!obj@meta.data$TLS_maturity %in%  c("Mature", "deviating","conforming"), ]$TLS <- "STRH05B_NA"

all_mab <- unique(tmp_bcr_cloneid$clone_subid)

tmp_bcr_cloneid <- tmp_bcr_cloneid[tmp_bcr_cloneid$type == "smrt" & tmp_bcr_cloneid$sample == cur_sample ,]
tmp_bcr_cloneid <- table(tmp_bcr_cloneid$Bin50, tmp_bcr_cloneid$clone_subid) %>% as.data.frame()
colnames(tmp_bcr_cloneid) <- c("BinID", "mab_name", "clone")
tmp_bcr_cloneid <- tmp_bcr_cloneid[tmp_bcr_cloneid$clone != 0,]
tmp_bcr_cloneid <- dplyr::distinct(tmp_bcr_cloneid, BinID, .keep_all = T)
# tmp_bcr_cloneid$BinID <- paste0(cur_sample, "_", tmp_bcr_cloneid$BinID)

obj@meta.data <- dplyr::left_join(obj@meta.data, tmp_bcr_cloneid, by = "BinID") %>%
  column_to_rownames(var="BinID")
obj@meta.data$bin_id <- rownames(obj@meta.data)
obj@meta.data$mab_name <- as.character(obj@meta.data$mab_name)
obj@meta.data$mab_name <- ifelse(is.na(obj@meta.data$mab_name), "other", obj@meta.data$mab_name)
obj@meta.data$mab_name <- ifelse(obj@meta.data$TLS != "STRH05B_NA" & obj@meta.data$mab_name == "other", "TLS", obj@meta.data$mab_name)
obj@meta.data$mab_name <- ifelse(obj@meta.data$mab_name == "other", obj@meta.data$Bin_Region, obj@meta.data$mab_name)
obj@meta.data$mab_name <- ifelse(obj@meta.data$mab_name %in% c("Tumor","Invasive_zone"), "Tumor", obj@meta.data$mab_name)

# obj@meta.data$mab_name[obj@meta.data$mab_name %in% c("Ig_10", "Ig_14", "Ig_23", "Ig_30", "Ig_43", "Ig_5", "Ig_55", "Ig_67")] <- "anti-IFI30_antibody"
# obj@meta.data$mab_name[obj@meta.data$mab_name %in% c("Ig_38","Ig_39","Ig_49","Ig_58","Ig_88")] <- "anti-TP53_antibody"

# mab_color <- c("#E13B43", "#369945", "#1979A3", "#663D8F",  "#E0507D", "#EA7B25", "#40BCE4", "#EFA7B9", "#DB896C")
# mab_color <-  c("#33a02c")
# mab_color <-  c("#EE3743")
mab_color <-  c("#EE3743", "#33a02c", "#0B7FAB", "#6a3d9a", "#ff7f00",  "#F6A9BD", "#fdbf6f", "#c51b7d",  
                "#6F6C9E", "#01665e", "#a6cee3", "#b2df8a", "#ffff99", "#bf812d",  "#999999", "#BB7DB2", "#824615", "#ffff33", "#86DBD4",
                "#BFE2E3","#A1CFFA","#78BDAD","#D45651","#397A7F","#F0918E","#EEE8DA","#1F5392","#A0BFAF",
                "#AE98D6","#ECCBDC","#54BAD3","#8b4a4b","#DB896C","#AABAC2","#ffae3b",
                '#CCCCCC','#B5B5B5','#A092C3','#DDA0DD','#03A4C6','#7AC5CD','#A3D9ED','#00E5EE','#F9EBDA','#F5DEB3',
                '#98689E','#E84115','#FFB5C5','#00FF7F','#F9D01C','#B03060','#00ABDC','#D2691E','#03A464','#FF7F00',
                '#8968CD','#1C5B75')
# SpatialDimPython(obj = obj,
#                  prefix = glue("{outdir}"),
#                  plot_item = "mab_name",
#                  plot_order = c("mAb41", "mAb42", "mAb32,mAb40", "mAb43", "mAb34",  "mAb35",  
#                                 "mAb36",  "mAb39", "mAb33",      
#                                  "TLS","Tumor","Paratumor"),
# tmp_color = c(mab_color,"#9CCCE8", "#eaeaea","#EAD2DF"),
#                  name = glue("merge_bin20"))
obj$BinID = obj$bin_id
table(obj$mab_name)
unique(obj$mab_name)

## 外扩一圈
# df2 <- obj@meta.data[obj@meta.data$mab_name == "anti-TP53_antibody",]
for (i in c("Ig_88", "Ig_58", "Ig_39", "Ig_49", "Ig_38")){
  # for (i in c("Ig_5", "Ig_10", "Ig_14", "Ig_23", "Ig_30", "Ig_43", "Ig_55", "Ig_67")){
  
  df2 <- obj@meta.data[obj@meta.data$mab_name %in% i,]
  
  df_location <- df2 %>%
    dplyr::mutate(dplyr::across(c(col, row), as.numeric)) %>%
    dplyr::distinct(row, col)
  
  df_expand <- df_location %>%
    st_expand_n(
      n = 1, 
      x = "col", 
      y = "row", 
      expand_n = "expand_n", 
      mode = 2
    )
  df_expand <- df_expand %>%
    dplyr::select(col, row, TLS_center = group.by, distance = expand_n)
  
  
  obj@meta.data <- obj@meta.data %>%
    dplyr::mutate(dplyr::across(c(col, row), as.numeric)) %>%
    dplyr::select(-dplyr::any_of(c("TLS_center", "distance"))) %>%
    dplyr::left_join(df_expand, by = c("row", "col"))
  
  obj@meta.data$mab_name[obj@meta.data$distance %in% c("0", "1")] <- i
}


rownames(obj@meta.data) <- obj$bin_id

## 2837B
rownames(obj@meta.data) <- obj$bin_id
SpatialDimPython(obj = obj,
                 prefix = glue("{outdir}"),
                 plot_item = "mab_name",
                 plot_order = c( "Ig_39", "Ig_49", "Ig_58",
                                 "TLS", "Tumor_tissue"),

                 tmp_color = c(mab_color[2:4],"#9CCCE8", "#eaeaea","#EAD2DF"),
                 name = glue("{cur_sample}_cloneid_bin50"))


#### Figure 3C ====
bcr_smrt_merge <- readRDS("~/45sample_bcr0626.rds")
changeo_merge <- readRDS("~/changeo/merge4.rds")



use_id <- fread("~/use_id.txt")

for(cur_tls in unique(use_id$tls)){
  cur_sample <- substr(cur_tls, 3, 7)
  outdir <- glue("~/tree_plot/{cur_tls}")
  dir.create(outdir, recursive = T)
  ab_result <- readRDS("~/ab_merge_clonefam_v2.rds")
  ab_result$mab_name <- ""
  ab_result$mab_name[grepl("Ig", ab_result$clone_subid)] <- ab_result$clone_subid[grepl("Ig", ab_result$clone_subid)]
  # ab_result$mab_name[grepl("ifi30", ab_result$clone_subid)] <- "ifi30"
  ab_result$mab_name[grepl(cur_sample, ab_result$clone_subid)] <- "other_ig"
  ab_result <- ab_result[ab_result$mab_name != ""]
  all_clone <- ab_result$clone_id[!grepl("other_ig", ab_result$mab_name)] %>% unique()
  for (x in all_clone) {
    tmp_bcr_cloneid <- subset(ab_result, ab_result$clone_id == x)
    mab_name <- tmp_bcr_cloneid$clone_subid[!grepl("other_ig", unique(tmp_bcr_cloneid$mab_name))] %>% unique()
    mab_name <- paste(mab_name, collapse = ",")
    ab_result$mab_name[ab_result$clone_id == x] <- mab_name
  }
  
  for(ID in use_id$id[use_id$tls == cur_tls]){
    test_graph_df <- fread(glue('~/tree_plot/{cur_tls}/{ID}_graph_df.txt'))
    test_graph_df$sequence_id <- gsub("seqID_", "", test_graph_df$to)
    
    cur_all_seq_df <- test_graph_df %>% tidyr::separate_rows(sequence_merge, sep = ",")
    cur_all_seq_df$sequence_merge <- gsub("seqID_", "", cur_all_seq_df$sequence_merge )
    cur_all_seq_df <- cur_all_seq_df[!grepl("Germline", cur_all_seq_df$sequence_merge),]
    cur_all_seq_df <- cur_all_seq_df[!grepl("Inferred", cur_all_seq_df$sequence_merge),]
    cur_all_seq_df <- as.data.frame(cur_all_seq_df)
    cur_all_seq_df$label <- str_replace_all(cur_all_seq_df$label, "\\|", "_")
    cur_changeo <- changeo_merge[changeo_merge$sequence_id %in% cur_all_seq_df$sequence_merge, ]
    cur_changeo <- cur_changeo  %>% 
      dplyr::left_join(bcr_smrt_merge[, c("sequence_id",
                                          "SpotLight_Anno", "CellSubType", "TLS_maturity",
                                          "TLS_raw", "Bin50", "Bin20")]) %>% 
      dplyr::left_join(cur_all_seq_df[, c("sequence_merge",
                                          "steps", "distance",
                                          "order", "weight", "label")], by = c("sequence_id" = "sequence_merge"))
    cur_changeo_smrt <- cur_changeo[cur_changeo$type == "Spatial",]
    
    obj <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
    obj@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                               sep = "\t",
                               header = T,
                               row.names = 1) %>% as.data.frame()
    obj@meta.data$bin_id <- rownames(obj@meta.data)
    obj@meta.data$TLS_maturity[is.na(obj@meta.data$TLS_maturity)] <- glue("ST{cur_sample}_NA")
    obj@meta.data$TLS_raw[is.na(obj@meta.data$TLS_raw)] <- glue("ST{cur_sample}_NA")
    obj@meta.data$Bin_Region[is.na(obj@meta.data$Bin_Region)] <- "None"
    obj@meta.data$CellSubType[is.na(obj@meta.data$CellSubType)] <- "None"
    obj@meta.data$CellSubType[is.na(obj@meta.data$CellSubType)] <- "None"
    
    # obj@meta.data$bin_id <- obj@meta.data[,1] 
    obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"
    obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Tumor_side_of_Margin_area"] <- "Tumor"
    obj@meta.data$Bin_Region[obj@meta.data$Bin_Region == "Paratumor_side_of_Margin_area"] <- "Paratumor"
    obj@meta.data$TLS <- obj@meta.data$TLS_raw
    obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Conforming"] <- "conforming"
    obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Deviating"] <- "deviating"
    obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "NotMature"] <- "deviating"
    obj@meta.data$TLS_maturity[obj@meta.data$TLS_maturity == "Rare_Naive_GCB"] <- glue("ST{cur_sample}_NA")
    obj@meta.data$TLS_maturity[is.na(obj@meta.data$TLS_maturity)] <- glue("ST{cur_sample}_NA")
    
    obj@meta.data[!obj@meta.data$TLS_maturity %in%  c("Mature", "deviating","conforming"), ]$TLS_raw <- glue("ST{cur_sample}_NA")
    
    
    obj@meta.data$TLS_raw[obj@meta.data$TLS_maturity == glue("ST{cur_sample}_NA")] <- glue("ST{cur_sample}_NA")
    obj@meta.data$TLS_maturity[obj@meta.data$TLS_raw == glue("ST{cur_sample}_NA")] <- glue("ST{cur_sample}_NA")
    cur_meta <- obj@meta.data
    
    
    
    cur_mab <- cur_all_seq_df$c_call[!grepl("G1", cur_all_seq_df$c_call)]%>% unique()
    
    tmp_bcr_cloneid <- cur_changeo
    tmp_bcr_cloneid <- tmp_bcr_cloneid[tmp_bcr_cloneid$type == "Spatial",] %>% dplyr::arrange(desc(order))
    
    tmp_bcr_cloneid_tab <- tmp_bcr_cloneid$Bin50 %>% table() %>% as.data.frame()
    colnames(tmp_bcr_cloneid_tab) <- c("bin_id","clone")
    tmp_bcr_cloneid_tab <- tmp_bcr_cloneid_tab %>%
      dplyr::left_join(tmp_bcr_cloneid[, c("Bin50", "order", "weight", "label")], by = c("bin_id" = "Bin50"))
    node_count_tab <- tmp_bcr_cloneid_tab %>%
      dplyr::distinct(bin_id, order, .keep_all = T) %>%
      as.data.frame() %>%
      dplyr::group_by(bin_id, order) %>%
      dplyr::summarise(count = n()) %>%
      as.data.frame()
    if(any(node_count_tab$count > 1)){
      print(cur_mab)
    }
    tmp_bcr_cloneid_tab <- dplyr::distinct(tmp_bcr_cloneid_tab, bin_id, .keep_all = T)
    seq_bin <- tmp_bcr_cloneid_tab$bin_id
    
    cur_obj <- obj
    
    cur_obj@meta.data <- dplyr::left_join(cur_obj@meta.data, tmp_bcr_cloneid_tab, by = "bin_id")
    cur_obj@meta.data$clone <- ifelse(is.na(cur_obj@meta.data$clone), "other", "other_ig")
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$TLS != glue("ST{cur_sample}_NA") & cur_obj@meta.data$clone == "other", "TLS", cur_obj@meta.data$clone)
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone == "other", cur_obj@meta.data$Bin_Region, cur_obj@meta.data$clone)
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone %in% c("Tumor", "Invasive_zone"), "Tumor", cur_obj@meta.data$clone)
    # b_naive_bin <- cur_obj@meta.data$bin_id[cur_obj@meta.data$TLS_raw == cur_tls& cur_obj@meta.data$celltype == "NaiveB"]
    gcb_bin <- cur_obj@meta.data$bin_id[cur_obj@meta.data$TLS_raw == cur_tls& cur_obj@meta.data$CellSubType == "GCB"]
    # cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone == "TLS" & cur_obj@meta.data$bin_id %in% b_naive_bin, "NaiveB", cur_obj@meta.data$clone)
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone == "TLS" & cur_obj@meta.data$bin_id %in% gcb_bin, "GCB", cur_obj@meta.data$clone)
    
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone == "other_ig" & cur_obj@meta.data$bin_id %in% seq_bin, cur_obj@meta.data$label, cur_obj@meta.data$clone)
    
    
    all_label <- unique(cur_obj@meta.data$label)
    all_label <- all_label[!is.na(all_label)] %>% sort()
    cur_mab_seq <- cur_changeo_smrt$sequence_id[!grepl(cur_sample, cur_changeo_smrt$label)]
    cur_mab_clone <- bcr_smrt_merge$clone_id[bcr_smrt_merge$sequence_id %in% cur_mab_seq]
    
    cur_mab_clone_bin <- bcr_smrt_merge[bcr_smrt_merge$clone_id %in% cur_mab_clone & bcr_smrt_merge$sample == cur_sample,]$Bin50
    cur_obj@meta.data$clone <- ifelse((!cur_obj@meta.data$clone %in%  all_label) & cur_obj@meta.data$bin_id %in% cur_mab_clone_bin,
                                      cur_mab, cur_obj@meta.data$clone)
    
    rownames(cur_obj@meta.data) <- cur_obj$bin_id
    
    label_col <- rep("#A2C4CF", length(all_label))
    label_col[which(grepl("Ig", all_label))] <- "#3A8339"
    
    SpatialDimPython(obj = cur_obj,
                     prefix = glue("{outdir}/"),
                     plot_item = "clone",
                     plot_order = c(all_label, cur_mab,"GCB", "TLS","Tumor","Paratumor", "None"),
                     tmp_color = c(label_col, "#3A8339", "black","#00ABDC", "#eaeaea","#EAD2DF", "white"),
                     name = glue("{ID}_spatial_all"))
    
    
    
    cur_obj <- obj
    
    cur_obj@meta.data <- dplyr::left_join(cur_obj@meta.data, tmp_bcr_cloneid_tab, by = "bin_id")
    cur_obj@meta.data$clone <- ifelse(is.na(cur_obj@meta.data$clone), "other", "other_ig")
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$TLS != glue("ST{cur_sample}_NA") & cur_obj@meta.data$clone == "other", "TLS", cur_obj@meta.data$clone)
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone == "other", cur_obj@meta.data$Bin_Region, cur_obj@meta.data$clone)
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone %in% c("Tumor", "Invasive_zone"), "Tumor", cur_obj@meta.data$clone)
    # b_naive_bin <- cur_obj@meta.data$bin_id[cur_obj@meta.data$TLS_raw == cur_tls& cur_obj@meta.data$celltype == "NaiveB"]
    gcb_bin <- cur_obj@meta.data$bin_id[cur_obj@meta.data$TLS_raw == cur_tls& cur_obj@meta.data$CellSubType == "GCB"]
    # cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone == "TLS" & cur_obj@meta.data$bin_id %in% b_naive_bin, "NaiveB", cur_obj@meta.data$clone)
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone == "TLS" & cur_obj@meta.data$bin_id %in% gcb_bin, "GCB", cur_obj@meta.data$clone)
    
    cur_obj@meta.data$clone <- ifelse(cur_obj@meta.data$clone == "other_ig" & cur_obj@meta.data$bin_id %in% seq_bin, cur_obj@meta.data$order, cur_obj@meta.data$clone)
    
    rownames(cur_obj@meta.data) <- cur_obj$bin_id
    all_order <- unique(cur_obj@meta.data$order)
    all_order <- all_order[!is.na(all_order)] %>% sort()
    
    SpatialDimPython(obj = subset(cur_obj, TLS_raw == cur_tls & clone %in% c(all_order,"GCB",  "TLS")),
                     prefix = glue("{outdir}/"),
                     plot_item = "clone",
                     plot_order = c(all_order, "GCB",  "TLS"),
                     tmp_color = c(mycolor_merge[1:length(all_order)],"black", "#eaeaea"),
                     name = glue("{ID}_tls"))
    
    
    SpatialDimPython(obj = cur_obj,
                     prefix = glue("{outdir}/"),
                     plot_item = "clone",
                     plot_order = c(all_order, "GCB", "TLS","Tumor","Paratumor", "None"),
                     tmp_color = c(mycolor_merge[1:length(all_order)],"black", "#00ABDC", "#eaeaea","#EAD2DF", "white"),
                     name = glue("{ID}_spatial"))
    
    
    tls_obj <- subset(cur_obj, TLS_raw == cur_tls)
    tls_x_min <- tls_obj$row %>% min() - 10
    tls_y_min <- tls_obj$col %>% min() - 10
    tls_x_max <- tls_obj$row %>% max() + 10
    tls_y_max <- tls_obj$col %>% max() + 10
    tls_expand <- subset(cur_obj, row %in% c(tls_x_min:tls_x_max) & col %in% c(tls_y_min:tls_y_max))
    
    all_mab_bin <- ab_result$Bin50[ab_result$mab_name == cur_mab]
    all_mab_bin <- all_mab_bin[!is.na(all_mab_bin)]
    
    tls_expand@meta.data$clone <- ifelse(tls_expand@meta.data$clone == 'TLS'& tls_expand@meta.data$TLS_raw != cur_tls, tls_expand@meta.data$Bin_Region, tls_expand@meta.data$clone)
    
    tls_expand@meta.data$clone <- ifelse(tls_expand@meta.data$clone %in% c("Tumor", "Invasive_zone", 'Paratumor', 'TLS') & tls_expand@meta.data$bin_id %in% all_mab_bin, cur_mab, tls_expand@meta.data$clone)
    
    tls_expand@meta.data$clone <- ifelse(tls_expand@meta.data$clone %in% all_order, tls_expand@meta.data$label, tls_expand@meta.data$clone)
    
    all_label <- tls_expand@meta.data$label %>% unique()
    all_label <- all_label[!is.na(all_label)]
    
    label_col <- rep("#A2C4CF", length(all_label))
    label_col[which(grepl("Ig", all_label))] <- "#3A8339"
    # c_plot_order <- c(tls_expand@meta.data$label %>% unique(), cur_mab, "GCB",  "TLS", "Tumor")
    
    
    SpatialDimPython(obj = tls_expand,
                     prefix = glue("{outdir}/"),
                     plot_item = "clone",
                     plot_order = c(all_label, cur_mab, "GCB",  "TLS", "Tumor", "Paratumor", "Invasive_zone"),
                     tmp_color = c(label_col, "#3A8339", "black", "#eaeaea", "white", "white", "white"),
                     name = glue("{ID}_tls_expand"))
    
    
  }
  
}






### check tree cell type====
source("~/expand_n_layer.R")
merge_bcr <- readRDS("~/merged_bcr.rds")
merge_bcr_rh05b <- merge_bcr[merge_bcr$sample == "RH05B",]
cur_sample <- "RH05B"
hcc_st <- readRDS(glue('~/{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds'))
hcc_st@meta.data <-read.table(glue("~/{cur_sample}.txt"),
                              sep = "\t",
                              header = T,
                              row.names = 1) %>% as.data.frame()
rownames(hcc_st@meta.data) <- str_replace(rownames(hcc_st@meta.data), "\\w{5}_", "")
hcc_st@meta.data$SpotLight_Anno[is.na(hcc_st@meta.data$SpotLight_Anno)] <- "None"
hcc_st@meta.data$CellSubType[is.na(hcc_st@meta.data$CellSubType)] <- "None"
hcc_st@meta.data$SpotLight_Anno[hcc_st@meta.data$SpotLight_Anno == "pDC"] <- "Myeloid"
hcc_st$BinID <- rownames(hcc_st@meta.data)

meta_data <- hcc_st@meta.data
meta_data$bin50 <- paste(meta_data$SampleID, meta_data$BinID, sep = "_")
meta_data$tree_label <- ""

## tree STRH05B_10_1072
tree_df <- fread("~/tree_plot/STRH05B_10/STRH05B_10_1072_graph_df.txt")
tree_df <- tidyr::separate_rows(tree_df, sequence_merge, sep = ",")
cur_seq <- tree_df$sequence_merge[tree_df$order == 3]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()

meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()

## tree STRH05B_21_1275
tree_df <- fread("~/tree_plot/STRH05B_21/STRH05B_21_1275_graph_df.txt")
tree_df <- tidyr::separate_rows(tree_df, sequence_merge, sep = ",")
cur_seq <- tree_df$sequence_merge[tree_df$order == 3]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()


cur_seq <- tree_df$sequence_merge[tree_df$order == 4]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()



## tree STRH05B_12_1173 ====
tree_df <- fread("~/tree_plot/STRH05B_12/STRH05B_12_1173_graph_df.txt")
tree_df <- tidyr::separate_rows(tree_df, sequence_merge, sep = ",")
cur_seq <- tree_df$sequence_merge[tree_df$order == 3]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()

meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()


cur_seq <- tree_df$sequence_merge[tree_df$order == 4]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()




## tree STRH05B_9_1072 ====
tree_df <- fread("~/tree_plot/STRH05B_9/STRH05B_9_1072_graph_df.txt")
tree_df <- tidyr::separate_rows(tree_df, sequence_merge, sep = ",")
cur_seq <- tree_df$sequence_merge[tree_df$order == 6]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()

meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()


cur_seq <- tree_df$sequence_merge[tree_df$order == 8]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()


## tree STRH05B_4_1173 ====
tree_df <- fread("~/tree_plot/STRH05B_4/STRH05B_4_1173_graph_df.txt")
tree_df <- tidyr::separate_rows(tree_df, sequence_merge, sep = ",")
cur_seq <- tree_df$sequence_merge[tree_df$order == 3]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()

meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()


cur_seq <- tree_df$sequence_merge[tree_df$order == 4]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()



cur_seq <- tree_df$sequence_merge[tree_df$order == 5]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()


## tree STRH05B_7_902 ====
tree_df <- fread("~/tree_plot/STRH05B_7/STRH05B_7_902_graph_df.txt")
tree_df <- tidyr::separate_rows(tree_df, sequence_merge, sep = ",")
cur_seq <- tree_df$sequence_merge[tree_df$order == 4]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()

meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()


cur_seq <- tree_df$sequence_merge[tree_df$order == 5]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()



cur_seq <- tree_df$sequence_merge[tree_df$order == 5]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()



## tree STRH05B_7_1173 ====
tree_df <- fread("~/tree_plot/STRH05B_7/STRH05B_7_1173_graph_df.txt")
tree_df <- tidyr::separate_rows(tree_df, sequence_merge, sep = ",")
cur_seq <- tree_df$sequence_merge[tree_df$order == 3]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()

meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()


cur_seq <- tree_df$sequence_merge[tree_df$order == 4]
cur_bin50 <- merge_bcr_rh05b$bin_id[merge_bcr_rh05b$sequence_id %in% cur_seq]
cur_cell_list <- meta_data$SpotLight_Anno[meta_data$bin50 %in% cur_bin50]
cur_cell_list %>% unique()
meta_data$tree_label <- ""
meta_data$tree_label[meta_data$bin50 %in% cur_bin50] <- "tree"
meta_data_e1 <-   meta_data %>%
  dplyr::filter(tree_label == "tree") %>% 
  st_expand_n(
    n = 1, 
    x = "col", 
    y = "row", 
    expand_n = "expand_n", 
    group.by = "tree_label", 
    mode = 2
  )
# cur_cell_sub_e <- meta_data$CellSubType[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
# cur_cell_sub_e %>% unique()
cur_cell_list_e <- meta_data$SpotLight_Anno[meta_data$row %in% meta_data_e1$row & meta_data$col %in% meta_data_e1$col]
cur_cell_list_e %>% unique()

