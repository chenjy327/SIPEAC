PROJECT_ROOT <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, 'data')
source('Code/Extended_Analyses/R/io.R')
source('Code/Extended_Analyses/R/plot.R')

#### 参数 ====

library(data.table)

library(dplyr)

library(glue)

library(stringr)

library(alakazam)

SC_BCR_ROOT <- file.path(dirs$raw, "scBCRseq")

DATA_BCR_CHG <- file.path(dirs$raw, "scRNAseq", "merge4.rds")

SIPEAC_FINAL <- file.path(dirs$sipeac, "final_version")

DATA_AB_125 <- file.path(dirs$metadata, "new_antibody", "final_genescript2.csv")

DATA_BCR_MERGED <- file.path(dirs$raw, "scRNAseq", "merged_bcr.rds")

DATA_SC_MAP <- file.path(dirs$metadata, "Pacbio_SingleCell.csv")

TBL_DIR <- dirs$tables

PROCESSED_DIR <- dirs$processed

dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

#### 预测结果与参考一致的重轻链配对 ====

changeo_output <- readRDS(DATA_BCR_CHG)

sample_list <- changeo_output$sample[changeo_output$type == "Spatial"] %>% unique()

changeo_output$clone_id <- as.character(changeo_output$clone_id)

changeo_output$clone_name <- paste(changeo_output$locus, changeo_output$clone_id, sep = "_")

allelic_filtered_results <- list()

binsize_list <- c(20, 50, 110)

for (current_sample in sample_list) {
    pred_matched_merge <- list()
    for (current_binsize in binsize_list) {
        result_file_path <- file.path(SIPEAC_FINAL, glue("bin{current_binsize}"), glue("{current_sample}_test_e100.tsv"))
        if (!file.exists(result_file_path)) {
            message(glue("{current_sample} file not exist"))
            next
        }
        pred_data <- fread(result_file_path)
        pred_matched <- pred_data[pred_data$Lclone == pred_data$Lclone_t, ]
        if (nrow(pred_matched) == 0) {
            message(glue("No matching predictions for {current_sample} bin{current_binsize}"))
            next
        }
        pred_matched_merge[[current_binsize]] <- pred_matched
    }
    if (length(pred_matched_merge) == 0) {
        message(glue("No predictions for {current_sample}"))
        next
    }
    pred_matched_merge <- dplyr::bind_rows(pred_matched_merge) %>% dplyr::distinct(Hclone, Lclone, .keep_all = TRUE)
    input_filtered <- changeo_output[changeo_output$clone_name %in% pred_matched_merge$Hclone, ]
    input_filtered <- dplyr::left_join(input_filtered, pred_matched_merge[, c("Hclone", "Lclone", "Lclone_t")], 
        by = c(clone_name = "Hclone"))
    input_filtered$clone_name <- paste(current_sample, input_filtered$clone_name, sep = "_")
    allelic_filtered_results[[current_sample]] <- input_filtered
}

final_result <- bind_rows(allelic_filtered_results)

write_rds(final_result, "pred_new_antibody.rds")

#### SHM、TLS与肿瘤区域筛选 ====

pred_sc <- final_result[final_result$type == "SingleCell", ]

ab_125 <- fread(DATA_AB_125)

pred_spatial <- final_result[final_result$type == "Spatial", ]

ab_125$sequence_id <- sub("-", "_", ab_125$raw_sc_id)

ab_125_igg <- ab_125[ab_125$isotype == "IGHG1", ]

pred_sc <- pred_sc %>% dplyr::left_join(ab_125_igg[, c("sequence_id", "IG_name")])

ab_clone_label <- pred_sc[!is.na(pred_sc$IG_name), ] %>% dplyr::distinct(IG_name, .keep_all = TRUE)

bcr_data <- readRDS(DATA_BCR_MERGED)

pred_spatial_new <- pred_spatial[!pred_spatial$clone_name %in% ab_clone_label$clone_name, ]

pred_spatial_new <- dplyr::left_join(pred_spatial_new, bcr_data[, c("sequence_id", "TLS_raw", "Bin_Region")])

pred_spatial_new_clone <- countClones(pred_spatial_new, clone = "clone_name")

pred_spatial_new <- dplyr::left_join(pred_spatial_new, pred_spatial_new_clone[, c("clone_name", "seq_count")])

pred_spatial_shm <- pred_spatial_new[pred_spatial_new$shm > 0.02, ]

pred_spatial_tls <- pred_spatial_shm[!grepl("NA", pred_spatial_shm$TLS_raw), ]

pred_spatial_tumor <- pred_spatial_tls[pred_spatial_tls$region_type == "Tumor", ]

pred_spatial_final <- dplyr::arrange(pred_spatial_tumor, desc(seq_count)) %>% dplyr::distinct(clone_name, 
    Lclone_t, .keep_all = TRUE) %>% dplyr::select(sample, c_call, TLS_raw, shm, region_type, seq_count, 
    clone_name, sequence_id, Lclone_t)

pred_spatial_final$c_call <- substr(pred_spatial_final$c_call, 1, 5)

pred_spatial_final <- dplyr::arrange(pred_spatial_final, desc(seq_count)) %>% dplyr::distinct(clone_name, 
    .keep_all = TRUE)

save_table(pred_spatial_final, file.path(TBL_DIR, "R4_new_antibody_candidates"), if_csv = TRUE)

#### 重轻链序列与合成信息 ====

pacbio_sc_info <- fread(DATA_SC_MAP)

cur_sc_path_list <- pacbio_sc_info$SC_Path %>% lapply(function(x) strsplit(x, ";")[[1]]) %>% unlist()

cur_sc_sample <- basename(cur_sc_path_list)

cur_sc_path_list <- file.path(SC_BCR_ROOT, basename(dirname(cur_sc_path_list)), cur_sc_sample)

sc_path_info <- data.frame(path = cur_sc_path_list, sample = cur_sc_sample)

pred_sc_new <- pred_sc[pred_sc$clone_name %in% pred_spatial_final$clone_name, ]

final_genescript <- list()

for (i in 1:nrow(pred_spatial_final)) {
    cur_clone <- pred_spatial_final$clone_name[i]
    cur_df <- pred_sc_new[pred_sc_new$clone_name == cur_clone, ]
    if (nrow(cur_df) == 0) 
        next
    seq <- cur_df$sequence[1]
    cur_lclone <- cur_df$Lclone_t[1]
    cur_lbarcode <- changeo_output[changeo_output$clone_name == cur_lclone & changeo_output$type == "SingleCell", 
        ]
    cur_sample <- pred_spatial_final$sample[i]
    cur_tls <- pred_spatial_final$TLS_raw[i]
    sc_i_s <- cur_df$sample[1]
    sc_i_seq <- cur_df$sequence_id[1]
    sc_i_path <- sc_path_info$path[sc_path_info$sample == sc_i_s] %>% unique()
    contig_csv <- glue("{sc_i_path}/filtered_contig_annotations.csv")
    if (!file.exists(contig_csv)) {
        message(glue("[{i}] 10x contig 不存在: {contig_csv}，跳过"))
        next
    }
    raw_10x <- fread(contig_csv)
    raw_10x$sequence_id <- paste(sc_i_s, raw_10x$contig_id, sep = "_")
    raw_10x_sel <- raw_10x[raw_10x$sequence_id == sc_i_seq, ]
    raw_10x_cell <- raw_10x[raw_10x$barcode == raw_10x_sel$barcode, ]
    raw_10x_h <- raw_10x_sel
    raw_10x_h <- raw_10x_h[raw_10x_h$length == max(raw_10x_h$length), ]
    seq_nt <- paste(raw_10x_h$fwr1_nt, raw_10x_h$cdr1_nt, raw_10x_h$fwr2_nt, raw_10x_h$cdr2_nt, raw_10x_h$fwr3_nt, 
        raw_10x_h$cdr3_nt, raw_10x_h$fwr4_nt, sep = "")
    seq_aa <- paste(raw_10x_h$fwr1, raw_10x_h$cdr1, raw_10x_h$fwr2, raw_10x_h$cdr2, raw_10x_h$fwr3, raw_10x_h$cdr3, 
        raw_10x_h$fwr4, sep = "")
    light_info <- raw_10x_cell[raw_10x_cell$chain != "IGH", ]
    light_info <- light_info[light_info$sequence_id %in% cur_lbarcode$sequence_id, ]
    if (nrow(light_info) == 0) 
        next
    if (nrow(light_info) > 1) {
        message(glue("[{i}] {sc_i_seq} multi light chains found"))
    }
    seq_nt_l <- paste(light_info$fwr1_nt, light_info$cdr1_nt, light_info$fwr2_nt, light_info$cdr2_nt, 
        light_info$fwr3_nt, light_info$cdr3_nt, light_info$fwr4_nt, sep = "")
    seq_aa_l <- paste(light_info$fwr1, light_info$cdr1, light_info$fwr2, light_info$cdr2, light_info$fwr3, 
        light_info$cdr3, light_info$fwr4, sep = "")
    type1 <- str_replace(raw_10x_h$c_gene, "GH", "g")
    type2 <- light_info$chain
    type <- paste(type1, type2, sep = ";")
    cur_result <- data.frame(type = type, IG_H_nt = seq_nt, IG_H_aa = seq_aa, IG_L_nt = seq_nt_l, IG_L_aa = seq_aa_l, 
        name = glue("newig_{i}"), tls = cur_tls, clone_name = cur_clone, cell_barcode = light_info$barcode, 
        sample_st = cur_sample, sample_sc = sc_i_s, seq_id = seq, stringsAsFactors = FALSE)
    final_genescript[[glue("{i}")]] <- cur_result
}

final_genescript <- dplyr::bind_rows(final_genescript)

save_table(final_genescript, file.path(TBL_DIR, "R4_new_antibody_synthesis"), if_csv = TRUE)
