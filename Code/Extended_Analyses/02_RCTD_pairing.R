PROJECT_ROOT <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, 'data')
source('Code/Extended_Analyses/utils_io.R')
source('Code/Extended_Analyses/utils_plot.R')
source('Code/Extended_Analyses/02_rctd_functions.R')

#### 参数 ====

library(Seurat)

library(spacexr)

library(Matrix)

library(dplyr)

library(tidyr)

library(data.table)

library(glue)

library(ggplot2)

library(ggpubr)

library(scales)

library(stringr)

library(SoupX)

library(DoubletFinder)

REFERENCE_LABEL_MODE <- "original"

RCTD_CELL_MIN_INSTANCE <- 1L

RCTD_DOUBLET_THRESHOLD <- 20

RCTD_MAX_CORES <- 30L

REFERENCE_TARGET_PER_TYPE <- 100L

REFERENCE_MAX_PER_TYPE <- 2000L

SAMPLE_FILTER <- character()

SKIP_EXISTING_RCTD <- TRUE

SC_DATA_ROOT <- file.path(dirs$raw, "scRNAseq", "singlemRNA")

SC_SAMPLE_66 <- file.path(dirs$raw, "scRNAseq", "mrna_tcr_66sample.rds")

SC_UMAP_RDS <- file.path(dirs$raw, "scRNAseq", "66_sc_umap_merged.rds")

SC_RAW_RDS_DIR <- file.path(dirs$processed, "02_RCTD_pairing", "sample_QC")

DATA_SC_REF <- file.path(dirs$processed, "02_RCTD_pairing", "66_sc_celltype_merged.rds")

DATA_SC_MAP <- file.path(dirs$metadata, "Pacbio_SingleCell.csv")

DATA_SC_LIST <- file.path(dirs$metadata, "sc_34.txt")

DATA_THRE <- file.path(dirs$metadata, "threshold.txt")

ST_ROOT <- dirs$st_root

ST_META_DIR <- dirs$st_meta

RCTD_OUT_ROOT <- file.path(dirs$stats_logs, "RCTD_Plots_new")

SIPEAC_ROOT <- dirs$sipeac

FIG_DIR <- file.path(dirs$figures, "R1_Q2b_RCTD_singlet反卷积")

TBL_DIR <- dirs$tables

dir.create(dirname(DATA_SC_REF), recursive = TRUE, showWarnings = FALSE)

dir.create(RCTD_OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

#### 单细胞质控与参考注释 ====

sample_list_66 <- readRDS(SC_SAMPLE_66)

dir.create(SC_RAW_RDS_DIR, recursive = TRUE, showWarnings = FALSE)

for (cur_s in sample_list_66) {
    out_rds <- file.path(SC_RAW_RDS_DIR, glue("{cur_s}.rds"))
    if (file.exists(out_rds)) 
        next
    df <- Read10X(glue("{SC_DATA_ROOT}/{cur_s}/filtered_feature_bc_matrix"))
    df <- CreateSeuratObject(counts = df, project = glue("ST{cur_s}PBMC"), min.cells = 3, min.features = 200)
    df <- run_preliminary_process(df)
    multiplet_rate <- if (ncol(df) < 10000) 
        0.075
    else 0.1
    df <- run_DoubletFinder(df, multiplet_rate)
    df <- RenameCells(df, add.cell.id = cur_s)
    saveRDS(df, out_rds)
}

rds_list <- lapply(sample_list_66, function(s) {
    obj <- readRDS(file.path(SC_RAW_RDS_DIR, glue("{s}.rds")))
    obj@meta.data$sampleid <- s
    RenameCells(obj, add.cell.id = s)
})

sample_66_merged <- merge(rds_list[[1]], rds_list[-1])

all_merged1 <- readRDS(SC_UMAP_RDS)

all_merged1@meta.data <- all_merged1@meta.data %>% dplyr::mutate(celltype = dplyr::case_when(seurat_clusters %in% 
    c(3, 15) ~ "B", seurat_clusters %in% c(20, 29) ~ "Plasma", seurat_clusters %in% c("0", "1", "2", 
    "4", "6", "8", "13", "14", "18", "19", "24", "27", "11") ~ "T", seurat_clusters %in% c(5, 9, 16) ~ 
    "Myeloid", seurat_clusters %in% c(6, 26) ~ "Endothelial", seurat_clusters %in% c(17) ~ "Fibroblast", 
    seurat_clusters %in% c(7, 10, 12, 21, 22, 23, 25, 28) ~ "Tumor&Hepatocyte"))

print(table(all_merged1@meta.data$celltype, useNA = "ifany"))

genes_to_check <- c("GPC3", "EPCAM", "KRT8", "KRT18", "KRT19", "APOA1", "APOA2", "ALB", "TTR", "APOB", 
    "PTPRC", "CD3D", "CD3E", "CD4", "CD8A", "KLRB1", "NCR1", "CD19", "CD79A", "MS4A1", "IGHG1", "MZB1", 
    "SDC1", "CD68", "CD163", "CD14", "C1QA", "C1QB", "S100A9", "S100A8", "MMP19", "KIT", "TPSAB1", "TPSB2", 
    "CLEC9A", "CD1C", "ITGAX", "LILRA4", "UGCG", "CLEC4C", "COL1A1", "COL3A1", "TAGLN", "RGS5", "MYH11", 
    "MCAM", "ACTA2", "PDGFRB", "PECAM1", "VWF", "ACKR1", "PLVAP", "MKI67", "TOP2A")

p_dot <- DotPlot(all_merged1, features = genes_to_check, group.by = "celltype") + RotatedAxis() + ggplot2:::coord_flip() + 
    scale_color_gradient2(low = "#005fbf", mid = "white", high = "#bf0000")

save_plot(p_dot, file.path(FIG_DIR, "celltype_dotplot"), width = 10, height = 10, if_png = TRUE, if_pdf = TRUE)

saveRDS(all_merged1, DATA_SC_REF)

#### RCTD doublet mode ====

all_merged1 <- readRDS(DATA_SC_REF)

all_merged1 <- apply_reference_label_mode(all_merged1, REFERENCE_LABEL_MODE)

sc_list <- fread(DATA_SC_LIST)

sc_list$sample_name <- gsub("ST", "", sc_list$sample_name)

if (length(SAMPLE_FILTER)) sc_list <- sc_list[sc_list$sample_name %in% SAMPLE_FILTER, , drop = FALSE]

sc_sample_list_new <- list(RH17B = c("RH17_T"), RH16B = c("RH16"), RH15B = "RH15", RH14B = c("RH14", 
    "RH14_A_P", "RH14_A_T"), RH10T = "RH10", RH10P = c("RH10_A_P"), RH10B = "RH10", RH08B = "RH08", RH06P = "RH06P", 
    RH06B = "RH06", RH05B = "RH05", RH03T = "RH03", `3163B` = "3163_T", `3154T` = "3154_T", `3154P` = "3154_P", 
    `3154B` = "3154_T", `3125T` = "3125_T", `3099T` = "3099T", `3099P` = "3099P", `3099B` = "3099T", 
    `3001B` = "3001T", `2976T` = "2976T", `2976B` = "2976T", `2934T` = "2934T", `2934B` = "2934T", `2931T` = "2931T", 
    `2931P` = "2931P", `2931B` = "2931B", `2907T` = "2907T", `2907B` = "2907T", `2896P` = "2896P", `2893T` = c("2893T", 
        "2893P"), `2893P` = "2893P", `2893B` = c("2893T", "2893P"), `2879T` = "2879T", `2879B` = "2879T", 
    `2846B` = "2846T", `2837B` = "2837T", `2740P` = "2740P", `2903B` = "2903P", `2903P` = "2903P", `2865T` = "2865T", 
    `2823T` = "2823T")

reference_composition <- list()

for (cur_sample in sc_list$sample_name) {
    sampleid <- cur_sample
    cur_sc_sample_list <- sc_sample_list_new[[cur_sample]]
    existing_result <- file.path(RCTD_OUT_ROOT, glue("RCTD_{sampleid}"), "result.rds")
    if (SKIP_EXISTING_RCTD && file.exists(existing_result)) {
        if (REFERENCE_LABEL_MODE == "refined") {
            reference_composition[[sampleid]] <- summarize_refined_reference(all_merged1, matched_sample_ids = cur_sc_sample_list, 
                target_per_type = REFERENCE_TARGET_PER_TYPE, max_per_type = REFERENCE_MAX_PER_TYPE) %>% 
                mutate(cur_sample = sampleid, .before = 1)
        }
        next
    }
    if (REFERENCE_LABEL_MODE == "refined") {
        refined_reference <- build_refined_reference(all_merged1, matched_sample_ids = cur_sc_sample_list, 
            target_per_type = REFERENCE_TARGET_PER_TYPE, max_per_type = REFERENCE_MAX_PER_TYPE, seed_key = sampleid)
        sc_ref <- refined_reference$object
        reference_composition[[sampleid]] <- refined_reference$composition %>% mutate(cur_sample = sampleid, 
            .before = 1)
    }
    else {
        sc_ref <- subset(all_merged1, sampleid %in% cur_sc_sample_list)
    }
    counts <- sc_ref@assays$RNA@counts
    cell_types <- as.factor(sc_ref@meta.data$celltype)
    names(cell_types) <- rownames(sc_ref@meta.data)
    nUMI <- colSums(counts)
    names(nUMI) <- rownames(sc_ref@meta.data)
    reference <- Reference(counts, cell_types, nUMI)
    st_obj <- readRDS(file.path(ST_ROOT, glue("{sampleid}_scRefST2934T_random100_marker20_spotlight_Spatial.rds")))
    st_obj@meta.data <- read.table(file.path(ST_META_DIR, glue("{sampleid}.txt")), sep = "\t", header = TRUE, 
        row.names = 1) %>% as.data.frame()
    st_obj@meta.data$bin_id <- rownames(st_obj@meta.data)
    rownames(st_obj@meta.data) <- str_remove_all(rownames(st_obj@meta.data), glue("{sampleid}_"))
    puck <- SpatialRNA(st_obj@images$slice1@coordinates[, c("row", "col")], st_obj@assays$Spatial@counts, 
        colSums(st_obj@assays$Spatial@counts))
    myRCTD <- create.RCTD(puck, reference, max_cores = RCTD_MAX_CORES, CELL_MIN_INSTANCE = RCTD_CELL_MIN_INSTANCE, 
        DOUBLET_THRESHOLD = RCTD_DOUBLET_THRESHOLD)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
    results <- myRCTD@results
    resultsdir <- file.path(RCTD_OUT_ROOT, glue("RCTD_{sampleid}"))
    dir.create(resultsdir, recursive = TRUE)
    norm_weights <- normalize_weights(results$weights)
    cell_type_names <- myRCTD@cell_type_info$info[[2]]
    spatialRNA <- myRCTD@spatialRNA
    plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
    plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
    plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df)
    plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
    plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
    saveRDS(results, file.path(resultsdir, "result.rds"))
    bcr_pie_df <- results$results_df %>% group_by(spot_class) %>% summarise(value = n(), .groups = "drop") %>% 
        mutate(ratio = value/sum(value), percent = percent(ratio, accuracy = 0.1), label = paste0(spot_class, 
            " (", percent, ")"))
    bcr_pie <- ggpie(bcr_pie_df, "value", label = bcr_pie_df$label, lab.pos = "in", lab.font = "white", 
        fill = "spot_class", color = "white", palette = mycolor_new)
    save_plot(bcr_pie, file.path(resultsdir, glue("{sampleid}_spot_class")), width = 4, height = 4, if_png = TRUE, 
        if_pdf = TRUE)
}

if (length(reference_composition) > 0) {
    composition_suffix <- if (length(SAMPLE_FILTER) == 0) {
        "all"
    }
    else {
        paste0("shard_", SAMPLE_FILTER[[1]], "_n", length(SAMPLE_FILTER))
    }
    write.csv(bind_rows(reference_composition), file.path(TBL_DIR, paste0("R1_Q2b_refined_reference_composition_", 
        composition_suffix, ".csv")), row.names = FALSE)
}

#### Singlet bin中的SIPEAC配对 ====

sc_list <- fread(DATA_SC_LIST)

sc_list$sample_name <- gsub("ST", "", sc_list$sample_name)

if (length(SAMPLE_FILTER)) sc_list <- sc_list[sc_list$sample_name %in% SAMPLE_FILTER, , drop = FALSE]

thre_list <- read.table(DATA_THRE, header = FALSE)$V1

binsize_list <- c(50)

pred_all <- list()

dir.create(file.path(RCTD_OUT_ROOT, "bin50"), recursive = TRUE, showWarnings = FALSE)

for (cur_sample in sc_list$sample_name) {
    cur_rctd_path <- file.path(RCTD_OUT_ROOT, glue("RCTD_{cur_sample}"), "result.rds")
    if (!file.exists(cur_rctd_path)) 
        next
    cur_rctd <- readRDS(cur_rctd_path)
    cur_singlet <- cur_rctd$results_df[cur_rctd$results_df$spot_class == "singlet", ]
    if (nrow(cur_singlet) == 0) 
        next
    singlet_spot <- rownames(cur_singlet)
    yang_score <- fread(file.path(SIPEAC_ROOT, "2025_fold0_epoch100", glue("{cur_sample}_cdr3aa.csv")))
    colnames(yang_score) <- c("IGH_cloneIds", "IGL_cloneIds", "score")
    pair_name_y <- paste(yang_score$IGH_cloneIds, yang_score$IGL_cloneIds, sep = "_")
    yang_score$IGH_cloneIds <- paste0("clone_", as.character(yang_score$IGH_cloneIds))
    yang_score$IGL_cloneIds <- paste0("clone_", as.character(yang_score$IGL_cloneIds))
    yang_score_mat <- tidyr::pivot_wider(yang_score, names_from = "IGL_cloneIds", values_from = "score") %>% 
        as.data.frame()
    yang_score$pair_name_y <- pair_name_y
    yang_score <- yang_score %>% group_by(IGH_cloneIds) %>% mutate(rank_s = n() - rank(score)) %>% as.data.table()
    all_igh <- yang_score_mat$IGH_cloneIds
    all_igl <- colnames(yang_score_mat[, -1]) %>% as.character()
    yang_score_mat <- as.matrix(yang_score_mat[, -1])
    rownames(yang_score_mat) <- all_igh
    colnames(yang_score_mat) <- all_igl
    for (bin_size in binsize_list) {
        input_path <- file.path(SIPEAC_ROOT, "input_allelic", glue("bin{bin_size}"), glue("{cur_sample}_info.txt"))
        if (!file.exists(input_path)) 
            next
        input <- read.table(input_path, sep = "\t", header = TRUE)
        input <- input[input$binid %in% singlet_spot, ]
        input$clone_id <- as.character(input$clone_id)
        input$clone_name <- paste(input$locus, input$clone_id, sep = "_")
        input_spatial <- input[grepl("Spatial", input$type), ]
        clone_num <- table(input_spatial$clone_name)
        df_clone_num <- input_spatial[, c("clone_name", "clone_num")] %>% distinct(clone_name, .keep_all = TRUE)
        input_spatial$clone_num <- as.numeric(clone_num[input_spatial$clone_name])
        df_bin_num <- input_spatial %>% distinct(binid, clone_name, .keep_all = TRUE) %>% count(clone_name, 
            name = "bin_num") %>% mutate(bin_num = as.numeric(bin_num))
        df_clone_num <- df_clone_num %>% left_join(df_bin_num, by = "clone_name")
        input_spatial <- input_spatial %>% left_join(df_bin_num, by = "clone_name")
        gt_path <- file.path(SIPEAC_ROOT, "grouth_truth", glue("bin{bin_size}"), glue("{cur_sample}.tsv"))
        if (!file.exists(gt_path)) 
            next
        test_mat <- read.table(gt_path, sep = "\t", header = TRUE, row.names = 1)
        test_mat <- test_mat[test_mat$Hclone %in% input_spatial$clone_name, ]
        if (nrow(test_mat) == 0) 
            next
        L_clone_t <- sapply(unique(test_mat$Hclone), function(h) paste(test_mat$Lclone[test_mat$Hclone == 
            h], collapse = ","))
        test_mat$Lclone_t <- L_clone_t[test_mat$Hclone]
        test_mat$pair_name <- paste(test_mat$Hclone, test_mat$Lclone, sep = "_")
        input_mat <- table(input_spatial$binid, input_spatial$clone_name)
        exclued_h <- c()
        exclued_l <- c()
        output_pair_list <- list()
        n_loop <- 1
        repeat {
            if (sum(!colnames(input_mat) %in% c(exclued_h, exclued_l)) == 1) 
                break
            cur_input <- input_mat[, !colnames(input_mat) %in% c(exclued_h, exclued_l), drop = FALSE]
            x_1 <- .pick_cols(cur_input, "IGH")
            x_2 <- .pick_cols(cur_input, "IGL|IGK")
            combined <- (rowSums(x_1) >= 1) & (rowSums(x_2) >= 1)
            if (sum(combined) == 0) 
                break
            cur_input_HL <- if (sum(combined) == 1) {
                out <- data.frame(matrix(cur_input[combined, ], nrow = 1))
                rownames(out) <- rownames(cur_input)[combined]
                colnames(out) <- colnames(cur_input)
                out
            }
            else cur_input[combined, , drop = FALSE]
            cur_IgH_mat <- .pick_cols(cur_input_HL, "IGH")
            cur_IgL_mat <- .pick_cols(cur_input_HL, "IGL|IGK")
            if (n_loop == 1) {
                oH <- c()
                oL <- c()
                ob <- c()
                ohc <- c()
                olc <- c()
                ohp <- c()
                olp <- c()
                for (cb in rownames(cur_IgH_mat)) {
                  hi <- which(cur_IgH_mat[cb, ] >= 1)
                  h1 <- cur_IgH_mat[cb, hi]
                  h1 <- h1/sum(h1)
                  h1 <- h1[h1 == max(h1)]
                  li <- which(cur_IgL_mat[cb, ] >= 1)
                  l1 <- cur_IgL_mat[cb, li]
                  l1 <- l1/sum(l1)
                  l1 <- l1[l1 == max(l1)]
                  for (h in names(h1)) for (l in names(l1)) {
                    oH <- c(oH, h)
                    oL <- c(oL, l)
                    ob <- c(ob, cb)
                    ohc <- c(ohc, cur_IgH_mat[cb, h])
                    olc <- c(olc, cur_IgL_mat[cb, l])
                    ohp <- c(ohp, h1[h])
                    olp <- c(olp, l1[l])
                  }
                }
                cop <- data.frame(Hclone = oH, Lclone = oL, hcount = ohc, lcount = olc, hp = ohp, lp = olp, 
                  binid = ob)
                if (nrow(cop) == 0) 
                  break
                cop$pair_name <- paste(cop$Hclone, cop$Lclone, sep = "_")
                cop <- cop %>% left_join(test_mat[, c("Hclone", "Lclone_t")])
                write.table(cop, file.path(RCTD_OUT_ROOT, glue("bin{bin_size}"), glue("{cur_sample}_pair_raw.tsv")), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
                exist_idx <- which(test_mat$pair_name %in% cop$pair_name)
                tm_e <- test_mat[exist_idx, ]
                sm_h <- cop %>% group_by(Hclone) %>% summarise(n_pair = n())
                tm1 <- left_join(tm_e, distinct(cop, pair_name, pair_count = 1L), by = "pair_name") %>% 
                  left_join(sm_h, by = "Hclone")
                tm1$IGH_cloneIds <- gsub("IGH_|IGL_|IGK_", "", tm1$Hclone)
                tm1$IGL_cloneIds <- gsub("IGH_|IGL_|IGK_", "", tm1$Lclone)
                tm1$pair_name_y <- paste(tm1$IGH_cloneIds, tm1$IGL_cloneIds, sep = "_")
                tm1 <- left_join(tm1, yang_score[, c("pair_name_y", "score", "rank_s")], by = "pair_name_y")
                write.table(tm1, file.path(RCTD_OUT_ROOT, glue("bin{bin_size}"), glue("{cur_sample}_test_summary.tsv")), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
            }
            oH <- c()
            oL <- c()
            ob <- c()
            pcl <- c()
            sl <- c()
            for (cb in rownames(cur_IgH_mat)) {
                mh_c <- max(cur_IgH_mat[cb, ])
                mh <- colnames(cur_IgH_mat)[which(cur_IgH_mat[cb, ] == mh_c)]
                ml_c <- max(cur_IgL_mat[cb, ])
                ml <- colnames(cur_IgL_mat)[which(cur_IgL_mat[cb, ] == ml_c)]
                for (h in mh) for (l in ml) {
                  oH <- c(oH, h)
                  oL <- c(oL, l)
                  ob <- c(ob, cb)
                  pcl <- c(pcl, mh_c * ml_c)
                  hn <- gsub("IGH_|IGK_|IGL_", "", h)
                  ln <- gsub("IGH_|IGK_|IGL_", "", l)
                  sl <- c(sl, yang_score_mat[glue("clone_{hn}"), glue("clone_{ln}")] %>% as.numeric())
                }
            }
            cop <- data.frame(Hclone = oH, Lclone = oL, binid = ob, pair_count_expand = pcl, score = sl)
            if (nrow(cop) == 0) 
                break
            cop$pair_name <- paste(cop$Hclone, cop$Lclone, sep = "_")
            cop <- cop %>% group_by(pair_name) %>% mutate(pair_count = n(), pair_count_expand = sum(pair_count_expand)) %>% 
                as.data.frame()
            cop$pc_dot_s <- cop$pair_count * cop$score
            cop <- cop %>% left_join(test_mat[, c("Hclone", "Lclone_t")]) %>% left_join(df_clone_num, 
                by = c(Hclone = "clone_name")) %>% arrange(desc(clone_num), desc(pair_count))
            cop2 <- cop %>% group_by(Lclone) %>% filter(pc_dot_s == max(pc_dot_s)) %>% distinct(pair_name, 
                .keep_all = TRUE)
            cop2$clone_num <- as.numeric(cop2$clone_num)
            h_count <- table(cop2$Hclone) %>% sort() %>% rev()
            save_pair <- c()
            for (h in names(h_count)) {
                sp <- cop2[cop2$Hclone == h, ]
                sm <- sp[sp$pc_dot_s == max(sp$pc_dot_s), ]
                if (nrow(sm) == 1) 
                  save_pair <- c(save_pair, sm$pair_name)
            }
            if (length(save_pair) == 0) {
                output_pair_list[[n_loop]] <- cop2
                break
            }
            cop3 <- cop2[cop2$pair_name %in% save_pair, ]
            output_pair_list[[n_loop]] <- cop3
            exclued_h <- c(exclued_h, cop3$Hclone)
            exclued_l <- c(exclued_l, cop3$Lclone)
            n_loop <- n_loop + 1
        }
        if (length(output_pair_list) == 0) 
            next
        opf <- bind_rows(output_pair_list) %>% mutate(pair_count_expand = pair_count_expand * pair_count) %>% 
            group_by(Hclone) %>% filter(pair_count_expand == max(pair_count_expand)) %>% as.data.frame()
        write.table(opf, file.path(RCTD_OUT_ROOT, glue("bin{bin_size}"), glue("{cur_sample}_final_e100.tsv")), 
            sep = "\t", quote = FALSE, row.names = FALSE)
        opf <- opf %>% distinct(pair_name, .keep_all = TRUE)
        opt <- opf[opf$Hclone %in% test_mat$Hclone, ]
        write.table(opt, file.path(RCTD_OUT_ROOT, glue("bin{bin_size}"), glue("{cur_sample}_test_e100.tsv")), 
            sep = "\t", quote = FALSE, row.names = FALSE)
        for (cur_thre in thre_list) {
            opf1 <- opf[opf$bin_num >= cur_thre, ]
            opt1 <- opf1[opf1$Hclone %in% test_mat$Hclone, ]
            dbn1 <- df_bin_num[df_bin_num$bin_num >= cur_thre, ]
            n_H <- nrow(dbn1[grepl("IGH", dbn1$clone_name), ])
            n_H_pred <- length(unique(opf1$Hclone))
            if (nrow(opt1) == 0) 
                next
            pred_all[[glue("{cur_sample}_bin{bin_size}_{cur_thre}")]] <- data.frame(cur_sample = cur_sample, 
                n_H_clone = n_H, n_L_clone = nrow(dbn1[!grepl("IGH", dbn1$clone_name), ]), n_H_clone_pred = n_H_pred, 
                eff = n_H_pred/n_H, n_test = nrow(test_mat[test_mat$Hclone %in% opf1$Hclone, ]), threshold = cur_thre, 
                bin_size = bin_size, acc = sum(opt1$pair_name %in% test_mat$pair_name)/nrow(opt1))
        }
    }
}

if (length(pred_all) == 0) {
    df_pred_all <- data.frame(cur_sample = character(), n_H_clone = numeric(), n_L_clone = numeric(), 
        n_H_clone_pred = numeric(), eff = numeric(), n_test = numeric(), threshold = numeric(), bin_size = numeric(), 
        acc = numeric())
    warning("No SIPEAC-on-singlet accuracy rows were generated for the selected samples.")
} else {
    df_pred_all <- bind_rows(pred_all) %>% arrange(desc(n_test))
}

write.table(df_pred_all, file.path(RCTD_OUT_ROOT, "sipeac_on_rctd_acc_summary.tsv"), sep = "\t", quote = FALSE, 
    row.names = FALSE)

#### 配对准确率与singlet比例 ====

rctd_acc_path <- file.path(RCTD_OUT_ROOT, "sipeac_on_rctd_acc_summary.tsv")

if (!file.exists(rctd_acc_path)) {
    stop("Missing RCTD-specific SIPEAC summary: ", rctd_acc_path, ". Run the SIPEAC pairing section first.")
}

sipeac_rctd_acc <- fread(rctd_acc_path)

sipeac_rctd_acc <- sipeac_rctd_acc[!is.na(sipeac_rctd_acc$acc) & threshold == 1 & bin_size == 50]

sipeac_acc <- fread(file.path(SIPEAC_ROOT, "allelic_yscore_e100", "summary_e100.tsv"))

sipeac_acc <- sipeac_acc[!is.na(sipeac_acc$acc) & threshold == 1 & bin_size == 50]

select_sample <- sipeac_rctd_acc$cur_sample %>% unique()

df_long <- bind_rows(sipeac_acc %>% select(cur_sample, n_H_clone_pred, acc) %>% mutate(condition = "sipeac"), 
    sipeac_rctd_acc %>% select(cur_sample, n_H_clone_pred, acc) %>% mutate(condition = "sipeac_rctd")) %>% 
    filter(cur_sample %in% select_sample)

sc_list <- fread(DATA_SC_LIST)

sc_list$sample_name <- gsub("ST", "", sc_list$sample_name)

if (length(SAMPLE_FILTER)) sc_list <- sc_list[sc_list$sample_name %in% SAMPLE_FILTER, , drop = FALSE]

singlet_df <- list()

for (s in sc_list$sample_name) {
    rp <- file.path(RCTD_OUT_ROOT, glue("RCTD_{s}"), "result.rds")
    if (!file.exists(rp)) 
        next
    rd <- readRDS(rp)$results_df
    singlet_df[[s]] <- data.frame(cur_sample = s, singlet_rate = sum(rd$spot_class == "singlet")/nrow(rd))
}

singlet_df <- bind_rows(singlet_df) %>% mutate(singlet_pct = singlet_rate * 100)

write.csv(singlet_df, file.path(RCTD_OUT_ROOT, "R1_Q2b_singlet_rate.csv"), row.names = FALSE)

write.csv(singlet_df, file.path(TBL_DIR, "R1_Q2b_singlet_rate.csv"), row.names = FALSE)

df_long$singlet_rate <- singlet_df$singlet_rate[match(df_long$cur_sample, singlet_df$cur_sample)]

df_wide <- df_long %>% select(-n_H_clone_pred, -singlet_rate) %>% pivot_wider(names_from = condition, 
    values_from = acc) %>% mutate(diff = sipeac - sipeac_rctd) %>% arrange(desc(diff))

print(df_wide)

tt <- t.test(df_wide$sipeac_rctd, df_wide$sipeac, paired = TRUE)

wt <- wilcox.test(df_wide$sipeac_rctd, df_wide$sipeac, paired = TRUE, exact = FALSE)

test_df <- data.frame(test = c("paired t-test", "paired Wilcoxon"), statistic = c(round(unname(tt$statistic), 
    4), round(unname(wt$statistic), 4)), p_value = c(tt$p.value, wt$p.value), n_pairs = nrow(df_wide), 
    mean_sipeac = round(mean(df_wide$sipeac), 4), mean_sipeac_rctd = round(mean(df_wide$sipeac_rctd), 
        4), mean_diff = round(mean(df_wide$diff), 4))

write.csv(test_df, file.path(TBL_DIR, "R1_Q2b_acc_paired_test.csv"), row.names = FALSE)

write.csv(df_wide, file.path(TBL_DIR, "R1_Q2b_acc_paired_wide.csv"), row.names = FALSE)

print(test_df)

sample_ord <- levels(reorder(df_long$cur_sample, df_long$acc))

df_long_plot <- df_long %>% mutate(cur_sample = factor(cur_sample, levels = sample_ord))

p_label_acc <- sprintf("paired t-test p = %.4f\npaired Wilcoxon p = %.4f", tt$p.value, wt$p.value)

n_wide <- df_long %>% select(cur_sample, condition, n_H_clone_pred) %>% pivot_wider(names_from = condition, 
    values_from = n_H_clone_pred) %>% filter(!is.na(sipeac), !is.na(sipeac_rctd))

n_tt <- t.test(n_wide$sipeac, n_wide$sipeac_rctd, paired = TRUE)

n_wt <- wilcox.test(n_wide$sipeac, n_wide$sipeac_rctd, paired = TRUE, exact = FALSE)

p_label_n <- sprintf("paired t-test p = %.4f\npaired Wilcoxon p = %.4f", n_tt$p.value, n_wt$p.value)

p1a <- ggplot(df_long_plot, aes(x = cur_sample, y = acc, group = condition, color = condition)) + geom_line(linewidth = 0.6) + 
    geom_point(size = 2) + scale_color_manual(values = c(sipeac_rctd = "#005fbf", sipeac = "#bf0000")) + 
    scale_y_continuous(name = "Accuracy", limits = c(0, 1.05)) + labs(x = "Sample", color = "Condition", 
    title = "Pairing accuracy (sipeac vs sipeac_rctd)", subtitle = p_label_acc) + theme_bw(base_size = 12) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", plot.subtitle = element_text(size = 9, 
        color = "grey20"))

save_plot(p1a, file.path(FIG_DIR, "R1_Q2b_acc_line"), width = 6, height = 4, if_png = TRUE, if_pdf = TRUE)

p1b <- ggplot(df_long_plot, aes(x = cur_sample, y = n_H_clone_pred, fill = condition)) + geom_col(position = position_dodge(0.8), 
    width = 0.7) + geom_text(aes(label = n_H_clone_pred), position = position_dodge(0.8), size = 2.5, 
    vjust = -0.3) + scale_fill_manual(values = c(sipeac_rctd = "#005fbf", sipeac = "#bf0000")) + labs(x = "Sample", 
    y = "Predicted H-clones (n)", fill = "Condition", title = "Predicted H-clones (n) per condition", 
    subtitle = p_label_n) + theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top", plot.subtitle = element_text(size = 9, color = "grey20"))

save_plot(p1b, file.path(FIG_DIR, "R1_Q2b_acc_n_bar"), width = 6, height = 4, if_png = TRUE, if_pdf = TRUE)

p2 <- ggplot(singlet_df, aes(x = reorder(cur_sample, singlet_rate), y = singlet_pct)) + geom_col(aes(fill = singlet_pct), 
    width = 0.7) + geom_text(aes(label = sprintf("%.3f%%", singlet_pct)), hjust = -0.1, size = 3) + scale_fill_gradient(low = "#9ecae1", 
    high = "#005fbf", name = "Singlet %") + coord_flip() + expand_limits(y = max(singlet_df$singlet_pct) * 
    1.15) + labs(x = "Sample", y = "Singlet rate (%)", title = "Singlet rate per sample (sorted)") + 
    theme_bw(base_size = 12) + theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"))

save_plot(p2, file.path(FIG_DIR, "R1_Q2b_singlet_rate_by_sample"), width = 5, height = 5, if_png = TRUE, 
    if_pdf = TRUE)

p_label <- sprintf("paired t-test p = %.4f\npaired Wilcoxon p = %.4f", tt$p.value, wt$p.value)

dw <- df_wide %>% filter(!is.na(diff)) %>% mutate(better = ifelse(diff > 0, "sipeac > rctd", ifelse(diff < 
    0, "sipeac < rctd", "tie")), cur_sample = factor(cur_sample, levels = cur_sample[order(diff)]))

p3 <- ggplot(dw, aes(x = cur_sample, y = diff, fill = better)) + geom_hline(yintercept = 0, linetype = "dashed", 
    color = "grey40") + geom_col(width = 0.7) + scale_fill_manual(values = c(`sipeac > rctd` = "#bf0000", 
    `sipeac < rctd` = "#005fbf", tie = "grey70")) + coord_flip() + labs(x = "Sample", y = "Accuracy difference (sipeac − sipeac_rctd)", 
    title = "Per-sample paired accuracy difference", subtitle = p_label) + theme_bw(base_size = 12) + 
    theme(legend.position = "top", plot.subtitle = element_text(size = 9, color = "grey20"))

save_plot(p3, file.path(FIG_DIR, "R1_Q2b_acc_paired_diff"), width = 6, height = 4, if_png = TRUE, if_pdf = TRUE)

singlet_tbl <- singlet_df

ss <- intersect(df_wide$cur_sample, singlet_tbl$cur_sample)

df_cor <- data.frame(sample = ss, singlet_rate = singlet_tbl$singlet_rate[match(ss, singlet_tbl$cur_sample)], 
    accuracy_sipeac = df_wide$sipeac[match(ss, df_wide$cur_sample)], accuracy_sipeac_rctd = df_wide$sipeac_rctd[match(ss, 
        df_wide$cur_sample)], n_clone = sipeac_rctd_acc$n_H_clone_pred[match(ss, sipeac_rctd_acc$cur_sample)]) %>% 
    mutate(singlet_pct = singlet_rate * 100)

out_cor <- rbind(run_cor("accuracy_sipeac", "sipeac"), run_cor("accuracy_sipeac_rctd", "sipeac_rctd"))

write.csv(out_cor, file.path(TBL_DIR, "R1_Q2b_acc_vs_singlet.csv"), row.names = FALSE)

p4 <- plot_scatter("accuracy_sipeac", "sipeac", "#bf0000")

p5 <- plot_scatter("accuracy_sipeac_rctd", "sipeac_rctd", "#005fbf")

save_plot(p4, file.path(FIG_DIR, "R1_Q2b_acc_vs_singlet_sipeac"), width = 6, height = 5, if_png = TRUE, 
    if_pdf = TRUE)

save_plot(p5, file.path(FIG_DIR, "R1_Q2b_acc_vs_singlet_rctd"), width = 6, height = 5, if_png = TRUE, 
    if_pdf = TRUE)
