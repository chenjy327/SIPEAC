apply_reference_label_mode <- function(object, mode) {
    if (mode == "original") 
        return(object)
    if (mode != "refined") 
        stop("Unknown reference label mode: ", mode)
    cluster <- as.character(object@meta.data$seurat_clusters)
    object@meta.data$celltype <- dplyr::case_when(cluster %in% c("3", "15") ~ "B", cluster %in% c("20", 
        "29") ~ "Plasma", cluster %in% c("0", "1", "2", "4", "8", "11", "13", "14", "18", "19", "24", 
        "27") ~ "T", cluster %in% c("5", "9", "16") ~ "Myeloid", cluster %in% c("6", "26") ~ "Endothelial", 
        cluster == "17" ~ "Fibroblast", cluster %in% c("7", "12", "22", "25") ~ "Tumor", cluster %in% 
            c("10", "21", "23", "28") ~ "Hepatocyte", TRUE ~ NA_character_)
    if (anyNA(object@meta.data$celltype)) {
        stop("Refined reference mapping produced NA cell types.")
    }
    object
}

build_refined_reference <- function(object, matched_sample_ids, target_per_type, max_per_type, seed_key) {
    metadata <- object@meta.data
    all_types <- sort(unique(metadata$celltype))
    selected <- character()
    composition <- list()
    set.seed(sum(utf8ToInt(seed_key)))
    for (cell_type in all_types) {
        all_cells <- rownames(metadata)[metadata$celltype == cell_type]
        matched_cells <- rownames(metadata)[metadata$celltype == cell_type & metadata$sampleid %in% matched_sample_ids]
        supplemented_cells <- character()
        if (length(matched_cells) < target_per_type) {
            supplement_pool <- setdiff(all_cells, matched_cells)
            supplement_n <- min(target_per_type - length(matched_cells), length(supplement_pool))
            if (supplement_n > 0) {
                supplemented_cells <- sample(supplement_pool, supplement_n)
            }
        }
        type_cells <- unique(c(matched_cells, supplemented_cells))
        if (length(type_cells) > max_per_type) {
            type_cells <- sample(type_cells, max_per_type)
        }
        selected <- c(selected, type_cells)
        composition[[cell_type]] <- data.frame(celltype = cell_type, n_matched_available = length(matched_cells), 
            n_supplemented = sum(type_cells %in% supplemented_cells), n_reference = length(type_cells))
    }
    list(object = subset(object, cells = selected), composition = bind_rows(composition))
}

summarize_refined_reference <- function(object, matched_sample_ids, target_per_type, max_per_type) {
    metadata <- object@meta.data
    all_types <- sort(unique(metadata$celltype))
    bind_rows(lapply(all_types, function(cell_type) {
        all_cells <- rownames(metadata)[metadata$celltype == cell_type]
        matched_cells <- rownames(metadata)[metadata$celltype == cell_type & metadata$sampleid %in% matched_sample_ids]
        supplement_n <- min(max(target_per_type - length(matched_cells), 0), length(setdiff(all_cells, 
            matched_cells)))
        data.frame(celltype = cell_type, n_matched_available = length(matched_cells), n_supplemented = supplement_n, 
            n_reference = min(length(matched_cells) + supplement_n, max_per_type))
    }))
}

run_preliminary_process <- function(sc_data) {
    sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
    sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 500)
    sc_data <- ScaleData(sc_data, features = rownames(sc_data))
    sc_data <- RunPCA(sc_data, features = VariableFeatures(sc_data))
    sc_data <- FindNeighbors(sc_data, dims = 1:20)
    sc_data <- FindClusters(sc_data, resolution = 0.1)
    sc_data <- RunUMAP(sc_data, dims = 1:20)
    sc_data
}

run_SoupX <- function(dataDir.soupx, contaminationRange, project.name) {
    Sc <- load10X(dataDir.soupx)
    Sc <- autoEstCont(Sc, contaminationRange = contaminationRange)
    Sc_adj <- adjustCounts(Sc, roundToInt = TRUE)
    CreateSeuratObject(Sc_adj, project = project.name, min.cells = 3, min.features = 200)
}

run_DoubletFinder <- function(seurat_obj, multiplet_rate) {
    sweep.res.list <- DoubletFinder::paramSweep_v3(seurat_obj, PCs = 1:20, sct = FALSE)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)
    nExp_poi <- round(multiplet_rate * nrow(seurat_obj@meta.data))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, 
        reuse.pANN = FALSE, sct = FALSE)
    colnames(seurat_obj@meta.data)[ncol(seurat_obj@meta.data)] <- "DoubletFinder"
    seurat_obj <- subset(seurat_obj, subset = DoubletFinder == "Singlet")
    seurat_obj@meta.data <- seurat_obj@meta.data[, c(1, 2, 3)]
    seurat_obj
}

run_cor <- function(av, label) {
    ct <- cor.test(df_cor$singlet_rate, df_cor[[av]], method = "spearman", exact = FALSE)
    fit <- lm(reformulate("singlet_rate", response = av), data = df_cor, weights = n_clone)
    s <- summary(fit)
    data.frame(accuracy_of = label, n = nrow(df_cor), spearman_rho = round(unname(ct$estimate), 4), spearman_p = signif(ct$p.value, 
        3), lm_slope = round(coef(fit)["singlet_rate"], 4), lm_slope_p = signif(coef(s)[, "Pr(>|t|)"]["singlet_rate"], 
        3), lm_r2 = round(s$r.squared, 4), n_clone_sum = sum(df_cor$n_clone))
}

plot_scatter <- function(av, label, col) {
    ct <- cor.test(df_cor$singlet_rate, df_cor[[av]], method = "spearman", exact = FALSE)
    p <- ggplot(df_cor, aes(x = singlet_pct, y = .data[[av]])) + geom_smooth(method = "lm", se = TRUE, 
        color = col, fill = col, alpha = 0.15) + geom_point(aes(size = n_clone), color = col, alpha = 0.8) + 
        scale_size_continuous(range = c(2, 8), name = "n_clone (predicted H)") + labs(x = "Singlet rate (%)", 
        y = glue("{label} accuracy"), title = glue("{label} accuracy vs singlet rate (n = {nrow(df_cor)})"), 
        subtitle = sprintf("Spearman rho = %.2f, p = %.3f", ct$estimate, ct$p.value)) + theme_bw(base_size = 12) + 
        theme(plot.subtitle = element_text(size = 9, color = "grey20"), legend.position = "right") + 
        coord_cartesian(ylim = c(-0.05, 1.05))
    if (requireNamespace("ggrepel", quietly = TRUE)) 
        p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 2.8, max.overlaps = 20)
    else p <- p + geom_text(aes(label = sample), vjust = -1.2, size = 2.8)
    p
}

.pick_cols <- function(df, pat) {
    hit <- grepl(pat, colnames(df))
    if (sum(hit) == 1) {
        out <- data.frame(V1 = df[, hit])
        colnames(out) <- colnames(df)[hit]
        out
    }
    else df[, hit]
}
