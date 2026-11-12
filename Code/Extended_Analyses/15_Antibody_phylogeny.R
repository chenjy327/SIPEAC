#### 新增抗体家族与TLS系统发育树 ====
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
library(Seurat)
library(alakazam)
library(igraph)
library(readxl)

bcr_smrt_merge <- readRDS(file.path(input_dir, "data/merged_bcr.rds"))

changeo_merge <- readRDS(file.path(input_dir, "phylogeny/changeo_merge_07.rds"))

changeo_merge$clone_family_07 <- paste(changeo_merge$patient, changeo_merge$clone_family_07, sep = "_")

ab_info <- readxl::read_xlsx(file.path(input_dir, "phylogeny/genes_57.xlsx")) %>% as.data.frame()

ab_ifi30_seq <- ab_info$raw_seq

ab_ifi30_clone <- changeo_merge[changeo_merge$sequence %in% ab_ifi30_seq, ]

dir.create(dirname(file.path(processed_dir, "phylogeny/new_57ab_ifi30_clone_batch2.rds")), recursive = TRUE, 
    showWarnings = FALSE)

saveRDS(ab_ifi30_clone, file.path(processed_dir, "phylogeny/new_57ab_ifi30_clone_batch2.rds"))

ab_ifi30_cloneid <- ab_ifi30_clone$clone_id %>% unique()

ab_ifi30_clone <- dplyr::left_join(ab_ifi30_clone, ab_info[, c("raw_seq", "IG_name")], by = c(sequence = "raw_seq"))

ab_merge_clonefam_id <- changeo_merge[changeo_merge$clone_id %in% ab_ifi30_cloneid, ]$clone_family_07 %>% 
    unique()

ab_merge_clonefam <- changeo_merge[changeo_merge$clone_family_07 %in% ab_merge_clonefam_id, ]

ab_merge_clonefam$isotype <- substr(ab_merge_clonefam$c_call, 1, 5)

ab_merge_clonefam$isotype[is.na(ab_merge_clonefam$isotype)] <- "IGHG1"

ab_ifi30_clone$clone_subid <- ab_ifi30_clone$IG_name

ab_merge_clone_subid <- ab_ifi30_clone %>% group_by(clone_id) %>% summarise(clone_subid = paste(sort(unique(IG_name)), 
    collapse = ","))

ab_merge_clonefam <- dplyr::left_join(ab_merge_clonefam, ab_merge_clone_subid)

ab_merge_clonefam$clone_subid[is.na(ab_merge_clonefam$clone_subid)] <- ab_merge_clonefam$clone_id[is.na(ab_merge_clonefam$clone_subid)]

ab_merge_clonefam$isotype_group <- substr(ab_merge_clonefam$isotype, 1, 4)

ab_merge_clonefam <- dplyr::left_join(ab_merge_clonefam, bcr_smrt_merge[, c("sequence_id", "cell_annotation", 
    "TLS_maturity", "TLS_raw", "cell_subtype")])

ab_merge_clonefam$TLS_maturity[is.na(ab_merge_clonefam$TLS_maturity)] <- ab_merge_clonefam$type[is.na(ab_merge_clonefam$TLS_maturity)]

ab_merge_clonefam$cell_annotation[is.na(ab_merge_clonefam$cell_annotation)] <- ab_merge_clonefam$type[is.na(ab_merge_clonefam$cell_annotation)]

ab_merge_clonefam$TLS_raw[is.na(ab_merge_clonefam$TLS_raw)] <- ab_merge_clonefam$type[is.na(ab_merge_clonefam$TLS_raw)]

ab_merge_clonefam <- dplyr::left_join(ab_merge_clonefam, bcr_smrt_merge[, c("sequence_id", "Bin50", 
    "Bin20")])

dir.create(dirname(file.path(processed_dir, "phylogeny/new_57ab_merge_clonefam_batch2.rds")), recursive = TRUE, 
    showWarnings = FALSE)

saveRDS(ab_merge_clonefam, file.path(processed_dir, "phylogeny/new_57ab_merge_clonefam_batch2.rds"))

createTreeLabel <- function(graph, changeo_clone, clone) {
    graph_df <- graph %>% igraph::as_data_frame()
    graph_df <- rbind(graph_df, c("Germline", "Germline", 0, ""))
    dfs_ordered_vertex <- graph %>% igraph::dfs(root = "Germline", order = T)
    dfs_ordered_vertex_table <- cbind(dfs_ordered_vertex$order$label, c(1:length(dfs_ordered_vertex$order)))
    colnames(dfs_ordered_vertex_table) <- c("to", "order")
    graph_df <- merge(graph_df, dfs_ordered_vertex_table, by = "to")
    graph_df_isotype <- changeo_clone@data %>% as_tibble() %>% dplyr::select(sequence_id, c_call)
    graph_df_isotype$c_call <- gsub("IGH", "", graph_df_isotype$c_call)
    graph_df <- merge(graph_df, graph_df_isotype, by.x = "to", by.y = "sequence_id", all.x = T)
    graph_df[is.na(graph_df$c_call), "c_call"] <- substr(graph_df[is.na(graph_df$c_call), "to"], 
        1, 2)
    clone_selected <- clone %>% ungroup() %>% dplyr::select(sequence_id, type) %>% dplyr::mutate(label_ab = dplyr::if_else(type != 
        "other_ig", "other_ig", type))
    graph_df <- merge(graph_df, clone_selected, by.x = "to", by.y = "sequence_id", all.x = T, all.y = F)
    graph_df$type[is.na(graph_df$type)] <- "NA"
    graph_df <- graph_df %>% mutate(label = paste(order, c_call, sep = "|"))
    graph_df$label <- gsub("\\*00", "", graph_df$label)
    graph_df$label <- gsub(",", ".", graph_df$label)
    return(graph_df)
}

changeo_merge <- readRDS(file.path(input_dir, "data/merge4.rds"))

ab_merge_clonefam_batch2 <- readRDS(file.path(processed_dir, "phylogeny/new_57ab_merge_clonefam_batch2.rds"))

merge_metadata <- read.table(file.path(input_dir, "data/Bin50_classic_celltype_metadata_merge.txt"), 
    header = TRUE)

phylip_exec <- Sys.which("dnapars")

ab_merge_clonefam <- ab_merge_clonefam_batch2 %>% filter(type == "smrt") %>% group_by(sequence_id) %>% 
    arrange(desc(grepl("^Ig", clone_subid)), .by_group = TRUE) %>% slice(1) %>% ungroup()

family_list <- c("RH05B_287", "RH05B_902", "RH05B_747", "RH05B_2318", "RH05B_1173", "RH05B_2183", 
    "RH05B_773", "RH05B_2082", "RH05B_1059", "RH05B_1893", "RH05B_276", "RH05B_76", "RH05B_1058", 
    "RH05B_2500", "RH05B_2771")

#### TLS内外家族树 ====
cur_sample <- "RH05B"

obj_ST_raw <- readRDS(glue(file.path(input_dir, "spatial_objects/{cur_sample}.rds")))

df <- merge_metadata %>% filter(SampleID == cur_sample)

obj_ST_base <- obj_ST_raw

obj_ST_base@meta.data <- obj_ST_base@meta.data %>% dplyr::select(BinID, bin) %>% left_join(df, by = "BinID") %>% 
    mutate(Location = case_when(Bin_Region %in% c("Paratumor_tissue", "Paratumor_side_of_Margin_area") ~ 
        "Paratumor", Bin_Region %in% c("Tumor_capsule") ~ "Invasive_zone", Bin_Region %in% c("Tumor_side_of_Margin_area", 
        "Tumor_tissue") ~ "Tumor", TRUE ~ "Other"))

obj_ST_base@meta.data$TLS_Final <- ifelse(obj_ST_base@meta.data$TLS_maturity %in% c("Rare_Naive_GCB", 
    "Filtered"), glue("ST{cur_sample}_NA"), obj_ST_base@meta.data$TLS_Final)

obj_ST_base@meta.data <- obj_ST_base@meta.data %>% distinct(bin, .keep_all = TRUE)

rownames(obj_ST_base@meta.data) <- obj_ST_base@meta.data$bin

tls_table <- table(obj_ST_base@meta.data$TLS_Final)

tls_list <- names(tls_table)[names(tls_table) != glue("ST{cur_sample}_NA")]

rh05b_data <- ab_merge_clonefam %>% filter(sample == cur_sample)

b_cell_types <- c("MemoryB", "NaiveB", "other_B", "Plasma", "GCB")

all_meta <- obj_ST_base@meta.data %>% select(bin, row, col, CellSubType)

obj_ST_base@meta.data$celltype_1quan <- NA

for (i in 1:nrow(all_meta)) {
    cur <- all_meta[i, ]
    bin_id <- cur$bin
    ct <- cur$CellSubType
    if (ct %in% b_cell_types) {
        obj_ST_base@meta.data[bin_id, "celltype_1quan"] <- ct
    }
    else {
        neighbor_coords <- data.frame(row = c(cur$row - 1, cur$row + 1, cur$row, cur$row), col = c(cur$col, 
            cur$col, cur$col - 1, cur$col + 1))
        neighbor_bins <- all_meta %>% semi_join(neighbor_coords, by = c("row", "col")) %>% pull(bin)
        if (length(neighbor_bins) > 0) {
            neighbor_cts <- all_meta %>% filter(bin %in% neighbor_bins) %>% pull(CellSubType)
            if (any(neighbor_cts == "Plasma")) {
                obj_ST_base@meta.data[bin_id, "celltype_1quan"] <- "Plasma"
            }
            else if (any(neighbor_cts %in% b_cell_types)) {
                obj_ST_base@meta.data[bin_id, "celltype_1quan"] <- "Bcell"
            }
            else {
                obj_ST_base@meta.data[bin_id, "celltype_1quan"] <- ct
            }
        }
        else {
            obj_ST_base@meta.data[bin_id, "celltype_1quan"] <- ct
        }
    }
}

bin_info_list <- list()

for (fam in family_list) {
    cat("\n===== Processing family:", fam, "=====\n")
    fam_data <- rh05b_data %>% filter(clone_family_07 == fam)
    if (nrow(fam_data) == 0) {
        cat("  No sequences for this family.\n")
        next
    }
    fam_bins <- unique(fam_data$Bin50)
    if (length(fam_bins) == 0) {
        cat("  No Bin50 info.\n")
        next
    }
    bin_tls_map <- obj_ST_base@meta.data %>% filter(bin %in% fam_bins) %>% select(bin, TLS_Final) %>% 
        distinct()
    if (nrow(bin_tls_map) == 0) {
        cat("  No spots found in ST object.\n")
        next
    }
    for (tls in tls_list) {
        tls_bins <- bin_tls_map %>% filter(TLS_Final == tls) %>% pull(bin)
        if (length(tls_bins) == 0) 
            next
        tls_coords <- obj_ST_base@meta.data %>% filter(bin %in% tls_bins) %>% select(row, col) %>% 
            distinct()
        all_spots <- obj_ST_base@meta.data %>% select(bin, row, col, CellSubType, TLS_Final, Location, 
            Bin_Region, celltype_1quan)
        all_spots$dist_to_tls <- apply(all_spots[, c("row", "col")], 1, function(pt) {
            min(abs(pt[1] - tls_coords$row) + abs(pt[2] - tls_coords$col))
        })
        extended_spots <- all_spots %>% filter(dist_to_tls <= 10)
        if (nrow(extended_spots) == 0) 
            next
        seq_in_extended <- fam_data %>% filter(Bin50 %in% extended_spots$bin)
        if (nrow(seq_in_extended) == 0) 
            next
        clone_count <- length(unique(seq_in_extended$clone_id))
        if (clone_count < 2) 
            next
        seq_in_tls <- seq_in_extended
        cat("  ✔ TLS:", tls, "扩展区域内共有", clone_count, "个克隆，", nrow(seq_in_extended), 
            "条序列。\n")
        seq_with_dist <- seq_in_extended %>% left_join(extended_spots %>% select(bin, dist_to_tls, 
            CellSubType, TLS_Final, Location, Bin_Region, celltype_1quan), by = c(Bin50 = "bin"))
        for (i in 1:nrow(seq_with_dist)) {
            row_data <- seq_with_dist[i, ]
            bin_info <- data.frame(clone_family_07 = fam, clone_id = row_data$clone_id, distance_to_nearest_TLS = row_data$dist_to_tls, 
                bin_id = row_data$Bin50, CellSubType = row_data$CellSubType, celltype_1quan = row_data$celltype_1quan, 
                TLS_Final = row_data$TLS_Final, Location = row_data$Location, Bin_Region = row_data$Bin_Region, 
                stringsAsFactors = FALSE)
            bin_info_list[[length(bin_info_list) + 1]] <- bin_info
        }
        obj_ST <- obj_ST_base
        df_location <- obj_ST@meta.data %>% filter(bin %in% tls_bins) %>% left_join(seq_in_tls %>% 
            select(Bin50, clone_id), by = c(bin = "Bin50")) %>% mutate(across(c(col, row), as.numeric)) %>% 
            distinct(row, col, clone_id)
        if (nrow(df_location) == 0) {
            cat("    No spots found for this family in this TLS\n")
            next
        }
        df_expand <- df_location %>% mutate(distance = 0, clone_label = clone_id) %>% dplyr::select(col, 
            row, clone_label, distance)
        obj_ST@meta.data <- obj_ST@meta.data %>% mutate(across(c(col, row), as.numeric)) %>% select(-any_of(c("clone_label", 
            "distance"))) %>% left_join(df_expand, by = c("row", "col"))
        obj_ST@meta.data <- obj_ST@meta.data %>% distinct(bin, .keep_all = TRUE)
        rownames(obj_ST@meta.data) <- obj_ST@meta.data$bin
        obj_ST@meta.data <- obj_ST@meta.data %>% mutate(label = case_when(!is.na(clone_label) ~ 
            as.character(clone_label), TLS_Final == glue("ST{cur_sample}_NA") ~ Location, TRUE ~ 
            "TLS"))
        all_labels <- unique(obj_ST@meta.data$label)
        background_levels <- c("TLS", "Tumor", "Invasive_zone", "Paratumor")
        clone_levels <- all_labels[!all_labels %in% background_levels]
        plot_order <- c(background_levels, clone_levels)
        bg_colors <- c("#9CCCE8", "#E9ECEB", "#FFE8BF", "#F9F7EE")
        if (length(clone_levels) > 0) {
            clone_colors <- rainbow(length(clone_levels))
            names(clone_colors) <- clone_levels
        }
        else {
            clone_colors <- NULL
        }
        all_colors <- c(bg_colors, clone_colors)
        if (length(all_colors) != length(plot_order)) {
            all_colors <- rep(all_colors, length.out = length(plot_order))
        }
        tmp_color_str <- paste(all_colors, collapse = " ")
        out_dir <- glue(file.path(figure_dir, "phylogeny/per_family_tls/{fam}/{tls}"))
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        if (exists("SpatialDimPython")) {
            dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)
            SpatialDimPython(obj = obj_ST, prefix = out_dir, plot_item = "label", plot_order = plot_order, 
                tmp_color = tmp_color_str, name = glue("{cur_sample}_{fam}_{tls}_clone_spatial"))
            cat("    Spatial plot saved to:", out_dir, "\n")
        }
        else {
            warning("SpatialDimPython not found, cannot generate spatial plot.")
        }
        cur_changeo <- seq_in_tls
        cur_changeo$clone_name <- paste(cur_changeo$TLS_raw, str_remove(cur_changeo$clone_family_07, 
            cur_sample), sep = "")
        cur_changeo$type <- cur_changeo$clone_subid
        cur_changeo$is_antibody <- grepl("^Ig", cur_changeo$clone_subid)
        selected_clone <- cur_changeo
        selected_clone$duplicate_count <- 1
        selected_clone$locus <- "IGH"
        selected_clone$c_call <- "IgG1"
        selected_clone$c_call[selected_clone$type != cur_sample] <- selected_clone$clone_subid[selected_clone$type != 
            cur_sample]
        selected_clone$clone_name <- paste(fam, tls, sep = "_")
        changeo_clone <- tryCatch({
            makeChangeoClone(selected_clone, text_fields = "c_call", clone = "clone_name", num_fields = "duplicate_count", 
                pad_end = TRUE, verbose = FALSE)
        }, error = function(e) NULL)
        if (!is.null(changeo_clone)) {
            changeo_clone@germline <- paste0(changeo_clone@germline, paste(rep("N", nchar(changeo_clone@data[1, 
                2]) - nchar(changeo_clone@germline)), collapse = ""))
            graph <- tryCatch({
                buildPhylipLineage(changeo_clone, phylip_exec, rm_temp = TRUE)
            }, error = function(e) NULL)
            if (!is.null(graph)) {
                graph_df <- createTreeLabel(graph, changeo_clone, selected_clone)
                for (i in 1:length(V(graph)$label)) {
                  if (V(graph)$label[i] %in% graph_df$to) {
                    V(graph)$label[i] <- graph_df[graph_df$to == V(graph)$label[i], "label"]
                  }
                }
                V(graph)$color <- "grey"
                V(graph)$color[V(graph)$name == "Germline"] <- "black"
                V(graph)$color[grepl("molecule", V(graph)$name)] <- "lightblue"
                anti_flag <- setNames(selected_clone$is_antibody, selected_clone$sequence_id)
                is_antibody_vertex <- anti_flag[V(graph)$name]
                is_antibody_vertex[is.na(is_antibody_vertex)] <- FALSE
                V(graph)$color[is_antibody_vertex] <- "#377537"
                size_table <- merge(V(graph)$name %>% matrix(ncol = 1), changeo_clone@data, by.x = "V1", 
                  by.y = "sequence_id", all = TRUE) %>% dplyr::select(V1, duplicate_count)
                size_table[is.na(size_table)] <- 1
                size_table <- size_table[match(V(graph)$name, size_table$V1), ]
                map_size <- function(x, min_out = 3, max_out = 9) {
                  x_log <- log1p(x)
                  x_min <- min(x_log)
                  x_max <- max(x_log)
                  if (x_min == x_max) 
                    return(rep((min_out + max_out)/2, length(x)))
                  scaled <- (x_log - x_min)/(x_max - x_min)
                  min_out + scaled * (max_out - min_out)
                }
                node_sizes <- map_size(size_table$duplicate_count)
                dir.create(dirname(glue("{out_dir}/tree.pdf")), recursive = TRUE, showWarnings = FALSE)
                pdf(glue("{out_dir}/tree.pdf"), width = 12, height = 10)
                lay <- layout.reingold.tilford(graph)
                plot(graph, layout = -lay[, 2:1], edge.arrow.mode = 0, vertex.frame.color = "black", 
                  vertex.label.color = "black", vertex.label.dist = 1.5, vertex.label.degree = 0, 
                  vertex.label = V(graph)$label, vertex.size = node_sizes, vertex.color = V(graph)$color)
                legend("topleft", legend = c("Germline", "Inferred", "Antibody", "other ig"), fill = c("black", 
                  "grey", "#377537", "lightblue"), cex = 0.8, bty = "n")
                dev.off()
                cat("    Tree saved to:", glue("{out_dir}/tree.pdf"), "\n")
            }
        }
        else {
            cat("    Failed to build tree.\n")
        }
    }
}

if (length(bin_info_list) > 0) {
    all_bin_info <- do.call(rbind, bin_info_list)
    all_bin_info <- all_bin_info %>% arrange(clone_family_07, distance_to_nearest_TLS)
    dir.create(dirname(file.path(table_dir, "phylogeny/bin_info_table.csv")), recursive = TRUE, 
        showWarnings = FALSE)
    write.csv(all_bin_info, file.path(table_dir, "phylogeny/bin_info_table.csv"), row.names = FALSE)
}

