suppressPackageStartupMessages({
  library(glue); library(dplyr)
})

source("Code/Extended_Analyses/R/expand_n_layer.R")

.pick_mat <- function(m, pat) {
  hit <- grepl(pat, colnames(m))
  if (sum(hit) == 1) {
    out <- matrix(
      m[, hit], ncol = 1,
      dimnames = list(rownames(m), colnames(m)[hit])
    )
    return(out)
  }
  m[, hit, drop = FALSE]
}

build_ground_truth <- function(sc_shared) {
  stopifnot(all(c("clone_name", "locus") %in% names(sc_shared)))
  if (!"IGH" %in% sc_shared$locus) return(data.frame())

  gt_mat <- table(sc_shared$cell_barcode, sc_shared$clone_name)
  selected_rows <- rowSums(gt_mat) >= 2
  if (sum(selected_rows) == 0) return(data.frame())

  gt_mat_filtered <- gt_mat[selected_rows, , drop = FALSE]
  gt_df <- list()
  for (cb in rownames(gt_mat_filtered)) {
    clones <- names(which(gt_mat_filtered[cb, ] > 0))
    if (any(grepl("IGH", clones)) && any(grepl("IGL|IGK", clones))) {
      gt_df[[cb]] <- expand.grid(
        Hclone = clones[grepl("IGH", clones)],
        Lclone = clones[grepl("IGL|IGK", clones)]
      )
    }
  }
  bind_rows(gt_df) %>%
    mutate(pair_name = paste(Hclone, Lclone, sep = "_")) %>%
    distinct(pair_name, .keep_all = TRUE)
}

prepare_sipeac_input <- function(cur_result, sipeac_root, cur_sample, bin_size = 1) {
  cur_result_st <- cur_result[!is.na(cur_result$binid), ]

  in_dir <- file.path(sipeac_root, "input_allelic", glue("bin{bin_size}"))
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  save_table(cur_result, file.path(in_dir, glue("{cur_sample}_info")), if_csv = FALSE)

  cur_result_clone <- cur_result %>%
    distinct(clone_id, locus, junction_aa, junction) %>%
    rename(cdr3 = junction, cdr3_aa = junction_aa)
  cdr_dir <- file.path(sipeac_root, "input_cdr3aa", glue("bin{bin_size}"))
  dir.create(cdr_dir, recursive = TRUE, showWarnings = FALSE)
  save_table(cur_result_clone, file.path(cdr_dir, glue("{cur_sample}_cdr3aa")), if_csv = FALSE)

  cur_result_sc <- cur_result[cur_result$type == "SingleCell" & cur_result$type_shared_main == "Shared", ]
  if ("IGH" %in% cur_result_sc$locus) {
    cur_result_sc$cell_barcode <- sapply(cur_result_sc$sequence_id, function(x) strsplit(x, "_contig")[[1]][1])
    gt_df <- build_ground_truth(cur_result_sc)
    if (nrow(gt_df) > 0) {
      gt_out <- file.path(sipeac_root, "grouth_truth", glue("bin{bin_size}"))
      dir.create(gt_out, recursive = TRUE, showWarnings = FALSE)
      write.table(gt_df, file.path(gt_out, glue("{cur_sample}.tsv")),
                  sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
    }
  }

  input_mat <- table(cur_result_st$binid, cur_result_st$clone_name)
  rep_dir <- file.path(sipeac_root, "input_repair", glue("bin{bin_size}"))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(input_mat, file.path(rep_dir, glue("{cur_sample}_mat.tsv")),
              row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
  invisible(NULL)
}

sipeac_pair_greedy <- function(input_mat, score_mat, test_mat, df_clone_num) {
  exclued_h <- c(); exclued_l <- c()
  output_pair_list <- list(); n_loop <- 1

  repeat {
    if (sum(!colnames(input_mat) %in% c(exclued_h, exclued_l)) == 1) break
    cur_input <- input_mat[, !colnames(input_mat) %in% c(exclued_h, exclued_l), drop = FALSE]

    x_1 <- .pick_mat(cur_input, "IGH"); x_2 <- .pick_mat(cur_input, "IGL|IGK")
    combined_condition <- (rowSums(x_1) >= 1) & (rowSums(x_2) >= 1)
    if (sum(combined_condition) == 0) break

    cur_input_HL <- if (sum(combined_condition) == 1) {
      out <- data.frame(matrix(cur_input[combined_condition, ], nrow = 1))
      rownames(out) <- rownames(cur_input)[combined_condition]
      colnames(out) <- colnames(cur_input); out
    } else {
      cur_input[combined_condition, , drop = FALSE]
    }
    cur_IgH_mat <- .pick_mat(cur_input_HL, "IGH")
    cur_IgL_mat <- .pick_mat(cur_input_HL, "IGL|IGK")

    output_pair_H <- c(); output_pair_L <- c(); output_pair_binid <- c()
    pair_count_list <- c(); score_list <- c()
    for (cur_bin in rownames(cur_IgH_mat)) {
      cur_max_IgH_count <- max(cur_IgH_mat[cur_bin, ])
      cur_max_IgH <- colnames(cur_IgH_mat)[which(cur_IgH_mat[cur_bin, ] == cur_max_IgH_count)]
      cur_max_IgL_count <- max(cur_IgL_mat[cur_bin, ])
      cur_max_IgL <- colnames(cur_IgL_mat)[which(cur_IgL_mat[cur_bin, ] == cur_max_IgL_count)]
      for (h in cur_max_IgH) for (l in cur_max_IgL) {
        output_pair_H <- c(output_pair_H, h); output_pair_L <- c(output_pair_L, l)
        output_pair_binid <- c(output_pair_binid, cur_bin)
        pair_count_list <- c(pair_count_list, cur_max_IgH_count * cur_max_IgL_count)
        h_n <- gsub("IGH_", "", h); l_n <- gsub("IGH_|IGK_|IGL_", "", l)
        score_list <- c(score_list, score_mat[glue("clone_{h_n}"), glue("clone_{l_n}")] %>% as.numeric())
      }
    }
    cur_output_pair <- data.frame(
      Hclone = output_pair_H, Lclone = output_pair_L, binid = output_pair_binid,
      pair_count_expand = pair_count_list, score = score_list
    )
    if (nrow(cur_output_pair) == 0) break

    cur_output_pair$pair_name <- paste(cur_output_pair$Hclone, cur_output_pair$Lclone, sep = "_")
    cur_output_pair <- cur_output_pair %>%
      group_by(pair_name) %>%
      mutate(pair_count = n(), pair_count_expand = sum(pair_count_expand)) %>%
      as.data.frame()
    cur_output_pair$pc_dot_s <- cur_output_pair$pair_count * cur_output_pair$score
    cur_output_pair <- cur_output_pair %>%
      left_join(distinct(test_mat, Hclone, Lclone_t), by = "Hclone")
    cur_output_pair <- cur_output_pair %>% left_join(df_clone_num, by = c("Hclone" = "clone_name"))
    cur_output_pair <- arrange(cur_output_pair, desc(clone_num), desc(pair_count))

    cur_output_pair_2 <- cur_output_pair %>%
      group_by(Lclone) %>%
      filter(pc_dot_s == max(pc_dot_s)) %>%
      as.data.table()
    cur_output_pair_2 <- distinct(cur_output_pair_2, pair_name, .keep_all = TRUE)
    cur_output_pair_2$clone_num <- as.numeric(cur_output_pair_2$clone_num)

    h_count <- table(cur_output_pair_2$Hclone) %>% sort() %>% rev()
    save_pair <- c()
    for (h in names(h_count)) {
      sub_pair <- cur_output_pair_2[cur_output_pair_2$Hclone == h, ]
      sub_max <- sub_pair[sub_pair$pc_dot_s == max(sub_pair$pc_dot_s), ]
      if (nrow(sub_max) == 1) save_pair <- c(save_pair, sub_max$pair_name)
    }
    if (length(save_pair) == 0) {
      output_pair_list[[n_loop]] <- cur_output_pair_2
      break
    }
    cur_output_pair_3 <- cur_output_pair_2[cur_output_pair_2$pair_name %in% save_pair, ]
    output_pair_list[[n_loop]] <- cur_output_pair_3
    exclued_h <- c(exclued_h, cur_output_pair_3$Hclone)
    exclued_l <- c(exclued_l, cur_output_pair_3$Lclone)
    n_loop <- n_loop + 1
  }

  if (length(output_pair_list) == 0) return(data.frame())
  output_pair_final <- bind_rows(output_pair_list) %>%
    mutate(pair_count_expand = pair_count_expand * pair_count) %>%
    group_by(Hclone) %>%
    filter(pair_count_expand == max(pair_count_expand)) %>%
    as.data.frame()
  output_pair_final
}

read_cdr3aa_scores <- function(cdr3aa_csv) {
  yang_score <- fread(cdr3aa_csv)
  colnames(yang_score) <- c("IGH_cloneIds", "IGL_cloneIds", "score")
  yang_score$IGH_cloneIds <- paste0("clone_", as.character(yang_score$IGH_cloneIds))
  yang_score$IGL_cloneIds <- paste0("clone_", as.character(yang_score$IGL_cloneIds))
  yang_score_mat <- tidyr::pivot_wider(yang_score, names_from = "IGL_cloneIds", values_from = "score") %>%
    as.data.frame()
  all_igh <- yang_score_mat$IGH_cloneIds
  all_igl <- colnames(yang_score_mat[, -1]) %>% as.character()
  score_mat <- as.matrix(yang_score_mat[, -1])
  rownames(score_mat) <- all_igh
  colnames(score_mat) <- all_igl
  list(score_mat = score_mat, score_long = yang_score)
}

load_ground_truth <- function(gt_path, clone_names = NULL) {
  test_mat <- read.table(gt_path, sep = "\t", header = TRUE, row.names = 1)
  if (!is.null(clone_names)) {
    test_mat <- test_mat[test_mat$Hclone %in% clone_names, ]
  }
  if (nrow(test_mat) == 0) return(data.frame())
  L_clone_t <- c()
  for (h in unique(test_mat$Hclone)) {
    L_clone_t <- c(L_clone_t, paste(test_mat$Lclone[test_mat$Hclone == h], collapse = ","))
  }
  names(L_clone_t) <- unique(test_mat$Hclone)
  test_mat$Lclone_t <- L_clone_t[test_mat$Hclone]
  test_mat$pair_name <- paste(test_mat$Hclone, test_mat$Lclone, sep = "_")
  test_mat
}

pairing_accuracy_bin1 <- function(bcr_sample_df, gt) {
  input <- bcr_sample_df
  input$clone_id <- as.character(input$clone_id)
  input$clone_name <- paste(input$locus, input$clone_id, sep = "_")

  input_mat <- table(input$bin1_id, input$clone_name)
  x_1 <- .pick_mat(input_mat, "IGH"); x_2 <- .pick_mat(input_mat, "IGL|IGK")
  combined_condition <- (rowSums(x_1) >= 1) & (rowSums(x_2) >= 1)
  if (sum(combined_condition) == 0) return(list(pairs = data.frame(), acc = NA, n = 0))

  cur_input_HL <- if (sum(combined_condition) == 1) {
    out <- data.frame(matrix(input_mat[combined_condition, ], nrow = 1))
    rownames(out) <- rownames(input_mat)[combined_condition]
    colnames(out) <- colnames(input_mat); out
  } else {
    input_mat[combined_condition, , drop = FALSE]
  }
  cur_IgH_mat <- .pick_mat(cur_input_HL, "IGH")
  cur_IgL_mat <- .pick_mat(cur_input_HL, "IGL|IGK")

  output_pair_H <- c(); output_pair_L <- c(); output_pair_binid <- c()
  for (cur_bin in rownames(cur_IgH_mat)) {
    cur_IgH_idx <- which(cur_IgH_mat[cur_bin, ] >= 1)
    cur_IgL_idx <- which(cur_IgL_mat[cur_bin, ] >= 1)
    for (h in names(cur_IgH_idx)) for (l in names(cur_IgL_idx)) {
      output_pair_H <- c(output_pair_H, h)
      output_pair_L <- c(output_pair_L, l)
      output_pair_binid <- c(output_pair_binid, cur_bin)
    }
  }
  cur_output_pair <- data.frame(Hclone = output_pair_H, Lclone = output_pair_L, binid = output_pair_binid)
  if (nrow(cur_output_pair) == 0) return(list(pairs = data.frame(), acc = NA, n = 0))

  cur_output_pair$pair_name <- paste(cur_output_pair$Hclone, cur_output_pair$Lclone, sep = "_")
  cur_output_pair <- distinct(cur_output_pair, pair_name, .keep_all = TRUE)
  cur_output_pair <- cur_output_pair %>% left_join(distinct(gt, Hclone, Lclone_t), by = "Hclone")

  cur_output_pair_test <- cur_output_pair[cur_output_pair$Hclone %in% gt$Hclone, ]
  if (nrow(cur_output_pair_test) == 0) return(list(pairs = data.frame(), acc = NA, n = 0))

  cur_output_pair_test <- cur_output_pair_test %>%
    mutate(
      bin1_correct = pair_name %in% gt$pair_name,
      if_true = ifelse(bin1_correct, "True", "False")
    )
  n_correct <- sum(cur_output_pair_test$bin1_correct)
  acc <- n_correct / nrow(cur_output_pair_test)
  list(
    pairs = cur_output_pair_test,
    acc = acc,
    n_correct = n_correct,
    n = nrow(cur_output_pair_test)
  )
}

compare_pair_confusion <- function(bin1_pairs, sipeac_pairs, gt) {
  gt$pair_name <- paste(gt$Hclone, gt$Lclone, sep = "_")
  h_intersect <- intersect(unique(bin1_pairs$Hclone), unique(sipeac_pairs$Hclone))
  gt <- gt[gt$Hclone %in% h_intersect, ]
  if (nrow(gt) == 0) return(data.frame())

  bin1_hit <- gt$pair_name %in% unique(bin1_pairs$pair_name[bin1_pairs$bin1_correct %in% TRUE])
  sip_hit  <- gt$pair_name %in% unique(sipeac_pairs$pair_name[sipeac_pairs$sipeac_correct %in% TRUE])
  data.frame(
    pair_name = gt$pair_name,
    bin1_correct = bin1_hit,
    sipeac_correct = sip_hit,
    n_bin1 = nrow(bin1_pairs),
    n_sip = nrow(sipeac_pairs),
    bin1_acc_all = round(sum(bin1_pairs$bin1_correct) / nrow(bin1_pairs), 4),
    sip_acc_all  = round(sum(sipeac_pairs$sipeac_correct) / nrow(sipeac_pairs), 4)
  )
}

sipeac_accuracy_summary <- function(output_pair_final, test_mat, df_bin_num,
                                     cur_sample, bin_size, thre_list = c(1)) {
  pred_all <- list()
  for (cur_thre in thre_list) {
    output_pair_final1 <- output_pair_final[output_pair_final$bin_num >= cur_thre, ]
    output_pair_test <- output_pair_final1[output_pair_final1$Hclone %in% test_mat$Hclone, ]
    df_bin_num1 <- df_bin_num[df_bin_num$bin_num >= cur_thre, ]
    n_H_clone <- nrow(df_bin_num1[grepl("IGH", df_bin_num1$clone_name), ])
    n_H_clone_pred <- length(unique(output_pair_final1$Hclone))
    n_test <- nrow(test_mat[test_mat$Hclone %in% output_pair_final1$Hclone, ])
    if (nrow(output_pair_test) == 0) next
    n_correct <- sum(output_pair_test$pair_name %in% test_mat$pair_name)
    n_evaluable <- nrow(output_pair_test)
    acc <- n_correct / n_evaluable
    eff <- n_H_clone_pred / n_H_clone
    pred_all[[glue("{cur_thre}")]] <- data.frame(
      cur_sample = cur_sample, n_H_clone = n_H_clone,
      n_L_clone = nrow(df_bin_num1[!grepl("IGH", df_bin_num1$clone_name), ]),
      n_H_clone_pred = n_H_clone_pred, eff = eff, n_test = n_test,
      n_evaluable_predictions = n_evaluable,
      n_correct_predictions = n_correct,
      threshold = cur_thre, bin_size = bin_size, acc = acc
    )
  }
  bind_rows(pred_all)
}

compare_sipeac_repair <- function(compare_df,
                                  bins = c(20, 50, 110),
                                  thresholds = c(1, 2, 5, 10, 15, 20)) {
  stopifnot(all(c("cur_sample","bin_size","threshold","acc","repair_acc") %in% names(compare_df)))
  sub <- compare_df %>%
    filter(bin_size %in% bins, threshold %in% thresholds,
           !is.na(acc), !is.na(repair_acc))
  bind_rows(
    sub %>% transmute(cur_sample, bin_size, threshold, method = "SIPEAC", value = acc),
    sub %>% transmute(cur_sample, bin_size, threshold, method = "repair", value = repair_acc)
  )
}

get_shared_info <- function(df, source = "clone_id", target_list = c()) {
  df[[source]] <- as.character(df[[source]])
  for (cur_target in target_list) {
    shared_type_info  <- c()
    shared_type_info2 <- c()
    for (cur.id in unique(df[[source]])) {
      sample_list <- unique(df[df[[source]] == cur.id, ][[cur_target]]) %>% sort()
      if (length(sample_list) > 1) {
        sl <- paste(sample_list, collapse = ", ")
        shared_type <- paste0("Shared_", sl)
        st2 <- "Shared"
      } else {
        shared_type <- paste0("Private_", sample_list)
        st2 <- "Private"
      }
      shared_type_info  <- c(shared_type_info, shared_type)
      shared_type_info2 <- c(shared_type_info2, st2)
    }
    names(shared_type_info)  <- unique(df[[source]])
    names(shared_type_info2) <- unique(df[[source]])
    df[[glue("{cur_target}_shared_sub")]]  <- shared_type_info[df[[source]]]
    df[[glue("{cur_target}_shared_main")]] <- shared_type_info2[df[[source]]]
  }
  clone_num <- table(df[[source]])
  df$clone_num <- clone_num[df[[source]]]
  return(df)
}

createTreeLabel <- function(graph, changeo_clone, clone) {
  graph_df <- graph %>% igraph::as_data_frame()

  graph_df <- rbind(graph_df, c("Germline", "Germline", 0, ""))

  dfs_ordered_vertex <- graph %>% dfs(root = "Germline", order = TRUE)
  dfs_ordered_vertex_table <- cbind(dfs_ordered_vertex$order$label, c(1:length(dfs_ordered_vertex$order)))
  colnames(dfs_ordered_vertex_table) <- c("to", "order")
  graph_df <- merge(graph_df, dfs_ordered_vertex_table, by = "to")

  graph_df_isotype <- changeo_clone@data %>% as_tibble() %>% dplyr::select(sequence_id, c_call)
  graph_df_isotype$c_call <- gsub("IGH", "", graph_df_isotype$c_call)
  graph_df <- merge(graph_df, graph_df_isotype, by.x = "to", by.y = "sequence_id", all.x = TRUE)
  graph_df[is.na(graph_df$c_call), "c_call"] <- substr(graph_df[is.na(graph_df$c_call), "to"], 1, 2)

  clone_selected <- clone %>%
    ungroup() %>%
    dplyr::select(sequence_id, type) %>%
    dplyr::mutate(label_ab = dplyr::if_else(type != "other_ig", "other_ig", type))

  graph_df <- merge(graph_df, clone_selected, by.x = "to", by.y = "sequence_id", all.x = TRUE, all.y = FALSE)
  graph_df$type[is.na(graph_df$type)] <- "NA"

  graph_df <- graph_df %>%
    mutate(label = paste(order, c_call, sep = "|"))
  graph_df$label <- gsub("\\*00", "", graph_df$label)
  graph_df$label <- gsub(",", ".", graph_df$label)

  return(graph_df)
}
