.pick_mat <- function(m, pattern) {
    hit <- grepl(pattern, colnames(m))
    if (sum(hit) == 1L) {
        out <- matrix(m[, hit], ncol = 1L)
        colnames(out) <- colnames(m)[hit]
        rownames(out) <- rownames(m)
        return(out)
    }
    m[, hit, drop = FALSE]
}

load_ground_truth <- function(path) {
    gt <- read.table(path, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    if (!all(c("Hclone", "Lclone") %in% names(gt))) 
        stop("Invalid ground truth: ", path)
    if (!nrow(gt)) {
        return(gt)
    }
    l_map <- vapply(unique(gt$Hclone), function(h) paste(gt$Lclone[gt$Hclone == h], collapse = ","), 
        character(1))
    gt$Lclone_t <- unname(l_map[gt$Hclone])
    gt$pair_name <- paste(gt$Hclone, gt$Lclone, sep = "_")
    gt
}

prepare_spatial_input <- function(raw) {
    x <- as.data.table(copy(raw))[!is.na(binid) & nzchar(binid) & !is.na(clone_name) & grepl("^IGH|^IGK|^IGL", 
        clone_name)]
    if (!nrow(x)) 
        stop("No spatial immunoglobulin observations")
    pair_counts <- x[, .(count = .N), by = .(binid, clone_name)]
    bin_levels <- sort(unique(pair_counts$binid))
    clone_levels <- sort(unique(pair_counts$clone_name))
    input_mat <- sparseMatrix(i = match(pair_counts$binid, bin_levels), j = match(pair_counts$clone_name, 
        clone_levels), x = pair_counts$count, dims = c(length(bin_levels), length(clone_levels)), dimnames = list(bin_levels, 
        clone_levels))
    clone_stats <- x[, .(clone_num = suppressWarnings(max(as.numeric(clone_num), na.rm = TRUE)), bin_num = uniqueN(binid)), 
        by = clone_name]
    clone_stats[!is.finite(clone_num), `:=`(clone_num, 0)]
    list(input_mat = input_mat, clone_stats = as.data.frame(clone_stats), filtered = x)
}

read_score_lookup <- function(path, counts, chunk_size = 1000000L) {
    h <- counts[chain == "H", .(binid, h_id = sub("^IGH_", "", clone_name))]
    l <- counts[chain == "L", .(binid, l_id = sub("^(IGL|IGK)_", "", clone_name))]
    possible <- unique(merge(h, l, by = "binid", allow.cartesian = TRUE)[, .(h_id, l_id)])
    possible_keys <- paste(possible$h_id, possible$l_id, sep = "\r")
    connection <- file(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    readLines(connection, n = 1L)
    retained <- list()
    index <- 1L
    repeat {
        lines <- readLines(connection, n = chunk_size)
        if (!length(lines)) 
            break
        chunk <- fread(text = paste(lines, collapse = "\n"), header = FALSE, col.names = c("h_id", "l_id", 
            "score"), showProgress = FALSE)
        chunk[, `:=`(h_id = as.character(h_id), l_id = as.character(l_id))]
        keep <- paste(chunk$h_id, chunk$l_id, sep = "\r") %chin% possible_keys
        if (any(keep)) {
            retained[[index]] <- chunk[keep]
            index <- index + 1L
        }
    }
    scores <- unique(rbindlist(retained))
    setkey(scores, h_id, l_id)
    scores
}

sipeac_pair_greedy <- function(input_mat, test_mat, clone_stats, score_lookup = NULL, tie_policy = c("legacy", 
    "spatial")) {
    tie_policy <- match.arg(tie_policy)
    excluded_h <- character()
    excluded_l <- character()
    output <- list()
    iteration <- 1L
    repeat {
        keep <- !colnames(input_mat) %in% c(excluded_h, excluded_l)
        if (sum(keep) <= 1L) 
            break
        current <- input_mat[, keep, drop = FALSE]
        h_all <- .pick_mat(current, "IGH")
        l_all <- .pick_mat(current, "IGL|IGK")
        if (!ncol(h_all) || !ncol(l_all)) 
            break
        has_both <- rowSums(h_all) >= 1 & rowSums(l_all) >= 1
        if (!any(has_both)) 
            break
        current_hl <- if (sum(has_both) == 1L) {
            out <- data.frame(matrix(current[has_both, ], nrow = 1L))
            rownames(out) <- rownames(current)[has_both]
            colnames(out) <- colnames(current)
            out
        }
        else {
            current[has_both, , drop = FALSE]
        }
        h_mat <- .pick_mat(current_hl, "IGH")
        l_mat <- .pick_mat(current_hl, "IGL|IGK")
        rows <- vector("list", nrow(h_mat))
        row_index <- 1L
        for (bin_id in rownames(h_mat)) {
            h_max <- max(h_mat[bin_id, ])
            l_max <- max(l_mat[bin_id, ])
            h_names <- colnames(h_mat)[which(h_mat[bin_id, ] == h_max)]
            l_names <- colnames(l_mat)[which(l_mat[bin_id, ] == l_max)]
            rows[[row_index]] <- expand.grid(Hclone = h_names, Lclone = l_names, stringsAsFactors = FALSE) %>% 
                mutate(binid = bin_id, pair_count_expand = h_max * l_max)
            row_index <- row_index + 1L
        }
        candidates <- bind_rows(rows)
        if (!nrow(candidates)) 
            break
        if (is.null(score_lookup)) {
            candidates$score <- 1
        }
        else {
            query <- data.table(h_id = sub("^IGH_", "", candidates$Hclone), l_id = sub("^(IGL|IGK)_", 
                "", candidates$Lclone))
            candidates$score <- score_lookup[query, on = .(h_id, l_id), x.score]
            if (anyNA(candidates$score)) 
                stop("Missing AI score for a candidate pair")
        }
        candidates <- candidates %>% mutate(pair_name = paste(Hclone, Lclone, sep = "_")) %>% group_by(pair_name) %>% 
            mutate(pair_count = n(), pair_count_expand = sum(pair_count_expand)) %>% ungroup() %>% mutate(pc_dot_s = pair_count * 
            score) %>% left_join(distinct(test_mat, Hclone, Lclone_t), by = "Hclone") %>% left_join(clone_stats, 
            by = c(Hclone = "clone_name")) %>% arrange(desc(clone_num), desc(pair_count))
        if (tie_policy == "spatial") {
            selected_l <- candidates %>% arrange(desc(pc_dot_s), desc(pair_count_expand), desc(pair_count), 
                desc(clone_num), pair_name) %>% group_by(Lclone) %>% slice_head(n = 1L) %>% ungroup() %>% 
                distinct(pair_name, .keep_all = TRUE)
            selected <- selected_l %>% arrange(desc(pc_dot_s), desc(pair_count_expand), desc(pair_count), 
                desc(clone_num), pair_name) %>% group_by(Hclone) %>% slice_head(n = 1L) %>% ungroup()
        }
        else {
            selected_l <- candidates %>% group_by(Lclone) %>% filter(pc_dot_s == max(pc_dot_s)) %>% ungroup() %>% 
                distinct(pair_name, .keep_all = TRUE)
            h_count <- sort(table(selected_l$Hclone), decreasing = TRUE)
            saved_pairs <- unlist(lapply(names(h_count), function(h) {
                current <- selected_l[selected_l$Hclone == h, , drop = FALSE]
                best <- current[current$pc_dot_s == max(current$pc_dot_s), , drop = FALSE]
                if (nrow(best) == 1L) 
                  best$pair_name
                else character()
            }), use.names = FALSE)
            if (!length(saved_pairs)) {
                output[[iteration]] <- selected_l
                break
            }
            selected <- selected_l %>% filter(pair_name %in% saved_pairs)
        }
        if (!nrow(selected)) 
            break
        output[[iteration]] <- selected
        excluded_h <- c(excluded_h, selected$Hclone)
        excluded_l <- c(excluded_l, selected$Lclone)
        iteration <- iteration + 1L
    }
    if (!length(output)) {
        return(data.frame())
    }
    final <- bind_rows(output) %>% mutate(pair_count_expand = pair_count_expand * pair_count)
    if (tie_policy == "spatial") {
        final <- final %>% arrange(desc(pair_count_expand), desc(pc_dot_s), desc(pair_count), desc(clone_num), 
            pair_name) %>% group_by(Hclone) %>% slice_head(n = 1L) %>% ungroup()
    }
    else {
        final <- final %>% group_by(Hclone) %>% filter(pair_count_expand == max(pair_count_expand)) %>% 
            ungroup()
    }
    as.data.frame(final)
}

spatial_clone_counts <- function(raw) {
    x <- as.data.table(copy(raw))[!is.na(binid) & nzchar(binid) & !is.na(clone_name) & grepl("^IGH|^IGK|^IGL", 
        clone_name)]
    x[, `:=`(chain, fifelse(grepl("^IGH", clone_name), "H", "L"))]
    out <- x[, .(clone_count = .N), by = .(binid, chain, clone_name)]
    out[, `:=`(chain_total = sum(clone_count), chain_max = max(clone_count), n_top = sum(clone_count == 
        max(clone_count))), by = .(binid, chain)]
    out[, `:=`(abundance = clone_count/chain_total, is_top = clone_count == chain_max)]
    setorder(out, binid, chain, clone_name)
    out[]
}

annotate_prediction_dominance <- function(pred, counts) {
    pred_dt <- as.data.table(as.data.frame(pred))
    if (!"pair_name" %in% names(pred_dt)) {
        pred_dt[, `:=`(pair_name, paste(Hclone, Lclone, sep = "_"))]
    }
    pred_dt[, `:=`(prediction_id, .I)]
    if (!nrow(pred_dt)) {
        return(as.data.frame(pred_dt[, `:=`(prediction_id, NULL)]))
    }
    h <- counts[chain == "H", .(binid, Hclone = clone_name, pH = abundance, H_is_top = is_top, n_top_H = n_top)]
    l <- counts[chain == "L", .(binid, Lclone = clone_name, pL = abundance, L_is_top = is_top, n_top_L = n_top)]
    h_obs <- merge(pred_dt[, .(prediction_id, Hclone, Lclone)], h, by = "Hclone", allow.cartesian = TRUE)
    hl <- merge(h_obs, l, by = c("binid", "Lclone"), allow.cartesian = TRUE)
    if (nrow(hl)) {
        hl[, `:=`(joint_abundance = pmin(pH, pL), jointly_top = H_is_top & L_is_top, jointly_unique_top = H_is_top & 
            L_is_top & n_top_H == 1L & n_top_L == 1L, both_100 = pH == 1 & pL == 1)]
        dominance <- hl[, {
            top_values <- joint_abundance[jointly_top]
            list(joint_dominance = if (length(top_values)) max(top_values) else 0, joint_dominance_all = max(joint_abundance), 
                dominance_state = if (any(both_100)) {
                  "Both 100%"
                } else if (any(jointly_unique_top)) {
                  "Unique dominant"
                } else if (any(jointly_top)) {
                  "Tied dominant"
                } else {
                  "Not jointly dominant"
                }, n_joint_bins = uniqueN(binid), n_joint_top_bins = uniqueN(binid[jointly_top]))
        }, by = prediction_id]
    }
    else {
        dominance <- data.table(prediction_id = integer(), joint_dominance = numeric(), joint_dominance_all = numeric(), 
            dominance_state = character(), n_joint_bins = integer(), n_joint_top_bins = integer())
    }
    out <- merge(pred_dt, dominance, by = "prediction_id", all.x = TRUE, sort = FALSE)
    out[is.na(joint_dominance), `:=`(joint_dominance = 0, joint_dominance_all = 0, dominance_state = "Not jointly dominant", 
        n_joint_bins = 0L, n_joint_top_bins = 0L)]
    setorder(out, prediction_id)
    out[, `:=`(prediction_id, NULL)]
    as.data.frame(out)
}

h_candidate_states <- function(counts) {
    h <- counts[chain == "H" & is_top, .(binid, Hclone = clone_name, pH = abundance, n_top_H = n_top)]
    l <- counts[chain == "L" & is_top, .(binid, Lclone = clone_name, pL = abundance, n_top_L = n_top)]
    candidates <- merge(h, l, by = "binid", allow.cartesian = TRUE)
    if (!nrow(candidates)) {
        return(data.frame())
    }
    candidates[, `:=`(supporting_top_tie, n_top_H > 1L | n_top_L > 1L)]
    out <- candidates[, .(n_candidate_L = uniqueN(Lclone), any_supporting_top_tie = any(supporting_top_tie), 
        both_100_support = any(pH == 1 & pL == 1)), by = Hclone]
    out[, `:=`(h_candidate_state, fifelse(n_candidate_L == 1L & !any_supporting_top_tie, "Single unambiguous top candidate", 
        "Tied or multiple top candidates"))]
    as.data.frame(out)
}

assemble_arm_predictions <- function(pred, gt, counts, method) {
    gt_pairs <- gt %>% mutate(pair_name = paste(Hclone, Lclone, sep = "_")) %>% distinct(Hclone, Lclone, 
        pair_name)
    out <- as.data.frame(pred) %>% mutate(pair_name = paste(Hclone, Lclone, sep = "_")) %>% distinct(pair_name, 
        .keep_all = TRUE) %>% filter(Hclone %in% gt_pairs$Hclone) %>% mutate(across(any_of(c("Hclone", 
        "Lclone", "Lclone_t", "pair_name")), as.character), correct = pair_name %in% gt_pairs$pair_name, 
        method = method)
    annotate_prediction_dominance(out, counts)
}

build_common_h_agreement <- function(full, no_ai, gt, h_states) {
    gt_pairs <- unique(paste(gt$Hclone, gt$Lclone, sep = "_"))
    arm_summary <- function(pred, prefix) {
        out <- unique(as.data.table(pred)[, .(Hclone, Lclone)])
        out[, `:=`(pair_name, paste(Hclone, Lclone, sep = "_"))]
        out[, `:=`(correct, pair_name %in% gt_pairs)]
        out <- out[, .(L_set = paste(sort(unique(Lclone)), collapse = "|"), n_predictions = uniqueN(Lclone), 
            correct = any(correct)), by = Hclone]
        setnames(out, c("L_set", "n_predictions", "correct"), paste0(prefix, c("_L_set", "_n_predictions", 
            "_correct")))
        out
    }
    out <- merge(arm_summary(full, "full"), arm_summary(no_ai, "no_ai"), by = "Hclone")
    out[, `:=`(set_agreement, full_L_set == no_ai_L_set)]
    out[, `:=`(transition, fcase(full_correct & no_ai_correct, "Both correct", full_correct & !no_ai_correct, 
        "AI correction", !full_correct & no_ai_correct, "AI harm", default = "Both incorrect"))]
    states <- as.data.table(h_states)
    if (nrow(states)) {
        out <- merge(out, states[, .(Hclone, n_candidate_L, h_candidate_state, both_100_support)], by = "Hclone", 
            all.x = TRUE, sort = FALSE)
    }
    if (!"h_candidate_state" %in% names(out)) 
        out[, `:=`(h_candidate_state, NA_character_)]
    if (!"both_100_support" %in% names(out)) 
        out[, `:=`(both_100_support, NA)]
    out[is.na(h_candidate_state), `:=`(h_candidate_state, "No raw joint-top candidate")]
    out[is.na(both_100_support), `:=`(both_100_support, FALSE)]
    as.data.frame(out)
}

dominance_bin <- function(x) {
    labels <- case_when(x == 0 ~ "0", x <= 0.25 ~ "(0,25%]", x <= 0.5 ~ "(25%,50%]", x <= 0.75 ~ "(50%,75%]", 
        x < 1 ~ "(75%,100%)", TRUE ~ "100%")
    factor(labels, levels = DOMINANCE_LEVELS, ordered = TRUE)
}
