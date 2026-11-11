PROJECT_ROOT <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, 'data')
source('Code/Extended_Analyses/R/io.R')
source('Code/Extended_Analyses/R/plot.R')
source('Code/Extended_Analyses/R/09_Bin50_functions.R')

#### 参数与输入 ====

library(data.table)

library(dplyr)

library(Matrix)

library(ggplot2)

library(patchwork)

library(scales)

raw_root <- dirs$sipeac

processed_dir <- file.path(dirs$processed, "R1_Q2c_SIPEAC_Bin50_selected_examples")

figure_dir <- file.path(dirs$figures, "R1_Q2c_SIPEAC_Bin50_selected_examples")

table_dir <- dirs$tables

log_dir <- file.path(dirs$stats_logs, "R1_Q2c_SIPEAC_Bin50_selected_examples")

METHOD_NO_AI <- "No AI / rule-based only"

METHOD_FULL <- "SIPEAC (with AI)"

METHOD_LEVELS <- c(METHOD_NO_AI, METHOD_FULL)

METHOD_COLOURS <- c(`No AI / rule-based only` = "#D84A5B", `SIPEAC (with AI)` = "#0B5EAA")

SAMPLES <- c("2931B", "RH16B", "2931T")

DOMINANCE_LEVELS <- c("0", "(0,25%]", "(25%,50%]", "(50%,75%]", "(75%,100%)", "100%")

DOMINANCE_LABELS <- c(`0` = "Never jointly\ndominant", `(0,25%]` = "Weak\n(0–25%]", `(25%,50%]` = "Moderate\n(25–50%]", 
    `(50%,75%]` = "Clear\n(50–75%]", `(75%,100%)` = "Strong\n(75–<100%)", `100%` = "Exclusive pair\n(100%)")

gt_dir_name <- if (dir.exists(file.path(raw_root, "grouth_truth"))) {
    "grouth_truth"
} else {
    "ground_truth"
}

input_paths <- data.frame(sample = SAMPLES, info = file.path(raw_root, "input_allelic", "bin50", paste0(SAMPLES, 
    "_info.txt")), ground_truth = file.path(raw_root, gt_dir_name, "bin50", paste0(SAMPLES, ".tsv")), 
    ai_score = file.path(raw_root, "2025_fold0_epoch100", paste0(SAMPLES, "_cdr3aa.csv")), stringsAsFactors = FALSE)

missing_inputs <- unlist(input_paths[-1], use.names = FALSE)

missing_inputs <- missing_inputs[!file.exists(missing_inputs)]

if (length(missing_inputs)) {
    stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"))
}

#### Bin50配对与优势度 ====

results <- vector("list", nrow(input_paths))

for (i in seq_len(nrow(input_paths))) {
    sample_id <- input_paths$sample[[i]]
    paths <- input_paths[i, , drop = FALSE]
    message("Reading ", sample_id, " raw Bin50 observations")
    raw <- fread(paths$info, select = c("clone_name", "clone_num", "binid"), showProgress = FALSE)
    spatial <- prepare_spatial_input(raw)
    counts <- spatial_clone_counts(spatial$filtered)
    gt <- load_ground_truth(paths$ground_truth)
    score_lookup <- read_score_lookup(paths$ai_score, counts)
    full_raw <- sipeac_pair_greedy(spatial$input_mat, gt, spatial$clone_stats, score_lookup = score_lookup, 
        tie_policy = "legacy")
    no_ai_raw <- sipeac_pair_greedy(spatial$input_mat, gt, spatial$clone_stats, tie_policy = "spatial")
    rm(score_lookup)
    full <- assemble_arm_predictions(full_raw, gt, counts, METHOD_FULL)
    no_ai <- assemble_arm_predictions(no_ai_raw, gt, counts, METHOD_NO_AI)
    if (anyDuplicated(full$Hclone) || anyDuplicated(no_ai$Hclone)) {
        stop("More than one final prediction per Hclone in ", sample_id)
    }
    agreement <- build_common_h_agreement(full, no_ai, gt, h_candidate_states(counts))
    common <- full %>% distinct(Hclone, .keep_all = TRUE) %>% transmute(Hclone, full_correct_prediction = correct, 
        joint_dominance, dominance_state) %>% inner_join(agreement, by = "Hclone")
    if (!all(common$full_correct_prediction == common$full_correct)) {
        stop("Full-arm correctness mismatch in ", sample_id)
    }
    analysis <- common %>% mutate(dominance = dominance_bin(joint_dominance), mechanism_eligible = case_when(joint_dominance < 
        1 & h_candidate_state == "Tied or multiple top candidates" ~ TRUE, joint_dominance == 1 & h_candidate_state == 
        "Single unambiguous top candidate" & both_100_support ~ TRUE, TRUE ~ FALSE)) %>% filter(mechanism_eligible)
    panel_counts <- analysis %>% group_by(dominance) %>% summarise(denominator = n(), no_ai_correct = sum(no_ai_correct), 
        sipeac_correct = sum(full_correct), .groups = "drop") %>% mutate(sample = sample_id, .before = 1)
    sample_result <- list(predictions = bind_rows(full, no_ai) %>% mutate(sample = sample_id, .before = 1), 
        agreement = agreement %>% mutate(sample = sample_id, .before = 1), analysis = analysis %>% mutate(sample = sample_id, 
            .before = 1), panel_counts = panel_counts, reference_n = nrow(full))
    results[[i]] <- sample_result
}

#### 共同重链与分组准确率 ====

predictions <- bind_rows(lapply(results, `[[`, "predictions"))

agreement <- bind_rows(lapply(results, `[[`, "agreement"))

analysis_data <- bind_rows(lapply(results, `[[`, "analysis"))

counts <- bind_rows(lapply(results, `[[`, "panel_counts")) %>% mutate(dominance = factor(dominance, levels = DOMINANCE_LEVELS, 
    ordered = TRUE))

sample_meta <- data.frame(sample = SAMPLES, reference_n = vapply(results, `[[`, integer(1), "reference_n"), 
    stringsAsFactors = FALSE)

sample_stats <- bind_rows(lapply(SAMPLES, function(sample_id) {
    d <- counts %>% filter(sample == sample_id)
    ambiguous <- d %>% filter(as.character(dominance) != "100%")
    data.frame(sample = sample_id, reference_n = sample_meta$reference_n[sample_meta$sample == sample_id], 
        plotted_n = sum(d$denominator), ambiguous_n = sum(ambiguous$denominator), no_ai_correct = sum(ambiguous$no_ai_correct), 
        sipeac_correct = sum(ambiguous$sipeac_correct), ambiguous_gain_pp = 100 * (sum(ambiguous$sipeac_correct) - 
            sum(ambiguous$no_ai_correct))/sum(ambiguous$denominator), stringsAsFactors = FALSE)
}))

plot_data <- bind_rows(counts %>% transmute(sample, dominance, denominator, method = METHOD_NO_AI, correct = no_ai_correct), 
    counts %>% transmute(sample, dominance, denominator, method = METHOD_FULL, correct = sipeac_correct)) %>% 
    mutate(method = factor(method, levels = METHOD_LEVELS), accuracy = correct/denominator)

#### 作图 ====

plots <- list()

for (sample_id in SAMPLES) {
    d <- plot_data %>% filter(sample == sample_id)
    wide <- counts %>% filter(sample == sample_id)
    equal_groups <- wide$dominance[wide$no_ai_correct == wide$sipeac_correct]
    label_data <- d %>% filter(!dominance %in% equal_groups)
    equal_data <- wide %>% filter(dominance %in% equal_groups)
    if (nrow(label_data)) {
        label_data <- label_data %>% group_by(dominance) %>% mutate(label = paste0(correct, "/", denominator), 
            label_y = if_else(accuracy == max(accuracy), pmin(accuracy + 0.065, 1.06), pmax(accuracy - 
                0.065, 0.015))) %>% ungroup()
    }
    equal_data <- equal_data %>% mutate(accuracy = sipeac_correct/denominator, label_y = pmin(accuracy + 
        0.065, 1.06), label = paste0("Both: ", sipeac_correct, "/", denominator))
    p <- ggplot(d, aes(dominance, accuracy, colour = method, group = method)) + geom_line(linewidth = 0.75) + 
        geom_point(size = 2.6) + geom_text(data = label_data, aes(y = label_y, label = label), size = 3, 
        show.legend = FALSE) + geom_text(data = equal_data, aes(dominance, label_y, label = label), inherit.aes = FALSE, 
        size = 3) + scale_colour_manual(values = METHOD_COLOURS, drop = FALSE) + scale_x_discrete(labels = DOMINANCE_LABELS, 
        drop = FALSE) + scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1.08), 
        breaks = c(0, 0.25, 0.5, 0.75, 1)) + labs(title = sample_id, x = "Joint IgH–IgL dominance", 
        y = "Accuracy", colour = NULL) + theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", 
        size = 12), axis.text.x = element_text(size = 7.5), legend.position = "top", legend.justification = "left")
    plots[[sample_id]] <- p
}

panel <- wrap_plots(plots, ncol = 1, guides = "collect") + plot_annotation(title = "Representative Bin50 samples") & 
    theme(legend.position = "top")

#### 输出 ====

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(predictions, file.path(processed_dir, "prediction_level_metrics.tsv"), sep = "\t")

fwrite(agreement, file.path(processed_dir, "common_H_agreement.tsv"), sep = "\t")

fwrite(analysis_data, file.path(processed_dir, "mechanism_analysis_data.tsv"), sep = "\t")

fwrite(counts %>% mutate(dominance = as.character(dominance)), file.path(table_dir, "R1_Q2c_SIPEAC_Bin50_selected_examples_counts.csv"))

fwrite(sample_stats, file.path(table_dir, "R1_Q2c_SIPEAC_Bin50_selected_examples_summary.csv"))

manifest <- bind_rows(input_paths %>% transmute(sample, role = "raw_spatial_info", path = info), input_paths %>% 
    transmute(sample, role = "ground_truth", path = ground_truth), input_paths %>% transmute(sample, 
    role = "frozen_ai_pair_scores", path = ai_score))

file_metadata <- file.info(manifest$path)

manifest$size_bytes <- file_metadata$size

manifest$modified <- format(file_metadata$mtime, "%Y-%m-%d %H:%M:%S %z")

fwrite(manifest, file.path(log_dir, "input_manifest.csv"))

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

figure_stem <- file.path(figure_dir, "SIPEAC_Bin50_selected_examples")

save_plot(panel, figure_stem, width = 7.2, height = 12, if_pdf = TRUE, if_png = TRUE)
ggsave(paste0(figure_stem, ".svg"), panel, device = svglite::svglite, width = 7.2, height = 12)

print(sample_stats)

message("Figures: ", figure_dir)

message("Tables: ", table_dir)

message("Processed data: ", processed_dir)
