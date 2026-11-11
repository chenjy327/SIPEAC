#### ssDNA与三分辨率汇总图 ====
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/R/io.R")
source("Code/Extended_Analyses/R/plot.R")
input_dir <- file.path(data_dir, "hzy")
processed_dir <- file.path(dirs$processed, "hzy")
table_dir <- file.path(dirs$tables, "hzy", "ssDNA")
figure_dir <- file.path(dirs$figures, "hzy")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(glue)
library(rstatix)

base_dir <- table_dir

bin_sizes <- c(20L, 50L, 110L)

df_old <- bind_rows(lapply(bin_sizes, function(bs) {
    f <- file.path(base_dir, glue("acc_eff_old_gold_bin{bs}.csv"))
    if (!file.exists(f)) {
        warning("缺少文件: ", f)
        return(NULL)
    }
    read.csv(f) %>% mutate(group = glue("old bin{bs}"))
}))

df_ssdna <- read.csv(file.path(base_dir, "acc_eff_ssDNA_vs_scRNA_gold_bin.csv")) %>% mutate(group = "ssDNA")

bin_labels_present <- unique(df_old$group)

group_levels <- c("ssDNA", bin_labels_present)

df_plot <- bind_rows(df_old, df_ssdna) %>% filter(!is.na(acc)) %>% mutate(group = factor(group, 
    levels = group_levels))

n_by_group <- df_plot %>% count(group)

fit <- aov(acc ~ group, data = df_plot)

anova_p <- summary(fit)[[1]][["Pr(>F)"]][1]

tukey <- tukey_hsd(df_plot, acc ~ group) %>% filter(group1 == "ssDNA" | group2 == "ssDNA") %>% mutate(group1 = factor(group1, 
    levels = group_levels), group2 = factor(group2, levels = group_levels), p.adj.format = format.pval(p.adj, 
    digits = 3, eps = 1e-04)) %>% arrange(group2)

y_max <- max(df_plot$acc, na.rm = TRUE)

y_rng <- diff(range(df_plot$acc, na.rm = TRUE))

if (is.na(y_rng) || y_rng == 0) y_rng <- 0.1

tukey <- tukey %>% mutate(y.position = y_max + y_rng * (0.1 + 0.1 * (seq_len(n()) - 1)))

y_top <- max(tukey$y.position) + y_rng * 0.18

bin_colors <- setNames(colorRampPalette(c("#FDDBC7", "#D6604D"))(length(bin_labels_present)), bin_labels_present)

p <- ggplot(df_plot, aes(x = group, y = acc, fill = group)) + geom_boxplot(width = 0.6, outlier.shape = NA, 
    alpha = 0.55, color = "grey25", linewidth = 0.5) + geom_jitter(width = 0.12, size = 1.8, alpha = 0.65, 
    color = "grey15", shape = 16) + stat_summary(fun = median, geom = "point", shape = 23, size = 3, 
    fill = "white", color = "black", stroke = 0.8) + stat_pvalue_manual(tukey, label = "p.adj.format", 
    xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0.02, size = 3.6, 
    inherit.aes = FALSE) + scale_fill_manual(values = c(ssDNA = "#92C5DE", bin_colors)) + scale_y_continuous(limits = c(NA, 
    y_top), expand = expansion(mult = c(0.02, 0))) + labs(x = NULL, y = "Accuracy (acc)", title = "BCR pairing accuracy: ssDNA vs old gold standard", 
    subtitle = paste0(paste(sprintf("%s (n=%d)", n_by_group$group, n_by_group$n), collapse = " | "), 
        "\nOne-way ANOVA p = ", format.pval(anova_p, digits = 3, eps = 1e-04), "   (post-hoc: Tukey HSD)")) + 
    theme_bw(base_size = 13) + theme(legend.position = "none", plot.title = element_text(face = "bold"), 
    plot.subtitle = element_text(color = "grey35", size = 10), axis.text.x = element_text(size = 11, 
        color = "black"), panel.grid.minor = element_blank())

print(p)

dir.create(dirname(file.path(figure_dir, "acc_boxplot_ssDNA_vs_bins_anova_tukey.png")), recursive = TRUE, 
    showWarnings = FALSE)

ggsave(file.path(figure_dir, "acc_boxplot_ssDNA_vs_bins_anova_tukey.png"), p, width = 7.5, height = 5.5, 
    dpi = 300)

dir.create(dirname(file.path(figure_dir, "acc_boxplot_ssDNA_vs_bins_anova_tukey.pdf")), recursive = TRUE, 
    showWarnings = FALSE)

ggsave(file.path(figure_dir, "acc_boxplot_ssDNA_vs_bins_anova_tukey.pdf"), p, width = 7.5, height = 5.5)

