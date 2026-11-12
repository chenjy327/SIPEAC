#### TCGA TLS与Tumor/Paratumor差异 ====
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/utils_io.R")
source("Code/Extended_Analyses/utils_plot.R")
input_dir <- file.path(data_dir, "hzy")
processed_dir <- file.path(dirs$processed, "hzy")
table_dir <- file.path(dirs$tables, "hzy")
figure_dir <- file.path(dirs$figures, "hzy")
library(dplyr)
library(limma)

#### TLS limma ====
exp <- read.csv(file.path(input_dir, "data/TCGA_LIHC_count_symbol.csv"), header = T, check.names = F, 
    row.names = 1)

group <- read.csv(file.path(input_dir, "data/TCGA_group.csv"), header = T, check.names = F)

df_HE <- group %>% filter(TissueType == "Tumor" & TLS != "none")

HE_sample <- pull(df_HE, samplename)

pro_ref <- read.csv(file.path(input_dir, "reference_windows/pro_coding_genes.csv"), header = T, 
    row.names = 1)

pro_coding_genes <- pull(pro_ref, gene_name)

exp_Tumor <- exp[intersect(pro_coding_genes, rownames(exp)), HE_sample]

exp_Tumor[1:5, 1:5]

exp_Tumor <- log2(exp_Tumor + 1)

group_list <- factor(df_HE[, 3], ordered = F)

design <- model.matrix(~0 + group_list)

colnames(design) <- levels(group_list)

rownames(design) <- colnames(exp_Tumor)

design

fit <- lmFit(exp_Tumor, design)

cont.matrix <- makeContrasts(contrasts = "positive-negative", levels = design)

cont.matrix

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

tmpOut <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "logFC")

dif <- tmpOut[tmpOut$P.Value <= 0.05 & abs(tmpOut$logFC) > 0, ]

dim(dif)

dir.create(dirname(file.path(table_dir, "TCGA/diff_TCGA_TLS_limma_log.csv")),
    recursive = TRUE, showWarnings = FALSE)

write.csv(tmpOut, file = file.path(table_dir, "TCGA/diff_TCGA_TLS_limma_log.csv"))

#### Tumor/Paratumor limma ====
exp <- read.csv(file.path(input_dir, "data/TCGA_LIHC_count_symbol.csv"), header = T, check.names = F, 
    row.names = 1)

group <- read.csv(file.path(input_dir, "data/TCGA_group.csv"), header = T, check.names = F)

df_TP <- group

TP_sample <- pull(df_TP, samplename)

pro_ref <- read.csv(file.path(input_dir, "data/pro_coding_genes.csv"), header = T, row.names = 1)

pro_coding_genes <- pull(pro_ref, gene_name)

exp_Tumor <- exp[intersect(pro_coding_genes, rownames(exp)), TP_sample]

exp_Tumor[1:5, 1:5]

exp_Tumor <- log2(exp_Tumor + 1)

group_list <- factor(df_TP[, 2], ordered = F)

design <- model.matrix(~0 + group_list)

colnames(design) <- levels(group_list)

rownames(design) <- colnames(exp_Tumor)

design

fit <- lmFit(exp_Tumor, design)

cont.matrix <- makeContrasts(contrasts = "Tumor-Paratumor", levels = design)

cont.matrix

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

tmpOut <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "logFC")

dif <- tmpOut[tmpOut$P.Value <= 0.05 & abs(tmpOut$logFC) > 0, ]

dim(dif)

dir.create(dirname(file.path(table_dir, "TCGA/diff_TCGA_TP_limma_log.csv")),
    recursive = TRUE, showWarnings = FALSE)

write.csv(tmpOut, file = file.path(table_dir, "TCGA/diff_TCGA_TP_limma_log.csv"))

#### TLS Wilcoxon ====
exp <- read.csv(file.path(input_dir, "data/TCGA_LIHC_count_symbol.csv"), header = T, check.names = F, 
    row.names = 1)

group <- read.csv(file.path(input_dir, "data/TCGA_group.csv"), header = T, check.names = F)

df_HE <- group %>% filter(TissueType == "Tumor" & TLS != "none")

HE_sample <- pull(df_HE, samplename)

pro_ref <- read.csv(file.path(input_dir, "data/pro_coding_genes.csv"), header = T, row.names = 1)

pro_coding_genes <- pull(pro_ref, gene_name)

exp_Tumor <- exp[intersect(pro_coding_genes, rownames(exp)), HE_sample]

exp_Tumor[1:5, 1:5]

gene <- rownames(exp_Tumor)

p_matrix <- matrix(nrow = length(gene), ncol = 3, dimnames = list(gene, c("p_TLS+_vs_TLS-", "mean_TLS+", 
    "mean_TLS-")))

df_gene <- df_HE

for (g in gene) {
    df_gene$gene <- as.numeric(exp_Tumor[g, HE_sample])
    res <- wilcox.test(gene ~ TLS, data = df_gene)
    p_matrix[g, "p_TLS+_vs_TLS-"] <- res$p.value
    df <- df_gene %>% group_by(TLS) %>% summarise(mean = mean(gene))
    rownames(df) <- df$TLS
    p_matrix[g, "mean_TLS+"] <- as.numeric(df["positive", "mean"])
    p_matrix[g, "mean_TLS-"] <- as.numeric(df["negative", "mean"])
}

df_TLS <- as.data.frame(p_matrix) %>% mutate(FC = `mean_TLS+`/`mean_TLS-`)

df_TLS_diff <- df_TLS %>% filter(`p_TLS+_vs_TLS-` < 0.05)

dim(df_TLS_diff)

dir.create(dirname(file.path(table_dir, "TCGA/diff_TCGA_TLS_wilcox.csv")), recursive = TRUE, showWarnings = FALSE)

write.csv(df_TLS, file = file.path(table_dir, "TCGA/diff_TCGA_TLS_wilcox.csv"))

#### Tumor/Paratumor Wilcoxon ====
exp <- read.csv(file.path(input_dir, "data/TCGA_LIHC_count_symbol.csv"), header = T, check.names = F, 
    row.names = 1)

group <- read.csv(file.path(input_dir, "data/TCGA_group.csv"), header = T, check.names = F)

df_TP <- group

TP_sample <- pull(df_TP, samplename)

pro_ref <- read.csv(file.path(input_dir, "reference_windows/pro_coding_genes.csv"), header = T, 
    row.names = 1)

pro_coding_genes <- pull(pro_ref, gene_name)

exp_Tumor <- exp[intersect(pro_coding_genes, rownames(exp)), TP_sample]

exp_Tumor[1:5, 1:5]

gene <- rownames(exp_Tumor)

p_matrix <- matrix(nrow = length(gene), ncol = 3, dimnames = list(gene, c("p_T_vs_P", "mean_T", 
    "mean_P")))

df_gene <- df_TP

for (g in gene) {
    df_gene$gene <- as.numeric(exp_Tumor[g, TP_sample])
    res <- wilcox.test(gene ~ TissueType, data = df_gene)
    p_matrix[g, "p_T_vs_P"] <- res$p.value
    df <- df_gene %>% group_by(TissueType) %>% summarise(mean = mean(gene))
    rownames(df) <- df$TissueType
    p_matrix[g, "mean_T"] <- as.numeric(df["Tumor", "mean"])
    p_matrix[g, "mean_P"] <- as.numeric(df["Paratumor", "mean"])
}

df_TP <- as.data.frame(p_matrix) %>% mutate(FC = mean_T/mean_P)

df_TP_diff <- df_TP %>% filter(p_T_vs_P < 0.05)

dim(df_TP_diff)

dir.create(dirname(file.path(table_dir, "TCGA/diff_TCGA_TP_wilcox.csv")),
    recursive = TRUE, showWarnings = FALSE)

write.csv(df_TP, file = file.path(table_dir, "TCGA/diff_TCGA_TP_wilcox.csv"))

