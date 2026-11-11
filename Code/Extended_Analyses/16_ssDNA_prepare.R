#### ssDNA h5ad转Seurat和坐标导出 ====
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(PROJECT_ROOT, "data")
source("Code/Extended_Analyses/R/io.R")
source("Code/Extended_Analyses/R/plot.R")
input_dir <- file.path(data_dir, "hzy")
processed_dir <- file.path(dirs$processed, "hzy")
table_dir <- file.path(dirs$tables, "hzy", "ssDNA")
figure_dir <- file.path(dirs$figures, "hzy")
library(glue)

library(reticulate)

library(Seurat)

use_condaenv("hotspot", required = TRUE)

py_config()

library(anndata)

library(Matrix)

adata <- read_h5ad(file.path(input_dir, "ssDNA/2879T.h5ad"))

counts_mat <- Matrix(adata$X, sparse = TRUE)

print(class(counts_mat))

rownames(counts_mat) <- adata$obs_names

colnames(counts_mat) <- adata$var_names

counts_mat <- Matrix::t(counts_mat)

obj <- CreateSeuratObject(counts = counts_mat, meta.data = as.data.frame(adata$obs))

dir.create(dirname(file.path(processed_dir, "ssDNA/2879T.rds")), recursive = TRUE, showWarnings = FALSE)

saveRDS(obj, file.path(processed_dir, "ssDNA/2879T.rds"))

cur_sample = "2879T"

seurat_obj <- readRDS(glue(file.path(processed_dir, "ssDNA/{cur_sample}.rds")))

cell_meta <- data.frame(cell_id = rownames(seurat_obj@meta.data), x = seurat_obj@meta.data$x, y = seurat_obj@meta.data$y, 
    stringsAsFactors = FALSE)

dir.create(dirname(glue(file.path(table_dir, "cell_meta.csv"))), recursive = TRUE, showWarnings = FALSE)

write.csv(cell_meta, glue(file.path(table_dir, "cell_meta.csv")), row.names = FALSE)

tb <- readRDS(file.path(input_dir, "ssDNA/merged_bcr.rds"))

dir.create(dirname(file.path(table_dir, "merged_bcr.csv")), recursive = TRUE, showWarnings = FALSE)

write.csv(tb, file.path(table_dir, "merged_bcr.csv"), row.names = FALSE)

