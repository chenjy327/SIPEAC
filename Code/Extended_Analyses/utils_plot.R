suppressPackageStartupMessages({
  library(ggplot2); library(glue); library(scales)
})

PYTHON_BIN <- "python3"
PLOT_CLUSTER_PY <- proj_path("Code", "Extended_Analyses", "utils_spatial_plot.py")

theme_csl <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text   = element_text(color = "black"),
      axis.line   = element_line(linewidth = 0.3, color = "black"),
      axis.ticks  = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.2, "pt"),
      plot.title  = element_text(face = "bold", size = base_size),
      legend.key.size    = unit(4, "pt"),
      legend.title       = element_text(size = base_size),
      legend.text        = element_text(size = base_size - 1),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold", size = base_size)
    )
}

pal_csl <- list(
  tissue = c(
    HCC            = "#C0392B",   # 肝细胞癌（主色）
    `Normal liver` = "#5DADE2",
    Cirrhosis      = "#F39C12",
    CRC            = "#27AE60",
    Hepatoblastoma = "#8E44AD"
  ),
  group = c(
    TLS_in   = "#C0392B",
    TLS_out  = "#2C3E50",
    `anti-IFI30+` = "#C0392B",
    `anti-IFI30-` = "#7F8C8D",
    HCC      = "#C0392B",
    Healthy  = "#5DADE2"
  ),
  patient = scales::hue_pal()   # 按需 palette，患者多时调用
)

mycolor_new <- c("#F7CCC1", "#BCD1EC", "#9688BF", "#A8DEE0","#F3BDD7",  "#1B97A8","#C3DDF3", "#F09AB7", "#E0E0DF", "#F9E0BD","#78B0DD","#FADBBF", "#8ACFE2",
                 "#C198C4", "#ECA579", "#E97BAC", "#C0DA82", "#92D3F1", "#CFB5D5", "#94CAED", "#639EBC", "#D47AAD", "#94C54B", "#D6BEDB",
                 "#E9866D", "#F1ADB4", "#7985C0", "#87BEE7",  "#F1AE7A", "#F2A289", "#EBCC83", "#C67FB2", "#C9D1E9",
                 "#F8C7A7", "#C6CDE8", "#E98D9B", "#F4D37A", "#6765AB", "#E585A4", "#F1A681", "#78C8D5", "#C58ADE", "#B1BEEA", "#90C4E9",
                 "#8095CE", "#FEA78C", "#FFA3A6", "#F583B4", "#CD69A7", "#ED7179", "#D6A3DC", "#F7DB70", "#EABEBF", "#75CCE8",
                 "#A5DEE5", "#80BEAF", "#B3DDD1", "#D1DCE2", "#F5B994", "#EE9C6A", "#FFF5E9", "#FFE5C0", "#EAD2DF", "#EDB5D2",
                 "#C8A8DA", "#FEB396", "#FFCCBF", "#AED4D5", "#F9CC88", "#F5CEC7", "#E79796", "#FFB284", "#FFC98B", "#C6C09C",
                 "#879E46", "#BBD5A6", "#FBCEB9", "#FEBD3D", "#E57B87", "#FFC4C8", "#FED154","#BEE2F7", "#FEB25E", "#6AC7E6", "#FB93AC",
                 "#9BC768", "#FDC7D4", "#FDD84C", "#BB9DCF", "#E49AAB", "#FFDFA2", "#BFE4FF", "#A3B5FD", "#C3C5F8",
                 "#60EFDB", "#BEF2E5", "#C5E7F1", "#79CEED", "#6F89A2", "#85CBCD",  "#F9E2AE", "#FBC78D", "#A7D676",
                 "#85CBCD", "#A8DEE0", "#F9E2AE", "#FBC78D", "#D2A8B2", "#EDE1E5", "#DAE6E4", "#EECF65", "#DA9C25")

mycolor_cell <- c("#238b45","#ce1256","#ce1256","#1f78b4","#6a3d9a","#fdbf6f","#df65b0",
                  "#eaeaea","#C9BEE8","#eaeaea","#C9BEE8",
                  "#6F89A2", "#eaeaea", "#FDD84C", "#A1CFFA", "#F7CCC1", "white")
mycell <- c('B', 'T/NK','T',  'Endothelial', 'Myeloid', 'Plasma', 'Fibroblast',
            'Hepatocyte', "Tumor", "Tumor&Hepatocyte","Hep/Tumor",
            "Mast", "Doublet", "Epithelial",'pDC', 'Neutrophil', "None")
names(mycolor_cell) <- mycell

cursubtype1 <- c("NaiveB", "MemoryB", "GCB", 'Plasma','B_Other',
                 'Tnaive', 'Tmemory', 'Tpex', 'Tex', 'Treg', 'Tfh','T_Other',
                 'Mph_CXCL10', 'Mph_FOLR2', 'Mono-like_FCGR3A', 'Mph_TREM2',
                 'cDC1_CLEC9A', 'cDC2_CD1C', 'cDC3_LAMP3',  'pDC','myeloid_others',
                 "vCAF","iCAF","myCAF", "CAF_others" ,
                 "Tumor&Hepatocyte","Hepatocyte",
                 "Vein" ,"LSEC" , "Lymphatic" ,  "HEV",
                 "PVC"  ,"Tip", "Stalk", "Capillary"  ,
                 "Endo_others", "Artery", "HLA_EC",
                 "None",
                 'T/NK','T',  'Endothelial', 'Myeloid',  'Fibroblast',
                 "Mast", "Doublet", "Epithelial",'pDC', 'Neutrophil')
mycolorcellsub <- c('#78c679', '#a1d99b', '#008000', '#eb5c2f', '#66c2a4',
                    '#FF0000', '#e45f56', '#FF4C3B', '#cd5554', '#d40e52', '#82052f', '#fcbba1',
                    '#6a51a3',  "#edf8fb", "#b3cde3", "#8c96c6",
                    "#810f7c",'#dd3497', '#996699', '#756bb1', '#c994c7',
                    "#f1eef6", "#d7b5d8", "#df65b0", "#ce1256",
                    "#eaeaea","#eaeaea",
                    "#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4",
                    "#00A0FF", "#2BDE73", "#FFCC00", "#9178Fe",
                    "#1d91c0", "#225ea8", "#0c2c84",
                    "#eaeaea",
                    "#ce1256","#ce1256","#1f78b4","#6a3d9a","#df65b0",
                    "#6F89A2", "#eaeaea", "#FDD84C", "#A1CFFA", "#F7CCC1")
names(mycolorcellsub) <- cursubtype1

tls_color <- c('#db3196', '#f57622', '#4cb1d2', '#a6bddb', '#bdbdbd', "#999999")
names(tls_color) <- c("Mature", "Developing", "NotMature", 'Rare_Naive_GCB', 'Filtered', "No_TLS")

tls_color_2 <- c('#db3196', '#f57622', '#4cb1d2', '#a6bddb', '#bdbdbd', "#999999")
names(tls_color_2) <- c("Mature", "deviating","conforming", 'Rare_Naive_GCB', 'Filtered', "Out_TLS")

mycolor16 <- c("#BFE2E3","#A1CFFA","#78BDAD","#D45651","#397A7F","#F0918E","#EEE8DA","#1F5392","#A0BFAF",
               "#AE98D6","#ECCBDC","#54BAD3","#8b4a4b","#DB896C","#AABAC2","#ffae3b")
mycolor22 <- c('#CCCCCC','#B5B5B5','#A092C3','#DDA0DD','#03A4C6','#7AC5CD','#A3D9ED','#00E5EE','#F9EBDA','#F5DEB3',
               '#98689E','#E84115','#FFB5C5','#00FF7F','#F9D01C','#B03060','#00ABDC','#D2691E','#03A464','#FF7F00',
               '#8968CD','#1C5B75')
mycolor21 <- c("#EE3743", "#33a02c", "#0B7FAB", "#6a3d9a", "#ff5084", "#ff7f00", "#01ceff", "#F6A9BD", "#fdbf6f", "#c51b7d",
               "#6F6C9E", "#01665e", "#a6cee3", "#b2df8a", "#ffff99", "#bf812d",  "#999999", "#BB7DB2", "#824615", "#ffff33", "#86DBD4")
mycolor_merge <- c(mycolor21, mycolor16, mycolor22)
color_paper <- c("#00A0FF", "#2BDE73", "#FFCC00", "#9178Fe", "#08D9D6", "#fbf579", "#eb94a2", "#FC4A1A", "#009A31","#2060ff", "#EEA282", "#db50B2")

colallmain <- c("#78BDAD", "#AE98D6","#8b4a4b","#DB896C","#ffae3b")
names(colallmain) <- c('IGHG','IGHA', 'IGHM', 'IGHD', 'IGHE')
isotype_sub <- c('IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHA1', 'IGHA2', 'IGHM', 'IGHD', 'IGHE')
colallsub <- c('#FF0000', '#e45f56', '#FF4C3B', '#cd5554',"#810f7c",'#dd3497',"#225ea8", "#0c2c84", '#008000')
names(colallsub) <- isotype_sub
colorsecret <- c("#5FBFFC","#FFAF3D")
names(colorsecret) <- c("Secret", "Membrane")

col_clone_expanded_5 <- c("#FF412A","#FFAF3D","#B8FFEA","#5FBFFC","#0042A4")
col_clone_expanded <- c('#ffae3b', '#1F5392')
names(col_clone_expanded) <- c("Expanded","Singleton")

save_plot <- function(gg, prefix, width, height, if_png = FALSE, if_pdf = FALSE) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  if (if_pdf) {
    ggsave(glue("{prefix}.pdf"), gg, width = width, height = height, bg = "white")
  }
  if (if_png) {
    ggsave(glue("{prefix}.png"), gg, width = width, height = height,
           dpi = 600, bg = "white", limitsize = FALSE)
  }
}

SpatialColors <- function(n) {
  x <- ramp(seq.int(0, 1, length.out = n))
  if (ncol(x) == 4L)
    rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
  else
    rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
}

spatial_dim_plot <- function(obj, target, colors, prefix, point.size = 0.5) {
  p <- ggplot(data = obj@meta.data, aes_string(colour = target)) +
    geom_point(aes(x = max(col) - col, y = max(row) - row),
               size = point.size, show.legend = FALSE) +
    scale_color_manual(values = colors) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(), panel.background = element_blank()
    ) +
    coord_fixed() + theme(legend.position = "bottom")
  save_plot(p, glue("{prefix}"), width = 8, height = 8, if_png = TRUE)
}

spatial_feature_plot <- function(obj, target, colors, prefix, point.size = 0.5, width = 8, height = 8) {
  p <- ggplot(data = obj@meta.data, aes_string(color = target)) +
    geom_point(aes(x = max(col) - col, y = max(row) - row),
               size = point.size, show.legend = TRUE) +
    scale_color_gradientn(colors = colors) +
    scale_alpha(range = c(1, 1)) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(), panel.background = element_blank(),
      plot.title = element_text(size = 70, color = "black", face = "bold", hjust = 0.5, vjust = -3)
    ) +
    coord_fixed() + theme(legend.position = "bottom")
  save_plot(p, glue("{prefix}"), width = width, height = height, if_png = TRUE)
}

SpatialFeaturePython <- function(obj, if_gene = FALSE, plot_item, prefix,
                                 if_tls_contour = "F", name = "",
                                 vmin = "", vmax = "", if_count = FALSE) {
  if (!dir.exists(prefix)) dir.create(prefix, recursive = TRUE)

  tmp_meta <- obj@meta.data
  row_n <- rownames(obj@meta.data)
  if (if_gene) {
    if (if_count) {
      tmp_meta[[plot_item]] <- obj@assays$Spatial@counts[plot_item, ]
    } else {
      tmp_meta[[plot_item]] <- obj@assays$Spatial@data[plot_item, ]
    }
  }
  if (if_tls_contour == "T") {
    tmp_meta <- tmp_meta[, c("row", "col", "TLS", plot_item)]
  } else {
    tmp_meta <- tmp_meta[, c("row", "col", plot_item)]
  }
  if (vmin == "" & vmax == "") {
    vmin <- min(tmp_meta[[plot_item]])
    vmax <- max(tmp_meta[[plot_item]])
  }
  name <- glue("{name}")
  rownames(tmp_meta) <- row_n
  write.table(obj@images$slice1@coordinates, glue('{prefix}/df_pos_tmp.txt'), quote = FALSE, sep = '\t')
  write.table(tmp_meta, glue('{prefix}/df_meta_tmp.txt'), quote = FALSE, sep = '\t')

  command <- glue("{PYTHON_BIN} {PLOT_CLUSTER_PY} -i {plot_item} --plot_type FeaturePlot -o {prefix} -n {name} --plot_tls_contour {if_tls_contour} --vmin {vmin} --vmax {vmax}")
  system(command)
  command <- glue("{PYTHON_BIN} {PLOT_CLUSTER_PY} -i {plot_item} --plot_type Colorbar -o {prefix} -n {name}")
  system(command)
  system(glue('rm {prefix}/df_pos_tmp.txt {prefix}/df_meta_tmp.txt'))
}

SpatialDimPython <- function(obj, plot_item, prefix, tmp_color, name,
                             plot_order = NULL, if_tls_contour = "F") {
  if (!dir.exists(prefix)) dir.create(prefix, recursive = TRUE)

  if (is.null(plot_order)) {
    plot_order <- unique(obj@meta.data[[plot_item]]) %>% as.character()
  }

  row_n <- rownames(obj@meta.data)
  tmp_meta <- obj@meta.data
  rownames(tmp_meta) <- tmp_meta$BinID
  if (if_tls_contour == "T") {
    tmp_meta <- tmp_meta[, c("row", "col", "TLS", plot_item)]
  } else {
    tmp_meta <- tmp_meta[, c("row", "col", plot_item)]
  }
  rownames(tmp_meta) <- row_n

  tmp_color <- str_c(tmp_color, collapse = ' ')
  pos_tmp <- obj@images$slice1@coordinates
  write.table(pos_tmp, glue('{prefix}/df_pos_tmp.txt'), quote = FALSE, sep = '\t')
  write.table(tmp_meta, glue('{prefix}/df_meta_tmp.txt'), quote = FALSE, sep = '\t')

  tmp <- str_c(plot_order, collapse = ' ')
  command <- glue("{PYTHON_BIN} {PLOT_CLUSTER_PY} -i {plot_item} --plot_type DimPlot -o {prefix} --plot_order '{tmp}' -n {name} --plot_color '{tmp_color}' --plot_tls_contour {if_tls_contour}")
  system(command)
  system(glue('rm {prefix}/df_pos_tmp.txt {prefix}/df_meta_tmp.txt'))
}

single_cell_chisq_heatmap <- function(single.cell.obj, column.to.chisq, compare.groups,
                                      text_use, file_name,
                                      heatmap_palette = colorRampPalette(c("#FFFFFF", '#FFCC66', "#CC3300"))(100),
                                      heatmap_path = "chisq_heatmap/",
                                      cellwidth = 60, cellheight = 32, fontsize = 24, width = 20, height = 12) {
  require(Seurat); require(tidyverse); require(circlize)
  if (!endsWith(heatmap_path, "/")) heatmap_path <- paste0(heatmap_path, "/")
  if (!dir.exists(heatmap_path)) dir.create(heatmap_path, recursive = TRUE)

  clusters <- sort(unique(single.cell.obj[, column.to.chisq]))
  meta.data <- single.cell.obj

  split_df <- split(meta.data, list(meta.data[, compare.groups], meta.data[, column.to.chisq]))
  counts <- sapply(split_df, nrow)
  a <- data.frame(do.call(rbind, strsplit(names(counts), "\\.")), a = counts, row.names = NULL)
  names(a) <- c(compare.groups, column.to.chisq, "a")

  split_df <- split(meta.data, list(meta.data[, compare.groups]))
  counts <- sapply(split_df, nrow)
  ac <- data.frame(do.call(rbind, strsplit(names(counts), "\\.")), ac = counts, row.names = NULL)
  names(ac) <- c(compare.groups, "ac")

  split_df <- split(meta.data, list(meta.data[, column.to.chisq]))
  counts <- sapply(split_df, nrow)
  ab <- data.frame(do.call(rbind, strsplit(names(counts), "\\.")), ab = counts, row.names = NULL)
  names(ab) <- c(column.to.chisq, "ab")

  d1 <- merge(a, ac); d1 <- merge(d1, ab)
  d1["b"] <- d1$ab - d1$a
  d1["c"] <- d1$ac - d1$a
  d1["d"] <- sum(d1$a) - d1$ab - d1$c

  d1[, c("p.value", "Ea", "Ec", "Eb", "Ed")] <- t(apply(d1, 1, function(x) {
    x <- as.numeric(x[c('a', 'c', 'b', 'd')])
    k.test <- chisq.test(matrix(x, nrow = 2, ncol = 2))
    return(c(k.test$p.value, as.vector(k.test$expected)))
  }))
  d1["Reo"] <- d1$a / d1$Ea
  d1$p.value.symbol <- cut(d1$p.value, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                           labels = c('****', '***', '**', '*', 'N.S.'))

  require(reshape2); require(pheatmap); require(ComplexHeatmap)
  dist.heatmap.data <- d1[, c(column.to.chisq, compare.groups, "Reo")]
  dist.heatmap <- reshape2::dcast(dist.heatmap.data, dist.heatmap.data[, 1] ~ dist.heatmap.data[, 2], value.var = "Reo")
  names(dist.heatmap)[1] <- column.to.chisq
  rownames(dist.heatmap) <- dist.heatmap[, column.to.chisq]
  dist.heatmap <- dist.heatmap[, -1]
  dist.heatmap1 <- round(dist.heatmap, digits = 2)
  dist.heatmap1 <- dist.heatmap1[order(factor(rownames(dist.heatmap1), levels = clusters)), ]

  dist.heatmap.pval.symbols <- d1[, c(column.to.chisq, compare.groups, "p.value.symbol")]
  dist.heatmap.pval.symbols <- reshape2::dcast(dist.heatmap.pval.symbols, dist.heatmap.pval.symbols[, 1] ~ dist.heatmap.pval.symbols[, 2], value.var = "p.value.symbol")
  names(dist.heatmap.pval.symbols)[1] <- column.to.chisq
  rownames(dist.heatmap.pval.symbols) <- dist.heatmap.pval.symbols[, column.to.chisq]
  dist.heatmap.pval.symbols <- dist.heatmap.pval.symbols[, -1]
  dist.heatmap.pval.symbols <- dist.heatmap.pval.symbols[order(factor(rownames(dist.heatmap.pval.symbols), levels = clusters)), ]

  label_df <- as.data.frame(matrix(ncol = ncol(dist.heatmap1), nrow = nrow(dist.heatmap1)))
  for (o in 1:nrow(dist.heatmap1)) {
    for (p in 1:ncol(dist.heatmap1)) {
      label_df[o, p] <- paste(dist.heatmap1[o, p], dist.heatmap.pval.symbols[o, p], sep = "\n")
    }
  }

  palette <- colorRamp2(c(0, max(dist.heatmap1) / 2, max(dist.heatmap1)), c("#FFFFFF", '#FFCC66', "#CC3300"))
  ht <- Heatmap(dist.heatmap1,
               name = "Ro/e", col = palette,
               show_row_names = TRUE, show_column_names = TRUE,
               row_title = column.to.chisq,
               column_title = paste0("Chi-square test for ", text_use, " in each cluster"),
               cluster_columns = FALSE, cluster_rows = FALSE,
               row_title_gp = gpar(fontsize = 0),
               row_names_gp = gpar(fontsize = fontsize),
               column_names_gp = gpar(fontsize = fontsize),
               column_title_gp = gpar(fontsize = fontsize),
               column_names_rot = 45,
               width = unit(cellwidth / 1.5, "mm") * ncol(dist.heatmap1),
               height = unit(cellheight / 2.5, "mm") * nrow(dist.heatmap1),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(label_df[i, j], x, y, gp = gpar(fontsize = fontsize))
               },
               heatmap_legend_param = list(
                 direction = "vertical",
                 title_gp = gpar(fontsize = fontsize),
                 labels_gp = gpar(fontsize = fontsize - 4),
                 legend_height = unit(4, "cm"),
                 grid_width = unit(.8, "cm")
               ))
  pdf(file = paste0(heatmap_path, file_name, "_chisq_test_heatmap.pdf"), width = width, height = height)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}

My_AddModuleScore <- function(data_obj, geneSet, geneSet_name, mycolor, crop, pt.size, cur_dir) {
  data_obj <- AddModuleScore(object = data_obj, features = geneSet, ctrl = 100,
                             name = glue('{geneSet_name}'))
  p1 <- VlnPlot(data_obj, features = glue('{geneSet_name}'), pt.size = 0) +
    scale_fill_manual(values = mycolor) + labs(title = glue('{geneSet_name}'))
  save_plot(p1, glue("{cur_dir}/{geneSet_name}1"), width = 8, height = 8, if_png = TRUE)

  names(data_obj@meta.data)[names(data_obj@meta.data) == glue('{geneSet_name}')] <- glue('{geneSet_name}(AddModuleScore)')
  p2 <- SpatialPlot(data_obj, image.alpha = 0, stroke = 0, crop = crop,
                    features = glue('{geneSet_name}(AddModuleScore)'), pt.size.factor = pt.size) +
    theme(legend.position = "bottom")
  save_plot(p2, glue("{cur_dir}/{geneSet_name}(AddModuleScore)"), width = 8, height = 8, if_png = TRUE)

  data_obj <- calculate_countSum(data_obj, geneSet)
  names(data_obj@meta.data)[names(data_obj@meta.data) == '{geneSet_name}_CountSum'] <- glue('{geneSet_name}(CountSum)')
  p3 <- SpatialPlot(data_obj, image.alpha = 0, stroke = 0, crop = crop,
                    features = glue('{geneSet_name}(CountSum)'), pt.size.factor = pt.size) +
    theme(legend.position = "bottom")
  save_plot(p3, glue("{cur_dir}/{geneSet_name}(CountSum)"), width = 8, height = 8, if_png = TRUE)
  return(data_obj)
}

save_fig <- function(plot, name, subdir = NULL, width = 90, height = 90) {
  out <- dirs$figures
  if (!is.null(subdir)) out <- file.path(out, subdir)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  prefix <- file.path(out, name)
  save_plot(plot, prefix, width = width / 25.4, height = height / 25.4, if_png = TRUE, if_pdf = TRUE)
}

harmonize_st_meta <- function(obj) {
  m <- obj@meta.data
  cursampleid <- obj$SampleID %>% unique()
  m$Bin_Region[m$Bin_Region == "Tumor_capsule"] <- "Invasive_zone"
  m$Bin_Region[m$Bin_Region %in% c("Tumor_side_of_Margin_area", "Tumor_tissue")] <- "Tumor"
  m$Bin_Region[m$Bin_Region %in% c("Paratumor_side_of_Margin_area", "Paratumor_tissue")] <- "Paratumor"
  m$TLS <- m$TLS_raw
  m$TLS_maturity[m$TLS_maturity == "Conforming"]     <- "conforming"
  m$TLS_maturity[m$TLS_maturity == "Deviating"]      <- "deviating"
  m$TLS_maturity[m$TLS_maturity == "NotMature"]      <- "deviating"
  m$TLS_maturity[m$TLS_maturity == "Rare_Naive_GCB"] <- glue("ST{cursampleid}_NA")
  m[!m$TLS_maturity %in% c("Mature", "deviating", "conforming"), ]$TLS <- glue("ST{cursampleid}_NA")
  obj@meta.data <- m
  obj
}

annotate_clone <- function(obj, IGH_bin, IGL_bin, pair1, IGHL_bin = intersect(IGH_bin, IGL_bin)) {
  cursampleid <- obj$SampleID %>% unique()
  obj@meta.data$mab_name <- NULL
  obj@meta.data$mab_name[obj@meta.data$BinID %in% IGL_bin]  <- pair1$IGL
  obj@meta.data$mab_name[obj@meta.data$BinID %in% IGH_bin]  <- pair1$IGH
  obj@meta.data$mab_name[obj@meta.data$BinID %in% IGHL_bin] <- paste(pair1$IGH, pair1$IGL, sep = "_")
  obj@meta.data$bin_id <- rownames(obj@meta.data)
  obj@meta.data$mab_name <- as.character(obj@meta.data$mab_name)
  obj@meta.data$mab_name <- ifelse(is.na(obj@meta.data$mab_name), "other", obj@meta.data$mab_name)
  obj@meta.data$mab_name <- ifelse(obj@meta.data$TLS != glue("ST{cursampleid}_NA") & obj@meta.data$mab_name == "other",
                                   "TLS", obj@meta.data$mab_name)
  obj@meta.data$mab_name <- ifelse(obj@meta.data$mab_name == "other",
                                   obj@meta.data$Bin_Region, obj@meta.data$mab_name)
  obj$TLS <- obj$Bin_Region
  obj$TLS[obj$TLS != "Invasive_zone"] <- "_NA"
  obj
}

expand_bin_ring <- function(bin_set, meta_rc) {
  rc <- meta_rc %>% dplyr::filter(BinID %in% bin_set) %>%
    dplyr::transmute(row = as.integer(row), col = as.integer(col), clone = "c")
  if (nrow(rc) == 0) return(bin_set)
  rc_exp <- st_expand_n(rc, n = 1, x = "col", y = "row",
                        expand_n = "expand_n", group.by = "clone", mode = 2)
  rc_all <- bind_rows(
    rc     %>% dplyr::transmute(row, col),
    rc_exp %>% dplyr::transmute(row = as.integer(row), col = as.integer(col))
  ) %>% dplyr::distinct()
  meta_rc %>%
    dplyr::mutate(row_i = as.integer(row), col_i = as.integer(col)) %>%
    dplyr::inner_join(rc_all, by = c("row_i" = "row", "col_i" = "col")) %>%
    dplyr::pull(BinID) %>% unique()
}

load_st_obj <- function(cur_sample, st_root = dirs$st_root, st_meta = dirs$st_meta) {
  obj_path <- file.path(st_root, glue("{cur_sample}_scRefST2934T_random100_marker20_spotlight_Spatial.rds"))
  meta_path <- file.path(st_meta, glue("{cur_sample}.txt"))
  if (!file.exists(obj_path) || !file.exists(meta_path)) return(NULL)
  obj <- readRDS(obj_path)
  obj@meta.data <- read.table(meta_path, header = TRUE) %>% as.data.frame()
  obj@meta.data$BinID <- substr(obj@meta.data$BinID, 7, nchar(obj@meta.data$BinID))
  rownames(obj@meta.data) <- obj@meta.data$BinID
  obj@meta.data$bin_id <- obj@meta.data[, 1]
  harmonize_st_meta(obj)
}
