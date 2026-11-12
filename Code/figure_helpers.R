# Figure1.R、Figure E1.R、Figure E2.R、Figure2,E3.R、Figure3.R 共用的函数和颜色
# 摘自原来的 各种包和颜色.R，只保留这几份脚本用到的部分
# 工作目录设为 SIPEAC 根目录后 source

suppressPackageStartupMessages({
  library(glue)
  library(stringr)
  library(dplyr)
  library(tibble)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(ggforce)
  library(scales)
  library(Seurat)
  library(alakazam)
  library(scuttle)
  library(ComplexHeatmap)
})

source("Code/Extended_Analyses/utils_spatial_expand.R")

# python 路径可用环境变量 SIPEAC_PYTHON 指定
PYTHON_BIN <- Sys.getenv("SIPEAC_PYTHON", "python3")
PLOT_CLUSTER_PY <- "Code/Extended_Analyses/utils_spatial_plot.py"

mycolor_merge <- c("#EE3743", "#33a02c", "#0B7FAB", "#6a3d9a", "#ff5084", "#ff7f00", "#01ceff", "#F6A9BD", "#fdbf6f", "#c51b7d",
                   "#6F6C9E", "#01665e", "#a6cee3", "#b2df8a", "#ffff99", "#bf812d",  "#999999", "#BB7DB2", "#824615", "#ffff33", "#86DBD4",
                   "#BFE2E3","#A1CFFA","#78BDAD","#D45651","#397A7F","#F0918E","#EEE8DA","#1F5392","#A0BFAF",
                   "#AE98D6","#ECCBDC","#54BAD3","#8b4a4b","#DB896C","#AABAC2","#ffae3b",
                   '#CCCCCC','#B5B5B5','#A092C3','#DDA0DD','#03A4C6','#7AC5CD','#A3D9ED','#00E5EE','#F9EBDA','#F5DEB3',
                   '#98689E','#E84115','#FFB5C5','#00FF7F','#F9D01C','#B03060','#00ABDC','#D2691E','#03A464','#FF7F00',
                   '#8968CD','#1C5B75')

col_clone_expanded <- c('#ffae3b', '#1F5392')
names(col_clone_expanded) <- c("Expanded", "Singleton")


SpatialFeaturePython <- function(obj,
                                 if_gene = F,
                                 plot_item, # plot的基因
                                 prefix, # 当前路径
                                 if_tls_contour = "F",
                                 name = "",
                                 vmin = "",
                                 vmax = "",
                                 if_count = F
){
  if(!dir.exists(prefix)){
    dir.create(prefix, recursive = T)
  }

  # 创建临时元数据tmp_meta，包括"row"、"col"和plot_item列
  tmp_meta <- obj@meta.data
  row_n <- rownames(obj@meta.data)
  if(if_gene){
    if(if_count){
      tmp_meta[[plot_item]] <- obj@assays$Spatial@counts[plot_item, ]
    }
    else{
      tmp_meta[[plot_item]] <- obj@assays$Spatial@data[plot_item, ]
    }
  }
  if(if_tls_contour == "T"){
    tmp_meta <- tmp_meta[, c("row", "col", "TLS", plot_item)]
  }else{
    tmp_meta <- tmp_meta[, c("row", "col", plot_item)]
  }
  if(vmin == "" & vmax == ""){
    vmin <- min(tmp_meta[[plot_item]])
    vmax <- max(tmp_meta[[plot_item]])
  }
  name <- glue("{name}")
  rownames(tmp_meta) <- row_n
  # 将hcc_st@images$slice1@coordinates写入到df_pos_tmp.txt文件中
  write.table(obj@images$slice1@coordinates, glue('{prefix}/df_pos_tmp.txt'), quote=F, sep='\t')
  # 将tmp_meta写入到df_meta_tmp.txt文件中
  write.table(tmp_meta, glue('{prefix}/df_meta_tmp.txt'), quote=F, sep='\t')

  command <- glue("{PYTHON_BIN} {PLOT_CLUSTER_PY} -i {plot_item} --plot_type FeaturePlot -o {prefix} -n {name} --plot_tls_contour {if_tls_contour} --vmin {vmin} --vmax {vmax}")
  system(command)

  # colorbar
  command <- glue("{PYTHON_BIN} {PLOT_CLUSTER_PY} -i {plot_item} --plot_type Colorbar -o {prefix} -n {name} --vmin {vmin} --vmax {vmax}")
  system(command)

  # 删除临时文件
  system(glue('rm {prefix}/df_pos_tmp.txt {prefix}/df_meta_tmp.txt'))
}

SpatialDimPython <- function(obj,
                             plot_item, # plot的基因
                             prefix, # 当前路径
                             tmp_color,
                             name,
                             plot_order = NULL,
                             if_tls_contour = "F"
){
  if(!dir.exists(prefix)){
    dir.create(prefix, recursive = T)
  }

  # 设置绘图顺序为当前单元格值的列表
  if(is.null(plot_order)){
    plot_order <- unique(obj@meta.data[[plot_item]]) %>% as.character()
  }

  # 创建临时元数据tmp_meta，包括"row","col"和plot_item列
  row_n <- rownames(obj@meta.data)
  tmp_meta <- obj@meta.data
  rownames(tmp_meta) <- tmp_meta$BinID
  if(if_tls_contour == "T"){
    tmp_meta <- tmp_meta[, c("row", "col", "TLS", plot_item)]
  }else{
    tmp_meta <- tmp_meta[, c("row", "col", plot_item)]
  }
  rownames(tmp_meta) <- row_n

  tmp_color <- str_c(tmp_color, collapse=' ')
  # 将hcc_st@images$slice1@coordinates写入到df_pos_tmp.txt文件中
  pos_tmp <- obj@images$slice1@coordinates
  write.table(pos_tmp, glue('{prefix}/df_pos_tmp.txt'), quote=F, sep='\t')

  # 将tmp_meta写入到df_meta_tmp.txt文件中
  write.table(tmp_meta, glue('{prefix}/df_meta_tmp.txt'), quote=F, sep='\t')

  tmp <- str_c(plot_order, collapse=' ')
  command <- glue("{PYTHON_BIN} {PLOT_CLUSTER_PY} -i {plot_item} --plot_type DimPlot -o {prefix} --plot_order '{tmp}' -n {name} --plot_color '{tmp_color}' --plot_tls_contour {if_tls_contour}")
  system(command)
  # 删除临时文件
  system(glue('rm {prefix}/df_pos_tmp.txt {prefix}/df_meta_tmp.txt'))
}

ST_nor_sca <- function(data_obj) {
  #先把seurat对象转回为Spatial,不然用SCT很容易出问题
  if (data_obj@active.assay == 'SCT') {
    print('Convert the default assay of seurat object to \'Spatial\'')
    DefaultAssay(data_obj) <- 'Spatial'
  }else{
    print('The current default assay of seurat object is \'Spatial\'')
  }

  #去掉在所有spot中表达量为0的基因
  genes_use <- row.names(data_obj)[rowSums(data_obj[['Spatial']]@counts)!=0]

  #去掉线粒体核糖体等基因
  genes_use <- genes_use[which(!grepl('^MT-', genes_use))]
  genes_use <- genes_use[which(!grepl('^RPS', genes_use))]
  genes_use <- genes_use[which(!grepl('^RPL', genes_use))]
  genes_use <- genes_use[which(!grepl('^GM[0-9]+$', genes_use))]
  genes_use <- genes_use[which(!grepl('^MRPS', genes_use))]

  data_obj <- subset(data_obj, features = genes_use)

  #去掉所有基因表达量都为0的spot
  spots_use <- colnames(data_obj)[colSums(data_obj[['Spatial']]@counts)!=0]

  data_obj <- subset(data_obj, cells = spots_use)

  #subset后重新Normalization，不做Scale
  data_obj <- NormalizeData(data_obj, verbose = FALSE)

  return(data_obj)
}

save_table <- function(df, prefix, row.names = F, if_csv = F, col.names = T){
  if(if_csv){
    write.csv(df,
              glue("{prefix}.csv"),
              row.names = row.names,
              quote = F,
              fileEncoding = "UTF-8")
  }else{
    write.table(df,
                glue("{prefix}.txt"),
                col.names = col.names,
                row.names = row.names,
                quote = F,
                sep = "\t",
                fileEncoding = "UTF-8")
  }
}

save_plot <- function(gg, prefix, width, height, if_png = F, if_pdf = F){
  if(if_pdf){
    ggsave(glue("{prefix}.pdf"),
           gg, width = width, height = height)
  }
  if(if_png){
    ggsave(glue("{prefix}.png"),
           gg, width = width, height = height, limitsize = FALSE)
  }
}

get_shared_info <- function(df, source  = "clone_id",target_list = c()){
  df[[source]] <- as.character(df[[source]])
  for(cur_target in target_list){
    shared_type_info <- c()
    shared_type_info2 <- c()
    for(cur.id in unique(df[[source]])){
      sample_list <- unique(df[df[[source]] == cur.id,][[cur_target]]) %>% sort()
      if(length(sample_list) > 1){
        sl <- paste(sample_list,collapse = ", ")
        shared_type <- paste0("Shared_", sl)
        st2 <- "Shared"
      }else{
        shared_type <- paste0("Private_",sample_list)
        st2 <- "Private"
      }
      shared_type_info <- c(shared_type_info, shared_type)
      shared_type_info2 <- c(shared_type_info2,st2)
    }
    names(shared_type_info) <- unique(df[[source]])
    df[[glue("{cur_target}_shared_sub")]] <- shared_type_info[df[[source]]]
    names(shared_type_info2) <- unique(df[[source]])
    df[[glue("{cur_target}_shared_main")]] <- shared_type_info2[df[[source]]]
  }
  clone_num <- table(df[[source]])
  df$clone_num <- clone_num[df[[source]]]
  return(df)
}
