library(argparse)
library(data.table)
library(glue)
library(dplyr)
library(dplyr)
library(data.table)
library(glue)
library(stringr)
library(shazam)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(rstatix)
library(ggsci)
library(ggsignif)
# 创建命令行参数解析器
parser <- ArgumentParser()
parser$add_argument('-s', '--sample', type='character', default=64, help='prefix')
parser$add_argument('-o', '--outdir', type='character', default=64, help='prefix')

# 解析命令行参数
args <- parser$parse_args()
sample <- args$sample
outdir <- args$outdir

# 读取数据框
df_changeo <- fread(glue("{outdir}/{sample}_clone-pass.tsv"))
df_changeo$clone_family_07 <- df_changeo$clone_id
write.table(df_changeo, 
            file=glue("{outdir}/{sample}_family.tsv"),
            sep='\t', row.names=FALSE, quote = F, col.names = T)

