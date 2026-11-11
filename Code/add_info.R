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
parser$add_argument('-b', '--b_outdir', type='character', default=64, help='prefix')  # run.changeo.sh 会传入，脚本中未使用

# 解析命令行参数
args <- parser$parse_args()
sample <- args$sample
outdir <- args$outdir
# b_outdir <- args$b_outdir

# 读取数据框
df_changeo <- fread(glue("{outdir}/{sample}_igblast_db-pass_germ-pass_parse-select.tsv"))
# df_isotype <- fread(glue("{b_outdir}/{sample}_isotype.info.txt"))
# df_isotype$sequence_id = df_isotype$molecule

# 提取所需列并进行处理
# df_changeo$molecule <- sapply(df_changeo$sequence_id, function(x)strsplit(x, " ")[[1]][1])
df_changeo <- df_changeo[!is.na(df_changeo$germline_alignment), ]
# df_changeo <- dplyr::left_join(df_changeo, df_isotype)
# df_changeo <- df_changeo[!is.na(df_changeo$isotype),]
# df_changeo <- dplyr::left_join(df_changeo, df_barcode[, c("id","binid")], by = c("molecule" = "id"))
# df_changeo <- df_changeo[!is.na(df_changeo$binid),]

# df_changeo <- df_changeo[df_changeo$binid %in% rownames(meta@meta.data),]

write.table(df_changeo,
            file=glue("{outdir}/{sample}.tsv"),
            sep='\t', row.names=FALSE, quote = F, col.names = T)



changeo_germ_pass <- df_changeo
mut_info <- changeo_germ_pass
## calculate by shazam
mutation_info <- observedMutations(
  mut_info,
  sequenceColumn = "sequence_alignment",
  germlineColumn = "germline_alignment_d_mask",
  regionDefinition=IMGT_V,
  mutationDefinition = CHARGE_MUTATIONS,
  frequency = T,
  combine = T,
  nproc = 40
)

mutation_info <- mutation_info %>% dplyr::rename(shm = mu_freq)
write.table(mutation_info, 
            file=glue("{outdir}/{sample}.tsv"),
            sep='\t', row.names=FALSE, quote = F, col.names = T)

shm_mut_rate_density <- ggplot(mutation_info,aes(x=shm,y=..density..))+
  geom_density(size  = 1.5, colour = "black")+
  labs(title= glue("{sample} ({nrow(mutation_info)} BCR)"),
       x ="SHM", y = "Density")+
  scale_x_continuous(labels = scales::percent)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(size = 10, color = "black", face = "bold",angle = 45,vjust = 1,hjust = 0.95),
    axis.line = element_line(size = 0.6),
    axis.title = element_text(size = 15, color = "black", face = "bold"),
    legend.text = element_text(size = 10, color = "black", face = "bold"),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"),
    plot.title = element_text(size = 20, color = "black", face = "bold",hjust = 0.5, vjust = 2),
  )

ggsave(filename = glue("{outdir}/{sample}_shm_density.pdf"),
       width = 8, height = 7,
       shm_mut_rate_density)

mutation_info <- mutation_info[mutation_info$locus == "IGH",]
dist_ham <- distToNearest(mutation_info, 
                          sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=40)

threshold <- 0.2

library(ggplot2)
p1 <- ggplot(subset(dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
ggsave(glue("{outdir}/{sample}_dist_ham.pdf"),
       p1, width = 4,height = 4)

# 
# if(nrow(mutation_info) < 500){
#   output <- findThreshold(dist_ham$dist_nearest, method="density")
#   threshold <- output@threshold
#   pdf(glue("{outdir}/{sample}_ham_density.pdf"), width = 4, height = 4)
#   plot(output, binwidth=0.02, title="Density Method")
#   dev.off()
# }else{
#   output <- findThreshold(dist_ham$dist_nearest, method="gmm", model="gamma-gamma")
#   threshold <- output@threshold
#   # Plot distance histogram, density estimate and optimum threshold
#   pdf(glue("{outdir}/{sample}_ham_gmm.pdf"), width = 4, height = 4)
#   plot(output, binwidth=0.02, title="GMM Method: gamma-gamma")
#   dev.off()
# }

# output <- findThreshold(dist_ham$dist_nearest, method="gmm", model="gamma-gamma")
# threshold <- output@threshold
# # Plot distance histogram, density estimate and optimum threshold
# pdf(glue("{outdir}/{sample}_ham_gmm.pdf"), width = 4, height = 4)
# plot(output, binwidth=0.02, title="GMM Method: gamma-gamma")
# dev.off()



if(!is.na(threshold)){
  write.table(data.frame(threshold),
              glue("{outdir}/{sample}_threshold.txt"),
              row.names = F, col.names = F, quote = F, sep = "\t")
}


