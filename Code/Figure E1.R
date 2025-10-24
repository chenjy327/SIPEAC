#### Figure E1B smrt igh ====
P1 <- fread('~/P1_all_ontissue_clone_list.csv') 
P1$cloneId <- paste0('P1_',P1$cloneId )
P1 <- P1[,c('cloneId','allCHitsWithScore','cloneCount')]
P2 <- fread('~/P2_all_ontissue_clone_list.csv') 
P2$cloneId <- paste0('P2_',P2$cloneId )
P2 <- P2[,c('cloneId','allCHitsWithScore','cloneCount')]
merge_vdj <- rbind(P1,P2)
merge_vdj$isotype <- substr(merge_vdj$allCHitsWithScore,1,5)
table(merge_vdj$isotype)

cur_meta <- subset(merge_vdj,isotype %in% c('IGHA1','IGHA2','IGHD*','IGHG1','IGHG2','IGHG3','IGHG4','IGHM*'))
pie_df <- cur_meta %>% dplyr::group_by(isotype) %>% 
  dplyr::summarise(value=sum(cloneCount))
pie_df[pie_df$isotype=='IGHD*',]$isotype <- 'IGHD'
pie_df[pie_df$isotype=='IGHM*',]$isotype <- 'IGHM'

pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE"))
fill_colors <- c("#fdc7d4",
                 "#ef8e89",
                 "#f6ba59",
                 "#9fc0d6",
                 "#aacb89",
                 "#5c9b4a",
                 "#7775b6",
                 "#caabd0",
                 "#225ea8")
names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "isotype", color = NULL,
           palette = fill_colors)
save_plot(p, '~/smrt_BCR_IGH', if_pdf = T,
          width = 4, height = 4)

#### Figure E1B smrt igl ====
merge_vdj$isotype <- substr(merge_vdj$allCHitsWithScore,1,4)
cur_meta <- subset(merge_vdj,isotype %in% c('IGLC','IGKC'))
pie_df <- cur_meta %>% dplyr::group_by(isotype) %>% 
  dplyr::summarise(value=sum(cloneCount))

pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4",
                 "#9fc0d6")
names(fill_colors) <- c("IGKC", "IGLC")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "isotype", color = NULL,
           palette = fill_colors)
save_plot(p, '~/Smrt_BCR_IGK_pie', if_pdf = T,
          width = 4, height = 4)



#### Figure E1B smrt tcr ====
merge_vdj$isotype <- substr(merge_vdj$allCHitsWithScore,1,4)
cur_meta <- subset(merge_vdj,isotype %in% c('TRBC','TRAC'))
pie_df <- cur_meta %>% dplyr::group_by(isotype) %>% 
  dplyr::summarise(value=sum(cloneCount))

pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
pie_df$isotype <- factor(pie_df$isotype, levels = c("TRBC", "TRAC"))
fill_colors <- c("#fdc7d4",
                 "#ef8e89")
names(fill_colors) <- c("TRBC", "TRAC")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "isotype", color = NULL,
           palette = fill_colors)
save_plot(p, '~/smrt_SC_TCR_isotype', if_pdf = T,
          width = 4, height = 4)
save_table(pie_df, "~/smrt_SC_TCR_isotype_df")

#### Figure E1B sc igh ====
sc_bcr <- fread('~/20220309_BCSA23_Bcell_filtered_contig_annotations.csv')
cur_meta <- subset(sc_bcr,chain=='IGH')
cur_meta <- subset(cur_meta,c_gene!='')
cur_meta$isotype <- cur_meta$c_gene
cur_meta$isotype %>% unique()
table(cur_meta$isotype)

pie_df <- cur_meta %>% dplyr::group_by(isotype) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
pie_df$isotype <- factor(pie_df$isotype, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE"))
fill_colors <- c("#fdc7d4",
                 "#ef8e89",
                 "#f6ba59",
                 "#9fc0d6",
                 "#aacb89",
                 "#5c9b4a",
                 "#7775b6",
                 "#caabd0",
                 "#225ea8")
names(fill_colors) <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

labs <- paste0(pie_df$isotype, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "isotype", color = NULL,
           palette = fill_colors)
save_plot(p, '~/SC_BCR_IGH', if_pdf = T,
          width = 4, height = 4)

#### Figure E1B sc igl ====

cur_meta <- subset(sc_bcr,chain!='IGH')
pie_df <- cur_meta %>% dplyr::group_by(chain) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)

fill_colors <- c("#fdc7d4",
                 "#9fc0d6")
names(fill_colors) <- c("IGK", "IGL")

labs <- paste0(pie_df$chain, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "chain", color = NULL,
           palette = fill_colors)
save_plot(p, "~/SC_BCR_IGK", if_pdf = T,
          width = 4, height = 4)

p <- ggplot() + geom_arc_bar(data=pie_df,
                             stat = "pie",
                             aes(x0=0,y0=0,r0=1,r=2,
                                 amount=ratio, fill=chain
                                 # explode=c(0.05,0.1,0.05,0.05,
                                 #           0.05,0.05,0.05,0.05,0.05,0.1,0.1)
                             )
) + scale_fill_manual(values = fill_colors)+ theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')
save_plot(p, "~/SC_BCR_IGK_circle", if_pdf = T,
          width = 4, height = 4)
print(sum(pie_df$value))
save_table(pie_df, "~/SC_BCR_IGK_circle_df")



#### Figure E1B sc tcr ====
sc_tcr <- fread('~/20220309_BCSA23_Tcel_filtered_contig_annotations.csv')
pie_df <- sc_tcr %>% dplyr::group_by(chain) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
pie_df$chain <- factor(pie_df$chain, levels = c("TRB", "TRA"))
fill_colors <- c("#fdc7d4",
                 "#ef8e89")
names(fill_colors) <- c("TRB", "TRA")

labs <- paste0(pie_df$chain, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "chain", color = NULL,
           palette = fill_colors)
save_plot(p, '~/BC_SC_TCR_isotype', if_pdf = T,
          width = 4, height = 4)
save_table(pie_df, "~/BC_SC_TCR_isotype_df")


#### Figure E1C smrt igh ####
sci_p1_clone <- readRDS("~/sci_p1_clone.rds")
sci_p1_clone$sample <- "p1"
sci_p1_clone <- sci_p1_clone %>% dplyr::select(-"locus")
sci_p2_clone <- readRDS("~/sci_p2_clone.rds")
sci_p2_clone$sample <- "p1"
sci_clone <- rbind(sci_p1_clone, sci_p2_clone)
sci_clone$clone_id <- paste(sci_clone$sample, sci_clone$cloneid, sep = "_")
sci_clone <- sci_clone[grepl("IGH", sci_clone$cloneid),]
clones <- dplyr::group_by(sci_clone, clone_id) %>% 
  dplyr::summarise(bin_count = n()) %>% as.data.frame()
clones <-  clones  %>%
  dplyr::mutate(expand_type =  dplyr::case_when(
    bin_count >= 10 ~ "Expanded >= 10",
    bin_count == 1 ~ "Singleton",
    bin_count > 1 & bin_count < 10 ~"Expanded < 10"
  ))
pie_df <- clones %>% dplyr::group_by(expand_type) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
fill_colors <- c("#946625", "#edb55a", "#65a3d1")
names(fill_colors) <- c("Expanded >= 10", "Expanded < 10", "Singleton")
labs <- paste0(pie_df$expand_type, " (", pie_df$percent, ")")




p <- ggplot() + geom_arc_bar(data=pie_df,
                             stat = "pie",
                             aes(x0=0,y0=0,r0=1,r=2,
                                 amount=ratio, fill=expand_type

                             )
) + scale_fill_manual(values = fill_colors)+ theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')
save_plot(p, "~/BC-smrt-IGH-circle", if_pdf = T,
          width = 6, height = 6)

#### Figure E1C smrt tcr ####
sci_p1_clone <- readRDS("~/sci_p1_clone.rds")
sci_p1_clone$sample <- "p1"
sci_p1_clone <- sci_p1_clone %>% dplyr::select(-"locus")
sci_p2_clone <- readRDS("~/sci_p2_clone.rds")
sci_p2_clone$sample <- "p1"
sci_clone <- rbind(sci_p1_clone, sci_p2_clone)
sci_clone$clone_id <- paste(sci_clone$sample, sci_clone$cloneid, sep = "_")
sci_clone <- sci_clone[grepl("TRB", sci_clone$cloneid),]
clones <- dplyr::group_by(sci_clone, clone_id) %>% 
  dplyr::summarise(bin_count = n()) %>% as.data.frame()
clones <-  clones  %>%
  dplyr::mutate(expand_type =  dplyr::case_when(
    bin_count >= 10 ~ "Expanded >= 10",
    bin_count == 1 ~ "Singleton",
    bin_count > 1 & bin_count < 10 ~"Expanded < 10"
  ))
pie_df <- clones %>% dplyr::group_by(expand_type) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
fill_colors <- c("#946625", "#edb55a", "#65a3d1")
names(fill_colors) <- c("Expanded >= 10", "Expanded < 10", "Singleton")
labs <- paste0(pie_df$expand_type, " (", pie_df$percent, ")")



p <- ggplot() + geom_arc_bar(data=pie_df,
                             stat = "pie",
                             aes(x0=0,y0=0,r0=1,r=2,
                                 amount=ratio, fill=expand_type
                                 # explode=c(0.05,0.1,0.05,0.05,
                                 #           0.05,0.05,0.05,0.05,0.05,0.1,0.1)
                             )
) + scale_fill_manual(values = fill_colors)+ theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')
save_plot(p, "~/BC-smrt-TRB-circle", if_pdf = T,
          width = 6, height = 6)



#### Figure E1C sc igh ####
IGH <- subset(sc_bcr,chain == 'IGH')
IGH <- subset(IGH,c_gene!='')
a <- IGH %>% 
  group_by(raw_clonotype_id) %>% 
  summarise(count=n())

Expanded10 <- a[a$count>=10,]$raw_clonotype_id
Expanded <- a[a$count>1 & a$count<10,]$raw_clonotype_id
Singleton <- a[a$count==1,]$raw_clonotype_id
a$expand_type <- ''
a[a$raw_clonotype_id %in% Expanded10,]$expand_type <- 'Expanded >= 10'
a[a$raw_clonotype_id %in% Expanded,]$expand_type <- 'Expanded < 10'
a[a$raw_clonotype_id %in% Singleton,]$expand_type <- 'Singleton'

pie_df <- a %>% dplyr::group_by(expand_type) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
fill_colors <- c( "Expanded >= 10"='#7E5E2B',"Expanded < 10" = "#f6ba59","Singleton" = "#65abde" )


p <- ggplot() + geom_arc_bar(data=pie_df,
                             stat = "pie",
                             aes(x0=0,y0=0,r0=1,r=2,
                                 amount=ratio, fill=expand_type
                                 # explode=c(0.05,0.1,0.05,0.05,
                                 #           0.05,0.05,0.05,0.05,0.05,0.1,0.1)
                             )
) + scale_fill_manual(values = fill_colors)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')+ggtitle(paste0("all_",sum(pie_df$value)))
save_plot(p, glue('~/SC_IGH_expand_circle'), if_pdf = T,
          width = 6, height = 6)
#### Figure E1C sc tcr ####

sc_TCR <- fread('~/20220309_BCSA23_Tcel_filtered_contig_annotations.csv')
TRB <- subset(sc_TCR,chain == 'TRB')
a <- TRB %>% 
  group_by(raw_clonotype_id) %>% 
  summarise(count=n())
Expanded10 <- a[a$count>=10,]$raw_clonotype_id
Expanded <- a[a$count>1 & a$count<10,]$raw_clonotype_id
Singleton <- a[a$count==1,]$raw_clonotype_id
a$expand_type <- ''
a[a$raw_clonotype_id %in% Expanded10,]$expand_type <- 'Expanded >= 10'
a[a$raw_clonotype_id %in% Expanded,]$expand_type <- 'Expanded < 10'
a[a$raw_clonotype_id %in% Singleton,]$expand_type <- 'Singleton'

pie_df <- a %>% dplyr::group_by(expand_type) %>% 
  dplyr::summarise(value = n())
pie_df$ratio <- pie_df$value/sum(pie_df$value)
pie_df$percent <- percent(pie_df$ratio, accuracy = 0.1)
fill_colors <- c( "Expanded >= 10"='#7E5E2B',"Expanded < 10" = "#f6ba59","Singleton" = "#65abde" )

labs <- paste0(pie_df$expand_type, " (", pie_df$percent, ")")

p <- ggpie(pie_df, "value", label = labs,
           lab.pos = "in", lab.font = "white",
           fill = "expand_type", color = NULL,
           palette = fill_colors)
save_plot(p, '~/SC_TRB_expand', if_pdf = T,
          width = 4, height = 4)
save_table(pie_df, '~/SC_TRB_expand_df')

p <- ggplot() + geom_arc_bar(data=pie_df,
                             stat = "pie",
                             aes(x0=0,y0=0,r0=1,r=2,
                                 amount=ratio, fill=expand_type
                                 # explode=c(0.05,0.1,0.05,0.05,
                                 #           0.05,0.05,0.05,0.05,0.05,0.1,0.1)
                             )
) + scale_fill_manual(values = fill_colors)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')+ggtitle(paste0("all_",sum(pie_df$value)))
save_plot(p, glue('~/T注释/SC_TRB_expand_circle'), if_pdf = T,
          width = 6, height = 6)

#### Figure E1D ####
bcr_merge <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_bcr.rds")
tcr_merge <- readRDS("/Users/mac/Desktop/Rproject/Pacbio/文章/RDS/merged_tcr.rds")

tmp_tcr <- tcr_merge[tcr_merge$sample == "2907T",]
tmp_bcr <- bcr_merge[bcr_merge$sample == "2907T",]

obj <- readRDS(glue("~/2907T_scRefST2934T_random100_marker20_spotlight_Spatial.rds"))
obj@meta.data <- read.table("~/2907T.txt",sep = "\t",header = T)
rownames(obj@meta.data) <- sapply(obj@meta.data$BinID, function(x)strsplit(x, "_")[[1]][2]) %>% as.character()
obj@meta.data$BinID <- rownames(obj@meta.data)
cur_gene <- c("IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2","IGHM","IGHD","IGLC2","IGKC","TRBC1", "TRBC2","TRAC", "TRA2B")
for (i in cur_gene) {
  SpatialFeaturePython(obj = obj,
                       prefix = glue("~/bcrgene_count"),
                       if_gene = T,
                       plot_item = i,
                       name = glue("2907T_{i}_count"),
                       if_count = T)
}



bcr_isotype <- c('IGHA1', 'IGHA2',  'IGHD', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4',  'IGHM',  'IGKC', 'IGLC')
tmp_bcr$isotype[tmp_bcr$isotype %in% c('IGLC1', 'IGLC2', 'IGLC3')] <- "IGLC"
for (i in bcr_isotype) {
  print(glue("{i}"))
  isotype <- subset(tmp_bcr, tmp_bcr$isotype == i)
  isotype <- isotype$Bin50 %>% table() %>% as.data.frame()
  colnames(isotype) <- c("BinID",i)
  isotype[,i] <- i
  obj_i <- obj
  obj_i@meta.data <- dplyr::left_join(obj_i@meta.data, isotype, by = "BinID") %>% column_to_rownames(var="BinID")
  obj_i@meta.data[,i] <- ifelse(is.na(obj_i@meta.data[,i]), "Other", obj_i@meta.data[,i])
  SpatialDimPython(obj = obj_i, 
                   prefix = "~/bcr_smrt",
                   plot_item = i,
                   plot_order = c(i,"Other"), 
                   tmp_color = c("red","#E9ECEB"),
                   name = glue("2907T_{i}"))
}




tcr_isotype <- c("TRA","TRB",  "TRG", "TRD")
for (i in tcr_isotype) {
  print(glue("{i}"))
  isotype <- subset(tmp_tcr, tmp_tcr$topChains == i)
  isotype <- isotype$Bin50_BinID %>% table() %>% as.data.frame()
  colnames(isotype) <- c("BinID",i)
  isotype[,i] <- i
  # isotype$BinID <- paste0(glue("2907T_"),isotype$BinID)
  obj_i <- obj
  
  obj_i@meta.data <- dplyr::left_join(obj_i@meta.data, isotype, by = "BinID") %>% column_to_rownames(var="BinID")
  obj_i@meta.data[,i] <- ifelse(is.na(obj_i@meta.data[,i]), "Other", obj_i@meta.data[,i])
  SpatialDimPython(obj = obj_i, 
                   prefix =  "~/tcr_smrt",
                   plot_item = i,
                   plot_order = c(i,"Other"), 
                   tmp_color = c("#00BF05","#E9ECEB"),
                   name = glue("2907T_{i}"))
}