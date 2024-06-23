library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsignif)
library(scales)
library(EnhancedVolcano)


source("utils.R")

my_colors5=c("#DAC5DF","#BBD5EB","#6C8DC1","#EAAB9E","#DD6E60")
my_colors7=c("#e2ad9f","#809b5c","#c0d6eb","#d37562","#d7c6de","#9b7aad","#90aed2")

my_colors2 <- c("#ef8775", "#89afd2")

my_colors2_2 <- c(my_colors5[4], my_colors5[2])
my_colors4 = my_colors7[c(1,3,2,4)]

# #create colors with any number
# cols = colorRampPalette(colors = c("#e2ad9f","#d37562","#d7c6de","#9b7aad",
#                                             "#90aed2","#c0d6eb"))(length(unique(tcell.batch.corrected$subtype1))))

# expression dotplot
  # scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")


packageVersion("Seurat")

packageVersion("SeuratObject")

meta=read.csv("processed_data/data_B2-19/metadata_all_ct_subtype_condition_noTLB.csv", row.names = 1)


table(meta$cell_type)

dim(meta)



GC <- readRDS("plots/data/GC_no_TLB_new_meta.rds")

# used
my_colors7=c("#e2ad9f","#809b5c","#c0d6eb","#d37562","#d7c6de","#9b7aad","#90aed2")

options(repr.plot.width=6, repr.plot.height=5)
p<-DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols=my_colors7)
plot(p)
pdf("plots/figures/fig1B_umap.pdf",width=6, height=5)
plot(p)
dev.off()

options(repr.plot.width=5.5, repr.plot.height=4.5)
p <- DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "patient",label=FALSE)
plot(p)
pdf("plots/figures/figS1B_umap.pdf",width=6, height=5)
plot(p)
dev.off()

options(repr.plot.width=5.5, repr.plot.height=4.5)
p <- DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE)
plot(p)
pdf("plots/figures/fig1A_umap.pdf",width=5.5, height=4.5)
plot(p)
dev.off()

options(repr.plot.width=6, repr.plot.height=5)
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols=c("#d7c6de","#c0d6eb","#738bc1","#e2ad9f","#d37562","#809b5c","#9b7aad"))

options(repr.plot.width=6, repr.plot.height=5)
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols=c("#e2ad9f","#d7c6de","#c0d6eb","#738bc1","#d37562","#809b5c","#6ea39f"))

user_markers<-c("CD2", "CD3E", "CD3D", #T cell
                "CD68", "CSF1R","CD14", #macro
                "TPSAB1","TPSB2",#mast cell
#                 "FCER1G","XCL2","XCL1","TYROBP","TRDC",#NK
                
                "ACTA2",  "COL1A1", "LUM", #"DCN", #fibro
                "EPCAM", "KRT18", "PGC",# "TFF1", "MUC5AC", "GIF", "CHGA", # Epithelial
                "PECAM1", "VWF", #endo
#                 "CXCL8","CSF3R", "S100A8", "S100A9", #neutrophil
                "CD79A"#, # Bcell
                
#                 "IGHG1", "IGKC", "SDC1", "CD38"#plasma
               )

options(repr.plot.width=9, repr.plot.height=5)
DotPlot(GC, features = user_markers, group.by = "cell_type") +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

pdf("plots/figures/fig1C_expression.pdf",width=9, height=5)
DotPlot(GC, features = user_markers, group.by = "cell_type") +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")
dev.off()




# GC_meta <- GC@meta.data
meta=read.csv("processed_data/data_B2-19/metadata_all_ct_subtype_condition_noTLB.csv", row.names = 1)


options(repr.plot.width=4, repr.plot.height=3.5)
my_colors7=c("#e2ad9f","#809b5c","#c0d6eb","#d37562","#d7c6de","#9b7aad","#90aed2")

cluster_group <- table(meta$cell_type, meta$LN_condition)
cluster_group <- melt(cluster_group, id.vars = "LN_condition")
colnames(cluster_group) <- c("cell_type", "LN_condition", "Freq")
cluster_group$cell_type <- as.factor(cluster_group$cell_type)

cluster_group$LN_condition <- factor(cluster_group$LN_condition,
       levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
# pdf("cell_typing/cell_type_composition.pdf")
p1d <- ggplot(cluster_group, aes(x = LN_condition, y = Freq, fill = cell_type)) +
  geom_bar(position = "fill", stat = "identity")+
  scale_fill_manual(values = my_colors7) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(panel.background = element_blank())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line())+
  theme(axis.text.y = element_text(color = "black"))+
  theme(axis.text.x = element_text(color = "black"))+
  ylab("Fraction among all cells")
plot(p1d)
# dev.off()

ggsave(filename = "plots/figures/fig1D_fraction.pdf", plot = p1d, width = 4, height = 3.5)

# try color
options(repr.plot.width=4, repr.plot.height=3.5)
my_colors7=c("#e2ad9f","#809b5c","#c0d6eb","#d37562","#d7c6de","#9b7aad","#90aed2")

cluster_group <- table(meta$cell_type, meta$LN_condition)
cluster_group <- melt(cluster_group, id.vars = "LN_condition")
colnames(cluster_group) <- c("cell_type", "LN_condition", "Freq")
cluster_group$cell_type <- as.factor(cluster_group$cell_type)

cluster_group$LN_condition <- factor(cluster_group$LN_condition,
       levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
# pdf("cell_typing/cell_type_composition.pdf")
p1d <- ggplot(cluster_group, aes(x = LN_condition, y = Freq, fill = cell_type)) +
  geom_bar(position = "fill", stat = "identity")+
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(panel.background = element_blank())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line())+
  theme(axis.text.y = element_text(color = "black"))+
  theme(axis.text.x = element_text(color = "black"))+
  ylab("Fraction among all cells")
plot(p1d)
# dev.off()

ggsave(filename = "plots/figures/fig1D_fraction.pdf", plot = p1d, width = 4, height = 3.5)



# not used
options(repr.plot.width=4, repr.plot.height=3.5)
my_colors7=c("#e2ad9f","#809b5c","#c0d6eb","#d37562","#d7c6de","#9b7aad","#90aed2")

cluster_group <- table(GC_meta$cell_type, GC_meta$LN_condition)
cluster_group <- melt(cluster_group, id.vars = "LN_condition")
colnames(cluster_group) <- c("cell_type", "LN_condition", "Freq")
cluster_group$cell_type <- as.factor(cluster_group$cell_type)

cluster_group$LN_condition <- factor(cluster_group$LN_condition,
       levels=c("N", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
# pdf("cell_typing/cell_type_composition.pdf")
ggplot(cluster_group, aes(x = LN_condition, y = Freq, fill = cell_type)) +
  geom_bar(position = "fill", stat = "identity")+
  scale_fill_manual(values = my_colors7) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(panel.background = element_blank())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line())+
  theme(axis.text.y = element_text(color = "black"))+
  theme(axis.text.x = element_text(color = "black"))+
  ylab("% among all cells")

# dev.off()

# ggplot(cluster_group, aes(x = cell_type, y = Freq, fill = LN_condition)) +
#   geom_bar(position = "fill", stat = "identity")+
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
#   theme(panel.background = element_blank())+
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
#   theme(axis.line = element_line())+
#   theme(axis.text.y = element_text(color = "black"))+
#   theme(axis.text.x = element_text(color = "black"))

# ggplot(cluster_group, aes(x = cell_type, y = Freq, fill = LN_condition)) +
#   geom_bar(stat = "identity")+
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
#   theme(panel.background = element_blank())+
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
#   theme(axis.line = element_line())+
#   theme(axis.text.y = element_text(color = "black"))+
#   theme(axis.text.x = element_text(color = "black"))

# ggplot(cluster_group, aes(x = LN_condition, y = Freq, fill = cell_type)) +
#   geom_bar(stat = "identity")+
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
#   theme(panel.background = element_blank())+
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
#   theme(axis.line = element_line())+
#   theme(axis.text.y = element_text(color = "black"))+
#   theme(axis.text.x = element_text(color = "black"))

# 

epithelial.batch.corrected <- readRDS("plots/data/epithelial_new_meta.rds")
Idents(epithelial.batch.corrected) <- "LN_condition"

epithelial.batch.corrected$cnv_score[which(is.na(epithelial.batch.corrected$cnv_score))]<-0
epithelial.batch.corrected$cnv_score[epithelial.batch.corrected$cnv_score>0.015]<-0.015
#为了方便展示调整了部分CNV的分数

# fig 1F used
options(repr.plot.width = 3, repr.plot.height = 3)
my_colors2 <- c("#ef8775", "#89afd2")

# 绘制小提琴图并应用自定义颜色和y轴刻度
vln_plot <- VlnPlot(epithelial.batch.corrected, 
                    features = "cnv_score",
                    idents = c("Pri.GC", "Met.LN"),
                    pt.size = 0) + #, y.max = 0.025
  geom_signif(comparisons = list(c("Pri.GC", "Met.LN")), 
              map_signif_level = TRUE, y_position = 0.018) +
  NoLegend() +   
  scale_fill_manual(values = c("Pri.GC" = "#ef8775", "Met.LN" = "#89afd2")) +
 # scale_fill_manual(values = c("Pri.GC" = my_colors5[2], "Met.LN" = my_colors5[4])) +
  scale_y_continuous(limits = c(0, 0.02), breaks = c(0, 0.01, 0.02)) + # 自定义y轴刻度
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))   # 取消x轴标签旋转
plot(vln_plot)

ggsave(filename = "plots/figures/fig1E_cnv_score.pdf", plot = vln_plot, width = 3, height = 3)

# # color not used备选2色方案
options(repr.plot.width = 3, repr.plot.height = 3)


# 绘制小提琴图并应用自定义颜色和y轴刻度
vln_plot <- VlnPlot(epithelial.batch.corrected, 
                    features = "cnv_score",
                    idents = c("T", "LN_meta"),
                    pt.size = 0) + #, y.max = 0.025
  geom_signif(comparisons = list(c("T", "LN_meta")), 
              map_signif_level = TRUE, y_position = 0.018) +
  NoLegend() +   
  scale_fill_manual(values = c("T" = my_colors5[2], "LN_meta" = my_colors5[4])) +
  scale_y_continuous(limits = c(0, 0.02), breaks = c(0, 0.01, 0.02)) + # 自定义y轴刻度
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1))   # 取消x轴标签旋转
plot(vln_plot)

meta=read.csv("processed_data/data_B2-19/metadata_all_ct_subtype_condition_noTLB.csv", row.names = 1)

sample_meta <- unique(meta[, c("sample", "region", "patient", 'LN_condition','LN_station','LN_condition.LN_station',
                              "patient_condition","patient_condition.LN_condition",
                               "Lauren","Lauren.LN_condition")])

# # 对sample进行normalize后的结果，把所有细胞类型作为一个barplot作图
# options(repr.plot.width=9, repr.plot.height=9)
# sample_ct_comp <- list()
# for(ct_i in c("all", "all_sub", unique(meta$cell_type1))){
#     if(ct_i == "all"){
#         sample_ct_comp[[ct_i]] <- table(meta$cell_type1, meta$sample)
#     }else if(ct_i == "all_sub"){
#         sample_ct_comp[[ct_i]] <- table(meta$subtype1, meta$sample)
# #     }else if(ct_i == "all_sub1"){
# #         sample_ct_comp[[ct_i]] <- table(meta$subtype1, meta$sample)
#     }else{
#         sample_ct_comp[[ct_i]] <- table(meta[meta$cell_type1==ct_i,]$subtype1, meta[meta$cell_type1==ct_i,]$sample)
#     }
#     sample_ct_comp[[ct_i]] <- melt(sample_ct_comp[[ct_i]], id.vars = "sample")
#     colnames(sample_ct_comp[[ct_i]]) <- c("cell_type1", "sample", "Fraction")
#     sample_ct_comp[[ct_i]]$cell_type1 <- as.factor(sample_ct_comp[[ct_i]]$cell_type1)
    
#     # normalize到1
#     for(condition in unique(sample_ct_comp[[ct_i]]$sample)){
#         mask=sample_ct_comp[[ct_i]]$sample==condition
#         sample_ct_comp[[ct_i]]$Fraction[mask]=sample_ct_comp[[ct_i]]$Fraction[mask]/sum(sample_ct_comp[[ct_i]]$Fraction[mask])
#     }

# }


# for(ct_i in names(sample_ct_comp)){
#     print(ct_i)
#     sample_ct_comp[[ct_i]]<-merge(sample_ct_comp[[ct_i]], sample_meta, by = "sample")
# }

# saveRDS(sample_ct_comp,"processed_data/data_B2-19/cell_type1_composition/sample_ct_comp.rds")

sample_ct_comp <- readRDS("processed_data/data_B2-19/cell_type1_composition/sample_ct_comp.rds")
# from cell_type_composition.ipynb

#几个subtype放在一张图上的示例
my_colors5=c("#DAC5DF","#BBD5EB","#6C8DC1","#EAAB9E","#DD6E60")
sample_ct_comp[["all"]]$LN_condition <- factor(sample_ct_comp[["all"]]$LN_condition,
                                                               levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
tmp <- sample_ct_comp[["all"]][sample_ct_comp[["all"]]$cell_type1 !="Remove",]

options(repr.plot.width=12, repr.plot.height=10)
p <- ggplot(tmp, aes(x=LN_condition, y=Fraction)) + 
        geom_boxplot(aes(fill=LN_condition)) +
        scale_fill_manual(values = my_colors5) +
              geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                          map_signif_level=TRUE) +
          geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  facet_wrap(. ~  cell_type1, ncol = 4) +
  labs(title="",x="LN_condition", y = "Fraction")+
  scale_y_continuous(limits = c(0, 1.))+
  theme_classic()
plot(p)   

ggsave(filename = "plots/figures/figS1A_fraction.pdf", plot = p, width = 12, height = 10)

meta_Temra <- meta[meta$subtype1 %in% c("CD8_GZMK+ Temra","CD8_GZMB+ Temra"),]

meta_Temra$LN_condition <- factor(meta_Temra$LN_condition,
                                                               levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
# fig
options(repr.plot.width=5, repr.plot.height=1.8)
cluster_group <- table(meta_Temra$subtype1, meta_Temra$LN_condition)
cluster_group <- melt(cluster_group, id.vars = "LN_condition")
colnames(cluster_group) <- c("subtype1", "LN_condition", "Fraction")
cluster_group$subtype1 <- as.factor(cluster_group$subtype1)

p1<-ggplot(cluster_group, aes(x = subtype1, y = Fraction, fill = LN_condition)) +
      geom_bar(position = "fill", stat = "identity")+
        scale_fill_manual(values = my_colors5) +

#       theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
      theme(panel.background = element_blank())+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(axis.line = element_line())+
      theme(axis.text.y = element_text(color = "black"))+
      theme(axis.text.x = element_text(color = "black"))+
coord_flip() 

plot(p1)
ggsave(filename = "plots/figures/fig3I_GZMB_fraction.pdf", plot = p1, width = 5, height = 1.8)




#特别关注几个亚型，可以单独作图
options(repr.plot.width=2.5, repr.plot.height=3.8)
ct_i="T-CD8"
subtype1 = "CD8_TIM3+ Trm"
print(subtype1)
    title=paste0(subtype1)
    sample_ct_comp[[ct_i]]$LN_condition <- factor(sample_ct_comp[[ct_i]]$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
    p <- ggplot(sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="none",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/fig2C_tim3_fraction.pdf", plot = p, width = 2.5, height = 3.8)
        

#特别关注几个亚型，可以单独作图
options(repr.plot.width=2.5, repr.plot.height=3.8)
ct_i="T-CD8"
subtype1 = "CD8_LAYN+ Tex"
print(subtype1)
    title=paste0(subtype1)
    sample_ct_comp[[ct_i]]$LN_condition <- factor(sample_ct_comp[[ct_i]]$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
    p <- ggplot(sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="none",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/fig2D_layn_fraction.pdf", plot = p, width = 2.5, height = 3.8)
        

#特别关注几个亚型，可以单独作图
options(repr.plot.width=2.5, repr.plot.height=3.8)
ct_i="T-CD8"
subtype1 = "CD8_CCL5+ Tex"
print(subtype1)
    title=paste0(subtype1)
    sample_ct_comp[[ct_i]]$LN_condition <- factor(sample_ct_comp[[ct_i]]$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
    p <- ggplot(sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="none",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/fig2D_ccl5_fraction.pdf", plot = p, width = 2.5, height = 3.8)
        

#特别关注几个亚型，可以单独作图
options(repr.plot.width=3, repr.plot.height=3.8)
ct_i="T-CD4"
subtype1 = "CD4_Tex"
print(subtype1)
    title=paste0(subtype1)
    sample_ct_comp[[ct_i]]$LN_condition <- factor(sample_ct_comp[[ct_i]]$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
    p <- ggplot(sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="none",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/figS3A_cd4_tex_fraction.pdf", plot = p, width = 3, height = 3.8)
        

#特别关注几个亚型，可以单独作图
options(repr.plot.width=3, repr.plot.height=4)
ct_i="T-CD8"
subtype1 = "CD8_HLA-DRA+ Tem"
print(subtype1)
    title=paste0(subtype1)
    sample_ct_comp[[ct_i]]$LN_condition <- factor(sample_ct_comp[[ct_i]]$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
    p <- ggplot(sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="none",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/figS5C_HLA_fraction.pdf", plot = p, width = 3, height = 4)
        

#特别关注几个亚型，可以单独作图
ct_i="DC"
ncol=1  
options(repr.plot.width=3.5, repr.plot.height=5)
tmp <- sample_ct_comp[[ct_i]]
tmp <- tmp[tmp$cell_type1 %in% c('FCER1A+ cDC2','NLRP3+ cDC2'),]

# reorder conditions in plots
tmp$LN_condition <- factor(tmp$LN_condition,
                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))
p <- ggplot(tmp, aes(x=LN_condition, y=Fraction)) + 
      geom_boxplot(aes(fill=LN_condition)) +
      scale_fill_manual(values = my_colors5) +
      geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                  map_signif_level=TRUE,y_position = 0.9) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_classic()+
facet_wrap(. ~  cell_type1, ncol = ncol) +
  theme(text = element_text(size=11),
        legend.position="none")+
   ylab("Fraction")+
   xlab("")
plot(p)
ggsave(filename = "plots/figures/fig4F_cDC2_fraction.pdf", plot = p, width = 3.5, height = 5)



#每个cell_type1的subtype做一个整的图
options(repr.plot.width=10, repr.plot.height=4)

ct_i ="NK"
ncol=3   
    
sample_ct_comp[[ct_i]]$LN_condition <- factor(sample_ct_comp[[ct_i]]$LN_condition,
                                              levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN", "PBMC"))

p <- ggplot(sample_ct_comp[[ct_i]], 
            aes(x=LN_condition, y=Fraction)) + 
        geom_boxplot(aes(fill=LN_condition)) +
      scale_fill_manual(values = my_colors5) +
              geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                          map_signif_level=TRUE,y_position = 0.9) +
          geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  facet_wrap(. ~  cell_type1, ncol = ncol) +
  labs(title="",x="LN_condition", y = "Fraction")+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+    
  scale_y_continuous(limits = c(0, 1.))+
  theme_classic()
plot(p)   

ggsave(filename = "plots/figures/figS4D_nk_fraction.pdf", plot = p, width = 10, height = 4)



#特别关注几个亚型，可以单独作图
# 这里的condition少了PBMC
options(repr.plot.width=4, repr.plot.height=4.2)
ct_i="T-CD8"
subtype1 = 'CD8_GZMB+ Temra'
print(subtype1)
title=paste0(subtype1)#, " in ", ct_i
    # reorder conditions in plots
tmp <- sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$LN_condition %in% c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"),]
    tmp$LN_condition <- factor(tmp$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"))
p <- ggplot(tmp[tmp$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5[1:4]) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="right",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/figS6B_GZMB_fraction.pdf", plot = p, width = 4, height = 4.2)
        


#特别关注几个亚型，可以单独作图
# 这里的condition少了PBMC
options(repr.plot.width=4, repr.plot.height=4.2)
ct_i="T-CD8"
subtype1 = 'CD8_GZMK+ Temra'
print(subtype1)
title=paste0(subtype1)#, " in ", ct_i
    # reorder conditions in plots
tmp <- sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$LN_condition %in% c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"),]
    tmp$LN_condition <- factor(tmp$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"))
p <- ggplot(tmp[tmp$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5[1:4]) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="right",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/figS6C_GZMK_fraction.pdf", plot = p, width = 4, height = 4.2)
        


#特别关注几个亚型，可以单独作图
# 这里的condition少了PBMC
options(repr.plot.width=4, repr.plot.height=4.2)
ct_i="T-CD8"
subtype1 = 'CD8_GZMK+ Temra'
print(subtype1)
title=paste0(subtype1)#, " in ", ct_i
    # reorder conditions in plots
tmp <- sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$LN_condition %in% c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"),]
    tmp$LN_condition <- factor(tmp$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"))
p <- ggplot(tmp[tmp$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5[1:4]) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="right",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/figS6C_GZMK_fraction.pdf", plot = p, width = 4, height = 4.2)
        


#特别关注几个亚型，可以单独作图
# 这里的condition少了PBMC
options(repr.plot.width=4, repr.plot.height=4.2)
ct_i="Plasma"
subtype1 = 'RRM2+ Plasma'
print(subtype1)
title=paste0(subtype1)
    # reorder conditions in plots
tmp <- sample_ct_comp[[ct_i]][sample_ct_comp[[ct_i]]$LN_condition %in% c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"),]
tmp$LN_condition <- factor(tmp$LN_condition,levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"))

p <- ggplot(tmp[tmp$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5[1:4]) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="right",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction among Plasma cells")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/figS10C_RRM2_fraction.pdf", plot = p, width = 4, height = 4.2)
        


unique(meta$subtype1[meta$cell_type=="B cell"])

meta$subtype0 <- meta$subtype1
meta$subtype0[meta$cell_type1=="Plasma"] <- "Plasma"

# 对sample进行normalize后的结果，把所有细胞类型作为一个barplot作图
options(repr.plot.width=9, repr.plot.height=9)
sample_ct_comp <- list()
ct_i = "B cell"

tmp <- table(meta[meta$cell_type==ct_i,]$subtype0, meta[meta$cell_type==ct_i,]$sample)
tmp <- melt(tmp, id.vars = "sample")
colnames(tmp) <- c("cell_type1", "sample", "Fraction")
tmp$cell_type1 <- as.factor(tmp$cell_type1)

# normalize到1
for(condition in unique(tmp$sample)){
    mask=tmp$sample==condition
    tmp$Fraction[mask]=tmp$Fraction[mask]/sum(tmp$Fraction[mask])
}
tmp<-merge(tmp, sample_meta, by = "sample")



#特别关注几个亚型，可以单独作图
# 这里的condition少了PBMC
options(repr.plot.width=4, repr.plot.height=4.2)
ct_i="B cell"
subtype1 = 'Plasma'
print(subtype1)
title=paste0(subtype1)#, " in ", ct_i
    # reorder conditions in plots
tmp <- tmp[tmp$LN_condition %in% c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"),]
tmp$LN_condition <- factor(tmp$LN_condition,
                                                           levels=c("Adj.Nor", "Pri.GC", "Met.LN", "Neg.LN"))
p <- ggplot(tmp[tmp$cell_type1==subtype1,], aes(x=LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=LN_condition)) +
          scale_fill_manual(values = my_colors5[1:4]) +
          geom_signif(comparisons = list(c("Pri.GC", "Met.LN"),c("Met.LN", "Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(title) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
            legend.position="right",
           plot.title = element_text(hjust = 0.5))+
       ylab("Fraction among all B cells")+
       xlab("") 
    plot(p)
ggsave(filename = "plots/figures/figS10D_plasma_fraction.pdf", plot = p, width = 4, height = 4.2)
        


unique(meta$patient_condition.LN_condition)

mask <- meta$patient_condition.LN_condition %in% c('Meta_Neg.LN','Meta_PBMC','Meta_Pri.GC',
                                                   'Neg_Neg.LN','Neg_PBMC','Neg_Pri.GC')
meta_patient_condition.LN_condition <- meta[mask,]


# 对sample进行normalize后的结果，做矩阵
patient_condition.LN_condition_sample_ct_comp <- list()
for(ct_i in c("all", "all_sub", 'T-CD8','T-CD4','B cell','NK','DC','Plasma','Macro','Tgd','Mono','Neutro')){
    if(ct_i == "all"){
        patient_condition.LN_condition_sample_ct_comp[[ct_i]] <- table(meta_patient_condition.LN_condition$cell_type1, meta_patient_condition.LN_condition$sample)
    }else if(ct_i == "all_sub"){
        patient_condition.LN_condition_sample_ct_comp[[ct_i]] <- table(meta_patient_condition.LN_condition$subtype1, meta_patient_condition.LN_condition$sample)
    }else{
        patient_condition.LN_condition_sample_ct_comp[[ct_i]] <- table(meta_patient_condition.LN_condition[meta_patient_condition.LN_condition$cell_type1==ct_i,]$subtype1, meta_patient_condition.LN_condition[meta_patient_condition.LN_condition$cell_type1==ct_i,]$sample)
    }
    patient_condition.LN_condition_sample_ct_comp[[ct_i]] <- melt(patient_condition.LN_condition_sample_ct_comp[[ct_i]], id.vars = "sample")
    colnames(patient_condition.LN_condition_sample_ct_comp[[ct_i]]) <- c("cell_type1", "sample", "Fraction")
    patient_condition.LN_condition_sample_ct_comp[[ct_i]]$cell_type1 <- as.factor(patient_condition.LN_condition_sample_ct_comp[[ct_i]]$cell_type1)
    
    # normalize到1
    for(condition in unique(patient_condition.LN_condition_sample_ct_comp[[ct_i]]$sample)){
        mask=patient_condition.LN_condition_sample_ct_comp[[ct_i]]$sample==condition
        patient_condition.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask]=patient_condition.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask]/sum(patient_condition.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask])
    }
    
    print(ct_i)
    patient_condition.LN_condition_sample_ct_comp[[ct_i]]<-merge(patient_condition.LN_condition_sample_ct_comp[[ct_i]], 
                                                                 sample_meta, by = "sample")
    
}


meta$subtype2 <- meta$subtype1
meta$subtype2[meta$subtype2 %in% c("MS4A1+ Plasmablast","STMN1+ Plasmablast")] <- "GCB"

mask <- meta$patient_condition.LN_condition %in% c('Meta_Neg.LN','Meta_PBMC','Meta_Pri.GC',
                                                   'Neg_Neg.LN','Neg_PBMC','Neg_Pri.GC')
meta_patient_condition.LN_condition <- meta[mask,]


unique(meta_patient_condition.LN_condition$subtype2)

# 对sample进行normalize后的结果，把所有细胞类型作为一个barplot作图
patient_condition.LN_condition_sample_ct_comp <- list()
for(ct_i in c('B cell')){
    patient_condition.LN_condition_sample_ct_comp[[ct_i]] <- table(meta_patient_condition.LN_condition[meta_patient_condition.LN_condition$cell_type==ct_i,]$subtype2, 
                                                                       meta_patient_condition.LN_condition[meta_patient_condition.LN_condition$cell_type==ct_i,]$sample)
    patient_condition.LN_condition_sample_ct_comp[[ct_i]] <- melt(patient_condition.LN_condition_sample_ct_comp[[ct_i]], id.vars = "sample")
    colnames(patient_condition.LN_condition_sample_ct_comp[[ct_i]]) <- c("cell_type", "sample", "Fraction")
    patient_condition.LN_condition_sample_ct_comp[[ct_i]]$cell_type <- as.factor(patient_condition.LN_condition_sample_ct_comp[[ct_i]]$cell_type)
    
    # normalize到1
    for(condition in unique(patient_condition.LN_condition_sample_ct_comp[[ct_i]]$sample)){
        mask=patient_condition.LN_condition_sample_ct_comp[[ct_i]]$sample==condition
        patient_condition.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask]=patient_condition.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask]/sum(patient_condition.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask])
    }
    
    print(ct_i)
    patient_condition.LN_condition_sample_ct_comp[[ct_i]]<-merge(patient_condition.LN_condition_sample_ct_comp[[ct_i]], 
                                                                 sample_meta, by = "sample")
    
}


tmp <- patient_condition.LN_condition_sample_ct_comp[["B cell"]]
tmp <- tmp[tmp$LN_condition=="Pri.GC",]

options(repr.plot.width=3.5, repr.plot.height=3)

for(subtype2 in c("GCB")){
    tmp$patient_condition.LN_condition <- factor(tmp$patient_condition.LN_condition,
                                                           levels=c("Meta_Pri.GC", "Neg_Pri.GC"))
    p <- ggplot(tmp[tmp$cell_type==subtype2,], aes(x=patient_condition.LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=patient_condition.LN_condition)) +
        scale_fill_manual(values = my_colors2) +
    
          # geom_signif(comparisons = list(c("Meta_Pri.GC", "Neg_Pri.GC")), 
          #             map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(subtype2) +
      theme_classic()+
      theme(axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5))+
    xlab("") +
    labs(fill = NULL)
    theme(legend.position="right")
    plot(p)
    print(subtype2)
}
ggsave(filename = "plots/figures/fig5B_GCB_fraction.pdf", plot = p, width = 3.5, height = 3)

tmp <- patient_condition.LN_condition_sample_ct_comp[["T-CD8"]]
tmp <- tmp[tmp$LN_condition=="Pri.GC",]

options(repr.plot.width=3.5, repr.plot.height=4)

for(subtype2 in c("CD8_HLA-DRA+ Tem")){
    tmp$patient_condition.LN_condition <- factor(tmp$patient_condition.LN_condition,
                                                           levels=c("Meta_Pri.GC", "Neg_Pri.GC"))
    p <- ggplot(tmp[tmp$cell_type==subtype2,], aes(x=patient_condition.LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=patient_condition.LN_condition)) +
        scale_fill_manual(values = my_colors2) +
    
          geom_signif(comparisons = list(c("Meta_Pri.GC", "Neg_Pri.GC")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(subtype2) +
      theme_classic()+
      theme(axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5))+
    xlab("") +
    labs(fill = NULL)
    theme(legend.position="right")
    plot(p)
    print(subtype2)
}
ggsave(filename = "plots/figures/fig3C_HLA_fraction.pdf", plot = p, width = 3.5, height = 4)

# for(ct_i in names(patient_condition.LN_condition_sample_ct_comp)){
ct_i="Neutro"
print(ct_i)
ncol=3   
options(repr.plot.width=9, repr.plot.height=3)
tmp_neutro <- patient_condition.LN_condition_sample_ct_comp[[ct_i]]
tmp_neutro<- tmp_neutro[tmp_neutro$LN_condition=="Neg.LN",]
tmp_neutro$patient_condition.LN_condition <- factor(tmp_neutro$patient_condition.LN_condition,
                                                       levels=c(#"Meta_Pri.GC", "Neg_Pri.GC",
                                     "Meta_Neg.LN", "Neg_Neg.LN"))
p <- ggplot(tmp_neutro, 
            aes(x=patient_condition.LN_condition, y=Fraction)) + 
      geom_boxplot(aes(fill=patient_condition.LN_condition)) +
        scale_fill_manual(values = my_colors2) +
#       geom_signif(comparisons = list(
#                                      c("Meta_Neg.LN", "Neg_Neg.LN")
#                                      ), 
#                   map_signif_level=TRUE,y_position = 1.02) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
scale_y_continuous(limits = c(0, 1.1))+
facet_wrap(. ~  cell_type1, ncol = ncol) +
  theme_classic()+
    xlab("") +
    labs(fill = NULL)+
  theme(axis.text.x=element_blank(),text = element_text(size=12))
plot(p)
ggsave(filename = "plots/figures/figS8A_neutro_fraction.pdf", plot = p, width = 9, height = 3)

# for(ct_i in names(patient_condition.LN_condition_sample_ct_comp)){
ct_i="Neutro"
print(ct_i)
ncol=3   
options(repr.plot.width=9, repr.plot.height=4)
tmp_neutro <- patient_condition.LN_condition_sample_ct_comp[[ct_i]]
tmp_neutro<- tmp_neutro[tmp_neutro$LN_condition=="Neg.LN",]
tmp_neutro$patient_condition.LN_condition <- factor(tmp_neutro$patient_condition.LN_condition,
                                                       levels=c(#"Meta_Pri.GC", "Neg_Pri.GC",
                                     "Meta_Neg.LN", "Neg_Neg.LN"))
p <- ggplot(tmp_neutro, 
            aes(x=patient_condition.LN_condition, y=Fraction)) + 
      geom_boxplot(aes(fill=patient_condition.LN_condition)) +
        scale_fill_manual(values = my_colors2_2) +
#       geom_signif(comparisons = list(
#                                      c("Meta_Neg.LN", "Neg_Neg.LN")
#                                      ), 
#                   map_signif_level=TRUE,y_position = 1.02) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
scale_y_continuous(limits = c(0, 1.1))+
facet_wrap(. ~  cell_type1, ncol = ncol) +
  theme_classic()+
  theme(axis.text.x=element_blank(),text = element_text(size=13))
plot(p)

tmp <- patient_condition.LN_condition_sample_ct_comp[["Neutro"]]
tmp <- tmp[tmp$LN_condition=="Pri.GC",]

options(repr.plot.width=3.5, repr.plot.height=3.5)

subtype2 = "MNDA+ Neutro"
tmp$patient_condition.LN_condition <- factor(tmp$patient_condition.LN_condition,
                                                       levels=c("Meta_Pri.GC", "Neg_Pri.GC"))
p <- ggplot(tmp[tmp$cell_type==subtype2,], aes(x=patient_condition.LN_condition, y=Fraction)) + 
      geom_boxplot(aes(fill=patient_condition.LN_condition)) +
    scale_fill_manual(values = my_colors2) +

      geom_signif(comparisons = list(c("Meta_Pri.GC", "Neg_Pri.GC")), 
                  map_signif_level=TRUE) + #,y_position = 0.12
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  ggtitle(subtype2) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))+
xlab("") +
labs(fill = NULL)
theme(legend.position="right")
plot(p)

ggsave(filename = "plots/figures/figS7C_MNDA_fraction.pdf", plot = p, width = 3.5, height = 3.5)



tmp <- patient_condition.LN_condition_sample_ct_comp[["Macro"]]
tmp <- tmp[tmp$LN_condition=="Neg.LN",]

options(repr.plot.width=3.5, repr.plot.height=3.5)

subtype2 = "FCN1+ Mph"
    tmp$patient_condition.LN_condition <- factor(tmp$patient_condition.LN_condition,
                                                           levels=c("Meta_Neg.LN", "Neg_Neg.LN"))
    p <- ggplot(tmp[tmp$cell_type==subtype2,], aes(x=patient_condition.LN_condition, y=Fraction)) + 
          geom_boxplot(aes(fill=patient_condition.LN_condition)) +
        scale_fill_manual(values = my_colors2) +
    
          geom_signif(comparisons = list(c("Meta_Neg.LN", "Neg_Neg.LN")), 
                      map_signif_level=TRUE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle(subtype2) +
      theme_classic()+
      theme(axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5))+
    xlab("") +
    labs(fill = NULL)
    theme(legend.position="right")
    plot(p)

ggsave(filename = "plots/figures/fig4A_FCN1_fraction.pdf", plot = p, width = 3.5, height = 3.5)




meta_Lauren.LN_condition <- meta[meta$Lauren.LN_condition %in% c('Diffuse_Met.LN','Diffuse_Pri.GC','Intestinal_Pri.GC','Intestinal_Met.LN'),]


# 对sample进行normalize后的结果
Lauren.LN_condition_sample_ct_comp <- list()
for(ct_i in c("all", "all_sub", 'T-CD8','T-CD4','B cell','NK','DC','Plasma','Macro','Tgd','Mono','Neutro')){
    if(ct_i == "all"){
        Lauren.LN_condition_sample_ct_comp[[ct_i]] <- table(meta_Lauren.LN_condition$cell_type1, meta_Lauren.LN_condition$sample)
    }else if(ct_i == "all_sub"){
        Lauren.LN_condition_sample_ct_comp[[ct_i]] <- table(meta_Lauren.LN_condition$subtype1, meta_Lauren.LN_condition$sample)
    }else{
        Lauren.LN_condition_sample_ct_comp[[ct_i]] <- table(meta_Lauren.LN_condition[meta_Lauren.LN_condition$cell_type1==ct_i,]$subtype1, meta_Lauren.LN_condition[meta_Lauren.LN_condition$cell_type1==ct_i,]$sample)
    }
    Lauren.LN_condition_sample_ct_comp[[ct_i]] <- melt(Lauren.LN_condition_sample_ct_comp[[ct_i]], id.vars = "sample")
    colnames(Lauren.LN_condition_sample_ct_comp[[ct_i]]) <- c("cell_type1", "sample", "Fraction")
    Lauren.LN_condition_sample_ct_comp[[ct_i]]$cell_type1 <- as.factor(Lauren.LN_condition_sample_ct_comp[[ct_i]]$cell_type1)
    
    # normalize到1
    for(condition in unique(Lauren.LN_condition_sample_ct_comp[[ct_i]]$sample)){
        mask=Lauren.LN_condition_sample_ct_comp[[ct_i]]$sample==condition
        Lauren.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask]=Lauren.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask]/sum(Lauren.LN_condition_sample_ct_comp[[ct_i]]$Fraction[mask])
    }
    
    Lauren.LN_condition_sample_ct_comp[[ct_i]]<-merge(Lauren.LN_condition_sample_ct_comp[[ct_i]], sample_meta, by = "sample")
    
}


#单独看几个亚型
options(repr.plot.width=2.8, repr.plot.height=4.2)
ct_i="DC"
# mask = Lauren.LN_condition_sample_ct_comp[[ct_i]]$"LN_condition"=="PBMC"
tmp <- Lauren.LN_condition_sample_ct_comp[[ct_i]]#[mask,]
dim(tmp)
subtype1 = "SPIB+ DC"

title=paste0(subtype1)#, " in ", ct_i
# reorder conditions in plots
tmp$Lauren.LN_condition <- factor(tmp$Lauren.LN_condition,levels=c('Intestinal_Pri.GC','Intestinal_Met.LN','Diffuse_Pri.GC','Diffuse_Met.LN'))
p <- ggplot(tmp[tmp$cell_type1==subtype1,], 
            aes(x=Lauren.LN_condition, y=Fraction)) + 
      geom_boxplot(aes(fill=Lauren.LN_condition)) +
scale_fill_manual(values = my_colors4) +
geom_signif(comparisons = list(c('Intestinal_Met.LN','Diffuse_Met.LN')),
        map_signif_level=TRUE,y_position = c(0.15)) +
#               geom_signif(comparisons = list(c("PBMC_Meta", "PBMC_Neg")), 
#                           map_signif_level=TRUE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  ggtitle(title) +
  theme_classic()+
  theme(text = element_text(size=13),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    strip.text.y =  element_text(size = 10),
        legend.position="none",
       plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(limits = c(0, 0.17))+
   ylab("Fraction")+
   xlab("") 
plot(p)
        
ggsave(filename = "plots/figures/figS9C_SPIB_fraction.pdf", plot = p, width = 2.8, height = 4.2)


meta_LN_condition.LN_station <- meta[meta$LN_condition.LN_station %in% c('Neg.LN_S2','Neg.LN_S1',
                                                                           'Met.LN_S1','Met.LN_S2'),]


# 产生做秩和检验的重要步骤，每次都需要重跑
# 对sample进行normalize后的结果，把所有细胞类型作为一个barplot作图
LN_condition.LN_station_sample_ct_comp <- list()
for(ct_i in c("all", "all_sub", 'T-CD8','T-CD4','B cell','NK','DC','Plasma','Macro','Tgd','Mono','Neutro')){
    if(ct_i == "all"){
        LN_condition.LN_station_sample_ct_comp[[ct_i]] <- table(meta_LN_condition.LN_station$cell_type1, meta_LN_condition.LN_station$sample)
    }else if(ct_i == "all_sub"){
        LN_condition.LN_station_sample_ct_comp[[ct_i]] <- table(meta_LN_condition.LN_station$subtype1, meta_LN_condition.LN_station$sample)
    }else{
        LN_condition.LN_station_sample_ct_comp[[ct_i]] <- table(meta_LN_condition.LN_station[meta_LN_condition.LN_station$cell_type1==ct_i,]$subtype1, meta_LN_condition.LN_station[meta_LN_condition.LN_station$cell_type1==ct_i,]$sample)
    }
    LN_condition.LN_station_sample_ct_comp[[ct_i]] <- melt(LN_condition.LN_station_sample_ct_comp[[ct_i]], id.vars = "sample")
    colnames(LN_condition.LN_station_sample_ct_comp[[ct_i]]) <- c("cell_type1", "sample", "Fraction")
    LN_condition.LN_station_sample_ct_comp[[ct_i]]$cell_type1 <- as.factor(LN_condition.LN_station_sample_ct_comp[[ct_i]]$cell_type1)
    
    # normalize到1
    for(condition in unique(LN_condition.LN_station_sample_ct_comp[[ct_i]]$sample)){
        mask=LN_condition.LN_station_sample_ct_comp[[ct_i]]$sample==condition
        LN_condition.LN_station_sample_ct_comp[[ct_i]]$Fraction[mask]=LN_condition.LN_station_sample_ct_comp[[ct_i]]$Fraction[mask]/sum(LN_condition.LN_station_sample_ct_comp[[ct_i]]$Fraction[mask])
    }
    
    LN_condition.LN_station_sample_ct_comp[[ct_i]]<-merge(LN_condition.LN_station_sample_ct_comp[[ct_i]], sample_meta, by = "sample")
 
    
}


#单独看几个亚型
options(repr.plot.width=3.5, repr.plot.height=3.5)

tmp <- LN_condition.LN_station_sample_ct_comp[["DC"]]
tmp <- tmp[tmp$LN_condition=="Met.LN",]

for(subtype1 in c("FCER1A+ cDC2")){
    title=paste0(subtype1)
    tmp$LN_condition.LN_station <- factor(tmp$LN_condition.LN_station,levels=c("Met.LN_S1",
                                                                    "Met.LN_S2"))
    p <- ggplot(tmp[tmp$cell_type1==subtype1,], 
                aes(x=LN_condition.LN_station, y=Fraction)) + 
                  geom_boxplot(aes(fill=LN_condition.LN_station)) +
        scale_fill_manual(values = my_colors2) +
    
            geom_signif(comparisons = list(c("Met.LN_S1", "Met.LN_S2")), 
                                      map_signif_level=TRUE,y_position = c(0.04)) +
                  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
                  ggtitle(title) +
            theme_classic()+
              theme(axis.text.x = element_blank(),
                    plot.title = element_text(hjust = 0.5))+
            xlab("") +
            labs(fill = NULL)+
            theme(legend.position="right")
    plot(p)
    print(subtype1)
        
}
ggsave(filename = "plots/figures/figS9E_cDC2_fraction.pdf", plot = p, width = 3.5, height = 3.5)

# library(openxlsx)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggsignif)

frac_TIM3_by_GSM <- read.csv("eACA_scan/Fraction_df/frac_TIM3_by_CD8.csv",row.names = 1)

tmp <- frac_TIM3_by_GSM[frac_TIM3_by_GSM$patient.condition=="Cancer",]
tmp$sample.condition[tmp$sample.condition=="tissue"]="Normal"
tmp$sample.condition[tmp$sample.condition=="lymph.node"]="Lymph node"
tmp$sample.condition[tmp$sample.condition=="LN_meta"]="Lymph node" #不分是否为转移淋巴结了
tmp$sample.condition[tmp$sample.condition=="tumor"]="Tumor"
tmp$sample.condition[tmp$sample.condition=="pbmc"]="PBMC"


table(tmp$sample.condition)

tmp_df <- aggregate(tmp$Fraction,list(tmp$sample.condition),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

tmp$sample.condition <- factor(tmp$sample.condition,levels=c("Normal", "Tumor",
                                                       "Lymph node", "PBMC"))

my_colors4 = my_colors7[c(1,3,2,4)]
options(repr.plot.width=5, repr.plot.height=4)
fig2f <- ggplot(tmp, aes(x=sample.condition, y=Fraction)) + 
    geom_boxplot(aes(fill=sample.condition)) +
          scale_fill_manual(values = my_colors4) +

    geom_signif(comparisons = list(
                                  c("Lymph node", "Tumor")), 
                map_signif_level=TRUE,
                y_position = c(0.4)) +
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
    theme_classic()+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    strip.text.y =  element_text(size = 10),
#     legend.position="none",
    plot.title = element_text(hjust = 0.5))+
    ggtitle("CD8_TIM3+ Trm") + #
   ylab("Fraction among CD8+ Tcells")+ #
   xlab("") 
plot(fig2f)
ggsave("plots/figures/fig2F_CWAS.pdf",fig2f,width=5, height=4)
# +coord_flip()

frac_HLA_by_GSM <- read.csv("eACA_scan/Fraction_df/frac_HLA_by_CD8.csv",row.names = 1)

tmp <- frac_HLA_by_GSM[frac_HLA_by_GSM$GSM.condition=="Cancer-tumor",]

as.data.frame(table(tmp$curated_cancer_type))

tmp_df <- aggregate(tmp$Fraction,list(tmp$curated_cancer_type),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

# 统计 curated_cancer_type 出现次数
type_counts <- tmp %>%
  count(curated_cancer_type) %>%
  filter(n > 3)  # 筛选出现次数大于 3 的类型

# 按照中位数从大到小对 curated_cancer_type 进行排序
sorted_types <- tmp %>%
  filter(!is.na(curated_cancer_type)) %>%  # 删除 NA 值
  filter(curated_cancer_type %in% type_counts$curated_cancer_type) %>%
  group_by(curated_cancer_type) %>%
  summarise(median_fraction = median(Fraction)) %>%
  arrange(desc(median_fraction)) %>%
  pull(curated_cancer_type)

# 使用 ggplot2 绘制箱线图
options(repr.plot.width=6, repr.plot.height=3.5)
ggplot(tmp %>% filter(curated_cancer_type %in% sorted_types), aes(x = factor(curated_cancer_type, levels = sorted_types), y = Fraction)) +
  geom_boxplot(aes(fill=curated_cancer_type)) +
  labs(x = "Cancer type", y = "Fraction") +  # 修改 x 轴标题和 y 轴范围
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
        legend.position="none",
        plot.title = element_text(hjust = 0.5)) +
    ggtitle("CD8_HLA-DRA+ Tem") + #
  ylim(0, 0.2)  # 修改 y 轴范围
ggsave("plots/figures/fig3F_CWAS.pdf",width=6, height=3.5)


tmp <- frac_HLA_by_GSM[frac_HLA_by_GSM$patient.condition=="Cancer",]
tmp$sample.condition[tmp$sample.condition=="tissue"]="Normal"
tmp$sample.condition[tmp$sample.condition=="lymph.node"]="Lymph node"
tmp$sample.condition[tmp$sample.condition=="LN_meta"]="Lymph node" #不分是否为转移淋巴结了
tmp$sample.condition[tmp$sample.condition=="tumor"]="Tumor"
tmp$sample.condition[tmp$sample.condition=="pbmc"]="PBMC"
# tmp$Fraction <- tmp$Fraction * 100#


table(tmp$sample.condition)

tmp_df <- aggregate(tmp$Fraction,list(tmp$sample.condition),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

tmp$sample.condition <- factor(tmp$sample.condition,levels=c("Normal", "Tumor",
                                                       "Lymph node", "PBMC"))

options(repr.plot.width=5, repr.plot.height=4)
ggplot(tmp, aes(x=sample.condition, y=Fraction)) + 
    geom_boxplot(aes(fill=sample.condition)) +
          scale_fill_manual(values = my_colors4) +
    geom_signif(comparisons = list(c("Normal", "Tumor"),
#                                    c("Tumor","LN_meta"),
#                                    c("LN_meta","Lymph node"),
                                   c("Lymph node", "PBMC"),
                                  c("Tumor", "PBMC"),
                                  c("Normal", "PBMC"),
#                                   c("Lymph node", "Normal"),
                                  c("Lymph node", "Tumor")), 
                map_signif_level=TRUE,
                y_position = c(0.14,0.16,0.23,0.26,0.2)) +
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
    theme_classic()+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    strip.text.y =  element_text(size = 10),
#     legend.position="none",
    plot.title = element_text(hjust = 0.5))+
   ylab("Fraction")+
   xlab("") + ylim(0, 0.25) +
ggtitle("CD8_HLA-DRA+ Tem")
# +coord_flip()
ggsave("plots/figures/figS5D_CWAS.pdf",width=5, height=4)




frac_GZMB_by_GSM <- read.csv("eACA_scan/Fraction_df/frac_GZMB_by_CD8.csv",row.names = 1)

tmp <- frac_GZMB_by_GSM[frac_GZMB_by_GSM$GSM.condition=="Cancer-pbmc",]
# tmp$Fraction <- tmp$Fraction * 100#

as.data.frame(table(tmp$curated_cancer_type))

tmp_df <- aggregate(tmp$Fraction,list(tmp$curated_cancer_type),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

# 统计 curated_cancer_type 出现次数
type_counts <- tmp %>%
  count(curated_cancer_type) %>%
  filter(n > 3)  # 筛选出现次数大于 3 的类型

# 按照中位数从大到小对 curated_cancer_type 进行排序
sorted_types <- tmp %>%
  filter(!is.na(curated_cancer_type)) %>%  # 删除 NA 值
  filter(curated_cancer_type %in% type_counts$curated_cancer_type) %>%
  group_by(curated_cancer_type) %>%
  summarise(median_fraction = median(Fraction)) %>%
  arrange(desc(median_fraction)) %>%
  pull(curated_cancer_type)

# 使用 ggplot2 绘制箱线图
options(repr.plot.width=5, repr.plot.height=3.5)
ggplot(tmp %>% filter(curated_cancer_type %in% sorted_types), aes(x = factor(curated_cancer_type, levels = sorted_types), y = Fraction)) +
  geom_boxplot(aes(fill=curated_cancer_type)) +
          scale_fill_manual(values = my_colors5) +

  labs(x = "Cancer type", y = "Fraction") +  # 修改 x 轴标题和 y 轴范围
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
        legend.position="none",
        plot.title = element_text(hjust = 0.5)) +
    ggtitle("CD8_GZMB+ Temra")
#   ylim(0, 0.2)  # 修改 y 轴范围
ggsave("plots/figures/figS6I_CWAS.pdf",width=5, height=3.5)




tmp <- frac_GZMB_by_GSM[frac_GZMB_by_GSM$GSM.condition %in% c("Cancer-pbmc","Health-pbmc"),]
# tmp$Fraction <- tmp$Fraction * 100#

tmp$curated_cancer_type[tmp$patient.condition=="Health"]="Health"
as.data.frame(table(tmp$curated_cancer_type))

tmp_df <- aggregate(tmp$Fraction,list(tmp$curated_cancer_type),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

# 统计 curated_cancer_type 出现次数
type_counts <- tmp %>%
  count(curated_cancer_type) %>%
  filter(n > 3)  # 筛选出现次数大于 3 的类型

# 按照中位数从大到小对 curated_cancer_type 进行排序
sorted_types <- tmp %>%
  filter(!is.na(curated_cancer_type)) %>%  # 删除 NA 值
  filter(curated_cancer_type %in% type_counts$curated_cancer_type) %>%
  group_by(curated_cancer_type) %>%
  summarise(median_fraction = median(Fraction)) %>%
  arrange(desc(median_fraction)) %>%
  pull(curated_cancer_type)

# 使用 ggplot2 绘制箱线图
options(repr.plot.width=6, repr.plot.height=5)
ggplot(tmp %>% filter(curated_cancer_type %in% sorted_types), aes(x = factor(curated_cancer_type, levels = sorted_types), y = Fraction)) +
  geom_boxplot(aes(fill=curated_cancer_type)) +
  labs(x = "Cancer type", y = "Fraction") +  # 修改 x 轴标题和 y 轴范围
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.text.y =  element_text(size = 10),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))# +
#   ylim(0, 0.2)  # 修改 y 轴范围

tmp <- frac_GZMB_by_GSM[frac_GZMB_by_GSM$sample.condition=="pbmc",]
# tmp$Fraction <- tmp$Fraction * 100#

table(tmp$patient.condition)

tmp_df <- aggregate(tmp$Fraction,list(tmp$patient.condition),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

tmp <- tmp[tmp$patient.condition %in% c("Cancer", "Health"),]
tmp$patient.condition <- factor(tmp$patient.condition,levels=c("Cancer", "Health"))

options(repr.plot.width=2.5, repr.plot.height=2.5)
ggplot(tmp, aes(x=patient.condition, y=Fraction)) + 
    geom_boxplot(aes(fill=patient.condition)) +
          scale_fill_manual(values = my_colors2) +

    geom_signif(comparisons = list(
                                   c("Cancer", "Health")), 
                map_signif_level=TRUE,
                y_position = c(0.7)) +
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
    theme_classic()+
  theme(axis.text.x = element_text(size=15), 
#         axis.text.y = element_text(size=15),
    strip.text.y =  element_text(size = 20),
    legend.position="none",
    plot.title = element_text(hjust = 0.5))+
   ylab("Fraction")+
   xlab("") 
# +coord_flip()
ggsave("plots/figures/figS6H_CWAS.pdf",width=2.5, height=2.5)


tmp <- frac_GZMB_by_GSM[frac_GZMB_by_GSM$patient.condition=="Cancer",]
tmp$sample.condition[tmp$sample.condition=="tissue"]="Normal"
tmp$sample.condition[tmp$sample.condition=="lymph.node"]="Lymph node"
tmp$sample.condition[tmp$sample.condition=="LN_meta"]="Lymph node" #不分是否为转移淋巴结了
tmp$sample.condition[tmp$sample.condition=="tumor"]="Tumor"
tmp$sample.condition[tmp$sample.condition=="pbmc"]="PBMC"
# tmp$Fraction <- tmp$Fraction * 100#

table(tmp$sample.condition)

tmp_df <- aggregate(tmp$Fraction,list(tmp$sample.condition),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

tmp$sample.condition <- factor(tmp$sample.condition,levels=c("Normal", "Tumor",
                                                       "Lymph node", "PBMC"))

options(repr.plot.width=5, repr.plot.height=3.5)
ggplot(tmp, aes(x=sample.condition, y=Fraction)) + 
    geom_boxplot(aes(fill=sample.condition)) +
          scale_fill_manual(values = my_colors4) +

    geom_signif(comparisons = list(#c("Normal", "Tumor"),
#                                    c("Tumor","LN_meta"),
#                                    c("LN_meta","Lymph node"),
                                   c("Lymph node", "PBMC"),
                                  c("Tumor", "PBMC"),
                                  c("Normal", "PBMC")), 
                map_signif_level=TRUE,
                y_position = c(0.6,0.7,0.8)) +
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
    theme_classic()+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    strip.text.y =  element_text(size = 10),
#     legend.position="none",
    plot.title = element_text(hjust = 0.5))+
   ylab("Fraction")+
   xlab("") 
# +coord_flip()
ggsave("plots/figures/figS6G_CWAS.pdf",width=5, height=3.5)


frac_CD4_Tex_by_GSM <- read.csv("eACA_scan/Fraction_df/frac_CD4_Tex_by_CD4.csv",row.names = 1)

tmp <- frac_CD4_Tex_by_GSM[frac_CD4_Tex_by_GSM$patient.condition=="Cancer",]
tmp$sample.condition[tmp$sample.condition=="tissue"]="Normal"
tmp$sample.condition[tmp$sample.condition=="lymph.node"]="LN_neg"
tmp$sample.condition[tmp$sample.condition=="LN_meta"]="LN_meta"
tmp$sample.condition[tmp$sample.condition=="tumor"]="Tumor"
tmp$sample.condition[tmp$sample.condition=="pbmc"]="PBMC"
# tmp$Fraction <- tmp$Fraction * 100#

table(tmp$sample.condition)
aggregate(tmp$Fraction,list(tmp$sample.condition),mean)

tmp$sample.condition <- factor(tmp$sample.condition,levels=c("Normal", "Tumor",
                                                       "LN_meta","LN_neg", "PBMC"))

options(repr.plot.width=4, repr.plot.height=3.5)
ggplot(tmp, aes(x=sample.condition, y=Fraction)) + 
    geom_boxplot(aes(fill=sample.condition)) +
          scale_fill_manual(values = my_colors5) +
    geom_signif(comparisons = list(
                                   c("Tumor","LN_meta"),
                                   c("LN_meta","LN_neg")), 
                map_signif_level=TRUE,
                y_position = c(0.3,0.3)) +
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
    theme_classic()+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    strip.text.y =  element_text(size = 10),
    legend.position="none",
    plot.title = element_text(hjust = 0.5))+
   ylab("Fraction")+
   xlab("") 
# +coord_flip()
ggsave("plots/figures/figS3E_CWAS.pdf",width=4, height=3.5)




frac_SPIB_by_GSM <- read.csv("eACA_scan/Fraction_df/frac_SPIB_by_DC.csv",row.names = 1)

tmp <- frac_SPIB_by_GSM[frac_SPIB_by_GSM$patient.condition=="Cancer",]
tmp$sample.condition[tmp$sample.condition=="tissue"]="Normal"
tmp$sample.condition[tmp$sample.condition=="lymph.node"]="Lymph node"
tmp$sample.condition[tmp$sample.condition=="LN_meta"]="Lymph node" #不分是否为转移淋巴结了
tmp$sample.condition[tmp$sample.condition=="tumor"]="Tumor"
tmp$sample.condition[tmp$sample.condition=="pbmc"]="PBMC"
# tmp$Fraction <- tmp$Fraction * 100#

table(tmp$sample.condition)

tmp_df <- aggregate(tmp$Fraction,list(tmp$sample.condition),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

tmp$sample.condition <- factor(tmp$sample.condition,levels=c("Normal", "Tumor",
                                                       "Lymph node", "PBMC"))

options(repr.plot.width=5.5, repr.plot.height=4)
ggplot(tmp, aes(x=sample.condition, y=Fraction)) + 
    geom_boxplot(aes(fill=sample.condition)) +
          scale_fill_manual(values = my_colors4) +
    geom_signif(comparisons = list(c("Normal", "Tumor")), 
                map_signif_level=TRUE,
                y_position = c(0.3)) +
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
    theme_classic()+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    strip.text.y =  element_text(size = 10),
#     legend.position="none",
    plot.title = element_text(hjust = 0.5))+
   ylab("Fraction")+
   xlab("") +scale_y_continuous(limits = c(0, 0.5))+
    ggtitle("SPIB+ DC")
# +coord_flip()
ggsave("plots/figures/figS9D_CWAS.pdf",width=5.5, height=4)




frac_cDC2_by_GSM <- read.csv("eACA_scan/Fraction_df/frac_cDC2_by_DC.csv",row.names = 1)

tmp <- frac_cDC2_by_GSM[frac_cDC2_by_GSM$sample.condition=="pbmc" & 
                        frac_cDC2_by_GSM$patient.condition=="Cancer",]
# tmp$Fraction <- tmp$Fraction * 100#

table(tmp$Subtype)

tmp_df <- aggregate(tmp$Fraction,list(tmp$Subtype),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

# tmp$Subtype <- factor(tmp$Subtype,levels=c("disease", "tumor", "health"))

options(repr.plot.width=3, repr.plot.height=4)
ggplot(tmp, aes(x=Subtype, y=Fraction)) + 
    geom_boxplot(aes(fill=Subtype)) +
          scale_fill_manual(values = my_colors2) +
    geom_signif(comparisons = list(c("FCER1A+ cDC2","NLRP3+ cDC2")), 
                map_signif_level=TRUE,
                y_position = c(0.9)) +
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
    theme_classic()+
#   theme(axis.text.x = element_text(size=15), 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 12),
        
#         axis.text.y = element_text(size=15),
    strip.text.y =  element_text(size = 20),
    legend.position="none",
    plot.title = element_text(hjust = 0.5))+
   ylab("Fraction")+
   xlab("") +
ggtitle("PBMC in cancer patients")
# +coord_flip()
ggsave("plots/figures/figS9G_CWAS.pdf",width=3, height=4)




tmp <- frac_cDC2_by_GSM[frac_cDC2_by_GSM$sample.condition=="pbmc" & 
                        frac_cDC2_by_GSM$patient.condition=="Health",]
# tmp$Fraction <- tmp$Fraction * 100#

table(tmp$Subtype)

tmp_df <- aggregate(tmp$Fraction,list(tmp$Subtype),mean)
tmp_df$x <- round(tmp_df$x,4)
tmp_df

# tmp$Subtype <- factor(tmp$Subtype,levels=c("disease", "tumor", "health"))

options(repr.plot.width=3, repr.plot.height=4)
ggplot(tmp, aes(x=Subtype, y=Fraction)) + 
    geom_boxplot(aes(fill=Subtype)) +
          scale_fill_manual(values = my_colors2) +
    geom_signif(comparisons = list(c("FCER1A+ cDC2","NLRP3+ cDC2")), 
                map_signif_level=TRUE,
                y_position = c(0.9)) +
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
    theme_classic()+
#   theme(axis.text.x = element_text(size=15), 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 12),
        
#         axis.text.y = element_text(size=15),
    strip.text.y =  element_text(size = 20),
    legend.position="none",
    plot.title = element_text(hjust = 0.5))+
   ylab("Fraction")+
   xlab("") +
ggtitle("PBMC in healthy donors")

# +coord_flip()
ggsave("plots/figures/figS9H_CWAS.pdf",width=3, height=4)


tcell.batch.corrected <- readRDS("plots/data/tcell_no_TLB_new_meta.rds")
Idents(tcell.batch.corrected) <- "subtype1"
DefaultAssay(tcell.batch.corrected)

options(repr.plot.width=12, repr.plot.height=7)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "subtype1",label=FALSE)

pdf("plots/figures/fig2A_umap.pdf",width=12, height=7)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "subtype1",label=FALSE)
dev.off()

# try color
options(repr.plot.width=12, repr.plot.height=7)

DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "subtype1",label=FALSE,
         cols = colorRampPalette(colors = c("#e2ad9f","#d37562","#d7c6de","#9b7aad",
                                            "#90aed2","#c0d6eb"))(length(unique(tcell.batch.corrected$subtype1))))
        # cols = DiscretePalette(length(unique(tcell.batch.corrected$subtype1)), palette = "parade", shuffle = FALSE))

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE)
pdf("plots/figures/figS2A_umap.pdf",width=6, height=6)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE)
dev.off()

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE,
         cols = colorRampPalette(colors = c("#e2ad9f","#d37562","#d7c6de","#9b7aad",
                                            "#90aed2","#c0d6eb"))(length(unique(tcell.batch.corrected$cell_type1))))

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE,
         cols = colorRampPalette(colors = my_colors5)(length(unique(tcell.batch.corrected$cell_type1))))

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE,
       cols="Set3")

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE,
       cols="Set1")

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE,
       cols="Paired")

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(tcell.batch.corrected, reduction = "umap", shuffle=TRUE, group.by = "cell_type1",label=FALSE,
       cols="Dark2")

tex <- subset(tcell.batch.corrected, subtype1 %in% c("CD8_TIM3+ Trm","CD8_LAYN+ Tex","CD8_CCL5+ Tex"))

Idents(tex) <- "subtype1"
Idents(tex) <- factor(Idents(tex), levels = c("CD8_LAYN+ Tex","CD8_CCL5+ Tex","CD8_TIM3+ Trm"))
unique(Idents(tex))

options(repr.plot.width=8, repr.plot.height=4)

vln_features=c("IL7R", "CCL5","CD44","GZMK","CD69",#canonical TMEM-associated markers
"TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN','TNFRSF9','CXCL13','IFNG')

p <- VlnPlot(tex, features = vln_features,
        stack=TRUE)+
    scale_fill_manual(values = colorRampPalette(colors = my_colors5)(length(vln_features))) +
  RotatedAxis()+ NoLegend()
plot(p)

pdf("plots/figures/fig2B_expression.pdf",width=8, height=4)
plot(p)
dev.off()

options(repr.plot.width=8, repr.plot.height=4)
VlnPlot(tex, features = c("IL7R", "CCL5","CD44","GZMK","CD69",#canonical TMEM-associated markers
"TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN','TNFRSF9','CXCL13','IFNG'),
        stack=TRUE)+
  RotatedAxis()+ NoLegend()

options(repr.plot.width=8, repr.plot.height=4)
vln_features=c("IL7R", "CCL5","CD44","GZMK","CD69",#canonical TMEM-associated markers
"TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN','TNFRSF9','CXCL13','IFNG')
VlnPlot(tex, features = vln_features,
        stack=TRUE)+
    scale_fill_manual(values = colorRampPalette(colors = my_colors7)(length(vln_features))) +
  RotatedAxis()+ NoLegend()

ESCA_CD8 <- read.table("public_LN/pancan_Tcell/GSE156728_ESCA_10X.CD8.counts.txt.gz")
esca_cd8 <- CreateSeuratObject(ESCA_CD8)

meta_pancan <- read.table("public_LN/pancan_Tcell/GSE156728_metadata.txt.gz",header = TRUE,row.names = 1)
meta_esca <- meta_pancan[colnames(esca_cd8),]

esca_cd8 <- AddMetaData(esca_cd8, meta_esca, col.name = NULL)

esca_cd8 <- NormalizeData(esca_cd8, normalization.method = "LogNormalize", scale.factor = 10000)

Idents(esca_cd8) <- "meta.cluster"
saveRDS(esca_cd8,"public_LN/pancan_Tcell/ESCA_CD8.rds")

esca_cd8 <- readRDS("public_LN/pancan_Tcell/ESCA_CD8.rds")
dim(esca_cd8)
unique(Idents(esca_cd8))


options(repr.plot.width=8, repr.plot.height=8)

vln_features=c("IL7R", "CCL5","CD44","GZMK","CD69",#canonical TMEM-associated markers
"TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN','TNFRSF9','CXCL13','IFNG')

VlnPlot(esca_cd8, features = vln_features,
        stack=TRUE)+
    scale_fill_manual(values = colorRampPalette(colors = my_colors5)(length(vln_features))) +
  RotatedAxis()+ NoLegend()

pdf("plots/figures/figS2B_expression.pdf",width=8, height=8)
VlnPlot(esca_cd8, features = vln_features,
        stack=TRUE)+
    scale_fill_manual(values = colorRampPalette(colors = my_colors5)(length(vln_features))) +
  RotatedAxis()+ NoLegend()
dev.off()



cd4_tcell.batch.corrected <- readRDS("plots/data/cd4_tcell_new_meta.rds")


cd4_tex <- subset(cd4_tcell.batch.corrected, subtype1=="CD4_Tex")
cd4_tex <- subset(cd4_tex, LN_condition!="PBMC")

options(repr.plot.width=6, repr.plot.height=4)
cd4_features <- c("TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN')

p <- VlnPlot(cd4_tex, features = cd4_features,
        group.by="LN_condition",
        stack=TRUE)+
    scale_fill_manual(values = colorRampPalette(colors = my_colors7)(length(cd4_features))) +
  RotatedAxis()+ NoLegend()
plot(p)
pdf("plots/figures/figS3B_expression.pdf",width=6, height=4)
plot(p)
dev.off()

options(repr.plot.width=6, repr.plot.height=4)
cd4_features <- c("TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN',"FOXP3")

VlnPlot(cd4_tex, features = cd4_features,
        group.by="LN_condition",
        stack=TRUE)+
    scale_fill_manual(values = colorRampPalette(colors = my_colors5)(length(cd4_features))) +
  RotatedAxis()+ NoLegend()

options(repr.plot.width=6, repr.plot.height=7)
cd4_features <- c("TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN',"FOXP3")

VlnPlot(cd4_tcell.batch.corrected, features = cd4_features,
        group.by="subtype1",
        stack=TRUE)+
    scale_fill_manual(values = colorRampPalette(colors = my_colors5)(length(cd4_features))) +
  RotatedAxis()+ NoLegend()

options(repr.plot.width=6, repr.plot.height=4)
VlnPlot(cd4_tex, features = c(#"IL7R", "CCL5","CD44","GZMK","CD69",#canonical TMEM-associated markers
        "TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN'),#'TNFRSF9','CXCL13','IFNG'),#idents="CD4+ Tex",
        group.by="LN_condition",
        stack=TRUE)+
  RotatedAxis()+ NoLegend()



NK <- subset(tcell.batch.corrected, subtype1 %in% c('NK1','NK2','NK3'))
Idents(NK) <- "subtype1"
Idents(NK) <- factor(Idents(NK), levels = c('NK1','NK2','NK3'))

options(repr.plot.width=10, repr.plot.height=4)
DotPlot(NK, features = c("FGFBP2", "GZMB", "GZMH", "PRF1","FCGR3A","S1PR5",
                        "CCL3", "CCL4", "XCL1", "XCL2", "GZMK","CCL3L1","AREG",
                         "LTB","SELL","CD44"
                        )
       )+RotatedAxis()+
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

pdf("plots/figures/figS4A_expression.pdf",width=10, height=4)
DotPlot(NK, features = c("FGFBP2", "GZMB", "GZMH", "PRF1","FCGR3A","S1PR5",
                        "CCL3", "CCL4", "XCL1", "XCL2", "GZMK","CCL3L1","AREG",
                         "LTB","SELL","CD44"
                        )
       )+
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")
dev.off()

options(repr.plot.width=9, repr.plot.height=3)
VlnPlot(NK, features = c("KLRB1", # NK cell CD161,i.e.NK1.1
                                                                 "NCAM1",#CD56
                                                                 "FCGR3A"#CD16
                        ),ncol=3,pt.size=0#,cols=c("#d7c6de", "#c0d6eb","#e2ad9f")
       )
pdf("plots/figures/fig2H_expression.pdf",width=9, height=3)
VlnPlot(NK, features = c("KLRB1", # NK cell CD161,i.e.NK1.1
                                                                 "NCAM1",#CD56
                                                                 "FCGR3A"#CD16
                        ),ncol=3,pt.size=0#,cols=c("#d7c6de", "#c0d6eb","#e2ad9f")
       )
dev.off()

options(repr.plot.width=9, repr.plot.height=3)
VlnPlot(NK, features = c("KLRB1", # NK cell CD161,i.e.NK1.1
                                                                 "NCAM1",#CD56
                                                                 "FCGR3A"#CD16
                        ),ncol=3,pt.size=0,cols=c("#d7c6de", "#c0d6eb","#e2ad9f")
       )


Idents(tcell.batch.corrected) <- "subtype1"

hla_tex <- subset(tcell.batch.corrected, subtype1 %in% c('CD8_MKI67+ Tex','CD8_LAYN+ Tex','CD8_CCL5+ Tex',"CD8_HLA-DRA+ Tem"))

options(repr.plot.width=8, repr.plot.height=4)
hla_features <- c("GZMK","GZMB",'IFNG',"PRF1",#canonical TMEM-associated markers
"TCF7","TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN')
p <- VlnPlot(hla_tex, features = hla_features,
        # idents=unique(Idents(tcell.batch.corrected)[tcell.batch.corrected$cell_type1=="T-CD8"]),
        stack=TRUE)+
    scale_fill_manual(values = colorRampPalette(colors = my_colors7)(length(hla_features))) +
  RotatedAxis()+ NoLegend()
plot(p)
pdf("plots/figures/fig3B_expression.pdf",width=8, height=4)
plot(p)
dev.off()

options(repr.plot.width=12, repr.plot.height=12)

VlnPlot(tcell.batch.corrected, features = c("CD3D","CD3E","CD8A","CD8B","HLA-DRA","HLA-DRB5","HAVCR2",
                                            "TCF7","TOX","PDCD1", "LAG3","CTLA4","TIGIT",'LAYN'),
        # idents=unique(Idents(tcell.batch.corrected)[tcell.batch.corrected$cell_type1=="T-CD8"]),
        stack=TRUE)

library(openxlsx)
surface_protein=read.xlsx("../GENOME/human_surface_protein/2018_PNAS_table_S3_surfaceome.xlsx", 
                          sheet = 2, startRow = 2)
surface_gene=surface_protein$UniProt.gene

c("CD3D","CD8A","HLA-DRA","HLA-DRB5","HAVCR2","CTLA4") %in% surface_gene

options(repr.plot.width=8, repr.plot.height=4)
VlnPlot(tcell.batch.corrected, idents = c("CD8_GZMB+ Temra","CD8_GZMK+ Temra"),
        features = c("FGFBP2","FCGR3A"),ncol=2,
        pt.size=0.00001#,stack=TRUE
       )
pdf("plots/figures/fig3H_expression.pdf",width=8, height=4)
VlnPlot(tcell.batch.corrected, idents = c("CD8_GZMB+ Temra","CD8_GZMK+ Temra"),
        features = c("FGFBP2","FCGR3A"),ncol=2,
        pt.size=0.00001#,stack=TRUE
       )
dev.off()

options(repr.plot.width=4, repr.plot.height=8)
VlnPlot(tcell.batch.corrected, idents = c("CD8_GZMB+ Temra","CD8_GZMK+ Temra"),
        features = c("GZMK", "GZMB"),ncol=1,
        pt.size=0.00001#,stack=TRUE
       )
pdf("plots/figures/figS6A_expression.pdf",width=4, height=8)
VlnPlot(tcell.batch.corrected, idents = c("CD8_GZMB+ Temra","CD8_GZMK+ Temra"),
        features = c("GZMK", "GZMB"),ncol=1,
        pt.size=0.00001#,stack=TRUE
       )
dev.off()

FunctionScoreMatrix <- readRDS("processed_data/data_B2-19/cell_typing/tcell/cd8_FunctionScoreMatrix.rds")
signatureType_row <- readRDS("processed_data/data_B2-19/cell_typing/tcell/cd8_signatureType_row.rds")


# orderC = c("CD8_c13", "CD8_c3", "CD8_c6", "CD8_c0", "CD8_c11", "CD8_c9", "CD8_c10", "CD8_c12", "CD8_c8", "CD8_c2", "CD8_c7", "CD8_c4", "CD8_c5", "CD8_c1")
# FunctionScoreMatrix <- FunctionScoreMatrix[,orderC]
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
    colorRampPalette(colors = c("#89afd2", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#ef8775"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)

options(repr.plot.width=8, repr.plot.height=5)

p <- pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(3, 14, 17),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8) #,
#          filename = file.path(figurePath, paste0("CD8_FunctionScore_heatmap.pdf")))
p
pdf("plots/figures/figS6D_heatmap.pdf", width = 8,height = 5)
p
dev.off()

bcell.batch.corrected <- readRDS("plots/data/bcell_no_TLB_new_meta.rds")
DefaultAssay(bcell.batch.corrected) <- "RNA"

options(repr.plot.width=11.5, repr.plot.height=8.5)
p <- DotPlot(bcell.batch.corrected, group.by="subtype1",cluster.idents = TRUE,
        features = c("BCL6", "AICDA","CD38","MS4A1","MKI67","JCHAIN","IGHA1","IGHG1","IGHM","IGHD")
       )+
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")
plot(p)

pdf("plots/figures/fig5A_expression.pdf",width=11.5, height=8.5)
plot(p)
dev.off()

options(repr.plot.width=12, repr.plot.height=9)
DotPlot(bcell.batch.corrected, group.by="subtype1",cluster.idents = TRUE,
        features = c("BCL6", "AICDA","CD38","MS4A1","MKI67","JCHAIN","IGHA1","IGHG1","IGHM","IGHD" ))#Germinal center B cell, TLS

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(bcell.batch.corrected, reduction = "umap", group.by = "subtype1",shuffle=TRUE, label = FALSE)
pdf("plots/figures/figS10A_umap.pdf",width=8, height=6)
DimPlot(bcell.batch.corrected, reduction = "umap", group.by = "subtype1",shuffle=TRUE, label = FALSE)
dev.off()

DimPlot(bcell.batch.corrected)

myeloid.batch.corrected <- readRDS("plots/data/myeloid_new_meta.rds")
Idents(myeloid.batch.corrected) <- "subtype1"

options(repr.plot.width=10, repr.plot.height=6)
DimPlot(myeloid.batch.corrected, reduction = "umap", group.by = "subtype1",shuffle=TRUE, label = FALSE)
pdf("plots/figures/figS7A_umap.pdf",width=10, height=6)
DimPlot(myeloid.batch.corrected, reduction = "umap", group.by = "subtype1",shuffle=TRUE, label = FALSE)
dev.off()

DC_so <- subset(myeloid.batch.corrected, 
                    subtype1 %in% grep("DC",unique(myeloid.batch.corrected$subtype1),value = TRUE))

Idents(DC_so) <- "subtype1"

unique(Idents(DC_so))

# fig
# myeloid subtype markers from cancer cell small cell lung cancer
DC_markers=c("HLA-DRA", "HLA-DRB1", # antigen presentation
                "CD1C",# cDC2
             "CLEC9A","CADM1","XCR1",#cDC1
                  "LAMP3","CD274", # LAMP3+ DC, PD-L1
                  "PPA1","LSP1","CSF2RA","HLA-DQB1","ID2",#cDC
                "GZMB","JCHAIN","IRF7","ITM2C","LILRA4","PLD4","TCF4"#pDC
               )


options(repr.plot.width=10, repr.plot.height=5)
DotPlot(DC_so, features = DC_markers, cluster.idents = TRUE)+
  RotatedAxis()+
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

pdf(file = "plots/figures/fig4E_dc_expression.pdf", width = 10, height = 5)
DotPlot(DC_so, features = DC_markers, cluster.idents = TRUE)+
  RotatedAxis()+
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

dev.off()

options(repr.plot.width=8, repr.plot.height=4)
#GO交集的基因
DotPlot(DC_so, features = c('AXL','NRARP','LYN','IL4I1','HAVCR2','RIPOR2','LILRB4','LST1'), cluster.idents = TRUE)+
  RotatedAxis()+
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

pdf(file = "plots/figures/figS9B_dc_expression.pdf", width = 8, height = 4)
DotPlot(DC_so, features = c('AXL','NRARP','LYN','IL4I1','HAVCR2','RIPOR2','LILRB4','LST1'), cluster.idents = TRUE)+
  RotatedAxis()+
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

dev.off()

options(repr.plot.width=8, repr.plot.height=4)
#GO交集的基因
DotPlot(DC_so, features = c('AXL','NRARP','LYN','IL4I1','HAVCR2','RIPOR2','LILRB4','LST1'), cluster.idents = TRUE)+
  RotatedAxis()#+
  # scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

# pdf(file = "plots/figures/figS9B_dc_expression.pdf", width = 8, height = 4)
# DotPlot(DC_so, features = c('AXL','NRARP','LYN','IL4I1','HAVCR2','RIPOR2','LILRB4','LST1'), cluster.idents = TRUE)+
#   RotatedAxis()+
#   scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

# dev.off()

neutro_so <- subset(myeloid.batch.corrected, subtype1 %in% c("CCL4+ Neutro", "S100A12+ Neutro","MNDA+ Neutro"))

# fig
options(repr.plot.width=10, repr.plot.height=4)
VlnPlot(neutro_so, 
#         stack = TRUE,
       features = c("OSM","TGFB1","MMP9"),
       pt.size=0)+RotatedAxis()+NoLegend()

pdf(file = "plots/figures/figS7D_neutro_expression.pdf", width = 10, height = 4)
VlnPlot(neutro_so, 
#         stack = TRUE,
       features = c("OSM","TGFB1","MMP9"),
       pt.size=0)+RotatedAxis()+NoLegend()

dev.off()

# DefaultAssay(neutro_so) <- "RNA"

# neutro.markers_1.5 <- FindAllMarkers(neutro_so, only.pos = TRUE, min.pct = 0.25, 
#                                       logfc.threshold = 0.25)
# neutro.markers_1.5 %>%
#     group_by(cluster) %>%
#     slice_max(n = 4, order_by = avg_log2FC)

# write.csv(neutro.markers_1.5,"processed_data/data_B2-19/cell_typing/Neutrophil/neutro_3ct_markers.csv")

neutro.markers_1.5 <- read.csv("processed_data/data_B2-19/cell_typing/Neutrophil/neutro_3ct_markers.csv",
                               row.names = 1)

neutro_top4_markers <- neutro.markers_1.5 %>%
    group_by(cluster) %>%
    slice_max(n = 7, order_by = avg_log2FC) %>%
    pull(gene)

neutro_top4_markers

options(repr.plot.width=10, repr.plot.height=4)
DotPlot(neutro_so, features = c(unique(neutro_top4_markers)),group.by = "subtype1")+RotatedAxis()+
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")


pdf(file = "plots/figures/figS7B_neutro_expression.pdf", width = 10, height = 4)
DotPlot(neutro_so, features = c(unique(neutro_top4_markers)),group.by = "subtype1")+RotatedAxis()+
  scale_color_gradient(low = "lightgrey", high = "#c45a48", name = "Expression Level")

dev.off()

options(repr.plot.width=10, repr.plot.height=4)
DotPlot(neutro_so, features = c(unique(neutro_top4_markers)),group.by = "subtype1")+RotatedAxis()

# from run_monocle
library(Seurat)
library(tidyverse)
library(magrittr)
library(monocle)

cds_plasma <- readRDS("processed_data/data_B2-19/monocle/plasma_self_markers_root.rds")

# fig
options(repr.plot.width=6, repr.plot.height=5)
my_genes <- c("MKI67","IGHA1","IGHG1","IGHM")
cds_plasma_subset <- cds_plasma[my_genes,]
p <- plot_genes_in_pseudotime(cds_plasma_subset, color_by = "cell_type",ncol = 2)
plot(p)
pdf("plots/figures/fig5E_monocle.pdf",width=6, height=5)
plot(p)
dev.off()

options(repr.plot.width=6, repr.plot.height=4)
p <- plot_cell_trajectory(cds_plasma, color_by = "cell_type")+ theme(legend.position="right")
plot(p)
pdf("plots/figures/fig5D_monocle.pdf",width=6, height=4)
plot(p)
dev.off()

cds_nk <- readRDS("processed_data/data_B2-19/monocle/nk_self_markers_root.rds")

# pdf("monocle/nk_all.pdf")
options(repr.plot.width=4, repr.plot.height=4)
p1 <- plot_cell_trajectory(cds_nk, color_by = "cell_type") 
plot(p1)
pdf("plots/figures/figS4B_monocle.pdf",width=4, height=4)
plot(p1)
dev.off()

p2 <- monocle::plot_cell_trajectory(cds_nk, color_by = "Pseudotime")
plot(p2)
pdf("plots/figures/figS4C_monocle.pdf",width=4, height=4)
plot(p2)
dev.off()

monocle::plot_cell_trajectory(cds_nk, color_by = "Pseudotime") +
  scale_color_gradient(low = "#89afd2", high = "#ef8775")
# dev.off()

# font_theme <- theme(
#   plot.title = element_text(size = 8, face = "bold", family = "Arial"), # 修改标题字体、字号、字体系列
#   axis.title.x = element_text(size = 7, family = "Arial"),  # 修改x轴标题的字号和字体系列
#   axis.title.y = element_text(size = 7, family = "Arial"),  # 修改y轴标题的字号和字体系列
#   axis.text.x = element_text(size = 7, family = "Arial"),   # 修改x轴标签的字号和字体系列
#   axis.text.y = element_text(size = 7, family = "Arial"),   # 修改y轴标签的字号和字体系列
#   legend.title = element_text(size = 7, family = "Arial"),  # 修改图例标题的字号和字体系列
#   legend.text = element_text(size = 7, family = "Arial")    # 修改图例标签的字号和字体系列
# )

font_theme <- theme(
  plot.title = element_text(size = 14, face = "bold", family = "Arial"), # 修改标题字体、字号、字体系列
  axis.title.x = element_text(size = 12, family = "Arial"),  # 修改x轴标题的字号和字体系列
  axis.title.y = element_text(size = 12, family = "Arial"),  # 修改y轴标题的字号和字体系列
  axis.text.x = element_text(size = 12, family = "Arial"),   # 修改x轴标签的字号和字体系列
  axis.text.y = element_text(size = 12, family = "Arial"),   # 修改y轴标签的字号和字体系列
  legend.title = element_text(size = 12, family = "Arial"),  # 修改图例标题的字号和字体系列
  legend.text = element_text(size = 12, family = "Arial")    # 修改图例标签的字号和字体系列
)

#from subtyping_epi
epi_go <- readRDS("processed_data/data_B2-19/cell_typing/epithelial/GO_all_conditions.rds")

options(repr.plot.width=7, repr.plot.height=8)
custom_colors <- rev(c("#4792c4", "#9cc2d7", "#f5f4e9", "#f7ddc6", "#e7a988", "#c86955"))
go_epi <- dotplot(epi_go,showCategory=3)+ 
  scale_color_gradientn(colors = custom_colors) +
#   scale_color_gradient(low = "#bb5a74", high = "#bfdaf3")+
    theme(axis.text.x = element_text(angle = 45, hjust=1))#+font_theme
plot(go_epi)
pdf("plots/figures/fig1F_epi_go.pdf", width = 7, height = 8)
plot(go_epi)
dev.off()


# try color
options(repr.plot.width=7, repr.plot.height=8)
custom_colors <- c("#60a5a0", "#99c7c9", "#b9dadb", "whitesmoke")
dotplot(epi_go,showCategory=3)+ 
scale_color_gradientn(colors = custom_colors) +
#   scale_color_gradient(low = "#bb5a74", high = "#bfdaf3")+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+font_theme

go_hla <- readRDS("processed_data/data_B2-19/cell_typing/tcell/HLA_DRA_Tem/GO_HLA_DRA_Tem_Tex.rds")
# data from tcell_functions.ipynb

options(repr.plot.width=6, repr.plot.height=5)
custom_colors <- rev(c("#4792c4", "#9cc2d7", "#f5f4e9", "#f7ddc6", "#e7a988", "#c86955"))

p <- dotplot(go_hla,showCategory=c("positive regulation of cell activation",
                              "negative regulation of cell activation",
                              "positive regulation of leukocyte activation",
                              "negative regulation of leukocyte activation",
                              "positive regulation of T cell activation",
                              "negative regulation of T cell activation"
                             ))+
    scale_color_gradientn(colors = custom_colors) #+  #,by="p.adjust"
#     theme(axis.text.x = element_text(angle = 45, hjust=1))
plot(p)
pdf("plots/figures/figS5F_go.pdf", width = 6, height = 5)
plot(p)
dev.off()


go_rrm2 <- readRDS("processed_data/data_B2-19/cell_typing/bcell/RRM2_Plasma/go_compare_RRM2_Plasma.RDS")

options(repr.plot.width=6, repr.plot.height=8)
custom_colors <- rev(c("#4792c4", "#9cc2d7", "#f5f4e9", "#f7ddc6", "#e7a988", "#c86955"))

p <- dotplot(go_rrm2,showCategory=5)+
    scale_color_gradientn(colors = custom_colors) +  #,by="p.adjust"
    theme(axis.text.x = element_text(angle = 45, hjust=1))
plot(p)
pdf("plots/figures/fig5C_go.pdf", width = 6, height = 8)
plot(p)
dev.off()

options(repr.plot.width=6, repr.plot.height=7)

dotplot(go_rrm2,showCategory=c("immunoglobulin complex",
                                    "endoplasmic reticulum protein-containing complex",
                              "IgG immunoglobulin complex",
                              "IgA immunoglobulin complex",
                              "immunoglobulin complex, circulating"
                             ))+  
    theme(axis.text.x = element_text(angle = 30, hjust=1))

GO_RRM2_Plasmablast <- readRDS("processed_data/data_B2-19/cell_typing/bcell/RRM2_Plasma/GO_RRM2_Plasma_Plasmablast.rds")

options(repr.plot.width=6, repr.plot.height=7)

p <- barplot(GO_RRM2_Plasmablast$up_BP, showCategory=10, title=paste0("RRM2+ Plasma")) +
  scale_fill_gradient(low = "#D75A4A", high = "#4A90D9")

    # scale_color_gradientn(colors = custom_colors)
plot(p)
pdf("plots/figures/figS10B_go.pdf", width = 6, height = 7)
plot(p)
dev.off()

GO_GZMB_Temra_GZMK_Temra <- readRDS("processed_data/data_B2-19/cell_typing/tcell/GZMB_Temra/GO_GZMB_Temra_GZMK_Temra.rds")

options(repr.plot.width=6, repr.plot.height=7)

p <- barplot(GO_GZMB_Temra_GZMK_Temra$up_BP, showCategory=10, title=paste0("GZMB+ Temra")) +
  scale_fill_gradient(low = "#D75A4A", high = "#4A90D9")
plot(p)
pdf("plots/figures/figS6E_go.pdf", width = 6, height = 7)
plot(p)
dev.off()

GO_cDC2_markers <- readRDS("processed_data/data_B2-19/cell_typing/DC/FCER1A_cDC2/GO_FCER1A_cDC2.rds")

options(repr.plot.width=6, repr.plot.height=7.3)

p <- dotplot(GO_cDC2_markers$down_BP, showCategory=14, title=paste0("NLRP3+ cDC2")) +
    scale_color_gradientn(colors = custom_colors)
plot(p)
pdf("plots/figures/figS9F_go.pdf", width = 6, height = 7.3)
plot(p)
dev.off()

barplot(GO_cDC2_markers$down_BP, showCategory=15, title="NLRP3+ cDC2") +
  scale_fill_gradient(low = "#D75A4A", high = "#4A90D9")

library(EnhancedVolcano)
library(ggplot2)

lamp3_de_markers <- read.csv("processed_data/data_B2-19/cell_typing/DC/LAMP3_DC/lamp3_lauren_de_markers.csv",
                         row.names = 1)

lamp3_de_markers[lamp3_de_markers$gene=="CD274",]

options(repr.plot.width=5, repr.plot.height=6)

p <- EnhancedVolcano(lamp3_de_markers, 
                lab = lamp3_de_markers$gene,
                x = "avg_log2FC", 
                y = "p_val_adj",
#                selectLab = "CD274",
                title = "LAMP3+ DC",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black') +
  theme_classic() +
  theme(legend.position = "top")
plot(p)
pdf("plots/figures/fig6A_DE.pdf", width = 5, height = 6)
plot(p)
dev.off()

plasma_de_markers <- read.csv("processed_data/data_B2-19/cell_typing/bcell/Plasma/plasma_LN_condition_de_markers.csv",
                         row.names = 1)

options(repr.plot.width=7, repr.plot.height=6)

p <- EnhancedVolcano(plasma_de_markers, 
                lab = plasma_de_markers$gene,
                x = "avg_log2FC", 
                y = "p_val_adj",
#                selectLab = "CD274",
                title = "Plasma",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black') +
  theme_classic() +
  theme(legend.position = "top")
plot(p)
pdf("plots/figures/figS10E_DE.pdf", width = 7, height = 6)
plot(p)
dev.off()

library(scRepertoire)
# suppressMessages(library(scRepertoire))
library(Seurat)
library(reshape2)
library(pheatmap)

seurat_t_tcr <- readRDS("processed_data/data_B2-19/TCR/tcell_with_tcr.rds")

seurat_t_all <- readRDS("processed_data/data_B2-19/TCR/tcell_all.rds")

combined_ct_tcr <- expression2List(seurat_t_tcr, split.by = "subtype1")
length(combined_ct_tcr)

options(repr.plot.width=4, repr.plot.height=5)
p <- compareClonotypes(combined_ct_tcr, 
                  numbers = 1000, 
                  samples = c('CD8_LAYN+ Tex','CD8_CCL5+ Tex'),
                    cloneCall="gene+nt", graph = "alluvial")+NoLegend()
# ggsave("plots/figures/figS6F_tcr1.pdf",width = 4,height = 2.8)


p

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
              "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

options(repr.plot.width=8, repr.plot.height=5)
slot(seurat_t_all, "meta.data")$cloneType <- factor(slot(seurat_t_all, "meta.data")$cloneType, 
                levels = c("Hyperexpanded (100 < X <= 500)", 
                           "Large (20 < X <= 100)", 
                            "Medium (5 < X <= 20)", 
                            "Small (1 < X <= 5)", 
                            "Single (0 < X <= 1)", NA))
p <- DimPlot(seurat_t_all, group.by = "cloneType") +
    scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())
plot(p)
pdf("plots/figures/fig3G_umap.pdf", width = 8,height = 5)
plot(p)
dev.off()

options(repr.plot.width=8, repr.plot.height=5)
slot(seurat_t_tcr, "meta.data")$cloneType <- factor(slot(seurat_t_tcr, "meta.data")$cloneType, 
                levels = c("Hyperexpanded (100 < X <= 500)", 
                           "Large (20 < X <= 100)", 
                            "Medium (5 < X <= 20)", 
                            "Small (1 < X <= 5)", 
                            "Single (0 < X <= 1)"))
DimPlot(seurat_t_tcr, group.by = "cloneType") +
    scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

options(repr.plot.width=14, repr.plot.height=6)
quantContig(combined_ct_tcr, cloneCall="gene+nt", scale = T)+RotatedAxis()
# B6T 的TCR只是warning，map比例有45.22%，另外valid barcode 83.27%,低于85%了，除了B6T以外，其他肿瘤样本均有一定程度扩增。



seurat_t_tcr_temra <- subset(seurat_t_tcr, 
                       subtype1 %in% c('CD8_MKI67+ Tex','CD8_LAYN+ Tex','CD8_CCL5+ Tex',"CD8_GZMB+ Temra"))
#                            c('CD8_LAYN+ Tex','CD8_CCL5+ Tex',"CD4_Tex"))
combined_ct_tcr_temra <- expression2List(seurat_t_tcr_temra, split.by = "subtype1")
length(combined_ct_tcr_temra)

options(repr.plot.width=5, repr.plot.height=4)
clonalOverlap(combined_ct_tcr_temra, cloneCall = "gene+nt", method = "overlap") + 
  scale_fill_gradient2(low = "lightgrey", high = "#D75A4A",na.value = "white") +
  RotatedAxis() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
ggsave("plots/figures/figS6F_tcr.pdf",width = 5,height = 4)

options(repr.plot.width=4, repr.plot.height=2.8)
compareClonotypes(combined_ct_tcr_temra, 
                  numbers = 1000, 
                  samples = c('CD8_LAYN+ Tex','CD8_GZMB+ Temra'),
                    cloneCall="gene+nt", graph = "alluvial")+NoLegend()
ggsave("plots/figures/figS6F_tcr1.pdf",width = 4,height = 2.8)


options(repr.plot.width=4, repr.plot.height=2.8)
compareClonotypes(combined_ct_tcr_temra, 
                  numbers = 1000, 
                  samples = c('CD8_CCL5+ Tex','CD8_GZMB+ Temra'),
                    cloneCall="gene+nt", graph = "alluvial")+NoLegend()
ggsave("plots/figures/figS6F_tcr2.pdf",width = 4,height = 2.8)


seurat_t_tcr_tim3 <- subset(seurat_t_tcr, 
                       subtype1 %in% c('CD8_LAYN+ Tex','CD8_CCL5+ Tex','CD8_TIM3+ Trm','CD8_MKI67+ Tex'))
combined_ct_tcr_tim3 <- expression2List(seurat_t_tcr_tim3, split.by = "subtype1")

options(repr.plot.width=5, repr.plot.height=4)
clonalOverlap(combined_ct_tcr_tim3, cloneCall = "gene+nt", method = "overlap") + 
  scale_fill_gradient2(low = "lightgrey", high = "#D75A4A",na.value = "white") +
  RotatedAxis() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
ggsave("plots/figures/figS2D_tcr.pdf",width = 5,height = 4)


# seurat_t_tcr$subtype1.sample <- paste(seurat_t_tcr$subtype1, seurat_t_tcr$sample, sep = "_")
# seurat_t_tcr$subtype1.LN_condition <- paste(seurat_t_tcr$subtype1, seurat_t_tcr$LN_condition, 
#                                              sep = "_")

seurat_t_tcr_HLA_Tex <- subset(seurat_t_tcr, 
                       subtype1 %in% c('CD8_MKI67+ Tex','CD8_LAYN+ Tex','CD8_CCL5+ Tex',"CD8_HLA-DRA+ Tem"))
combined_ct_tcr_tex <- expression2List(seurat_t_tcr_HLA_Tex, split.by = "subtype1")

options(repr.plot.width=5, repr.plot.height=4)
clonalOverlap(combined_ct_tcr_tex, cloneCall = "gene+nt", method = "overlap") + 
  scale_fill_gradient2(low = "lightgrey", high = "#D75A4A",na.value = "white") +
  RotatedAxis() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
ggsave("plots/figures/fig3A_tcr.pdf",width = 5,height = 4)


# mask <- seurat_t_tcr_HLA_Tex$CTstrict %in% unique(seurat_t_tcr_HLA_Tex$CTstrict[seurat_t_tcr_HLA_Tex$Frequency>5])
HLA_Tex_expansion_meta <- seurat_t_tcr_HLA_Tex@meta.data


df <- data.frame(table(HLA_Tex_expansion_meta$subtype1, HLA_Tex_expansion_meta$CTstrict)) 
colnames(df) <- c("subtype","TCR","count")

tcr_sum <- aggregate(df$count,list(df$TCR),sum)
colnames(tcr_sum) <- c("TCR","count_sum")

mask_2 <- df$TCR %in% unique(tcr_sum$TCR[tcr_sum$count_sum>3])
df <- df[mask_2,]
# tcr in hla and tex cells >3

df_mat<-dcast(df,subtype~TCR,value.var = 'count')
rownames(df_mat)<- df_mat$subtype
df_mat$subtype <- NULL
dim(df_mat)
df_mat[,1:2]

df_inter_tcr <- df_mat[,colSums(df_mat>1)>1 & df_mat["CD8_HLA-DRA+ Tem",]>0]
#至少两个类型共享，且HLA亚型中有至少1个该TCR的cell
dim(df_inter_tcr)
# 将数值限制在最大值为5
df_inter_tcr[df_inter_tcr > 5] <- 5

options(repr.plot.width=9, repr.plot.height=3)
# 自定义 breaks 和 labels
breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
labels <- c("0", "1", "2", "3", "4", ">=5")

# 指定离散颜色
# colors <- c("#4575b4", "#9cc6df", "#ebf7e3", "#feeda3", "#fc9e64", "#d6281e")
colors <- c("lightgrey",colorRampPalette(c("#feddc5", "#d6281e"))(length(breaks) - 1))

# 绘制热图并自定义颜色和图例标签
p <- pheatmap(
  df_inter_tcr,
  show_colnames = FALSE, 
  cluster_rows = F,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  breaks = breaks,
  color = colors,
  legend_labels = labels,
  legend_breaks = c(0, 1, 2, 3, 4, 5)
)
p
pdf("plots/figures/figS5A_tcr.pdf",width = 9,height = 3)
p
dev.off()

options(repr.plot.width=10, repr.plot.height=3)
# 自定义 breaks 和 labels
breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
labels <- c("0", "1", "2", "3", "4", ">=5")

# 指定离散颜色
colors <- c("#4575b4", "#9cc6df", "#ebf7e3", "#feeda3", "#fc9e64", "#d6281e")

# 绘制热图并自定义颜色和图例标签
pheatmap(
  df_inter_tcr,
  show_colnames = FALSE, 
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  breaks = breaks,
  color = colors,
  legend_labels = labels,
  legend_breaks = c(0, 1, 2, 3, 4, 5)
)



seurat_t_tcr_HLA <- readRDS("processed_data/data_B2-19/TCR/seurat_t_tcr_HLA.rds")

slot(seurat_t_tcr_HLA, "meta.data")$cloneType <- factor(slot(seurat_t_tcr_HLA, "meta.data")$cloneType, 
                levels = c("Hyperexpanded (100 < X <= 500)", 
                           "Large (20 < X <= 100)", 
                            "Medium (5 < X <= 20)", 
                            "Small (1 < X <= 5)", 
                            "Single (0 < X <= 1)"))

options(repr.plot.width=8, repr.plot.height=4)
p <- VlnPlot(seurat_t_tcr_HLA, features = c("TOX","PDCD1","HAVCR2", "LAG3","CTLA4","TIGIT",'LAYN'),
        stack=TRUE, group.by="cloneType")+
  RotatedAxis()+ NoLegend()
plot(p)

pdf("plots/figures/figS5B_expression.pdf",width=8, height=4)
plot(p)
dev.off()


library(ggplot2)
library(dplyr)
library(Seurat)
library(GSVA)

library(pheatmap)
library(patchwork)
library(msigdbr)


kegg <- readRDS("processed_data/data_B2-19/cell_typing/myeloid/subtype1_GSVA.rds")

# by cell_type1
options(repr.plot.width=12, repr.plot.height=10)
p_cell_type1<-pheatmap(kegg, show_rownames=1, show_colnames=T)
ggsave(p_cell_type1,
       filename = "plots/figures/figS8B_GSVA_myeloid.pdf",width=12,height=10)

# by cell_type1
my.colors <- c(
    colorRampPalette(colors = c("#89afd2", "white", "#ef8775"))(100))

options(repr.plot.width=12, repr.plot.height=10)
p_cell_type1<-pheatmap(kegg, show_rownames=1, show_colnames=T,color = my.colors)
ggsave(p_cell_type1,
       filename = "plots/figures/figS8B_color_2_GSVA_myeloid.pdf",width=12,height=10)



# epithelial.batch.corrected <- readRDS("plots/data/epithelial_new_meta.rds")

kegg_patient_condition_mean <- readRDS("processed_data/data_B2-19/cell_typing/epithelial/patient_condition_mean_GSVA.rds")

# by patient condition
options(repr.plot.width=5, repr.plot.height=2.5)

p_patient_condition_mean<-pheatmap(rescale(kegg_patient_condition_mean[c("HALLMARK_ANGIOGENESIS",
                                                                             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS"),], 
                                             to=c(-1, 1)),
                                   show_rownames=1, show_colnames=T,
                                   cluster_cols = FALSE,cluster_rows = FALSE)

ggsave(p_patient_condition_mean,
       filename = "plots/figures/figS1C_GSVA_epi_patient_condition.pdf",width=6,height=3)

kegg_Lauren.LN_condition_mean <- readRDS("processed_data/data_B2-19/cell_typing/epithelial/Lauren.LN_condition_mean_GSVA.rds")

# by patient condition
options(repr.plot.width=5, repr.plot.height=2.5)
p_Lauren.LN_condition_mean<-pheatmap(rescale(kegg_Lauren.LN_condition_mean[c("HALLMARK_ANGIOGENESIS",
                                                                             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS"),], 
                                             to=c(-1, 1)),
                                   show_rownames=1, show_colnames=T,
                                   cluster_cols = FALSE,cluster_rows = FALSE)

ggsave(p_Lauren.LN_condition_mean,
       filename = "plots/figures/figS1C_GSVA_epi_lauren.pdf",width=6,height=3)

library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)

options(stringsAsFactors = FALSE)

packageVersion("CellChat")

# cellchat.patient_meta_T <- readRDS("processed_data/data_B2-19/cellchat_v1/cellchat_all_subtype/cellchat.patient_meta_T.rds")
cellchat.patient_neg_T <- readRDS("processed_data/data_B2-19/cellchat_v1/cellchat_all_subtype/cellchat.patient_neg_T.rds")

# object.list <- list(patient_meta_T = cellchat.patient_meta_T, patient_neg_T = cellchat.patient_neg_T)
# cellchat_patient_meta_T_patient_neg_T <- mergeCellChat(object.list, add.names = names(object.list))

# cellchat_patient_meta_T_patient_neg_T <- liftCellChat(cellchat_patient_meta_T_patient_neg_T, 
#                                                                 group.new = levels(cellchat_patient_meta_T_patient_neg_T@idents$joint) )


# figure
options(repr.plot.width=9, repr.plot.height=5)
netVisual_bubble(cellchat.patient_neg_T, targets.use = which(rownames(cellchat.patient_neg_T@net$count) %in% c("CD8_HLA-DRA+ Tem",
                                                                                                               "CD8_LAYN+ Tex",
                                                                                                               "CD8_CCL5+ Tex"
                                                                                                              )),
                 #targets.use = myeloid_group,
                 signaling=c("CD80","PD-L1","CDH1"),
                 remove.isolate = TRUE,thresh=0.01)#targets.use = c(5:11),
ggsave(filename = "plots/figures/fig3D_cellchat.pdf",width=9, height=5)

cellchat.patient_neg_LN_neg <- readRDS("processed_data/data_B2-19/cellchat_v1/cellchat_all_subtype/cellchat.patient_neg_LN_neg.rds")
cellchat.patient_meta_LN_neg <- readRDS("processed_data/data_B2-19/cellchat_v1/cellchat_all_subtype/cellchat.patient_meta_LN_neg.rds")

meta <- read.csv("processed_data/data_B2-19/metadata_all_ct_subtype_condition_noTLB.csv",row.names = 1)

grep("Tgd",unique(meta$subtype1),value=TRUE)

# object.list <- list(patient_neg_LN_neg = cellchat.patient_neg_LN_neg, patient_meta_LN_neg = cellchat.patient_meta_LN_neg)
# cellchat_patient_neg_LN_neg_patient_meta_LN_neg <- mergeCellChat(object.list, add.names = names(object.list))

# cellchat_patient_neg_LN_neg_patient_meta_LN_neg <- liftCellChat(cellchat_patient_neg_LN_neg_patient_meta_LN_neg, 
#                                                                 group.new = levels( cellchat_patient_neg_LN_neg_patient_meta_LN_neg@idents$joint) )


dc_ct <- unique(meta$subtype1[meta$cell_type1=="DC"])
dc_group <- which(rownames(cellchat.patient_meta_LN_neg@net$count) %in% dc_ct)

options(repr.plot.width=4, repr.plot.height=4)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 52, targets.use = dc_group,
                 signaling="IL1",
                 remove.isolate = TRUE,thresh=0.01,return.data = TRUE)
ggsave(filename = "plots/figures/figS8D_cellchat.pdf",width=4, height=4)


options(repr.plot.width=9, repr.plot.height=5)
patient_neg_LN_neg_VEGF <- netVisual_bubble(cellchat.patient_neg_LN_neg, targets.use = 47,#targets.use = myeloid_group,
                 signaling="VEGF",
                 remove.isolate = TRUE,thresh=0.01,return.data = TRUE)#targets.use = c(5:11),
patient_neg_LN_neg_VEGF$gg.obj
ggsave(filename = "plots/figures/figS8C_cellchat.pdf",width=9, height=5)
write.csv(patient_neg_LN_neg_VEGF$communication,"plots/data/cellchat/patient_neg_LN_neg_VEGF.csv")

options(repr.plot.width=10, repr.plot.height=5)
patient_meta_LN_neg_VEGF <- netVisual_bubble(cellchat.patient_meta_LN_neg, targets.use = 47,#targets.use = myeloid_group,
                 signaling="VEGF",
                 remove.isolate = TRUE,thresh=0.01,return.data = TRUE)#targets.use = c(5:11),
patient_meta_LN_neg_VEGF$gg.obj
ggsave(filename = "plots/figures/fig4D_cellchat.pdf",width=10, height=5)
write.csv(patient_meta_LN_neg_VEGF$communication,"plots/data/cellchat/patient_metga_LN_neg_VEGF.csv")

max(patient_neg_LN_neg_VEGF$communication$prob[patient_neg_LN_neg_VEGF$communication$ligand=="VEGFA"])

max(patient_meta_LN_neg_VEGF$communication$prob[patient_meta_LN_neg_VEGF$communication$ligand=="VEGFA"])

neutro_ct <- unique(meta$subtype1[meta$cell_type1=="Neutro"])
neutro_group <- which(rownames(cellchat.patient_neg_LN_neg@net$count) %in% neutro_ct)

options(repr.plot.width=4, repr.plot.height=5)
netVisual_bubble(cellchat.patient_neg_LN_neg, sources.use = 52, targets.use = neutro_group,
                 signaling="CXCL",
                 remove.isolate = TRUE,thresh=0.01,return.data = TRUE)

ggsave(filename = "plots/figures/fig4C_cellchat.pdf",width=4, height=5)


neutro_ct <- unique(meta$subtype1[meta$cell_type1=="Neutro"])
neutro_group <- which(rownames(cellchat.patient_meta_LN_neg@net$count) %in% neutro_ct)

options(repr.plot.width=4, repr.plot.height=5)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 52, targets.use = neutro_group,
                 signaling="CXCL",
                 remove.isolate = FALSE,thresh=0.01,return.data = TRUE)
ggsave(filename = "plots/figures/fig4B_cellchat.pdf",width=4, height=5)



cellchat.LN_meta_diffuse <- readRDS("processed_data/data_B2-19/cellchat_v1/cellchat.LN_meta_diffuse.rds")
cellchat.LN_meta_intestinal <- readRDS("processed_data/data_B2-19/cellchat_v1/cellchat.LN_meta_intestinal.rds")

object.list <- list(Intestianl_Met.LN = cellchat.LN_meta_intestinal, 
                    Diffuse_Met.LN = cellchat.LN_meta_diffuse)
cellchat_LN_meta_intestinal_LN_meta_diffuse <- mergeCellChat(object.list, 
                                                             add.names = names(object.list))

pathways.show <- c("PD-L1") 
options(repr.plot.width=16, repr.plot.height=6)
netAnalysis_signalingRole_network(cellchat.LN_meta_diffuse, signaling = pathways.show, 
                                  width = 36, height = 10, font.size = 10)#+font_theme

pdf("plots/figures/fig6C_pathway.pdf",width=16, height=6)
netAnalysis_signalingRole_network(cellchat.LN_meta_diffuse, signaling = pathways.show, 
                                  width = 36, height = 10, font.size = 10)#+font_theme
dev.off()

pathways.show <- c("IL10") 
options(repr.plot.width=16, repr.plot.height=6)
netAnalysis_signalingRole_network(cellchat.LN_meta_diffuse, signaling = pathways.show, 
                                  width = 36, height = 10, font.size = 10)#+font_theme

pdf("plots/figures/figS11D_pathway.pdf",width=16, height=6)
netAnalysis_signalingRole_network(cellchat.LN_meta_diffuse, signaling = pathways.show, 
                                  width = 36, height = 10, font.size = 10)#+font_theme
dev.off()

pathways.show <- c("IL10") 
options(repr.plot.width=16, repr.plot.height=6)
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.LN_meta_diffuse, signaling = pathways.show, 
                                  width = 36, height = 10, font.size = 10)#+font_theme

macro_group <- which(rownames(cellchat.LN_meta_diffuse@net$count) %in% c("SPP1+ Mph","C1QC+ Mph","MKI67+ Mph",
                                                                         # "CD8_TIM3+MKI67+ Trm",
                                                                         "CLEC10A+ Mph"
                                                                           ))

lamp3 <- which(rownames(cellchat.LN_meta_diffuse@net$count) %in% c("LAMP3+ DC"))

options(repr.plot.width=7, repr.plot.height=6)
netVisual_bubble(cellchat_LN_meta_intestinal_LN_meta_diffuse, sources.use = macro_group, targets.use = lamp3,  
                 comparison = c(1, 2), 
                 signaling="IL10")

ggsave("plots/figures/fig6B_bigger_cellchat.pdf",width=7, height=6)


options(repr.plot.width=8, repr.plot.height=5)
netVisual_bubble(cellchat_LN_meta_intestinal_LN_meta_diffuse, sources.use = macro_group, targets.use = lamp3,  
                 comparison = c(1, 2), 
                 signaling="IL10",
                 angle.x = 45)
ggsave("plots/figures/fig6B_cellchat.pdf",width=8, height=5)

# macro_group <- which(rownames(cellchat.LN_meta_diffuse@net$count) %in% c("SPP1+ Mph","C1QC+ Mph","MKI67+ Mph","CLEC10A+ Mph"))

options(repr.plot.width=7, repr.plot.height=4)
netVisual_bubble(cellchat.LN_meta_diffuse, sources.use = macro_group, targets.use = lamp3,
                 signaling="IL10",
                 remove.isolate = TRUE,return.data = FALSE)

options(repr.plot.width=7, repr.plot.height=14)

gg1 <- rankNet(cellchat_LN_meta_intestinal_LN_meta_diffuse, mode = "comparison", 
               stacked = T, do.stat = TRUE,
               ylim = c(0, 10))
               # , segments = list(c(11, 14),c(16, 28)), 
               # tick_width = c(5,2,5), rel_heights = c(0.8,0,0.1,0,0.1))
gg1
# ggsave("plots/figures/figS11A_cellchat.pdf",width=7, height=14)




options(repr.plot.width=7, repr.plot.height=14)

gg1 <- rankNet(cellchat_LN_meta_intestinal_LN_meta_diffuse, mode = "comparison", 
               stacked = T, do.stat = TRUE)
gg1
ggsave("plots/figures/figS11A_cellchat.pdf",width=7, height=14)


options(repr.plot.width=6, repr.plot.height=8)

gg1 <- rankNet(cellchat_LN_meta_intestinal_LN_meta_diffuse, mode = "comparison", 
               targets.use=c("LAMP3+ DC"),
               stacked = T, do.stat = TRUE)
gg1
ggsave("plots/figures/figS11B_cellchat.pdf",width=6, height=8)


# # 加载必要的包
# library(dplyr)
# library(tidyr)

# tmp <- gg1$signaling.contribution[,c("name","contribution","group","pvalues")]

# df_wide <- tmp %>%
#   pivot_wider(names_from = group, values_from = contribution, names_prefix = "contribution_")

# P_plot <- df_wide$name[df_wide$contribution_LN_meta_intestinal < df_wide$contribution_LN_meta_diffuse & df_wide$pvalues < 0.01]
# # try to plot only target pathways, failed

# str(cellchat_LN_meta_intestinal_LN_meta_diffuse1@netP)



library(ggplot2)
library(dplyr)
library(Seurat)
library(GSVA)

library(pheatmap)
library(patchwork)
library(msigdbr)

library(nichenetr)
library(tidyverse)

vis_ligand_target <- readRDS("plots/NicheNet/IL10_CD274.rds")
# from run_nichenet.ipynb

options(repr.plot.width=4, repr.plot.height=3)
vis_ligand_target %>% 
make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                    color = "purple",legend_position = "right", 
                    x_axis_position = "top",legend_title = "Regulatory potential")  + 
theme(axis.text.x = element_text(face = "italic")) + 
scale_fill_gradient2(low = "#c98df8", high = "#9930ef")
ggsave("plots/figures/figS11C_nichenet.pdf",width=4, height=3)



library(ggplot2)
library(dplyr)
library(Seurat)
library(GSVA)

library(pheatmap)
library(patchwork)
library(msigdbr)


myeloid.batch.corrected <- readRDS("processed_data/data_B2-19/cell_typing/myeloid/myeloid.batch.corrected_ct_subtype.rds")
DefaultAssay(myeloid.batch.corrected)<-"RNA"


subtype_mean_myeloid=data.frame(row.names = rownames(myeloid.batch.corrected))
for(i in unique(myeloid.batch.corrected$subtype1)){
#     print(i)
    sub_i=subset(myeloid.batch.corrected, subset=subtype1==i)
    subtype_mean_myeloid[,as.character(i)]=apply(sub_i@assays$RNA@data, 1, mean)
}

# find differential pathway -----------------------------------------------
m_df = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)

dim(subtype_mean_myeloid)

# packageVersion("matrixStats")

# kegg <- gsva(as.matrix(subtype_mean_myeloid), msigdbr_list, kcdf="Gaussian",method = "gsva",parallel.sz=10)


# saveRDS(kegg,"processed_data/data_B2-19/cell_typing/myeloid/subtype1_GSVA.rds")

kegg <- readRDS("processed_data/data_B2-19/cell_typing/myeloid/subtype1_GSVA.rds")

# by cell_type1
options(repr.plot.width=12, repr.plot.height=10)
p_cell_type1<-pheatmap(kegg, show_rownames=1, show_colnames=T)

epithelial.batch.corrected <- readRDS("plots/data/epithelial_new_meta.rds")

epithelial.batch.corrected@meta.data[,c('LN_condition.LN_station','patient_condition.region','patient_condition.LN_condition',
       'LN_condition.Lauren','condition')] <- NULL

DefaultAssay(epithelial.batch.corrected)

epithelial.batch.corrected$patient_condition.LN_condition <- paste(epithelial.batch.corrected$patient_condition, 
                                                                   epithelial.batch.corrected$LN_condition, sep = "_")
table(epithelial.batch.corrected$patient_condition.LN_condition)

patient_condition_mean=data.frame(row.names = rownames(epithelial.batch.corrected))
for(i in c("Meta_Pri.GC", "Neg_Pri.GC")){
    sub_i=subset(epithelial.batch.corrected, subset=patient_condition.LN_condition==i)
    patient_condition_mean[,as.character(i)]=apply(sub_i@assays$RNA@data, 1, mean)
}

# find differential pathway -----------------------------------------------
m_df = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)

kegg_patient_condition_mean <- gsva(as.matrix(patient_condition_mean), msigdbr_list, 
                                    kcdf="Gaussian",method = "gsva",parallel.sz=10)


# saveRDS(kegg_patient_condition_mean,"processed_data/data_B2-19/cell_typing/epithelial/patient_condition_mean_GSVA.rds")

kegg_patient_condition_mean <- readRDS("processed_data/data_B2-19/cell_typing/epithelial/patient_condition_mean_GSVA.rds")

# by patient condition
options(repr.plot.width=6, repr.plot.height=3)
p_patient_condition_mean<-pheatmap(kegg_patient_condition_mean[c("HALLMARK_ANGIOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS"),], 
                                   show_rownames=1, show_colnames=T,
                                   cluster_cols = FALSE,cluster_rows = FALSE)

epithelial.batch.corrected$Lauren.LN_condition <- paste(epithelial.batch.corrected$Lauren, epithelial.batch.corrected$LN_condition, sep = "_")
table(epithelial.batch.corrected$Lauren.LN_condition)

Lauren.LN_condition_mean=data.frame(row.names = rownames(epithelial.batch.corrected))
for(i in c("Diffuse_Pri.GC", "Intestinal_Pri.GC")){
    sub_i=subset(epithelial.batch.corrected, subset=Lauren.LN_condition==i)
    Lauren.LN_condition_mean[,as.character(i)]=apply(sub_i@assays$RNA@data, 1, mean)
}

kegg_Lauren.LN_condition_mean <- gsva(as.matrix(Lauren.LN_condition_mean), msigdbr_list, kcdf="Gaussian",method = "gsva",parallel.sz=10)


# saveRDS(kegg_Lauren.LN_condition_mean,"processed_data/data_B2-19/cell_typing/epithelial/Lauren.LN_condition_mean_GSVA.rds")

kegg_Lauren.LN_condition_mean <- readRDS("processed_data/data_B2-19/cell_typing/epithelial/Lauren.LN_condition_mean_GSVA.rds")

head(kegg_Lauren.LN_condition_mean)

custom_colors <- rev(c("#4792c4", "#9cc2d7", "#f5f4e9", "#f7ddc6", "#e7a988", "#c86955"))
dotplot(epi_go,showCategory=3)+ 
  scale_color_gradientn(colors = custom_colors) +

head(kegg_Lauren.LN_condition_mean)

# 缩放数据到 -1 到 1 范围
kegg_Lauren.LN_condition_mean <- apply(kegg_Lauren.LN_condition_mean, 2, function(x) pmax(pmin(x, 1), -1))

# 自定义颜色断点和颜色
my.breaks <- c(seq(-1, 0, by = 0.1), seq(0.1, 1, by = 0.1))
# my.colors <- c(
#     colorRampPalette(colors = c("#89afd2", "white"))(length(my.breaks) / 2),
#     colorRampPalette(colors = c("white", "#ef8775"))(length(my.breaks) / 2)
# )
custom_colors <- colorRampPalette(colors = c("#4792c4", "#fdf8b4", "#c86955"))(20)

# 调整图表宽度和高度
options(repr.plot.width = 6, repr.plot.height = 3)

# 绘制热图
p_Lauren.LN_condition_mean <- pheatmap(
    rescale(kegg_Lauren.LN_condition_mean, to=c(-1, 1))[c("HALLMARK_ANGIOGENESIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                                    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS"), ], 
    show_rownames = 1, show_colnames = T,
    cluster_cols = FALSE, cluster_rows = FALSE,
    breaks = my.breaks,
color = custom_colors)
  
    
    # color = my.colors


# by patient condition
options(repr.plot.width=6, repr.plot.height=3)
p_Lauren.LN_condition_mean<-pheatmap(kegg_Lauren.LN_condition_mean[c("HALLMARK_ANGIOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS"),], 
                                   show_rownames=1, show_colnames=T,
                                   cluster_cols = FALSE,cluster_rows = FALSE)

library(scales)
kegg_Lauren.LN_condition_mean <- t(apply(kegg_Lauren.LN_condition_mean, 1, rescale, to=c(-1, 1)))

rescale(kegg_Lauren.LN_condition_mean, to=c(-1, 1))

library(scales)
kegg_Lauren.LN_condition_mean <- t(apply(kegg_Lauren.LN_condition_mean, 1, rescale, to=c(-1, 1)))

# by patient condition
options(repr.plot.width=6, repr.plot.height=3)
p_Lauren.LN_condition_mean<-pheatmap(kegg_Lauren.LN_condition_mean[c("HALLMARK_ANGIOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS"),], 
                                   show_rownames=1, show_colnames=T,
                                   cluster_cols = FALSE,cluster_rows = FALSE)

# by patient condition
options(repr.plot.width=6, repr.plot.height=3)
p_Lauren.LN_condition_mean<-pheatmap(kegg_Lauren.LN_condition_mean[c("HALLMARK_ANGIOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS"),], 
                                   show_rownames=1, show_colnames=T,
                                   cluster_cols = FALSE,cluster_rows = FALSE)

p_Lauren.LN_condition_mean

# ggsave(p_seurat_clusters,
#        filename = "pathway_results/GSVA_H_pathways_for_seurat_clusters_malign.pdf",width=12,height=9)



myeloid_tmp <- myeloid.batch.corrected

DefaultAssay(myeloid_tmp) <- "integrated"

# subset big subtype to 500 cells
set.seed=0
cell_use <- NULL
for(subtype in unique(myeloid_tmp$subtype1)){
    print(subtype)
    # print(sum(myeloid_tmp$subtype1==subtype))
    # if(sum(myeloid_tmp$subtype1==subtype) > 500){
        rowname <- sample(colnames(myeloid_tmp)[myeloid_tmp$subtype1==subtype], 10, replace = FALSE)
    # }else{
    #     rowname <- colnames(myeloid_tmp)[myeloid_tmp$subtype1==subtype]
    # }
    cell_use <- c(cell_use,rowname)
    # print(length(cell_use))
}

myeloid_tmp_use <- subset(myeloid_tmp,cells=cell_use)
dim(myeloid_tmp_use)

tmp_myeloid_markers <- FindAllMarkers(myeloid_tmp_use, only.pos = TRUE)

tmp_myeloid_markers <- tmp_myeloid_markers %>%
    group_by(cluster) %>%
    slice_max(n = 4, order_by = avg_log2FC)

dim(tmp_myeloid_markers)

options(repr.plot.width=8, repr.plot.height=3)
my.colors <- c(
    colorRampPalette(colors = c("#89afd2", "white"))(10),
    colorRampPalette(colors = c("white", "#ef8775"))(10))

DoHeatmap(object = myeloid_tmp_use,
          features = tmp_myeloid_markers$gene,
          group.bar = F,draw.lines = F)+
  scale_fill_gradientn(colors = c("#89afd2", "white", "#ef8775"))

options(repr.plot.width=6, repr.plot.height=5)
n=length(unique(GC$cell_type))
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = DiscretePalette(n, palette = "alphabet", shuffle = FALSE))
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
n=length(unique(GC$cell_type))
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = DiscretePalette(n, palette = "alphabet2", shuffle = FALSE))
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
n=length(unique(GC$cell_type))
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = DiscretePalette(n, palette = "glasbey", shuffle = FALSE))
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
n=length(unique(GC$cell_type))
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = DiscretePalette(n, palette = "polychrome", shuffle = FALSE))
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
n=length(unique(GC$cell_type))
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = DiscretePalette(n, palette = "stepped", shuffle = FALSE))
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
n=length(unique(GC$cell_type))
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = DiscretePalette(n, palette = "parade", shuffle = FALSE))
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = "Set3")
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = "Set2")
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = "Paired")
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 

options(repr.plot.width=6, repr.plot.height=5)
DimPlot(GC, reduction = "umap", shuffle=TRUE, group.by = "cell_type",label=F,
        cols = "Dark2")
# "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade". 



# # 安装并加载必要的包
# install.packages("ggplot2")
# install.packages("scales")

library(ggplot2)
library(scales)

# 创建示例数据集
set.seed(42)
data <- data.frame(
  Value = rnorm(100)
)

# 自定义渐变色函数
low_high_gradient <- scale_color_gradient2(low = "#89afd2", mid = "yellow", high = "#ef8775", midpoint = 0)

# 绘制散点图
ggplot(data, aes(x = Value, y = rnorm(100), color = Value)) +
  geom_point(size = 3) +
  low_high_gradient +
  ggtitle("Scatter Plot with Custom Low Saturation Gradient") +
  theme_minimal()


# # 安装并加载必要的包
# install.packages("ggplot2")
# install.packages("scales")

library(ggplot2)
library(scales)

# 创建示例数据集
set.seed(42)
data <- data.frame(
  Value = rnorm(100)
)

# 自定义渐变色函数
low_high_gradient <- scale_color_gradient2(low = "#89afd2", high = "#ef8775")

# 绘制散点图
ggplot(data, aes(x = Value, y = rnorm(100), color = Value)) +
  geom_point(size = 3) +
  low_high_gradient +
  ggtitle("Scatter Plot with Custom Low Saturation Gradient") +
  theme_minimal()

