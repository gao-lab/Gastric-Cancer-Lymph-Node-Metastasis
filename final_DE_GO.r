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

packageVersion("Seurat")

packageVersion("SeuratObject")

meta=read.csv("processed_data/data_B2-19/metadata_all_ct_subtype_condition_noTLB.csv", row.names = 1)


tcell.batch.corrected <- readRDS("plots/data/tcell_no_TLB_new_meta.rds")
Idents(tcell.batch.corrected) <- "subtype1"
DefaultAssay(tcell.batch.corrected)

myeloid.batch.corrected <- readRDS("plots/data/myeloid_new_meta.rds")
Idents(myeloid.batch.corrected) <- "subtype1"
DefaultAssay(myeloid.batch.corrected) <- "RNA"


bcell.batch.corrected <- readRDS("plots/data/bcell_no_TLB_new_meta.rds")
DefaultAssay(bcell.batch.corrected) <- "RNA"
Idents(bcell.batch.corrected) <- "subtype1"


epithelial.batch.corrected <- readRDS("plots/data/epithelial_new_meta.rds")



up_GO_enrichment <- function(markers){
    up_markers <- markers %>%
          filter(avg_log2FC>0)
    if(dim(up_markers)[1]>10){ # number of DE markers >10, then do GO enrichment
        up_transID=bitr(up_markers$gene,
                         fromType="SYMBOL",
                         toType=c("ENSEMBL", "ENTREZID"),
                         OrgDb="org.Hs.eg.db"
        )
        up_BP=enrichGO(up_transID$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP"
        )
        up_MF=enrichGO(up_transID$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF"
        )
        up_CC=enrichGO(up_transID$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC"
        )
        up_KEGG=enrichKEGG(up_transID$ENTREZID,
                        "hsa"
        )}
    return(list(up_BP = up_BP, 
                up_MF = up_MF, 
                up_CC = up_CC,
                up_KEGG = up_KEGG
               ))
}

down_GO_enrichment <- function(markers){       
    down_markers <- markers %>%
          filter(avg_log2FC<0)
    if(dim(down_markers)[1]>10){
        down_transID=bitr(down_markers$gene,
                         fromType="SYMBOL",
                         toType=c("ENSEMBL", "ENTREZID"),
                         OrgDb="org.Hs.eg.db"
        )
        down_BP=enrichGO(down_transID$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP"
        )
        down_MF=enrichGO(down_transID$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF"
        )
        down_CC=enrichGO(down_transID$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC"
        )
        down_KEGG=enrichKEGG(down_transID$ENTREZID,
                        "hsa"
        )
    }
    return(list(down_BP = down_BP, 
                down_MF = down_MF, 
                down_CC = down_CC,
                down_KEGG = down_KEGG
               ))
}

tcell.batch.corrected$cell_type_tmp <- tcell.batch.corrected$subtype1
tcell.batch.corrected$cell_type_tmp[tcell.batch.corrected$cell_type_tmp %in% c("CD8_CCL5+ Tex",
                                                                                       "CD8_LAYN+ Tex")] <- "Tex"
tcell.batch.corrected$cell_type_tmp[tcell.batch.corrected$cell_type_tmp == "CD8_HLA-DRA+ Tem"] <- "HLA-DRA+ Tem"

Idents(tcell.batch.corrected) <- "cell_type_tmp"


table(Idents(tcell.batch.corrected)[tcell.batch.corrected$cell_type_tmp %in% c("Tex", "HLA-DRA+ Tem")])

markers_HLA_DRA_Tem_Tex=FindMarkers(tcell.batch.corrected, 
                                                               ident.1 = "HLA-DRA+ Tem", 
                                                               ident.2 = "Tex", 
                                                               verbose = FALSE)
markers_HLA_DRA_Tem_Tex$gene <- rownames(markers_HLA_DRA_Tem_Tex)

markers_HLA_DRA_Tem_Tex_0 <- markers_HLA_DRA_Tem_Tex

markers_HLA_DRA_Tem_Tex <- markers_HLA_DRA_Tem_Tex_0[markers_HLA_DRA_Tem_Tex_0$p_val_adj < 0.05,]
write.csv(markers_HLA_DRA_Tem_Tex,
          "processed_data/data_B2-19/cell_typing/tcell/HLA_DRA_Tem/markers_HLA_DRA_Tem_Tex.csv")

mask <- markers_HLA_DRA_Tem_Tex_0$p_val_adj < 0.05 & abs(markers_HLA_DRA_Tem_Tex_0$avg_log2FC) > 0.5
table(mask)
markers_HLA_DRA_Tem_Tex <- markers_HLA_DRA_Tem_Tex_0[mask,]
table(markers_HLA_DRA_Tem_Tex$avg_log2FC>0)

# markers_HLA_DRA_Tem_Tex <- read.csv("processed_data/data_B2-19/cell_typing/tcell/HLA_DRA_Tem/markers_HLA_DRA_Tem_Tex.csv",
#                                    row.names=1)

dim(markers_HLA_DRA_Tem_Tex)

epi_markers_symbol <- list()
up_markers <- markers_HLA_DRA_Tem_Tex %>%
          filter(avg_log2FC>0) %>% pull(gene)
epi_markers_symbol[["HLA-DRA+ Tem"]] <- up_markers
down_markers <- markers_HLA_DRA_Tem_Tex %>%
          filter(avg_log2FC<0) %>% pull(gene)
epi_markers_symbol[["Tex"]] <- down_markers


epi_markers_id <- list()
for(condition in names(epi_markers_symbol)){
     tran_id <- bitr(epi_markers_symbol[[condition]],
                 fromType="SYMBOL",
                 toType=c("ENSEMBL", "ENTREZID"),
                 OrgDb="org.Hs.eg.db"
    )
    epi_markers_id[[condition]] <- tran_id$ENTREZID
}

ck_all <- compareCluster(geneCluster = epi_markers_id, fun = enrichGO, OrgDb='org.Hs.eg.db',ont = "BP")
ck_all <- setReadable(ck_all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# head(ck_all) 

saveRDS(ck_all,"processed_data/data_B2-19/cell_typing/tcell/HLA_DRA_Tem/GO_HLA_DRA_Tem_Tex.rds")

options(repr.plot.width=7, repr.plot.height=6)

dotplot(ck_all,showCategory=c("positive regulation of cell activation",
                              "negative regulation of cell activation",
                              "positive regulation of leukocyte activation",
                              "negative regulation of leukocyte activation",
                              "positive regulation of T cell activation",
                              "negative regulation of T cell activation"
                             ))#+  #,by="p.adjust"
#     theme(axis.text.x = element_text(angle = 45, hjust=1))

epithelial.batch.corrected$patient_condition.LN_condition <- paste(epithelial.batch.corrected$patient_condition, 
                                                                   epithelial.batch.corrected$LN_condition, sep = "_")
table(meta$patient_condition.LN_condition)

epithelial.batch.corrected$Lauren.LN_condition <- paste(epithelial.batch.corrected$Lauren, epithelial.batch.corrected$LN_condition, sep = "_")
table(epithelial.batch.corrected$Lauren.LN_condition)

table(epithelial.batch.corrected$LN_condition)

Idents(epithelial.batch.corrected) <- "LN_condition"

markers_epi_LN_condition=FindMarkers(epithelial.batch.corrected, 
                                                               ident.1 = "Met.LN", 
                                                               ident.2 = "Pri.GC", 
                                                               verbose = FALSE)
markers_epi_LN_condition$gene <- rownames(markers_epi_LN_condition)

markers_epi_LN_condition_0 <- markers_epi_LN_condition

Idents(epithelial.batch.corrected) <- "patient_condition.LN_condition"

markers_epi_patient_condition.LN_condition=FindMarkers(epithelial.batch.corrected, 
                                                               ident.1 = "Meta_Pri.GC", 
                                                               ident.2 = "Neg_Pri.GC", 
                                                               verbose = FALSE)
markers_epi_patient_condition.LN_condition$gene <- rownames(markers_epi_patient_condition.LN_condition)

markers_epi_patient_condition.LN_condition_0 <- markers_epi_patient_condition.LN_condition

Idents(epithelial.batch.corrected) <- "Lauren.LN_condition"

markers_epi_Lauren.LN_condition=FindMarkers(epithelial.batch.corrected, 
                                                               ident.1 = "Diffuse_Pri.GC", 
                                                               ident.2 = "Intestinal_Pri.GC", 
                                                               verbose = FALSE)
markers_epi_Lauren.LN_condition$gene <- rownames(markers_epi_Lauren.LN_condition)

markers_epi_Lauren.LN_condition_0 <- markers_epi_Lauren.LN_condition

mask <- markers_epi_LN_condition_0$p_val_adj < 0.05 & abs(markers_epi_LN_condition_0$avg_log2FC) > 0.5
table(mask)
markers_epi_LN_condition <- markers_epi_LN_condition_0[mask,]
table(markers_epi_LN_condition$avg_log2FC>0)

mask <- markers_epi_patient_condition.LN_condition_0$p_val_adj < 0.05 & abs(markers_epi_patient_condition.LN_condition_0$avg_log2FC) > 0.5
table(mask)
markers_epi_patient_condition.LN_condition <- markers_epi_patient_condition.LN_condition_0[mask,]
table(markers_epi_patient_condition.LN_condition$avg_log2FC>0)

mask <- markers_epi_Lauren.LN_condition_0$p_val_adj < 0.05 & abs(markers_epi_Lauren.LN_condition_0$avg_log2FC) > 0.5
table(mask)
markers_epi_Lauren.LN_condition <- markers_epi_Lauren.LN_condition_0[mask,]
table(markers_epi_Lauren.LN_condition$avg_log2FC>0)

epi_de_markers <- list()

epi_de_markers[["T.patient_condition"]] <- markers_epi_patient_condition.LN_condition
epi_de_markers[["T.Lauren"]] <- markers_epi_Lauren.LN_condition
epi_de_markers[["T_LN_meta"]] <- markers_epi_LN_condition


saveRDS(epi_de_markers,
        "processed_data/data_B2-19/cell_typing/epithelial/epi_de_markers_list.rds")



# epi_de_markers <- readRDS("processed_data/data_B2-19/cell_typing/epithelial/epi_de_markers_list.rds")

epi_markers_symbol <- list()
up_markers <- epi_de_markers[["T.patient_condition"]] %>%
          filter(avg_log2FC>0) %>% pull(gene)
epi_markers_symbol$Meta_Pri.GC <- up_markers
down_markers <- epi_de_markers[["T.patient_condition"]] %>%
          filter(avg_log2FC<0) %>% pull(gene)
epi_markers_symbol$Neg_Pri.GC <- down_markers

up_markers <- epi_de_markers[["T_LN_meta"]] %>%
          filter(avg_log2FC>0) %>% pull(gene)
epi_markers_symbol$Met.LN <- up_markers
down_markers <- epi_de_markers[["T_LN_meta"]] %>%
          filter(avg_log2FC<0) %>% pull(gene)
epi_markers_symbol$Pri.GC <- down_markers

up_markers <- epi_de_markers[["T.Lauren"]] %>%
          filter(avg_log2FC>0) %>% pull(gene)
epi_markers_symbol$Diffuse_Pri.GC <- up_markers

down_markers <- epi_de_markers[["T.Lauren"]] %>%
          filter(avg_log2FC<0) %>% pull(gene)
epi_markers_symbol$Intestinal_Pri.GC <- down_markers


names(epi_markers_symbol)

epi_markers_id <- list()
for(condition in names(epi_markers_symbol)){
     tran_id <- bitr(epi_markers_symbol[[condition]],
                 fromType="SYMBOL",
                 toType=c("ENSEMBL", "ENTREZID"),
                 OrgDb="org.Hs.eg.db"
    )
    epi_markers_id[[condition]] <- tran_id$ENTREZID
}

ck_all <- compareCluster(geneCluster = epi_markers_id, fun = enrichGO, OrgDb='org.Hs.eg.db',ont = "BP")
ck_all <- setReadable(ck_all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# head(ck_all) 

saveRDS(ck_all,
        "processed_data/data_B2-19/cell_typing/epithelial/GO_all_conditions.rds")

options(repr.plot.width=7, repr.plot.height=8)
custom_colors <- rev(c("#4792c4", "#9cc2d7", "#f5f4e9", "#f7ddc6", "#e7a988", "#c86955"))
go_epi <- dotplot(ck_all,showCategory=3)+ 
  scale_color_gradientn(colors = custom_colors) +
#   scale_color_gradient(low = "#bb5a74", high = "#bfdaf3")+
    theme(axis.text.x = element_text(angle = 45, hjust=1))#+font_theme
plot(go_epi)
# pdf("plots/figures/fig1F_epi_go.pdf", width = 7, height = 8)
# plot(go_epi)
# dev.off()


Idents(tcell.batch.corrected) <- "subtype1"
unique(Idents(tcell.batch.corrected))

markers_GZMB_Temra_GZMK_Temra=FindMarkers(tcell.batch.corrected, 
                                                               ident.1 = "CD8_GZMB+ Temra", 
                                                               ident.2 = "CD8_GZMK+ Temra", 
                                                               verbose = FALSE)
markers_GZMB_Temra_GZMK_Temra$gene <- rownames(markers_GZMB_Temra_GZMK_Temra)
markers_GZMB_Temra_GZMK_Temra_0 <- markers_GZMB_Temra_GZMK_Temra


markers_GZMB_Temra_GZMK_Temra <- markers_GZMB_Temra_GZMK_Temra[markers_GZMB_Temra_GZMK_Temra$p_val_adj < 0.05,]
write.csv(markers_GZMB_Temra_GZMK_Temra,
          "processed_data/data_B2-19/cell_typing/tcell/GZMB_Temra/markers_GZMB_Temra_GZMK_Temra.csv")

mask <- markers_GZMB_Temra_GZMK_Temra_0$p_val_adj < 0.05 & abs(markers_GZMB_Temra_GZMK_Temra_0$avg_log2FC) > 0.5
table(mask)
markers_GZMB_Temra_GZMK_Temra <- markers_GZMB_Temra_GZMK_Temra_0[mask,]
table(markers_GZMB_Temra_GZMK_Temra$avg_log2FC>0)

GO_GZMB_Temra_GZMK_Temra <- GO_enrichment(markers_GZMB_Temra_GZMK_Temra)

options(repr.plot.width=6, repr.plot.height=6)
plot_GO(GO_GZMB_Temra_GZMK_Temra,"GZMB+ Temra to GZMK+ Temra",10)

saveRDS(GO_GZMB_Temra_GZMK_Temra,"processed_data/data_B2-19/cell_typing/tcell/GZMB_Temra/GO_GZMB_Temra_GZMK_Temra.rds")

DC_so <- subset(myeloid.batch.corrected, 
                    subtype1 %in% grep("DC",unique(myeloid.batch.corrected$subtype1),value = TRUE))
DefaultAssay(DC_so)

table(Idents(DC_so))

markers_SPIB_DC <- FindMarkers(DC_so, `ident.1` = "SPIB+ DC")

markers_SPIB_DC$gene <- rownames(markers_SPIB_DC)
markers_SPIB_DC_0 <- markers_SPIB_DC


markers_SPIB_DC <- markers_SPIB_DC_0[markers_SPIB_DC_0$p_val_adj < 0.05,]
write.csv(markers_SPIB_DC,
          "processed_data/data_B2-19/cell_typing/DC/SPIB_DC/markers_SPIB_DC.csv")

mask <- markers_SPIB_DC_0$p_val_adj < 0.05 & abs(markers_SPIB_DC_0$avg_log2FC) > 0.5 & markers_SPIB_DC_0$`pct.1`>0.25
table(mask)
markers_SPIB_DC <- markers_SPIB_DC_0[mask,]
table(markers_SPIB_DC$avg_log2FC>0)

GO_SPIB_DC <- GO_enrichment(markers_SPIB_DC)

clusterProfiler::dotplot(GO_SPIB_DC$up_BP,showCategory=8)

options(repr.plot.width=6, repr.plot.height=6)
plot_GO(GO_SPIB_DC,"SPIB+ DC to DC",10)

saveRDS(GO_SPIB_DC, "processed_data/data_B2-19/cell_typing/DC/SPIB_DC/GO_SPIB_DC_DC.rds")

DC_so <- subset(myeloid.batch.corrected, 
                    subtype1 %in% grep("DC",unique(myeloid.batch.corrected$subtype1),value = TRUE))
DefaultAssay(DC_so)
Idents(DC_so) <- "subtype1"


markers_cDC2_markers <- FindMarkers(DC_so, ident.1 = "FCER1A+ cDC2", ident.2 = "NLRP3+ cDC2", verbose = FALSE)

markers_cDC2_markers$gene <- rownames(markers_cDC2_markers)
markers_cDC2_markers_0 <- markers_cDC2_markers


markers_cDC2_markers <- markers_cDC2_markers_0[markers_cDC2_markers_0$p_val_adj < 0.05,]
write.csv(markers_cDC2_markers,"processed_data/data_B2-19/cell_typing/DC/FCER1A_cDC2/FCER1A_cDC2_markers.csv")


mask <- markers_cDC2_markers_0$p_val_adj < 0.05 & abs(markers_cDC2_markers_0$avg_log2FC) > 0.5
table(mask)
markers_cDC2_markers <- markers_cDC2_markers_0[mask,]
table(markers_cDC2_markers$avg_log2FC>0)



GO_cDC2_markers <- GO_enrichment(markers_cDC2_markers)

options(repr.plot.width=6, repr.plot.height=6)
plot_GO(GO_cDC2_markers,"FCER1A+ cDC2 to NLRP3+ cDC2",15)

saveRDS(GO_cDC2_markers, "processed_data/data_B2-19/cell_typing/DC/FCER1A_cDC2/GO_FCER1A_cDC2.rds")

plasma.batch.corrected <- readRDS("processed_data/data_B2-19/cell_typing/bcell/plasma.batch.corrected_ct_0.5.rds")

table((plasma.batch.corrected$cell_type))

plasma.batch.corrected$subtype0 <- "Plasma"
plasma.batch.corrected$subtype0[plasma.batch.corrected$cell_type %in% c("STMN1+ Plasmablast","MS4A1+ Plasmablast")] <- "Plasmablast"
plasma.batch.corrected$subtype0[plasma.batch.corrected$cell_type == "RRM2+ Plasma"] <- "RRM2+ Plasma"


Idents(plasma.batch.corrected) <- "subtype0"

table(Idents(plasma.batch.corrected))

markers_RRM2_Plasma=FindMarkers(plasma.batch.corrected, 
                                                               ident.1 = "RRM2+ Plasma", 
                                                               ident.2 = "Plasma", 
                                                               verbose = FALSE)
markers_RRM2_Plasma$gene <- rownames(markers_RRM2_Plasma)
markers_RRM2_Plasma_0 <- markers_RRM2_Plasma


markers_RRM2_Plasma <- markers_RRM2_Plasma_0[markers_RRM2_Plasma_0$p_val_adj < 0.05,]
write.csv(markers_RRM2_Plasma,"processed_data/data_B2-19/cell_typing/bcell/RRM2_Plasma/markers_RRM2_Plasma_Plasma.csv")


mask <- markers_RRM2_Plasma_0$p_val_adj < 0.05 & abs(markers_RRM2_Plasma_0$avg_log2FC) > 0.5
table(mask)
markers_RRM2_Plasma <- markers_RRM2_Plasma_0[mask,]
table(markers_RRM2_Plasma$avg_log2FC>0)

GO_RRM2_Plasma <- GO_enrichment(markers_RRM2_Plasma)

options(repr.plot.width=6, repr.plot.height=6)
plot_GO(GO_RRM2_Plasma,"RRM2+ Plasma to Plasma",10)

saveRDS(GO_RRM2_Plasma, "processed_data/data_B2-19/cell_typing/bcell/RRM2_Plasma/GO_RRM2_Plasma_Plasma.rds")

markers_RRM2_Plasmablast=FindMarkers(plasma.batch.corrected, 
                                                               ident.1 = "RRM2+ Plasma", 
                                                               ident.2 = "Plasmablast", 
                                                               verbose = FALSE)
markers_RRM2_Plasmablast$gene <- rownames(markers_RRM2_Plasmablast)
markers_RRM2_Plasmablast_0 <- markers_RRM2_Plasmablast


markers_RRM2_Plasmablast <- markers_RRM2_Plasmablast_0[markers_RRM2_Plasmablast_0$p_val_adj < 0.05,]
write.csv(markers_RRM2_Plasma,"processed_data/data_B2-19/cell_typing/bcell/RRM2_Plasma/markers_RRM2_Plasma_Plasmablast.csv")

mask <- markers_RRM2_Plasmablast_0$p_val_adj < 0.05 & abs(markers_RRM2_Plasmablast_0$avg_log2FC) > 0.5
table(mask)
markers_RRM2_Plasmablast <- markers_RRM2_Plasmablast_0[mask,]
table(markers_RRM2_Plasmablast$avg_log2FC>0)

GO_RRM2_Plasmablast <- GO_enrichment(markers_RRM2_Plasmablast)

options(repr.plot.width=6, repr.plot.height=5)
p <- barplot(GO_RRM2_Plasmablast$up_BP, showCategory=9, title="RRM2+ Plasma to Plasmablast")
plot(p)

saveRDS(GO_RRM2_Plasmablast, "processed_data/data_B2-19/cell_typing/bcell/RRM2_Plasma/GO_RRM2_Plasma_Plasmablast.rds")

epi_markers_symbol <- list()
up_markers <- markers_RRM2_Plasmablast %>%
          filter(avg_log2FC>0) %>% pull(gene)
epi_markers_symbol[["RRM2+ Plasma to Plasmablast"]] <- up_markers
down_markers <- markers_RRM2_Plasma %>%
          filter(avg_log2FC<0) %>% pull(gene)
epi_markers_symbol[["Plasma to RRM2+ Plasma"]] <- down_markers



epi_markers_symbol

epi_markers_id <- list()
for(condition in names(epi_markers_symbol)){
     tran_id <- bitr(epi_markers_symbol[[condition]],
                 fromType="SYMBOL",
                 toType=c("ENSEMBL", "ENTREZID"),
                 OrgDb="org.Hs.eg.db"
    )
    epi_markers_id[[condition]] <- tran_id$ENTREZID
}

ck_all_0.5 <- compareCluster(geneCluster = epi_markers_id, fun = enrichGO, OrgDb='org.Hs.eg.db',ont = "CC")
#不说明ont则默认应用MF
ck_all_0.5 <- setReadable(ck_all_0.5, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# head(ck_all_0.5) 

saveRDS(ck_all_0.5,"processed_data/data_B2-19/cell_typing/bcell/RRM2_Plasma/go_compare_RRM2_Plasma.RDS")

options(repr.plot.width=6, repr.plot.height=7)

dotplot(`ck_all_0.5`,showCategory=c("immunoglobulin complex",
                                    "endoplasmic reticulum protein-containing complex",
                              "IgG immunoglobulin complex",
                              "IgA immunoglobulin complex",
                              "immunoglobulin complex, circulating"
                             ))+  
    theme(axis.text.x = element_text(angle = 30, hjust=1))

options(repr.plot.width=7, repr.plot.height=9)

dotplot(ck_all_0.5,showCategory=5)+  #,by="p.adjust"
    theme(axis.text.x = element_text(angle = 45, hjust=1))



lamp3_dc <- subset(myeloid.batch.corrected, subtype1=="LAMP3+ DC")

lamp3_dc$Lauren.LN_condition <- paste(lamp3_dc$Lauren, lamp3_dc$LN_condition, sep = "_")
table(lamp3_dc$Lauren.LN_condition)
Idents(lamp3_dc) <- "Lauren.LN_condition"


lamp3_de_markers=FindMarkers(lamp3_dc, 
                                    ident.1 = "Diffuse_Met.LN", 
                                    ident.2 = "Intestinal_Met.LN", 
                                    verbose = FALSE)
lamp3_de_markers$gene <- rownames(lamp3_de_markers)

lamp3_de_markers_0 <- lamp3_de_markers

# lamp3_de_markers <- lamp3_de_markers_0[lamp3_de_markers_0$p_val_adj < 0.05,]
write.csv(lamp3_de_markers_0,
          "processed_data/data_B2-19/cell_typing/DC/LAMP3_DC/lamp3_lauren_de_markers.csv")

meta_plasma <- meta[colnames(plasma.batch.corrected),]

all(colnames(plasma.batch.corrected)==rownames(meta_plasma))

plasma.batch.corrected$LN_condition <- meta_plasma$LN_condition
plasma.batch.corrected$patient_condition <- meta_plasma$patient_condition
plasma.batch.corrected$Lauren <- meta_plasma$Lauren

Idents(plasma.batch.corrected) <- "LN_condition"
table(Idents(plasma.batch.corrected))

plasma_de_markers=FindMarkers(plasma.batch.corrected, 
                                    ident.1 = "Pri.GC", 
                                    ident.2 = "Met.LN", 
                                    verbose = FALSE)
plasma_de_markers$gene <- rownames(plasma_de_markers)

plasma_de_markers_0 <- plasma_de_markers

# plasma_de_markers <- plasma_de_markers_0[plasma_de_markers_0$p_val_adj < 0.05,]
write.csv(plasma_de_markers_0,
          "processed_data/data_B2-19/cell_typing/bcell/Plasma/plasma_LN_condition_de_markers.csv")




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



library(tidyverse)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
# library(wesanderson)
library(viridis)
library(ggsci)
# library(tidytext)
library(ggpubr)
library(cowplot)
# library(facetscales)
# library(latex2exp)
# library(ggstatsplot)
library(scales)

cd8_tcell.batch.corrected <- readRDS("processed_data/data_B2-19/cell_typing/tcell/cd8_tcell.batch_corrected.ct.rds")


unique(cd8_tcell.batch.corrected$subtype1)

library(openxlsx)
signatures <- read.xlsx("public_markers/t cell signatures.xlsx",sheet =5,startRow=2,colNames = TRUE,rowNames=FALSE)

marker.list <- list()
for(pathway in colnames(signatures)){
    print(pathway)
    marker.list[[pathway]] <- na.omit(as.vector(signatures[,pathway]))
    print(length(marker.list[[pathway]]))
}

CD8_Obj <- cd8_tcell.batch.corrected
CD8_Obj <-  AddModuleScore(cd8_tcell.batch.corrected,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")

Differentiation <- c("Naive", "Activation:Effector function", "Exhaustion")
Function <- c("TCR Signaling", "Cytotoxicity", "Cytokine:Cytokine receptor",
              "Chemokine:Chemokine receptor", "Senescence", "Anergy",
              "NFKB Signaling", "Stress response", "MAPK Signaling", "Adhesion",
              "IFN Response")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)


MarkerNameVector
names(marker.list)

for(i in 1:length(marker.list)){
    colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- MarkerNameVector[i]
}
Idents(CD8_Obj) <- CD8_Obj$subtype1

FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(CD8_Obj$subtype1)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- unique(CD8_Obj$subtype1)
rownames(FunctionScoreMatrix) <- MarkerNameVector

for(ci in colnames(FunctionScoreMatrix)){
    for(ri in 1:nrow(FunctionScoreMatrix)){
        FunctionVec <- as_tibble(CD8_Obj@meta.data) %>% pull(MarkerNameVector[ri])
        print(head(FunctionVec))
        fv <- mean(FunctionVec[CD8_Obj$subtype1 == ci],na.rm = TRUE)
        print(fv)
        FunctionScoreMatrix[ri, ci] <- fv
    }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))


saveRDS(FunctionScoreMatrix,"processed_data/data_B2-19/cell_typing/tcell/cd8_FunctionScoreMatrix.rds")


signatureType_row <- data.frame(Signature.type = c(
                                    rep("Differentiation", length(Differentiation)),
                                    rep("Function", length(Function)),
                                    rep("Metabolism", length(Metabolism)),
                                    rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector

saveRDS(signatureType_row,"processed_data/data_B2-19/cell_typing/tcell/cd8_signatureType_row.rds")


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

options(repr.plot.width=8, repr.plot.height=5)
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(3, 14, 17),
         cluster_rows = F,
         cluster_cols = TRUE,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8)

# orderC = c("CD8_c13", "CD8_c3", "CD8_c6", "CD8_c0", "CD8_c11", "CD8_c9", "CD8_c10", "CD8_c12", "CD8_c8", "CD8_c2", "CD8_c7", "CD8_c4", "CD8_c5", "CD8_c1")
# FunctionScoreMatrix <- FunctionScoreMatrix[,orderC]
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)
signatureType_row <- data.frame(Signature.type = c(
                                    rep("Differentiation", length(Differentiation)),
                                    rep("Function", length(Function)),
                                    rep("Metabolism", length(Metabolism)),
                                    rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector
pheatmap(FunctionScoreMatrix,
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

options(repr.plot.width=8, repr.plot.height=5)
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(3, 14, 17),
         cluster_rows = F,
         cluster_cols = TRUE,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8)
