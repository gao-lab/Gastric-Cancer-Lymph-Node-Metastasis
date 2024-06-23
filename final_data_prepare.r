library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsignif)
library(scales)


packageVersion("Seurat")

packageVersion("SeuratObject")

meta=read.csv("processed_data/data_B2-19/metadata_all_ct_subtype_condition_noTLB.csv", row.names = 1)


tcell.batch.corrected <- readRDS("processed_data/data_B2-19/cell_typing/tcell/tcell.batch.corrected_ct_subtype.rds")


unique(tcell.batch.corrected$subtype1)

mask <- tcell.batch.corrected$subtype1 %in% c("TLB",grep("Low",unique(tcell.batch.corrected$subtype1),value = TRUE))
table(mask)
tcell.batch.corrected <- tcell.batch.corrected[,!mask]
dim(tcell.batch.corrected)

table(tcell.batch.corrected$subtype1[tcell.batch.corrected$cell_type1=="Tgd"])

meta_tcell <- meta[colnames(tcell.batch.corrected),]

all(colnames(tcell.batch.corrected)==rownames(meta_tcell))

mask <- tcell.batch.corrected$subtype1 != meta_tcell$subtype1
tcell.batch.corrected$subtype1[mask] <- meta_tcell$subtype1[mask]

tcell.batch.corrected$LN_condition <- meta_tcell$LN_condition
tcell.batch.corrected$patient_condition <- meta_tcell$patient_condition
tcell.batch.corrected$Lauren <- meta_tcell$Lauren

Idents(tcell.batch.corrected) <- "subtype1"

saveRDS(tcell.batch.corrected,"plots/data/tcell_no_TLB_new_meta.rds")



cd4_tcell.batch.corrected <- readRDS("processed_data/data_B2-19/cell_typing/tcell/cd4_tcell.batch_corrected.ct.rds")

meta_cd4_tcell <- meta[colnames(cd4_tcell.batch.corrected),]
all(colnames(cd4_tcell.batch.corrected)==rownames(meta_cd4_tcell))

cd4_tcell.batch.corrected$LN_condition <- meta_cd4_tcell$LN_condition
cd4_tcell.batch.corrected$patient_condition <- meta_cd4_tcell$patient_condition
cd4_tcell.batch.corrected$Lauren <- meta_cd4_tcell$Lauren

saveRDS(cd4_tcell.batch.corrected,"plots/data/cd4_tcell_new_meta.rds")

myeloid.batch.corrected <- readRDS("processed_data/data_B2-19/cell_typing/myeloid/myeloid.batch.corrected_ct_subtype.rds")


mask <- myeloid.batch.corrected$subtype1 %in% c("TLB",
                           #"CD3E+CXCL8+ Mph","CD3E+GNLY+ Mph", "CD3E+HLA-DRA+ Mph", "CD79A+ Mph",
  "Remove",grep("Low",unique(myeloid.batch.corrected$subtype1),value = TRUE))
table(mask)
myeloid.batch.corrected <- myeloid.batch.corrected[,!mask]
dim(myeloid.batch.corrected)

meta_myeloid <- meta[colnames(myeloid.batch.corrected),]

all(colnames(myeloid.batch.corrected)==rownames(meta_myeloid))

myeloid.batch.corrected$LN_condition <- meta_myeloid$LN_condition
myeloid.batch.corrected$patient_condition <- meta_myeloid$patient_condition
myeloid.batch.corrected$Lauren <- meta_myeloid$Lauren

saveRDS(myeloid.batch.corrected,"plots/data/myeloid_new_meta.rds")

# bcell.batch.corrected <- readRDS("processed_data/data_B2-19/cell_typing/bcell/bcell.batch.corrected_ct_subtype.rds")

mask <- bcell.batch.corrected$subtype1 %in% c("TLB",
                           #"CD3E+CXCL8+ Mph","CD3E+GNLY+ Mph", "CD3E+HLA-DRA+ Mph", "CD79A+ Mph",
  "Remove",grep("Low",unique(bcell.batch.corrected$subtype1),value = TRUE))
table(mask)
bcell.batch.corrected <- bcell.batch.corrected[,!mask]
dim(bcell.batch.corrected)

meta_bcell <- meta[colnames(bcell.batch.corrected),]

all(colnames(bcell.batch.corrected)==rownames(meta_bcell))

bcell.batch.corrected$LN_condition <- meta_bcell$LN_condition
bcell.batch.corrected$patient_condition <- meta_bcell$patient_condition
bcell.batch.corrected$Lauren <- meta_bcell$Lauren

# saveRDS(bcell.batch.corrected,"plots/data/bcell_no_TLB_new_meta.rds")

epithelial.batch.corrected <- readRDS("processed_data/data_B2-19/cell_typing/epithelial.batch.corrected_cnv.rds")
# Idents(epithelial.batch.corrected) <- "LN_condition"

meta_epithelial <- meta[colnames(epithelial.batch.corrected),]

all(colnames(epithelial.batch.corrected)==rownames(meta_epithelial))

epithelial.batch.corrected$LN_condition <- meta_epithelial$LN_condition
epithelial.batch.corrected$patient_condition <- meta_epithelial$patient_condition
epithelial.batch.corrected$Lauren <- meta_epithelial$Lauren

saveRDS(epithelial.batch.corrected,"plots/data/epithelial_new_meta.rds")

ls()

GC <- readRDS("processed_data/data_B2-19/GC_qc500_mt0.2_ct_subtype.rds")

dim(GC)

mask <- GC$subtype1 %in% c("TLB",
                           #"CD3E+CXCL8+ Mph","CD3E+GNLY+ Mph", "CD3E+HLA-DRA+ Mph", "CD79A+ Mph",
  "Remove",grep("Low",unique(GC$subtype1),value = TRUE))
table(mask)
GC <- GC[,!mask]
dim(GC)

all(colnames(GC)==rownames(meta))

mask <- GC$subtype1 != meta$subtype1
GC$subtype1[mask] <- meta$subtype1[mask]


GC$LN_condition <- meta$LN_condition
GC$patient_condition <- meta$patient_condition
GC$Lauren <- meta$Lauren

GC@meta.data[,c('LN_condition.LN_station','patient_condition.region','patient_condition.LN_condition',
       'LN_condition.Lauren','condition')] <- NULL

table(GC$patient_condition)

table(GC$Lauren)

GC

saveRDS(GC,"plots/data/GC_no_TLB_new_meta.rds")

library(Matrix)

GC <- readRDS("plots/data/GC_no_TLB_new_meta.rds")

expr <- GC@assays$RNA@counts

dim(expr)

expr[1:10,1:10]

writeMM(expr,"plots/data/10X.GC.counts.mtx")

colnames(GC@meta.data)

meta_data <- GC@meta.data[,c("subtype1"), drop=FALSE]

dim(meta_data)
head(meta_data)

write.table(meta_data, "plots/data/10X.GC.metadata.txt")


