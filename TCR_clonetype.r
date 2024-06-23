library(Seurat)
library(reshape2)
library(pheatmap)
library(scRepertoire)
# suppressMessages(library(scRepertoire))
library(reshape2)
library(pheatmap)

samples=dir("02_cr_results_raw/scRNA/")
length(samples)

contig_t_list <- list()
for(sample in samples){ 
    print(sample)
    contig_t_list[[sample]] <- read.csv(paste0("02_cr_results_raw/scRNA/", sample, "/outs/per_sample_outs/",
                              sample, "/vdj_t/filtered_contig_annotations.csv"))
}

names(contig_t_list)[grep("16",names(contig_t_list))] <- paste0("B",names(contig_t_list)[grep("16",names(contig_t_list))])
names(contig_t_list)[grep("B2",names(contig_t_list))] <- gsub("-MULTI","",names(contig_t_list)[grep("B2",names(contig_t_list))])
names(contig_t_list)[grep("G",names(contig_t_list))] <- gsub("G","B",names(contig_t_list)[grep("G",names(contig_t_list))])
names(contig_t_list)[grep("B13L3",names(contig_t_list))] <- "B13L3"

names(contig_t_list)

for(sample in names(contig_t_list)){
    cat(sample, unique(contig_t_list[[sample]]$chain), "\n")
}


combined <- combineTCR(contig_t_list, 
                samples = names(contig_t_list), 
                cells ="T-AB")

saveRDS(combined, "processed_data/data_B2-19/TCR/combined.rds")

combined <- readRDS("processed_data/data_B2-19/TCR/combined.rds")

tcell.batch.corrected <- readRDS("plots/data/tcell_no_TLB_new_meta.rds")

seurat_t_all <- combineExpression(combined, tcell.batch.corrected, 
            cloneCall="gene", group.by = "sample", proportion = FALSE, 
            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

saveRDS(seurat_t_all, "processed_data/data_B2-19/TCR/tcell_all.rds")

seurat_t_all

seurat_t_tcr <- seurat_t_all[, !is.na(seurat_t_all$barcode)]
# seurat_t_tcr$subtype1.sample <- paste(seurat_t_tcr$subtype1, seurat_t_tcr$sample, sep = "_")

seurat_t_tcr

saveRDS(seurat_t_tcr, "processed_data/data_B2-19/TCR/tcell_with_tcr.rds")

dim(seurat_t_tcr)

table(seurat_t_tcr$sample)

combined <- readRDS("processed_data/data_B2-19/TCR/combined.rds")

combined2 <- combined
for(sample_name in names(combined2)){
#     print(sample_name)
    pattern <- "B\\d*"
    m <- regexpr(pattern, sample_name)
    combined2[[sample_name]]$patient <- regmatches(sample_name, m)
}

hla_tem <- subset(tcell.batch.corrected, subtype1=='CD8_HLA-DRA+ Tem')

# clono types expansion computed by only CD8_HLA-DRA+ Tem in patient level

seurat_t_HLA <- combineExpression(combined2, hla_tem, 
            cloneCall="gene", group.by = "patient", proportion = FALSE, #根据patient计算每个clone的frequency，
            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))#cloneTypes是基于上面计算的frequency的

seurat_t_tcr_HLA <- seurat_t_HLA[, !is.na(seurat_t_HLA$barcode)]


seurat_t_tcr_HLA

saveRDS(seurat_t_tcr_HLA,"processed_data/data_B2-19/TCR/seurat_t_tcr_HLA.rds")
