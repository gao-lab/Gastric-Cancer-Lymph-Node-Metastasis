library(Seurat)
library(infercnv)

GC <- readRDS("processed_data/data_B2-19/GC_qc500_mt0.2_ct_subtype.rds")
# prepare matrix file
raw_matrix <- GC@assays$RNA@counts

annotation <- GC@meta.data


annotation$cell_type1 <- as.character(annotation$cell_type1)

GC

raw_matrix0 <- raw_matrix
annotation0 <- annotation

mask_ct=annotation0$cell_type1 %in% c("Mast", "Endothelial", "Epithelial")
raw_matrix <- raw_matrix0[, mask_ct]
annotation <- annotation0[mask_ct,]
annotation$cell_type1 <- as.character(annotation$cell_type1)

dim(raw_matrix)
dim(annotation)

table(annotation$cell_type1, annotation$sample)

table(annotation$cell_type1)

mask_epi=annotation$cell_type1 == "Epithelial"
annotation$cell_type1[mask_epi] <- annotation$sample[mask_epi]
annotation <- annotation[,"cell_type1", drop=FALSE]
colnames(annotation) <- c("cell_type")

head(annotation)

tmp <- data.frame(table(annotation$cell_type))
sample_use <- as.character(tmp[tmp$Freq>5,1])
sample_use

mask_sample <- annotation$cell_type %in% sample_use
raw_matrix <- raw_matrix[, mask_sample]
annotation <- annotation[mask_sample, ,drop=FALSE]

dim(annotation)

# all Epi
write.table(annotation, "processed_data/data_B2-19/inferCNV_new/anno_file.txt", sep="\t", col.names = FALSE)

annotation<- read.table("processed_data/data_B2-19/inferCNV_new/anno_file.txt", sep="\t")

table(annotation$V2)

# create infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_matrix,
                                    annotations_file="processed_data/data_B2-19/inferCNV_new/anno_file.txt",
                                    delim="\t",
                                    gene_order_file="../GENOME/inferCNV/gencode_v21_gen_pos.complete_symbol.txt",
                                    ref_group_names=c("Endothelial", "Mast"))


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="processed_data/data_B2-19/inferCNV_new/results/", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)




