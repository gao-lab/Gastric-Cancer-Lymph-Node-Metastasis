library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)

packageVersion("CellChat")

options(stringsAsFactors = FALSE)

subsample_so <- function(seurat_obj){       
    # subset big subtype1 to 1000 cells
    set.seed=0
    cell_use <- NULL
    for(subtype1 in unique(seurat_obj$subtype1)){
        print(subtype1)
        print(sum(seurat_obj$subtype1==subtype1))
        if(sum(seurat_obj$subtype1==subtype1) > 1000){
            rowname <- sample(colnames(seurat_obj)[seurat_obj$subtype1==subtype1], 1000, replace = FALSE)
        }else{
            rowname <- colnames(seurat_obj)[seurat_obj$subtype1==subtype1]
        }
        cell_use <- c(cell_use,rowname)
        print(length(cell_use))
    }
    seurat_obj <- subset(seurat_obj,cells=cell_use)
}

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

run_CellChat <- function(seurat_obj){   
    cellchat.obj <- createCellChat(seurat_obj, group.by = "ident", assay = "RNA")
    #below is universe
    # CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    # # showDatabaseCategory(CellChatDB)
    # dplyr::glimpse(CellChatDB$interaction)

    # # use a subset of CellChatDB for cell-cell communication analysis
    # # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
    # # use all CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- CellChatDB # simply use the default CellChatDB

    # set the used database in the object

    cellchat.obj@DB <- CellChatDB.use

    options(future.globals.maxSize = 8000 * 1024^2)

    # subset the expression data of signaling genes for saving computation cost
    cellchat.obj <- subsetData(cellchat.obj) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 4) # do parallel
    #> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
    #> explicitly specify either 'multisession' or 'multicore'. In the current R
    #> session, 'multiprocess' equals 'multisession'.
    #> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
    #> processing ('multicore') is not supported when running R from RStudio
    #> because it is considered unstable. For more details, how to control forked
    #> processing or not, and how to silence this warning in future R sessions, see ?
    #> parallelly::supportsMulticore
    cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)


    # project gene expression data onto PPI (Optional: when running it, USER should 
    # set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
    # cellchat.obj <- projectData(cellchat.obj, PPI.human)

    cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)

    cellchat.obj <- computeCommunProb(cellchat.obj)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat.obj <- filterCommunication(cellchat.obj, min.cells = 10)

    cellchat.obj <- computeCommunProbPathway(cellchat.obj) #Infer the cell-cell communication at a signaling pathway level

    cellchat.obj <- aggregateNet(cellchat.obj) #Calculate the aggregated cell-cell communication network

    cellchat.obj <- netAnalysis_computeCentrality(cellchat.obj, slot.name = "netP") 
    return(cellchat.obj)    
} 

GC <- readRDS("plots/data/GC_no_TLB_new_meta.rds")

table(GC$LN_condition)

Idents(GC) <- "subtype1"

GC_all <- GC

table(GC$subtype1[GC$cell_type1=="Macro"])

mask <- GC$subtype1 %in% c("TLB","CD3E+CXCL8+ Mph","CD3E+GNLY+ Mph", "CD3E+HLA-DRA+ Mph", "CD79A+ Mph",
  "Remove",grep("Low",unique(GC$subtype1),value = TRUE))
table(mask)
GC <- GC[,!mask]
GC

length(unique(GC$subtype1))

grep("TLB",unique(GC$subtype1))

table(GC$patient_condition.LN_condition)

t1<-Sys.time()

GC_patient_meta_LN_neg <- subset(GC, patient_condition.LN_condition=="patient_meta_LN_neg")
GC_patient_meta_LN_neg <- subsample_so(GC_patient_meta_LN_neg)

Sys.time()

cellchat.patient_meta_LN_neg <- run_CellChat(GC_patient_meta_LN_neg)

Sys.time()

saveRDS(cellchat.patient_meta_LN_neg, "processed_data/data_B2-19/cellchat/cellchat.patient_meta_LN_neg.rds")

t2<-Sys.time()
t2-t1



table(GC$patient_condition.LN_condition)

t1<-Sys.time()

GC_patient_neg_LN_neg <- subset(GC, patient_condition.LN_condition=="patient_neg_LN_neg")
GC_patient_neg_LN_neg <- subsample_so(GC_patient_neg_LN_neg)

Sys.time()

cellchat.patient_neg_LN_neg <- run_CellChat(GC_patient_neg_LN_neg)

Sys.time()

saveRDS(cellchat.patient_neg_LN_neg, "processed_data/data_B2-19/cellchat/cellchat.patient_neg_LN_neg.rds")

t2<-Sys.time()
t2-t1

table(GC$patient_condition.LN_condition)

t1<-Sys.time()

GC_patient_meta_T <- subset(GC, patient_condition.LN_condition=="patient_meta_T")
GC_patient_meta_T <- subsample_so(GC_patient_meta_T)

Sys.time()

cellchat.patient_meta_T <- run_CellChat(GC_patient_meta_T)

Sys.time()

saveRDS(cellchat.patient_meta_T, "processed_data/data_B2-19/cellchat/cellchat.patient_meta_T.rds")

t2<-Sys.time()
t2-t1

table(GC$patient_condition.LN_condition)

t1<-Sys.time()

GC_patient_neg_T <- subset(GC, patient_condition.LN_condition=="patient_neg_T")
GC_patient_neg_T <- subsample_so(GC_patient_neg_T)

Sys.time()

cellchat.patient_neg_T <- run_CellChat(GC_patient_neg_T)

Sys.time()

saveRDS(cellchat.patient_neg_T, "processed_data/data_B2-19/cellchat/cellchat.patient_neg_T.rds")

t2<-Sys.time()
t2-t1

table(GC$LN_condition)

t1<-Sys.time()

GC_LN_neg <- subset(GC, LN_condition=="LN_neg")
GC_LN_neg <- subsample_so(GC_LN_neg)

Sys.time()

cellchat.LN_neg <- run_CellChat(GC_LN_neg)

Sys.time()

saveRDS(cellchat.LN_neg, "processed_data/data_B2-19/cellchat/cellchat.LN_neg.rds")

t2<-Sys.time()
t2-t1

table(GC$LN_condition.Lauren)

t1<-Sys.time()

GC_LN_meta_diffuse <- subset(GC, LN_condition.Lauren=="LN_meta_diffuse")
GC_LN_meta_diffuse <- subsample_so(GC_LN_meta_diffuse)

Sys.time()

cellchat.LN_meta_diffuse <- run_CellChat(GC_LN_meta_diffuse)

Sys.time()

saveRDS(cellchat.LN_meta_diffuse, "processed_data/data_B2-19/cellchat_v1/cellchat.LN_meta_diffuse.rds")

t2<-Sys.time()
t2-t1

t1<-Sys.time()

GC_LN_meta_intestinal <- subset(GC, LN_condition.Lauren=="LN_meta_intestinal")
GC_LN_meta_intestinal <- subsample_so(GC_LN_meta_intestinal)

Sys.time()

cellchat.LN_meta_intestinal <- run_CellChat(GC_LN_meta_intestinal)

Sys.time()

saveRDS(cellchat.LN_meta_intestinal, "processed_data/data_B2-19/cellchat_v1/cellchat.LN_meta_intestinal.rds")

t2<-Sys.time()
t2-t1

table(GC$LN_condition.Lauren)

t1<-Sys.time()

GC_T_diffuse <- subset(GC, LN_condition.Lauren=="T_diffuse")
GC_T_diffuse <- subsample_so(GC_T_diffuse)

Sys.time()

cellchat.T_diffuse <- run_CellChat(GC_T_diffuse)

Sys.time()

saveRDS(cellchat.T_diffuse, "processed_data/data_B2-19/cellchat/cellchat.T_diffuse.rds")

t2<-Sys.time()
t2-t1

t1<-Sys.time()

GC_T_intestinal <- subset(GC, LN_condition.Lauren=="T_intestinal")
GC_T_intestinal <- subsample_so(GC_T_intestinal)

Sys.time()

cellchat.T_intestinal <- run_CellChat(GC_T_intestinal)

Sys.time()

saveRDS(cellchat.T_intestinal, "processed_data/data_B2-19/cellchat/cellchat.T_intestinal.rds")

t2<-Sys.time()
t2-t1





GC_LN_meta <- subset(GC, LN_condition=="LN_meta")
# GC.condition.list <- SplitObject(GC, split.by = "LN_condition")

GC_LN_meta <- subsample_so(GC_LN_meta)

cellchat.LN_meta <- run_CellChat(GC_LN_meta)

saveRDS(cellchat.LN_meta, "processed_data/data_B2-19/cellchat/cellchat.LN_meta.rds")



GC_PT <- subset(GC, LN_condition=="T")
GC_PT <- subsample_so(GC_PT)

cellchat.PT <- run_CellChat(GC_PT)

saveRDS(cellchat.PT, "processed_data/data_B2-19/cellchat/cellchat.PT.rds")



GC_all <- subsample_so(GC)

cellchat.all <- run_CellChat(GC_all)

saveRDS(cellchat.all, "processed_data/data_B2-19/cellchat/cellchat.all.rds")
