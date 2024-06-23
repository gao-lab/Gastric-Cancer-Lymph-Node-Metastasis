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

GC <- readRDS("processed_data/data_B2-19/GC_qc500_mt0.2_ct_subtype.rds")


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

ls()

which(rownames(cellchat.LN_meta_diffuse@net$count)=="SPIB+ DC")

cd8_tcell <- which(rownames(cellchat.LN_meta_diffuse@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type1=='T-CD4']))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.LN_meta_diffuse, sources.use = 73,targets.use = cd8_tcell,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


cd8_tcell <- which(rownames(cellchat.LN_meta_diffuse@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type1=='T-CD8']))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.LN_meta_diffuse, sources.use = 73,targets.use = cd8_tcell,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


cd8_tcell <- which(rownames(cellchat.LN_meta_diffuse@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type1=='T-CD8']))

options(repr.plot.width=11, repr.plot.height=15)
netVisual_bubble(cellchat.LN_meta_diffuse, sources.use = 73,#targets.use = cd8_tcell,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


which(rownames(cellchat.LN_meta_diffuse@net$count)=="LAMP3+ DC")

options(repr.plot.width=11, repr.plot.height=15)
netVisual_bubble(cellchat.LN_meta_diffuse, sources.use = 35,#targets.use = cd8_tcell,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


which(rownames(cellchat.LN_meta_diffuse@net$count)=="LILRA4+ pDC")

options(repr.plot.width=11, repr.plot.height=15)
netVisual_bubble(cellchat.LN_meta_diffuse, sources.use = 43,#targets.use = cd8_tcell,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),




# cellchat.LN_neg <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_neg.rds")
# cellchat.LN_meta <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_meta.rds")
# cellchat.PT <- readRDS("processed_data/data_B2-19/cellchat/cellchat.PT.rds")

cellchat.patient_meta_T <- readRDS("processed_data/data_B2-19/cellchat/cellchat.patient_meta_T.rds")
cellchat.patient_neg_T <- readRDS("processed_data/data_B2-19/cellchat/cellchat.patient_neg_T.rds")
cellchat.patient_meta_LN_neg <- readRDS("processed_data/data_B2-19/cellchat/cellchat.patient_meta_LN_neg.rds")
cellchat.patient_neg_LN_neg <- readRDS("processed_data/data_B2-19/cellchat/cellchat.patient_neg_LN_neg.rds")

# cellchat.T_diffuse <- readRDS("processed_data/data_B2-19/cellchat/cellchat.T_diffuse.rds")
# cellchat.T_intestinal <- readRDS("processed_data/data_B2-19/cellchat/cellchat.T_intestinal.rds")
# cellchat.LN_meta_diffuse <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_meta_diffuse.rds")
# cellchat.LN_meta_intestinal <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_meta_intestinal.rds")

# par(mfrow=c(1,1))
options(repr.plot.width=12, repr.plot.height=12)
netVisual_heatmap(cellchat.patient_meta_LN_neg, color.heatmap = "Reds")

meta_condition <- read.csv("processed_data/data_B2-19/metadata_all_ct_subtype_condition.csv", row.names = 1)


unique(meta_condition$cell_type1)

table(meta_condition$cell_type1,meta_condition$cell_type)

unique(meta_condition$subtype1[meta_condition$cell_type1=='T-CD8'])

unique(meta_condition$subtype1[meta_condition$cell_type1=='T-CD4'])





Myeloid <- which(rownames(cellchat.patient_meta_LN_neg@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type %in% c("Endothelial","Epithelial", "Fibroblast")]))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51,targets.use = Myeloid,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


Myeloid <- which(rownames(cellchat.patient_meta_LN_neg@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type1 %in% c("NK","NKT")]))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51,targets.use = Myeloid,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


Myeloid <- which(rownames(cellchat.patient_meta_LN_neg@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type==c('Myeloid',"Mast")]))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51,targets.use = Myeloid,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


Myeloid <- which(rownames(cellchat.patient_neg_LN_neg@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type=='Myeloid']))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_neg_LN_neg, sources.use = 51,targets.use = Myeloid,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


cd4_tcell <- which(rownames(cellchat.patient_meta_LN_neg@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type1=='T-CD4']))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51,targets.use = cd4_tcell,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


cd8_tcell <- which(rownames(cellchat.patient_meta_LN_neg@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type1=='T-CD8']))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51,targets.use = cd8_tcell,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


cd8_tcell <- which(rownames(cellchat.patient_meta_LN_neg@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type=='B cell']))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51,targets.use = cd8_tcell,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


which(rownames(cellchat.patient_meta_LN_neg@net$count)=="Endothelial")

which(rownames(cellchat.patient_meta_LN_neg@net$count)=="FCN1+ Mph")

options(repr.plot.width=11, repr.plot.height=16)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51, remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51,targets.use = 46,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),



options(repr.plot.width=12, repr.plot.height=9)
par(mfrow=c(1,2))
netVisual_bubble(cellchat.patient_meta_LN_neg, sources.use = 51,targets.use = 46,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),
netVisual_bubble(cellchat.patient_neg_LN_neg, sources.use = 51,targets.use = 46,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_neg_LN_neg, sources.use = 51,targets.use = 46,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),




cellchat.patient_neg_LN_neg <- readRDS("processed_data/data_B2-19/cellchat/cellchat.patient_neg_LN_neg.rds")
# cellchat.LN_meta_intestinal <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_meta_intestinal.rds")
# cellchat.T_patient_meta <- readRDS("processed_data/data_B2-19/cellchat/cellchat.T_patient_meta.rds")
# cellchat.T_patient_neg <- readRDS("processed_data/data_B2-19/cellchat/cellchat.T_patient_neg.rds")
# cellchat.LN_neg <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_neg.rds")
# cellchat.LN_meta <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_meta.rds")
# cellchat.PT <- readRDS("processed_data/data_B2-19/cellchat/cellchat.PT.rds")

which(rownames(cellchat.patient_neg_LN_neg@net$count)=="Endothelial")

which(rownames(cellchat.patient_neg_LN_neg@net$count)=="CD8_GZMB+ Temra")

which(rownames(cellchat.patient_neg_LN_neg@net$count)=="FCN1+ Mph")

Myeloid <- which(rownames(cellchat.patient_neg_LN_neg@net$count) %in% 
      unique(meta_condition$subtype1[meta_condition$cell_type %in% c("Endothelial","Epithelial", "Fibroblast")]))

options(repr.plot.width=9, repr.plot.height=9)
netVisual_bubble(cellchat.patient_neg_LN_neg, sources.use = 51,targets.use = Myeloid,
                 remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),


options(repr.plot.width=11, repr.plot.height=16)
netVisual_bubble(cellchat.patient_neg_LN_neg, sources.use = 51, remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),




groupSize <- as.numeric(table(cellchat.patient_meta_LN_neg@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.patient_meta_LN_neg@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.patient_meta_LN_neg@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

head(mat)

rownames(mat)

mat <- cellchat.patient_meta_LN_neg@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dim(mat2[,colSums(mat2)>0])
dim(mat2)

options(repr.plot.width=11, repr.plot.height=11)
mat <- cellchat.patient_meta_LN_neg@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
for (i in c(51)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

options(repr.plot.width=11, repr.plot.height=11)
mat <- cellchat.patient_meta_LN_neg@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
for (i in c(18)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[, i] <- mat[, i]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}# max(mat) can change to other max value

options(repr.plot.width=11, repr.plot.height=11)
mat <- cellchat.patient_meta_LN_neg@net$count
# par(mfrow = c(3,4), xpd=TRUE)
for (i in c(18)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[, i] <- mat[, i]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

par(mfrow=c(1,1))
netVisual_heatmap(cellchat.patient_neg_LN_neg, color.heatmap = "Reds")

ls()

which(rownames(mat)=="FCN1+ Mph")

options(repr.plot.width=11, repr.plot.height=16)
netVisual_bubble(cellchat.patient_neg_LN_neg, sources.use = 51,  remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),

options(repr.plot.width=11, repr.plot.height=11)
mat <- cellchat.LN_neg@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
for (i in c(18)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

options(repr.plot.width=11, repr.plot.height=11)
mat <- cellchat.LN_neg@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
for (i in c(18)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[, i] <- mat[, i]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}# max(mat) can change to other max value

options(repr.plot.width=11, repr.plot.height=11)
mat <- cellchat.LN_neg@net$count
# par(mfrow = c(3,4), xpd=TRUE)
for (i in c(18)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[, i] <- mat[, i]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

par(mfrow=c(1,1))
netVisual_heatmap(cellchat.LN_neg, color.heatmap = "Reds")

which(rownames(mat)=="Epithelial")

options(repr.plot.width=11, repr.plot.height=16)
netVisual_bubble(cellchat.LN_neg, sources.use = 18,  remove.isolate = FALSE,thresh=0.01)#targets.use = c(5:11),

netVisual_bubble(cellchat.LN_neg, sources.use = 18, targets.use = 37, remove.isolate = FALSE)#

ls()

cellchat.LN_meta_diffuse <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_meta_diffuse.rds")
cellchat.LN_meta_intestinal <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_meta_intestinal.rds")
cellchat.T_patient_meta <- readRDS("processed_data/data_B2-19/cellchat/cellchat.T_patient_meta.rds")
cellchat.T_patient_neg <- readRDS("processed_data/data_B2-19/cellchat/cellchat.T_patient_neg.rds")
cellchat.LN_neg <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_neg.rds")
cellchat.LN_meta <- readRDS("processed_data/data_B2-19/cellchat/cellchat.LN_meta.rds")
cellchat.PT <- readRDS("processed_data/data_B2-19/cellchat/cellchat.PT.rds")

object.list <- list(LN_meta_diffuse = cellchat.LN_meta_diffuse, 
                    LN_meta_intestinal = cellchat.LN_meta_intestinal,
                    T_diffuse = cellchat.T_diffuse,
                    T_intestinal = cellchat.T_intestinal,
                    T_patient_meta = cellchat.T_patient_meta,
                    T_patient_neg = cellchat.T_patient_neg,
                    PT = cellchat.PT,
                    LN_meta = cellchat.LN_meta,
                    LN_neg = cellchat.LN_neg                    
                   )

names(object.list)

# condition <- "PT"
for(condition in names(object.list)){
    pdf(paste0("processed_data/data_B2-19/cellchat/LAMP3_DC/",condition,".pdf"), width = 12,height=16)
    plot(netVisual_heatmap(object.list[[condition]], color.heatmap = "Reds",
                      title.name = paste0("Number of interactions in ", condition),
                     font.size.title = 20))
    plot(netVisual_bubble(object.list[[condition]], sources.use = 18, #targets.use = myeloid_group, 
                     remove.isolate = FALSE,thresh=0.01,
                    title.name = paste0("Signaling from LAMP3+ DC in ", condition), 
                     font.size.title = 20))

    plot(netVisual_bubble(object.list[[condition]], targets.use = 18, #targets.use = myeloid_group, 
                     remove.isolate = FALSE,thresh=0.01,
                    title.name = paste0("Signaling to LAMP3+ DC in ", condition), 
                     font.size.title = 20))
    dev.off()
}

# condition <- "PT"
pdf(paste0("processed_data/data_B2-19/cellchat/LAMP3_DC/heatmap_per_condition.pdf"), width = 10,height=10)
for(condition in names(object.list)){
    p1<-plot(netVisual_heatmap(object.list[[condition]], color.heatmap = "Reds",
                      title.name = paste0("Number of interactions in ", condition),
                     font.size.title = 20))
    p1
#     plot(netVisual_bubble(object.list[[condition]], sources.use = 18, #targets.use = myeloid_group, 
#                      remove.isolate = FALSE,thresh=0.01,
#                     title.name = paste0("Signaling from LAMP3+ DC in ", condition), 
#                      font.size.title = 20))

#     plot(netVisual_bubble(object.list[[condition]], targets.use = 18, #targets.use = myeloid_group, 
#                      remove.isolate = FALSE,thresh=0.01,
#                     title.name = paste0("Signaling to LAMP3+ DC in ", condition), 
#                      font.size.title = 20))   
}
dev.off()

# only show the last heatmap, not ok!

condition <- "PT"
pdf("processed_data/data_B2-19/cellchat/LAMP3_DC/condition.pdf", width = 12,height=16)
# groupSize <- as.numeric(table(object.list[[condition]]@idents))
# options(repr.plot.width=12, repr.plot.height=12)

# mat <- cellchat.T_patient_meta@net$count
# i=18
# par(mfrow = c(2,1), xpd=TRUE)
# for (i in c(18)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
#                    edge.weight.max = max(mat[i,],mat[,i]), 
#                    title.name = paste0("outgoing signaling for LAMP3+ DC in ", condition))
# }
# for (i in c(18)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[, i] <- mat[, i]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
#                    edge.weight.max = max(mat[i,],mat[,i]), 
#                    title.name = paste0("incoming signaling for LAMP3+ DC in ", condition))
# par(mfrow=c(1,1))
netVisual_heatmap(object.list[[condition]], color.heatmap = "Reds",
                  title.name = paste0("Number of interactions in ", condition),
                 font.size.title = 20)

# options(repr.plot.width=12, repr.plot.height=16)
netVisual_bubble(object.list[[condition]], sources.use = 18, #targets.use = myeloid_group, 
                 remove.isolate = FALSE,thresh=0.01,
                title.name = "Signaling for LAMP3+ DC and T cell in T_intestinal", 
                 font.size.title = 20)

}
dev.off()

par(mfrow=c(1,1))
netVisual_heatmap(object.list[[condition]], color.heatmap = "Reds",
                  title.name = paste0("Number of interactions in ", condition),
                 font.size.title = 20)

options(repr.plot.width=11, repr.plot.height=16)
netVisual_bubble(object.list[[condition]], sources.use = 18, targets.use = myeloid_group, 
                 remove.isolate = FALSE,thresh=0.01,
                title.name = "Signaling for LAMP3+ DC and T cell in T_intestinal", 
                 font.size.title = 20)

options(repr.plot.width=11, repr.plot.height=16)
netVisual_bubble(object.list[[condition]], sources.use = 18, targets.use = t_group, 
                 remove.isolate = FALSE,thresh=0.01,
                title.name = "Signaling for LAMP3+ DC and T cell in T_intestinal", 
                 font.size.title = 20)

options(repr.plot.width=12, repr.plot.height=16)
netVisual_bubble(object.list[[condition]], sources.use = 18, #targets.use = myeloid_group, 
                 remove.isolate = FALSE,thresh=0.01,
                title.name = "Signaling for LAMP3+ DC and T cell in T_intestinal", 
                 font.size.title = 20)


