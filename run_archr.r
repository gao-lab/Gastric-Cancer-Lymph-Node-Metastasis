library(ArchR)
library(ggplot2)
library(hexbin)
library(dplyr)

library(BSgenome.Hsapiens.UCSC.hg38)

packageVersion("ggplot2")

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

addArchRThreads(threads = 16) 

addArchRGenome("hg38")

out_dir <- "archr_results"

inputFiles<-c("../02_cr_results/scATAC/B2L1-ATAC/outs/fragments.tsv.gz", 
             "../02_cr_results/scATAC/B2L2-ATAC/outs/fragments.tsv.gz",
    "../02_cr_results/scATAC/B2T-ATAC/outs/fragments.tsv.gz", 
              "../02_cr_results/scATAC/B8T-ATAC/outs/fragments.tsv.gz", 
              "../02_cr_results/scATAC/B8L6-ATAC/outs/fragments.tsv.gz", 
             "../02_cr_results/scATAC/B8L11-ATAC/outs/fragments.tsv.gz",
              "../02_cr_results/scATAC/B17T-ATAC/outs/fragments.tsv.gz"
             )
names(inputFiles)<-c("B2L1", "B2L2", "B2T", "B8T", "B8L6", "B8L11","B17T")

inputFiles

addGeneScoreMatrix

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

###之后可通过如下方法进一步筛选
# idxPass <- which(projHeme1$TSSEnrichment >= 8)
# cellsPass <- projHeme1$cellNames[idxPass]
# projHeme1[cellsPass, ]

ArrowFiles

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

out_dir

proj <- ArchRProject(
  ArrowFiles = c('B8L6.arrow', 'B2L1.arrow','B8L11.arrow', 'B2L2.arrow','B2T.arrow','B8T.arrow','B17T.arrow'), 
  outputDirectory = out_dir,
  copyArrows = TRUE,  #This is recommened so that you maintain an unaltered copy for later usage.
)

proj

proj <- filterDoublets(proj)

cellMeta <- getCellColData(proj)
mask=(cellMeta$TSSEnrichment>=8)

table(mask, cellMeta$Sample)

# QC
proj<-proj[rownames(cellMeta)[mask], ]

pattern <- "B\\d*"
m <- regexpr(pattern, proj$Sample)
proj$patient <- regmatches(proj$Sample, m)
unique(proj$patient)

pattern <- "[[:alpha:]]+"
m <- gregexpr(pattern, proj$Sample)
proj$region <- as.character(lapply(regmatches(proj$Sample, m), function(x) x[[2]]))
unique(proj$region)

unique(proj$Sample)

# gene_score <- getMatrixFromProject(proj, "GeneScoreMatrix")

# tile_matrix <- getMatrixFromProject(proj, "TileMatrix",binarize = TRUE)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI",force=TRUE)

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force=TRUE)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

options(repr.plot.width=12, repr.plot.height=8)
ggAlignPlots(p1, p2, type = "h")

options(repr.plot.width=8, repr.plot.height=8)
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "patient", embedding = "UMAP")

options(repr.plot.width=8, repr.plot.height=8)
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

options(repr.plot.width=8, repr.plot.height=8)
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type1_glue", embedding = "UMAP")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)



saveArchRProject(ArchRProj = proj, load = TRUE)

# batch effect removal if needed
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "patient"
)

proj <- addClusters(input = proj, reducedDims = "Harmony", force=TRUE)

proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony", force=TRUE)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

options(repr.plot.width=12, repr.plot.height=8)
ggAlignPlots(p1, p2, type = "h")

# plotPDF(p1,p2, name = "Plot-Harmony-UMAP-Sample-Clusters.pdf",
#         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

options(repr.plot.width=8, repr.plot.height=8)
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

options(repr.plot.width=8, repr.plot.height=8)
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type1_glue", embedding = "UMAP")

options(repr.plot.width=8, repr.plot.height=8)
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "patient", embedding = "UMAP")

options(repr.plot.width=8, repr.plot.height=8)
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "region", embedding = "UMAP")

saveArchRProject(ArchRProj = proj, load = TRUE)

proj <- addImputeWeights(proj)

proj<-loadArchRProject(path = "archr_results")

packageVersion("hexbin")
# plotAs = "points"

markerGenes  <- c("CD79A", "CHGA", "ENG", "KRT18", "COL1A2", "CD14", "CPA3", "ATP4B", "CD3D")

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

p$CD3D

# #Rearrange for grid plotting
# p2 <- lapply(p, function(x){
#     x + guides(color = FALSE, fill = FALSE) + 
#     theme_ArchR(baseSize = 6.5) +
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#     theme(
#         axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(), 
#         axis.text.y=element_blank(), 
#         axis.ticks.y=element_blank()
#     )
# })
# do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

user_markers<-c("EPCAM", "KRT18", "PGC", "TFF1", "MUC5AC", "GIF", "CHGA",#epithelial
                "CD68", "CSF1R","CD14", #macro
                "PECAM1", "VWF", #endo
                "ACTA2",  "COL1A1", "LUM", "DCN", #fibro
                "TPSAB1","TPSB2",#mast cell
                "CXCL8", #neutrophil
                "CD79A", # Bcell
                "CD2", "CD3E", "CD3D", "CD4", "CD8A", #T cell
                
                "CD27", #memory B cell
                "TCL1A", "MME", #naive B
               "CD69",# activate B
               "SDC1","CD38",#"IGHG1", "IGKC",#plasma
                
                "LTB", "SPOCK2","SOCS3","PBXIP1",#CD4 Tconv
               "FOXP3","TNFRSF4","CARD16","IL2RA",#"TBCD4",#Treg
               "CST7","GZMA","GZMH",#"CCL4",#CD8 T cell
               "CTLA4","LAG3","PTMS","PDCD1","TIGIT",#Exhaustion
                "CD44","GZMK","IL7R","CD69","CD27"#Memory T
               )
p_more <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = user_markers, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p_more, 
    name = "Plot-UMAP-Detailed-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

seRNA_all <- readRDS("../plots/data/GC_no_TLB_new_meta.rds")#only paired samples

table(seRNA_all$sample)

table(proj$Sample)

table(seRNA_all$sample[seRNA_all$sample %in% c(unique(proj$Sample),"B14T","B16T")])

seRNA <- subset(seRNA_all, sample %in% c(unique(proj$Sample),"B14T","B16T"))
#two intestinal PT with metastasis

dim(seRNA)

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "cell_type",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

pal <- paletteDiscrete(values = seRNA$cell_type)

p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal
)

plotPDF(p1, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, load = TRUE)

# proj<-loadArchRProject(path = "Save-proj2")

getAvailableMatrices(proj)

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    groupRNA = "cell_type",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force=TRUE
)

getAvailableMatrices(proj)

table(proj$Clusters)

colnames(proj@cellColData)

saveArchRProject(ArchRProj = proj, load = TRUE)

proj$cell_type="T cell"
proj$cell_type[proj$Clusters %in% paste0("C", c(3))]="Endothelial"
proj$cell_type[proj$Clusters %in% paste0("C", c(18:23))]="B cell"
proj$cell_type[proj$Clusters %in% paste0("C", c(2))]="Myeloid"
proj$cell_type[proj$Clusters == paste0("C", 1)]="Epithelial"
proj$cell_type[proj$Clusters == paste0("C", 15)]="Plasma"

p3<-plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

plotPDF(p3, name = "Plot-UMAP-cell_type-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p3

table(proj$cell_type, proj$Sample)

cluster_group <- table(proj$cell_type, proj$Sample)
cluster_group <- melt(cluster_group, id.vars = "sample")
colnames(cluster_group) <- c("Cell_type", "Sample", "Freq")
cluster_group$Cell_type <- as.factor(cluster_group$Cell_type)

# pdf("cell_typing/ATAC_cell_type_composition.pdf")
ggplot(cluster_group, aes(x = Cell_type, y = Freq, fill = Sample)) +
  geom_bar(position = "fill", stat = "identity")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(panel.background = element_blank())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line())+
  theme(axis.text.y = element_text(color = "black"))+
  theme(axis.text.x = element_text(color = "black"))
ggplot(cluster_group, aes(x = Sample, y = Freq, fill = Cell_type)) +
  geom_bar(position = "fill", stat = "identity")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(panel.background = element_blank())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line = element_line())+
  theme(axis.text.y = element_text(color = "black"))+
  theme(axis.text.x = element_text(color = "black"))
# dev.off()

# proj$Project <- "all"

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "cell_type")
# assing cell type based on gene score and integration with scRNA, 
# to call peaks with the pseudobulk divided by cell type

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "cell_type", 
    pathToMacs2 = "/home/weil/software/miniconda3/envs/macs2/bin/macs2"
)

getPeakSet(proj)

saveArchRProject(ArchRProj = proj, load = TRUE)

proj <- addPeakMatrix(proj)

getAvailableMatrices(proj)

peak_matrix <- getMatrixFromProject(proj, "PeakMatrix")

proj

peak_matrix

dim(peak_matrix)

assay(peak_matrix)[1:3,1:3]

peaks <- rowRanges(peak_matrix)
rownames(peak_matrix) <- paste(as.character(seqnames(peaks)), start(peaks), end(peaks), sep = "_")

peaks[12345]

cellMeta <- getCellColData(proj)
all(colnames(peak_matrix)==rownames(cellMeta))

table(cellMeta$Sample)

cellMeta=cellMeta[colnames(peak_matrix),]
all(colnames(peak_matrix)==rownames(cellMeta))

# write out for GLUE
writeMM(assay(peak_matrix), file = paste0(out_dir, "/peak_matrix/peak_matrix.mtx"))
write.table(rownames(peak_matrix), file = paste0(out_dir, "/peak_matrix/peakinfo.txt"), sep = "\t", 
            row.names = F, col.names = F, quote = F)

# cellMeta <- getCellColData(proj)
write.csv(cellMeta, paste0(out_dir, "/peak_matrix/cellname.txt"))



table(seRNA$Fresh)

table(seRNA$sample)

# write out for GLUE
library(Matrix)
writeMM(seRNA@assays$RNA@counts, file = "seRNA/expr_counts.mtx")
write.table(data.frame(row.names = rownames(seRNA@assays$RNA@counts)), file = "seRNA/genes.csv", sep = ",", 
            row.names = T)

write.table(seRNA@meta.data, file = "seRNA/metadata.csv", sep = ",", 
            row.names = T)

dim(seRNA@assays$RNA@counts)

dim(data.frame(row.names = rownames(seRNA@assays$RNA@counts)))



table(proj_cd8_t$subtype1_glue,proj_cd8_t$cell_type1_glue)

proj_cd8_t_ct <- proj[proj$cell_type1_glue == "T-CD8",]

peak_matrix_cd8 <- getMatrixFromProject(proj_cd8_t_ct, "PeakMatrix")

peaks_cd8 <- rowRanges(peak_matrix_cd8)
rownames(peak_matrix_cd8) <- paste(as.character(seqnames(peaks_cd8)), start(peaks_cd8), end(peaks_cd8), sep = "_")

peak_matrix_cd8[1:3,1:3]
dim(peak_matrix_cd8)

cellMeta <- getCellColData(proj_cd8_t_ct)
all(colnames(peak_matrix_cd8)==rownames(cellMeta))

dim(cellMeta)

table(cellMeta$Sample)

cellMeta=cellMeta[colnames(peak_matrix_cd8),]
all(colnames(peak_matrix_cd8)==rownames(cellMeta))

# write out for GLUE
writeMM(assay(peak_matrix_cd8), file = paste0("GLUE_cd8", "/peak_matrix_cd8/peak_matrix_cd8.mtx"))
write.table(rownames(peak_matrix_cd8), file = paste0("GLUE_cd8", "/peak_matrix_cd8/peakinfo.txt"), sep = "\t", 
            row.names = F, col.names = F, quote = F)

# cellMeta <- getCellColData(proj_cd8_t_ct)
write.csv(cellMeta, paste0("GLUE_cd8", "/peak_matrix_cd8/cellname.txt"))



seRNA_cd8 <- subset(seRNA, cell_type1=="T-CD8")
#two intestinal PT with metastasis

table(seRNA_cd8$sample)

# write out for GLUE
library(Matrix)
writeMM(seRNA_cd8@assays$RNA@counts, file = "seRNA_cd8/expr_counts.mtx")
write.table(data.frame(row.names = rownames(seRNA_cd8@assays$RNA@counts)), file = "seRNA_cd8/genes.csv", sep = ",", 
            row.names = T)

write.table(seRNA_cd8@meta.data, file = "seRNA_cd8/metadata.csv", sep = ",", 
            row.names = T)

dim(seRNA_cd8@assays$RNA@counts)

proj<-loadArchRProject(path = "archr_results")

getAvailableMatrices(proj)

colnames(proj@cellColData)

cellMeta <- getCellColData(proj)

all(rownames(cellMeta)==rownames(proj@cellColData))

cell_subtype=read.csv("GLUE/atac_ct_glue.csv", row.names = 1)

table(cell_subtype$subtype1)

colnames(cell_subtype)

all(rownames(proj@cellColData)==rownames(cell_subtype))

cell_subtype=cell_subtype[rownames(proj@cellColData),]
all(rownames(proj@cellColData)==rownames(cell_subtype))

proj$cell_type1_glue=cell_subtype$cell_type1
proj$cell_type_glue=cell_subtype$cell_type
proj$subtype1_glue=cell_subtype$subtype1

table(proj$subtype1_glue)



p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type1_glue", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "subtype1_glue", embedding = "UMAP")

p1

p2

options(repr.plot.width=12, repr.plot.height=8)
ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-GLUE-celltypes.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
# need to downgrade specifically ggplot fixed this. This worked for me:
# conda install conda-forge::r-ggplot2=3.3



colnames(proj@cellColData)

proj$LN_condition <- proj$region 
proj$LN_condition[proj$region=="L"] <- "Met.LN"
proj$LN_condition[proj$region=="T"] <- "Pri.GC"
proj$LN_condition[proj$Sample=="B8L11"] <- "Neg.LN"
table(proj$LN_condition,proj$Sample)

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type_glue", embedding = "UMAP")

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

# add motif presence
proj <- addMotifAnnotations(ArchRProj = proj,
                            motifSet = "cisbp",
                            name = "Motif",
                            species = "homo sapiens",
                            force = TRUE)

if("Motif" %ni% names(proj@peakAnnotation)){
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}

proj <- addBgdPeaks(proj)

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

plotVarDev

proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
)

proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
)

saveArchRProject(ArchRProj = proj, load = TRUE)



proj_cd4_tex <- proj[proj$subtype1_glue=="CD4_Tex", ]


proj_cd4_tex

p_CTLA4 <- plotBrowserTrack(
    ArchRProj = proj_cd4_tex, 
#     region=tmp,
    groupBy = "LN_condition", 
#     useGroups = c("CD4_DUSP4+ Treg","CD4_BATF+ Treg","CD4_MKI67+ Treg","CD4_Tex","CD4_MKI67+ Tex"),
    plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
    geneSymbol = "CTLA4",
#     features =  getMarkers(markerTest, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR = TRUE),
    loops = getPeak2GeneLinks(proj),
    
    upstream = 50000,
    downstream = 50000
)

plotPDF(plotList = p_CTLA4, 
    name = "Plot-Tracks-CD4_Tex-CTLA4.pdf", 
    ArchRProj = proj_cd4_tex, 
    addDOC = FALSE, width = 5, height = 5)

p_PDCD1 <- plotBrowserTrack(
    ArchRProj = proj_cd4_tex, 
#     region=tmp,
    groupBy = "LN_condition", 
#     useGroups = c("CD4_DUSP4+ Treg","CD4_BATF+ Treg","CD4_MKI67+ Treg","CD4_Tex","CD4_MKI67+ Tex"),
    plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
    geneSymbol = "PDCD1",
#     features =  getMarkers(markerTest, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR = TRUE),
    loops = getPeak2GeneLinks(proj),
    
    upstream = 50000,
    downstream = 50000
)

plotPDF(plotList = p_PDCD1, 
    name = "Plot-Tracks-CD4_Tex-PDCD1.pdf", 
    ArchRProj = proj_cd4_tex, 
    addDOC = FALSE, width = 5, height = 5)
