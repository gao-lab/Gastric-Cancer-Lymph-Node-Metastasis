# -*- coding: utf-8 -*-
cross_knn <- function(sce.query, sce.ref, hvgs, cosineNorm = T, npc=50, k=30,seed = 816)
{
  coldata.query <- as.data.frame(colData(sce.query))
  coldata.ref <- as.data.frame(colData(sce.ref))
  mt = cbind(logcounts(sce.query), logcounts(sce.ref))
  cellnames = colnames(mt)
  if(cosineNorm) {
    mt <- batchelor::cosineNorm(mt[hvgs,], mode = "matrix")
  }
  #--- PCA
  library(BiocNeighbors)
  set.seed(seed)
  cat("[INFO] Calculate PCA ...\n")
  pc <- irlba::prcomp_irlba(t(mt), n = npc)$x
  rownames(pc) <- cellnames
  pc.query <- pc[colnames(sce.query),]
  pc.ref <- pc[colnames(sce.ref),]
  cat("[INFO] Query KNN ...\n")
  N12 <- queryKNN(pc.ref, query=pc.query, k=k, get.distance=T, BPPARAM= BiocParallel::MulticoreParam(workers = 8))
  # Use ascites as reference
  cat("[INFO] Reverse Query KNN ...\n")
  N21 <- queryKNN(pc.query, query=pc.ref, k=k, get.distance=T, BPPARAM= BiocParallel::MulticoreParam(workers = 8))
  cat("[INFO] Find MNN ...\n")
  mnn.sets <- batchelor::findMutualNN(pc.ref, pc.query, k1 = k, k2 = k, BPPARAM= BiocParallel::MulticoreParam(workers = 8))
  N12$cellnames <- rownames(pc.query)
  N21$cellnames <- rownames(pc.ref)
  return(list(pca.ref = pc.ref, pca.query = pc.query, queryToRef = N12, refToQuery = N21, mnnSets = mnn.sets, coldata.ref = coldata.ref, coldata.query = coldata.query))
}

construct_crossknn_result <- function(sce.query, sce.ref, knn.result) {
  coldata.query <- as.data.frame(colData(sce.query))
  coldata.query$cell <- rownames(coldata.query)
  coldata.ref <- as.data.frame(colData(sce.ref))
  coldata.ref$cell <- rownames(coldata.ref)
  pc.ref.index <- knn.result$mnnSets$first
  pc.query.index <- knn.result$mnnSets$second
  query.mnn <- unique(pc.query.index)
  query.mnn.ordered <- sort(query.mnn)
  for(i in 1:30){
      coldata.query[, paste0("nn.cellId",i)] <- coldata.ref[knn.result$queryToRef$index[,i], "cell"]
  }
#   coldata.query[, "nn.cellId1"] <- coldata.ref[knn.result$queryToRef$index[,1], "cell"]
#   coldata.query[, "nn.cellId2"] <- coldata.ref[knn.result$queryToRef$index[,2], "cell"]
#   coldata.query[, "nn.cellId3"] <- coldata.ref[knn.result$queryToRef$index[,3], "cell"]
#   coldata.query[, "nn.cellId4"] <- coldata.ref[knn.result$queryToRef$index[,4], "cell"]
#   coldata.query[, "nn.cellId5"] <- coldata.ref[knn.result$queryToRef$index[,5], "cell"]
#   coldata.query[, "nn.cellId6"] <- coldata.ref[knn.result$queryToRef$index[,6], "cell"]
  coldata.query[, "nn.celltype"] <- coldata.ref[knn.result$queryToRef$index[,1], "celltype_sub"]
  coldata.query[, "nn.tissue"] <- coldata.ref[knn.result$queryToRef$index[,1], "tissue"]
  coldata.query[, "nn.donor"] <- coldata.ref[knn.result$queryToRef$index[,1], "patient"]
  coldata.query[, "nn.dist"] <- knn.result$queryToRef$distance[,1]
  coldata.query$cell.anchor <- ifelse(rownames(coldata.query) %in% rownames(knn.result$pca.query)[query.mnn.ordered], "anchor", "others")
  coldata.query <- coldata.query %>% 
    dplyr::mutate(nn.label = paste0(nn.celltype,"_", nn.tissue)) 
  return(coldata.query)
}

PhenoAligner_KNN <- function(sce_all, markers, cell_type="Myeloid", meta_tissue = "Met.LN",seed = 0)
{    
#     sce_all need to contain LN_condition meta, and subtypes in subtype1 column
    if(!file.exists(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type))){
        dir.create(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type))
        print("Creat directory ok!")
    }
    ### prepare input data (SingleCellExperiments object)
    sce_all$tissue <- sce_all$LN_condition
    table(sce_all$tissue)
    sce_all$celltype_sub <- sce_all$subtype1
    sce_all$celltype_major <- cell_type

    Myeloid_counts <- as.matrix(GetAssayData(sce_all, slot = "data"))
    Myeloid_meta.data <- sce_all@meta.data
    Myeloid_meta.data$barcode <- rownames(Myeloid_meta.data)
    
    cat("[INFO] Create SingleCellExperiment Object ...\n")
    sce.10x <- SingleCellExperiment(assays=list(counts = Myeloid_counts,logcounts = Myeloid_counts),
                                    colData = Myeloid_meta.data)
    #seperate data into query data and reference data
    sce.mt <- sce.10x[markers, sce.10x$tissue == meta_tissue]
    sce.others <- sce.10x[markers, sce.10x$tissue != meta_tissue]

    ###########################==================performe KNN method based on PCA result=================
    knn.list <- cross_knn(sce.mt, sce.others, markers,cosineNorm=T, npc=50, k=30,seed = seed)
    coldata.mt <- construct_crossknn_result(sce.query = sce.mt, sce.ref = sce.others, knn.result = knn.list)
    #--- Saving out results
    saveRDS(knn.list, paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/sce_all_results.mt2others.rds"))
    readr::write_tsv(coldata.mt, paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/sce_all_knn_results.mt2others.tsv"))
    return(list(sce.others = sce.others, sce.mt = sce.mt, coldata.mt = coldata.mt))
}

PhenoAligner_index <- function(sce.others, sce.mt, coldata.mt, 
                               cell_type="Myeloid", cor_cutoff=0.6)
{    
    ## filter KNN by corr and gap
    #calculate pearson correlation between each cell pairs from KNN result===============#################
    sce.ref <- sce.others
    sce.query <- sce.mt
    sce.mt<-assay(sce.mt,"logcounts")
    sce.mt<-sce.mt+0.1
    sce.others<-assay(sce.others,"logcounts")
    sce.others<-sce.others+0.1
    select_col <- c()
    for (i in 1:30){select_col <- append(select_col,(paste0("nn.cellId",i)))}

    cat("[INFO] Calculate Pearson correlation ...\n")
    result <- list()
    for(i in 1:ncol(sce.mt)){
      cell_ref<-coldata.mt[i,]$barcode
      y<-sce.mt[,cell_ref]
      cell_query<-as.character(coldata.mt[i,select_col]) #30 KNN cellID 
      query_mat<-sce.others[,cell_query]
      cor.mat<-Hmisc::rcorr(as.matrix(y),query_mat,type="pearson")
      cor.result<-data.frame(cor.mat$r[1,2:31])
      colnames(cor.result)<-"cor"
      cor.result$barcode <-rownames(cor.result)
      cor.result <-cor.result[order(cor.result$cor,decreasing = T),]
      pvalue<-data.frame(cor.mat$P[1,2:31])
      colnames(pvalue)<-"pvalue"
      pvalue$barcode <- rownames(pvalue)
      myresult <- plyr::join(cor.result,pvalue,by="barcode")
      result[[i]] <-myresult
    }
    print(length(result))
    saveRDS(result, paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/sce_all_COR_MAT.rds"))


    #filter correlation < cor_cutoff cell pairs & gap<0.01 for filter out cells---================#############################
    cat("[INFO] Filter by correlation and GAP ...\n")
    sce_all_coldata<- readRDS(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/sce_all_COR_MAT.rds"))
    cor_mat <- as.data.frame(do.call(rbind,sce_all_coldata))
    cor_mat$BH <- p.adjust(cor_mat$pvalue,method = "BH",n=length(cor_mat$pvalue))
    cor_mat$origin<-rep(colnames(sce.mt),each=30)#weil: origin is query barcode，此时还没有filter掉query
    cor_mat_keep <-cor_mat[which(cor_mat$cor>cor_cutoff),]
    print(dim(cor_mat))
    print(dim(cor_mat_keep))
    print(length(unique(cor_mat_keep$origin)))
    cor_mat_keep <-cor_mat_keep[which(cor_mat_keep$BH<0.05),]
    #这里是所有的query的结果一起filter，按照correlation filter

    # save only the nearest neighbor as nn.tissue, and add meta for nearest hit
    sce_all_phenoaligner <-  cor_mat_keep %>% 
      group_by(origin) %>% 
      dplyr::mutate(cor_max=max(cor),cor_max2=sort(cor,decreasing = TRUE)[2]) %>%
      dplyr::filter(cor==max(cor))
    map_data <-colData(sce.ref)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
    map_data <- data.frame(map_data)
    sce_all_phenoaligner <-left_join(sce_all_phenoaligner,map_data,by="barcode")
    colnames(sce_all_phenoaligner) <- c("cor","barcode","pvalue","BH","origin","cor_max","cor_max2",
                                    "nn.celltype_sub","nn.celltype_major","nn.tissue","nn.sample","nn.patient")
    #nn指的是reference中对应query的nearest neighbor，这里只保留了correlation最大的hit

    map_data <-colData(sce.query)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
    #注意我的数据中的meta列名，cell_type改为celltype_sub，celltype_major加上，标记为他的大类，加上tissue，内容为LN_condition

    map_data <- data.frame(map_data)
    colnames(map_data)[1]<-"origin"
    sce_all_phenoaligner <-left_join(sce_all_phenoaligner,map_data,by="origin")#这里加上的是query的meta信息
    head(sce_all_phenoaligner)
    write.csv(sce_all_phenoaligner,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_PheonAligner.csv"))

    #############================add gap<0.01 as filter cutoff
    map_data <-colData(sce.ref)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
    map_data <- data.frame(map_data)
    different_tissue <-left_join(cor_mat_keep,map_data,by="barcode")
    colnames(different_tissue) <- c("cor","barcode","pvalue","BH","origin","nn.celltype_sub",
                                    "nn.celltype_major","nn.tissue","nn.sample","nn.patient")
    sce_all_tissue_max_cor<-different_tissue %>% 
      group_by(origin, nn.tissue) %>% 
      summarise(tissue_max = max(cor),cells=n())#计算对每个query来自每个tissue的各自的最大的correlation，以及该tissue出现次数。
    sce_all_tissue_max_diff<-sce_all_tissue_max_cor %>% 
      group_by(origin) %>% 
      dplyr::mutate(top_one=max(tissue_max),top_two=sort(tissue_max,decreasing = TRUE)[2])
    #给每个query，给出各个tissue中correlation最大值最大的两个tissue，用来计算gap
    sce_all_tissue_max_diff<-data.frame(sce_all_tissue_max_diff)
    sce_all_tissue_max_diff$gap <- sce_all_tissue_max_diff$top_one-sce_all_tissue_max_diff$top_two
    write.csv(sce_all_tissue_max_diff,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_tissue_max_gap.csv"))
    #这里加上了meta，还没有用gap进行filter
    sce_all_tissue_max_diff <- sce_all_tissue_max_diff[!duplicated(sce_all_tissue_max_diff$origin),] 
    print(dim(sce_all_tissue_max_diff))
    ##[1] 5324    7
    sce_all_tissue_max_diff<-sce_all_tissue_max_diff[-which(sce_all_tissue_max_diff$gap<0.01),]#按照gap filter
    print(dim(sce_all_tissue_max_diff))
    ##[1] 3582    7
    # coldata.mt<- read.csv(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_PheonAligner.csv",row.names = 1)
    coldata.mt <- sce_all_phenoaligner
    coldata.mt <- coldata.mt[which(coldata.mt$origin%in%sce_all_tissue_max_diff$origin),]
    #按照gap filter，此时是filter之后的query，每个query对应一个correlation最大的hit
    write.csv(coldata.mt,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_PhenoAligner_with_gap.csv"))


    # barplot of NN tissue only by the first KNN cell

    ################==========================generate barplot=============##########################
#     cells.count <- as.data.frame(colData(sce.10x)) %>%
#       group_by(tissue, celltype_sub) %>% 
#       dplyr::summarise(cells = n())
#     cells.count$celltype_sub <- as.character(cells.count$celltype_sub)
#     df <- inner_join(coldata.mt, cells.count , by = c("nn.tissue" = "tissue", "nn.celltype_sub" = "celltype_sub"))
#     plotdata  <- df %>% 
    plotdata  <- coldata.mt %>% 
#       dplyr::mutate(nn.celltype = nn.celltype_sub)  %>% 
      group_by(nn.tissue, celltype_sub) %>% 
      dplyr::summarise(cells = n())  %>% 
      group_by(celltype_sub) %>% 
      dplyr::mutate(frac = cells/sum(cells))#每个亚型返回的最近邻组织中各个组织的占比

    p <- ggplot(plotdata, aes(x = nn.tissue, y = frac)) +
      geom_bar(aes(fill = nn.tissue), stat = "identity") +
      geom_text(aes(y = ifelse(frac <= .1, 0.1,  frac-.1) , label = paste0("N=", cells)), size = 3) +
      scale_y_continuous(limits = c(0, 1), labels =  scales::percent) +
      #scale_fill_manual(values = tissue.cols) +
      facet_wrap(. ~  celltype_sub, ncol = 4) +
      theme_bw(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            strip.text.y =  element_text(size = 10))+
      ylab("Fraction") +
      xlab("") + ggtitle("Merged Donors")

    #"#E64B35B2" "#4DBBD5B2" "#00A087B2" "#3C5488B2" "#F39B7FB2" "#8491B4B2" #"#91D1C2B2" "#sce_all0000B2" "#7E6148B2"
    #"#1f77b4"  "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b" "#e377c2" "#7f7f7f" "#bcbd22" "#17becf"
    #"#aec7e8" "#ffbb78" "#98df8a" "#ff9896" "#c5b0d5" "#c49c94" "#f7b6d2" "#c7c7c7" "#dbdb8d" "#9edae5"
    human_tissue_color_panel <- c("PBMC"="#1f77b4","Neg.LN"="#ff7f0e","Adj.Nor"="#2ca02c","Pri.GC"="#9467bd")
    #                               ,"metastasis normal"="#9467bd","metastasis tumor"="#8491B4B2")
    pdf(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_result.pdf"), height = 12, width = 12)
    plot(p + scale_fill_manual(values = human_tissue_color_panel))
    dev.off()
    # #没有被注释到的细胞（被filter掉的）
    # test <- colData(sce.query)
    # cells <- setdiff(test$barcode, coldata.mt$origin)
    # test <- test[which(test$barcode%in%cells),]
    # dim(test)
    # table(test$celltype_sub)



    ##########################--- mean top 3 cells correlation to difine nearest tissue---###########################
    ############################---test mean top 3 cells correlation as gap 
    cat("[INFO] Calculate mean top 3 cells correlation ...\n")
    sce_all_coldata<- readRDS(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/sce_all_COR_MAT.rds"))

    cor_mat <- as.data.frame(do.call(rbind,sce_all_coldata))
    cor_mat$BH <- p.adjust(cor_mat$pvalue,method = "BH",n=length(cor_mat$pvalue))
    cor_mat$origin<-rep(colnames(sce.mt),each=30)
    cor_mat_keep <-cor_mat[cor_mat$origin %in% coldata.mt$origin,]
    # use corr and gap filtered results

    #############================calc top 3 mean corr
    map_data <-colData(sce.ref)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
    map_data <- data.frame(map_data)
    different_tissue <-left_join(cor_mat_keep,map_data,by="barcode")
    colnames(different_tissue) <- c("cor","barcode","pvalue","BH","origin","nn.celltype_sub",
                                    "nn.celltype_major","nn.tissue","nn.sample","nn.patient")
    #sce_all_tissue_max_cor<-different_tissue %>% group_by(origin, nn.tissue) %>% summarise(tissue_max = max(cor),cells=n())
    sce_all_tissue_max_cor<-different_tissue %>% 
      group_by(origin, nn.tissue) %>% 
      dplyr::mutate(cells=n(),top_one=max(cor),
                    top_two=sort(cor,decreasing = TRUE)[2],top_three=sort(cor,decreasing = TRUE)[3])
    sce_all_tissue_max_cor$mean <- (sce_all_tissue_max_cor$top_one+sce_all_tissue_max_cor$top_two+sce_all_tissue_max_cor$top_three)/3
    sce_all_tissue_max_cor$nn.lable <- paste(sce_all_tissue_max_cor$origin,sce_all_tissue_max_cor$nn.tissue,sce_all_tissue_max_cor$mean,
                                         sep='_')
    sce_all_tissue_max_cor <- data.frame(sce_all_tissue_max_cor)
    dim(sce_all_tissue_max_cor)
    sce_all_tissue_max_cor <- sce_all_tissue_max_cor[!duplicated(sce_all_tissue_max_cor$nn.lable),]
    #每个query每个tissue都取前三个corr的均值
    dim(sce_all_tissue_max_cor)
    #[1] 11938    17
    sce_all_tissue_max_mean<-sce_all_tissue_max_cor %>% group_by(origin) %>% dplyr::mutate(mean_max=max(mean,na.rm = TRUE))
    ###每个query每个tissue都取前三个corr的均值，然后用均值最大的组织作为结果
    dim(sce_all_tissue_max_mean)
    # 524317 17
    # write.csv(sce_all_tissue_max_mean,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/test_B.csv"))
    sce_all_tissue_max_diff<-sce_all_tissue_max_mean[which(sce_all_tissue_max_mean$mean==sce_all_tissue_max_mean$mean_max),]
    print(dim(sce_all_tissue_max_diff))
    # write.csv(sce_all_tissue_max_diff,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/test_B_diff.csv"))

    cells <- setdiff(unique(sce_all_tissue_max_mean$origin),unique(sce_all_tissue_max_diff$origin))
    ### 81
    #这里我理解的是保留了mean最大的组织，其他的组织丢掉，但这里是想说有些origin也会丢掉？结果应该是0吧。
    #之前师姐给的代码这里只用了corr filter，没有用gap filter，那样的时候，在这里的diff有81个结果，但不太明白是怎么来的diff，
    #而且结合她文章里面的说法，应该是corr和gap都filter之后，再计算top3平均最高的作为nn.tissue
    length(cells)

    #coldata.mt<- read.csv(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/sce_all_PheonAligner.csv",row.names = 1)
    #coldata.mt <- coldata.mt[which(coldata.mt$origin%in%sce_all_tissue_max_diff$origin),]
    map_data <-colData(sce.query)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
    map_data <- data.frame(map_data)
    colnames(map_data)[1]<-"origin"
    sce_all_tissue_max_diff <-left_join(sce_all_tissue_max_diff,map_data,by="origin")
    head(sce_all_tissue_max_diff)
    coldata.mt <- sce_all_tissue_max_diff
    write.csv(sce_all_tissue_max_diff,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_mean_top3_result.csv"))
    #最后把结果导出，然后看每个亚型的所有细胞，分别被标记为各个组织的占比。
    #现在这里的nn.tissue就是根据前三个均值最大得到的tissue
    #注意他文章里说的是标记到每个组织的数目占该亚型总数目的占比，但实际上是总数，还是filter之后的有标记结果的总数还不好说。

    ## barplot of nearest tissue by the top 3 KNN cell
    ################==========================generate barplot=============##########################
#     cells.count <- as.data.frame(colData(sce.10x)) %>%
#       group_by(tissue, celltype_sub) %>% 
#       dplyr::summarise(cells = n())
#     cells.count$celltype_sub <- as.character(cells.count$celltype_sub)
#     df <- inner_join(coldata.mt, cells.count , by = c("nn.tissue" = "tissue", "nn.celltype_sub" = "celltype_sub"))
#     plotdata  <- df %>% 
    plotdata  <- coldata.mt %>% 
#       dplyr::mutate(nn.celltype = nn.celltype_sub)  %>% 
      group_by(nn.tissue, celltype_sub) %>% 
      dplyr::summarise(cells = n())  %>% 
      group_by(celltype_sub) %>% 
      dplyr::mutate(frac = cells/sum(cells))#每个亚型返回的最近邻组织中各个组织的占比
    write.csv(plotdata, paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/PhenoAligner_index_final.csv"))
    
    
    p <- ggplot(plotdata, aes(x = nn.tissue, y = frac)) +
      geom_bar(aes(fill = nn.tissue), stat = "identity") +
      geom_text(aes(y = ifelse(frac <= .1, 0.1,  frac-.1) , label = paste0("N=", cells)), size = 3) +
      scale_y_continuous(limits = c(0, 1), labels =  scales::percent) +
      #scale_fill_manual(values = tissue.cols) +
      facet_wrap(. ~  celltype_sub, ncol = 4) +
      theme_bw(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            strip.text.y =  element_text(size = 10))+
      ylab("Fraction") +
      xlab("") + ggtitle("Merged Donors")

    human_tissue_color_panel <- c("PBMC"="#1f77b4","Neg.LN"="#ff7f0e","Adj.Nor"="#2ca02c","Pri.GC"="#9467bd")
    #                               ,"metastasis normal"="#9467bd","metastasis tumor"="#8491B4B2")
    pdf(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_top3_result.pdf"), height = 12, width = 12)
    plot(p + scale_fill_manual(values = human_tissue_color_panel))
    dev.off()

    return(plotdata)
}


# shuffle label for null distribution
# PhenoAligner_shuffle <- NULL

PhenoAligner_shuffle <- function(sce.others, sce.mt, 
                               cell_type="Myeloid", cor_cutoff=0.6, threads=10)
{   
    #filter correlation < cor_cutoff cell pairs & gap<0.01 for filter out cells---================#############################
    cat("[INFO] Filter by correlation ...\n")
    sce_all_coldata<- readRDS(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/sce_all_COR_MAT.rds"))
    cor_mat <- as.data.frame(do.call(rbind,sce_all_coldata))
    cor_mat$BH <- p.adjust(cor_mat$pvalue,method = "BH",n=length(cor_mat$pvalue))
    cor_mat$origin<-rep(colnames(sce.mt),each=30)#weil: origin is query barcode，此时还没有filter掉query
    cor_mat_keep_0 <-cor_mat[which(cor_mat$cor>cor_cutoff),]
    print(dim(cor_mat))
    print(dim(cor_mat_keep_0))
    print(length(unique(cor_mat_keep_0$origin)))
    cor_mat_keep_0 <-cor_mat_keep_0[which(cor_mat_keep_0$BH<0.05),]
    #这里是所有的query的结果一起filter，按照correlation filter
    # cor_mat, cor_mat_keep_0 are not changed during the for loop
    
    run_shuffle <- function(random_seed) {
#     for(random_seed in 1:1000){# only random_seed changes during the for loop
        print(random_seed)
        sce.ref <- sce.others
        sce.query <- sce.mt
        set.seed(random_seed) # shuffle reference label   
        sce.ref$tissue <- sample(sce.ref$tissue)
        
        # save only the nearest neighbor as nn.tissue, and add meta for nearest hit
        sce_all_phenoaligner <-  cor_mat_keep_0 %>% 
          group_by(origin) %>% 
          dplyr::mutate(cor_max=max(cor),cor_max2=sort(cor,decreasing = TRUE)[2]) %>%
          dplyr::filter(cor==max(cor))
        map_data <-colData(sce.ref)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
        map_data <- data.frame(map_data)
        sce_all_phenoaligner <-left_join(sce_all_phenoaligner,map_data,by="barcode")
        colnames(sce_all_phenoaligner) <- c("cor","barcode","pvalue","BH","origin","cor_max","cor_max2",
                                        "nn.celltype_sub","nn.celltype_major","nn.tissue","nn.sample","nn.patient")
        #nn指的是reference中对应query的nearest neighbor，这里只保留了correlation最大的hit

        map_data <-colData(sce.query)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
        map_data <- data.frame(map_data)
        colnames(map_data)[1]<-"origin"
        sce_all_phenoaligner <-left_join(sce_all_phenoaligner,map_data,by="origin")#这里加上的是query的meta信息
        # write.csv(sce_all_phenoaligner,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_PheonAligner.csv"))

        #############================add gap<0.01 as filter cutoff
        cat("[INFO] Filter by GAP ...\n")
        map_data <-colData(sce.ref)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
        map_data <- data.frame(map_data)
        different_tissue <-left_join(cor_mat_keep_0,map_data,by="barcode")
        colnames(different_tissue) <- c("cor","barcode","pvalue","BH","origin","nn.celltype_sub",
                                        "nn.celltype_major","nn.tissue","nn.sample","nn.patient")
        sce_all_tissue_max_cor<-different_tissue %>% 
          group_by(origin, nn.tissue) %>% 
          summarise(tissue_max = max(cor),cells=n())#计算对每个query来自每个tissue的各自的最大的correlation，以及该tissue出现次数。
        sce_all_tissue_max_diff<-sce_all_tissue_max_cor %>% 
          group_by(origin) %>% 
          dplyr::mutate(top_one=max(tissue_max),top_two=sort(tissue_max,decreasing = TRUE)[2])
        #给每个query，给出各个tissue中correlation最大值最大的两个tissue，用来计算gap
        sce_all_tissue_max_diff<-data.frame(sce_all_tissue_max_diff)
        sce_all_tissue_max_diff$gap <- sce_all_tissue_max_diff$top_one-sce_all_tissue_max_diff$top_two
        # write.csv(sce_all_tissue_max_diff,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_tissue_max_gap.csv"))
        #这里加上了meta，还没有用gap进行filter
        sce_all_tissue_max_diff <- sce_all_tissue_max_diff[!duplicated(sce_all_tissue_max_diff$origin),] 
#         dim(sce_all_tissue_max_diff)
        sce_all_tissue_max_diff<-sce_all_tissue_max_diff[-which(sce_all_tissue_max_diff$gap<0.01),]#按照gap filter
#         print(dim(sce_all_tissue_max_diff)[1])
        # coldata.mt<- read.csv(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_PheonAligner.csv",row.names = 1)
        coldata.mt <- sce_all_phenoaligner
        coldata.mt <- coldata.mt[which(coldata.mt$origin%in%sce_all_tissue_max_diff$origin),]
        #按照gap filter，此时是filter之后的query，每个query对应一个correlation最大的hit
        # write.csv(coldata.mt,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/global_PhenoAligner_with_gap.csv"))


        ##########################--- mean top 3 cells correlation to difine nearest tissue---###########################
        cat("[INFO] Calculate mean top 3 cells correlation ...\n")
        cor_mat_keep <-cor_mat[cor_mat$origin %in% coldata.mt$origin,]# use corr and gap filtered results

        #############================calc top 3 mean corr
        map_data <-colData(sce.ref)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
        map_data <- data.frame(map_data)
        different_tissue <-left_join(cor_mat_keep,map_data,by="barcode")
        colnames(different_tissue) <- c("cor","barcode","pvalue","BH","origin","nn.celltype_sub",
                                        "nn.celltype_major","nn.tissue","nn.sample","nn.patient")
        #sce_all_tissue_max_cor<-different_tissue %>% group_by(origin, nn.tissue) %>% summarise(tissue_max = max(cor),cells=n())
        sce_all_tissue_max_cor<-different_tissue %>% 
          group_by(origin, nn.tissue) %>% 
          dplyr::mutate(cells=n(),top_one=max(cor),
                        top_two=sort(cor,decreasing = TRUE)[2],top_three=sort(cor,decreasing = TRUE)[3])
        sce_all_tissue_max_cor$mean <- (sce_all_tissue_max_cor$top_one+sce_all_tissue_max_cor$top_two+sce_all_tissue_max_cor$top_three)/3
        sce_all_tissue_max_cor$nn.lable <- paste(sce_all_tissue_max_cor$origin,sce_all_tissue_max_cor$nn.tissue,sce_all_tissue_max_cor$mean,
                                             sep='_')
        sce_all_tissue_max_cor <- data.frame(sce_all_tissue_max_cor)
#         dim(sce_all_tissue_max_cor)
        sce_all_tissue_max_cor <- sce_all_tissue_max_cor[!duplicated(sce_all_tissue_max_cor$nn.lable),]
        #每个query每个tissue都取前三个corr的均值
#         dim(sce_all_tissue_max_cor)
        sce_all_tissue_max_mean<-sce_all_tissue_max_cor %>% group_by(origin) %>% dplyr::mutate(mean_max=max(mean,na.rm = TRUE))
        ###每个query每个tissue都取前三个corr的均值，然后用均值最大的组织作为结果
#         dim(sce_all_tissue_max_mean)
        # # write.csv(sce_all_tissue_max_mean,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/test_B.csv"))
        sce_all_tissue_max_diff<-sce_all_tissue_max_mean[which(sce_all_tissue_max_mean$mean==sce_all_tissue_max_mean$mean_max),]
#         print(dim(sce_all_tissue_max_diff))
        # # write.csv(sce_all_tissue_max_diff,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/test_B_diff.csv"))

        map_data <-colData(sce.query)[,c("barcode","celltype_sub","celltype_major","tissue","sample","patient")]
        map_data <- data.frame(map_data)
        colnames(map_data)[1]<-"origin"
        sce_all_tissue_max_diff <-left_join(sce_all_tissue_max_diff,map_data,by="origin")
        coldata.mt <- sce_all_tissue_max_diff
        coldata.mt$nn.tissue <- factor(coldata.mt$nn.tissue, levels = c("Neg.LN", "Adj.Nor", "PBMC", "Pri.GC"))
        
        ## barplot of nearest tissue by the top 3 KNN cell
        plotdata  <- coldata.mt %>% 
          group_by(nn.tissue, celltype_sub) %>% 
          dplyr::summarise(cells = n())  %>% 
          group_by(celltype_sub) %>% 
          dplyr::mutate(frac = cells/sum(cells))#每个亚型返回的最近邻组织中各个组织的占比
        # write.csv(plotdata, paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/PhenoAligner_index_final.csv"))
      
#         pheno_index_list[[random_seed]] <- plotdata
        return(plotdata)
    }

    
    time_start <- Sys.time()
    bpparam <- BiocParallel::MulticoreParam(threads, log = FALSE, stop.on.error = TRUE)
    pheno_index_list <- BiocParallel::bplapply(1:1000, run_shuffle, BPPARAM=bpparam) 
    time_end <- Sys.time()
    time_use <- time_end - time_start
    print(time_use)
    saveRDS(pheno_index_list, paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/pheno_shuffle.list.rds"))
    
    return(pheno_index_list) 
}

plot_pheno <- function(pheno_result,pheno_shuffle.list,cell_type="mph"){    
    # merge shuffle results
    tmp_list <- pheno_shuffle.list
    for(i in 1:length(tmp_list)){
#         print(dim(tmp_list[[i]]))
        tmp_list[[i]] <- tmp_list[[i]][,c(1,2,4)]
        colnames(tmp_list[[i]])[3] <- paste0("frac_",i)
    }
    index_df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("nn.tissue", "celltype_sub"), all = TRUE), tmp_list)
    index_df[is.na(index_df)] <- 0
    dim(index_df)
                   
#     pdf(paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/I_phenoaligner_by_subtype.pdf"), height = 4, width = 4)
    plot_df <- list()
    for(subtype in unique(pheno_result$celltype_sub)){
        # shuffle index
        shuffle_index <- t(index_df[index_df$celltype_sub==subtype,])
        colnames(shuffle_index) <- shuffle_index["nn.tissue",]
        shuffle_index <- shuffle_index[3:dim(shuffle_index)[1],]

        shuffle_index <- melt(shuffle_index)
        colnames(shuffle_index) <- c("random","Tissue","frac_shuffle")
        shuffle_index$frac_shuffle <- as.numeric(shuffle_index$frac_shuffle)

        # final index
        final_index <- t(pheno_result[pheno_result$celltype_sub==subtype,])
        colnames(final_index) <- final_index["nn.tissue",]

        final_index_2 <- final_index["frac",,drop=FALSE]
        final_index_2 <- melt(final_index_2[])
        colnames(final_index_2) <- c("random","Tissue","I_phenoaligner")
        final_index_2$Tissue <- gsub("TRUE","Pri.GC",final_index_2$Tissue)
        final_index_2$cells <- as.numeric(final_index["cells",])
        final_index_2$I_phenoaligner <- as.numeric(final_index_2$I_phenoaligner)
        final_index_2$Tissue <- as.character(final_index_2$Tissue)

        # calculate p-value
        for(tissue in c("Neg.LN", "Adj.Nor", "PBMC", "Pri.GC")){
            if(sum(final_index_2$Tissue==tissue)==0){
                next()
            }else{
                test_p <- pnorm(final_index_2$I_phenoaligner[final_index_2$Tissue==tissue], 
                                mean=mean(shuffle_index$frac_shuffle[shuffle_index$Tissue==tissue]), 
                               sd=sd(shuffle_index$frac_shuffle[shuffle_index$Tissue==tissue]),
                               lower.tail=FALSE)
                final_index_2[final_index_2$Tissue==tissue, "p.value"] <- test_p
                if(test_p < 0.001){
                    sig <- "***"
                }else if(test_p < 0.01){
                    sig <- "**"
                }else if(test_p < 0.05){
                    sig <- "*"
                }else{
                    sig <- NA
                }
                final_index_2[final_index_2$Tissue==tissue, "significance"] <- sig
            }
        }

        plot_df[[subtype]] <- merge(shuffle_index, final_index_2,by="Tissue",all.x=TRUE)

        # options(repr.plot.width=4, repr.plot.height=4)
#         p <- ggplot(plot_df[[subtype]], aes(x=Tissue, y=frac_shuffle,fill=Tissue)) + 
#           geom_violin(trim=FALSE,na.rm = TRUE, scale = "width")+ 
#           geom_point(aes(y = I_phenoaligner), size = 2)+
#           geom_text(aes(y = I_phenoaligner+0.05 , label = paste0("N=", cells)), size = 3)+
#           geom_text(aes(y = I_phenoaligner-0.05, label = significance), size = 4) +
#           labs(title=subtype,x="Tissue", y = "I_phenoaligner")+
#           scale_y_continuous(limits = c(-0.05, 1.05))+
#           theme_classic()
#         plot(p)    
    }
#     dev.off()
    saveRDS(plot_df,paste0("processed_data/data_B2-19/PhenoAligner/", cell_type, "/plots.list.rds"))
    return(plot_df)
}
