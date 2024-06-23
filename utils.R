# +
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

# GO enrichment function
GO_enrichment <- function(markers){
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
        )
    }

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
    return(list(up_BP = up_BP, 
                up_MF = up_MF, 
                up_CC = up_CC,
                up_KEGG = up_KEGG,
                down_BP = down_BP, 
                down_MF = down_MF, 
                down_CC = down_CC,
                down_KEGG = down_KEGG
               ))
}

plot_try <- function(x){
    tryCatch(
        expr = {
            plot(x)
#             message("Successfully executed the plot(x) call.")
        },
        error = function(e){
            message('Caught an error!')
#             print(e)
        },
        warning = function(w){
            message('Caught an warning!')
#             print(w)
        }
    )    
}
plot_GO <- function(enrich.list, condition_name,show_category=12){ #condition_name is only used as pdf name here
    up_BP = enrich.list$up_BP 
    up_MF = enrich.list$up_MF 
    up_CC = enrich.list$up_CC
    up_KEGG = enrich.list$up_KEGG
    down_BP = enrich.list$down_BP 
    down_MF = enrich.list$down_MF 
    down_CC = enrich.list$down_CC
    down_KEGG = enrich.list$down_KEGG
    
    if(!is.null(up_BP)){
        p1<-barplot(up_BP, showCategory=show_category, title=paste0(condition_name, " up BP"))
        plot_try(p1)
    }
    if(!is.null(up_MF)){
        p2<-barplot(up_MF, showCategory=show_category, title=paste0(condition_name, " up MF"))
        plot_try(p2)
    }
    if(!is.null(up_CC)){
        p3<-barplot(up_CC, showCategory=show_category, title=paste0(condition_name, " up CC"))  
        plot_try(p3)
    }
    if(!is.null(up_KEGG)){
        p4<-barplot(up_KEGG, showCategory=show_category, title=paste0(condition_name, " up KEGG"))  
        plot_try(p4)
    }
    if(!is.null(down_BP)){
        p5<-barplot(down_BP, showCategory=show_category, title=paste0(condition_name, " down BP"))
        plot_try(p5)
    }
    if(!is.null(down_MF)){
        p6<-barplot(down_MF, showCategory=show_category, title=paste0(condition_name, " down MF"))
        plot_try(p6)
    }
    if(!is.null(down_CC)){
        p7<-barplot(down_CC, showCategory=show_category, title=paste0(condition_name, " down CC"))
        plot_try(p7)
    }
    if(!is.null(down_KEGG)){
        p8<-barplot(down_KEGG, showCategory=show_category, title=paste0(condition_name, " down KEGG"))
        plot_try(p8)
    }
}
