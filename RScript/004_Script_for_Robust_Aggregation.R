# global options
rm(list = ls())
options(stringsAsFactors = FALSE)

# library
library(RobustRankAggreg)
library(pheatmap)
library(tidyverse)

# load data
loads <- function(...){
  file <- list.files(...)
  
  for (i in 1:length(file)) {
    load(file[i], envir = parent.frame(1))
  }
}

loads(path = "./RData/", pattern = "result_B.*", full.names = T)

# isolate DEGs
isolateDEGS <- function(dif){
  require(clusterProfiler)
  
  n <- (nrow(dif) * 0.025) %>% ceiling()
  
  dif <- na.omit(dif)
  colnames(dif)[grep("pvalue", colnames(dif))] <- "pValue"
  colnames(dif)[grep("P.Value", colnames(dif))] <- "pValue"
  colnames(dif)[grep("log2FoldChange", colnames(dif))] <- "log2FC"
  colnames(dif)[grep("logFC", colnames(dif))] <- "log2FC"
  
  dif <- dif[dif$pValue <= 0.05,]
  dif <- arrange(dif, desc(log2FC))
  
  if(sum(dif$log2FC > 0) >= n){
    dif_up <- dif[dif$log2FC > 0,] %>% head(n) %>% magrittr::extract2("symbol")
  }else{
    dif_up <- dif[dif$log2FC > 0,] %>% head(sum(dif$log2FC > 0)) %>% magrittr::extract2("symbol")
  }
  
  if(sum(dif$log2FC < 0) >= n){
    dif_down <- dif[dif$log2FC < 0,] %>% tail(n) %>% arrange(log2FC) %>% magrittr::extract2("symbol")
  }else{
    dif_down <- dif[dif$log2FC < 0,] %>% tail(sum(dif$log2FC < 0)) %>% arrange(log2FC) %>% magrittr::extract2("symbol")
  }
  
  #dif_up <- bitr(dif_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  #dif_down <- bitr(dif_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") 
  
  #dif_up <- dif_up$ENTREZID
  #dif_down <- dif_down$ENTREZID
  
  difList <- list(up = dif_up, down = dif_down)
  
  difList
  
}

difList_cell_1 <- isolateDEGS(dif = res_BC1_GSE56633$dif)
difList_cell_2 <- isolateDEGS(dif = res_BC2_GSE57896$dif)
difList_cell_3 <- isolateDEGS(dif = res_BC3_GSE71293$dif)
difList_cell_4 <- isolateDEGS(dif = res_BC4_GSE125331$dif)
difList_tissue_1 <- isolateDEGS(dif = res_BT1_GSE113764$dif)
difList_tissue_2 <- isolateDEGS(dif = res_BT2_GSE122721$dif)

glist_up <- list(difList_cell_1$up,
                 difList_cell_2$up,
                 difList_cell_3$up,
                 difList_cell_4$up,
                 difList_tissue_1$up,
                 difList_tissue_2$up)

glist_down <- list(difList_cell_1$down,
                   difList_cell_2$down,
                   difList_cell_3$down,
                   difList_cell_4$down,
                   difList_tissue_1$down,
                   difList_tissue_2$down)

ups <- aggregateRanks(glist = glist_up, N = length(unique(unlist(glist_up))))

tmp <- as.data.frame(table(unlist(glist_up)))
ups$Freq <- tmp[match(ups$Name,tmp[,1]),2]

ups <- ups[ups$Score < 0.05 & ups$Freq >=4,]$Name

downs <- aggregateRanks(glist = glist_down, N = length(unique(unlist(glist_down))))

tmp <- as.data.frame(table(unlist(glist_down)))
downs$Freq <- tmp[match(downs$Name,tmp[,1]),2]

downs <- downs[downs$Score < 0.05 & downs$Freq >=4,] %>% arrange(desc(Score)) %>% magrittr::extract2("Name")

mf <- c(ups,downs)

# heatmap
isolateDEGS <- function(dif){
  require(clusterProfiler)
  
  n <- (nrow(dif) * 0.025) %>% ceiling()
  
  dif <- na.omit(dif)
  colnames(dif)[grep("pvalue", colnames(dif))] <- "pValue"
  colnames(dif)[grep("P.Value", colnames(dif))] <- "pValue"
  colnames(dif)[grep("log2FoldChange", colnames(dif))] <- "log2FC"
  colnames(dif)[grep("logFC", colnames(dif))] <- "log2FC"
  
  dif <- dif[dif$pValue <= 0.05,]
  dif <- arrange(dif, desc(log2FC))
  
  if(sum(dif$log2FC > 0) >= n){
    dif_up <- dif[dif$log2FC > 0,] %>% head(n) %>% magrittr::extract2("symbol")
  }else{
    dif_up <- dif[dif$log2FC > 0,] %>% head(sum(dif$log2FC > 0)) %>% magrittr::extract2("symbol")
  }
  
  if(sum(dif$log2FC < 0) >= n){
    dif_down <- dif[dif$log2FC < 0,] %>% tail(n) %>% arrange(log2FC) %>% magrittr::extract2("symbol")
  }else{
    dif_down <- dif[dif$log2FC < 0,] %>% tail(sum(dif$log2FC < 0)) %>% arrange(log2FC) %>% magrittr::extract2("symbol")
  }
  
  #dif_up <- bitr(dif_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  #dif_down <- bitr(dif_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") 
  
  #dif_up <- dif_up$ENTREZID
  #dif_down <- dif_down$ENTREZID
  
  difList <- list(up = dif_up, down = dif_down)
  
  difList
  
}

mfList <- data.frame(BC1 = res_BC1_GSE56633$dif[mf,"logFC"],
                     BC2 = res_BC2_GSE57896$dif[mf,"log2FoldChange"],
                     BC3 = res_BC3_GSE71293$dif[mf,"logFC"],
                     BC4 = res_BC4_GSE125331$dif[mf,"log2FoldChange"],
                     BT1 = res_BT1_GSE113764$dif[mf,"logFC"],
                     BT2 = res_BT2_GSE122721$dif[mf,"logFC"],
                     row.names = mf)


bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

pdf(file = "./plot/heatmap_robust.pdf", width = 6, height = 6)
pheatmap(mfList,
         display_numbers = T,
         scale = "none",
         cluster_rows = F, cluster_cols = F,
         border_color = "black", number_color = "black",
         color = c(colorRampPalette(colors = c("#2376b7","white"))(length(bk)/2),colorRampPalette(colors = c("white","#f0a1a8"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk)
dev.off()






