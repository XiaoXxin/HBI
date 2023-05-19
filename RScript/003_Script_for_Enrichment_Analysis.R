# global options
rm(list = ls())
options(stringsAsFactors = FALSE)

# library
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(UpSetR)
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
  
  dif_up <- bitr(dif_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  dif_down <- bitr(dif_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") 
  
  dif_up <- dif_up$ENTREZID
  dif_down <- dif_down$ENTREZID
  
  difList <- list(up = dif_up, down = dif_down)
  
  difList
  
}

difList_cell_1 <- isolateDEGS(dif = res_BC1_GSE56633$dif)
difList_cell_2 <- isolateDEGS(dif = res_BC2_GSE57896$dif)
difList_cell_3 <- isolateDEGS(dif = res_BC3_GSE71293$dif)
difList_cell_4 <- isolateDEGS(dif = res_BC4_GSE125331$dif)
difList_tissue_1 <- isolateDEGS(dif = res_BT1_GSE113764$dif)
difList_tissue_2 <- isolateDEGS(dif = res_BT2_GSE122721$dif)

# up-regulated genes
geneList.up <- list(BC1 = difList_cell_1$up,
                    BC2 = difList_cell_2$up,
                    BC3 = difList_cell_3$up,
                    BC4 = difList_cell_4$up,
                    BT1 = difList_tissue_1$up,
                    BT2 = difList_tissue_2$up) 

# enrichment analysis
# note that this rscript was run in 24 May, 2022
# the database should be downloaded in 20220524 version for reproducibility
go.up <- compareCluster(geneList.up, 
                         fun = "enrichGO", 
                         OrgDb = "org.Hs.eg.db",
                         ont = "BP",
                         pAdjustMethod = "none",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.99,
                         minGSSize = 5,
                         maxGSSize = 1000
)

kegg.up <- compareCluster(geneList.up, 
                          fun = "enrichKEGG", 
                          organism = "hsa",
                          pAdjustMethod = "none",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.99,
                          minGSSize = 5,
                          maxGSSize = 1000
)

wp.up <- compareCluster(geneList.up, 
                          fun = "enrichWP", 
                        organism = "Homo sapiens",
                          pAdjustMethod = "none",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.99,
                          minGSSize = 5,
                          maxGSSize = 1000
)

color <- c("#0095C4", "#B44946", "#559575", "#6163A3", "#DB802D","#374E55")

go.up2 <- pairwise_termsim(go.up)
kegg.up2 <- pairwise_termsim(kegg.up)
wp.up2 <- pairwise_termsim(wp.up)

set.seed(2333)
p1 <- emapplot(go.up2, 
               showCategory = 10, 
               legend_n = 2, 
               pie = "count", 
               cex_category = 2,
               repel = T,
               cex_label_category = 0.7,
               node_label = "none")+
  scale_fill_manual(values = color)+
  labs(title = "GO enrichment analysis")+
  theme(plot.title = element_text(hjust = 0.5))

set.seed(2333)
p2 <- emapplot(kegg.up2, 
               showCategory = 10, 
               legend_n = 2, 
               pie = "count", 
               cex_category = 2,
               repel = T,
               cex_label_category = 0.7)+
  scale_fill_manual(values = color)+
  labs(title = "KEGG enrichment analysis")+
  theme(plot.title = element_text(hjust = 0.5))


pdf(file = "./plot/enrichment_up1.pdf", width = 5, height = 6)
p1
dev.off()

pdf(file = "./plot/enrichment_up2.pdf", width = 5, height = 6)
p2
dev.off()


