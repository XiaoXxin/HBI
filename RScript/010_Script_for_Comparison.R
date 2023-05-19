# library
library(GEOquery)
library(limma)
library(oligo)
library(GSVA)
library(R.utils)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(paletteer)
library(FactoMineR)
library(factoextra)
library(sva)
library(tidyverse)

#####################
#### Cell cohort ####
#####################

# global options
rm(list = ls())
options(stringsAsFactors = FALSE)

# download data
eset <- getGEO("GSE71293", destdir="./rawData/GSE71293", AnnotGPL = F, getGPL = F)
eset <- eset[[1]]
phenoData <- pData(eset)
eset <- exprs(eset)
eset[is.na(eset)] = 0

dim(eset)
keep1 <- apply(eset, 1, function(x) sum(is.na(x))) <= 6
eset <- eset[keep1,]
dim(eset)

# prepare sample annotations
phenoData <- data.frame(sample = phenoData$geo_accession, group = phenoData$title)

phenoData$group <- gsub("Control.*", "White", phenoData$group)
phenoData$group <- gsub("Rosiglitazone.*", "Beige", phenoData$group)
phenoData$group <- gsub("GW7647.*", "Beige", phenoData$group)

# prepare gene annotations
library(HsAgilentDesign026652.db)
IDs <- toTable(HsAgilentDesign026652ENSEMBL)

# remove duplicated genes
dim(eset)
eset <- eset[rownames(eset) %in% IDs$probe_id, ]
dim(eset)
IDs <- IDs[match(rownames(eset), IDs$probe_id), ]
eset <- avereps(eset, ID = IDs$ensembl_id)
dim(eset)

write.table(eset, file = "BC3.txt", sep = "\t", col.names = T, row.names = T, quote = F)

# analysis using ProFAT
res <- read.csv(file = "./prediction_BC3.csv", header = T, row.names = 1)
res$Group <- factor(phenoData$group, levels = c("White", "Beige"))

p <- ggplot(res, aes(x = Group, y = BAT, fill = Group))+
  geom_boxplot()+
  theme_classic()+
  labs(x = NULL, y = "Prediction")+
  ggtitle("Cell cohort (BC3)")+
  scale_fill_manual(values = c("#559575", "#6163A3"))+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), 
                     method = "t.test",  
                     method.args = list(p.adjust.method = "BH"))+
  theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_line(size = 0.3),
        plot.title = element_text(hjust = 0.5))


pdf(file = "./plot/method_comparison1.pdf", width = 4, height = 3)
p
dev.off()

#######################
#### Tissue cohort ####
#######################

# global options
rm(list = ls())

# load  raw data
celFiles <-  list.files(path = "./rawData/GSE122721", pattern = "\\.gz$", full.names = T)
rawData <- oligo::read.celfiles(celFiles)

# rename
sampleNames(rawData)
sampleNames(rawData) <- gsub("_.*", "", sampleNames(rawData))
sampleNames(rawData)

# call rma function
rawData.rma <- oligo::rma(rawData)

# isolate matrix
eset <- exprs(rawData.rma)

# prepare sample annotations
phenoData <- data.frame(sample = colnames(eset), group = rep(c("WAT", "BAT"), 3))

# prepare gene annotations
library(hugene10sttranscriptcluster.db)
IDs <- toTable(hugene10sttranscriptclusterENSEMBL)


# remove duplicated genes
dim(eset)
eset <- eset[rownames(eset) %in% IDs$probe_id, ]
dim(eset)
IDs <- IDs[match(rownames(eset), IDs$probe_id), ]
eset <- avereps(eset, ID = IDs$ensembl_id)
dim(eset)

write.table(eset, file = "BT2.txt", sep = "\t", col.names = T, row.names = T, quote = F)

# analysis using ProFAT
res <- read.csv(file = "prediction_BT2.csv", header = T, row.names = 1)
res$Group <- factor(phenoData$group, levels = c("WAT", "BAT"))

p <- ggplot(res, aes(x = Group, y = BAT, fill = Group))+
  geom_boxplot()+
  theme_classic()+
  labs(x = NULL, y = "Prediction")+
  ggtitle("Tissue cohort (BT2)")+
  scale_fill_manual(values = c("#559575", "#6163A3"))+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), 
                     method = "t.test",  
                     method.args = list(p.adjust.method = "BH"))+
  theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_line(size = 0.3),
        plot.title = element_text(hjust = 0.5))


pdf(file = "./plot/method_comparison2.pdf", width = 4, height = 3)
p
dev.off()

####################
#### WAT cohort ####
####################

# global options
rm(list = ls())
options(stringsAsFactors = FALSE)

# load data
load("./RData/cohort_GTEx_AT.RData")

# isolate SAT samples
phenoData <- phenoData[phenoData$group == "SAT",]
counts <- counts[,phenoData$sample]

rownames(counts) <- gsub("[.].*", "", rownames(counts))

set.seed(123)
index <- createFolds(colnames(counts), k = 10, list = F) 
table(index)

for (i in 1:10) {
  write.table(counts[,index == i], file = paste0("GTEX_SAT_", i, ".txt"), sep = "\t", col.names = T, row.names = T, quote = F)
}

# analysis using ProFAT
res <- data.frame()
for (i in 1:10) {
  resi <- read.csv(file = paste0("prediction_", i, ".csv"), header = T, row.names = 1)
  res <- rbind(res, resi)
}

rownames(res) <- gsub("_R00", "", rownames(res))
rownames(res) <- gsub("[.]", "-", rownames(res))
res_ProFAT <- res[phenoData$sample,]

save(res_ProFAT, file = "./RData/result_ProFAT.RData")

rm(list = ls())

# load data
load("./RData/result_GTEx_SAT.RData")
load("./RData/result_score.RData")
load("./RData/result_ProFAT.RData")

# HBI
meanNor <- function(x){(x-mean(x))/(max(x)-min(x))}

matrix <- exprSet[c("DHRS11", "PREX1", "REEP6", "STX11"),] %>% t() %>% apply(2, meanNor) %>% as.data.frame()
res <- predict(fit, matrix, interval = "prediction",level = 0.95) %>% as.data.frame() %>% magrittr::set_colnames("fit") %>% t()


# absHBI
load("./RData/result_abs_models.RData")
load("./RData/result_abs_background.RData")

absHBIpre <- function(exprSet, background){
  
  locAll <- c()
  for (i in names(background)) {
    bg <- c(background[[i]], i) %>% unique()
    exp <- exprSet[rownames(exprSet) %in% bg,]
    size <- nrow(exp)
    locList <- c()
    
    for (j in 1:ncol(exp)) {
      exprSelect <- exp[,j] %>% 
        as.data.frame() %>% 
        magrittr::set_colnames("sample") %>% 
        arrange(sample) %>% 
        rownames_to_column("symbol")
      
      geneloc <- grep(i, exprSelect$symbol)
      geneloc <- (geneloc/size)
      locList <- c(locList, geneloc)
      
    }
    locAll <- cbind(locAll, locList)%>% as.data.frame() %>% remove_rownames()
  }
  
  rownames(locAll) <- colnames(exprSet)
  colnames(locAll) <- names(background)
  
  locAll <- as.data.frame(locAll)
  
  model_preds <- predict(models[[5]], newdata = locAll)
  model_preds
}

res_abs <- absHBIpre(exprSet, background = genePre) %>% as.data.frame() %>% magrittr::set_colnames("fit") 

# isolate KEGG term for ssGSEA
buildGseaList <- function(species, selected) {
  require(clusterProfiler)
  require(tidyverse)
  
  R.utils::setOption("clusterProfiler.download.method",'auto')
  keggList <- download_KEGG(species = species)
  
  keggRef <- keggList$KEGGPATHID2NAME %>% filter(to %in% selected)
  
  keggPath <- keggList$KEGGPATHID2EXTID %>% filter(from %in% keggRef$from)
  
  keggPath$from <- keggRef$to[match(keggPath$from, keggRef$from)]
  
  keggPath <- split(keggPath$to, keggPath$from)
  
  keggPath
  
}

geneList <- buildGseaList(species = "hsa", 
                          selected = c("Thermogenesis",
                                       "Fatty acid biosynthesis", 
                                       "Fatty acid degradation",
                                       "PPAR signaling pathway",
                                       "AMPK signaling pathway"))

# rename gene names
IDs <- bitr(rownames(exprSet), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrezMatrix <- exprSet[IDs$SYMBOL,]
rownames(entrezMatrix) <- IDs$ENTREZID

# ssGSEA
beigeSig <- gsva(entrezMatrix, 
                 geneList,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)  %>% t() %>% scale() %>% as.data.frame()

colnames(beigeSig) <- c("AMPK", "biosynthesis", "degradation", "PPAR", "Thermogenesis")

groupEXP <- function(exprSet, spliteBy, spliteWho, myTitle){
  
  markerEXP <- exprSet[spliteBy, ]
  
  spliteWho$Marker <- ifelse(markerEXP >= mean(markerEXP), "High", "Low")
  
  spliteWho <- gather(spliteWho, process, score, -Marker) %>% 
    mutate(Marker = factor(Marker, levels = c("Low", "High"))) %>% 
    mutate(process = factor(process, levels = c("Thermogenesis", "biosynthesis", "degradation", "PPAR", "AMPK" )))
  
  ggplot(spliteWho, aes(x = process, y = score, fill = Marker))+ 
    geom_boxplot(outlier.size = 0.5)+
    ylim(-4,6)+
    scale_fill_manual(values = c("#559575", "#6163A3"))+
    theme_classic()+
    labs(x = NULL, y = "Z-scored relative indices")+
    scale_x_discrete(labels= c("Thermogenesis", "FA biosynthesis", "FA degradation", "PPAR pathway", "AMPK pathway"))+
    ggtitle(myTitle)+
    theme(axis.line.x = element_line(size = 0.3),
          axis.line.y = element_line(size = 0.3),
          axis.ticks.x = element_line(size = 0.3),
          axis.ticks.y = element_line(size = 0.3),
          plot.title = element_text(hjust = 0.5))+
    stat_compare_means(aes(label = paste0("p=", ..p.format..)), method = "t.test", label.y = 4)+ coord_flip()
}

p1 <- groupEXP(exprSet = res, spliteBy = "fit", spliteWho = beigeSig, myTitle = "HBI classification")
p2 <- groupEXP(exprSet = t(res_abs), spliteBy = "fit", spliteWho = beigeSig, myTitle = "absHBI classification")
p3 <- groupEXP(exprSet = t(res_ProFAT), spliteBy = "BAT", spliteWho = beigeSig, myTitle = "ProFAT classification")

pdf(file = "./plot/method_comparison3.pdf", width = 9, height = 4)
p1+p2+p3+plot_layout(ncol = 3, guides = "collect")
dev.off()


