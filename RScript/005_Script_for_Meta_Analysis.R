# library
library(meta)
library(patchwork)
library(tidyverse)

# global options
rm(list = ls())
options(stringsAsFactors = FALSE)

# load data
loads <- function(...){
  file <- list.files(...)
  
  for (i in 1:length(file)) {
    load(file[i], envir = parent.frame(1))
  }
}

loads(path = "./RData/", pattern = "result_B.*", full.names = T)

# human brown-like adipocyte markers
feaList <- c("OASL", "UCP1", "IL1B", "FAR2", "SORL1", "CD96", "STX11", "CA5B", "REEP6", "DHRS11", "BANK1")
traList <- c("ADIPOQ", "CA4", "CD36", "CIDEA", "CITED1", "DIO2", "ELOVL3", "EPSTI1", "FABP4", "HOXC9", "IRF4", "KCNK3", "LHX8", "MTUS1", "MYF5", "PDK4", "PLIN1", "PPARG", "PPARGC1A", "PRDM16", "SHOX2", "SLC2A4", "TMEM26", "TNFRSF9", "UCP1", "ZIC1")
geneList <- c(feaList, traList) %>% unique()

cohortList <- list(res_BC2_GSE57896_A,
                   res_BC2_GSE57896_B,
                   res_BC3_GSE71293_GW7647,
                   res_BC3_GSE71293_Rosi,
                   res_BC5_GSE115421,
                   res_BT1_GSE113764, 
                   res_BT2_GSE122721,
                   res_BT3_GSE13070_IR,
                   res_BT3_GSE13070_IS)

auto_meta <- function(cohortList, geneList){
  require(meta)
  require(tidyverse)
  
  # merge phenoData
  phenoList <- lapply(cohortList, function(x) x$phenoData %>% add_column(cohort = x$gse))
  phenoList <- lapply(phenoList, function(x) x[,c("sample", "group", "cohort")])
  phenoData <- Reduce(rbind, phenoList)
  phenoData$group <- gsub("Deep neck progenitors", "Beige", phenoData$group)
  phenoData$group <- gsub("Subcutaneous progenitors", "White", phenoData$group)
  phenoData$group <- gsub("BAT", "Beige", phenoData$group)
  phenoData$group <- gsub("WAT", "White", phenoData$group)
  
  # merge gene expression data
  expList <- lapply(cohortList, function(x) x$expMatrix) %>% 
    lapply(function(x) x[rownames(x) %in% geneList,]) %>%
    lapply(function(x) t(x) %>% scale() %>% as.data.frame())
  
  expMatrix <- dplyr::bind_rows(expList)
  
  resList <- c()
  for (i in 1:ncol(expMatrix)) {
    dat <- phenoData %>% 
      add_column(gene = expMatrix[,i]) %>% 
      na.omit() %>%
      group_by(cohort, group) %>%
      summarise(mean = mean(gene),
                sd = sd(gene),
                n = n())
    
    dat_e <- dat[dat$group == "Beige",]
    colnames(dat_e)[2:ncol(dat_e)] <- paste0(colnames(dat_e)[2:ncol(dat_e)], ".e")
    dat_c <- dat[dat$group == "White",]
    colnames(dat_c)[2:ncol(dat_c)] <- paste0(colnames(dat_c)[2:ncol(dat_c)], ".c")
    dat_s <- merge(dat_e,dat_c, by = "cohort") %>% na.omit()
    
    res <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
                    data = dat_s,
                    studlab = cohort,
                    sm = "SMD",
                    fixed = F, random = T, overall = T)
    
    resList <- rbind(resList, c(res$TE.random, res$pval.random, nrow(dat_s)))
    
  }
  rownames(resList) <- colnames(expMatrix)
  colnames(resList) <- c("SMD", "pvalue", "n")
  resList
}

res <- auto_meta(cohortList, geneList)

res %>% as.data.frame() %>% rownames_to_column("symbol") %>% arrange(pvalue) %>% write.csv(file = "./RData/result_meta_gene.csv", row.names = F, quote = F)

data <- res %>% 
  as.data.frame() %>%
  subset(.$SMD > 0) %>%  
  mutate(color = ifelse(pvalue < 0.05, "significant", "non-significant")) %>%
  mutate(pvalue = -log10(pvalue)) %>% 
  arrange(pvalue) %>% 
  rownames_to_column("symbol") %>%
  mutate(symbol = factor(symbol, levels = symbol))



# plot
p <- ggplot(data, aes(x = symbol, y = pvalue, color = color, fill = color)) +
  geom_segment(aes(x = symbol, xend = symbol, y = 0, yend = pvalue),size = 2) +
  geom_point(size = 4, shape = 21, stroke = 2) +
  scale_fill_manual(values = c("#2376b7", "#f0a1a8"))+
  scale_color_manual(values = c("#2376b766", "#f0a1a866"))+
  labs(x = NULL, y = "-log10(Overall p value)")+
  geom_hline(yintercept = -log10(0.05))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35,vjust = 0.5,hjust = 0.5),
        legend.title=element_blank(),
        legend.justification = c(-0.2,1.5), 
        legend.position =  c(0, 1), 
        legend.direction = 'horizontal',
        legend.background = element_rect(color = "black", size = 0.5, linetype="solid"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5))


pdf(file = "./plot/pop_meta.pdf", width = 9, height = 3)
p
dev.off()

# plot the meta forest
datList <- lapply(cohortList, function(x) x$expMatrix %>% t %>% scale %>% as.data.frame) %>% Reduce(dplyr::bind_rows,.)

phenoList <- lapply(cohortList, function(x) x$phenoData %>% add_column(cohort = x$gse))
phenoList <- lapply(phenoList, function(x) x[,c("sample", "group", "cohort")])
phenoList <- Reduce(rbind, phenoList)
phenoList$group <- gsub("Deep neck progenitors", "Beige", phenoList$group)
phenoList$group <- gsub("Subcutaneous progenitors", "White", phenoList$group)
phenoList$group <- gsub("BAT", "Beige", phenoList$group)
phenoList$group <- gsub("WAT", "White", phenoList$group)

phenoList$SORL1 <- datList[,"SORL1"]
phenoList$DHRS11 <- datList[,"DHRS11"]
phenoList$CD96 <- datList[,"CD96"]
phenoList$FAR2 <- datList[,"FAR2"]
phenoList$REEP6 <- datList[,"REEP6"]
phenoList$OASL <- datList[,"OASL"]
phenoList$UCP1 <- datList[,"UCP1"]
phenoList$CA4 <- datList[,"CA4"]
phenoList$BANK1 <- datList[,"BANK1"]
phenoList$IL1B <- datList[,"IL1B"]
phenoList$STX11 <- datList[,"STX11"]
phenoList$PRDM16 <- datList[,"PRDM16"]
phenoList$ELOVL3 <- datList[,"ELOVL3"]
phenoList$KCNK3 <- datList[,"KCNK3"]

plotMeta <- function(phenoList, studlab, xlim){
  require(meta)
  require(tidyverse)
  
  colnames(phenoList)[colnames(phenoList) == studlab] <- "gene"
  
  dat <- phenoList %>%
    group_by(cohort, group) %>%
    summarise(mean = mean(gene),
              sd = sd(gene),
              n = n()) %>% na.omit()
  
  dat_e <- dat[dat$group == "Beige",]
  colnames(dat_e)[2:ncol(dat_e)] <- paste0(colnames(dat_e)[2:ncol(dat_e)], ".e")
  dat_c <- dat[dat$group == "White",]
  colnames(dat_c)[2:ncol(dat_c)] <- paste0(colnames(dat_c)[2:ncol(dat_c)], ".c")
  dat_s <- merge(dat_e,dat_c, by = "cohort") %>% na.omit()
  
  res <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
                  data = dat_s,
                  studlab = cohort,
                  sm = "SMD",
                  fixed = F, random = T, overall = T)
  
  pdf(file = paste0("./plot/meta_", studlab, ".pdf"), width = 7, height = 4)
  forest(res,
         test.overall = T,layout = "meta",leftcols = "studlab",
         print.I2 = F, print.tau2 = F, print.pval.Q = F,xlim = xlim,
         smlab = "White     Brown",
         label.test.overall.random = "Overall  ",
         fontsize = 11,colgap.forest.left = "0.5cm",
         leftlabs = paste0("Index ", "(", studlab, ")"))
  dev.off()
}


plotMeta(phenoList, studlab = "SORL1", xlim = c(-6,6))
plotMeta(phenoList, studlab = "DHRS11", xlim = c(-6,6))
plotMeta(phenoList, studlab = "CD96", xlim = c(-6,6))
plotMeta(phenoList, studlab = "FAR2", xlim = c(-6,6))
plotMeta(phenoList, studlab = "REEP6", xlim = c(-6,6))
plotMeta(phenoList, studlab = "OASL", xlim = c(-6,6))
plotMeta(phenoList, studlab = "UCP1", xlim = c(-6,6))
plotMeta(phenoList, studlab = "CA4", xlim = c(-6,6))
plotMeta(phenoList, studlab = "BANK1", xlim = c(-6,6))
plotMeta(phenoList, studlab = "IL1B", xlim = c(-6,6))
plotMeta(phenoList, studlab = "STX11", xlim = c(-6,6))
plotMeta(phenoList, studlab = "PRDM16", xlim = c(-6,6))
plotMeta(phenoList, studlab = "ELOVL3", xlim = c(-6,6))
plotMeta(phenoList, studlab = "KCNK3", xlim = c(-6,6))


# plot gene expression
phenoList$group <- factor(phenoList$group, levels = c("White", "Beige"))

plotMetaEXp <- function(phenoList, gene){
  p <- ggplot(phenoList, aes(x = cohort, y = phenoList[,gene], fill = group))+
    geom_boxplot()+
    theme_classic()+
    labs(x = NULL, y = gene,fill = "Group")+
    scale_fill_manual(values = c("#2376b7", "#f0a1a8"),labels=c("White", "Brown-like"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
}

p1 <- plotMetaEXp(phenoList, gene = "SORL1")
p2 <- plotMetaEXp(phenoList, gene = "DHRS11")
p3 <- plotMetaEXp(phenoList, gene = "CD96")
p4 <- plotMetaEXp(phenoList, gene = "FAR2")
p5 <- plotMetaEXp(phenoList, gene = "REEP6")
p6 <- plotMetaEXp(phenoList, gene = "OASL")
p7 <- plotMetaEXp(phenoList, gene = "UCP1")
p8 <- plotMetaEXp(phenoList, gene = "CA4")
p9 <- plotMetaEXp(phenoList, gene = "BANK1")
p10 <- plotMetaEXp(phenoList, gene = "IL1B")
p11 <- plotMetaEXp(phenoList, gene = "STX11")
p12 <- plotMetaEXp(phenoList, gene = "PRDM16")
p13 <- plotMetaEXp(phenoList, gene = "ELOVL3")
p14 <- plotMetaEXp(phenoList, gene = "KCNK3")

pdf(file = "./plot/meta_gene_exp.pdf", width = 15, height = 15)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+plot_layout(ncol = 3, guides = "collect")
dev.off()








