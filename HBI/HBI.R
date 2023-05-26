# library
if (!require("caret", quietly = TRUE)) install.packages("caret")
if (!require("caretEnsemble", quietly = TRUE)) install.packages("caretEnsemble")
if (!require("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library(caret)
library(caretEnsemble)
library(tidyverse)

# load data
load("./hbi_models.RData")
load("./abs_models.RData")
load("./abs_background.RData")

# HBI function
HBIpre <- function(exprSet){
  meanNor <- function(x){(x-mean(x))/(max(x)-min(x))}
  test <- exprSet %>% t %>% apply(2, meanNor)  %>% as.data.frame
  res <- predict(fit, test, interval = "prediction", level = 0.95)
  res
}

# absHBI function
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
