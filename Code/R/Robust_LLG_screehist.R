rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/Dissertation/Code/R")

# dataNameVec <- c("m2g", "ndmg2")
# dataNameDisplayVec <- c("m2g", "ndmg2")
dataNameVec <- c("m2g")
dataNameDisplayVec <- c("ndmg1")

source("mylibrary.R")
require(Matrix)
require(ggplot2)
eigenResult <- list()
pp_scree <- list()
pp_hist <- list()
for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  tmpList <- ReadData(dataName)
  AList <- tmpList[[1]]
  n <- tmpList[[2]]
  M <- tmpList[[3]]
  rm(tmpList)
  Abar <- add(AList)/M
  D0 <- Diagonal(n, x = rowSums(Abar)/(n-1))
  eigenResult[[iData]] <- eigen(Abar + D0)$values
  df <- data.frame(eval=eigenResult[[iData]], k=1:n)
  
  pp_scree[[iData]] <- ggplot(df,aes(x=k,y=eval))+
    geom_line()+
    scale_linetype_manual(name="",values=c("longdash","dotted","dotdash"))+
    xlab("order in algebraic") + ylab("eigenvalue")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(dataNameDisplayVec[[iData]])
  
  ggsave(paste0("../../Figure/screeplot_", dataName, ".pdf"),
         pp_scree[[iData]]+theme(text=element_text(size=10,family="Times")),
         # pp_scree[[iData]]+theme(text=element_text(size=10,family="CM Roman")),
         width=2.5, height=2)
  
  pp_hist[[iData]] <- ggplot(df, aes(x=eval))+
    geom_histogram()+
    xlab("eigenvalue") + ylab("count")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(dataNameDisplayVec[[iData]])
  
  ggsave(paste0("../../Figure/hist_", dataName, ".pdf"),
         pp_hist[[iData]]+theme(text=element_text(size=10,family="Times")),
         # pp_scree[[iData]]+theme(text=element_text(size=10,family="CM Roman")),
         width=2.5, height=2)
}

