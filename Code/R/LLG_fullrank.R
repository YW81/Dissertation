rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/Dissertation/Code/R")

lSize = .8
legendSize = 1.5

# dataName = "CPAC200"
dataName = "desikan"
# dataName = "JHU"

mVec = c(1, 5, 10)

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

fileName = paste("../../Result/result_", dataName, "_brute_fullrank1_",
                 "_eig_dim.RData", sep="")
load(fileName)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  
  fileName = paste("../../Result/result_", dataName, "_brute_fullrank1_",
                   "m_", m, "_eig.RData", sep="")
  load(fileName)
  
  if (iM == 1) {
    errorAbarMean = array(0, dim=c(length(mVec), 1))
    errorAbarLower = array(0, dim=c(length(mVec), 1))
    errorAbarUpper = array(0, dim=c(length(mVec), 1))
    errorPhatEIGMean = array(0, dim=c(length(mVec), n))
    errorPhatEIGLower = array(0, dim=c(length(mVec), n))
    errorPhatEIGUpper = array(0, dim=c(length(mVec), n))
  }
  
  error_A_bar = error_A_bar^2*n*(n-1)
  error_P_hat = error_P_hat^2*n*(n-1)
  
  errorAbarMean[iM] = mean(error_A_bar)
  errorAbarLower[iM] = errorAbarMean[iM] - sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
  errorAbarUpper[iM] = errorAbarMean[iM] + sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
  errorPhatEIGMean[iM,] = rowMeans(error_P_hat)
  errorPhatEIGLower[iM,] = errorPhatEIGMean[iM,] - 
    sqrt(apply(error_P_hat, 1, var))/sqrt(dim(error_P_hat)[2])*1.96
  errorPhatEIGUpper[iM,] = errorPhatEIGMean[iM,] + 
    sqrt(apply(error_P_hat, 1, var))/sqrt(dim(error_P_hat)[2])*1.96
}






errorPhatUSVT = rep(0, length(mVec))
errorPhatZG = rep(0, length(mVec))
for (iM in 1:length(mVec)) {
  x = dUSVTMean[iM]
  x1 = floor(x)
  y1 = errorPhatEIGMean[iM, x1]
  x2 = ceiling(x)
  y2 = errorPhatEIGMean[iM, x2]
  errorPhatUSVT[iM] = (y2-y1)/(x2-x1)*(x-x1)+y1
  
  x = dZGMean[iM]
  x1 = floor(x)
  y1 = errorPhatEIGMean[iM, x1]
  x2 = ceiling(x)
  y2 = errorPhatEIGMean[iM, x2]
  errorPhatZG[iM] = (y2-y1)/(x2-x1)*(x-x1)+y1
}




error_by_dim_df <- rbind(
  data.frame(mse=errorAbarMean,lci=errorAbarLower,uci=errorAbarUpper,
             which="Abar",m=mVec,d=n),
  data.frame(mse=errorAbarMean,lci=errorAbarLower,uci=errorAbarUpper,
             which="Abar",m=mVec,d=2),
  data.frame(mse=c(errorPhatEIGMean),lci=c(errorPhatEIGLower),uci=c(errorPhatEIGUpper),
             which="Phat",m=rep(mVec,n),d=rep(1:n,each=3))) %>%
  mutate(m=factor(paste0("m=",m),c("m=1","m=5","m=10")))

dim_selection_df <- rbind(
  data.frame(mse=errorPhatZG,lci=errorPhatZG,uci=errorPhatZG,
             which="ZG 3rd",m=mVec,d=dZGMean),
  data.frame(mse=errorPhatUSVT,lci=errorPhatUSVT,uci=errorPhatUSVT,
             which="USVT c=0.7",m=mVec,d=dUSVTMean)) %>%
  mutate(m=factor(paste0("m=",m),c("m=1","m=5","m=10")))

label_y <- with(error_by_dim_df, .75*max(mse)+.25*min(mse))

gg <- ggplot(error_by_dim_df,aes(x=d,y=mse,linetype=factor(which),shape=factor(which)))+
  facet_wrap(~m, nrow=1)+
  geom_point(data=dim_selection_df,size=1.5)+
  scale_linetype_manual(name="",values=c(1,2,0,0),
                        labels=c(expression(bar(A)), expression(hat(P)), "USVT c=0.7", "ZG 3rd"))+
  scale_shape_manual(name="",values=c(-1,-1,15,17),
                     labels=c(expression(bar(A)), expression(hat(P)), "USVT c=0.7", "ZG 3rd"))+
  geom_line()+
  geom_linerange(aes(ymin=lci,ymax=uci),linetype=1,alpha=.5,size=.5)+
  xlab("dimension")+ylab("MSE")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="bottom")+
  ggtitle(paste0(dataName, ", n=", n))


gg <- ggplot(error_by_dim_df[error_by_dim_df$d>1,],
             aes(x=d,y=mse,linetype=factor(which),shape=factor(which)))+
  facet_wrap(~m)+
  geom_point(data=dim_selection_df,size=1.5)+
  scale_linetype_manual(name="",values=c(1,2,0,0),
                        labels=c(expression(bar(A)), expression(hat(P)), "USVT c=0.7", "ZG 3rd"))+
  scale_shape_manual(name="",values=c(-1,-1,15,17),
                     labels=c(expression(bar(A)), expression(hat(P)), "USVT c=0.7", "ZG 3rd"))+
  geom_line()+
  # geom_linerange(aes(ymin=lci,ymax=uci),alpha=.5,size=1)+
  xlab("dimension")+ylab("MSE")+
  # theme(strip.text.x = element_text(size=20,face="bold"))+
  # theme(axis.text=element_text(size=15),
  #       axis.title=element_text(size=20,face="bold"))+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  # theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
# ggtitle(paste0("Simulation based on ", dataName, ", N=", n, ", M=", m))+
# theme(legend.key.size=unit(legendSize,"line"))+
# theme(plot.title=element_text(lineheight=.8,size=20,face="bold"))
print(gg)



ggsave("../../Figure/sim_desikan.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=5.5,height=2.5)


