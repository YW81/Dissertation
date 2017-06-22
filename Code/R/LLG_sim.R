rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/Dissertation/Code/R")

lSize = 0.8

###### Fix m ######
maxIter = 10000
m = 100
nVec = c(30, 50, 100, 250, 500, 1000)
isSVD = 0
iModel = 1
d = 2

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
source("function_collection.R")

pp = list()

for (iN in 1:length(nVec)) {
  n = nVec[iN]
  
  if (isSVD) {
    fileName = paste("../../Result/result_sim_", iModel, "_d_", d,
                     "_n_", n, "_m_", m, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_sim_", iModel, "_d_", d,
                     "_n_", n, "_m_", m, "_eig.RData", sep="")
  }
  
  load(fileName)
  
  if (iN == 1) {
    mse1M = array(NaN, dim=c(length(nVec), maxIter))
    mse2M = array(NaN, dim=c(length(nVec), maxIter))
    mse3M = array(NaN, dim=c(length(nVec), maxIter))
    
    re1M = array(NaN, dim=c(length(nVec), maxIter))
    re2M = array(NaN, dim=c(length(nVec), maxIter))
    re3M = array(NaN, dim=c(length(nVec), maxIter))
  }
  
  mse1M[iN,] = error_P_hat[1,1:maxIter]*n
  mse2M[iN,] = error_P_hat[2,1:maxIter]*n
  mse3M[iN,] = error_P_hat[3,1:maxIter]*n
  
  re1M[iN,] = mse1M[iN,]/error_A_bar[1,1:maxIter]
  re2M[iN,] = mse2M[iN,]/error_A_bar[2,1:maxIter]
  re3M[iN,] = mse3M[iN,]/error_A_bar[3,1:maxIter]
  
}

mse1ExpectM = B[1,1]*(1-B[1,1])*(1/rho[1]+1/rho[1])/m
mse2ExpectM = B[2,2]*(1-B[2,2])*(1/rho[2]+1/rho[2])/m
mse3ExpectM = B[1,2]*(1-B[1,2])*(1/rho[1]+1/rho[2])/m

reExpect = (1/rho[1]+1/rho[2])

mse1MeanM = rowMeans(mse1M)
mse1LowerM = mse1MeanM - 
  sqrt(apply(mse1M, 1, var))/sqrt(dim(mse1M)[2])*1.96
mse1UpperM = mse1MeanM + 
  sqrt(apply(mse1M, 1, var))/sqrt(dim(mse1M)[2])*1.96

mse2MeanM = rowMeans(mse2M)
mse2LowerM = mse2MeanM - 
  sqrt(apply(mse2M, 1, var))/sqrt(dim(mse2M)[2])*1.96
mse2UpperM = mse2MeanM + 
  sqrt(apply(mse2M, 1, var))/sqrt(dim(mse2M)[2])*1.96

mse3MeanM = rowMeans(mse3M)
mse3LowerM = mse3MeanM - 
  sqrt(apply(mse3M, 1, var))/sqrt(dim(mse3M)[2])*1.96
mse3UpperM = mse3MeanM + 
  sqrt(apply(mse3M, 1, var))/sqrt(dim(mse3M)[2])*1.96

re1MeanM = rowMeans(re1M)
re1LowerM = re1MeanM - 
  sqrt(apply(re1M, 1, var))/sqrt(dim(re1M)[2])*1.96
re1UpperM = re1MeanM + 
  sqrt(apply(re1M, 1, var))/sqrt(dim(re1M)[2])*1.96

re2MeanM = rowMeans(re2M)
re2LowerM = re2MeanM - 
  sqrt(apply(re2M, 1, var))/sqrt(dim(re2M)[2])*1.96
re2UpperM = re2MeanM + 
  sqrt(apply(re2M, 1, var))/sqrt(dim(re2M)[2])*1.96

re3MeanM = rowMeans(re3M)
re3LowerM = re3MeanM - 
  sqrt(apply(re3M, 1, var))/sqrt(dim(re3M)[2])*1.96
re3UpperM = re3MeanM + 
  sqrt(apply(re3M, 1, var))/sqrt(dim(re3M)[2])*1.96


df <- rbind(
  data.frame(re=reExpect,lci=reExpect,uci=reExpect,
             which="Theoretical",m=m,n=max(nVec)),
  data.frame(re=reExpect,lci=reExpect,uci=reExpect,
             which="Theoretical",m=m,n=min(nVec)),
  data.frame(re=re1MeanM,lci=re1LowerM,uci=re1UpperM,
             which="B11",m=m,n=nVec),
  data.frame(re=re2MeanM,lci=re2LowerM,uci=re2UpperM,
             which="B22",m=m,n=nVec),
  data.frame(re=re3MeanM,lci=re3LowerM,uci=re3UpperM,
             which="B12",m=m,n=nVec)) %>%
  mutate(m=factor(paste0("Simulation based on SBM, M=",m),
                  c("Simulation based on SBM, M=100")))

label_y <- with(df, .75*max(re)+.25*min(re))




p <- ggplot(df,aes(x=n,y=re,linetype=factor(which),alpha=factor(which)))+
  # facet_wrap(~m)+
  #   geom_point()+
  geom_line()+
  #   geom_linerange(aes(ymin=lci,ymax=uci),alpha=.5,size=1)+
  #   geom_vline(data=dim_selection_df,
  #              aes(xintercept=value,color=which,linetyRpe=variable))+
  # scale_linetype_manual(name="",values=c("solid","longdash","dotted","dotdash"))+
  scale_linetype_manual(name="",values=c("solid","longdash","dotted","dotdash"),
                        labels=c("Theoretical", expression(B[11]), expression(B[22]), expression(B[12])))+
  # scale_shape_manual(name="",values=c(-1,-1,-1,-1),
  #                    labels=c("Theoretical", expression(B[11]), expression(B[22]), expression(B[12])))+
  scale_alpha_manual(name="",values=c(.5,1,1,1),
                     labels=c("Theoretical", expression(B[11]), expression(B[22]), expression(B[12])))+
  #   scale_color_discrete(guide=FALSE)+
  xlab("number of vertices") + ylab("scaled relative efficiency")+
  # theme(strip.text.x = element_text(size=20,face="bold"))+
  # theme(axis.text=element_text(size=20),
  # axis.title=element_text(size=20,face="bold"))+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  # theme(legend.text=element_text(size=20,face="bold"))+
  # theme(legend.key.size=unit(2,"line"))+
  theme(legend.position="bottom")
print(p)

ggsave("../../Figure/scaled_RE.pdf",
       p+theme(text=element_text(size=10,family="Times")),
       # p+theme(text=element_text(size=10,family="CM Roman")),
       width=5, height=3)


p <- ggplot(df %>% filter(which!="Theoretical"),aes(x=n,y=re/n,linetype=factor(which)))+
  geom_line()+
  scale_linetype_manual(name="",values=c("longdash","dotted","dotdash"))+
  xlab("number of vertices") + ylab("relative efficiency")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="none")
print(p)

ggsave("../../Figure/RE.pdf",
       p+theme(text=element_text(size=10,family="Times")),
       # p+theme(text=element_text(size=10,family="CM Roman")),
       width=5, height=2)


# ggplot(df,aes(x=n,y=n*(re-4),linetype=factor(which),alpha=factor(which)))+
  # geom_line()





###### Change rho ######
maxIter = 1000
m = 100
n = 500
isSVD = 0
iModel = 1
d = 2

rhoVec = (1:100)/100
nRho = length(rhoVec)

re1Expectrho = rep(0,nRho)
re2Expectrho = rep(0,nRho)
re3Expectrho = rep(0,nRho)

for (iRho in 1:nRho) {
  rho1 = rhoVec[iRho]
  rho = c(rho1, 1-rho1)
  re1Expectrho[iRho] = (1/rho[1]+1/rho[1])
  re2Expectrho[iRho] = (1/rho[2]+1/rho[2])
  re3Expectrho[iRho] = (1/rho[1]+1/rho[2])
}

df <- rbind(
  data.frame(re=re1Expectrho,which="B11",m=m,rho=rhoVec),
  data.frame(re=re3Expectrho,which="B12",m=m,rho=rhoVec))
# data.frame(re=re2Expectrho,which="B22",m=m,rho=rhoVec))

label_y <- with(df, .75*max(re)+.25*min(re))

p <- ggplot(df,aes(x=rho,y=re,linetype=factor(which)))+
  # facet_wrap(~m)+
  #   geom_point()+
  geom_line(alpha=1)+
  #   geom_linerange(aes(ymin=lci,ymax=uci),alpha=.5,size=1)+
  #   geom_vline(data=dim_selection_df,
  #              aes(xintercept=value,color=which,linetype=variable))+
  scale_linetype_manual(name="",values=c("solid","longdash","dotted"),
                        labels=c(expression(B[11]), expression(B[12]), expression(B[22])))+
  #   scale_color_discrete(guide=FALSE)+
  #   xlab("N") + ylab("Normalized Relative Efficiency")+
  scale_x_continuous(name="proportion of vertices in block one")+
  scale_y_continuous(name="scaled relative efficiency")+
  coord_trans(limy=(c(0,40)))+
  # theme(strip.text.x = element_text(size=20,face="bold"))+
  # theme(axis.text=element_text(size=20),
  # axis.title=element_text(size=20,face="bold"))+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70')) +
  # theme(legend.text=element_text(size=20,face="bold"))+
  # theme(legend.key.size=unit(2,"line"))+
  theme(legend.position=c(.82,.75))
print(p)

ggsave("../../Figure/rho.pdf",
       p+theme(text=element_text(size=10,family="Times")),
       # p+theme(text=element_text(size=10,family="CM Roman")),
       width=5, height=2.5)






