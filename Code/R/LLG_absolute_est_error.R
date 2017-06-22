rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/Dissertation/Code/R")

dataName = "desikan"

set.seed(12345)

m = 5
isSVD = 0

require(Matrix)
library(latex2exp)
source("function_collection.R")
tmpList = read_data(dataName, threshold=0, DA=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)

# sum(add(A_all)/M == 1)/(n^2)
# sum(add(A_all)/M == 0)/(n^2)

library(lattice)
new.palette=colorRampPalette(c("white","black"),space="rgb")
# new.palette=colorRampPalette(c("black","white"),space="rgb")

nColor <- 100
myAt <- seq(0, 1, length.out=nColor)
myCkey <- list(at=myAt)


add <- function(x) Reduce("+", x)
P = add(A_all)/M

mse1 <- c()
mse2 <- c()
abs1 <- c()
abs2 <- c()
# maxMSE <- 0
maxAbs <- 0
sampleBest <- c()
source("getElbows.R")
nElb = 3
dMax = ceiling(n*3/5)

sampleVec = sample.int(M, m)
A_bar = add(A_all[sampleVec])/m


####### Estimate dimensions ######
source("getElbows.R")
nElb = 3
dMax = ceiling(n*3/5)
evalVec = ase(A_bar, dMax, isSVD)[[1]]
dHat = getElbows(evalVec, n=nElb, plot=F)[[nElb]]

####### Calculate Phat ######
A.ase = ase(diag_aug(A_bar), dHat, isSVD)
if (dHat == 1) {
  Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
} else {
  Ahat <- A.ase[[3]][,1:dHat] %*% diag(A.ase[[1]][1:dHat]) %*% t(A.ase[[2]][,1:dHat])
}
P_hat = regularize(Ahat)

####### Plot the difference ######
Diff_A_bar = abs(A_bar - P)
Diff_P_hat = abs(P_hat - P)
Diff_Between = abs(A_bar - P_hat)
valLow = 0.4

nv = (Diff_A_bar<valLow)
nv[upper.tri(nv,diag=T)] = FALSE
Diff_A_bar[nv] = 0
pdf("../../Figure/Diff2_desikan_m5.pdf", family="Times", width=4, height=3.5)
levelplot(Diff_A_bar[1:n,n:1],col.regions=new.palette(nColor),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$|\\bar{A} - P|$ for Desikan with m=5')),
          at=myAt, colorkey=FALSE)
dev.off()

nv = (Diff_P_hat<valLow)
nv[upper.tri(nv,diag=T)] = FALSE
Diff_P_hat[nv] = 0
pdf("../../Figure/Diff3_desikan_m5.pdf", family="Times", width=4.5, height=3.5)
levelplot(Diff_P_hat[1:n,n:1],col.regions=new.palette(nColor),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$|\\hat{P} - P|$ for Desikan with m=5')),
          at=myAt, colorkey=FALSE)
dev.off()

nv = (Diff_Between<valLow)
nv[upper.tri(nv,diag=T)] = FALSE
Diff_Between[nv] = 0
pdf("../../Figure/Diff1_desikan_m5.pdf", family="Times", width=4, height=3.5)
levelplot(Diff_Between[1:n,n:1],col.regions=new.palette(nColor),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$|\\bar{A} - \\hat{P}|$ for Desikan with m=5')),
          at=myAt, colorkey=myCkey)
dev.off()




###### Write to files, for MATLAB use ######

Diff_between = abs(A_bar - P_hat)
Diff_A_bar = abs(A_bar - P)
Diff_P_hat = abs(P_hat - P)


valLow = sort(Diff_A_bar-Diff_P_hat, decreasing=T)[50]
nv = ((Diff_A_bar-Diff_P_hat)<valLow)
s = ""
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    if (nv[i,j]==F) {
      s = paste0(s,",",i,",",j,",",P_hat[i,j]-P[i,j])
    }
  }
}
s = substr(s,2,nchar(s))
write(s,file="../../Result/Edge_Diff_between_desikan.csv")




rowSumDiffBetween = rowSums(Diff_A_bar-Diff_P_hat)

valLow = sort(rowSumDiffBetween, decreasing=T)[5]
nv = (rowSumDiffBetween<valLow)
s = ""
for (i in 1:n) {
  if (nv[i]==F) {
    s = paste0(s,",",i)
  }
}
s = substr(s,2,nchar(s))
write(s,file="../../Result/Vertex_Diff_Between_desikan.csv")

