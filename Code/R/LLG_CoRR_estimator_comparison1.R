rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/Dissertation/Code/R")
dataName = "desikan"

set.seed(12345)

m = 1
isSVD = 0

require(Matrix)
library(latex2exp)
source("function_collection.R")
tmpList = read_data(dataName, threshold=0, DA=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)

library(lattice)
new.palette=colorRampPalette(c("white","black"),space="rgb")

nColor <- 100
myAt <- seq(0, 1, length.out=nColor)
myCkey <- list(at=myAt)


add <- function(x) Reduce("+", x)
P = add(A_all)/M


pdf("../../Figure/P_desikan.pdf", family="Times", width=4, height=3.5)
levelplot(P[1:n,n:1],col.regions=new.palette(nColor),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$P$ for Desikan')),
          at=myAt, colorkey=F, lwd=0)
dev.off()

mse1 <- c()
mse2 <- c()
abs1 <- c()
abs2 <- c()
maxAbs <- 0
sampleBest <- c()
source("getElbows.R")
nElb = 3
dMax = ceiling(n*3/5)

sampleVec = sample.int(M, m)
A_bar = add(A_all[sampleVec])/m

valLow = 0.4
Diff_A_bar = abs(A_bar - P)
nv1 = (Diff_A_bar < valLow)
sum((Diff_A_bar >= valLow))/2
nv2 = lower.tri(A_bar, diag = T)
A_bar_combine = A_bar
A_bar_combine[nv1 & nv2] = 0

pdf(paste0("../../Figure/Abar_desikan_m", m, ".pdf"), family="Times", width=4, height=3.5)
levelplot(A_bar_combine[1:n,n:1],col.regions=new.palette(nColor),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$\\bar{A}$ for Desikan with m=1')),
          at=myAt, colorkey=F, lwd=0)
dev.off()

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

Diff_P_hat = abs(P_hat - P)
nv1 = (Diff_P_hat < valLow)
sum((Diff_P_hat >= valLow))/2
nv2 = lower.tri(P_hat, diag = T)
P_hat_combine = P_hat
P_hat_combine[nv1 & nv2] = 0


pdf(paste0("../../Figure/Phat_desikan_m", m, ".pdf"), family="Times", width=4.5, height=3.5)
levelplot(P_hat_combine[1:n,n:1],col.regions=new.palette(nColor),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX('$\\hat{P}$ for Desikan with m=1')),
          at=myAt, colorkey=myCkey, lwd=0)
dev.off()
print(dHat)

MSE_Abar = norm(A_bar - P, "F")^2/n/(n-1)
MSE_Phat = norm(P_hat - P, "F")^2/n/(n-1)
(MSE_Abar - MSE_Phat)/MSE_Abar

