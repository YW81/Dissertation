rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/Dissertation/Code/R")
library(extrafont)
set.seed(12345)

require(Matrix)
require(igraph)
library(lattice)
library(latex2exp)
new.palette=colorRampPalette(c("white","black"),space="rgb")

source("function_collection.R")

B = array(c(0.9,0.27,0.05,0.1,0.3,
            0.27,0.67,0.02,0.26,0.14,
            0.05,0.02,0.44,0.25,0.33,
            0.1,0.26,0.25,0.7,0.18,
            0.3,0.14,0.33,0.18,0.58), dim=c(5,5))

isSVD <- 0


rho = c(0.22,0.39,0.05,0.16,0.18)
K <- length(rho)
n = 200

tau = rep(1:5,round(n*rho))
P = B[tau,tau]
diag(P) = 0

myAt <- seq(0, 1, length.out=20)
myCkey <- list(at=myAt)

pdf("../../Figure/SBM_P.pdf", family="Times", width=4, height=4.4)
# levelplot(P[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
#           ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
#           main=list(label="Probability matrix"),
#           at=myAt, colorkey=F, lwd=0)
levelplot(P[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX(sprintf('Probability matrix $P$'))),
          at=myAt, colorkey=F, lwd=0)
dev.off()


g = sample_sbm(n, B, round(n*rho), directed=F, loops=F)
A = as_adj(g, type="both", sparse=FALSE)

pdf("../../Figure/SBM_A.pdf", family="Times", width=4.53, height=4.4)
# levelplot(as.matrix(A[1:n,n:1]),col.regions=new.palette(20),xlab=list(cex=0),
#           ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
#           main=list(label="Adjacency matrix"),
#           at=myAt, colorkey=myCkey, lwd=0)
levelplot(as.matrix(A[1:n,n:1]),col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX(sprintf('Adjacency matrix $A$'))),
          at=myAt, colorkey=myCkey, lwd=0)
dev.off()

A2 = sample_sbm(n, B, round(n*rho), directed=F, loops=F)[]
A3 = sample_sbm(n, B, round(n*rho), directed=F, loops=F)[]

Abar <- (A+A2+A3)/3


####### Estimate dimensions ######
source("getElbows.R")
nElb = 2
dMax = ceiling(n*3/5)
evalVec = ase(Abar, dMax, isSVD)[[1]]
dHat = getElbows(evalVec, n=nElb, plot=F)[[nElb]]

dHat = 5;

####### Calculate Phat ######
A.ase = ase(diag_aug(Abar), dHat, isSVD)
if (dHat == 1) {
  Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
} else {
  Ahat <- A.ase[[3]][,1:dHat] %*% diag(A.ase[[1]][1:dHat]) %*% t(A.ase[[2]][,1:dHat])
}
P_hat = regularize(Ahat)


Phat <- as.matrix(with(eigen(Abar), vectors[,1:K]%*%diag(values[1:K])%*%t(vectors[,1:K])))
Phat[Phat<0]<- 0
Phat[Phat>1]<- 1

pdf("../../Figure/SBM_Abar.pdf", family="Times", width=4.53, height=4.4)
# levelplot(as.matrix(Abar[1:n,n:1]),col.regions=new.palette(20),xlab=list(cex=0),
#           ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
#           main=list(label="Element-wise mean (m=3)"),
#           at=myAt, colorkey=myCkey, lwd=0)
levelplot(as.matrix(Abar[1:n,n:1]),col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX(sprintf('Element-wise mean $\\bar{A}$ ($m=3$)'))),
          at=myAt, colorkey=myCkey, lwd=0)
dev.off()


pdf("../../Figure/SBM_Phat.pdf", family="Times", width=4, height=4.4)
# levelplot(P_hat[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
#           ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
#           main=list(label="Rank-5 approximation (m=3)"),
#           at=myAt, colorkey=F, lwd=0)
levelplot(P_hat[1:n,n:1],col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=TeX(sprintf('Rank-5 approximation $\\hat{P}$ ($m=3$)'))),
          at=myAt, colorkey=F, lwd=0)
dev.off()
