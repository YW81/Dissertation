rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/Dissertation/Code/R")

dataName = "desikan"

source("function_collection.R")
require(ggplot2)

tmpList = read_data(dataName, DA=F, newGraph=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)
P = add(A_all)/M
diag(P) = NA
pVec = c(P)
pVec = pVec[!is.na(pVec)]
pVec
# hist(pVec)

(sum(pVec == 0))/n/(n-1)
(sum(pVec == 1))/n/(n-1)

df <- data.frame(p=pVec)
# label_y <- with(df, .75*yMax+.25*yMin)

gg <- ggplot(df,aes(x=p))+
  geom_histogram()+
  xlab("P") + ylab("count")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="none")+
  ggtitle("Histogram of P for Desikan")

ggsave("../../Figure/P_hist_desikan.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=4.5,height=3.5)
