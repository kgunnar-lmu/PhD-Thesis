options(stringsAsFactors = F)
library(ggplot2)
library(tidyverse)
library(ggbeeswarm)
library(ggrepel)
library(RColorBrewer)

data=read.csv("/Home Office/Thesis/Figures/sc meetodite tabel.csv")
head(data)

data$order2=as.numeric(paste0(data$Year,".",data$Month))

cc1 <- c(brewer.pal(name = "YlOrRd",n=9))[6]#red
cc2 <- c(brewer.pal(name = "YlOrBr",n=9))[5]#or
cc3 <- c(brewer.pal(name = "YlOrRd",n=9))[3]#yellow
cc4 <- c(brewer.pal(name = "YlGn",n=9))[6]#green
cc5 <- c(brewer.pal(name = "YlGnBu",n=9))[5]#lightblue
cc6 <- c(brewer.pal(name = "Blues",n=9))[6]#blue
cc7 <- c(brewer.pal(name = "Purples",n=9))[6]#viol
cc8 <- c(brewer.pal(name = "RdPu",n=9))[7]#purp
cc=c(cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8)


ggplot(data,aes(x=order2,y=Cells,color=Type,label=Name))+ggrepel::geom_text_repel(color="black", min.segment.length = unit(0, 'lines'))+
  geom_point(shape=1,size=3)+
  scale_y_log10(breaks=c(1,10,100,1000,10000,100000,1000000))+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  scale_x_continuous(breaks = c(2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018))+scale_color_manual(values = cc[c(4,6,1)])



