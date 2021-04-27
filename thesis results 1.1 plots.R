library(ggplot2)
library(gtools)
library(ggrepel)
library(ggbeeswarm)
library(RColorBrewer)
options(stringsAsFactors = F)
######################################################################################################
# Mac vs Scrb
######################################################################################################
data=read.csv("/Home Office/Thesis Data/Tecan data/2017/Used for thesis/mcplot.csv")
head(data)
plot_order=data.frame(name=unique(data$Input),order=1:6)
head(plot_order)
cc=c(brewer.pal(name = "Dark2",n=8))
#######################
ggplot(data,aes(x=Position,y=Value,color=factor(Method)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
  scale_x_continuous(label = plot_order$name,breaks = c(1,2,3,4,5,6))+scale_color_manual(values = cc[])+xlab("")+
  labs(color = "Protocol")+ggtitle("SCRB-seq VS mcSCRB-seq")

######################################################################################################
# Superscript vs maxima now! 
######################################################################################################

enz_data=read.csv("/Home Office/Thesis Data/Tecan data/2017/Used for thesis/enzymeplot.csv")
head(enz_data)
plot_order=data.frame(name=unique(enz_data$Input),order=1:5)
head(plot_order)
#######################
ggplot(enz_data,aes(x=Order,y=Value,color=factor(Enzyme)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
  scale_x_continuous(label = plot_order$name,breaks=plot_order$order)+scale_color_manual(values = cc)+ stat_smooth()+xlab("")+
  labs(color = "Protocol")+ylab("ng/ul")+ggtitle("MAXIMA H- VS SSIV RT")
######################################################################################################
ggplot(enz_data,aes(x=Order,y=Value,color=factor(Enzyme)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
  scale_x_continuous(label = plot_order$name,breaks=plot_order$order)+scale_color_manual(values = cc)+xlab("")+
  labs(color = "Protocol")+ylab("ng/ul")+ggtitle("MAXIMA H- VS SSIV RT")
######################################################################################################
ggplot(enz_data,aes(x=Order,y=Value,color=factor(Enzyme)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
  scale_x_continuous(label = plot_order$name,breaks=plot_order$order)+scale_color_manual(values = cc)+ stat_smooth()+xlab("")+
  labs(color = "Protocol")+ylab("ng/ul") + facet_grid(cols = vars(Enzyme))+ggtitle("MAXIMA H- VS SSIV RT")

######################################################################################################
# Old vs new primer 
######################################################################################################

