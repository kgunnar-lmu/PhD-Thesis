library(ggplot2)
library(gtools)
library(ggrepel)
library(ggbeeswarm)
library(RColorBrewer)
options(stringsAsFactors = F)
cc=c(brewer.pal(name = "Dark2",n=8))
######################################################################################################
# Mac vs Scrb
######################################################################################################
data=read.csv("/Home Office/Thesis Data/Tecan data/2017/Used for thesis/mcplot.csv")
head(data)

# plot_order=data.frame(name=unique(data$Input),order=1:6)
# plot_order=data.frame(name=unique(data$Input),order=c(1,2,5,3,4,6))
# head(plot_order)
# data$Input=factor(data$Input,levels = unique(data$Input)[c(1,2,5,3,4,6)])
# data$Position=factor(data$Position,levels = unique(data$Position)[c(1,2,5,3,4,6)])
# cc=c(brewer.pal(name = "Dark2",n=8))
#######################
ggplot(data,aes(x=Order,y=Value,color=factor(Method)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+ 
  scale_x_continuous(label = unique(data$Input)[c(6,4,1,5,3,2)],breaks = c(1,2,3,4,5,6))+scale_color_manual(values = cc[])+xlab("")+
  labs(color = "Protocol")+ggtitle("SCRB-seq VS mcSCRB-seq")+ylab("ng/ul")+ stat_smooth()

plot4

scale_x_continuous(label = plot_order$name,breaks = c(1,2,3,4,5,6))
######################################################################################################
# Superscript vs maxima now! 
######################################################################################################

enz_data=read.csv("/Home Office/Thesis Data/Tecan data/2017/Used for thesis/enzymeplot.csv")
head(enz_data)
plot_order2=data.frame(name=unique(enz_data$Input),order=1:5)
head(plot_order2)
#######################
plot3=ggplot(enz_data,aes(x=Order,y=Value,color=factor(Enzyme)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
  scale_x_continuous(label = plot_order2$name,breaks=plot_order2$order)+scale_color_manual(values = cc)+ stat_smooth()+xlab("")+
  labs(color = "RT Enzyme")+ylab("ng/ul")+ggtitle("MAXIMA H- VS SSIV RT")
######################################################################################################
# ggplot(enz_data,aes(x=Order,y=Value,color=factor(Enzyme)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
#   scale_x_continuous(label = plot_order2$name,breaks=plot_order2$order)+scale_color_manual(values = cc)+xlab("")+
#   labs(color = "RT Enzyme")+ylab("ng/ul")+ggtitle("MAXIMA H- VS SSIV RT")
######################################################################################################
# ggplot(enz_data,aes(x=Order,y=Value,color=factor(Enzyme)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
#   scale_x_continuous(label = plot_order2$name,breaks=plot_order2$order)+scale_color_manual(values = cc)+ stat_smooth()+xlab("")+
#   labs(color = "RT Enzyme")+ylab("ng/ul") + facet_grid(cols = vars(Enzyme))+ggtitle("MAXIMA H- VS SSIV RT")

######################################################################################################
# Old vs new primer 
######################################################################################################

prim_data=read.csv("/Home Office/Thesis Data/Tecan data/2017/Used for thesis/primerplot.csv")
head(prim_data)
prim_data$Primer[prim_data$Primer=="New"]="New primer"
plot_order3=data.frame(name=unique(prim_data$Input),order=1:5)
head(plot_order3)

#######################
plot2=ggplot(prim_data,aes(x=Order,y=Value,color=factor(Primer)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
  scale_x_continuous(label = plot_order3$name,breaks=plot_order3$order)+scale_color_manual(values = cc)+ stat_smooth()+xlab("")+
  labs(color = "Primer")+ylab("ng/ul")+ggtitle("Old VS New Oligo(dT) primer")
######################################################################################################
# Stat tests
head(prim_data)
atw=aov(Value~Primer + Input, data = prim_data)
summary(atw)
prim_data
prim_data$In=c(0,0,1,1,1,1,1,10,10,10,10,10,100,100,100,1000,1000,1000,0,0,1,1,1,1,1,10,10,10,10,10,100,100,100,1000,1000,1000)
res.man=manova(cbind(Value, In) ~ Primer, data = prim_data)
summary(res.man)
# summary(manova(cbind(prim_data$Value, prim_data$In)~prim_data$Primer))
# summary(aov(cbind(prim_data$Value, prim_data$In)~prim_data$Primer))
# ?aov()

######################################################################################################
# ggplot(prim_data,aes(x=Order,y=Value,color=factor(Primer)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
#   scale_x_continuous(label = plot_order3$name,breaks=plot_order3$order)+scale_color_manual(values = cc)+xlab("")+
#   labs(color = "Primer")+ylab("ng/ul")+ggtitle("Old VS New Oligo(dT) primer")
######################################################################################################
# ggplot(prim_data,aes(x=Order,y=Value,color=factor(Primer)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
#   scale_x_continuous(label = plot_order3$name,breaks=plot_order3$order)+scale_color_manual(values = cc)+ stat_smooth()+xlab("")+
#   labs(color = "Primer")+ylab("ng/ul")+ggtitle("Old VS New Oligo(dT) primer")+facet_grid(cols=vars(Primer))

######################################################################################################
# Old vs new primer 
######################################################################################################

comp_data=read.csv("/Home Office/Thesis Data/Tecan data/2020/Used for thesis/old new mac.csv")
head(comp_data)
plot_order4=data.frame(name=unique(comp_data$Input),order=c(2,1,3:5))
head(plot_order4)
comp_data$Value[47]=mean(comp_data$Value[c(45,46,48)])
comp_data$Value[47]=NA
comp_data=comp_data[comp_data$Method!="McSCRB-seq",]
#######################
plot1=ggplot(comp_data,aes(x=Order,y=Value,color=factor(Method)))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
  scale_x_continuous(label = plot_order4$name,breaks=plot_order4$order)+scale_color_manual(values = cc)+xlab("")+
  labs(color = "Method")+ylab("ng/ul")+ggtitle("SC Methods Comparisson")+ stat_smooth(method = "loess")
plot1

##########################################################################################################################
library(cowplot)
plot_grid(plot1,plot2,plot3,plot4)
##########################################################################################################################
##########################################################################################################################

# Mapping percentage
##########################################################################################################################
tbrb=read.csv("/Home Office/Thesis/Figures/thomas brb-seq mapping percentage.csv")
tbrb$Length=factor(tbrb$Length,levels = unique(tbrb$Length))
head(tbrb)
cc=c(brewer.pal(name = "Set1",n=8))

ggplot(tbrb, aes(x=Length, y=Percent,fill="")) + geom_bar(stat = "identity", position = "dodge", width=0.9,alpha=0.8)+
  geom_errorbar(data = tbrb,aes(ymin = Percent-SD, ymax = Percent+SD),width=0.9,colour="black", alpha=1, size=0.3,position = position_dodge())+
  ggtitle("Uniquely Mapped sequences")+xlab("")+ylab('Percent mapped')+theme_classic()+scale_fill_manual(values=cc[2]) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="")

# Cell input 
##########################################################################################################################
input=read.csv("/Home Office/Thesis/Figures/brb testing.csv")
head(input)
# plot_order5=data.frame(name=unique(comp_data$Input),order=c(2,1,3:5))
# head(plot_order4)
# input$Value[47]=mean(comp_data$Value[c(45,46,48)])
# comp_data$Value[47]=NA
# comp_data=comp_data[comp_data$Method!="McSCRB-seq",]
#######################
ggplot(input,aes(x=order,y=value,color="moo"))+geom_beeswarm(dodge.width=0.5,size=2,cex=2)+theme_classic()+
  scale_x_continuous(label = input$type,breaks=input$order)+
  scale_color_manual(values = cc[3])+xlab("")+labs(color = "Input # of Cells")+ylab("ng/ul")+ggtitle("")+ stat_smooth(method = "loess")

