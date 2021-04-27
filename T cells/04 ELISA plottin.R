

library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(gtools)
options(stringsAsFactors = F)
####
data=read.csv("/Home Office/CSF2 Elisa project/elisa data frame.csv")
head(data)
data=data[data$Cytokine%in%c("IFN-g","IL-4","IL-17","GM-CSF"),]
data=data[data$Stimulation=="PMA+Iono",]
data7=data[data$Day=="Day-7",]
data14=data[data$Day=="Day-14",]
head(data14)
cc=c(brewer.pal(name = "Set1",n=8))

SEM <- function(x) sd(x)/sqrt(length(x))
SEM <- function(x) sd(x)
head(data14)
summarised_14=data14 %>% group_by(Condition,Cytokine)%>% 
  summarise(Value = mean(value), sem = SEM(value))
head(summarised_14)
summarised_14$Condition=factor(summarised_14$Condition,levels = mixedsort(unique(summarised_14$Condition)))

# Ifng
ifng=summarised_14[summarised_14$Cytokine=="IFN-g",]
ifngp=ggplot(ifng, aes(x=Condition, y=Value, fill=Condition)) + 
  geom_bar(stat="identity") + geom_errorbar(aes(ymin=Value-sem, ymax=Value+sem), width=.5) + 
  ggtitle(paste0("IFN-g"))+xlab("")+ylab(paste0("pg/ml"))+theme_classic()+scale_fill_manual(values = cc[c(4,1,2,3)])+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="",shape="")
ifngp
###
# IL4
il4=summarised_14[summarised_14$Cytokine=="IL-4",]
il4p=ggplot(il4, aes(x=Condition, y=Value, fill=Condition)) + 
  geom_bar(stat="identity") + geom_errorbar(aes(ymin=Value-sem, ymax=Value+sem), width=.5) + 
  ggtitle(paste0("IL-4"))+xlab("")+ylab(paste0("pg/ml"))+theme_classic()+scale_fill_manual(values = cc[c(4,1,2,3)])+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="",shape="")
il4p
plot_grid(ifngp,il4p)

library(ggpubr)

head(data)
neliteist=data[data$Day=="Day-14",]
head(neliteist)
neliteist=neliteist[neliteist$Cytokine%in%c("IFN-g","IL-4"),]
head(neliteist)
neliteist$Condition=factor(neliteist$Condition,levels = unique(neliteist$Condition))
ifng=neliteist[neliteist$Cytokine=="IFN-g",]
head(ifng)
ifng=ifng %>% group_by(Condition,Donor)%>% 
  summarise(Value = mean(value))
x1=ggbarplot(ifng, x = "Condition", y = "Value",fill = "Condition",add = c("mean_se", "jitter"))
il4=neliteist[neliteist$Cytokine=="IL-4",]
il4=il4 %>% group_by(Condition,Donor)%>% 
  summarise(Value = mean(value))
x2=ggbarplot(il4, x = "Condition", y = "Value",fill = "Condition",add = c("mean_se", "jitter"))
plot_grid(x1,x2)







cc=c(cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8)

pal <- c("blue", "black", "purple","red")


# XY plot 
dataa=read.csv("/Home Office/Data/2021.01.15.BRB_Elisas/3 plates together complete.csv")
head(dataa)
dataa$Plate=as.character(dataa$Plate)
dataa=dataa[dataa$Plate==1,]

dataa=dataa[dataa$Well%in%dds$Origin_well,]
dataa=dataa[dataa$Well%in%xyagain$well,]
prot=ggplot(dataa,aes(x=IFNG, y = IL13,color=KO,label=Well))+ggrepel::geom_text_repel(color="black")+geom_point(size=2)+
  theme_classic()+scale_color_manual(values = cc[c(6,1,3)])
prot
dataa$KO=factor(dataa$KO,levels = unique(dataa$KO)[c(3,1,2)])
head(dataa)

b1=ggbarplot(dataa, x = "KO", y = "IFNG",add = c("mean_se", "jitter"),fill = "KO")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("IFN-G")+ylab("")
b2=ggbarplot(dataa, x = "KO", y = "TNF",add = c("mean_se", "jitter"),fill = "KO")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("TNF")+ylab("")
b3=ggbarplot(dataa, x = "KO", y = "IL4",add = c("mean_se", "jitter"),fill = "KO")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("IL-4")+ylab("")
b4=ggbarplot(dataa, x = "KO", y = "IL13",add = c("mean_se", "jitter"),fill = "KO")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("IL-13")+ylab("")
plot_grid(b1,b2,b3,b4,ncol = 4)

