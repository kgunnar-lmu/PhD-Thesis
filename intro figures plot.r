library(ggplot2)

trends=read.csv("/Home Office/Thesis/Figures/pubmed publication trends vol 2.csv")
head(trends)
trends=trends[trends$Year>1980,]
ggplot(trends,aes(x=Year,y=Count,col=Type))+geom_line(size=1.5)+theme_light()


seq=read.csv("/Home Office/Thesis/Figures/seq cost.csv")
head(seq)
ggplot(seq,aes(x=Order,y=Cost,color=Type))+geom_line()+scale_x_log10()

seq2=read.csv("/Home Office/Thesis/Figures/seq cost2.csv")
head(seq2)
ggplot(seq2,aes(x=Order))+
  geom_line(aes(y=log2(MB)),color="red",size=1.5)+
  geom_line(aes(y=log2(Genome)),color="blue",size=1.5)+theme_minimal()+scale_y_continuous(breaks = log2(c(100000000,10000000,1000000,100000,10000,1000,100,10,1,0.1,0.01,0.001)))+
  scale_x_continuous(breaks = 1:71)+theme(axis.text.x = element_text(angle = 45))

