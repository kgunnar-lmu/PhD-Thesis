# # Purpose of the script:
# Seurat Analysis of iPSCs LPS stim
# Load libraries
options(stringsAsFactors = F)
library(biomaRt)
library(cowplot)
library(data.table)
library(dplyr)
library(GGally)
library(ggplot2)
library(ggsci)
library(gtools)
library(Hmisc)
library(knitr)
library(MAST)
library(Matrix)
library(mclust)
library(NMF)
library(pheatmap)
library(plotly)
library(RColorBrewer)
library(reshape2)
library(rsvd)
library(Rtsne)
library(scater)
library(scran)
library(sctransform)
library(Seurat)
library(stringr)
library(tidyverse)

# Load Data
Thomas_brb_dge=readRDS("/Home Office/Data/2018_06_01._Thomas-brb/L100_D_D/R data files mine/Thomas_DGE_filtered.rds")
dim(Thomas_brb_dge)
cell_annot=readRDS("/Home Office/Data/2018_06_01._Thomas-brb/L100_D_D/R data files mine/Thomas_BRB_annot_Filtered.rds")
dim(cell_annot)
########################
#########################
for (i in seq_along(cell_annot$Type)) {
  cell_annot$KO[i]=strsplit(cell_annot$Type,"[.]")[[i]][1]}

mixedsort(unique(cell_annot$KO))
# cell_annot=cell_annot[cell_annot$KO%in%c("Rela","WT","IKKa","IKKb","IKKa_IKKb"),]
# cell_annot=cell_annot[cell_annot$Stim.%in%c("LPS","Unstim."),]
# head(cell_annot)

cell_annot$KOness=ifelse(cell_annot$KO=="WT","WT","KO")

Thomas_brb_dge=Thomas_brb_dge[,cell_annot$Full_BC]
dim(Thomas_brb_dge)
dim(cell_annot)
cc=c(brewer.pal(name = "Set1",n=9))
# Rename ens to sym
Gene_data=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv",col.names =c("gene_id","gene_name","description") )
Gene_data=Gene_data[Gene_data$gene_id%in%rownames(Thomas_brb_dge),]
paste("You have ", sum(duplicated(Gene_data$gene_name)),"duplicated gene names in your data")
Gene_data=Gene_data[!duplicated(Gene_data$gene_name),]
Thomas_brb_dge=Thomas_brb_dge[rownames(Thomas_brb_dge)%in%Gene_data$gene_id,]
paste("You have ", sum(duplicated(Gene_data$gene_name)),"duplicated gene names in your data")
rownames(Thomas_brb_dge)=Gene_data[match(rownames(Thomas_brb_dge),Gene_data$gene_id),"gene_name"]

# Make Seurat
Thomas_dots=CreateSeuratObject(counts = Thomas_brb_dge, project = "BRB")
Thomas_dots@meta.data<- cbind(Thomas_dots@meta.data,cell_annot[match(colnames(Thomas_brb_dge), 
                                                                     cell_annot$Full_BC),])

Thomas_dots[["percent.mt"]]=PercentageFeatureSet(Thomas_dots, pattern = "^MT-")
VlnPlot(Thomas_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Thomas_dots <- subset(Thomas_dots, subset = percent.mt <  10 & percent.mt >  5 & nFeature_RNA > 11000)

VlnPlot(Thomas_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dim(Thomas_dots)
##########################################################################################################
"Transform and Normalizing"
##########################################################################################################
Thomas_dots <- NormalizeData(Thomas_dots, normalization.method = "LogNormalize",)
Thomas_dots <- FindVariableFeatures(Thomas_dots, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Thomas_dots), 10)
plot1 <- VariableFeaturePlot(Thomas_dots)+  theme(legend.position="none")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)+  theme(legend.position="none")
plot2

all.genes <- rownames(Thomas_dots)
head(Thomas_dots@meta.data)
Thomas_dots <- ScaleData(Thomas_dots, features = all.genes)
# Thomas_dots <- ScaleData(Thomas_dots, features = all.genes,vars.to.regress = "UMIs")
# Thomas_dots <- ScaleData(Thomas_dots, vars.to.regress = "Set") # regressing out plates messes up the thing
# Thomas_dots=SCTransform(Thomas_dots,vars.to.regress = "Set")
# Thomas_dots=SCTransform(Thomas_dots,vars.to.regress = "UMIs")
Thomas_dots <- RunPCA(Thomas_dots, features = VariableFeatures(object = Thomas_dots),npcs = 50)
DimPlot(Thomas_dots, reduction = "pca",group.by = "Stim.",pt.size = 2)
DimPlot(Thomas_dots, reduction = "pca",group.by = "KOness",pt.size = 2)

############## Finding clusters:
Thomas_dots <- FindNeighbors(Thomas_dots, dims = 1:2)
Thomas_dots <- FindClusters(Thomas_dots, resolution = 0.1)

pca = Thomas_dots@reductions$pca
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
varExplained=round(varExplained*100,digits = 2)

pca_data=Thomas_dots@reductions$pca@cell.embeddings[,1:4]
pca_data=cbind(pca_data,Thomas_dots@meta.data)
head(pca_data)

ggplot(pca_data, aes(x = PC_1, y = PC_2, col = Stim.,label=KO)) + geom_point(size = 2,alpha=0.9) + xlab(paste("PC1 (", varExplained[1], "%)"))+theme_classic()+
  ylab(paste("PC2 (",varExplained[2], "%)"))+ ggtitle("PCA")+scale_colour_manual(values = cc,name="Treatment")

ggplot(pca_data, aes(x = PC_1, y = PC_2, col = seurat_clusters ,label=KO)) + geom_point(size = 2,alpha=0.9) + xlab(paste("PC1 (", varExplained[1], "%)"))+theme_classic()+
  ylab(paste("PC2 (",varExplained[2], "%)"))+ ggtitle("PCA")+scale_colour_manual(values = cc[c(3,1,2)],name="Clusters")

pca_data$nimed[pca_data$KO%in%c("WT","Rela","IRF3_7","IKKa_IKKb")]
ggplot(pca_data, aes(x = PC_1, y = PC_2, col = seurat_clusters ,label=KO)) + geom_point(size = 2,alpha=0.9) + xlab(paste("PC1 (", varExplained[1], "%)"))+theme_classic()+
  ylab(paste("PC2 (",varExplained[2], "%)"))+ ggtitle("PCA")+scale_colour_manual(values = cc[c(3,1,2)],name="Clusters")+geom_text_repel()

FeaturePlot(Thomas_dots,features = c("IL1B","OASL"))

nii=ggplot(pca_data, aes(x = PC_1, y = PC_2, col = Stim. ,label=KO)) + geom_point(size = 2,alpha=0.9) + xlab(paste("PC1 (", varExplained[1], "%)"))+theme_classic()+
  ylab(paste("PC2 (",varExplained[2], "%)"))+ ggtitle("PCA")+scale_colour_manual(values = cc,name="Clusters")+geom_text_repel()
ggplotly(nii)

# 0-unstim
# 1-DNA
# 2-LPS
LPS_markers <- FindMarkers(Thomas_dots, ident.1 = 2, ident.2 = 0, min.pct = 0.25)
# LPS_markers <- FindMarkers(Thomas_dots, ident.1 = 2, ident.2 = c(0,1), min.pct = 0.25)
head(LPS_markers, n = 11)
top_lps=LPS_markers[1:100,]
top_lps=top_lps[order(top_lps$avg_logFC,decreasing = T),][1:10,]
top_lps
# IL1B and CCL4
DNA_markers <- FindMarkers(Thomas_dots, ident.1 = 1, ident.2 = 0,min.pct = 0.25)
# DNA_markers <- FindMarkers(Thomas_dots, ident.1 = 1, ident.2 = c(0,2),min.pct = 0.25)
head(DNA_markers, n = 11)
top_dna=DNA_markers[1:100,]
top_dna=top_dna[order(top_dna$avg_logFC,decreasing = T),][1:10,]
top_dna
# CXCL10 and OASL
info=as.data.frame(Thomas_dots@meta.data)
info=info[info$Stim.=="HT-DNA",]

FeaturePlot(Thomas_dots,features = c(rownames(top_dna)[1:3],rownames(top_lps)[1:3]),reduction = "pca")
FeaturePlot(Thomas_dots,features = c("NFKB1","OASL","CXCL10","IL6"),reduction = "pca")
pca_data=Thomas_dots@reductions$pca@cell.embeddings[,1:5]
pca_data=cbind(pca_data,Thomas_dots@meta.data)
pg=ggplot(pca_data, aes(x = PC_1, y = PC_2, col = Stim.,label=KO)) + geom_point(size = 1,alpha=0.9) + xlab(paste("PC1 (", varExp[1], "%)"))+theme_classic()+
  ylab(paste("PC2 (",varExp[2], "%)"))+ ggtitle("PCA of top 2000 HVGs")+scale_colour_manual(values = cc,name="Treatment")
fp=FeaturePlot(Thomas_dots,features = c("IL1B"),reduction = "pca")
plot_grid(pg,p1,fp)
ggplotly(pg)

# Manual PCA
# brb=as.data.frame(as.matrix(Thomas_dots@assays$RNA@data))
# mat=t(brb[VariableFeatures(Thomas_dots),])
# mat=t(brb)
# pc <- prcomp(mat, scale = T)
# pc.sum <- summary(pc)$importance
# df <- data.frame(pcs = 1:50, variance = pc.sum[2,][1:50], rule1 = pc.sum[3, ][1:50] > 0.9)
# varExp <- round(pc.sum[2, ] * 100, 2)
# pcs <- data.frame(pc$x)
# pcs <- cbind(pcs,Thomas_dots@meta.data)
# head(pcs)
# 
# pg=ggplot(pcs, aes(x = PC1, y = PC2, col = Stim.,label=KO)) + geom_point(size = 2,alpha=0.9) + xlab(paste("PC1 (", varExp[1], "%)"))+theme_classic()+
#   ylab(paste("PC2 (",varExp[2], "%)"))+ ggtitle("PCA of top 2000 HVGs")+scale_colour_manual(values = cc,name="Treatment")
# pg
# ggplotly(pg)

# pcs1=pcs[pcs$Stim.=="Unstim.",]
# p1=ggplot(pcs1, aes(x = PC1, y = PC2, col = KO,label=KO)) + geom_point(size = 2,alpha=0.9)+theme_classic()+
#   ggtitle("PCA of top 2000 HVGs Unstimulated samples only")+scale_colour_discrete(name="KO")
# ggplotly(p1)
#################################
# pcs2=pcs[pcs$Stim.=="LPS",]
# p2=ggplot(pcs2, aes(x = PC1, y = PC2, col = KO,label=KO)) + geom_point(size = 2,alpha=0.9)+theme_classic()+
#   ggtitle("PCA of top 2000 HVGs LPS stim. samples only")+scale_colour_discrete(name="KO")
# ggplotly(p2)
#################################
# pcs3=pcs[pcs$Stim.=="HT-DNA",]
# p3=ggplot(pcs3, aes(x = PC1, y = PC2, col = KO,label=KO)) + geom_point(size = 2,alpha=0.9)+theme_classic()+
#   ggtitle("PCA of top 2000 HVGs DNA stim. samples only")+scale_colour_discrete(name="KO")
# ggplotly(p3)
#################################

head(cell_annot)

PCA_data=as.data.frame(Thomas_dots@reductions$pca@cell.embeddings)
head(PCA_data)
head(Thomas_dots@meta.data)
PCA_data=cbind(PCA_data,Thomas_dots@meta.data)
" TOO MUCH INFORMATION"
# ggplot(PCA_data,aes(x=PC_1,y=PC_2,color=Stim.,label=Type))+geom_point(size=2)+ggrepel::geom_text_repel(size=3)+scale_color_manual(values = cc)+theme_classic()

# tSNE
Thomas_dots=RunTSNE(Thomas_dots,perplexity=15)
# Umap
Thomas_dots=RunUMAP(Thomas_dots, dims = 1:50, n.neighbors = 15,umap.method = "umap-learn",metric = "correlation",min.dist = 0.3,verbose = FALSE)
p1=DimPlot(Thomas_dots,reduction = "pca",group.by = "Stim.",pt.size = 2)+ggtitle("PCA")+xlab("")+ylab("")+scale_color_manual(values = cc)+xlab("")+ylab("")
p2=DimPlot(Thomas_dots,reduction = "tsne",group.by = "Stim.",pt.size = 2)+ggtitle("tSNE")+xlab("")+ylab("")+scale_color_manual(values = cc)+xlab("")+ylab("")
p3=DimPlot(Thomas_dots,reduction = "umap",group.by = "Stim.",pt.size = 2)+ggtitle("UMAP")+xlab("")+ylab("")+scale_color_manual(values = cc)+xlab("")+ylab("")
plot_grid(p1,p2,p3)

DimPlot(Thomas_dots,reduction = "umap",group.by = "seurat_clusters",pt.size = 2)+ggtitle("UMAP")+xlab("")+ylab("")+scale_color_manual(values = cc)+xlab("")+ylab("")

###################
Tsne_data=as.data.frame(Thomas_dots@reductions$tsne@cell.embeddings)
head(Tsne_data)
Tsne_data=cbind(Tsne_data,Thomas_dots@meta.data)
# ggplot(Tsne_data,aes(x=tSNE_1,y=tSNE_2,color=Stim.,label=Type))+geom_point(size=2)+ggrepel::geom_text_repel(size=3)+scale_color_manual(values = cc)+theme_classic()
###################
UMAP_data=as.data.frame(Thomas_dots@reductions$umap@cell.embeddings)
head(UMAP_data)
UMAP_data=cbind(UMAP_data,Thomas_dots@meta.data)
# ggplot(UMAP_data,aes(x=UMAP_1,y=UMAP_2,color=Stim.,label=Type))+geom_point(size=2)+ggrepel::geom_text_repel(size=3)+scale_color_manual(values = cc)+theme_classic()
###################
sel=c("IL1B","IL6","TNF","CCL4","OASL","NLRP3")
VlnPlot(Thomas_dots, features = sel, slot = "counts", log = TRUE)
VlnPlot(Thomas_dots, features = sel,group.by = "Stim.")
# VlnPlot(Thomas_dots, features = sel,group.by = "KO")
### HALF Violins would be awesome!
FeaturePlot(Thomas_dots, features = sel,pt.size = 2,reduction = "pca")

#############################################################################################################################################

scale_data=as.data.frame(as.matrix(Thomas_dots@assays$RNA@data))
info_data=Thomas_dots@meta.data
head(scale_data)

library('gplots')
info_data$info=paste(info_data$Info,info_data$KO,info_data$Well,sep = "_")
colnames(scale_data)
colnames(scale_data)=info_data[match(colnames(scale_data),info_data$Full_BC),"info"]
# heatmap.2(cor(scale_data,method = "spearman"), trace='none', main='Sample correlations (raw)',density.info = "none")
heatmap.2(cor(scale_data,method = "pearson"), trace='none', main='Sample correlations (raw)',density.info = "none" )

library("PoiClaClu")
# poisd=PoissonDistance(scale_data)
# samplePoisDistMatrix=as.matrix(poisd$dd)
# rownames(samplePoisDistMatrix)=info_data$info
# colnames(samplePoisDistMatrix)=NULL
# c3=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(samplePoisDistMatrix,
#          clustering_distance_rows = poisd$dd,
#          clustering_distance_cols = poisd$dd,
#          col = c3)

Top_genes=VariableFeatures(Thomas_dots)
Top_genes=Top_genes[!grepl("^AC[0-9]",Top_genes)]
Top_genes=Top_genes[!grepl("^AL[0-9]",Top_genes)]
Top_genes=Top_genes[!grepl("^AP[0-9]",Top_genes)]
head(Top_genes,20)
Top_genes=Top_genes[!grepl("orf",Top_genes)]
head(Top_genes,20)
Top_genes=Top_genes[!grepl("MIR",Top_genes)]
head(Top_genes,20)
Top_genes=Top_genes[1:200]
Top_genes=unique(c(Top_genes,"CH25H","IFIT5","SAMD9L","IFI44L","OAS2","HELZ2","MX2","IFIT3","ISG15","IFIT1","IFIT2","OASL","HERC6","IFI44"))
matrix=as.matrix(scale_data[Top_genes,])
mat_col=data.frame(Stim=info_data$Stim.,KO=info_data$KO)
rownames(mat_col)=colnames(matrix)
library("viridis")
pheatmap(matrix, border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col,
         drop_levels=T, main = "Heatmap top 200HVGs",color = rev(redgreen(75)))
# heatmap.2(matrix, trace='none', main="Gene clustering", density.info = "none" ,)

###############################################################################################################################################
###############################################################################################################################################
# FOR VEIT AND FIONAN!
###############################################################################################################################################
###############################################################################################################################################
# sel=c("IL1B","IL6","TNF","CCL4","OASL","NLRP3")
# # sel=sel[sel%in%rownames(Thomas_dots)]
# brb_data=as.data.frame(as.matrix(Thomas_dots@assays$RNA@counts))
# sel_bar=brb_data[sel,]
# sel_bar$Gene=rownames(sel_bar)
# # head(sel_bar)
# sel_bar=gather(sel_bar,key = XC,value = Exprs,-Gene)
# # head(sel_bar)
# sel_bar$Stim=Thomas_dots@meta.data$Stim.[match(sel_bar$XC,Thomas_dots@meta.data$Full_BC)]
# sel_bar$KO=Thomas_dots@meta.data$KO[match(sel_bar$XC,Thomas_dots@meta.data$Full_BC)]
# sel_bar$Name=Thomas_dots@meta.data$Type[match(sel_bar$XC,Thomas_dots@meta.data$Full_BC)]
# # head(sel_bar)
# SEM <- function(x) sd(x)/sqrt(length(x))
# sel_bar$Exprs=2^sel_bar$Exprs
# 
# mean_sel_data=sel_bar %>% group_by(Gene, Stim, KO)%>%
#   summarise(Value = mean(Exprs),count=length(Exprs),sem = SEM(Exprs))
# 
# head(mean_sel_data)
# mean_sel_data$Val=mean_sel_data$Value
# # mean_sel_data$Val=2^mean_sel_data$Value
# head(mean_sel_data)
# max=max(mean_sel_data$Val)
# max=max+(0.1*max)
# # WT
# mean_sel_data1=mean_sel_data[mean_sel_data$KO=="WT",]
# p1=ggplot(mean_sel_data1, aes(x=Stim, y=Val, fill=Stim)) + geom_bar(stat = "identity", position = "dodge", width=0.9,alpha=0.8)+
#   geom_errorbar(data = mean_sel_data1,aes(ymin = Val-sem, ymax = Val+sem),width=0.9,colour="black", alpha=1, size=0.3,position = position_dodge())+
#   ggtitle("Thomas Data - WT")+xlab("")+ylab('RAW UMI COUNTS')+theme_classic()+scale_fill_manual(values=cc) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="")+facet_grid(~Gene)+ylim(c(0,max))
# # rela
# mean_sel_data2=mean_sel_data[mean_sel_data$KO=="Rela",]
# p2=ggplot(mean_sel_data2, aes(x=Stim, y=Val, fill=Stim)) + geom_bar(stat = "identity", position = "dodge", width=0.9,alpha=0.8)+
#   geom_errorbar(data = mean_sel_data2,aes(ymin = Val-sem, ymax = Val+sem),width=0.9,colour="black", alpha=1, size=0.3,position = position_dodge())+
#   ggtitle("Thomas Data - Rela KO")+xlab("")+ylab('RAW UMI COUNTS')+theme_classic()+scale_fill_manual(values=cc) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="")+facet_grid(~Gene)+ylim(c(0,max))
# # Ikka
# mean_sel_data3=mean_sel_data[mean_sel_data$KO=="IKKa",]
# p3=ggplot(mean_sel_data3, aes(x=Stim, y=Val, fill=Stim)) + geom_bar(stat = "identity", position = "dodge", width=0.9,alpha=0.8)+
#   geom_errorbar(data = mean_sel_data3,aes(ymin = Val-sem, ymax = Val+sem),width=0.9,colour="black", alpha=1, size=0.3,position = position_dodge())+
#   ggtitle("Thomas Data - IKKA KO")+xlab("")+ylab('RAW UMI COUNTS')+theme_classic()+scale_fill_manual(values=cc) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="")+facet_grid(~Gene)+ylim(c(0,max))
# # IKKb
# mean_sel_data4=mean_sel_data[mean_sel_data$KO=="IKKb",]
# p4=ggplot(mean_sel_data4, aes(x=Stim, y=Val, fill=Stim)) + geom_bar(stat = "identity", position = "dodge", width=0.9,alpha=0.8)+
#   geom_errorbar(data = mean_sel_data4,aes(ymin = Val-sem, ymax = Val+sem),width=0.9,colour="black", alpha=1, size=0.4,position = position_dodge())+
#   ggtitle("Thomas Data - IKKB KO")+xlab("")+ylab('RAW UMI COUNTS')+theme_classic()+scale_fill_manual(values=cc) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="")+facet_grid(~Gene)+ylim(c(0,max))
# # IKKab
# mean_sel_data5=mean_sel_data[mean_sel_data$KO=="IKKa_IKKb",]
# p5=ggplot(mean_sel_data5, aes(x=Stim, y=Val, fill=Stim)) + geom_bar(stat = "identity", position = "dodge", width=0.9,alpha=0.8)+
#   geom_errorbar(data = mean_sel_data5,aes(ymin = Val-sem, ymax = Val+sem),width=0.9,colour="black", alpha=1, size=0.3,position = position_dodge())+
#   ggtitle("Thomas Data - IKKA/B double KO")+xlab("")+ylab('RAW UMI COUNTS')+theme_classic()+scale_fill_manual(values=cc) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="")+facet_grid(~Gene)+ylim(c(0,max))
# plot_grid(p1,p2,p3,p4,p5,ncol = 1)
# write.csv(mean_sel_data,"/Home Office/Thomas_data_RAW_do_not_use.csv")

brb_data=as.data.frame(as.matrix(Thomas_dots@assays$RNA@data))
brb_data$Gene=rownames(brb_data)
long_brb=gather(brb_data,key = XC,value = Value,-Gene)

long_brb$Stim=Thomas_dots@meta.data$Stim.[match(long_brb$XC,Thomas_dots@meta.data$Full_BC)]
long_brb$KO=Thomas_dots@meta.data$KO[match(long_brb$XC,Thomas_dots@meta.data$Full_BC)]
long_brb$Name=Thomas_dots@meta.data$Type[match(long_brb$XC,Thomas_dots@meta.data$Full_BC)]
head(long_brb)
mean_brb=long_brb %>% group_by(Gene, Stim, KO)%>%
  summarise(Value = mean(Value))
mean_brb$Sample=paste(mean_brb$KO,mean_brb$Stim,sep = "_")
head(mean_brb)
wide_brb=mean_brb %>% spread(Gene, Value)
wide_brb=as.data.frame(wide_brb)
wide_brb[1:5,1:5]
rownames(wide_brb)=wide_brb$Sample
wide_brb=wide_brb[,-c(1:3)]
dim(wide_brb)
wide_brb[1:5,1:5]
wide_brb=t(wide_brb)
head(wide_brb)

#top genes
# Top_genes=VariableFeatures(Thomas_dots)
# Top_genes=Top_genes[!grepl("^AC[0-9]",Top_genes)]
# Top_genes=Top_genes[!grepl("^AL[0-9]",Top_genes)]
# Top_genes=Top_genes[!grepl("^AP[0-9]",Top_genes)]
# Top_genes=Top_genes[!grepl("orf",Top_genes)]
# Top_genes=Top_genes[!grepl("MIR",Top_genes)]
# I throw away a few weird genes - based on their expression pattern
Top_genes=Top_genes[!Top_genes%in%c("NEBL","INO80D","ZC3HAV1","ZCCHC7")]
Top_genes=Top_genes
# [1:200]

heat_mat=wide_brb[Top_genes,]
mat_col=data.frame(Stim=colnames(heat_mat),Name=colnames(heat_mat))
mat_col$Stim[grepl("LPS",mat_col$Name)]="LPS"
mat_col$Stim[grepl("DNA",mat_col$Name)]="HT-DNA"
mat_col$Stim[grepl("Unstim",mat_col$Name)]="Unstimulated"
rownames(mat_col)=colnames(heat_mat)
mat_col=mat_col[-2]
pheatmap(heat_mat, border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col,
         drop_levels=T, main = "Heatmap top 200HVGs",color = rev(redgreen(75)),cutree_rows = 3,cutree_cols = 3)

head(heat_mat)
colnames(heat_mat)
write.csv(colnames(heat_mat),"/Home Office/temp/ssss.csv")
xxx=read.csv("/Home Office/temp/ssss.csv")

xxx$x[match(colnames(heat_mat),xxx$x)]
colnames(heat_mat)[match(colnames(heat_mat),xxx$x)]
colnames(heat_mat)%in%xxx$x
heat_matx=heat_mat[,match(colnames(heat_mat),xxx$x)]
heat_matx=heat_mat[,match(xxx$x,colnames(heat_mat))]
# heat_mat[,match(xxx$x,colnames(heat_mat))]
# xxx$x[match(xxx$x,colnames(heat_mat))]
pheatmap(heat_matx, border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col,
         drop_levels=T, main = "Heatmap top 200HVGs",color = rev(redgreen(75)),cutree_rows = 3,cutree_cols = 3,cluster_cols = F)




heat_mat2=heat_mat[c("IL1B","CCL4","CXCL10","OASL"),]
pheatmap(heat_mat2, border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col,
         drop_levels=T, main = "Heatmap top 200HVGs",color = rev(redgreen(75)),cutree_cols = 3)

########## LINE PLOTS
c("IL1B","CCL4","CXCL10","OASL")

heat_mat2=as.data.frame(heat_mat2)
heat_mat2$Gene=rownames(heat_mat2)
head(heat_mat2)
heat_mat3=gather(heat_mat2,"Sample","Value",-Gene)
head(heat_mat3)
# heat_mat3=heat_mat2 %>% group_by(info,Gene)%>%
#   summarise(Means = mean(Value))
# head(plotting_data2)
##### LINE PLOT
nii=as.data.frame(Thomas_dots@meta.data)
head(nii)
moo=data.frame(Full=paste0(nii$KO,"_",nii$Stim.),KO=nii$KO,Stim=nii$Stim.)
moo=moo[!duplicated(moo$Full),]
moo[order(moo$Stim[order(moo$KO)]),]
moo$KO=as.character(moo$KO)
moo$Stim=as.character(moo$Stim)
moo[with(moo, order(Stim, KO)), ]
moo=moo[order(moo[,3], moo[,2]), ]

nii=unique(heat_mat3$Sample)
nii=data.frame(full=nii,stim=nii,name=nii)

heat_mat3$Sample=factor(heat_mat3$Sample,levels = moo$Full)

ggplot(heat_mat3, aes(x=Sample, y=Value, color=Gene,group=Gene)) + geom_line(size=1.2)+ theme_light() + 
  ylab("Norm. expression")+xlab("")+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+scale_color_aaas()



mat_cluster_cols <- hclust(dist(t(heat_mat)))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
trees=cutree(mat_cluster_cols,k = 3)

trees[trees==1]
trees[trees==2]
trees[trees==3]



# Group number for a gene
average_heatmap=pheatmap(heat_mat, border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col,
                         drop_levels=T, main = "Heatmap top 200HVGs",color = rev(redgreen(75)),cutree_rows = 3,cutree_cols = 3,silent = TRUE)
pheatmap(heat_mat, border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col,
                         drop_levels=T, main = "Heatmap top 200HVGs",color = rev(redgreen(75)),cutree_rows = 3,cutree_cols = 3)

avg_hm_data=data.frame(names=names(cutree(average_heatmap$tree_row,k = 3)),ass=cutree(average_heatmap$tree_row,k = 3))
#################################################################
# Full matrix
heat_mat2=as.data.frame(as.matrix(Thomas_dots@assays$RNA@data))
heat_mat2=heat_mat2[Top_genes,]
mat_col2=data.frame(Stim=Thomas_dots@meta.data$Stim.,KO=Thomas_dots@meta.data$Type)
head(mat_col2)

rownames(mat_col2)=mat_col2$KO
colnames(heat_mat2)=mat_col2$KO

pheatmap(heat_mat2, annotation_legend = F,border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col2,
         drop_levels=T, main = "Heatmap top 200HVGs",color = rev(redgreen(75)),cutree_rows = 5,cutree_cols = 3)

# Group number for a gene
all_heatmap=pheatmap(heat_mat2, annotation_legend = F,border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col2,
                     drop_levels=T, main = "Heatmap top 200HVGs",color = rev(redgreen(75)),cutree_rows = 5,cutree_cols = 3,silent = T)


# all_hm_data=data.frame(names=names(cutree(all_heatmap$tree_row,k = 5)),ass=cutree(all_heatmap$tree_row,k = 5))
# avg_hm_data$Cluster="Not_sure"
# avg_hm_data$ass[avg_hm_data$names=="IL6"]
# avg_hm_data$ass[avg_hm_data$names=="CXCL8"]
# avg_hm_data$Cluster[avg_hm_data$ass==3]="LPS"
# ################################################
# avg_hm_data$ass[avg_hm_data$names=="CXCL10"]
# avg_hm_data$ass[avg_hm_data$names=="WARS"]
# avg_hm_data$ass[avg_hm_data$names=="IFIT1"]
# avg_hm_data$ass[avg_hm_data$names=="ISG15"]
# avg_hm_data$ass[avg_hm_data$names=="IFIT3"]
# avg_hm_data$ass[avg_hm_data$names=="IFIT2"]
# avg_hm_data$ass[avg_hm_data$names=="OASL"]
# avg_hm_data$Cluster[avg_hm_data$ass==1]="HTDNA"
# ###############################################
# all_hm_data$Cluster="Not_sure"
# all_hm_data$ass[all_hm_data$names=="IL6"]
# all_hm_data$ass[all_hm_data$names=="CXCL8"]
# all_hm_data$Cluster[all_hm_data$ass==5]="LPS"
# all_hm_data$Cluster[all_hm_data$ass==4]="LPS"
# all_hm_data$ass[all_hm_data$names=="CXCL10"]
# all_hm_data$ass[all_hm_data$names=="WARS"]
# all_hm_data$ass[all_hm_data$names=="IFIT1"]
# all_hm_data$ass[all_hm_data$names=="ISG15"]
# all_hm_data$ass[all_hm_data$names=="IFIT3"]
# all_hm_data$ass[all_hm_data$names=="IFIT2"]
# all_hm_data$ass[all_hm_data$names=="OASL"]
# all_hm_data$ass[all_hm_data$names=="CCL2"]
# all_hm_data$ass[all_hm_data$names=="CCL8"]
# all_hm_data$Cluster[all_hm_data$ass==1]="HTDNA"
# all_hm_data$Cluster[all_hm_data$ass==3]="HTDNA"
# tot_data=cbind(all_hm_data,avg_hm_data)

# VIOLIN PLOTS OF OSAL AND CXCL MAybe do it based on the 3 groups on heatmaps? 
nii=data.frame(Full_BC=Thomas_dots$Full_BC,Cluster=cutree(all_heatmap$tree_col,k = 3))
nii$Clust_Name="Unstimulated Cluster"
nii$Clust_Name[nii$Cluster==2]="LPS Cluster"
nii$Clust_Name[nii$Cluster==3]="HTDNA Cluster"
head(nii)
tail(nii)

Viol_matr=as.data.frame(as.matrix(Thomas_dots@assays$RNA@data))
Viol_matr=Viol_matr[c("OASL","CXCL10","NLRP3","IL1B"),]
Viol_matr$Gene=rownames(Viol_matr)
Viol_matr=Viol_matr %>% gather(Full_BC,Value,-Gene)
Viol_matr$Cluster=nii$Clust_Name[match(Viol_matr$Full_BC,nii$Full_BC)]
# Viol_matr$Cluster=nii$Clust_Name[match(Viol_matr$Full_BC,nii$Full_BC)]

source("/Home Office/My_awesome_RNA-seq_Project/data_viz/Geom_Splitviolin.R")
# Viol_matr$Value=2^Viol_matr$Value
Viol_matr$Gene=factor(Viol_matr$Gene,levels = unique(Viol_matr$Gene)[c(2,1,4,3)])
Viol_matr2=Viol_matr[Viol_matr$Cluster!="LPS Cluster",]
head(Viol_matr2)
v1=ggplot(Viol_matr2, aes(x = Gene, y = Value, fill = Cluster,color=Cluster)) + geom_split_violin()+theme_classic()+
  scale_fill_manual(values = cc)+ggtitle("HT-DNA stimulation")+scale_color_manual(values = cc)+
  ylim(0,4)
Viol_matr3=Viol_matr[Viol_matr$Cluster!="HTDNA Cluster",]
v2=ggplot(Viol_matr3, aes(x = Gene, y = Value, fill = Cluster,color=Cluster)) + geom_split_violin()+theme_classic()+
  scale_fill_manual(values = cc)+ggtitle("LPS stimulation")+scale_color_manual(values = cc)+
  ylim(0,4)
plot_grid(v1,v2)

Viol_matr$Cluster=factor(Viol_matr$Cluster,levels = unique(Viol_matr$Cluster)[c(2,1,3)])
library(ggbeeswarm)
ggplot(Viol_matr, aes(x = Gene, y = Value, fill = Cluster)) + geom_boxplot()+#geom_beeswarm(dodge.width = 0.7,cex = 0.7)+
  theme_light()+scale_fill_manual(values = cc)+ggtitle("HT-DNA vs. LPS vs. Unstimulated")+scale_color_manual(values = cc)


############### MANUAL PCA on average data
head(wide_brb)
mat=wide_brb
rowVar <- apply(mat, 1, var)
mv2000 <- order(rowVar, decreasing = T)[1:500]
mat <- t(mat[mv2000,])
pc <- prcomp(mat, scale = T)
pc.sum <- summary(pc)$importance
pc.sum[, 1:5]
df <- data.frame(pcs = 1:dim(pc.sum)[2], variance = pc.sum[2,], rule1 = pc.sum[3, ] > 0.8)
varExp=round(pc.sum[2, ] * 100, 2)
pcs=data.frame(pc$x, stringsAsFactors = F)
write.csv(pcs,"/Users/Ashto/Desktop/moooooo.csv")
nii=read.csv("/Users/Ashto/Desktop/moooooo.csv")
nii=nii[-1]
rownames(nii)=nii$info
head(nii)

res=ggplot(nii, aes(PC1, PC2, color=stim,label=ko))+
  geom_point(size = 2) + theme_classic()+ggtitle("PCA based on top 500 HVGs and averaged samples")+
  xlab(paste("PC1 (", varExp[1], "%)"))+ ylab(paste("PC2 (",varExp[2], "%)"))+
  scale_color_manual(values =cc)
ggplotly(res)

geenid=c("CCL4L2","CCL4","IFIT2","OASL","UCK2")

hmmm=data.frame(names=colnames(wide_brb),CCL4L2=wide_brb["CCL4L2",],CCL4=wide_brb["CCL4",],IFIT2=wide_brb["IFIT2",],OASL=wide_brb["OASL",],UCK2=wide_brb["UCK2",])
head(hmmm)
head(nii)
nii=cbind(nii,hmmm[match(nii$info,hmmm$names),])
head(hmmm)

p1=ggplot(nii, aes(PC1, PC2, color=CCL4L2,label=ko))+
  geom_point(size = 2) + theme_classic()+ggtitle("PCA based on top 500 HVGs and averaged samples")+
  xlab(paste("PC1 (", varExp[1], "%)"))+ ylab(paste("PC2 (",varExp[2], "%)"))

p2=ggplot(nii, aes(PC1, PC2, color=CCL4,label=ko))+
  geom_point(size = 2) + theme_classic()+ggtitle("PCA based on top 500 HVGs and averaged samples")+
  xlab(paste("PC1 (", varExp[1], "%)"))+ ylab(paste("PC2 (",varExp[2], "%)"))

p3=ggplot(nii, aes(PC1, PC2, color=IFIT2,label=ko))+
  geom_point(size = 2) + theme_classic()+ggtitle("PCA based on top 500 HVGs and averaged samples")+
  xlab(paste("PC1 (", varExp[1], "%)"))+ ylab(paste("PC2 (",varExp[2], "%)"))

p4=ggplot(nii, aes(PC1, PC2, color=OASL,label=ko))+
  geom_point(size = 2) + theme_classic()+ggtitle("PCA based on top 500 HVGs and averaged samples")+
  xlab(paste("PC1 (", varExp[1], "%)"))+ ylab(paste("PC2 (",varExp[2], "%)"))

p5=ggplot(nii, aes(PC1, PC2, color=UCK2,label=ko))+
  geom_point(size = 2) + theme_classic()+ggtitle("PCA based on top 500 HVGs and averaged samples")+
  xlab(paste("PC1 (", varExp[1], "%)"))+ ylab(paste("PC2 (",varExp[2], "%)"))

plot_grid(p1,p2,p3,p4,p5)
##################################################
# DE between groups

mat_cluster_cols <- hclust(dist(t(heat_mat)))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
trees=cutree(mat_cluster_cols,k = 3)

trees[trees==1]
trees[trees==2]
trees[trees==3]
head(nii)
naa=nii[,c(39:42)]
head(naa)
naa$clusters="Unstimulated"
naa$clusters[naa$info%in%names(trees[trees==1])]="HT-DNA"
naa$clusters[naa$info%in%names(trees[trees==3])]="LPS"
naa

library(DESeq2)
colnames(wide_brb)
wide_brb2=wide_brb[,order(colnames(wide_brb))]
# wide_brb2=wide_brb2*10
wide_brb2=round(2^wide_brb2)

wide_brb2=as.matrix(wide_brb2)
wide_brb2=as.data.frame(wide_brb2)
sum(is.na(wide_brb2))
wide_brb2[1:5,1:5]
write.csv(wide_brb2,"/Users/Ashto/Desktop/naasdka.csv")
again=read.csv("/Users/Ashto/Desktop/naasdka.csv")
again=again[-1]
again[1:10,1:10]
again=as.matrix(again)

write.csv(naa,"/Users/Ashto/Desktop/sanpledata-naasdka.csv")
naa=read.csv("/Users/Ashto/Desktop/sanpledata-naasdka.csv")

dds= DESeqDataSetFromMatrix(countData = wide_brb2, colData = naa ,design = ~ clusters)

# dds= DESeqDataSetFromMatrix(countData = again, colData = naa ,design = ~ clusters)

dds
dds <- DESeq(dds)
#################################################################################################
resultsNames(dds)
# CONTRASTS
# WTs
wt_vs_lps=as.data.frame(results(dds, contrast=c("clusters", "LPS", "Unstimulated")))
wt_vs_lps=wt_vs_lps[order(wt_vs_lps$pvalue),]
wt_vs_lps=na.omit(wt_vs_lps)
head(wt_vs_lps)

wt_vs_DNA=as.data.frame(results(dds, contrast=c("clusters", "HT-DNA", "Unstimulated")))
wt_vs_DNA=wt_vs_DNA[order(wt_vs_DNA$pvalue),]
wt_vs_DNA=na.omit(wt_vs_DNA)
head(wt_vs_DNA)

# add a grouping column; default value is "not significant"
volcano_data=wt_vs_lps
volcano_data["group"] <- "NotSignificant"
volcano_data[base::which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) < 1.5 ),"group"] <- "Significant"
volcano_data[which(volcano_data['padj'] > 0.05 & abs(volcano_data['log2FoldChange']) > 1.5 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) > 1.5 ),"group"] <- "Significant&FoldChange"
# to manualy set colors
pal <- c("blue", "black", "purple","red")
volcano_data$Symbol=rownames(volcano_data)
plot_ly(volcano_data, x=~log2FoldChange, y=~-log10(padj),
        color=~group,colors=pal,
        text=~Symbol,type="scatter",size = ~-log10(padj),
        mode="markers") %>%
  layout(title = 'Unstimulated vs. cGAMP + CD3/CD28',
         yaxis = list(zeroline = T),
         xaxis = list(zeroline = T))


# add a grouping column; default value is "not significant"
volcano_data=wt_vs_DNA
volcano_data["group"] <- "NotSignificant"
volcano_data[base::which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) < 1.5 ),"group"] <- "Significant"
volcano_data[which(volcano_data['padj'] > 0.05 & abs(volcano_data['log2FoldChange']) > 1.5 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) > 1.5 ),"group"] <- "Significant&FoldChange"
# to manualy set colors
pal <- c("blue", "black", "purple","red")
volcano_data$Symbol=rownames(volcano_data)
plot_ly(volcano_data, x=~log2FoldChange, y=~-log10(padj),
        color=~group,colors=pal,
        text=~Symbol,type="scatter",size = ~-log10(padj),
        mode="markers") %>%
  layout(title = 'Unstimulated vs. cGAMP + CD3/CD28',
         yaxis = list(zeroline = T),
         xaxis = list(zeroline = T))
