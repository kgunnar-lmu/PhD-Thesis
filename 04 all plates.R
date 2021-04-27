# Generate stats for further filtering and analysis. 
# Filtering is the next script
library(tidyverse)
library(mclust)
library(stringr)
library(Seurat)
library(sctransform)
library(scran)
library(scater)
library(Rtsne)
library(rsvd)
library(Hmisc)
library(reshape2)
library(plotly)
library(RColorBrewer)
library(pheatmap)
library(NMF)
library(MAST)
library(limma)
library(knitr)
library(gtools)
library(GSEABase)
library(ggsci)
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)
library(cowplot)
library(biomaRt)
options(stringsAsFactors = F)
cc=c(brewer.pal(name = "Set1",n=9))

df1=first_working_plate
dim(first_working_plate)
df2=ribo_plate
dim(ribo_plate)
df3=PK_plate
dim(PK_plate)

colnames(df1)=paste0("A","_",colnames(df1))
colnames(df2)=paste0("B","_",colnames(df2))
colnames(df3)=paste0("C","_",colnames(df3))
head(df3)
n1=rownames(df1)
n2=rownames(df2)
n3=rownames(df3)
#
n1_2=n1[n1%in%n2]
n_all=n1_2[n1_2%in%n3]
#
length(n_all)

df1=df1[n_all,]
df2=df2[n_all,]
df3=df3[n_all,]

df_all=cbind(df1,df2,df3)
df_all[1:5,1:5]
head(df_all)

# Seurat
all_dots=CreateSeuratObject(counts = df_all, project = "all plate")
all_dots[["percent.mt"]]=PercentageFeatureSet(all_dots, pattern = "^MT-")
VlnPlot(all_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Mito_cut=10
RNA_cut=2500
minRNA=500
all_dots=subset(all_dots, subset = nFeature_RNA < RNA_cut & nFeature_RNA > minRNA & percent.mt< Mito_cut)

v1=VlnPlot(all_dots, features = c("nFeature_RNA"))+scale_fill_manual(values = cc)
v2=VlnPlot(all_dots, features = c("nCount_RNA"))+scale_fill_manual(values = cc)
v3=VlnPlot(all_dots, features = c("percent.mt"))+scale_fill_manual(values = cc)

plot_grid(v1,v2,v3,ncol = 3)
plot_grid(v1,v2,v3,ncol = 3)

FeatureScatter(all_dots, feature1 = "nFeature_RNA", feature2 = "percent.mt")
dim(all_dots)
# Info
all_dots$plate=substr(colnames(all_dots),1,1)
all_dots$stim="Unstim."
all_dots$well=substr(colnames(all_dots),3,6)
all_dots$stim[grepl("7|8|9|10|11|12",all_dots$well)]="LPS"
all_dots$name=colnames(all_dots)
all_dots$UMIs=colSums(all_dots@assays$RNA@counts)
head(all_dots@meta.data)
# Normalizing with SCtransform
# all_dots=SCTransform(all_dots)
# cof_var=c("UMIs")
# all_dots=SCTransform(all_dots,vars.to.regress = cof_var)
# all_dots=all_dots[,all_dots$plate!="A"]
cof_var=c("plate","UMIs")
all_dots=SCTransform(all_dots,vars.to.regress = cof_var,return.only.var.genes = F)
# PCA
all_dots=RunPCA(all_dots, npcs = 50, verbose = FALSE)
p1=DimPlot(all_dots,reduction = "pca",group.by = "stim",pt.size = 2)+scale_color_manual(values = cc)
p2=DimPlot(all_dots,reduction = "pca",group.by = "plate",pt.size = 2)+scale_color_manual(values = cc)
plot_grid(p1,p2)
# Umap
all_dots=RunUMAP(all_dots, dims = 1:50, n.neighbors = 30,umap.method = "umap-learn", 
                  metric = "correlation",min.dist = 0.2,verbose = FALSE)

DimPlot(all_dots,reduction = "umap",group.by = "stim",pt.size = 2)+ggtitle("all 3 plates")+
  xlab("")+ylab("")+scale_color_manual(values = cc)+xlab("")+ylab("")
DimPlot(all_dots,reduction = "umap",group.by = "plate",pt.size = 2)+ggtitle("all 3 plates")+
  xlab("")+ylab("")+scale_color_manual(values = cc)+xlab("")+ylab("")

# all_dots=FindNeighbors(all_dots,dims = 1:50,compute.SNN = T)
# all_dots=FindClusters(all_dots,resolution = 0.6)
# DimPlot(all_dots,reduction = "umap",pt.size = 2) + xlab("")+ylab("")+scale_color_manual(values = cc)

# all_dots$kmeans=as.character(kmeans(all_dots@reductions$umap@cell.embeddings, centers = 2)$cluster)
# cc=c(brewer.pal(name = "Dark2",n=8))

# DimPlot(all_dots,reduction = "pca",group.by = "kmeans",pt.size = 2)+scale_color_manual(values = cc)

# u1=DimPlot(all_dots,reduction = "umap",pt.size = 2,group.by = "stim") + xlab("")+ylab("")+scale_color_manual(values = cc)+labs(color="Stimulus")+theme(legend.position = "none")
# u2=DimPlot(all_dots,reduction = "umap",pt.size = 2,group.by = "plate") + xlab("")+ylab("")+scale_color_manual(values = cc)+labs(color="Batch")+theme(legend.position = "none")
# u3=DimPlot(all_dots,reduction = "umap",pt.size = 2,group.by = "kmeans") + xlab("")+ylab("")+scale_color_manual(values = cc)+labs(color="K-means")+theme(legend.position = "none")
# 
# plot_grid(p1,p2,p3)
# plot_grid(u1,u3,u2,ncol = 3)

u1=DimPlot(all_dots,reduction = "umap",pt.size = 2,group.by = "plate") + xlab("")+ylab("")+scale_color_manual(values = cc)+labs(color="PlateXXXXX")
u2=DimPlot(all_dots,reduction = "umap",pt.size = 2,group.by = "stim") + xlab("")+ylab("")+scale_color_manual(values = cc)+labs(color="stimuXXXXX")
# u3=DimPlot(all_dots,reduction = "umap",pt.size = 2,group.by = "kmeans") + xlab("")+ylab("")+scale_color_manual(values = cc)+labs(color="K-meansXXX")
all_dots=FindNeighbors(all_dots,dims = 1:40,compute.SNN = T)
all_dots=FindClusters(all_dots,resolution = 0.6,algorithm = 2 ,method = "igraph")
u4=DimPlot(all_dots,reduction = "umap",pt.size = 2) + xlab("")+ylab("")+scale_color_manual(values = cc[c(2,1,3,4,5,6,7,8)])+labs(color="SNN-Louvain")
u4
# plot_grid(u1,u2,u3,ncol = 3)
plot_grid(u1,u2,u4,ncol = 3)
FeaturePlot(all_dots,features = c("CCL4","IL1B","SOD2","CCL3L1"),reduction = "umap",ncol = 2,order = T,pt.size = 2)
############################################################

# Idents(all_dots)=all_dots$kmeans
dim(all_dots@assays$RNA@counts)
dim(all_dots@assays$RNA@data)


all_genes=FindMarkers(object = all_dots, ident.1 = 1,test.use = "DESeq2")
dim(all_genes)
volcano_data=all_genes
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-15]=volcano_data$Symbol[volcano_data$p_val_adj<1e-15]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
library(ggrepel)
pal <- c("blue", "black", "purple","red")
ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("LPS vs control")+xlim(-2,2) + geom_vline(xintercept = 0)
#

# Keep only up genes in middle and FDR 0.05% or less
marker_genes=all_genes[all_genes$p_val_adj<0.05,]
marker_genes=marker_genes[marker_genes$avg_logFC>0,]
marker_genes=marker_genes[order(marker_genes$p_val_adj),]
dim(marker_genes)
head(marker_genes)
found_genes=rownames(marker_genes)
# VlnPlot(all_dots,group.by = "stim",features = found_genes[1:9],ncol = 3)
# VlnPlot(all_dots,group.by = "plate",features = found_genes[1:9],ncol = 3)
# VlnPlot(all_dots,group.by = "kmeans",features = found_genes[1:9],ncol = 3)
###########
volcano=data.frame(Gene=rownames(all_genes),FDR=all_genes$p_val_adj,LFC=all_genes$avg_logFC)
head(volcano)
sum(is.na(volcano$Pval))
volcano$group="NotSignificant"
volcano[which(volcano$FDR< 0.05),"group"]="Significant"
volcano[which(abs(volcano$LFC) > 0.5 ),"group"]="FoldChange"
volcano[which(volcano$FDR < 0.05 & abs(volcano$LFC) > 0.5 ),"group"]="Significant&FoldChange"
table(volcano$group)
volcano$top=NA
volcano$top[volcano$FDR<1e-10]=volcano$Gene[volcano$FDR<1e-10]
volcano$top[abs(volcano$LFC)<1]=NA

ggplot(volcano,aes(x=LFC,y=-log10(FDR),label=top,color=group))+geom_point(aes(size=-log2(FDR)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("LPS vs control")+xlim(-2,2)

plot_ly(volcano, x=~LFC, y=~-log10(FDR),
        color=~group,colors=cc[c(2,3,1)],
        text=~Gene,type="scatter",size = ~-log10(FDR),
        mode="markers") %>%
  layout(title = 'Volcano',
         yaxis = list(zeroline = FALSE),
         xaxis = list(zeroline = FALSE))  


# Heatmap
mat=(df_all)
marker_genes
pmat=mat[rownames(marker_genes),]
dim(pmat)
pmat=pmat[,all_dots$name]
dim(pmat)
rownames(pmat)
colnames(pmat)
all_dots@meta.data$name

annotation=data.frame(cluster=all_dots@meta.data$kmeans)
rownames(annotation)=all_dots@meta.data$name
head(pmat)
pheatmap(log2(pmat+1),annotation_col = annotation, clustering_distance_rows="euclidean",
         cutree_row = 1, main="DE genes",cluster_rows = F)


pmat=as.matrix(all_dots@assays$SCT@data)
pmat=pmat[rownames(marker_genes),]
dim(pmat)
rownames(pmat)
colnames(pmat)
all_dots@meta.data$name
annotation=data.frame(cluster=all_dots@meta.data$kmeans)
rownames(annotation)=all_dots@meta.data$name
pheatmap(log2(pmat+1),annotation_col = annotation, clustering_distance_rows="euclidean",
         cutree_row = 1, main="DE genes",cluster_rows = F)

leitud=all_genes
leitud=leitud[leitud$avg_logFC>0.5,]
leitud=leitud[leitud$p_val_adj<0.01,]
leitud
write.csv("/Home Office/Thesis/Figures/aall_plate_new.csv",x = leitud)
# koik=rownames(all_dots)
# write.csv("/Home Office/Thesis/Plotting/all_plate_bckgrnd.csv",x = koik)
# GO PLOT
GO=read.csv("/Home Office/Thesis/Figures/all_plate_go_vol2.csv")
go_plot=GO
go_plot$order=nrow(go_plot):1
go_plot$log10fdr=-log10(go_plot$FDR)
library(Hmisc)
go_plot$Description=capitalize(go_plot$Description)  
go_plot$col=rev(go_plot$log10fdr)
ggplot(go_plot,aes(x=order,y=Enrichment, color=log10fdr,size=Genes,label=Description)) +geom_point(stat = "identity") + coord_flip()+
  theme_classic() + scale_color_gradient(low="blue", high="red") + xlab("") + ylab("Enrichment") + 
  scale_x_continuous(label = go_plot$Description, breaks=go_plot$order)+
  theme(axis.text.y = element_text(color="#000000"))+
  ggtitle("TOP 20 GO Hits")+labs(size="Number of Genes",color="-log10(FDR)")+
  ylim(-0.5,max(go_plot$Enrichment))+geom_hline(yintercept = 0,linetype='dashed')
# #  #    #       #            #                   #                               #                                         #
plot_grid(g1,g2)

# MAKE SCE
head(df_all)

read_matrix=as.matrix(df_all)
par(mfrow=c(1,1))
barplot(colSums(read_matrix),border="dark green")
sum(colSums(read_matrix)<500)
sum(colSums(read_matrix)>15000)

read_matrix=read_matrix[,colSums(read_matrix)>500]
read_matrix=read_matrix[,colSums(read_matrix)<15000]

sce=SingleCellExperiment(list(counts=read_matrix))
clusters=quickCluster(sce, method = c("igraph"), min.size = 15)
sce=computeSumFactors(sce, clusters=clusters)
sce=logNormCounts(sce)
stats=modelGeneVar(sce) 
TopHVGs=getTopHVGs(stats,var.field = "bio")

# before and after norm
barplot(colSums(read_matrix),border="dark green")
barplot(colSums(logcounts(sce)),col="dark green",border = "dark green")
# create matrix of normalized expression values
sce$Treatment="Unstim"
sce$Treatment[grep("7|8|9|10|11|12",colnames(sce))]="LPS"
sce$Plate=substr(colnames(sce),1,1)
sce$Name=colnames(sce)

sce=runPCA(sce)
plotPCA(sce)+geom_point(aes(color=sce$Treatment),size=2,alpha=1)+scale_color_manual(values=cc)
plotPCA(sce)+geom_point(aes(color=sce$Plate),size=2,alpha=1)+scale_color_manual(values=cc)
plotPCA(sce)+geom_point(aes(color=sce$Treatment),size=2,alpha=1)+scale_color_manual(values=cc)
library(batchelor)
firstMNN=fastMNN(sce, batch=sce$Plate, subset.row=TopHVGs)
plotReducedDim(firstMNN,dimred="corrected")+ geom_point(aes(color=firstMNN$batch),alpha=1,size=1.5)+
  scale_color_manual(values=cc)+ggtitle("Batchelor - fastMNN - PCA") + labs(color = "Plate")
firstMNN$Treatment=sce$Treatment[match(colnames(firstMNN),colnames(sce))]
plotReducedDim(firstMNN,dimred="corrected")+ geom_point(aes(color=firstMNN$Treatment),alpha=1,size=1.5)+
  scale_color_manual(values=cc)+ggtitle("Batchelor - fastMNN - PCA") + labs(color = "Plate")


library(umap)
norm=reducedDim(firstMNN)
sce$kmeans=as.character(kmeans(norm[,1:2], centers = 2)$cluster)

custom.config = umap.defaults
custom.config$n_neighbors=30
custom.config$min_dist=0.2
custom.config$metric="pearson"
Muuumap=umap(norm,custom.config)
colnames(Muuumap$layout)=c("UC1","UC2")
muumap=cbind(Muuumap$layout,colData(sce))
muumap=as.data.frame(muumap)

ggplot(muumap, aes(x = UC1, y = UC2, color =Plate))+  geom_point(size =1.5,alpha=1)+
  xlab("Umap 1")+ ylab("Umap 2")+  ggtitle("Batchelor - fastMNN - UMAP")+labs(color = "Time")+
  theme_classic()+scale_color_manual(values=cc)
ggplot(muumap, aes(x = UC1, y = UC2, color =Treatment)) +  geom_point(size =1.5,alpha=1)+
  xlab("Umap 1")+ ylab("Umap 2")+  ggtitle("Batchelor - fastMNN - UMAP")+labs(color = "Time")+
  theme_classic()+scale_color_manual(values=cc)
ggplot(muumap, aes(x = UC1, y = UC2, color =kmeans)) +  geom_point(size =1.5,alpha=1)+
  xlab("Umap 1")+ ylab("Umap 2")+  ggtitle("Batchelor - fastMNN - UMAP")+labs(color = "Time")+
  theme_classic()+scale_color_manual(values=cc)

plotPCA(sce)+geom_point(aes(color=sce$kmeans),size=2,alpha=1)+scale_color_manual(values=cc)

########################################################################
sce$UMIs=colSums(assay(sce))
sce$Genes=colSums(assay(sce)>1)
sce$wellKey=colnames(sce)
sca = SceToSingleCellAssay(sce)
# sca$Treatment=factor(sca$Treatment)
# sca$Treatment=relevel(sca$Treatment,"Unstim")

nii=c("Unstim","LPS")
sca$clusters=nii[factor(sca$kmeans)]
sca$clusters=factor(sca$clusters)
sca$clusters=relevel(sca$clusters,"Unstim")
################################################
#
#
#

plotPCA(sca)

# Recalculating cellular detection rate:
cdr2=colSums(assay(sca)>0)
# qplot(x=cdr2, y=colData(sca)$Genes) + xlab('New CDR') + ylab('Old CDR')
colData(sca)$cngeneson=scale(cdr2)
plotPCA(sca)

scaSample=sca[sample(which(freq(sca)>.1), 20),]
flat=as(scaSample, 'data.table')
flat=as.data.frame(flat)
flat=flat[,!duplicated(colnames(flat))]
ggplot(flat, aes(x=logcounts))+geom_density() +facet_wrap(~primerid, scale='free_y')

thres=thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
par(mfrow=c(5,4))
plot(thres)

#####
# zlmCond=zlm(~Treatment + cngeneson, sca)
zlmCond=zlm(~clusters + cngeneson, sca)
# zlm.output=zlm(formula = ~Type+Plate,sca =  sca)
show(zlmCond)
summaryCond=summary(zlmCond, doLRT='clustersLPS')$datatable
# summaryCond=summary(zlmCond, doLRT='TreatmentLPS')$datatable
summaryCond
table(summaryCond$contrast)
#############################################################################################
fcHurdle <- merge(summaryCond[contrast=='clustersLPS' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryCond[contrast=='clustersLPS' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
#############################################################################################
fcHurdle$qval=p.adjust(fcHurdle$`Pr(>Chisq)`,method ="BH", length(fcHurdle$`Pr(>Chisq)`))
# Cutoff - qvalue and LFC
head(fcHurdle)
# fcHurdleSig=fcHurdle[qval<.05 & coef > 0.5]
# fcHSigBcell=fcHurdle[qval<.05 & coef < -0.5]
fcHurdleSig_all=fcHurdle[qval<.05 & coef > (0.5)]
dim(fcHurdleSig_all)
head(fcHurdleSig_all)
setorder(fcHurdleSig_all, qval)
head(fcHurdleSig_all)
#
sig_genes=fcHurdleSig_all$primerid
#
volcano=data.frame(Gene=fcHurdle$primerid,FDR=fcHurdle$qval,LFC=fcHurdle$coef)
head(volcano)
sum(is.na(volcano$Pval))
volcano=volcano[!volcano$LFC=="NaN",]



# volcano$LFC[volcano$LFC=="NaN"]=0
volcano$group="NotSignificant"
volcano[which(volcano$FDR< 0.05),"group"]="Significant"
volcano[which(abs(volcano$LFC) > 0.6 ),"group"]="FoldChange"
volcano[which(volcano$FDR < 0.05 & abs(volcano$LFC) > 0.6 ),"group"]="Significant&FoldChange"
table(volcano$group)

plot_ly(volcano, x=~LFC, y=~-log10(FDR),
        color=~group,colors=cc[c(2,3,1)],
        text=~Gene,type="scatter",size = ~-log10(FDR),
        mode="markers") %>%
  layout(title = 'Volcano',
         yaxis = list(zeroline = FALSE),
         xaxis = list(zeroline = FALSE))  
volcano$text=""
volcano$text[volcano$FDR<1e-24]=volcano$Gene[volcano$FDR<1e-24]
volcano$text[volcano$LFC<(-0.7)]=volcano$Gene[volcano$LFC<(-0.7)]

cc1=c(brewer.pal(name = "Set1",n=9))
cc2=c(brewer.pal(name = "Set1",n=9))

sig=volcano[volcano$group!="NotSignificant",]
unsi=volcano[volcano$group=="NotSignificant",]
unsi2=unsi[sample(1:nrow(unsi),size = 500),]

volcano=rbind(sig,unsi2)

ggplot(volcano,aes(x=LFC,y=-log10(FDR),color=group,label=text))+geom_point()+
  ggrepel::geom_text_repel()+scale_color_manual(values = cc1[c(3,5,1)])+labs(color="")+theme_classic()+
  xlim(-3,3)
  


  
theme(legend.position = "none")
par(mfrow=c(1,1))
display.brewer.all()

# Heatmap
mat=logcounts(sce)
pmat=mat[sig_genes,]
annotation=data.frame(Cluster=sca$clusters)
rownames(annotation)=colnames(sce)
head(pmat)
pheatmap(log2(pmat+1),annotation_col = annotation, clustering_distance_rows="euclidean",
         cutree_row = 1, main="DE genes",cluster_rows = F)

# int_Genes=rownames(pmat)[13:length(rownames(pmat))]
# write.csv("/Home Office/Thesis/Plotting/all_plates_int_vol1.csv",x = sig_genes)
# koik=rownames(sca)
# write.csv("/Home Office/Thesis/Plotting/all_plate_bckgrnd.csv_vol2",x = koik)
#############
bulk=read.csv("/Home Office/bulksig.csv")[,2]
sc=read.csv("/Home Office/Thesis/Plotting/all_plates_int_vol1.csv")[,2]
eulerplot=data.frame(Genes=unique(c(bulk,sc)),Bulk=FALSE,SC=FALSE)
eulerplot$Bulk[eulerplot$Genes%in%bulk]=TRUE
eulerplot$SC[eulerplot$Genes%in%sc]=TRUE
cc=c(brewer.pal(name = "Set3",n=8))
plot(eulerr::euler(eulerplot[, 2:3]),col=cc[2],fill=cc[c(3,1,4)],labels = list(font = 4,col="black"),
     quantities = list(font = 4,col="black"),main=list(col="black", label="Overlapping DE genes"))

# 
# GO1=read.csv("/Home Office/Thesis/Plotting/all_plate_go_vol1.csv")
# GO2=read.csv("/Home Office/Thesis/Plotting/bulk go.csv")
# sum(GO1$GO.term%in%GO2$GO.term)
# sum(GO1$Description%in%GO2$Description)
# 
# GO3=read.csv("/Home Office/Thesis/Plotting/comp GO.csv")
# sum(GO3$SC[1:100]%in%GO3$Bulk[1:100])


# GO PLOT
GO=read.csv("/Home Office/Thesis/Plotting/all_plate_go_vol1.csv")
GO$log10fdr=-log10(GO$FDR.q.value)
head(GO,10)
GO=GO[order(GO$log10fdr,decreasing = T),]
GO$order=nrow(GO):1
GO$Description=capitalize(GO$Description)  
GO$col=rev(GO$log10fdr)
head(GO,10)


ggplot(GO,aes(x=order,y=log10fdr, fill=log10fdr)) +geom_bar(stat = "identity") + coord_flip() + 
  theme_classic() + scale_fill_continuous(low="blue", high="red") + xlab("") + ylab("-log10(FDR)") + 
  scale_x_continuous(label = GO$Description, breaks=GO$order)+
  geom_label(label=GO$GO.term,aes(y = 1.5,color=col),label.padding = unit(0.15, "lines"))+
  theme(legend.position = "none", axis.text.y = element_text(face="bold", color="#000000"))+
  ggtitle("Top 20 GO hits all plates together")

# NEW GO
GO=read.csv("/Home Office/Thesis/Plotting/all plates new go.csv")
head(GO)
GO$log10fdr=-log10(GO$FDR.q.value)
head(GO,10)
GO=GO[order(GO$log10fdr,decreasing = T),]
GO$order=nrow(GO):1
GO$Description=capitalize(GO$Description)  
GO$col=rev(GO$log10fdr)
head(GO,10)
GO=GO[1:11,]

ggplot(GO,aes(x=order,y=Enrichment, color=log10fdr,size=no.of.hits)) +geom_point(stat = "identity") + coord_flip()+
  theme_classic() + scale_color_continuous(low="#5D00A9", high="red") + xlab("") + ylab("Enrichment") + 
  scale_x_continuous(label = GO$Description, breaks=GO$order)+
  theme(axis.text.y = element_text(color="#000000"))+
  ggtitle("Top 20 pathway enrichment - GOrilla")+labs(size="Number of Genes",color="FDR")
