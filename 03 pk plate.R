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
library(ggsci)
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)
library(cowplot)
library(biomaRt)
options(stringsAsFactors = F)
cc=c(brewer.pal(name = "Set1",n=9))

# Read in ribo scrb
rds_file="/Home Office/Data/2017-18_Early_seq_results/04 PK_Test 2018/PK_Test.dgecounts.rds"
rds_here=readRDS(rds_file)
reads=rds_here$exons$downsampled$downsampled_37614$umicounts_downsampled
reads=as.data.frame(as.matrix(reads))
summary(colSums(reads)) # UMI
summary(colSums(reads>0)) # Gene
# Top expr genes
sort(rowMeans(reads),decreasing = T)[1:50]
keenid=names(sort(rowMeans(reads),decreasing = T)[1:50])
genes=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv")
genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
genes[match(keenid,genes$ensembl_gene_id),"description"]
# Top var genes
keenid2=names(sort(apply(reads, 1,var),decreasing = T)[1:50])
genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
#
summary(colSums(reads))
new_bcs=read.csv("/Home Office/PhD/Experiments/SCRP-Seq_Barcodes/My_1st_set/1st_set.csv",header = F)
lugendid=reads
sum(colnames(lugendid)%in%new_bcs$V3)
colnames(lugendid)=new_bcs[match(colnames(lugendid),new_bcs$V3),"V1"]
lugendid=lugendid[,!is.na(colnames(lugendid))]
lugendid=lugendid[,mixedorder(colnames(lugendid))]
dim(lugendid)
colnames(lugendid)
# Renaming genes
Gene_data=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv",col.names =c("gene_id","gene_name","description") )
Gene_data=Gene_data[Gene_data$gene_id%in%rownames(lugendid),]
paste("You have ", sum(duplicated(Gene_data$gene_name)),"duplicated gene names in your data")
Gene_data=Gene_data[!duplicated(Gene_data$gene_name),]
lugendid=lugendid[rownames(lugendid)%in%Gene_data$gene_id,]
paste("You have ", sum(duplicated(Gene_data$gene_name)),"duplicated gene names in your data")
# Order
rownames(lugendid)=Gene_data[match(rownames(lugendid),Gene_data$gene_id),"gene_name"]
PK_plate=lugendid
PK_plate[1:5,1:5]
# Seurat
colSums(PK_plate)
barplot(colSums(PK_plate))
summary(colSums(PK_plate))


PK_dots=CreateSeuratObject(counts = PK_plate, project = "ribo plate")
#
PK_dots[["percent.mt"]]=PercentageFeatureSet(PK_dots, pattern = "^MT-")
pk_violin=VlnPlot(PK_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pk_violin
Mito_cut=10
RNA_cut=1000
PK_dots=subset(PK_dots, subset = nFeature_RNA > RNA_cut & percent.mt< Mito_cut)
VlnPlot(PK_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(PK_dots, feature1 = "nFeature_RNA", feature2 = "percent.mt")
dim(PK_dots)
# Normalizing with SCtransform
PK_dots=SCTransform(PK_dots)
PK_dots$stim="Unstim."
PK_dots$well=rownames(PK_dots@meta.data)
PK_dots$stim[grepl("7|8|9|10|11|12",PK_dots$well)]="LPS"
# PCA
PK_dots=RunPCA(PK_dots, npcs = 30, verbose = FALSE)
pk_p=DimPlot(PK_dots,reduction = "pca",group.by = "stim",pt.size = 2)+ggtitle("PK Test")
pk_p
# Umap
PK_dots=RunUMAP(PK_dots, dims = 1:15, n.neighbors = 25,umap.method = "umap-learn", 
                  metric = "correlation",min.dist = 0.04,verbose = FALSE)

DimPlot(PK_dots,reduction = "umap",group.by = "stim",pt.size = 2)+ggtitle("PK. test")+
  xlab("")+ylab("")+scale_color_manual(values = cc)+xlab("")+ylab("")

# PK_dots=FindNeighbors(PK_dots)
# PK_dots=FindClusters(PK_dots)
# DimPlot(PK_dots,reduction = "umap",pt.size = 2) + xlab("")+ylab("")+scale_color_manual(values = cc)

PK_dots$kmeans=as.character(kmeans(PK_dots@reductions$umap@cell.embeddings, centers = 2)$cluster)
cc=c(brewer.pal(name = "Dark2",n=8))

u3=DimPlot(PK_dots,reduction = "umap",pt.size = 2,group.by = "stim")+
  theme(legend.position = "none")+ scale_color_manual(values = cc)+xlab("")+ylab("")
u3.2=DimPlot(PK_dots,reduction = "umap",pt.size = 2,group.by = "kmeans")+
  theme(legend.position = "none")+ scale_color_manual(values = cc)+xlab("")+ylab("")

u3b=DimPlot(PK_dots,reduction = "umap",pt.size = 2,group.by = "stim")+
  scale_color_manual(values = cc)+ggtitle("ribo")
u3.2b=DimPlot(PK_dots,reduction = "umap",pt.size = 2,group.by = "kmeans")+
  scale_color_manual(values = cc)

# sort(rowMeans(PK_dots),decreasing = T)[1:30]
# pk_violin
# pk_p
# pk_p$data[which(pk_p$data$PC_2==max(pk_p$data$PC_2)),]


# MAKE SCE
read_matrix=as.matrix(PK_plate)
par(mfrow=c(1,1))
barplot(colSums(read_matrix))
sum(colSums(read_matrix)<1000)
read_matrix=read_matrix[,colSums(read_matrix)>1000]
read_matrix=read_matrix[,colSums(read_matrix)<6000]
sce=SingleCellExperiment(list(counts=read_matrix))
clusters=quickCluster(sce, method = c("igraph"), min.size = 5)
sce=computeSumFactors(sce, clusters=clusters)
sce=logNormCounts(sce)
stats=modelGeneVar(sce) 
TopHVGs=getTopHVGs(stats,var.field = "bio")
# before and after norm
# barplot(colSums(first_working_plate),border="dark green")
# barplot(colSums(logcounts(sce)),col="dark green",border = "dark green")
# create matrix of normalized expression values
sce$Treatment="Unstim"
sce$Treatment[grep("7|8|9|10|11|12",colnames(sce))]="LPS"
# info=as.data.frame(first_working_dots@meta.data)
# sum(!colnames(sce)%in%rownames(info))
# sce=sce[,rownames(info)]
# sce$clusters=info[match(colnames(sce),info$well),"seurat_clusters"]
sce=runPCA(sce)
plotPCA(sce)+geom_point(aes(color=sce$Treatment),size=2,alpha=1)+scale_color_manual(values=cc)
sce_pca=reducedDim(sce)
sce$kmeans=as.character(kmeans(sce_pca[,1:2], centers = 2)$cluster)
plotPCA(sce)+geom_point(aes(color=sce$kmeans),size=2,alpha=1)+scale_color_manual(values=cc)

sce=runUMAP(sce)
# sce_umap=reducedDims(sce)$UMAP
# sce$kmeans_umap=as.character(kmeans(sce_umap, centers = 2)$cluster)
plotUMAP(sce)+geom_point(aes(color=sce$Treatment),size=2,alpha=1)+scale_color_manual(values=cc)
plotUMAP(sce)+geom_point(aes(color=sce$kmeans),size=2,alpha=1)+scale_color_manual(values=cc)

########################################################################
sce$UMIs=colSums(assay(sce))
sce$Genes=colSums(assay(sce)>1)
sce$wellKey=colnames(sce)
sca = SceToSingleCellAssay(sce)
sca$Treatment=factor(sca$Treatment)
sca$Treatment=relevel(sca$Treatment,"Unstim")

nii=c("LPS","Unstim")
sca$clusters=nii[factor(sca$kmeans)]
sca$clusters=factor(sca$clusters)
sca$clusters=relevel(sca$clusters,"Unstim")
################################################
plotPCA <- function(sca_obj){
  projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
  colnames(projection)=c("PC1","PC2","PC3","PC4")
  pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
  print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'UMIs','Genes'),
                mapping=aes(color=Treatment), upper=list(continuous='blank')))
  invisible(pca)
}
plotPCA(sca)

# Recalculating cellular detection rate:
cdr2=colSums(assay(sca)>0)
qplot(x=cdr2, y=colData(sca)$Genes) + xlab('New CDR') + ylab('Old CDR')
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
volcano[which(abs(volcano$LFC) > 1 ),"group"]="FoldChange"
volcano[which(volcano$FDR < 0.05 & abs(volcano$LFC) > 1 ),"group"]="Significant&FoldChange"
table(volcano$group)

# plot_ly(volcano, x=~LFC, y=~-log10(FDR),
#         color=~group,colors=cc[c(2,3,1)],
#         text=~Gene,type="scatter",size = ~-log10(FDR),
#         mode="markers") %>%
#   layout(title = 'Volcano',
#          yaxis = list(zeroline = FALSE),
#          xaxis = list(zeroline = FALSE))  

volcano$text=""
head(volcano)
volcano$logfdr=-log10(volcano$FDR)
# volcano$text[abs(volcano$LFC)>(1.5)]=volcano$Gene[abs(volcano$LFC)>(1.5)]
# volcano$text[volcano$logfdr<5]=""
# volcano$Gene[volcano$logfdr>5]
# volcano[order(volcano$LFC),]$text[1:10]=volcano[order(volcano$LFC),]$Gene[1:10]
volcano[order(volcano$logfdr,decreasing = T),]$text[1:10]=volcano[order(volcano$logfdr,decreasing = T),]$Gene[1:10]

cc1=c(brewer.pal(name = "Set1",n=9))

sig=volcano[volcano$group!="NotSignificant",]
unsi=volcano[volcano$group=="NotSignificant",]
unsi2=unsi[sample(1:nrow(unsi),size = 1000),]
volcano=rbind(sig,unsi2)

ggplot(volcano,aes(x=LFC,y=-log10(FDR),color=group,label=text))+geom_point()+
  ggrepel::geom_text_repel()+scale_color_manual(values = cc1[c(3,5,1)])+labs(color="")+theme_classic()+
  xlim(-3,4)

#

#

#

#

# Heatmap
mat=logcounts(sce)
pmat=mat[sig_genes,]
annotation=data.frame(Cluster=sca$clusters)
rownames(annotation)=colnames(sce)
head(pmat)
pheatmap(log2(pmat+1),annotation_col = annotation, clustering_distance_rows="euclidean",
         cutree_row = 1, main="DE genes",cluster_rows = F)

# int_Genes=rownames(pmat)[13:length(rownames(pmat))]
int_Genes=rownames(pmat)
int_Genes
write.csv("/Home Office/Thesis/Plotting/pk_plate.csv",x = int_Genes)
koik=rownames(sce)
write.csv("/Home Office/Thesis/Plotting/pk_plate_bckgrnd.csv",x = koik)

# GO PLOT
GO=read.csv("/Home Office/Thesis/Plotting/PK_plate GO results.csv")
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
  ggtitle("Top 20 GO hits first working plate")


go_essa=read.csv("/Home Office/Thesis/Plotting/essa plate Go results.csv")
go_ribbo=read.csv("/Home Office/Thesis/Plotting/ribo_plate GO results.csv")
go_pk=read.csv("/Home Office/Thesis/Plotting/PK_plate GO results.csv")

Reduce(intersect, list(go_essa$GO.term,go_ribbo$GO.term,go_pk$GO.term))
Reduce(intersect, list(go_essa$Description,go_ribbo$Description,go_pk$Description))


