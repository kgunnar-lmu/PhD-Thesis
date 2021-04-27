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
# First plate ever

rds_file="/Home Office/Data/2017-18_Early_seq_results/02 Scrb working 2017/zUMI_scrb_seq.dgecounts.rds"
rds_here=readRDS(rds_file)
reads=rds_here$exons$downsampled$downsampled_118012$umicounts_downsampled
reads=as.data.frame(as.matrix(reads))
summary(colSums(reads)) # Per UMI
summary(colSums(reads>0)) # Per Gene
# Top expr genes
sort(rowMeans(reads),decreasing = T)[1:50]
keenid=names(sort(rowMeans(reads),decreasing = T)[1:50])
genes=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv")
genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
genes[match(keenid,genes$ensembl_gene_id),"description"]
# Top var genes
keenid2=names(sort(apply(reads, 1,var),decreasing = T)[1:50])
genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
# plots
barplot(colSums(reads))
summary(colSums(reads))
dim(reads)
lugendid=reads
colnames(lugendid)
new_bcs=read.csv("/Home Office/PhD/Experiments/SCRP-Seq_Barcodes/My_1st_set/1st_set.csv",header = F)
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
first_working_plate=lugendid
first_working_plate[1:5,1:5]

first_working_dots=CreateSeuratObject(counts = first_working_plate, project = "first working")
#
first_working_dots[["percent.mt"]]=PercentageFeatureSet(first_working_dots, pattern = "^MT-")
VlnPlot(first_working_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Mito_cut=9
max_RNA=3000
first_working_dots=subset(first_working_dots, subset = nFeature_RNA < max_RNA & percent.mt < Mito_cut)
VlnPlot(first_working_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dim(first_working_dots)
#
first_working_dots$stim="Unstim."
first_working_dots$well=rownames(first_working_dots@meta.data)
first_working_dots$stim[grepl("7|8|9|10|11|12",first_working_dots$well)]="LPS"
###
first_working_dots=SCTransform(first_working_dots)
first_working_dots=RunPCA(first_working_dots, npcs = 30, verbose = FALSE)
DimPlot(first_working_dots,reduction = "pca",group.by = "stim",pt.size = 2)+ggtitle("SCRB-seq 2017")+xlab("")+ylab("")

first_working_dots=RunUMAP(first_working_dots, dims = 1:30, n.neighbors = 15,umap.method = "umap-learn", 
                           metric = "correlation",min.dist = 0.5,verbose = FALSE)# Read in first working mcscrb

cc=c(brewer.pal(name = "Set1",n=9))

DimPlot(first_working_dots,reduction = "umap",group.by = "stim",pt.size = 2)+ggtitle("SCRB-seq 2017")+
  theme(legend.position = "none")+xlab("")+ylab("")+scale_color_manual(values = cc)

# first_working_dots=FindNeighbors(first_working_dots)
# first_working_dots=FindClusters(first_working_dots)
# DimPlot(first_working_dots,reduction = "umap",pt.size = 2)+
#   theme(legend.position = "none")+scale_color_manual(values = cc)

first_working_dots$kmeans=as.character(kmeans(first_working_dots@reductions$umap@cell.embeddings, centers = 2)$cluster)
cc=c(brewer.pal(name = "Dark2",n=8))

u1=DimPlot(first_working_dots,reduction = "umap",pt.size = 2,group.by = "stim")+
  theme(legend.position = "none")+ scale_color_manual(values = cc)+xlab("")+ylab("")
u1.2=DimPlot(first_working_dots,reduction = "umap",pt.size = 2,group.by = "kmeans")+
  theme(legend.position = "none")+ scale_color_manual(values = cc)+xlab("")+ylab("")

u1b=DimPlot(first_working_dots,reduction = "umap",pt.size = 2,group.by = "stim")+
  scale_color_manual(values = cc)+ggtitle("first ever")
u1.2b=DimPlot(first_working_dots,reduction = "umap",pt.size = 2,group.by = "kmeans")+
  scale_color_manual(values = cc)

"#####################################################################################################################"
"#####################################################################################################################"
"#####################################################################################################################"
'Idents(first_working_dots)=first_working_dots$kmeans
first_genes=FindMarkers(object=first_working_dots,ident.1=1)
# Keep only up genes in middle and FDR 0.05% or less
marker_genes=first_genes[first_genes$p_val_adj<0.05,]
marker_genes=marker_genes[marker_genes$avg_logFC>0,]
marker_genes=marker_genes[order(marker_genes$p_val_adj),]
dim(marker_genes)
head(marker_genes)
found_genes=rownames(marker_genes)
VlnPlot(first_working_dots,group.by = "stim",features = found_genes[1:9],ncol = 3)
VlnPlot(first_working_dots,group.by = "kmeans",features = found_genes[1:9],ncol = 3)
###########
volcano=data.frame(Gene=rownames(first_genes),FDR=first_genes$p_val_adj,LFC=first_genes$avg_logFC)
head(volcano)
sum(is.na(volcano$Pval))
volcano$group="NotSignificant"
volcano[which(volcano$FDR< 0.05),"group"]="Significant"
volcano[which(abs(volcano$LFC) > 0.5 ),"group"]="FoldChange"
volcano[which(volcano$FDR < 0.05 & abs(volcano$LFC) > 0.5 ),"group"]="Significant&FoldChange"
table(volcano$group)

plot_ly(volcano, x=~LFC, y=~-log10(FDR),
        color=~group,colors=cc[c(2,3,1)],
        text=~Gene,type="scatter",size = ~-log10(FDR),
        mode="markers") %>%
  layout(title = "Volcano",
         yaxis = list(zeroline = FALSE),
         xaxis = list(zeroline = FALSE))  

# Heatmap
mat=(first_working_plate)
marker_genes
pmat=mat[rownames(marker_genes),]
dim(pmat)
pmat=pmat[,colnames(first_working_dots)]
dim(pmat)
rownames(pmat)
colnames(pmat)
annotation=data.frame(cluster=first_working_dots@meta.data$stim)
rownames(annotation)=colnames(first_working_dots)
head(pmat)
pheatmap(log2(pmat+1),annotation_col = annotation, clustering_distance_rows="euclidean",
         cutree_row = 1, main="DE genes",cluster_rows = F)


nii=read.csv("/Home Office/Thesis/Plotting/essa_plate.csv")
length(rownames(marker_genes))
length(nii$x)
sum(nii$x%in%rownames(marker_genes))
write.csv(rownames(marker_genes),"/Home Office/Thesis/Plotting/essa_plate_vol2.csv")
write.csv(rownames(first_working_dots),"/Home Office/Thesis/Plotting/essa_plate_bckgrnd_vol2.csv")'
"#####################################################################################################################"
"#####################################################################################################################"
"#####################################################################################################################"

# MAKE SCE
read_matrix=as.matrix(first_working_plate)
sum(colSums(read_matrix)<1000)
read_matrix=read_matrix[,colSums(read_matrix)>1000]
read_matrix=read_matrix[,colSums(read_matrix)<15000]
sce=SingleCellExperiment(list(counts=read_matrix))
clusters=quickCluster(sce, method = c("igraph"), min.size = 5)
sce=computeSumFactors(sce, clusters=clusters)
sce=logNormCounts(sce)
stats=modelGeneVar(sce) 
TopHVGs=getTopHVGs(stats,var.field = "bio")
# before and after norm
barplot(colSums(first_working_plate),border="dark green")
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

nii=c("Unstim","LPS")
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
fcHurdleSig_all=fcHurdle[qval<.05 & coef > 0.5]
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

#### GGPLOT Volcano
volcano$text=""
head(volcano)
volcano$logfdr=-log10(volcano$FDR)
# volcano$text[abs(volcano$LFC)>(1.5)]=volcano$Gene[abs(volcano$LFC)>(1.5)]
# volcano$text[volcano$logfdr<5]=""
# volcano$Gene[volcano$logfdr>5]
volcano[order(volcano$LFC),]$text[1:10]=volcano[order(volcano$LFC),]$Gene[1:10]
volcano[order(volcano$logfdr,decreasing = T),]$text[1:10]=volcano[order(volcano$logfdr,decreasing = T),]$Gene[1:10]

cc1=c(brewer.pal(name = "Set1",n=9))

sig=volcano[volcano$group!="NotSignificant",]
unsi=volcano[volcano$group=="NotSignificant",]
unsi2=unsi[sample(1:nrow(unsi),size = 4000),]
volcano=rbind(sig,unsi2)

ggplot(volcano,aes(x=LFC,y=-log10(FDR),color=group,label=text))+geom_point()+
  ggrepel::geom_text_repel()+scale_color_manual(values = cc1[c(2,3,5,1)])+labs(color="")+theme_classic()+
  xlim(-3,4)
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

int_Genes=rownames(pmat)[13:length(rownames(pmat))]
int_Genes=rownames(pmat)
# write.csv("/Home Office/Thesis/Plotting/essa_plate.csv",x = int_Genes)
koik=rownames(sce)
# write.csv("/Home Office/Thesis/Plotting/essa_plate_bckgrnd.csv",x = koik)

# GO PLOT
GO=read.csv("/Home Office/Thesis/Plotting/essa plate Go results.csv")
GO$order=nrow(GO):1
GO$log10fdr=-log10(GO$FDR.q.value)
GO$Description=capitalize(GO$Description)  
GO$col=rev(GO$log10fdr)
head(GO,10)

ggplot(GO,aes(x=order,y=log10fdr, fill=log10fdr)) +geom_bar(stat = "identity") + coord_flip() + 
  theme_classic() + scale_fill_continuous(low="blue", high="red") + xlab("") + ylab("-log10(FDR)") + 
  scale_x_continuous(label = GO$Description, breaks=GO$order)+
  geom_label(label=GO$GO.term,aes(y = 1.5,color=col),label.padding = unit(0.15, "lines"))+
  theme(legend.position = "none", axis.text.y = element_text(face="bold", color="#000000"))+
  ggtitle("Top 20 GO hits first working plate")





