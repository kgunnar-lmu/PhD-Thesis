# Generate stats for further filtering and analysis. 
# Filtering is the next script
library(tidyverse)
library(stringr)
library(Seurat)
library(sctransform)
library(scran)
library(scater)
library(Rtsne)
library(rsvd)
library(reshape2)
library(plotly)
library(Hmisc)
library(RColorBrewer)
library(pheatmap)
library(NMF)
library(MAST)
library(limma)
library(knitr)
library(gtools)
# library(GSEABase)
library(ggsci)
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)
library(cowplot)
library(biomaRt)
options(stringsAsFactors = F)
cc=c(brewer.pal(name = "Set1",n=9))


# Read in first working scrb
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
p1=DimPlot(first_working_dots,reduction = "pca",group.by = "stim",pt.size = 2)+ggtitle("SCRB-seq 2017")+xlab("")+ylab("")

first_working_dots=RunUMAP(first_working_dots, dims = 1:30, n.neighbors = 15,umap.method = "umap-learn", 
                        metric = "correlation",min.dist = 0.5,verbose = FALSE)# Read in first working mcscrb
u1=DimPlot(first_working_dots,reduction = "umap",group.by = "stim",pt.size = 2)+ggtitle("SCRB-seq 2017")+
  theme(legend.position = "none")+xlab("")+ylab("")+scale_color_manual(values = cc)
dev.off()
# display.brewer.all()
first_working_dots=FindNeighbors(first_working_dots)
first_working_dots=FindClusters(first_working_dots)
u1.2=DimPlot(first_working_dots,reduction = "umap",pt.size = 2)+
  theme(legend.position = "none")+scale_color_manual(values = cc)

DimPlot(first_working_dots,reduction = "umap",pt.size = 2)+
  scale_color_manual(values = cc)



#################################################################################################################
#################################################################################################################
#################################################################################################################
# Read in new mcscrb
rds_file="/Home Office/Data/2017-18_Early_seq_results/03 mcscrb 2017-18/mc_scrb_seq_2018.dgecounts.rds"
rds_here=readRDS(rds_file)
reads=rds_here$exons$downsampled$downsampled_241164$umicounts_downsampled
reads=as.data.frame(as.matrix(reads))
summary(colSums(reads)) # UMI
summary(colSums(reads>0)) # Gene
# Top expr genes
sort(rowMeans(reads),decreasing = T)[1:50]
keenid=names(sort(rowMeans(reads),decreasing = T)[1:50])
# genes=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv")
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
# new_bcs=read.csv("/Home Office/PhD/Experiments/SCRP-Seq_Barcodes/My_1st_set/1st_set.csv",header = F)
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
mcscrb_plate=lugendid
mcscrb_plate[1:5,1:5]
mcscrb_dots=CreateSeuratObject(counts=mcscrb_plate, project = "mcscrb plate")
#
mcscrb_dots[["percent.mt"]]=PercentageFeatureSet(mcscrb_dots, pattern = "^MT-")
VlnPlot(mcscrb_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Mito_cut=6
max_RNA=3500
mcscrb_dots=subset(mcscrb_dots, subset = nFeature_RNA < max_RNA & percent.mt < Mito_cut)
VlnPlot(mcscrb_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(mcscrb_dots, feature1 = "nFeature_RNA", feature2 = "percent.mt")
dim(mcscrb_dots)
#
mcscrb_dots$stim="Unstim."
mcscrb_dots$well=rownames(mcscrb_dots@meta.data)
mcscrb_dots$stim[grepl("7|8|9|10|11|12",mcscrb_dots$well)]="LPS"
#
mcscrb_dots=SCTransform(mcscrb_dots)
mcscrb_dots=RunPCA(mcscrb_dots, npcs = 30, verbose = FALSE)
DimPlot(mcscrb_dots,reduction = "pca",group.by = "stim",pt.size = 2)+ggtitle("mcSCRB 2018")
mcscrb_dots=RunUMAP(mcscrb_dots, dims = 1:30, n.neighbors = 15,umap.method = "umap-learn", 
                        metric = "correlation",min.dist = 0.2,verbose = FALSE)

cc=c(brewer.pal(name = "Dark2",n=8))

mcscrb_dots$kmeans=as.character(kmeans(mcscrb_dots@reductions$umap@cell.embeddings, centers = 2)$cluster)

u4=DimPlot(mcscrb_dots,reduction = "umap",group.by = "stim",pt.size = 2)+scale_color_manual(values = cc)+
  theme(legend.position = "none")+xlab("")+ylab("")
u4.2=DimPlot(mcscrb_dots,reduction = "umap",group.by = "kmeans",pt.size = 2)+scale_color_manual(values = cc)+
  theme(legend.position = "none")+xlab("")+ylab("")

u4b=DimPlot(mcscrb_dots,reduction = "umap",group.by = "stim",pt.size = 2)+ggtitle("mcSCRB 2018")+scale_color_manual(values = cc)
u4.2b=DimPlot(mcscrb_dots,reduction = "umap",group.by = "kmeans",pt.size = 2)+scale_color_manual(values = cc)



##############################################################################################################
##############################################################################################################
##############################################################################################################
# Read in ribo scrb
rds_file="/Home Office/Data/2017-18_Early_seq_results/04 Ribo_test 2018/Ribo_test.dgecounts.rds"
rds_here=readRDS(rds_file)
reads=rds_here$exons$downsampled$downsampled_130445$umicounts_downsampled
reads=as.data.frame(as.matrix(reads))
summary(colSums(reads)) # UMI
summary(colSums(reads>0)) # Gene
# Top expr genes
sort(rowMeans(reads),decreasing = T)[1:50]
keenid=names(sort(rowMeans(reads),decreasing = T)[1:50])
# genes=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv")
genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
genes[match(keenid,genes$ensembl_gene_id),"description"]
# Top var genes
keenid2=names(sort(apply(reads, 1,var),decreasing = T)[1:50])
genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
# plots
barplot(colSums(reads))
summary(colSums(reads))
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
ribo_plate=lugendid
ribo_plate[1:5,1:5]
# Seurat
ribo_dots=CreateSeuratObject(counts = ribo_plate, project = "ribo plate")
#
ribo_dots[["percent.mt"]]=PercentageFeatureSet(ribo_dots, pattern = "^MT-")
VlnPlot(ribo_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Mito_cut=10
RNA_cut=1000
ribo_dots=subset(ribo_dots, subset = nFeature_RNA > RNA_cut & percent.mt< Mito_cut)
VlnPlot(ribo_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(ribo_dots, feature1 = "nFeature_RNA", feature2 = "percent.mt")
dim(ribo_dots)
# Normalizing with SCtransform
ribo_dots=SCTransform(ribo_dots)
ribo_dots$stim="Unstim."
ribo_dots$well=rownames(ribo_dots@meta.data)
ribo_dots$stim[grepl("7|8|9|10|11|12",ribo_dots$well)]="LPS"
# PCA
ribo_dots=RunPCA(ribo_dots, npcs = 30, verbose = FALSE)
p2=DimPlot(ribo_dots,reduction = "pca",group.by = "stim",pt.size = 2)+ggtitle("Ribo Plate")

# Umap
ribo_dots=RunUMAP(ribo_dots, dims = 1:15, n.neighbors = 25,umap.method = "umap-learn", 
                          metric = "correlation",min.dist = 0.04,verbose = FALSE)
u2=DimPlot(ribo_dots,reduction = "umap",group.by = "stim",pt.size = 2)+ggtitle("Ribonuk. test")+
  theme(legend.position = "none")+xlab("")+ylab("")+scale_color_manual(values = cc)+xlab("")+ylab("")+scale_color_manual(values = cc)
ribo_dots=FindNeighbors(ribo_dots)
ribo_dots=FindClusters(ribo_dots)
u2.2=DimPlot(ribo_dots,reduction = "umap",pt.size = 2)+
  theme(legend.position = "none")+xlab("")+ylab("")+scale_color_manual(values = cc)

###################################################################################################
###################################################################################################
# Read in prot.k scrb
rds_file="/Home Office/Data/2017-18_Early_seq_results/04 PK_Test 2018/PK_Test.dgecounts.rds"
rds_here=readRDS(rds_file)
reads=rds_here$exons$downsampled$downsampled_37614$umicounts_downsampled
reads=as.data.frame(as.matrix(reads))
summary(colSums(reads)) # UMIs
summary(colSums(reads>0)) # Genes
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
PK_plate=lugendid
PK_plate[1:5,1:5]
# Seurat
# PK_dots=CreateSeuratObject(counts = PK_plate, project = "ribo plate",min.cells = 10,min.features = 200)
PK_dots=CreateSeuratObject(counts = PK_plate, project = "PK plate")
#
PK_dots[["percent.mt"]]=PercentageFeatureSet(PK_dots, pattern = "^MT-")
VlnPlot(PK_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Mito_cut=6
RNA_cut=600
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
p3=DimPlot(PK_dots,reduction = "pca",group.by = "stim",pt.size = 2)+ggtitle("Prot.K. test")
# Umap
PK_dots=RunUMAP(PK_dots, dims = 1:30, n.neighbors = 25,umap.method = "umap-learn", 
                  metric = "correlation",min.dist = 0.04,verbose = FALSE)
u3=DimPlot(PK_dots,reduction = "umap",group.by = "stim",pt.size = 2)+ggtitle("Prot.K. test")+xlab("")+ylab("")+
  scale_color_manual(values = cc)+labs(color = "Stimulus")
PK_dots=FindNeighbors(PK_dots)
PK_dots=FindClusters(PK_dots)
u3.2=DimPlot(PK_dots,reduction = "umap",pt.size = 2)+xlab("")+ylab("")+scale_color_manual(values = cc)+labs(color = "Clusters")

###################################################################################################
plot_grid(p1,p2,p3,ncol = 3)

plot_grid(u1,u2,u3,u1.2,u2.2,u3.2,ncol = 3)
# dim(u1$data)
# dim(u2$data)
# dim(u3$data)

###################################################################################################
first_working_plate
read_matrix=as.matrix(first_working_plate)
read_matrix=read_matrix[,which(colSums(read_matrix>0)>500)]
# read_matrix=read_matrix[rowMeans(read_matrix)>0.05,]
sce=SingleCellExperiment(list(counts=read_matrix))
# sce <- calculateQCMetrics(sce)
clusters=quickCluster(sce, method = c("igraph"), min.size = 5)
sce=computeSumFactors(sce, clusters=clusters)
sce=logNormCounts(sce)
stats=modelGeneVar(sce) 
TopHVGs=getTopHVGs(stats,var.field = "bio")
barplot(colSums(first_working_plate),border="dark green")
barplot(colSums(logcounts(sce)),col="dark green",border = "dark green")
# create matrix of normalized expression values
sce$Treatment="Unstim"
sce$Treatment[grep("7|8|9|10|11|12",colnames(sce))]="LPS"

cc=c(brewer.pal(name = "Set1",n=9))

info=as.data.frame(first_working_dots@meta.data)
sum(!colnames(sce)%in%rownames(info))
sce=sce[,rownames(info)]
sce$clusters=info[match(colnames(sce),info$well),"seurat_clusters"]

sce=runPCA(sce)
plotPCA(sce)+geom_point(aes(color=sce$Treatment),size=2,alpha=1)+scale_color_manual(values=cc)
plotPCA(sce)+geom_point(aes(color=sce$clusters),size=2,alpha=1)+scale_color_manual(values=cc)

sce=runUMAP(sce)
plotUMAP(sce)+geom_point(aes(color=sce$Treatment),size=2,alpha=1)+scale_color_manual(values=cc)
plotUMAP(sce)+geom_point(aes(color=sce$clusters),size=2,alpha=1)+scale_color_manual(values=cc)

"Heatmap"

pmat=mat[,1:100]
pheatmap(pmat,annotation_row = annotation, clustering_distance_rows="euclidean",
         cutree_row = 1, main="clustering top 100 genes")
############### 
"t-SNE Plot"
############### 
## Executing the algorithm on curated data and the top 500 most variable genes
tsne15 <- Rtsne(mat, perplexity=15, verbose = TRUE)
data.frame(tsne15$Y, batch = pcs$batch,sample=pcs$labels)
ggplot(df, aes(x = X1, y = X2, col = batch,label=sample)) + geom_point(size=0.1,alpha=0.1) +
  ggtitle("perplexity 15")+  geom_text(hjust = 1.5,position = position_jitter(0.5,0.5))

" Mean vs CV"
counts_matrix=round(as.matrix(scrb.scran))
mean_per_gene=rowMeans(scrb.scran)
njaa=as.data.frame(mean_per_gene)
ens_genes=rownames(scrb.scran)
mean_per_gene=njaa$mean_per_gene
counts_matrix=as.matrix(scrb.scran)
SD=transform(counts_matrix, SD=rowSds(counts_matrix, na.rm=F))$SD
CV=SD/mean_per_gene
plotting_df=data.frame(names=ens_genes,mean=mean_per_gene,CV2=CV^2)
head(plotting_df)
#plotting
plot(x=log10(plotting_df$mean),y=plotting_df$CV2,pch=1,cex=1,col="grey50",lwd=0.5,
     ylab = "average squared coeffcient of variations across samples",
     xlab = "gene log10 mean across samples",
     main = "nii")
########################
"Mean vs drop-out"
sum(scrb.scran[1,]==0)/length(scrb.scran)*100
(drop_out=apply(scrb.scran,1,function(x) sum((x==0)/length(scrb.scran)*100)))
drop_out=as.data.frame(drop_out)
plotting_df$drop_out=drop_out$drop_out
plot(x=log10(plotting_df$mean),y=plotting_df$drop_out,pch=1,cex=1,col="grey50",lwd=0.5,
     ylab = "gene drop-out rate in % across samples",
     xlab = "gene mean expression value in log10 across samples",
     main= "nii")
########################################################################

########################################################################
sce$UMIs=colSums(assay(sce))
sce$Genes=colSums(assay(sce)>1)
sce$wellKey=colnames(sce)
sca = SceToSingleCellAssay(sce)
sca$Treatment=factor(sca$Treatment)
sca$Treatment=relevel(sca$Treatment,"Unstim")

sca$clusters=nii[sca$clusters]
sca$clusters=factor(sca$clusters)
sca$clusters=relevel(sca$clusters,"Unstim")
################################################
set.seed(123)
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
fcHurdleSig_all=fcHurdle[qval<.05 & abs(coef) > 0.5]
dim(fcHurdleSig_all)
head(fcHurdleSig_all)
setorder(fcHurdleSig_all, qval)
head(fcHurdleSig_all)
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

plot_ly(volcano, x=~LFC, y=~-log10(FDR),
        color=~group,colors=cc[c(2,3,1)],
        text=~Gene,type="scatter",size = ~-log10(FDR),
        mode="markers") %>%
  layout(title = 'Volcano',
         yaxis = list(zeroline = FALSE),
         xaxis = list(zeroline = FALSE))  



################### Extra scripts
# BC_names=colnames(reads_counts)
# BC_names=as.data.frame(BC_names)
# colnames(BC_names)="V1"
# Barcodes=Barcodes[c("V3","V2","V1")]
# colnames(Barcodes)=c("V1","V2","V3")
###merging 2 DFs based on the common "V1" and keeping all elements from the first DF
# Merged_names=merge(BC_names,Barcodes,by="V1",all.x=T)
# Merged_names2=Merged_names
# sum(is.na(Merged_names$V3))
# Merged_names2$V2 <- ifelse(is.na(Merged_names2$V2)==1,paste0(rownames(Merged_names2),"X"),Merged_names2$V2)
# Merged_names2$V3 <- ifelse(is.na(Merged_names2$V3)==1,paste0("X",rownames(Merged_names2)),Merged_names2$V3)
### assigning names to count matrix
# reads_counts_named=reads_counts
# colnames(reads_counts_named)=Merged_names2$V3
# colSums(reads_counts_named)

# reads_counts_named=reads_counts_named[, -grep("X", colnames(reads_counts_named))]
# reads_counts_named=reads_counts_named[!(colSums(reads_counts_named)<3000)]
#### using scran package
# matrix=as.matrix(read_counts)
# matrix=matrix[,which(colSums(matrix>0)>1500)]
# matrix=matrix[rowMeans(matrix)>0.1,]
# sce <- SingleCellExperiment(list(counts=matrix))
# scran <- computeSumFactors(sce,sizes=c(20),positive=T)
# norm <- convertTo(scran, type="DESeq2")
# norm_counts=counts(norm,normalized=T)
# norm_counts=norm_counts[,order(colnames(norm_counts))]
# coludata=as.data.frame(colnames(norm_counts))
# coludata$Treatment = "Unstim"
# names(coludata)=c("Barcode","Treatment")
# coludata$Treatment[grep("7|8|9|10|11|12",coludata$Barcode)]="LPS"
# nrow(coludata)
# norm_round_counts=round(norm_counts)
# ncol(norm_round_counts)
# data_dds <- DESeqDataSetFromMatrix(countData = norm_round_counts,
#                                    colData = coludata,
#                                    design = ~ Treatment)
# sizeFactors(data_dds)=c(rep(1,length(coludata$Barcode)))
# data_dds$Treatment <- factor(data_dds$Treatment, levels = c("Unstim","LPS"))
# data_dds=DESeq(data_dds)
# test1=varianceStabilizingTransformation(data_dds)


# display.brewer.all() 
# colors=brewer.pal(9,"Blues")[3:9]
# 
# Heatmap(clus_mat, clustering_distance_columns = "euclidean",col = colors,
#         clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
#         cluster_rows = T, cluster_columns = T, show_row_dend = T,
#         split = batch, gap = unit(7, "mm"), show_row_names = T)


#



#


#

#


#


#


#


#

#


# # Read in first one?
# rds_file="/Home Office/Data/2017-18_Early_seq_results/together with brb/SCRB-seq.dgecounts.rds"
# rds_here=readRDS(rds_file)
# reads=rds_here$exons$downsampled$downsampled_5$umicounts_downsampled
# reads=as.data.frame(as.matrix(reads))
# mean(colSums(reads>0))
# median(colSums(reads>0))
# # Top expr genes
# sort(rowMeans(reads),decreasing = T)[1:50]
# keenid=names(sort(rowMeans(reads),decreasing = T)[1:50])
# genes=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv")
# genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
# genes[match(keenid,genes$ensembl_gene_id),"description"]
# # Top var genes
# keenid2=names(sort(apply(reads, 1,var),decreasing = T)[1:50])
# genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
# # plots
# barplot(colSums(reads))
# summary(colSums(reads))
# lugendid=reads
# colnames(lugendid)
# new_bcs=read.csv("/Home Office/PhD/Experiments/SCRP-Seq_Barcodes/My_1st_set/1st_set.csv",header = F)
# sum(colnames(lugendid)%in%new_bcs$V3)
# colnames(lugendid)=new_bcs[match(colnames(lugendid),new_bcs$V3),"V1"]
# lugendid=lugendid[,!is.na(colnames(lugendid))]
# lugendid=lugendid[,mixedorder(colnames(lugendid))]
# dim(lugendid)
# colnames(lugendid)
# barplot(colSums(lugendid))
# summary(colSums(lugendid))
# # Renaming genes
# Gene_data=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv",col.names =c("gene_id","gene_name","description") )
# Gene_data=Gene_data[Gene_data$gene_id%in%rownames(lugendid),]
# paste("You have ", sum(duplicated(Gene_data$gene_name)),"duplicated gene names in your data")
# Gene_data=Gene_data[!duplicated(Gene_data$gene_name),]
# lugendid=lugendid[rownames(lugendid)%in%Gene_data$gene_id,]
# paste("You have ", sum(duplicated(Gene_data$gene_name)),"duplicated gene names in your data")
# # Order
# rownames(lugendid)=Gene_data[match(rownames(lugendid),Gene_data$gene_id),"gene_name"]
# messed_up_plate=lugendid
# messed_up_plate[1:5,1:5]
# # Seurat
# messed_dots=CreateSeuratObject(counts = messed_up_plate, project = "mess?")
# messed_dots[["percent.mt"]]=PercentageFeatureSet(messed_dots, pattern = "^MT-")
# VlnPlot(messed_dots, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Mito_cut=30
# Feature_hig=40
# messed_dots=subset(messed_dots, subset = nFeature_RNA < Feature_hig & percent.mt< Mito_cut)
# FeatureScatter(messed_dots, feature1 = "nFeature_RNA", feature2 = "percent.mt")
# Feature_low=10
# messed_dots=subset(messed_dots, subset = nFeature_RNA > Feature_low)
# dim(messed_dots)
# # Normalizing with SCtransform
# messed_dots=SCTransform(messed_dots)
# messed_dots$stim="Unstim."
# messed_dots$well=rownames(messed_dots@meta.data)
# messed_dots$stim[grepl("7|8|9|10|11|12",messed_dots$well)]="LPS"
# # PCA
# messed_dots=RunPCA(messed_dots, npcs = 30, verbose = FALSE)
# DimPlot(messed_dots,reduction = "pca",group.by = "stim",pt.size = 2)+ggtitle("First Plate")
# # Umap
# messed_dots=RunUMAP(messed_dots, dims = 1:30, n.neighbors = 30,umap.method = "umap-learn", 
#                 metric = "correlation",min.dist = 0.01,verbose = FALSE)
# DimPlot(messed_dots,reduction = "umap",group.by = "stim",pt.size = 2)+ggtitle("First Plate")

###############################
###############################
#### FIRST PLATE THAT DIDNT WORK!
# rds_file="/Home Office/Data/2017-18_Early_seq_results/01 First seq 2016/Plate1_2016.dgecounts.rds"
# rds_here=readRDS(rds_file)
# reads=rds_here$exons$downsampled$downsampled_477719$umicounts_downsampled
# reads=rds_here$intron.exon$umicounts
# head(reads)
# reads=as.data.frame(as.matrix(reads))
# head(reads)
# mean(colSums(reads>0))
# median(colSums(reads>0))
# Top genes
# sort(rowMeans(reads),decreasing = T)[1:50]
# keenid=names(sort(rowMeans(reads),decreasing = T)[1:50])
# genes=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv")
# genes[match(keenid,genes$ensembl_gene_id),"external_gene_name"]
# genes[match(keenid,genes$ensembl_gene_id),"description"]
# barplot(colSums(reads))
# summary(colSums(reads))
# reads2=reads[colSums(reads)>1000]
# dim(reads2)
# barplot(colSums(reads2))
# lugendid=reads
# old_bcs=read.csv("/Home Office/PhD/Experiments/SCRP-Seq_Barcodes/SET1.csv",header = F)
# head(old_bcs)
# sum(colnames(lugendid)%in%old_bcs$V3)
# colnames(lugendid)=old_bcs[match(colnames(lugendid),old_bcs$V3),"V1"]
# dim(lugendid)
# lugendid=lugendid[,!is.na(colnames(lugendid))]
# colSums(lugendid)
# lugendid=lugendid[,mixedorder(colnames(lugendid))]
# leer=lugendid[,grepl("9|10|11|12",colnames(lugendid))]
# dim(leer)
# head(leer)
# summary(colSums(leer))
# Blaer1=lugendid[,grepl("5|6|7|8",colnames(lugendid))]
# dim(Blaer1)
# head(Blaer1)
# summary(colSums(Blaer1))
# THP=lugendid[,!grepl("5|6|7|8|9|10|11|12",colnames(lugendid))]
# dim(THP)
# head(THP)
# summary(colSums(THP))
# write.csv(colSums(lugendid),"/Home Office/nii.csv")
# kosti=as.data.frame(as.matrix(restimT1$umicount$inex$downsampling$downsampled_))
# kosti[1:5,1:5]
# bisiis=read.csv("/Home Office/PhD/Experiments/SCRP-Seq_Barcodes/My_1st_set/1st_set.csv",header = F)
# head(bisiis)
# colnames(kosti)=bisiis[match(substr(colnames(kosti),9,16),bisiis$V3),"V1"]
# colSums(kosti[,colnames(kosti)=="H6"])
# colSums(kosti[,colnames(kosti)=="H8"])
########################################################################################################
########################################################################################################

###############################
###############################

