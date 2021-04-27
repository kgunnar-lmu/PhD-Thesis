options(stringsAsFactors = F)
library(Matrix)
library(ggplot2)
library(dplyr)
library(gtools)
library(tidyr)
library(cowplot)
library(ggsci)
library(data.table)
library(Matrix)
library(tidyverse)
library(plotly)
library(ggrepel)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggsci)
library(data.table)
library(Seurat)
library(biomaRt)
library(ggsci)
library(RColorBrewer)

# samples=read.delim("/Home Office/Data/2019.07.18.T-cells/Sample_table.tsv")
# head(samples)
# # Reading in zUMIs output
# (zDir<-"/Home Office/Data/2019.07.18.T-cells/T-cells_lane1_binning/")
# (Tcell_macs=list.files(paste0(zDir,"expression/")))
# Tcell_macs=readRDS(paste0(zDir,"expression/",Tcell_macs))
# # Ex In-Ex downS
# summary(colSums(as.matrix(Tcell_macs$umicount$exon$all))/colSums(as.matrix(Tcell_macs$umicount$inex$all)))
# summary(colSums(as.matrix(Tcell_macs$umicount$exon$downsampling$downsampled_))/colSums(as.matrix(Tcell_macs$umicount$inex$downsampling$downsampled_)))
# #
# Tcell_macs_ds_ex=as.data.frame(as.matrix(Tcell_macs$umicount$inex$downsampling$downsampled_))
# # STATS
# stats_macs=samples
# stats_macs=stats_macs[stats_macs$Experiment==1,]
# genecounts=read.delim("/Home Office/Data/2019.07.18.T-cells/T-cells_lane1_binning/stats/All_together_binning.genecounts.txt")
# genecounts=genecounts[genecounts$type=="Intron+Exon",]
# UMIcounts=read.delim("/Home Office/Data/2019.07.18.T-cells/T-cells_lane1_binning/stats/All_together_binning.UMIcounts.txt")
# UMIcounts=UMIcounts[UMIcounts$type=="Intron+Exon",]
# stats_macs$UMIs=UMIcounts$Count[match(stats_macs$XC,UMIcounts$SampleID)]
# stats_macs$Genes=genecounts$Count[match(stats_macs$XC,genecounts$SampleID)]
# stats_macs=stats_macs[!stats_macs$Well%in%c("H6","H7"),]
# head(stats_macs)
# Tcell_macs_ds_ex=Tcell_macs_ds_ex[colnames(Tcell_macs_ds_ex)%in%stats_macs$XC]
# Tcell_macs_ds_ex=Tcell_macs_ds_ex[order(colnames(Tcell_macs_ds_ex))]
# stats_macs=stats_macs[stats_macs$XC%in%colnames(Tcell_macs_ds_ex),]
# stats_macs=stats_macs[order(stats_macs$XC),]
######################################################################################################################################################
# # Reading in zUMIs output
# (zDir<-"/Home Office/Data/2019.07.18.T-cells/T-cells_presorted_binning/")
# (Tcell_facs=list.files(paste0(zDir,"expression/")))
# Tcell_facs=readRDS(paste0(zDir,"expression/",Tcell_facs))
# # Ex In-Ex downS
# summary(colSums(as.matrix(Tcell_facs$umicount$exon$all))/colSums(as.matrix(Tcell_facs$umicount$inex$all)))
# summary(colSums(as.matrix(Tcell_facs$umicount$exon$downsampling$downsampled_))/colSums(as.matrix(Tcell_facs$umicount$inex$downsampling$downsampled_)))
# #
# Tcell_facs_ds_ex=as.data.frame(as.matrix(Tcell_facs$umicount$inex$downsampling$downsampled_))
# # STATS
# stats_facs=samples
# stats_facs=stats_facs[stats_facs$Experiment==2,]
# genecounts=read.delim("/Home Office/Data/2019.07.18.T-cells/T-cells_presorted_binning/stats/All_pre_binning.genecounts.txt")
# genecounts=genecounts[genecounts$type=="Intron+Exon",]
# UMIcounts=read.delim("/Home Office/Data/2019.07.18.T-cells/T-cells_presorted_binning/stats/All_pre_binning.UMIcounts.txt")
# UMIcounts=UMIcounts[UMIcounts$type=="Intron+Exon",]
# stats_facs$UMIs=UMIcounts$Count[match(stats_facs$XC,UMIcounts$SampleID)]
# stats_facs$Genes=genecounts$Count[match(stats_facs$XC,genecounts$SampleID)]
# stats_facs=stats_facs[!stats_facs$Well%in%c("H6","H7"),]
# head(stats_facs)
# Tcell_facs_ds_ex=Tcell_facs_ds_ex[colnames(Tcell_facs_ds_ex)%in%stats_facs$XC]
# Tcell_facs_ds_ex=Tcell_facs_ds_ex[order(colnames(Tcell_facs_ds_ex))]
# stats_facs=stats_facs[stats_facs$XC%in%colnames(Tcell_facs_ds_ex),]
# stats_facs=stats_facs[order(stats_facs$XC),]
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
# stats_together=rbind(stats_facs,stats_macs)
# stats_together=na.omit(stats_together)
# ###
# Tcell_facs_ds_ex[1:5,1:5]
# Tcell_macs_ds_ex[1:5,1:5]
# Tcell_facs_ds_ex$ens=rownames(Tcell_facs_ds_ex)
# Tcell_macs_ds_ex$ens=rownames(Tcell_macs_ds_ex)
# Tcells_together=full_join(Tcell_macs_ds_ex,Tcell_facs_ds_ex,by="ens")
# rownames(Tcells_together)=Tcells_together$ens
# Tcells_together=Tcells_together[colnames(Tcells_together)!="ens"]
# dim(Tcells_together)
# Tcells_together[is.na(Tcells_together)]=0
# Tcells_together[1:5,1:5]
# 
# saveRDS(Tcells_together,"/Home Office/Data/2019.07.18.T-cells/Analysis/New analysis/together_reads.rds")
# saveRDS(stats_together, "/Home Office/Data/2019.07.18.T-cells/Analysis/New analysis/together_annotation.rds")


#######################################################################################################################
# SEURAT
#######################################################################################################################


pal <- c("blue", "black", "purple","red")
# Load Data
DGE<- readRDS("/Home Office/Data/2019.07.18.T-cells/Analysis/New analysis/together_reads.rds")
dim(DGE)
cell_annotation<-readRDS("/Home Office/Data/2019.07.18.T-cells/Analysis/New analysis/together_annotation.rds")
dim(cell_annotation)
# GIVING SYMBOL NAMES TO GENES
# Adding my ERCCs to the complete gene list
ERCCs=read.table("/Home Office/Data/ERCCs/ERCC_sequences.txt",header = T)
ERCCs=data.frame(gene_id=ERCCs$ERCC_ID,gene_name=ERCCs$ERCC_ID,description="ERCC, internal control")
head(ERCCs)
Complete_gene_names=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv",col.names =c("gene_id","gene_name","description") )
head(Complete_gene_names)
Complete_gene_names=rbind(Complete_gene_names,ERCCs)
dim(Complete_gene_names)
Complete_gene_names=Complete_gene_names[Complete_gene_names$gene_id%in%rownames(DGE),]
dim(Complete_gene_names)
Complete_gene_names=Complete_gene_names[!duplicated(Complete_gene_names$gene_name),]
dim(Complete_gene_names)
nrow(DGE)
DGE=DGE[rownames(DGE)%in%Complete_gene_names$gene_id,]
nrow(DGE)
rownames(DGE)=Complete_gene_names$gene_name[match(rownames(DGE),Complete_gene_names$gene_id)]

# Removing mito? and ribo?
# (mito=Complete_gene_names[grep("Mt", Complete_gene_names$gene_biotype),]$gene_name)
# (ribo=Complete_gene_names[grep("rRNA*", Complete_gene_names$gene_biotype),]$gene_name)
# dim(DGE_filtered)
# DGE_filtered=DGE_filtered[!(rownames(DGE_filtered) %in% mito),]
# dim(DGE_filtered)
# DGE_filtered=DGE_filtered[!(rownames(DGE_filtered) %in% ribo),]
# dim(DGE_filtered)

# Create Seurat object
Dotto <- CreateSeuratObject(counts = DGE, min.cells = 10, min.features = 200, project = "T-cells")
Dotto@meta.data<- cbind(Dotto@meta.data,cell_annotation[match(colnames(DGE), cell_annotation$XC),])
# QC and selecting cells for further analysis
Dotto[["percent.mt"]] <- PercentageFeatureSet(Dotto, pattern = "^MT-")
Dotto[["percent.ercc"]] <- PercentageFeatureSet(Dotto, pattern = "^ERCC-")
VlnPlot(Dotto, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ercc"), ncol = 4)
# Filtering
Dotto <- subset(Dotto, subset =UMIs< 80000 & nFeature_RNA< 10000 & percent.mt< 11 & percent.ercc< 20)
VlnPlot(Dotto, features = c("nFeature_RNA","UMIs", "nCount_RNA", "percent.mt","percent.ercc"), ncol = 5)

cc1 <- c(brewer.pal(name = "YlOrRd",n=9))[6]#red
cc2 <- c(brewer.pal(name = "YlOrBr",n=9))[5]#or
cc3 <- c(brewer.pal(name = "YlOrRd",n=9))[3]#yellow
cc4 <- c(brewer.pal(name = "YlGn",n=9))[6]#green
cc5 <- c(brewer.pal(name = "YlGnBu",n=9))[5]#lightblue
cc6 <- c(brewer.pal(name = "Blues",n=9))[6]#blue
cc7 <- c(brewer.pal(name = "Purples",n=9))[6]#viol
cc8 <- c(brewer.pal(name = "RdPu",n=9))[7]#purp

cc=c(cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8)


#################################################################################################################################
"First with no regress and standard workflow"
T_seu2 <- NormalizeData(Dotto, normalization.method = "LogNormalize", scale.factor = 10000)
T_seu2 <- FindVariableFeatures(T_seu2, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(T_seu2)
T_seu2=ScaleData(T_seu2,features = all.genes)
T_seu2=RunPCA(T_seu2, npcs = 100, verbose = FALSE)
T_seu2=RunUMAP(T_seu2, dims = 1:100, n.neighbors = 30,umap.method = "umap-learn", metric = "correlation",min.dist = 0.3,verbose = FALSE)
U1=DimPlot(T_seu2, label = TRUE,group.by = "Time",pt.size = 0.5)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
U2=DimPlot(T_seu2, label = TRUE,group.by = "Experiment",pt.size = 0.5)+scale_color_manual(values = cc[c(1,5)])
plot_grid(U1,U2)
####################################################################
# T_seux=T_seu2[,T_seu2$Time%in%c(72,144)]
# T_seux <- NormalizeData(T_seux, normalization.method = "LogNormalize", scale.factor = 10000)
# T_seux <- FindVariableFeatures(T_seux, selection.method = "vst", nfeatures = 5000)
# all.genes <- rownames(T_seux)
# T_seux=ScaleData(T_seux,features = all.genes)
# T_seux=ScaleData(T_seux,vars.to.regress = "Set")
# T_seux=RunPCA(T_seux, npcs = 100, verbose = FALSE)
# T_seux=RunUMAP(T_seux, dims = 1:100, n.neighbors = 30,umap.method = "umap-learn", metric = "correlation",min.dist = 0.3,verbose = FALSE)
# Px1=DimPlot(T_seux, label = TRUE,group.by = "Time",pt.size = 2,reduction = "pca",label.size = 0)+scale_color_manual(values = cc[c(1,5,7)])
# Px2=DimPlot(T_seux, label = TRUE,group.by = "Experiment",pt.size = 2,reduction = "pca",label.size = 0)+scale_color_manual(values = cc[c(1,5)])
# plot_grid(Px1,Px2)
# FeaturePlot(T_seux,features = c("IL13","IFNG","CSF2","IL2RA"),order = T,reduction = "pca")
# Ux1=DimPlot(T_seux, label = TRUE,group.by = "Time",pt.size = 2)+scale_color_manual(values = cc[c(1,5,7)])
# Ux2=DimPlot(T_seux, label = TRUE,group.by = "Experiment",pt.size = 2)+scale_color_manual(values = cc[c(1,5)])
# Ux3=DimPlot(T_seux, label = TRUE,group.by = "Set",pt.size = 2)+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
# plot_grid(Ux1,Ux2,Ux3)
# FeaturePlot(T_seux,features = c("IL13","IFNG","CSF2","IL2RA"),order = T,reduction = "umap")

####################################################################
# trying to mitigate the batch effect with regressing out unwanted variation
# conf_var=c("UMIs","percent.mt","percent.ercc","Batch","Experiment")
# conf_var=c("Set","Experiment","UMIs")
# T_seu2=ScaleData(T_seu2,vars.to.regress = conf_var)
# T_seu2=ScaleData(T_seu2,vars.to.regress = conf_var,features = top100)
# T_seu2=RunPCA(T_seu2, npcs = 100, verbose = FALSE)
# T_seu2=RunUMAP(T_seu2, dims = 1:100, n.neighbors = 20,umap.method = "umap-learn", metric = "correlation",min.dist = 0.6,verbose = FALSE)
# U3=DimPlot(T_seu2, label = TRUE,group.by = "Time",pt.size = 0.5)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
# U4=DimPlot(T_seu2, label = TRUE,group.by = "Experiment",pt.size = 0.5)+scale_color_manual(values = cc[c(1,5)])
# plot_grid(U3,U4)
# 
# FeaturePlot(T_seu2,features = c("MAOA","RORA","IL17F","CTLA4"))
# FeaturePlot(T_seu2,features = c("IL13","ONECUT3","IFNG","CSF2","GZMB","S100A4"))
# FeaturePlot(T_seu2,features =c("HLA-DPB1","LGALS1","IL2RA","GZMB"))
############################################
##
# conf_var=c("Batch")
# T_seu2=ScaleData(T_seu2,vars.to.regress = conf_var)
# T_seu2=RunPCA(T_seu2, npcs = 100, verbose = FALSE)
# T_seu2=RunUMAP(T_seu2, dims = 1:100, n.neighbors = 20,umap.method = "umap-learn", metric = "correlation",min.dist = 0.3,verbose = FALSE)
# U3=DimPlot(T_seu2, label = TRUE,group.by = "Time",pt.size = 0.5)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
# U4=DimPlot(T_seu2, label = TRUE,group.by = "Experiment",pt.size = 0.5)+scale_color_manual(values = cc[c(1,5)])
#################################################################################################################################
# conf_var=c("Set")
# T_seu2=ScaleData(T_seu2,vars.to.regress = conf_var)
# T_seu2=RunPCA(T_seu2, npcs = 100, verbose = FALSE)
# T_seu2=RunUMAP(T_seu2, dims = 1:100, n.neighbors = 20,umap.method = "umap-learn", metric = "correlation",min.dist = 0.3,verbose = FALSE)
# U5=DimPlot(T_seu2, label = TRUE,group.by = "Time",pt.size = 0.5)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
# U6=DimPlot(T_seu2, label = TRUE,group.by = "Experiment",pt.size = 0.5)+scale_color_manual(values = cc[c(1,5)])
# plot_grid(U1,U2,U3,U4,U5,U6,ncol=2)
#################################################################################################################################
"Now with sctransform"
#################################################################################################################################
head(Dotto@meta.data)
# conf_var=c("UMIs","percent.mt","percent.ercc","Batch","Experiment")
# conf_var=c("Set","UMIs","percent.mt","percent.ercc","Batch","Experiment")
conf_var=c("Set","UMIs")
# conf_var=c("UMIs","Set","Experiment")
# here=FindVariableFeatures(Dotto)
# top500=VariableFeatures(here)[1:1000]
# top500=top500[!grepl("ERCC|AC[0-9]|LINC|AL[0-9]|AF[0-9]|FP[0-9]|AP[0-9]",top500)]
# top500
# top100=top500[1:100]
# top200=top500[1:200]
# top300=top500[1:300]
# top400=top500[1:400]
# 
# T_seu=SCTransform(Dotto,return.only.var.genes = F,vars.to.regress = conf_var)
# T_seu=SCTransform(Dotto,vars.to.regress = conf_var,variable.features.n = 50)
# T_seu=SCTransform(Dotto,variable.features.n = 50)
T_seu=RunPCA(T_seu, npcs = 100, verbose = FALSE)
T_seu=RunUMAP(T_seu, dims = 1:100, n.neighbors = 10,umap.method = "umap-learn", metric = "correlation",min.dist = 0.6,verbose = FALSE)
# DimPlot(T_seu, label = TRUE,group.by = "Time",pt.size = 0.5,reduction = "pca")+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
U3=DimPlot(T_seu, label = TRUE,group.by = "Time",pt.size = 0.5)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
U4=DimPlot(T_seu, label = TRUE,group.by = "Experiment",pt.size = 0.5)+scale_color_manual(values = cc[c(1,5)])
plot_grid(U3,U4)

conf_var=c("Batch","UMIs")
T_seu=SCTransform(Dotto,return.only.var.genes = F,vars.to.regress = conf_var)
T_seu=RunPCA(T_seu, npcs = 100, verbose = FALSE)
T_seu=RunUMAP(T_seu, dims = 1:100, n.neighbors = 10,umap.method = "umap-learn", metric = "correlation",min.dist = 0.6,verbose = FALSE)
U5=DimPlot(T_seu, label = TRUE,group.by = "Time",pt.size = 0.5)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
U6=DimPlot(T_seu, label = TRUE,group.by = "Experiment",pt.size = 0.5)+scale_color_manual(values = cc[c(1,5)])
plot_grid(U5,U6)

conf_var=c("Set","UMIs")
T_seu=SCTransform(Dotto,return.only.var.genes = F,vars.to.regress = conf_var)
T_seu=RunPCA(T_seu, npcs = 100, verbose = FALSE)
T_seu=RunUMAP(T_seu, dims = 1:100, n.neighbors = 10,umap.method = "umap-learn", metric = "correlation",min.dist = 0.6,verbose = FALSE)
U7=DimPlot(T_seu, label = TRUE,group.by = "Time",pt.size = 0.5,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
U8=DimPlot(T_seu, label = TRUE,group.by = "Experiment",pt.size = 0.5,label.size = 0)+scale_color_manual(values = cc[c(1,5)])
plot_grid(U7,U8)

plot_grid(U1,U3,U5,U7,U2,U4,U6,U8,ncol = 4)

FeaturePlot(T_seu,features = c("IFNG","IL4","TBX21","GATA3","TGFB1","IL13"),ncol = 2,order = T,pt.size = 1.5)
VlnPlot(T_seu,features = c("IFNG","IL4","TBX21","GATA3","TGFB1","IL13"),group.by = "Time",ncol = 2)
VlnPlot(T_seu,features = c("IFNG","IL4","TBX21","GATA3","TGFB1","IL13"),group.by = "Time",ncol = 2)

#################
# FeaturePlot(T_seu,features = c("BCAR3","PLEKHG1","AMER1","VPS9D1","CHAF1A","ZNF80B","FSTL3"),order = T,pt.size = 1)
# VlnPlot(T_seu,group.by = "Time",features = c("BCAR3","PLEKHG1","AMER1","VPS9D1","CHAF1A","ZNF80B","FSTL3"),ncol = 4)
# #################
# T_seu1=T_seu[,T_seu$Experiment==1]
# T_seu2=T_seu[,T_seu$Experiment==2]
# FeaturePlot(T_seu1,features = c("FSTL3"),order = T,pt.size = 1)
# FeaturePlot(T_seu2,features = c("FSTL3"),order = T,pt.size = 1)
# VlnPlot(T_seu1,group.by = "Time",features = c("FSTL3"),ncol = 4)



FeaturePlot(T_seu,features = rownames(T_seu)[grep("IL17",rownames(T_seu))],ncol = 2,order = T)

FeaturePlot(T_seu,features = c("IFNG","IL13","GATA3","IL4","TBX21","GZMB","LTA"),pt.size = 1.5,order = T)

 Ovs144_down=Ovs144[Ovs144$avg_logFC<(-1),]
Ovs144_down=Ovs144_down[Ovs144_down$p_val_adj<0.01,]

write.csv(rownames(Ovs144_down),"/Home Office/Data/xdsa.csv")

xxx=read.csv("/Home Office/Data/All_Gene_Names/list_of_combined_TFs.csv")
head(xxx)
xxx[xxx$x%in%rownames(Ovs144_down),]
# FeaturePlot(T_seu,features = c("IL13","ONECUT3","IFNG","CSF2","GZMB","S100A4"))
# FeaturePlot(T_seu,features = c("IL13","ONECUT3","IFNG","CSF2","GZMB","S100A4"),reduction = "pca")
# FeaturePlot(T_seu,features =c("HLA-DPB1","LGALS1","IL2RA","GZMB"))
# FeaturePlot(T_seu,features = c("IL13","IFNG"),pt.size = 1.2)
#######################################################################################
#######################################################################################
#######################################################################################
# Dotto$Experiment=as.character(Dotto$Experiment)
# T_seu1=SCTransform(Dotto)
# T_seu1=RunPCA(T_seu1, npcs = 100, verbose = FALSE)
# T_seu1=RunUMAP(T_seu1, dims = 1:100, n.neighbors = 10,umap.method = "umap-learn", metric = "correlation",min.dist = 0.6,verbose = FALSE)
# U1=DimPlot(T_seu1, label = TRUE,group.by = "Time",pt.size = 0.5,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
# U1.2=DimPlot(T_seu1, label = TRUE,group.by = "Set",pt.size = 0.5)
# U1.2
# U1
T_seu2=SCTransform(Dotto,vars.to.regress = c("Set","UMIs"))
T_seu2=RunPCA(T_seu2, npcs = 100, verbose = FALSE)
T_seu2=RunUMAP(T_seu2, dims = 1:100, n.neighbors = 10,umap.method = "umap-learn", metric = "correlation",min.dist = 0.6,verbose = FALSE)
# U2=DimPlot(T_seu2, label = TRUE,group.by = "Time",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
# U3=DimPlot(T_seu2, label = TRUE,group.by = "Experiment",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,5)])
# plot_grid(U2,U3)

# T_seu3=SCTransform(Dotto,vars.to.regress = "Batch")
# T_seu3=RunPCA(T_seu3, npcs = 100, verbose = FALSE)
# T_seu3=RunUMAP(T_seu3, dims = 1:100, n.neighbors = 10,umap.method = "umap-learn", metric = "correlation",min.dist = 0.6,verbose = FALSE)
# U3=DimPlot(T_seu3, label = TRUE,group.by = "Time",pt.size = 0.5,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
# U3
# T_seu4=SCTransform(Dotto,vars.to.regress = "Experiment")
# T_seu4=RunPCA(T_seu4, npcs = 100, verbose = FALSE)
# T_seu4=RunUMAP(T_seu4, dims = 1:100, n.neighbors = 10,umap.method = "umap-learn", metric = "correlation",min.dist = 0.6,verbose = FALSE)
# U4=DimPlot(T_seu4, label = TRUE,group.by = "Time",pt.size = 0.5,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,7)])
# U4
# plot_grid(U1,U4,U3,U2,)
# U2.2=DimPlot(T_seu2, label = TRUE,group.by = "Experiment",pt.size = 0.5,label.size = 0)+scale_color_manual(values = cc[c(1,5)])
# U3.2=DimPlot(T_seu3, label = TRUE,group.by = "Experiment",pt.size = 0.5,label.size = 0)+scale_color_manual(values = cc[c(1,5)])
# U4.2=DimPlot(T_seu4, label = TRUE,group.by = "Experiment",pt.size = 0.5,label.size = 0)+scale_color_manual(values = cc[c(1,5)])
# plot_grid(U1,U4,U3,U3.2,U2,U2.2,ncol = 2)
# plot_grid(U4,U4.2,U3,U3.2,U2,U2.2,ncol = 2)
# plot_grid(U4,U4.2,U3,U3.2,U2,U2.2,ncol = 2)

T_seu <- FindNeighbors(T_seu, dims = 1:10)
# T_seu <- FindClusters(T_seu, resolution = 0.2)
T_seu <- FindClusters(T_seu, resolution = 0.3,method = "igraph")
u1=DimPlot(T_seu, label = TRUE,pt.size = 1,label.size = 0)#+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
# u1=DimPlot(T_seu2, label = TRUE,pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
u2=DimPlot(T_seu, label = TRUE,group.by = "Time",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
u3=DimPlot(T_seu, label = TRUE,group.by = "Experiment",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,5)])
plot_grid(u1,u2,u3)

markers_late_early=FindMarkers(T_seu, ident.1 = c("1"), ident.2 = c("0"))
markers_late_early=markers_late_early[markers_late_early$avg_logFC>0.5,]
markers_late_early=markers_late_early[markers_late_early$p_val_adj<0.01,]
head(markers_late_early,20)

markers_late_standyouty=FindMarkers(T_seu, ident.1 = c("1"), ident.2 = c("2"))
markers_late_standyouty=markers_late_standyouty[markers_late_standyouty$avg_logFC>0.5,]
markers_late_standyouty=markers_late_standyouty[markers_late_standyouty$p_val_adj<0.01,]
head(markers_late_standyouty,20)

markers_early_vs_standyouty=FindMarkers(T_seu, ident.1 = c("0"), ident.2 = c("2"))
markers_early_vs_standyouty=markers_early_vs_standyouty[markers_early_vs_standyouty$avg_logFC>0.5,]
markers_early_vs_standyouty=markers_early_vs_standyouty[markers_early_vs_standyouty$p_val_adj<0.01,]
head(markers_early_vs_standyouty,20)



ggplotly(u1)

Idents(T_seu2)=T_seu2$Time
markers_1vs_3=FindMarkers(T_seu2, ident.1 = c("4","8"), ident.2 = c("0","24","72","144"))
markers_1vs_3a=markers_1vs_3[markers_1vs_3$avg_logFC>1,]
markers_1vs_3a=markers_1vs_3a[order(markers_1vs_3a$p_val_adj),]
head(markers_1vs_3a,25)


# pc1=DimPlot(T_seu2, label = TRUE,pt.size = 1,label.size = 0,reduction = "pca")+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
# pc2=DimPlot(T_seu2, label = TRUE,group.by = "Time",pt.size = 1,label.size = 0,reduction = "pca")+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
# pc3=DimPlot(T_seu2, label = TRUE,group.by = "Experiment",pt.size = 1,label.size = 0,reduction = "pca")+scale_color_manual(values = cc[c(1,5)])
# plot_grid(pc1,pc2,pc3)


# set.seed(123)
# muumap=as.data.frame(T_seu2@reductions$umap@cell.embeddings)
# muumap=cbind(muumap,as.data.frame(T_seu2@meta.data))
# head(muumap)
# muumap$Time=factor(muumap$Time,levels = mixedsort(unique(muumap$Time)))
# ux=ggplot(muumap,aes(x=UMAP_1,y=UMAP_2,col=Time))+geom_point()+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
# plot_grid(u2,ux)
# 
# library(mclust)
# library(factoextra)
# library(pheatmap)
# fviz_nbclust(muumap[,1:2], kmeans, method = c("silhouette"))
# fviz_nbclust(muumap[,1:2], kmeans, method = c("wss"))
# fviz_nbclust(muumap[,1:2], kmeans, method = c("wss"))
# 
# cl2 <- kmeans(muumap[,1:2], centers = 4)$cluster
# muumap$Kmeans=as.character(cl2)
# uy=ggplot(muumap,aes(x=UMAP_1,y=UMAP_2,col=Kmeans))+geom_point()+scale_color_manual(values = cc[c(1,4,5,6,7,8)])
# uy
# 
# keenid=VariableFeatures(T_seu)#[1:500]
top10 <- head(VariableFeatures(T_seu), 25)
plot1 <- VariableFeaturePlot(T_seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# keenid=VariableFeatures(T_seu2)[1:500]
# dadaaa=as.data.frame(T_seu2@assays$SCT@counts)
# dadaaa=dadaaa[keenid,]
# dadaaa[1:5,1:5]
# annot=data.frame(Time=as.character(T_seu2@meta.data$Time),Donor=as.character(T_seu2@meta.data$Experiment))
# rownames(annot)=colnames(dadaaa)
# pheatmap(dadaaa,cluster_cols = T,cluster_rows = T,annotation_col = annot,show_colnames = F)


###
FeaturePlot(T_seu2,reduction = "pca",features = c("IL13","IFNG"),pt.size = 1.2,order = T)
FeaturePlot(T_seu2,reduction = "pca",features = c("GZMB","GZMA","CSF2","GNLY"),pt.size = 1, order = T)

FeaturePlot(T_seu2,features = c("IFNG","IL4","IL13","GZMB","CSF2","GNLY"),pt.size = 1,order = T)


keenid=VariableFeatures(T_seu)[1:200]
keenid=keenid[!grepl("ERCC",keenid)]
keenid=keenid[!grepl("AC[0-9]",keenid)]
keenid=keenid[!grepl("AL[0-9]",keenid)]
keenid=keenid[!grepl("AF[0-9]",keenid)]
keenid=keenid[1:50]
dim(T_seu)
dadaaa=as.data.frame(T_seu@assays$SCT@data)
dadaaa=dadaaa[keenid,]
dadaaa[1:5,1:5]
annot=data.frame(Time=as.character(T_seu@meta.data$Time),Donor=as.character(T_seu@meta.data$Experiment))
annot=data.frame(Time=as.character(T_seu@meta.data$Time))
rownames(annot)=colnames(dadaaa)
library(pheatmap)
pheatmap(dadaaa,cluster_cols = T,cluster_rows = T,annotation_col = annot,show_colnames = F)

keenid=VariableFeatures(T_seu2)[1:200]
head(keenid,50)

VlnPlot(T_seu,features = c("GZMB","IFNG","CSF2","GNLY","CCL4","IL13","CD69","S100A4"),group.by = "Time",ncol = 4)
VlnPlot(T_seu2,features = c("GZMB","IFNG","CSF2","GNLY","CCL4","IL13","CD69","TRPM3"),group.by = "Time",ncol = 4)

VlnPlot(T_seu,features = c("ELOF1","SYNE1"),group.by = "Time",ncol = 4)

VlnPlot(T_seu2,features = c("GZMB","IFNG","CSF2","GZMH","GNLY","CCL4"),group.by = "Time")

FeaturePlot(T_seu2,features = c("GZMB","GZMA","CSF2","GNLY"),pt.size = 1, order = T)

FeaturePlot(T_seu2,features = c("CSF2"),pt.size = 1.2)
VlnPlot(T_seu2,features = c("CSF2"),pt.size = 1.2,group.by = "Time")

# p1=DimPlot(T_seu2, reduction = "pca", label = TRUE,pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
# p2=DimPlot(T_seu2, reduction = "pca", label = TRUE,group.by = "Time",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
# p3=DimPlot(T_seu2, reduction = "pca", label = TRUE,group.by = "Experiment",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,2,3,4,5,6,7,8)])
# plot_grid(p1,p2,p3)

Idents(T_seu2)=T_seu2$Time

Ovs144=FindMarkers(T_seu2, ident.1 = "144", ident.2 = "0")
Ovs72=FindMarkers(T_seu2, ident.1 = "72", ident.2 = "0")
Ovs24=FindMarkers(T_seu2, ident.1 = "24", ident.2 = "0")
Ovs4=FindMarkers(T_seu2, ident.1 = "4", ident.2 = "0")
Ovs8=FindMarkers(T_seu2, ident.1 = "8", ident.2 = "0")
loppvs72=FindMarkers(T_seu2, ident.1 = "144", ident.2 = "72")
earlyvlate=FindMarkers(T_seu2, ident.1 = c("72","144"), ident.2 = c("4","8"))
earlyvlate=FindMarkers(T_seu2, ident.1 = c("72","144"), ident.2 = c("0","4","8"))
ovsearly=FindMarkers(T_seu2, ident.1 = c("4","8"), ident.2 = c("0"))

activation_markers=FindMarkers(T_seu2, ident.1 = c("4","8"), ident.2 = c("0","24","72","144"))
diff_markers=FindMarkers(T_seu2, ident.1 = c("72","144"), ident.2 = c("0","4","8","24"))

head(earlyvlate,15)
head(ovsearly,15)
################################################################################################################
volcano_data=Ovs144
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-25]=volcano_data$Symbol[volcano_data$p_val_adj<1e-25]
volcano_data$top[abs(volcano_data$avg_logFC)<1.5]=NA
v144=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
    geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(2,3,4)],0.7))+theme(legend.position = "none")+
    ggtitle("0h vs 144h")
v144


# FOR GO
int=volcano_data
int=int[int$avg_logFC>1,]
int=int[int$p_val_adj<0.01,]
int$Symbol
# ovs144_genes_less=int$Symbol
# ovs144_genes_more=int$Symbol
write.csv(int$Symbol,"/Home Office/Data/2019.07.18.T-cells/Analysis/DE144vs0.csv")

#######################################################################################
options(stringsAsFactors = F)
GO1=read.csv("/Home Office/Data/2019.07.18.T-cells/Analysis/GO144vs0.csv")
head(GO1)
go_plot=GO1
go_plot$order=nrow(go_plot):1
go_plot$log10fdr=-log10(go_plot$FDR)
library(Hmisc)
go_plot$Description=capitalize(go_plot$Description)  
go_plot$col=rev(go_plot$log10fdr)
g1=ggplot(go_plot,aes(x=order,y=Enrichment, color=log10fdr,size=Genes,label=Description)) +geom_point(stat = "identity") + coord_flip()+
  theme_classic() + scale_color_gradient(low="blue", high="red") + xlab("") + ylab("Enrichment") + 
  scale_x_continuous(label = go_plot$Description, breaks=go_plot$order)+
  theme(axis.text.y = element_text(color="#000000"))+
  ggtitle("TOP 20 GO Hits")+labs(size="Number of Genes",color="-log10(FDR)")+
  ylim(-0.5,max(go_plot$Enrichment))+geom_hline(yintercept = 0,linetype='dashed')
g1
#######################################################################################
#######################################################################################
options(stringsAsFactors = F)
GO2=read.csv("/Home Office/Data/2019.07.18.T-cells/Analysis/GO4vs0.csv")
head(GO2)
go_plot=GO2
go_plot$order=nrow(go_plot):1
go_plot$log10fdr=-log10(go_plot$FDR)
go_plot$Description=capitalize(go_plot$Description)  
go_plot$col=rev(go_plot$log10fdr)
head(go_plot)
g2=ggplot(go_plot,aes(x=order,y=Enrichment, color=log10fdr,size=Genes,label=Description)) +geom_point(stat = "identity") + coord_flip()+
  theme_classic() + scale_color_gradient(low="blue", high="red") + xlab("") + ylab("Enrichment") + 
  scale_x_continuous(label = go_plot$Description, breaks=go_plot$order)+
  theme(axis.text.y = element_text(color="#000000"))+
  ggtitle("TOP 20 GO Hits")+labs(size="Number of Genes",color="-log10(FDR)")+
  ylim(-0.5,max(go_plot$Enrichment))+geom_hline(yintercept = 0,linetype='dashed')
g2
plot_grid(g2,g1)
plot_grid(early,v144)

GO1[GO1$Description%in%GO2$Description,]
#######################################################################################
#######################################################################################
volcano_data=Ovs72
dim(volcano_data)
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
dim(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1.5 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1.5 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1.5 ),"group"] <- "Significant&FoldChange"

volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]

volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-25]=volcano_data$Symbol[volcano_data$p_val_adj<1e-25]
volcano_data$top[abs(volcano_data$avg_logFC)<1.5]=NA
head(volcano_data)

v72=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("0h vs 72h")
v72
pp=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=Symbol,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("0h vs 72h")
ggplotly(pp)
##############################################################################################################################################################################
volcano_data=Ovs24
dim(volcano_data)
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
dim(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < .7 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > .7 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > .7 ),"group"] <- "Significant&FoldChange"

volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]

volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-15]=volcano_data$Symbol[volcano_data$p_val_adj<1e-15]
volcano_data$top[abs(volcano_data$avg_logFC)<.7]=NA
head(volcano_data)

v24=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("0h vs 24h")
v24
#######################################################################################
#######################################################################################
plot_grid(v144,v72,v24)
#######################################################################################
VlnPlot(T_seu2,features = c("GZMB","CSF2","ACTB","IL2RA"),group.by = "Time")

#######################################################################################
volcano_data=Ovs8
dim(volcano_data)
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
dim(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 0.5 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "Significant&FoldChange"

volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]

volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-25]=volcano_data$Symbol[volcano_data$p_val_adj<1e-25]
volcano_data$top[abs(volcano_data$avg_logFC)<0.5]=NA
head(volcano_data)

ovs8volc=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group,size=-log2(p_val_adj)))+geom_point()+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("0 vs 8h")
ovs8volc

up_Ovs4=Ovs4[Ovs4$avg_logFC>0.4,]
up_Ovs4=up_Ovs4[up_Ovs4$p_val_adj<0.01,]

up_Ovs8=Ovs8[Ovs8$avg_logFC>0.4,]
up_Ovs8=up_Ovs8[up_Ovs8$p_val_adj<0.01,]

length(rownames(up_Ovs4))
length(rownames(up_Ovs8))

rownames(up_Ovs4)[rownames(up_Ovs4)%in%rownames(up_Ovs8)]
##############################################################
volcano_data=Ovs4
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 0.5 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-10]=volcano_data$Symbol[volcano_data$p_val_adj<1e-10]
volcano_data$top[abs(volcano_data$avg_logFC)<0.5]=NA
early=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("4h vs 0h")
early
#
plot_grid(early,v144)
plot_grid(g1,g2)

int=volcano_data
int=int[int$avg_logFC>(0.5),]
int=int[int$p_val_adj<0.01,]
int$Symbol
ovs4_genes=int$Symbol


write.csv(int$Symbol,"/Home Office/Data/2019.07.18.T-cells/Analysis/0vs4.csv")

ovs144_genes[ovs144_genes%in%ovs4_genes]
# similar genes
ovs4_genes[ovs4_genes%in%ovs144_genes_more]
# unique genes
ovs4_genes[!ovs4_genes%in%ovs144_genes_more]
####################################################################


volcano_data=Ovs8
dim(volcano_data)
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
dim(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 0.5 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-15]=volcano_data$Symbol[volcano_data$p_val_adj<1e-15]
volcano_data$top[abs(volcano_data$avg_logFC)<0.5]=NA
head(volcano_data)
mid_early=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group,size=-log2(p_val_adj)))+geom_point()+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("8h vs 0h")
mid_early
plot_grid(early,mid_early)

int=volcano_data
int=int[int$avg_logFC>(0.5),]
int=int[int$p_val_adj<0.01,]
int$Symbol
mid_early_genes=int$Symbol

early_genes
mid_early_genes
late

yyy=early_genes[!early_genes%in%late_genes]
mid_early_genes[mid_early_genes%in%late_genes]
early_genes[early_genes%in%mid_early_genes]
xxx=xxx[xxx%in%mid_early_genes]
yyy=yyy[!yyy%in%mid_early_genes]
VlnPlot(T_seu2,features = xxx,group.by = "Time")
VlnPlot(T_seu2,features = yyy[c(2,4,6,8,10,11,12,13,14,15,16,18)],group.by = "Time")
##################################
volcano_data=earlyvlate
dim(volcano_data)
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
dim(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-90]=volcano_data$Symbol[volcano_data$p_val_adj<1e-90]
volcano_data$top[abs(volcano_data$avg_logFC)<1.2]=NA
head(volcano_data)
erlate=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("late vs early")
# erlate
int=volcano_data
int=int[int$avg_logFC>1,]
int=int[int$p_val_adj<0.01,]
int$Symbol
early_to_late_genes2=int$Symbol
early_to_late_genes2[!early_to_late_genes2%in%early_to_late_genes]
##############################################################
volcano_data=ovsearly
dim(volcano_data)
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
dim(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 0.5 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-15]=volcano_data$Symbol[volcano_data$p_val_adj<1e-15]
volcano_data$top[abs(volcano_data$avg_logFC)<0.5]=NA
head(volcano_data)
ovsearlyplot=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group,size=-log2(p_val_adj)))+geom_point()+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("Early vs 0h")
plot_grid(ovsearlyplot,erlate)
int=volcano_data
int=int[int$avg_logFC>(0.5),]
int=int[int$p_val_adj<0.01,]
int$Symbol
ovsearlyingenes=int$Symbol
early_to_late_genes[early_to_late_genes%in%ovsearlyingenes]
#
write.csv(ovsearlyingenes,"/Home Office/Data/2019.07.18.T-cells/Analysis/0vsEarly.csv")
write.csv(early_to_late_genes2,"/Home Office/Data/2019.07.18.T-cells/Analysis/early vs late2.csv")

#
volcano_data=activation_markers
dim(volcano_data)
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
dim(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 0.5 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 0.5 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-40]=volcano_data$Symbol[volcano_data$p_val_adj<1e-40]
volcano_data$top[abs(volcano_data$avg_logFC)<0.7]=NA
head(volcano_data)
actvolc=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group,size=-log2(p_val_adj)))+geom_point()+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("activation volcano")
actvolc

int=volcano_data
int=int[int$avg_logFC>(0.5),]
int=int[int$p_val_adj<0.01,]
int$Symbol
act_genes=int$Symbol
#
volcano_data=diff_markers
dim(volcano_data)
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
dim(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-50]=volcano_data$Symbol[volcano_data$p_val_adj<1e-50]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
head(volcano_data)
diffvolc=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("activation volcano")
diffvolc

int=volcano_data
int=int[int$avg_logFC>1,]
int=int[int$p_val_adj<0.01,]
int$Symbol
diff_genes=int$Symbol
diff_genes
write.csv(diff_genes,"/Home Office/Data/2019.07.18.T-cells/Analysis/diff_genes.csv")
act_genes
write.csv(act_genes,"/Home Office/Data/2019.07.18.T-cells/Analysis/act_genes.csv")

#
T_seu2=FindVariableFeatures(T_seu2, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10=head(VariableFeatures(T_seu2), 30)
# plot variable features with and without labels
plot1=VariableFeaturePlot(T_seu2)
plot2=LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

FeaturePlot(T_seu2,features = c("GZMB","IFNG","CSF2","GNLY"),pt.size = 1.5)
FeaturePlot(T_seu2,features = c("IFNG","CSF2"),pt.size = 1.5)
VlnPlot(T_seu2,,features = c("GZMB","CSF2"),group.by = "Time")
VlnPlot(T_seu2,,features = c("GZMB","CSF2","GZMH","GZMA","PRF1"),group.by = "Time",ncol = 5)
VlnPlot(T_seu2,,features = c("GNLY","PRF1"),group.by = "Time")
all=Ovs144
all=all[all$avg_logFC>0.7,]
all=all[all$p_val_adj<0.01,]
all$names=rownames(all)
all[grepl("GZM",all$names),]
all[grepl("PR",all$names),]
#

#

#



FeaturePlot(T_seu2,features = c("IL13","IFNG"),reduction = "umap",pt.size = 2)





################
plot_grid(U1,U2,U3,U4,U5,U6,ncol = 2)
plot_grid(U1,U2,U3,U4,U5,U6,U7,U8,ncol = 2)
plot_grid(U1,U2,U3,U4,U7,U8,ncol = 2)

#################################################################################################################################
FeaturePlot(T_seu,features = c("MAOA","RORA","IL17F","CTLA4"))
#
FeaturePlot(T_seu,features = c("IL13","ONECUT3","GNLY","IFNG","GZMB","CSF2"))
FeaturePlot(T_seu2,features = c("IL13","ONECUT3","GNLY","IFNG","GZMB","CSF2"))
FeaturePlot(T_seu,features = c("IL13","ONECUT3","GNLY","GZMB","CCR7","S100A4"))
FeaturePlot(T_seu,features = c("IL13","ONECUT3","GNLY","GZMB","CCR7","S100A4"))
FeaturePlot(T_seu,features = c("IL13","IL","GNLY","GZMB","CCR7","S100A4"))
###############################################################################
VlnPlot(T_seu, features = c("IL13","ONECUT3","GNLY","IFNG","GZMB","CSF2","IL22","CCR7","S100A4"), slot = "counts", log = TRUE,group.by = "Time")
############################################################################
top20 <- head(VariableFeatures(T_seu2), 20)
plot1 <- VariableFeaturePlot(T_seu2)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
head(VariableFeatures(T_seu2), 50)
############################################################################
top10 <- head(VariableFeatures(T_seu), 20)
plot1 <- VariableFeaturePlot(T_seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
############################################################################
T_seu <- FindNeighbors(T_seu, dims = 1:50)
T_seu <- FindClusters(T_seu, resolution = 0.9)
UX=DimPlot(T_seu, label = TRUE)
plot_grid(U5,U6,UX,ncol=3)
##### DE 
clust5 <- FindMarkers(T_seu, ident.1 = 5, min.pct = 0.25)
head(clust5, n = 10)
clust1 <- FindMarkers(T_seu, ident.1 = 1, min.pct = 0.25)
head(clust1, n = 10)
#######################################################################################
Activation_markers=c("IL2RA","GZMB","CSF2","GNLY","S100A4")
deact_feat=c("EOMES","ONECUT3")

##############################################################################
##############################################################################
TH1_Features=c("IFNG","TBX21","CXCR3","CCR4","CCR6")
TH2_Features=c("IL4","GATA3","IL5","IL13","CCR4","CCR6")
TH17_Features=c("RORC","IL17RA","IL22","CCR4","CCR6")
TH22_Features=c("IL22","CCR10")
TH_reg_Features=c("FOXP3","IL2RA","IL7R")
T_cytotox_Features=c("GZMB","PRF1")

cclow <- c(brewer.pal(name = "YlGnBu",n=9))[8]
cchi <- c(brewer.pal(name = "YlOrRd",n=9))[6]
ccheat=c(cclow,cchi)

nii=c("GATA3","IFNG","TBX21","IL4","IL5","IL13","IL17RA","IL22","CCR4","IL22","FOXP3","IL2RA")
# nii=c("GATA3","IFNG","TBX21","IL4","IL5","IL13","RORC","IL17RA","IL22","CCR4","CCR6","IL22","FOXP3","IL2","IL2RA")
FeaturePlot(T_seu,features = nii,pt.size = 1,cols=ccheat,ncol = 2)

naa=c("ONECUT3","IFNG","CSF2","GZMB","CCR7","S100A4")
FeaturePlot(T_seu,features = naa,pt.size = 1,cols=ccheat)

##############################################################################
##############################################################################


##############################################################################
"DE 0 vs 144"
##############################################################################
Samples<- as.factor(T_seu@meta.data$Time)
names(Samples)<- T_seu@meta.data$XC
Samples[Samples=="72"]=144
Samples[Samples=="0"]=0
Samples[Samples=="4"]=0
Samples[Samples=="8"]=0
Samples[Samples=="24"]=0
Samples=droplevels(Samples)
T_seu@active.ident<- Samples

DE_dotto <- FindMarkers(object = T_seu, ident.1 = "144", ident.2 = "0", test.use = "MAST",logfc.threshold = 0)
head(DE_dotto,20)

clust144 <- FindMarkers(T_seu, ident.1 = 144, min.pct = 0.25)
head(clustX,20)

clust0 <- FindMarkers(T_seu, ident.1 = 0, min.pct = 0.25)
head(clust0[order(clust0$avg_logFC,decreasing = T),],20)
#############################
DE_dotto$Gene<- rownames(DE_dotto)
DE_dotto$p_val_adj[DE_dotto$p_val_adj==0]=(1e-310)
ggplot()+
  geom_point(data=DE_dotto, aes(x=avg_logFC, y=-log10(p_val_adj)), alpha=0.5)+
  geom_point(data=DE_dotto[DE_dotto$p_val_adj <0.05 & DE_dotto$avg_logFC >1.5 | DE_dotto$avg_logFC < -1.5,], 
             aes(x=avg_logFC, y=-log10(p_val_adj)), colour="red")+
  geom_text(data=DE_dotto[DE_dotto$p_val_adj <0.05 & DE_dotto$avg_logFC >1.5 |DE_dotto$avg_logFC < -1.5,], 
            aes(x=avg_logFC, y=-log10(p_val_adj), label=Gene), vjust=-0.2,colour="red")+
  ggtitle("DE  \n 0h  vs. 144h")+theme_classic()
head(DE_dotto)
ppp=ggplot(DE_dotto,aes(label=Gene))+
  geom_point(data=DE_dotto, aes(x=avg_logFC, y=-log10(p_val_adj)), alpha=0.5)+
  geom_point(data=DE_dotto[DE_dotto$p_val_adj <0.05 & DE_dotto$avg_logFC >1.5 | DE_dotto$avg_logFC < -1.5,], 
             aes(x=avg_logFC, y=-log10(p_val_adj)), colour="red")+ ggtitle("DE  \n 0h  vs. 144h")+theme_classic()

library(plotly)
ggplotly(ppp)







