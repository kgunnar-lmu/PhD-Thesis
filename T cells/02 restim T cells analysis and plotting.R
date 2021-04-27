options(stringsAsFactors = F)
library(Matrix)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(gtools)
library(tidyr)
library(cowplot)
library(ggsci)
library(data.table)
library(Matrix)
library(tidyverse)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(ggsci)
library(data.table)
library(Seurat)
library(biomaRt)
library(slingshot)

# # Sample table
# samples=read.csv("/Home Office/Data/2020.01.27.re-Stim_T-cells/Info/sample_table.csv")
# # Reading in zUMIs output
# (zDir<-"/Home Office/Data/2020.01.27.re-Stim_T-cells/lane_1/")
# (Tcell_L1=list.files(paste0(zDir,"expression/")))
# Tcell_L1=readRDS(paste0(zDir,"expression/",Tcell_L1))
# # Ex In-Ex downS
# summary(colSums(as.matrix(Tcell_L1$umicount$exon$all))/colSums(as.matrix(Tcell_L1$umicount$inex$all)))
# summary(colSums(as.matrix(Tcell_L1$umicount$exon$downsampling$downsampled_))/colSums(as.matrix(Tcell_L1$umicount$inex$downsampling$downsampled_)))
# # DGE
# Tcell_L1_ds_ex=as.data.frame(as.matrix(Tcell_L1$umicount$inex$downsampling$downsampled_))
# dim(Tcell_L1_ds_ex)
# # STATS
# stats_L1=samples
# head(samples)
# genecounts=read.delim("/Home Office/Data/2020.01.27.re-Stim_T-cells/lane_1/stats/restim_tcells.genecounts.txt")
# genecounts=genecounts[genecounts$type=="Intron+Exon",]
# UMIcounts=read.delim("/Home Office/Data/2020.01.27.re-Stim_T-cells/lane_1/stats/restim_tcells.UMIcounts.txt")
# UMIcounts=UMIcounts[UMIcounts$type=="Intron+Exon",]
# stats_L1$UMIs=UMIcounts$Count[match(stats_L1$XC,UMIcounts$SampleID)]
# stats_L1$Genes=genecounts$Count[match(stats_L1$XC,genecounts$SampleID)]
# stats_L1=stats_L1[!stats_L1$Well%in%c("H6","H7"),]
# head(stats_L1)
# Tcell_L1_ds_ex=Tcell_L1_ds_ex[colnames(Tcell_L1_ds_ex)%in%stats_L1$XC]
# dim(Tcell_L1_ds_ex)
# stats_L1=stats_L1[stats_L1$XC%in%colnames(Tcell_L1_ds_ex),]
# ######################################################################################################################################################
# # Reading in zUMIs output
# (zDir<-"/Home Office/Data/2020.01.27.re-Stim_T-cells/lane_2/")
# (Tcell_L2=list.files(paste0(zDir,"expression/")))
# Tcell_L2=readRDS(paste0(zDir,"expression/",Tcell_L2))
# # Ex In-Ex downS
# summary(colSums(as.matrix(Tcell_L2$umicount$exon$all))/colSums(as.matrix(Tcell_L2$umicount$inex$all)))
# summary(colSums(as.matrix(Tcell_L2$umicount$exon$downsampling$downsampled_))/colSums(as.matrix(Tcell_L2$umicount$inex$downsampling$downsampled_)))
# # DGE
# Tcell_L2_ds_ex=as.data.frame(as.matrix(Tcell_L2$umicount$inex$downsampling$downsampled_))
# # STATS
# stats_L2=samples
# # stats_L2=stats_L2[stats_L2$Experiment=="Do_87",]
# genecounts=read.delim("/Home Office/Data/2020.01.27.re-Stim_T-cells/lane_2/stats/restim_tcells_lane_2.genecounts.txt")
# genecounts=genecounts[genecounts$type=="Intron+Exon",]
# UMIcounts=read.delim("/Home Office/Data/2020.01.27.re-Stim_T-cells/lane_2/stats/restim_tcells_lane_2.UMIcounts.txt")
# UMIcounts=UMIcounts[UMIcounts$type=="Intron+Exon",]
# stats_L2$UMIs=UMIcounts$Count[match(stats_L2$XC,UMIcounts$SampleID)]
# stats_L2$Genes=genecounts$Count[match(stats_L2$XC,genecounts$SampleID)]
# stats_L2=stats_L2[!stats_L2$Well%in%c("H6","H7"),]
# head(stats_L2)
# Tcell_L2_ds_ex=Tcell_L2_ds_ex[colnames(Tcell_L2_ds_ex)%in%stats_L2$XC]
# dim(Tcell_L2_ds_ex)
# Tcell_L2_ds_ex=Tcell_L2_ds_ex[colnames(Tcell_L2_ds_ex)%in%colnames(Tcell_L1_ds_ex)]
# stats_L2=stats_L2[stats_L2$XC%in%colnames(Tcell_L2_ds_ex),]
# 
# hist(colSums(Tcell_L2_ds_ex),breaks = 30)
# hist(colSums(Tcell_L1_ds_ex),breaks = 30)
# 
# ######################################################################################################################################################
# stats_L1=stats_L1[order(stats_L1$XC),]
# stats_L2=stats_L2[order(stats_L2$XC),]
# stats_together=stats_L1
# head(stats_together)
# stats_together$UMIs=stats_L1$UMIs+stats_L2$UMIs
# stats_together$Genes=stats_L1$Genes+stats_L2$Genes
# head(stats_together)
# ### Merging stuff
# Tcells_stim_together=bind_rows(Tcell_L1_ds_ex %>% add_rownames(), 
#           Tcell_L2_ds_ex %>% add_rownames()) %>% 
#   group_by(rowname) %>% 
#   summarise_all(sum)
# 
# Tcells_stim_together=as.data.frame(Tcells_stim_together)
# rownames(Tcells_stim_together)=Tcells_stim_together$rowname
# Tcells_stim_together=Tcells_stim_together[-1]
# Tcells_stim_together[1:5,1:5]
# 
# saveRDS(Tcells_stim_together,"/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/New_DGE_2020.12.rds")
# saveRDS(stats_together, "/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/New_annot_2020.12.rds")

#######################################################################################################################
# SEURAT
#######################################################################################################################
cc1 <- c(brewer.pal(name = "YlOrRd",n=9))[6]#red
cc2 <- c(brewer.pal(name = "YlOrBr",n=9))[5]#or
cc3 <- c(brewer.pal(name = "YlOrRd",n=9))[3]#yellow
cc4 <- c(brewer.pal(name = "YlGn",n=9))[6]#green
cc5 <- c(brewer.pal(name = "YlGnBu",n=9))[5]#lightblue
cc6 <- c(brewer.pal(name = "Blues",n=9))[6]#blue
cc7 <- c(brewer.pal(name = "Purples",n=9))[6]#viol
cc8 <- c(brewer.pal(name = "RdPu",n=9))[7]#purp

cc=c(cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8)

pal <- c("blue", "black", "purple","red")

# Load Data
DGE<- readRDS("/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/New_DGE_2020.12.rds")
dim(DGE)
cell_annotation<-readRDS("/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/New_annot_2020.12.rds")
head(cell_annotation)
cell_annotation$Stimulus="Resting"
filter=cell_annotation$Rep==1&as.numeric(substr(cell_annotation$Well,2,4))>6
cell_annotation$Stimulus[filter]="Re-stim"
filter=cell_annotation$Rep==2&as.numeric(substr(cell_annotation$Well,2,4))<7
cell_annotation$Stimulus[filter]="Re-stim"
table(cell_annotation$Stimulus)
dim(cell_annotation)
# VARS exp plot
DGE<- readRDS("/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/New_DGE_2020.12.rds")

# GIVING SYMBOL NAMES TO GENES
# Adding my ERCCs to the complete gene list
ERCCs=read.table("/Home Office/Data/ERCCs/ERCC_sequences.txt",header = T)
head(ERCCs,2)
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
# write.csv(cell_annotation,"/Home Office/moo.csv")
cell_annotation=read.csv("/Home Office/moo.csv")
head(cell_annotation)
# Create Seurat object
Dotto <- CreateSeuratObject(counts = DGE, min.cells = 10, min.features = 200, project = "T-cells")
Dotto@meta.data<- cbind(Dotto@meta.data,cell_annotation[match(colnames(DGE), cell_annotation$XC),])
# QC and selecting cells for further analysis
Dotto[["percent.mt"]] <- PercentageFeatureSet(Dotto, pattern = "^MT-")
Dotto[["percent.ercc"]] <- PercentageFeatureSet(Dotto, pattern = "^ERCC-")
Idents(Dotto)=Dotto$Stimulus
VlnPlot(Dotto, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ercc","UMIs"), ncol = 5)
# Filtering
Dotto <- subset(Dotto, subset =UMIs< 100000 & nFeature_RNA< 6000 & percent.mt< 10 & percent.ercc< 37)
VlnPlot(Dotto, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ercc","UMIs"), ncol = 5)
#################################################################################################################################
"First with no regress and standard workflow"
T_seu2 <- NormalizeData(Dotto, normalization.method = "LogNormalize", scale.factor = 10000)
T_seu2 <- FindVariableFeatures(T_seu2, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(T_seu2)
T_seu2=ScaleData(T_seu2,features = all.genes)
T_seu2=RunPCA(T_seu2, npcs = 100, verbose = FALSE)
T_seu2=RunUMAP(T_seu2, dims = 1:100, n.neighbors = 30,umap.method = "umap-learn", metric = "correlation",min.dist = 0.3,verbose = FALSE)
# head(T_seu2@meta.data)
# 
# P1=DimPlot(T_seu2,reduction = "pca", label = TRUE,group.by = "Timepoint")
# P2=DimPlot(T_seu2,reduction = "pca", label = TRUE,group.by = "Stimulus")
# 
U1=DimPlot(T_seu2, label = TRUE,group.by = "Timepoint")
U2=DimPlot(T_seu2, label = TRUE,group.by = "Stimulus")
plot_grid(U1,U2)
# ####################################################################
# 
# FeaturePlot(T_seu,features = c("BCAR3","PLEKHG1","AMER1","VPS9D1","CHAF1A","ZNF80B","FSTL3"),order = T,pt.size = 1)
# VlnPlot(T_seu,group.by = "Timepoint",features = c("BCAR3","PLEKHG1","AMER1","VPS9D1","CHAF1A","ZNF80B","FSTL3"),ncol = 4)

plot_grid(U5,U6,U7,ncol = 3)

# # trying to mitigate the batch effect with regressing out unwanted variation
# # conf_var=c("UMIs","percent.mt","percent.ercc","Batch","Experiment")
# conf_var=c("Set")
# T_seu2=ScaleData(T_seu2,vars.to.regress = conf_var)
# T_seu2=RunPCA(T_seu2, npcs = 100, verbose = FALSE)
# T_seu2=RunUMAP(T_seu2, dims = 1:100, n.neighbors = 30,umap.method = "umap-learn", metric = "correlation",min.dist = 0.3,verbose = FALSE)
# 
# P1=DimPlot(T_seu2,reduction = "pca", label = TRUE,group.by = "Timepoint")
# P2=DimPlot(T_seu2,reduction = "pca", label = TRUE,group.by = "Stimulus")
# plot_grid(P1,P2)
# 
# U3=DimPlot(T_seu2, label = TRUE,group.by = "Timepoint")
# U4=DimPlot(T_seu2, label = TRUE,group.by = "Stimulus")
# plot_grid(U3,U4)
saveRDS(T_seu,"/Home Office/restim_tcells_jan_2020.rds")
#################################################################################################################################
"Now with sctransform"
#################################################################################################################################
head(Dotto@meta.data)
# conf_var=c("UMIs","percent.mt","percent.ercc","Batch","Experiment")
conf_var=c("Set","UMIs")
conf_var=c("Batch","Experiment","UMIs")
T_seu=SCTransform(Dotto,return.only.var.genes = F,vars.to.regress = conf_var,variable.features.n = 5000)
T_seu=SCTransform(Dotto,return.only.var.genes = F,variable.features.n = 5000)
T_seu=RunPCA(T_seu, npcs = 100, verbose = FALSE)
T_seu=RunUMAP(T_seu, dims = 1:100, n.neighbors = 100,umap.method = "umap-learn", metric = "correlation",min.dist = 0.01,verbose = FALSE)
T_seu=RunUMAP(T_seu, dims = 1:100, n.neighbors = 5,umap.method = "umap-learn", metric = "correlation",min.dist = 0.01,verbose = FALSE)
U5=DimPlot(T_seu, label = TRUE,group.by = "Timepoint",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,2,4,5,6,7,8)])#+theme(legend.position = "none")
U6=DimPlot(T_seu, label = TRUE,group.by = "Stimulus",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,5)])+theme(legend.position = "none")
U7=DimPlot(T_seu, label = TRUE,group.by = "Experiment",pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,5)])+theme(legend.position = "none")
plot_grid(U5,U6,U7,ncol = 3)
################
DimPlot(T_seu, label = TRUE,group.by = "Timepoint",pt.size = 1,label.size = 0,reduction = "pca")+
  scale_color_manual(values = cc[c(1,2,4,5,6,7,8)])
################
T_seu <- FindNeighbors(T_seu, dims = 1:10)
T_seu <- FindClusters(T_seu, resolution = 0.6)
U7=DimPlot(T_seu, label = TRUE,pt.size = 1,label.size = 0)+scale_color_manual(values = cc[c(1,2,4,5,6,7,8)])
U7
theme(legend.position = "none")
plot_grid(U5,U7)
################
FeaturePlot(T_seu,features = c("DUSP2","PIM3","MYC","IL2","CD69","TNF","CYCS"),order = T,pt.size = 2)
FeaturePlot(T_seu,features = c("GZMB","CSF2","IFNG","IL2","CCL4","IL13"),order = T,pt.size = 1)
# plot_grid(U1,U2,U3,U4,U5,U6,ncol = 2)
############################################
# top20 <- head(VariableFeatures(T_seu2), 20)
# plot1 <- VariableFeaturePlot(T_seu2)
# plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
# plot2
# head(VariableFeatures(T_seu2), 50)
############################################################################
top10 <- head(VariableFeatures(T_seu), 20)
# all <- VariableFeatures(T_seu)
plot1 <- VariableFeaturePlot(T_seu)
plot2 <- LabelPoints(plot = plot1, points = top10)
ggplotly(plot2)
#################################################################################################################################
# FeaturePlot(T_seu,features = c("CCL4","S100A4","CCR7","IFNG","GZMB","CSF2"))
T_seu_stim=T_seu[,T_seu$Stimulus=="Re-stim"]
T_seu_rest=T_seu[,T_seu$Stimulus!="Re-stim"]
VlnPlot(T_seu_rest,features = c("IL13","IFNG","CSF2","GZMB"),group.by = "Timepoint",ncol = 4)+scale_fill_manual(values = cc[c(1,2,4,5,6,7,8)])
VlnPlot(T_seu_stim,features = c("IL13","IFNG","CSF2","GZMB"),group.by = "Timepoint",ncol = 4)


FeaturePlot(T_seu_rest,features = c("GZMB","CSF2","GZMA","GZMH","PRF1"))
FeaturePlot(T_seu_rest,features = c("CCR7","CD3","CD4","CD45RA"))

FeaturePlot(T_seu,features = c("IL13","IL4","IFNG","IL17RA"),pt.size = 2)
FeaturePlot(T_seu,features = c("IL13","IFNG","IL4","IL22","GATA3","TBX21"),pt.size = 2,order = T,ncol = 2)
FeaturePlot(T_seu,features = c("IL13","IFNG","ACTB","IL22","GATA3","TBX21","TUBB","GAPDH"),pt.size = 2,order = T,ncol = 2)
FeaturePlot(T_seu,features = c("IL6ST","ASIP"),pt.size = 2,order = T,ncol = 2)
VlnPlot(T_seu,features = c("GAPDH","ACTB","TUBB","GATA3","IL13"),group.by = "Timepoint",ncol = 4)
VlnPlot(T_seu,features = c("ASIP","IL6ST"),group.by = "Timepoint",ncol = 4)


xxx=as.data.frame(T_seu@assays$SCT@data)
xxx_var= apply(xxx, 1, var)
xxx_var=xxx_var[order(xxx_var)]
head(xxx_var,40)
xxx_var[names(xxx_var)%in%"IFNG"]
xxx_var_500=xxx_var[1:1000]

yyy=as.data.frame(T_seu@assays$SCT@data)
yyy_sum= rowSums(yyy)
yyy_sum=yyy_sum[order(yyy_sum,decreasing = T)]
head(yyy_sum,40)
yyy_sum_500=yyy_sum[1:1000]
xxx_var_500[names(xxx_var_500)%in%names(yyy_sum_500)]
yyy_sum[names(yyy_sum)%in%names(xxx_var[1:5])]
yyy_sum[names(yyy_sum)%in%"IFNG"]


FeaturePlot(T_seu,features = c("IL13","IFNG","GZMB","CSF2","IL4","TNF","GZMA","GNLY"),pt.size = 2,order = T,ncol = 4)


VlnPlot(T_seu,features = c("IL13","IFNG","CSF2","TNF"),group.by = "Timepoint",ncol = 4)
VlnPlot(T_seu_rest,features = c("IL13","IFNG","CSF2","TNF"),group.by = "Timepoint",ncol = 4)
VlnPlot(T_seu_stim,features = c("IL13","IFNG","CSF2","TNF"),group.by = "Timepoint",ncol = 4)


VlnPlot(T_seu_rest, features = c("GZMB","CSF2","GZMA","GZMH","PRF1"),group.by = c("Timepoint"))
VlnPlot(T_seu_stim, features = c("GZMB","CSF2","GZMA","GZMH","PRF1"),group.by = c("Timepoint"))
#################
FeaturePlot(T_seu_stim,features = c("GZMB","CSF2","GZMA","GZMH","PRF1"))
VlnPlot(T_seu_stim, features = c("GZMB","CSF2","GZMA","GZMH","PRF1"),group.by = c("Timepoint"))
##################################################################################################
thing=as.data.frame(as.matrix(T_seu@assays$SCT@data))
thing=as.data.frame(as.matrix(T_seu@assays$SCT@counts))
thing=as.data.frame(as.matrix(T_seu@assays$RNA@counts))
t_thing=t(thing)
t_thing=as.data.frame(t_thing)
t_thing[1:5,1:5]

ggscatter(t_thing, x = "IL13", y = "IFNG",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
          xlab = "IL13", ylab = "IFNG")

FeaturePlot(T_seu,features = c("IL13","IFNG"),pt.size = 2,slot = "counts",order = T)
FeaturePlot(T_seu,features = c("IL13","IFNG"),pt.size = 2,slot = "data",order = T)

Idents(T_seu)=T_seu$Stimulus
ResVstiM=FindMarkers(T_seu, ident.1 = "Re-stim", ident.2 = "Resting")
head(ResVstiM,10)
volcano_data=ResVstiM
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-5]=volcano_data$Symbol[volcano_data$p_val_adj<1e-5]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("Resting vs restimulated")
#
int=volcano_data
int=int[int$avg_logFC>1,]
int=int[int$p_val_adj<0.01,]
int$Symbol
stim_Genes=int$Symbol
# write.csv(int$Symbol,"/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/stimvsres.csv")
# RESTING GO
GO=read.csv("/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/GOstimvsres.csv")
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

#

Idents(T_seu_stim)=T_seu_stim$Timepoint
ovs144_stim=FindMarkers(T_seu_stim, ident.1 ="7", ident.2 = "0")
head(ovs144_stim,10)
volcano_data=ovs144_stim
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-25]=volcano_data$Symbol[volcano_data$p_val_adj<1e-25]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("0h vs 144h re-stimulated")
#
int=volcano_data
int=int[int$avg_logFC>1,]
int=int[int$p_val_adj<0.01,]
int$Symbol
restim0vs144_genes=int$Symbol
# write.csv(int$Symbol,"/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/stimvsres.csv")
Idents(T_seu_stim)=T_seu_stim$Timepoint
ovs336_stim=FindMarkers(T_seu_stim, ident.1 ="14", ident.2 = "0")
head(ovs336_stim,10)
volcano_data=ovs336_stim
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-25]=volcano_data$Symbol[volcano_data$p_val_adj<1e-25]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
resti=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("0d vs 14d re-stimulated")
#
#################################################################################################################################
Idents(T_seu_rest)=T_seu_rest$Timepoint
ovs336_rest=FindMarkers(T_seu_rest, ident.1 ="14", ident.2 = "0")
head(ovs336_rest,10)
volcano_data=ovs336_rest
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-25]=volcano_data$Symbol[volcano_data$p_val_adj<1e-25]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
puhka=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("0d vs 14d resting")
##################################################################################################################################
plot_grid(resti,puhka)
ovs336_rest_up=ovs336_rest[ovs336_rest$p_val_adj<0.01,]
ovs336_rest_up=ovs336_rest[ovs336_rest$avg_logFC>0.5,]
ovs336_rest_up=ovs336_rest_up[order(ovs336_rest_up$p_val_adj,decreasing = F),]
head(ovs336_rest_up)
top100=rownames(ovs336_rest_up)[1:100]
ovs336_stim_up=ovs336_stim[ovs336_stim$p_val_adj<0.01,]
ovs336_stim_up=ovs336_stim[ovs336_stim$avg_logF>0.5,]
ovs336_stim_up=ovs336_stim_up[order(ovs336_stim_up$p_val_adj,decreasing = F),]
head(ovs336_stim_up)
top100rest=rownames(ovs336_stim_up)[1:100]
sum(top100rest%in%top100)
top100rest[top100rest%in%top100]

top100rest[top100rest%in%top100]
top100rest[!top100rest%in%top100]
top100[!top100%in%top100rest]
##################################################################################################################################

int=volcano_data
int=int[int$avg_logFC>1,]
int=int[int$p_val_adj<0.01,]
int$Symbol
restim0vs336_genes=int$Symbol
# write.csv(int$Symbol,"/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/stimvsres.csv")
restim0vs144_genes
restim0vs144_genes[!restim0vs144_genes%in%stim_Genes]
write.csv(restim0vs144_genes,"/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/0vs144stim.csv")
###
restim0vs336_genes
restim0vs336_genes[!restim0vs336_genes%in%stim_Genes]
write.csv(restim0vs336_genes,"/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/0vs336stim.csv")
####################################################################################################
T_seu$tot=paste0(T_seu$Stimulus,"_",T_seu$Timepoint)
T_seu$tot
Idents(T_seu)=T_seu$tot
early_late_resting=FindMarkers(T_seu, ident.1 ="Resting_14", ident.2 = "Resting_0")
early_late_restim=FindMarkers(T_seu, ident.1 ="Re-stim_14", ident.2 = "Re-stim_0")
Idents(T_seu)=T_seu$Timepoint
early_late=FindMarkers(T_seu, ident.1 ="14", ident.2 = "0")
head(early_late_resting,13)
head(early_late_restim,13)
head(early_late,13)


volcano_data=early_late_resting
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-25]=volcano_data$Symbol[volcano_data$p_val_adj<1e-25]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
v1=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("Early late resting")
v1
####################################
volcano_data=early_late_restim
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-25]=volcano_data$Symbol[volcano_data$p_val_adj<1e-25]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
v2=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("Early late restimulated")
v2
plot_grid(v1,v2)
plot_grid(v2,vx)
t1=early_late_restim[1:50,]
t1$nimed=rownames(t1)
t2=early_late_resting[1:50,]
t2$nimed=rownames(t2)
sum(rownames(t1)%in%rownames(t2))
rownames(t1)[!rownames(t1)%in%rownames(t2)]

200-100
149-x
3200/50
35800/500
#########################################
volcano_data=early_late
volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
volcano_data["group"] <- "NotSignificant"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data$Symbol=rownames(volcano_data)
volcano_data=volcano_data[order(volcano_data$p_val_adj),]
volcano_data$top=NA
volcano_data$top[volcano_data$p_val_adj<1e-50]=volcano_data$Symbol[volcano_data$p_val_adj<1e-50]
volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
v3=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
  geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(2,3,4)],0.7))+theme(legend.position = "none")+
  ggtitle("Early late")
v3
################################
plot_grid(v1,v2,v3)
VlnPlot(T_seu,features = c("IL13","IFNG"),group.by = "Timepoint")

VlnPlot(T_seu_rest,features = c("CSF2"),group.by = "Timepoint")
VlnPlot(T_seu_stim,features = c("CSF2"),group.by = "Timepoint")


###############################################################################
VlnPlot(T_seu, features = c("S100A4","IFNG","GZMB","CCR7","IL2","CSF2"), slot = "counts", log = TRUE,group.by = "Timepoint")
VlnPlot(T_seu, features = c("S100A4","IFNG","GZMB","CCR7","IL2","CSF2"), slot = "counts", log = TRUE,group.by = "Stimulus")
VlnPlot(T_seu, features = c("CCL4"),group.by = c("Stimulus","Timepoint"))
## Looking at them separately
T_seu_resting=subset(T_seu, subset = Stimulus == "Resting")
T_seu_stim=subset(T_seu, subset = Stimulus == "Re-stim")
VlnPlot(T_seu_resting, features = c("S100A4","IFNG","GZMB","CCR7","IL2","CSF2"), slot = "counts", log = TRUE,group.by = "Timepoint")
VlnPlot(T_seu_stim, features = c("S100A4","IFNG","GZMB","CCR7","IL2","CSF2"), slot = "counts", log = TRUE,group.by = "Timepoint")
###
T_seu_stim=RunPCA(T_seu_stim, npcs = 100, verbose = FALSE)
T_seu_stim=RunUMAP(T_seu_stim, dims = 1:100, n.neighbors = 30,umap.method = "umap-learn", metric = "correlation",min.dist = 0.3,verbose = FALSE)
DimPlot(T_seu_stim, label = TRUE,group.by = "Timepoint")

FeaturePlot(T_seu_stim,features = c("GZMB","IFNG","IL2","CSF2"))

FeaturePlot(T_seu,features = c("IFNG","IL12","IL4","IL13","IL5"),pt.size = 2,order = T)
T_seu <- FindNeighbors(T_seu, dims = 1:50)
T_seu <- FindClusters(T_seu, resolution = 1.5)
DimPlot(T_seu,pt.size = 2,reduction = "umap")

T_seu_stim_14=subset(T_seu_stim, Timepoint > 7)

FeaturePlot(T_seu_stim_14,features = c("GZMB","IFNG","IL2","CSF2"))
FeaturePlot(T_seu_stim_14,features = c("GZMB","IFNG","IL2","CSF2"),reduction = "pca")

T_seu_stim_14 <- FindNeighbors(T_seu_stim_14, dims = 1:50)
T_seu_stim_14 <- FindClusters(T_seu_stim_14, resolution = 0.9)
DimPlot(T_seu_stim_14,reduction = "pca")
DimPlot(T_seu_stim_14,reduction = "umap")

High_low <- FindMarkers(T_seu_stim_14, ident.1 = 1, min.pct = 0.25)
head(High_low, n = 15)
High_low=High_low[High_low$p_val_adj<0.05,]

DE_High_low <- FindMarkers(object = T_seu_stim_14, ident.1 = 1, ident.2 = 0, test.use = "MAST",logfc.threshold = 0)
head(DE_High_low,15)
DE_High_low=DE_High_low[DE_High_low$p_val_adj<0.05,]
common=rownames(DE_High_low)[rownames(DE_High_low)%in%rownames(High_low)]

list_of_TFs=read.csv("/Home Office/Data/All_Gene_Names/list_of_combined_TFs.csv")
head(list_of_TFs)
common[common%in%list_of_TFs$x]


FeaturePlot(T_seu_stim_14,features = c("GZMB","IFNG","IL2","CSF2","ZBED2"))

FeaturePlot(T_seu_rest,features = c("GZMB","IFNG","IL2","CSF2","ZBED2"))

FeaturePlot(T_seu,features = c("GZMB","IFNG","IL2","CSF2","ZBED2"))
# FOUND !!!!!!
"ZBED2"
# goes up in where other stuff does too -.-
############################################################################
T_seu=RunPCA(T_seu, npcs = 100, verbose = FALSE)
T_seu=RunUMAP(T_seu, dims = 1:100, n.neighbors = 30,umap.method = "umap-learn", metric = "correlation",min.dist = 0.3,verbose = FALSE)
TH1_Features=c("IFNG","TBX21","CXCR3","CCR4","CCR6")
TH2_Features=c("IL4","GATA3","IL5","IL13","CCR4","CCR6")
TH17_Features=c("RORC","IL17RA","IL22","CCR4","CCR6")
TH22_Features=c("IL22","CCR10")
TH_reg_Features=c("FOXP3","IL2RA","IL7R")
T_cytotox_Features=c("GZMB","PRF1")
##############################################################################
DimPlot(T_seu,group.by = "Timepoint")
FeaturePlot(T_seu,features = TH1_Features)
FeaturePlot(T_seu,features = TH2_Features)
FeaturePlot(T_seu,features = TH17_Features)
FeaturePlot(T_seu,features = TH22_Features)
FeaturePlot(T_seu,features = TH_reg_Features)
FeaturePlot(T_seu,features = T_cytotox_Features)
#############################################################################

metamata=T_seu@meta.data
head(metamata)
metamata$Timepoint=factor(metamata$Timepoint,levels=mixedsort(unique(metamata$Timepoint)))
umi=ggplot(metamata,aes(x=Timepoint,y=UMIs,fill=Stimulus))+geom_bar(stat = "identity")+facet_grid(~Stimulus)+theme_classic()+scale_fill_manual(values = cc[c(1,5)])
gene=ggplot(metamata,aes(x=Timepoint,y=Genes,fill=Stimulus))+geom_bar(stat = "identity")+facet_grid(~Stimulus)+theme_classic()+scale_fill_manual(values = cc[c(1,5)])
plot_grid(umi,gene)

#############################################################################
T_seu_stim=T_seu[,T_seu$Stimulus=="Re-stim"] 
T_Cell_SCE=SingleCellExperiment(assays=list(counts=T_seu_stim@assays$SCT@counts,logcounts=T_seu_stim@assays$SCT@data),
                                reducedDims=SimpleList(pca=T_seu_stim@reductions$pca@cell.embeddings[,1:10]),
                                colData=T_seu_stim@meta.data)

T_SL <- slingshot(T_Cell_SCE, reducedDim = "pca") # SL needs an SCE object otherways it doesnt know how to make stuff up
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

par(mfrow = c(1, 1))
plot(reducedDims(T_SL)$pca, 
     col = colors[cut(T_SL$slingPseudotime_1,breaks=100)], 
     pch=16, asp = 1, main="Slingshot")
lines(SlingshotDataSet(T_SL), lwd=2)
DimPlot(T_seu_stim,reduction = "pca",group.by = "Timepoint")+scale_color_manual(values = cc[c(1,2,4,5,6,7,8)])

xx=T_seu_stim@reductions$pca@cell.embeddings
head(xx)
xx=as.data.frame(xx)
head(xx)
xx$time=as.character(T_seu_stim@meta.data$Timepoint[match(rownames(xx),T_seu_stim@meta.data$XC)])
# xx$time=T_seu_stim@meta.data$Timepoint[match(rownames(xx),T_seu_stim@meta.data$XC)]
ggplot(xx,aes(x=PC_1,y=PC_2,color=time))+geom_point()+theme_minimal()+scale_color_manual(values = cc[c(1,2,4,5,6,7,8)])
ggplot(xx,aes(x=PC_1,y=PC_3,color=time))+geom_point()+theme_minimal()+scale_color_manual(values = cc[c(1,2,4,5,6,7,8)])

library(plotly)
plot_ly(xx, x= ~PC_1, y=~PC_2, z=~PC_3,color=~time)

#############################################################################
T_seu_stim=T_seu[,T_seu$Stimulus=="Re-stim"] 

library(scran)
T_Cell_SCE=SingleCellExperiment(assays=list(counts=T_seu_stim@assays$SCT@counts,logcounts=T_seu_stim@assays$SCT@data),
                                reducedDims=SimpleList(umap=T_seu_stim@reductions$umap@cell.embeddings),
                                colData=T_seu_stim@meta.data)

umap=as.data.frame(reducedDim(T_Cell_SCE, "umap"))
T_Cell_SCE$U1=umap$UMAP_1
T_Cell_SCE$U2=umap$UMAP_2
T_Cell_SCE=T_Cell_SCE[,T_Cell_SCE$U1<0]
T_Cell_SCE$pseudotime_U <- rank(T_Cell_SCE$U1)  # rank cells by their Umap score


T_SL <- slingshot(T_Cell_SCE, reducedDim = "umap") # SL needs an SCE object otherways it doesnt know how to make stuff up
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

par(mfrow = c(1, 1))
plot(reducedDims(T_SL)$umap, 
     col = colors[cut(T_SL$slingPseudotime_1,breaks=100)], 
     pch=16, asp = 1, main="Slingshot")
lines(SlingshotDataSet(T_SL), lwd=2)
T_Cell_SCE$Timepoint

# T_SL$Timepoint=factor(T_SL$Timepoint,levels = mixedsort(unique(T_SL$Timepoint)))
# plot(reducedDims(T_SL)$umap, 
#      col = cc[T_SL$Timepoint], 
#      pch=16, asp = 1, main="Slingshot")
T_seu_stim=T_seu_stim[,T_seu_stim$XC%in%T_Cell_SCE$XC]
DimPlot(T_seu_stim,reduction = "umap",group.by = "Timepoint")+scale_color_manual(values = cc[c(1,2,4,5,6,7,8)])

# T_SL <- slingshot(T_Cell_SCE, reducedDim = "umap",clusterLabels = 'Timepoint',start.clus = 0)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(T_SL$slingPseudotime_1, breaks=100)]
# plot(reducedDims(T_SL)$umap, 
#      col = plotcol, 
#      pch=16, asp = 1, main="Slingshot with starting point")
# lines(SlingshotDataSet(T_SL), lwd=2, col='black')



Y <- log2(counts(T_SL) + 1)
Top_genes=FindVariableFeatures(T_seu_stim, selection.method = "vst", nfeatures = 3000)
var1K <- VariableFeatures(Top_genes)
Y <- Y[var1K, ]  # only counts for variable genes
t <- T_SL$slingPseudotime_1
library(gam)
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p})

important_genes=as.data.frame(gam.pval)
important_genes$Genes=rownames(important_genes)
important_genes=important_genes[order(important_genes$gam.pval),]
important_genes= important_genes[c(2,1)]
head(important_genes,30)
# write.csv(important_genes,"/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/GAM pseudo.csv")
important_genes=important_genes[!grepl("ERCC",rownames(important_genes)),]

topgenes <- rownames(important_genes)[1:100]
pst.ord <- order(T_SL$slingPseudotime_1, na.last = NA)
heatdata <- log1p(as.matrix(assays(T_SL)$counts[topgenes, pst.ord])+1)
heatdata[1:5,1:5]
heatclus <- (T_SL$Timepoint[pst.ord]+1)

heatmap(heatdata, Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

tfs=read.csv("/Home Office/Data/All_Gene_Names/list_of_combined_TFs.csv")

important_genes_cut=important_genes[important_genes$gam.pval<0.001,]
important_genes_cut=important_genes_cut[rownames(important_genes_cut)%in%tfs$x,]
important_genes_cut
write.csv(important_genes_cut,"/Home Office/Data/2020.01.27.re-Stim_T-cells/Analysis/GAM pseudo TFs olnly.csv")
#####################################################
# PHEATMAP!!!
keenid=VariableFeatures(T_seu)[1:200]
keenid=keenid[!grepl("ERCC",keenid)]
keenid=keenid[!grepl("AC[0-9]",keenid)]
keenid=keenid[!grepl("AL[0-9]",keenid)]
keenid=keenid[!grepl("AF[0-9]",keenid)]
keenid=keenid[1:50]

dadaaa=as.data.frame(T_seu@assays$SCT@data)
dadaaa=dadaaa[keenid,]
dadaaa[1:5,1:5]
annot=data.frame(Time=as.character(T_seu@meta.data$Time),Stimulus=as.character(T_seu@meta.data$Stimulus))
# annot=data.frame(Time=as.character(T_seu@meta.data$Time))
rownames(annot)=colnames(dadaaa)
library(pheatmap)
pheatmap(dadaaa,cluster_cols = T,cluster_rows = T,annotation_col = annot,show_colnames = F)

# 
# Idents(T_seu_stim)=T_seu_stim$Timepoint
# xxx2=FindMarkers(T_seu_stim, ident.1 ="14", ident.2 = "0",test.use =  "DESeq2")
# head(xxx,13)
# head(xxx2,13)
# volcano_data=xxx2
# volcano_data=volcano_data[!is.na(volcano_data$p_val_adj),]
# volcano_data["group"] <- "NotSignificant"
# volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) < 1 ),"group"] <- "Significant"
# volcano_data[which(volcano_data$p_val_adj > 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "FoldChange"
# volcano_data[which(volcano_data$p_val_adj < 0.01 & abs(volcano_data$avg_logFC) > 1 ),"group"] <- "Significant&FoldChange"
# volcano_data$Symbol=rownames(volcano_data)
# volcano_data=volcano_data[order(volcano_data$p_val_adj),]
# volcano_data$top=NA
# volcano_data$top[volcano_data$p_val_adj<1e-50]=volcano_data$Symbol[volcano_data$p_val_adj<1e-50]
# volcano_data$top[abs(volcano_data$avg_logFC)<1]=NA
# v2=ggplot(volcano_data,aes(x=avg_logFC,y=-log10(p_val_adj),label=top,color=group))+geom_point(aes(size=-log2(p_val_adj)))+
#   geom_text_repel(color="black")+theme_classic()+scale_color_manual(values = ggplot2::alpha(pal[c(1,2,3,4)],0.7))+theme(legend.position = "none")+
#   ggtitle("deseq")+geom_vline(xintercept = 0)
# plot_grid(v1,v2)
# v2
# v1
