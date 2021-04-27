## Tbrb test stats and analysis:
options(stringsAsFactors = F)
library(cowplot)
library(data.table)
library(gplots)
library(Hmisc)
library(mclust)
library(pheatmap)
library(sctransform)
library(Seurat)
library(tidyverse)
library(stringr)
library(ggrepel)
library("Rsamtools")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("apeglm")
library("gtools")
library("IHW")
library("ashr")
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library(ggbeeswarm)
library("plotly")
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
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

# Filtering and stats
Test_DGE=readRDS("/Home Office/Data/2020.10.10. current t cell experiment/2021.02.01.TEST-seq/P2/zUMIs_output/expression/Tcell_final_test_p2.dgecounts.rds")
head(Test_DGE)
Test_DGE=as.data.frame(as.matrix(Test_DGE$umicount$exon$all))
Test_annot=read.csv("/Home Office/Data/2020.10.10. current t cell experiment/2021.02.01.TEST-seq/Plate2_sample_table.csv")
head(Test_annot)
barcodes=read.csv("/Home Office/PhD/Experiments/SCRP-Seq_Barcodes/My_1st_set/1st_set.csv",header = F)
colnames(Test_DGE)=barcodes$V1[match(colnames(Test_DGE),barcodes$V3)]

# Consider throwing away ribos and MTs
Complete_gene_names=read.delim("/Home Office/Data/All_Gene_Names/Complete_gene_names_ENS.tsv",col.names =c("gene_id","gene_name","description") )
head(Complete_gene_names)
# Complete_gene_names[Complete_gene_names$gene_name=="MTATP6P1",]
(mito1=Complete_gene_names[grep("mitochondria", Complete_gene_names$description,ignore.case = T),]$gene_name)
(mito2=Complete_gene_names[grep("^MT-", Complete_gene_names$gene_name,ignore.case = T),]$gene_name)
(mito3=Complete_gene_names[grep("pseudogene", Complete_gene_names$description,ignore.case = T),]$gene_name)
# mito=unique(c(mito1,mito2))
mito=unique(c(mito1,mito2,mito3))
# mito=mito[order(mito)]
# mito[mito=="MTATP6P1"]
(ribo=Complete_gene_names[grep("ribosomal", Complete_gene_names$description,ignore.case = T),]$gene_name)
#
symbol_names=Complete_gene_names
symbol_names=symbol_names[symbol_names$gene_id%in%rownames(Test_DGE),]
symbol_names[duplicated(symbol_names$gene_name),]
symbol_names=symbol_names[!duplicated(symbol_names$gene_name),]
Test_DGE=Test_DGE[rownames(Test_DGE)%in%symbol_names$gene_id,]
rownames(Test_DGE)=symbol_names$gene_name[match(rownames(Test_DGE),symbol_names$gene_id)]
head(Test_DGE)

"many_mito"
sum(rownames(Test_DGE)%in%mito)
mitosums1=colSums(Test_DGE[rownames(Test_DGE)%in%mito1,])
mitosums2=colSums(Test_DGE[rownames(Test_DGE)%in%mito2,])
mitosums3=colSums(Test_DGE[rownames(Test_DGE)%in%mito,])
"many_ribo"
sum(rownames(Test_DGE)%in%ribo)
ribosums=colSums(Test_DGE[rownames(Test_DGE)%in%ribo,])
##############
allsums=colSums(Test_DGE)
dim(Test_DGE)
Test_DGE=Test_DGE[!(rownames(Test_DGE) %in% mito),]
Test_DGE=Test_DGE[!(rownames(Test_DGE) %in% ribo),]
dim(Test_DGE)
rownames(Test_DGE)[grepl("^AC[0-9]|^AL[0-9]",rownames(Test_DGE))]

sumsums=colSums(Test_DGE)
data=data.frame(all=allsums,mito=mitosums3,ribo=ribosums,left=sumsums)
head(data)
data$nimed="Plate2"

data$stim=Test_annot$Stimulus[match(rownames(data),Test_annot$BC_well)]
p1=ggplot(data,aes(y=all,x=nimed,color=stim))+geom_violin(color="black")+geom_beeswarm(cex = 2)+theme_classic()+labs(color="")+ggtitle("All UMIs")+scale_y_continuous(limits = c(0,150000))+xlab("")+ylab("")
p2=ggplot(data,aes(y=mito,x=nimed,color=stim))+geom_violin(color="black")+geom_beeswarm(cex = 2)+theme_classic()+labs(color="")+ggtitle("Mitochondiral UMIs")+scale_y_continuous(limits = c(0,150000))+xlab("")+ylab("")
p3=ggplot(data,aes(y=ribo,x=nimed,color=stim))+geom_violin(color="black")+geom_beeswarm(cex = 2)+theme_classic()+labs(color="")+ggtitle("Ribosomal UMIs")+scale_y_continuous(limits = c(0,150000))+xlab("")+ylab("")
p4=ggplot(data,aes(y=left,x=nimed,color=stim))+geom_violin(color="black")+geom_beeswarm(cex = 2)+theme_classic()+labs(color="")+ggtitle("Left after removal")+scale_y_continuous(limits = c(0,150000))+xlab("")+ylab("")
plot_grid(p1,p2,p3,p4,ncol = 4)

################################################################################################
colnames(Test_DGE)%in%Test_annot$BC_well
Test_annot$BC_well%in%colnames(Test_DGE)
Test_annot=Test_annot[Test_annot$BC_well%in%colnames(Test_DGE),]
Test_annot=Test_annot[mixedorder(Test_annot$BC_well),]
Test_DGE=Test_DGE[,mixedorder(colnames(Test_DGE))]
Test_annot$UMIs=colSums(Test_DGE)
Test_annot$Info=Test_annot$KO
Test_annot$Info[grepl("WT",Test_annot$Info)]="WT"
################################################################################################
################################################################################################
# Create Seurat object

Test_Dotto <- CreateSeuratObject(counts = Test_DGE, min.cells = 5, min.features = 100, project = "Test-Tcells")
Test_Dotto@meta.data<- cbind(Test_Dotto@meta.data,Test_annot[match(colnames(Test_DGE), Test_annot$BC_well),])
# QC and selecting cells for further analysis
# Test_Dotto[["percent.mt"]] <- PercentageFeatureSet(Test_Dotto, pattern = "^MT-")
VlnPlot(Test_Dotto, features = c("nFeature_RNA", "nCount_RNA","UMIs"), ncol = 3)
# Filtering
# Test_Dotto <- subset(Test_Dotto, subset =UMIs< 100000 & nFeature_RNA< 6000 & percent.mt< 10 & percent.ercc< 37)
# VlnPlot(Test_Dotto, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ercc","UMIs"), ncol = 5)
#################################################################################################################################
###############################################################################################################
###############################################################################################################
table(paste0(Test_Dotto$Stimulus,Test_Dotto$Info))
# DESEQ2 for fewer samples
Test_norm_deseq=Test_Dotto
# Test_norm_deseq=Test_norm_deseq[,Test_norm_deseq$Plate=="P2"]
Test_norm_deseq=Test_norm_deseq[,grepl("GATA3|TBX21|WT",Test_norm_deseq$KO)]
# Test_norm_deseq=Test_norm_deseq[,grepl("GATA3|TBX21",Test_norm_deseq$KO)]
Test_DGE=as.data.frame(as.matrix(Test_norm_deseq@assays$RNA@counts))
annot=as.data.frame(Test_norm_deseq@meta.data)
head(annot)
dds <- DESeqDataSetFromMatrix(countData = Test_DGE, colData = annot ,design = ~ Stimulus + Info)
##########################################################################################################
keep <- rowSums(counts(dds)) >= 10
sum(keep)
paste("number of genes kept -",sum(keep))
dim(dds)
dds <- dds[keep,]
# Doing the DESEQ
dds <- DESeq(dds)
#################################################################################################
# norm
barplot(colSums(counts(dds,normalized=F)))
barplot(colSums(counts(dds,normalized=T)))
# 
vsd <- vst(dds)
#
# vsd=varianceStabilizingTransformation(dds)

# Heatmap
################
vsd2=vsd[,grepl("GATA|TBX|WT",vsd$KO)]
################
select <- order(rowVars(assay(vsd,normalized=TRUE)),decreasing=TRUE)[1:500]
df <- as.data.frame(colData(vsd)[,c("Stimulus","Info","Origin_well")])
nii=assay(vsd)[select,]
# nii=nii[!grepl("^MT",rownames(nii)),]
xxx=rownames(nii)[grepl("GNLY|PRF|GZM",rownames(nii))]
top100=rownames(nii)[1:100]
x=c(xxx,top100)
nii=nii[x,]
# nii=nii[order(rowMeans(nii))[1:100],]
rownames(df)==colnames(nii)
annot$data=paste(annot$Stimulus,annot$Info,annot$Origin_well,sep = "_")
colnames(nii)=annot$data[match(colnames(nii),annot$BC_well)]
rownames(df)=colnames(nii)
pheatmap(nii, cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col=df)
###################################################### 
### PCA !!!
vsd1=vsd[,vsd$Origin_well!="H9"]
mat=assay(vsd1)
# mat=assay(vsd)
rowVar <- apply(mat, 1, var)
mv2000 <- order(rowVar, decreasing = T)[1:1000]
# mat <- t(mat)
mat <- t(mat[mv2000,])
pc <- prcomp(mat, scale = T)
pc.sum <- summary(pc)$importance
pc.sum[, 1:5]
df <- data.frame(pcs = 1:dim(pc.sum)[2], variance = pc.sum[2,], rule1 = pc.sum[3, ] > 0.8)
# ggplot(df, aes(x = pcs, y = variance, alpha = rule1)) + geom_bar(stat = "identity", col = 1) + 
#   geom_hline(yintercept = mean(df$variance),col = "red") + 
#   geom_segment(aes(x = 4, y = 0.08, xend = 4, yend = 0.04), arrow = arrow(), col = "darkred") + 
#   scale_x_discrete(breaks = 1:dim(pc.sum)[2])+  theme(legend.position = "None")
# always put the variances on the axes
varExp=round(pc.sum[2, ] * 100, 2)
pcs=data.frame(pc$x, stringsAsFactors = F)
pcs=cbind(pcs,annot[match(rownames(pcs),annot$BC_well),])
head(pcs)
away=pcs$Origin_well[pcs$PC2<(-5)]
away=c(away,"H9")
# pcs=pcs[pcs$Origin_well!="H9",]
x=ggplot(pcs, aes(PC1, PC2, color=Info,shape=Stimulus,label=Origin_well))+
  geom_point(size=4) + theme_classic()+ggtitle("PCA based on top 1000 HVGs")+
  xlab(paste("PC1 (", varExp[1], "%)"))+ ylab(paste("PC2 (",varExp[2], "%)"))+
    coord_fixed()+scale_color_manual(values = cc[c(6,1,3)])
x
# pcs$dataaa=paste0(pcs$KO,"_",pcs$Origin_well)
# x=ggplot(pcs, aes(PC1, PC2, color=Info, shape=Stimulus, label=dataaa))+
#   geom_point(size=4) + theme_classic()+ggtitle("PCA based on top 1000 HVGs")+
#   xlab(paste("PC1 (", varExp[1], "%)"))+ ylab(paste("PC2 (",varExp[2], "%)"))+
#   coord_fixed()

ggplotly(x)


table(paste0(pcs$Stimulus,pcs$Info))

########################################################################################
# CONTRAST
###############
# Gata vs Tbx everything
########################################################################################
dds$Stimulus
dds=dds[,!dds$Origin_well%in%away]
dds=DESeq(dds)
##############
res=as.data.frame(results(dds, contrast=c("Stimulus","Stimulated","Unstimulated")))
res=res[order(res$pvalue),]
res=na.omit(res)
head(res)
volcano_data=res
head(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[base::which(volcano_data['padj'] < 0.01 & abs(volcano_data['log2FoldChange']) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data['padj'] > 0.01 & abs(volcano_data['log2FoldChange']) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data['padj'] < 0.01 & abs(volcano_data['log2FoldChange']) > 1 ),"group"] <- "Significant&FoldChange"
pal <- c("blue", "black", "purple","red")
volcano_data$Symbol=rownames(volcano_data)
plot_ly(volcano_data, x=~log2FoldChange, y=~-log10(padj),
        color=~group,colors=pal,
        text=~Symbol,type="scatter",size = ~-log10(padj),
        mode="markers") %>%
  layout(title = 'Unstimulated vs. cGAMP + CD3/CD28',
         yaxis = list(zeroline = T),
         xaxis = list(zeroline = T))

volcano_data=as.data.frame(volcano_data)
volcano_data$txt=NA
volcano_data$txt[volcano_data$padj<1e-20]=volcano_data$Symbol[volcano_data$padj<1e-20]

v_stim=ggplot(volcano_data,aes(x=log2FoldChange,y=-log10(padj),color=group,label=txt))+geom_point(aes(size=-log10(padj)))+geom_text_repel(color="black")+theme_classic()+
  scale_color_manual(values = pal[c(1,2,3,4)])+geom_vline(xintercept = 0)
v_stim

xxxx=volcano_data
xxxx=xxxx[xxxx$log2FoldChange>2,]
xxxx=xxxx[xxxx$padj<0.01,]
xxxx
write.csv(xxxx,"/Home Office/moooooooooooooooooooo.csv")
###############
# CONTRAST
###############
# Gata vs Tbx everything
########################################################################################
Test_norm_deseq=Test_Dotto
Test_norm_deseq=Test_norm_deseq[,grepl("GATA3|TBX21|WT",Test_norm_deseq$KO)]
#
Test_norm_deseq=Test_norm_deseq[,grepl("Stimulated",Test_norm_deseq$Stimulus)]
Test_norm_deseq=Test_norm_deseq[,Test_norm_deseq$Origin_well!="C6"]
#
Test_norm_deseq=Test_norm_deseq[,grepl("Unstimulated",Test_norm_deseq$Stimulus)]
#
Test_DGE=as.data.frame(as.matrix(Test_norm_deseq@assays$RNA@counts))
annot=as.data.frame(Test_norm_deseq@meta.data)
dds2 <- DESeqDataSetFromMatrix(countData = Test_DGE, colData = annot ,design = ~ Info)
dds2 <- DESeq(dds2)
GATAvsTH1=as.data.frame(results(dds2, contrast=c("Info", "GATA3_KO", "TBX21_KO")))
GATAvsTH1=GATAvsTH1[order(GATAvsTH1$pvalue),]
GATAvsTH1=na.omit(GATAvsTH1)
head(GATAvsTH1)
volcano_data=GATAvsTH1
head(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[base::which(volcano_data['padj'] < 0.01 & abs(volcano_data['log2FoldChange']) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data['padj'] > 0.01 & abs(volcano_data['log2FoldChange']) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data['padj'] < 0.01 & abs(volcano_data['log2FoldChange']) > 1 ),"group"] <- "Significant&FoldChange"
pal <- c("blue", "black", "purple","red")
volcano_data$Symbol=rownames(volcano_data)
plot_ly(volcano_data, x=~log2FoldChange, y=~-log10(padj),
        color=~group,colors=pal,
        text=~Symbol,type="scatter",size = ~-log10(padj),
        mode="markers") %>%
  layout(title = 'Unstimulated vs. cGAMP + CD3/CD28',
         yaxis = list(zeroline = T),
         xaxis = list(zeroline = T))

volcano_data[volcano_data$Symbol=="IFNG",]
volcano_data=as.data.frame(volcano_data)
volcano_data$txt=NA
volcano_data$txt[volcano_data$padj<0.01]=volcano_data$Symbol[volcano_data$padj<0.01]

v_stim=ggplot(volcano_data,aes(x=log2FoldChange,y=-log10(padj),color=group,label=txt))+geom_point()+geom_text_repel()+theme_classic()+
  scale_color_manual(values = pal[c(3,2,4)])+geom_vline(xintercept = 0)
v_stim
plot_grid(v_stim,v_unstim)
###############

# Gata vs Tbx p2 only and only stim
########################################################################################
# dds2=dds[,dds$Plate=="P2"]
dds2=dds[,dds$Stimulus=="Stimulated"]
dds2=dds[,dds$Stimulus=="Unstimulated"]
design(dds2)= ~ Info
dds2=dds2[,dds2$Origin_well!="C6"]
dds2=dds2[,dds2$KO%in%c("GATA3_KO","TBX21_KO")]
dds2$Info=droplevels(dds2$Info)
dds2 <- DESeq(dds2)
resultsNames(dds2)
GATAvsTBX=as.data.frame(results(dds2, contrast=c("Info", "GATA3_KO", "TBX21_KO")))
# GATAvsTH1=as.data.frame(results(dds2, contrast=c("Stimulus", "Stimulated", "Unstimulated")))
GATAvsTBX=GATAvsTBX[order(GATAvsTBX$pvalue),]
GATAvsTBX=na.omit(GATAvsTBX)
head(GATAvsTBX)
volcano_data=GATAvsTBX
head(volcano_data)
volcano_data["group"] <- "NotSignificant"
volcano_data[base::which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) < 1 ),"group"] <- "Significant"
volcano_data[which(volcano_data['padj'] > 0.05 & abs(volcano_data['log2FoldChange']) > 1 ),"group"] <- "FoldChange"
volcano_data[which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) > 1 ),"group"] <- "Significant&FoldChange"
pal <- c("blue", "black", "purple","red")
volcano_data$Symbol=rownames(volcano_data)
plot_ly(volcano_data, x=~log2FoldChange, y=~-log10(padj),
        color=~group,colors=pal,
        text=~Symbol,type="scatter",size = ~-log10(padj),
        mode="markers") %>%
  layout(title = 'Unstimulated vs. cGAMP + CD3/CD28',
         yaxis = list(zeroline = T),
         xaxis = list(zeroline = T))

volcano_data[volcano_data$Symbol%in%c("IL4","IL5","IL6","IL10","IL13"),]


volcano_data=as.data.frame(volcano_data)
volcano_data$txt=NA
volcano_data$txt[volcano_data$padj<0.05]=volcano_data$Symbol[volcano_data$padj<0.05]
v2=ggplot(volcano_data,aes(x=log2FoldChange,y=-log10(padj),color=group,label=txt))+geom_point()+geom_text_repel(color="black")+theme_classic()+
  scale_color_manual(values = pal[c(3,2,4)])+geom_vline(xintercept = 0)+xlim(-6,6)
v2
plot_grid(v2,v1)
############################################################################
############################################################################
xyagain=assay(vsd)
# xyagain=counts(dds,normalized=T)
# xyagain=counts(dds,normalized=T)
xyagain=t(xyagain)
xyagain=xyagain[,c("IL13","IFNG","IL4","IL5","TNF","IL9")]
xyagain=as.data.frame(xyagain)
xyagain$name=rownames(xyagain)
xyagain$stim=annot$Stimulus[match(xyagain$name,annot$BC_well)]
xyagain$info=annot$Info[match(xyagain$name,annot$BC_well)]
xyagain$well=annot$Origin_well[match(xyagain$name,annot$BC_well)]
xyagain$plate=annot$Plate[match(xyagain$name,annot$BC_well)]
xyagain$ko=annot$KO[match(xyagain$name,annot$BC_well)]
xyagain=xyagain[!xyagain$info=="IRF2_HET",]
xyagain=xyagain[xyagain$stim=="Stimulated",]
xyagain=xyagain[!xyagain$well=="C6",]
xyagain=xyagain[grepl("GATA|TBX|WT",xyagain$info),]
xyagain$info=factor(xyagain$info,levels = sort(unique(xyagain$info)))
head(xyagain)

xyagain$IL13=2^xyagain$IL13
xyagain$IFNG=2^xyagain$IFNG
xyagain$IL4=2^xyagain$IL4
xyagain$IL5=2^xyagain$IL5
xyagain$TNF=2^xyagain$TNF
xyagain$IL9=2^xyagain$IL9

library(ggpubr)
p1=ggbarplot(xyagain, x = "info", y = "IFNG",add = c("mean_se", "jitter"),fill = "info")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("IFN-G")+ylab("")
p2=ggbarplot(xyagain, x = "info", y = "TNF",add = c("mean_se", "jitter"),fill = "info")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("TNF")+ylab("")
p3=ggbarplot(xyagain, x = "info", y = "IL4",add = c("mean_se", "jitter"),fill = "info")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("IL-4")+ylab("")
# p3=ggbarplot(xyagain, x = "info", y = "IL5",add = c("mean_se", "jitter"),fill = "info")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
#   xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("IL-5")+ylab("")
p4=ggbarplot(xyagain, x = "info", y = "IL13",add = c("mean_se", "jitter"),fill = "info")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("IL-13")+ylab("")
# p6=ggbarplot(xyagain, x = "info", y = "IL9",add = c("mean_se", "jitter"),fill = "info")+theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
#   xlab("")+scale_fill_manual(values = cc[c(6,1,3)])+ggtitle("IL-9")+ylab("")

plot_grid(p1,p2,p3,p4,ncol = 4)

# XY PLOT IF NEED BE
rna=ggplot(xyagain,aes(x=log2(IFNG), y = IL13,color=info,label=well))+ggrepel::geom_text_repel(color="black")+geom_point(size=2)+
  theme_classic()+scale_color_manual(values = cc[c(6,1,3)])
plot_grid(rna,prot)

head(xyagain)
pca_data=xyagain[c(1,2,3,5)]
head(pca_data)
rownames(pca_data)=xyagain$name
pca_data=prcomp(pca_data,scale. = T)
pc.sum=summary(pca_data)$importance
pc.sum
df <- data.frame(pcs = 1:dim(pc.sum)[2], variance = pc.sum[2,], rule1 = pc.sum[3, ] > 0.8)
varExp <- round(pc.sum[2, ] * 100, 2)
pcs <- data.frame(pca_data$x, stringsAsFactors = F)
pcs=cbind(pcs,xyagain[7:10])
head(pcs)
ggplot(pcs, aes(x = PC1, y = PC2,color=info,label=well )) + geom_point(size = 2,alpha=1)+theme_classic()+ 
  xlab(paste("PC1 (", varExp[1], "%)"))+  ylab(paste("PC2 (",varExp[2], "%)"))+scale_color_manual(values = cc[c(6,1,3)])+
  ggrepel::geom_text_repel(color="black")

