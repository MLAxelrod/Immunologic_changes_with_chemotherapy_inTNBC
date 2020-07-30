##############################################################################
## Changes in peripheral & local tumor immunity following NAC in BC         ##
## Data analysis by Justin M. Balko, Pharm.D., Ph.D. & Margaret L. Axelrod  ##
##                                                                          ##
##############################################################################

###################################################
############scRNAseq on PBMCs #####################
###################################################

#Analysis of scRNAseq on PD-1Hi CD8+ T cells and whole PBMCs
#Figures 5, S11, S12


#####################
##   Libraries     ##
#####################
library(circlize)
library(colorspace)
library(ComplexHeatmap)
library(cowplot)
library(plyr)
library(dplyr)
library(EnhancedVolcano)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggfortify)
library(gplots)
library(gtools)
library(heatmap3)
library(RColorBrewer)
library(rstatix)
library(Seurat)
library(survival)
library(survminer)
library(utils)
library(harmony)
library(hdf5r)
library(reticulate)
library(umap)
library(SeuratWrappers)
library(patchwork)
library(SingleR)

############ PD-1 Hi CD8+ T cell data #######

# Data input #Data available on request
sample.1.in <- Read10X_h5('Path/Sample1PD1hi/filtered_gene_bc_matrices_h5.h5')
sample.1 <- CreateSeuratObject(sample.1.in, project = 'Sample 1', min.cells = 5, min.features = 100)
sample.2.in <- Read10X_h5('Path/Sample2PD1hi/filtered_gene_bc_matrices_h5.h5')
sample.2 <- CreateSeuratObject(sample.2.in, project = 'Sample 2', min.cells = 5, min.features = 100)
merge.samples <- merge(sample.1, sample.2, add.cell.ids = c('S1', 'S2'), project="TNBC")

# Normalization and correction
merge.samples <- CellCycleScoring(merge.samples, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
merge.samples <- PercentageFeatureSet(merge.samples, pattern = "^MT-", col.name = "percent.mt")
merge.samples <- SCTransform(merge.samples, vars.to.regress = c('Phase', 'percent.mt'))
merge.samples <- RunPCA(merge.samples)
merge.samples <- JackStraw(merge.samples, reduction = 'pca')
merge.samples <- ScoreJackStraw(merge.samples)

JackStrawPlot(merge.samples)
Idents(merge.samples)<-"orig.ident"

merge.samples <- RunHarmony(merge.samples, group.by.vars = c('orig.ident'), max.iter.harmony = 20)

# Clustering, reduction for visualization, imputation
merge.samples <- FindNeighbors(merge.samples, reduction = 'harmony')
merge.samples <- FindClusters(merge.samples, resolution=0.3, graph.name = 'SCT_nn')
merge.samples <- RunUMAP(merge.samples, dims=1:5)#, reduction='harmony')
merge.samples <- RunALRA(merge.samples)
markers<-FindAllMarkers(merge.samples, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top50<-markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top15<-markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#Plots
FeaturePlot(merge.samples, features = c("CD8A", "CD4", "PDCD1"),split.by="orig.ident", cols=c("grey","red"))
FeaturePlot(merge.samples, features = c("GNLY", "GZMB", "FGFBP2"),split.by="orig.ident", cols=c("grey","red"))
FeaturePlot(merge.samples, features = c("SERPINB9"),split.by="orig.ident", cols=c("grey","red"))
VlnPlot(merge.samples, 'GNLY', split.by = 'orig.ident')
table(merge.samples@active.ident)

##Heatmaps##
merge.samples <- ScaleData(merge.samples, assay="alra")

DoHeatmap(merge.samples, features=top10$gene,slot="scale.data",cells = NULL, assay="alra", label=FALSE)+
  scale_fill_gradientn(colors = c("blue", "white", "red")) + labs(color="Cluster")+
  theme(axis.text.y = element_text(color = "grey20", size = 9, angle = 0, hjust = 1, vjust = 0, face = "plain"))

DoHeatmap(merge.samples, features=top50$gene,slot="scale.data",cells = NULL, group.by = "ident", assay="alra")+
  scale_fill_gradientn(colors = c("blue", "white", "red"))
DoHeatmap(merge.samples, features=markers$gene[markers$cluster==0][1:100],slot="scale.data",cells = NULL, group.by = "ident", assay="alra")+
  scale_fill_gradientn(colors = c("blue", "white", "red"))+
  theme(axis.text.y = element_text(color = "grey20", size = 7, angle = 0, hjust = 1, vjust = 0, face = "plain"))

c("TBX21","GZMB","GNLY","FASLG","IFNG","HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB5", "PDCD1", "LAG3","FCRL6")->genes
DoHeatmap(merge.samples, features=genes , group.by = "ident", assay="alra", label=FALSE)+
  scale_fill_gradientn(colors = c("blue", "white", "red"))+ labs(color="Cluster")+
  theme(axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"))


###############Volcano Plots for PD-1Hi CD8+ T cells##########################
Idents(merge)<-"orig.ident"
diff_sample <- FindMarkers(merge, ident.1 = "Sample 2", ident.2 = "Sample 1", assay= "alra")
Idents(merge)<-"seurat_clusters"
diff_cluster <- FindMarkers(merge, ident.1 = 0, assay= "alra")

#nanostring genes
ns<-c("FGFBP2", "HLA-DRB5", "HLA-G", "LAG3", "PDCD1", "NKG7", "GNLY", "GZMB", "GZMH")

EnhancedVolcano(diff_sample, lab = rownames(diff_sample), x= 'avg_logFC', y= 'p_val_adj',
                FCcutoff = 0.5, ylab = bquote(~-Log[10]~adjusted~italic(P)),
                legendPosition = 'top', legendLabSize = 10, legendIconSize = 2,
                drawConnectors = TRUE, selectLab = ns, transcriptLabSize = 4, axisLabSize =11,
                title="", subtitle = "Sample 1 vs. Sample 2", gridlines.major = FALSE, gridlines.minor = FALSE, caption="")

EnhancedVolcano(diff_cluster, lab = rownames(diff_cluster), x= 'avg_logFC', y= 'p_val_adj',
                FCcutoff = 0.5, ylab = bquote(~-Log[10]~adjusted~italic(P)),
                legendPosition = 'top', legendLabSize = 10, legendIconSize = 2,
                drawConnectors = TRUE, selectLab = ns, transcriptLabSize = 4, axisLabSize =11,
                title="", subtitle = "Enriched in Cluster 0", gridlines.major = FALSE, gridlines.minor = FALSE, caption="")


########## Whole Blood Single Cell RNA Seq ##############################
#Pt1 = pt 3399; metaplastic TNBC with pCR
#Pt2 = pt 2020-01; TNBC with pCR, newly collected

#Import data
s1<-Read10X("Path/pt3399/filtered_feature_bc_matrix/")
s2<-Read10X("Path/pt202001/filtered_feature_bc_matrix/")
pt1<- CreateSeuratObject(s1, project = 'pt1', min.cells = 5, min.features = 100)
pt2<- CreateSeuratObject(s2, project = 'pt2', min.cells = 5, min.features = 100)
#QC
pt1[["percent.mt"]] <- PercentageFeatureSet(pt1, pattern = "^MT-")
pt2[["percent.mt"]] <- PercentageFeatureSet(pt2, pattern = "^MT-")
VlnPlot(pt1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dim(pt2)
pt1<- subset(pt1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
pt2<- subset(pt2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
pt1<-NormalizeData(pt1)
pt2<-NormalizeData(pt2)
pt1 <- FindVariableFeatures(pt1, selection.method = "vst", nfeatures = 2000)
pt2 <- FindVariableFeatures(pt2, selection.method = "vst", nfeatures = 2000)
#Integrate
anchors<-FindIntegrationAnchors(object.list= list(pt1, pt2), dims=1:20)
combo<-IntegrateData(anchorset = anchors, dims=1:20)
DefaultAssay(combo)<-"integrated"
#Transform Data
combo<-ScaleData(combo)
combo <- CellCycleScoring(combo, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
combo <- SCTransform(combo, vars.to.regress = c('Phase', 'percent.mt'))
combo<-RunPCA(combo, npcs=30)

#Chose # of PCs
combo<- JackStraw(combo, num.replicate = 100, assay="SCT")
combo<-ScoreJackStraw(combo, dims = 1:20)
JackStrawPlot(combo, dims = 1:20)
ElbowPlot(combo) #Set to 15

combo<- RunHarmony(combo, group.by.vars = c('orig.ident'), max.iter.harmony = 20, assay.use = "SCT")

combo<-RunALRA(combo, assay="SCT")

#UMAP and Cluster
combo<-RunUMAP(combo, reduction="harmony", dims=1:15,
               n.neighbors = 15, seed.use = 42, n.epochs = 500, min.dist = 0.5)
combo<-FindNeighbors(combo, reduction = "harmony", dims=1:15)
combo<-FindClusters(combo, resolution = 0.2, graph.name = "SCT_nn")
#Plot
DimPlot(combo, reduction="umap")
Idents(combo)<-"seurat_clusters"
DimPlot(combo, reduction="umap", split.by = "orig.ident")
DimPlot(combo, reduction="umap", group.by = "orig.ident")

#SINGLER
hc<-GetAssayData(combo, slot= "data", assay="alra")
myref<-BlueprintEncodeData()

pred<-SingleR(test=hc, ref=myref, method="single", labels = myref$label.main)

plotScoreHeatmap(pred)
to.remove<-pruneScores(pred)
plotScoreDistribution(pred,show.nmads = 3)

labs<-pred$labels
labs<-revalue(labs, c("Erythrocytes"="Other", "DC"="Other", "Eosinophils"="Other", "HSC"="Other", "Neutrophils"="Other"))
#labs[is.na(labs)]<-"Other"

combo[["CellTypes"]]<-labs
Idents(combo)<-"CellTypes"
DimPlot(combo, reduction = "umap", label=TRUE)

plt<-function(gene){
  VlnPlot(combo, features = gene, split.by = "orig.ident", group.by="CellTypes", pt.size = 0, sort=FALSE, log=FALSE)+
    scale_x_discrete(limits=c("NK cells", "CD8+ T-cells", "CD4+ T-cells", "B-cells", "Monocytes", "Other"))+
    ggtitle(gene)+xlab("")
}

plt("CD3E")
plt("CD8A")
plt("CD4")
plt("CD19")
plt("CD14")

plt("FGFBP2")
plt("GNLY")
plt("GZMB")
plt("GZMH")
plt("NKG7")
plt("HLA-DRB5")
plt("HLA-G")
plt("LAG3")
plt("PDCD1")

###Gene Set Score###
gene.set<-c("FGFBP2", "GNLY", "GZMB", "GZMH", "NKG7","LAG3", "PDCD1")

# Get mean expression of genes of interest per cell
mean.exp <- colMeans(x = combo[gene.set, ], na.rm = TRUE)-colMeans(x = combo["HLA-G", ], na.rm = TRUE)

# Add mean expression values in 'object@meta.data$gene.set.score'
if (all(names(x = mean.exp) == rownames(x = combo@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'combo@meta.data':\n", 
      "adding gene set mean expression values in 'combo@meta.data$gene.set.score'")
  combo@meta.data$gene.set.score <- mean.exp
}

# Plot mean expression using Seurat::FeaturePlot()
FeaturePlot(object = combo, features = "gene.set.score", cols=c("darkgrey", "grey","orangered","red")) + ggtitle("8 Gene Score")

plt("gene.set.score")+ ggtitle("8 Gene Score")

VlnPlot(combo, features = c("gene.set.score"),split.by="orig.ident", group.by="CellTypes",pt.size =0, sort=FALSE, log=FALSE)+
  scale_x_discrete(limits=c("NK cells", "CD8+ T-cells", "CD4+ T-cells", "B-cells", "Monocytes", "Other"))+ ggtitle("8 Gene Score")+xlab("")+ 
  geom_boxplot(width=0.1, outlier.shape = NA)+scale_fill_discrete(labels=c("Pt 4","Pt 5"))

