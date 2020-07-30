##############################################################################
## Changes in peripheral & local tumor immunity following NAC in BC         ##
## Data analysis by Justin M. Balko, Pharm.D., Ph.D. & Margaret L. Axelrod  ##
##                                                                          ##
##############################################################################

################################################################
######   Peripheral Blood NanoString Data Analysis   ###########
################################################################

#Figures 6, S13

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
library(rstatix)
library(Seurat)
library(survival)
library(survminer)
library(utils)
library(patchwork)

theme_set(theme_classic())


######Vanderbilt Cohort
#Load data ; normalized to positive/negative controls & endogenous housekeepers using nanoString nCounter software
read.delim("Normalizedcounts.txt", row.names=1, check.names=FALSE)->norm #Table S8
read.delim("metadata.txt", row.names=1)->metadata #TableS7
all(rownames(metadata)==colnames(norm)) #should be all true - data to metadata matching
metadata[!rownames(metadata) %in% c("31", "34"),]->foo
table(foo$NAC) 

#Transform Data
log(norm,2)->log.dataset
#log.dataset<-scaleRow(log.dataset) # z-score for sum 8 gene
t<-data.frame(t(log.dataset))
all(rownames(t)==rownames(metadata)) 
data<-cbind(t,metadata)
###Exlclude some samples: 31 has no recur data, 34 has low follow up time (200 days vs 1000 days for others)
data[!rownames(data) %in% c("31", "34"),]->data

#Calculate 8 gene score
data$combo<-((data$FGFBP2+data$GNLY+data$GZMB+data$LAG3+data$GZMH+data$NKG7+data$PDCD1)-data$HLA.G)

##HEATMAP
#Sort data in order by response.recur then by TNBC and then sub set
#Manually assigned patient order
read.delim("temp.txt", row.names=1, check.names=FALSE)->sorted
tsorted<-data.frame(t(sorted), check.names = FALSE)
tsorted[rownames(tsorted) %in% c("HLA.G","FGFBP2", "GNLY", "GZMB", "GZMH", "LAG3", "NKG7", "PDCD1", "combo"),]->sub
rownames(sub)<-c("FGFBP2", "GNLY", "GZMB", "GZMH", "HLA-G", "LAG3", "NKG7", "PDCD1", "8 Gene Score")

diverge_hcl(100,h=c(260,0),c=100,l = c(55, 100),power=1)->cols
colnames(sub)==colnames(tsorted)
colnames(tsorted)==rownames(sorted)
data.frame(sorted[["RR"]], sorted[["TNBC"]])->annData
colnames(annData)<-c("Response","TNBC")

#scale by row
scaleRow <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

foo<-scaleRow(data.matrix(sub))

f<-HeatmapAnnotation(df= annData, col = list(Response=c("No NAC" = "dodgerblue4","RD YES" = "deepskyblue", "RD NO" = "darkorchid3", "pCR" ="slateblue1" ),
                                             TNBC= c("NO" = "black", "YES"= "purple")),
                     annotation_legend_param = list(Response= list(labels=c("No NAC", "RD recur", "RD not recur", "pCR"))))

Heatmap(foo, col=cols, cluster_columns = FALSE, name = "",
        column_title = "Patient ID", column_title_side = "bottom", column_title_gp = gpar(fontsize= 10),
        column_split= factor(rep(annData$Response), levels= c("No NAC", "RD YES", "RD NO", "pCR")), cluster_column_slices = FALSE,
        top_annotation = f, row_order = c("8 Gene Score", "FGFBP2", "GNLY", "GZMB", "GZMH", "NKG7", "LAG3", "PDCD1", "HLA-G"),
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))


###########Plots for 8 genes by Response/Recur, with pCR as only one group, including stats############
#stats using Kruskal and post-hoc Dunn test
#individually add calculated kw to plot & locate significance brackets

maggieplt<-function(gene, ylab, p.pos, kw.pos){
  ggplot(data, aes(RR, gene)) +geom_boxplot(outlier.shape = NA)+
    geom_jitter(aes(colour=Response),width = 0.05, height = 0) +
    ylab(ylab)+ xlab("")+theme_classic() + 
    theme(legend.position = "none")+
    scale_x_discrete(limits= c("RD YES", "RD NO", "pCR", "No NAC"), labels = c("RD recur", "RD not recur", "pCR","No NAC"))+
    geom_vline(xintercept = 3.5, linetype="dashed")+
    stat_pvalue_manual(dunn.sub, y.position = c(p.pos))+
    annotate("text",x=1.8, y=kw.pos, label=paste("Kruskal-Wallis, p = ", kw$p))
}
maggieplt2<-function(gene, ylab,kw.pos){
  ggplot(data, aes(RR, gene)) +geom_boxplot(outlier.shape = NA)+
    geom_jitter(aes(colour=Response),width = 0.05, height = 0) +
    ylab(ylab)+ xlab("")+theme_classic() + 
    theme(legend.position = "none")+
    scale_x_discrete(limits= c("RD YES", "RD NO", "pCR", "No NAC"), labels = c("RD recur", "RD not recur", "pCR","No NAC"))+
    geom_vline(xintercept = 3.5, linetype="dashed")+
    annotate("text",x=1.8, y=kw.pos, label=paste("Kruskal-Wallis, p = ", kw$p))
} #no post-hoc dunn for those without significant KW

#8 Gene Score
kw<-kruskal_test(formula = combo~RR, data=data )
dunn<-dunn_test(data = data, formula = combo ~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]
maggieplt(data$combo,expression(paste("Log Normalized  8 Gene Score")), 80, 38)
maggieplt(data$combo,expression(paste("8 Gene Score (sum Z-scores)")), 17, -18)

#FGFBP2
kw<-kruskal_test(formula = FGFBP2~RR, data=data )
dunn<-dunn_test(data = data, formula = FGFBP2~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]
maggieplt(data$FGFBP2,expression(paste("Log Normalized  ", italic(FGFBP2), " Counts")), 13, 5.9)

#GNLY
kw<-kruskal_test(formula = GNLY~RR, data=data )
dunn<-dunn_test(data = data, formula = GNLY~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]
maggieplt(data$GNLY, expression(paste("Log Normalized  ", italic(GNLY), " Counts")), 16, 11)

#GZMB
kw<-kruskal_test(formula = GZMB~RR, data=data )
dunn<-dunn_test(data = data, formula = GZMB~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]
maggieplt(data$GZMB, expression(paste("Log Normalized  ", italic(GZMB), " Counts")), 14.2, 7.8)
maggieplt2(data$GZMB, expression(paste("Log Normalized  ", italic(GZMB), " Counts")), 7.8)


#GZMH
kw<-kruskal_test(formula = GZMH~RR, data=data )
dunn<-dunn_test(data = data, formula = GZMH~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]
maggieplt(data$GZMH, expression(paste("Log Normalized  ", italic(GZMH), " Counts")), 14, 6.7)
maggieplt2(data$GZMH, expression(paste("Log Normalized  ", italic(GZMH), " Counts")), 6.7)

#LAG3
kw<-kruskal_test(formula = LAG3~RR, data=data )
dunn<-dunn_test(data = data, formula = LAG3~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]
maggieplt(data$LAG3, expression(paste("Log Normalized  ", italic(LAG3), " Counts")), 11, 5.5)


#HLA-G
kw<-kruskal_test(formula = HLA.G~RR, data=data )
dunn<-dunn_test(data = data, formula = HLA.G~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]
maggieplt(data$HLA.G, expression(paste("Log Normalized  ", italic(HLA-G), " Counts")), 14, 10)
maggieplt2(data$HLA.G, expression(paste("Log Normalized  ", italic(HLA-G), " Counts")), 10)

#NKG7
kw<-kruskal_test(formula = NKG7~RR, data=data )
dunn<-dunn_test(data = data, formula = NKG7~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]

#PDCD1
kw<-kruskal_test(formula = PDCD1~RR, data=data )
dunn<-dunn_test(data = data, formula = PDCD1~RR)
dunn.sub<-dunn[rownames(dunn) %in% c("5"),]
maggieplt(data$PDCD1, expression(paste("Log Normalized  ", italic(PDCD1), " Counts")), 9, 2)



############DFCI Cohort#################
#Load data ; normalized to positive/negative controls & endogenous housekeepers using nanoString nCounter software
#Data available on request
read.delim("normalizedcounts.txt", row.names=1, check.names=FALSE)->norm
read.delim("metadata.txt", row.names=1)->metadata
all(rownames(metadata)==colnames(norm))

#Transform
log(norm,2)->log.dataset
z<-scaleRow(log.dataset)
t<-data.frame(t(log.dataset))
z.t<-data.frame(t(z))
all(rownames(t)==rownames(metadata))
all(rownames(z.t)==rownames(metadata))
data<-cbind(t,metadata)
z.data<-cbind(z.t, metadata)

#Calculate 8 gene score
data$combo<-((data$FGFBP2+data$GNLY+data$GZMB+data$LAG3+data$GZMH+data$NKG7+data$PDCD1)-data$HLA.G)
z.data$combo<-((z.data$FGFBP2+z.data$GNLY+z.data$GZMB+z.data$LAG3+z.data$GZMH+z.data$NKG7+z.data$PDCD1)-z.data$HLA.G)

#Combine categories
z.data$RR<-paste(z.data$RCB, z.data$Recur)
data$RR<-paste(data$RCB, data$Recur)

z.data$RR<-revalue(z.data$RR, c("III YES"= "III", "III NO"= "III", "III NA"= "III"))
data$RR<-revalue(data$RR, c("III YES"= "III", "III NO"= "III", "III NA"= "III"))

z.data$RRcombine<-revalue(z.data$RR, c("pCR NO"= "0/I/II NO", "I NO"="0/I/II NO", "II NO" ="0/I/II NO", 
                                       "I YES" = "0/I/II YES", "II YES" = "0/I/II YES",
                                       "I NA" = "NA", "II NA"= "NA", "III NA" = "NA", "NA NA" = "NA", "NA NO"= "NA", "pCR NA"= "NA"))

data$RRcombine<-revalue(data$RR, c("pCR NO"= "0/I/II NO", "I NO"="0/I/II NO", "II NO" ="0/I/II NO", 
                                   "I YES" = "0/I/II YES", "II YES" = "0/I/II YES",
                                   "I NA" = "NA", "II NA"= "NA", "III NA" = "NA", "NA NA" = "NA", "NA NO"= "NA", "pCR NA"= "NA"))


z.data[z.data$Timepoint=="P-S",]->ps.z
data[data$Timepoint=="P-S",]->ps
ps.z.clean<-na.omit(ps.z)
ps.clean<-na.omit(ps)

ggplot(ps.z, aes(RCB, combo)) +geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.05, height = 0)+ ylab("8 Gene Score")


kw<-kruskal_test(data=ps.z.clean, formula= combo~RRcombine)
stat<-dunn_test(data=ps.z.clean, formula= combo~RRcombine)
sub<-stat[stat$p.adj<0.05,]

ggplot(ps.z.clean, aes(RRcombine, combo)) +geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.05, height = 0)+ ylab("8 Gene Score (sum Z-scores)")+
  scale_x_discrete(labels=c("RCB0/I/II not recur", "RCB0/I/II recur", "RCBIII"))+
  xlab("")+theme(axis.text.x=element_text(angle=90, size=12, color="black"),
                 axis.title.y= element_text(size=13))+
  annotate("text", label= paste("Kruskal-Wallis, p=", kw$p ), x=2, y=-9)+
  stat_pvalue_manual(sub, y.position = 16)


# Individual gene plots for RCB w/ outcome data for DFCI SAMPLES#
dfci<-function(gene, y, p.y){
  ggplot(ps.clean, aes(RRcombine, gene)) +geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.05, height = 0)+ ylab(y)+
    scale_x_discrete(labels=c("RCB0/I/II not recur", "RCB0/I/II recur", "RCBIII"))+
    xlab("")+theme(axis.text.x=element_text(angle=90, size=12, color="black"),
                   axis.title.y= element_text(size=13))+
    annotate("text", label= paste("Kruskal-Wallis, p=", kw$p ), x=2, y=p.y)
}

kw<-kruskal_test(data=ps.clean, formula= FGFBP2~RRcombine)
stat<-dunn_test(data=ps.clean, formula= FGFBP2~RRcombine)
sub<-stat[stat$p.adj<0.05,]
dfci(ps.clean$FGFBP2, expression(paste("Log Normalized  ", italic(FGFBP2), " Counts")) ,4)

kw<-kruskal_test(data=ps.clean, formula= GNLY~RRcombine)
stat<-dunn_test(data=ps.clean, formula= GNLY~RRcombine)
sub<-stat[stat$p.adj<0.05,]
dfci(ps.clean$GNLY, expression(paste("Log Normalized  ", italic(GNLY), " Counts")) ,7)

kw<-kruskal_test(data=ps.clean, formula= GZMB~RRcombine)
stat<-dunn_test(data=ps.clean, formula= GZMB~RRcombine)
sub<-stat[stat$p.adj<0.05,]
dfci(ps.clean$GZMB, expression(paste("Log Normalized  ", italic(GZMB), " Counts")) ,5)

kw<-kruskal_test(data=ps.clean, formula= GZMH~RRcombine)
stat<-dunn_test(data=ps.clean, formula= GZMH~RRcombine)
sub<-stat[stat$p.adj<0.05,]
dfci(ps.clean$GZMH, expression(paste("Log Normalized  ", italic(GZMH), " Counts")) ,4.4)+  
  stat_pvalue_manual(sub, y.position = 10.1)

kw<-kruskal_test(data=ps.clean, formula= HLA.G~RRcombine)
stat<-dunn_test(data=ps.clean, formula= HLA.G~RRcombine)
sub<-stat[stat$p.adj<0.05,]
dfci(ps.clean$HLA.G, expression(paste("Log Normalized  ", italic(HLA-G), " Counts")) ,6)  

kw<-kruskal_test(data=ps.clean, formula= LAG3~RRcombine)
stat<-dunn_test(data=ps.clean, formula= LAG3~RRcombine)
sub<-stat[stat$p.adj<0.05,]
dfci(ps.clean$LAG3, expression(paste("Log Normalized  ", italic(LAG3), " Counts")) ,4.8)  

kw<-kruskal_test(data=ps.clean, formula= NKG7~RRcombine)
stat<-dunn_test(data=ps.clean, formula= NKG7~RRcombine)
sub<-stat[stat$p.adj<0.05,]
dfci(ps.clean$NKG7, expression(paste("Log Normalized  ", italic(NKG7), " Counts")) ,9)  

kw<-kruskal_test(data=ps.clean, formula= PDCD1~RRcombine)
stat<-dunn_test(data=ps.clean, formula= PDCD1~RRcombine)
sub<-stat[stat$p.adj<0.05,]
dfci(ps.clean$PDCD1, expression(paste("Log Normalized  ", italic(PDCD1), " Counts")) ,1)  
