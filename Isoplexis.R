##############################################################################
## Changes in peripheral & local tumor immunity following NAC in BC         ##
## Data analysis by Justin M. Balko, Pharm.D., Ph.D. & Margaret L. Axelrod  ##
##                                                                          ##
##############################################################################

###########################################################
###    Single cell cytokine profiling with ISOPLEXIS    ###
###########################################################

#Figures 4, S9

#####################
##   Libraries     ##
#####################
library(circlize)
library(colorspace)
library(ComplexHeatmap)
library(cowplot)
library(plyr)
library(dplyr)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggfortify)
library(gplots)
library(gtools)
library(RColorBrewer)
library(rstatix)
library(patchwork)

##Pt 1 = 3388; pt2 = 3349 ; pt3 = 3400; pt4=3399

theme_set(theme_classic())

#Raw data available on request

##Polyfunctionality Stacked Bar Plots
setwd("C:/Users/axelroml/Desktop/Isoplexis/CD8")
read.delim("pf_stacked.csv", row.names=1, check.names=TRUE, sep = ",")->pfcd8
setwd("C:/Users/axelroml/Desktop/Isoplexis/CD4")
read.delim("CD4pf_2stacked.csv", row.names=1, check.names=TRUE, sep = ",")->pfcd4
##Just 2+ Cd8+
ggplot(pfcd8, aes(External.Donor))+
  geom_col(aes(External.Donor, Percent, fill=factor(Analytes, levels=c("5+ Analytes", "4 Analytes", "3 Analytes", "2 Analytes"))), position="stack")+
  labs(fill="Cytokines \nSecreted")+ xlab(bquote('PD-1'^'Hi'~ 'CD8+'))+ylab("Polyfunctionality (% of sample)")+
  scale_x_discrete(limits=c("3388_Pre", "3388_Post", "3349_Pre", "3349_Post", "3400_Pre", "3400_Post", "3399_Pre", "3399_Post"),
                   labels=c("Pt1 Pre", "Pt1 Post", "Pt2 Pre", "Pt2 Post", "Pt3 Pre", "Pt3 Post", "Pt4 Pre", "Pt4 Post"))+
  annotate("text", x=1.5, y=38, label="TNBC \nRCB-II", size=5)+ geom_vline(xintercept=2.5, linetype="dashed")+
  annotate("text", x=3.5, y=38, label="ER+ \nRCB-II", size=5)+ geom_vline(xintercept=4.5, linetype="dashed")+
  annotate("text", x=5.5, y=38, label="ER+HER2+ \nRCB-III", size=5)+ geom_vline(xintercept=6.5, linetype="dashed")+
  annotate("text", x=7.5, y=38, label="TNBC \nRCB-0", size=5)+ theme(axis.text.x= element_text(size=14, angle=45, vjust=0.6),
                                                                     axis.title= element_text(size=15),
                                                                     legend.text = element_text(size=13),
                                                                     legend.title = element_text(size=15))+
  scale_fill_npg(labels = c("5+", "4", "3", "2"))

##Just 2+ Cd4+
ggplot(pfcd4, aes(External.Donor))+
  geom_col(aes(External.Donor, Percent, fill=factor(Analytes, levels=c("5+ Analytes", "4 Analytes", "3 Analytes", "2 Analytes"))), position="stack")+
  labs(fill="Cytokines \nSecereted")+ xlab(bquote('PD-1'^'Hi'~ 'CD4+'))+ylab("Polyfunctionality (% of sample)")+
  scale_x_discrete(limits=c("3388_Pre", "3388_Post", "3349_Pre", "3349_Post", "3400_Pre", "3400_Post", "3399_Pre", "3399_Post"),
                   labels=c("Pt1 Pre", "Pt1 Post", "Pt2 Pre", "Pt2 Post", "Pt3 Pre", "Pt3 Post", "Pt4 Pre", "Pt4 Post"))+
  annotate("text", x=1.5, y=38, label="TNBC \nRCB-II", size=5)+ geom_vline(xintercept=2.5, linetype="dashed")+
  annotate("text", x=3.5, y=38, label="ER+ \nRCB-II", size=5)+ geom_vline(xintercept=4.5, linetype="dashed")+
  annotate("text", x=5.5, y=38, label="ER+HER2+ \nRCB-III", size=5)+ geom_vline(xintercept=6.5, linetype="dashed")+
  annotate("text", x=7.5, y=38, label="TNBC \nRCB-0", size=5)+ theme(axis.text.x= element_text(size=14, angle=45, vjust=0.6),
                                                                     axis.title= element_text(size=15),
                                                                     legend.text = element_text(size=13),
                                                                     legend.title = element_text(size=15))+
  scale_fill_npg(labels = c("5+", "4", "3", "2"))


#Load data for heatmaps of cytokines per single cell
setwd("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/Isoplexis/CD8")
read.delim("PD-1+T cell pre,post_CD8+_raw_data.csv", row.names=1, check.names=TRUE, sep = ",")->cd8
setwd("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/Isoplexis/CD4")
read.delim("PD-1+T cell pre,post_CD4+_raw_data.csv", row.names=1, check.names=TRUE, sep = ",")->cd4

#Subset Patients & log transform
#Pt 1 = 3388; pt2 = 3349 ; pt3 = 3400; pt4=3399
#CD8s
cd8[cd8$External.Donor=="3388_Post" | cd8$External.Donor=="3388_Pre",]->pt1.cd8
pt1.cd8labels<- pt1.cd8[,"External.Donor"]
pt1.cd8[,!colnames(pt1.cd8) %in% c("Cell.Subset", "Stimulation", "External.Donor")] -> pt1.cd8
log(pt1.cd8+1,2)->log.pt1.cd8 #need plus 1 to keep zeros

cd8[cd8$External.Donor=="3349_Post" | cd8$External.Donor=="3349_Pre",]->pt2.cd8
pt2.cd8labels<- pt2.cd8[,"External.Donor"]
pt2.cd8[,!colnames(pt2.cd8) %in% c("Cell.Subset", "Stimulation", "External.Donor")] -> pt2.cd8
log(pt2.cd8+1,2)->log.pt2.cd8 #need plus 1 to keep zeros

cd8[cd8$External.Donor=="3400_Post" | cd8$External.Donor=="3400_Pre",]->pt3.cd8
pt3.cd8labels<- pt3.cd8[,"External.Donor"]
pt3.cd8[,!colnames(pt3.cd8) %in% c("Cell.Subset", "Stimulation", "External.Donor")] -> pt3.cd8
log(pt3.cd8+1,2)->log.pt3.cd8 #need plus 1 to keep zeros

cd8[cd8$External.Donor=="3399_Post" | cd8$External.Donor=="3399_Pre",]->pt4.cd8
pt4.cd8labels<- pt4.cd8[,"External.Donor"]
pt4.cd8[,!colnames(pt4.cd8) %in% c("Cell.Subset", "Stimulation", "External.Donor")] -> pt4.cd8
log(pt4.cd8+1,2)->log.pt4.cd8 #need plus 1 to keep zeros

#CD4s
cd4[cd4$External.Donor=="3388_Post" | cd4$External.Donor=="3388_Pre",]->pt1.cd4
pt1.cd4labels<- pt1.cd4[,"External.Donor"]
pt1.cd4[,!colnames(pt1.cd4) %in% c("Cell.Subset", "Stimulation", "External.Donor")] -> pt1.cd4
log(pt1.cd4+1,2)->log.pt1.cd4

cd4[cd4$External.Donor=="3349_Post" | cd4$External.Donor=="3349_Pre",]->pt2.cd4
pt2.cd4labels<- pt2.cd4[,"External.Donor"]
pt2.cd4[,!colnames(pt2.cd4) %in% c("Cell.Subset", "Stimulation", "External.Donor")] -> pt2.cd4
log(pt2.cd4+1,2)->log.pt2.cd4

cd4[cd4$External.Donor=="3400_Post" | cd4$External.Donor=="3400_Pre",]->pt3.cd4
pt3.cd4labels<- pt3.cd4[,"External.Donor"]
pt3.cd4[,!colnames(pt3.cd4) %in% c("Cell.Subset", "Stimulation", "External.Donor")] -> pt3.cd4
log(pt3.cd4+1,2)->log.pt3.cd4

cd4[cd4$External.Donor=="3399_Post" | cd4$External.Donor=="3399_Pre",]->pt4.cd4
pt4.cd4labels<- pt4.cd4[,"External.Donor"]
pt4.cd4[,!colnames(pt4.cd4) %in% c("Cell.Subset", "Stimulation", "External.Donor")] -> pt4.cd4
log(pt4.cd4+1,2)->log.pt4.cd4


#Heatmaps
col_fun= colorRamp2(c(0,.1, 15), c("white", "chocolate1", "red"))

#Pt1 CD8
hm.pt1.cd8<-data.matrix(log.pt1.cd8)
labs.pt1cd8<-data.matrix(pt1.cd8labels)

Heatmap(hm.pt1.cd8, cluster_columns= FALSE, cluster_rows = FALSE,col=col_fun, border= TRUE, show_row_names = FALSE,
        column_names_gp = gpar(fontsize=10), column_labels = c("GzmB", "IFNg", "MIP1a", "MIP1b", "TNFa", "TNFb", "sCD137"),
        row_split = factor(rep(labs.pt1cd8), levels = (c("3388_Pre", "3388_Post"))), row_title = c("Pt1 Post", "Pt1 Pre"), row_title_gp = gpar(fontsize=12),
        heatmap_legend_param = list(title="Log \nIntensity"))

#Pt2 CD8
hm.pt2.cd8<-data.matrix(log.pt2.cd8)
labs.pt2cd8<-data.matrix(pt2.cd8labels)

Heatmap(hm.pt2.cd8, cluster_columns= FALSE, cluster_rows = FALSE,  col=col_fun, border= TRUE, show_row_names = FALSE,
        column_labels = c("GzmB", "IFNg", "MIP1a", "MIP1b", "TNFa", "TNFb", "sCD137"), column_names_gp = gpar(fontsize=10),
        row_split = factor(rep(labs.pt2cd8), levels= c("3349_Pre", "3349_Post")), row_title = c("Pt2 Post", "Pt2 Pre"), row_title_gp = gpar(fontsize=12),
        heatmap_legend_param = list(title="Log \nIntensity"))

#Pt3 CD8
hm.pt3.cd8<-data.matrix(log.pt3.cd8)
labs.pt3cd8<-data.matrix(pt3.cd8labels)

Heatmap(hm.pt3.cd8, cluster_columns= FALSE, cluster_rows = FALSE, col=col_fun, border= TRUE, show_row_names = FALSE,
        column_labels = c("GzmB", "IFNg", "MIP1a", "MIP1b", "TNFa", "TNFb", "sCD137"), column_names_gp = gpar(fontsize=10),
        row_split = factor(rep(labs.pt3cd8), levels = c("3400_Pre","3400_Post")),row_title = c("Pt3 Post", "Pt3 Pre"), row_title_gp = gpar(fontsize=12),
        heatmap_legend_param = list(title="Log \nIntensity"))

#Pt4 CD8
hm.pt4.cd8<-data.matrix(log.pt4.cd8)
labs.pt4cd8<-data.matrix(pt4.cd8labels)

setwd("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/Isoplexis/CD8")
pdf("test2.pdf", width=3, height=20)
Heatmap(hm.pt4.cd8, cluster_columns= FALSE, cluster_rows = FALSE, col=col_fun, border= TRUE, show_row_names = FALSE, show_column_names = FALSE,
        column_labels = c("GzmB", "IFNg", "MIP1a", "MIP1b", "TNFa", "TNFb", "sCD137"), column_names_gp = gpar(fontsize=10),
        row_split = factor(rep(labs.pt4cd8), levels = c("3399_Pre", "3399_Post")),row_title = c("Pt4 Pre", "Pt4 Post"), row_title_gp = gpar(fontsize=12),
        show_heatmap_legend = FALSE)
dev.off()



##CD4 heatmaps
#Pt1 CD4
hm.pt1.cd4<-t(data.matrix(log.pt1.cd4))
labs.pt1cd4<-t(data.matrix(pt1.cd4labels))

Heatmap(hm.pt1.cd4, cluster_columns= FALSE, cluster_rows = TRUE, col=col_fun, border= TRUE, show_column_names = FALSE,
        row_labels = c("GzmB", "IFNg", "IL-10", "MIP1a", "MIP1b", "TNFa", "TNFb", "sCD137"), row_names_gp = gpar(fontsize=10),
        column_split = factor(labs.pt1cd4),column_title = c("Pt1 Post", "Pt1 Pre"), column_title_gp = gpar(fontsize=12),
        heatmap_legend_param = list(title="Log Intensity", title_position= "lefttop-rot"))


#Pt2 CD4
hm.pt2.cd4<-t(data.matrix(log.pt2.cd4))
labs.pt2cd4<-t(data.matrix(pt2.cd4labels))

Heatmap(hm.pt2.cd4, cluster_columns= FALSE, cluster_rows = TRUE, col=col_fun, border= TRUE, show_column_names = FALSE,
        row_labels = c("GzmB", "IFNg", "IL-10", "MIP1a", "MIP1b", "TNFa", "TNFb", "sCD137"), row_names_gp = gpar(fontsize=10),
        column_split = factor(labs.pt2cd4),column_title = c("Pt2 Post", "Pt2 Pre"), column_title_gp = gpar(fontsize=12),
        heatmap_legend_param = list(title="Log Intensity", title_position= "lefttop-rot"))

#Pt3 CD4
hm.pt3.cd4<-t(data.matrix(log.pt3.cd4))
labs.pt3cd4<-t(data.matrix(pt3.cd4labels))

Heatmap(hm.pt3.cd4, cluster_columns= FALSE, cluster_rows = TRUE, col=col_fun, border= TRUE, show_column_names = FALSE,
        row_labels = c("GzmB", "IFNg", "IL-10", "MIP1a", "MIP1b", "TNFa", "TNFb", "sCD137"), row_names_gp = gpar(fontsize=10),
        column_split = factor(labs.pt3cd4),column_title = c("Pt3 Post", "Pt3 Pre"), column_title_gp = gpar(fontsize=12),
        heatmap_legend_param = list(title="Log Intensity", title_position= "lefttop-rot"))

#Pt4 CD4
hm.pt4.cd4<-t(data.matrix(log.pt4.cd4))
labs.pt4cd4<-t(data.matrix(pt4.cd4labels))

Heatmap(hm.pt4.cd4, cluster_columns= FALSE, cluster_rows = TRUE, col=col_fun, border= TRUE, show_column_names = FALSE,
        row_labels = c("GzmB", "IFNg", "IL-10", "MIP1a", "MIP1b", "TNFa", "TNFb", "sCD137"), row_names_gp = gpar(fontsize=10),
        column_split = factor(labs.pt4cd4),column_title = c("Pt4 Post", "Pt4 Pre"), column_title_gp = gpar(fontsize=12),
        heatmap_legend_param = list(title="Log Intensity", title_position= "lefttop-rot"))