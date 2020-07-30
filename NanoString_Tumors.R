##############################################################################
## Changes in peripheral & local tumor immunity following NAC in BC         ##
## Data analysis by Justin M. Balko, Pharm.D., Ph.D. & Margaret L. Axelrod  ##
##                                                                          ##
##############################################################################

##########################
##  NanoString Tumors   ##
##########################

#PCA, gene and gene sets associated with outcome, heatmap 
#Figures 1B, 1C, 2, 3; S5, S6


#####################
##   Libraries     ##
#####################
library(colorspace)
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
library(survival)
library(survminer)
library(utils)


setwd("path_to_your_files")
read.delim("dataset.txt", row.names=1)->dataset #Table S2
read.delim("categories.txt", row.names=1)->categories #Table S3
read.delim("metadata.txt", row.names=1)->metadata #Table S1
all(rownames(metadata)==colnames(dataset)) #should be true

###Log transform
log(dataset,2)->log.data
log.data[,metadata$TIME=="post"]-log.data[,metadata$TIME=="pre"]->log.delta
log.data[,metadata$TIME=="post"]->log.post
log.data[,metadata$TIME=="pre"]->log.pre

#Subset data
metadata[metadata$TIME=="post",]->post.metadata
post.metadata[post.metadata$TNBC==TRUE,]->post.metadata.tnbc
post.metadata.tnbc[post.metadata.tnbc$TAXANE==TRUE,]->post.metadata.tnbc.t
post.metadata.tnbc[post.metadata.tnbc$TAXANE==FALSE,]->post.metadata.tnbc.not
post.metadata[post.metadata$TNBC==FALSE,]->post.metadata.nontnbc
post.metadata[post.metadata$ER_BIN==TRUE,]->post.metadata.er
post.metadata[post.metadata$HER2_BIN==TRUE,]->post.metadata.her2
log.delta[,post.metadata$TNBC==TRUE]->log.delta.tnbc
log.delta.tnbc[,post.metadata.tnbc$TAXANE==TRUE]->log.delta.tnbc.t
log.delta.tnbc[,post.metadata.tnbc$TAXANE==FALSE]->log.delta.tnbc.not
log.delta[,post.metadata$TNBC==FALSE]->log.delta.nontnbc
log.delta[,post.metadata$ER_BIN==TRUE]->log.delta.er
log.delta[,post.metadata$HER2_BIN==TRUE]->log.delta.her2
log.post[,post.metadata$TNBC==TRUE]->log.post.tnbc
log.pre[,post.metadata$TNBC==TRUE]->log.pre.tnbc
metadata[metadata$TIME=="pre",]->pre.metadata
pre.metadata[pre.metadata$TNBC==TRUE,]->pre.metadata.tnbc
pre.metadata.tnbc$PAM50->post.metadata.tnbc$PAM50

#PCA
pcadata<- cbind(t(log.data), metadata)
pca<-prcomp(t(log.data),scale=TRUE)
autoplot(pca, data= pcadata, colour= 'TIME') +theme_classic()

all(rownames(post.metadata)==colnames(log.delta))
pcadata.delta<- cbind(t(log.delta), post.metadata)
pca.delta<-prcomp(t(log.delta),scale=TRUE)
autoplot(pca.delta, data= pcadata.delta, colour= 'COHORT') +theme_classic()


#####################################################
#### Genes and gene sets associated with outcome ####
#####################################################


result<-NULL
for(i in 1:nrow(log.delta.tnbc)){
  log.delta.tnbc[i,]->gene
  coxph(Surv(post.metadata.tnbc$RFS_month,post.metadata.tnbc$RFS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="BH"))->result
rownames(log.delta.tnbc)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->delta.rfs

delta.rfs<-delta.rfs[order(delta.rfs$qval),]

delta.rfs2 <- mutate(delta.rfs, Significance=ifelse(delta.rfs$pval<0.05, "p-value<0.05", "Not Significant"))
rownames(delta.rfs)->rownames(delta.rfs2)
mutatedelta.rfs<- mutate(delta.rfs2, Significance=ifelse(delta.rfs2$qval<0.10, "q-value<0.10", Significance))

mutatedelta.rfs<- cbind(pathway=rownames(delta.rfs2), mutatedelta.rfs)
volc = ggplot(mutatedelta.rfs, aes(Coefficient, -log10(pval))) +
  geom_point(aes(col=Significance)) + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune genes vs. RFS")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatedelta.rfs[mutatedelta.rfs$qval<0.10,], 100),
                       aes(label=pathway))+
  scale_x_continuous(limits = c(-1.5,1.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))

#ggsave("RFS_genes_volcano.jpeg", device="jpeg")

#OS
result<-NULL
for(i in 1:nrow(log.delta.tnbc)){
  log.delta.tnbc[i,]->gene
  coxph(Surv(post.metadata.tnbc$OS_month,post.metadata.tnbc$OS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="BH"))->result
rownames(log.delta.tnbc)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->delta.os

delta.os<-delta.os[order(abs(delta.os$Coefficient),decreasing=TRUE),]

delta.os2 <- mutate(delta.os, Significance=ifelse(delta.os$pval<0.05, "p-value<0.05", "Not Significant"))
rownames(delta.os)->rownames(delta.os2)
mutatedelta.os<- mutate(delta.os2, Significance=ifelse(delta.os2$qval<0.10, "q-value<0.10", Significance))

mutatedelta.os<- cbind(pathway=rownames(delta.os2), mutatedelta.os)


volc = ggplot(mutatedelta.os, aes(Coefficient, -log10(pval))) +
  geom_point(aes(col=Significance)) + #add points colored by significance
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune genes vs. OS")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatedelta.os[mutatedelta.os$qval<0.10,], 50),
                       aes(label=pathway))+
  scale_x_continuous(limits = c(-1.5,1.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))

colnames(mutatedelta.rfs)<-paste(colnames(mutatedelta.rfs), ".rfs")
colnames(mutatedelta.os)<-paste(colnames(mutatedelta.os), ".os")
foo<-mutatedelta.os[match(mutatedelta.rfs$`pathway .rfs .rfs`, mutatedelta.os$`pathway .os`),]
genes.outcomes<-cbind(mutatedelta.rfs, foo)
#write.table(genes.outcomes, "C:/Users/marga/Desktop/TNBCgenesoutcomes.txt", quote=FALSE, sep="\t")

#Check genes for trends on KM plots

quantcut(as.numeric(log.delta.tnbc["CDH1",]),seq(0,1,by=0.33), na.rm=TRUE)->strat

data.frame("surv.time"=post.metadata.tnbc$RFS_month, "surv.event"=post.metadata.tnbc$RFS_event, "strat"=strat)->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ strat,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE, conf.int=TRUE, pval.coord=c(60,0.80), legend.labs=c("low", "intermediate",  "high"),
           xlab="Time after surgery (months)", ylab="Recurrence-free survival")

data.frame("surv.time"=post.metadata.tnbc$OS_month, "surv.event"=post.metadata.tnbc$OS_event, "strat"=strat)->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ strat,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE, pval.coord=c(75,0.50),legend.labs=c("low", "intermediate",  "high"),
           xlab="Time after surgery (months)", ylab="Overall survival")


##################Repeat for non-TNBC

#RFS
result<-NULL
for(i in 1:nrow(log.delta.nontnbc)){
  log.delta.nontnbc[i,]->gene
  coxph(Surv(post.metadata.nontnbc$RFS_month,post.metadata.nontnbc$RFS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="BH"))->result
rownames(log.delta.nontnbc)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->delta.rfs

delta.rfs<-delta.rfs[order(delta.rfs$qval),]

delta.rfs2 <- mutate(delta.rfs, Significance=ifelse(delta.rfs$pval<0.05, "p-value<0.05", "Not Significant"))
rownames(delta.rfs)->rownames(delta.rfs2)
mutatedelta.rfs<- mutate(delta.rfs2, Significance=ifelse(delta.rfs2$qval<0.10, "q-value<0.10", Significance))

mutatedelta.rfs<- cbind(pathway=rownames(delta.rfs2), mutatedelta.rfs)
volc = ggplot(mutatedelta.rfs, aes(Coefficient, -log10(pval))) +
  geom_point(aes(col=Significance)) + #add points colored by significance
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune genes vs. RFS")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatedelta.rfs[mutatedelta.rfs$qval<0.10,], 100),
                       aes(label=pathway))+
  scale_x_continuous(limits = c(-1.5,1.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))
#ggsave("RFS_genes_volcano_nonTNBC.jpeg", device="jpeg") #In case you want to easily save to disk


#OS
result<-NULL
for(i in 1:nrow(log.delta.nontnbc)){
  log.delta.nontnbc[i,]->gene
  coxph(Surv(post.metadata.nontnbc$OS_month,post.metadata.nontnbc$OS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="BH"))->result
rownames(log.delta.nontnbc)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->delta.os

delta.os<-delta.os[order(abs(delta.os$Coefficient),decreasing=TRUE),]

delta.os2 <- mutate(delta.os, Significance=ifelse(delta.os$pval<0.05, "p-value<0.05", "Not Significant"))
rownames(delta.os)->rownames(delta.os2)
mutatedelta.os<- mutate(delta.os2, Significance=ifelse(delta.os2$qval<0.10, "q-value<0.10", Significance))

mutatedelta.os<- cbind(pathway=rownames(delta.os2), mutatedelta.os)


volc = ggplot(mutatedelta.os, aes(Coefficient, -log10(pval))) +
  geom_point(aes(col=Significance)) + #add points colored by significance
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune genes vs. OS")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatedelta.os[mutatedelta.os$qval<0.10,], 50),
                       aes(label=pathway))+
  scale_x_continuous(limits = c(-1.5,1.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))



###Genes for ER+
#RFS
result<-NULL
for(i in 1:nrow(log.delta.er)){
  log.delta.er[i,]->gene
  coxph(Surv(post.metadata.er$RFS_month,post.metadata.er$RFS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="BH"))->result
rownames(log.delta.er)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->delta.rfs

delta.rfs<-delta.rfs[order(delta.rfs$qval),]

delta.rfs2 <- mutate(delta.rfs, Significance=ifelse(delta.rfs$pval<0.05, "p-value<0.05", "Not Significant"))
rownames(delta.rfs)->rownames(delta.rfs2)
mutatedelta.rfs<- mutate(delta.rfs2, Significance=ifelse(delta.rfs2$qval<0.10, "q-value<0.10", Significance))

mutatedelta.rfs<- cbind(pathway=rownames(delta.rfs2), mutatedelta.rfs)
volc = ggplot(mutatedelta.rfs, aes(Coefficient, -log10(pval))) +
  geom_point(aes(col=Significance)) + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune genes vs. RFS in ER+")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatedelta.rfs[mutatedelta.rfs$qval<0.10,], 100),
                       aes(label=pathway))+
  scale_x_continuous(limits = c(-1.5,1.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20)) +theme(plot.title= element_text(size=15), axis.title = element_text(size=15), 
                                             legend.title = element_text(size=14), legend.text = element_text(size=14))


#OS ER+
result<-NULL
for(i in 1:nrow(log.delta.er)){
  log.delta.er[i,]->gene
  coxph(Surv(post.metadata.er$OS_month,post.metadata.er$OS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="BH"))->result
rownames(log.delta.er)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->delta.os

delta.os<-delta.os[order(abs(delta.os$Coefficient),decreasing=TRUE),]

delta.os2 <- mutate(delta.os, Significance=ifelse(delta.os$pval<0.05, "p-value<0.05", "Not Significant"))
rownames(delta.os)->rownames(delta.os2)
mutatedelta.os<- mutate(delta.os2, Significance=ifelse(delta.os2$qval<0.10, "q-value<0.10", Significance))

mutatedelta.os<- cbind(pathway=rownames(delta.os2), mutatedelta.os)

volc = ggplot(mutatedelta.os, aes(Coefficient, -log10(pval))) +
  geom_point(aes(col=Significance)) + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune genes vs. OS in ER+")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatedelta.os[mutatedelta.os$qval<0.10,], 50),
                       aes(label=pathway))+
  scale_x_continuous(limits = c(-1.5,1.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))+theme(plot.title= element_text(size=15), axis.title = element_text(size=15), 
                                            legend.title = element_text(size=14), legend.text = element_text(size=14))

###Genes for HER2+
#RFS
result<-NULL
for(i in 1:nrow(log.delta.her2)){
  log.delta.her2[i,]->gene
  coxph(Surv(post.metadata.her2$RFS_month,post.metadata.her2$RFS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="BH"))->result
rownames(log.delta.her2)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->delta.rfs

delta.rfs<-delta.rfs[order(delta.rfs$qval),]

delta.rfs2 <- mutate(delta.rfs, Significance=ifelse(delta.rfs$pval<0.05, "p-value<0.05", "Not Significant"))
rownames(delta.rfs)->rownames(delta.rfs2)
mutatedelta.rfs<- mutate(delta.rfs2, Significance=ifelse(delta.rfs2$qval<0.10, "q-value<0.10", Significance))

mutatedelta.rfs<- cbind(pathway=rownames(delta.rfs2), mutatedelta.rfs)
volc = ggplot(mutatedelta.rfs, aes(Coefficient, -log10(pval))) +
  geom_point(aes(col=Significance)) + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune genes vs. RFS in HER2+")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatedelta.rfs[mutatedelta.rfs$qval<0.10,], 100),
                       aes(label=pathway))+
  scale_x_continuous(limits = c(-1.5,1.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20)) +theme(plot.title= element_text(size=15), axis.title = element_text(size=15), 
                                             legend.title = element_text(size=14), legend.text = element_text(size=14))


#OS HER2+
result<-NULL
for(i in 1:nrow(log.delta.her2)){
  log.delta.her2[i,]->gene
  coxph(Surv(post.metadata.her2$OS_month,post.metadata.her2$OS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="BH"))->result
rownames(log.delta.her2)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->delta.os

delta.os<-delta.os[order(abs(delta.os$Coefficient),decreasing=TRUE),]

delta.os2 <- mutate(delta.os, Significance=ifelse(delta.os$pval<0.05, "p-value<0.05", "Not Significant"))
rownames(delta.os)->rownames(delta.os2)
mutatedelta.os<- mutate(delta.os2, Significance=ifelse(delta.os2$qval<0.10, "q-value<0.10", Significance))

mutatedelta.os<- cbind(pathway=rownames(delta.os2), mutatedelta.os)

volc = ggplot(mutatedelta.os, aes(Coefficient, -log10(pval))) +
  geom_point(aes(col=Significance)) + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune genes vs. OS in HER2+")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatedelta.os[mutatedelta.os$qval<0.10,], 50),
                       aes(label=pathway))+
  scale_x_continuous(limits = c(-1.5,1.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))+theme(plot.title= element_text(size=15), axis.title = element_text(size=15), 
                                            legend.title = element_text(size=14), legend.text = element_text(size=14))

##############################Signature loop########################################################

#for all patients
result<-NULL
for(i in 1:ncol(categories)){
  rownames(categories)[categories[,i]==1]->sig
  print(length(which(rownames(log.data)%in%sig))==length(sig))
  apply(log.data[rownames(log.data)%in%sig,],2,mean,na.rm=TRUE)->sig.score
  cbind(result,sig.score)->result
}
colnames(categories)->colnames(result)
result->sig.data
t(sig.data)->sig.data
sig.data[,metadata$TIME=="post"]-sig.data[,metadata$TIME=="pre"]->sig.delta
sig.delta[,post.metadata$TNBC==TRUE]->sig.delta.tnbc
sig.delta.tnbc[,post.metadata.tnbc$TAXANE==TRUE]->sig.delta.tnbc.t
sig.delta.tnbc[,post.metadata.tnbc$TAXANE==FALSE]->sig.delta.tnbc.not

sig.delta[,post.metadata$TNBC==FALSE]->sig.delta.nontnbc
sig.delta[,post.metadata$ER_BIN==TRUE]->sig.delta.er
sig.delta[,post.metadata$HER2_BIN==TRUE]->sig.delta.her2

#Genesets TNBC RFS
result<-NULL
for(i in 1:nrow(sig.delta.tnbc)){
  sig.delta.tnbc[i,]->gene
  coxph(Surv(post.metadata.tnbc$RFS_month,post.metadata.tnbc$RFS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="fdr"))->result
rownames(sig.delta.tnbc)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->sig.rfs
sig.rfs.perm<-sig.rfs
sig.rfs<-sig.rfs[order(sig.rfs$pval),]
sig.rfs2 <- mutate(sig.rfs, Significance=ifelse(sig.rfs$pval<0.05, "p-value<0.05", "Not Significant"))
mutatesig.rfs <- mutate(sig.rfs2, Significance=ifelse(sig.rfs2$qval<0.10, "q-value<0.10", Significance))
mutatesig.rfs<- cbind(pathway=rownames(sig.rfs), mutatesig.rfs)
volc = ggplot(mutatesig.rfs, aes(Coefficient, -log10(pval))) +
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  geom_point(aes(col=Significance)) + #add points colored by significance
  scale_color_manual(values=c("black", "green")) +
  ggtitle("Changes in immune signatures vs. RFS")+
  ylab("-Log10(p-value)")
volc+geom_text_repel(data=head(mutatesig.rfs[mutatesig.rfs$qval<0.10,], 8), aes(label=pathway),
                     point.padding=.3)+
  scale_x_continuous(limits = c(-2,0.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))
#ggsave("RFS_geneset_volcano.jpeg", device="jpeg") #In case you want to easily save to disk


#OS - TNBC signatures
result<-NULL
for(i in 1:nrow(sig.delta.tnbc)){
  sig.delta.tnbc[i,]->gene
  coxph(Surv(post.metadata.tnbc$OS_month,post.metadata.tnbc$OS_event)~as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="fdr"))->result
rownames(sig.data)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->sig.os
sig.os.perm<-sig.os
sig.os<-sig.os[order(sig.os$pval),]
sig.os2 <- mutate(sig.os, Significance=ifelse(sig.os$pval<0.05, "p-value<0.05", "Not Significant"))
mutatesig.os <- mutate(sig.os2, Significance=ifelse(sig.os2$qval<0.10, "q-value<0.10", Significance))
mutatesig.os<- cbind(pathway=rownames(sig.os), mutatesig.os)
volc = ggplot(mutatesig.os, aes(Coefficient, -log10(pval))) +
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  geom_point(aes(col=Significance)) + #add points colored by significance
  scale_color_manual(values=c("black", "green")) +
  ggtitle("Changes in immune signatures vs. OS")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatesig.os[mutatesig.os$qval<0.10,], 8), aes(label=pathway),
                       point.padding = 0.35)+
  scale_x_continuous(limits = c(-2,0.5), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))
#ggsave("OS_geneset_volcano.jpeg", device="jpeg") #In case you want to easily save to disk

colnames(mutatesig.rfs)<-paste(colnames(mutatesig.rfs), ".rfs")
colnames(mutatesig.os)<-paste(colnames(mutatesig.os), ".os")
foo<-mutatesig.os[match(mutatesig.rfs$`pathway .rfs`, mutatesig.os$`pathway .os`),]
sig.outcomes<-cbind(mutatesig.rfs, foo)
#write.table(sig.outcomes, "C:/Users/marga/Desktop/TNBCsigoutcomes.txt", quote=FALSE, sep="\t")


#Sig check
rainbow_hcl(5,c=260,l=70, fixup=TRUE)->cols
#quantcut(as.numeric(sig.delta.tnbc["T.cell_activation",]),seq(0,1,by=0.33), na.rm=TRUE)->strat
#quantcut(as.numeric(sig.delta.tnbc["Phagyocytosis_signal_transduction",]),seq(0,1,by=0.33), na.rm=TRUE)->strat
quantcut(as.numeric(sig.delta.tnbc["NK_Cell_Functions",]),seq(0,1,by=0.33), na.rm=TRUE)->strat


data.frame("surv.time"=post.metadata.tnbc$RFS_month, "surv.event"=post.metadata.tnbc$RFS_event, "strat"=strat)->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ strat,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE, pval.coord=c(60,0.80), legend.labs=c("down", "equivocal", "up"),
           xlab="Time after surgery (months)", ylab="Recurrence-free survival")

data.frame("surv.time"=post.metadata.tnbc$OS_month, "surv.event"=post.metadata.tnbc$OS_event, "strat"=strat)->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ strat,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE, pval.coord=c(75,0.50),legend.labs=c("down", "equivocal", "up"),
           xlab="Time after surgery (months)", ylab="Overall survival")




##Non TNBC
sig.delta[,post.metadata$TNBC==TRUE]->sig.delta.tnbc
sig.delta[,post.metadata$TNBC==FALSE]->sig.delta.nontnbc


#RFS nonTNBC gene sets
result<-NULL
for(i in 1:nrow(sig.delta.nontnbc)){
  sig.delta.nontnbc[i,]->gene
  coxph(Surv(post.metadata.nontnbc$RFS_month,post.metadata.nontnbc$RFS_event) ~ as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="fdr"))->result
rownames(sig.delta.nontnbc)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->sig.rfs
sig.rfs<-sig.rfs[order(sig.rfs$pval),]
sig.rfs2 <- mutate(sig.rfs, Significance=ifelse(sig.rfs$pval<0.05, "p-value<0.05", "Not Significant"))
mutatesig.rfs <- mutate(sig.rfs2, Significance=ifelse(sig.rfs2$qval<0.10, "q-value<0.10", Significance))
mutatesig.rfs<- cbind(pathway=rownames(sig.rfs), mutatesig.rfs)
volc = ggplot(mutatesig.rfs, aes(Coefficient, -log10(pval))) +
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  geom_point(aes(col=Significance)) + #add points colored by significance
  scale_color_manual(values=c("black", "red","green")) +
  ggtitle("Changes in immune signatures vs. RFS")+
  ylab("-Log10(p-value)")
volc+geom_text_repel(data=head(mutatesig.rfs[mutatesig.rfs$qval<0.10,], 10), aes(label=pathway))+
  scale_x_continuous(limits = c(-2,2), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))

#ggsave("RFS_geneset_volcano.jpeg", device="jpeg") #In case you want to easily save to disk

#OS nonTNBC genesets
result<-NULL
for(i in 1:nrow(sig.delta.nontnbc)){
  sig.delta.nontnbc[i,]->gene
  coxph(Surv(post.metadata.nontnbc$OS_month,post.metadata.nontnbc$OS_event)~as.numeric(gene))->tmp
  summary(tmp)$coeff[1,5]->p.val
  summary(tmp)$coeff[1,1]->coeff
  rbind(result,c(coeff,p.val))->result
}
cbind(result,p.adjust(result[,2],method="fdr"))->result
rownames(sig.data)->rownames(result)
colnames(result)<-c("Coefficient","pval","qval")
as.data.frame(result)->sig.os
sig.os<-sig.os[order(sig.os$pval),]
sig.os2 <- mutate(sig.os, Significance=ifelse(sig.os$pval<0.05, "p-value<0.05", "Not Significant"))
mutatesig.os <- mutate(sig.os2, Significance=ifelse(sig.os2$qval<0.10, "q-value<0.10", Significance))
mutatesig.os<- cbind(pathway=rownames(sig.os), mutatesig.os)
volc = ggplot(mutatesig.os, aes(Coefficient, -log10(pval))) +
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=1.5)+
  geom_point(aes(col=Significance)) + #add points colored by significance
  scale_color_manual(values=c("black", "red")) +
  ggtitle("Changes in immune signatures vs. OS")+
  ylab("-Log10(p-value)")
volc + geom_text_repel(data=head(mutatesig.os[mutatesig.os$pval<0.05,], 5), aes(label=pathway))+
  scale_x_continuous(limits = c(-2,2), expand = c(0,0) )+
  theme_set(theme_bw(base_size = 20))

#####Heatmap of nanoString TUMORS###########
log.data->heatmap.data
metadata->tmp.meta

colclus <- hclust(as.dist(1-cor(heatmap.data,method="spearman")), method="ward.D")
rowclus <- hclust(as.dist(1-cor(t(heatmap.data),method="spearman")), method="ward.D")
diverge_hcl(100,h=c(260,0),c=100,l = c(55, 100),power=1)->cols

#Create color side bars
levels(tmp.meta$TIME)<-c("darkred","pink")
time_colors=as.character(tmp.meta$TIME)
levels(tmp.meta$COHORT)<-c("red","green","darkorchid")
cohort_colors=as.character(tmp.meta$COHORT)
levels(tmp.meta$PAM50)<-c("red","pink","cyan","blue","green")
subtype_colors=as.character(tmp.meta$PAM50)
as.factor(tmp.meta$HER2_BIN)->tmp
levels(tmp)<-c("black","grey","white")
HER2colors=as.character(tmp)
addNA(as.factor(tmp.meta$PR_BIN))->tmp
levels(tmp)<-c("black","grey","white")
PRcolors=as.character(tmp)
addNA(as.factor(tmp.meta$ER_BIN))->tmp
levels(tmp)<-c("black","grey","white")
ERcolors=as.character(tmp)
addNA(tmp.meta$TP53_type)->tmp.meta$TP53_type
levels(tmp.meta$TP53_type)<-c("red","yellow","purple","green","grey")
TP53colors=as.character(tmp.meta$TP53_type)

clab=cbind(time_colors, cohort_colors, subtype_colors,HER2colors,PRcolors,ERcolors,TP53colors)
colnames(clab)=c("Time","Cohort","Subtype","HER2","PR","ER","TP53")

heatmap3(heatmap.data, na.rm = TRUE, scale="row", margins=c(6,10),
         Rowv = NULL, Colv= as.dendrogram(colclus),ColSideColors=clab,
         labCol=FALSE, labRow=FALSE, cexRow=1, col=cols)
legend("topright",legend=c("Time","Pre","Post","","Cohort","DART","VICC","PERU","","Subtype","Basal","Her2","LumA","LumB","Normal","","IHC","Positive","Negative","","TP53","Del/Trunc", "Other","Point","WT","NA"),
       fill=c("white","pink","darkred","white","white","red","green","darkorchid","white","white","red","pink","cyan","blue","green","white","white","black","grey","white","white","red","yellow","purple","green","grey"),
       border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

#DELTA Dataset
log.delta->heatmap.data
metadata[metadata$TIME=="post",]->tmp.meta

#Create color side bars
levels(tmp.meta$COHORT)<-c("red","green","darkorchid")
cohort_colors=as.character(tmp.meta$COHORT)
levels(tmp.meta$PAM50)<-c("red","pink","cyan","blue","green")
subtype_colors=as.character(tmp.meta$PAM50)
as.factor(tmp.meta$HER2_BIN)->tmp
levels(tmp)<-c("black","grey","white")
HER2colors=as.character(tmp)
addNA(as.factor(tmp.meta$PR_BIN))->tmp
levels(tmp)<-c("black","grey","white")
PRcolors=as.character(tmp)
addNA(as.factor(tmp.meta$ER_BIN))->tmp
levels(tmp)<-c("black","grey","white")
ERcolors=as.character(tmp)
addNA(tmp.meta$TP53_type)->tmp.meta$TP53_type
levels(tmp.meta$TP53_type)<-c("red","yellow","purple","green","grey")
TP53colors=as.character(tmp.meta$TP53_type)

clab=cbind(cohort_colors, subtype_colors,HER2colors,PRcolors,ERcolors,TP53colors)
colnames(clab)=c("Cohort","Subtype","HER2","PR","ER","TP53")

colclus <- hclust(as.dist(1-cor(heatmap.data,method="spearman")), method="ward.D")
rowclus <- hclust(as.dist(1-cor(t(heatmap.data),method="spearman")), method="ward.D")
diverge_hcl(100,h=c(260,0),c=100,l = c(55, 100),power=0.7)->cols

heatmap3(heatmap.data, na.rm = TRUE, scale="row", margins=c(6,20),
         Rowv = NULL, Colv= as.dendrogram(colclus),ColSideColors=clab,
         labCol=FALSE, labRow=FALSE, cexRow=1, col=cols)
legend("topright",legend=c("Time","Pre","Post","","Cohort","DART","VICC","PERU","","Subtype","Basal","Her2","LumA","LumB","Normal","","IHC","Positive","Negative","","TP53","Del/Trunc", "Other","Point","WT","NA"),
       fill=c("white","pink","darkred","white","white","red","green","darkorchid","white","white","red","pink","cyan","blue","green","white","white","black","grey","white","white","red","yellow","purple","green","grey"),
       border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)



#DELTA Dataset SIGS
sig.delta.tnbc[sig.rfs.perm$pval<0.05,]->heatmap.data
apply(heatmap.data,2,sum)->tmp
heatmap.data[,order(tmp)]->heatmap.data
post.metadata.tnbc[order(tmp),]->tmp.meta
#heatmap.data[,order(post.metadata.tnbc$TP53_type)]->heatmap.data #for reodering by p53
#post.metadata.tnbc[order(addNA(post.metadata.tnbc$TP53_type)),]->tmp.meta
#Create color side bars

addNA(tmp.meta$PAM50)->tmp.meta$PAM50
levels(tmp.meta$PAM50)<-c("red","pink","cyan","blue","green","grey")
subtype_colors=as.character(tmp.meta$PAM50)
addNA(tmp.meta$TP53_type)->tmp.meta$TP53_type
levels(tmp.meta$TP53_type)<-c("red","yellow","purple","green","grey")
TP53colors=as.character(tmp.meta$TP53_type)

clab=cbind(subtype_colors,TP53colors)
colnames(clab)=c("Subtype","TP53")

colclus <- hclust(as.dist(1-cor(heatmap.data,method="spearman")), method="ward.D")
rowclus <- hclust(as.dist(1-cor(t(heatmap.data),method="spearman")), method="ward.D")
diverge_hcl(100,h=c(260,0),c=100,l = c(55, 100),power=0.7)->cols

heatmap3(heatmap.data, na.rm = TRUE, scale="row", margins=c(10,20),
         Rowv = as.dendrogram(rowclus), Colv= NA,ColSideColors=clab,labRow=rownames(heatmap.data),
         cexRow=0.6, col=cols)
legend("topright",legend=c("Subtype","Basal","Her2","LumA","LumB","Normal","NA","","TP53","Del/Trunc", "Other","Point","WT","NA"),
       fill=c("white","red","pink","cyan","blue","green","grey","white","white","red","yellow","purple","green","grey"),
       border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

##rearrange by p53 status
heatmap.data[,order(tmp.meta$TP53_type)]->heatmap.data
tmp.meta[order(addNA(tmp.meta$TP53_type)),]->tmp.meta
addNA(tmp.meta$PAM50)->tmp.meta$PAM50
levels(tmp.meta$PAM50)<-c("red","pink","cyan","blue","green","grey")
subtype_colors=as.character(tmp.meta$PAM50)
addNA(tmp.meta$TP53_type)->tmp.meta$TP53_type
levels(tmp.meta$TP53_type)<-c("red","yellow","purple","green","grey")
TP53colors=as.character(tmp.meta$TP53_type)

clab=cbind(subtype_colors,TP53colors)
colnames(clab)=c("Subtype","TP53")

colclus <- hclust(as.dist(1-cor(heatmap.data,method="spearman")), method="ward.D")
rowclus <- hclust(as.dist(1-cor(t(heatmap.data),method="spearman")), method="ward.D")
diverge_hcl(100,h=c(260,0),c=100,l = c(55, 100),power=2)->cols

heatmap3(heatmap.data, na.rm = TRUE, scale="row", margins=c(10,20),
         Rowv = as.dendrogram(rowclus), Colv= NA,ColSideColors=clab,labRow=rownames(heatmap.data),
         cexRow=0.6, col=cols)
legend("topright",legend=c("Subtype","Basal","Her2","LumA","LumB","Normal","NA","","TP53","Del/Trunc", "Other","Point","WT","NA"),
       fill=c("white","red","pink","cyan","blue","green","grey","white","white","red","yellow","purple","green","grey"),
       border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
