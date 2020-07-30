##############################################################################
## Changes in peripheral & local tumor immunity following NAC in BC         ##
## Data analysis by Justin M. Balko, Pharm.D., Ph.D. & Margaret L. Axelrod  ##
##                                                                          ##
##############################################################################

############################
#### Basic TILS and KM #####
############################

#Figures 1A, S2, S3, S4
#Metadata = Table S1


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

###Paired TILS plot
metadata$uniquept<-paste(metadata$COHORT,metadata$PATIENT)
metadata.post<-metadata[metadata$TIME=="post",]
metadata.pre<-metadata[metadata$TIME=="pre",]
#Can only combine like this because metadata perfectly in order with post, pre for each patient sequential
foo<-data.frame(PRE=metadata.pre$TILS, POST=metadata.post$TILS, PAM50=metadata.pre$PAM50, TNBC=metadata.pre$TNBC, ER=metadata.pre$ER_BIN, HER2=metadata.pre$HER2_BIN, ID=metadata.pre$uniquept)
ggpaired(foo, cond1= "PRE", cond2 = "POST", fill= "condition", alpha=0.2)+ylab("TILs") +xlab("Time") +theme(legend.position = "none")
ggpaired(foo, cond1= "PRE", cond2 = "POST", facet.by = "PAM50")+ylab("TILs") +xlab("Time")
ggpaired(foo, cond1= "PRE", cond2 = "POST", facet.by = "TNBC")+ylab("TILs") +xlab("Time") + ggtitle("TNBC")
ggpaired(foo, cond1= "PRE", cond2 = "POST", facet.by = "ER")+ylab("TILs") +xlab("Time") + ggtitle("ER+")
ggpaired(foo, cond1= "PRE", cond2 = "POST", facet.by = "HER2")+ylab("TILs") +xlab("Time") + ggtitle("HER2+")

foo$diff<-foo$POST-foo$PRE
mean<-round(mean(foo[["diff"]], na.rm=TRUE), digits=2)
ggplot(metadata, aes(TIME, TILS))+ geom_point()+ geom_line(aes(group=uniquept))+
  scale_x_discrete(limits=c("pre", "post"))+
  ggtitle("All Patients")+
  annotate("text", x=1.5, y=92, label=paste("Mean Difference=", mean))+
  theme_classic()
#do for TNBC, ER+, HER2+
meta.tnbc<-metadata[metadata$TNBC==TRUE,]
foo.tnbc<-foo[foo$TNBC==TRUE,]
mean<-round(mean(foo.tnbc[["diff"]], na.rm=TRUE), digits=2)
ggplot(meta.tnbc, aes(TIME, TILS))+ geom_point()+ geom_line(aes(group=uniquept))+
  scale_x_discrete(limits=c("pre", "post"))+
  ggtitle("TNBC Only")+
  annotate("text", x=1.5, y=90, label=paste("Mean Difference=", mean))+
  theme_classic()
#ER
meta.er<-metadata[metadata$ER_BIN==TRUE,]
foo.er<-foo[foo$ER==TRUE,]
mean<-round(mean(foo.er[["diff"]], na.rm=TRUE), digits=2)
ggplot(meta.er, aes(TIME, TILS))+ geom_point()+ geom_line(aes(group=uniquept))+
  scale_x_discrete(limits=c("pre", "post"))+
  ggtitle("ER+ Only")+
  annotate("text", x=1.5, y=90, label=paste("Mean Difference=", mean))+
  theme_classic()
#HER2
meta.her2<-metadata[metadata$HER2_BIN==TRUE,]
foo.her2<-foo[foo$HER2==TRUE,]
mean<-round(mean(foo.her2[["diff"]], na.rm=TRUE), digits=2)
ggplot(meta.her2, aes(TIME, TILS))+ geom_point()+ geom_line(aes(group=uniquept))+
  scale_x_discrete(limits=c("pre", "post"))+
  ggtitle("HER2+ Only")+
  annotate("text", x=1.5, y=90, label=paste("Mean Difference=", mean))+
  theme_classic()



#############################
##  Basic TILs KM analysis ##
#############################
attach(metadata)
coxph(Surv(RFS_month[TIME=="pre"],RFS_event[TIME=="pre"]) ~ TILS[TIME=="pre"])
coxph(Surv(OS_month[TIME=="pre"],OS_event[TIME=="pre"]) ~ TILS[TIME=="pre"])
coxph(Surv(RFS_month[TIME=="post"],RFS_event[TIME=="post"]) ~ TILS[TIME=="post"])
coxph(Surv(OS_month[TIME=="post"],OS_event[TIME=="post"]) ~ TILS[TIME=="post"])

##Only in TNBC
coxph(Surv(RFS_month[TNBC=="TRUE" & TIME=="pre"],RFS_event[TNBC=="TRUE" & TIME=="pre"]) ~ TILS[TNBC=="TRUE" & TIME=="pre"])
coxph(Surv(OS_month[TNBC=="TRUE" & TIME=="pre"],OS_event[TNBC=="TRUE" & TIME=="pre"]) ~ TILS[TNBC=="TRUE" & TIME=="pre"])
coxph(Surv(RFS_month[TNBC=="TRUE" & TIME=="post"],RFS_event[TNBC=="TRUE" & TIME=="post"]) ~ TILS[TNBC=="TRUE" & TIME=="post"])
coxph(Surv(OS_month[TNBC=="TRUE" & TIME=="post"],OS_event[TNBC=="TRUE" & TIME=="post"]) ~ TILS[TNBC=="TRUE" & TIME=="post"])

#TILs delta in TNBC
TILS[TNBC=="TRUE"&TIME=="post"]-TILS[TNBC=="TRUE"&TIME=="pre"]->tils.delta
TILS[TIME=="post"]-TILS[TIME=="pre"]->tils.delta.all
coxph(Surv(RFS_month[TNBC=="TRUE"&TIME=="pre"],RFS_event[TNBC=="TRUE"&TIME=="pre"]) ~ tils.delta)
coxph(Surv(OS_month[TNBC=="TRUE"&TIME=="pre"],OS_event[TNBC=="TRUE"&TIME=="pre"]) ~ tils.delta)


#### sTILs KM ####
as.factor(TILS>29)->TILs # factor TILs by a 29% cutoff
c("sTILs<30","sTILs>30")-> levels(TILs)

#ALL PATIENTS
#Pre, all patients, RFS
data.frame("surv.time"=RFS_month[TIME=="pre"], "surv.event"=RFS_event[TIME=="pre"], "TILs"=TILs[TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE, xlab="Time after surgery (months)", ylab="Recurrence-free survival" )
#Pre, all patients, OS
data.frame("surv.time"=OS_month[TIME=="pre"], "surv.event"=OS_event[TIME=="pre"], "TILs"=TILs[TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival" )
#Post, all patients, RFS
data.frame("surv.time"=RFS_month[TIME=="pre"], "surv.event"=RFS_event[TIME=="pre"], "TILs"=TILs[TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival" )
#Post, all patients, OS
data.frame("surv.time"=OS_month[TIME=="pre"], "surv.event"=OS_event[TIME=="pre"], "TILs"=TILs[TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival" )
####TNBC+/-
#Pre, TNBC, RFS
data.frame("surv.time"=RFS_month[TNBC=="TRUE"&TIME=="pre"], "surv.event"=RFS_event[TNBC=="TRUE"&TIME=="pre"], "TILs"=TILs[TNBC=="TRUE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival" )
#Pre, TNBC, OS
data.frame("surv.time"=OS_month[TNBC=="TRUE"&TIME=="pre"], "surv.event"=OS_event[TNBC=="TRUE"&TIME=="pre"], "TILs"=TILs[TNBC=="TRUE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival" )
#Post, TNBC, RFS
data.frame("surv.time"=RFS_month[TNBC=="TRUE"&TIME=="pre"], "surv.event"=RFS_event[TNBC=="TRUE"&TIME=="pre"], "TILs"=TILs[TNBC=="TRUE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", pval.coord=c(60,0.9))
#Post, TNBC, OS
data.frame("surv.time"=OS_month[TNBC=="TRUE"&TIME=="pre"], "surv.event"=OS_event[TNBC=="TRUE"&TIME=="pre"], "TILs"=TILs[TNBC=="TRUE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival" )
#Post, Non-TNBC, RFS
data.frame("surv.time"=RFS_month[TNBC=="FALSE"&TIME=="pre"], "surv.event"=RFS_event[TNBC=="FALSE"&TIME=="pre"], "TILs"=TILs[TNBC=="FALSE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival")
#Post, non-TNBC, OS
data.frame("surv.time"=OS_month[TNBC=="FALSE"&TIME=="pre"], "surv.event"=OS_event[TNBC=="FALSE"&TIME=="pre"], "TILs"=TILs[TNBC=="FALSE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival")
###TNBC +/- TAXANE
#Pre, TNBC, w/ taxane, RFS
data.frame("surv.time"=RFS_month[TNBC=="TRUE"&TIME=="pre"&TAXANE=="TRUE"], "surv.event"=RFS_event[TNBC=="TRUE"&TIME=="pre"&TAXANE=="TRUE"], "TILs"=TILs[TNBC=="TRUE"&TIME=="pre"&TAXANE=="TRUE"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Pre-NAC sTILs in TNBC w/ Taxane")
#Pre, TNBC, w/o taxane, RFS
data.frame("surv.time"=RFS_month[TNBC=="TRUE"&TIME=="pre"&TAXANE=="FALSE"], "surv.event"=RFS_event[TNBC=="TRUE"&TIME=="pre"&TAXANE=="FALSE"], "TILs"=TILs[TNBC=="TRUE"&TIME=="pre"&TAXANE=="FALSE"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Pre-NAC sTILs in TNBC w/o Taxane")
#Pre, TNBC, OS, w/ TAXANE
data.frame("surv.time"=OS_month[TNBC=="TRUE"&TIME=="pre"&TAXANE=="TRUE"], "surv.event"=OS_event[TNBC=="TRUE"&TIME=="pre"&TAXANE=="TRUE"], "TILs"=TILs[TNBC=="TRUE"&TIME=="pre"&TAXANE=="TRUE"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title= "Pre-NAC sTILS in TNBC w/ Taxane")
#Pre, TNBC, OS, w/o TAXANE
data.frame("surv.time"=OS_month[TNBC=="TRUE"&TIME=="pre"&TAXANE=="FALSE"], "surv.event"=OS_event[TNBC=="TRUE"&TIME=="pre"&TAXANE=="FALSE"], "TILs"=TILs[TNBC=="TRUE"&TIME=="pre"&TAXANE=="FALSE"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title= "Pre-NAC sTILS in TNBC w/o Taxane")
#Post, TNBC, w/ taxane, RFS
data.frame("surv.time"=RFS_month[TNBC=="TRUE"&TIME=="post"&TAXANE=="TRUE"], "surv.event"=RFS_event[TNBC=="TRUE"&TIME=="post"&TAXANE=="TRUE"], "TILs"=TILs[TNBC=="TRUE"&TIME=="post"&TAXANE=="TRUE"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Post-NAC sTILs in TNBC w/ Taxane")
#Post, TNBC, w/o taxane, RFS
data.frame("surv.time"=RFS_month[TNBC=="TRUE"&TIME=="post"&TAXANE=="FALSE"], "surv.event"=RFS_event[TNBC=="TRUE"&TIME=="post"&TAXANE=="FALSE"], "TILs"=TILs[TNBC=="TRUE"&TIME=="post"&TAXANE=="FALSE"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Post-NAC sTILs in TNBC w/o Taxane")
#Post, TNBC, OS, w/ TAXANE
data.frame("surv.time"=OS_month[TNBC=="TRUE"&TIME=="post"&TAXANE=="TRUE"], "surv.event"=OS_event[TNBC=="TRUE"&TIME=="post"&TAXANE=="TRUE"], "TILs"=TILs[TNBC=="TRUE"&TIME=="post"&TAXANE=="TRUE"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title= "Post-NAC sTILS in TNBC w/ Taxane")
#Post, TNBC, OS, w/o TAXANE
data.frame("surv.time"=OS_month[TNBC=="TRUE"&TIME=="post"&TAXANE=="FALSE"], "surv.event"=OS_event[TNBC=="TRUE"&TIME=="post"&TAXANE=="FALSE"], "TILs"=TILs[TNBC=="TRUE"&TIME=="post"&TAXANE=="FALSE"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title= "Post-NAC sTILS in TNBC w/o Taxane")

##Additional KMs
#####ER POSITIVE
#Pre, ER+, RFS
data.frame("surv.time"=RFS_month[ER_BIN=="TRUE"&TIME=="pre"], "surv.event"=RFS_event[ER_BIN=="TRUE"&TIME=="pre"], "TILs"=TILs[ER_BIN=="TRUE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Pre-NAC sTILs in ER+")
#Pre, ER+, OS
data.frame("surv.time"=OS_month[ER_BIN=="TRUE"&TIME=="pre"], "surv.event"=OS_event[ER_BIN=="TRUE"&TIME=="pre"], "TILs"=TILs[ER_BIN=="TRUE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title= "Pre-NAC sTILs in ER+")
#Post, ER+, RFS
data.frame("surv.time"=RFS_month[ER_BIN=="TRUE"&TIME=="pre"], "surv.event"=RFS_event[ER_BIN=="TRUE"&TIME=="pre"], "TILs"=TILs[ER_BIN=="TRUE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Post-NAC sTILs in ER+")
#Post, ER+, OS
data.frame("surv.time"=OS_month[ER_BIN=="TRUE"&TIME=="pre"], "surv.event"=OS_event[ER_BIN=="TRUE"&TIME=="pre"], "TILs"=TILs[ER_BIN=="TRUE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title="Post-NAC sTILs in ER+")
#####ER NEGATIVE
#Pre, ALL ER-, RFS
data.frame("surv.time"=RFS_month[ER_BIN=="FALSE"&TIME=="pre"], "surv.event"=RFS_event[ER_BIN=="FALSE"&TIME=="pre"], "TILs"=TILs[ER_BIN=="FALSE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Pre-NAC sTILs in ER-")
#Pre, ER-, OS
data.frame("surv.time"=OS_month[ER_BIN=="FALSE"&TIME=="pre"], "surv.event"=OS_event[ER_BIN=="FALSE"&TIME=="pre"], "TILs"=TILs[ER_BIN=="FALSE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title= "Pre-NAC sTILs in ER-")
#Post, ER-, RFS
data.frame("surv.time"=RFS_month[ER_BIN=="FALSE"&TIME=="pre"], "surv.event"=RFS_event[ER_BIN=="FALSE"&TIME=="pre"], "TILs"=TILs[ER_BIN=="FALSE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Post-NAC sTILs in ER-")
#Post, ER-, OS
data.frame("surv.time"=OS_month[ER_BIN=="FALSE"&TIME=="pre"], "surv.event"=OS_event[ER_BIN=="FALSE"&TIME=="pre"], "TILs"=TILs[ER_BIN=="FALSE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title="Post-NAC sTILs in ER-")
#####HER2 POSITIVE
#Pre, HER2+, RFS
data.frame("surv.time"=RFS_month[HER2_BIN=="TRUE"&TIME=="pre"], "surv.event"=RFS_event[HER2_BIN=="TRUE"&TIME=="pre"], "TILs"=TILs[HER2_BIN=="TRUE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Pre-NAC sTILs in HER2+")
#Pre, HER2+, OS
data.frame("surv.time"=OS_month[HER2_BIN=="TRUE"&TIME=="pre"], "surv.event"=OS_event[HER2_BIN=="TRUE"&TIME=="pre"], "TILs"=TILs[HER2_BIN=="TRUE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title= "Pre-NAC sTILs in HER2+")
#Post, HER2+, RFS
data.frame("surv.time"=RFS_month[HER2_BIN=="TRUE"&TIME=="pre"], "surv.event"=RFS_event[HER2_BIN=="TRUE"&TIME=="pre"], "TILs"=TILs[HER2_BIN=="TRUE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Post-NAC sTILs in HER2+")
#Post, HER2+, OS
data.frame("surv.time"=OS_month[HER2_BIN=="TRUE"&TIME=="pre"], "surv.event"=OS_event[HER2_BIN=="TRUE"&TIME=="pre"], "TILs"=TILs[HER2_BIN=="TRUE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title="Post-NAC sTILs in HER2+")
#####HER2 NEGATIVE
#Pre, ALL HER2-, RFS
data.frame("surv.time"=RFS_month[HER2_BIN=="FALSE"&TIME=="pre"], "surv.event"=RFS_event[HER2_BIN=="FALSE"&TIME=="pre"], "TILs"=TILs[HER2_BIN=="FALSE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Pre-NAC sTILs in HER2-")
#Pre, HER2-, OS
data.frame("surv.time"=OS_month[HER2_BIN=="FALSE"&TIME=="pre"], "surv.event"=OS_event[HER2_BIN=="FALSE"&TIME=="pre"], "TILs"=TILs[HER2_BIN=="FALSE"&TIME=="pre"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title= "Pre-NAC sTILs in HER2-")
#Post, HER2-, RFS
data.frame("surv.time"=RFS_month[HER2_BIN=="FALSE"&TIME=="pre"], "surv.event"=RFS_event[HER2_BIN=="FALSE"&TIME=="pre"], "TILs"=TILs[HER2_BIN=="FALSE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", title="Post-NAC sTILs in HER2-")
#Post, HER2-, OS
data.frame("surv.time"=OS_month[HER2_BIN=="FALSE"&TIME=="pre"], "surv.event"=OS_event[HER2_BIN=="FALSE"&TIME=="pre"], "TILs"=TILs[HER2_BIN=="FALSE"&TIME=="post"])->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ TILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival", title="Post-NAC sTILs in HER2-")

###Change in TILS KM Plots
#all patients
TILS[TIME=="post"]-TILS[TIME=="pre"]->tils.delta.all
as.factor(tils.delta.all>0)->sTILs
c("down/equal","increased")-> levels(sTILs)
#RFS
data.frame("surv.time"=RFS_month[TIME=="pre"], "surv.event"=RFS_event[TIME=="pre"], "TILs"=sTILs)->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ sTILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", pval.coord=c(60,0.9))
#OS
data.frame("surv.time"=OS_month[TIME=="pre"], "surv.event"=OS_event[TIME=="pre"], "TILs"=sTILs)->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ sTILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival")

#delta TILs ONLY TNBC
TILS[TNBC=="TRUE"&TIME=="post"]-TILS[TNBC=="TRUE"&TIME=="pre"]->tils.delta
as.factor(tils.delta>0)->sTILs
c("down/equal","increased")-> levels(sTILs)
#delta TILs TNBC RFS
data.frame("surv.time"=RFS_month[TNBC=="TRUE"&TIME=="pre"], "surv.event"=RFS_event[TNBC=="TRUE"&TIME=="pre"], "TILs"=sTILs)->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ sTILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Recurrence-free survival", pval.coord=c(60,0.9))
#delta TILs TNBC OS
data.frame("surv.time"=OS_month[TNBC=="TRUE"&TIME=="pre"], "surv.event"=OS_event[TNBC=="TRUE"&TIME=="pre"], "TILs"=sTILs)->km.data
fit <- survfit(Surv(surv.time, surv.event) ~ sTILs,data = km.data)
ggsurvplot(fit, data = km.data, risk.table = TRUE, pval=TRUE,xlab="Time after surgery (months)", ylab="Overall survival")
##
