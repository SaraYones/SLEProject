#Figures Manuscript
#Before and after pruning PCA
library(ggfortify)
library(grid)
library(gridExtra)
rownames(DA13)=DA13$V1
DA13pca=SLEclassifier[,SLE_Johnson$MCFSFeatures]
temp=decisionSLE13
temp[which(temp==1)]="DA1"
temp[which(temp==3)]="DA3"

temp=removedVisits$Removal
temp[which(temp=="no")]="Not Removed"
temp[which(temp=="yes")]="Removed"
DA13pca$decision=temp
DA13pca$decision=as.factor(DA13pca$decision)
temp=removedVisits$treatment
DA13pca$treatment=temp
DA13pca=as.data.frame(DA13pca)
pca_res <- prcomp(DA13pca[,1:dim(DA13pca)[2]-1], scale. = TRUE)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p<-autoplot(pca_res, data = DA13pca, colour = "decision",shape = "treatment")
p<-p+theme



DA13pca2=SLEclassifier[rownames(dt_40_remove1),SLE_Johnson$MCFSFeatures]
indextemp=which(rownames(SLEclassifier) %in% rownames(dt_40_remove1))
DA13pca2$decision=temp[indextemp]
DA13pca2$decision=as.factor(DA13pca2$decision)
DA13pca2=as.data.frame(DA13pca2)
pca_res2 <- prcomp(DA13pca2[,1:dim(DA13pca2)[2]-1], scale. = TRUE)
p2<-autoplot(pca_res2, data = DA13pca2, colour = "decision")
p2<-p2+theme

cowplot::plot_grid(p, p2,labels = "AUTO")

#----------------------------------------------------------------------------
tempMetadata13= metadata13[rownames(dt_40_remove1),c("subject","visit")]
tempMetadata13$decision=temp[indextemp]
tempMetadata13$subject=gsub("subject: (.*)","\\1",tempMetadata13$subject)
theme2<-theme(legend.title = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p3 <- ggplot(data = tempMetadata13, aes(x = visit, y = subject, group = decision))
p3 + geom_point(aes(colour = factor(decision)))+theme2
#----------------------------------------------------------------------
#Removed objects table (Supplementary)

removedVisits=(cbind(rownames(metadata13),!(rownames(metadata13) %in% rownames(dt_40_remove1))))
colnames(removedVisits)=c("Visit_ID","Removal")
removedVisits=as.data.frame(removedVisits)
removedVisits[] <- lapply(removedVisits, as.character)
removedVisits$Removal[which(removedVisits$Removal=="TRUE")]="yes"
removedVisits$Removal[which(removedVisits$Removal=="FALSE")]="no"
decisionRemoved=dt_40_remove1[match(removedVisits$Visit_ID,rownames(dt_40_remove1)),"decisionSLE"]
treatment=metadata13[match(removedVisits$Visit_ID,rownames(metadata13)),"treatment"]
days_since_diagnosis=metadata13[match(removedVisits$Visit_ID,rownames(metadata13)),"daysSinceDiagnosis"]
nephritis_class=metadata13[match(removedVisits$Visit_ID,rownames(metadata13)),"nephritis_class:.1"]

SLEDAI=metadata13[match(removedVisits$Visit_ID,rownames(metadata13)),"sledai"]
SLEDAI=gsub("sledai: (.*)","\\1",SLEDAI)
SLEDAI=as.numeric(SLEDAI)
treatment=gsub("treatment: (.*)","\\1",treatment)
nephritis_class=gsub("nephritis_class: (.*)","\\1",nephritis_class)
decisionRemoved=as.character(decisionRemoved)
removedVisits=as.data.frame(cbind(removedVisits,decisionRemoved))
removedVisits=as.data.frame(cbind(removedVisits,treatment))
removedVisits=as.data.frame(cbind(removedVisits,days_since_diagnosis))
removedVisits=as.data.frame(cbind(removedVisits,SLEDAI))
removedVisits=as.data.frame(cbind(removedVisits,nephritis_treatment))
#---------------nephritis class------------------------------

temp=metadata13[which(rownames(metadata13)%in% rownames(dt_40_remove1)),]
#change name of coloumn to nephritisClass1
colnames(temp)[86] <- "nephritisClass1"
#distribution for nephritis class
ggplot(temp, aes(x=reorder(nephritisClass_1, nephritisClass1, function(x)-length(x))))+     geom_bar(fill='red') +  labs(x='nephritisClass1')
#How many non removed objects had DA3
temp2=removedVisits[which((removedVisits$decisionRemoved==3)&(removedVisits$Removal=="no")),]
samplesNephritis=rownames(temp[which(!(grepl("nephritis_class: NoLN|nephritis_class: Data Not Available|nephritis_class: Not Applicable",temp$nephritisClass1))),])
length(which(temp2$Visit_ID %in% samplesNephritis))

#------------------Creating PCA plots for the removed and non removed based on logistic regression results-----------------------------

temp=nephritis_class
DA13pca$treatment=temp
DA13pca=as.data.frame(DA13pca)
pca_res <- prcomp(DA13pca[,1:dim(DA13pca)[2]-1], scale. = TRUE)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p<-autoplot(pca_res, data = DA13pca, colour = "decision",shape = "treatment")
p<-p+theme


temp=nephritis_treatment
DA13pca$treatment=temp
DA13pca=as.data.frame(DA13pca)
pca_res <- prcomp(DA13pca[,1:dim(DA13pca)[2]-1], scale. = TRUE)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p<-autoplot(pca_res, data = DA13pca, colour = "decision",shape = "treatment")
p<-p+theme


days_since_diagnosis[which(days_since_diagnosis<median(days_since_diagnosis,na.rm = TRUE))]=0
days_since_diagnosis[which(days_since_diagnosis>median(days_since_diagnosis,na.rm = TRUE))]=1
days_since_diagnosis[which(days_since_diagnosis==0)]="short"
days_since_diagnosis[which(days_since_diagnosis==1)]="long"
days_since_diagnosis[which(is.na(days_since_diagnosis))]="Not Available"

temp=days_since_diagnosis


DA13pca$treatment=temp
DA13pca=as.data.frame(DA13pca)
pca_res <- prcomp(DA13pca[,1:dim(DA13pca)[2]-1], scale. = TRUE)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p<-autoplot(pca_res, data = DA13pca, colour = "decision",shape = "treatment")
p<-p+theme

days_since_diagnosis[which(days_since_diagnosis<median(days_since_diagnosis,na.rm = TRUE))]=0
days_since_diagnosis[which(days_since_diagnosis>median(days_since_diagnosis,na.rm = TRUE))]=1
days_since_diagnosis[which(days_since_diagnosis==0)]="short"
days_since_diagnosis[which(days_since_diagnosis==1)]="long"
days_since_diagnosis[which(is.na(days_since_diagnosis))]="Not Available"


SLEDAI[which(SLEDAI<=2)]=0
SLEDAI[which(SLEDAI>2)]=1
SLEDAI[which(SLEDAI==0)]="low"
SLEDAI[which(SLEDAI==1)]="high"

temp=SLEDAI


DA13pca$treatment=temp
DA13pca=as.data.frame(DA13pca)
pca_res <- prcomp(DA13pca[,1:dim(DA13pca)[2]-1], scale. = TRUE)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p<-autoplot(pca_res, data = DA13pca, colour = "decision",shape = "treatment")
p<-p+theme
#------------------Distribution of Logistic regression results-------------------
tempdistribution=removedVisits[which(removedVisits$Removal=="yes"),]
ggplot(tempdistribution, aes(x=SLEDAI))+     geom_histogram(binwidth=1, fill="#FF9999", color="#e9ecef", alpha=0.9) +  labs(x='SLEDAI')

tempdistribution=removedVisits[which(removedVisits$Removal=="yes"),]
ggplot(tempdistribution, aes(x=days_since_diagnosis))+     geom_histogram(fbinwidth=1, fill="#FF9999", color="#e9ecef", alpha=0.9) +  labs(x='days_since_diagnosis')

ggplot(tempdistribution, aes(x=treatment))+ geom_bar(binwidth=1, fill="#FF9999", color="#e9ecef", alpha=0.9) +  labs(x='histogram')

#---------------------------------------------------------------

