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

DA13pca$decision=temp
DA13pca$decision=as.factor(DA13pca$decision)
DA13pca=as.data.frame(DA13pca)
pca_res <- prcomp(DA13pca[,1:dim(DA13pca)[2]-1], scale. = TRUE)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p<-autoplot(pca_res, data = DA13pca, colour = "decision")
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


