source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(GEOquery)
library(data.table)
BiocManager::install("MethylMix")
library(MethylMix)
library(ggpubr)
library(cowplot)
library(ggplot2)
library("RColorBrewer")
library(ggplotify)
#For ggarrange
library("gridExtra")
library("ggpubr")
library(forcats)
library("easyGgplot2")
library(xlsx)
#Used with frequency and difference between cohorts plot
library(grid)
library(gridExtra)
#For legends
library(lemon)




methylationSLE=fread("GSE118144_Normalized_data.txt")
methylationSLE=as.data.frame(methylationSLE)
rownames(methylationSLE) <- methylationSLE[,1]
methylationSLE[,1] <- NULL
cpgSites=rownames(methylationSLE)
metaMethylation=getGEO(filename="GSE118144_series_matrix.txt")
metaMethylation=as.data.frame(metaMethylation)
ageMethylation=fread("AgeMethylationStudy.txt",header = TRUE)
indicesWB=which(grepl("SLE[0-9][0-9][0-9]_WB\\.AVG_Beta|BUC[0-9][0-9][0-9]_WB\\.AVG_Beta",colnames(methylationSLE)))
methylationSLE=methylationSLE[,indicesWB]
rownames(methylationSLE)<-cpgSites
annotation.table = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
PediatricCases=ageMethylation$`Sample ID`[which(ageMethylation$`Age at recruitment`<20)]
methylationSLEPediatric=methylationSLE[,which(colnames(methylationSLE)%in%paste(PediatricCases,"_WB.AVG_Beta",sep=""))]

compareMethylationProfile=function(methylationSLEPediatric,geneName,annotation.table)
{
  geneName="OTOF"
  cpgAnnotations=annotation.table[annotation.table$UCSC_RefGene_Name==geneName,]
  methylationProfile=methylationSLEPediatric[rownames(methylationSLEPediatric)%in% rownames(cpgAnnotations),]
  cases=methylationProfile[grepl("SLE.*",rownames(methylationProfile)),]
  controls=methylationProfile[grepl("BUC.*",rownames(methylationProfile)),]
  classType=append(rep("controls",dim(as.data.frame(controls))[1]),rep("cases",dim(as.data.frame(cases))[1]))
  methylationProfile=t(methylationProfile)
  methylationProfile=as.data.frame(methylationProfile)
  methylationProfile$classType=as.character(classType)
  df.m <- melt(methylationProfile, id.var = "classType")
  classes=unique(as.character(df.m$classType))
  pvalues=calculatePvalues(df.m,"Not",NULL)

  df.mfilter=df.m[df.m$variable %in% names(pvalues)[which(pvalues<=0.05)],]
  View(annotation.table[which(annotation.table$UCSC_RefGene_Name==geneName&&rownames(annotation.table) %in% names(pvalues)[which(pvalues<=0.05)]),])
  
  filterAnnotation=annotation.table[which(annotation.table$UCSC_RefGene_Name==geneName),]
  filterAnnotation=filterAnnotation[which(rownames(filterAnnotation) %in% names(pvalues)[which(pvalues<=0.05)]),]
  
  p1<-ggplot(data = df.mfilter,aes(x = fct_reorder(variable, value,.desc=TRUE), y = value,fill=classType))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+theme_bw()+ geom_point(aes(y=value, group=classType),alpha=1,size=0.2, position = position_dodge(width=0.75))
  p1<-p1+ylim(0, 1)+geom_boxplot(inherit.aes = TRUE,aes(fill=classType),alpha=0.3,outlier.size=0,lwd=0.1,stat = "boxplot") +theme(axis.text=element_text(size=8), axis.title=element_text(size=12),legend.title=element_text(size=10),legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "AVG beta")+stat_compare_means(aes(group = classType,label=sprintf("p = %5.4f", as.numeric(..p.format..))),label.y = 105,paired = FALSE ,label.x.npc="middle",size = 2)
  
  p2<-ggplot(as.data.frame(filterAnnotation), aes(x=as.factor(unique(filterAnnotation$UCSC_RefGene_Name)),fill=factor(filterAnnotation$UCSC_RefGene_Group)))+geom_bar()+theme_classic()+theme(legend.title = element_blank())+labs(x = "")
  
  pcombined=plot_grid( as.grob(as.ggplot(p1)), as.grob(as.ggplot(p2)), labels="AUTO",label_size = 10)
}


  
  
  
  
calculatePvalues=function(TiTvMatrix,flagPaired,orderNames)
  {
    TiTv=unique(as.character(TiTvMatrix$variable))
    names(TiTv)=TiTv
    TiTv=TiTv[order(factor(names(TiTv), levels = orderNames))]
    classes=unique(as.character(TiTvMatrix$classType))
    pvalues=NULL  
    for (i in 1:length(TiTv))
    {
      temp=TiTvMatrix[TiTvMatrix$variable==TiTv[i],]
      #Normality test
      result1=shapiro.test(temp[temp$classType==classes[1],"value"])
      result2=shapiro.test(temp[temp$classType==classes[2],"value"])
      #t-test
      if(result1$p.value>0.05 && result2$p.value>0.05)
      { print(paste("using t-test for",TiTv[i]))
        #print(resultDiagnosis)
        if(flagPaired=="Paired")
        {
          result=compare_means(value ~ classType, data = temp, paired = TRUE,method = "t.test")
          pvalues=append(pvalues,result$p.adj)
        }else{
          result=compare_means(value ~ classType, data = temp, paired = FALSE,method = "t.test")
          pvalues=append(pvalues,result$p.adj)
          
        }
      }else{#Wilxcon any of them is significantly differnt than normal 
        print(paste("using Wilcoxon-test for",TiTv[i]))
        if(flagPaired=="Paired")
        {
          result=compare_means(value ~ classType, data = temp, paired = TRUE,method = "wilcox.test")
          pvalues=append(pvalues,result$p.adj)
          
        }else{
          
          result=compare_means(value ~ classType, data = temp, paired = FALSE,method = "wilcox.test")
          pvalues=append(pvalues,result$p.adj)
        }
        
      }
    }
    names(pvalues)<-TiTv
    return(pvalues)
  }
  
#----------Tria with MethylMix package
cancerSite <- "OV"
targetDirectory <- paste0(getwd(), "/")

METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory)
# Processing methylation data
METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)

