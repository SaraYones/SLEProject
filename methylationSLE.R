source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(GEOquery)
library(data.table)
#BiocManager::install("MethylMix")
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
#For showing the cpg sites on tracks
#BiocManager::install("coMET")
library(coMET)
library(corrr)
#install.packages("remotes")
#remotes::install_github("tidymodels/corrr")

methylationSLE=fread("GSE118144_Normalized_data.txt")
methylationSLE=as.data.frame(methylationSLE)
rownames(methylationSLE) <- methylationSLE[,1]
methylationSLE[,1] <- NULL
cpgSites=rownames(methylationSLE)
metaMethylation=getGEO(filename="GSE118144_series_matrix.txt")
metaMethylation=as.data.frame(metaMethylation)
ageMethylation=fread("AgeMethylationStudy.txt",header = TRUE)
indicesWB=which(grepl("SLE[0-9][0-9][0-9]_WB\\.AVG_Beta|BUC[0-9][0-9][0-9]_WB\\.AVG_Beta",colnames(methylationSLE)))
indicesPV=which(grepl("SLE[0-9][0-9][0-9]_WB\\.Detection Pval|BUC[0-9][0-9][0-9]_WB\\Detection Pval",colnames(methylationSLE)))
methylationSLEPval=methylationSLE[,indicesPV]
methylationSLEtemp=methylationSLE
methylationSLE=methylationSLE[,indicesWB]
rownames(methylationSLE)<-cpgSites
annotation.table = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
PediatricCases=ageMethylation$`Sample ID`[which(ageMethylation$`Age at recruitment`<20)]
methylationSLEPediatric=methylationSLE[,which(colnames(methylationSLE)%in%paste(PediatricCases,"_WB.AVG_Beta",sep=""))]

compareMethylationProfile(methylationSLEPediatric,"CKAP4",annotation.table)

compareMethylationProfile=function(methylationSLEPediatric,geneName,annotation.table)
{
  
  geneName="ZNFX1"
  #cpgAnnotations=annotation.table[annotation.table$UCSC_RefGene_Name %in% geneName,]
 # cpgAnnotations=annotation.table[which(grepl(paste(".*",geneName,";?.*",sep=""),annotation.table$UCSC_RefGene_Name )),]
  cpgAnnotations1=annotation.table[annotation.table$UCSC_RefGene_Name %in% geneName,]
  cpgAnnotations2=annotation.table[which(grepl(paste(".*",geneName,";.*",sep=""),annotation.table$UCSC_RefGene_Name )),]
  cpgAnnotations=annotation.table[which(rownames(annotation.table) %in% unique(union(rownames(cpgAnnotations1),rownames(cpgAnnotations2)))),]
  
   methylationProfile=methylationSLEPediatric[rownames(methylationSLEPediatric)%in% rownames(cpgAnnotations),]
  #cases=methylationProfile[grepl("SLE.*",colnames(methylationProfile)),]
  cases=methylationProfile[,which(grepl("SLE.*",colnames(methylationProfile)))]
  #controls=methylationProfile[grepl("BUC.*",colnames(methylationProfile)),]
  controls=methylationProfile[,which(grepl("BUC.*",colnames(methylationProfile)))]
  classType=append(rep("controls",dim(as.data.frame(controls))[2]),rep("cases",dim(as.data.frame(cases))[2]))
  methylationProfile=t(methylationProfile)
  methylationProfile=as.data.frame(methylationProfile)
  methylationProfile$classType=as.character(classType)
  df.m <- melt(methylationProfile, id.var = "classType")
  classes=unique(as.character(df.m$classType))
  View(df.m)
#  pvalues=calculatePvalues(df.m,"Not",NULL)
#print(pvalues)
 # df.mfilter=df.m[df.m$variable %in% names(pvalues)[which(pvalues<=0.05)],]
  df.mfilter=df.m
 # View(annotation.table[which(annotation.table$UCSC_RefGene_Name==geneName&&rownames(annotation.table) %in% names(pvalues)[which(pvalues<=0.05)]),])
  
  filterAnnotation=annotation.table[which(annotation.table$UCSC_RefGene_Name==geneName),]
  #filterAnnotation2=annotation.table[which(grepl(paste(".*",geneName,";.*",sep=""),annotation.table$UCSC_RefGene_Name )),]
  #filterAnnotation=filterAnnotation[which(rownames(filterAnnotation) %in% names(pvalues)[which(pvalues<=0.05)]),]
  
  p1<-ggplot(data = df.mfilter,aes(x = fct_reorder(variable, value,.desc=TRUE), y = value,fill=classType))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+theme_bw()+ geom_point(aes(y=value, group=classType),alpha=1,size=0.2, position = position_dodge(width=0.75))
  p1<-p1+ylim(0, 1)+geom_boxplot(inherit.aes = TRUE,aes(fill=classType),alpha=0.3,outlier.size=0,lwd=0.1,stat = "boxplot") +theme(axis.text=element_text(size=8,angle = 90), axis.title=element_text(size=12),legend.title=element_text(size=10),legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "AVG beta")#+stat_compare_means(aes(group = classType),label = "p.signif")
  

  p2<-ggplot(as.data.frame(cpgAnnotations), aes(x=as.factor(unique(cpgAnnotations$UCSC_RefGene_Name)),fill=factor(cpgAnnotations$UCSC_RefGene_Group)))+geom_bar()+theme_classic()+theme(legend.title = element_blank())+labs(x = "")
  
  
#  p2<-ggplot(as.data.frame(cpgAnnotations), aes(x=strsplit((cpgAnnotations$UCSC_RefGene_Name),";"),fill=strsplit(cpgAnnotations$UCSC_RefGene_Group,";")))+geom_bar()+theme_classic()+theme(legend.title = element_blank())+labs(x = "")
 # p2<-ggplot(as.data.frame(cpgAnnotations), aes(x=strsplit((cpgAnnotations$UCSC_RefGene_Name),";"),fill=strsplit(cpgAnnotations$UCSC_RefGene_Group,";")))+geom_bar()+theme_classic()+theme(legend.title = element_blank())+labs(x = "")
  
  temp=as.data.frame(cbind(unlist(strsplit((cpgAnnotations$UCSC_RefGene_Name),";")),unlist(strsplit(cpgAnnotations$UCSC_RefGene_Group,";"))))
  ggplot(temp, aes(V2))
  p2<-ggplot(temp, aes(x=as.factor(V1),fill=as.factor(V2)))+geom_bar()+theme_classic()+theme(legend.title = element_blank() ,axis.text.y=element_blank(), axis.ticks.y=element_blank())+labs(x = "")
  
  plot_grid( as.grob(as.ggplot(p1)), as.grob(as.ggplot(p2)), labels="AUTO",label_size = 10)
}

plotGeneMethylationTrack=function(methylationSLEPediatric,methylationSLEPval,annotation.table,geneName)
{
  
  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
  configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver_Grch38.txt")
  data_config <-read.csv(configfile, quote = "", sep="\t", header=FALSE)
  geneName="ZNFX1"
  cpgAnnotations1=annotation.table[annotation.table$UCSC_RefGene_Name %in% geneName,]
  cpgAnnotations2=annotation.table[which(grepl(paste(".*",geneName,";.*",sep=""),annotation.table$UCSC_RefGene_Name )),]
  cpgAnnotations=annotation.table[which(rownames(annotation.table) %in% unique(union(rownames(cpgAnnotations1),rownames(cpgAnnotations2)))),]
  methylationProfile=methylationSLEPediatric[rownames(methylationSLEPediatric)%in% rownames(cpgAnnotations),]
  methylationPval=methylationSLEPval[rownames(methylationSLEPval)%in% rownames(cpgAnnotations),]
  methylationPval=methylationPval[match(rownames(cpgAnnotations),rownames(methylationPval)),]
  TargetDF=cbind(rownames(cpgAnnotations),cpgAnnotations$chr,cpgAnnotations$pos,methylationPval$`SLE014_WB.Detection Pval`)
  colnames(TargetDF)<-c("TargetID","CHR","MAPINFO","Pval")
  TargetDF=as.data.frame(TargetDF)
  write.xlsx(TargetDF, file = "CPGLocations.xlsx", sheetName=geneName, append=TRUE)
  
  TargetDF$CHR=as.integer(gsub("(chr([0-9]*))","\\2",as.character(TargetDF$CHR)))
  TargetDF$MAPINFO=as.integer(as.character(TargetDF$MAPINFO))
  TargetDF$TargetID=as.character(TargetDF$TargetID)
  TargetDF$Pval=as.numeric(as.character(TargetDF$Pval))
  TargetDF$Pval<-runif(dim(TargetDF)[1],min=0.1,max = 0.11)
  methylationProfile[is.na(methylationProfile)]<-0.1
  
  
  TargetDF=TargetDF[order(TargetDF$MAPINFO),]
 temp=methylationProfile
  methylationProfile=methylationProfile[match(TargetDF$TargetID,rownames(methylationProfile)),]
  data_config <-read.csv(configfile, quote = "", sep="\t", header=FALSE)
   
  comet(config.file=configfile, mydata.file=as.data.frame(TargetDF),mydata.type="dataframe",cormatrix.file=cormatrix.data.raw ,cormatrix.type="dataframe",
            print.image=FALSE,verbose=FALSE)
  
  matrix.dnamethylation <- read.delim(myinfofile, header=TRUE, sep="\t", as.is=TRUE,
                                      blank.lines.skip = TRUE, fill=TRUE)
  matrix.expression <- read.delim(myexpressfile, header=TRUE, sep="\t", as.is=TRUE,
                                  blank.lines.skip = TRUE, fill=TRUE)
  cormatrix.data.raw <- read.delim(mycorrelation, sep="\t", header=TRUE, as.is=TRUE,
                                   blank.lines.skip = TRUE, fill=TRUE)

  listmatrix.expression <- list(matrix.expression)
  listcormatrix.data.raw <- list(cormatrix.data.raw)
   
  
  listcormatrix.data.raw <- list(cormatrix.data.raw)
  
  chrom <- paste("chr",unique(TargetDF$CHR),sep="")
  start <- min(TargetDF$MAPINFO)
  end <- max(TargetDF$MAPINFO)
  gen <- "hg19"
  strand <- "*"
  BROWSER.SESSION="UCSC"
  mySession <- browserSession(BROWSER.SESSION)
  genome(mySession) <- gen
  
  genetrack <-genes_ENSEMBL(gen,chrom,start,end-200,showId=TRUE)
 # snptrack <- snpBiomart_ENSEMBL(gen, chrom, start, end,
  #                               dataset="hsapiens_snp_som",showId=FALSE,title = "SNPs ENSEMBL")
  #strutrack <- structureBiomart_ENSEMBL(chrom, start, end,
   #                                     strand, dataset="hsapiens_structvar_som")
  #clinVariant<-ClinVarMain_UCSC(gen,chrom,start,end,showId = TRUE)
  #clinCNV<-ClinVarCnv_UCSC(gen,chrom,start,end,showId = TRUE)
 # gwastrack <-GWAScatalog_UCSC(gen,chrom,start,end,showId=TRUE)
  regtrack<-regulatoryFeaturesBiomart_ENSEMBL(gen, chrom, start, end, featureDisplay = "all",
                                              datasetEnsembl = "hsapiens_regulatory_feature",
                                              title="Regulatory Features ENSEMBL")
  #geneRtrack <-GeneReviews_UCSC(gen,chrom,start,end,showId=TRUE)
  listgviz <- list(regtrack,genetrack)
  comet(config.file=configfile,mydata.file=TargetDF,
        mydata.type="dataframe",cormatrix.file= list(t(methylationProfile)),cormatrix.method="spearman",cormatrix.format="raw",cormatrix.type="listdataframe",fontsize.gviz =12, font.factor=3,tracks.gviz=listgviz,verbose=FALSE, print.image=FALSE)
 
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

