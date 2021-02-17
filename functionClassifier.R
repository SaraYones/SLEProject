
Classifier <- setRefClass("Classifier",
                          #          fields = list(classifier = "data.frame",MCFSFeatures = "list", Accuracies = "list",classifierMCFSgenes="data.frame",
                          #             path="string",Accuracies_Genetic="list",classifierMCFSgenes="list",Classifier_MCFS="as.data.frame",resultRosettaMCFSGenetic="as.data.frame"
                          # ),
                          fields = list(classifier = "data.frame",flagAccuracy="character",MCFSFeatures = "character",rank="numeric", Accuracies = "numeric",AccuraciesGenetic="numeric", path="character",classifierMCFSgenes="character",Classifier_MCFS="data.frame",annotation="character",resultRosettaMCFSGenetic="list",resultRosettaMCFSJohnson="list",resultRosettaMCFSMerged="list",recalculatedResultRosettaMCFSJohnson="data.frame",recalculatedResultRosettaMCFSGenetic="data.frame",enrichment="list",kegg="gseaResult",ontology="character",numberOfFeatures="numeric",keyType="character",flagEnrichment="logical",underSample="logical"),
                          methods =  list(findAccuracies = function() {
                            classifier <<- classifier
                            MCFSFeatures<<-MCFSFeatures
                            flagAccuracy<<-flagAccuracy
                            underSample<<-underSample
                            #ontology<<--ontology
                            numberOfFeatures<<--numberOfFeatures
                            #  keyType<<--keyType
                            # kegg<<--kegg
                            # enrichment<<--enrichment
                            
                            if(flagAccuracy=="Genetic")
                            {
                              if(underSample==TRUE)
                                AccuraciesGenetic<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy,underSample = TRUE)
                              else
                                AccuraciesGenetic<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy,underSample = FALSE)
                              
                            }else if(flagAccuracy=="Johnson")
                            {
                              Sys.sleep(4)
                              if(underSample==TRUE)
                                Accuracies<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy,TRUE)
                              else
                                Accuracies<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy,FALSE)
                              
                              
                            }
                            
                            
                            
                          },
                          createAnnotatedClassifier = function() {
                            
                            Accuracies<<-Accuracies
                            AccuraciesGenetic<<-AccuraciesGenetic
                            classifierMCFSgenes<<-classifierMCFSgenes
                            Classifier_MCFS<<-Classifier_MCFS
                            resultRosettaMCFSJohnson<<- resultRosettaMCFSJohnson
                            resultRosettaMCFSGenetic<<- resultRosettaMCFSGenetic
                            recalculatedResultRosettaMCFSJohnson<<-recalculatedResultRosettaMCFSJohnson
                            recalculatedResultRosettaMCFSGenetic<<-recalculatedResultRosettaMCFSGenetic
                            underSample<<-underSample
                            #   recalculatedResultRosettaMCFS<<-recalculatedResultRosettaMCFS
                            classifier<<-classifier
                            flagAccuracy<<-flagAccuracy
                            resultRosettaMCFSMerged<<-resultRosettaMCFSMerged
                            genes<<-genes
                            if(flagAccuracy=="Johnson")
                            {
                              if(length(which(Accuracies==max(Accuracies))>1))
                                
                              {
                                temp=which(Accuracies==max(Accuracies))
                                maxAccuracy<-temp[length(temp)]*10
                              }else{
                                maxAccuracy<-which(Accuracies==max(Accuracies))*10
                              }
                            }else if(flagAccuracy=="Genetic") {
                              if(length(which(AccuraciesGenetic==max(AccuraciesGenetic))>1))
                                
                              {
                                temp=which(AccuraciesGenetic==max(AccuraciesGenetic))
                                maxAccuracy<-temp[length(temp)]*10
                              }else{
                                maxAccuracy<-which(AccuraciesGenetic==max(AccuraciesGenetic))*10
                              }
                              
                            }
                            
                            filtered_MCFS<<-getAnnotatedGenes(genes,MCFSFeatures[1:maxAccuracy])
                            #All Genes MCFS
                            
                            classifierMCFSgenes<<-append(filtered_MCFS[[1]]$geneID,filtered_MCFS[[2]])
                            
                            
                            #Ordering according to MCFS Features
                            classifierMCFSgenes<<-classifierMCFSgenes[!is.na(match(MCFSFeatures, classifierMCFSgenes))]
                            
                            
                            annotation<<-annotateInOrder(filtered_MCFS,classifierMCFSgenes)
                            
                            annotation<<-make.unique(annotation)
                            
                            
                            Classifier_MCFS<<-classifier[,append(classifierMCFSgenes,"decision")]
                            
                            
                            colnames(Classifier_MCFS)<<-append(annotation,"decision")
                            
                            if(underSample==TRUE)
                              resultRosettaMCFSGenetic<<-rosetta(Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5,reducer="Genetic",underSample =T,ruleFiltration=TRUE,ruleFiltrSupport=c(1,3) , discreteParam=3)
                            else
                              resultRosettaMCFSGenetic<<-rosetta(Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5,reducer="Genetic",underSample =F,ruleFiltration=TRUE,ruleFiltrSupport=c(1,3) , discreteParam=3)
                            
                            Sys.sleep(10)
                            # resultRosettaMCFSJohnson<<-NULL
                            if(underSample==TRUE)
                              resultRosettaMCFSJohnson<<-rosetta(Classifier_MCFS,cvNum=5,underSample=T)
                            else
                              resultRosettaMCFSJohnson<<-rosetta(Classifier_MCFS,underSample=F)
                            resultRosettaMCFSMerged<<-list("Johnson"=resultRosettaMCFSJohnson,"Genetic"=resultRosettaMCFSGenetic)
                            
                            recalculatedResultRosettaMCFSJohnson<<-recalculateRules(Classifier_MCFS,resultRosettaMCFSJohnson$main)
                            Sys.sleep(10)
                            # recalculatedResultRosettaMCFSGenetic<<-NULL
                            recalculatedResultRosettaMCFSGenetic<<-recalculateRules(Classifier_MCFS,resultRosettaMCFSGenetic$main)
                            resultRosettaMCFSMerged<<-list("JohnsonResult"=resultRosettaMCFSJohnson,"GeneticResult"=resultRosettaMCFSGenetic,"recalculatedResultJohnson"=recalculatedResultRosettaMCFSJohnson,"recalculatedResultGenetic"=recalculatedResultRosettaMCFSGenetic)
                            plotDistributionCategory(list(filtered_MCFS[[1]]$gene_type),list(c("MCFS")),paste("Distribution_Classes_FS","_",flagAccuracy,sep=""),paste(getwd(),path,sep=""),"")
                            
                          },
                          
                          
                          clusterRulesandWriteoutput = function() {
                            #recalculatedResultRosettaMCFS<<-recalculatedResultRosettaMCFS
                            resultRosettaMCFSGenetic<<-resultRosettaMCFSGenetic
                            recalculatedResultRosettaMCFSGenetic<<-recalculatedResultRosettaMCFSGenetic
                            resultRosettaMCFSJohnson<<-resultRosettaMCFSJohnson
                            recalculatedResultRosettaMCFSJohnson<<-recalculatedResultRosettaMCFSJohnson
                            path<<-path
                            enrichment<<-enrichment
                            flagEnrichment<<-flagEnrichment
                            if(flagEnrichment==FALSE)
                            {
                              tempJohnson=recalculatedResultRosettaMCFSJohnson[(which(!(recalculatedResultRosettaMCFSJohnson$supportLHS==0))),]
                              tempGenetic=recalculatedResultRosettaMCFSGenetic[(which(!(recalculatedResultRosettaMCFSGenetic$supportLHS==0))),]
                              enrichment<<-enrichment
                              
                              clusteredRulesMCFS<<-clusterRules(tempJohnson,rownames(Classifier_MCFS))
                              
                              clusteredRulesMCFSGenetic<<-clusterRules(tempGenetic,rownames(Classifier_MCFS))
                              enrichment<<-computeEnrichment()
                              
                              writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Genetic",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFSGenetic,tempGenetic,resultRosettaMCFSGenetic$main,enrichment,"enrichment",FALSE,FALSE)
                              
                              writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Johnson",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFS, tempJohnson,resultRosettaMCFSJohnson$main,enrichment ,"enrichment",FALSE,FALSE)
                            }
                            else
                            {
                              print(length(enrichment))
                              writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Genetic",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFSGenetic,tempGenetic,resultRosettaMCFSGenetic$main,enrichment,"enrichment",TRUE,FALSE)
                              
                              
                              
                            }
                            
                          }, computeEnrichment = function() {
                            
                            
                            MCFSFeatures<<-MCFSFeatures
                            classifier<<-classifier
                            enrichment<<-enrichment
                            ontology<<-ontology
                            numberOfFeatures<<-numberOfFeatures
                            keyType<<-keyType
                            classifierMCFSgenes<<-classifierMCFSgenes
                            optimalNumber<<-length(classifierMCFSgenes)
                            enrichment<<-enrichment
                            
                            enrichment1=GOenrichment(MCFSFeatures[1:numberOfFeatures],colnames(classifier)[1:length(classifier)-1],keyType,ontology)
                            enrichment2=GOenrichment(MCFSFeatures[1:optimalNumber],colnames(classifier)[1:length(classifier)-1],keyType,ontology)
                            enrichment3=GOenrichment(MCFSFeatures,colnames(classifier)[1:length(classifier)-1],keyType,ontology)
                            enrichment=list(enrichment1=enrichment1,enrichment2=enrichment2,enrichment3=enrichment3)
                            # plotEnrichment(enrichment,"Enrichment  All Genes",paste(paste(getwd(),path,sep = ""),Sys.Date(),"/","enrichment.csv",sep = ""))
                            
                          },
                          computeKEGG = function() {
                            kegg<<-kegg
                            keyType<<-keyType
                            MCFSFeatures<<-MCFSFeatures
                            numberOfFeatures<<-numberOfFeatures
                            rank<<-rank
                            geneList=NULL
                            geneList=1:numberOfFeatures
                            names(geneList)=MCFSFeatures[1:numberOfFeatures]
                            
                            eg = bitr(MCFSFeatures[1:numberOfFeatures], fromType=keyType,toType="ENTREZID", OrgDb="org.Hs.eg.db")
                            # geneList[match(eg$ENSEMBL, names(geneList))]
                            # temp=geneList[match(eg$ENSEMBL, names(geneList))]
                            if(keyType=="ENSMBL")
                            {
                              temp=eg[match(names(geneList),eg$ENSEMBL),]
                              temp=temp[which(!is.na(temp)),]
                              temp=temp[which(!is.na(temp$ENSEMBL)),]
                              temp=unique(temp)
                              
                              geneList=geneList[names(geneList)%in% temp$ENSEMBL]
                              View(geneList)
                            } else if(keyType=="SYMBOL")
                              
                            {
                              temp=eg[match(names(geneList),eg$SYMBOL),]
                              temp=temp[which(!is.na(temp)),]
                              temp=temp[which(!is.na(temp$SYMBOL)),]
                              temp=unique(temp)
                              
                              geneList=geneList[names(geneList)%in% temp$SYMBOL]
                              View(geneList)
                              
                              
                            }
                            
                            names(geneList)=temp$ENTREZID
                            
                            kegg <- gseKEGG(geneList     = sort(geneList, decreasing = TRUE),
                                            organism     = 'hsa',
                                            nPerm        = 1000,
                                            minGSSize    = 1000,
                                            pvalueCutoff = 0.8,
                                            verbose      = FALSE)
                            
                            
                          }
                          
                          )
                          
)




#-----------------Functions Classifier-------------------------------------------------
extractFeaturesBoruta=function(boruta.train)
{
  confirmed=rownames(as.data.frame(boruta.train$finalDecision[which(boruta.train$finalDecision=="Confirmed")]))
  tentative=rownames(as.data.frame(boruta.train$finalDecision[which(boruta.train$finalDecision=="Tentative")]))
  BorutaFeatures=append(confirmed,tentative)
  return(BorutaFeatures)
}

FilterFeatures=function(file,numberOfFeatures)
{
  
  features=fread(file)
  features=features[seq(1,numberOfFeatures),c("attribute","RI_norm")]
  return(list(features$attribute,features$RI_norm))
  
}
compareAccuracies=function(Decisiontable,limit,features,flag,underSample)
{
  Accuracies=NULL
  resultRosetta=NULL
  
  for(i in seq(10, limit, by = 10)){
    print(i)
    
    # resultRosettaWithUSMod7=rosetta(DecisiontableBinary,classifier="StandardVoter",discrete=TRUE,underSample = TRUE)
    #For SLE 
    #resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discreteMask=TRUE,discrete = TRUE,ruleFiltration=TRUE,ruleFiltrSupport=c(1,3))
    #FOR AML
    #Johnson
    # View(Decisiontable)
    #print(class(Decisiontable))
    #print(limit)
    #print(features)
    if(flag =="Johnson")
    { if(underSample==TRUE)
      resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter", discreteMethod="EqualFrequency",cvNum = 5,underSample=T,discreteParam=3)
    else
      resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter", discreteMethod="EqualFrequency",underSample=F,discreteParam=3)
    
    } else if(flag =="Genetic")
    {
      #resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5, discreteParam=3)
      # print("hiii")
      #    print(str(Decisiontable[,append(features[1:i],"decision")]))
      #        #Genetic
      #      View(Decisiontable[,append(features[1:i],"decision")])
      #      print(str(Decisiontable[,append(features[1:i],"decision")]))
      if(underSample==TRUE)
        resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],cvNum = 5,reducer="Genetic",underSample=T,ruleFiltration=TRUE)
      else{
        resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],cvNum = 10,reducer="Genetic",underSample=F,ruleFiltration=TRUE)
        
      }
      
    }# print(resultRosetta)
    Accuracies=append(Accuracies,resultRosetta$quality$accuracyMean)
    resultRosetta=NULL
    Sys.sleep(4)
    # print(Accuracies)
  }
  plot(Accuracies,type='l',xaxt="n",main="Accuracies")
  axis(1,at=seq(1,as.numeric(limit/10), by = 1),labels=as.character(seq(10, limit, by = 10)))
  return(Accuracies)
  
}
#For Gene expression Decision tables
prepareDT=function(DT,values)
{
  temp=NULL;
  temp=DT[,1:dim(DT)[2]-1]
  temp=as.matrix(temp)
  temp[temp=="up"]<-values[3]
  temp[temp=="down"]<-values[1]
  temp[temp=="normal"]<-values[2]
  temp=apply(temp,2,as.numeric)
  temp=as.data.frame(temp)
  decisiontemp=unlist(lapply(DT$decisionSLE,as.character))
  temp=cbind.data.frame(temp,decisiontemp)
  #temp=as.data.frame(apply(temp,2,as.factor))
  names(temp)[names(temp)=="decisiontemp"]="decisionSLE"
  
  #temp$decision=as.factor(temp$decision)
  #temp$decisionSLE=as.integer(temp$decisionSLE)
  rownames(temp)<-rownames(DT)
  return (temp)
}
clusterRules=function(result,objects)
{
  #result=filterResultRosetta13
  #objects=rownames(DescretizedDF13)
  resultMatrix=matrix(0,nrow=dim(result)[1],ncol=length(objects))
  colnames(resultMatrix)<-objects
  for(i in 1:dim(result)[1])
  {
    SUPP_SET_LHS=unlist(as.list(strsplit(as.character(result$supportSetLHS[[i]]), ",")))
    #SUPP_SET_RHS=unlist(as.list(strsplit(as.character(result$SUPP_SET_RHS[[i]]), ",")))
    #SUPP_SET=intersect(SUPP_SET_LHS,SUPP_SET_RHS)
    resultMatrix[i,which(colnames(resultMatrix) %in% SUPP_SET_LHS)]=1
    
  }
  resultMatrix=resultMatrix[,which(colSums(resultMatrix)!=0)]
  #resultMatrix=resultMatrix[,which(colSums(resultMatrix)!=-dim(resultMatrix)[1])]
  return(resultMatrix)
}

associateRulesClusters=function(clusterRules1,clusters1)
{
matrixCluster=clusterRules1
colors=unique(clusters1[,"color_cluster"])
for(i in 1:length(colors))
{
x=rownames(clusters1)[which(clusters1[,"color_cluster"] %in% colors[i])]
index=which(colnames(clusterRules1) %in% x)
print(index)
print(colors[i])
matrixCluster[,index][matrixCluster[,index] == 1]=i+1
}

return(matrixCluster)

}

rulesInCluster <- function(recal, samples,percent){
  clust <- array(0, NROW(recal))
  

    for(k in 1:NROW(recal)){
      a <- unlist(as.list(strsplit(as.character(recal$supportSetRHS[[k]]), ",")))
      a<-as.character(a)
     # print(a)
   if((length((which(a %in% samples))))>=((length(a))*percent))
   {
     
     clust[k] <- clust[k] + 1 
   }
      }
    
  
  clust <- cbind.data.frame(recal$features, recal$levels, recal$decision, clust)
  colnames(clust) <- c("FEATURES", "LEVELS", "DECISION", "n")
  
  return(clust)
}

#The distribution of Rules of the model in each of the cluster based on the length of the support set compared to the support of rules of each cluster
returnDistribution<-function(rules,clusters)
{
  colors=unique(clusters$color_cluster)
  distributionClusters=matrix(0,nrow=10,ncol=length(colors))
  distribution=NULL
  dfplot=NULL
  dfplotall=NULL
  allplots=NULL
  filt <- rules[rules$pValue < 0.05,]

    for (i in 1:length(colors))
    {
      for(j in seq(10, 100, by = 10))
        {
        color<-rulesInCluster(filt,rownames(clusters[clusters$color == as.character(colors[i]),]),(j/100))
            distribution=append(distribution,dim(color[color$n==1,])[1])
            
      }
      print(distribution)
      distributionClusters[,i]=distribution
      distribution=NULL
      dfplot$percent=seq(10, 100, by = 10)
      dfplot$frequency=distributionClusters[,i]
      dfplot$portion=dfplot$frequency/(dim(rules)[1])
      dfplotall=append(dfplotall,dfplot)
      p<-ggplot(as.data.frame(dfplot), aes(x = percent, y = frequency)) +
        geom_bar(fill = "red", stat = "identity")+labs(title=colors[i])

      allplots=append(allplots, p)
      
      
    }


 # plot_grid(allplots)
 # ggarrange(plotlist=as.list(allplots))
  return(list(distributionClusters=distributionClusters,plots=allplots,dfplotall=dfplotall))
}

getmax<-function(thelist)
{
  
  maximum=lapply(thelist, function(x) x[which.max(abs(x))])
  return(maximum)
}

metaInformation<-function(color,clusters,index,phenotype)
  
{
  print(paste("Total number of samples in cluster:",length(rownames(clusters[which(clusters$color_cluster==color),]))))
  
  print(paste("Number of samples in cluster that have phenotype value:",
  
  length(rownames(clusters)[which(rownames(clusters[which(clusters$color_cluster==color),]) %in% rownames(phenotype[index,]))])))
 
  #print(rownames(clusters)[which(rownames(clusters[which(clusters$color_cluster==color),]) %in% rownames(phenotype[index,]))])
  
  print(paste("Total Number of samples that have value of phenotype: ",length(index)))

  
}
#We are not using this variant of the function

rulesInCluster <- function(recal, samples,percent){
  clust <- array(0, NROW(recal))
  
  for(i in 1:length(samples)){
    
    for(k in 1:NROW(recal)){
      a <- unlist(as.list(strsplit(as.character(recal$supportSetRHS[[k]]), ",")))
      a<-as.character(a)
      for(j in 1:length(a)){
        
        if(samples[i] == a[j]){
          clust[k] <- clust[k] + 1 
        }
      }
    }
  }
  clust <- cbind.data.frame(recal$features, recal$levels, recal$decision, clust)
  colnames(clust) <- c("FEATURES", "LEVELS", "DECISION", "n")
  
  return(list(All=as.data.frame(clust),filt=as.data.frame(clust[clust$n>=(max(clust$n)*percent),])))
}

#phenotype=metadata13[which(rownames(metadata13) %in% rownames(matrixCluster)),]

sigRulesTableNC=createSigRulesTable(indexNC,66,phenotype$neutrophil_count,phenotype,matrixCluster)
sigRulesTableNP=createSigRulesTable(indexNP,66,phenotype$neutrophil_percent,phenotype,matrixCluster)
sigRulesTableLC=createSigRulesTable(indexLC,66,phenotype$lymphocyte_count,phenotype,matrixCluster)
sigRulesTableLP=createSigRulesTable(indexLP,66,phenotype$lymphocyte_percent,phenotype,matrixCluster)
sigRulesTablec3=createSigRulesTable(indexc3,66,phenotype$c3,phenotype,matrixCluster)
sigRulesTablec4=createSigRulesTable(indexc4,66,phenotype$c4,phenotype,matrixCluster)

write.xlsx(sigRulesTablec3,"Phenotype-SignificantRulesClusters.xlsx", sheetName="c3",append = TRUE)
write.xlsx(sigRulesTablec4,"Phenotype-SignificantRulesClusters.xlsx", sheetName="c4",append = TRUE)
write.xlsx(sigRulesTableNC,"Phenotype-SignificantRulesClusters.xlsx", sheetName="Neutrophil count",append = TRUE)
write.xlsx(sigRulesTableNP,"Phenotype-SignificantRulesClusters.xlsx", sheetName="Neutrophil Percentage",append = TRUE)
write.xlsx(sigRulesTableLC,"Phenotype-SignificantRulesClusters.xlsx", sheetName="Lymphocyte count",append = TRUE)
write.xlsx(sigRulesTableLP,"Phenotype-SignificantRulesClusters.xlsx", sheetName="Lymphocyte Percentage",append = TRUE)
write.xlsx(phenotype[,c("c3","c4","neutrophil_count","neutrophil_percent","lymphocyte_count","lymphocyte_percent")],"Phenotype-SignificantRulesClusters.xlsx", sheetName="Phenotypes",append = TRUE)

#Plot distribution of Rules 
#Use binary clusters

plotDistribution=NULL
distributionClusters=returnDistribution(filt,clusters)
plotDistribution$clusters=unique(clusters$color_cluster)
plotDistribution$clusters=as.factor(plotDistribution$clusters)
plotDistribution$distribution=paste(getmax(resultdistribution$dfplotall[3]),getmax(resultdistribution$dfplotall[6]),getmax(resultdistribution$dfplotall[9]),getmax(resultdistribution$dfplotall[12]),getmax(resultdistribution$dfplotall[15]))

plotDistribution$distribution=as.numeric(unlist(str_split(plotDistribution$distribution," ")))
ggplot(as.data.frame(plotDistribution), aes(x = clusters, y = distribution,fill = factor(clusters)))+geom_bar(stat = "identity")+scale_fill_manual("", values = c("red" = "#a30019ff","blue"="#3a5187ff","orange" = "#d2722eff", "green" = "#12664cff","purple"="#877194ff"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),text = element_text(size=12),axis.line = element_line(colour = "grey"),panel.background = element_blank())+geom_hline(yintercept=0.25, linetype="dotted",size = 1)




allclusterssplot=data.frame(distributionClusters$dfplotall[2],distributionClusters$dfplotall[5],distributionClusters$dfplotall[8],distributionClusters$dfplotall[11],distributionClusters$dfplotall[14])
colnames(allclusterssplot)<-unique(clusters$color_cluster)
rownames(allclusterssplot)=as.character(str_split(unlist(distributionClusters$dfplotall[1])," "))
trial=melt(allclusterssplot, id.var = "classtype")
colnames(trial)<-c("classtype","clusters","value")
trial$classtype <- factor(trial$classtype,levels = c("10", "20", "30", "40","50","60","70","80","90","100"))
p1<-ggplot2.barplot(data=as.data.frame(trial),xName='classtype', yName='value',groupName='clusters', groupColors=c('red','green','orange','purple','blue'),position=position_dodge(),backgroundColor="white", color="black", xtitle="percent", ytitle="number of rules", mainTitle="",removePanelGrid=TRUE,removePanelBorder=TRUE,axisLine=c(0.5, "solid", "black")) 


#Correlation between common values


correlationIndices=cor(as.matrix(apply(phenotype[intersectIndices,c("c3","c4","neutrophil_count","neutrophil_percent","lymphocyte_count","lymphocyte_percent","wbc","MC","MP","esr","hgb","hct","mcv","mch","mchc","rdw","mpv","cr","alb","ds_dna","ast","ald","alt","ldh","sledai","age")],2,as.numeric)))
phenotypeCorrelationPlot=corrplot(cor(as.matrix(apply(phenotype[intersectIndices,c("c3","c4","neutrophil_count","neutrophil_percent","lymphocyte_count","lymphocyte_percent","wbc","MC","MP","esr","hgb","hct","mcv","mch","mchc","rdw","mpv","cr","alb","ds_dna","ast","ald","alt","ldh","sledai","age")],2,as.numeric))))
write.xlsx(correlationIndices, paste("ClustersHeatmaps/","CorrelationPhenotypes.xlsx",sep=""), sheetName="CorrelationPhenotypes", append=TRUE)
write.xlsx(intersectIndices, paste("ClustersHeatmaps/","CorrelationPhenotypes.xlsx",sep=""), sheetName="IndicesofPhenotypes", append=TRUE)
write.xlsx(apply(phenotype[intersectIndices,c("c3","c4","neutrophil_count","neutrophil_percent","lymphocyte_count","lymphocyte_percent","wbc","monocyte_count","monocyte_percent","esr","hgb","hct","mcv","mch","mchc","rdw","mpv","cr","alb","ds_dna","ast","ald","alt","ldh","sledai","age")],2,as.numeric), paste("ClustersHeatmaps/","CorrelationPhenotypes.xlsx",sep=""), sheetName="PhenotyeMatrix", append=TRUE)

#Dominance
DominanceTable=createDominanceTable(consolidated,colors=c("red","blue","green","orange","purple"),indicesPhenotypes=indicesPhenotypes,indicesPhenotypesCategorical=indicesPhenotypesCategorical)
DominanceTable=DominanceTable[,-dim(DominanceTable)[2]]
rowSums(DominanceTable)
DominanceTablePercentage=DominanceTable/rowSums(DominanceTable)[row(DominanceTable)]

DominanceTablePercentage=DominanceTablePercentage*100
DominanceTablePercentage$Clusters=rownames(DominanceTablePercentage)



dominant=read_excel("October2020SLE/ComparisonsRules/dominantBinary.xlsx",sheet = "dominanceContinous")
rows=dominant$Clusters
dominant=dominant[,2:dim(dominant)[2]]
rownames(dominant)=rows
dominant=as.data.frame(dominant)
Clustersdominant=dominant$Clusters
dominant=dominant[,1:(dim(dominant)[2]-1)]
dominant=dominant*100
dominant$Clusters=Clustersdominant
#trial=melt(DominanceTablePercentage, id.var = "Clusters")
trial=melt(dominant, id.var = "Clusters")
trial$Clusters=as.factor(trial$Clusters)

p1<-ggplot(trial, aes(x = variable, y = value ,fill = factor(Clusters))) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = scales::percent_format())+coord_flip()+guides(fill=guide_legend(title=""))+
 # scale_fill_manual("", values = c("red" = "#661100","blue"="#0072B2","orange" = "#D55E00", "green" = "#117733","purple"="#CC79A7"))+xlab("")+ylab("")
 scale_fill_manual("legend", values = c("red" = "#a30019ff","blue"="#3a5187ff","orange" = "#d2722eff", "green" = "#12664cff","purple"="#877194ff"))+xlab("")+ylab("")


dominantCategorical=read_excel("October2020SLE/ComparisonsRules/dominantBinary.xlsx",sheet = "dominanceCategorical")
rows=dominantCategorical$Clusters
dominantCategorical=dominantCategorical[,2:dim(dominantCategorical)[2]]
rownames(dominantCategorical)=rows
dominantCategorical=as.data.frame(dominantCategorical)
Clustersdominant=dominantCategorical$Clusters
dominantCategorical=dominantCategorical[,1:(dim(dominantCategorical)[2]-1)]
dominantCategorical=dominantCategorical*100
dominantCategorical$Clusters=Clustersdominant
#trial=melt(DominanceTablePercentage, id.var = "Clusters")
trial=melt(dominantCategorical, id.var = "Clusters")
trial$Clusters=as.factor(trial$Clusters)
p2<-ggplot(trial, aes(x = variable, y = value ,fill = factor(Clusters))) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = scales::percent_format())+coord_flip()+guides(fill=guide_legend(title=""))+
 # scale_fill_manual("", values = c("red" = "#661100","blue"="#0072B2","orange" = "#D55E00", "green" = "#117733","purple"="#CC79A7"))+xlab("")+ylab("")
  scale_fill_manual("legend", values = c("red" = "#a30019ff","blue"="#3a5187ff","orange" = "#d2722eff", "green" = "#12664cff","purple"="#877194ff"))+xlab("")+ylab("")
cowplot::plot_grid(p1, p2,labels = "AUTO")
  #trial$value=trial$value*100
colnames(trial)<-c("clusters","classtype","value")
ggplot(data=trial, aes(x=value,y=classtype,fill=factor(clusters))) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual("legend", values = c("red" = "darkred","blue"="darkblue","orange" = "darkorange", "green" = "darkgreen","purple"="purple"))+xlab("dominance")+ylab("")


write.xlsx(DominanceTablePercentage,"dominantBinary.xlsx", sheetName="dominance",append = TRUE)

#---------------------------Permutation Test--------------------------
clusterHeatmapPer=as.character(clusters$color_cluster)
names(clusterHeatmapPer)=rownames(clusters)
heatmapPermutations(recal_John40_remove1,dt_40_remove1,rownames(dt_40_remove1),clusterHeatmapPer)

heatmapPermutations<-function(rules,decisionSystem,samples,clusters)
{  #decisionSys1=decisionSystem
  pvalue=0
  for(i in 1:100)
  {
    
  index=shuffle(rules$decision)
  newdecision=rules$decision[index]
  rules$decision=newdecision
 # decisionSys1$decision=decisionSystem$decision[shuffle(decisionSystem$decision)]
 #df <- decisionSys1 #[recal_John40_remove1$DECISION == 1,] # & recal_John40_remove1$DECISION == 3,]
   clusterRules <- clusterRule(rules, 
                              samples)
  
  
  #df <- df[rownames(df) %in% colnames(clusterRules),]
 
  #clusterRules <- clusterRules[, order(colnames(clusterRules))]
  # View(clusterRules)
 # print("I am here")
  #print(colnames(clusterRules))
  #rownames(df) =colnames(clusterRules)
  rule_cluster2=heatmap.F(t(clusterRules),
                          colors = c("white", "blue"), 
                          cutoffmethod = "number",
                          cutoff = 5,
                          distmethod='binary')
  #order rule_Cluster2 based on 1
  
  rule_cluster2=rule_cluster2[order(factor(names(rule_cluster2), levels = names(clusters)))]
  
  #compute percentage of difference between binary and pearson correlation
  View(rule_cluster2)
  #View(clusters)
comparison=length(which((clusters==rule_cluster2)))/length(clusters)
print(paste("comparison is",as.character(comparison)))


  if(comparison<0.5)
    pvalue=pvalue+1
 # print("Hello")
 # decisionSys1=decisionSystem
print( paste("pvalue for",as.character(i),"is:",as.character(pvalue)))
  
  }
print(pvalue/100)
}

ModelPermutations<-function(rules,decisionSystem,samples,clusters)
{  #decisionSys1=decisionSystem
  pvalue=0
  comparison=0
  a=0
  s=0
  m=0
  for(i in 1:1000)
  {
    
    index=shuffle(decisionSystem$decision)
    newdecision=decisionSystem$decision[index]
    decisionSystem$decision=newdecision
    
    result=rosetta(decisionSystem,classifier="StandardVoter",discrete = TRUE)
    #View(clusters)
    comparison= append(comparison,result$quality$accuracyMean)
   # print(paste("comparison is",as.character(comparison)))
    print(i)
    
    
    
    
   
  }
  
  a <- mean(comparison)
   s <-sd(comparison)
   n <- 500
   error <- qnorm(0.975)*s/sqrt(n)
   left <- a-error
  right <- a+error
  
#  if(sd(comparison)>0.810357 )
   # valuesleft=comparison<left
   # valuesright=comparison>right
 # pvalue=length(valuesleft)+length(valuesright)
  if(0.810357< left | 0.810357 > right)
    pvalue=0.05
 
    #pvalue=pvalue/500
  # print("Hello")
  # decisionSys1=decisionSystem
  print( paste("pvalue for",as.character(i),"is:",as.character(pvalue)))
  
  
  print(pvalue)
  return(comparison)
}


valuesClusterPhenotypes<-function(phenotype1,phenotypeName1,indicesPhenotypes1,clusters1)
{
  colors=unique(clusters1$color_cluster)
  for (i in 1:length(colors))
 #   for (i in 1:1)
  {
    
    tempRowscluster=rownames(clusters1)[clusters1$color_cluster==as.character(colors[i])]
    for(j in 1:length(indicesPhenotypes1))
   # for(j in 1:6)
    {
      print(phenotypeName1[j])
                tempphenotype=as.data.frame(phenotype1[,as.character(phenotypeName1[j])])
              rownames(tempphenotype)<-rownames(phenotype1)
              colnames(tempphenotype)<-phenotypeName1[j]
View(tempphenotype)
tempindices=as.list(indicesPhenotypes1[[j]])
print(tempindices)
tempindices=tempindices[[1]]
#print(tempindices)
#print(tempindices[[1]])
tempPhenoRows=rownames(tempphenotype)[tempindices]
tempPhenoRows=tempRowscluster[which(tempRowscluster%in%tempPhenoRows)]

tempphenotype=as.data.frame(tempphenotype[tempPhenoRows,1])
rownames(tempphenotype)<-tempPhenoRows
colnames(tempphenotype)<-phenotypeName1[j]

               # tempPhenoRows=rownames(tempphenotype)[which(rownames(tempphenotype)[tempindices] %in% tempRowscluster)]
               
#View(tempphenotype)
#tempPhenoRows=rownames(tempphenotype)[which(tempPhenoRows%in%tempRowscluster)]
            #print(tempPhenoRows)
                #tempValues=tempphenotype[which(rownames(tempphenotype)[tempindices] %in% tempPhenoRows),1]
       #        tempValues=tempphenotype[tempPhenoRows,1]
               # print(tempValues)
              #  temp=as.data.frame(cbind(tempPhenoRows,as.character(tempValues)))
    #            View(tempphenotype)
           write.xlsx(tempphenotype, paste("ClustersHeatmaps/","ValuesPhenoPerCluster.xlsx",sep=""), sheetName=paste(colors[i],"-",phenotypeName1[j],sep=""), append=TRUE)
                
    }
    
    
  }
  
  
  
}


valuesClusterPhenotypes<-function(phenotype1,phenotypeName1,indicesPhenotypes1,clusters1)
{
  colors=unique(clusters1$color_cluster)
  for (i in 1:length(colors))
    #   for (i in 1:1)
  {
    
    tempRowscluster=rownames(clusters1)[clusters1$color_cluster==as.character(colors[i])]
    for(j in 1:length(indicesPhenotypes1))
      # for(j in 1:6)
    {
      print(phenotypeName1[j])
      tempphenotype=as.data.frame(phenotype1[,as.character(phenotypeName1[j])])
      rownames(tempphenotype)<-rownames(phenotype1)
      colnames(tempphenotype)<-phenotypeName1[j]
      View(tempphenotype)
      tempindices=as.list(indicesPhenotypes1[[j]])
      print(tempindices)
      tempindices=tempindices[[1]]
      #print(tempindices)
      #print(tempindices[[1]])
      tempPhenoRows=rownames(tempphenotype)[tempindices]
      tempPhenoRows=tempRowscluster[which(tempRowscluster%in%tempPhenoRows)]
      
      tempphenotype=as.data.frame(tempphenotype[tempPhenoRows,1])
      rownames(tempphenotype)<-tempPhenoRows
      colnames(tempphenotype)<-phenotypeName1[j]
      
      # tempPhenoRows=rownames(tempphenotype)[which(rownames(tempphenotype)[tempindices] %in% tempRowscluster)]
      
      #View(tempphenotype)
      #tempPhenoRows=rownames(tempphenotype)[which(tempPhenoRows%in%tempRowscluster)]
      #print(tempPhenoRows)
      #tempValues=tempphenotype[which(rownames(tempphenotype)[tempindices] %in% tempPhenoRows),1]
      #        tempValues=tempphenotype[tempPhenoRows,1]
      # print(tempValues)
      #  temp=as.data.frame(cbind(tempPhenoRows,as.character(tempValues)))
      #            View(tempphenotype)
      write.xlsx(tempphenotype, paste("ClustersHeatmaps/","ValuesPhenoPerCluster.xlsx",sep=""), sheetName=paste(colors[i],"-",phenotypeName1[j],sep=""), append=TRUE)
      
    }
    
    
  }
  
  
  
}

checkRedundantPhenotypes<-function(pathSigRules,indicesPhenotypes,clusters,filt,phenotypeTable,phenotypeName)
{
  
  colors=unique(clusters$color_cluster)
  namesPheno=names(indicesPhenotypes)
  
for (i in 1:length(colors))
     # for (i in 1:1)
  {
    sigPhenoRules<- read_excel(paste(pathSigRules,"/SigPhenoRules",colors[i],"-Binary.xlsx",sep=""), sheet =colors[i] )
  rules=unique(sigPhenoRules$'Rule Number')
  
  for(j in 1:length(rules))
  # for(j in 1:1)
    {
     print(paste("Rule",rules[j]))
     # print(j)
    # 
 
    temp=sigPhenoRules[which(sigPhenoRules$'Rule Number'==rules[j]),]
     # print(temp)
     # 
    
   phenotypes=temp$phenotype
   supportSet=unlist(as.list(strsplit(as.character(filt$supportSetRHS[[as.numeric(rules[j])]]), ",")))
   # print("SupportSet")
   #  print(length(supportSet))
   # print(supportSet)

   supported=supportSet[which(supportSet %in%  rownames(clusters[clusters$color_cluster==colors[i],]))]
   # print("Supported")
   #  print(length(supported))
   # print(supported)
   notSupported=supportSet[which(!(supportSet %in% rownames(clusters[clusters$color_cluster==colors[i],])))]
   
   dependant=append(rep(1,(length(supported))),rep(0,(length(notSupported))))
   
 #model=matrix(NA,nrow = length(dependant),ncol=length(indicesPhenotypes))
   model=matrix(NA,nrow = length(dependant),ncol=length(phenotypes))
 
 #colnames(model)=names(indicesPhenotypes)
   colnames(model)=phenotypes
 rownames(model)=append(supported,notSupported)
 # print(rownames(model))
 # print(colnames(model))
 model=phenotypeTable[rownames(model),phenotypes]
 #model=phenotypeTable[rownames(model),names(indicesPhenotypes)]
 tempmodel=as.data.frame(model)
 model=as.data.frame(model)
 #View(model)
 model=apply(model,2,as.numeric)
 
 model=as.data.frame(model)
# View(model)
 model$decision=as.factor(dependant)
 #Removing random effects
# model$subject=as.factor(c(1:length(rownames(model))))
 #model$decision <- relevel(model$decision, "0")
 model=as.data.frame(model)

# print(class(model))
#View(model)

 #model <- as.data.frame(model[!apply(model, 1, function(x) {any(x == NA)}),])
# print("I am here")
model=as.data.frame(na.omit(model))
# View(model)
if(dim(model)[1]>1)
{
if(dim(tempmodel)[2]>1)

{
  #tempmodel=apply(model[,1:(((dim(model)[2])-1))],2,function(x) scale(x))
  tempmodel=model[,1:(((dim(model)[2])-1))]
#tempmodel=apply(model[,1:(((dim(model)[2])-1))],2,function(x) scale(x))
model[,1:(((dim(model)[2])-1))]=tempmodel
#model[,1:(((dim(model)[2])-1))]=tempmodel
model=as.data.frame(model)
tempmodel=as.data.frame(tempmodel)

#Using logistic regression from here
# 
#  #print(class(model$hct))
#  #print(class(model$rdw))
# # View(model)
# 
# #View(model)
# 
#  dimension=dim(tempmodel)[2]
#  print(dimension)
#  
#  
# 
#  
#  #tempmodel=apply(model[,1:(((dim(model)[2])-1)-1)],2,function(x) scale(x))
#  #model[,1:(((dim(model)[2])-1)-1)]=tempmodel
#  #model=as.data.frame(model)
#  #tempmodel=as.data.frame(tempmodel)
#  
#   f <- as.formula(
#      paste("decision", paste(paste(append(colnames(tempmodel[,1:dimension]),"(1|subject)"),sep = ""), collapse = " + "), sep = " ~ "))
#    
# # f <- as.formula(paste("decision", paste(colnames(tempmodel[,1:dimension]), collapse = " + "),
#  #        sep = " ~ "))
# 
#  print(f)
#  
#  #View(model)
#  
#  x=as.data.frame(prop.table(table(model$decision)))
#  print(class(x[1,"Freq"]))
#  print(x)
#  value=0
# 
# if(x[1,"Freq"]<0.2 | x[2,"Freq"]<0.2 )
# {
# 
#     y=as.data.frame(table(model$decision))
# 
#     #print(class(y$'0'))
#   if(y[1,"Freq"]>y[2,"Freq"])
#   { print("smthg")
# 
# value = y[1,"Freq"]
# modelOver= ROSE::ovun.sample(decision ~., data = get("model", sys.frame(1)), method = "over",N = get("value", sys.frame(1))*2)$data
# 
#   }
#   else
#   {
# 
#   value=y[2,"Freq"]
#   modelOver= ROSE::ovun.sample(decision ~., data = get("model", sys.frame(1)), method = "over",N = get("value", sys.frame(1))*2)$data
#   
# 
#   }
#   #  modelOver=ovun.sample(decision ~ ., data = model, method = "over",N = value*2)$data
# 
# }
#  View(model)
#  
#    if(value==0)
#  {
#     modelOver=model
#   }
# View(modelOver)
#  #glm.fit <- glmer(f, data =modelOver, family = binomial, control=glmerControl(optimizer="bobyqa", tolPwrss=1e-10,optCtrl = list(maxfun = 10000000)))#,control = glmerControl(optimizer = "bobyqa"),
#                # nAGQ = 10)
#  glm.fit <- glmer(f, data =modelOver, family = binomial, control=glmerControl(tolPwrss=1000,optimizer = "bobyqa"))#control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
#  
#  print(check_collinearity( glm.fit))
#  x<-check_collinearity( glm.fit)
#  print(x$VIF)
 
#To there
tempprint=model[,!(colnames(model)%in% c("decision"))]
rownames(tempprint)=append(supported,notSupported)[as.numeric(rownames(model[,!(colnames(model)%in% c("decision","subject"))]))]
#print(tempprint)
print(model)
print("Number of relevant samples are:")
print(dim(tempprint)[1])


 correlation=cor(as.matrix(model[,!(colnames(model)%in% c("decision"))]))
 print(correlation)
 corrplot(correlation)
 average=NULL;
 standardDV=NULL;
 minimum=NULL;
 maximum=NULL;
 
 x=(model[which(model$decision==1),])
 y=(model[which(model$decision==0),])
 
 average$supported=apply(x[,1:dim(x)[2]-1], 2, mean)

 average$Notsupported=apply(y[,1:dim(y)[2]-1], 2, mean)

 standardDV$supported=apply(x[,1:dim(x)[2]-1], 2, sd)

 standardDV$Notsupported=apply(y[,1:dim(y)[2]-1], 2, sd)
 
 minimum$supported=apply(x[,1:dim(x)[2]-1],2,min)
 minimum$Notsupported=apply(y[,1:dim(y)[2]-1],2,min)
 
 maximum$supported=apply(x[,1:dim(x)[2]-1],2,max)
 maximum$Notsupported=apply(y[,1:dim(y)[2]-1],2,max)
 
  print("Average")
  print(average)
  print("Standard Deviation")
  print(standardDV)
  print("Minimum")
  print(minimum)
  print("Maximum")
  print(maximum)
  

 df.m=reshape2::melt(model,id.var="decision")
 View(df.m)
print(p<-ggplot(as.data.frame(df.m) ,aes(x=variable, y=value, fill=decision)) +
  geom_boxplot())
  
 #boxplot(df.m$value)

}else {
  

  print(paste("there is only 1 phenotype that is significant to rule",rules[j]))
}
}else{
print(paste("Correlation cannot be computed because the samples do not align for rule",rules[j]))
 
 #print(summary(glm.fit))
 
}
 
 
 
   }
 }
  

  
}

phenotypeComparison<-function(percent,clusters,filt,phenotype,index,phenotypeTable,phenotypeName,ClusterColor,flag,flagCont)
  {
    #rulesInCluster(filt,rownames(clusters[clusters$color == as.character(colors[i]),]),(j/100))
   # colors=unique(clusters$color_cluster)
    colors=ClusterColor
    distribution=NULL
    pvalue=NULL
    indicesRules=NULL
    indicesRulesSig=NULL
    pvalues=NULL
    tempSheet=NULL
    
   print(phenotypeName)
   
   for (i in 1:length(colors))
     # for (i in 1:1)
    {
     if(flag!=TRUE)
     {
       color<-rulesInCluster(filt,rownames(clusters[clusters$color_cluster == as.character(colors[i]),]),(percent/100))
       rulesCluster=as.data.frame(color[color$n==1,])
     
      # pvalues=matrix(0,nrow=1,ncol=length(rownames(clusters[clusters$color_cluster==colors[i],])))
       pvalues=matrix(0,nrow=1,ncol=length(rownames(rulesCluster)))
       
       colnames(pvalues)<-rownames(rulesCluster)
    #   pvalues=as.data.frame(pvalues)
       # View(pvalues)
       View(rulesCluster)
      
      for(j in 1:dim(rulesCluster)[1])
     #  for(j in 1:3)
          {
       
        ruleNumber=as.numeric(rownames(rulesCluster)[j])
       
        # print(ruleNumber)
        supportSet=unlist(as.list(strsplit(as.character(filt$supportSetRHS[[ruleNumber]]), ",")))
       # print("SupportSet")
        #  print(length(supportSet))
        # print(supportSet)
       
        supported=supportSet[which(supportSet %in%  rownames(clusters[clusters$color_cluster==colors[i],]))]
        # print("Supported")
        #  print(length(supported))
       #  print(supported)
        notSupported=supportSet[which(!(supportSet %in% rownames(clusters[clusters$color_cluster==colors[i],])))]
        # print("notSupported")
        # print(length(notSupported))
         #print(notSupported)
        supportedValues=phenotype[which(rownames(phenotypeTable[index,]) %in% supported)]
        
        #print(supportedValues)
        notSupportedValues=phenotype[which(rownames(phenotypeTable[index,]) %in% notSupported)]
        #result1=shapiro.test(supportedValues)
        #result2=shapiro.test(notSupportedValues)
     if(length(supportedValues)!=0 && length(notSupportedValues)!=0)
      {
        classtype=append(rep("supported",length(supportedValues)),rep("notsupported",length(notSupportedValues)))
        values=append(supportedValues,notSupportedValues)
      #  print("Values=")
      #  print(values)
        df1=as.data.frame(cbind(classtype,values))
        # View(df1)
        df1$values=as.numeric(as.character(df1$values))
        df1$classtype=as.factor(df1$classtype)
        # print(class(df1$values))
        # print(class(df1$classtype))
        # print(class(df1))
        # 
        colnames(df1)<-c("classtype","values")
     
    #print(df1)
        if (flagCont==TRUE)
           {
          result=compare_means(values ~ classtype, data = df1, paired = FALSE,method = "wilcox.test")
          pvalues[,which(colnames(pvalues)==as.character(ruleNumber))]=result$p.adj
         # print("Continous")
          
          p=result$p.adj
        # print(p)
           
        }else{
          # print("Categorical")
          # print(supportedValues)
          # print(notSupportedValues)
         # print(table(as.factor(supportedValues), as.factor(notSupportedValues)))
         # result=fisher.test(df1)
         #For Debugging here
           #View(df1)
      #  print(table(df1))
          # print( pairwise.chisq.test(df1$classtype, df1$values, p.adjust.method = "bonferroni"))
       if(length(unique(df1$values))!=1)
       {
           result=fisher.test(table(df1))
      #  print(result)
          pvalues[,which(colnames(pvalues)==as.character(ruleNumber))]=result$p.value
          p=result$p.value
         # print(class(p))
       }else{
         
         pvalues[,which(colnames(pvalues)==as.character(ruleNumber))]=1
         p=1
       }
          
           
        }
            # print(result)
            # print(ruleNumber)
   if(flag!=TRUE)
   {
            if(p<=0.05)
            {indicesRulesSig=append(indicesRulesSig,ruleNumber)
        #    print("I am significant")
             
            indicesRules=append(indicesRulesSig,ruleNumber)
            }
     }else if(length(notSupported)==0)
     {
       
       indicesRulesSig=append(indicesRulesSig,ruleNumber)
       pvalues[,which(colnames(pvalues)==as.character(ruleNumber))]=0
      indicesRules=append(indicesRulesSig,ruleNumber)
     }
   else{
     pvalues[,which(colnames(pvalues)==as.character(ruleNumber))]=1
          
   }
     } else{
          
     indicesRules=append(indicesRules,ruleNumber)
        }
            classtype=NULL
            values=NULL
           #pvalues=NULL
      }
       
      }
     #  print("I am here")
       
       #Heatmap
      # print(rownames(clusters[clusters$color_cluster==colors[i],]))
      # print(length(indicesRules))
      # View(pvalues)
 
    
   if(flag==TRUE)
    {
      indicesRules=rownames(filt)
   }
     # print("indicesSigRules")
     # print(indicesRules)
      # heatmap4=matrix(0,nrow=length(indicesRulesSig),ncol=length(rownames(clusters[clusters$color_cluster==colors[i],])))
     if(length(indicesRules)!=0)
     {
      
      heatmap4=matrix(0,nrow=length(indicesRules),ncol=length(rownames(clusters[clusters$color_cluster==colors[i],])))
       
           # View(heatmap4)
       heatpmap4=as.data.frame(heatmap4)
       #heatmap4=apply(heatmap4,2,as.character)
     # rownames(heatmap4)<-indicesRulesSig
      rownames(heatmap4)<-indicesRules
         colnames(heatmap4)<-rownames(clusters[clusters$color_cluster==colors[i],])
       
     
       #for(j in 1:length(indicesRulesSig))
         
        
        outputSigRules=NULL;
         
         # print("This is indicesRules")
        #print(indicesRules)
        
        for(j in 1:length(indicesRules))
       {
         
       
         # print("SIGNIFICANT RULE NUMBER ")
        ruleNumber=as.numeric(indicesRules[j])
        # ruleNumber=as.numeric(indicesRulesSig[j])
      #  print("this is the ruleNumber")
      #print(ruleNumber)
         supportSet=unlist(as.list(strsplit(as.character(filt$supportSetRHS[[ruleNumber]]), ",")))
         # print(supportSet)
         supported=supportSet[which(supportSet %in% rownames(clusters[clusters$color_cluster==colors[i],]))]
        # print(supported)
          supportedValues=phenotype[which(rownames(phenotypeTable[index,]) %in% supported)]
         # print(supportedValues)
         
         #Check here
         SupportRowNames=supported[which(supported %in% rownames(phenotypeTable[index,]))]
         #For Debugging here
       #  supportedupdated= supported[which(supported %in% SupportRowNames)]
         ## print(SupportRowNames)
          ##print(length(SupportRowNames))
        
      #   View(heatmap4)
         #print(SupportRowNames %in% colnames(heatmap4))
         #print(length(heatmap4[j,supportedupdated]))
       ##  print(length(supportedValues))
        # heatmap4[j,supportedupdated]=supportedValues
      
           
         heatmap4[j,SupportRowNames]=supportedValues
         
       
  
         
       }
   
        
   #tiff(paste("ClustersHeatmaps/",percent,"Sig","-",phenotypeName,"-",colors[i], Sys.Date(),".tiff",sep=""), units="in", width=10, height=5, res=400)
   
#QWhen wanted to save as image     tiff(paste("ClustersHeatmaps/",percent,"-",phenotypeName,"-",colors[i], Sys.Date(),".tiff",sep=""), units="in", width=10, height=5, res=400)
       # View(heatmap4) 
      
    #   paste(rownames(heatmap4),pvalues[match(rownames(heatmap4),rownames(pvalues))])
    
 
       # print( paste(rownames(heatmap4),"-",pvalues[,match(rownames(heatmap4),colnames(pvalues))],sep=""))
         if (flag!=TRUE)
         {
           
         tempRowNames<-paste(rownames(heatmap4),"-",pvalues[,match(rownames(heatmap4),colnames(pvalues))],sep="")
       #  print(filt$features[as.numeric(rownames(heatmap4))])
        # print(filt$features[as.numeric(as.character(rownames(heatmap4)))])
         repeatPhenotypeName=rep(phenotypeName,length(rownames(heatmap4)))
        outputSigRules=as.data.frame(cbind(repeatPhenotypeName,rownames(heatmap4),as.character(filt$features[as.numeric(rownames(heatmap4))]),as.character(filt$levels[as.numeric(rownames(heatmap4))]),
                             as.character(filt$decision[as.numeric(rownames(heatmap4))]),as.character(pvalues[,match(rownames(heatmap4),colnames(pvalues))])))
       colnames(outputSigRules)<-c("phenotype","Rule Number","Genes","Gene Values","decision","pvalue")
       #For Debugging here
       ##View(outputSigRules)
       #print(class(outputSigRules$pvalue))
       #print(as.numeric(as.character(outputSigRules$pvalue)))
       #print(class(as.numeric(as.character(outputSigRules$pvalue))))
       outputSigRules=outputSigRules[which(as.numeric(as.character(outputSigRules$pvalue))<=0.05),]
      if(!dim(outputSigRules)[1]<1)
       write.xlsx(outputSigRules, paste("ClustersHeatmaps/BinaryResults/",ClusterColor,"-SigRules.xlsx",sep=""),sheetName = phenotypeName, append=TRUE)
      
         #  write.csv(outputSigRules,paste("ClustersHeatmaps/","SigRules-",phenotypeName,"-",ClusterColor,".csv",sep=""))
      #tempSheet=rbind(tempSheet,outputSigRules)
             }
         else
         {
           tempRowNames<-rownames(heatmap4)
        #For Debugging here
           #   print(tempRowNames)
         }
       # print(class(heatmap4))
      #  heatmap4=as.data.frame(heatmap4)
        
       # heatmap4 %<>% mutate_if(is.character,as.numeric)
        #heatmap4=as.numeric(heatmap4)
        #print("dim before apply")
        #print(dim(heatmap4))
       # heatmap4=as.data.frame(heatmap4)
      #Dealing with only 1 rule in the heatmap
        if(dim(heatmap4)[1]==1)
        {
          heatmap4=apply(heatmap4,1,as.numeric)
        }else{
          heatmap4=apply(heatmap4,2,as.numeric)
        }
        # if(dim(heatmap4)[2]==1) 
      #For Debugging here
        #  View(heatmap4)
        #  print("Hello i am dimensions")
         # print(dim(heatmap4))
           if(dim(heatmap4)[2]==1) 
            {
              heatmap4=t(heatmap4)
              heatmap4=as.matrix(heatmap4)
            }   
       rownames(heatmap4)<-tempRowNames
       #Dealing with only 1 rule in the heatmap
       heatmap4=as.matrix(heatmap4)
     if(dim(heatmap4)[1]==1|dim(heatmap4)[2]==1)
     {
       print("There is only 1 significant rule. No heatmap to draw")
     }else if(sum(heatmap4)==0)
       {
         
        print("The heatmap values are all 0s won't be plotted")
       }
       
       
       else{
         heatmap.2(heatmap4,dendrogram='none', Rowv=FALSE, Colv=TRUE,trace='none',srtCol=45)
  #When wanted to save as image     dev.off()
     }
       indicesRulesSig=NULL
       indicesRules=NULL
       pvalues=NULL
       
     }
   else{
     print("There were no rules significant to this cluster")
   }
   }
  #  if(!is.null(tempSheet))
  # {
  #    tempSheet=as.data.frame(tempSheet)
  #  write.xlsx(tempSheet, paste("ClustersHeatmaps/BinaryResults/",ClusterColor,"-SigRules.xlsx",sep=""),sheetName = "SigRules-Phenotypes", append=TRUE)
  #  }
  # return(tempSheet)
   
}

# phenotypeComparisonCluster<-function(clusters,indices,phenotypeTable,ClusterColor,flagCont,confounderCheck)
# {
# 
#   namesIndices=names(indices)
#   #print(namesIndices)
# 
#  # clusters=subset(clusters, select=-c(rule_cluster))
#   tempBoxPlot=NULL
# 
# 
# 
#  # print(phenotypeName)
# 
# 
#   for(i in 1:length(indices))
#  # for(i in 1:1)
#   {
# 
#     if(flagCont==TRUE)
#     {
# 
# 
#       tempPhenotypeTable=phenotypeTable[as.numeric(unlist(as.list(indices[i]))),]
#       tempClusters=clusters[which(rownames(clusters)%in% rownames(tempPhenotypeTable)),]
#      # View(tempPhenotypeTable)
#       #View(tempClusters)
#       temp=tempPhenotypeTable[rownames(tempClusters),]
#       View(temp)
#       temp=as.data.frame(temp)
#       tempClusters=as.data.frame(tempClusters)
#       tempPhenotypeTable=as.data.frame(tempPhenotypeTable)
#       View(temp)
#       View(tempClusters)
#      # print(dim(temp))
#       # print(rownames(temp))
#       # print(rownames(tempClusters))
#       # print(factor(rownames(temp), levels = rownames(tempClusters)))
#       # print(class(tempClusters))
#       # print(class(temp))
# 
#      # temp=temp[order(factor(rownames(temp), levels = rownames(tempClusters))),]
#      temp=temp[match(rownames(tempClusters), (rownames(temp))), ]
#      temp=as.data.frame(temp)
#     # View(temp)
#     #  View(tempClusters)
#     # print("passed")
#      # tempBoxPlot$colors=tempClusters$color_cluster
#       #tempBoxPlot$values=temp$namesIndices[i]
#      # print(length(tempClusters$color_cluster))
#      #print(temp[,namesIndices[i]])
#      print(namesIndices[i])
#       tempBoxPlot=cbind.data.frame(as.character(tempClusters$color_cluster),as.character(temp[,namesIndices[i]]))
# 
#       #View(tempBoxPlot)
#       rownames(tempBoxPlot)=rownames(tempClusters)
# 
#       tempBoxPlot=as.data.frame(tempBoxPlot)
# 
#        print(tempBoxPlot)
#        colnames(tempBoxPlot)<-c("color_cluster","values")
#        tempBoxPlot$values=as.numeric(as.character(tempBoxPlot$values))
#        tempBoxPlot$color_cluster=as.character(tempBoxPlot$color_cluster)
#        print(tempBoxPlot)
#      #  print(class(tempBoxPlot$values))
#       p1<-ggplot(data = tempBoxPlot,aes(x = color_cluster, y = values,fill=color_cluster))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+theme_bw()+ geom_point(aes(y=values, group=color_cluster),alpha=1,size=0.2, position = position_dodge(width=0.75))
#       p1<-p1+geom_boxplot(inherit.aes = TRUE,aes(fill=color_cluster),alpha=0.3,outlier.size=0,lwd=0.1,stat = "boxplot")+
#        theme(axis.text=element_text(size=8), axis.title=element_text(size=8),legend.title=element_text(size=10),legend.text=element_text(size=9))+scale_fill_manual(values=c("blue","green","orange","purple", "red"))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "Frequency of substitutions")+stat_compare_means(aes(group = color_cluster))
#       print(p1)
#       statsClusters=tempBoxPlot
#       statsClusters$color_cluster=as.factor(statsClusters$color_cluster)
#       color_cluster=levels(statsClusters)
#       print(pairwise.t.test(statsClusters$values, statsClusters$color_cluster, p.adj = "holm"))
# 
#     print(  TukeyHSD(aov(statsClusters$values  ~  statsClusters$color_cluster)) )
# 
# }
# 
# 
#   else{
# 
#   #  print(as.numeric(unlist(as.list(indices[28]))))
#     tempPhenotypeTable=phenotypeTable[as.numeric(unlist(as.list(indices[i]))),]
#     tempClusters=clusters[which(rownames(clusters)%in% rownames(tempPhenotypeTable)),]
#     # View(tempPhenotypeTable)
#     #View(tempClusters)
# #    View(temp)
#    # View(tempPhenotypeTable)
#     #View(tempClusters)
# 
#     temp=tempPhenotypeTable[rownames(tempClusters),]
#     #View(temp)
#     temp=as.data.frame(temp)
#     tempClusters=as.data.frame(tempClusters)
#     tempPhenotypeTable=as.data.frame(tempPhenotypeTable)
# 
# 
# 
#     # print(dim(temp))
#     # print(rownames(temp))
#     # print(rownames(tempClusters))
#     # print(factor(rownames(temp), levels = rownames(tempClusters)))
#     # print(class(tempClusters))
#     # print(class(temp))
# 
#     # temp=temp[order(factor(rownames(temp), levels = rownames(tempClusters))),]
#     temp=temp[match(rownames(tempClusters), (rownames(temp))), ]
#     temp=as.data.frame(temp)
# 
#     # print("passed")
#     # tempBoxPlot$colors=tempClusters$color_cluster
#     #tempBoxPlot$values=temp$namesIndices[i]
#     # print(length(tempClusters$color_cluster))
#     #print(temp[,namesIndices[i]])
#     print(namesIndices[i])
#    # print(temp[,namesIndices[28]])
#     tempBoxPlot=cbind.data.frame(as.character(tempClusters$color_cluster),as.character(temp[,namesIndices[28]]))
# 
#     #View(tempBoxPlot)
#     rownames(tempBoxPlot)=rownames(tempClusters)
# 
#     tempBoxPlot=as.data.frame(tempBoxPlot)
# 
#      # print(tempBoxPlot)
#     colnames(tempBoxPlot)<-c("color_cluster","values")
#     tempBoxPlot$values=as.numeric(as.character(tempBoxPlot$values))
#     tempBoxPlot$color_cluster=as.character(tempBoxPlot$color_cluster)
#     print(tempBoxPlot)
#    # colnames(tempClusters)<-c("color_cluster","values")
#     #tempClusters$values=as.numeric(tempClusters$values)
# 
#    print(table(tempBoxPlot))
#    if(dim(table(tempBoxPlot))[2]>1)
#     print(fisher.test(table(tempBoxPlot), workspace = 2e8))
#    else
#      paste("There are no values for",as.character(namesIndices[[i]]),"for all clusters")
# 
#    # result=chisq.test(table(tempClusters))
#     #  print(result)
#       balloonplot(table(tempBoxPlot), main =namesIndices[[i]], xlab ="", ylab="",
#                   label = FALSE, show.margins = FALSE)
#    # pvalues[,which(colnames(pvalues)==as.character(ruleNumber))]=result$p.value
#     #p=result$p.value
# 
#    # print(result)
#     # p1<-ggplot(data=tempClusters, aes(x=color_cluster,y=values)) +
#     #   geom_bar(position="dodge",stat="identity") +scale_fill_manual("legend", values = c("red" = "darkred","blue"="darkblue","orange" = "darkorange", "green" = "darkgreen","purple"="purple"))+xlab("")+ylab("")
#     # print(p1)
#     #
# 
#       }
# 
#   tempBoxPlot=NULL
#   }
# return(temp)
#   }
phenotypeComparisonCluster<-function(clusters,indices,phenotypeTable,ClusterColor,flagCont,confounderCheck)
{
  
  
     namesIndices=names(indices)
  #   #print(namesIndices)
  #   
  #  # clusters=subset(clusters, select=-c(rule_cluster)) 
  #   tempBoxPlot=NULL
  # 
  # 
  #   
  #  # print(phenotypeName)
  # 
     FinalPlot=NULL;
     tempClusters=NULL
  
  for(i in 1:length(indices))
   # for(i in 1:1)
  {
    
    if(flagCont==TRUE)
    {    
      
      
      x=unlist(indices[[i]])
      names(x)=NULL
      print(namesIndices[[i]])
      
      tempPhenotype=phenotypeTable[x,]
     
      
      tempClusters=clusters[which(rownames(clusters) %in% rownames(tempPhenotype)),]
      
      y=which(rownames(clusters) %in% rownames(tempPhenotype))
      
      tempClusters=as.data.frame(tempClusters)
      rownames(tempClusters)<-rownames(clusters)[y]
      temp=tempPhenotype[(rownames(tempClusters)),]
      temp=temp[match(rownames(tempClusters), (rownames(temp))), ]
      View(temp)
      print(rownames(tempClusters))
      print(rownames(phenotypeTable))
      
  print(which(rownames(phenotypeTable) %in% rownames(tempClusters)))
      View(phenotypeTable[which(rownames(phenotypeTable) %in% rownames(tempClusters)),namesIndices[[i]]])
      
      tempClusters$values=temp[which(rownames(temp) %in% rownames(tempClusters)),namesIndices[[i]]]
      
      print(tempClusters)
      
      tempClusters$values=as.numeric(as.character(tempClusters$values))
      tempClusters$phenotype=as.character(rep(as.character(namesIndices[i]),length(tempClusters$values)))
        
        colnames(tempClusters)<-c("color_cluster","values","phenotype")
        tempClusters$values=as.numeric(tempClusters$values)
      FinalPlot=rbind.data.frame(FinalPlot,tempClusters)
      p1<-ggplot(data = tempClusters,aes(x = color_cluster, y = values,fill=color_cluster))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+theme_bw()+ geom_point(aes(y=values, group=color_cluster),alpha=1,size=0.2, position = position_dodge(width=0.75))
      p1<-p1+geom_boxplot(inherit.aes = TRUE,aes(fill=color_cluster),alpha=0.3,outlier.size=0,lwd=0.1,stat = "boxplot")+
        theme(axis.text=element_text(size=8), axis.title=element_text(size=8),legend.title=element_text(size=10),legend.text=element_text(size=9))+scale_fill_manual(values=c("red" = "#a30019ff","blue"="#3a5187ff","orange" = "#d2722eff", "green" = "#12664cff","purple"="#877194ff"))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "Frequency of substitutions")+stat_compare_means(aes(group = color_cluster))
      print(p1)
      statsClusters=tempClusters
      statsClusters$color_cluster=as.factor(statsClusters$color_cluster)
      color_cluster=levels(statsClusters)
      print(pairwise.t.test(statsClusters$values, statsClusters$color_cluster, p.adj = "holm"))
      
      print(  TukeyHSD(aov(statsClusters$values  ~  statsClusters$color_cluster)) )
      #  df.m=reshape2::melt(model,id.var="decision")
      # View(df.m)
      # print(p<-ggplot(as.data.frame(df.m) ,aes(x=variable, y=value, fill=decision)) +
      #        geom_boxplot())
      # View(model)
      # print(p<-ggplot(model ,aes(x=classtype, y=values)) +
      #   geom_boxplot())
      # result=compare_means(values ~ classtype, data = model, paired = FALSE,method = "wilcox.test")
      #   # print("Continous")
      
      # p=result$p.adj
      # print ("pvalue is =")
      # print(p)
      #boxplot(df.m$value)
    }
    #  values=NULL
    # classtype=NULL
    
    else{
      
      x=unlist(indices[[i]])
      names(x)=NULL
     print(namesIndices[[i]])
      
      tempPhenotype=phenotypeTable[x,]
      
      tempClusters=clusters[which(rownames(clusters) %in% rownames(tempPhenotype)),]
      
      y=which(rownames(clusters) %in% rownames(tempPhenotype))
      
      tempClusters=as.data.frame(tempClusters)
      rownames(tempClusters)<-rownames(clusters)[y]
      
      
      
      tempClusters$values=phenotypeTable[which(rownames(phenotypeTable) %in% rownames(tempClusters)),namesIndices[[i]]]
      tempClusters=tempClusters[,c("color_cluster","values")]
     # print(tempClusters)
      
      colnames(tempClusters)<-c("color_cluster","values")
    
      #View(tempClusters)
      #print("I want to check the data type of values")
      #print(class(tempClusters$values))
    #  tempClusters$values=as.numeric(tempClusters$values)
      tempClusters$color_cluster=as.character(tempClusters$color_cluster)
      #If you want to check confounders per cluster compared to all other clusters
      if(!is.na(confounderCheck))
      {
      #  print("I am here horrayyyyyy")
       # print(tempClusters$color_cluster)
       #print( tempClusters$color_cluster[tempClusters$color_cluster!=confounderCheck])
      tempClusters$color_cluster[which(tempClusters$color_cluster!=confounderCheck)]="other"
      
    #  print(tempClusters$color_cluster)
      }
   # View(tempClusters)
      if(dim(table(tempClusters))[2]<= 15 &&dim(table(tempClusters))[1]<= 3 &&dim(table(tempClusters))[2]>1)
     print(fisher.test(table(tempClusters), workspace = 2e8))
      else if(dim(table(tempClusters))[2]< 1 )
        paste("There are no values for",as.character(namesIndices[[i]]),"for all clusters")
      else
        paste("There are too many different classes and levels to perform fisher's exact test")
      # result=chisq.test(table(tempClusters))
      #  print(result)
      balloonplot(table(tempClusters), main =namesIndices[[i]], xlab ="", ylab="",
                  label = FALSE, show.margins = FALSE)
      # pvalues[,which(colnames(pvalues)==as.character(ruleNumber))]=result$p.value
      #p=result$p.value
      
      # print(result)
      # p1<-ggplot(data=tempClusters, aes(x=color_cluster,y=values)) +
      #   geom_bar(position="dodge",stat="identity") +scale_fill_manual("legend", values = c("red" = "darkred","blue"="darkblue","orange" = "darkorange", "green" = "darkgreen","purple"="purple"))+xlab("")+ylab("")
      # print(p1)
      # 
      
    }
    
    
  }
     return(FinalPlot)
}




createSigRulesTable=function(index1,ruleNumber,phenotype1,phenotypeTable,matrixCluster1)

{
  range=NULL
  average=NULL
  supportlist=NULL
  numberpheno=NULL
  print(index1)
x=apply(matrixCluster1[index1,],2,function (x) as.integer(x))
#phenotype c3

y=as.data.frame(cbind(as.numeric(phenotype1[index1]),x))


colnames(y)=append("decision",paste("V",c(1:ruleNumber), sep=""))
View(y)
# y=x
# y<-lapply(y,as.integer)
# #Draw histograms of all explanatory variables to check the values are distributed kind of the same accross classes
# hist(as.numeric(x$V5))
# #Check Multicollinearity
# j=cor(as.matrix(y))
# corrplot(cor(y))
# 
# 
# lm3 <- lm(decision ~  factor(V1) + factor(V2) + factor(V3) , data = x)
# j=cor(as.matrix(y))


# specifications of how to model,
# coming from somewhere else
outcome <- "decision"
variables <-paste("V",c(1:ruleNumber), sep="")

# our modeling effort, 
# fully parameterized!
f <- as.formula(
  paste(outcome, 
        paste(paste("factor(",variables,")",sep = ""), collapse = " + "), 
        sep = " ~ "))

# f <- as.formula(
#   paste(outcome, 
#         paste(variables, collapse = " + "), 
#         sep = " ~ "))
print(f)
# mpg ~ cyl + disp + hp + carb
lm3 <- lm(f , data = y)
result=summary(lm3)
print(result)
# values=names(which(result$coefficients[,4]<=0.1))
values=names(result$coefficients[,4])
print("here")
print(values)
#values=gsub("factor\\(.*\\)([0-3])","\\1",values)
#print(values)
#indexother=which(result$coefficients[,4]<=0.1)
indexother=1:dim(result$coefficients)[1]
indexother=indexother[2:length(indexother)]
indexother=as.integer(indexother)
#print(indexother)
#print(gsub("factor\\(V(.*)\\)[0-6]","\\1",values[2:length(values)]))
sigRulesIndex=gsub("factor\\(V(.*)\\)[0-6]","\\1",values[2:length(values)])
print(sigRulesIndex)
print(values)
sigclusters=gsub("factor\\(V.*\\)([0-6])","\\1",values[2:length(values)])
print(sigclusters)
sigRulesIndex=gsub("V(.*)","\\1",sigRulesIndex)
sigRulesIndex=as.integer(sigRulesIndex)
print(sigRulesIndex)
#colors=colors[as.integer(sigclusters)-1]



sigRulesTable=cbind(sigRulesIndex,result_John40_remove1$main$features[sigRulesIndex],result_John40_remove1$main$decision[sigRulesIndex],result_John40_remove1$main$levels[sigRulesIndex],colors[as.integer(sigclusters)-1],result$coefficients[,4][indexother])
 sigRulesTable=as.data.frame(sigRulesTable)
colnames(sigRulesTable)<-c("Rule-Number","Rule-genes","Rule-decision","Rule-ValuesGenes","Cluster","Pvalue")

#print(dim(sigRulesTable)[1])
#View(sigRulesTable)
heatmap4=matrix(0,nrow=(dim(sigRulesTable[as.numeric(as.character(sigRulesTable$Pvalue))<=0.05,])[1]),ncol=length(rownames(matrixCluster1)))
heatmap4=apply(heatmap4,2,as.character)
#heatmap4=sapply(heatmap4,as.numeric)
colnames(heatmap4)<-rownames(matrixCluster1)
rownames(heatmap4)<-(sigRulesTable[as.numeric(as.character(sigRulesTable$Pvalue))<=0.05,"Rule-Number"])
View(heatmap4)
for(i in 1:dim(sigRulesTable)[1])
{

  indexRule=intersect(which(as.character(recal_John40_remove1$features) %in% as.character(sigRulesTable$`Rule-genes`[i]) ),
                      which(as.character(recal_John40_remove1$levels) %in% as.character(sigRulesTable$`Rule-ValuesGenes`[i])))
  print(indexRule)
  SUPP_SET_LHS=unlist(as.list(strsplit(as.character(recal_John40_remove1$supportSetLHS[[indexRule]]), ",")))
  print(SUPP_SET_LHS)
  supportlist=append(supportlist,paste(as.character(SUPP_SET_LHS),",",collapse=''))
  temp=phenotype1[index1]
  valuepheno=temp[which(rownames(phenotypeTable)[index1] %in% SUPP_SET_LHS)]
  valuepheno=as.numeric(valuepheno)
  #heatmap4[i,rownames(phenotypeTable[which(rownames(phenotypeTable)[index1] %in% SUPP_SET_LHS),])]=valuepheno
  range=append(range,paste(paste(as.character(min(valuepheno)),"-",sep=""),as.character(max(valuepheno)),sep=""))
  numberpheno=append(numberpheno,length(valuepheno))
  average=append(average,as.character(round(sum(valuepheno)/length(valuepheno),2)))

}
sigRulesTable=cbind(sigRulesTable,supportlist,range,numberpheno,average)
colnames(sigRulesTable)<-c("Rule-Number","Rule-genes","Rule-decision","Rule-ValuesGenes","Cluster","Pvalue","Support list","Range","Number of Clinical values","Average")


significantTable=sigRulesTable[as.numeric(as.character(sigRulesTable$Pvalue))<=0.05,]
View(significantTable)
View(sigRulesTable)
supportlist=NULL
print("significant table rule number")
print(significantTable$`Rule-Number`)
print(dim(as.data.frame(significantTable)))
for(i in 1:dim(significantTable)[1])
{

  # indexRule=intersect(which(as.character(recal_John40_remove1$features) %in% as.character(significantTable$`Rule-genes`[i]) ),
  #                     which(as.character(recal_John40_remove1$levels) %in% as.character(significantTable$`Rule-ValuesGenes`[i])))
#  print(indexRule)
  print("SIGNIFICANT RULE NUMBER ")
  print(as.numeric(as.character(significantTable$`Rule-Number`[i])))
  SUPP_SET_LHS=unlist(as.list(strsplit(as.character(recal_John40_remove1$supportSetLHS[as.numeric(as.character(significantTable$`Rule-Number`[i]))]), ",")))
 print("Support set")
  print(SUPP_SET_LHS)
 # supportlist=append(supportlist,as.character(SUPP_SET_LHS[which(SUPP_SET_LHS %in% rownames(phenotypeTable)[index1])]))
  temp=phenotype1[index1]
  valuepheno=temp[which(rownames(phenotypeTable)[index1] %in% SUPP_SET_LHS)]
  SupportRowNames=rownames(phenotypeTable)[index1]
  #print("first")
  #print(SupportRowNames)
  SupportRowNames=SupportRowNames[which(as.character(SupportRowNames) %in% as.character(SUPP_SET_LHS))]
print("Second")
print(SupportRowNames)
  supportlist=unique(append(supportlist,as.character(SupportRowNames)))
  # print("value pheno")
  # print(length(valuepheno))
  # print("support list")
  # print(SUPP_SET_LHS)
  # print(length(SUPP_SET_LHS))
  # print(isUnique(SUPP_SET_LHS))
  # valuepheno=as.numeric(valuepheno)
  #print("here is the i")
  #print(i)
 # print(SupportRowNames)
  print("valuePheno")
  print(valuepheno)
  print(class(SupportRowNames))
  heatmap4[i,SupportRowNames]=valuepheno
  print(heatmap4[i,SupportRowNames])
 # heatmap4[i,SUPP_SET_LHS]=valuepheno
  
}
View(heatmap4)
#apply(heatmap4,2,function(x) x[is.na(x)] <- 0)
View(significantTable)
print("i am here")
#View(heatmap4)
#apply(as.matrix(heatmap4),2,function(x) x[is.na(x)] <- 0)
#print(supportlist)
#print(as.character(unlist(unique(supportlist))))

#View(heatmap4[,as.character(unique(supportlist))]) 
# #rownames(sigRulesTable)<-NULL
# #row.names(sigRulesTable)<-as.character(sigRulesIndex)
# 
# 
# View(temp)



#colnames(sigRulesTable)<-c("Rule-Number","Rule-genes","Rule-decision","Rule-ValuesGenes","Cluster","Pvalue")
#print(length(unique(supportlist)))
print(dim(clusters)[1])

return(list(sigRulesTable=as.data.frame(sigRulesTable),heatmap4=as.data.frame(heatmap4),supportlist=supportlist))
#return(list(sigRulesTable=as.data.frame(sigRulesTable),heatmap4=as.data.frame(heatmap4)))

}


#------------ CUSTOMIZE VISUNET -----------------------------------------

#filt <- result_John40_remove1$main[which(result_John40_remove1$main$pValue <= 0.05),]
filt<- recal_John40_remove1[which(result_John40_remove1$main$pValue <= 0.05),]
vis_out <- visunet(filt, type = "RDF")

red <- rulesInCluster(filt, rownames(clusters[clusters$color == "red",]),0.20)
green<- rulesInCluster(filt, rownames(clusters[clusters$color == "green",]),0.20)
blue <- rulesInCluster(filt, rownames(clusters[clusters$color == "blue",]),0.20)
purple <- rulesInCluster(filt, rownames(clusters[clusters$color == "purple",]),0.20)
orange <- rulesInCluster(filt, rownames(clusters[clusters$color == "orange",]),0.20)

id_green <- getId(green[green$n==1,])
id_red <- getId(red[red$n==1,])
id_blue <- getId(blue[blue$n==1,])
id_purple<- getId(blue[purple$n==1,])
id_orange<- getId(orange[orange$n==1,])

cust <- vis_out$all$nodes
cust$color <- "grey"
cust$color[which(as.character(cust$id) %in% id_blue)] <- "blue"
cust$color[which(as.character(cust$id) %in% id_green)] <- "green"
cust$color[which(as.character(cust$id) %in% id_red)] <- "red"
cust$color[which(as.character(cust$id) %in% id_purple)] <- "purple"
cust$color[which(as.character(cust$id) %in% id_orange)] <- "orange"
custList <- list(nodes = cust, CustCol =  c("color"))
#vis_out <- visunet(result_John40_remove1$main, type = "RDF", CustObjectNodes = custList)
vis_out <- visunet(recal_John40_remove1, type = "RDF", CustObjectNodes = custList)


getId <- function(rules){
  id_tot <- array()
  
  for(i in 1:NROW(rules)){
    features <- unlist(strsplit(as.character(rules$FEATURES[i]), ","))
    values <- unlist(strsplit(as.character(rules$LEVELS[i]), ","))
    id <- array("", length(features))
    
    for(j in 1:length(features)){
      id[j] <- paste(features[j], "=", values[j], sep = "")
    }
    id_tot <- c(id_tot, id)
  }
  
  return(unique(id_tot)[-1])
}

#-------------------------Gene enrichment analysis based on the visunet---------------------------------------------------------
#orange cluster
#First customize for each color visunet before running this 
orangegenes=custList$nodes[which(custList$nodes$color=="orange"),]$label
redgenes=custList$nodes[which(custList$nodes$color=="red"),]$label
greengenes=custList$nodes[which(custList$nodes$color=="green"),]$label
bluegenes=custList$nodes[which(custList$nodes$color=="blue"),]$label
purplegenes=custList$nodes[which(custList$nodes$color=="purple"),]$label
D1genes=unique(vis_out$'1'$nodes$label)
D3genes=unique(vis_out$'3'$nodes$label)
DA13genes=unique(append(D1genes,D3genes))

bpred=GOenrichment(redgenes,colnames(DA13)[1:dim(DA13)[2]-1],"SYMBOL","BP")
bpgreen=GOenrichment(greengenes,colnames(DA13)[1:dim(DA13)[2]-1],"SYMBOL","BP")
bpblue=GOenrichment(bluegenes,colnames(DA13)[1:dim(DA13)[2]-1],"SYMBOL","BP")
bppurple=GOenrichment(purplegenes,colnames(DA13)[1:dim(DA13)[2]-1],"SYMBOL","CC")
bporange=GOenrichment(orangegenes,colnames(DA13)[1:dim(DA13)[2]-1],"SYMBOL","BP")
bpD1=GOenrichment(D1genes,colnames(DA13)[1:dim(DA13)[2]-1],"SYMBOL","BP")
bpD3=GOenrichment(D3genes,colnames(DA13)[1:dim(DA13)[2]-1],"SYMBOL","BP")
bpDA13=GOenrichment(DA13genes,colnames(DA13)[1:dim(DA13)[2]-1],"SYMBOL","BP")

colorGSE=append(append(append(append(append(append(rep("red",dim(bpred@result)[1]),rep("green",dim(bpgreen@result)[1])),rep("blue",dim(bpblue@result)[1])),
    rep("orange",dim(bporange@result)[1])),rep("D1",dim(bpD1@result)[1])),rep("D3",dim(bpD3@result)[1])),
rep("DA13",dim(bpDA13@result)[1]))

finalGSE=cbind(rbind(bpred@result,bpgreen@result,bpblue@result,bporange@result,bpD1@result,bpD3@result,bpDA13@result),colorGSE)

write.xlsx(finalGSE,"GEA.xlsx",sheetName="sheet1")

p1=cnetplot(bpred, node_label="all") + ggtitle("Red")+ theme(plot.title = element_text(hjust = 0.5))
p2=cnetplot(bpblue, node_label="all") + ggtitle("Blue")+ theme(plot.title = element_text(hjust = 0.5))
p3=cnetplot(bpgreen, node_label="all") + ggtitle("Green")+ theme(plot.title = element_text(hjust = 0.5))
#p4=cnetplot(bppurple, node_label="all") 
p5=cnetplot(bporange, node_label="all") + ggtitle("Orange")+ theme(plot.title = element_text(hjust = 0.5))
cowplot::plot_grid(p1, p2, p3, p5,labels = "AUTO")

pDA1<-cnetplot(bpD1, node_label="all") + ggtitle("DA1")+ theme(plot.title = element_text(hjust = 0.5))
pDA3<-cnetplot(bpD3, node_label="all") + ggtitle("DA3")+ theme(plot.title = element_text(hjust = 0.5))

pDA13<-cnetplot(bpDA13, node_label="all") + ggtitle("DA13")+ theme(plot.title = element_text(hjust = 0.5))

cowplot::plot_grid(pDA1, pDA3, pDA13,labels = "AUTO")


d1 <- dotplot(bpred, showCategory=30) + ggtitle("Red")
d2 <- dotplot(bpgreen, showCategory=30) + ggtitle("Green")
d3 <- dotplot(bpblue, showCategory=30) + ggtitle("Blue")
d5 <- dotplot(bporange, showCategory=30) + ggtitle("Orange")
cowplot::plot_grid(d1, d2, d3, d5,labels = "AUTO")
dDA1<-dotplot(bpD1, showCategory=30) + ggtitle("DA1")
dDA3<-dotplot(bpD3, showCategory=30) + ggtitle("DA3")
dDA13<-dotplot(bpDA13, showCategory=30) + ggtitle("DA13")

cowplot::plot_grid(dDA1, dDA3,labels = "AUTO")

  


qplot(as.character(RosettaEnrichment$features), as.factor(as.character(BP)),xlab = "Feature", ylab = "Biological Process")+theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1),plot.background=element_rect(fill="white", colour=NA))#+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

#----------------------Create Dominance Table--------------------------------------------------

redfile=read_excel_allsheets("/Users/saryo614/Desktop/Projects/SLEProject/ClustersHeatmaps/BinaryResults/red-SigRules.xlsx")
bluefile=read_excel_allsheets("/Users/saryo614/Desktop/Projects/SLEProject/ClustersHeatmaps/BinaryResults/blue-SigRules.xlsx")
greenfile=read_excel_allsheets("/Users/saryo614/Desktop/Projects/SLEProject/ClustersHeatmaps/BinaryResults/green-SigRules.xlsx")
orangefile=read_excel_allsheets("/Users/saryo614/Desktop/Projects/SLEProject/ClustersHeatmaps/BinaryResults/orange-SigRules.xlsx")
purplefile=read_excel_allsheets("/Users/saryo614/Desktop/Projects/SLEProject/ClustersHeatmaps/BinaryResults/purple-SigRules.xlsx")

#write.xlsx(redfile, paste("ClustersHeatmaps/BinaryResults/","red","-consolidatedUpdated","-SigRules.xlsx",sep=""))
redfile$color= rep("red",dim(redfile)[1])
redfile=redfile[,-c(1,2)]
orangefile$color= rep("orange",dim(orangefile)[1])
orangefile=orangefile[,-c(1,2)]
greenfile$color= rep("green",dim(greenfile)[1])
greenfile=greenfile[,-c(1,2)]
bluefile$color= rep("blue",dim(bluefile)[1])
bluefile=bluefile[,-c(1,2)]
purplefile$color= rep("purple",dim(purplefile)[1])
purplefile=purplefile[,-c(1,2)]
consolidated=rbind(redfile,greenfile,bluefile,orangefile,purplefile)
write.xlsx(consolidated, paste("ClustersHeatmaps/BinaryResults/","consolidated","-SigRules.xlsx",sep=""))
consolidated=read.xlsx(paste("ClustersHeatmaps/BinaryResults/","consolidatedUpdated","-SigRules.xlsx",sep=""),sheetIndex = 1)


read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x=ldply(x, rbind)
  return(x)
}

createDominanceTable<-function(consolidated,colors,indicesPhenotypes,indicesPhenotypesCategorical)
{
  sigVarFrequnecy=matrix(0,nrow=length(colors),ncol=length(append(names(indicesPhenotypes),names(indicesPhenotypesCategorical))))
  colnames( sigVarFrequnecy)=append(names(indicesPhenotypes),names(indicesPhenotypesCategorical))
  rownames( sigVarFrequnecy)=colors
  sigVarFrequnecy=as.data.frame( sigVarFrequnecy)
  View(sigVarFrequnecy)
  for(i in 1:length(colors))
  {
   
    
      
      temp= consolidated[which(consolidated$color==colors[i]),]
      View(temp)
      temp=as.data.frame(temp)
      print(colnames(temp))
      temp1=unique(temp$phenotype)
     #for(j in 1:length(unique(temp$phenotype)))
        
      for(j in 1:length(temp1))
        
      { 
     #   temp1=unique(temp$phenotype)
        print(temp1[j])
        print(temp$phenotype)
        print(dim(temp[which(temp$phenotype==temp1[j]),])[1]) 
       # sigVarFrequnecy[i,which(colnames( sigVarFrequnecy)==as.character(temp1[j]))]=dim(temp[which(temp$phenotype==temp1[j]),])[1]
       print(temp1[j])
           print( grep(temp1[j],colnames( sigVarFrequnecy),ignore.case=TRUE,value=FALSE))
         sigVarFrequnecy[i,grep(paste("^",temp1[j],"$", sep=""),colnames( sigVarFrequnecy),ignore.case=TRUE,value=FALSE)]=dim(temp[which(temp$phenotype==temp1[j]),])[1]
     #   sigVarFrequnecy[i,j]=dim(temp[which(temp$phenotype==temp1[j]),])[1]
        
             }
      
 
  }

  return(sigVarFrequnecy)
  
}
#------------------------Values of phenotypes for each rule---------------------------------------------

returnRulesPhenoValues<-function(clusteredRules,filt,samples,phenotypeTable,phenotype)
  
{
# write.xlsx(clusteredRules,paste("ClustersHeatmaps/","RulesPhenoValues.xlsx",sep=""),  sheetName="SamplesSupportedPerRule")
  #write.csv(clusterRules,paste("ClustersHeatmaps/","SamplesSupportedPerRule.csv",sep=""))
  
  # namesPheno=names(indicesPhenotypes)
  # 
  # 
  #   sigPhenoRules<- filt
  #   rules=rep(1:dim(filt)[1])
  #   
  #   
  #  # for(j in 1:length(namesPheno))
  #     for(j in 1:1)
  #   {
  #     print(paste("Rule",rules[j]))
  #     # print(j)
  #     # 
  #     
  #     temp=sigPhenoRules[j,]
  #     # print(temp)
  #     # 
  #     
  #     phenotypes=temp$phenotype
  #     supportSet=unlist(as.list(strsplit(as.character(filt$supportSetRHS[[j]]), ",")))
  #     # print("SupportSet")
  #     #  print(length(supportSet))
  #     # print(supportSet)
  #     
  #     supported=supportSet
  #     # print("Supported")
  #     #  print(length(supported))
  #     # print(supported)
  #     notSupported=supportSet[which(!(supportSet %in% samples))]
  #     
  #     dependant=append(rep(1,(length(supported))),rep(0,(length(notSupported))))
  #     
  #     #model=matrix(NA,nrow = length(dependant),ncol=length(indicesPhenotypes))
  #     model=matrix(NA,nrow = length(dependant),ncol=2)
  #     
  #     #colnames(model)=names(indicesPhenotypes)
  #     colnames(model)=c(namesPheno[j],"decision")
  #     rownames(model)=append(supported,notSupported)
  #     # print(rownames(model))
  #     # print(colnames(model))
  #     model[,1]=phenotypeTable[rownames(model),namesPheno[j]]
  #     model[,2]=dependant
  #     #model=phenotypeTable[rownames(model),names(indicesPhenotypes)]
  #     
  #     model=as.data.frame(model)
  #     #View(model)
  #     model=apply(model,2,as.numeric)
  #     
  #     model=as.data.frame(model)
  #     # View(model)
  #    # model$decision=as.factor(dependant)
  #     #Removing random effects
  #     # model$subject=as.factor(c(1:length(rownames(model))))
  #     #model$decision <- relevel(model$decision, "0")
  #     model=as.data.frame(na.omit(model))
  #     # print(class(model))
  #     View(model)
  #     
  #     #model <- as.data.frame(model[!apply(model, 1, function(x) {any(x == NA)}),])
  #     # print("I am here")
  #     
  #     # View(model)
  #    }
  
  
  # heatmap4=matrix(0,nrow=dim(filt)[1],ncol=length(samples))
  # 
  # # View(heatmap4)
  # heatpmap4=as.data.frame(heatmap4)
  # #heatmap4=apply(heatmap4,2,as.character)
  # # rownames(heatmap4)<-indicesRulesSig
  # rownames(heatmap4)<-rep(1:dim(filt)[1])
  # colnames(heatmap4)<-samples
  # 

  #for(j in 1:length(indicesRulesSig))
  
  
  outputSigRules=NULL;
  
  # print("This is indicesRules")
  #print(indicesRules)
  final=NULL
  temp=NULL
  rulesNumber=NULL
  for(i in 1:length(phenotype))
{
    temp=append(temp,as.character(rep(as.character(phenotype[i]),dim(filt)[1])))
    rulesNumber=append(rulesNumber,rep(1:dim(filt)[1]))
      
    
    heatmap4=matrix(0,nrow=dim(filt)[1],ncol=length(samples))
    
    # View(heatmap4)
    heatpmap4=as.data.frame(heatmap4)
    #heatmap4=apply(heatmap4,2,as.character)
    # rownames(heatmap4)<-indicesRulesSig
    rownames(heatmap4)<-rep(1:dim(filt)[1])
    colnames(heatmap4)<-samples
    
  for(j in 1:dim(filt)[1])
  {
    
    
    # print("SIGNIFICANT RULE NUMBER ")
    ruleNumber=j
    # ruleNumber=as.numeric(indicesRulesSig[j])
    #  print("this is the ruleNumber")
    print(ruleNumber)
    supportSet=unlist(as.list(strsplit(as.character(filt$supportSetRHS[ruleNumber]), ",")))
    
     print(supportSet)
    supported=supportSet
    # print(supported)
    supportedValues=phenotypeTable[which(rownames(phenotypeTable) %in% supported),phenotype[i]]
    print(supportedValues)
    
    #Check here
    
    #For Debugging here
    #  supportedupdated= supported[which(supported %in% SupportRowNames)]
    ## print(SupportRowNames)
    ##print(length(SupportRowNames))
    
    #   View(heatmap4)
    #print(SupportRowNames %in% colnames(heatmap4))
    #print(length(heatmap4[j,supportedupdated]))
    ##  print(length(supportedValues))
    # heatmap4[j,supportedupdated]=supportedValues
    
    print(dim(heatmap4))
    heatmap4[j,supported]=supportedValues
  
    
    
    
  }
    
     
  View(heatmap4)
x=as.data.frame(heatmap4)
View(x)
print(temp)
#heatmap4$phenotype=as.character(temp)

View(heatmap4)
#temp=phenotype[i]
if(i==1)
  final=as.data.frame(heatmap4)
else
  
final=as.data.frame(rbind(final,heatmap4))

#View(final)
   # write.xls(x, file=paste("ClustersHeatmaps/","RulesPhenoValues.xls",sep=""), sh.names =as.character(temp))
  
  
  }
View(temp)
final$phenotype=temp
final$ruleNumber=rulesNumber
final=as.data.frame(final)

write.csv(final, file=paste("ClustersHeatmaps/","RulesPhenoValues.csv",sep=""))

return(final)
}
final=returnRulesPhenoValues(clusters,filt,rownames(dt_40_remove1),phenotype,append(names(indicesPhenotypes),names(indicesPhenotypesCategorical)))
redValuesPheno=final[,rownames(clusters[which(clusters$color_cluster=="red"),])]
greenValuesPheno=final[,rownames(clusters[which(clusters$color_cluster=="green"),])]
blueValuesPheno=final[,rownames(clusters[which(clusters$color_cluster=="blue"),])]
orangeValuesPheno=final[,rownames(clusters[which(clusters$color_cluster=="orange"),])]
purpleValuesPheno=final[,rownames(clusters[which(clusters$color_cluster=="purple"),])]




write.xlsx(clusterRules,paste("ClustersHeatmaps/","RulesPhenoValues.xlsx",sep=""),  sheetName="SamplesSupportedPerRule",append = TRUE)
write.xlsx(final2,paste("ClustersHeatmaps/","RulesPhenoValues.xlsx",sep="") ,sheetName="RulesPhenoValues",append = TRUE)
write.xlsx(redValuesPheno,paste("ClustersHeatmaps/","RulesPhenoValues.xlsx",sep=""),  sheetName="RedSamplesSupportedPerRule",append = TRUE)
write.xlsx(greenValuesPheno,paste("ClustersHeatmaps/","RulesPhenoValues.xlsx",sep=""),  sheetName="GreenRedSamplesSupportedPerRule",append = TRUE)
write.xlsx(blueValuesPheno,paste("ClustersHeatmaps/","RulesPhenoValues.xlsx",sep=""),  sheetName="BlueSamplesSupportedPerRule",append = TRUE)
write.xlsx(purpleValuesPheno,paste("ClustersHeatmaps/","RulesPhenoValues.xlsx",sep=""),  sheetName="PurpleSamplesSupportedPerRule",append = TRUE)
write.xlsx(orangeValuesPheno,paste("ClustersHeatmaps/","RulesPhenoValues.xlsx",sep=""),  sheetName="OrangeSamplesSupportedPerRule",append = TRUE)


#----------------------------------------------------------------------------------------

final2=apply(final,2,as.numeric)
final2[final2==0] <- NA
final2=as.data.frame(apply(final2,2,as.character))
redValuesPheno2=returnWithMeansSD(redValuesPheno)
redValuesPheno2=as.data.frame(apply(redValuesPheno2,2,as.character))
blueValuesPheno2=returnWithMeansSD(blueValuesPheno)
orangeValuesPheno2=returnWithMeansSD(orangeValuesPheno)
greenValuesPheno2=returnWithMeansSD(greenValuesPheno)
purpleValuesPheno2=returnWithMeansSD(purpleValuesPheno)

returnWithMeansSD<-function(matrixvalues)
{
  
  matrixvalues=apply(matrixvalues,2,as.numeric)
  matrixvalues[matrixvalues==0] <- NA
meansFinal=rowMeans(matrixvalues, na.rm = TRUE)
sdFinal=rowSds(matrixvalues,na.rm=TRUE)
matrixvalues=as.data.frame(matrixvalues)
matrixvalues$mean=meansFinal
matrixvalues$SD=sdFinal
matrixvalues=as.data.frame(apply(matrixvalues,2,as.character))
return(as.data.frame(matrixvalues))

}




write.xlsx(clusterRules,paste("ClustersHeatmaps/","RulesPhenoValuesNAs.xlsx",sep=""),  sheetName="SamplesSupportedPerRule",append = TRUE)
write.xlsx(final2,paste("ClustersHeatmaps/","RulesPhenoValuesNAs.xlsx",sep="") ,sheetName="RulesPhenoValues",append = TRUE)
write.xlsx(redValuesPheno2,paste("ClustersHeatmaps/","RulesPhenoValuesNAs.xlsx",sep=""),  sheetName="Red",append = TRUE)
write.xlsx(greenValuesPheno2,paste("ClustersHeatmaps/","RulesPhenoValuesNAs.xlsx",sep=""),  sheetName="Green",append = TRUE)
write.xlsx(blueValuesPheno2,paste("ClustersHeatmaps/","RulesPhenoValuesNAs.xlsx",sep=""),  sheetName="Blue",append = TRUE)
write.xlsx(purpleValuesPheno2,paste("ClustersHeatmaps/","RulesPhenoValuesNAs.xlsx",sep=""),  sheetName="Purple",append = TRUE)
write.xlsx(orangeValuesPheno2,paste("ClustersHeatmaps/","RulesPhenoValuesNAs.xlsx",sep=""),  sheetName="Orange",append = TRUE)





 