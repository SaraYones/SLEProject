library(ggbiplot)
library(arules)
library('Hmisc')
library("ggpubr")
library(sva)
library(limma)
library(plyr)
library("caret")
library('doParallel')
library('data.table')
cl <- makeCluster(4)
registerDoParallel(cl)
require("VennDiagram")
library(devtools)
require(ggplot2)

install_github("GabrielHoffman/variancePartition")
library(variancePartition)
library('variancePartition')

preprocessMetaData=function(metadata)
{
  metadata[,"visit:"]=gsub("(visit: (([0-9]*[.])?[0-9]+))","\\2",metadata[,"visit:"])
  metadata[,"visit_count:"]=gsub("(visit_count: (([0-9]*[.])?[0-9]+))","\\2",metadata[,"visit_count:"])
  metadata[,"cumulative_time:"]=gsub("(cumulative_time: (([0-9]*[.])?[0-9]+))","\\2",metadata[,"cumulative_time:"])
  metadata[,"days_since_diagnosis:"]=gsub("(days_since_diagnosis: (([0-9]*[.])?[0-9]+))","\\2",metadata[,"days_since_diagnosis:"])
  metadata[,"days_since_last_visit:"]=gsub("(days_since_last_visit: (([0-9]*[.])?[0-9]+))","\\2",metadata[,"days_since_last_visit:"])
  metadata[,"days_between_diagnosis_and_last_visit:"]=gsub("(days_between_diagnosis_and_last_visit: (([0-9]*[.])?[0-9]+))","\\2",metadata[,"days_between_diagnosis_and_last_visit:"])
  metadata[,"gender:"]=gsub("(gender: ([.])*)","\\2",metadata[,"gender:"])
  metadata[,"race:"]=gsub("(race: ([.])*)","\\2",metadata[,"race:"])
  metadata[,"age:"]=gsub("(age: (([0-9]*[.])?[0-9]+))","\\2",metadata[,"age:"])
  metadata[,"biopsy_history:"]=gsub("(biopsy_history: ([.])*)","\\2",metadata[,"biopsy_history:"])
  
  #rename(metadata, c("visit:"="visit", "visit_count:"="visit_count","cumulative_time:"="cumulative_time","days_since_diagnosis:"="days_since_diagnosis","days_since_last_visit:"="days_since_last_visit","days_between_diagnosis_and_last_visit:"="days_between_diagnosis_and_last_visit",
   #                  "gender:"="gender","race:"="race","age:"="age","biopsy_history:"="biopsy_history"))
  names(metadata)[names(metadata)=="visit:"]="visit"
  names(metadata)[names(metadata)=="visit_count:"]="visitCount"
  names(metadata)[names(metadata)=="cumulative_time:"]="cumulativeTime"
  names(metadata)[names(metadata)=="days_since_diagnosis:"]="daysSinceDiagnosis"
  names(metadata)[names(metadata)=="days_since_last_visit:"]="daysSinceLastVisit"
  names(metadata)[names(metadata)=="days_between_diagnosis_and_last_visit:"]="daysBetweenDiagnosisAndLastVisit"
  names(metadata)[names(metadata)=="gender:"]="gender"
  names(metadata)[names(metadata)=="race:"]="race"
  names(metadata)[names(metadata)=="age:"]="age"
  names(metadata)[names(metadata)=="biopsy_history:"]="biopsyHistory"
  names(metadata)[names(metadata)=="subject:"]="subject"
  
  metadata[,"visit"]<-as.numeric(metadata[,"visit"])
  metadata[,"visitCount"]=as.numeric(metadata[,"visitCount"])
  metadata[,"cumulativeTime"]=as.numeric(metadata[,"cumulativeTime"])
  metadata[,"daysSinceDiagnosis"]=as.numeric(metadata[,"daysSinceDiagnosis"])
  metadata[,"daysSinceLastVisit"]=as.numeric(metadata[,"daysSinceLastVisit"])
  metadata[,"daysBetweenDiagnosisAndLastVisit"]= as.numeric(metadata[,"daysBetweenDiagnosisAndLastVisit"])
  metadata[,"gender"]=as.character(metadata[,"gender"])
  metadata[,"race"]= as.character(metadata[,"race"])
  metadata[,"age"]=as.numeric(metadata[,"age"])
  metadata[,"biopsyHistory:"]=as.character(metadata[,"biopsyHistory"])
  

  
  return(metadata)
}
#rcorr(decisionSLE12,metadata12[,"gender:"])

#shapiro.test()
#ggscatter(as.data.frame(cbind(trial,decisionSLE12)), x = "trial", y = "decisionSLE12", 
#         add = "reg.line", conf.int = TRUE, 
#        cor.coef = TRUE, cor.method = "spearman",
#       xlab = "gender", ylab = "decision")



calculateCorrelation=function(datatype,x,y)
{
  if(datatype=="categorical")
  {
    chi2 = chisq.test(table(cbind(trial,decisionSLE12)), correct=F)
    c(chi2$statistic, chi2$p.value)
    sqrt(chi2$statistic / sum(table(cbind(trial,decisionSLE12))))
  }else{
    
    
  }
  spine(as.factor(trial), as.factor(decisionSLE12))
}



plotPCAmeta=function(GE,metadata,variable,filepath,discretizemethod)
{
  
  
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    
    
    temp=metadata[[i]]
    temp=temp[,variable]
    
    if(variable=="age:"|variable=="Age.at.Diagnosis.in.Days")
    {
      temp=discretizeAge(temp,discretizemethod)
      
      
    }
    
    
    
    #print(GE[i])
    #print(temp[,variable])
    my.plots[[i]]=plotPCA(filepath,GE[[i]],temp,variable)
    
  }
  
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}

savePDF=function(myplots,variable,filepath,discretizemethod="")
{
  
  
  
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}

discretizeAge=function(values,method)
{
  
  #if using SLE
 # age=as.numeric(gsub("(age: (([0-9]*[.])?[0-9]+))","\\2",values))
  #if using AML data
  age=values
  print(age)
  discretizedAgeFrequency=discretize(age,method = method, breaks = 4,  labels = c("young", "youth","adult","senior"))
  return(discretizedAgeFrequency)
}


removeBatchEffect= function(logGeneExpressionMatrixBatch,pheno,requiredPheno,flag){
  
  #Ask Mateusz how to deal with this because I dont want an exact name for the factor i want to substitute the value of requiredPheno
  #pheno2=as.data.frame(pheno) 
  
  #Creating a model for sva 
  # mod = model.matrix(~cyto_risk_group, data=pheno)
  # mod0 = model.matrix(~1,as.data.frame(pheno))
  #   n.sv = num.sv(t(logGeneExpressionMatrixBatch),mod,method="leek")
  #  svobj = sva(t(logGeneExpressionMatrixBatch),mod,mod0,n.sv=n.sv)
  #   modSv = cbind(mod,svobj$sv)
  #  mod0Sv = cbind(mod0,svobj$sv)
  #fit = lmFit(t(logGeneExpressionMatrixBatch),modSv)
  #Batches known
  if(flag==1)
  #Creating a model for COMBAT with known batches
  {
  modcombat = model.matrix(~1, data=as.data.frame(seq(1,length(pheno))))
  print(modcombat)
  combat_edata = ComBat(dat=as.matrix(t(logGeneExpressionMatrixBatch)), batch=pheno, mod=modcombat)
  }

  
  
  return(combat_edata)
}

pheno = pData(bladderEset)


plotPCA= function(filepath,GeneExpressionMatrixlocal,Groups,variable){
  
  myplots=NULL
  #pdf(filepath)
  # log transform 
  # apply PCA - scale. = TRUE is highly 
  # advisable, but default is FALSE. 
  
  #ir.pca <- prcomp(GeneExpressionMatrixlocal,
  #                center = TRUE,
  #               scale. = TRUE) 
  ir.pca <- prcomp(GeneExpressionMatrixlocal,
                   center = TRUE) 
  
  
  g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
                groups =Groups, ellipse = FALSE, 
                circle = FALSE,var.axes = FALSE)
  g <- g + scale_color_discrete(name = variable)
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  #print(g)
  #dev.off()
  plot(g)
  myplots=recordPlot()
  return(myplots)
}

removeColsWithZeroVar=function(GeneExpression)
{
  for(i in 1:length(GeneExpression))
  {
    remove_cols=nearZeroVar(as.data.frame(GeneExpression[i]))
    temp=as.data.frame(GeneExpression[i])
   
    if(length(remove_cols)!=0)
    {
      GeneExpression[[i]]=as.data.frame(temp[,-remove_cols])
    }else{
      GeneExpression[[i]]=as.data.frame(temp)
      
    }
   
  }
return(GeneExpression)
  }

checkVariableEffects=function(GeneExpr,metadata,form,filepath,variable,discretizemethod="")
{ my.plots <- vector(length(GeneExpr), mode='list')
  myplots=NULL
  for(i in 1:length(GeneExpr)){
  varPart <- fitExtractVarPartModel( t(as.data.frame(GeneExpr[[i]])), form, as.data.frame(metadata[[i]]))
  vp <- sortCols(varPart)
  print(plotVarPart(vp))
  myplots=recordPlot()
  my.plots[[i]]=myplots
  }
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots){
    replayPlot(my.plot)
  }
  graphics.off()

}

coorelationVariables=function(metadata,decisionSLE,form,filepath,variable,discretizemethod="")
  
{

  
  my.plots <- vector(length(metadata), mode='list')
  myplots=NULL
  for(i in 1:length(metadata)){
    print("hello")
    metadata[[i]]=as.data.frame(cbind(metadata[[i]],decisionSLE[[i]]))
    names(metadata[[i]])[names(metadata[[i]])=="decisionSLE[[i]]"]="decisionSLE"
    View(metadata[[i]])
    C = canCorPairs( form, metadata[[i]])
    print( plotCorrMatrix( C ))
    myplots=recordPlot()
    my.plots[[i]]=myplots
  }
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots){
    replayPlot(my.plot)
  }
  graphics.off()
}

exploreSamples=function(GeneExpression,titles,form,filepath,variable,discretizemethod="")
{
  
  par(mar=c(7,5,1,1))
  my.plots <- vector(length(GeneExpression), mode='list')
  myplots=NULL
  for(i in 1:length(GeneExpression)){
    #print("hello")
    
    sample <- rownames(as.data.frame(GeneExpression[[i]]))
    d.f <- data.frame(sample,as.data.frame(GeneExpression[[i]]))
    d.f2 <- melt(d.f, id.vars = "sample")
    
    
    
    boxplot(value~sample,
            data=d.f2,
            main=titles[[i]],
            xlab="samples",
            ylab="gene expr",
            col="orange",
            border="brown",las=2)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
  }
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots){
    replayPlot(my.plot)
  }
  graphics.off()
  
  
}

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
  features=features[seq(1,numberOfFeatures),"attribute"]
  return(features$attribute)
  
}
compareAccuracies=function(Decisiontable,limit,features)
{
  Accuracies=NULL
  
  for(i in seq(10, limit, by = 10)){
    print(i)
      
     # resultRosettaWithUSMod7=rosetta(DecisiontableBinary,classifier="StandardVoter",discrete=TRUE,underSample = TRUE)
     #For SLE 
    #resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discreteMask=TRUE,discrete = TRUE,ruleFiltration=TRUE,ruleFiltrSupport=c(1,3))
     #FOR AML
    resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3)
    
     Accuracies=append(Accuracies,resultRosetta$quality$Accuracy.Mean)
  }
    plot(Accuracies,type='l',xaxt="n",main="Accuracies")
    axis(1,at=seq(1,as.numeric(limit/10), by = 1),labels=as.character(seq(10, limit, by = 10)))
    return(Accuracies)
    
  }
  
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
    SUPP_SET_LHS=unlist(as.list(strsplit(as.character(result$SUPP_SET_LHS[[i]]), ",")))
    #SUPP_SET_RHS=unlist(as.list(strsplit(as.character(result$SUPP_SET_RHS[[i]]), ",")))
    #SUPP_SET=intersect(SUPP_SET_LHS,SUPP_SET_RHS)
    resultMatrix[i,which(colnames(resultMatrix) %in% SUPP_SET_LHS)]=1
    
  }
  resultMatrix=resultMatrix[,which(colSums(resultMatrix)!=0)]
  #resultMatrix=resultMatrix[,which(colSums(resultMatrix)!=-dim(resultMatrix)[1])]
  return(resultMatrix)
}

#-----------------Plotting Functions--------------------------------
plotVenn=function(geneLS,parameters)
{
  
  VENN.LIST <- geneLS
  venn.plot <- venn.diagram(VENN.LIST , NULL, fill=parameters[[1]], alpha=parameters[[2]], cex = as.numeric(parameters[[3]]), cat.fontface=4, category.names=parameters[[5]], main="Gene Lists")
  
  # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
  grid.draw(venn.plot)
  
  # To get the list of gene present in each Venn compartment we can use the gplots package
  require("gplots")
  
  a <- venn(VENN.LIST, show.plot=FALSE)
  
  # You can inspect the contents of this object with the str() function
  inters <- attr(a,"intersections")
  return(inters)
  
  
}

