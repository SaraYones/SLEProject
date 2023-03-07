library(ggbiplot)
library(arules)
library('Hmisc')
library("ggpubr")
library(sva)
library(limma)
library(plyr)
library("caret")
library(ggplot2)
require(reshape2)
library('doParallel')
cl <- makeCluster(4)
registerDoParallel(cl)
library(Boruta)
library(R.ROSETTA)
library(plyr)
library(gplots) 
#install.packages("RVenn")
library(purrr)
library(Rvenn)
library(corrplot)
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("Ringo", version = "3.8")
library(Ringo)
options(stringsAsFactors = FALSE)
options(numericAsFactors=FALSE)

#Prerequiste you have to run first the correct GE script

#For running test to remove multicolleniarity in a logistic regression
library(performance)

#For oversampling 
#install.packages("ROSE")
library(ROSE)


#install_github("kassambara/easyGgplot2")

# Load easyGgplot2
library(easyGgplot2)


#for permutation test
library(permute)
library(stringr)

#--------------------------------Discretizations-----------------------------------------
#Here we want to extract 1 healthy individual per each technical replicate
metadata=preprocessMetaData(metadata)
#make rownames to a coloumn using dyplr
metadata=tibble::rownames_to_column(metadata)
#aggregate according to rowname to see what are the subjects that share the same rowname
aggMetadata=aggregate(subject ~ rowname, as.matrix(metadata),unlist)
#Filter only healthy
aggMetadata=aggMetadata[which(grepl("subject: BAY-([.])*",aggMetadata[,"subject"])),]
#check physically what subjects share the same rowname
rows=apply(aggMetadata,1,function(x) which(aggMetadata[,"subject"]==x[2]) )
#Pick up a subject from each replicate
rows=unique(aggMetadata[unlist(lapply(rows,function(x) x[[1]])),"rowname"])
#print(ggplot(d.f2, aes(x=as.factor(sample), y=value)) + geom_boxplot(fill="slateblue", alpha=0.2)+xlab("samples"))
#make sure that all the picked are healthy individuals "They have to be all starting with "BAY"
metadata[which(metadata[,"rowname"] %in% rows),"subject"]
#get Gene expressions of healthy 
controls=logDAsPatientWithoutBatch[which(rownames(logDAsPatientWithoutBatch) %in% rows),]

#Calculate Means and standard deviation of controls
means=apply(controls,2,function(x) mean(x))
sd=apply(controls,2,function(x) sd(x) )

#compare to mean and sd for each gene
trial=mapply(function(x,y,z) sapply(x,function(x) if(x>y+2*z) x="up" else if (x<y-2*z) x="down" else x="normal")  ,logDAsPatientWithoutBatch,means,sd)

#Check that it actually works the correct way
trial1=mapply(function(x,y,z) if(x>y+2*z) x="up" else if (x<y-2*z) x="down" else x="normal",logDAsPatientWithoutBatch[,1],means[1],sd[1])

DescretizedDF<-as.data.frame(unlist(trial), nrow = dim(logDAsPatientWithoutBatch)[1], ncol = dim(logDAsPatientWithoutBatch)[2])
#If you want fast replacement transform into matrix instead of DF
#trialDF<-as.matrix(trialDF)

#trialDF[trialDF==TRUE]="High"
DescretizedDF$decisionSLE=decisionSLE
#---------------------------------Feature Selection using Boruta-------------------------
#I ran it on ULAM
boruta.train13=Boruta(decision~., data=as.data.frame(logDA13PatientWithoutBatch),holdHistory=FALSE)

borutaFeatures13=extractFeaturesBoruta(boruta.train13)
MCFSFeatures13=FilterFeatures("output13/output13__RI.csv",200)
inters=plotVenn(list(MCFSGenes=MCFSFeatures13,BorutaGenes=borutaFeatures13),list(c("darkmagenta", "darkblue"),c(0.5,0.5),2,2, c("MCFSGenes","BorutaGenes")))
Features=MCFSFeatures13[which(MCFSFeatures13 %in% inters$`MCFSGenes:BorutaGenes`)]
DescretizedDF13=DescretizedDF[which(decisionSLE==1 | decisionSLE==3),]
DescretizedDF13$decisionSLE=as.character(DescretizedDF13$decisionSLE)
#resultRosetta13=rosetta(DescretizedDF13[,append(borutaFeatures13,"decisionSLE")],classifier="StandardVoter",)
rownames(DescretizedDF13)<-rownames(logDA13PatientWithoutBatch)

tempfoundation=prepareDT(DescretizedDF13[,append(borutaFeatures13[1:40],"decisionSLE")],c(1,2,3))

SLEclassifier=fread("DA13.csv") 
SLEclassifier=as.data.frame(SLEclassifier)
row.names(SLEclassifier) <- SLEclassifier$V1
SLEclassifier=SLEclassifier[,2:dim(SLEclassifier)[2]]


SLE_Johnson<- Classifier(classifier = dt_40_remove1,flagAccuracy="Johnson",path="/newRun",
                                          MCFSFeatures=FilterFeatures("/Users/hp/Desktop/out_remove1/output13_remove1_RI.csv",1000)[[1]],
                                          
                         ontology="BP",numberOfFeatures=34,keyType="SYMBOL",underSample=FALSE)

SLE_Johnson$computeEnrichment()
resultRosetta13=rosetta(temp,classifier="StandardVoter",discrete = TRUE)
>>>>>>> aa77e9918b51e4c48910f07bae5bfd662841554a

SLE_Johnson<- Classifier(classifier = discDA13_remove1,flagAccuracy="Johnson",path="/newRun",
                         MCFSFeatures=FilterFeatures("/Users/hp/Desktop/out_remove1/output13_remove1_RI.csv",1000),
                         
                         ontology="BP",numberOfFeatures=40,keyType="SYMBOL",underSample=FALSE)


#For computing the Feature boosting step
temporaryTable <- data.frame(lapply(discDA13_remove1,as.character), stringsAsFactors=FALSE)

SLE_Johnson<- Classifier(classifier = temporaryTable,flagAccuracy="Johnson",path="/newRun",
                         MCFSFeatures=FilterFeatures("out_remove1/output13_remove1_RI.csv",1000),
                         
                         ontology="BP",numberOfFeatures=40,keyType="SYMBOL",underSample=FALSE)

x=as.character(seq(10, 500, by = 10))


figure=cbind(x,SLE_Johnson$Accuracies)
figure=as.data.frame(figure)
colnames(figure)=c("x","Accuracies")
figure$x=as.numeric(x)
figure$Accuracies=as.numeric(SLE_Johnson$Accuracies)

ggplot(figure, aes(x=x, y = Accuracies)) + geom_line()+theme_classic() +xlab("Features count")+ylab("Accuracies")+geom_vline(xintercept=c(50,60), linetype="dotted")




SLE_Johnson$computeEnrichment()
#resultRosetta13=rosetta(tempfoundation,classifier="StandardVoter",discrete = TRUE)
resultRosetta13_foundation=rosetta(tempfoundation,classifier="StandardVoter",discrete = TRUE)
#temp=as.matrix(temp)
#temp=apply(as.matrix(temp),1,function(x) lapply(x,function(y) print(y)))
#temp$decisionSLE=as.factor(temp$decisionSLE)

#recalculatedResultRosetta13=recalculateRules(temp,resultRosetta13$main,discrete = TRUE)
recalculatedResultRosetta13_foundation=recalculateRules(tempfoundation,resultRosetta13_foundation$main,discrete = TRUE)
#Filter according to pval of the rules

filt_foundation<- resultRosetta13_foundation$main[which(resultRosetta13_foundation$main$pValue <= 0.05),]
vis_out_foundation<- visunet(filt, type = "RDF")
#recalculatedResultRosetta13_foundation=recalculatedResultRosetta13_foundation[which(recalculatedResultRosetta13$PVAL<=0.08),]

ruleHeatmap(temp, recalculatedResultRosetta13, discrete =TRUE,ind=1)

saveLineByLine(resultRosetta13$main, "resultRosetta13.txt", discrete=TRUE, filterByPval=FALSE, pval=0.01)
#Get the best number of features to use
Accuracies=compareAccuracies(DescretizedDF13[,append(borutaFeatures13,"decisionSLE")],120,borutaFeatures13)
resultRosetta13Genetic=rosetta(DescretizedDF13[,append(borutaFeatures13,"decisionSLE")],classifier="StandardVoter",reducer="Genetic",ruleFiltration=TRUE,ruleFiltrSupport=c(1,7),discrete = TRUE)
#Didnt outpreform johnson in this case
resultRosetta13Genetic60=readRDS("resultRosetta13Genetic60")
resultRosetta13Genetic40=readRDS("resultRosetta13Genetic40")

filterResultRosetta13=recalculatedResultRosetta13[recalculatedResultRosetta13$DECISION=="3",]

clusteredRules=clusterRules(filterResultRosetta13,rownames(DescretizedDF13))


clusteredRules=clusterRules(recalculatedResultRosetta13,rownames(DescretizedDF13))



#heatmap(clusteredRules, scale = "none", Rowv = NA, Colv = NA, col = cm.colors(2), main = "HeatMap Example") 
heatmap(clusteredRules, scale = "none",  col = cm.colors(2), main = "HeatMap Example")


clusters=heatmap.F(t(clusteredRules),distmethod='pearson')
f.FeatureHeatmap(resultRosetta13$main)
