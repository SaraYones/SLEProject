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



#load('/Users/sarayones/Desktop/Projects/SLEProject/SLEWorkspace1.RData')
#Designing data structures------------------------------------
metadata=readRDS('/Users/sarayones/Desktop/Projects/SLEProject/data/metadata.rds')
DAsPatients=readRDS('/Users/sarayones/Desktop/Projects/SLEProject/data/DAsPatients.rds')

SampleIDs=DAsPatients[,1]
rownames(DAsPatients)<-SampleIDs
DAsPatients=DAsPatients[order(rownames(DAsPatients)), ] 
decisionSLE=DAsPatients$decision
DAsPatients=DAsPatients[,3:length(DAsPatients)-1]
logDAsPatients=DAsPatients+0.00001
logDAsPatients=log2(logDAsPatients)
logDAsPatients=logDAsPatients[order(rownames(logDAsPatients)), ] 
metadata=metadata[order(rownames(metadata)), ] 

logDA12Patients=logDAsPatients[which(decisionSLE==1 | decisionSLE==2),]
decisionSLE12=decisionSLE[which(decisionSLE==1 | decisionSLE==2)]
metadata12=metadata[which(decisionSLE==1 | decisionSLE==2),]

logDA13Patients=logDAsPatients[which(decisionSLE==1 | decisionSLE==3),]
decisionSLE13=decisionSLE[which(decisionSLE==1 | decisionSLE==3)]
metadata13=metadata[which(decisionSLE==1 | decisionSLE==3),]

logDA23Patients=logDAsPatients[which(decisionSLE==2 | decisionSLE==3),]
decisionSLE23=decisionSLE[which(decisionSLE==2 | decisionSLE==3)]
metadata23=metadata[which(decisionSLE==2 | decisionSLE==3),]

exploreSamples(list(as.data.frame(logDA12Patients),as.data.frame(logDA13Patients),as.data.frame(logDA23Patients)),list("DA12","DA13","DA23"),form,filepath,"ExploreSamples")

#--------------------------Discretize decision--------------------------------------------------------------------------------------
table(discretize(as.numeric(metadata12[,"age:"]), breaks = 3))
#-----------------Plotting-------------------------------------
#----------Samples with batch-----------------------------------------------------------------------------------------------------
filepath=paste('/Users/saryo614/Desktop/Projects/SLEProject/plots/')
plotPCAmeta(list(as.data.frame(logDA12Patients),as.data.frame(logDA13Patients),as.data.frame(logDA23Patients)),list(as.data.frame(metadata12),as.data.frame(metadata13),as.data.frame(metadata23)),"batch:",filepath)
plotPCAmeta(list(as.data.frame(logDA12Patients),as.data.frame(logDA13Patients),as.data.frame(logDA23Patients)),list(as.data.frame(metadata12),as.data.frame(metadata13),as.data.frame(metadata23)),"gender:",filepath)
plotPCAmeta(list(as.data.frame(logDA12Patients),as.data.frame(logDA13Patients),as.data.frame(logDA23Patients)),list(as.data.frame(metadata12),as.data.frame(metadata13),as.data.frame(metadata23)),"age:",filepath,"cluster")
plotPCAmeta(list(as.data.frame(logDA12Patients),as.data.frame(logDA13Patients),as.data.frame(logDA23Patients)),list(as.data.frame(metadata12),as.data.frame(metadata13),as.data.frame(metadata23)),"age:",filepath,"frequency")
#With batch
plotPCAmeta(list(as.data.frame(logDAsPatients)),list(as.data.frame(metadata)),"batch:",filepath,"")
#Without  batch
plotPCAmeta(list(as.data.frame(logDAsPatientWithoutBatch)),list(as.data.frame(metadata)),"batch:",filepath,"")


#------------Samples After Removing batches---------------------------------------------------------------------------------------
logDA12PatientWithoutBatch=as.data.frame(cbind(as.data.frame(t(removeBatchEffect(logDA12Patients,as.numeric(gsub("(batch: (([0-9]*[.])?[0-9]+))","\\2",metadata12[,"batch:"])),NULL))),decisionSLE12))
logDA13PatientWithoutBatch=as.data.frame(cbind(as.data.frame(t(removeBatchEffect(logDA13Patients,as.numeric(gsub("(batch: (([0-9]*[.])?[0-9]+))","\\2",metadata13[,"batch:"])),NULL))),decisionSLE13))
logDA23PatientWithoutBatch=as.data.frame(cbind(as.data.frame(t(removeBatchEffect(logDA23Patients,as.numeric(gsub("(batch: (([0-9]*[.])?[0-9]+))","\\2",metadata23[,"batch:"])),NULL))),decisionSLE23))
logDAsPatientWithoutBatch=as.data.frame(t(removeBatchEffect(logDAsPatient,as.numeric(gsub("(batch: (([0-9]*[.])?[0-9]+))","\\2",metadata[,"batch:"])),NULL)))


temp=removeColsWithZeroVar(list(as.data.frame(logDAsPatientWithoutBatch)))

logDA12PatientWithoutBatch=temp[[1]][which(decisionSLE==1 | decisionSLE==2),] 
logDA13PatientWithoutBatch=temp[[1]][which(decisionSLE==1 | decisionSLE==3),]
logDA23PatientWithoutBatch=temp[[1]][which(decisionSLE==2 | decisionSLE==3),]

#temp=removeColsWithZeroVar(list(as.data.frame(logDA12PatientWithoutBatch),as.data.frame(logDA13PatientWithoutBatch),as.data.frame(logDA23PatientWithoutBatch)))
#logDA12PatientWithoutBatch=temp[[1]]
#logDA13PatientWithoutBatch=temp[[2]]
#logDA23PatientWithoutBatch=temp[[3]]



logDA12PatientWithoutBatch=as.data.frame(cbind(logDA12PatientWithoutBatch,decisionSLE12))
logDA13PatientWithoutBatch=as.data.frame(cbind(logDA13PatientWithoutBatch,decisionSLE13))
logDA23PatientWithoutBatch=as.data.frame(cbind(logDA23PatientWithoutBatch,decisionSLE23))




#Plotting after Removing batches----------------------------------------------------------------
plotPCAmeta(list(as.data.frame(logDA12PatientWithoutBatch),as.data.frame(logDA13PatientWithoutBatch),as.data.frame(logDA23PatientWithoutBatch)),list(as.data.frame(metadata12),as.data.frame(metadata13),as.data.frame(metadata23)),"batch:",filepath,"WithoutBatch")
plotPCAmeta(list(as.data.frame(logDAsPatientWithoutBatch)),list(as.data.frame(metadata)),"batch:",filepath,"")

#----------------Find coorelation between variables------------------------------------------------
metadata12=preprocessMetaData(metadata12)
metadata13=preprocessMetaData(metadata13)
metadata23=preprocessMetaData(metadata23)
#------remove Zero Variance--------------------------------------------------------------------
temp=removeColsWithZeroVar(list(as.data.frame(logDA12Patients),as.data.frame(logDA13Patients),as.data.frame(logDA23Patients)))
logDA12Patients=temp[[1]]
logDA13Patients=temp[[2]]
logDA23Patients=temp[[3]]

#--------------------------------------------------------------------------------------------
temp=NULL
trial=sapply(metadata12[,"gender:"], function(x) if(x=="gender: F") append(temp,"F") else append(temp,"M"))

#--------------------------Correlation between variables-----------------------------
form <- ~ visit + visitCount + cumulativeTime + daysSinceDiagnosis + daysSinceLastVisit + daysBetweenDiagnosisAndLastVisit + gender + race + age + biopsyHistory + decisionSLE
coorelationVariables(list(as.data.frame(metadata12),as.data.frame(metadata13),as.data.frame(metadata23)),list(decisionSLE=decisionSLE12,decisionSLE=decisionSLE13,decisionSLE=decisionSLE23),form,filepath,"variableCorrelations")

#---------------------------Check the effect of different variables----------------------------------
form <- ~ visit + visitCount + cumulativeTime + daysSinceDiagnosis + daysSinceLastVisit + daysBetweenDiagnosisAndLastVisit + (1|gender) + (1|race) + age + (1|biopsyHistory)
checkVariableEffects(list(as.data.frame(logDA12Patients),as.data.frame(logDA13Patients),as.data.frame(logDA23Patients)),list(as.data.frame(metadata12),as.data.frame(metadata13),as.data.frame(metadata23)),form,filepath,"variableEffects")
#-----------------------------Save the tables ---------------------------------------------------------
write.table(logDA12PatientWithoutBatch, file ="DA12.csv",row.names = TRUE)
write.csv(logDA13PatientWithoutBatch, file ="DA13.csv")
write.table(logDA23PatientWithoutBatch, file ="DA23.csv",row.names = TRUE)

form <- ~ visit + visitCount + cumulativeTime + daysSinceDiagnosis + daysSinceLastVisit +(1|treatment)+ (1|gender) + (1|race) + age + (1|biopsyHistory)
checkVariableEffects(list(as.data.frame(logDA13PatientWithoutBatch)),list(as.data.frame(metadata13)),form,filepath,"variableEffects")



