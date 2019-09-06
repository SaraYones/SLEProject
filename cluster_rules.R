
install.packages("dynamicTreeCut")
install.packages("WGCNA")
install.packages("preprocessCore")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
BiocManager::install("preprocessCore")

{
  library(gplots)
  library(dynamicTreeCut)
  library(WGCNA)
  library(impute)
  library(preprocessCore)
}

#------------ DATA ----------------------------------------------------
result_John40_remove1 <- readRDS("/Users/hp/Desktop/mapp/result_John40_remove1.rds")
discDA13_remove1 <- read.csv("/Users/hp/Desktop/mapp/discDA13_remove1.csv")
features13_remove1 <- FilterFeatures("/Users/hp/Desktop/out_remove1/output13_remove1_RI.csv", 1000)
rownames(discDA13_remove1) <- discDA13_remove1$X
dt_40_remove1=prepareDT(discDA13_remove1[,append(features13_remove1[1:40],"decisionSLE")],c(1,2,3))


#--------- CLUSTER RULES ---------------------------------------------

recal_John40_remove1 <- recalculateRules(dt_40_remove1, 
                                         result_John40_remove1$main[result_John40_remove1$main$ACC_RHS >= 0.8 & result_John40_remove1$main$PVAL < 0.05,], 
                                         discrete = TRUE)

#cluster rules 
df <- recal_John40_remove1 #[recal_John40_remove1$DECISION == 1,] # & recal_John40_remove1$DECISION == 3,]
rowname <- rownames(dt_40_remove1) #[DA13_remove1$decisionSLE == 1,]) #& DA13_remove1$decisionSLE == 3,])

clusterRules <- clusterRule(df, 
                            rowname)

df <- dt_40_remove1
df <- df[rownames(df) %in% colnames(clusterRules),]
clusterRules <- clusterRules[, order(colnames(clusterRules))]

rownames(df) == colnames(clusterRules)

#build heat map 
rule_cluster=heatmap.F((clusterRules),
                       colors = c("white", "blue"), 
                       cutoffmethod = "number",
                       cutoff = 5,
                       distmethod='pearson')

rule_cluster <- as.data.frame(rule_cluster)
rule_clusterDA3 <- rule_clusterDA3[order(rownames(rule_clusterDA3)), , drop = FALSE]


#------------ FUNCTIONS ---------------------------

f.FeatureHeatmap = function(rule.table){
  fe.x = strsplit(rule.table[,"FEATURES"], split=",")
  con.x = strsplit(rule.table[,"CUTS_COND"], split=",")
  acc.x = rule.table[,"ACC_RHS"]
  supp.x = rule.table[,"SUPP_RHS"]
  
  fx = function(i){return(apply(FUN=paste, X=cbind(fe.x[[i]], con.x[[i]]), MARGIN=1, collapse="_"))}
  rule.list = lapply(FUN=fx, X=1:length(fe.x))
  feature.list = unique(unlist(rule.list))
  
  ti.x = c("W0D0","W0D1","W0D3","W0D7","W18D0","W18D1","W18D3","W18D7","W88D0")
  Mx = matrix(0, nrow=length(feature.list), ncol=length(ti.x))
  colnames(Mx) = ti.x
  for(i in 1:length(feature.list)){
    x1 = unlist(strsplit(feature.list[i], split="_"))
    Mx[i,][ti.x==strsplit(x1[2], split=".", fixed=T)[[1]][1]] = as.numeric(x1[3])
    Mx[i,][ti.x==strsplit(x1[2], split=".", fixed=T)[[1]][2]] = -1
  }
  Mx[Mx==1] = -2; Mx[Mx==2] = 1; Mx[Mx==3] = 2
  View(Mx)
  #heatmap.F(Mx, clusterdim="row", distmethod = "pearson", colors=c("blue","grey","white","red","red"))
}


heatmap.F = function(dataM, 
                     colors=c('darkblue', 'mediumblue', 'dodgerblue', 'white', 'orange', 'red', 'darkred'),
                     colorsaturation=0.25,
                     distmethod='bicor',
                     clustermethod='ward.D2',
                     clusterdim='both',
                     exclude=NULL,
                     cutoffmethod='depth',
                     cutoff=1,
                     main=NULL)
{
  if(length(colnames(dataM))==0){colnames(dataM) = as.character(1:ncol(dataM))}
  fx = function(v){return(!all(is.na(v)))}
  dataM = dataM[apply(FUN=fx, X=dataM, MARGIN=1),apply(FUN=fx, X=dataM, MARGIN=2)]
  distmethod = match.arg(arg=distmethod, choices=c('bicor', 'pearson', 'spearman', 'direct', 'euclidean'))
  clusterdim = match.arg(arg=clusterdim, choices=c('both', 'row', 'column', 'none'))
  cutoffmethod = match.arg(arg=cutoffmethod, choices=c('depth', 'height', 'number', 'none'))
  Rowv=T
  Colv=T
  dendrogram='both'
  if(clusterdim=='row'){Colv=F; dendrogram='row'}
  if(clusterdim=='column'){Rowv=F; dendrogram='column'}
  if(clusterdim=='none'){Rowv=F; Colv=F; dendrogram='none'}
  # Removing excluded samples:
  dataM = dataM[,which(colnames(dataM)%in%exclude==F)]
  # Color scale:
  neg.col = rev(colors[1:ceiling(length(colors)/2)])
  pos.col = colors[ceiling(length(colors)/2):length(colors)]
  bins = diff(round(quantile(abs(dataM), seq(from=0, to=1, length.out=length(pos.col))^colorsaturation, na.rm=T)/max(abs(range(dataM, na.rm=T)), na.rm=T), digits=3)*1000)
  a1 = list()
  a2 = list()
  for(i in 1:length(bins)){
    a1[[i]] = t(colorRamp(neg.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
    a2[[i]] = t(colorRamp(pos.col[c(i, i+1)])(seq(from=0, to=1, length.out=bins[i])))
  }
  a1 = matrix(unlist(a1), ncol=3, byrow=T); a2 = matrix(unlist(a2), ncol=3, byrow=T)
  
  b1 = rgb(red=a1[,1], green=a1[,2], blue=a1[,3], max=255)
  b2 = rgb(red=a2[,1], green=a2[,2], blue=a2[,3], max=255)
  c = abs(range(dataM))
  if(c[1]>c[2]){
    
    b2 = b2[1:round(length(b2)*c[2]/c[1])]
  } else {
    
    b1 = b1[1:round(length(b1)*c[1]/c[2])]
  }
  color.vector = c(rev(b1), b2)
  # Distance metric
  if(distmethod=='bicor'){fd = function(x){return(as.dist(1-bicor(t(x), use = 'pairwise.complete.obs')))}}
  if(distmethod=='pearson'){fd = function(x){return(as.dist(1-cor(t(x), method='pearson', use = 'pairwise.complete.obs')))}}
  if(distmethod=='spearman'){fd = function(x){return(as.dist(1-cor(t(x), method='spearman', use = 'pairwise.complete.obs')))}}
  if(distmethod=='direct'){fd = function(x){return(as.dist(x))}}
  if(distmethod=='euclidean'){fd = function(x){return(dist(x))}}
  
  # Clustering method
  fh = function(x){return(stats::hclust(x,method=clustermethod))}
  
  # Rowside colors
  
  if(cutoffmethod=='depth'){fc = function(M){return(cutreeHybrid(dendro=fh(fd(M)), distM=as.matrix(fd(M)), deepSplit=cutoff, verbose=0, minClusterSize=1)$labels)}}
  if(cutoffmethod=='height'){fc = function(M){return(cutree(fh(fd(M)), h=cutoff))}}
  if(cutoffmethod=='number'){fc = function(M){return(cutree(fh(fd(M)), k=cutoff))}}
  if(cutoffmethod=='none'){fc = function(M){return(rep(1, nrow(M)))}}
  
  if(dendrogram%in%c('none','column')){rowcol=rep('grey70', nrow(dataM))} else {
    rowcol = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")[suppressWarnings(fc(dataM))]
    if(length(rowcol)!=nrow(dataM)){rowcol=rep("gray", nrow(dataM))}
  }
  
  col <- as.character(recal_John40_remove1$DECISION)
  col[col == "1"] <- "black"
  col[col == "3"] <- "darkgrey"
  
  hm = heatmap.2(dataM,
                 col= c("#F5F5F5", "blue"),
                 hclustfun=fh,
                 distfun=fd,
                 trace='none',
                 Rowv=Rowv,
                 Colv=Colv,
                 dendrogram=dendrogram,
                 RowSideColors=rowcol,
                 #ColSideColors = col,
                 symbreaks=T, main=main,cexRow=0.6,
                 cexCol=0.8, 
                 breaks = c(0, 0.5, 1),
                 labRow = FALSE, 
                 labCol = FALSE, 
                 key = FALSE)
  
  names(rowcol) = rownames(dataM)
  #return(rowcol)
  return(rev(rowcol[hm$rowInd]))
}

