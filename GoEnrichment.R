GOenrichment=function(features,background,flag,ontology)
{
  
  if (flag=="ENSEMBL")
  {
    
    ego <- enrichGO(gene          = features,
                    universe      = background,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = ontology,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

    bp2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)

   # gosters=gost(features,organism = "hsapiens",significant = FALSE,domain_scope = "custom",custom_bg=background,evcodes = TRUE)
    
  }else{


     ego <- enrichGO(gene          = features,
                   universe      = background,
                   OrgDb         = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                 ont           = ontology,
                pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
    #
    
    # ego<-gseGO(geneList     = features,
    #       OrgDb        = org.Hs.eg.db,
    #       ont          = ontology,
    #       nPerm        = 1000,
    #       minGSSize    = 100,
    #       maxGSSize    = 500,
    #       pvalueCutoff = 0.05,
    #       verbose      = FALSE)
    # 
    bp2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    
    #gosters=gost(features,organism = "hsapiens",significant = FALSE,domain_scope = "custom",custom_bg=background,evcodes = TRUE)
    
  }
 # return(gosters)
  return(bp2)
  
  
}

findReleventBP=function(enrichment,features,decisiontable,pval)
{
  #  enrichment=enrichmentTARGET
  # features=append(MCFSFeatures[1:20],"decision")
  #  decisiontable=TARGET_GE_Classifier
  result=rosetta(decisiontable[,features],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",reducer="Genetic",ruleFiltration=TRUE, discreteParam=3)
  result=recalculateRules(decisiontable[,features],result$main)
  features=getFeatures(decisiontable[,features],result,1.0)
  #based that we have only 2 decisions can be extended for multi class
  features=unique(append(levels(features[[1]]),levels(features[[2]])))
  
  print(class(features))
  print(features)
  BP=vector(length(features), mode='list')
  print(length(BP))
  print(length(BP))
  print(length(features))
  for(i in 1:length(features))
  {
    print(features[i])
    print("Iam here")
    for(j in 1:dim(enrichment)[1])
    {
      
      temp=strsplit(enrichment[i]$geneID, "/")[[1]]
      if(features[i] %in% temp)
      {BP[i]=enrichment[i]$Description
      print('yes')
      }else{BP[i]="None"}
      
    }
  }
  
  return(as.data.frame(cbind(features,BP)))
}

#DEPRECATED CHECK script function_classifier
#----------------------------------------------------------------
#findReleventBP(enrichmentTARGET,append(MCFSFeatures[1:20],"decision"),TARGET_GE_Classifier)
options(java.parameters = "-Xmx2048m")
library("xlsx")
 
SLE_Annotation=read.xlsx("SLE-genes-Annotation2.xlsx", sheetName="WhichGeneWhichCluster", header=TRUE)
Orange=SLE_Annotation[which(SLE_Annotation[,"ORANGE"]=="y"),"All.Genes"]
Blue=SLE_Annotation[which(SLE_Annotation[,"BLUE"]=="y"),"All.Genes"]
Red=SLE_Annotation[which(SLE_Annotation[,"RED"]=="y"),"All.Genes"]
Purple=SLE_Annotation[which(SLE_Annotation[,"PURPLE"]=="y"),"All.Genes"]
Green=SLE_Annotation[which(SLE_Annotation[,"GREEN"]=="y"),"All.Genes"]
DA1Genes=SLE_Annotation[which(SLE_Annotation[,"DA1"]=="y"),"All.Genes"]
DA3Genes=SLE_Annotation[which(SLE_Annotation[,"DA3"]=="y"),"All.Genes"]
DA13=fread("DA13.csv")
background=colnames(DA13)
SLE_Johnson$enrichment[[1]]=GOenrichment(SLE_Johnson$MCFSFeatures[1:SLE_Johnson$numberOfFeatures],colnames(SLE_Johnson$classifier)[1:length(SLE_Johnson$classifier)-1],SLE_Johnson$keyType,SLE_Johnson$ontology)
x=as.data.frame(cbind(SLE_Johnson$enrichment[[1]]$result$term_id,SLE_Johnson$enrichment[[1]]$result$term_name,SLE_Johnson$enrichment[[1]]$result$term_size,SLE_Johnson$enrichment[[1]]$result$p_value,SLE_Johnson$enrichment[[1]]$result$intersection))
write.xlsx(x,"EnrichedTerms.xlsx", sheetName="DA13",append = TRUE)
SLE_Johnson$enrichment[[2]]=GOenrichment(DA1Genes,colnames(SLE_Johnson$classifier)[1:length(SLE_Johnson$classifier)-1],SLE_Johnson$keyType,SLE_Johnson$ontology)
x=as.data.frame(cbind(SLE_Johnson$enrichment[[2]]$result$term_id,SLE_Johnson$enrichment[[2]]$result$term_name,SLE_Johnson$enrichment[[2]]$result$term_size,SLE_Johnson$enrichment[[2]]$result$p_value,SLE_Johnson$enrichment[[2]]$result$intersection))
write.xlsx(x,"EnrichedTerms.xlsx", sheetName="DA1",append = TRUE)





SLE_Johnson$enrichment[[3]]=GOenrichment(DA3Genes,colnames(SLE_Johnson$classifier)[1:length(SLE_Johnson$classifier)-1],SLE_Johnson$keyType,SLE_Johnson$ontology)
x=as.data.frame(cbind(SLE_Johnson$enrichment[[3]]$result$term_id,SLE_Johnson$enrichment[[3]]$result$term_name,SLE_Johnson$enrichment[[3]]$result$term_size,SLE_Johnson$enrichment[[3]]$result$p_value,SLE_Johnson$enrichment[[3]]$result$intersection))
write.xlsx(x,"EnrichedTerms.xlsx", sheetName="DA3",append = TRUE)
p1<-dotplot(SLE_Johnson$enrichment[[2]],showCategory=20,font.size=10.5)+ ggtitle("DA1")
p2<-dotplot(SLE_Johnson$enrichment[[3]],showCategory=20,font.size=10.5)+ ggtitle("DA3")

prow <- cowplot::plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 10))
)

# add the legend to the row we made earlier. Give it one-third of
# the width of one plot (via rel_widths).
cowplot::plot_grid(prow, legend, rel_widths = c(4, .3))



#Orange
SLE_Johnson$enrichment[[4]]=GOenrichment(Orange,colnames(SLE_Johnson$classifier)[1:length(SLE_Johnson$classifier)-1],SLE_Johnson$keyType,SLE_Johnson$ontology)
x=as.data.frame(cbind(SLE_Johnson$enrichment[[4]]$result$term_id,SLE_Johnson$enrichment[[4]]$result$term_name,SLE_Johnson$enrichment[[4]]$result$term_size,SLE_Johnson$enrichment[[4]]$result$p_value,SLE_Johnson$enrichment[[4]]$result$intersection))
write.xlsx(x,"EnrichedTerms.xlsx", sheetName="Orange",append = TRUE)

#Blue
SLE_Johnson$enrichment[[5]]=GOenrichment(Blue,colnames(SLE_Johnson$classifier)[1:length(SLE_Johnson$classifier)-1],SLE_Johnson$keyType,SLE_Johnson$ontology)
x=as.data.frame(cbind(SLE_Johnson$enrichment[[5]]$result$term_id,SLE_Johnson$enrichment[[5]]$result$term_name,SLE_Johnson$enrichment[[5]]$result$term_size,SLE_Johnson$enrichment[[5]]$result$p_value,SLE_Johnson$enrichment[[5]]$result$intersection))
write.xlsx(x,"EnrichedTerms.xlsx", sheetName="Blue",append = TRUE)
#Red
SLE_Johnson$enrichment[[6]]=GOenrichment(Red,colnames(SLE_Johnson$classifier)[1:length(SLE_Johnson$classifier)-1],SLE_Johnson$keyType,SLE_Johnson$ontology)
x=as.data.frame(cbind(SLE_Johnson$enrichment[[6]]$result$term_id,SLE_Johnson$enrichment[[6]]$result$term_name,SLE_Johnson$enrichment[[6]]$result$term_size,SLE_Johnson$enrichment[[6]]$result$p_value,SLE_Johnson$enrichment[[6]]$result$intersection))
write.xlsx(x,"EnrichedTerms.xlsx", sheetName="Red",append = TRUE)
#Purple
SLE_Johnson$enrichment[[7]]=GOenrichment(Purple,colnames(SLE_Johnson$classifier)[1:length(SLE_Johnson$classifier)-1],SLE_Johnson$keyType,SLE_Johnson$ontology)
x=as.data.frame(cbind(SLE_Johnson$enrichment[[7]]$result$term_id,SLE_Johnson$enrichment[[7]]$result$term_name,SLE_Johnson$enrichment[[7]]$result$term_size,SLE_Johnson$enrichment[[7]]$result$p_value,SLE_Johnson$enrichment[[7]]$result$intersection))
write.xlsx(x,"EnrichedTerms.xlsx", sheetName="Purple",append = TRUE)
#Green
SLE_Johnson$enrichment[[8]]=GOenrichment(Green,colnames(SLE_Johnson$classifier)[1:length(SLE_Johnson$classifier)-1],SLE_Johnson$keyType,SLE_Johnson$ontology)
x=as.data.frame(cbind(SLE_Johnson$enrichment[[8]]$result$term_id,SLE_Johnson$enrichment[[8]]$result$term_name,SLE_Johnson$enrichment[[8]]$result$term_size,SLE_Johnson$enrichment[[8]]$result$p_value,SLE_Johnson$enrichment[[8]]$result$intersection))
write.xlsx(x,"EnrichedTerms.xlsx", sheetName="Green",append = TRUE)

p4<-dotplot(SLE_Johnson$enrichment[[4]],showCategory=20,font.size=10.5)+ ggtitle("Orange")
p5<-dotplot(SLE_Johnson$enrichment[[5]],showCategory=10,font.size=10.5)+ ggtitle("Blue")
p6<-dotplot(SLE_Johnson$enrichment[[6]],showCategory=20,font.size=10.5)+ ggtitle("Red")
p7<-dotplot(SLE_Johnson$enrichment[[7]],showCategory=20,font.size=10.5)+ ggtitle("Purple")
p8<-dotplot(SLE_Johnson$enrichment[[8]],showCategory=20,font.size=10.5)+ ggtitle("Green")

prow2 <- cowplot::plot_grid(
  p4+ theme(legend.position="none"),
  p5 + theme(legend.position="none"),
  p6 + theme(legend.position="none"),
  p7 + theme(legend.position="none"),
  p8 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B","C","D","E"),
  hjust = -1,
  nrow = 2
)
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p4 + theme(legend.box.margin = margin(0, 0, 0, 10))
)

# add the legend to the row we made earlier. Give it one-third of
# the width of one plot (via rel_widths).
cowplot::plot_grid(prow2, legend, rel_widths = c(4, .3))



gostplot(SLE_Johnson$enrichment[[1]], capped = TRUE, interactive = TRUE)
gostplot(SLE_Johnson$enrichment[[2]], capped = TRUE, interactive = TRUE)
gostplot(SLE_Johnson$enrichment[[3]], capped = TRUE, interactive = TRUE)

                         