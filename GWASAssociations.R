#GWAS associationlibrary(BiocManager)
library(BiocManager)
#BiocManager::install("gwascat")
#https://rcc-uchicago.github.io/R-intro/demo_gwas_catalog.html
library(gwascat)
makeCurrentGwascat()
data(ebicat37) 

#https://rmagno.eu/gwasrapidd/articles/gwasrapidd.html
library('gwasrapidd')
# install.packages("remotes")
remotes::install_github("ramiromagno/gwasrapidd")

#Visualization of associated variants
install.packages("qqman")
library(qqman)

#SLE EFO_ID='EFO_0002690'

#Extracted Genes from rules ranked according to mean accuracy * mean support using Alva script
genesDA1=fread('gene_rank_DA1.csv')
genesDA3=fread('gene_rank_DA3.csv')
EFO_ID='EFO_0002690'
extractedGenes=unique(append(genesDA1$label,genesDA3$label))

associatedGenes=associateGenesWithGWASCat(EFO_ID,extractedGenes)
associateGenesWithGWASCat=function(EFO_ID,extractedGenes)
{
  
  #You could have also used get_associations(efo_trait = 'autoimmune disease')
  my_associations <- get_associations(efo_id = EFO_ID)
 associatedGenes=extractedGenes[extractedGenes %in% my_associations@ensembl_ids$gene_name]
 #associatedGenes=genesDA3$label[genesDA3$label %in% my_associations@ensembl_ids$gene_name]
 #get relevant association IDs for genes
 relevantAssociationIDs=my_associations@ensembl_ids[which(my_associations@ensembl_ids$gene_name %in% genesDA1$label),]
 #get P-Values of these associations
 relevantAssociationPvalues=my_associations@associations[my_associations@associations$association_id %in% relevantAssociationIDs$association_id,]
 #
 variant_id=NULL
 chromosome_name=NULL
 chromosome_position=NULL
 for(i in 1:length(relevantAssociationIDs$association_id))
 {
  
   variantsGWAS=get_variants(association_id = relevantAssociationIDs$association_id[i])
   
   variant_id=append(variant_id,variantsGWAS@variants$variant_id)
   chromosome_name=append(chromosome_name,variantsGWAS@variants$chromosome_name)
   chromosome_position=append(chromosome_position,variantsGWAS@variants$chromosome_position)
   
 }
 
 GWASTable=cbind(variant_id,chromosome_name,chromosome_position,relevantAssociationPvalues$pvalue)
 GWASTable=as.data.frame( GWASTable)                   
 colnames(GWASTable)=c("SNP","CHR","BP","P")
 GWASTable$CHR=as.numeric(as.character(GWASTable$CHR))
 GWASTable$BP=as.numeric(as.character(GWASTable$BP))
 GWASTable$P=as.numeric(as.character(GWASTable$P))
 
 plotGWAS(GWASTable,"GWAS",getwd,"")
 return(associatedGenes)
}