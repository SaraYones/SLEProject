library(devtools)

install_github("GabrielHoffman/variancePartition")
library('variancePartition')
library('variancePartition')

library('doParallel')
cl <- makeCluster(4)
registerDoParallel(cl)

data(varPartData)

form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)

varPart <- fitExtractVarPartModel( geneExpr, form, info )
vp <- sortCols( varPart )

plotPercentBars( vp[1:10,] )
plotVarPart( vp )


i <- which.max( varPart$Tissue )
GE <- data.frame( Expression = geneExpr[i,], Tissue = info$Tissue)

plotStratify( Expression ~ Tissue, GE, main=rownames(geneExpr)[i]) 

i <- which.max( varPart$Individual )
GE <- data.frame( Expression = geneExpr[i,],
                  Individual = info$Individual)

label <- paste("Individual:", format(varPart$Individual[i]*100,
                                     digits=3), "%")
main <- rownames(geneExpr)[i]
plotStratify(  Expression ~ Individual, GE, colorBy=NULL,
               text=label, main=main)
