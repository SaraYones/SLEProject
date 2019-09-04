
#-------------- DATA ---------------------
recal_John40_remove1 <- readRDS('C:/Users/hp/Desktop/mapp/recal_John40_remove1.rds')
result_John40_remove1 <- readRDS('C:/Users/hp/Desktop/mapp/result_John40_remove1.rds')
dt_40_remove1 <- readRDS('C:/Users/hp/Desktop/mapp/dt_John40_remove1.rds')

genes_model <- getFeatures(dt_40_remove1, recal)

#--------------- DA3 -------------------------------------

#Make two data frame 
genes_model_DA3 <- as.data.frame(genes_model[[2]])

#get sum of p-values and number of rules 
gene_model_DA3 <- rank_features(recal, genes_model_DA3$gene)

#calculate avrage p-value 
gene_model_DA3$av_p_value <- (gene_model_DA3$sum_p_val/gene_model_DA3$n_rules)

#Order based on avrage p-value 
gene_model_DA3 <- gene_model_DA3[order(gene_model_DA3$av_p_value),]

#Add decision 
gene_model_DA3$decision <- "DA3"

#--------------- DA1 -------------------------------------

#Make two data frame 
genes_model_DA1 <- as.data.frame(genes_model[[1]])

#get sum of p-values and number of rules 
gene_model_DA1 <- rank_features(recal, genes_model_DA1$gene)

#calculate avrage p-value 
gene_model_DA1$av_p_value <- (gene_model_DA1$sum_p_val/gene_model_DA1$n_rules)

#Order based on avrage p-value 
gene_model_DA1 <- gene_model_DA1[order(gene_model_DA1$av_p_value),]

#Add decision 
gene_model_DA1$decision <- "DA1"

#---------- FUNCTIONS ------------------------------------

rank_features <- function(recal, features){
   
   sum_p_val <- array(0, length(features))
   n_rules <- array(0, length(features))
   expression_level <- array(0, length(features))
   
   
   for(i in 1:length(features)){
       
     for(j in 1:NROW(recal)){
        s <- strsplit(as.character(recal[j,1]), ",")
             
          for(k in 1:length(s[[1]])){
                 
            if(s[[1]][k] == features[i]){
              sum_p_val[i] <- sum_p_val[i] + recal$PVAL[j]
              n_rules[i] <- n_rules[i] + 1  
              expression_level[i] <- strsplit(as.character(recal$DISC_CLASSES[j]), ",")[[1]][k]
            }
          }
        }
       }
  return(cbind.data.frame(features, expression_level, sum_p_val, n_rules))
}
