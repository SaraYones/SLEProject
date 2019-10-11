
library(ggplot2)

#---------- WEB SLE ----------------------

webSLE_upDA3 <- getGenesWebSLE("webSLE_upDA3.csv")
webSLE_downDA3 <- getGenesWebSLE("webSLE_downDA3.csv")

gene_rank_DA3 <- read.csv("gene_rank_DA3.csv")
gene_rank_DA1 <- read.csv("gene_rank_DA1.csv")

gene_rank_DA3 <- compareWebSLE(gene_rank_DA3, webSLE_upDA3, webSLE_downDA3)
gene_rank_DA1 <- compareWebSLE(gene_rank_DA1, webSLE_downDA3, webSLE_upDA3)

write.csv(gene_rank_DA1, "gene_rank_DA1.csv")
write.csv(gene_rank_DA3, "gene_rank_DA3.csv")

#---------- PLOT ------------------------------------------------------

#Plot line plot DA3 
df_DA3 <- GeneLinePlotDA3(gene_rank_DA3)

ggplot(df_DA3[df_DA3$DiscState != "NA",], 
       aes(x = label, 
           y = DiscState, 
           fill = source, 
           group = source)) +
  geom_line(size = 0.8) + 
  geom_point(aes(shape = source), size = 3, stroke = 1.5) +
  scale_shape_manual(values = c(21, 25)) +
  scale_fill_manual(values = alpha(c("red", "blue"), 0.6)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  xlab("") +
  ylab("Expression Level")

#Plot line pot DA1 
df_DA1 <- GeneLinePlotDA1(gene_rank_DA1)

ggplot(df_DA1[df_DA1$DiscState != "NA",], 
       aes(x = label, 
           y = DiscState, 
           fill = source, 
           group = source)) +
  geom_line(size = 0.8) + 
  geom_point(aes(shape = source), size = 3, stroke = 1.5) +
  scale_shape_manual(values = c(21, 25)) +
  scale_fill_manual(values = alpha(c("red", "blue"), 0.6)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  xlab("") +
  ylab("Expression Level") + 


#------- FUNCTION ---------------------------------------------------

compareWebSLE <- function(gene_rank, webSLEgenesUp, webSLEgenesDown){
  webSLE_expression <- array("NA", NROW(gene_rank))
  for(i in 1:NROW(gene_rank)){
    if(gene_rank$label[i] %in% webSLEgenesUp){
       webSLE_expression[i] <- "3"
    }
    if(gene_rank$label[i] %in% webSLEgenesDown){
      webSLE_expression[i] <- "1"
    }
  }
  gene_rank <- cbind(gene_rank, webSLE_expression)
  return(gene_rank)
}

getGenesWebSLE <- function(file){
  webSLE <- read.csv(file)
  features <- webSLE$SYMBOL
  return(features)
}

GeneLinePlotDA3 <- function(gene_rank_DA3){
  df.1 <- gene_rank_DA3[,4:5]
  df.1$source <- "Network DA3"
  
  df.2 <- gene_rank_DA3[,c(4,10)]
  df.2$source <- "WebSLE"
  
  colnames(df.2) <- c("label", "DiscState", "source")
  
  df_DA3 <- rbind(df.1, df.2)
  lev <- c(as.character(unique(df_DA3[!df_DA3$label %in% df_DA3[df_DA3$DiscState == "NA",]$label,]$label)),
           as.character(df_DA3[df_DA3$DiscState == "NA",]$label))
  
  df_DA3$label <- factor(df_DA3$label, 
                     levels = lev[1:30])
  
  return(df_DA3)
}


GeneLinePlotDA1 <- function(gene_rank_DA1){
  df.3 <- gene_rank_DA1[,4:5]
  df.3$source <- "Network DA1"
  
  df.4 <- gene_rank_DA1[,c(4,10)]
  df.4$source <- "WebSLE"
  
  colnames(df.4) <- c("label", "DiscState", "source")
  
  df_DA1 <- rbind(df.3, df.4)
  lev <- c(as.character(unique(df_DA1[!df_DA1$label %in% df_DA1[df_DA1$DiscState == "NA",]$label,]$label)),
           as.character(df_DA1[df_DA1$DiscState == "NA",]$label))
  
  df_DA1$label <- factor(df_DA1$label, 
                         levels = lev[1:20])
  
  return(df_DA1)
}
