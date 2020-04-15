require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
library(GGally)
library(rmcfs)
library(R.ROSETTA)
library(arules)
library(VisuNet)
library(tidyverse)
library(broom)
library(lme4)

rule_cluster <- readRDS("/Users/hp/Desktop/mapp/rule_cluster.rds")

meta <- metadata[which(decisionSLE == 3),]
meta <- meta[!rownames(meta) %in% rows_remove1$intersect..2..,]
meta <- preProcessMetaData(meta)

#---------- DATA ------------------------------

df <-  cbind.data.frame(meta$subject, 
                        meta$daysSinceDiagnosis, 
                        meta$age, 
                        meta$visit, 
                        meta$daysSinceLastVisit, 
                        meta$gender, 
                        meta$race,
                        meta$biopsyHistory,
                        meta$treatment,
                        meta$sledai, 
                        meta$`nephritis_class:.1`,
                        meta$`sledai_component_class:`)

colnames(df) <- c("subject", 
                  "daysSinceDiagnosis", 
                  "age", 
                  "visit", 
                  "daysSinceLastVisit",
                  "gender",
                  "race", 
                  "biopsyHistory",
                  "treatment",
                  "sledai", 
                  "nephritis_class",
                  "sledai_component_class")  

df <- cbind.data.frame(df, 
                       as.data.frame(unclass(meta[,57:80])))

#df <- as.data.frame(unclass(meta[,57:80]))
rownames(df) <- rownames(meta)
                       
# check if continous variables ar correlated --> no 
ggpairs(df[, c("daysSinceDiagnosis", 
               "age", 
               "visit", 
               "daysSinceLastVisit",
               "sledai")])

# add decision (rule cluster)
df <- merge(df, rule_clusterDA3, by = "row.names")
rownames(df) <- df$Row.names
df <- df[,2:38]

dfDA1 <- df
df1 <- df[df$rule_cluster == "blue" | df$rule_cluster == "red",]
df2 <- df[df$rule_cluster == "blue" | df$rule_cluster == "orange",]
df3 <- df[df$rule_cluster == "red" | df$rule_cluster == "orange",]
#--- MCFS -----------------------------------
result_mcfs3 <- mcfs(rule_cluster~., 
                     df3)
plot(result_mcfs3, type = "features", size = 15)

g <- build.idgraph(result_mcfsDA1, size = 10, size_ID = 12, orphan_nodes = TRUE)
plot(g)

ggplot(df, aes(x = rule_cluster, fill = sledai_component_class)) + 
  geom_bar()

ggplot(df, aes(x = visit, y = subject, color = rule_cluster)) + 
      geom_point()

#------- ROSETTA ----------------------------

df_DA1 <- df[,colnames(df) %in% result_mcfsDA1$RI$attribute[1:5]]
df_DA1 <- cbind(df_DA1, df$rule_cluster)
colnames(df_DA1)[6] <- "rule_cluster"

resultDA1 <- rosetta(df_DA1,
                     discrete = TRUE,
                     underSample = TRUE, 
                     underSampleNum = 3)

resultDA1$quality
resultDA1$rules

recalDA1 <- recalculateRules(df_DA1, resultDA1$main, discrete = TRUE)

View(resultDA1$main)

#------- REGRESSION -----------------------

model1 <- glmer(rule_cluster ~
                  #gender +
                  #race +
                  #leukopenia. +
                  #biopsyHistory +
                  #daysSinceDiagnosis +
                  #new_rash +
                  #vasculitis +
                  #daysSinceLastVisit +
                  (1 | subject), 
                data = df1, 
                family = binomial, 
                control = glmerControl(optimizer = "bobyqa"),
                nAGQ = 10) 

summary(model1)

model2 <- glmer(rule_cluster ~ 
                  daysSinceDiagnosis +
                  #age +
                  #visit +
                  daysSinceLastVisit + 
                  gender + 
                  race + 
                  biopsyHistory +
                  treatment +
                  sledai + 
                  nephritis_class + 
                  new_rash +
                  #low_complement. +
                  increased_dna_binding. +
                  fever. +
                  leukopenia. +
                  (1 | subject), 
                data = dfDA1, 
                family = binomial, 
                control = glmerControl(optimizer = "bobyqa"),
                nAGQ = 10) 
summary(model2)



anova(model1, model2)

plot_model(model2, type = "pred")


ggplot(dfDA1, aes(x = biopsyHistory, fill = rule_cluster)) + 
  geom_bar()

ggplot(dfDA1, aes(x = daysSinceDiagnosis, 
                y = subject, 
                color = rule_cluster)) +
  geom_point()


#--------- GLM -------------------------------

dfDA1 <- dfDA1[, sapply(dfDA1, nlevels) != 1]
dfDA1 <- na.omit(dfDA1)

result_mcfsDA1 <- mcfs(rule_cluster~., 
                     dfDA1)

plot(result_mcfsDA1, type = "features", size = 10)

g <- build.idgraph(result_mcfsDA1, size = 10, size_ID = 12, orphan_nodes = TRUE)
plot(g)

#create test and training sets 
s <- NULL 
k <- sample(1:107, 70)

for(i in 1:length(k)){
  s <- c(s, rownames(df1)[k[i]])
}

df_train <- df1[rownames(df1) %in% s,]
df_test <- df1[!rownames(df1) %in% s,]

model1 <- glmer(rule_cluster ~ 
                  #gender +
                  race +
                  leukopenia. +
                 #daysSinceDiagnosis +
                 #new_rash +
                 #vasculitis +
                 #age +
                 #visit +
                 #biopsyHistory +
                 #biopsyHistory +
                 #daysSinceLastVisit + 
                 #fever. +
                 #thrombocytopenia. +
                 #urinary_casts +
                 #sledai +
               (1 | subject), 
             data = dfDA1, 
             family = binomial, 
             control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)

model2 <- glmer(rule_cluster ~ 
                  #gender +
                  race +
                  #leukopenia. +
                 #mucosal_ulcers. +
                 #age +
                 #visit +
                 #daysSinceDiagnosis +
                 #new_rash +
                 #biopsyHistory +
                 #thrombocytopenia. +
                 #leukopenia. + 
                 #daysSinceLastVisit + 
                 #fever. +
                 #thrombocytopenia. +
                 #urinary_casts +
                 #sledai +
                 (1 | subject), 
               data = dfDA1, 
               family = binomial, 
               control = glmerControl(optimizer = "bobyqa"),
               nAGQ = 10)

anova(model1, model2)

summary(model2)

prob <- predict(model1, 
                newdata = dfDA1,
                allow.new.levels = TRUE,
                type = "response")

pred_test <- ifelse(prob > 0.5, "red", "blue")

pred_test <- as.data.frame(pred_test)

c <- cbind.data.frame(as.vector(pred_test), as.vector(dfDA1$rule_cluster))
colnames(c) <- c("predicted", "true")
as.data.frame(table(c))

ggplot(df2, aes(x = daysSinceDiagnosis, y = subject, color = rule_cluster)) + 
  geom_point() + 
  scale_color_manual(values = c("blue", "orange"))

ggplot(dfDA1, aes(x = race, fill = rule_cluster)) + 
  geom_bar() + 
  scale_fill_manual(values = c("blue", "red"))

plot_model(model2, type = "pred")

predicted.classes <- ifelse(probabilities > 0.5, "red", "blue")
head(predicted.classes) 

# Select only numeric predictors

probabilities <- predict(model, 
                         type = "response")

mydata <- df_train %>%
  dplyr::select_if(is.numeric) 
predictors <- colnames(mydata)
mydata$row <- rownames(mydata)

# Bind the logit and tidying the data for plot
mydata <- mydata %>%
  mutate(logit = log(probabilities/(1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit, -row)

#inspect linearity 
theme_set(theme_classic())
ggplot(mydata, aes(logit, predictor.value))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() + 
  facet_wrap(~predictors, scales = "free_y")

car::vif(model)

plot(model, which = 4, id.n = 3)

# Extract model results
model.data <- augment(model) %>% 
  mutate(index = 1:n()) 

model.data %>% top_n(3, .cooksd)

ggplot(model.data, aes(index, .std.resid)) + 
  geom_point(aes(color = rule_cluster), alpha = .5) +
  theme_bw()

model.data %>% 
  filter(abs(.std.resid) > 3)

#--------------------------------------------------------------- 

df <- 
df <- df[,which(colSums(df)!=0)]
df <- ifelse(df != 0, 1, df)
df <- data.matrix(df)

heatmap.2(df)
