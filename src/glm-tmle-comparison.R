
#glm and static tmle comparison

rm(list=ls())
library(ltmle)
library(here)
library(tlverse)
library(tidyverse)
library(tmle)
library(data.table)
# load("./data/clean/finaldata.RData")
data <- readRDS(here("simulated data/novo_registry_simulated.RDS"))

str(data)
names(data)
head(data)
table(data$dementia)
table(data$A)
table(data$T1)

table(data$A, data$any_dementia)
data$A <- factor(data$A, levels = c("sglt2","glp1"))

d <- data %>% mutate(obese=ifelse(obese_1+obese_2+obese_3+obese_4+obese_5+obese_6+obese_7+obese_8+obese_9>0,1,0),  
                     kidney=ifelse(kidney_1+kidney_2+kidney_3+kidney_4+kidney_5+kidney_6+kidney_7+kidney_8+kidney_9>0,1,0)) %>%
                    subset(., select = c(A, age, sex, stroke, any_dementia, edu, income, diabduration, obese, kidney))

res_unadj <- glm(any_dementia~A, data=d, family=binomial(link="log"))
summary(res_unadj)
exp(res_unadj$coefficients)


res_adj <- glm(any_dementia~A + ., data=d, family=binomial(link="log"))
summary(res_adj)
exp(res_adj$coefficients)


d$A <- ifelse(d$A=="glp1",1,0)
d$W1 <- d$W2 <- 1
res_tmle_unadj <- tmle(Y=d$any_dementia, A=d$A, W=subset(d, select= c("W1","W2")), family="binomial", Q.SL.library = c("SL.glm"),g.SL.library = c("SL.glm"), prescreenW.g=F)
res_tmle_unadj$estimates$RR

res_tmle_adj_glm <- tmle(Y=d$any_dementia, A=d$A, W=subset(d, select= -c(any_dementia,A)), family="binomial",Q.SL.library = c("SL.glm"),g.SL.library = c("SL.glm"))
res_tmle_adj_glm$estimates$RR

lib = c("SL.mean","SL.glm","SL.glmnet","SL.gam")
res_tmle_adj_sl <- tmle(Y=d$any_dementia, A=d$A, W=subset(d, select= -c(any_dementia,A)), family="binomial",Q.SL.library = lib,
                         g.SL.library = lib)
res_tmle_adj_sl$estimates$RR





