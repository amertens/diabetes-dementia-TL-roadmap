rm(list=ls())
library(ltmle)
library(here)
library(dplyr)
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
table(data$A, data$dementia)
table(data$A, data$death)




# need to expand out the data to get distinct A and Y nodes 

#choose outcome and create outcome var
clean_outcome <- function(outvar){
for(i in 1:9){
  data[[outvar]] <-as.numeric(data[[outvar]])
  #head(data[[outvar]])
  #table(data[[outvar]])
  #table(is.na(data[[outvar]]))
  data[,paste0("Y",i)]<- as.numeric(ifelse(data[[outvar]]<=i,1,0))
  # table(data[,paste0("Y",i)],data[[outvar]])
  table(is.na(data[,paste0("Y",i)]),is.na(data[[outvar]]))
  class(data[,paste0("Y",i)])
  data[is.na(data[paste0("Y",i)]),paste0("Y",i)] <-0
  names(data)[names(data) == paste0("T",i)] <- paste0("A",i)
}
  return(data)
}

cleandata <- clean_outcome(outvar="dementia")
names(cleandata)
table(cleandata$Y1,cleandata$Y2)

#rename the time-varying covariates
#BMI and kidney disease:
cleandata <- data.table(cleandata)
setnames(cleandata, old = c("obese_1","obese_2","obese_3","obese_4","obese_5","obese_6",
                    "obese_7","obese_8","obese_9","kidney_1","kidney_2","kidney_3",    
                    "kidney_4","kidney_5","kidney_6","kidney_7","kidney_8",    
                     "kidney_9" ), new = c("L1a","L2a","L3a","L4a","L5a","L6a",
                        "L7a","L8a","L9a","L1b","L2b","L3b","L4b","L5b","L6b",
                        "L7b","L8b","L9b"))
setnames(cleandata, old = c("sex","stroke","age","income","diabduration"),new=c("L1c","L1d","L1f","L1g","L1h"))
names(cleandata)
Anodes <- c("A1","A2","A3","A4","A5","A6","A7","A8","A9")
Ynodes <- c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9")
Cnodes <-NULL
Lnodes <- c("L1a","L1b","L1c","L1d","L1f","L1g","L1h",
             "L2a","L2b",
             "L3a","L3b",
             "L4a","L4b",
             "L5a","L5b",
             "L6a","L6b",
             "L7a","L7b",
             "L8a","L8b",
             "L9a","L9b")

# subset <- cleandata %>% select(
#                                "L1a","L1b","L1c","L1d","L1f","A1","Y1",
#                                "L2a","L2b","A2","Y2",
#                                "L3a","L3b","A3","Y3",
#                                "L4a","L4b","A4","Y4",
#                                "L5a","L5b","A5","Y5",
#                                "L6a","L6b","A6","Y6",
#                                "L7a","L7b","A7","Y7",
#                                "L8a","L8b","A8","Y8",
#                                "L9a","L9b","A9","Y9")
subset <- cleandata %>% select(
  "L1a","L1b","L1c","L1d","L1f","L1g","L1h","A1",
  "L2a","L2b","A2",
  "L3a","L3b","A3",
  "L4a","L4b","A4",
  "L5a","L5b","A5",
  "L6a","L6b","A6",
  "L7a","L7b","A7",
  "L8a","L8b","A8",
  "L9a","L9b","A9",Ynodes)

names(subset)
#SL.library<- c( "SL.step", "SL.mean")#, "SL.ranger","SL.nnet", "SL.biglasso")#"SL.glm","SL.xgboost",
SL.library<- c( "SL.glm", "SL.mean")
#SL.library<- c(  "SL.mean","SL.glm","SL.glmnet", "SL.randomForest")

set.seed(12345)
abar <- list(a=rep(1,(length(Anodes))), b=rep(0,(length(Anodes))))
# result <- ltmle(subset, Anodes = Anodes, Ynodes = Ynodes, 
#                 Cnodes=Cnodes, Lnodes=Lnodes, abar = abar,
#                 survivalOutcome=F,SL.library=SL.library,variance.method = "ic")
# 
# summary(result)

get_estimates <- function(gcomp_output,est_type){
  #CRUDE ESTIMATE
  data_crude <- cleandata %>% select(Anodes, Ynodes)
  
  out <-  ltmle(data_crude, Anodes = Anodes, Ynodes = Ynodes, 
                Cnodes=Cnodes, Lnodes=NULL, abar = abar,
                survivalOutcome=F,SL.library=SL.library,gcomp=gcomp_output,
                variance.method = "ic")
  #ADJUSTED ESTIMATE
  out2 <-  ltmle(subset, Anodes = Anodes, Ynodes = Ynodes, 
                 Cnodes=Cnodes, Lnodes=Lnodes, abar = abar,
                 survivalOutcome=F,SL.library=SL.library,gcomp=gcomp_output,
                 variance.method = "ic")
  
  Crude <- summary(out,est_type)$effect.measures$ATE$estimate
  Adjusted <- summary(out2,est_type)$effect.measures$ATE$estimate
  CI_crude <- summary(out,est_type)$effect.measures$ATE$CI
  CI_adj <- summary(out2,est_type)$effect.measures$ATE$CI
  
  return(list(res=data.frame(est_type,Crude,CI_crude,Adjusted,CI_adj), fit_unadj=out, fit_adj=out2))
}


system.time((est_tmle <- get_estimates(gcomp_output=F,est_type="tmle")))
summary(est_tmle$fit_unadj)$effect.measures$RR
summary(est_tmle$fit_adj)$effect.measures$RR

rio::export(est_tmle$res, file="outcome_estimates_tmle.csv")

system.time((est_iptw <- get_estimates(gcomp_output=F,est_type="iptw")))
rio::export(est_iptw$res, file="outcome_estimates_iptw.csv")

system.time((est_gcomp <- get_estimates(gcomp_output=T,est_type="gcomp")))
rio::export(est_gcomp$res, file="outcome_estimates_gcomp.csv")


estimates <-rbind(est_tmle,est_iptw,est_gcomp)
rio::export(estimates, file="outcome_estimates_ltmle.csv")


save(list = ls(all.names = TRUE), file = "res.RData", envir = .GlobalEnv)


#plot of results

# library(ggplot2)
# p = ggplot(data=allest,
#            aes(x = type,y = est, ymin = lowerCI, ymax = upperCI ))+
#   geom_pointrange(aes(col=type))+
#   geom_hline(aes(fill=type),yintercept =0, linetype=2)+
#   xlab('')+ ylab("Average Treatment Effect (95% CI)")+coord_flip()+
#   geom_errorbar(aes(ymin=lowerCI, ymax=upperCI,col=type),width=0.5,cex=1)+
#   facet_wrap(~Parameter,strip.position="left",nrow=9,scales = "free_y") +
#   geom_text(label=allest$label,nudge_x = 0.25, 
#              check_overlap = T) 
# 
# ggsave("./forestplot.png",  width = 7.5,
#        height = 4.26,)