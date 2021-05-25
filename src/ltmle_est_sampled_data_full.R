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

table(data$A, data$any_dementia)
table(data$A, data$dementia)
table(data$A, data$death)

table(data$Aglp_1, data$dementia)
table(data$Asglt_1, data$dementia)
table(data$Aglp_2, data$dementia)
table(data$Asglt_2, data$dementia)
table(data$Aglp_3, data$dementia)
table(data$Asglt_3, data$dementia)

# data$death[data$death>3] <- 3
# table(data$A, data$death)
# 
# data$dementia[data$dementia>3] <- 3
# table(data$A, data$dementia)

# need to expand out the data to get distinct A and Y nodes 

#choose outcome and create outcome var
clean_outcome <- function(outvar){
for(i in 1:9){
  data[[outvar]] <-as.numeric(data[[outvar]])
  #head(data[[outvar]])
  #table(data[[outvar]])
  #table(is.na(data[[outvar]]))
  data[,paste0("Y",i)]<- as.numeric(ifelse(data[[outvar]]<=i & data[[outvar]]>0,1,0))
  # table(data[,paste0("Y",i)],data[[outvar]])
  table(is.na(data[,paste0("Y",i)]),is.na(data[[outvar]]))
  class(data[,paste0("Y",i)])
  data[is.na(data[paste0("Y",i)]),paste0("Y",i)] <-0
}
  return(data)
}

cleandata <- clean_outcome(outvar="dementia")

table(cleandata$Aglp_1, cleandata$Y3)
table(cleandata$Asglt_1 , cleandata$Y3)
table(cleandata$Aglp_2, cleandata$Y3)
table(cleandata$Asglt_2 , cleandata$Y3)
table(cleandata$Aglp_3, cleandata$Y3)
table(cleandata$Asglt_3 , cleandata$Y3)


#rename the time-varying covariates
#BMI and kidney disease:
cleandata <- data.table(cleandata)
setnames(cleandata, old = c("obese_1","obese_2","obese_3","obese_4","obese_5","obese_6",
                    "obese_7","obese_8","obese_9","kidney_1","kidney_2","kidney_3",    
                    "kidney_4","kidney_5","kidney_6","kidney_7","kidney_8",    
                     "kidney_9" ), new = c("L1a","L2a","L3a","L4a","L5a","L6a",
                        "L7a","L8a","L9a","L1b","L2b","L3b","L4b","L5b","L6b",
                        "L7b","L8b","L9b"))
setnames(cleandata, old = c("sex","stroke","age","income", "diabduration"),new=c("L1c","L1d","L1f","L1g","L1h"))
names(cleandata)
Anodes <- c("Aglp_1","Asglt_1","Aglp_2","Asglt_2","Aglp_3","Asglt_3",
            "Aglp_4","Asglt_4","Aglp_5","Asglt_5","Aglp_6","Asglt_6",
            "Aglp_7","Asglt_7","Aglp_8","Asglt_8","Aglp_9","Asglt_9")
Ynodes <- c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9")
Cnodes <-NULL
Lnodes <- c("L1a","L1b","L1c","L1d","L1f","L1g","L1h",
             "L2a","L2b","L3a","L3b","L4a","L4b","L5a","L5b",
            "L6a","L6b","L7a","L7b","L8a","L8b","L9a","L9b")

subset <- cleandata %>% select(
  "L1a","L1b","L1c","L1d","L1f","L1g","L1h","Aglp_1","Asglt_1",
  "L2a","L2b","Aglp_2","Asglt_2",
  "L3a","L3b","Aglp_3","Asglt_3",
  "L4a","L4b","Aglp_4","Asglt_4",
  "L5a","L5b","Aglp_5","Asglt_5",
  "L6a","L6b","Aglp_6","Asglt_6",
  "L7a","L7b","Aglp_7","Asglt_7",
  "L8a","L8b","Aglp_8","Asglt_8",
  "L9a","L9b","Aglp_9","Asglt_9",
  Ynodes)

names(subset)
#SL.library<- c( "SL.step", "SL.mean")#, "SL.ranger","SL.nnet", "SL.biglasso")#"SL.glm","SL.xgboost",
SL.library<- c( "SL.glm")
#SL.library<- c(  "SL.mean","SL.glm","SL.glmnet", "SL.randomForest")

set.seed(12345)
abar <- list(a=rep(c(1,0),(length(Anodes)/2)), b=rep(c(0,1),(length(Anodes)/2)))
# result <- ltmle(subset, Anodes = Anodes, Ynodes = Ynodes,
#                 Cnodes=Cnodes, Lnodes=Lnodes, abar = abar,
#                 survivalOutcome=F,SL.library=SL.library,variance.method = "ic")
# 
# summary(result)
# 
# 
# #CRUDE ESTIMATE
# data_crude <- subset %>% select(Anodes, Ynodes)
# out <-  ltmle(data_crude, Anodes = Anodes, Ynodes = Ynodes, 
#               Cnodes=Cnodes, Lnodes=NULL, abar = abar,
#               survivalOutcome=F,SL.library=SL.library,gcomp=F,
#               variance.method = "ic")
# summary(out)
# 
# 
# gcomp_output=F
# est_type="tmle"

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

system.time((est_iptw <- get_estimates(gcomp_output=F,est_type="iptw")))

system.time((est_gcomp <- get_estimates(gcomp_output=T,est_type="gcomp")))

save(est_tmle, est_iptw, est_gcomp, file=here("results/sampled_data_estimates.R"))

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