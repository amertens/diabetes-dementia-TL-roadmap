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

cleandata <- clean_outcome(outvar="death")
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

names(cleandata)
Anodes <- c("A1","A2","A3","A4","A5","A6","A7","A8","A9")
Ynodes <- c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9")
Cnodes <-NULL
Lnodes <- c("sex","stroke","age","L1a","L2a","L3a","L4a","L5a","L6a",
            "L7a","L8a","L9a","L1b","L2b","L3b","L4b","L5b","L6b",
            "L7b","L8b","L9b")
subset <- cleandata %>% select(Anodes,Lnodes,Ynodes)
names(subset)


abar <- list(a=rep(1,(length(Anodes))), b=rep(0,(length(Anodes))))
result <- ltmle(subset, Anodes = Anodes, Ynodes = Ynodes, 
                Cnodes=Cnodes, Lnodes=Lnodes, abar = abar,
                survivalOutcome=F)

summary(result)

#need to account for censoring?
