rm(list=ls())
library(ltmle)
library(dplyr)
# load("./data/clean/finaldata.RData")
data <- rio::import("./diabetes-dementia-TL-roadmap/simulated data/novo_registry_simulated.RDS")

str(data)
head(data)
table(data$death)
table(data$A)
table(data$T1)

# need to expand out the data to get distinct A and Y nodes 
data$death <-as.numeric(data$death)
head(data$death)
table(data$death)
table(is.na(data$death))

#looking at death first:
for(i in 1:9){
  data[,paste0("Y",i)]<- ifelse(data$death>i,0,1)
  table(data[,paste0("Y",i)],data$death)
  table(is.na(data[,paste0("Y",i)]),is.na(data$death))

  data[is.na(data[,paste0("Y",i)]),] <-0
  names(data)[names(data) == paste0("T",i)] <- paste0("A",i)
}
names(data)
table(data$Y1,data$Y2)

Anodes <- c("A1","A2","A3","A4","A5","A6","A7","A8","A9")
Ynodes <- c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9")
Cnodes <-NULL
Lnodes <- NULL #c("sex","stroke","age")
subset <- data %>% select(Anodes,Ynodes)
head(subset)

abar <- matrix(1,nrow=nrow(data),ncol=length(Anodes))

result <- ltmle(subset, Anodes = Anodes, Ynodes = Ynodes, 
                Cnodes=Cnodes, Lnodes=Lnodes, abar = abar,
                survivalOutcome=TRUE)

summary(result)


