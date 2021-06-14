

#TODO:
#Add censoring. Make censoring associated with drug use.


rm(list=ls())
library(here)
library(tidyverse)
library(simstudy)
library(simcausal)
library(tmle)
library(ltmle)
#library(knitr)
#https://cran.r-project.org/web/packages/simhelpers/vignettes/MCSE.html
library(simhelpers)

#function for Truncated continuous distribution
rnorm_trunc <- function(n, mean, sd, minval = 17){
     out <- rnorm(n = n, mean = mean, sd = sd)
     minval <- minval[1]
     out[out < minval] <- minval
     out
}




#SEM for complicated longitudinal data with time-dependent confounders La and Lb
D <- DAG.empty() + 
  node("CVD", distr="rcat.b1", probs = c(0.5, 0.25, 0.25)) +
  node("U.A", distr = "rnorm", mean = 0, sd = 1) +
  node("U.Y", distr = "rnorm", mean = 0, sd = 5) +
  node("A", t=0, distr="rbern", prob=plogis(.08*CVD -.4)) + 
  node("La", t=0, distr="rbern", prob=plogis(-2 - 0.3*CVD - 0.5*A[t])) +
  node("A", t=1:2, distr="rbern", prob=plogis(A[t-1]  - 0.5*La[t-1])) +
  node("La", t=1:2, distr="rbern", prob=plogis(-3 - 0.3*CVD +  -0.5*A[t] + 1.5*La[t-1])) +
  node("Y", t=1:3, distr="rbern", prob=plogis( -3 - 1.2*La[t-1] + 0.1*CVD - A[t-1]), EFU=TRUE)

#Set DAG and save
D <- set.DAG(D)
saveRDS(D, file=here("simulated data/simulated_dag_RD.RDS"))

act_t0_theta <- node("A",t=0:2, distr="rbern", prob=ifelse(theta==1,1,0))
D <- D + action("A_th0", nodes=c(act_t0_theta), theta=0)
D <- D + action("A_th1", nodes=c(act_t0_theta), theta=1)
D <- set.targetE(D, outcome="Y", t=1:3, param="A_th1 - A_th0")
eval <- eval.target(D, n=500000, rndseed=12345)
eval$res

trueRD <- eval$res[3]
trueRD



# set up simulation parameters
N_sim<-100 #100 datasets to estimate on

#Get a list of different dataset sizes
set.seed(12345)
ndata<- ceiling(runif(N_sim)*10)^4 + ceiling(runif(N_sim)*100)^2 +100
summary(ndata)


#-------------------------------------
#Unadjusted regression
#-------------------------------------

sim_res_glm <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res=try(glm(Y_3~A_0, data=dat, family = "gaussian"))
  sim_res_glm$est[i] <- try((res$coefficients)[2])
  sim_res_glm$var[i] <- try(vcov(res)[2,2])
}
sim_res_glm
saveRDS(sim_res_glm,  file=here("results/simulation_results_glm_RD.RDS"))


#-------------------------------------
#Adjusted regression
#-------------------------------------  

sim_res_glm_adj <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res = try(glm(Y_3~A_0 + CVD + La_0, data=dat, family = "gaussian"))
  sim_res_glm_adj$est[i] <- try((res$coefficients)[2])
  sim_res_glm_adj$var[i] <- try(vcov(res)[2,2])
  }
sim_res_glm_adj
saveRDS(sim_res_glm_adj,  file=here("results/simulation_results_glm_adj_RD.RDS"))


#-------------------------------------
# TMLE fit with regression
#-------------------------------------

sim_res_tmle_glm <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res <- try(tmle(Y=dat$Y_3, A=dat$A_0, W=subset(dat, select= c(CVD, La_0,La_1,La_2)), family="gaussian",Q.SL.library = c("SL.glm"),g.SL.library = c("SL.glm")))
  sim_res_tmle_glm$est[i] <-  try(res$estimates$ATE$psi)
  sim_res_tmle_glm$var[i] <-  try(res$estimates$ATE$var.psi)
}
sim_res_tmle_glm
saveRDS(sim_res_tmle_glm,  file=here("results/simulation_results_tmle_glm_RD.RDS"))


#-------------------------------------
# TMLE fit with SL
#-------------------------------------

sim_res_tmle <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
#lib = c("SL.glm","SL.glmnet","SL.glm.interaction","SL.randomForest")
lib = c("SL.glm","SL.glmnet","SL.step")
glib = c("SL.glm","SL.glmnet")
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res <- try(tmle(Y=dat$Y_3, A=dat$A_0, W=subset(dat, select= c(CVD,  La_0,La_1,La_2)), family="gaussian",Q.SL.library = lib,g.SL.library = glib))
  sim_res_tmle$est[i] <- try(res$estimates$ATE$psi)
  sim_res_tmle$var[i] <- try(res$estimates$ATE$var.psi)
}
sim_res_tmle
saveRDS(sim_res_tmle,  file=here("results/simulation_results_tmle_RD.RDS"))


#-------------------------------------
# LTMLE 
#-------------------------------------
Anodes <- c("A_0","A_1","A_2")
Ynodes <- c("Y_1","Y_2","Y_3")
Cnodes <-NULL
Lnodes <- c("La_0","L0c","La_1","La_2")


lib = c("SL.glm","SL.glmnet","SL.step")
abar <- list(a=rep(1,(length(Anodes))), b=rep(0,(length(Anodes))))

sim_res_ltmle <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  subset <- dat %>% rename(
    L0c=CVD
  ) %>%
    select(
      "A_0","La_0","L0c",
      "La_1","A_1",
      "La_2","A_2",Ynodes)
  res=NULL
  
  set.seed(12345)
  abar <- list(a=rep(1,(length(Anodes))), b=rep(0,(length(Anodes))))
  result <- ltmle(subset, Anodes = Anodes, Ynodes = Ynodes,
                  Cnodes=Cnodes, Lnodes=Lnodes, abar = abar,
                  survivalOutcome=F,SL.library=lib,variance.method = "ic")
  
  res<-summary(result)
  sim_res_ltmle$est[i] <- try(res$effect.measures$ATE$estimate)
  sim_res_ltmle$var[i] <- try(res$effect.measures$ATE$std.dev^2)
}
sim_res_ltmle
saveRDS(sim_res_ltmle,  file=here("results/simulation_results_tmle_RD.RDS"))


save(
  trueRD,
  sim_res_glm,  
  sim_res_glm_adj, 
  sim_res_tmle_glm, 
  sim_res_tmle, 
  sim_res_ltmle,
  file=here("results/simulation_results_RD.Rdata")
)



