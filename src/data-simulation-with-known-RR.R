

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

#Truncated continious
rnorm_trunc <- function(n, mean, sd, minval = 17){
     out <- rnorm(n = n, mean = mean, sd = sd)
     minval <- minval[1]
     out[out < minval] <- minval
     out
}

#calc mse and bootstrapped CI

#The relative bias was calculated as the average of estimated RR in log scale minus the log of the true RR divided by the log of the true RR. 
#Standard error was defined as the empirical standard error of the estimated RR in log scale over all 1,000 simulations. MSE was calculated by taking 
#the sum of the squared bias in log scale and the variances.
calc_mse <- function(sim, trueRR, model){
  #variance
  var<-mean(sim$var)
  var.ci <- quantile(sim$var, c(0.025, 0.975))
  #bias
  bias <- (mean(log(sim$est))-log(trueRR))/log(trueRR)
  bias.ci <- quantile(((log(sim$est))-log(trueRR))/log(trueRR), c(0.025, 0.975))
  
  #mse
  mse_ests <- (log(sim$est)-log(trueRR))^2 + sim$var
  mse <- mean(mse_ests)
  mse.ci <- quantile(mse_ests, c(0.025, 0.975))
  res <- bind_rows(
    data.frame(model=model, metric="Variance", est=var, est.lb=var.ci[1], est.ub=var.ci[2]), 
    data.frame(model=model, metric="Bias",est=bias, est.lb=bias.ci[1], est.ub=bias.ci[2]), 
    data.frame(model=model, metric="MSE",est=mse, est.lb=mse.ci[1], est.ub=mse.ci[2]))
  rownames(res) <- NULL
  return(res)
}


#Replicate sampled dataset:
data <- readRDS(here("simulated data/clean_data.RDS"))
head(data)

D <- DAG.empty() + 
  node("CVD", distr="rcat.b1", probs = c(0.5, 0.25, 0.25)) +
  node("educ", distr="rcat.b1", probs = c(0.52,0.308,0.094,0.078)) +
  node("diab.dur", distr="rcat.b1", probs = c(0,0,0,0,0,.3,.2,.15,.1,.1,.05,.05,.03,.02)) +
  node("income", distr = "rnorm_trunc", mean = 100, sd = 25) +
  node("bmi", distr = "rnorm_trunc", mean = 25, sd = 10) +
  node("U.W", distr = "rnorm", mean = 0, sd = 1) +
  node("U.A", distr = "rnorm", mean = 0, sd = 1) +
  node("U.Y", distr = "rnorm", mean = 0, sd = 5) +
  node("U.C", distr = "rnorm", mean = 0, sd = 1) +
  node("A", t=0, distr="rbern", prob=plogis(.08*CVD + .005*diab.dur + 0.1*educ + .0025*bmi + U.A -.4)) + 
  node("La", t=0, distr="rbern", prob=plogis(-2 - 0.3*CVD - 0.5*A[t])) +
  node("Lb", t=0, distr="rbern", prob=plogis(-2 - 0.3*CVD - 0.5*A[t])) +
  node("A", t=1:9, distr="rbern", prob=plogis(A[t-1] )) +
  node("La", t=1:9, distr="rbern", prob=plogis(-3 - 0.3*CVD + diab.dur +bmi + educ + -0.5*A[t] + 1.5*La[t-1])) +
  node("Lb", t=1:9, distr="rbern", prob=plogis(-3 - 0.3*CVD + diab.dur +bmi + educ + -0.5*A[t] + 1.5*Lb[t-1])) +
  node("Y", t=0:9, distr="rbern", prob=plogis( -3 - 1.2*La[t] - 1.2*Lb[t] - 0.001*income*educ + diab.dur/10 * 0.1*CVD - 0.4*A[t] +0.01*bmi -0.02*A[t]*bmi + U.Y), EFU=TRUE)


D <- set.DAG(D)
saveRDS(D, file=here("simulated data/simulated_dag.RDS"))

dat <- sim(D,n=5000, LTCF="Y", verbose=T, rndseed=12345) #Do I need to fill in censoring with LTCF argument?
head(dat)

summary(dat$bmi)
table(dat$CVD)
table(dat$diab.dur)
table(dat$A_0)
table(dat$A_9)
table(dat$La_0)
table(dat$La_9)
table(dat$Lb_0)
table(dat$Lb_9)

table(dat$Y_0)
table(dat$Y_9)
table(dat$A_9==1, dat$Y_9)
tab<- table(dat$A_9==1, dat$Y_9)
(tab[1,1]*tab[2,2])/(tab[1,2]*tab[2,1])




# act_t0_theta <- node("A",t=0:9, distr="rbern", prob=ifelse(theta==1,1,0))
# D <- D + action("A_th0", nodes=c(act_t0_theta), theta=0)
# D <- D + action("A_th1", nodes=c(act_t0_theta), theta=1)
# D <- set.targetE(D, outcome="Y", t=0:9, param="A_th1 / A_th0")
# eval.target(D, n=5000)
# 
# D <- set.targetE(D, outcome="Y", t=9, param="A_th1 / A_th0")
# trueRR <- eval.target(D, n=5000)$res



#simulation parameters
N_sim<-100
#ndata<-5000 #could vary this with runif()
set.seed(12345)
ndata<- ceiling(runif(N_sim)*10)^4 + ceiling(runif(N_sim)*100)^2 +100
summary(ndata)

sim_res_glm <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res=try(glm(Y_9~A_0, data=dat, family = binomial(link="log")))
  sim_res_glm$est[i] <- try(exp(res$coefficients)[2])
  sim_res_glm$var[i] <- try(vcov(res)[2,2])
}
sim_res_glm
saveRDS(sim_res_glm,  file=here("results/simulation_results_glm.RDS"))

  

sim_res_glm_adj <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res = try(glm(Y_9~A_0 + CVD + diab.dur +bmi + educ + income, data=dat, family = binomial(link="log")))
  sim_res_glm_adj$est[i] <- try(exp(res$coefficients)[2])
  sim_res_glm_adj$var[i] <- try(vcov(res)[2,2])
  }
sim_res_glm_adj
saveRDS(sim_res_glm_adj,  file=here("results/simulation_results_glm_adj.RDS"))

sim_res_tmle_glm <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res <- try(tmle(Y=dat$Y_9, A=dat$A_0, W=subset(dat, select= c(CVD:bmi)), family="binomial",Q.SL.library = c("SL.glm"),g.SL.library = c("SL.glm")))
  sim_res_tmle_glm$est[i] <- try(res$estimates$RR$psi)
  sim_res_tmle_glm$var[i] <- try(res$estimates$RR$var.log.psi)
}
sim_res_tmle_glm
saveRDS(sim_res_tmle_glm,  file=here("results/simulation_results_tmle_glm.RDS"))

sim_res_tmle <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
lib = c("SL.glm","SL.glmnet","SL.glm.interaction","SL.randomForest")
glib = c("SL.glm","SL.glmnet")
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res <- try(tmle(Y=dat$Y_9, A=dat$A_0, W=subset(dat, select= c(CVD:bmi)), family="binomial",Q.SL.library = lib,g.SL.library = glib))
  sim_res_tmle$est[i] <- try(res$estimates$RR$psi)
  sim_res_tmle$var[i] <- try(res$estimates$RR$var.log.psi)
}
sim_res_tmle
saveRDS(sim_res_tmle,  file=here("results/simulation_results_tmle.RDS"))



i<-1
dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
colnames(dat)
Anodes <- c("A_1","A_2","A_3","A_4","A_5","A_6","A_7","A_8","A_9")
Ynodes <- c("Y_1","Y_2","Y_3","Y_4","Y_5","Y_6","Y_7","Y_8","Y_9")
Cnodes <-NULL
Lnodes <- c("La_1","Lb_1","L1c","L1d","L1f","L1g","L1h",
            "La_2","Lb_2",
            "La_3","Lb_3",
            "La_4","Lb_4",
            "La_5","Lb_5",
            "La_6","Lb_6",
            "La_7","Lb_7",
            "La_8","Lb_8",
            "La_9","Lb_9")




lib = c("SL.glm")
abar <- list(a=rep(1,(length(Anodes))), b=rep(0,(length(Anodes))))

sim_res_ltmle <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  subset <- dat %>% rename(
    L1c=CVD, L1d=educ, L1f=diab.dur, L1g=income,  L1h=bmi,
  ) %>%
    select(
      "La_1","Lb_1","L1c","L1d","L1f","L1g","L1h","A_1",
      "La_2","Lb_2","A_2",
      "La_3","Lb_3","A_3",
      "La_4","Lb_4","A_4",
      "La_5","Lb_5","A_5",
      "La_6","Lb_6","A_6",
      "La_7","Lb_7","A_7",
      "La_8","Lb_8","A_8",
      "La_9","Lb_9","A_9",Ynodes)
  res=NULL
  
  set.seed(12345)
  abar <- list(a=rep(1,(length(Anodes))), b=rep(0,(length(Anodes))))
  result <- ltmle(subset, Anodes = Anodes, Ynodes = Ynodes,
                  Cnodes=Cnodes, Lnodes=Lnodes, abar = abar,
                  survivalOutcome=F,SL.library=lib,variance.method = "ic")
  
  res<-summary(result)
  sim_res_ltmle$est[i] <- try(res$effect.measures$RR$estimate)
  sim_res_ltmle$var[i] <- try(res$effect.measures$RR$std.dev^2)
}
sim_res_ltmle
saveRDS(sim_res_ltmle,  file=here("results/simulation_results_tmle.RDS"))


save(
  sim_res_glm,  
  sim_res_glm_adj, 
  sim_res_tmle_glm, 
  sim_res_tmle, 
  sim_res_ltmle,
  file=here("results/simulation_results.Rdata")
)


