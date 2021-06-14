

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

#function for Truncated continious distribution
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

#example simulation of the data and crosstabs
dat <- sim(D,n=5000, LTCF="Y", verbose=T, rndseed=12345) #Do I need to fill in censoring with LTCF argument?
head(dat)

table(dat$CVD)
table(dat$A_0)
table(dat$A_2)
table(dat$La_0)
table(dat$La_2)


table(dat$Y_1)
table(dat$Y_3)
table(dat$A_2==1, dat$Y_3)
tab<- table(dat$A_2==1, dat$Y_3)

#Crude OR and RD
(tab[1,1]*tab[2,2])/(tab[1,2]*tab[2,1])
(tab[1,1]/(tab[1,1]+tab[1,2])) - (tab[2,1]/(tab[2,1]+tab[2,2]))

  
#Nerissa: here is the code to estimate the true RD, but I fear I am setting it up 
#wrong. I set actions, but the eval.target command gives me the warning
#"no actions specified, sampling full data for ALL actions from the DAG"
# and the RD seems very far off the simulated results. 
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

#-------------------------------------
# TMLE fit with SL
#-------------------------------------

sim_res_tmle <- data.frame(est=rep(NA, N_sim), var=rep(NA, N_sim))
#lib = c("SL.glm","SL.glmnet","SL.glm.interaction","SL.randomForest")
lib = c("SL.glm","SL.glmnet")
glib = c("SL.glm","SL.glmnet")
for(i in 1:N_sim){
  dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
  res=NULL
  res <- try(tmle(Y=dat$Y_3, A=dat$A_0, W=subset(dat, select= c(CVD,La_0,La_1, La_2)), family="gaussian",Q.SL.library = lib,g.SL.library = glib))
  sim_res_tmle$est[i] <- try(res$estimates$ATE$psi)
  sim_res_tmle$var[i] <- try(res$estimates$ATE$var.psi)
}
sim_res_tmle


#-------------------------------------
# LTMLE 
#-------------------------------------
Anodes <- c("A_0","A_1","A_2")
Ynodes <- c("Y_1","Y_2","Y_3")
Cnodes <-NULL
Lnodes <- c("La_0","L0c","La_1","La_2")



lib = c("SL.glm","SL.glmnet")
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
  sim_res_ltmle$var[i] <- try((res$effect.measures$ATE$std.dev^2)/nrow(subset))
}
sim_res_ltmle



sim_res_glm <- sim_res_glm %>% mutate(est=as.numeric(est), var=as.numeric(var))
sim_res_tmle <- sim_res_tmle %>% mutate(est=as.numeric(est), var=as.numeric(var))
sim_res_ltmle <- sim_res_ltmle %>% mutate(est=as.numeric(est), var=as.numeric(var))


#Get true RD
act_t0_theta <- node("A",t=0:2, distr="rbern", prob=ifelse(theta==1,1,0))
D <- D + action("A_th0", nodes=c(act_t0_theta), theta=0)
D <- D + action("A_th1", nodes=c(act_t0_theta), theta=1)
D <- set.targetE(D, outcome="Y", t=0:3, param="A_th1 - A_th0")
eval <- eval.target(D, n=5000, rndseed=12345)
eval$res

D <- set.targetE(D, outcome="Y", t=3, param="A_th1 - A_th0")
eval <- eval.target(D, n=5000, rndseed=12345)
eval$res

# trueRD <- eval$res[3]
# trueRD

mean(sim_res_glm$est)
mean(sim_res_tmle$est)
mean(sim_res_ltmle$est)
mean(sim_res_glm$var)
mean(sim_res_tmle$var)
mean(sim_res_ltmle$var)


#combine simulation results
sim_res <- bind_rows(
  data.frame(sim_res_glm, trueRD=trueRD, model="Unadjusted GLM"),
  data.frame(sim_res_tmle, trueRD=trueRD, model="TMLE fit with SL"),
  data.frame(sim_res_ltmle, trueRD=trueRD, model="longitudinal TMLE")
) 

#use simhelpers package to estimate simulation performances
sim_res <- sim_res %>% mutate(
  model=factor(model, levels=rev(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL", "longitudinal TMLE")))#,
  #metric=factor(metric, levels=c("Variance", "Bias", "MSE"))
) %>% arrange(model)

sim_res %>%
  group_by(model) %>% # grouping 
  do(calc_absolute(., estimates = est, true_param = trueRD)) %>%
  mutate(model=factor(model, levels=(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL", "longitudinal TMLE")))) %>%
  arrange(model) %>%
  knitr::kable()

sim_res %>%
  group_by(model) %>% # grouping 
  do(calc_relative(., estimates = est, true_param = trueRD)) %>%
  mutate(model=factor(model, levels=(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL")))) %>%
  arrange(model) %>%
  knitr::kable()


plot_res <- sim_res %>%
  group_by(model) %>% # grouping 
  do(calc_absolute(., estimates = est, true_param = trueRD)) %>%
  gather(bias:rmse_mcse,  key = metric, value = est) %>%
  filter(metric!="rmse" & metric!="rmse_mcse")
plot_res_est <- plot_res %>% filter(!grepl("_mcse",metric))
plot_res_var <- plot_res %>% filter(grepl("_mcse",metric)) %>% mutate(metric=gsub("_mcse","",metric)) %>% rename(se=est)
plot_res <- merge(plot_res_est, plot_res_var, by=c("model","metric","K")) %>%
  mutate(est.lb=est-1.96*se, est.ub=est+1.96*se)
plot_res

#Plot performance metrics
ggplot(plot_res, aes(x=model, y=est)) + 
  geom_point() + geom_linerange(aes(ymin=est.lb, ymax=est.ub)) +
  geom_hline(yintercept=0) +
  facet_wrap(~metric, scales="free") + coord_flip() +
  theme_bw()

save(sim_res, plot_res, file=here("results/perf_metrics1.Rdata"))

#to do:
#make sure adjusted glm is better than glm
#then try tmle and make sure it is better than adj glm
#then try ltmle
#compare true RR from checking just at Y9 vs Y0:9. 
#maybe just make a table of metrics and Se's












# lib = c("SL.glm")
# abar <- list(a=rep(1,(length(Anodes))), b=rep(0,(length(Anodes))))
# i<-1
# 
#   dat <- sim(D,n=ndata[i], LTCF="Y", verbose=F, rndseed=i) 
#   subset <- dat %>% rename(
#     L1c=CVD, L1d=educ, L1f=diab.dur, L1g=income,  L1h=bmi,
#   ) %>%
#     select(
#       "La_1","Lb_1","L1c","L1d","L1f","L1g","L1h","A_1",
#       "La_2","Lb_2","A_2",
#       "La_3","Lb_3","A_3",Ynodes)
#   res=NULL
#   
#   set.seed(12345)
#   abar <- list(a=rep(1,(length(Anodes))), b=rep(0,(length(Anodes))))
#   result <- ltmle(subset, Anodes = Anodes, Ynodes = Ynodes,
#                   Cnodes=Cnodes, Lnodes=Lnodes, abar = abar,
#                   survivalOutcome=T,SL.library=lib,variance.method = "ic")
#   
#   res<-summary(result)
#   res

