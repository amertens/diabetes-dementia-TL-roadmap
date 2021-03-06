

rm(list=ls())
library(here)
library(tidyverse)
library(simcausal)
library(simhelpers)

#Truncated continious
rnorm_trunc <- function(n, mean, sd, minval = 17){
  out <- rnorm(n = n, mean = mean, sd = sd)
  minval <- minval[1]
  out[out < minval] <- minval
  out
}


#Load  DAG and simulation results
load(here("results/simulation_results_RD.Rdata"))
D <- readRDS(here("simulated data/simulated_dag_RD.RDS"))

sim_res_glm <- sim_res_glm %>% mutate(est=as.numeric(est), var=as.numeric(var))
sim_res_glm_adj  <- sim_res_glm_adj %>% mutate(est=as.numeric(est), var=as.numeric(var))
sim_res_tmle_glm <- sim_res_tmle_glm %>% mutate(est=as.numeric(est), var=as.numeric(var))
sim_res_tmle <- sim_res_tmle %>% mutate(est=as.numeric(est), var=as.numeric(var))
sim_res_ltmle <- sim_res_ltmle %>% mutate(est=as.numeric(est), var=as.numeric(var))


#Get true RD
act_t0_theta <- node("A",t=0:3, distr="rbern", prob=ifelse(theta==1,1,0))
D <- D + action("A_th0", nodes=c(act_t0_theta), theta=0)
D <- D + action("A_th1", nodes=c(act_t0_theta), theta=1)
D <- set.targetE(D, outcome="Y", t=0:3, param="A_th1 - A_th0")
eval <- eval.target(D, n=500000, rndseed=12345)
eval$res

trueRD <- eval$res[4]
trueRD


#combine simulation results
sim_res <- bind_rows(
  data.frame(sim_res_glm, trueRD=trueRD, model="Unadjusted GLM"),
  data.frame(sim_res_glm_adj, trueRD=trueRD, model="Adjusted GLM"),
  data.frame(sim_res_tmle_glm, trueRD=trueRD, model="TMLE fit with GLM"),
  data.frame(sim_res_tmle, trueRD=trueRD, model="TMLE fit with SL"),
  data.frame(sim_res_ltmle, trueRD=trueRD, model="longitudinal TMLE")
)

#use simhelpers package to estimate simulation performances
sim_res <- sim_res %>% mutate(
  model=factor(model, levels=rev(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL", "longitudinal TMLE"))),
  metric=factor(metric, levels=c("Variance", "Bias", "MSE"))
) %>% arrange(model, metric)

sim_res %>%
  group_by(model) %>% # grouping 
  do(calc_absolute(., estimates = est, true_param = trueRD)) %>%
  mutate(model=factor(model, levels=(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL")))) %>%
  arrange(model) %>%
  knitr::kable()

sim_res %>%
  group_by(model) %>% # grouping 
  do(calc_relative(., estimates = est, true_param = trueRD)) %>%
  mutate(model=factor(model, levels=(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL")))) %>%
  arrange(model) %>%
  knitr::kable()


sim_res <- sim_res %>%
  group_by(model) %>% # grouping 
  do(calc_absolute(., estimates = est, true_param = trueRD)) %>%
  gather(bias:rmse_mcse,  key = metric, value = est) %>%
  filter(metric!="rmse" & metric!="rmse_mcse")
sim_res_est <- sim_res %>% filter(!grepl("_mcse",metric))
sim_res_var <- sim_res %>% filter(grepl("_mcse",metric)) %>% mutate(metric=gsub("_mcse","",metric)) %>% rename(se=est)
sim_res <- merge(sim_res_est, sim_res_var, by=c("model","metric","K")) %>%
  mutate(est.lb=est-1.96*se, est.ub=est+1.96*se)
sim_res

#Plot performance metrics
ggplot(sim_res, aes(x=model, y=est)) + 
  geom_point() + geom_linerange(aes(ymin=est.lb, ymax=est.ub)) +
  geom_hline(yintercept=0) +
  facet_wrap(~metric, scales="free") + coord_flip() +
  theme_bw()



#to do:
#make sure adjusted glm is better than glm
#then try tmle and make sure it is better than adj glm
#then try ltmle
#compare true RR from checking just at Y9 vs Y0:9. 
#maybe just make a table of metrics and Se's







