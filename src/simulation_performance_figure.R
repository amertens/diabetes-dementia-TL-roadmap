
rm(list=ls())
library(ltmle)
library(here)
library(tidyverse)
library(cowplot)
load(here("results/perf_metrics1.Rdata"))

# 
# load(here("results/simulation_results_RD.Rdata"))
# 
# 
# 
# sim_res_glm <- sim_res_glm %>% mutate(est=as.numeric(est), var=as.numeric(var))
# sim_res_glm_adj <- sim_res_glm_adj %>% mutate(est=as.numeric(est), var=as.numeric(var))
# sim_res_tmle_glm <- sim_res_tmle_glm %>% mutate(est=as.numeric(est), var=as.numeric(var))
# sim_res_tmle <- sim_res_tmle %>% mutate(est=as.numeric(est), var=as.numeric(var))
# sim_res_ltmle <- sim_res_ltmle %>% mutate(est=as.numeric(est), var=as.numeric(var))
# 
# 
# 
# #combine simulation results
# sim_res <- bind_rows(
#   data.frame(sim_res_glm, trueRD=trueRD, model="Unadjusted GLM"),
#   data.frame(sim_res_glm_adj, trueRD=trueRD, model="Adjusted GLM"),
#   data.frame(sim_res_tmle_glm, trueRD=trueRD, model="UTMLE fit with GLM"),
#   data.frame(sim_res_tmle, trueRD=trueRD, model="TMLE fit with SL"),
#   data.frame(sim_res_ltmle, trueRD=trueRD, model="longitudinal TMLE")
# ) 
# 
# #use simhelpers package to estimate simulation performances
# sim_res <- sim_res %>% mutate(
#   model=factor(model, levels=rev(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL", "longitudinal TMLE")))#,
#   #metric=factor(metric, levels=c("Variance", "Bias", "MSE"))
# ) %>% arrange(model)
# 
# sim_res %>%
#   group_by(model) %>% # grouping 
#   do(calc_absolute(., estimates = est, true_param = trueRD)) %>%
#   mutate(model=factor(model, levels=(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL", "longitudinal TMLE")))) %>%
#   arrange(model) %>%
#   knitr::kable()
# 
# sim_res %>%
#   group_by(model) %>% # grouping 
#   do(calc_relative(., estimates = est, true_param = trueRD)) %>%
#   mutate(model=factor(model, levels=(c("Unadjusted GLM", "Adjusted GLM", "TMLE fit with GLM", "TMLE fit with SL")))) %>%
#   arrange(model) %>%
#   knitr::kable()
# 
# 
# plot_res <- sim_res %>%
#   group_by(model) %>% # grouping 
#   do(calc_absolute(., estimates = est, true_param = trueRD)) %>%
#   gather(bias:rmse_mcse,  key = metric, value = est) %>%
#   filter(metric!="rmse" & metric!="rmse_mcse")
# plot_res_est <- plot_res %>% filter(!grepl("_mcse",metric))
# plot_res_var <- plot_res %>% filter(grepl("_mcse",metric)) %>% mutate(metric=gsub("_mcse","",metric)) %>% rename(se=est)
# plot_res <- merge(plot_res_est, plot_res_var, by=c("model","metric","K")) %>%
#   mutate(est.lb=est-1.96*se, est.ub=est+1.96*se)
# plot_res
# 
# #Plot performance metrics
# ggplot(plot_res, aes(x=model, y=est)) + 
#   geom_point() + geom_linerange(aes(ymin=est.lb, ymax=est.ub)) +
#   geom_hline(yintercept=0) +
#   facet_wrap(~metric, scales="free") + coord_flip() +
#   theme_bw()




theme_ki <- function() {
  theme_bw() %+replace%
    theme(
      strip.background = element_blank(),
      legend.position="none",
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size=14),
      axis.title = element_text(size=12),
      axis.text.y = element_text(size=10),
      axis.text.x = element_text(size=10, angle = 0, hjust = 0.5, vjust=.1)
    )
}

#hbgdki pallets
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
tableau11 <- c("Black","#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

theme_set(theme_ki())


sim_res
plot_res
sim_res %>% group_by(model) %>%
  summarise(bias=mean(est-trueRD), var=mean(var)) %>%
  mutate(mse=bias^2-var)

#Clean up for plotting
plot_res <- plot_res %>%
  mutate(
    metric=case_when(metric=="bias" ~"Bias",
                     metric=="var" ~"Variance",
                     metric=="mse" ~"MSE"),
    metric=factor(metric, levels=c("Bias","Variance","MSE"))
  )


p <-ggplot(plot_res, aes(x=model, y=est, color=model)) + 
  geom_point() + geom_linerange(aes(ymin=est.lb, ymax=est.ub)) +
  scale_color_manual(values = tableau10) + 
  geom_hline(yintercept=0) + ylab("") + xlab("") +
  facet_wrap(~metric, scales="free") + coord_flip() 
p




ggsave(p, file=here("figures/simulation_performance.png"),  width = 10,     height = 4)