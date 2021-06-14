
rm(list=ls())
library(ltmle)
library(here)
library(dplyr)
library(cowplot)
load(here("results/sampled_data_estimates.Rdata"))

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

res=est_tmle
param="ATE"
get_est <- function(res, param="ATE", type="tmle"){
  est_unadj= summary(res$fit_unadj,type)$effect.measures[[param]]$estimate
  ci1_unadj= summary(res$fit_unadj,type)$effect.measures[[param]]$CI[1]
  ci2_unadj= summary(res$fit_unadj,type)$effect.measures[[param]]$CI[2]
  est_adj= summary(res$fit_adj,type)$effect.measures[[param]]$estimate
  ci1_adj= summary(res$fit_adj,type)$effect.measures[[param]]$CI[1]
  ci2_adj= summary(res$fit_adj,type)$effect.measures[[param]]$CI[2]
  res_unadj = data.frame(type="unadj", model=type, est=est_unadj, ci1=ci1_unadj, ci2=ci2_unadj, param=param)
  res_adj = data.frame(type="adj",  model=type, est=est_adj, ci1=ci1_adj, ci2=ci2_adj, param=param)
  res=bind_rows(res_unadj, res_adj)
  return(res)
}

res <- bind_rows(
  get_est(est_gcomp, type="gcomp"),
  get_est(est_iptw, type="iptw"),
  get_est(est_tmle),
  get_est(est_gcomp, param="RR", type="gcomp"),
  get_est(est_iptw, param="RR", type="iptw"),
  get_est(est_tmle, param="RR"))

#Clean up for plotting
res <- res %>%
  mutate(
    model=case_when(model=="gcomp" ~"G-comp.",
                    model=="iptw" ~"IPTW",
                    model=="tmle" ~"TMLE"),
    model=factor(model, levels=c("TMLE","IPTW","G-comp."))
  )
  


pATE = ggplot(data=res[res$param=="ATE"& res$type=="adj",],
           aes(x = model,y = est, ymin = ci1, ymax = ci2 ))+
  geom_pointrange(aes(col=model))+
  geom_hline(yintercept =0, linetype=2)+
  xlab('Estimator')+ ylab("Average Treatment Effect (95% CI)")+coord_flip()+
  geom_errorbar(aes(col=model),width=0.5,cex=1) +
  #facet_wrap(~model,nrow=9,scales = "free_y") +
  scale_color_manual(values = tableau10)
pATE

pRR = ggplot(data=res[res$param=="RR" & res$type=="adj",],
           aes(x = model,y = est, ymin = ci1, ymax = ci2 ))+
  geom_pointrange(aes(col=model))+
  geom_hline(yintercept =1, linetype=2)+
  xlab('')+ ylab("Relative risk (95% CI)")+coord_flip()+
  scale_y_continuous(trans='log10')+
  geom_errorbar(aes(col=model),width=0.5,cex=1) +
  scale_color_manual(values = tableau10) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pRR

p <- plot_grid(pATE, pRR, rel_widths=c(0.57,0.4))
p

ggsave(p, file=here("figures/synthetic_data_ests.png"),  width = 7.5,     height = 4.26,)