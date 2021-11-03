getwd()  # under folder diabetes-dementia-TL-roadmap
pkg_dir <- "./simulation_ZW/DK_trip_2021/pkgs/ltmle/cvSL_snow_skipspeedglm_20211005/ltmle-master/R/"  # load ltmle with snow parallel SL; skip speedglm runs
# pkg_dir <- "./simulation_ZW/DK_trip_2021/pkgs/ltmle/cvSL_skipspeedglm_20210915/ltmle-master/R/"  # load ltmle with multicore parallel SL; skip speedglm runs
data_path <- "./simulation_ZW/DK_trip_2021/dt_use_backup_MSM.rds"  # example data


#-------------------------------------------------------------------
# Simulate data
#-------------------------------------------------------------------

library(dplyr)
library(data.table)
library(pryr)

# load simulated data
{
  logit <- function(x) log(x/(1-x))
  expit <- function(x) 1/(1+exp(-x))
  coef_intercept <- 3.5
  coef_a <- -0.3  # marginal time coefficient for survival prob
  coef_b <- 0.15  # marginal dose coefficient for survival prob
  
  dt_use_backup <- readRDS(data_path)
  
  set.seed(123)
  dt_use <- dt_use_backup[sample(nrow(dt_use_backup), 20000, T), ]  # target sample size
  K <- 10  # total time points
  dt_use[, first_date_2nd_line := NULL]  # remove index dates; equivalence with censoring process
  node_names <- names(dt_use)
  
  # fake dual traetment nodes
  set.seed(123)
  for (i in 1:10) {
    dt_use[, paste0("A2_", i) := get(paste0("A1_", i))[sample(nrow(dt_use))]]
    setcolorder(dt_use, c(colnames(dt_use)[1:grep(paste0("^A1_", i, "$"), colnames(dt_use))], 
                          paste0("A2_", i))
    )
  }
  
  # fake competing risk; only happens when uncensored and dementia-free with 0.05 probability
  # once dead, all future Y2_t become 1 and Y_t become 0
  set.seed(123)
  for (i in 1:11) {
    # transform all non-dementia to death
    temp_Y <- dt_use[, paste0("Y_", i), with = F]
    temp_Y <- temp_Y + 1
    temp_Y[temp_Y==2] <- 0
    # only keep 0.05 death
    temp_Y[sample(which(temp_Y == 1), round(length(which(temp_Y == 1)) * 0.95))] <- 0
    
    dt_use[, paste0("Y2_", i) := temp_Y]
    setcolorder(dt_use, c(colnames(dt_use)[1:grep(paste0("^Y_", i, "$"), colnames(dt_use))], 
                          paste0("Y2_", i))
    )
  }
  for (i in 1:11) {
    temp_Y <- dt_use[, paste0("Y2_", i), with = F]
    for (j in i:11) {
      dt_use[which(temp_Y == 1), (paste0("Y2_", j)) := 1] 
    }
  }
  # if any dementia, remove death
  dt_use[which(apply(dt_use[, paste0("Y_", 1:11)] == 1, 1, any)), 
         (paste0("Y2_", 1:11)) := ifelse(get((paste0("Y2_", 1:11))) %>% is.na, NA, 0)
  ]
  # # if dementia, remain that way
  # for (i in 1:11) {
  #   temp_Y <- dt_use[, paste0("Y2_", i), with = F]
  #   for (j in i:11) {
  #     dt_use[which(temp_Y == 1), (paste0("Y2_", j)) := 1] 
  #   }
  # }
  
  
  # if censored, remove event
  dt_use[which(apply(dt_use[, names(dt_use)[grep("^C", names(dt_use))], with = F] == 0, 1, any)),
         (paste0("Y_", 1:11)) := 0
  ]
  dt_use[which(apply(dt_use[, names(dt_use)[grep("^C", names(dt_use))], with = F] == 0, 1, any)),
         (paste0("Y2_", 1:11)) := 0
  ]
  
  
  C.nodes.index <- grep("^C_", names(dt_use))
  for (i in 1:10) {
    temp_C <- dt_use[, paste0("C_", i), with = F]
    cur.C.index <- grep(paste0("^C_", i, "$"), names(dt_use))
    for (j in (cur.C.index+1):ncol(dt_use)) {
      if (j %in% C.nodes.index) {} else {
        dt_use[which(temp_C == 0), (names(dt_use)[j]) := NA] 
      }
    }
  }
  dt_use[which(apply(dt_use[, paste0("Y_", 1:11)] == 1, 1, any)), 
         (paste0("Y2_", 1:11)) := ifelse(get((paste0("Y2_", 1:11))) %>% is.na, NA, 0)
  ]
  
  
  d <- as.data.frame(dt_use)
}


for(i in 1:N_time){
  d[[paste0(censor_var,i)]] <- ifelse(d[[paste0(censor_var,i)]]==1,"uncensored","censored")
}

dfull<-d
d<-dfull

#-------------------------------------------------------------------
#get summary tables 
#-------------------------------------------------------------------

head(d)
colnames(d)

#Set to number of N-nodes to calc over
N_time <- 10

#Need to set vars we want to tabulate to NA's for censored variables
group_var <- "A1_"
censor_var <- "C_"
for(i in 1:N_time){
   d[[paste0(group_var,i)]][d[[paste0(censor_var,i)]] == "censored"] <- NA
}
head(d)


tab_treatment_regime_glp1 <- d %>% 
  mutate(Total_N = n()) %>% 
  #group by all longitudinal nodes
  group_by(A1_1, A1_2, A1_3, A1_4, A1_5, A1_6, A1_7, A1_8, A1_9, A1_10) %>% 
  #count occurrences and number of events by the end of followup (update Y_10 to last outcome node)
  summarise(N=n(), Total_N=Total_N[1]) %>% ungroup() %>%
  #make sure no 1's or 2's
  mutate(N=case_when(N==1 ~ 3, N==2 ~ 3, (N==0|N>2) ~ as.numeric(N)),
         weights=N/Total_N)


write.csv(tab_treatment_regime_glp1, file=here::here("/reports/tab_treatment_regime_glp1.csv"))


#Need to set vars we want to tabulate to NA's for censored variables
group_var <- "Y_"
censor_var <- "C_"
for(i in 1:N_time){
  d[[paste0(group_var,i)]][d[[paste0(censor_var,i)]] == "censored"] <- NA
}
head(d)

tab_treatment_regime_dementia <- d %>% 
  mutate(Total_N = n()) %>% 
  #group by all longitudinal nodes
  group_by(Y_1, Y_2, Y_3, Y_4, Y_5, Y_6, Y_7, Y_8, Y_9, Y_10) %>% 
  #count occurrences and number of events by the end of followup (update Y_10 to last outcome node)
  summarise(N=n(), Total_N=Total_N[1]) %>% ungroup() %>%
  #make sure no 1's or 2's
  mutate(N=case_when(N==1 ~ 3, N==2 ~ 3, (N==0|N>2) ~ as.numeric(N)),
         weights=N/Total_N)


write.csv(tab_treatment_regime_glp1, file=here::here("/reports/tab_treatment_regime_dementia.csv"))


