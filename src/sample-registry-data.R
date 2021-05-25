
rm(list=ls())
library(here)
library(tidyverse)
library(holodeck)
rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))

#-------------------------------- 
# Load data summaries
#-------------------------------- 
glp1 <- read.csv(here::here("data summaries/glp1-regimes-npar-table.csv"), header = FALSE)
sglt2 <- read.csv(here::here("data summaries/sglt2-regimes-npar-table.csv"), header = FALSE)

colnames(glp1) <- c("age","sex","Aglp_1","Aglp_2","Aglp_3","Aglp_4","Aglp_5","Aglp_6","Aglp_7","Aglp_8","Aglp_9","death","dementia","stroke","ami","weight")
colnames(sglt2) <- c("age","sex","Asglt_1","Asglt_2","Asglt_3","Asglt_4","Asglt_5","Asglt_6","Asglt_7","Asglt_8","Asglt_9","death","dementia","stroke","ami","weight")

glp1$A <- "glp1"
sglt2$A <- "sglt2"

df <- bind_rows(glp1,sglt2) %>% mutate(
  dementia=as.numeric(dementia),
  death=as.numeric(death),
  stroke=as.numeric(stroke),
  ami=as.numeric(ami),
  any_event = ifelse(!is.na(death) | !is.na(dementia) | !is.na(stroke) | !is.na(ami), 1, 0),
  CVOT = ifelse( !is.na(stroke) | !is.na(ami), 1, 0),
  any_death=ifelse(is.na(death),0,1),
  any_dementia=ifelse(is.na(dementia),0,1)
)
df <- df %>% mutate(tot_weight=sum(weight))


#Add SGLT2 dementia based on frequency in 
#file:///C:/Users/andre/Dropbox/CTL/Joint%20Initiative%20for%20Causal%20Inference%20Novo%20UCB%20UCPH/Active%20projects/Dementia/FinalReport-2019-04-29.html
62 /(40572+4557)
df <- bind_rows(
    data.frame(age=4, sex=0, any_dementia=1, dementia=3, 
               #weight=62 /(40572+4557) * df$tot_weight[1],
               weight=3 /(40572+4557) * df$tot_weight[1],
    Asglt_1=1, Asglt_2=1, Asglt_3=1, Asglt_4=1, Asglt_5=1, Asglt_6=1,  Asglt_7=1, Asglt_8=1, Asglt_9=1,
    A="sglt2", any_event=1, CVOT=0, total_usage=1, any_drug=1), df
    )

df[is.na(df)] <- 0
table(is.na(df))

df <- df %>% mutate(num_sglt=Asglt_1 +Asglt_2 +Asglt_3 +Asglt_4 +Asglt_5 +Asglt_6 +Asglt_7 +Asglt_8 +Asglt_9,
                    num_glp=Aglp_1 +Aglp_2 +Aglp_3 +Aglp_4 +Aglp_5 +Aglp_6 +Aglp_7 +Aglp_8 +Aglp_9)

table(df$num_sglt, df$dementia)
table(df$num_glp, df$dementia)
df %>% filter(dementia!=0, num_sglt!=0|num_glp!=0) %>% group_by(A) %>% 
  summarise(sum(weight))

#-------------------------------- 
# Make dataset by sampling relative 
# frequency of each covariate/treatment pattern
#-------------------------------- 

#Dataset size
ndf <- 45000
#Sample with replacement with the probability of the observed frequency
set.seed(12)
d <- df[sample(nrow(df), ndf, replace = TRUE, prob = df$weight), ]

table(d$num_sglt, d$dementia)
table(d$num_glp, d$dementia)

table(d$A, d$dementia)
table(d$A, d$any_dementia)
table(d$A, d$any_death)
table(d$A, d$CVOT)
table(d$A, d$any_event)


d$any_dementia[is.na(d$any_dementia)] <- 0
d$any_death[is.na(d$any_death)] <- 0
d$CVOT[is.na(d$CVOT)] <- 0
d$any_event[is.na(d$any_event)] <- 0
d$stroke[is.na(d$stroke)] <- 0
d$ami[is.na(d$ami)] <- 0

# #recode T1=T9 as 1==GLP1
# d <- d %>% mutate(
#   T1 =ifelse(A=="glp1",1,0),
#   T2 =ifelse(A=="glp1",1,0),
#   T3 =ifelse(A=="glp1",1,0),
#   T4 =ifelse(A=="glp1",1,0),
#   T5 =ifelse(A=="glp1",1,0),
#   T6 =ifelse(A=="glp1",1,0),
#   T7 =ifelse(A=="glp1",1,0),
#   T8 =ifelse(A=="glp1",1,0),
#   T9 =ifelse(A=="glp1",1,0)
# )



#Add baseline confounders by outcome status
# https://cran.r-project.org/web/packages/holodeck/vignettes/simulating-data.html

set.seed(12345)
d$edu <- NA
d$edu[d$any_dementia==0] <- sample(rep( c("Basic","Medium","Advanced","Unknown"), nrow(d)*c(0.52,0.308,0.094,0.078)), length(d$edu[d$any_dementia==0]))
d$edu[d$any_dementia==1] <- sample(rep( c("Basic","Medium","Advanced","Unknown"), nrow(d)*c(0.512,0.309,0.085,0.095)), length(d$edu[d$any_dementia==1]))
d$edu <- factor(d$edu)
table(d$edu)

d <- d %>%
  group_by(any_dementia) %>% 
  sim_discr(n_vars = 3, var = 1, cov = 0.5, group_means = c(10, 15),  seed=12345) %>%
  rename(income=V1, diabduration=V2, start_year=V3) %>%
  mutate(start_year=ntile(start_year,6))


# #simulate missingness
# df2 <-
#   df1 %>% 
#   sim_missing(prop = 0.1)
# df2

#Add LTMLE confounder- BMI, comorbidity diagnoses, comorbidity

d <- d %>% mutate(
  baseline_prob=0.5 * any_dementia + any_death + 0.1*income + 0.1*diabduration + 0.3*ifelse(edu=="Basic",1,0), 
  normalized_baseline_prob = (baseline_prob-min(baseline_prob))/(max(baseline_prob)-min(baseline_prob)))


summary(plogis(d$normalized_baseline_prob-1))

set.seed(12345)
d <- d %>%
  mutate(kidney_1= rexpit(normalized_baseline_prob-6),
         kidney_2=ifelse(kidney_1==1,1,rexpit(normalized_baseline_prob-6)),
         kidney_3=ifelse(kidney_2==1,1,rexpit(normalized_baseline_prob-6)),
         kidney_4=ifelse(kidney_3==1,1,rexpit(normalized_baseline_prob-6)),
         kidney_5=ifelse(kidney_4==1,1,rexpit(normalized_baseline_prob-6)),
         kidney_6=ifelse(kidney_5==1,1,rexpit(normalized_baseline_prob-6)),
         kidney_7=ifelse(kidney_6==1,1,rexpit(normalized_baseline_prob-6)),
         kidney_8=ifelse(kidney_7==1,1,rexpit(normalized_baseline_prob-6)),
         kidney_9=ifelse(kidney_8==1,1,rexpit(normalized_baseline_prob-6)),
         
         obese_1= rexpit(normalized_baseline_prob-1),
         obese_2= rexpit((normalized_baseline_prob+0.1*obese_1-1)),
         obese_3= rexpit((normalized_baseline_prob+0.1*obese_2-1)),
         obese_4= rexpit((normalized_baseline_prob+0.1*obese_3-1)),
         obese_5= rexpit((normalized_baseline_prob+0.1*obese_4-1)),
         obese_6= rexpit((normalized_baseline_prob+0.1*obese_5-1)),
         obese_7= rexpit((normalized_baseline_prob+0.1*obese_6-1)),
         obese_8= rexpit((normalized_baseline_prob+0.1*obese_7-1)),
         obese_9= rexpit((normalized_baseline_prob+0.1*obese_8-1)))

         


table(d$kidney_1)
table(d$kidney_9)
table(d$obese_1)
table(d$obese_9)

table(is.na(d))
table(is.na(d$T1))
colnames(d)
#d <- d %>% subset(., select=-c(weight,drop,baseline_prob,normalized_baseline_prob))

saveRDS(d, file=here("simulated data/novo_registry_simulated.RDS"))







