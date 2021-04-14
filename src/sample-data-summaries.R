
library(here)
library(tidyverse)
library(holodeck)
rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))

#-------------------------------- 
# Load data summaries
#-------------------------------- 
glp1 <- read.csv(here::here("data summaries/glp1-regimes-npar-table.csv"), header = FALSE)
sglt2 <- read.csv(here::here("data summaries/sglt2-regimes-npar-table.csv"), header = FALSE)
head(glp1)
head(sglt2)

#-------------------------------- 
# Give variable names to datasets
#-------------------------------- 

# in the glp1 and sglt2 tables I have added age, sex stroke, ami, dementia
# and death information but did not bother to make the wide format. so,
# columns are:
#   
#   c1: agegroups 1=40-50, 2=50-60, 3=60-70, 4=70-80, 5=80+
#   c2: sex: 0 = male, 1= female
# c3-c11: drug regime in 9 periods (each period is 6 months long)
# c12: period where death occurs
# c13: period where dementia occurs
# c14: period where stroke occurs
# c15: period where ami occurs
# c16: relative frequency of the whole pattern
colnames(glp1) <- c("age","sex","T1","T2","T3","T4","T5","T6","T7","T8","T9","death","dementia","stroke","ami","weight")
colnames(sglt2) <- c("age","sex","T1","T2","T3","T4","T5","T6","T7","T8","T9","death","dementia","stroke","ami","weight")


#-------------------------------- 
# Make dataset by sampling relative 
# frequency of each covariate/treatment pattern
#-------------------------------- 

#Dataset size
ndf <- 550000
#Sample with replacement with the probability of the observed frequency
set.seed(12345)
glp1_sim <- glp1[sample(nrow(glp1), ndf, replace = TRUE, prob = glp1$weight), ]
sglt2_sim <- sglt2[sample(nrow(sglt2), ndf, replace = TRUE, prob = sglt2$weight), ]
glp1_sim$A <- "glp1"
sglt2_sim$A <- "sglt2"

d <- bind_rows(glp1_sim, sglt2_sim)

#Drop gaps in treatment
dim(d)
d <- d %>% mutate(drop=
        case_when(T1==0 & T2==0 ~ 1,
                  T2==0 & T3==0 ~ 1,
                  T3==0 & T4==0 ~ 1,
                  T4==0 & T5==0 ~ 1,
                  T5==0 & T6==0 ~ 1,
                  T6==0 & T7==0 ~ 1,
                  T7==0 & T8==0 ~ 1,
                  T8==0 & T9==0 ~ 1)) %>%
  filter(is.na(drop))

dim(d)
head(d)
table(d$death)
table(d$dementia)

d <- d %>% mutate(
  age=factor(age),
  any_death = ifelse(is.na(as.numeric(death)),0,1),
  any_dementia= ifelse(is.na(as.numeric(dementia)),0,1)
)

table(d$any_death)
table(d$any_dementia)

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

colnames(d)
d <- d %>% subset(., select=-c(weight,drop,baseline_prob,normalized_baseline_prob))

saveRDS(d, file=here("simulated data/novo_registry_simulated.RDS"))







