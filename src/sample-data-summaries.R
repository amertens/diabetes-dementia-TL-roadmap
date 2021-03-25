
library(here)
library(tidyverse)

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
ndf <- 10000
#Sample with replacement with the probability of the observed frequency
glp1_sim <- glp1[sample(nrow(glp1), ndf, replace = TRUE, prob = glp1$weight), ]
sglt2_sim <- sglt2[sample(nrow(sglt2), ndf, replace = TRUE, prob = sglt2$weight), ]
glp1_sim$A <- "glp1"
sglt2_sim$A <- "sglt2"

d <- bind_rows(glp1_sim, sglt2_sim)
saveRDS(d, file=here("simulated data/novo_registry_simulated.RDS"))





