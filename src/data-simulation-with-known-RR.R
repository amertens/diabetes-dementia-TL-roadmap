
library(here)
library(tidyverse)
library(simstudy)
library(simcausal)

#Replicate sampled dataset:
data <- readRDS(here("simulated data/novo_registry_simulated.RDS"))
head(data)

D <- DAG.empty() + 
  node("CVD", distr="rcat.b1", probs = c(0.5, 0.25, 0.25)) +
  node("A1C", t=0, distr="rnorm", mean=5 + (CVD > 1)*10 + (CVD > 2)*5) + 
  node("T", t=0, distr="rbern", prob=plogis(-5 - 0.3*CVD + 0.5*A1C[t])) +
  
  node("A1C", t=1:9, distr="rnorm", mean=-T[t-1]*10 + 5 + (CVD > 1)*10 + (CVD > 2)*5) +
  node("T", t=1:9, distr="rbern", prob=plogis(-5 - 0.3*CVD + 0.5*A1C[t] + 1.5*T[t-1])) +
  node("Y", t=0:9, distr="rbern", prob=plogis(-6 - 1.2*T[t] + 0.1*CVD + 0.3*A1C[t]), EFU=TRUE)
D <- set.DAG(D)
dat <- sim(D,n=200, LTCF="Y", verbose=T, rndseed=12345) #Do I need to fill in censoring with LTCF argument?

head(dat)
