#############################################
#############################################
# Project: Novo Nordisk GLP-1RA Project
# Purpose: Cohort creation file
# Input: pop, lmdb, lpr
# Output: cohort_file
# Author (date): Nerissa Nance (July 2021)
# Edits (author, date):

#############################################
#############################################
# SIMULATE source data
library(heaven)
library(data.table)
library(rio)
set.seed(05021992)
N <- 1000
pop <- simPop(N)
lpr <- simAdmissionData(N*100, startDate = "2012-01-01")
#import atc codes for diabetes meds
diabmeds <- rio::import("./cohort_creation/data/reference/diab_med_atcs.csv")
atcs <- rep(list(c(200, 30)), length = nrow(dementiacodes[dementiacodes$codetype=="atc"]))
names(atcs) <- c(diabmeds$ATC)
lmdb <-
  data.table(simPrescriptionData(N, packages = atcs, startDate = "2012-01-01"))
head(lmdb)
table(lmdb$atc)

load("./cohort_creation/data/output/base_cohort.Rdata")


#EXCLUSION for insulin before index date
insulin_use <- 
  data.table(lmdb[lmdb$atc%in%diabmeds$ATC[diabmeds$Type=="insulins"]])
table(insulin_use$atc)

library(sqldf)

cohort_exclusions <- data.table(sqldf("select distinct a.*, 
max(case when b.pnr is not null then 1 else 0 end) as exclude_insulin
from metformin_users a
         left join insulin_use b on a.pnr=b.pnr
              and b.eksd<a.index_dt
          group by a.pnr"))
table(test$exclude_insulin)
#look at some examples
cohort_exclusions[pnr=="1"]
insulin_use[pnr=="1"]


#exclude for dementia diagnosis or dementia med prior to index:
(dementiacodes <- rio::import("./cohort_creation/data/reference/dementia_codes.csv"))
 dementia_meds <-
  data.table(lmdb[lmdb$atc%in%dementiacodes$code[dementiacodes$codetype=="atc"]])
table(dementia_meds$atc)
# dementia_inpatient <- lpr
head(lpr)
#subset to encounters with a dx for dementia
dementia_encounters <-  
  data.table(lpr[lpr$diag%in%dementiacodes$code[dementiacodes$codetype=="icd10"]])
table(dementia_encounters$diag)
#if the inpatient encounter started anytime before index, then count it
cohort_exclusions2 <- data.table(sqldf("select distinct a.*, 
max(case when b.pnr is not null then 1 else 0 end) as excl_dementia_rx
from cohort_exclusions a
         left join dementia_encounters b on a.pnr=b.pnr
              and b.inddto<a.index_dt
          group by a.pnr"))
table(cohort_exclusions2$excl_dementia_rx)

# QUESTIONS--
# - diabetes inclusion criteria
# - what does index date mean in the admission data?


#save as exclusions file 

