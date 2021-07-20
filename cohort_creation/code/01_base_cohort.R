#############################################
#############################################
# Project: Novo Nordisk GLP-1RA Project
# Purpose: base cohort creation file
# Input: pop, lmdb, lpr
# Output: cohort_file
# Author (date): Nerissa Nance (July 2021)
# Edits (author, date, description):

#############################################
#############################################
# SIMULATE source data
library(heaven)
library(data.table)
set.seed(05021992)
N <- 1000
pop <- simPop(N)
lpr <- simAdmissionData(N)
#import atc codes for diabetes meds
diabmeds <- rio::import("./diab_med_atcs.csv")
diabmeds$ATC <- stringr::str_trim(diabmeds$ATC)
atcs <- rep(list(c(200, 30)), length = nrow(diabmeds))
names(atcs) <- c(diabmeds$ATC)
lmdb <-
  data.table(simPrescriptionData(N, packages = atcs, startDate = "2012-01-01"))
head(lmdb)
table(lmdb$atc)

# Step 0(?): identify diabetes patients--require a dx and/or rx? Years?



#Step 1: find first instance of second line medication
# must be AFTER 2015 (grab date=> index)
secondline <-
  data.table(diabmeds[diabmeds$Type %in% c("glp1", "dpp4", "sglt2"),])
lmdb_secondline <-
  merge(
    lmdb,
    secondline,
    by.x = "atc",
    by.y = "ATC",
    all.x = F,
    all.y = F
  )

table(lmdb_secondline$Type)
table(lubridate::year(lmdb$eksd))
#sort first then take first element
pat_secondline <-
  lmdb[lubridate::year(eksd) >= 2015, .SD[which.min(eksd)], by = pnr]
setnames(pat_secondline, c("eksd", "atc"), c("index_dt", "index_atc"))
head(pat_secondline)

#Step 2: look back for prior metformin use
# just do 6 months for now
metformin <-
  data.table(lmdb[lmdb$atc %in% diabmeds$ATC[diabmeds$Name == "metformin"]])
table(metformin$atc)

metformin_usage <-
  merge(
    pat_secondline[, c("pnr", "index_dt", "index_atc")],
    metformin,
    by = 'pnr',
    all.x = T,
    all.y = F
  )

metformin_users <-
  metformin_usage[between(eksd, index_dt, index_dt + 180),
                  .SD[which.min(eksd)], by = pnr]
#need to add in some sanity checks
#question: what about metformin, then something else (sulfoneryas?), then glp

finaldata <- metformin_users[, .(pnr,index_dt,index_atc)]
head(finaldata)

#save base cohort dataset
objects <- c("finaldata")

save(list = objects, file = "./cohort_creation/data/output/base_cohort.Rdata")

