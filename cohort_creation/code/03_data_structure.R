# data structure based on 1) original data tables; 2) selected pnr/IDs after exclusions
# Note that we do note apply exclusion rules in this file
# At is the (t-1)-th time interval drug usage; to avoid reverse causality; last time interval prescriptions are dropped
# each prescription assumes 6mos/180d use
# TODO
# current simPrescriptionData() are too frequent; severe 2nd line drug overlap, everyone developed dementia at some point
# if single prescription coverage is not 180d, it should be more careful when deciding At
# add insulin
# add baseline/history hypertension; maybe hypertension as categorical
# Questions
# count 9mos/270d interruption after the 180d coverage or not? 3mos no coverage vs 9mos (current) no actual use; but 3mos is essentially always 1
# does exclusion rule also exclude patients without continuous 2nd line drug use? 
# years between first metformin use and first 2nd line drug use & accumulated metformin prescriptions, as baseline covariates? 
# Any other covariates to be generated from diseasecode? 

library(dplyr)
library(data.table)
library(heaven)

# atc/icd codes from group document, used in the last version; checked they are identical with the info below
glp1 <- c("A10BJ01", "A10BJ02", "A10BJ03", "A10BJ04", "A10BJ05", "A10BJ06", "A10BJ07")  
dpp4 <- c("A10BH01", "A10BH02", "A10BH03", "A10BH04", "A10BH05", "A10BH06", "A10BH07", "A10BH08", "A10BH51", "A10BH52")
SGLT2i <- c("A10BK01", "A10BK02", "A10BK03", "A10BK04", "A10BK05", "A10BK06", "A10BK07")
dementia <- c("N06DA01", "N06DA02", "N06DA03", "N06DA04", "N06DA05", "N06DA52", "N06DA53")
# diabete drugs; glp1, dpp4, sglt2 are used for intervention nodes; 
diabmeds <- rio::import("./cohort_creation/data/reference/diab_med_atcs.csv")
dementiacodes <- rio::import("./cohort_creation/data/reference/dementia_codes.csv")
hypertension_list <- hypertensionATC

# generate "raw" data tables
N <- 10000
# atc codes include diabete drugs, dementia, and hypertension
names_atcs <- c(diabmeds$ATC, 
                dementiacodes[dementiacodes$codetype == "atc", "code"], 
                unlist(hypertension_list)
)
atcs <- rep(list(c(200, 30)), length = length(names_atcs))
names(atcs) <- names_atcs
set.seed(123)
lmdb <- data.table(simPrescriptionData(N, max.prescriptions = 5, 
                                       packages = atcs, startDate = "2012-01-01"))
set.seed(123)
pop <- simPop(N)
set.seed(123)
lpr <- simAdmissionData(N, startDate = "2012-01-01", 
                        diagnoses = c("DN162D", "DV1180", "DN982", "DT698", 
                                      "DJ343", "DP389C", "DD484", "DB741", "DO721A", "DQ728D", 
                                      "DK254E", "DT635", "DB601A", "DD239E", "DQ794A", 
                                      "DO010", "DL923B", "DD223Z", "DF0122", "DZ237", "DE519", 
                                      "DG461", "DO472", "DK265D", "DN330", "DM92", "DUA19", 
                                      dementiacodes[dementiacodes$codetype == "icd10", "code"]
                                      )  # make sure there are dementia records
                        )
# make sure the dates are ordered before atc types
setkey(lmdb, pnr, eksd, atc)

# assume that we have a list of selected subject IDs
set.seed(123)
selected_IDs <- sample(N, round(N*4/5)) %>% sort


# generate structured data for the input subject list
dt_tmle <- data.table(pnr = selected_IDs)
# index date
dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Type %in% c("glp1", "dpp4", "sglt2"), "ATC"], .(date = eksd %>% first), pnr], 
        on = "pnr", first_date_2nd_line := i.date]  # search for the dates of first 2nd line drug use (index date)
# age
dt_tmle[pop, on = "pnr", age := as.numeric(first_date_2nd_line - i.birthdate)/365.25]
# sex
dt_tmle[pop, sex := (i.sex == "female")*1, on = "pnr"]
# any previous metformin use
dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Name %in% c("metformin"), "ATC"], .(date = eksd %>% first), pnr], 
        on = "pnr", first_date_metformin := i.date]  # search for the dates of first 2nd line drug use (index date)

# decide intervention groups
{
  # find people with <9mos gap of glp1 use (can be transformed with 6 mos cover of each prescription) among target population
  max.date <- as.numeric(lmdb$eksd) %>% max  # last available date; maybe the end date
  dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Type %in% c("glp1"), "ATC"], .(test = all(c(eksd, max.date)[-1] - (eksd+180) < 270)), by = "pnr"], 
          "if_target_glp1" := i.test,
          on = "pnr"]
  # allow starting another 2nd line drug first, but the first "interruption" between the index date to the glp1 start date should also <270
  dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Type %in% c("glp1"), "ATC"], .(date = eksd %>% first), pnr], 
          on = "pnr", first_date_glp1 := i.date]  # search for the dates of first 2nd line drug use (index date)
  dt_tmle[, if_target_glp1 := if_target_glp1 & (first_date_glp1 - first_date_2nd_line) < 270]
  
  # find people with <9mos gap of dpp4 use (can be transformed with 6 mos cover of each prescription) among target population
  dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Type %in% c("dpp4"), "ATC"], .(test = all(c(eksd, max.date)[-1] - (eksd+180) < 270)), by = "pnr"], 
          "if_target_dpp4" := i.test,
          on = "pnr"]
  # allow starting another 2nd line drug first, but the first "interruption" between the index date to the glp1 start date should also <270
  dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Type %in% c("glp1"), "ATC"], .(date = eksd %>% first), pnr], 
          on = "pnr", first_date_dpp4 := i.date]  # search for the dates of first 2nd line drug use (index date)
  dt_tmle[, if_target_dpp4 := if_target_dpp4 & (first_date_dpp4 - first_date_2nd_line) < 270]
  
  # find people with <9mos gap of sglt2 use (can be transformed with 6 mos cover of each prescription) among target population
  dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Type %in% c("sglt2"), "ATC"], .(test = all(c(eksd, max.date)[-1] - (eksd+180) < 270)), by = "pnr"], 
          "if_target_sglt2" := i.test,
          on = "pnr"]
  # allow starting another 2nd line drug first, but the first "interruption" between the index date to the glp1 start date should also <270
  dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Type %in% c("glp1"), "ATC"], .(date = eksd %>% first), pnr], 
          on = "pnr", first_date_sglt2 := i.date]  # search for the dates of first 2nd line drug use (index date)
  dt_tmle[, if_target_sglt2 := if_target_sglt2 & (first_date_sglt2 - first_date_2nd_line) < 270]
  
  table(dt_tmle$if_target_glp1 + dt_tmle$if_target_dpp4 + dt_tmle$if_target_sglt2)  # severe 2nd line drugs overlap; should be avoided in simulation
}





# define wide format A1t nodes for glp1 use; 1/2 year time interval; index date is t=0; (index date, index date + interval] is t=1
date_cuts <- as.numeric(paste0(2015:2020, "-01-01") %>% sapply(as.Date))
date_to_t <- function(date, start) {
  ((date - start) %/% (365.25/2)) - ((date - start) %% (365.25/2) == 0) * 1 + 1
}
date_to_t(c(NA, date_cuts), date_cuts[1])
date_to_t(c(0, 
            1, 365.25/2, 
            365.25/2 + 0.01, 365, 365.25, 
            367), 0)

# max should be the maximum time point
fill_zero <- function(vec, max = 10) {
  temp <- rep(0, max)
  vec <- vec[vec <= max]  # caps at max; might be used when ignoring prescriptions in the last interval
  temp[vec] <- 1
  return(temp)
}
event_free <- function(i, max = 10) {
  temp <- rep(1, max)
  if (!is.na(i)) temp[i:max] <- 0
  return(temp)
}



# find all time points with glp1 prescription
# consider 180d<365.25 coverage for each prescription; one prescription covers at most two time interval
for (x in paste0("A1_", 1:10)) dt_tmle[, (x) := 0]
dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% diabmeds[diabmeds$Type %in% c("glp1"), "ATC"], 
             c(date_to_t(as.numeric(eksd), min(as.numeric(eksd))) + 1, date_to_t(as.numeric(eksd) + 180, min(as.numeric(eksd))) + 1), 
             by = "pnr"][, V1 %>% as.vector %>% fill_zero %>% as.list, by = "pnr"], 
        ':=' (paste0("A1_", 1:10), mget(paste0("i.V", 1:10))), 
        on = "pnr"]
dt_tmle[which(if_target_glp1), ]
dt_tmle[if_target_glp1 == F, ]
dt_tmle %>% head(10)




# for target population, collect hypertension prescription dates
# make ht nodes; baseline/history hypertension is generated elsewhere
for (x in paste0("L_", 1:10)) dt_tmle[, (x) := 0]
dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% unlist(hypertension_list), 
             date_to_t(as.numeric(eksd), min(as.numeric(eksd))), 
             by = "pnr"][V1 > 0, ][, V1 %>% as.vector %>% fill_zero %>% as.list, by = "pnr"] , 
        ':=' (paste0("L_", 1:10), mget(paste0("i.V", 1:10))), 
        on = "pnr"]


# make event process
dt_tmle[lmdb[pnr %in% selected_IDs & atc %in% dementiacodes[dementiacodes$codetype == "atc", "code"], .(date = eksd %>% first), pnr], 
        on = "pnr", first_date_dementia := i.date]  # search for the dates of first 2nd line drug use (index date)
dt_tmle[lpr[pnr %in% selected_IDs & diag %in% dementiacodes[dementiacodes$codetype == "icd10", "code"], .(date = inddto %>% first), pnr], 
        on = "pnr", first_date_dementia_ad := i.date]
dt_tmle[, first_date_dementia:=ifelse(!is.na(first_date_dementia_ad) & first_date_dementia_ad < first_date_dementia, first_date_dementia_ad, first_date_dementia) ]
dt_tmle[, first_date_dementia:=as.Date(first_date_dementia, "1970-01-01")]

dt_tmle[dt_tmle[, date_to_t(first_date_dementia %>% as.numeric, first_date_2nd_line %>% as.numeric), by = "pnr"][, V1 %>% event_free %>% as.list, by = "pnr"], 
        ':=' (paste0("Y_", 1:10), mget(paste0("i.V", 1:10))), 
        on = "pnr"]

# censoring process; study ends
censoring_free <- function(start, last_date = max.date, event_t, max = 10) {
  if (is.na(event_t)) {
    return(((start + 365.25/2 * (1:max)) <= last_date) * 1)
  } else {
    temp <- rep(0, max)
    temp[1:event_t] <- 1
    return(temp)
  }
}

dt_tmle[, event_t := which(mget(paste0("Y_", 1:10)) == 0) %>% first, by = "pnr"]
dt_tmle[, paste0("C_", 1:10) := (first_date_2nd_line %>% as.numeric) %>% censoring_free(event_t = event_t) %>% as.list, by = "pnr"]
dt_tmle[, event_t := NULL]
