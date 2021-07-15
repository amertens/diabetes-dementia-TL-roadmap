# completely fake distribution; to incorporate summary statistics info (train lkd objs with fake+intermediate samples) later
# fully ATC prescription based; ICD-8/10 or other info later

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

library(dplyr)
library(data.table)
library(heaven)

source("./simulation_ZW/simPrescriptionData_atc.R")
glp1 <- c("A10BJ01", "A10BJ02", "A10BJ03", "A10BJ04", "A10BJ05", "A10BJ06", "A10BJ07")  
dpp4 <- c("A10BH01", "A10BH02", "A10BH03", "A10BH04", "A10BH05", "A10BH06", "A10BH07", "A10BH08", "A10BH51", "A10BH52")
SGLT2i <- c("A10BK01", "A10BK02", "A10BK03", "A10BK04", "A10BK05", "A10BK06", "A10BK07")
hypertension <- unlist(hypertensionATC)
dementia <- c("N06DA01", "N06DA02", "N06DA03", "N06DA04", "N06DA05", "N06DA52", "N06DA53")

set.seed(123)
n_total_population <- 100000
dt_atc <- simPrescriptionData_atc(n = n_total_population, startDate = "2015-01-01")
dt_subject <- simPop(n = n_total_population)

# target population
# start any second line drugs after 2015; record the date
dt_target_population <- dt_atc[, .(if_any_2nd_line = any(atc %in% c(glp1, dpp4, SGLT2i))), pnr]
dt_target_population[dt_atc[, .(date = eksd[(which(atc %in% c(glp1, dpp4, SGLT2i)))] %>% sort %>% first), pnr], on = "pnr", first_date_2nd_line := i.date]
dt_target_population[dt_atc[, .(date = eksd[(which(atc %in% c(glp1)))] %>% sort %>% first), pnr], on = "pnr", first_date_glp1 := i.date]

# exclude prior dementia
dt_target_population[dt_atc[, .(date = eksd[(which(atc %in% c(dementia)))] %>% sort %>% first), pnr], on = "pnr", first_date_dementia := i.date]
dt_target_population[, "if_prior_dementia" := (first_date_dementia < max(first_date_2nd_line, first_date_glp1, na.rm = T)), by = "pnr"]
dt_target_population[, "if_still_in" := if_any_2nd_line & (is.na(if_prior_dementia) | !if_prior_dementia)]
# exclude age under 50
dt_target_population[dt_subject, on = "pnr", age := as.numeric(first_date_2nd_line - i.birthdate)/365.25]
dt_target_population[, if_still_in := if_still_in & (age >= 50)]

# find people with <12mos gap of glp1 use among target population
max.date <- as.numeric(dt_atc$eksd) %>% max
dt_atc[dt_target_population, target_glp1_dates := ifelse(i.if_still_in & (atc %in% glp1), eksd, NA), on = "pnr"]  # for target population, record glp1 dates
dt_atc[dt_target_population, target_2nd_line_dates := ifelse(i.if_still_in & (atc %in% c(glp1, dpp4, SGLT2i)), eksd, NA), on = "pnr"]  # for target population, record glp1 dates
dt_target_population[dt_atc[, .(sort(target_glp1_dates)), by = "pnr"][!is.na(V1), ][, .(test = all(diff(c(V1, max.date)) < 270)), by = "pnr"], 
                     "if_target_glp1" := i.test,
                     on = "pnr"]
# if there is no glp1 use at all, but in target population, label as not glp1 group
dt_target_population[, if_target_glp1 := ifelse(is.na(if_target_glp1), ifelse(if_still_in, F, NA), if_target_glp1)]
dt_target_population$if_target_glp1 %>% table(useNA = "a")
dt_target_population$if_still_in %>% table(useNA = "a")

# define id and baseline nodes
dt_tmle <- dt_target_population[which(if_still_in), c("pnr","age")]
dt_tmle[dt_subject, sex := (i.sex == "female")*1, on = "pnr"]

# define wide format A1t nodes for glp1 use
dt_tmle[dt_target_population, if_glp1_group := i.if_target_glp1, on = "pnr"]
dt_tmle[dt_target_population, first_date_glp1 := i.first_date_glp1, on = "pnr"]
dt_tmle[dt_target_population, first_date_2nd_line := i.first_date_2nd_line, on = "pnr"]

date_cuts <- as.numeric(paste0(2015:2020, "-01-01") %>% sapply(as.Date))
date_to_t <- function(date, start) {
  (date - start) %/% 365.25 + 1
}
date_to_t(c(NA, date_cuts), date_cuts[1])

# max should be the maximum time point
fill_zero <- function(vec, max = 5) {
  temp <- rep(0, max)
  temp[vec] <- 1
  return(temp)
}
event_free <- function(i, max = 5) {
  temp <- rep(1, max)
  if (!is.na(i)) temp[i:max] <- 0
  return(temp)
}

# find all time points with glp1 prescription
for (x in paste0("A1_", 1:5)) dt_tmle[, (x) := 0]
dt_tmle[dt_atc[pnr %in% dt_tmle$pnr, ifelse(is.na(target_glp1_dates), 0,
                                            date_to_t(target_glp1_dates, min(as.numeric(target_glp1_dates), na.rm = T))
), by = "pnr"][V1 != 0, ][, V1 %>% as.vector %>% fill_zero %>% as.list, by = "pnr"] 
, ':=' (paste0("A1_", 1:5), mget(paste0("i.V", 1:5))), on = "pnr"]
dt_tmle[which(if_glp1_group), ]
dt_tmle %>% head(10)

# for target population, collect hypertension prescription dates
dt_atc[dt_target_population, target_ht_dates := ifelse(i.if_still_in & (atc %in% hypertension), eksd, NA), on = "pnr"]  # for target population, record glp1 dates
# make ht nodes
for (x in paste0("L_", 1:5)) dt_tmle[, (x) := 0]
dt_tmle[dt_atc[pnr %in% dt_tmle$pnr, ifelse(is.na(target_ht_dates), 0,
                                            date_to_t(target_ht_dates, ifelse(any(!is.na(target_glp1_dates)), min(as.numeric(target_glp1_dates), na.rm = T), min(as.numeric(target_2nd_line_dates), na.rm = T)) )
), by = "pnr"][V1 > 0, ][, V1 %>% as.vector %>% fill_zero %>% as.list, by = "pnr"] 
, ':=' (paste0("L_", 1:5), mget(paste0("i.V", 1:5))), on = "pnr"]

# make event process
dt_tmle[dt_target_population[pnr %in% dt_tmle$pnr, date_to_t(first_date_dementia %>% as.numeric, max(first_date_2nd_line, first_date_glp1, na.rm = T) %>% as.numeric  ), by = "pnr"]
        [, V1 %>% event_free %>% as.list, by = "pnr"] 
,   ':=' (paste0("Y_", 1:5), mget(paste0("i.V", 1:5))), on = "pnr"]

# censoring process; study ends
last_date <- date_cuts %>% last
censoring_free <- function(start, last_date = date_cuts %>% last, event_t, max = 5) {
  if (is.na(event_t)) {
    return(((start + 365.25 * (1:5)) <= last_date) * 1)
  } else {
    temp <- rep(0, max)
    temp[1:event_t] <- 1
  }
}

dt_tmle[, event_t := which(mget(paste0("Y_", 1:5)) == 0) %>% first, by = "pnr"]
dt_tmle[, paste0("C_", 1:5) := (max(first_date_2nd_line, first_date_glp1, na.rm = T) %>% as.numeric) %>% censoring_free(event_t = event_t) %>% as.list, by = "pnr"]
dt_tmle[, event_t := NULL]


saveRDS(dt_tmle, file = "./simulation_ZW/dt_tmle_20210712.rds")

