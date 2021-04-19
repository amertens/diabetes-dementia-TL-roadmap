# https://github.com/romainkp/LtAtStructuR
# library(remotes)
# install_github("romainkp/LtAtStructuR")
library(LtAtStructuR)
#runShiny()
library(data.table)
library(lubridate)
library(heaven)
#library(future) # optional (for parallel processing)
#plan(multiprocess) # optional (for parallel processing)

#Step 1: simulate the raw data formats
set.seed(05021992)
N<- 1000
cohort_raw <- simPop(N)
lmdb=simPrescriptionData(N,packages=packs)
lpr=simAdmissionData(N)
rio::export(cohort_raw, "./data_raw/cohort_data.csv")


#Step 2: preprocessing of med and admission data
packs = list("R03AK11"=list(c(750,75),c(500,200),c(400,200)),
             "R03AL03"=list(c(750,75),c(500,200),c(400,200)),
             "C01CA01"=list(c(200,100),c(750,30)))
lpr <- getAdmLimits(lpr,collapse=TRUE)
drug1 = list(atc=c("R03AK11","R03AL03","R03AC02","R03AC04","R03AC19",
                   "R03AL02","R03AA01","R03AC18","R03AL01"),
             maxdepot=4000,
             period=as.Date(c("1995-01-01", "2012-12-31")),
             prescriptionwindow=2,
             doses=list(value=c(750,500,400,200,75),
                        min = c(250,200,200,100,25),
                        max = c(1000,600,800,600,100),
                        def = c(750,500,400,200,75)))
drug2=list(atc=c("C01CA01","C01AA05"),
           maxdepot=4000,
           period=as.Date(c("1995-01-01", "2012-12-31")),
           prescriptionwindow=2,
           doses=list(value=c(200, 400, 500,750),
                      min = c(100, 200, 250,750),
                      max = c(300, 800, 1000,750),
                      def = c(200, 400, 500,750)))
med_data <- medicinMacro(drugs=list("drug1"=drug1,"drug2"=drug2),drugdb=lmdb,
                         admdb=lpr)
med_data_1 <-cbind(med_data$drug1, type="drug1")
head(med_data_1)
med_data_2 <-cbind(med_data$drug2, type="drug2")
head(med_data_2)

med_data_all <- rbind(med_data_1,med_data_2)
head(med_data_all)


med_data_1[,lag_last := shift(lastday,1,"lag"),by=pnr]
med_data_1[,lead_last := shift(lastday,1,"lead"),by=pnr]
med_data_1[,lead_first := shift(firstday,1,"lead"),by=pnr]

med_data_1$firstday_cp <- mapply( function(x,y){
  if(!is.na(y) & x==y){
  as.Date(x+1)
  }else{
    as.Date(x)
}}, med_data_1$firstday, med_data_1$lag_last)
med_data_1$firstday_cp <- as.Date(med_data_1$firstday_cp, origin='1970-01-01')
head(med_data_1)
rio::export(med_data_1, "exposure_data.csv")

cohortDT <- fread(file='cohort_data.csv') 
cohortDT[,pnr:=as.character(pnr)] 
cohortDT[,birthdate:=as_date(parse_date_time(birthdate, c('Ymd','Ydm','mdY','mYd','dmY','dYm')))] 
cohortDT[,doddate:=as_date(parse_date_time(doddate, c('Ymd','Ydm','mdY','mYd','dmY','dYm')))] 
cohortDT[,status:=as.character(status)] 
col.to.num <- names(cohortDT)[sapply(cohortDT,class)%in%'integer'] 
cohortDT[,eval(col.to.num) := as.numeric(get(col.to.num))] 
cohort <- setCohort(data = cohortDT, IDvar = 'pnr', index_date = 'birthdate', EOF_date = 'doddate', EOF_type = 'status', Y_name = '1', L0 = c('sex'), L0_timeIndep = list('sex'=list('categorical'=TRUE,'impute'='mode','impute_default_level'=NA))) 

exposureDT <- fread(file='exposure_data.csv') 
exposureDT[,pnr:=as.character(pnr)]
exposureDT[,firstday_cp:=as_date(parse_date_time(firstday_cp, c('Ymd','Ydm','mdY','mYd','dmY','dYm')))]
exposureDT[,lastday:=as_date(parse_date_time(lastday, c('Ymd','Ydm','mdY','mYd','dmY','dYm')))]
exposure <- setExposure(data = exposureDT, IDvar = 'pnr', start_date = 'firstday_cp', end_date = 'lastday', exp_level = NA, exp_ref = NA) 

LtAt.specification <- cohort + exposure 
system.time(LtAt.data <- construct(LtAt.specification, time_unit = 100, 
                                first_exp_rule = 1, exp_threshold = 0.5))
head(LtAt.data)


fwrite(LtAt.data, file = './data/clean/final_dataset.csv', row.names = FALSE)
save(LtAt.data, file="./data/clean/finaldata.RData")
