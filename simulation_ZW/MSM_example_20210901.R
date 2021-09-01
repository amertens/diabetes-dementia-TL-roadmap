library(dplyr)
library(data.table)
# library(ltmle)
# load modified version; mcSuperLearner, deterministic censoring, etc. 
{
  library(Matrix)
  library(matrixStats)
  library(speedglm)
  temp <- list.files("./simulation_ZW/ltmle-master/R", full.names = T)
  for (i in temp) source(i)
}
library(hal9001)
library(pryr)

###
# demonstration of simulated outcome distribution
###

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))
coef_intercept <- 3.5
coef_a <- -0.3
coef_b <- 0.1
expit(coef_intercept + coef_a*(10) + coef_b*(0:10)) %>% plot  # endpoint survival prob by accumulated exposure
expit(coef_intercept + coef_a*(1:10) + coef_b*(0)) %>% plot(type = "l")  # survival curve with 0 exposure
expit(coef_intercept + coef_a*(1:10) + coef_b*(1:10)) %>% lines(col = "red")  # survival curve with full exposure

###
# load fully random data; make positivity issue less severe; incorporate the outcome distribution above
###

dt_tmle <- readRDS("./simulation_ZW/dt_tmle_202108.rds")
set.seed(200)
# add more "treated" subjects to make positivity issue less severe for now
dt_tmle <- dt_tmle[c(sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) > 5), 100000, T), sample(nrow(dt_tmle), 100000, T)), ]

K <- 10 
node_names <- c("age", "sex", "L_0", 
                "first_date_2nd_line",
                expand.grid(c("L", "Y", "A1", "C"), as.character(1:(K+1))) %>% apply(1, function(row) paste0(row, collapse = "_"))
)
node_names <- node_names[!node_names %in% c(paste0(c("C_", "A1_"), K+1))]
# create random L_0
dt_tmle[, "L_0" := sample(0:1, nrow(dt_tmle), replace = T)]
dt_use <- dt_tmle[, ..node_names]
# max.date <- max(dt_use$first_date_2nd_line)
max.date <- as.Date("2016-12-31")  # the last day of possible follow-up

# ltmle use 1 as event
for (x in grep("Y_", node_names)) dt_use[, (node_names[x]) := 1 - get(node_names[x])]

# noninformative censoring with 0.05 prob; but will be coded back into the fake index date, so there is no need to fit C models
set.seed(123)
# censoring is the within interval change
prob_censoring <- 0.05
dt_use[, C_1 := rbinom(nrow(dt_use), 1, 1-prob_censoring)]
for (i in 2:10) {
  risk_set <- dt_use[[paste0("C_", i-1)]] == 1
  temp_input <- sapply(rep(1-prob_censoring, sum(risk_set)), function(each_p) rbinom(1, 1, each_p))
  dt_use[risk_set, ':=' (paste0("C_", i), temp_input)]
  dt_use[!risk_set, ':=' (paste0("C_", i), 0)]
}

# regenerate outcomes according to target distribution
# Yt is the start of interval status; always start with event-free
# dt_use[, Y_1 := 1-rbinom(nrow(dt_use), 1, expit(3))]
dt_use[, Y_1 := 0]
for (i in 2:11) {
  risk_set <- !is.na(dt_use[[paste0("Y_", i-1)]]) & dt_use[[paste0("Y_", i-1)]] == 0 & dt_use[[paste0("C_", i-1)]] == 1 
  if (i==2) prev_accumulated_time <- rep(0, nrow(dt_use)) else prev_accumulated_time <- dt_use[, paste0("A1_", 1:(i-2)), with=F] %>% rowSums
  accumulated_time <- dt_use[, paste0("A1_", 1:(i-1)), with=F] %>% rowSums
  if (i > 2) {
    prob_ratio <- expit(coef_intercept + coef_a*(i-1) + coef_b*accumulated_time)/expit(coef_intercept + coef_a*(i-2) + coef_b*prev_accumulated_time)
    # /(1-prob_censoring)
  } else {
    prob_ratio <- expit(coef_intercept + coef_a*(i-1) + coef_b*accumulated_time)
    # /(1-prob_censoring)
  }
  prob_ratio <- ifelse(prob_ratio <=1, prob_ratio, 1)
  prob_ratio <- ifelse(prob_ratio >=0, prob_ratio, 0)
  temp_input <- sapply(
    # rep(0.95, sum(risk_set))
    # (1 - 0.03*accumulated_time)[which(risk_set)]
    prob_ratio[which(risk_set)]
    , function(each_p) rbinom(1, 1, each_p))
  dt_use[risk_set, ':=' (paste0("Y_", i), 1-temp_input)]
  dt_use[!risk_set, ':=' (paste0("Y_", i), ifelse(dt_use[!risk_set, paste0("Y_", i-1), with=F] == 1, 1, NA))]
}


# code censoring into a regenerated index date
# first time bin t where censoring status changes
censored_t <- dt_use[, paste0("C_", 1:10)] %>% apply(1, function(eachRow) first(which(eachRow == 0))) %>% lapply(function(x) ifelse(length(x) == 0, 0, x)) %>% unlist
censored_t[censored_t==0] <- 11
fake_index_dates <- max.date - (censored_t - 1) * (365.25/2) - sample(10:170, nrow(dt_use), T)

# check fake index date inferred Ct process matches observed censoring
temp <- floor((max.date - fake_index_dates) / (365.25/2)) + 1 # decide max.date is in which interval, then it means censoring status changes within this interval
table(temp)
censored_t %>% table

dt_use[, first_date_2nd_line := fake_index_dates]

# keep it as NP truth for different sample sizes
dt_use_backup <- dt_use %>% copy


###
# generate 10k sample, and fit ltmleMSM
###

set.seed(123)
dt_use <- dt_use_backup[sample(nrow(dt_use_backup), 10000, T), ]

# deterministic censoring
administrativeCensoring <- function(data, current.node, nodes, last_day = max.date, which_baseline_date = 4) {
  Cnodes <- nodes$C
  if (!(current.node %in% Cnodes)) 
    return(NULL)
  current_t <- which(Cnodes == current.node)  # this is the censored t
  temp_t <- floor((last_day - data[[which_baseline_date]]) / (365.25/2)) + 1
  is.deterministic <- rep(T, nrow(data))
  return(list(is.deterministic = is.deterministic, prob1 = ifelse(temp_t <= current_t, 0, 1)))
}

# code summary covaraites in MSM models; now include accumulated exposure and time
n <- nrow(dt_use)
time.points <- 10
n_pool <- 7
regime.matrix <- as.matrix(expand.grid(rep(list(0:1), time.points)))
dim(regime.matrix)
num.regimes <- 2^time.points
regimes <- array(dim = c(n, time.points, num.regimes)) #n x numAnodes x numRegimes = n x time.points x 2^time.points
# summary.measures <- array(dim = c(num.regimes, 1, 1)) #numRegimes x num.summary.measures x num.final.Ynodes = 2^time.points x 1 x 1
summary.measures <- array(dim = c(num.regimes, 2, n_pool)) #numRegimes x num.summary.measures x num.final.Ynodes
test.treated <- array(dim = c(num.regimes, 1, 1))
for (i in 1:num.regimes) {
  regimes[, , i] <- matrix(regime.matrix[i, ], byrow = TRUE, nrow = n, ncol = time.points)
  for (j in 1:n_pool) {
    summary.measures[i, , j] <- c(sum(regime.matrix[i, j]), j)  # accumulated exposure, time
  }
  # test.treated[i, 1, 1] <- !any(diff(which(regime.matrix[i, ] == 0)) == 1)
  test.treated[i, 1, 1] <- all(diff(c(0, which(regime.matrix[i, ] == 1), 11)) <= 3)  # max interruption: this number -1
}

# colnames(summary.measures) <- "time.on.treatment"
colnames(test.treated) <- "if.in.group"
test.treated[, 1, 1] %>% table

summary.measures <- summary.measures[test.treated[, 1, 1], , ]
# %>% array(dim = c(sum(test.treated[, 1, 1]), 1, 1))
# colnames(summary.measures) <- "time.on.treatment"
for (j in 1:time.points) {
  colnames(summary.measures[, , j]) <- c("time.on.treatment", "time")
}
colnames(summary.measures) <- c("time.on.treatment", "time")

{
  start_time <- Sys.time()
  ss <- mem_change(
    test <- ltmleMSM(dt_use, Anodes = grep("^A1_", node_names), Lnodes = paste0("L_", 1:11), Ynodes = grep("^Y_", node_names), Cnodes = grep("^C_", node_names), survivalOutcome = T, 
                     # SL.library = c("SL.glm"),
                     # SL.library = c("SL.glm", c("SL.glm","screen.corP")),
                     regimes = regimes[, , test.treated[, 1, 1]], 
                     summary.measures = summary.measures, 
                     working.msm = "Y~time.on.treatment + time",
                     # working.msm = "Y~1",
                     variance.method = "ic",  # direct EIC plug in; might be underestimated with positivity and rare outcomes
                     # final.Ynodes = "Y_11",
                     final.Ynodes = paste0("Y_", seq(to = K+1, length.out = n_pool, by = 1)),   # to pool across these nodes
                     msm.weights = "empirical",  # h weights by obs data
                     # SL.cvControl = list(V = 8)  # control CV fold numbers; might be used in paralleled version
                     deterministic.g.function = administrativeCensoring,  # incorporate deterministic censoring based on index dates
                     estimate.time = F  # do not run on a n=50 subsample for predicting computation time
                     
    )
  )
  end_time <- Sys.time()
}
ss  # MEM change in MB
test  # coef est
difftime(end_time, start_time, units = "mins")  # time in min

# function of checking np truth
np_truth <- function(t, ss) {
  vec_rowsums <- regime.matrix[as.vector(test.treated), 1:t] %>% rowSums
  target_dt <- regime.matrix[as.vector(test.treated), 1:t] %>% as.data.table
  target_dt <- target_dt[vec_rowsums == ss, 1:t] %>% unique
  setnames(target_dt, paste0("A1_", 1:t))
  to_match <- dt_use_backup %>% copy()
  setkeyv(to_match, paste0("A1_", 1:t))
  setkeyv(target_dt, paste0("A1_", 1:t))
  temp_vec <- to_match[.(target_dt)][!is.na(age) & get(paste0("C_", t)) == 1, ][[paste0("Y_", t+1)]] %>% table(useNA = "always")
  temp_vec["0"]/sum(temp_vec[1:2])
}

png("./simulation_ZW/MSM_example.png", width = 12, height = 8, units = "in", res = 300)
ll <- 10
expit(test$beta[1] + test$beta[3] * (6:ll) + test$beta[2]*c((6:ll))) %>% plot(x = 6:ll, type = "l", ylim = c(0.05, 0.35), ylab = "Risk", xlab = "t")
expit(test$beta[1] + test$beta[3] * (6:ll) + test$beta[2]*c((6:ll) - 1)) %>% lines(x = 6:ll, col = "blue")
expit(test$beta[1] + test$beta[3] * (6:ll) + test$beta[2]*c((6:ll) - 2)) %>% lines(x = 6:ll, col = "red")
# (1-expit(coef_intercept +coef_a*(6:ll) + coef_b*c((6:ll)))) %>% lines(x = 6:ll, col = "black", lty = 2)
# (1-expit(coef_intercept +coef_a*(6:ll) + coef_b*c((6:ll)-1))) %>% lines(x = 6:ll, col = "blue", lty = 2)
# (1-expit(coef_intercept +coef_a*(6:ll) + coef_b*c((6:ll)-2))) %>% lines(x = 6:ll, col = "red", lty = 2)
(1-sapply(6:ll, function(u) np_truth(u, u) )) %>% lines(x = 6:ll, col = "black", lty = 2)
(1-sapply(6:ll, function(u) np_truth(u, u-1) )) %>% lines(x = 6:ll, col = "blue", lty = 2)
(1-sapply(6:ll, function(u) np_truth(u, u-2) )) %>% lines(x = 6:ll, col = "red", lty = 2)

legend(x = "topleft",          # Position
       legend = c("Estimate", "Target", "Missing 2", "Missing 1", "Full exposure"),  # Legend texts
       lty = c(1, 2, 1, 1, 1),           # Line types
       col = c("black", "black", "red", "blue", "black")
       )

dev.off()
