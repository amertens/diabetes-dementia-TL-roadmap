getwd()  # under folder diabetes-dementia-TL-roadmap
pkg_dir <- "./simulation_ZW/mcSuperLearner_example/pkgs/ltmle/cvSL_skipspeedglm_20210915/ltmle-master/R"  
dt_tmle_path <- "./simulation_ZW/parallel/dt_tmle_202108.rds"

ncores <- 8

library(dplyr)
library(data.table)

# library(ltmle)

library(Matrix)
library(matrixStats)
library(speedglm)
temp <- list.files(pkg_dir, full.names = T)
for (i in temp) source(i)

library(hal9001)
library(pryr)
# library(nnls)

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))
# logit(p_t) = 100 - t - accum - 0.2 t*accum
# p_t = \prod_{j=1}^t s_j
coef_intercept <- 3.5
coef_a <- -0.3
coef_b <- 0.1

expit(coef_intercept + coef_a*(10) + coef_b*(0:10)) %>% plot
expit(coef_intercept + coef_a*(1:10) + coef_b*(0)) %>% plot(type = "l")
expit(coef_intercept + coef_a*(1:10) + coef_b*(10)) %>% lines(col = "red")

dt_tmle <- readRDS(dt_tmle_path)
set.seed(200)
dt_tmle <- dt_tmle[c(sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 10), 100000, T), sample(nrow(dt_tmle), 100000, T)), ]

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
max.date <- as.Date("2016-12-31")


# ltmle use 1 as event
for (x in grep("Y_", node_names)) dt_use[, (node_names[x]) := 1 - get(node_names[x])]


set.seed(123)
# censoring is the within interval change
prob_censoring <- 0.1
dt_use[, C_1 := rbinom(nrow(dt_use), 1, 1-prob_censoring)]
for (i in 2:10) {
  risk_set <- dt_use[[paste0("C_", i-1)]] == 1
  temp_input <- sapply(rep(1-prob_censoring, sum(risk_set)), function(each_p) rbinom(1, 1, each_p))
  dt_use[risk_set, ':=' (paste0("C_", i), temp_input)]
  dt_use[!risk_set, ':=' (paste0("C_", i), 0)]
}
# 0.95 
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
# make Lt nodes fully random too
for(i in 1:(K+1)) {
  dt_use[, paste0("L_", i) := sample(0:1, nrow(dt_use), T)]
}



# first time bin t where censoring status changes
censored_t <- dt_use[, paste0("C_", 1:10)] %>% apply(1, function(eachRow) first(which(eachRow == 0))) %>% lapply(function(x) ifelse(length(x) == 0, 0, 
                                                                                                                                    ifelse(is.na(x), 0, x))) %>% unlist
censored_t[censored_t==0] <- 11
fake_index_dates <- max.date - (censored_t - 1) * (365.25/2) - sample(10:170, nrow(dt_use), T)

# check fake index date inferred Ct process
temp <- floor((max.date - fake_index_dates) / (365.25/2)) + 1 # decide max.date is in which interval, then it means censoring status changes within this interval
table(temp)
censored_t %>% table

dt_use[, first_date_2nd_line := fake_index_dates]

expit(3.5 + coef_a*K + coef_b*K)
temp_count <- dt_use$Y_11[rowSums(dt_use[, paste0("A1_", 1:10)]) == 10 & dt_use$C_10 == 1] %>% table(useNA = "always")
temp_count[1]/sum(temp_count)

rowSums(dt_use[, paste0("A1_", 1:10)]) %>% table

dt_use_backup <- dt_use %>% copy


set.seed(000)
dt_use <- dt_use_backup[sample(nrow(dt_use_backup), 1000, T), ]  # target sample size
# rm(dt_use_backup)
# rm(dt_tmle)

# if(any(dt_use$age > 1)) {
#   dt_use[, age := age/100]  
# }
dt_use[, first_date_2nd_line := NULL]
node_names <- node_names[!(node_names %in% c("first_date_2nd_line"))]

for (i in 1:10) {
  dt_use[, ':='(paste0("C_", i), BinaryToCensoring(is.uncensored = paste0("C_", i) %>% get))]
}


# K <- 5
# uu <- grep(paste0("^Y_", K+1, "$"), node_names)
# node_names <- node_names[1:uu]
# dt_use <- dt_use[, 1:uu]

SL.hal9001 <- function(Y,
                       X,
                       newX = NULL,
                       family = stats::binomial(),
                       obsWeights = rep(1, length(Y)),
                       id = NULL,
                       max_degree = ifelse(ncol(X) >= 20, 2, 3),
                       smoothness_orders = 1,
                       num_knots = ifelse(smoothness_orders >= 1, 25, 50),
                       reduce_basis = 1 / sqrt(length(Y)),
                       lambda = NULL,
                       ...) {
  
  # create matrix version of X and newX for use with hal9001::fit_hal
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.null(newX) & !is.matrix(newX)) newX <- as.matrix(newX)
  
  # fit hal
  hal_fit <- fit_hal(
    X = X, Y = Y, family = family$family,
    fit_control = list(weights = obsWeights), id = id, max_degree = max_degree,
    smoothness_orders = smoothness_orders, num_knots = num_knots, reduce_basis
    = reduce_basis, lambda = lambda
  )
  
  # compute predictions based on `newX` or input `X`
  if (!is.null(newX)) {
    pred <- stats::predict(hal_fit, new_data = newX)
  } else {
    pred <- stats::predict(hal_fit, new_data = X)
  }
  
  # build output object
  fit <- list(object = hal_fit)
  class(fit) <- "SL.hal9001"
  out <- list(pred = pred, fit = fit)
  return(out)
}



# deterministic censoring
SI_function <- function(data, current.node, nodes) {
  # if not A node, skip
  # if there is no previous node, then treat as usual
  Anodes <- nodes$A
  if (!(current.node %in% Anodes)) return(NULL) 
  if (!(any(Anodes < current.node))) return(NULL) 

  # if t-1 is treated, then deterministically set to be the observed value; (when observed value is not NA)
  # prob1 is the prob of being 1, equal to obs value here: observed 1, set a^d_t = obs value with prob 1; observed 0, set a^d_t = obs value with prob 0
  # if t-1 is not treated, then treat as usual (not deterministic)
  prev.a.node <- max(Anodes[Anodes < current.node])
  is.deterministic <- ifelse(!is.na(data[, current.node]), data[, prev.a.node] == 1, F)
  prob1 <- data[, current.node][is.deterministic]
  # prob1 <- prob1[!is.na(prob1)]  # match length(which(det.list$is.deterministic); det.list is the returned output of this function
  return(list(is.deterministic=is.deterministic, prob1=prob1))  
}

options(mc.cores = ncores)

Anodes = grep("^A1_", node_names)
Lnodes = paste0("L_", 1:(K+1))
Ynodes = grep("^Y_", node_names)
Cnodes = grep("^C_", node_names)
all.nodes <- CreateNodes(dt_use, Anodes, Cnodes, Lnodes, Ynodes)
data <- CleanData(dt_use %>% as.data.frame, all.nodes, NULL, T)

Anodes <- paste0("A1_", 1:K)
abar <- data[, Anodes]
for (k in 2:K) {
  abar[,k] = ifelse(abar[, k-1] == 1, abar[, k], 1)
}
abar <- as.matrix(abar)

{
  start_time <- Sys.time()
  ss <- mem_change(
    test <- ltmle(data, Anodes = grep("^A1_", node_names), Lnodes = paste0("L_", 1:(K+1)), Ynodes = grep("^Y_", node_names), Cnodes = grep("^C_", node_names), survivalOutcome = T, 
                  abar = abar, 
                  deterministic.g.function = SI_function, 
                  SL.library = c(
                    "SL.hal9001",
                                 "SL.mean", "SL.glm"),
                  variance.method = "ic",
                  SL.cvControl = list(V = ncores), 
                  estimate.time = F
    )
  )
  end_time <- Sys.time()
}
ss
test
difftime(end_time, start_time, units = "mins")

temp_result <- list(mem = ss, 
                    time = difftime(end_time, start_time, units = "mins"), 
                    est = test)

# temp_output
# temp_output

# file end
