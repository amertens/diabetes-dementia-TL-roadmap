

rm(list=ls())

#devtools::install_github('osofr/stremr')
library(ltmle)
library(stremr)



# Simulated data example
# Load the data:
  
  require("magrittr")
#> Loading required package: magrittr
require("data.table")
#> Loading required package: data.table
require("stremr")
#> Loading required package: stremr

data(OdataNoCENS)
datDT <- as.data.table(OdataNoCENS, key=c("ID", "t"))
#Define some summaries (lags):
datDT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
datDT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]
#Define counterfactual exposures. In this example we define one intervention as always treated and another as never treated. Such intervention can be defined conditionally on other variables (dynamic intervention). Similarly, one can define the intervention as a probability that the counterfactual exposure is 1 at each time-point t (for stochastic interventions).

datDT[, ("TI.set1") := 1L]
datDT[, ("TI.set0") := 0L]
#Import input data into stremr object DataStorageClass and define relevant covariates:
  
  OData <- importData(datDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "N.tminus1"), CENS = "C", TRT = "TI", OUTCOME = "Y.tplus1")
#Once the data has been imported, it is still possible to inspect it and modify it, as shown in this example:
  
get_data(OData)[, ("TI.set0") := 1L]
get_data(OData)[, ("TI.set0") := 0L]
#Regressions for modeling the propensity scores for censoring (CENS) and exposure (TRT). By default, each of these propensity scores is fit with a common model that pools across all available time points (smoothing over all time-points).

gform_CENS <- "C ~ highA1c + lastNat1"
gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
#Stratification, that is, fitting separate models for different time-points, is enabled with logical expressions in arguments stratify_... (see ?fitPropensity). For example, the logical expression below states that we want to fit the censoring mechanism with a separate model for time point 16, while pooling with a common model fit over time-points 0 to 15. Any logical expression can be used to define such stratified modeling. This can be similarly applied to modeling the exposure mechanism (stratify_TRT) and the monitoring mechanism (stratify_MONITOR).

stratify_CENS <- list(C=c("t < 16", "t == 16"))
#Fit the propensity scores for censoring, exposure and monitoring:
  
  OData <- fitPropensity(OData,
                         gform_CENS = gform_CENS,
                         gform_TRT = gform_TRT,
                         stratify_CENS = stratify_CENS)
#Estimate survival based on non-parametric/saturated IPW-MSM (IPTW-ADJUSTED KM):
  
  AKME.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
  survNPMSM(OData) %$%
  estimates
#The result is a data.table that contains the estimates of the counterfactual survival for each time-point, for the treatment regimen TI.set1. In this particular case, the column St.NPMSM contains the survival estimates for IPW-NPMSM and the first row represents the estimated proportion alive at the end of the first cycle / time-point. Note that the column St.KM contains the unadjusted/crude estimates of survival (should be equivalent to standard Kaplan-Meier estimates for most cases).

head(AKME.St.1[],2)
#>    est_name time sum_Y_IPW sum_all_IPAW   ht.NPMSM  St.NPMSM      ht.KM
#> 1:    NPMSM    0 1.6610718     38.13840 0.04355379 0.9564462 0.04733728
#> 2:    NPMSM    1 0.8070748     48.10323 0.01677797 0.9403990 0.01863354
#>        St.KM rule.name
#> 1: 0.9526627   TI.set1
#> 2: 0.9349112   TI.set1
#Estimate survival with bounded IPW:
  
  IPW.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
  directIPW(OData) %$%
  estimates
#As before, the result is a data.table with estimates of the counterfactual survival for each time-point, for the treatment regimen TI.set1, located in column St.directIPW.

head(IPW.St.1[],2)
#>     est_name time sum_Y_IPW  sum_IPW St.directIPW rule.name
#> 1: directIPW    0  9.828827 225.6710    0.9564462   TI.set1
#> 2: directIPW    1 14.841714 308.6067    0.9519073   TI.set1
#Estimate hazard with IPW-MSM then map into survival estimate. Using two regimens and smoothing over two intervals of time-points:
  
wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1")
wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "TI.set0", rule_name = "TI0")
survMSM_res <- survMSM(list(wts.DT.1, wts.DT.0), OData, tbreaks = c(1:8,12,16)-1,)
#In this particular case the output is a little different, with separate survival tables for each regimen. The output of survMSM is hence a list, with one item for each counterfactual treatment regimen considered during the estimation. The actual estimates of survival are located in the column(s) St.MSM. Note that survMSM output also contains the standard error estimates of survival at each time-point in column(s) SE.MSM. Finally, the output table also contains the subject-specific estimates of the influence-curve (influence-function) in column(s) IC.St. These influence function estimates can be used for constructing the confidence intervals of the counterfactual risk-differences for two contrasting treatments (see help for get_RDs function for more information).

head(survMSM_res[["TI0"]][["estimates"]],2)
#>    est_name time      ht.MSM    St.MSM      SE.MSM rule.name
#> 1:      MSM    0 0.004214338 0.9957857 0.002105970       TI0
#> 2:      MSM    1 0.013068730 0.9827720 0.004100295       TI0
#>                                                                       IC.St
#> 1: 0.004543242,0.004543242,0.004543242,0.004543242,0.004543242,0.004543242,
#> 2: 0.004483868,0.016683119,0.016797415,0.017900770,0.017900770,0.017900770,
head(survMSM_res[["TI1"]][["estimates"]],2)
#>    est_name time     ht.MSM    St.MSM     SE.MSM rule.name        IC.St
#> 1:      MSM    0 0.04355379 0.9564462 0.01521910       TI1 0,0,0,0,0,0,
#> 2:      MSM    1 0.01677797 0.9403990 0.01786105       TI1 0,0,0,0,0,0,
#Longitudinal GCOMP (G-formula) and TMLE
#Define time-points of interest, regression formulas and software to be used for fitting the sequential outcome models:
  
  tvals <- c(0:10)
Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(tvals)+1))
#To run iterative means substitution estimator (G-Computation), where all at risk observations are pooled for fitting each outcome regression (Q-regression):
  
  gcomp_est <- fit_GCOMP(OData, tvals = tvals, intervened_TRT = "TI.set1", Qforms = Qforms)
#The output table of fit_GCOMP contains the following information, with the column St.GCOMP containing the survival estimates for each time period:
  
  head(gcomp_est$estimates[],2)
#>    est_name time  St.GCOMP St.TMLE   type    cum.inc              IC.St
#> 1:    GCOMP    0 0.9837583      NA pooled 0.01624168 NA,NA,NA,NA,NA,NA,
#> 2:    GCOMP    1 0.9699022      NA pooled 0.03009778 NA,NA,NA,NA,NA,NA,
#>    fW_fit rule.name
#> 1:   NULL   TI.set1
#> 2:   NULL   TI.set1
#To run the longitudinal long format Targeted Minimum-Loss Estimation (TMLE), stratified by rule-followers for fitting each outcome regression (Q-regression):
  
  tmle_est <- fit_TMLE(OData, tvals = tvals, intervened_TRT = "TI.set1", Qforms = Qforms)
#> GLM TMLE update cannot be performed since the outcomes (Y) are either all 0 or all 1, setting epsilon to 0
#> GLM TMLE update cannot be performed since the outcomes (Y) are either all 0 or all 1, setting epsilon to 0

#The output table of fit_TMLE contains the following information, with the column St.TMLE containing the survival estimates for each time period. In addition, the column SE.TMLE contains the standard error estimates and the column and the column IC.St contains the subject-specific estimates of the efficient influence curve. The letter estimates are useful for constructing the confidence intervals of risk differences for two contrasting treatments (see help for get_RDs function for more information).

head(tmle_est$estimates[],2)
#>    est_name time St.GCOMP   St.TMLE   type    cum.inc     SE.TMLE
#> 1:     TMLE    0       NA 0.9839271 pooled 0.01607286 0.003449949
#> 2:     TMLE    1       NA 0.9707676 pooled 0.02923243 0.004492235
#>                                                                             IC.St
#> 1: -0.007292922,-0.007292922,-0.010190141,-0.007292922,-0.007292922,-0.007292922,
#> 2: -0.009469707,-0.009469707,-0.010503891,-0.009469707,-0.009469707,-0.009469707,
#>    fW_fit rule.name
#> 1:   NULL   TI.set1
#> 2:   NULL   TI.set1
#To parallelize estimation over several time-points (tvals) for either GCOMP or TMLE use argument parallel = TRUE:
  
  require("doParallel")
registerDoParallel(cores = parallel::detectCores())
tmle_est <- fit_TMLE(OData, tvals = tvals, intervened_TRT = "TI.set1", Qforms = Qforms, parallel = TRUE)
#Data-adaptive estimation, cross-validation and Super Learning
#Nuisance parameters can be modeled with any machine learning algorithm supported by sl3 R package. For example, for GLMs use learner Lrnr_glm_fast, for xgboost use learner Lrnr_xgboost, for h2o GLM learner use Lrnr_h2o_glm, for any other ML algorithm implemented in h2o use Lrnr_h2o_grid$new(algorithm = "algo_name"), for glmnet use learner Lrnr_glmnet. All together, these learners provide access to a wide variety of ML algorithms. To name a few: GLM, Regularized GLM, Distributed Random Forest (RF), Extreme Gradient Boosting (GBM) and Deep Neural Nets.

#Model selection can be performed via V-fold cross-validation or random validation splits and model stacking and Super Learner combination can be accomplished by using the learner Lrnr_sl and specifying the meta-learner (e.g., Lrnr_solnp). In the example below we define a Super Learner ensemble consisting of several learning algorithms.
#First, we define sl3 learners for for xgboost, two types of GLMs and glmnet. Then we will stack these learners into a single learner called Stack:
  
  library("sl3")
lrn_xgb <- Lrnr_xgboost$new(nrounds = 5)
lrn_glm <- Lrnr_glm_fast$new()
lrn_glm2 <- Lrnr_glm_fast$new(covariates = c("CVD"))
lrn_glmnet <- Lrnr_glmnet$new(nlambda = 5, family = "binomial")
## Stack the above candidates:
lrn_stack <- Stack$new(lrn_xgb, lrn_glm, lrn_glm2, lrn_glmnet)
#Next, we will define a Super Learner on the above defined stack, by feeding the stack into the Lrnr_sl object and then specifying the meta-learner that will find the optimal convex combination of the learners in a stack (Lrnr_solnp):
  
  lrn_sl <- Lrnr_sl$new(learners = lrn_stack, metalearner = Lrnr_solnp$new())
#We will now use stremr to estimate the exposure / treatment propensity model with the above defined Super Learner (lrn_sl):
  
  OData <- fitPropensity(OData,
                         gform_CENS = gform_CENS,
                         gform_TRT = gform_TRT,
                         models_TRT = lrn_sl,
                         stratify_CENS = stratify_CENS)




library(dplyr)
library(data.table)
library(Matrix)
library(matrixStats)
library(speedglm)
temp <- list.files("/home/leo42k/Downloads/ltmle-master/R", full.names = T)
for (i in temp) source(i)

library(hal9001)
library(pryr)

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

# saveRDS(dt_tmle, "./simulation_ZW/dt_tmle_202108.rds")
dt_tmle <- readRDS("./simulation_ZW/dt_tmle_202108.rds")
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
prob_censoring <- 0.05
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



# first time bin t where censoring status changes
censored_t <- dt_use[, paste0("C_", 1:10)] %>% apply(1, function(eachRow) first(which(eachRow == 0))) %>% lapply(function(x) ifelse(length(x) == 0, 0, x)) %>% unlist
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

# # first time bin t where event status changes
# event_t <- dt_use[, paste0("Y_", 1:11)] %>% apply(1, function(eachRow) first(which(eachRow == 1))) %>% lapply(function(x) ifelse(length(x) == 0, 0, x)) %>% unlist
# 
# 
# if_censored <- censored_t > 0
# if_event <- event_t > 0
# # table(if_censored&if_event)  # we used Yt=1 after C_{t-1} censored (not in the risk set)
# censored_t[!(if_censored | if_event)] <- 11  # everyone with no events is censored by hypothetical trial administration
# if_censored <- censored_t > 0
# # one can be simulated as both censored and having event; but ltmle only takes which ever comes in first
# # when deciding the fake index date, use the first one
# max.date <- dt_use$first_date_2nd_line %>% max
# fake_index_dates <- (sapply(1:nrow(dt_use), function(x) {
#   if (if_censored[x] & if_event[x]) {
#     if (cencored_t[x] < event_t[x]) {
#       round(max.date - (cencored_t[x]-1) * 365.25/2) %>% return
#     } else 
#       round(max.date - 100 - (event_t[x]-1) * 365.25/2) %>% return
#   } else if (if_censored[x]) {
#     round(max.date - (cencored_t[x]-1) * 365.25/2) %>% return
#   } else {
#     round(max.date - 100 - (event_t[x]-1) * 365.25/2) %>% return
#   }
# }) + sample(30:70, nrow(dt_use), T)) %>% as.Date(origin = "1970-01-01")  # allow actual dates to be somewhere in between, so long as the status (event/censored) at the edge of an interval won't change
# dt_use[, first_date_2nd_line := fake_index_dates]


# get_censoring <- function(day_first, day_last = max.date, tau = K) {
#   temp_t <- floor((day_last - day_first - 0.0001) / (365.25/2)) + 1 # decide max.date is in which (, ] interval, then it means censoring status changes at the end of this interval
#   if (temp_t > tau) return(rep(1, tau)) else return(c(rep(1, temp_t), rep(0, tau - temp_t)))
# }
# temp_mat <- dt_use$first_date_2nd_line %>% sapply(get_censoring)
# temp_mat <- t(temp_mat)
# 
# 
# temp <- floor((max.date - dt_use$first_date_2nd_line - 0.0001) / (365.25/2)) + 1 # decide max.date is in which (, ] interval, then it means censoring status changes at the end of this interval
# identical(temp[if_censored & !if_event], censored_t[if_censored & !if_event])
# which(temp[if_censored & !if_event] != censored_t[if_censored & !if_event])
# which(temp[if_censored] != censored_t[if_censored])
# temp[if_censored][6]
# temp[if_censored & !if_event][4]
# censored_t[if_censored & !if_event][4]
# 
# censored_t[if_censored][6]
# event_t[if_censored][6]
# 
# (max.date - dt_use$first_date_2nd_line - 0.0001) 




dt_use_backup <- dt_use %>% copy


# dt_use <- dt_use_backup[sample(nrow(dt_use_backup), 2000, T), ]
# 
# start_time <- Sys.time()
# ss <- mem_change(
#   test <- ltmle(dt_use, Anodes = grep("^A1_", node_names), Lnodes = grep("^L_", node_names), Ynodes = grep("^Y_", node_names), Cnodes = grep("^C_", node_names), survivalOutcome = T, 
#                 abar = rep(1, 10))
# )
# end_time <- Sys.time()
# ss
# test
# difftime(end_time, start_time, units = "mins")



set.seed(135)
dt_use <- dt_use_backup[sample(nrow(dt_use_backup), 10000, T), ]




#like seq, but returns integer(0) if from > to   (always increments by 1)
sseq <- function(from, to) {
  if (from > to) return(integer(0))
  seq(from, to)
}

SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}
NodeToIndex <- function(data, node) {
  if (is.numeric(node) || is.null(node)) return(node)
  if (! is.character(node)) stop("nodes must be numeric, character, or NULL")
  index <- match(node, names(data))
  if (anyNA(index)) {
    stop(paste("\nnamed node(s) not found:", node[is.na(index)]))
  }
  return(index)
}
CreateNodes <- function(data, Anodes, Cnodes, Lnodes, Ynodes) {  
  Anodes <- NodeToIndex(data, Anodes)
  Cnodes <- NodeToIndex(data, Cnodes)
  Lnodes <- NodeToIndex(data, Lnodes)
  Ynodes <- NodeToIndex(data, Ynodes)
  nodes <- SuppressGivenWarnings(list(A=Anodes, C=Cnodes, L=Lnodes, Y=Ynodes, AC=sort(c(Anodes, Cnodes))), "is.na() applied to non-(list or vector) of type 'NULL'") #suppress warnings if no A/C nodes
  nodes$baseline <- sseq(1, min(c(nodes$A, nodes$L, nodes$C, nodes$Y)) - 1)
  nodes$LY <- CreateLYNodes(data, nodes, check.Qform=FALSE)
  return(nodes)
}

# Get the LY nodes but don't include "blocks" of L/Y nodes uninterrupted by A/C nodes
CreateLYNodes <- function(data, nodes, check.Qform, Qform) {
  LYnodes <- sort(c(nodes$L, nodes$Y))
  SuppressGivenWarnings(nodes.to.remove <- LYnodes[LYnodes < min(nodes$AC)], "no non-missing arguments to min; returning Inf") #no warning if no AC nodes
  
  #if there are no A/C nodes between two or more LY nodes, only the first LY node in the block is considered an LY node
  if (length(LYnodes) > 1) {
    for (i in 1:(length(LYnodes) - 1)) {
      cur.node <- LYnodes[i]
      next.node <- LYnodes[i + 1]
      if (! any(cur.node:next.node %in% nodes$AC)) {
        nodes.to.remove <- c(nodes.to.remove, next.node)
      }
    }
  }
  new.LYnodes <- setdiff(LYnodes, nodes.to.remove)
  if (check.Qform) {
    removed.Qform.index <- NULL
    for (i in nodes.to.remove) {
      index <- which(names(Qform) == names(data)[i])
      if (length(index) > 0) {
        removed.Qform.index <- c(removed.Qform.index, index)
      }
    }
    if (!is.null(removed.Qform.index)) {
      message("L/Y nodes (after removing blocks)  : ", paste(names(data)[new.LYnodes], collapse = " "), "\n")
      message("Qform names                        : ", paste(names(Qform), collapse = " "), "\n")
      message(paste("The following nodes are not being considered as L/Y nodes because they are part of a block\nof L/Y nodes. They are being dropped from Qform:\n"), paste(names(Qform)[removed.Qform.index], "\n", collapse=" "))
      Qform <- Qform[-removed.Qform.index]
    }
    return(list(LYnodes=new.LYnodes, Qform=Qform))
  }
  return(new.LYnodes)
}
all.nodes <- CreateNodes(dt_use, Anodes = paste0("A1_", 1:10), Cnodes = paste0("C_", 1:10), Lnodes = paste0("L_", 1:11), Ynodes = paste0("Y_", 1:11))
Cnodes <- all.nodes$C
current.node <- 20
current_t <- which(Cnodes == current.node)  # this is the censored t
temp_t <- floor((max.date - dt_use[[4]]) / (365.25/2)) + 1

data.frame(p = ifelse(temp_t <= current_t, 1, 0), 
           obs = (dt_use[["C_4"]] == "censored")*1
) %>% apply(1, function(u) u[1] == u[2]) %>% table



administrativeCensoring <- function(data, current.node, nodes, last_day = max.date, which_baseline_date = 4) {
  Cnodes <- nodes$C
  if (!(current.node %in% Cnodes)) 
    return(NULL)
  current_t <- which(Cnodes == current.node)  # this is the censored t
  temp_t <- floor((last_day - data[[which_baseline_date]]) / (365.25/2)) + 1
  is.deterministic <- rep(T, nrow(data))
  return(list(is.deterministic = is.deterministic, prob1 = ifelse(temp_t <= current_t, 0, 1)))
}

# try fake Y_11 ~ accumulated time relationship
table(dt_use$Y_10)
table(dt_use$Y_11)
dt_use$C_11

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
  test.treated[i, 1, 1] <- all(diff(c(0, which(regime.matrix[i, ] == 1), 11)) <= 2)  # max this number -1 interruption
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

# any(diff(which(c(0, 0, 0, 0, 0) == 0)) == 1)
# any(diff(which(c(1, 0, 1, 0, 0) == 0)) == 1)
# any(diff(which(c(1, 0, 1, 1, 1) == 0)) == 1)
# any(diff(which(c(1, 0, 1, 0, 1) == 0)) == 1)

for (i in 1:10) {
  dt_use[, ':='(paste0("C_", i), BinaryToCensoring(is.uncensored = paste0("C_", i) %>% get))]
}


# options(mc.cores = 8)

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
                     variance.method = "ic",
                     # final.Ynodes = "Y_11",
                     final.Ynodes = paste0("Y_", seq(to = K+1, length.out = n_pool, by = 1)), 
                     msm.weights = "empirical", 
                     # SL.cvControl = list(V = 8)
                     deterministic.g.function = administrativeCensoring,
                     estimate.time = F
                     
    )
  )
  end_time <- Sys.time()
}
ss
test
difftime(end_time, start_time, units = "mins")

ll <- 10
data.frame(truth = 1-expit(coef_intercept +coef_a*(ll) + coef_b*c(6:ll)), 
           est = expit(test$beta[1] + test$beta[3] * ll + test$beta[2]*c(6:ll))
)

test %>% summary



dt_use_backup$Y_11[rowSums(dt_use[, paste0("A1_", 1:10)]) == 10 & dt_use$C_10 == 1] %>% table(useNA = "always")


x <- seq(5, 10, by = 0.1)
plot(x, (1 - expit(coef_intercept + coef_a * 10 + coef_b * x)), type = "l", ylab = "Risk", xlab = "Accumulated Exposure")
lines(x, expit(summary(test)[[1]][1, 1] + summary(test)[[1]][2, 3] * x), lty = 2)
lines(x, expit(summary(test)[[1]][1, 1] + summary(test)[[1]][2, 4] * x), lty = 2)
# lines(x, expit(summary(test)[[1]][1, 4] + summary(test)[[1]][2, 4] * x), col = "red")
# lines(x, expit(summary(test)[[1]][1, 4] + summary(test)[[1]][2, 3] * x), col = "blue")

t = 10
ss = 10

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

sapply(5:10, function(u) np_truth(10, u))



ll <- 10
data.frame(truth = 1-expit(coef_intercept +coef_a*(ll) + coef_b*c(6:ll)), 
           est = expit(test$beta[1] + test$beta[3] * ll + test$beta[2]*c(6:ll)), 
           np = 1 - sapply(6:ll, function(u) np_truth(ll, u))
)
