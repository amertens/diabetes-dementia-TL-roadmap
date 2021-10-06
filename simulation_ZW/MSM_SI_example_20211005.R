getwd()  # under folder diabetes-dementia-TL-roadmap
pkg_dir <- "./simulation_ZW/DK_trip_2021/pkgs/ltmle/cvSL_snow_skipspeedglm_20211005/ltmle-master/R/"  # load ltmle with snow parallel SL; skip speedglm runs
# pkg_dir <- "./simulation_ZW/DK_trip_2021/pkgs/ltmle/cvSL_skipspeedglm_20210915/ltmle-master/R/"  # load ltmle with multicore parallel SL; skip speedglm runs
data_path <- "./simulation_ZW/DK_trip_2021/dt_use_backup_MSM.rds"  # example data

ncores <- 8  # parallel across folds
options(snow.cores = ncores)

library(dplyr)
library(data.table)
# load ltmle
{
  library(Matrix)
  library(matrixStats)
  library(speedglm)
  temp <- list.files(pkg_dir, full.names = T)
  for (i in temp) source(i)
}
library(hal9001)
library(pryr)


# load simulated data
{
  logit <- function(x) log(x/(1-x))
  expit <- function(x) 1/(1+exp(-x))
  coef_intercept <- 3.5
  coef_a <- -0.3  # marginal time coefficient for survival prob
  coef_b <- 0.15  # marginal dose coefficient for survival prob
  
  dt_use_backup <- readRDS(data_path)
  
  set.seed(123)
  dt_use <- dt_use_backup[sample(nrow(dt_use_backup), 20000, T), ]  # target sample size
  K <- 10  # total time points
  dt_use[, first_date_2nd_line := NULL]  # remove index dates; equivalence with censoring process
  node_names <- names(dt_use)
  
  # Anodes = grep("^A1_", node_names)
  # Lnodes = paste0("L_", 1:(K+1))
  # Ynodes = grep("^Y_", node_names)
  # Cnodes = grep("^C_", node_names)
  # all.nodes <- CreateNodes(dt_use, Anodes, Cnodes, Lnodes, Ynodes)
  # data <- CleanData(dt_use %>% as.data.frame, all.nodes, NULL, T)
  data <- as.data.frame(dt_use)
}

# tier 0, no super learner, logistic regression full sample
SL.library <- "glm"
# tier 1, only learners being logistic and mean
SL.library <- c("SL.glm", "SL.mean")
# tier 2, main term lasso
SL.library <- c("SL.glmnet", "SL.glm", "SL.mean")
# tier 3, faster HAL
SL.hal9001.flexible <- function(Y,
                                X,
                                newX = NULL,
                                family = stats::binomial(),
                                obsWeights = rep(1, length(Y)),
                                id = NULL,
                                # max_degree = ifelse(ncol(X) >= 20, 2, 3),
                                max_degree = ifelse(ncol(X) >= 5, 2, 3),
                                smoothness_orders = 1,
                                # num_knots = ifelse(smoothness_orders >= 1, 25, 50),
                                num_knots = if (max_degree == 2) c(25, 10) else c(25, 10, 5),  # ZW: use decreased number of knots for interactions; recommended fast setting in ?fit_hal
                                reduce_basis = 1 / sqrt(length(Y)),
                                lambda = NULL,
                                # bounds = c(0.005, 0.995),
                                ...) {
  
  if (inherits(family, "family")) {
    
  } else if (all(Y %in% c(0, 1, NA))) {
    family <- stats::binomial()
  } else {
    family <- stats::gaussian()
  }
  
  # create matrix version of X and newX for use with hal9001::fit_hal
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.null(newX) & !is.matrix(newX)) newX <- as.matrix(newX)
  
  # fit hal
  hal_fit <- hal9001::fit_hal(  # ZW: use side loaded dev version, with family as obj not character
    X = X, Y = Y, family = family$family,
    fit_control = list(weights = obsWeights), id = id, max_degree = max_degree,
    smoothness_orders = smoothness_orders, 
    num_knots = num_knots,  
    reduce_basis
    = reduce_basis, lambda = lambda
  )
  
  # compute predictions based on `newX` or input `X`
  if (!is.null(newX)) {
    pred <- stats::predict(hal_fit, new_data = newX)
  } else {
    pred <- stats::predict(hal_fit, new_data = X)
  }
  # pred[pred < bounds[1]] <- bounds[1]
  # pred[pred > bounds[2]] <- bounds[2]
  
  # build output object
  fit <- list(object = hal_fit)
  class(fit) <- "SL.hal9001"
  out <- list(pred = pred, fit = fit)
  return(out)
}
environment(SL.hal9001.flexible) <- asNamespace("SuperLearner") 
SL.library <- c("SL.hal9001.flexible", "SL.glm", "SL.mean")
# tier 4, faster HAL with lasso
SL.library <- c("SL.hal9001.flexible", "SL.glmnet", "SL.glm", "SL.mean")
# tier 5, slower HAL with lasso
SL.hal9001.flexible <- function(Y,
                                X,
                                newX = NULL,
                                family = stats::binomial(),
                                obsWeights = rep(1, length(Y)),
                                id = NULL,
                                max_degree = ifelse(ncol(X) >= 20, 2, 3),
                                smoothness_orders = 1,
                                # num_knots = ifelse(smoothness_orders >= 1, 25, 50),
                                num_knots = if (max_degree == 2) c(50, 25) else c(50, 25, 15),  # ZW: use decreased number of knots for interactions; recommended fast setting in ?fit_hal
                                reduce_basis = 1 / sqrt(length(Y)),
                                lambda = NULL,
                                # bounds = c(0.005, 0.995),
                                ...) {
  
  if (inherits(family, "family")) {
    
  } else if (all(Y %in% c(0, 1, NA))) {
    family <- stats::binomial()
  } else {
    family <- stats::gaussian()
  }
  
  # create matrix version of X and newX for use with hal9001::fit_hal
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.null(newX) & !is.matrix(newX)) newX <- as.matrix(newX)
  
  # fit hal
  hal_fit <- hal9001::fit_hal(  # ZW: use side loaded dev version, with family as obj not character
    X = X, Y = Y, family = family$family,
    fit_control = list(weights = obsWeights), id = id, max_degree = max_degree,
    smoothness_orders = smoothness_orders, 
    num_knots = num_knots,  
    reduce_basis
    = reduce_basis, lambda = lambda
  )
  
  # compute predictions based on `newX` or input `X`
  if (!is.null(newX)) {
    pred <- stats::predict(hal_fit, new_data = newX)
  } else {
    pred <- stats::predict(hal_fit, new_data = X)
  }
  # pred[pred < bounds[1]] <- bounds[1]
  # pred[pred > bounds[2]] <- bounds[2]
  
  # build output object
  fit <- list(object = hal_fit)
  class(fit) <- "SL.hal9001"
  out <- list(pred = pred, fit = fit)
  return(out)
}
environment(SL.hal9001.flexible) <- asNamespace("SuperLearner") 
SL.library <- c("SL.hal9001.flexible", "SL.glmnet", "SL.glm", "SL.mean")


# code summary covaraites in MSM models; now include accumulated exposure and time
{
  n <- nrow(data)
  time.points <- 10
  n_pool <- 10
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
      summary.measures[i, , j] <- c(sum(regime.matrix[i, 1:(time.points - n_pool + j)]), time.points - n_pool + j)  # accumulated exposure, time
    }
    # test.treated[i, 1, 1] <- !any(diff(which(regime.matrix[i, ] == 0)) == 1)
    
    ###
    # specify maximum lengths of each interruption and total interruption here
    ###
    test.treated[i, 1, 1] <- all(diff(c(0, which(regime.matrix[i, ] == 1), 11)) <= 2) &  # max length of one interruption: this first number -1; for example, 2-1 = 1 max each interruption
      sum(regime.matrix[i, ] == 0) <=5  # max total interruption: this second number; for example, max total interruption 2
  }
  colnames(test.treated) <- "if.in.group"
  test.treated[, 1, 1] %>% table
  
  summary.measures <- summary.measures[test.treated[, 1, 1], , ]
  for (j in 1:n_pool) {
    colnames(summary.measures[, , j]) <- c("time.on.treatment", "time")
  }
  colnames(summary.measures) <- c("time.on.treatment", "time")
}


SL.library <- "glm"
# ltmleMSM call
options(snow.cores = ncores)
test <- ltmleMSM(data, Anodes = grep("^A1_", node_names), 
                 Lnodes = paste0("L_", 1:11), 
                 Ynodes = grep("^Y_", node_names), 
                 Cnodes = grep("^C_", node_names), 
                 survivalOutcome = T, 
                 SL.library = SL.library, 
                 regimes = regimes[, , test.treated[, 1, 1]],  # only the candidate regimes
                 summary.measures = summary.measures,   # corresponding summary measures for the candidate regimes
                 working.msm = "Y~time.on.treatment + time",
                 variance.method = "ic",  # direct EIC plug in; might be underestimated with positivity and rare outcomes
                 final.Ynodes = paste0("Y_", seq(to = K+1, length.out = n_pool, by = 1)),   # to pool across these nodes
                 msm.weights = "empirical",  # h weights by obs data
                 SL.cvControl = list(V = ncores),  # control CV fold numbers; might be used in paralleled version
                 # gbounds = c(0.05, 0.95),
                 estimate.time = F  # do not run on a n=50 subsample for predicting computation time
)



# stochastic intervention
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
# specified the intervened treatment
Anodes = grep("^A1_", node_names)
abar <- data[, Anodes]
for (k in 2:K) {
  abar[,k] = ifelse(abar[, k-1] == 1, abar[, k], 1)
}
abar <- as.matrix(abar)
setattr(SL.library, "return.fit", F)
test_si <- ltmle(data,  # can use raw, or the processed one like this
                 Anodes = grep("^A1_", node_names), Lnodes = paste0("L_", 1:(K+1)), Ynodes = grep("^Y_", node_names), Cnodes = grep("^C_", node_names), survivalOutcome = T, 
                 abar = abar,  # the trt of interest under stochastic intervention
                 deterministic.g.function = SI_function, 
                 SL.library = SL.library,  # use one of the libraries above
                 variance.method = "ic",  # no targeting of gn
                 SL.cvControl = list(V = ncores),  # V folds cross validation in super learners
                 estimate.time = F,  # skip runtime est
)








png("./simulation_ZW/MSM_SI_example.png", width = 12, height = 8, units = "in", res = 300)
ll <- 10
expit(test$beta[1] + test$beta[3] * (1:ll) + test$beta[2]*c((1:ll))) %>% plot(x = 1:ll, type = "l", ylim = c(0.00, 0.3), ylab = "Risk", xlab = "t")
(1-expit(coef_intercept +coef_a*(1:ll) + coef_b*c((1:ll)))) %>% lines(x = 1:ll, col = "black", lty = 2)
temp_sec <- (1:ll) - 1
temp_sec[temp_sec < 0] <- 0
expit(test$beta[1] + test$beta[3] * (1:ll) + test$beta[2]*temp_sec) %>% lines(x = 1:ll, col = "blue")
(1-expit(coef_intercept +coef_a*(1:ll) + coef_b*temp_sec)) %>% lines(x = 1:ll, col = "blue", lty = 2)
temp_sec <- (1:ll) - 3
temp_sec[temp_sec < 0] <- 0
expit(test$beta[1] + test$beta[3] * (1:ll) + test$beta[2]*temp_sec) %>% lines(x = 1:ll, col = "red")
(1-expit(coef_intercept +coef_a*(1:ll) + coef_b*temp_sec)) %>% lines(x = 1:ll, col = "red", lty = 2)
temp_sec <- (1:ll) - 5
temp_sec[temp_sec < 0] <- 0
expit(test$beta[1] + test$beta[3] * (1:ll) + test$beta[2]*temp_sec) %>% lines(x = 1:ll, col = "green")
(1-expit(coef_intercept +coef_a*(1:ll) + coef_b*temp_sec)) %>% lines(x = 1:ll, col = "green", lty = 2)

legend(x = "topleft",          # Position
       legend = c("Estimate", "Truth", "Skipped 5", "Skipped 3", "Skipped 1", "Full exposure"),  # Legend texts
       lty = c(1, 2, 1, 1, 1),           # Line types
       col = c("black", "black", "green", "red", "blue", "black")
)
points(ll, test_si$estimates[1], pch = 8)
dev.off()