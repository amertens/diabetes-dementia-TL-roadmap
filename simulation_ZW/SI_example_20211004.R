getwd()  # under folder diabetes-dementia-TL-roadmap
pkg_dir <- "./simulation_ZW/DK_trip_2021/pkgs/ltmle/cvSL_snow_skipspeedglm_20211005/ltmle-master/R/"  
data_path <- "./simulation_ZW/DK_trip_2021/dt_use_backup.rds"

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

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))
coef_intercept <- 3.5
coef_a <- -0.3
coef_b <- 0.1
# expit(coef_intercept + coef_a*(10) + coef_b*(0:10)) %>% plot
# expit(coef_intercept + coef_a*(1:10) + coef_b*(0)) %>% plot(type = "l")
# expit(coef_intercept + coef_a*(1:10) + coef_b*(10)) %>% lines(col = "red")

dt_use_backup <- readRDS(data_path)

set.seed(000)
dt_use <- dt_use_backup[sample(nrow(dt_use_backup), 1000, T), ]  # target sample size
K <- 10

# if(any(dt_use$age > 1)) {
#   dt_use[, age := age/100]  
# }
dt_use[, first_date_2nd_line := NULL]
node_names <- names(dt_use)
node_names <- node_names[!(node_names %in% c("first_date_2nd_line"))]

for (i in 1:10) {
  dt_use[, ':='(paste0("C_", i), BinaryToCensoring(is.uncensored = paste0("C_", i) %>% get))]
}



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

environment(SL.hal9001.flexible) <- asNamespace("SuperLearner")  # 

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

options(snow.cores = ncores)
SL.library <- c(
  "SL.hal9001.flexible"
  , "SL.mean"
  , "SL.glm"
)
setattr(SL.library, "return.fit", F)
attr(SL.library, "return.fit") == TRUE
{
  start_time <- Sys.time()
  ss <- mem_change(
    test <- ltmle(data, Anodes = grep("^A1_", node_names), Lnodes = paste0("L_", 1:(K+1)), Ynodes = grep("^Y_", node_names), Cnodes = grep("^C_", node_names), survivalOutcome = T, 
                  abar = abar, 
                  deterministic.g.function = SI_function, 
                  SL.library = SL.library,
                  variance.method = "ic",
                  SL.cvControl = list(V = ncores), 
                  estimate.time = F,
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
