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
  
  # fake dual traetment nodes
  set.seed(123)
  for (i in 1:10) {
    dt_use[, paste0("A2_", i) := get(paste0("A1_", i))[sample(nrow(dt_use))]]
    setcolorder(dt_use, c(colnames(dt_use)[1:grep(paste0("^A1_", i, "$"), colnames(dt_use))], 
                          paste0("A2_", i))
    )
  }
  
  # fake competing risk; only happens when uncensored and dementia-free with 0.05 probability
  # once dead, all future Y2_t become 1 and Y_t become 0
  set.seed(123)
  for (i in 1:11) {
    # transform all non-dementia to death
    temp_Y <- dt_use[, paste0("Y_", i), with = F]
    temp_Y <- temp_Y + 1
    temp_Y[temp_Y==2] <- 0
    # only keep 0.05 death
    temp_Y[sample(which(temp_Y == 1), round(length(which(temp_Y == 1)) * 0.95))] <- 0

    dt_use[, paste0("Y2_", i) := temp_Y]
    setcolorder(dt_use, c(colnames(dt_use)[1:grep(paste0("^Y_", i, "$"), colnames(dt_use))], 
                          paste0("Y2_", i))
    )
  }
  for (i in 1:11) {
    temp_Y <- dt_use[, paste0("Y2_", i), with = F]
    for (j in i:11) {
      dt_use[which(temp_Y == 1), (paste0("Y2_", j)) := 1] 
    }
  }
  # if any dementia, remove death
  dt_use[which(apply(dt_use[, paste0("Y_", 1:11)] == 1, 1, any)), 
         (paste0("Y2_", 1:11)) := ifelse(get((paste0("Y2_", 1:11))) %>% is.na, NA, 0)
         ]
  # # if dementia, remain that way
  # for (i in 1:11) {
  #   temp_Y <- dt_use[, paste0("Y2_", i), with = F]
  #   for (j in i:11) {
  #     dt_use[which(temp_Y == 1), (paste0("Y2_", j)) := 1] 
  #   }
  # }
  

  # if censored, remove event
  dt_use[which(apply(dt_use[, names(dt_use)[grep("^C", names(dt_use))], with = F] == 0, 1, any)),
         (paste0("Y_", 1:11)) := 0
         ]
  dt_use[which(apply(dt_use[, names(dt_use)[grep("^C", names(dt_use))], with = F] == 0, 1, any)),
         (paste0("Y2_", 1:11)) := 0
         ]
  
  
  C.nodes.index <- grep("^C_", names(dt_use))
  for (i in 1:10) {
    temp_C <- dt_use[, paste0("C_", i), with = F]
    cur.C.index <- grep(paste0("^C_", i, "$"), names(dt_use))
    for (j in (cur.C.index+1):ncol(dt_use)) {
      if (j %in% C.nodes.index) {} else {
        dt_use[which(temp_C == 0), (names(dt_use)[j]) := NA] 
      }
    }
  }
  dt_use[which(apply(dt_use[, paste0("Y_", 1:11)] == 1, 1, any)), 
         (paste0("Y2_", 1:11)) := ifelse(get((paste0("Y2_", 1:11))) %>% is.na, NA, 0)
         ]
  

  data <- as.data.frame(dt_use)
}


det.Q.function <- function(data, current.node, nodes, called.from.estimate.g) {
  Y2.index <- grep("^Y2_", names(data))
  if (length(Y2.index) == 0) stop("ZW: no Y2 node found")
  history.Y2.index <- Y2.index[Y2.index < current.node]
  if (length(history.Y2.index) == 0) return(NULL)
  if (length(history.Y2.index) == 1) {
    is.deterministic <- data[, history.Y2.index] == 1 
  } else {
    is.deterministic <- apply(data[, history.Y2.index] == 1, 1, any) 
  }  # if death before, then remove from fitting
  is.deterministic[is.na(is.deterministic)] <- F
  return(list(is.deterministic=is.deterministic, Q.value=0))
}


# stochastic intervention with dual treatment nodes
SI_function <- function(data, current.node, nodes) {
  Anodes <- nodes$A  # Anodes can include A1_t and A2_t in order
  if (!(current.node %in% Anodes)) return(NULL)   # if not A node, skip
  if (sum(Anodes < current.node) < 2) return(NULL)  # for A1_1 and A2_1, skip
  # if (!(any(Anodes < current.node))) return(NULL)  # if there is no previous A node, then treat as usual
  
  # if t-1 is treated, then deterministically set to be the observed value; (when observed value is not NA)
  # prob1 is the prob of being 1, equal to obs value here: observed 1, set a^d_t = obs value with prob 1; observed 0, set a^d_t = obs value with prob 0
  # if t-1 is not treated, then treat as usual (not deterministic)
  prev.a.node <- max(Anodes[Anodes < current.node])
  prev.prev.a.node <- max(Anodes[Anodes < prev.a.node])
  is.deterministic <- ifelse(!is.na(data[, current.node]), data[, prev.prev.a.node] == 1, F)
  prob1 <- data[, current.node][is.deterministic]
  # prob1 <- prob1[!is.na(prob1)]  # match length(which(det.list$is.deterministic); det.list is the returned output of this function
  return(list(is.deterministic=is.deterministic, prob1=prob1))  
}

# specified the intervened treatment
Anodes = grep("^A[1-2]", node_names)
abar1 <- abar2 <- data[, Anodes]
abar1[, 2] <- abar2[, 1] <- 0
for (k in 2:K) {
  abar1[,k*2-1] = ifelse(abar1[, k*2-3] == 1, abar1[, k*2-1], 1)
  abar1[,k*2] = 0
  abar2[,k*2] = ifelse(abar2[, k*2-2] == 1, abar2[, k*2], 1)
  abar2[,k*2-1] = 0
}
abar1 <- as.matrix(abar1)
abar2 <- as.matrix(abar2)




SL.library <- "glm"
# SL.library <- c("SL.mean", "SL.glm")
# ltmle call
setattr(SL.library, "return.fit", F)
temp_Lnodes <- paste0(rep(c("L_", "Y2_"), K+1), rep(1:(K+1), rep(2, K+1)))
temp_Lnodes <- head(temp_Lnodes, -1)
data <- data[, -ncol(data)]
node_names <- names(data)
test_si <- ltmle(data,  # can use raw, or the processed one like this
                 Anodes = grep("^A[1-2]_", node_names), 
                 Lnodes = temp_Lnodes, 
                 Ynodes = grep("^Y_", node_names), Cnodes = grep("^C_", node_names), 
                 survivalOutcome = T, 
                 abar = list(abar1, abar2),  # the trt of interest under stochastic intervention
                 deterministic.g.function = SI_function, 
                 deterministic.Q.function = det.Q.function, 
                 SL.library = SL.library,  # use one of the libraries above
                 variance.method = "ic",  # no targeting of gn
                 SL.cvControl = list(V = ncores),  # V folds cross validation in super learners
                 estimate.time = F,  # skip runtime est
)
# summary(test_si)
