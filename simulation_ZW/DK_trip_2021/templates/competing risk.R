getwd()  # under folder diabetes-dementia-TL-roadmap
pkg_dir <- "./simulation_ZW/DK_trip_2021/pkgs/ltmle/cvSL_snow_skipspeedglm_20211005/ltmle-master/R/"  # load ltmle with snow parallel SL; skip speedglm runs
# pkg_dir <- "./simulation_ZW/DK_trip_2021/pkgs/ltmle/cvSL_skipspeedglm_20210915/ltmle-master/R/"  # load ltmle with multicore parallel SL; skip speedglm runs
data_path <- "./simulation_ZW/DK_trip_2021/dt_use_backup_MSM.rds"  # example data

ncores <- 8  # parallel across folds
options(snow.cores = ncores)

library(dplyr)
library(data.table)
# load ltmle functions
{
  library(Matrix)
  library(matrixStats)
  library(speedglm)
  temp <- list.files(pkg_dir, full.names = T)
  for (i in temp) source(i)
}
library(ltmle)
library(hal9001)
library(pryr)

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))
coef_intercept <- 3.5
coef_a <- -0.3  # marginal time coefficient for survival prob
coef_b <- 0.15  # marginal dose coefficient for survival prob


{
  dt_tmle <- readRDS("./simulation_ZW/dt_tmle_202108.rds")
  
  set.seed(200)
  dt_tmle <- dt_tmle[sample(nrow(dt_tmle), 100000, T), ]
  for (i in 1:10) {
    dt_tmle[, paste0("A1_", i) := sample(0:1, nrow(dt_tmle), replace = T)]
  }
  # add more "treated" subjects to make positivity issue less severe for now
  dt_tmle <- dt_tmle[c(sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 10), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 9), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 8), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 7), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 6), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 5), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 4), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 3), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 2), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 1), 10000, T), 
                       sample(which((dt_tmle[, paste0("A1_", 1:10)] %>% rowSums) == 0), 10000, T), 
                       sample(nrow(dt_tmle), 10000, T)), ]
  
  
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
  prob_censoring <- 0.01
  dt_use[, C_1 := rbinom(nrow(dt_use), 1, 1-prob_censoring)]
  dt_use[which(dt_use$C_1 == 0), (paste0("C_", 2:10)) := 0]
  for (i in 2:10) {
    risk_set <- dt_use[[paste0("C_", i-1)]] == 1
    temp_input <- sapply(rep(1-prob_censoring, sum(risk_set)), function(each_p) rbinom(1, 1, each_p))
    dt_use[risk_set, ':=' (paste0("C_", i), temp_input)]
    dt_use[!risk_set, ':=' (paste0("C_", i), 0)]
    dt_use[which(dt_use[[paste0("C_", i)]] == 0), (paste0("C_", i:10)) := 0]
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
}

# load simulated data
{
  
  # dt_use_backup <- readRDS(data_path)
  
  # set.seed(123)
  # dt_use <- dt_use_backup[sample(nrow(dt_use_backup), 200000, T), ]  # target sample size
  K <- 10  # total time points
  dt_use[, first_date_2nd_line := NULL]  # remove index dates; equivalence with censoring process
  node_names <- names(dt_use)
  
  # fake dual treatment nodes
  if_any_A1 <- apply(dt_use[, paste0("A1_", 1:10), with = F] == 1, 1, any)
  set.seed(123)
  for (i in 1:10) {
    temp_A <- dt_use[, paste0("A1_", i), with = F][[1]]
    temp_A2 <- rep(0, nrow(dt_use))
    temp_A2[temp_A == 0] <- rbinom(sum(temp_A == 0, na.rm = T), 1, 0.5)
    temp_A2[if_any_A1] <- 0
    dt_use[, paste0("A2_", i) := temp_A2
             # rbinom(nrow(dt_use), 1, 0.05)
           ]
    setcolorder(dt_use, c(colnames(dt_use)[1:grep(paste0("^A1_", i, "$"), colnames(dt_use))], 
                          paste0("A2_", i))
    )
  }
  
  # specified the intervened treatment
  Anodes_backup = grep("^A[1-2]", names(dt_use))
  abar1_backup <- abar2_backup <- dt_use[, Anodes_backup, with = F] %>% as.data.frame
  abar1_backup[, 2] <- abar2_backup[, 1] <- 0
  for (k in 2:K) {
    abar1_backup[,k*2-1] = ifelse(abar1_backup[, k*2-3] == 1, abar1_backup[, k*2-1], 1)
    abar1_backup[,k*2] = 0
    abar2_backup[,k*2] = ifelse(abar2_backup[, k*2-2] == 1, abar2_backup[, k*2], 1)
    abar2_backup[,k*2-1] = 0
  }
  abar1_backup <- as.matrix(abar1_backup)
  abar2_backup <- as.matrix(abar2_backup)
  s_at_K <- abar1_backup[, paste0("A1_", 1:10)] %>% rowSums()
  1 - mean(expit(coef_intercept+coef_a*10+coef_b*s_at_K))
  
  
  # fake competing risk; only happens when uncensored and completely dementia-free with 0.05 probability
  # once dead, all future Y2_t become 1 and Y_t become 0
  dementia_process <- dt_use[, paste0("Y_", 1:11), with = F]
  if_any_dementia <- apply(dementia_process == 1, 1, any)
  if_any_dementia[is.na(if_any_dementia)] <- F
  set.seed(123)
  for (i in 1:11) {
    # transform all non-dementia to death
    temp_Y <- dt_use[, paste0("Y_", i), with = F]
    temp_Y <- temp_Y + 1
    temp_Y[temp_Y==2] <- 0
    # only keep 0.05 or less death
    temp_Y[sample(which(temp_Y == 1), round(length(which(temp_Y == 1)) * 0.95))] <- 0
    # remove death with later dementia
    temp_Y[if_any_dementia] <- 0

    dt_use[, paste0("Y2_", i) := temp_Y]
    setcolorder(dt_use, c(colnames(dt_use)[1:grep(paste0("^Y_", i, "$"), colnames(dt_use))], 
                          paste0("Y2_", i))
    )
  }
  # if dementia or death, remain that way
  for (i in 1:11) {
    temp_Y <- dt_use[, paste0("Y2_", i), with = F]
    for (j in i:11) {
      dt_use[which(temp_Y == 1), (paste0("Y2_", j)) := 1] 
    }
  }
  for (i in 1:11) {
    temp_Y <- dt_use[, paste0("Y_", i), with = F]
    for (j in i:11) {
      dt_use[which(temp_Y == 1), (paste0("Y_", j)) := 1] 
    }
  }
  # # if any dementia, remove death
  # dt_use[which(apply(dt_use[, paste0("Y_", 1:11)] == 1, 1, any)), 
  #        (paste0("Y2_", 1:11)) := ifelse(get((paste0("Y2_", 1:11))) %>% is.na, NA, 0)
  #        ]


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
  # dt_use[which(apply(dt_use[, paste0("Y_", 1:11)] == 1, 1, any)), 
  #        (paste0("Y2_", 1:11)) := ifelse(get((paste0("Y2_", 1:11))) %>% is.na, NA, 0)
  #        ]
  
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




RNGkind("L'Ecuyer-CMRG")
set.seed(123)

library(parallel)
temp_result <- mclapply(X = 1:500, mc.cores = 8, FUN = function(aaa) {
  data <- as.data.frame(dt_use[sample(nrow(dt_use), 50000, replace = T), ])
  # specified the intervened treatment
  Anodes = grep("^A[1-2]", colnames(data))
  abar1 <- abar2 <- data[, Anodes]
  abar1[, 2] <- abar2[, 1] <- 0
  for (k in 2:K) {
    abar1[,k*2-1] = ifelse(abar1[, k*2-3] == 1, abar1[, k*2-1], 1)
    abar1[,k*2] = 0
    abar2[,k*2] = ifelse(abar2[, k*2-2] == 1, abar2[, k*2], 1)
    abar2[,k*2-1] = 0
  }
  for (k in 1:K){
    abar1[is.na(abar1[, 2*k-1]), 2*k] <- NA
    abar2[is.na(abar2[, 2*k]), 2*k-1] <- NA
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
  
  CheckData(data)
  
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
  summary(test_si)[[2]][[1]]$estimate %>% return
})

temp_result %>% unlist %>% mean
(temp_result %>% unlist %>% mean) - 0.1701289
temp_result %>% unlist %>% sd

temp_result %>% unlist %>% hist(xlab = paste0("Bias ", round((temp_result %>% unlist %>% mean) - 0.1701289, 5), ", SD ", temp_result %>% unlist %>% sd %>% round(5)))
abline(v = 0.1701289, col = "red", lty = 2)


