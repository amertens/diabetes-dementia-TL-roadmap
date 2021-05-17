
# -------------------------------------------------------------------------
# R code for replication of results in the main article.
# -------------------------------------------------------------------------
library("simcausal")
options(simcausal.verbose = FALSE)
options(width = 90)
set.seed(1121)

# -------------------------------------------------------------------------
# Functions for plotting MSM survival results (Figure 5)
# -------------------------------------------------------------------------
# Obtain the MSM survival predictions from the \code{data.table} in long format by time (\code{t\_vec}), 
# and by the MSM term (\code{MSMtermName}).
survbyMSMterm <- function(MSMres, t_vec, MSMtermName, 
                          use_actions = NULL, est.msm = NULL) {
  library("data.table")
  if (!is.null(MSMres$S.msm.map)) {
    mapS_exprs <- as.character(MSMres$S.msm.map[, "S_exprs_vec"])
    XMSMterms <- as.character(MSMres$S.msm.map[, "XMSMterms"])
    map_idx <- which(mapS_exprs %in% MSMtermName)
    XMSMtermName <- XMSMterms[map_idx]
    if (!is.null(XMSMtermName) && length(XMSMtermName) > 0) {
      MSMtermName <- XMSMtermName
    }
  }
  print("MSMtermName used")
  print(MSMtermName)
  t_dt <- data.table(t = as.integer(t_vec))
  setkey(t_dt, t)
  get_predict <- function(actname) {
    setkeyv(MSMres$df_long[[actname]], c("t", MSMtermName))
    MSMterm_vals <- as.numeric(MSMres$df_long[[actname]][t_dt, 
                                                         mult = "last"][[MSMtermName]])
    newdata <- data.frame(t = t_vec, MSMterm_vals = MSMterm_vals)
    colnames(newdata) <- c("t", MSMtermName)
    if (!is.null(est.msm)) {
      pred <- predict(est.msm, newdata = newdata, type = "response")
    } else {
      pred <- predict(MSMres$m, newdata = newdata, type = "response")
    }
    return(data.frame(t = t_vec, pred = pred))
  }
  action_names <- names(MSMres$df_long)
  if (!is.null(use_actions)) {
    action_names <- action_names[action_names%in%use_actions]
  }
  surv <- lapply(action_names, function(actname) {
    res <- get_predict(actname)
    if (MSMres$hazard) {
      res$surv <- cumprod(1 - res$pred)
    } else {
      res$surv <- 1 - res$pred
    }
    res$pred <- NULL
    res$action <- actname
    res
  })
  names(surv) <- names(MSMres$df_long)
  surv_melt <- do.call("rbind", surv)
  surv_melt$action <- factor(surv_melt$action, 
                             levels = unique(surv_melt$action), 
                             ordered = TRUE)
  surv_melt
}
plotsurvbyMSMterm <- function(surv_melt_dat) {
  library("ggplot2")
  f_ggplot_surv_wS <- ggplot(data= surv_melt_dat, aes(x=t, y=surv)) +
    geom_line(aes(group = action, color = action), 
              size = .4, linetype = "dashed") +
    theme_bw()
  
}
plotsurvbyMSMterm_facet <- function(surv_melt_dat1, surv_melt_dat2, msm_names = NULL) {
  library("ggplot2")
  if (is.null(msm_names)) {
    msm_names <- c("MSM1", "MSM2")
  }
  surv_melt_dat1$MSM <- msm_names[1]
  surv_melt_dat2$MSM <- msm_names[2]
  
  surv_melt_dat <- rbind(surv_melt_dat1, surv_melt_dat2)
  f_ggplot_surv_wS <- ggplot(data= surv_melt_dat, aes(x = t, y = surv)) +
    geom_line(aes(group = action, color = action), size=.4, 
              linetype = "dashed") +
    theme_bw() +
    facet_wrap( ~ MSM)
}


# #########################################################################
# Replication R code for Section 3
# 3. Simulation study with single time point interventions
# #########################################################################

# -------------------------------------------------------------------------
# 3.1. Specifying parametric structural equation models in simcausal
# -------------------------------------------------------------------------
library("simcausal")
D <- DAG.empty()
D <- D + 
  node("I", distr = "rcategor.int", 
       probs = c(0.1, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1)) + 
  node("W1", distr = "rnorm", 
       mean = ifelse(I == 1, 0, ifelse(I == 2, 3, 10)) + 0.6 * I, sd = 1) + 
  node("W2", distr = "runif", 
       min = 0.025*I, max = 0.7*I) +
  node("W3", distr = "rbern", 
       prob = plogis(-0.5 + 0.7*W1 + 0.3*W2 - 0.2*I)) +
  node("A", distr = "rbern", 
       prob = plogis(+4.2 - 0.5*W1 + 0.1*W2 + 0.2*W3))

D <- D + node("U.Y", distr = "rnorm", mean = 0, sd = 1)

# Example 1:
D <- D + node("Y", distr = "rconst", 
              const = -0.5 + 1.2*A + 0.2*I + 0.1*W1 + 0.3*W2 + 0.2*W3 + U.Y)
Dset1 <- set.DAG(D, latent.v = c("I", "U.Y"))

str(Dset1)

# Example 2:
D <- D + node("Y", distr = "rconst", 
              const = -0.5 + 1.2*A - 0.5*(A * W3) + 0.2*I + 0.2*(W1 + W2 + W3) + U.Y)
Dset2 <- set.DAG(D, latent.v = c("I", "U.Y"))

# Example 3:
D <- D + node("Y", distr = "rconst", 
              const = 
                +1.2*A + 0.05*(W1^2 + W2^3 / 10 + W3) + 0.7*abs(U.Y) + 0.002*I^2 +
                +0.02*abs(1 / sin(U.Y * W2 + A)) * (abs(1/sin(U.Y * W2)) <= 10) +
                +5*(abs(1/sin(U.Y * W2)) > 10))
Dset3 <- set.DAG(D, latent.v = c("I", "U.Y"))

# -------------------------------------------------------------------------
# 3.1. Specifying the structural equation model
# Figure 1: Graphical representation of the structural equation model using a DAG, where the latent nodes \\code{I} and \\code{U.Y} are enclosed in circles.
# -------------------------------------------------------------------------
plotDAG(Dset1, xjitter = 0.2, yjitter = 0.85, 
        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8), 
        vertex_attrs = list(size = 15, label.cex = 1.4))

# -------------------------------------------------------------------------
# 3.2. Simulating observed data (sim)
# -------------------------------------------------------------------------
Odat <- sim(DAG = Dset3, n = 10000, rndseed = 123)
Odat[1, ]

# -------------------------------------------------------------------------
# 3.3. Specifying interventions (+ action)
# -------------------------------------------------------------------------
A1 <- node("A", distr = "rbern", prob = 1)
Dset3 <- Dset3 + action("A1", nodes = A1)
A0 <- node("A", distr = "rbern", prob = 0)
Dset3 <- Dset3 + action("A0", nodes = A0)

names(A(Dset3))
class(A(Dset3)[["A0"]])

# -------------------------------------------------------------------------
# 3.4. Simulating full data (sim)
# -------------------------------------------------------------------------
Xdat1 <- sim(DAG = Dset3, actions = c("A1", "A0"), n = 100000, rndseed = 123)
names(Xdat1)
nrow(Xdat1[["A1"]])
nrow(Xdat1[["A0"]])

Xdat1[["A1"]][1, ]
Xdat1[["A0"]][1, ]

# -------------------------------------------------------------------------
# 3.5. Defining and evaluating various causal target parameters
# -------------------------------------------------------------------------
# Causal parameters defined with set.targetE
Dset3 <- set.targetE(Dset3, outcome = "Y", param = "A1")

eval.target(Dset3, data = Xdat1)$res

eval.target(Dset3, n = 100000, rndseed = 123)$res

Dset3 <- set.targetE(Dset3, outcome = "Y", param = "A1 - A0")
eval.target(Dset3, data = Xdat1)$res

Dset3 <- set.targetE(Dset3, outcome = "Y", param = "A1 / A0")
eval.target(Dset3, data = Xdat1)$res


# Causal parameters defined with set.targetMSM
newA <- node("A", distr = "rbern", prob = d)
Dset3 <- Dset3 + action("A1", nodes = newA, d = 1)
Dset3 <- Dset3 + action("A0", nodes = newA, d = 0)

msm.form <- "Y ~ d"
Dset3 <- set.targetMSM(Dset3, outcome = "Y", form = msm.form, 
                       family = "gaussian")
msm.res <- eval.target(Dset3, n = 100000, rndseed = 123)
msm.res$coef

# -------------------------------------------------------------------------
# 3.6. Defining node distributions
# -------------------------------------------------------------------------
distr.list()

rbern

rnorm_trunc <- function(n, mean, sd, minval = 0) {
  out <- rnorm(n = n, mean = mean, sd = sd)
  minval <- minval[1]
  out[out < minval] <- minval
  out
}

Dmin0 <- DAG.empty()
Dmin0 <- Dmin0 +
  node("W", distr = "rbern", 
       prob = plogis(-0.5)) +
  node("A", distr = "rbern", 
       prob = plogis(-0.5 - 0.3 * W)) +
  node("Y", distr = "rnorm_trunc", 
       mean = -0.1 + 1.2 * A + 0.3 * W, 
       sd = 10)
Dmin0set <- set.DAG(Dmin0)

Dmin0 <- Dmin0 +
  node("Y", distr = "rnorm_trunc", 
       mean = -0.1 + 1.2 * A + 0.3 * W, 
       sd = 10, 
       minval = 10)
Dmin10set <- set.DAG(Dmin0)

Dmin0 <- Dmin0 +
  node("Y", distr = "rnorm_trunc", 
       mean = -0.1 + 1.2 * A + 0.3 * W, 
       sd = 10, 
       minval = ifelse(A == 0, 5, 10))
Dminset <- set.DAG(Dmin0)


# #########################################################################
# Replication R code for Section 4
# 4. Simulation study with multiple time point interventions
# #########################################################################

# -------------------------------------------------------------------------
# 4.1. Specifying the structural equation model
# -------------------------------------------------------------------------
library("simcausal")
D <- DAG.empty()
D <- D +
  node("L2", t = 0, distr = "rbern", 
       prob = 0.05) +
  node("L1", t = 0, distr = "rbern", 
       prob = ifelse(L2[0] == 1, 0.5, 0.1)) +
  node("A1", t = 0, distr = "rbern", 
       prob =
         ifelse(L1[0] == 1 & L2[0] == 0, 0.5, 
                ifelse(L1[0] == 0 & L2[0] == 0, 0.1, 
                       ifelse(L1[0] == 1 & L2[0] == 1, 0.9, 0.5))))

t.end <- 16
D <- D +
  node("Y", t = 1:t.end, distr = "rbern", 
       prob =
         plogis(-6.5 + L1[0] + 4 * L2[t-1] +
                  0.05 * sum(I(L2[0:(t-1)] == rep(0, t)))), 
       EFU = TRUE) +
  node("L2", t = 1:t.end, distr = "rbern", 
       prob =
         ifelse(A1[t-1] == 1, 0.1, 
                ifelse(L2[t-1] == 1, 0.9, min(1, 0.1 + t / 16)))) +
  node("A1", t = 1:t.end, distr = "rbern", 
       prob =
         ifelse(A1[t-1] == 1, 1, 
                ifelse(L1[0] == 1 & L2[t] == 0, 0.3, 
                       ifelse(L1[0] == 0 & L2[t] == 0, 0.1, 
                              ifelse(L1[0] == 1 & L2[t] == 1, 0.7, 0.5)))))
lDAG <- set.DAG(D)

# -------------------------------------------------------------------------
# 4.1. Specifying the structural equation model
# Figure 2: Graphical representation of a portion of the structural equation model using a DAG.
# -------------------------------------------------------------------------
plotDAG(lDAG, tmax = 3, xjitter = 0.3, yjitter = 0.01, 
        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8), 
        vertex_attrs = list(size = 12, label.cex = 0.8))

# -------------------------------------------------------------------------
# 4.2. Simulating observed data (sim)
# -------------------------------------------------------------------------
Odat <- sim(DAG = lDAG, n = 10000, rndseed = 123)
Odat[1, ]

# -------------------------------------------------------------------------
# 4.3. Specifying interventions (+ action)
# -------------------------------------------------------------------------
act_theta <-c(
  node("A1", t = 0, distr = "rbern", 
       prob = ifelse(L2[0] >= theta , 1, 0)), 
  node("A1", t = 1:(t.end), distr = "rbern", 
       prob = ifelse(A1[t-1] == 1, 1, ifelse(L2[t] >= theta, 1, 0))))

Ddyn <- lDAG
Ddyn <- Ddyn + action("A1_th0", nodes = act_theta, theta = 0)
Ddyn <- Ddyn + action("A1_th1", nodes = act_theta, theta = 1)

class(A(Ddyn)[["A1_th0"]])
A(Ddyn)[["A1_th0"]]

A(Ddyn)[["A1_th0"]]$A1_0
Ddyntry <- Ddyn +
  action("A1_th0", nodes = node("A1", t = 0, distr = "rbern", prob = 0))
A(Ddyntry)[["A1_th0"]]$A1_0

A(Ddyntry)[["A1_th0"]]
Ddyntry <- Ddyntry +
  action("A1_th0", nodes = act_theta, theta = 1, newparam = 100)
A(Ddyntry)[["A1_th0"]]

`%+%` <- function(a, b) paste0(a, b)
Dstat <- lDAG
act_A1_tswitch <- node("A1", t = 0:(t.end), distr = "rbern", 
                       prob = ifelse(t >= tswitch, 1, 0))

tswitch_vec <- (0:t.end)
for (tswitch_i in tswitch_vec) {
  abar <- rep(0, length(tswitch_vec))
  abar[which(tswitch_vec >= tswitch_i)] <- 1
  Dstat <- Dstat + action("A1_ts"%+%tswitch_i, 
                          nodes = act_A1_tswitch, 
                          tswitch = tswitch_i, 
                          abar = abar)
}

A(Dstat)[["A1_ts3"]]

# -------------------------------------------------------------------------
# 4.4. Simulating full data (sim)
# -------------------------------------------------------------------------
Xdyn <- sim(Ddyn, actions = c("A1_th0", "A1_th1"), n = 200000, rndseed = 123)

Xdyn[["A1_th0"]][1, ]
Xdyn[["A1_th1"]][1, ]

# -------------------------------------------------------------------------
# 4.5. Converting a dataset from wide to long format (DF.to.long)
# -------------------------------------------------------------------------
Odat.wide <- sim(DAG = lDAG, n = 1000, wide = TRUE, rndseed = 123)
Odat.wide[1:2, 1:16]
Odat.long <- sim(DAG = lDAG, n = 1000, wide = FALSE, rndseed = 123)
Odat.long[1:7, ]

# -------------------------------------------------------------------------
# 4.6. Defining and evaluating various causal target parameters
# Causal parameters defined with set.targetE:
# -------------------------------------------------------------------------
Ddyn <- set.targetE(Ddyn, outcome = "Y", t = 1:16, param = "A1_th1")
surv_th1 <- 1 - eval.target(Ddyn, data = Xdyn)$res
Ddyn <- set.targetE(Ddyn, outcome = "Y", t = 1:16, param = "A1_th0")
surv_th0 <- 1 - eval.target(Ddyn, data = Xdyn)$res

# -------------------------------------------------------------------------
# 4.6. Defining and evaluating various causal target parameters
# Figure 3: Estimates of the true survival curves under the two dynamic interventions
# -------------------------------------------------------------------------
plotSurvEst(surv = list(d_theta1 = surv_th1, d_theta0 = surv_th0), 
            xindx = 1:17, 
            ylab = "Counterfactual survival for each intervention", 
            ylim = c(0.75, 1.0))

Ddyn <- set.targetE(Ddyn, outcome = "Y", t = 12, param = "A1_th1 - A1_th0")
(psi <- round(eval.target(Ddyn, data = Xdyn)$res, 3))

# -------------------------------------------------------------------------
# 4.6. Defining and evaluating various causal target parameters
# Causal parameters defined with set.targetMSM:
# -------------------------------------------------------------------------
# Example 1. Working dynamic MSM for survival probabilities over time.
msm.form <- "Y ~ theta + t + I(theta * t)"
Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = 1:16, form = msm.form, 
                      family = "binomial", hazard = FALSE)
MSMres1 <- eval.target(Ddyn, n = 10000, rndseed = 123)
MSMres1$coef

msm.form <- "Y ~ theta + as.factor(t) + as.factor(t):theta"
Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = 1:16, formula = msm.form, 
                      family = "binomial", hazard = FALSE)
MSMres2 <- eval.target(Ddyn, n = 200000, rndseed = 123)
MSMres2$coef

# -------------------------------------------------------------------------
# 4.6. Defining and evaluating various causal target parameters
# Figure 4: Survival curve estimates evaluated based on working MSM 1 (left) and saturated MSM 2 (right)
# -------------------------------------------------------------------------
par(mfrow = c(1, 2))
surv_MSM1_th1 <- 1 - predict(MSMres1$m, newdata = data.frame(theta = rep(1, 16), t = 1:16), type = "response")
surv_MSM1_th0 <- 1 - predict(MSMres1$m, newdata = data.frame(theta = rep(0, 16), t = 1:16), type = "response")
plotSurvEst(surv = list(MSM_theta1 = surv_MSM1_th1, MSM_theta0 = surv_MSM1_th0), 
            xindx = 1:16, 
            ylab = "MSM Survival, P(T>t)", 
            ylim = c(0.75, 1.0))
surv_MSM2_th1 <- 1 - predict(MSMres2$m, newdata = data.frame(theta = rep(1, 16), t = 1:16), type = "response")
surv_MSM2_th0 <- 1 - predict(MSMres2$m, newdata = data.frame(theta = rep(0, 16), t = 1:16), type = "response")
plotSurvEst(surv = list(MSM_theta1 = surv_MSM2_th1, MSM_theta0 = surv_MSM2_th0), 
            xindx = 1:16, 
            ylab = "MSM Survival, P(T>t)", 
            ylim = c(0.75, 1.0))
abline(v=11, lty = "dashed")
psi0 <- round((surv_MSM2_th0 - surv_MSM2_th1)[12], 3)
text(x=11.5, y=0.87, expression(psi))
text(x=12.3, y=0.871, "=" %+% psi0)

# -------------------------------------------------------------------------
# 4.6. Defining and evaluating various causal target parameters
# Example 2. Working static MSM for discrete-time hazards over time.
# -------------------------------------------------------------------------
Xts <- sim(Dstat, actions = names(A(Dstat)), n = 1000, rndseed = 123)
msm.form_1 <- "Y ~  t + S(mean(abar[0:(t-1)])) + I(t * S(mean(abar[0:(t-1)])))"
Dstat <- set.targetMSM(Dstat, outcome = "Y", t = 1:16, form = msm.form_1, 
                       family = "binomial", hazard = TRUE)
MSMres <- eval.target(Dstat, data = Xts)
MSMres$coef

names(MSMres)
MSMres$S.msm.map
names(MSMres$df_long)
MSMres$df_long[["A1_ts2"]]

# -------------------------------------------------------------------------
# 4.6. Defining and evaluating various causal target parameters
# Figure 5: Survival curve estimates evaluated based on working MSM 2
# -------------------------------------------------------------------------
survMSMh_wS <- survbyMSMterm(MSMres = MSMres, t_vec = 1:16, 
                             MSMtermName = "mean(abar[0:(t - 1)])")
print(plotsurvbyMSMterm(survMSMh_wS))



# #########################################################################
# Replication R code for Section 5
# Replication of Lefbvre et al. simulation study
# #########################################################################

`%+%` <- function(a, b) paste0(a, b)
library("simcausal")
options(simcausal.verbose=FALSE)
options(width=180)

# simulation parameters for Scenario 1:
set.seed(1121)
nsamp <- c(300, 1000, 10000)
trueA <- c(A_0 = -0.294, A_1 = -0.370)
nsims <- 10000

# -------------------------------------------------------------------------
# Lefebvre et al. Tab 2:
# -------------------------------------------------------------------------
covnmT2 <- c(
  c("\\emph{Lefebvre et al.}: Confounder(s) only", rep("", 2)), 
  c("\\emph{Lefebvre et al.}: Confounder(s) &", "risk factors", rep("", 1)))
lefebvreT2 <- data.frame(
  covnm = covnmT2, 
  N = rep(nsamp, 2), 
  A0Bias10 = sprintf("%.3f", c(0.768, 0.265, 0.057, 0.757, 0.283, 0.056)), 
  A0MSE10 = sprintf("%.3f", c(1.761, 0.761, 0.146, 1.642, 0.718, 0.139)), 
  A1Bias10 = sprintf("%.3f", c(0.889, 0.312, 0.086, 0.836, 0.330, 0.081)), 
  A1MSE10 = sprintf("%.3f", c(1.728, 0.723, 0.120, 1.505, 0.638, 0.114)), 
  stringsAsFactors = FALSE)

# -------------------------------------------------------------------------
# Lefebvre et al. Tab 4:
# -------------------------------------------------------------------------
covnmT4 <- c(
  c("\\emph{Lefebvre et al.}: Confounder(s) only", rep("", 2)), 
  c("\\emph{Lefebvre et al.}: Confounder(s) &", "risk factors", ""), 
  c("\\emph{Lefebvre et al.}: Confounder(s) &", "IVs", ""), 
  c("\\emph{Lefebvre et al.}: Confounder(s), ", "IVs & risk factors", ""), 
  c("\\emph{Lefebvre et al.}: Mis-specified", rep("", 2)), 
  c("\\emph{Lefebvre et al.}: Full Model", rep("", 2)))

lefebvreT4 <- data.frame(
  covnm = covnmT4, 
  N = rep(nsamp, 6), 
  A0Bias10 = sprintf("%.3f", c(
    -0.080, -0.371, -0.368, -0.110, -0.330, -0.378, 1.611, 
    0.824, 0.241, 1.600, 0.867, 0.235, 3.146, 2.460, 2.364, 
    1.524, 0.878, 0.240)), 
  A0MSE10 = sprintf("%.3f", c(
    1.170, 0.385, 0.056, 1.092, 0.340, 0.051, 3.538, 2.063, 
    0.684, 3.477, 2.053, 0.676, 3.326, 1.700, 0.832, 3.648, 
    2.099, 0.679)), 
  A1Bias10 = sprintf("%.3f", c(
    0.099, -0.035, -0.203, 0.112, -0.108, -0.207, 2.069, 1.245, 
    0.379, 2.143, 1.170, 0.372, 5.591, 5.258, 4.943, 2.221, 1.185, 
    0.377)), 
  A1MSE10 = sprintf("%.3f", c(
    1.155, 0.331, 0.043, 0.865, 0.245, 0.037, 3.841, 2.188, 0.622, 
    3.598, 2.043, 0.625, 5.494, 3.851, 2.705, 3.907, 2.099, 0.630)), 
  stringsAsFactors = FALSE)

col1name <- "Covariates in $P(A|L)$"
colnames(lefebvreT2)[1] <- colnames(lefebvreT4)[1] <- col1name
col36names <- c(
  "\\specialcell[t]{A(0)\\\\ Bias*10}", 
  "\\specialcell[t]{A(0)\\\\ MSE*10}", 
  "\\specialcell[t]{A(1)\\\\ Bias*10}", 
  "\\specialcell[t]{A(1)\\\\ MSE*10}")
colnames(lefebvreT2)[3:6] <- colnames(lefebvreT4)[3:6] <- col36names

# -------------------------------------------------------------------------
# Lefebvre et al. Scenario 1
# -------------------------------------------------------------------------
# (1) Specify the DAG object with true data generating distribution;
# (2) Set the MSM target parameter;
# (3) Evaluate the true values of the MSM parameters using full data samples of 500K, repeating 50 times;
# (4) Save the true MSM values under trueMSMreps.sc1.Rdata
# -------------------------------------------------------------------------
rbivNorm <- function(n, whichbiv, norms, mu, var1 = 1, var2 = 1, rho = 0.7) {
  whichbiv <- whichbiv[1]; var1 <- var1[1]; var2 <- var2[1]; rho <- rho[1]
  sigma <- matrix(c(var1, rho, rho, var2), nrow = 2)
  Scol <- chol(sigma)[, whichbiv]
  bivX <- (Scol[1] * norms[, 1] + Scol[2] * norms[, 2]) + mu
  bivX
}

Lnames <- c("LO1", "LO2", "LO3", "LC1")
D <- DAG.empty()
for (Lname in Lnames) {
  D <- D +
    node(Lname%+%".norm1", distr = "rnorm", mean = 0, sd = 1) +
    node(Lname%+%".norm2", distr = "rnorm", mean = 0, sd = 1)
}

D <- D +
  node("LO1", t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
       norms = c(LO1.norm1, LO1.norm2), 
       mu = 0) +
  node("LO2", t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
       norms = c(LO2.norm1, LO2.norm2), 
       mu = 0) +
  node("LO3", t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
       norms = c(LO3.norm1, LO3.norm2), 
       mu = 0) +
  node("LC1", t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
       norms = c(LC1.norm1, LC1.norm2), 
       mu = {if (t == 0) {0} else {-0.30 * A[t-1]}}) +
  node("alpha", t = 0:1, distr = "rconst", 
       const = {if(t == 0) {log(0.6)} else {log(1.0)}}) +
  node("A", t = 0:1, distr = "rbern", 
       prob =
         plogis(alpha[t] +
                  log(5)*LC1[t] + {if(t == 0) {0} else {log(5)*A[t-1]}})) +
  node("Y", t = 1, distr = "rnorm", 
       mean = (0.98 * LO1[t] + 0.58 * LO2[t] + 0.33 * LO3[t] +
                 0.98 * LC1[t] - 0.37 * A[t]), 
       sd = 1)
DAGO.sc1 <- set.DAG(D)

defAct <- function (Dact) {
  act.At <- node("A", t = 0:1, distr = "rbern", prob = abar[t])
  Dact <- Dact +
    action("A00", nodes = act.At, abar = c(0, 0)) +
    action("A10", nodes = act.At, abar = c(1, 0)) +
    action("A01", nodes = act.At, abar = c(0, 1)) +
    action("A11", nodes = act.At, abar = c(1, 1))
  return(Dact)
}

Dact.sc1 <- defAct(DAGO.sc1)
msm.form <- "Y ~ S(abar[0]) + S(abar[1])"
Dact.sc1 <- set.targetMSM(Dact.sc1, outcome = "Y", t = 1, 
                          form = msm.form, family = "gaussian")

repstudy2.sc1.truetarget <- function() {
  trueMSMreps.sc1 <- NULL
  reptrue <- 50
  for (i in (1:reptrue)) {
    res.sc1.i <- eval.target(Dact.sc1, n = 500000)$coef
    trueMSMreps.sc1 <- rbind(trueMSMreps.sc1, res.sc1.i)
  }
  return(trueMSMreps.sc1)
}

# NOTE 1: Whenever the file "trueMSMreps.sc1.Rdata" is not found, the following code
# will automatically run repstudy2.sc1.truetarget() to evalute the value of the true causal parameter
# NOTE 2: Evaluating true causal parameter value takes a long time to run
f1name <- "trueMSMreps.sc1.Rdata"
if (file.exists(f1name)) {
  load(f1name)
} else {
  trueMSMreps.sc1 <- repstudy2.sc1.truetarget()
  save(list = "trueMSMreps.sc1", file = f1name)
}
(trueMSM.sc1 <- apply(trueMSMreps.sc1, 2, mean))
# -------------------------------------------------------------------------
# Lefebvre et al. Scenario 3
# -------------------------------------------------------------------------
# (1) Specify the DAG object with true data generating distribution;
# (2) Set the MSM target parameter;
# (3) Evaluate the true values of the MSM parameters using full data samples of 500K, repeating 50 times;
# (4) Save the true MSM values under trueMSMreps.sc3.Rdata
# -------------------------------------------------------------------------
Lnames <- c("LO1", "LO2", "LO3", "LE1", "LE2", "LE3", "LC1", "LC2", "LC3")
D <- DAG.empty()

for (Lname in Lnames) {
  D <- D +
    node(Lname%+%".norm1", distr = "rnorm") +
    node(Lname%+%".norm2", distr = "rnorm")
}

coefAi <- c(-0.10, -0.20, -0.30)
sdLNi <- c(sqrt(1), sqrt(5), sqrt(10))

for (i in (1:3)) {
  D <- D +
    node("LO"%+%i, t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
         mu = 0, 
         params = list(norms = "c(LO"%+%i%+%".norm1, LO"%+%i%+%".norm2)")) +
    node("LE"%+%i, t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
         mu = 0, var1 = 1, var2 = 1, rho = 0.7, 
         params = list(norms = "c(LE"%+%i%+%".norm1, LE"%+%i%+%".norm2)")) +
    node("LC"%+%i, t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
         mu = {if (t == 0) {0} else {.(coefAi[i]) * A[t-1]}}, 
         params = list(norms = "c(LC"%+%i%+%".norm1, LC"%+%i%+%".norm2)")) +
    node("LN"%+%i, t = 0:1, distr = "rnorm", 
         mean = 0, sd = .(sdLNi[i]))
}

D <- D +
  node("alpha", t = 0:1, distr = "rconst", 
       const = {if(t == 0) {log(0.6)} else {log(1.0)}}) +
  node("A", t = 0:1, distr = "rbern", 
       prob = plogis(alpha[t] +
                       log(5) * LC1[t] + log(2) * LC2[t] + log(1.5) * LC3[t] +
                       log(5) * LE1[t] + log(2) * LE2[t] + log(1.5) * LE3[t] +
                       {if (t == 0) {0} else {log(5) * A[t-1]}})) +
  node("Y", t = 1, distr = "rnorm", 
       mean = 0.98 * LO1[t] + 0.58 * LO2[t] + 0.33 * LO3[t] +
         0.98 * LC1[t] + 0.58 * LC2[t] + 0.33 * LC3[t] - 0.39 * A[t], 
       sd = 1)
DAGO.sc3 <- set.DAG(D)

Dact.sc3 <- defAct(DAGO.sc3)
msm.form <- "Y ~ S(abar[0]) + S(abar[1])"
Dact.sc3 <- set.targetMSM(Dact.sc3, outcome = "Y", t = 1, 
                          form = msm.form, family = "gaussian")

repstudy2.sc3.truetarget <- function() {
  trueMSMreps.sc3 <- NULL
  reptrue <- 50
  for (i in (1:reptrue)) {
    res.sc3.i <- eval.target(Dact.sc3, n = 500000)$coef
    trueMSMreps.sc3 <- rbind(trueMSMreps.sc3, res.sc3.i)
  }
  return(trueMSMreps.sc3)
}

# NOTE 1: Whenever the file "trueMSMreps.sc1.Rdata" is not found, the following code
# will automatically run repstudy2.sc1.truetarget() to evalute the value of the true causal parameter
# NOTE 2: Evaluating true causal parameter value takes a long time to run
f2name <- "trueMSMreps.sc3.Rdata"
if (file.exists(f2name)) {
  load(f2name)
} else {
  trueMSMreps.sc3 <- repstudy2.sc3.truetarget()
  save(list = "trueMSMreps.sc3", file = f2name)
}
(trueMSM.sc3 <- apply(trueMSMreps.sc3, 2, mean))
# -------------------------------------------------------------------------
# Runs the simulation for IPTW MSM estimator
# 1) take the observed data DAG (DAGO);
# 2) sample one observed dataset of size nsamp;
# 3) estimate the IPTW MSM using covariates in Lnames;
# 4) repeat the process nsims times;
# 5) evaluate bias/MSE using the true causal MSM parameters in trueA;
# -------------------------------------------------------------------------
runMSMsw <- function(DAGO, Lnames, trueA, nsamp, nsims) {
  Lnames_0 <- Lnames%+%"_0"
  Lnames_1 <- Lnames%+%"_1"
  gforms <- c("A_0 ~ "%+%paste(Lnames_0, collapse = " + "), 
              "A_1 ~ A_0 + "%+%paste(Lnames_1, collapse = " + "))
  res_sw <- NULL
  for (sims in (1:nsims)) {
    datO <- sim(DAGO, n = nsamp)
    glmA_0 <- glm(datO[, c("A_0", Lnames_0)], formula = gforms[1], family = "binomial")
    glmA_1 <- glm(datO[, c("A_1", "A_0", Lnames_0, Lnames_1)], formula = gforms[2], family = "binomial")
    probA0_1 <- predict(glmA_0,  type = "response")
    weight_t0 <- 1 / (probA0_1^(datO$A_0) * (1-probA0_1)^(1-datO$A_0))
    probA1_1 <- predict(glmA_1,  type = "response")
    weight_t1 <- 1 / (probA1_1^(datO$A_1) * (1-probA1_1)^(1-datO$A_1))
    sw1 <- weight_t0*weight_t1
    emp.pA1cA0 <- table(datO$A_1, datO$A_0)/nrow(datO)
    empPA1 <- data.frame(A_0 = c(0, 0, 1, 1), A_1 = c(0, 1, 1, 0))
    empPA1$empPA_1_cA_0 <- apply(empPA1, 1, function(rowA)
      emp.pA1cA0[as.character(rowA["A_1"]), 
                 as.character(rowA["A_0"])])
    empPA1 <- merge(datO[, c("ID", "A_0", "A_1")], empPA1, sort = FALSE)
    empPA1 <- empPA1[order(empPA1$ID), ]
    swts <- empPA1$empPA_1_cA_0 * (weight_t0 * weight_t1)
    datO$swts <- swts
    MSMres_sw <- glm(datO, formula = "Y_1 ~ A_0 + A_1", weights = swts, family = "gaussian")
    res_sw <- rbind(res_sw, coef(MSMres_sw))
  }
  
  meanres <- apply(res_sw, 2, mean)
  Varres <- apply(res_sw, 2, var)
  bias <- c(meanres["A_0"]-trueA["A_0"], meanres["A_1"]-trueA["A_1"])
  MSE <- c(bias^2+Varres[c("A_0", "A_1")])
  bias10 <- sprintf("%.3f", bias*10)
  MSE10 <- sprintf("%.3f", MSE*10)
  resrow <- c(bias10[1], MSE10[1], bias10[2], MSE10[2])
  col36names <- c("\\specialcell[t]{A(0)\\\\ Bias*10}", 
                  "\\specialcell[t]{A(0)\\\\ MSE*10}", 
                  "\\specialcell[t]{A(1)\\\\ Bias*10}", 
                  "\\specialcell[t]{A(1)\\\\ MSE*10}")
  names(resrow) <- col36names
  return(resrow)
}

# -------------------------------------------------------------------------
# Run the simulation for IPTW based on data-generating distribution in Scenario 1
# Save results under "restabSc1_all_1Ksims.Rdata"
# NOTE: The simulation takes several hours to run on a single core
# -------------------------------------------------------------------------
restab <- NULL
runsim <- function(Lnames, DAGO) {
  for (n in nsamp) {
    resSc <- runMSMsw(DAGO = DAGO, Lnames = Lnames, trueA = trueA, 
                      nsamp = n, nsims = nsims)
    restab <- rbind(restab, c(N = n, resSc))
  }
  restab
}

Lnames <- c("LC1")
covnm <- c("Confounder(s) only", rep("", 2))
restab_1 <- cbind(covnm, runsim(Lnames, DAGO.sc1))
Lnames <- c("LC1", "LO1", "LO2", "LO3")
covnm <- c("Confounder(s) &", "risk factors", rep("", 1))
restab_2 <- cbind(covnm, runsim(Lnames, DAGO.sc1))
restab <- rbind(restab_1, restab_2)
col1name <- "Covariates in $P(A|L)$"
colnames(restab)[1] <- col1name
save(list = "restab", file = "restabSc1_all_1Ksims.Rdata")
restab

# -------------------------------------------------------------------------
# Replication R code for Tables 1 & 2 (Lefebvre Scenario 1)
# -------------------------------------------------------------------------
library("Hmisc")
load(file = "restabSc1_all_1Ksims.Rdata")
cat("\n")
latex(restab, file = "", where = "!htpb", caption.loc = "bottom", 
      caption = "Replication of the simulation results
              from \\citet{lefebvre2008} for Scenario 1.", 
      label = "tab2Lefebvre", booktabs = TRUE, rowname = NULL, landscape = FALSE, 
      col.just = c("l", rep("r", 5)), size = "small")

cat("\n")
latex(lefebvreT2, file = "", where = "!htpb", caption.loc = "bottom", 
      caption = "Simulation results for Scenario 1 as reported in
              Table II of \\citet{lefebvre2008}.", 
      label = "origtab2Lefebvre", booktabs = TRUE, rowname = NULL, landscape = FALSE, 
      col.just = c("l", rep("r", 5)), size = "small")


# -------------------------------------------------------------------------
# Run the simulation for IPTW based on the data-generating distribution in Scenario 3
# Save results under "restabSc3_all_1Ksims.Rdata"
# NOTE: The simulation takes several hours to run on a single core
# -------------------------------------------------------------------------
# simulation parameters for Scenario 1:
set.seed(1121)
nsamp <- c(300, 1000, 10000)
trueA <- c(A_0 = -0.316, A_1 = -0.390)
nsims <- 10000
restab <- NULL

runsim <- function(Lnames, DAGO) {
  for (n in nsamp) {
    resSc <- runMSMsw(DAGO = DAGO, Lnames = Lnames, trueA = trueA, 
                      nsamp = n, nsims = nsims)
    restab <- rbind(restab, c(N = n, resSc))
  }
  restab
}

Lnames <- c("LC1", "LC2", "LC3")
covnm <- c("Confounder(s) only", rep("", 2))
restab_1 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
Lnames <- c("LO1", "LO2", "LO3", "LC1", "LC2", "LC3")
covnm <- c("Confounder(s) &", "risk factors", "")
restab_2 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
Lnames <- c("LE1", "LE2", "LE3", "LC1", "LC2", "LC3")
covnm <- c("Confounder(s) &", "IVs", "")
restab_3 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
Lnames <- c("LO1", "LO2", "LO3", "LE1", "LE2", "LE3", "LC1", "LC2", "LC3")
covnm <- c("Confounder(s), ", "IVs & risk factors", "")
restab_4 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
Lnames <- c("LE1", "LE2", "LE3", "LC1")
covnm <- c("Mis-specified", rep("", 2))
restab_5 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
Lnames <- c("LO1", "LO2", "LO3", "LE1", "LE2", "LE3", "LC1", "LC2", "LC3", 
            "LN1", "LN2", "LN3")
covnm <- c("Full Model", rep("", 2))
restab_6 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
restab <- rbind(restab_1, restab_2, restab_3, restab_4, restab_5, restab_6)
col1name <- "Covariates in $P(A|L)$"
colnames(restab)[1] <- col1name
save(list = "restab", file = "restabSc3_all_1Ksims.Rdata")
restab

# -------------------------------------------------------------------------
# Replication R code for Tables 3 & 4 (Lefebvre Scenario 3)
# -------------------------------------------------------------------------
library("Hmisc")
load(file = "restabSc3_all_1Ksims.Rdata")
cat("\n")
latex(restab, file = "", where = "!htpb", caption.loc = "bottom", 
      caption = "Replication of the simulation results from \\citet{lefebvre2008} for Scenario 3.", 
      label = "tab4Lefebvre", booktabs = TRUE, rowname = NULL, landscape = FALSE, 
      col.just = c("l", rep("r", 5)), size = "small")

cat("\n")
latex(lefebvreT4, file = "", where = "!htpb", caption.loc = "bottom", 
      caption = "Simulation results for Scenario 3 as reported in Table IV of \\citet{lefebvre2008}.", 
      label = "origtab4Lefebvre", booktabs = TRUE, rowname = NULL, landscape = FALSE, 
      col.just = c("l", rep("r", 5)), size = "small")

