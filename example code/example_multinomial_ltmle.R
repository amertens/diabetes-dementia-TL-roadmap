
library(ltmle)
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))

GenerateData <- function(n, abar) {
  W <- rnorm(n)
  if (is.null(abar)) {
    A <- rexpit(W) + rexpit(W) #0, 1 or 2
  } else {
    A <- rep(abar, n)
  }
  Y <- rexpit(W + A)
  return(data.frame(W, A, Y))
}

psi0 <- mean(GenerateData(n=1e6, abar=2)$Y)
d <- GenerateData(n = 1000, abar = NULL)
dd <- data.frame(W = d$W, Ais0 = as.numeric(d$A == 0), Ais1 = as.numeric(d$A == 1), Y = d$Y)
head(dd)

# Basic method - I think this is probably fine
r <- ltmle(data = dd, Anodes = c("Ais0", "Ais1"), Ynodes = "Y", abar = c(0, 0))
print(summary(r))
print(r$fit$g)

# More complicated method - include the fact that if Ais0 is 1 then Ais1 is deterministically 0.   
# Avoids a regression with very large negative coefficient on Ais0, which might lead to weird results.   
# But in this case and others I've looked at, it gives the same result as above.   
# SuperLearner should also be able to avoid any problems.
det.g.fun <- function (data, current.node, nodes)
{
  if (names(data)[current.node] == "Ais1") {
    is.deterministic <- data[, "Ais0"] == 1
    prob1 <- 0 #if Ais0 is 1 then P(Ais1 = 1) = 0
    return(list(is.deterministic = is.deterministic, prob1 = prob1))
  } else {
    return(NULL)
  }
}

abar <- list(a=c(0,1,0,1,0,1,0), b=c(1,0,1,0,1,0))
abar <- list(a=c(0,1,0,1,0,1,0), b=c(0,0,0,0,0,0))

r2 <- ltmle(data = dd, Anodes = c("A1is0", "A1is1","A2is0", "A2is1"), Ynodes = "Y", abar = abar, deterministic.g.function = det.g.fun, gform = c("Ais0 ~ W", "Ais1 ~ W"))
print(summary(r2))
print(r2$fit$g)