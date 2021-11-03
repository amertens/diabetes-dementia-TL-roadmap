

if (!require(testthatsomemore)) {
  if (!require(devtools)) install.packages('devtools'); require(devtools)
  install_github('robertzk/testthatsomemore')
}

library(ltmle)
library(SuperLearner)

rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
n <- 1000
W1 <- rnorm(n)
W2 <- rnorm(n)
W3 <- rnorm(n)
A <- rexpit(-1 + 2 * W1 - 3*W2 + W3)
Y <- rexpit(W1 + W2 + W3 + A)
data <- data.frame(W1, W2, W3, A, Y)

result1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = F, SL.library = "SL.glmnet")
summary(result1)

SuperLearner_override <- function (Y, X, newX = NULL, family = gaussian(), SL.library,
                                   method = "method.NNLS", id = NULL, verbose = FALSE, control = list(),
                                   cvControl = list(), obsWeights = NULL, env = parent.frame()) {
  stopifnot(identical(SL.library, "SL.glmnet"))
  list(SL.predict = SL.glmnet(Y, X, newX, family, obsWeights, id)$pred)
}

testthatsomemore::package_stub("SuperLearner", "SuperLearner", SuperLearner_override, {
  result2 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = F, SL.library = "SL.glmnet")
}
)
summary(result2)


