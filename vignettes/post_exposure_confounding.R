## -----------------------------------------------------------------------------
set.seed(1)
expit <- function(x) exp(x)/(1+exp(x))
n <- 10000
C1 <- rnorm(n, mean = 1, sd = 0.1)
C2 <- rbinom(n, 1, 0.6)
A <- rbinom(n, 1, expit(0.2 + 0.5*C1 + 0.1*C2))
L <- rnorm(n, mean = 1 + A - C1 - 0.5*C2, sd = 0.5)
M <- rbinom(n, 1, expit(1 + 2*A - L + 1.5*C1 + 0.8*C2))
Y <- rbinom(n, 1, expit(-3 - 0.4*A - 1.2*M + 0.5*A*M - 0.5*L + 0.3*C1 - 0.6*C2))
data <- data.frame(A, M, Y, C1, C2, L)

## -----------------------------------------------------------------------------
library(CMAverse)
cmdag(outcome = "Y", exposure = "A", mediator = "M",
      basec = c("C1", "C2"), postc = "L", node = TRUE, text_col = "white")

## ----message=F,warning=F,results='hide'---------------------------------------
res_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                 mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                 ereg = "logistic", yreg = "logistic", mreg = list("logistic"),
                 wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                 astar = 0, a = 1, mval = list(1), 
                 estimation = "imputation", inference = "bootstrap", nboot = 2)

## ----message=F,warning=F------------------------------------------------------
summary(res_msm)

## ----message=F,warning=F,results='hide'---------------------------------------
res_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                      mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                      mreg = list("logistic"), yreg = "logistic", postcreg = list("linear"),
                      astar = 0, a = 1, mval = list(1), 
                      estimation = "imputation", inference = "bootstrap", nboot = 2)

## ----message=F,warning=F------------------------------------------------------
summary(res_gformula)

