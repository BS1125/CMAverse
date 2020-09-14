## -----------------------------------------------------------------------------
set.seed(1)
# data simulation
expit <- function(x) exp(x)/(1+exp(x))
n <- 1000000
C1 <- rnorm(n, mean = 1, sd = 0.1)
C2 <- rbinom(n, 1, 0.6)
A <- rbinom(n, 1, expit(0.2 + 0.5*C1 + 0.1*C2))
M <- rbinom(n, 1, expit(1 + 2*A + 1.5*C1 + 0.8*C2))
Y <- rbinom(n, 1, expit(-5 + 0.8*A - 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2))
yprevalence <- sum(Y)/n
data <- data.frame(A, M, Y, C1, C2)
case_indice <- sample(which(data$Y == 1), 2000, replace = FALSE)
control_indice <- sample(which(data$Y == 0), 2000, replace = FALSE)
data <- data[c(case_indice, control_indice), ]

## -----------------------------------------------------------------------------
library(CMAverse)
cmdag(outcome = "Y", exposure = "A", mediator = "M",
      basec = c("C1", "C2"), postc = NULL, node = TRUE, text_col = "white")

## ----message=F,warning=F,results='hide'---------------------------------------
res_yprevelence <- cmest(data = data, model = "rb", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         mreg = list("logistic"), yreg = "logistic",
                         astar = 0, a = 1, mval = list(1), yref = 1,
                         estimation = "paramfunc", inference = "delta")

## ----message=F,warning=F------------------------------------------------------
summary(res_yprevelence)

## ----message=F,warning=F,results='hide'---------------------------------------
res_yrare <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                   outcome = "Y", exposure = "A",
                   mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                   mreg = list("logistic"), yreg = "logistic",
                   astar = 0, a = 1, mval = list(1), yref = 1,
                   estimation = "paramfunc", inference = "delta")

## ----message=F,warning=F------------------------------------------------------
summary(res_yrare)

